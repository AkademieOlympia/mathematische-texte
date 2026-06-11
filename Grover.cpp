/*
 * Grover II – OpenMP C++ Port
 * Grover-inspirierte Inhaltssuche mit Amplitudenverstärkung.
 * Kompilieren: clang++ -std=c++17 -O3 -fopenmp -o Grover Grover.cpp
 * Mit libomp (macOS): clang++ -std=c++17 -O3 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp -o Grover Grover.cpp
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

static std::string to_lower(const std::string& s) {
    std::string r = s;
    for (char& c : r) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return r;
}

static size_t count_substring(const std::string& haystack, const std::string& needle) {
    size_t cnt = 0;
    size_t pos = 0;
    while ((pos = haystack.find(needle, pos)) != std::string::npos) {
        ++cnt;
        pos += needle.size();
    }
    return cnt;
}

static std::string read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    if (!f) return "";
    std::ostringstream buf;
    buf << f.rdbuf();
    return buf.str();
}

static void find_md_files(const fs::path& root, std::vector<fs::path>& out) {
    if (!fs::exists(root)) return;
    try {
        for (const auto& e : fs::recursive_directory_iterator(root, fs::directory_options::skip_permission_denied)) {
            if (e.is_regular_file()) {
                std::string name = e.path().filename().string();
                if (name.size() >= 3 && name.compare(name.size() - 3, 3, ".md") == 0 &&
                    (name.empty() || name[0] != '.')) {
                    out.push_back(e.path());
                }
            }
        }
    } catch (const fs::filesystem_error&) {}
}

/* Radix-2 FFT (in-place) */
static void fft_impl(std::vector<std::complex<double>>& a, bool inverse) {
    const size_t n = a.size();
    if (n <= 1) return;
    size_t j = 0;
    for (size_t i = 0; i < n; ++i) {
        if (i < j) std::swap(a[i], a[j]);
        size_t k = n >> 1;
        while (k > 0 && j >= k) { j -= k; k >>= 1; }
        j += k;
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * M_PI / len * (inverse ? 1 : -1);
        std::complex<double> wlen(std::cos(ang), std::sin(ang));
        for (size_t i = 0; i < n; i += len) {
            std::complex<double> w(1);
            for (size_t j = 0; j < len / 2; ++j) {
                std::complex<double> u = a[i + j], v = w * a[i + j + len / 2];
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (inverse)
        for (size_t i = 0; i < n; ++i) a[i] /= static_cast<double>(n);
}

static size_t next_pow2(size_t n) {
    size_t p = 1;
    while (p < n) p *= 2;
    return p;
}

static void fft_amplify(std::vector<double>& scores) {
    const size_t n = scores.size();
    if (n == 0) return;
    size_t N = next_pow2(n);
    std::vector<std::complex<double>> spec(N, 0.0);
    for (size_t i = 0; i < n; ++i) spec[i] = scores[i];
    fft_impl(spec, false);
    double mx = 0;
    for (size_t i = 0; i < N; ++i) {
        double ab = std::abs(spec[i]);
        if (ab > mx) mx = ab;
    }
    mx += 1e-12;
    for (size_t i = 0; i < N; ++i)
        spec[i] *= (1.0 + 0.5 * std::abs(spec[i]) / mx);
    fft_impl(spec, true);
    for (size_t i = 0; i < n; ++i) scores[i] = std::abs(spec[i].real());
}

struct GroverSearch {
    std::vector<fs::path> files;
    fs::path vault_root;
    std::vector<std::string> content_cache;

    void load_content() {
        content_cache.resize(files.size());
        const int n = static_cast<int>(files.size());
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            try {
                content_cache[i] = read_file(files[static_cast<size_t>(i)]);
            } catch (...) {
                content_cache[i] = "";
            }
        }
    }

    void relevance_scores(const std::string& term, bool case_sensitive, std::vector<double>& scores) {
        scores.assign(files.size(), 0.0);
        std::string term_check = case_sensitive ? term : to_lower(term);

        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < files.size(); ++i) {
            std::string text = case_sensitive ? content_cache[i] : to_lower(content_cache[i]);
            size_t cnt = count_substring(text, term_check);
            scores[i] = (cnt > 0) ? std::log1p(static_cast<double>(cnt)) : 0.0;
        }

        double mx = *std::max_element(scores.begin(), scores.end());
        if (mx > 0)
            for (double& s : scores) s /= mx;
    }

    void grover_amplify(std::vector<double>& scores, int iterations) {
        const size_t n = scores.size();
        if (n == 0 || iterations < 1) return;

        std::vector<double> amp(n);
        for (size_t i = 0; i < n; ++i)
            amp[i] = (scores[i] > 0) ? -scores[i] : scores[i];

        for (int it = 0; it < iterations - 1; ++it) {
            double mean = std::accumulate(amp.begin(), amp.end(), 0.0) / static_cast<double>(n);
            #pragma omp parallel for simd
            for (size_t i = 0; i < n; ++i)
                amp[i] = 2.0 * mean - amp[i];
        }

        double mean = std::accumulate(amp.begin(), amp.end(), 0.0) / static_cast<double>(n);
        for (size_t i = 0; i < n; ++i)
            scores[i] = std::abs(2.0 * mean - amp[i]);
    }

    void search(const std::string& term, bool case_sensitive, int top_k, int grover_iter, bool use_fft,
                std::vector<std::pair<size_t, double>>& result_indices) {
        std::vector<double> scores;
        relevance_scores(term, case_sensitive, scores);
        std::vector<bool> has_term(scores.size());
        for (size_t i = 0; i < scores.size(); ++i) has_term[i] = (scores[i] > 0);

        if (grover_iter > 0) grover_amplify(scores, grover_iter);
        if (use_fft) fft_amplify(scores);

        std::vector<size_t> idx(scores.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&scores](size_t a, size_t b) { return scores[a] > scores[b]; });

        result_indices.clear();
        for (size_t k = 0; k < idx.size() && static_cast<int>(result_indices.size()) < top_k; ++k) {
            size_t i = idx[k];
            if (!has_term[i]) continue;
            result_indices.emplace_back(i, scores[i]);
        }
    }
};

static void open_file(const fs::path& p, const std::string& obsidian_url) {
    std::string path = p.string();
    std::string escaped_path;
    for (char c : path) {
        if (c == '\\' || c == '"' || c == '$' || c == '`') escaped_path += '\\';
        escaped_path += c;
    }
#ifdef __APPLE__
    if (!obsidian_url.empty()) {
        std::string escaped_url;
        for (char c : obsidian_url) {
            if (c == '\\' || c == '"' || c == '$' || c == '`') escaped_url += '\\';
            escaped_url += c;
        }
        std::string cmd = "open \"" + escaped_url + "\"";
        if (std::system(cmd.c_str()) == 0) return;
    }
    std::string cmd = "open -a Obsidian \"" + escaped_path + "\" 2>/dev/null || open \"" + escaped_path + "\"";
    (void)std::system(cmd.c_str());
#elif defined(_WIN32)
    (void)obsidian_url;
    std::string cmd = "start \"\" \"" + escaped_path + "\"";
    (void)std::system(cmd.c_str());
#else
    (void)obsidian_url;
    std::string cmd = "xdg-open \"" + escaped_path + "\"";
    (void)std::system(cmd.c_str());
#endif
}

/* OSC 8: Klickbarer Hyperlink (iTerm2, Terminal.app, WezTerm, …) */
static void print_hyperlink(const std::string& url, const std::string& label) {
    std::cout << "\033]8;;" << url << "\033\\" << label << "\033]8;;\033\\";
}

static std::string url_encode(const std::string& s) {
    std::ostringstream out;
    for (unsigned char c : s) {
        if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '-' || c == '_' || c == '.' || c == '~')
            out << c;
        else if (c == ' ')
            out << "%20";
        else if (c == '/')
            out << "%2F";
        else if (c == '#')
            out << "%23";
        else {
            char buf[4];
            std::snprintf(buf, sizeof(buf), "%%%02X", c);
            out << buf;
        }
    }
    return out.str();
}

static void print_matches(const std::string& content, const std::string& term, bool case_sensitive,
                          int max_lines, std::vector<std::pair<int, std::string>>& matches) {
    matches.clear();
    std::string term_l = case_sensitive ? term : to_lower(term);
    std::istringstream is(content);
    std::string line;
    int ln = 0;
    while (std::getline(is, line) && static_cast<int>(matches.size()) < max_lines) {
        ++ln;
        std::string check = case_sensitive ? line : to_lower(line);
        if (check.find(term_l) != std::string::npos) {
            std::string snip = line;
            if (snip.size() > 150) snip.resize(150);
            matches.emplace_back(ln, snip);
        }
    }
}

int main(int argc, char** argv) {
    std::string term;
    std::string vault_str;
    int top_k = 10;
    bool case_sensitive = false;
    bool use_fft = true;
    int grover_iter = 2;
    int qubits = 0;  // 0 = Standard, sonst: Iterationen ~ π/4 * 2^(qubits/2)

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-v" && i + 1 < argc) { vault_str = argv[++i]; continue; }
        if (a == "-k" && i + 1 < argc) { top_k = std::atoi(argv[++i]); continue; }
        if (a == "-i") { case_sensitive = true; continue; }
        if (a == "--no-fft") { use_fft = false; continue; }
        if (a == "-g" && i + 1 < argc) { grover_iter = std::atoi(argv[++i]); continue; }
        if ((a == "-q" || a == "--qubits") && i + 1 < argc) { qubits = std::atoi(argv[++i]); continue; }
        if (a[0] != '-') {
            if (!term.empty()) term += ' ';
            term += a;
        }
    }

    if (vault_str.empty()) {
        const char* home = std::getenv("HOME");
        if (home)
            vault_str = std::string(home) + "/Library/Mobile Documents/iCloud~md~obsidian/Documents/Zettel";
        else
            vault_str = ".";
    }

    fs::path vault(vault_str);
    if (!fs::exists(vault)) {
        std::cerr << "Vault nicht gefunden: " << vault << "\n";
        return 1;
    }

    std::vector<fs::path> files;
    find_md_files(vault, files);
    std::cout << "Grover II – OpenMP Suche (" << files.size() << " Notizen)\n" << std::flush;

    if (term.empty()) {
        std::cout << "Verwendung: " << (argc ? argv[0] : "Grover") << " <Suchbegriff> [-v Vault] [-k Top] [-i] [-q Qubits] [--no-fft]\n";
        std::cout << "  -v Pfad   Vault-Verzeichnis\n";
        std::cout << "  -k N      Top-N Treffer (default 10)\n";
        std::cout << "  -i        Case-sensitive\n";
        std::cout << "  -q N      Qubits (setzt Grover-Iterationen ~ π/4·2^(N/2), z.B. -q 10)\n";
        std::cout << "  --no-fft  FFT-Verstärkung deaktivieren\n";
        return 0;
    }

    GroverSearch gs;
    gs.files = std::move(files);
    gs.vault_root = vault;

    if (qubits > 0) {
        grover_iter = std::max(1, static_cast<int>(M_PI / 4.0 * std::pow(2.0, qubits / 2.0)));
    }
    int equiv_qubits = (grover_iter <= 0) ? 0 : std::max(1, static_cast<int>(2.0 * std::log2(grover_iter * 4.0 / M_PI + 0.5)));
    std::cout << "Grover: " << grover_iter << " Iterationen (~" << equiv_qubits << " Qubits)\n" << std::flush;

    auto t0 = std::chrono::high_resolution_clock::now();
    gs.load_content();
    std::vector<std::pair<size_t, double>> results;
    gs.search(term, case_sensitive, top_k, grover_iter, use_fft, results);
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

#ifdef _OPENMP
    std::cout << "OpenMP Threads: " << omp_get_max_threads() << "\n";
#endif
    std::cout << "Suchbegriff: \"" << term << "\" (" << std::fixed << std::setprecision(1) << ms << " ms)\n\n" << std::flush;

    if (results.empty()) {
        std::cout << "Keine Treffer.\n" << std::flush;
        return 0;
    }

    std::cout << results.size() << " Treffer:\n\n" << std::flush;
    std::string term_l = case_sensitive ? term : to_lower(term);
    std::vector<std::pair<int, std::string>> matches;
    std::string vault_name = gs.vault_root.filename().string();
    if (vault_name.empty()) vault_name = gs.vault_root.parent_path().filename().string();
    if (vault_name.empty()) vault_name = "Vault";

    for (size_t r = 0; r < results.size(); ++r) {
        size_t idx = results[r].first;
        double score = results[r].second;
        fs::path fp = gs.files[idx];
        std::string rel;
        try {
            rel = fs::relative(fp, gs.vault_root).string();
        } catch (...) {
            rel = fp.filename().string();
        }

        std::string obsidian_link = "obsidian://open?vault=" + url_encode(vault_name) + "&file=" + url_encode(rel);
        std::string file_path = fp.string();
        std::string file_link = "file://";
        for (char c : file_path) {
            if (c == ' ') file_link += "%20";
            else if (c == '#') file_link += "%23";
            else if ((unsigned char)c >= 128) { char buf[4]; std::snprintf(buf, sizeof(buf), "%%%02X", (unsigned char)c); file_link += buf; }
            else file_link += c;
        }

        print_matches(gs.content_cache[idx], term, case_sensitive, 5, matches);

        std::cout << "[" << (r + 1) << "] --- " << rel << " (Score: " << std::fixed << std::setprecision(4) << score << ") ---\n";
        std::cout << "  Obsidian: ";
        print_hyperlink(obsidian_link, obsidian_link);
        std::cout << "\n  Datei:    ";
        print_hyperlink(file_link, file_link);
        std::cout << "\n";
        if (matches.empty()) {
            std::cout << "  (Zeilenkontext nicht ermittelt)\n";
        } else {
            for (const auto& [ln, snip] : matches)
                std::cout << "  Z" << ln << ": " << (snip.size() > 80 ? snip.substr(0, 80) + "..." : snip) << "\n";
        }
        std::cout << "\n" << std::flush;
    }

    std::cout << "Zahl 1–" << results.size() << " zum Öffnen in Obsidian, Enter = Ende: " << std::flush;
    std::string line;
    if (std::getline(std::cin, line) && !line.empty()) {
        int choice = std::atoi(line.c_str());
        if (choice >= 1 && choice <= static_cast<int>(results.size())) {
            size_t idx = results[static_cast<size_t>(choice - 1)].first;
            fs::path fp = gs.files[idx];
            std::string rel;
            try { rel = fs::relative(fp, gs.vault_root).string(); } catch (...) { rel = fp.filename().string(); }
            std::string obsidian_url = "obsidian://open?vault=" + url_encode(vault_name) + "&file=" + url_encode(rel);
            open_file(fp, obsidian_url);
            std::cout << "→ Geöffnet: " << fp.filename().string() << "\n" << std::flush;
        }
    }

    return 0;
}
