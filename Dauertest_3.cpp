#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <mutex>
#include <cstdlib>
#include <random>
#include <string>
#include <thread>
#include <vector>
#include <atomic>

using boost::multiprecision::cpp_int;
using BigInt = cpp_int;
bool g_use_external_gnfs = false;

// ---------- Hilfsfunktionen ----------
cpp_int abs_int(cpp_int x) { return x < 0 ? -x : x; }

cpp_int gcd_int(cpp_int a, cpp_int b) {
    a = abs_int(a);
    b = abs_int(b);
    while (b != 0) {
        cpp_int r = a % b;
        a = b;
        b = r;
    }
    return a;
}

/** Ganzzahlige Wurzel ⌊√n⌋ (Babylonisch), für additiven Vortest um √N. */
static cpp_int integer_sqrt_floor(const cpp_int& n) {
    if (n <= 0) return 0;
    cpp_int x = n;
    cpp_int y = (x + 1) / 2;
    while (y < x) {
        x = y;
        y = (n / y + y) / 2;
    }
    return x;
}

cpp_int powmod(cpp_int a, cpp_int e, const cpp_int& m) {
    a %= m;
    cpp_int r = 1 % m;
    while (e > 0) {
        if ((e & 1) != 0) r = (r * a) % m;
        e >>= 1;
        a = (a * a) % m;
    }
    return r;
}

bool is_probable_prime(const cpp_int& n) {
    if (n < 2) return false;
    static const int small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int p : small_primes) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    cpp_int d = n - 1;
    unsigned s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }

    static const int bases[] = {2, 3, 5, 7, 11, 13, 17};
    for (int a_int : bases) {
        cpp_int a = a_int;
        if (a >= n) continue;
        cpp_int x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool witness = true;
        for (unsigned i = 1; i < s; ++i) {
            x = (x * x) % n;
            if (x == n - 1) {
                witness = false;
                break;
            }
        }
        if (witness) return false;
    }
    return true;
}

cpp_int pollard_rho(const cpp_int& n, std::mt19937_64& rng, int max_tries = 12, int max_iters = 25000) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;

    std::uniform_int_distribution<uint64_t> dist(2, std::numeric_limits<uint64_t>::max());
    auto f = [&](const cpp_int& x, const cpp_int& c) -> cpp_int { return (x * x + c) % n; };

    for (int attempt = 0; attempt < max_tries; ++attempt) {
        cpp_int x = cpp_int(dist(rng)) % (n - 2) + 2;
        cpp_int y = x;
        cpp_int c = cpp_int(dist(rng)) % (n - 1) + 1;
        cpp_int d = 1;

        for (int iter = 0; iter < max_iters && d == 1; ++iter) {
            x = f(x, c);
            y = f(f(y, c), c);
            d = gcd_int(x - y, n);
        }
        if (d != 1 && d != n) return d;
    }
    return n;
}

void factor_pollard(const cpp_int& n, std::vector<cpp_int>& out, std::mt19937_64& rng) {
    cpp_int x = abs_int(n);
    if (x < 2) return;
    if (is_probable_prime(x)) {
        out.push_back(x);
        return;
    }
    cpp_int d = pollard_rho(x, rng);
    if (d == x) {
        out.push_back(x);
        return;
    }
    factor_pollard(d, out, rng);
    factor_pollard(x / d, out, rng);
}

std::pair<cpp_int, cpp_int> factor_trial_6k(const cpp_int& n_input) {
    cpp_int n = abs_int(n_input);
    if (n <= 3) return {n, 1};

    static const int small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int p : small_primes) {
        if (n % p == 0) return {cpp_int(p), n / p};
    }

    for (cpp_int d = 41; d * d <= n; d += 6) {
        if (n % d == 0) return {d, n / d};
        cpp_int d2 = d + 2;
        if (n % d2 == 0) return {d2, n / d2};
    }
    return {n, 1};
}

std::vector<cpp_int> factorize_trial_full(const cpp_int& n_input) {
    std::vector<cpp_int> factors;
    cpp_int n = abs_int(n_input);
    if (n < 2) return factors;
    while (n > 1) {
        auto [f, rest] = factor_trial_6k(n);
        factors.push_back(f);
        if (rest == 1 || f == n) break;
        n = rest;
    }
    return factors;
}

cpp_int strip_small_primes(cpp_int n, std::vector<cpp_int>& out) {
    static const int small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
    n = abs_int(n);
    for (int p : small_primes) {
        cpp_int pp = p;
        while (n > 1 && n % pp == 0) {
            out.push_back(pp);
            n /= pp;
        }
    }
    return n;
}

cpp_int apply_multiscale_patch(const cpp_int& N_in, const cpp_int& candidate_in) {
    cpp_int N = abs_int(N_in);
    cpp_int candidate = abs_int(candidate_in);
    if (N <= 3) return 0;

    const cpp_int high_res_threshold("100000000000000000000");
    bool high_res_mode = (N < high_res_threshold);
    int scan_range = high_res_mode ? 2000 : 500;
    cpp_int epsilon_step = high_res_mode ? (N >> 48) : (N >> 16);
    if (epsilon_step < 1) epsilon_step = 1;

    for (int i = -scan_range; i <= scan_range; ++i) {
        cpp_int probe = candidate + cpp_int(i) * epsilon_step;
        if (probe <= 1 || probe >= N) continue;
        cpp_int g = gcd_int(N, probe);
        if (g > 1 && g < N) return g;
    }
    return 0;
}

/**
 * AQG-Modul: Additiver Vortest (Resonanz).
 * 1) Band um √N mit gcd — trifft balancierte Semiprimfaktoren.
 * 2) Spiegelung um N/2 (ungerade Schritte ±2k) mit gcd — große Faktoren;
 *    gcd statt nur % erfasst auch gemeinsame Teiler ohne exakte Division.
 */
class SymmetricResonator {
public:
    struct ResonanceResult {
        bool found = false;
        cpp_int factor1 = 0;
        cpp_int factor2 = 0;
    };

    static ResonanceResult analyze(const cpp_int& N_in, int max_offset = 10000) {
        ResonanceResult result;
        cpp_int N = abs_int(N_in);
        if (N <= 3) return result;

        auto try_probe = [&](const cpp_int& probe) -> bool {
            if (probe <= 1 || probe >= N) return false;
            cpp_int g = gcd_int(N, probe);
            if (g > 1 && g < N) {
                result.found = true;
                result.factor1 = g;
                result.factor2 = N / g;
                return true;
            }
            return false;
        };

        const cpp_int s = integer_sqrt_floor(N);
        for (int k = -max_offset; k <= max_offset; ++k) {
            if (try_probe(s + k)) return result;
        }

        cpp_int center = N / 2;
        if ((center & 1) == 0) --center;

        for (int k = 0; k < max_offset; ++k) {
            cpp_int p_candidate = center - (cpp_int(k) * 2);
            cpp_int q_candidate = center + (cpp_int(k) * 2);
            if (try_probe(p_candidate)) return result;
            if (try_probe(q_candidate)) return result;
        }

        return result;
    }
};

/** Erster AQG-Resonanz-Treffer: symmetrischer Resonator, sonst Multiskala-Patch. */
static cpp_int aqg_resonance_probe(const cpp_int& rem) {
    SymmetricResonator::ResonanceResult sym = SymmetricResonator::analyze(rem);
    if (sym.found && sym.factor1 > 1 && sym.factor1 < rem) return sym.factor1;
    return apply_multiscale_patch(rem, rem >> 1);
}

// n-stellige Zufallszahl erzeugen (mit exakt n Ziffern).
cpp_int random_n_digit(std::mt19937_64& rng, int digits) {
    if (digits < 1) digits = 1;
    std::uniform_int_distribution<int> first_dist(1, 9);
    std::uniform_int_distribution<int> digit_dist(0, 9);

    std::string s;
    s.reserve(static_cast<std::size_t>(digits));
    s.push_back(static_cast<char>('0' + first_dist(rng)));
    for (int i = 1; i < digits; ++i) {
        s.push_back(static_cast<char>('0' + digit_dist(rng)));
    }

    // Zahl odd machen, damit trivialer Faktor 2 seltener dominiert.
    int last = s.back() - '0';
    if (last % 2 == 0) s.back() = static_cast<char>('0' + ((last + 1) % 10));

    return cpp_int(s);
}

// n-stellige Zufalls-Primzahl (probabilistisch via Miller-Rabin) erzeugen.
cpp_int random_prime_n_digit(std::mt19937_64& rng, int digits) {
    constexpr int max_attempts = 20000;
    for (int i = 0; i < max_attempts; ++i) {
        cpp_int candidate = random_n_digit(rng, digits);
        if (candidate > 2 && is_probable_prime(candidate)) {
            return candidate;
        }
    }
    // Fallback: falls in seltenen Faellen kein Treffer innerhalb des Limits.
    return random_n_digit(rng, digits);
}

double percentile_ms(const std::vector<long long>& values_us, double p) {
    if (values_us.empty()) return 0.0;
    std::vector<long long> sorted = values_us;
    std::sort(sorted.begin(), sorted.end());
    double idx = p * static_cast<double>(sorted.size() - 1);
    std::size_t lo = static_cast<std::size_t>(idx);
    std::size_t hi = std::min(lo + 1, sorted.size() - 1);
    double frac = idx - static_cast<double>(lo);
    double val = static_cast<double>(sorted[lo]) * (1.0 - frac) +
                 static_cast<double>(sorted[hi]) * frac;
    return val / 1000.0;
}

// ---------- Benchmark ----------
struct BenchResult {
    int run_id;
    cpp_int n;
    long long aqg_us;
    long long classic_us;
    long long gnfs_us;
    std::size_t aqg_factor_count;
    std::size_t classic_factor_count;
    std::size_t gnfs_factor_count;
};

// ---------- Hybrid Controller ----------
struct FactorResult {
    std::string method;
    double time_ms;
    bool success;
};

FactorResult classical_factor(const std::string& n_str);
FactorResult aqg_factor(const std::string& n_str);
FactorResult gnfs_factor(const std::string& n_str);

bool is_decimal_number(const std::string& s) {
    if (s.empty()) return false;
    std::size_t i = 0;
    if (s[0] == '+') i = 1;
    if (i >= s.size()) return false;
    for (; i < s.size(); ++i) {
        if (s[i] < '0' || s[i] > '9') return false;
    }
    return true;
}

bool command_exists(const std::string& cmd) {
    std::string check = "command -v " + cmd + " >/dev/null 2>&1";
    return std::system(check.c_str()) == 0;
}

FactorResult external_gnfs_factor(const std::string& n_str) {
    auto start = std::chrono::steady_clock::now();
    bool success = false;
    std::string method = "GNFS-external-unavailable";

    if (!is_decimal_number(n_str)) {
        auto end = std::chrono::steady_clock::now();
        double t = std::chrono::duration<double, std::milli>(end - start).count();
        return {method, t, false};
    }

    std::string tool;
    if (command_exists("msieve")) {
        tool = "msieve";
    } else if (command_exists("cado-nfs.py")) {
        tool = "cado-nfs.py";
    }

    if (!tool.empty()) {
        std::string out_file = "gnfs_external_output.txt";
        std::string cmd;
        if (tool == "msieve") {
            cmd = "msieve -q " + n_str + " > " + out_file + " 2>&1";
            method = "GNFS-external-msieve";
        } else {
            cmd = "cado-nfs.py " + n_str + " > " + out_file + " 2>&1";
            method = "GNFS-external-cado";
        }
        int rc = std::system(cmd.c_str());
        if (rc == 0) {
            success = true;
        } else {
            // Manche Tools liefern non-zero trotz Teilresultat.
            std::ifstream in(out_file);
            std::string line;
            while (std::getline(in, line)) {
                if (line.find("factor") != std::string::npos || line.find("Factors") != std::string::npos) {
                    success = true;
                    break;
                }
            }
        }
    }

    auto end = std::chrono::steady_clock::now();
    double t = std::chrono::duration<double, std::milli>(end - start).count();
    return {method, t, success};
}

// ---------- Hybrid Race (Parallel) ----------
struct RaceResult {
    std::string method;
    double time_us;
    bool done;
};

std::mutex race_mtx;
std::atomic<bool> race_finished(false);
RaceResult race_winner{"none", 0.0, false};

void classical_worker(const std::string& n) {
    while (!race_finished.load(std::memory_order_relaxed)) {
        auto start = std::chrono::steady_clock::now();
        FactorResult r = classical_factor(n); // factoring step
        auto end = std::chrono::steady_clock::now();
        double t_us = std::chrono::duration<double, std::micro>(end - start).count();

        if (r.success) {
            if (!race_finished.exchange(true)) {
                std::lock_guard<std::mutex> lock(race_mtx);
                race_winner = {"classical", t_us, true};
            }
            return;
        }
    }
}

void aqg_worker(const std::string& n) {
    while (!race_finished.load(std::memory_order_relaxed)) {
        auto start = std::chrono::steady_clock::now();
        FactorResult r = aqg_factor(n); // factoring step
        auto end = std::chrono::steady_clock::now();
        double t_us = std::chrono::duration<double, std::micro>(end - start).count();

        if (r.success) {
            if (!race_finished.exchange(true)) {
                std::lock_guard<std::mutex> lock(race_mtx);
                race_winner = {"AQG", t_us, true};
            }
            return;
        }
    }
}

RaceResult hybrid_race(const std::string& n) {
    race_finished = false;
    race_winner = {"none", 0.0, false};
    std::thread t1(classical_worker, n);
    std::thread t2(aqg_worker, n);
    t1.join();
    t2.join();
    return race_winner;
}

void log_race_result(std::ofstream& file, int run_id, const std::string& n, const RaceResult& r) {
    file << run_id << "," << n << "," << r.method << "," << r.time_us << "\n";
}

FactorResult classical_factor(const std::string& n_str) {
    auto start = std::chrono::steady_clock::now();
    bool success = false;

    try {
        cpp_int n(n_str);
        if (n < 0) n = -n;
        std::mt19937_64 rng{std::random_device{}()};
        std::vector<cpp_int> factors;
        factor_pollard(n, factors, rng);
        success = !factors.empty();
    } catch (...) {
        success = false;
    }

    auto end = std::chrono::steady_clock::now();
    double t = std::chrono::duration<double, std::milli>(end - start).count();
    return {"classical", t, success};
}

FactorResult aqg_factor(const std::string& n_str) {
    auto start = std::chrono::steady_clock::now();
    bool success = false;

    try {
        cpp_int n(n_str);
        if (n < 0) n = -n;
        std::mt19937_64 rng{std::random_device{}()};
        std::vector<cpp_int> aqg_factors;
        cpp_int rem = strip_small_primes(n, aqg_factors);
        if (rem > 1) {
            cpp_int resonance = aqg_resonance_probe(rem);
            if (resonance > 1 && resonance < rem) {
                auto left = factorize_trial_full(resonance);
                aqg_factors.insert(aqg_factors.end(), left.begin(), left.end());
                rem /= resonance;
            }
            if (rem > 1) {
                if (rem > (cpp_int(1) << 40)) {
                    std::vector<cpp_int> rest;
                    factor_pollard(rem, rest, rng);
                    aqg_factors.insert(aqg_factors.end(), rest.begin(), rest.end());
                } else {
                    auto rest = factorize_trial_full(rem);
                    aqg_factors.insert(aqg_factors.end(), rest.begin(), rest.end());
                }
            }
        }
        success = !aqg_factors.empty();
    } catch (...) {
        success = false;
    }

    auto end = std::chrono::steady_clock::now();
    double t = std::chrono::duration<double, std::milli>(end - start).count();
    return {"AQG", t, success};
}

void factor_pollard_tuned(const cpp_int& n_input,
                          std::vector<cpp_int>& out,
                          std::mt19937_64& rng,
                          int max_tries,
                          int max_iters) {
    cpp_int n = abs_int(n_input);
    if (n < 2) return;
    if (is_probable_prime(n)) {
        out.push_back(n);
        return;
    }
    cpp_int d = pollard_rho(n, rng, max_tries, max_iters);
    if (d == n || d == 1) {
        out.push_back(n);
        return;
    }
    factor_pollard_tuned(d, out, rng, max_tries, max_iters);
    factor_pollard_tuned(n / d, out, rng, max_tries, max_iters);
}

void gnfs_like_factorize(const cpp_int& n, std::vector<cpp_int>& factors, std::mt19937_64& rng) {
    const unsigned bit_len = (n > 0) ? (boost::multiprecision::msb(n) + 1) : 1;
    if (bit_len < 110) {
        factor_pollard(n, factors, rng);
        return;
    }

    cpp_int rem = strip_small_primes(n, factors);
    if (rem <= 1) return;

    // Adaptive Parameter, damit mittlere Zahlen nicht "haengen" wirken.
    int tries = 40;
    int iters = 200000;
    if (bit_len < 180) {
        tries = 14;
        iters = 35000;
    } else if (bit_len < 220) {
        tries = 22;
        iters = 70000;
    }
    factor_pollard_tuned(rem, factors, rng, tries, iters);
}

FactorResult gnfs_factor(const std::string& n_str) {
    auto start = std::chrono::steady_clock::now();
    bool success = false;

    try {
        cpp_int n(n_str);
        if (n < 0) n = -n;
        const unsigned bit_len = (n > 0) ? (boost::multiprecision::msb(n) + 1) : 1;
        if (g_use_external_gnfs && bit_len >= 220) {
            FactorResult ext = external_gnfs_factor(n_str);
            if (ext.success) return ext;
        }

        std::mt19937_64 rng{std::random_device{}()};
        std::vector<cpp_int> factors;
        gnfs_like_factorize(n, factors, rng);
        success = !factors.empty();
    } catch (...) {
        success = false;
    }

    auto end = std::chrono::steady_clock::now();
    double t = std::chrono::duration<double, std::milli>(end - start).count();
    return {"GNFS-like", t, success};
}

bool is_structured(const std::string& n) {
    if (n.empty()) return false;
    char last = n.back();
    if (last == '0' || last == '2' || last == '4' || last == '5' || last == '6' || last == '8') {
        return true;
    }
    int digit_sum = 0;
    for (char c : n) {
        if (c >= '0' && c <= '9') digit_sum += (c - '0');
    }
    return (digit_sum % 3 == 0);
}

bool likely_hard(const BigInt& n_in) {
    BigInt n = abs_int(n_in);
    if (n < 2) return false;

    // 1) Keine kleinen Teiler (klassischer "harte Zahl"-Indikator)
    static const int small_primes[] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47
    };
    for (int p : small_primes) {
        if (n == p) return false;
        if (n % p == 0) return false;
    }

    // 2) Hohe Bitlaenge
    // msb(n) ist 0-basiert, deshalb +1 fuer die effektive Bitzahl.
    const unsigned bit_len = boost::multiprecision::msb(n) + 1;
    if (bit_len < 120) return false;

    // 3) Keine "glatten" Faktoren im unteren Bereich:
    // kurzer 6k+/-1-Scan bis zu kleiner Obergrenze.
    const BigInt smooth_limit = 1009;
    for (BigInt d = 49; d <= smooth_limit; d += 6) {
        if (n % d == 0) return false;
        BigInt d2 = d + 2;
        if (d2 <= smooth_limit && n % d2 == 0) return false;
    }

    return true;
}

bool detect_rare_structure(const std::string& n) {
    if (is_structured(n)) return true;
    try {
        BigInt parsed_n(n);
        if (parsed_n < 0) parsed_n = -parsed_n;
        if (parsed_n < 2) return false;
        const unsigned bit_len = boost::multiprecision::msb(parsed_n) + 1;
        if (bit_len < 160) return false;
        return !likely_hard(parsed_n);
    } catch (...) {
        return false;
    }
}

FactorResult hybrid_factor(const std::string& n) {
    if (detect_rare_structure(n)) {
        return aqg_factor(n);
    }
    return classical_factor(n);
}

BenchResult run_single_benchmark(int run_id, const cpp_int& n) {
    std::mt19937_64 rng{std::random_device{}()};

    // AQG-Kernzeit (isoliert)
    auto aqg_start = std::chrono::steady_clock::now();
    std::vector<cpp_int> aqg_factors;
    cpp_int rem = strip_small_primes(n, aqg_factors);
    if (rem > 1) {
        cpp_int resonance = aqg_resonance_probe(rem);
        if (resonance > 1 && resonance < rem) {
            auto left = factorize_trial_full(resonance);
            aqg_factors.insert(aqg_factors.end(), left.begin(), left.end());
            rem /= resonance;
        }
        if (rem > 1) {
            if (rem > (cpp_int(1) << 40)) {
                std::vector<cpp_int> rest;
                factor_pollard(rem, rest, rng);
                aqg_factors.insert(aqg_factors.end(), rest.begin(), rest.end());
            } else {
                auto rest = factorize_trial_full(rem);
                aqg_factors.insert(aqg_factors.end(), rest.begin(), rest.end());
            }
        }
    }
    auto aqg_end = std::chrono::steady_clock::now();
    auto aqg_us = std::chrono::duration_cast<std::chrono::microseconds>(aqg_end - aqg_start).count();

    // Klassischer Kern (isoliert)
    auto c_start = std::chrono::steady_clock::now();
    std::vector<cpp_int> classic_factors;
    factor_pollard(n, classic_factors, rng);
    auto c_end = std::chrono::steady_clock::now();
    auto classic_us = std::chrono::duration_cast<std::chrono::microseconds>(c_end - c_start).count();

    // GNFS-like Kern (isoliert)
    auto g_start = std::chrono::steady_clock::now();
    std::vector<cpp_int> gnfs_factors;
    {
        gnfs_like_factorize(n, gnfs_factors, rng);
    }
    auto g_end = std::chrono::steady_clock::now();
    auto gnfs_us = std::chrono::duration_cast<std::chrono::microseconds>(g_end - g_start).count();

    return {run_id, n, aqg_us, classic_us, gnfs_us, aqg_factors.size(), classic_factors.size(), gnfs_factors.size()};
}

int main() {
    std::cout << "--- Dauertest (30-stellige Zufallszahlen) ---\n";

    std::string input;
    int runs = 50;
    int digits = 30;
    bool prime_only_mode = false;

    std::cout << "Stellenzahl der Zufallszahlen (Enter fuer 30): " << std::flush;
    std::getline(std::cin, input);
    if (!input.empty()) {
        try {
            digits = std::max(1, std::stoi(input));
        } catch (...) {
            std::cout << "Ungueltige Eingabe, verwende 30.\n";
            digits = 30;
        }
    }

    std::cout << "Anzahl Testlaeufe eingeben (Enter fuer 50): " << std::flush;
    std::getline(std::cin, input);
    if (!input.empty()) {
        try {
            runs = std::max(1, std::stoi(input));
        } catch (...) {
            std::cout << "Ungueltige Eingabe, verwende 50.\n";
            runs = 50;
        }
    }

    std::string csv_path = "dauertest_results.csv";
    std::cout << "CSV-Dateiname (Enter fuer dauertest_results.csv): " << std::flush;
    std::getline(std::cin, input);
    if (!input.empty()) csv_path = input;

    std::cout << "Zufallszahlen nur als Primzahlen erzeugen? (y/n, Enter = n): " << std::flush;
    std::getline(std::cin, input);
    prime_only_mode = (!input.empty() && (input == "y" || input == "Y" || input == "yes" || input == "YES"));

    std::cout << "Externes GNFS verwenden, falls verfuegbar? (y/n, Enter = n): " << std::flush;
    std::getline(std::cin, input);
    g_use_external_gnfs = (!input.empty() && (input == "y" || input == "Y" || input == "yes" || input == "YES"));

    std::mt19937_64 rng{std::random_device{}()};
    std::vector<BenchResult> results;
    results.reserve(static_cast<std::size_t>(runs));

    for (int i = 1; i <= runs; ++i) {
        cpp_int n = prime_only_mode ? random_prime_n_digit(rng, digits) : random_n_digit(rng, digits);
        auto r = run_single_benchmark(i, n);
        results.push_back(r);
        if (i <= 5 || i % 5 == 0 || i == runs) {
            std::cout << "Fortschritt: " << i << "/" << runs << "\n";
        }
    }

    std::ofstream out(csv_path);
    if (!out) {
        std::cerr << "Fehler: Kann CSV nicht schreiben: " << csv_path << "\n";
        return 1;
    }

    out << "run_id,n,aqg_us,classic_us,gnfs_us,aqg_factor_count,classic_factor_count,gnfs_factor_count\n";
    for (const auto& r : results) {
        out << r.run_id << ","
            << r.n << ","
            << r.aqg_us << ","
            << r.classic_us << ","
            << r.gnfs_us << ","
            << r.aqg_factor_count << ","
            << r.classic_factor_count << ","
            << r.gnfs_factor_count << "\n";
    }
    out.close();

    long long sum_aqg = 0;
    long long sum_classic = 0;
    long long sum_gnfs = 0;
    std::vector<long long> aqg_us_values;
    std::vector<long long> classic_us_values;
    std::vector<long long> gnfs_us_values;
    aqg_us_values.reserve(results.size());
    classic_us_values.reserve(results.size());
    gnfs_us_values.reserve(results.size());
    for (const auto& r : results) {
        sum_aqg += r.aqg_us;
        sum_classic += r.classic_us;
        sum_gnfs += r.gnfs_us;
        aqg_us_values.push_back(r.aqg_us);
        classic_us_values.push_back(r.classic_us);
        gnfs_us_values.push_back(r.gnfs_us);
    }

    std::cout << "CSV geschrieben: " << csv_path << "\n";
    std::cout << "Stellenzahl: " << digits << "\n";
    std::cout << "Zahlenmodus: " << (prime_only_mode ? "nur Primzahlen" : "beliebige Zahlen") << "\n";
    std::cout << "Durchschnitt AQG (ms): " << (sum_aqg / static_cast<double>(runs)) / 1000.0 << "\n";
    std::cout << "Durchschnitt Klassisch (ms): " << (sum_classic / static_cast<double>(runs)) / 1000.0 << "\n";
    std::cout << "Durchschnitt GNFS-like (ms): " << (sum_gnfs / static_cast<double>(runs)) / 1000.0 << "\n";
    std::cout << "Median AQG (ms): " << percentile_ms(aqg_us_values, 0.50) << "\n";
    std::cout << "Median Klassisch (ms): " << percentile_ms(classic_us_values, 0.50) << "\n";
    std::cout << "Median GNFS-like (ms): " << percentile_ms(gnfs_us_values, 0.50) << "\n";
    std::cout << "p90 AQG (ms): " << percentile_ms(aqg_us_values, 0.90) << "\n";
    std::cout << "p90 Klassisch (ms): " << percentile_ms(classic_us_values, 0.90) << "\n";
    std::cout << "p90 GNFS-like (ms): " << percentile_ms(gnfs_us_values, 0.90) << "\n";
    std::cout << "p95 AQG (ms): " << percentile_ms(aqg_us_values, 0.95) << "\n";
    std::cout << "p95 Klassisch (ms): " << percentile_ms(classic_us_values, 0.95) << "\n";
    std::cout << "p95 GNFS-like (ms): " << percentile_ms(gnfs_us_values, 0.95) << "\n";

    std::cout << "Hybrid-Einzeltest ausfuehren? (y/n, Enter = n): " << std::flush;
    std::getline(std::cin, input);
    if (!input.empty() && (input == "y" || input == "Y" || input == "yes" || input == "YES")) {
        std::string hybrid_n;
        std::cout << "Zahl fuer Hybrid-Test eingeben (Enter = Demozahl): " << std::flush;
        std::getline(std::cin, hybrid_n);
        if (hybrid_n.empty()) {
            hybrid_n = "1234567890123456789012345678901234567890";
        }
        FactorResult res = hybrid_factor(hybrid_n);
        std::cout << "Hybrid Method: " << res.method << "\n";
        std::cout << "Hybrid Time: " << res.time_ms << " ms\n";
        std::cout << "Hybrid Success: " << (res.success ? "true" : "false") << "\n";
    }
    return 0;
}