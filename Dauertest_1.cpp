#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <thread>
#include <vector>

using boost::multiprecision::cpp_int;

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
    std::size_t aqg_factor_count;
    std::size_t classic_factor_count;
};

// ---------- Hybrid Controller ----------
struct FactorResult {
    std::string method;
    double time_ms;
    bool success;
};

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
            cpp_int resonance = apply_multiscale_patch(rem, rem >> 1);
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

double quick_probe(std::function<FactorResult(const std::string&)> /*algo*/,
                   const std::string& /*n*/,
                   int probe_ms = 50) {
    auto start = std::chrono::steady_clock::now();
    std::this_thread::sleep_for(std::chrono::milliseconds(probe_ms));
    auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

FactorResult hybrid_factor(const std::string& n) {
    if (is_structured(n)) {
        return aqg_factor(n);
    }
    double t_classical = quick_probe(classical_factor, n);
    double t_aqg = quick_probe(aqg_factor, n);
    if (t_aqg < t_classical) {
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
        cpp_int resonance = apply_multiscale_patch(rem, rem >> 1);
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

    return {run_id, n, aqg_us, classic_us, aqg_factors.size(), classic_factors.size()};
}

int main() {
    std::cout << "--- Dauertest (30-stellige Zufallszahlen) ---\n";

    std::string input;
    int runs = 50;
    int digits = 30;

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

    std::mt19937_64 rng{std::random_device{}()};
    std::vector<BenchResult> results;
    results.reserve(static_cast<std::size_t>(runs));

    for (int i = 1; i <= runs; ++i) {
        cpp_int n = random_n_digit(rng, digits);
        auto r = run_single_benchmark(i, n);
        results.push_back(r);
        if (i % 10 == 0 || i == runs) {
            std::cout << "Fortschritt: " << i << "/" << runs << "\n";
        }
    }

    std::ofstream out(csv_path);
    if (!out) {
        std::cerr << "Fehler: Kann CSV nicht schreiben: " << csv_path << "\n";
        return 1;
    }

    out << "run_id,n,aqg_us,classic_us,aqg_factor_count,classic_factor_count\n";
    for (const auto& r : results) {
        out << r.run_id << ","
            << r.n << ","
            << r.aqg_us << ","
            << r.classic_us << ","
            << r.aqg_factor_count << ","
            << r.classic_factor_count << "\n";
    }
    out.close();

    long long sum_aqg = 0;
    long long sum_classic = 0;
    std::vector<long long> aqg_us_values;
    std::vector<long long> classic_us_values;
    aqg_us_values.reserve(results.size());
    classic_us_values.reserve(results.size());
    for (const auto& r : results) {
        sum_aqg += r.aqg_us;
        sum_classic += r.classic_us;
        aqg_us_values.push_back(r.aqg_us);
        classic_us_values.push_back(r.classic_us);
    }

    std::cout << "CSV geschrieben: " << csv_path << "\n";
    std::cout << "Stellenzahl: " << digits << "\n";
    std::cout << "Durchschnitt AQG (ms): " << (sum_aqg / static_cast<double>(runs)) / 1000.0 << "\n";
    std::cout << "Durchschnitt Klassisch (ms): " << (sum_classic / static_cast<double>(runs)) / 1000.0 << "\n";
    std::cout << "Median AQG (ms): " << percentile_ms(aqg_us_values, 0.50) << "\n";
    std::cout << "Median Klassisch (ms): " << percentile_ms(classic_us_values, 0.50) << "\n";
    std::cout << "p90 AQG (ms): " << percentile_ms(aqg_us_values, 0.90) << "\n";
    std::cout << "p90 Klassisch (ms): " << percentile_ms(classic_us_values, 0.90) << "\n";
    std::cout << "p95 AQG (ms): " << percentile_ms(aqg_us_values, 0.95) << "\n";
    std::cout << "p95 Klassisch (ms): " << percentile_ms(classic_us_values, 0.95) << "\n";

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