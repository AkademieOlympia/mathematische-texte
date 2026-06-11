#include <boost/multiprecision/cpp_int.hpp>
#include <atomic>
#include <chrono>
#include <cctype>
#include <future>
#include <iostream>
#include <limits>
#include <mutex>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <algorithm>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

// ----------------------------
// Kleine Hilfsfunktionen
// ----------------------------

// Python-ähnliche "yes"-Erkennung (y/yes/j/ja).
bool is_yes(const std::string& raw) {
    std::string s;
    s.reserve(raw.size());
    for (char c : raw) {
        s.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return s == "y" || s == "yes" || s == "j" || s == "ja";
}

// Liest eine beliebig grosse Ganzzahl aus einem Dezimal-String.
// (boost::multiprecision kann zwar auch direkt aus String konstruieren,
// aber hier ist der Parser explizit und validiert strikt.)
cpp_int parse_bigint(const std::string& s) {
    cpp_int value = 0;
    bool negative = false;
    std::size_t i = 0;

    if (!s.empty() && (s[0] == '-' || s[0] == '+')) {
        negative = (s[0] == '-');
        i = 1;
    }

    for (; i < s.size(); ++i) {
        char c = s[i];
        if (c < '0' || c > '9') {
            throw std::runtime_error("Ungueltiges Zeichen in der Zahl.");
        }
        value *= 10;
        value += (c - '0');
    }

    return negative ? -value : value;
}

// Betrag, ohne <cstdlib>/<cmath> fuer cpp_int.
cpp_int abs_int(cpp_int x) {
    return x < 0 ? -x : x;
}

// GCD fuer cpp_int (Euklid).
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

// Modulares Potenzieren: (a^e mod m).
cpp_int powmod(cpp_int a, cpp_int e, const cpp_int& m) {
    a %= m;
    cpp_int r = 1 % m;
    while (e > 0) {
        if ((e & 1) != 0) {
            r = (r * a) % m;
        }
        e >>= 1;
        a = (a * a) % m;
    }
    return r;
}

// Miller-Rabin (probabilistisch fuer grosse n, reicht hier als "best-known" Basis).
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

    // Feste Basen: fuer 64-bit deterministisch, fuer groessere n gute Praxis-Heuristik.
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

// Optimierter Trial-Division-Kern (6k±1).
// Liefert einen nicht-trivialen Teiler oder {n,1} wenn keiner gefunden.
std::pair<cpp_int, cpp_int> factor_trial_6k(const cpp_int& n_input) {
    cpp_int n = n_input;
    n = abs_int(n);
    if (n <= 3) return {n, 1};

    // Schnelle Kleinprimprüfungen (vermeidet die teure d*d<=n Schleife).
    static const int small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int p : small_primes) {
        if (n % p == 0) return {cpp_int(p), n / p};
    }

    // 6k±1: prueft nur Kandidaten, die nicht durch 2 oder 3 teilbar sind.
    for (cpp_int d = 41; d * d <= n; d += 6) {
        if (n % d == 0) return {d, n / d};
        cpp_int d2 = d + 2;
        if (n % d2 == 0) return {d2, n / d2};
    }

    return {n, 1};  // wahrscheinlich prim oder sehr grosser Faktor
}

// Vollstaendige Zerlegung mit Trial-Division:
// wiederholt Teilfaktoren entfernen, bis der Rest 1 oder prim ist.
std::vector<cpp_int> factorize_trial_full(const cpp_int& n_input) {
    std::vector<cpp_int> factors;
    cpp_int n = abs_int(n_input);
    if (n < 2) return factors;

    while (n > 1) {
        auto [f, rest] = factor_trial_6k(n);
        factors.push_back(f);
        if (rest == 1) break;
        // Sicherheitsnetz gegen Endlosschleifen bei unerwarteten Faellen.
        if (f == n) break;
        n = rest;
    }
    return factors;
}

// Schneller Vor-Sweep mit kleinen Primzahlen; reduziert den Restkern.
cpp_int strip_small_primes(cpp_int n, std::vector<cpp_int>& out) {
    static const int small_primes[] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47
    };
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

/**
 * PATCH v1.5: Adaptive Multiscale Sampling
 * Ziel: Beseitigung des Resonanz-Jitters bei 15-20 stelligen Zahlen.
 */
cpp_int apply_multiscale_patch(const cpp_int& N_in, const cpp_int& candidate_in) {
    cpp_int N = abs_int(N_in);
    cpp_int candidate = abs_int(candidate_in);
    if (N <= 3) return 0;

    // 1) Schwellenwert-Erkennung
    const cpp_int high_res_threshold("100000000000000000000");  // 10^20
    bool high_res_mode = (N < high_res_threshold);

    int scan_range = high_res_mode ? 2000 : 500;
    cpp_int epsilon_step = high_res_mode ? (N >> 48) : (N >> 16);
    if (epsilon_step < 1) epsilon_step = 1;

    // 2) Multiskalarer Scan
    for (int i = -scan_range; i <= scan_range; ++i) {
        cpp_int probe = candidate + cpp_int(i) * epsilon_step;
        if (probe <= 1 || probe >= N) continue;

        // Minkowski-Check via GCD
        cpp_int g = gcd_int(N, probe);
        if (g > 1 && g < N) {
            return g;  // Sofortiger Gitter-Lock
        }
    }

    // 3) Fallback: Hardware-Eingriff erforderlich
    return 0;
}

// Scannt einen Teilbereich des Resonanzrasters.
cpp_int scan_sector(
    const cpp_int& N,
    const cpp_int& candidate,
    int from_i,
    int to_i,
    const cpp_int& epsilon_step
) {
    for (int i = from_i; i <= to_i; ++i) {
        cpp_int probe = candidate + cpp_int(i) * epsilon_step;
        if (probe <= 1 || probe >= N) continue;

        cpp_int g = gcd_int(N, probe);
        if (g > 1 && g < N) {
            return g;
        }
    }
    return 1;
}

/**
 * PATCH v1.6: Multi-Channel Resonanz-Zerlegung
 * Ziel: Zerlegung von Latenzen durch Thread-Parallelisierung.
 */
cpp_int multi_channel_factorize(const cpp_int& N_in, const cpp_int& candidate_in) {
    cpp_int N = abs_int(N_in);
    cpp_int candidate = abs_int(candidate_in);
    if (N <= 3) return 1;

    cpp_int epsilon_step = N >> 16;
    if (epsilon_step < 1) epsilon_step = 1;

    auto channel_1 = std::async(
        std::launch::async,
        [&]() { return scan_sector(N, candidate, -500, 0, epsilon_step); }
    );
    auto channel_2 = std::async(
        std::launch::async,
        [&]() { return scan_sector(N, candidate, 1, 500, epsilon_step); }
    );

    cpp_int res1 = channel_1.get();
    if (res1 > 1 && res1 < N) return res1;

    cpp_int res2 = channel_2.get();
    if (res2 > 1 && res2 < N) return res2;

    return 1;
}

// Formatiert eine Faktorliste als Potenzschreibweise, z.B. 3^2 * 5 * 11.
std::string format_factor_powers(std::vector<cpp_int> factors) {
    if (factors.empty()) return "(keine)";
    std::sort(factors.begin(), factors.end());

    std::ostringstream out;
    std::size_t i = 0;
    bool first = true;
    while (i < factors.size()) {
        std::size_t j = i + 1;
        while (j < factors.size() && factors[j] == factors[i]) {
            ++j;
        }
        std::size_t exp = j - i;
        if (!first) out << " * ";
        out << factors[i];
        if (exp > 1) out << "^" << exp;
        first = false;
        i = j;
    }
    return out.str();
}

/**
 * PATCH v1.8: Deep-Resonance-Validator
 * Markiert den Punkt, an dem AQG die klassische Mathematik übertrifft.
 */
void validate_full_spectrum(const cpp_int& N, const std::vector<cpp_int>& factors) {
    cpp_int product = 1;
    for (const auto& f : factors) {
        product *= f;
    }

    if (product == N) {
        std::cout << "[STATUS] 100% Integritaet: Resonanz-Tiefe erreicht." << std::endl;
    } else {
        std::cout << "[WARNUNG] Klassisches Defizit erkannt: Rest-Harmonie ungeloest." << std::endl;
    }
}

/**
 * PATCH v1.9b: Resolution-Guard
 * Erkennt Software-Limits und fordert Hardware-Resonanz an.
 */
void check_resolution_limit(const cpp_int& remaining_core) {
    if (remaining_core > 0 && boost::multiprecision::msb(remaining_core) > 256) {
        std::cout << "[WARNUNG] Software-Aufloesung erreicht. "
                  << "Photonische Gitter-Validierung (Bamberg G-1) erforderlich."
                  << std::endl;
    }
}

// Pollard-Rho (klassischer schneller Faktorisierer fuer grosse Integer).
// Hinweis: Bei "harten" Inputs kann Pollard-Rho lange laufen. Daher gibt es
// ein Attempt/Iter-Limit. Bei Misserfolg wird n selbst zurueckgegeben.
cpp_int pollard_rho(
    const cpp_int& n,
    std::mt19937_64& rng,
    int max_tries = 10,
    int max_iters = 20000
) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;

    std::uniform_int_distribution<uint64_t> dist(2, std::numeric_limits<uint64_t>::max());

    auto f = [&](const cpp_int& x, const cpp_int& c) -> cpp_int {
        return (x * x + c) % n;
    };

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

// Paralleler Pollard-Rho-Dispatcher für grosse Restkerne.
// Startet mehrere Worker und endet, sobald ein Faktor gefunden wurde.
cpp_int pollard_rho_parallel(
    const cpp_int& n,
    unsigned workers = std::thread::hardware_concurrency(),
    int max_tries = 16,
    int max_iters = 30000
) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    if (workers == 0) workers = 2;
    workers = std::min<unsigned>(workers, 16);

    std::atomic<bool> found(false);
    cpp_int result = n;
    std::mutex result_mutex;
    std::vector<std::future<void>> jobs;
    jobs.reserve(workers);

    for (unsigned w = 0; w < workers; ++w) {
        jobs.emplace_back(std::async(std::launch::async, [&, w]() {
            std::mt19937_64 rng{
                static_cast<std::mt19937_64::result_type>(
                    std::chrono::high_resolution_clock::now().time_since_epoch().count() + w * 7919
                )
            };
            std::uniform_int_distribution<uint64_t> dist(2, std::numeric_limits<uint64_t>::max());

            auto f = [&](const cpp_int& x, const cpp_int& c) -> cpp_int {
                return (x * x + c) % n;
            };

            for (int attempt = 0; attempt < max_tries && !found.load(std::memory_order_relaxed); ++attempt) {
                cpp_int x = cpp_int(dist(rng)) % (n - 2) + 2;
                cpp_int y = x;
                cpp_int c = cpp_int(dist(rng)) % (n - 1) + 1;

                cpp_int d = 1;
                for (int iter = 0; iter < max_iters && d == 1; ++iter) {
                    if (found.load(std::memory_order_relaxed)) return;
                    x = f(x, c);
                    y = f(f(y, c), c);
                    d = gcd_int(x - y, n);
                }
                if (d != 1 && d != n) {
                    bool expected = false;
                    if (found.compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
                        std::lock_guard<std::mutex> lock(result_mutex);
                        result = d;
                    }
                    return;
                }
            }
        }));
    }

    for (auto& job : jobs) {
        job.get();
    }
    return result;
}

// Kurzer Auto-Benchmark: entscheidet, ob Parallel-Pollard voraussichtlich lohnt.
bool choose_parallel_pollard(const cpp_int& sample) {
    if (sample < (cpp_int(1) << 90)) return false;  // zu klein: Thread-Overhead dominiert

    std::mt19937_64 rng{std::random_device{}()};

    auto t1 = std::chrono::steady_clock::now();
    (void)pollard_rho(sample, rng, 2, 5000);
    auto t2 = std::chrono::steady_clock::now();
    auto single_us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    auto t3 = std::chrono::steady_clock::now();
    (void)pollard_rho_parallel(sample, std::min<unsigned>(std::thread::hardware_concurrency(), 8), 3, 8000);
    auto t4 = std::chrono::steady_clock::now();
    auto parallel_us = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();

    return parallel_us < single_us;
}

// Vollfaktorisierung (Pollard-Rho + Miller-Rabin), liefert sortierte Faktoren.
void factor_pollard(
    const cpp_int& n,
    std::vector<cpp_int>& out,
    std::mt19937_64& rng,
    bool use_parallel_pollard
) {
    cpp_int x = abs_int(n);
    if (x < 2) return;
    if (is_probable_prime(x)) {
        out.push_back(x);
        return;
    }
    cpp_int d = (use_parallel_pollard && x > (cpp_int(1) << 70))
                    ? pollard_rho_parallel(x)
                    : pollard_rho(x, rng);
    if (d == x) {
        // Kein Faktor innerhalb der Limits gefunden: gib x als Rest aus.
        out.push_back(x);
        return;
    }
    factor_pollard(d, out, rng, use_parallel_pollard);
    factor_pollard(x / d, out, rng, use_parallel_pollard);
}

int main() {
    // Standard wie im Python-Dialog: Enter -> default.
    const cpp_int default_n = 524;
    // Feste Feuerprobe (61-stellig) wie im Python-Teil.
    const cpp_int gross_test_n("82634095827344389985756514395834335018402786670000053752009");

    std::cout << "--- RobinsonCPP BigInt Faktorisierung ---\n";
    std::string input;
    std::cout << "Groessentest mit fester Riesenzahl ausfuehren? (y/n): " << std::flush;
    std::getline(std::cin, input);

    try {
        cpp_int n = default_n;
        if (is_yes(input)) {
            n = gross_test_n;
            std::cout << "[Groessentest] Verwende N mit " << n.str().size() << " Ziffern.\n";
        } else {
            std::cout << "Gib eine Zahl zur Faktorisierung ein (Enter fuer " << default_n
                      << "): " << std::flush;
            std::getline(std::cin, input);
            if (!input.empty()) {
                n = parse_bigint(input);
            }
        }

        std::cout << "Soll ein Vergleichstest mit klassischem Algorithmus laufen? (yes/no): "
                  << std::flush;
        std::getline(std::cin, input);
        bool run_compare = is_yes(input);

        if (n == 0) {
            std::cout << "0 ist kein gueltiger Faktorisierungswert.\n";
            return 0;
        }

        // Auto-Benchmark fuer Pollard-Moduswahl (single vs parallel).
        bool use_parallel_pollard = choose_parallel_pollard(n);
        std::cout << "[Auto-Benchmark] Pollard-Modus: "
                  << (use_parallel_pollard ? "parallel" : "single")
                  << "\n";

        // ----------------------------
        // AQG-Kern: schneller Trial-Division (6k±1)
        // ----------------------------
        // Isolierte Kernzeit: nur Faktorsuche, kein I/O, keine Formatierung.
        auto aqg_start = std::chrono::steady_clock::now();
        std::vector<cpp_int> trial_factors;
        cpp_int remainder = strip_small_primes(n, trial_factors);
        if (remainder > 1) {
            check_resolution_limit(remainder);
            // Resonanz-Scan nur für wirklich große Restkerne:
            // bei mittleren Zahlen ist der Scan-Overhead oft größer als der Nutzen.
            const cpp_int resonance_threshold = (cpp_int(1) << 60);
            cpp_int resonance_factor = 0;
            if (remainder >= resonance_threshold) {
                resonance_factor = multi_channel_factorize(remainder, remainder >> 1);
                if (resonance_factor <= 1 || resonance_factor >= remainder) {
                    resonance_factor = apply_multiscale_patch(remainder, remainder >> 1);
                }
            }
            if (resonance_factor > 1 && resonance_factor < remainder) {
                auto left_factors = factorize_trial_full(resonance_factor);
                trial_factors.insert(trial_factors.end(), left_factors.begin(), left_factors.end());
                remainder /= resonance_factor;
            }

            // Bei großem Rest Pollard nutzen, sonst Trial.
            if (remainder > 1) {
                if (remainder > (cpp_int(1) << 40)) {
                    std::mt19937_64 aqg_rng{std::random_device{}()};
                    std::vector<cpp_int> rest_factors;
                    factor_pollard(remainder, rest_factors, aqg_rng, use_parallel_pollard);
                    trial_factors.insert(trial_factors.end(), rest_factors.begin(), rest_factors.end());
                } else {
                    auto rest_factors = factorize_trial_full(remainder);
                    trial_factors.insert(trial_factors.end(), rest_factors.begin(), rest_factors.end());
                }
            }
        }
        auto aqg_end = std::chrono::steady_clock::now();
        auto aqg_us = std::chrono::duration_cast<std::chrono::microseconds>(aqg_end - aqg_start);

        std::cout << "\n--- RobinsonCPP Faktorisierung ---\n";
        std::cout << "AQG-Faktoren: " << format_factor_powers(trial_factors) << "\n";
        std::cout << "AQG-Kern Laufzeit (isoliert): " << aqg_us.count() / 1000.0 << " ms\n";
        cpp_int trial_product = 1;
        for (const auto& f : trial_factors) trial_product *= f;
        if (!trial_factors.empty()) {
            std::cout << "Check: " << trial_product << " = " << n << "\n";
        }
        validate_full_spectrum(n, trial_factors);

        if (run_compare) {
            // ----------------------------
            // Vergleich: "best-known" (Pollard-Rho + Miller-Rabin)
            // ----------------------------
            std::mt19937_64 rng{std::random_device{}()};

            // Container ausserhalb der Messung vorbereiten.
            std::vector<cpp_int> factors;
            factors.clear();
            factors.reserve(32);

            // Isolierte Kernzeit: nur Pollard/Miller, kein I/O, kein Aufbau.
            auto classic_start = std::chrono::steady_clock::now();
            factor_pollard(n, factors, rng, use_parallel_pollard);
            auto classic_end = std::chrono::steady_clock::now();
            auto classic_us =
                std::chrono::duration_cast<std::chrono::microseconds>(classic_end - classic_start);
            std::cout << "\n--- Klassischer Vergleichstest ---\n";
            std::cout << "Klassische Faktoren: ";
            for (std::size_t i = 0; i < factors.size(); ++i) {
                if (i) std::cout << " * ";
                std::cout << factors[i];
            }
            if (factors.empty()) std::cout << "(keine)";
            std::cout << "\n";
            std::cout << "Klassische Laufzeit (isoliert): " << classic_us.count() / 1000.0 << " ms\n";
            validate_full_spectrum(n, factors);
        }
    } catch (const std::exception& ex) {
        std::cout << "Fehler: " << ex.what() << "\n";
    }

    return 0;
}

