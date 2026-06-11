#include <gmpxx.h>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

bool is_yes(const std::string& raw) {
    std::string s;
    s.reserve(raw.size());
    for (char c : raw) {
        s.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return s == "y" || s == "yes" || s == "j" || s == "ja";
}

mpz_class parse_bigint(const std::string& s) {
    mpz_class value;
    if (s.empty()) return 0;
    if (mpz_set_str(value.get_mpz_t(), s.c_str(), 10) != 0) {
        throw std::runtime_error("Ungueltiges Integer-Format.");
    }
    return value;
}

mpz_class abs_int(mpz_class n) {
    return (n < 0) ? -n : n;
}

mpz_class gcd_int(const mpz_class& a, const mpz_class& b) {
    mpz_class g;
    mpz_gcd(g.get_mpz_t(), abs_int(a).get_mpz_t(), abs_int(b).get_mpz_t());
    return g;
}

static mpz_class integer_sqrt_floor(const mpz_class& n) {
    if (n <= 0) return 0;
    mpz_class r;
    mpz_sqrt(r.get_mpz_t(), n.get_mpz_t());
    return r;
}

mpz_class powmod(mpz_class a, mpz_class e, const mpz_class& m) {
    a %= m;
    mpz_class r = 1 % m;
    while (e > 0) {
        if ((e & 1) != 0) {
            r = (r * a) % m;
        }
        e >>= 1;
        a = (a * a) % m;
    }
    return r;
}

bool is_probable_prime(const mpz_class& n) {
    if (n < 2) return false;
    static const int small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int p : small_primes) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    mpz_class d = n - 1;
    unsigned s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }

    static const int bases[] = {2, 3, 5, 7, 11, 13, 17};
    for (int a_int : bases) {
        mpz_class a = a_int;
        if (a >= n) continue;
        mpz_class x = powmod(a, d, n);
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

std::pair<mpz_class, mpz_class> factor_trial_6k(mpz_class n) {
    if (n < 0) n = -n;
    if (n <= 3) return {n, 1};

    static const int small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int p : small_primes) {
        if (n % p == 0) return {mpz_class(p), n / p};
    }

    for (mpz_class d = 41; d * d <= n; d += 6) {
        if (n % d == 0) return {d, n / d};
        mpz_class d2 = d + 2;
        if (n % d2 == 0) return {d2, n / d2};
    }

    return {n, 1};
}

std::vector<mpz_class> factorize_trial_full(const mpz_class& n_input) {
    std::vector<mpz_class> factors;
    mpz_class n = abs_int(n_input);
    if (n < 2) return factors;

    while (n > 1) {
        auto [f, rest] = factor_trial_6k(n);
        factors.push_back(f);
        if (rest == 1) break;
        if (f == n) break;
        n = rest;
    }
    return factors;
}

std::string format_factor_powers(std::vector<mpz_class> factors) {
    if (factors.empty()) return "(keine)";
    std::sort(factors.begin(), factors.end());

    std::ostringstream out;
    std::size_t i = 0;
    bool first = true;
    while (i < factors.size()) {
        std::size_t j = i + 1;
        while (j < factors.size() && factors[j] == factors[i]) ++j;
        std::size_t exp = j - i;
        if (!first) out << " * ";
        out << factors[i];
        if (exp > 1) out << "^" << exp;
        first = false;
        i = j;
    }
    return out.str();
}

// PATCH v1.5 als GMP-Variante
mpz_class apply_multiscale_patch(const mpz_class& N_in, const mpz_class& candidate_in) {
    mpz_class N = abs_int(N_in);
    mpz_class candidate = abs_int(candidate_in);
    if (N <= 3) return 0;

    const mpz_class high_res_threshold("100000000000000000000");  // 10^20
    bool high_res_mode = (N < high_res_threshold);
    int scan_range = high_res_mode ? 2000 : 500;
    mpz_class epsilon_step = high_res_mode ? (N >> 48) : (N >> 16);
    if (epsilon_step < 1) epsilon_step = 1;

    for (int i = -scan_range; i <= scan_range; ++i) {
        mpz_class probe = candidate + mpz_class(i) * epsilon_step;
        if (probe <= 1 || probe >= N) continue;
        mpz_class g = gcd_int(N, probe);
        if (g > 1 && g < N) return g;
    }
    return 0;
}

/**
 * AQG-Modul: Additiver Vortest — Band um √N (gcd), dann Spiegel um N/2 (gcd).
 */
class SymmetricResonator {
public:
    struct ResonanceResult {
        bool found = false;
        mpz_class factor1 = 0;
        mpz_class factor2 = 0;
    };

    static ResonanceResult analyze(const mpz_class& N_in, int max_offset = 10000) {
        ResonanceResult result;
        mpz_class N = abs_int(N_in);
        if (N <= 3) return result;

        auto try_probe = [&](const mpz_class& probe) -> bool {
            if (probe <= 1 || probe >= N) return false;
            mpz_class g = gcd_int(N, probe);
            if (g > 1 && g < N) {
                result.found = true;
                result.factor1 = g;
                result.factor2 = N / g;
                return true;
            }
            return false;
        };

        const mpz_class s = integer_sqrt_floor(N);
        for (int k = -max_offset; k <= max_offset; ++k) {
            if (try_probe(s + k)) return result;
        }

        mpz_class center = N / 2;
        if ((center % 2) == 0) --center;

        for (int k = 0; k < max_offset; ++k) {
            mpz_class step = mpz_class(k) * 2;
            mpz_class p_candidate = center - step;
            mpz_class q_candidate = center + step;
            if (try_probe(p_candidate)) return result;
            if (try_probe(q_candidate)) return result;
        }

        return result;
    }
};

static mpz_class aqg_resonance_probe(const mpz_class& rem) {
    SymmetricResonator::ResonanceResult sym = SymmetricResonator::analyze(rem);
    if (sym.found && sym.factor1 > 1 && sym.factor1 < rem) return sym.factor1;
    return apply_multiscale_patch(rem, rem >> 1);
}

mpz_class pollard_rho(const mpz_class& n, std::mt19937_64& rng, int max_tries = 10, int max_iters = 20000) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;

    std::uniform_int_distribution<uint64_t> dist(2, std::numeric_limits<uint64_t>::max());
    auto f = [&](const mpz_class& x, const mpz_class& c) -> mpz_class { return (x * x + c) % n; };

    for (int attempt = 0; attempt < max_tries; ++attempt) {
        mpz_class x = mpz_class(static_cast<unsigned long>(dist(rng))) % (n - 2) + 2;
        mpz_class y = x;
        mpz_class c = mpz_class(static_cast<unsigned long>(dist(rng))) % (n - 1) + 1;

        mpz_class d = 1;
        for (int iter = 0; iter < max_iters && d == 1; ++iter) {
            x = f(x, c);
            y = f(f(y, c), c);
            d = gcd_int(x - y, n);
        }
        if (d != 1 && d != n) return d;
    }
    return n;
}

void factor_pollard(const mpz_class& n, std::vector<mpz_class>& out, std::mt19937_64& rng) {
    mpz_class x = abs_int(n);
    if (x < 2) return;
    if (is_probable_prime(x)) {
        out.push_back(x);
        return;
    }
    mpz_class d = pollard_rho(x, rng);
    if (d == x) {
        out.push_back(x);
        return;
    }
    factor_pollard(d, out, rng);
    factor_pollard(x / d, out, rng);
}

int main() {
    const mpz_class default_n = 524;
    const mpz_class gross_test_n("9582734438919857565143953563990834335018402786670000053752009");

    std::cout << "--- RobinsonGMP BigInt Faktorisierung ---\n";
    std::string input;

    std::cout << "Groessentest mit fester Riesenzahl ausfuehren? (y/n): " << std::flush;
    std::getline(std::cin, input);

    try {
        mpz_class n = default_n;
        if (is_yes(input)) {
            n = gross_test_n;
            std::cout << "[Groessentest] Verwende N mit " << n.get_str().size() << " Ziffern.\n";
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

        // AQG-Kern (isoliert): Multiscale-Resonanz + Vollzerlegung
        auto aqg_start = std::chrono::steady_clock::now();
        std::vector<mpz_class> trial_factors;
        mpz_class resonance_factor = aqg_resonance_probe(n);
        if (resonance_factor > 1 && resonance_factor < n) {
            auto left_factors = factorize_trial_full(resonance_factor);
            auto rest_factors = factorize_trial_full(n / resonance_factor);
            trial_factors.insert(trial_factors.end(), left_factors.begin(), left_factors.end());
            trial_factors.insert(trial_factors.end(), rest_factors.begin(), rest_factors.end());
        } else {
            trial_factors = factorize_trial_full(n);
        }
        auto aqg_end = std::chrono::steady_clock::now();
        auto aqg_us = std::chrono::duration_cast<std::chrono::microseconds>(aqg_end - aqg_start);

        std::cout << "\n--- RobinsonGMP Faktorisierung ---\n";
        std::cout << "AQG-Faktoren: " << format_factor_powers(trial_factors) << "\n";
        std::cout << "AQG-Kern Laufzeit (isoliert): " << aqg_us.count() / 1000.0 << " ms\n";
        mpz_class product = 1;
        for (const auto& f : trial_factors) product *= f;
        if (!trial_factors.empty()) {
            std::cout << "Check: " << product << " = " << n << "\n";
        }

        if (run_compare) {
            std::mt19937_64 rng{std::random_device{}()};
            std::vector<mpz_class> factors;
            factors.reserve(32);

            auto classic_start = std::chrono::steady_clock::now();
            factor_pollard(n, factors, rng);
            auto classic_end = std::chrono::steady_clock::now();
            auto classic_us = std::chrono::duration_cast<std::chrono::microseconds>(
                classic_end - classic_start
            );

            std::cout << "\n--- Klassischer Vergleichstest ---\n";
            std::cout << "Klassische Faktoren: ";
            for (std::size_t i = 0; i < factors.size(); ++i) {
                if (i) std::cout << " * ";
                std::cout << factors[i];
            }
            if (factors.empty()) std::cout << "(keine)";
            std::cout << "\n";
            std::cout << "Klassische Laufzeit (isoliert): " << classic_us.count() / 1000.0 << " ms\n";
        }
    } catch (const std::exception& ex) {
        std::cout << "Fehler: " << ex.what() << "\n";
    }

    return 0;
}

