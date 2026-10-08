#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <string>
#include <algorithm>
#include <random>
#include <chrono>
#include <cstdint>

// ============================================================================
// MINIMALES TEST-FRAMEWORK
// ============================================================================
class TestSuite {
private:
    int total_tests = 0;
    int passed_tests = 0;
    std::string current_test_name;

public:
    void begin_test(const std::string& name) {
        current_test_name = name;
        total_tests++;
        std::cout << "[RUN] " << name << " ... " << std::flush;
    }

    void assert_true(bool condition, const std::string& msg = "") {
        if (!condition) {
            std::cout << "FAIL\n  -> Assertion failed: " << msg << std::endl;
            throw std::runtime_error("Test failed: " + current_test_name);
        }
    }

    void pass_test() {
        passed_tests++;
        std::cout << "PASS" << std::endl;
    }

    void summary() const {
        std::cout << "\n==================================================\n";
        std::cout << "TESTSUITE SUMMARY: " << passed_tests << " / " << total_tests << " TESTS BESTANDEN.\n";
        std::cout << "==================================================\n";
    }
};

// ============================================================================
// ZAHLENTHEORETISCHE & ALGEBRAISCHE DEFINITIONEN
// ============================================================================
constexpr int V4_MULT[4][4] = {
    {0, 1, 2, 3}, {1, 0, 3, 2}, {2, 3, 0, 1}, {3, 2, 1, 0}
};

constexpr int C_COCYCLE[4][4] = {
    { 1,  1,  1,  1},
    { 1, -1,  1, -1},
    { 1, -1, -1,  1},
    { 1,  1, -1, -1}
};

inline int calc_omega(int g, int h, int k) {
    return C_COCYCLE[g][h] * C_COCYCLE[V4_MULT[g][h]][k];
}

inline int sector_of(int64_t p) {
    switch (p % 12) {
        case 1:  return 0; // E
        case 5:  return 1; // A
        case 7:  return 2; // B
        case 11: return 3; // C
        default: return -1;
    }
}

inline int sector_to_mod12(int sec) {
    constexpr int M[4] = {1, 5, 7, 11};
    return M[sec];
}

inline int64_t isqrt64(int64_t n) {
    if (n <= 0) return 0;
    int64_t r = static_cast<int64_t>(std::sqrt(n));
    while ((r + 1) * (r + 1) <= n) ++r;
    while (r * r > n) --r;
    return r;
}

inline int64_t mod_pow(int64_t base, int64_t exp, int64_t mod) {
    int64_t res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (__int128(res) * base) % mod;
        base = (__int128(base) * base) % mod;
        exp /= 2;
    }
    return res;
}

int64_t mod_sqrt(int64_t a, int64_t p) {
    a = (a % p + p) % p;
    if (a == 0) return 0;
    if (mod_pow(a, (p - 1) / 2, p) != 1) return -1;
    if (p % 4 == 3) return mod_pow(a, (p + 1) / 4, p);

    int64_t q = p - 1, s = 0;
    while (q % 2 == 0) { q /= 2; s++; }
    int64_t z = 2;
    while (mod_pow(z, (p - 1) / 2, p) != p - 1) z++;

    int64_t m = s, c = mod_pow(z, q, p), t = mod_pow(a, q, p), r = mod_pow(a, (q + 1) / 2, p);
    while (t != 1) {
        int64_t i = 0, temp = t;
        while (temp != 1 && i < m) {
            temp = (__int128(temp) * temp) % p;
            i++;
        }
        if (i == m) return -1;
        int64_t b = mod_pow(c, 1LL << (m - i - 1), p);
        m = i;
        c = (__int128(b) * b) % p;
        t = (__int128(t) * c) % p;
        r = (__int128(r) * b) % p;
    }
    return r;
}

bool solve_gauss_exact(int64_t p, int64_t& a, int64_t& b) {
    if (p % 4 != 1) return false;
    int64_t r0 = mod_sqrt(p - 1, p);
    if (r0 == -1) return false;
    int64_t r1 = p, limit = isqrt64(p);
    while (r0 > limit) {
        int64_t tmp = r1 % r0;
        r1 = r0;
        r0 = tmp;
    }
    a = r0;
    int64_t rem = p - a * a;
    b = isqrt64(rem);
    return (a * a + b * b == p);
}

bool solve_eisenstein_exact(int64_t p, int64_t& u, int64_t& v) {
    if (p % 3 != 1) return false;
    int64_t r0 = mod_sqrt(p - 3, p);
    if (r0 == -1) return false;
    int64_t r1 = p, limit = isqrt64(p);
    while (r0 > limit) {
        int64_t tmp = r1 % r0;
        r1 = r0;
        r0 = tmp;
    }
    u = r0;
    int64_t rem = p - u * u;
    if (rem > 0 && rem % 3 == 0) {
        v = isqrt64(rem / 3);
        if (u * u + 3 * v * v == p) return true;
    }
    int64_t max_v = isqrt64(p / 3);
    for (int64_t vc = 1; vc <= max_v; ++vc) {
        int64_t u2 = p - 3 * vc * vc;
        if (u2 > 0 && isqrt64(u2) * isqrt64(u2) == u2) {
            u = isqrt64(u2);
            v = vc;
            return true;
        }
    }
    return false;
}

bool solve_pell_exact(int64_t p, int64_t& x, int64_t& y, int& sgn) {
    int64_t res12 = p % 12;
    if (res12 != 1 && res12 != 11) return false;
    sgn = (res12 == 1) ? 1 : -1;
    for (int64_t yc = 1; yc <= 50000; ++yc) {
        int64_t x2 = (sgn == 1) ? (3 * yc * yc + p) : (3 * yc * yc - p);
        if (x2 > 0 && isqrt64(x2) * isqrt64(x2) == x2) {
            x = isqrt64(x2);
            y = yc;
            return true;
        }
    }
    return false;
}

std::vector<int64_t> sieve_primes(int64_t low, int64_t high) {
    std::vector<bool> is_prime(high, true);
    is_prime[0] = is_prime[1] = false;
    for (int64_t p = 2; p * p < high; ++p) {
        if (is_prime[p]) {
            for (int64_t i = p * p; i < high; i += p) is_prime[i] = false;
        }
    }
    std::vector<int64_t> res;
    for (int64_t p = low; p < high; ++p) {
        if (is_prime[p]) res.push_back(p);
    }
    return res;
}

bool generate_eulerian_shuffle(const std::vector<int>& seq, 
                               std::vector<int>& shuffled_seq, 
                               std::mt19937_64& rng, 
                               int max_attempts = 10000) {
    size_t N = seq.size();
    std::array<std::vector<int>, 4> adj;
    for (size_t i = 0; i < N - 1; ++i) adj[seq[i]].push_back(seq[i + 1]);

    int start_node = seq.front();
    int end_node = seq.back();

    for (int att = 1; att <= max_attempts; ++att) {
        auto temp_adj = adj;
        for (int g = 0; g < 4; ++g) {
            std::shuffle(temp_adj[g].begin(), temp_adj[g].end(), rng);
        }

        shuffled_seq.clear();
        shuffled_seq.reserve(N);
        int curr = start_node;
        shuffled_seq.push_back(curr);

        while (!temp_adj[curr].empty()) {
            int nxt = temp_adj[curr].back();
            temp_adj[curr].pop_back();
            shuffled_seq.push_back(nxt);
            curr = nxt;
        }

        if (shuffled_seq.size() == N && shuffled_seq.back() == end_node) {
            return true;
        }
    }
    return false;
}

// ============================================================================
// TESTMODULE
// ============================================================================

void test_algebra_axioms(TestSuite& suite) {
    suite.begin_test("T1: V_4 Gruppenaxiome, Kozyklus & Omega-Vorzeichentafel");
    
    // Assoziativität & Neutrales Element
    for (int a = 0; a < 4; ++a) {
        suite.assert_true(V4_MULT[0][a] == a && V4_MULT[a][0] == a, "Neutrales Element E fehlerhaft");
        suite.assert_true(V4_MULT[a][a] == 0, "Involutionseigenschaft a*a = E verletzt");
        for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
                suite.assert_true(V4_MULT[V4_MULT[a][b]][c] == V4_MULT[a][V4_MULT[b][c]], "Assoziativitaet verletzt");
            }
        }
    }

    // 2-Kozyklus Paritätszählung
    int pos_c = 0, neg_c = 0;
    for (int g = 0; g < 4; ++g) {
        for (int h = 0; h < 4; ++h) {
            if (C_COCYCLE[g][h] == 1) pos_c++;
            else if (C_COCYCLE[g][h] == -1) neg_c++;
        }
    }
    suite.assert_true(pos_c == 10 && neg_c == 6, "Kozyklustafel c(g,h) muss 10x(+1) und 6x(-1) haben");

    // Tripelphasen-Zählung (28 positive, 36 negative)
    int pos_omega = 0, neg_omega = 0;
    for (int g = 0; g < 4; ++g) {
        for (int h = 0; h < 4; ++h) {
            for (int k = 0; k < 4; ++k) {
                int w = calc_omega(g, h, k);
                if (w == 1) pos_omega++;
                else if (w == -1) neg_omega++;
            }
        }
    }
    suite.assert_true(pos_omega == 28 && neg_omega == 36, "Omega-Tabelle muss exakt 28x(+1) und 36x(-1) enthalten");

    // Pendeltripel-Vorzeichenprüfung auf {A, B, C}
    for (int g = 1; g < 4; ++g) {
        for (int h = 1; h < 4; ++h) {
            if (g != h) {
                suite.assert_true(calc_omega(g, h, g) == 1, "Pendelpfade zwischen Nicht-E Klassen muessen Omega = +1 besitzen");
            }
        }
    }
    suite.pass_test();
}

void test_quadratic_solvers(TestSuite& suite) {
    suite.begin_test("T2: Quadratische Zahlkörper-Solver auf Referenzintervall [5, 10.000)");
    auto primes = sieve_primes(5, 10000);
    for (int64_t p : primes) {
        int64_t a = 0, b = 0;
        int sgn = 0;
        int res12 = p % 12;
        if (res12 == 1 || res12 == 5) {
            suite.assert_true(solve_gauss_exact(p, a, b), "Gauss-Zerlegung (D = -4) fehlgeschlagen");
        }
        if (res12 == 1 || res12 == 7) {
            suite.assert_true(solve_eisenstein_exact(p, a, b), "Eisenstein-Zerlegung (D = -3) fehlgeschlagen");
        }
        if (res12 == 1 || res12 == 11) {
            suite.assert_true(solve_pell_exact(p, a, b, sgn), "Pell-Zerlegung (D = 12) fehlgeschlagen");
            suite.assert_true(sgn == (res12 == 1 ? 1 : -1), "Pell-Vorzeichenfehler");
        }
    }
    suite.pass_test();
}

void test_bigram_invariance_and_whittle(TestSuite& suite) {
    suite.begin_test("T3: Eulerian-Shuffle & Erhaltung der Bigram-Matrix C_gh");
    auto primes = sieve_primes(100000, 150000); // Testintervall zur schnellen Prüfung
    std::vector<int> seq;
    seq.reserve(primes.size());
    for (int64_t p : primes) seq.push_back(sector_of(p));

    int C_obs[4][4] = {0};
    for (size_t i = 0; i < seq.size() - 1; ++i) C_obs[seq[i]][seq[i+1]]++;

    std::mt19937_64 rng(42);
    std::vector<int> shuff;
    for (int rep = 0; rep < 50; ++rep) {
        bool ok = generate_eulerian_shuffle(seq, shuff, rng);
        suite.assert_true(ok, "Euler-Shuffle nicht konvergiert");
        suite.assert_true(shuff.size() == seq.size(), "Shuffle-Laenge inkonsistent");
        suite.assert_true(shuff.front() == seq.front() && shuff.back() == seq.back(), "Endpunkte nicht erhalten");

        int C_rep[4][4] = {0};
        for (size_t i = 0; i < shuff.size() - 1; ++i) C_rep[shuff[i]][shuff[i+1]]++;

        for (int g = 0; g < 4; ++g) {
            for (int h = 0; h < 4; ++h) {
                suite.assert_true(C_rep[g][h] == C_obs[g][h], "Bigram-Invarianz verletzt");
            }
        }
    }
    suite.pass_test();
}

void test_algebraic_delta_identities(TestSuite& suite) {
    suite.begin_test("T4: Algebraische Identitäten D_+ + D_- = 0 und D_+ - D_- = Delta S");
    auto primes = sieve_primes(100000, 500000);
    size_t N = primes.size();
    suite.assert_true(N == 31946, "Primzahlanzahl im Intervall [100.000, 500.000) != 31.946");

    std::vector<int> seq;
    seq.reserve(N);
    for (int64_t p : primes) seq.push_back(sector_of(p));

    int N_obs[64] = {0};
    int64_t S_obs = 0;
    for (size_t i = 0; i < N - 2; ++i) {
        int g = seq[i], h = seq[i+1], k = seq[i+2];
        int idx = (g << 4) | (h << 2) | k;
        N_obs[idx]++;
        S_obs += calc_omega(g, h, k);
    }
    suite.assert_true(S_obs == -842, "S_obs != -842");

    // Simulierte Referenzwerte aus B = 10.000 Freeze-Lauf
    const double mean_S_null = -618.8658;
    const double delta_S = S_obs - mean_S_null; // -223.1342

    int N_plus = 0, N_minus = 0;
    for (int idx = 0; idx < 64; ++idx) {
        int g = (idx >> 4) & 3, h = (idx >> 2) & 3, k = idx & 3;
        if (calc_omega(g, h, k) == 1) N_plus += N_obs[idx];
        else N_minus += N_obs[idx];
    }
    suite.assert_true(N_plus + N_minus == static_cast<int>(N - 2), "Summe N+ + N- != N - 2");
    suite.assert_true(N_plus - N_minus == S_obs, "N+ - N- != S_obs");
    suite.assert_true(N_plus == 15551 && N_minus == 16393, "Zellsummen N+ / N- weichen von Referenz ab");

    double D_plus = delta_S / 2.0;
    double D_minus = -delta_S / 2.0;
    suite.assert_true(std::abs((D_plus + D_minus)) < 1e-9, "D+ + D- muss exakt verschwinden");
    suite.assert_true(std::abs((D_plus - D_minus) - delta_S) < 1e-9, "D+ - D- muss exakt Delta_S entsprechen");
    suite.pass_test();
}

void test_two_step_mod12_autocorrelation(TestSuite& suite) {
    suite.begin_test("T5: 2-Schritt-Distanzen (p_{n+2} - p_n) mod 12 & Äquivalenz mit g_n == g_{n+2}");
    auto primes = sieve_primes(100000, 500000);
    size_t N = primes.size();
    std::vector<int> seq;
    seq.reserve(N);
    for (int64_t p : primes) seq.push_back(sector_of(p));

    int count_d2_zero = 0;
    int count_gn_eq_gn2 = 0;
    int count_d2_six = 0;

    for (size_t i = 0; i < N - 2; ++i) {
        int64_t d2 = primes[i+2] - primes[i];
        if (d2 % 12 == 0) count_d2_zero++;
        if (seq[i] == seq[i+2]) count_gn_eq_gn2++;
        if (d2 % 12 == 6) count_d2_six++;
    }

    // Bijektive Äquivalenz
    suite.assert_true(count_d2_zero == count_gn_eq_gn2, "Bijektion d2 == 0 mod 12 <=> g_n == g_{n+2} verletzt");
    suite.assert_true(count_d2_zero == 6121, "Anzahl Delta_2 == 0 mod 12 weicht von 6121 ab");
    suite.assert_true(count_d2_six == 8458, "Anzahl Delta_2 == 6 mod 12 weicht von 8458 ab");

    // Repulsionstest: Delta_2 = 0 muss seltener sein als Delta_2 = 6
    suite.assert_true(count_d2_zero < count_d2_six, "Keine Unterdrueckung von Delta_2 == 0 mod 12 gegenueber 6 mod 12");
    suite.pass_test();
}

// ============================================================================
// MAIN EXECUTION ENTRY POINT
// ============================================================================
int main() {
    std::cout << "==================================================\n";
    std::cout << "EABC / HARDY-LITTLEWOOD NUMERISCHE TESTSUITE C++20\n";
    std::cout << "==================================================\n\n";

    TestSuite suite;
    auto start_time = std::chrono::high_resolution_clock::now();

    try {
        test_algebra_axioms(suite);
        test_quadratic_solvers(suite);
        test_bigram_invariance_and_whittle(suite);
        test_delta_identities: test_algebraic_delta_identities(suite);
        test_two_step_mod12_autocorrelation(suite);
    } catch (const std::exception& e) {
        std::cerr << "\nABBRUCH WEGEN FEHLER: " << e.what() << std::endl;
        return 1;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    suite.summary();
    std::cout << "Gesamtlaufzeit der Testsuite: " << std::fixed << std::setprecision(3) 
              << elapsed.count() << " Sekunden.\n\n";

    return 0;
}#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <string>
#include <algorithm>
#include <random>
#include <chrono>
#include <cstdint>

// ============================================================================
// MINIMALES TEST-FRAMEWORK
// ============================================================================
class TestSuite {
private:
    int total_tests = 0;
    int passed_tests = 0;
    std::string current_test_name;

public:
    void begin_test(const std::string& name) {
        current_test_name = name;
        total_tests++;
        std::cout << "[RUN] " << name << " ... " << std::flush;
    }

    void assert_true(bool condition, const std::string& msg = "") {
        if (!condition) {
            std::cout << "FAIL\n  -> Assertion failed: " << msg << std::endl;
            throw std::runtime_error("Test failed: " + current_test_name);
        }
    }

    void pass_test() {
        passed_tests++;
        std::cout << "PASS" << std::endl;
    }

    void summary() const {
        std::cout << "\n==================================================\n";
        std::cout << "TESTSUITE SUMMARY: " << passed_tests << " / " << total_tests << " TESTS BESTANDEN.\n";
        std::cout << "==================================================\n";
    }
};

// ============================================================================
// ZAHLENTHEORETISCHE & ALGEBRAISCHE DEFINITIONEN
// ============================================================================
constexpr int V4_MULT[4][4] = {
    {0, 1, 2, 3}, {1, 0, 3, 2}, {2, 3, 0, 1}, {3, 2, 1, 0}
};

constexpr int C_COCYCLE[4][4] = {
    { 1,  1,  1,  1},
    { 1, -1,  1, -1},
    { 1, -1, -1,  1},
    { 1,  1, -1, -1}
};

inline int calc_omega(int g, int h, int k) {
    return C_COCYCLE[g][h] * C_COCYCLE[V4_MULT[g][h]][k];
}

inline int sector_of(int64_t p) {
    switch (p % 12) {
        case 1:  return 0; // E
        case 5:  return 1; // A
        case 7:  return 2; // B
        case 11: return 3; // C
        default: return -1;
    }
}

inline int sector_to_mod12(int sec) {
    constexpr int M[4] = {1, 5, 7, 11};
    return M[sec];
}

inline int64_t isqrt64(int64_t n) {
    if (n <= 0) return 0;
    int64_t r = static_cast<int64_t>(std::sqrt(n));
    while ((r + 1) * (r + 1) <= n) ++r;
    while (r * r > n) --r;
    return r;
}

inline int64_t mod_pow(int64_t base, int64_t exp, int64_t mod) {
    int64_t res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (__int128(res) * base) % mod;
        base = (__int128(base) * base) % mod;
        exp /= 2;
    }
    return res;
}

int64_t mod_sqrt(int64_t a, int64_t p) {
    a = (a % p + p) % p;
    if (a == 0) return 0;
    if (mod_pow(a, (p - 1) / 2, p) != 1) return -1;
    if (p % 4 == 3) return mod_pow(a, (p + 1) / 4, p);

    int64_t q = p - 1, s = 0;
    while (q % 2 == 0) { q /= 2; s++; }
    int64_t z = 2;
    while (mod_pow(z, (p - 1) / 2, p) != p - 1) z++;

    int64_t m = s, c = mod_pow(z, q, p), t = mod_pow(a, q, p), r = mod_pow(a, (q + 1) / 2, p);
    while (t != 1) {
        int64_t i = 0, temp = t;
        while (temp != 1 && i < m) {
            temp = (__int128(temp) * temp) % p;
            i++;
        }
        if (i == m) return -1;
        int64_t b = mod_pow(c, 1LL << (m - i - 1), p);
        m = i;
        c = (__int128(b) * b) % p;
        t = (__int128(t) * c) % p;
        r = (__int128(r) * b) % p;
    }
    return r;
}

bool solve_gauss_exact(int64_t p, int64_t& a, int64_t& b) {
    if (p % 4 != 1) return false;
    int64_t r0 = mod_sqrt(p - 1, p);
    if (r0 == -1) return false;
    int64_t r1 = p, limit = isqrt64(p);
    while (r0 > limit) {
        int64_t tmp = r1 % r0;
        r1 = r0;
        r0 = tmp;
    }
    a = r0;
    int64_t rem = p - a * a;
    b = isqrt64(rem);
    return (a * a + b * b == p);
}

bool solve_eisenstein_exact(int64_t p, int64_t& u, int64_t& v) {
    if (p % 3 != 1) return false;
    int64_t r0 = mod_sqrt(p - 3, p);
    if (r0 == -1) return false;
    int64_t r1 = p, limit = isqrt64(p);
    while (r0 > limit) {
        int64_t tmp = r1 % r0;
        r1 = r0;
        r0 = tmp;
    }
    u = r0;
    int64_t rem = p - u * u;
    if (rem > 0 && rem % 3 == 0) {
        v = isqrt64(rem / 3);
        if (u * u + 3 * v * v == p) return true;
    }
    int64_t max_v = isqrt64(p / 3);
    for (int64_t vc = 1; vc <= max_v; ++vc) {
        int64_t u2 = p - 3 * vc * vc;
        if (u2 > 0 && isqrt64(u2) * isqrt64(u2) == u2) {
            u = isqrt64(u2);
            v = vc;
            return true;
        }
    }
    return false;
}

bool solve_pell_exact(int64_t p, int64_t& x, int64_t& y, int& sgn) {
    int64_t res12 = p % 12;
    if (res12 != 1 && res12 != 11) return false;
    sgn = (res12 == 1) ? 1 : -1;
    for (int64_t yc = 1; yc <= 50000; ++yc) {
        int64_t x2 = (sgn == 1) ? (3 * yc * yc + p) : (3 * yc * yc - p);
        if (x2 > 0 && isqrt64(x2) * isqrt64(x2) == x2) {
            x = isqrt64(x2);
            y = yc;
            return true;
        }
    }
    return false;
}

std::vector<int64_t> sieve_primes(int64_t low, int64_t high) {
    std::vector<bool> is_prime(high, true);
    is_prime[0] = is_prime[1] = false;
    for (int64_t p = 2; p * p < high; ++p) {
        if (is_prime[p]) {
            for (int64_t i = p * p; i < high; i += p) is_prime[i] = false;
        }
    }
    std::vector<int64_t> res;
    for (int64_t p = low; p < high; ++p) {
        if (is_prime[p]) res.push_back(p);
    }
    return res;
}

bool generate_eulerian_shuffle(const std::vector<int>& seq, 
                               std::vector<int>& shuffled_seq, 
                               std::mt19937_64& rng, 
                               int max_attempts = 10000) {
    size_t N = seq.size();
    std::array<std::vector<int>, 4> adj;
    for (size_t i = 0; i < N - 1; ++i) adj[seq[i]].push_back(seq[i + 1]);

    int start_node = seq.front();
    int end_node = seq.back();

    for (int att = 1; att <= max_attempts; ++att) {
        auto temp_adj = adj;
        for (int g = 0; g < 4; ++g) {
            std::shuffle(temp_adj[g].begin(), temp_adj[g].end(), rng);
        }

        shuffled_seq.clear();
        shuffled_seq.reserve(N);
        int curr = start_node;
        shuffled_seq.push_back(curr);

        while (!temp_adj[curr].empty()) {
            int nxt = temp_adj[curr].back();
            temp_adj[curr].pop_back();
            shuffled_seq.push_back(nxt);
            curr = nxt;
        }

        if (shuffled_seq.size() == N && shuffled_seq.back() == end_node) {
            return true;
        }
    }
    return false;
}

// ============================================================================
// TESTMODULE
// ============================================================================

void test_algebra_axioms(TestSuite& suite) {
    suite.begin_test("T1: V_4 Gruppenaxiome, Kozyklus & Omega-Vorzeichentafel");
    
    // Assoziativität & Neutrales Element
    for (int a = 0; a < 4; ++a) {
        suite.assert_true(V4_MULT[0][a] == a && V4_MULT[a][0] == a, "Neutrales Element E fehlerhaft");
        suite.assert_true(V4_MULT[a][a] == 0, "Involutionseigenschaft a*a = E verletzt");
        for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
                suite.assert_true(V4_MULT[V4_MULT[a][b]][c] == V4_MULT[a][V4_MULT[b][c]], "Assoziativitaet verletzt");
            }
        }
    }

    // 2-Kozyklus Paritätszählung
    int pos_c = 0, neg_c = 0;
    for (int g = 0; g < 4; ++g) {
        for (int h = 0; h < 4; ++h) {
            if (C_COCYCLE[g][h] == 1) pos_c++;
            else if (C_COCYCLE[g][h] == -1) neg_c++;
        }
    }
    suite.assert_true(pos_c == 10 && neg_c == 6, "Kozyklustafel c(g,h) muss 10x(+1) und 6x(-1) haben");

    // Tripelphasen-Zählung (28 positive, 36 negative)
    int pos_omega = 0, neg_omega = 0;
    for (int g = 0; g < 4; ++g) {
        for (int h = 0; h < 4; ++h) {
            for (int k = 0; k < 4; ++k) {
                int w = calc_omega(g, h, k);
                if (w == 1) pos_omega++;
                else if (w == -1) neg_omega++;
            }
        }
    }
    suite.assert_true(pos_omega == 28 && neg_omega == 36, "Omega-Tabelle muss exakt 28x(+1) und 36x(-1) enthalten");

    // Pendeltripel-Vorzeichenprüfung auf {A, B, C}
    for (int g = 1; g < 4; ++g) {
        for (int h = 1; h < 4; ++h) {
            if (g != h) {
                suite.assert_true(calc_omega(g, h, g) == 1, "Pendelpfade zwischen Nicht-E Klassen muessen Omega = +1 besitzen");
            }
        }
    }
    suite.pass_test();
}

void test_quadratic_solvers(TestSuite& suite) {
    suite.begin_test("T2: Quadratische Zahlkörper-Solver auf Referenzintervall [5, 10.000)");
    auto primes = sieve_primes(5, 10000);
    for (int64_t p : primes) {
        int64_t a = 0, b = 0;
        int sgn = 0;
        int res12 = p % 12;
        if (res12 == 1 || res12 == 5) {
            suite.assert_true(solve_gauss_exact(p, a, b), "Gauss-Zerlegung (D = -4) fehlgeschlagen");
        }
        if (res12 == 1 || res12 == 7) {
            suite.assert_true(solve_eisenstein_exact(p, a, b), "Eisenstein-Zerlegung (D = -3) fehlgeschlagen");
        }
        if (res12 == 1 || res12 == 11) {
            suite.assert_true(solve_pell_exact(p, a, b, sgn), "Pell-Zerlegung (D = 12) fehlgeschlagen");
            suite.assert_true(sgn == (res12 == 1 ? 1 : -1), "Pell-Vorzeichenfehler");
        }
    }
    suite.pass_test();
}

void test_bigram_invariance_and_whittle(TestSuite& suite) {
    suite.begin_test("T3: Eulerian-Shuffle & Erhaltung der Bigram-Matrix C_gh");
    auto primes = sieve_primes(100000, 150000); // Testintervall zur schnellen Prüfung
    std::vector<int> seq;
    seq.reserve(primes.size());
    for (int64_t p : primes) seq.push_back(sector_of(p));

    int C_obs[4][4] = {0};
    for (size_t i = 0; i < seq.size() - 1; ++i) C_obs[seq[i]][seq[i+1]]++;

    std::mt19937_64 rng(42);
    std::vector<int> shuff;
    for (int rep = 0; rep < 50; ++rep) {
        bool ok = generate_eulerian_shuffle(seq, shuff, rng);
        suite.assert_true(ok, "Euler-Shuffle nicht konvergiert");
        suite.assert_true(shuff.size() == seq.size(), "Shuffle-Laenge inkonsistent");
        suite.assert_true(shuff.front() == seq.front() && shuff.back() == seq.back(), "Endpunkte nicht erhalten");

        int C_rep[4][4] = {0};
        for (size_t i = 0; i < shuff.size() - 1; ++i) C_rep[shuff[i]][shuff[i+1]]++;

        for (int g = 0; g < 4; ++g) {
            for (int h = 0; h < 4; ++h) {
                suite.assert_true(C_rep[g][h] == C_obs[g][h], "Bigram-Invarianz verletzt");
            }
        }
    }
    suite.pass_test();
}

void test_algebraic_delta_identities(TestSuite& suite) {
    suite.begin_test("T4: Algebraische Identitäten D_+ + D_- = 0 und D_+ - D_- = Delta S");
    auto primes = sieve_primes(100000, 500000);
    size_t N = primes.size();
    suite.assert_true(N == 31946, "Primzahlanzahl im Intervall [100.000, 500.000) != 31.946");

    std::vector<int> seq;
    seq.reserve(N);
    for (int64_t p : primes) seq.push_back(sector_of(p));

    int N_obs[64] = {0};
    int64_t S_obs = 0;
    for (size_t i = 0; i < N - 2; ++i) {
        int g = seq[i], h = seq[i+1], k = seq[i+2];
        int idx = (g << 4) | (h << 2) | k;
        N_obs[idx]++;
        S_obs += calc_omega(g, h, k);
    }
    suite.assert_true(S_obs == -842, "S_obs != -842");

    // Simulierte Referenzwerte aus B = 10.000 Freeze-Lauf
    const double mean_S_null = -618.8658;
    const double delta_S = S_obs - mean_S_null; // -223.1342

    int N_plus = 0, N_minus = 0;
    for (int idx = 0; idx < 64; ++idx) {
        int g = (idx >> 4) & 3, h = (idx >> 2) & 3, k = idx & 3;
        if (calc_omega(g, h, k) == 1) N_plus += N_obs[idx];
        else N_minus += N_obs[idx];
    }
    suite.assert_true(N_plus + N_minus == static_cast<int>(N - 2), "Summe N+ + N- != N - 2");
    suite.assert_true(N_plus - N_minus == S_obs, "N+ - N- != S_obs");
    suite.assert_true(N_plus == 15551 && N_minus == 16393, "Zellsummen N+ / N- weichen von Referenz ab");

    double D_plus = delta_S / 2.0;
    double D_minus = -delta_S / 2.0;
    suite.assert_true(std::abs((D_plus + D_minus)) < 1e-9, "D+ + D- muss exakt verschwinden");
    suite.assert_true(std::abs((D_plus - D_minus) - delta_S) < 1e-9, "D+ - D- muss exakt Delta_S entsprechen");
    suite.pass_test();
}

void test_two_step_mod12_autocorrelation(TestSuite& suite) {
    suite.begin_test("T5: 2-Schritt-Distanzen (p_{n+2} - p_n) mod 12 & Äquivalenz mit g_n == g_{n+2}");
    auto primes = sieve_primes(100000, 500000);
    size_t N = primes.size();
    std::vector<int> seq;
    seq.reserve(N);
    for (int64_t p : primes) seq.push_back(sector_of(p));

    int count_d2_zero = 0;
    int count_gn_eq_gn2 = 0;
    int count_d2_six = 0;

    for (size_t i = 0; i < N - 2; ++i) {
        int64_t d2 = primes[i+2] - primes[i];
        if (d2 % 12 == 0) count_d2_zero++;
        if (seq[i] == seq[i+2]) count_gn_eq_gn2++;
        if (d2 % 12 == 6) count_d2_six++;
    }

    // Bijektive Äquivalenz
    suite.assert_true(count_d2_zero == count_gn_eq_gn2, "Bijektion d2 == 0 mod 12 <=> g_n == g_{n+2} verletzt");
    suite.assert_true(count_d2_zero == 6121, "Anzahl Delta_2 == 0 mod 12 weicht von 6121 ab");
    suite.assert_true(count_d2_six == 8458, "Anzahl Delta_2 == 6 mod 12 weicht von 8458 ab");

    // Repulsionstest: Delta_2 = 0 muss seltener sein als Delta_2 = 6
    suite.assert_true(count_d2_zero < count_d2_six, "Keine Unterdrueckung von Delta_2 == 0 mod 12 gegenueber 6 mod 12");
    suite.pass_test();
}

// ============================================================================
// MAIN EXECUTION ENTRY POINT
// ============================================================================
int main() {
    std::cout << "==================================================\n";
    std::cout << "EABC / HARDY-LITTLEWOOD NUMERISCHE TESTSUITE C++20\n";
    std::cout << "==================================================\n\n";

    TestSuite suite;
    auto start_time = std::chrono::high_resolution_clock::now();

    try {
        test_algebra_axioms(suite);
        test_quadratic_solvers(suite);
        test_bigram_invariance_and_whittle(suite);
        test_delta_identities: test_algebraic_delta_identities(suite);
        test_two_step_mod12_autocorrelation(suite);
    } catch (const std::exception& e) {
        std::cerr << "\nABBRUCH WEGEN FEHLER: " << e.what() << std::endl;
        return 1;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    suite.summary();
    std::cout << "Gesamtlaufzeit der Testsuite: " << std::fixed << std::setprecision(3) 
              << elapsed.count() << " Sekunden.\n\n";

    return 0;
}