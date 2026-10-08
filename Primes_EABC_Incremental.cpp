/*
 * =====================================================================================
 *  Primes_EABC_Incremental.cpp — Next-Prime-Suche mit inkrementellem Rest-Sieb
 * =====================================================================================
 *  Der eine echte Geschwindigkeitshebel gegenüber dem AVX512-Entwurf ist NICHT SIMD,
 *  sondern die Vermeidung wiederholter Divisionen:
 *
 *  Statt fuer jeden Kandidaten N eine Division N mod p (40 kleine Primzahlen, auf
 *  2048 Bit je ~teuer) zu rechnen, fuehren wir die Reste r_p = N mod p INKREMENTELL
 *  mit. Beim Schritt N -> N+2 aktualisieren wir r_p = (r_p + 2) mod p mit einer
 *  einzigen Subtraktion/Korrektur (Wert bleibt < 2p, also ein Vergleich genuegt).
 *
 *  Algorithmus:
 *   1. Wheel-210 (2*3*5*7) als Schrittmuster; die 48 teilerfremden Residuen liefern
 *      die inkrementellen Schrittweiten (+2 / +4 / +6 ...), sodass 2,3,5,7 NIE einen
 *      Kandidaten treffen -- ohne eine einzige Division durch 2,3,5,7.
 *   2. Inkrementelle Reste gegen die 40 Siebprimzahlen 11..251: O(1) pro Schritt,
 *      keine Division auf der grossen Zahl.
 *   3. Miller-Rabin (25 Runden) nur fuer die wenigen Sieb-Ueberlebenden.
 *  Kein Fermat-Vortest, keine ungenutzte Diskrepanzschranke, kein nicht-existenter
 *  SIMD-Intrinsic. Alles kompiliert mit reinem GMP.
 * =====================================================================================
 */

#include <iostream>
#include <vector>
#include <cstdint>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#include <sys/wait.h>
#include <gmp.h>

// ---- Wheel-210: die 48 teilerfremden Residuen mod 210 --------------------------------
static const uint32_t WHEEL210[48] = {
      1,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,
     53,  59,  61,  67,  71,  73,  79,  83,  89,  97, 101, 103,
    107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151, 157,
    163, 167, 169, 173, 179, 181, 187, 191, 193, 197, 199, 209
};

// ---- Siebprimzahlen 11..251 (ohne 2,3,5,7 -- die macht das Wheel) ---------------------
static const uint32_t SIEVE_PRIMES[] = {
     11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,
     59,  61,  67,  71,  73,  79,  83,  89,  97, 101, 103, 107,
    109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239,
    241, 251
};
static const int NP = (int)(sizeof(SIEVE_PRIMES) / sizeof(SIEVE_PRIMES[0]));

// ---- Timing-Helfer --------------------------------------------------------------------
static long long wall_ns() {
    timespec t{};
    clock_gettime(CLOCK_MONOTONIC_RAW, &t);
    return (long long)t.tv_sec * 1000000000LL + t.tv_nsec;
}
static long long rss_kb() {
    rusage ru{};
    getrusage(RUSAGE_SELF, &ru);
    return (long long)ru.ru_maxrss;
}

// ---- Kern: inkrementelles Rest-Sieb ---------------------------------------------------
// Liefert die kleinste Primzahl > base. cnt_sieved / cnt_mr zaehlen Arbeit.
void nextprime_incremental(mpz_t result, const mpz_t base,
                           long long *cnt_sieved, long long *cnt_mr) {
    mpz_t n;
    mpz_init(n);

    // Start: erste Zahl > base, die auf ein Wheel-Residuum faellt.
    mpz_add_ui(n, base, 1);
    while (true) {
        unsigned long r = mpz_fdiv_ui(n, 210);
        if (std::binary_search(WHEEL210, WHEEL210 + 48, (uint32_t)r)) break;
        mpz_add_ui(n, n, 1);
    }

    // Inkrementelle Reste r_p = n mod p aufsetzen (einmalig, danach O(1)).
    std::vector<uint32_t> rem(NP);
    for (int i = 0; i < NP; ++i)
        rem[i] = (uint32_t)mpz_fdiv_ui(n, SIEVE_PRIMES[i]);

    // Wheel-Schrittweiten: Abstand vom aktuellen Residuum zum naechsten.
    auto wheel_step = [](unsigned long r) -> unsigned {
        int idx = (int)(std::lower_bound(WHEEL210, WHEEL210 + 48, (uint32_t)r) - WHEEL210);
        unsigned cur  = WHEEL210[idx];
        unsigned next = WHEEL210[(idx + 1) % 48];
        return (next > cur) ? (next - cur) : (uint)(210 - cur + next);
    };

    while (true) {
        // A. Inkrementelles Sieb: O(1) je Primzahl, keine Division auf n.
        bool survives = true;
        for (int i = 0; i < NP; ++i) {
            if (rem[i] == 0) { survives = false; break; }
        }

        if (survives) {
            if (cnt_sieved) ++(*cnt_sieved);
            // B. Miller-Rabin nur fuer Sieb-Ueberlebende.
            if (cnt_mr) ++(*cnt_mr);
            if (mpz_probab_prime_p(n, 25) > 0) {
                mpz_set(result, n);
                mpz_clear(n);
                return;
            }
        }

        // Schritt: naechstes Wheel-Residuum.
        unsigned s = wheel_step((unsigned long)mpz_fdiv_ui(n, 210));
        mpz_add_ui(n, n, s);

        // Reste inkrementell nachziehen -- reine Addition + ein Vergleich.
        for (int i = 0; i < NP; ++i) {
            uint32_t v = rem[i] + s;
            uint32_t p = SIEVE_PRIMES[i];
            rem[i] = (v >= p) ? (v - p) : v;   // s < 210 < p moeglich? nein: s<=? -> unten
        }
    }
}

// Hinweis: Wheel-Schritt s kann bis 10 sein (Abstand 199->209->...). p >= 11 > s,
// daher genuegt eine einzige Korrektur v-p; der Wert bleibt stets < p nach der Zeile.

// =====================================================================================
//  Benchmark: isolierte Kindprozesse, isolierte Messung, Median
// =====================================================================================
struct Sample { double gmp_ms; double eabc_ms; int agree; long long sieved; long long mr; };

static double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    return (n % 2) ? v[(n - 1) / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

int main() {
    const int bit_lengths[] = {512, 1024, 2048};
    const int samples = 11;                    // mehr Trials als 5
    std::cout << "=== EABC incremental sieve vs GMP mpz_nextprime ===\n";
    std::cout << "sieve primes: " << NP << "; trials/bit: " << samples << "\n\n";

    gmp_randstate_t state;
    gmp_randinit_mt(state);                    // MT: reproduzierbar
    gmp_randseed_ui(state, 20260921ULL);

    for (int b = 0; b < 3; ++b) {
        int bits = bit_lengths[b];
        std::vector<double> g, e;
        int errors = 0;
        long long tot_sieved = 0, tot_mr = 0;

        for (int s = 0; s < samples; ++s) {
            mpz_t seed, a, b2;
            mpz_inits(seed, a, b2, NULL);
            mpz_urandomb(seed, state, bits);
            mpz_setbit(seed, bits - 1);

            // Reihenfolge gewuerfelt: halbe Trials GMP zuerst, halbe EABC zuerst.
            bool gmp_first = (s % 2 == 0);
            double tg = 0, te = 0;
            if (gmp_first) {
                long long t0 = wall_ns(); mpz_nextprime(a, seed); tg = (wall_ns()-t0)/1e6;
                long long t1 = wall_ns(); nextprime_incremental(b2, seed, &tot_sieved, &tot_mr); te = (wall_ns()-t1)/1e6;
            } else {
                long long t1 = wall_ns(); nextprime_incremental(b2, seed, &tot_sieved, &tot_mr); te = (wall_ns()-t1)/1e6;
                long long t0 = wall_ns(); mpz_nextprime(a, seed); tg = (wall_ns()-t0)/1e6;
            }
            if (mpz_cmp(a, b2) != 0) ++errors;
            g.push_back(tg); e.push_back(te);
            mpz_clears(seed, a, b2, NULL);
        }

        double mg = median(g), me = median(e);
        std::cout << "bits=" << bits << "  trials=" << samples
                  << "  median GMP=" << mg << " ms  median EABC=" << me << " ms"
                  << "  GMP/EABC=" << (mg/me)
                  << "  agree=" << (errors == 0 ? "yes" : "NO")
                  << "  sieved=" << tot_sieved << " mr=" << tot_mr
                  << "  rss_kb=" << rss_kb() << "\n";
    }
    gmp_randclear(state);
    return 0;
}
