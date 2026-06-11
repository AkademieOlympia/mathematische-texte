#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using Clock = chrono::high_resolution_clock;

struct OddPrimeTable {
    int64_t limit;
    vector<uint8_t> bits; // nur ungerade Zahlen >= 3

    explicit OddPrimeTable(int64_t lim) : limit(lim) {
        if (lim >= 3) {
            size_t sz = static_cast<size_t>(((lim - 3) / 2) + 1);
            bits.assign(sz, 0);
        }
    }

    static inline size_t odd_index(int64_t n) {
        return static_cast<size_t>((n - 3) / 2);
    }

    inline void set_prime(int64_t n) {
        if (n >= 3 && (n & 1)) {
            bits[odd_index(n)] = 1;
        }
    }

    inline bool is_prime(int64_t n) const {
        if (n == 2) return true;
        if (n < 2 || n > limit || (n % 2 == 0)) return false;
        return bits[odd_index(n)] != 0;
    }
};

vector<int64_t> simple_sieve_odd(int64_t limit) {
    if (limit < 2) return {};
    if (limit == 2) return {2};

    size_t size = static_cast<size_t>(((limit - 3) / 2) + 1);
    vector<uint8_t> is_prime(size, 1);

    int64_t max_p = static_cast<int64_t>(sqrt(static_cast<double>(limit)));
    for (size_t i = 0; 2 * i + 3 <= max_p; ++i) {
        if (is_prime[i]) {
            int64_t p = 2 * i + 3;
            int64_t start = p * p;
            for (int64_t n = start; n <= limit; n += 2 * p) {
                is_prime[static_cast<size_t>((n - 3) / 2)] = 0;
            }
        }
    }

    vector<int64_t> primes;
    primes.push_back(2);
    for (size_t i = 0; i < size; ++i) {
        if (is_prime[i]) primes.push_back(2 * i + 3);
    }
    return primes;
}

pair<OddPrimeTable, vector<int64_t>> segmented_sieve_odd(int64_t limit, int64_t segment_size = 1'000'000) {
    OddPrimeTable table(limit);
    vector<int64_t> primes;
    if (limit < 2) return {std::move(table), primes};

    primes.push_back(2);
    vector<int64_t> base_primes = simple_sieve_odd(static_cast<int64_t>(sqrt(static_cast<double>(limit))));

    int64_t low = 3;
    if ((low & 1) == 0) ++low;

    while (low <= limit) {
        int64_t high = min(low + segment_size - 1, limit);
        if ((high & 1) == 0) --high;
        if (high < low) break;

        size_t seg_len = static_cast<size_t>(((high - low) / 2) + 1);
        vector<uint8_t> segment(seg_len, 1);

        for (int64_t p : base_primes) {
            if (p == 2) continue;
            int64_t p2 = p * p;
            if (p2 > high) break;

            int64_t start = max(p2, ((low + p - 1) / p) * p);
            if ((start & 1) == 0) start += p;

            for (int64_t n = start; n <= high; n += 2 * p) {
                segment[static_cast<size_t>((n - low) / 2)] = 0;
            }
        }

        for (size_t i = 0; i < seg_len; ++i) {
            if (segment[i]) {
                int64_t n = low + 2 * i;
                table.set_prime(n);
                primes.push_back(n);
            }
        }

        low = high + 2;
    }

    return {std::move(table), std::move(primes)};
}

char familie(int64_t p) {
    int64_t r = p % 12;
    if (r == 1) return 'E';
    if (r == 5) return 'A';
    if (r == 7) return 'B';
    if (r == 11) return 'C';
    return '?';
}

struct FamilyLists {
    vector<int64_t> E, A, B, C;
};

FamilyLists build_family_lists(const vector<int64_t>& primes) {
    FamilyLists fam;
    fam.E.reserve(primes.size() / 4);
    fam.A.reserve(primes.size() / 4);
    fam.B.reserve(primes.size() / 4);
    fam.C.reserve(primes.size() / 4);

    for (int64_t p : primes) {
        if (p <= 3) continue;
        switch (familie(p)) {
            case 'E': fam.E.push_back(p); break;
            case 'A': fam.A.push_back(p); break;
            case 'B': fam.B.push_back(p); break;
            case 'C': fam.C.push_back(p); break;
        }
    }
    return fam;
}

struct BenchResult {
    double search_time = 0.0;
    uint64_t twin_candidates = 0;
    uint64_t drill_i_candidates = 0;
    uint64_t drill_ii_candidates = 0;
    uint64_t quad_candidates = 0;
    uint64_t twins = 0;
    uint64_t drill_i = 0;
    uint64_t drill_ii = 0;
    uint64_t quads = 0;  // Vierlinge Typ I: (p,p+2,p+6,p+8)
};

BenchResult benchmark_direct_parallel(const OddPrimeTable& table, const vector<int64_t>& primes) {
    BenchResult res;
    auto t0 = Clock::now();

    uint64_t twin_candidates = 0, drill_i_candidates = 0, drill_ii_candidates = 0, quad_candidates = 0;
    uint64_t twins = 0, drill_i = 0, drill_ii = 0, quads = 0;

    const int64_t n = static_cast<int64_t>(primes.size());
    constexpr int64_t CHUNK_MAX = 1'000'000'000;  // OpenMP-safe: < INT_MAX

    for (int64_t chunk_start = 0; chunk_start < n; chunk_start += CHUNK_MAX) {
        const int64_t chunk_end = min(chunk_start + CHUNK_MAX, n);
        const int chunk_len = static_cast<int>(chunk_end - chunk_start);

        uint64_t l_tc = 0, l_di = 0, l_dii = 0, l_qc = 0;
        uint64_t l_tw = 0, l_d1 = 0, l_d2 = 0, l_qu = 0;

        #pragma omp parallel for reduction(+:l_tc,l_di,l_dii,l_qc,l_tw,l_d1,l_d2,l_qu) schedule(static)
        for (int i = 0; i < chunk_len; ++i) {
            int64_t p = primes[chunk_start + i];
            if (p <= 3) continue;

            if (p + 2 <= table.limit) {
                ++l_tc;
                if (table.is_prime(p + 2)) ++l_tw;
            }

            if (p + 6 <= table.limit) {
                ++l_di;
                if (table.is_prime(p + 2) && table.is_prime(p + 6)) ++l_d1;

                ++l_dii;
                if (table.is_prime(p + 4) && table.is_prime(p + 6)) ++l_d2;
            }

            if (p + 8 <= table.limit) {
                ++l_qc;
                if (table.is_prime(p + 2) && table.is_prime(p + 6) && table.is_prime(p + 8)) ++l_qu;
            }
        }

        twin_candidates += l_tc;
        drill_i_candidates += l_di;
        drill_ii_candidates += l_dii;
        quad_candidates += l_qc;
        twins += l_tw;
        drill_i += l_d1;
        drill_ii += l_d2;
        quads += l_qu;
    }

    auto t1 = Clock::now();
    res.search_time = chrono::duration<double>(t1 - t0).count();
    res.twin_candidates = twin_candidates;
    res.drill_i_candidates = drill_i_candidates;
    res.drill_ii_candidates = drill_ii_candidates;
    res.quad_candidates = quad_candidates;
    res.twins = twins;
    res.drill_i = drill_i;
    res.drill_ii = drill_ii;
    res.quads = quads;
    return res;
}

static void eabc_loop_ac(const OddPrimeTable& table, const vector<int64_t>& vec,
                         uint64_t& tc, uint64_t& di, uint64_t& qc, uint64_t& tw, uint64_t& d1, uint64_t& qu) {
    constexpr int64_t CHUNK_MAX = 1'000'000'000;
    const int64_t n = static_cast<int64_t>(vec.size());
    for (int64_t cs = 0; cs < n; cs += CHUNK_MAX) {
        int64_t ce = min(cs + CHUNK_MAX, n);
        int cl = static_cast<int>(ce - cs);
        uint64_t l_tc = 0, l_di = 0, l_qc = 0, l_tw = 0, l_d1 = 0, l_qu = 0;
        #pragma omp parallel for reduction(+:l_tc,l_di,l_qc,l_tw,l_d1,l_qu) schedule(static)
        for (int i = 0; i < cl; ++i) {
            int64_t p = vec[cs + i];
            if (p + 2 <= table.limit) { ++l_tc; if (table.is_prime(p + 2)) ++l_tw; }
            if (p + 6 <= table.limit) { ++l_di; if (table.is_prime(p + 2) && table.is_prime(p + 6)) ++l_d1; }
            if (p + 8 <= table.limit) { ++l_qc; if (table.is_prime(p + 2) && table.is_prime(p + 6) && table.is_prime(p + 8)) ++l_qu; }
        }
        tc += l_tc; di += l_di; qc += l_qc; tw += l_tw; d1 += l_d1; qu += l_qu;
    }
}

static void eabc_loop_eb(const OddPrimeTable& table, const vector<int64_t>& vec, uint64_t& dii, uint64_t& d2) {
    constexpr int64_t CHUNK_MAX = 1'000'000'000;
    const int64_t n = static_cast<int64_t>(vec.size());
    for (int64_t cs = 0; cs < n; cs += CHUNK_MAX) {
        int64_t ce = min(cs + CHUNK_MAX, n);
        int cl = static_cast<int>(ce - cs);
        uint64_t l_dii = 0, l_d2 = 0;
        #pragma omp parallel for reduction(+:l_dii,l_d2) schedule(static)
        for (int i = 0; i < cl; ++i) {
            int64_t p = vec[cs + i];
            if (p + 6 <= table.limit) { ++l_dii; if (table.is_prime(p + 4) && table.is_prime(p + 6)) ++l_d2; }
        }
        dii += l_dii; d2 += l_d2;
    }
}

BenchResult benchmark_eabc_parallel(const OddPrimeTable& table, const FamilyLists& fam) {
    BenchResult res;
    auto t0 = Clock::now();

    uint64_t twin_candidates = 0, drill_i_candidates = 0, drill_ii_candidates = 0, quad_candidates = 0;
    uint64_t twins = 0, drill_i = 0, drill_ii = 0, quads = 0;

    eabc_loop_ac(table, fam.A, twin_candidates, drill_i_candidates, quad_candidates, twins, drill_i, quads);
    eabc_loop_ac(table, fam.C, twin_candidates, drill_i_candidates, quad_candidates, twins, drill_i, quads);
    eabc_loop_eb(table, fam.E, drill_ii_candidates, drill_ii);
    eabc_loop_eb(table, fam.B, drill_ii_candidates, drill_ii);

    auto t1 = Clock::now();
    res.search_time = chrono::duration<double>(t1 - t0).count();
    res.twin_candidates = twin_candidates;
    res.drill_i_candidates = drill_i_candidates;
    res.drill_ii_candidates = drill_ii_candidates;
    res.quad_candidates = quad_candidates;
    res.twins = twins;
    res.drill_i = drill_i;
    res.drill_ii = drill_ii;
    res.quads = quads;
    return res;
}

double pct(uint64_t hit, uint64_t cand) {
    return cand ? 100.0 * static_cast<double>(hit) / static_cast<double>(cand) : 0.0;
}

void print_cluster_line(const string& name, uint64_t dc, uint64_t dh, uint64_t fc, uint64_t fh) {
    cout << name << "\n";
    cout << "  direkt: Kandidaten=" << dc
         << ", Treffer=" << dh
         << ", Quote=" << fixed << setprecision(4) << pct(dh, dc) << "%\n";
    cout << "  EABC:   Kandidaten=" << fc
         << ", Treffer=" << fh
         << ", Quote=" << fixed << setprecision(4) << pct(fh, fc) << "%\n\n";
}

void print_report(int64_t limit, int64_t segment_size = 1'000'000) {
    auto t0 = Clock::now();
    auto sieve_result = segmented_sieve_odd(limit, segment_size);
    auto t1 = Clock::now();

    OddPrimeTable table = std::move(sieve_result.first);
    vector<int64_t> primes = std::move(sieve_result.second);

    auto t2 = Clock::now();
    FamilyLists fam = build_family_lists(primes);
    auto t3 = Clock::now();

    BenchResult direct = benchmark_direct_parallel(table, primes);
    BenchResult filt = benchmark_eabc_parallel(table, fam);
    auto t4 = Clock::now();

    if (filt.twin_candidates > direct.twin_candidates)
        cerr << "FEHLER: EABC-Zwillingskandidaten > direkte Kandidaten\n";
    if (filt.drill_i_candidates > direct.drill_i_candidates)
        cerr << "FEHLER: EABC-Drilling-I-Kandidaten > direkte Kandidaten\n";
    if (filt.drill_ii_candidates > direct.drill_ii_candidates)
        cerr << "FEHLER: EABC-Drilling-II-Kandidaten > direkte Kandidaten\n";
    if (filt.quad_candidates > direct.quad_candidates)
        cerr << "FEHLER: EABC-Vierlingskandidaten > direkte Kandidaten\n";
    if (filt.twins > direct.twins || filt.drill_i > direct.drill_i ||
        filt.drill_ii > direct.drill_ii || filt.quads > direct.quads)
        cerr << "FEHLER: EABC-Treffer > direkte Treffer\n";

    double sieve_time = chrono::duration<double>(t1 - t0).count();
    double fam_time = chrono::duration<double>(t3 - t2).count();
    double total_time = chrono::duration<double>(t4 - t0).count();

    cout << string(86, '=') << "\n";
    cout << "C++ OpenMP odd-only EABC-Benchmark bis " << limit << "\n";
    cout << string(86, '=') << "\n";
#ifdef _OPENMP
    cout << "OpenMP Threads:        " << omp_get_max_threads() << "\n";
#else
    cout << "OpenMP Threads:        nicht aktiviert\n";
#endif
    cout << "Anzahl Primzahlen:     " << primes.size() << "\n";
    cout << "Odd-only Segmentsieb:  " << fixed << setprecision(6) << sieve_time << " s\n";
    cout << "Familienaufbau:        " << fam_time << " s\n";
    cout << "Direkte Suchphase:     " << direct.search_time << " s\n";
    cout << "EABC-Suchphase:        " << filt.search_time << " s\n";
    cout << "Gesamt inkl. Aufbau:   " << total_time << " s\n\n";

    cout << "Familiengrößen (ohne 2,3):\n";
    cout << "  E: " << fam.E.size() << "\n";
    cout << "  A: " << fam.A.size() << "\n";
    cout << "  B: " << fam.B.size() << "\n";
    cout << "  C: " << fam.C.size() << "\n\n";

    cout << "Kandidaten- und Treffervergleich\n";
    cout << string(86, '-') << "\n";
    print_cluster_line("reguläre Primzwillinge",
                       direct.twin_candidates, direct.twins,
                       filt.twin_candidates, filt.twins);

    print_cluster_line("primitive Drillinge Typ I (p,p+2,p+6) Abstand 6",
                       direct.drill_i_candidates, direct.drill_i,
                       filt.drill_i_candidates, filt.drill_i);

    print_cluster_line("primitive Drillinge Typ II (p,p+4,p+6) Abstand 6",
                       direct.drill_ii_candidates, direct.drill_ii,
                       filt.drill_ii_candidates, filt.drill_ii);

    print_cluster_line("echte Primvierlinge (p,p+2,p+6,p+8) Abstand 8",
                       direct.quad_candidates, direct.quads,
                       filt.quad_candidates, filt.quads);

    cout << "Geschwindigkeitsfaktor der reinen Suchphase\n";
    cout << string(86, '-') << "\n";
    if (filt.search_time > 0.0) {
        cout << "direkt / EABC = " << fixed << setprecision(3)
             << (direct.search_time / filt.search_time) << "\n";
    } else {
        cout << "direkt / EABC = n/a\n";
    }
}

int main() {
    int x;
    cout << "Exponent x fuer Limit 10^x eingeben (1-20, z.B. 6=1 Mio, 9=1 Mrd, 12=1 Bio): ";
    if (!(cin >> x) || x < 1 || x > 20) {
        cerr << "Ungueltig. Verwende x=8 (100 Mio) als Standard.\n";
        x = 8;
    }

    int64_t limit = static_cast<int64_t>(pow(10.0, x));
    int64_t segment_size = max(100'000LL, min(limit / 50, 10'000'000LL));

    cout << "Limit: " << limit << ", Segmentgroesse: " << segment_size << "\n\n";

    print_report(limit, segment_size);
    return 0;
}