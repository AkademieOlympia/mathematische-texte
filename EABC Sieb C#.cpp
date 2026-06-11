#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using Clock = chrono::high_resolution_clock;

struct OddPrimeTable {
    int limit;
    vector<uint8_t> bits; // nur ungerade Zahlen >= 3

    explicit OddPrimeTable(int lim) : limit(lim) {
        if (lim >= 3) {
            bits.assign(((lim - 3) / 2) + 1, 0);
        }
    }

    static inline int odd_index(int n) {
        return (n - 3) / 2;
    }

    inline void set_prime(int n) {
        if (n >= 3 && (n & 1)) {
            bits[odd_index(n)] = 1;
        }
    }

    inline bool is_prime(int n) const {
        if (n == 2) return true;
        if (n < 2 || n > limit || (n % 2 == 0)) return false;
        return bits[odd_index(n)] != 0;
    }
};

vector<int> simple_sieve_odd(int limit) {
    if (limit < 2) return {};
    if (limit == 2) return {2};

    int size = ((limit - 3) / 2) + 1;
    vector<uint8_t> is_prime(size, 1);

    int max_p = static_cast<int>(sqrt(limit));
    for (int i = 0; 2 * i + 3 <= max_p; ++i) {
        if (is_prime[i]) {
            int p = 2 * i + 3;
            long long start = 1LL * p * p;
            for (long long n = start; n <= limit; n += 2LL * p) {
                is_prime[(static_cast<int>(n) - 3) / 2] = 0;
            }
        }
    }

    vector<int> primes;
    primes.push_back(2);
    for (int i = 0; i < size; ++i) {
        if (is_prime[i]) primes.push_back(2 * i + 3);
    }
    return primes;
}

pair<OddPrimeTable, vector<int>> segmented_sieve_odd(int limit, int segment_size = 1'000'000) {
    OddPrimeTable table(limit);
    vector<int> primes;
    if (limit < 2) return {std::move(table), primes};

    primes.push_back(2);
    vector<int> base_primes = simple_sieve_odd(static_cast<int>(sqrt(limit)));

    int low = 3;
    if ((low & 1) == 0) ++low;

    while (low <= limit) {
        int high = min(low + segment_size - 1, limit);
        if ((high & 1) == 0) --high;
        if (high < low) break;

        int seg_len = ((high - low) / 2) + 1;
        vector<uint8_t> segment(seg_len, 1);

        for (int p : base_primes) {
            if (p == 2) continue;
            long long p2 = 1LL * p * p;
            if (p2 > high) break;

            long long start = max(p2, 1LL * ((low + p - 1) / p) * p);
            if ((start & 1) == 0) start += p;

            for (long long n = start; n <= high; n += 2LL * p) {
                segment[(static_cast<int>(n) - low) / 2] = 0;
            }
        }

        for (int i = 0; i < seg_len; ++i) {
            if (segment[i]) {
                int n = low + 2 * i;
                table.set_prime(n);
                primes.push_back(n);
            }
        }

        low = high + 2;
    }

    return {std::move(table), std::move(primes)};
}

char familie(int p) {
    int r = p % 12;
    if (r == 1) return 'E';
    if (r == 5) return 'A';
    if (r == 7) return 'B';
    if (r == 11) return 'C';
    return '?';
}

struct FamilyLists {
    vector<int> E, A, B, C;
};

FamilyLists build_family_lists(const vector<int>& primes) {
    FamilyLists fam;
    for (int p : primes) {
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
    uint64_t quads = 0;
};

BenchResult benchmark_direct(const OddPrimeTable& table, const vector<int>& primes) {
    BenchResult res;
    auto t0 = Clock::now();

    for (int p : primes) {
        if (p <= 3) continue;

        if (p + 2 <= table.limit) {
            ++res.twin_candidates;
            if (table.is_prime(p + 2)) ++res.twins;
        }

        if (p + 6 <= table.limit) {
            ++res.drill_i_candidates;
            if (table.is_prime(p + 2) && table.is_prime(p + 6)) ++res.drill_i;

            ++res.drill_ii_candidates;
            if (table.is_prime(p + 4) && table.is_prime(p + 6)) ++res.drill_ii;
        }

        if (p + 8 <= table.limit) {
            ++res.quad_candidates;
            if (table.is_prime(p + 2) && table.is_prime(p + 6) && table.is_prime(p + 8)) {
                ++res.quads;
            }
        }
    }

    auto t1 = Clock::now();
    res.search_time = chrono::duration<double>(t1 - t0).count();
    return res;
}

BenchResult benchmark_eabc_filtered(const OddPrimeTable& table, const FamilyLists& fam) {
    BenchResult res;
    auto t0 = Clock::now();

    for (int p : fam.A) {
        if (p + 2 <= table.limit) {
            ++res.twin_candidates;
            if (table.is_prime(p + 2)) ++res.twins;
        }
        if (p + 6 <= table.limit) {
            ++res.drill_i_candidates;
            if (table.is_prime(p + 2) && table.is_prime(p + 6)) ++res.drill_i;
        }
        if (p + 8 <= table.limit) {
            ++res.quad_candidates;
            if (table.is_prime(p + 2) && table.is_prime(p + 6) && table.is_prime(p + 8)) {
                ++res.quads;
            }
        }
    }

    for (int p : fam.C) {
        if (p + 2 <= table.limit) {
            ++res.twin_candidates;
            if (table.is_prime(p + 2)) ++res.twins;
        }
        if (p + 6 <= table.limit) {
            ++res.drill_i_candidates;
            if (table.is_prime(p + 2) && table.is_prime(p + 6)) ++res.drill_i;
        }
        if (p + 8 <= table.limit) {
            ++res.quad_candidates;
            if (table.is_prime(p + 2) && table.is_prime(p + 6) && table.is_prime(p + 8)) {
                ++res.quads;
            }
        }
    }

    for (int p : fam.E) {
        if (p + 6 <= table.limit) {
            ++res.drill_ii_candidates;
            if (table.is_prime(p + 4) && table.is_prime(p + 6)) ++res.drill_ii;
        }
    }

    for (int p : fam.B) {
        if (p + 6 <= table.limit) {
            ++res.drill_ii_candidates;
            if (table.is_prime(p + 4) && table.is_prime(p + 6)) ++res.drill_ii;
        }
    }

    auto t1 = Clock::now();
    res.search_time = chrono::duration<double>(t1 - t0).count();
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

void print_report(int limit, int segment_size = 1'000'000) {
    auto t0 = Clock::now();
    auto sieve_result = segmented_sieve_odd(limit, segment_size);
    auto t1 = Clock::now();

    OddPrimeTable table = std::move(sieve_result.first);
    vector<int> primes = std::move(sieve_result.second);

    auto t2 = Clock::now();
    FamilyLists fam = build_family_lists(primes);
    auto t3 = Clock::now();

    BenchResult direct = benchmark_direct(table, primes);
    BenchResult filt = benchmark_eabc_filtered(table, fam);
    auto t4 = Clock::now();

    double sieve_time = chrono::duration<double>(t1 - t0).count();
    double fam_time = chrono::duration<double>(t3 - t2).count();
    double total_time = chrono::duration<double>(t4 - t0).count();

    cout << string(84, '=') << "\n";
    cout << "C++ odd-only segmentierter EABC-Benchmark bis " << limit << "\n";
    cout << string(84, '=') << "\n";
    cout << "Anzahl Primzahlen:      " << primes.size() << "\n";
    cout << "Odd-only Segmentsieb:   " << fixed << setprecision(6) << sieve_time << " s\n";
    cout << "Familienaufbau:         " << fam_time << " s\n";
    cout << "Direkte Suchphase:      " << direct.search_time << " s\n";
    cout << "EABC-Suchphase:         " << filt.search_time << " s\n";
    cout << "Gesamt inkl. Aufbau:    " << total_time << " s\n\n";

    cout << "Familiengrößen (ohne 2,3):\n";
    cout << "  E: " << fam.E.size() << "\n";
    cout << "  A: " << fam.A.size() << "\n";
    cout << "  B: " << fam.B.size() << "\n";
    cout << "  C: " << fam.C.size() << "\n\n";

    cout << "Kandidaten- und Treffervergleich\n";
    cout << string(84, '-') << "\n";
    print_cluster_line("reguläre Primzwillinge",
                       direct.twin_candidates, direct.twins,
                       filt.twin_candidates, filt.twins);

    print_cluster_line("primitive Drillinge Typ I",
                       direct.drill_i_candidates, direct.drill_i,
                       filt.drill_i_candidates, filt.drill_i);

    print_cluster_line("primitive Drillinge Typ II",
                       direct.drill_ii_candidates, direct.drill_ii,
                       filt.drill_ii_candidates, filt.drill_ii);

    print_cluster_line("echte Primvierlinge",
                       direct.quad_candidates, direct.quads,
                       filt.quad_candidates, filt.quads);

    cout << "Geschwindigkeitsfaktor der reinen Suchphase\n";
    cout << string(84, '-') << "\n";
    if (filt.search_time > 0.0) {
        cout << "direkt / EABC = " << fixed << setprecision(3)
             << (direct.search_time / filt.search_time) << "\n";
    } else {
        cout << "direkt / EABC = n/a\n";
    }
}

int main() {
    // print_report(5'000'000, 500'000);

    // Für größere Läufe:
    print_report(10'000'000, 1'000'000);
    print_report(50'000'000, 2'000'000);
    print_report(100'000'000, 5'000'000);  // 10^8

    return 0;
}