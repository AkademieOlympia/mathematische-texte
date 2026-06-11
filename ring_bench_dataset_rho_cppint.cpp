
/*
 * ring_bench_dataset_rho_cppint.cpp
 *
 * Prefilter ist vollständig deaktiviert.
 *
 * Fokus: robuste Testgrundlage + Benchmark für
 *   - echte (wahrscheinlich) Primzahlen (Miller-Rabin via Boost)
 *   - Semiprimes (2 Faktoren)
 *   - Triples (3 Faktoren)
 *   - Quadrupel (4 Faktoren)
 *   - optional: Bamberg-Quadrupel (Faktoren in Residuenklassen 1,5,7,11 mod 12)
 *
 * Bigint:
 *   - Boost.Multiprecision cpp_int (header-only), problemlos bis 256 Bit (und darüber).
 *   - Miller-Rabin: boost::multiprecision::miller_rabin_test (schnell + zuverlässig in Boost).
 *
 * Build:
 *   g++ -std=c++17 -O3 -march=native -o ring_bench_dataset_rho ring_bench_dataset_rho_cppint.cpp
 *
 * Beispiel:
 *   ./ring_bench_dataset_rho --bits 128 --count 5 --mr-rounds 16 --seed 42 --bamberg on
 *   ./ring_bench_dataset_rho --bits 256 --count 3 --mr-rounds 20 --seed 7 --csv out.csv --json out.json
 */

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;
using Clock = std::chrono::steady_clock;

static inline double sec_since(const Clock::time_point& t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

static inline std::string to_dec(const cpp_int& x) {
    return x.convert_to<std::string>();
}

static inline cpp_int gcd_cpp(cpp_int a, cpp_int b) {
    while (b != 0) { cpp_int t = a % b; a = b; b = t; }
    return a;
}

// Small prime trial division to speed up prime generation and MR
static const uint32_t SMALL_PRIMES[] = {
    2u,3u,5u,7u,11u,13u,17u,19u,23u,29u,31u,37u,41u,43u,47u,53u,59u,61u,67u,71u,
    73u,79u,83u,89u,97u,101u,103u,107u,109u,113u,127u,131u,137u,139u,149u,151u,
    157u,163u,167u,173u,179u,181u,191u,193u,197u,199u,211u,223u,227u,229u,233u,
    239u,241u,251u,257u,263u,269u,271u,277u,281u,283u,293u,307u,311u,313u,317u,
    331u,337u,347u,349u,353u,359u,367u,373u,379u,383u,389u,397u,401u,409u,419u,
    421u,431u,433u,439u,443u,449u,457u,461u,463u,467u,479u,487u,491u,499u,503u,
    509u,521u,523u,541u,547u,557u,563u,569u,571u,577u,587u,593u,599u,601u,607u,
    613u,617u,619u,631u,641u,643u,647u,653u,659u,661u,673u,677u,683u,691u,701u,
    709u,719u,727u,733u,739u,743u,751u,757u,761u,769u,773u,787u,797u,809u,811u,
    821u,823u,827u,829u,839u,853u,857u,859u,863u,877u,881u,883u,887u,907u,911u,
    919u,929u,937u,941u,947u,953u,967u,971u,977u,983u,991u,997u
};

static inline bool divisible_by_small_primes(const cpp_int& n) {
    for (uint32_t p : SMALL_PRIMES) {
        if (n == p) return false;
        if ((n % p) == 0) return true;
    }
    return false;
}

// Random cpp_int with exactly "bits" bits (top bit set)
static cpp_int random_bits(unsigned bits, std::mt19937_64& rng) {
    std::uniform_int_distribution<uint64_t> dist64;
    cpp_int x = 0;
    unsigned produced = 0;
    while (produced < bits) {
        x <<= 64;
        x += dist64(rng);
        produced += 64;
    }
    if (produced > bits) {
        unsigned extra = produced - bits;
        x >>= extra;
    }
    // set top bit
    x |= (cpp_int(1) << (bits - 1));
    return x;
}

static cpp_int random_odd_candidate(unsigned bits, std::mt19937_64& rng) {
    cpp_int x = random_bits(bits, rng);
    x |= 1;
    return x;
}

static bool is_probable_prime_boost(const cpp_int& n, int rounds, std::mt19937_64& rng) {
    if (n < 2) return false;
    for (uint32_t p : SMALL_PRIMES) {
        if (n == p) return true;
        if ((n % p) == 0) return false;
    }
    if ((n & 1) == 0) return false;
    return boost::multiprecision::miller_rabin_test(n, rounds, rng);
}

// Generate probable prime with optional residue constraint mod 12 (for odd primes > 3)
static cpp_int random_probable_prime(unsigned bits, int rounds, std::mt19937_64& rng, int want_mod12 = -1) {
    if (bits < 4) bits = 4;
    for (;;) {
        cpp_int x = random_odd_candidate(bits, rng);

        // quick filters
        if (want_mod12 != -1) {
            // enforce x % 12 == want_mod12 by stepping by 2
            for (int tries = 0; tries < 4096; ++tries) {
                int r12 = int((x % 12).convert_to<int>());
                if (r12 < 0) r12 += 12;
                if (r12 == want_mod12 && (x % 3) != 0) break;
                x += 2;
            }
            int r12 = int((x % 12).convert_to<int>());
            if (r12 < 0) r12 += 12;
            if (r12 != want_mod12 || (x % 3) == 0) continue;
        }

        if (divisible_by_small_primes(x)) continue;
        if (is_probable_prime_boost(x, rounds, rng)) return x;
    }
}

// ---------------- Pollard-Rho (Brent) for cpp_int ----------------

static inline cpp_int rho_f(const cpp_int& x, const cpp_int& c, const cpp_int& n) {
    return (x * x + c) % n;
}

static cpp_int pollard_rho_brent(const cpp_int& n, std::mt19937_64& rng) {
    if ((n & 1) == 0) return 2;
    for (uint32_t p : SMALL_PRIMES) {
        if (n % p == 0) return cpp_int(p);
    }

    std::uniform_int_distribution<uint64_t> dist64;

    // Safety: bounded restarts to avoid "no output" / infinite loops
    const int MAX_RESTARTS = 64;

    for (int restart = 0; restart < MAX_RESTARTS; ++restart) {
        cpp_int y = random_bits(64, rng) % (n - 1) + 1;
        cpp_int c = random_bits(64, rng) % (n - 1) + 1;

        cpp_int m = 128;
        cpp_int g = 1, r = 1, q = 1;
        cpp_int x = 0, ys = 0;

        while (g == 1) {
            x = y;
            for (cpp_int i = 0; i < r; ++i) y = rho_f(y, c, n);

            cpp_int k = 0;
            while (k < r && g == 1) {
                ys = y;
                cpp_int rk = r - k;
                cpp_int lim = (m < rk ? m : rk);
                for (cpp_int i = 0; i < lim; ++i) {
                    y = rho_f(y, c, n);
                    cpp_int diff = (x > y) ? (x - y) : (y - x);
                    q = (q * diff) % n;
                }
                g = gcd_cpp(q, n);
                k += lim;
            }
            r <<= 1;
            if (r > (cpp_int(1) << 20)) break; // extra safety
        }

        if (g == n) {
            do {
                ys = rho_f(ys, c, n);
                cpp_int diff = (x > ys) ? (x - ys) : (ys - x);
                g = gcd_cpp(diff, n);
            } while (g == 1);
        }

        if (g != 1 && g != n) return g;
    }

    // If we failed repeatedly, fall back: likely prime-ish or hard case.
    // Caller will handle by re-testing primality with more rounds.
    return 0;
}

static void factorize_rho(cpp_int n, int mr_rounds, std::mt19937_64& rng, std::vector<cpp_int>& out) {
    if (n == 1) return;

    // If n is probable prime, accept.
    if (is_probable_prime_boost(n, mr_rounds, rng)) {
        out.push_back(n);
        return;
    }

    cpp_int d = pollard_rho_brent(n, rng);
    if (d == 0) {
        // Hard/rare: increase MR rounds; if still composite, restart recursion with different RNG state
        if (is_probable_prime_boost(n, mr_rounds + 8, rng)) {
            out.push_back(n);
            return;
        }
        // last resort: try a small random perturbation and retry
        cpp_int bump = (cpp_int(random_bits(32, rng)) | 1);
        d = gcd_cpp(n, bump);
        if (d == 1 || d == n) {
            // give up gracefully: return n as "unknown" factor
            out.push_back(n);
            return;
        }
    }

    factorize_rho(d, mr_rounds, rng, out);
    factorize_rho(n / d, mr_rounds, rng, out);
}

// ---------------- Dataset types ----------------

enum class Kind { PRIME, SEMIPRIME, TRIPLE, QUADRUPLE, BAMBERG_QUAD };

static inline const char* kind_name(Kind k) {
    switch (k) {
        case Kind::PRIME: return "prime";
        case Kind::SEMIPRIME: return "semiprime";
        case Kind::TRIPLE: return "triple";
        case Kind::QUADRUPLE: return "quadruple";
        case Kind::BAMBERG_QUAD: return "bamberg_quad";
        default: return "?";
    }
}

struct Sample {
    Kind kind;
    cpp_int n;
    std::vector<cpp_int> true_factors;
    std::vector<cpp_int> found_factors;
    double t_mr_ms = 0.0;
    double t_rho_ms = 0.0;
    bool ok = false;
};

static std::string factors_join(const std::vector<cpp_int>& v, int max_show = 6) {
    if (v.empty()) return "-";
    if ((int)v.size() <= max_show) {
        std::string s;
        for (size_t i=0;i<v.size();++i) {
            if (i) s += "*";
            s += to_dec(v[i]);
        }
        return s;
    }
    return to_dec(v.front()) + "*...*" + to_dec(v.back());
}

static Sample make_sample(Kind kind, unsigned bits, int mr_rounds, std::mt19937_64& rng) {
    Sample s;
    s.kind = kind;

    auto make_product = [&](int k, const int* want_mod12) {
        std::vector<cpp_int> facs;
        facs.reserve(k);
        unsigned each = std::max(16u, bits / (unsigned)k); // ensure factors aren't tiny
        for (int i=0;i<k;++i) {
            int want = want_mod12 ? want_mod12[i] : -1;
            facs.push_back(random_probable_prime(each, mr_rounds, rng, want));
        }
        std::sort(facs.begin(), facs.end());
        cpp_int n = 1;
        for (auto& f : facs) n *= f;
        s.n = n;
        s.true_factors = facs;
    };

    if (kind == Kind::PRIME) {
        s.n = random_probable_prime(bits, mr_rounds, rng, -1);
        s.true_factors = {s.n};
    } else if (kind == Kind::SEMIPRIME) {
        make_product(2, nullptr);
    } else if (kind == Kind::TRIPLE) {
        make_product(3, nullptr);
    } else if (kind == Kind::QUADRUPLE) {
        make_product(4, nullptr);
    } else if (kind == Kind::BAMBERG_QUAD) {
        int want[4] = {1,5,7,11};
        make_product(4, want);
    }

    return s;
}

static void run_sample(Sample& s, int mr_rounds, std::mt19937_64& rng) {
    // MR timing (primality)
    auto t0 = Clock::now();
    bool mr = is_probable_prime_boost(s.n, mr_rounds, rng);
    (void)mr;
    s.t_mr_ms = sec_since(t0) * 1000.0;

    // Factorization timing (Pollard-Rho)
    auto t1 = Clock::now();
    std::vector<cpp_int> facs;
    factorize_rho(s.n, mr_rounds, rng, facs);
    std::sort(facs.begin(), facs.end());
    s.t_rho_ms = sec_since(t1) * 1000.0;
    s.found_factors = facs;

    s.ok = (s.found_factors == s.true_factors);
}

// ---------------- IO ----------------

static void print_header() {
    std::cout
      << "i  kind          bits  ok  t_MR(ms)  t_RHO(ms)  #f  n  factors_true  factors_found\n"
      << "-------------------------------------------------------------------------------------------------------------\n";
}

static void print_row(int idx, unsigned bits, const Sample& s) {
    std::cout << std::setw(2) << idx << "  "
              << std::setw(12) << kind_name(s.kind) << "  "
              << std::setw(4) << bits << "  "
              << (s.ok ? "Y " : "N ")
              << std::fixed << std::setprecision(3)
              << std::setw(8) << s.t_mr_ms << "  "
              << std::setw(9) << s.t_rho_ms << "  "
              << std::setw(2) << s.true_factors.size() << "  "
              << to_dec(s.n) << "  "
              << factors_join(s.true_factors) << "  "
              << factors_join(s.found_factors)
              << "\n";
}

static void write_csv(const std::string& path, unsigned bits, const std::vector<Sample>& samples) {
    std::ofstream f(path);
    f << "kind,bits,n,ok,t_mr_ms,t_rho_ms,true_factors,found_factors\n";
    for (auto& s : samples) {
        f << kind_name(s.kind) << "," << bits << "," << to_dec(s.n) << "," << (s.ok?1:0) << ","
          << s.t_mr_ms << "," << s.t_rho_ms << ",\""
          << factors_join(s.true_factors, 1000) << "\",\""
          << factors_join(s.found_factors, 1000) << "\"\n";
    }
}

static void write_json(const std::string& path, unsigned bits, const std::vector<Sample>& samples) {
    std::ofstream f(path);
    f << "{\n  \"bits\": " << bits << ",\n  \"samples\": [\n";
    for (size_t i=0;i<samples.size();++i) {
        auto& s = samples[i];
        f << "    {\n";
        f << "      \"kind\": \"" << kind_name(s.kind) << "\",\n";
        f << "      \"n\": \"" << to_dec(s.n) << "\",\n";
        f << "      \"ok\": " << (s.ok ? "true" : "false") << ",\n";
        f << "      \"t_mr_ms\": " << s.t_mr_ms << ",\n";
        f << "      \"t_rho_ms\": " << s.t_rho_ms << ",\n";
        f << "      \"true_factors\": [";
        for (size_t j=0;j<s.true_factors.size();++j) {
            if (j) f << ", ";
            f << "\"" << to_dec(s.true_factors[j]) << "\"";
        }
        f << "],\n";
        f << "      \"found_factors\": [";
        for (size_t j=0;j<s.found_factors.size();++j) {
            if (j) f << ", ";
            f << "\"" << to_dec(s.found_factors[j]) << "\"";
        }
        f << "]\n";
        f << "    }" << (i+1<samples.size() ? "," : "") << "\n";
    }
    f << "  ]\n}\n";
}

// ---------------- Args ----------------

struct Args {
    unsigned bits = 128;   // typical: 128 or 256 (cpp_int supports larger too)
    int count = 5;         // per kind
    int mr_rounds = 16;
    uint64_t seed = 42;
    bool bamberg = true;
    std::string csv_path;
    std::string json_path;
};

static Args parse(int argc, char** argv) {
    Args a;
    for (int i=1;i<argc;++i) {
        std::string k = argv[i];
        auto need = [&](const std::string& name)->std::string {
            if (i+1>=argc) throw std::runtime_error("missing value for " + name);
            return argv[++i];
        };
        if (k=="--bits") a.bits = (unsigned)std::stoul(need(k));
        else if (k=="--count") a.count = std::stoi(need(k));
        else if (k=="--mr-rounds") a.mr_rounds = std::stoi(need(k));
        else if (k=="--seed") a.seed = std::stoull(need(k));
        else if (k=="--bamberg") a.bamberg = (need(k)=="on");
        else if (k=="--csv") a.csv_path = need(k);
        else if (k=="--json") a.json_path = need(k);
        else if (k=="--help" || k=="-h") {
            std::cout <<
                "Usage: ./ring_bench_dataset_rho [options]\n"
                "  --bits N            (e.g. 128 or 256)\n"
                "  --count N           (samples per kind)\n"
                "  --mr-rounds R       (MR rounds)\n"
                "  --seed SEED\n"
                "  --bamberg on|off    (include bamberg_quad)\n"
                "  --csv out.csv\n"
                "  --json out.json\n";
            std::exit(0);
        }
    }
    if (a.bits < 64) throw std::runtime_error("--bits should be >= 64");
    if (a.count <= 0) throw std::runtime_error("--count must be > 0");
    if (a.mr_rounds < 1) throw std::runtime_error("--mr-rounds must be >= 1");
    return a;
}

int main(int argc, char** argv) {
    try {
        Args args = parse(argc, argv);
        std::mt19937_64 rng(args.seed);

        std::vector<Kind> kinds = {Kind::PRIME, Kind::SEMIPRIME, Kind::TRIPLE, Kind::QUADRUPLE};
        if (args.bamberg) kinds.push_back(Kind::BAMBERG_QUAD);

        std::cout << "\n=== Dataset + Pollard-Rho Benchmark (prefilter DISABLED) ===\n";
        std::cout << "bits=" << args.bits
                  << "  count_per_kind=" << args.count
                  << "  mr_rounds=" << args.mr_rounds
                  << "  seed=" << args.seed
                  << "  bamberg_quad=" << (args.bamberg ? "on" : "off")
                  << "\n\n" << std::flush;

        print_header();

        std::vector<Sample> all;
        all.reserve((size_t)args.count * kinds.size());

        double sum_mr=0.0, sum_rho=0.0;
        int ok_cnt=0, total=0;

        int idx = 1;
        for (auto kind : kinds) {
            for (int i=0;i<args.count;++i) {
                // heartbeat BEFORE heavy generation, so user always sees progress:
                std::cout << "# progress: about to generate sample " << idx
                          << " / " << (args.count*(int)kinds.size())
                          << " (" << kind_name(kind) << ")\n" << std::flush;

                Sample s = make_sample(kind, args.bits, args.mr_rounds, rng);
                run_sample(s, args.mr_rounds, rng);

                print_row(idx, args.bits, s);

                sum_mr += s.t_mr_ms;
                sum_rho += s.t_rho_ms;
                ok_cnt += (s.ok ? 1 : 0);
                total += 1;

                all.push_back(std::move(s));
                idx++;
            }
        }

        std::cout << "\nSummary:\n";
        std::cout << "  total_samples   = " << total << "\n";
        std::cout << "  ok              = " << ok_cnt << "/" << total << "\n";
        std::cout << "  mean_t_MR(ms)   = " << std::fixed << std::setprecision(3) << (sum_mr/total) << "\n";
        std::cout << "  mean_t_RHO(ms)  = " << std::fixed << std::setprecision(3) << (sum_rho/total) << "\n\n";

        if (!args.csv_path.empty()) { write_csv(args.csv_path, args.bits, all); std::cout << "Wrote CSV: " << args.csv_path << "\n"; }
        if (!args.json_path.empty()) { write_json(args.json_path, args.bits, all); std::cout << "Wrote JSON: " << args.json_path << "\n"; }

        return (ok_cnt == total) ? 0 : 2;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        std::cerr << "Try: --help\n";
        return 1;
    }
}
