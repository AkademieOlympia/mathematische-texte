
/*
 * ring_bench_dataset_rho_boostfixed.cpp
 *
 * Optimierte Benchmark/Testbasis:
 *  - Prefilter komplett deaktiviert
 *  - Dataset: prime / semiprime / triple / quadruple / optional bamberg_quad (1,5,7,11 mod 12)
 *  - Benchmark: Boost Miller-Rabin (Primality) + Pollard-Rho (Brent)
 *
 * Speed-Optimierung:
 *  - Für bits<=128: boost::multiprecision::uint128_t mit uint256_t als "wider" Produkt
 *  - Für bits<=256: boost::multiprecision::uint256_t mit uint512_t als "wider" Produkt
 *  - Für bits>256: fallback auf boost::multiprecision::cpp_int
 *
 * Das reduziert dynamische Allokationen gegenüber cpp_int massiv und macht Rho auf 128/256-bit
 * typischerweise dramatisch schneller und stabiler.
 *
 * Build:
 *   g++ -std=c++17 -O3 -march=native -o ring_bench_dataset_rho_fast ring_bench_dataset_rho_boostfixed.cpp
 *
 * Beispiel:
 *   ./ring_bench_dataset_rho_fast --bits 128 --count 5 --mr-rounds 16 --seed 42 --bamberg on
 *   ./ring_bench_dataset_rho_fast --bits 256 --count 2 --mr-rounds 20 --seed 7 --bamberg on
 *
 * Optional:
 *   --rho-restarts N   (default 128)
 *   --rho-m M          (default 128)  block size in Brent
 *   --rho-timeout-ms T (default 0 = off) hartes Timeout pro Sample
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

namespace bmp = boost::multiprecision;
using Clock = std::chrono::steady_clock;

static inline double sec_since(const Clock::time_point& t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

// ---------- small primes ----------
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

// ---------- generic helpers ----------
template<class Int>
static inline Int gcd_int(Int a, Int b) {
    while (b != 0) { Int t = a % b; a = b; b = t; }
    return a;
}

template<class Int>
static inline bool divisible_by_small_primes(const Int& n) {
    for (uint32_t p : SMALL_PRIMES) {
        if (n == p) return false;
        if ((n % p) == 0) return true;
    }
    return false;
}

template<class Int>
static inline std::string to_dec(const Int& x) {
    return x.template convert_to<std::string>();
}

template<class Int>
static inline Int rand_bits(unsigned bits, std::mt19937_64& rng) {
    std::uniform_int_distribution<uint64_t> dist64;
    Int x = 0;
    unsigned produced = 0;
    while (produced < bits) {
        x <<= 64;
        x += dist64(rng);
        produced += 64;
    }
    if (produced > bits) x >>= (produced - bits);
    if (bits > 0) x |= (Int(1) << (bits - 1)); // top bit set
    return x;
}

template<class Int>
static inline Int rand_odd_candidate(unsigned bits, std::mt19937_64& rng) {
    Int x = rand_bits<Int>(bits, rng);
    x |= 1;
    return x;
}

template<class Int>
static inline bool is_probable_prime(const Int& n, int rounds, std::mt19937_64& rng) {
    if (n < 2) return false;
    for (uint32_t p : SMALL_PRIMES) {
        if (n == p) return true;
        if ((n % p) == 0) return false;
    }
    if ((n & 1) == 0) return false;
    return bmp::miller_rabin_test(n, rounds, rng);
}

// enforce residue class mod 12 (odd primes > 3 are in {1,5,7,11})
template<class Int>
static Int random_probable_prime(unsigned bits, int rounds, std::mt19937_64& rng, int want_mod12 = -1) {
    if (bits < 8) bits = 8;
    for (;;) {
        Int x = rand_odd_candidate<Int>(bits, rng);

        if (want_mod12 != -1) {
            // adjust by steps of 2 until residue fits (bounded)
            for (int tries = 0; tries < 8192; ++tries) {
                int r12 = int((x % 12).template convert_to<int>());
                if (r12 < 0) r12 += 12;
                if (r12 == want_mod12 && (x % 3) != 0) break;
                x += 2;
            }
            int r12 = int((x % 12).template convert_to<int>());
            if (r12 < 0) r12 += 12;
            if (r12 != want_mod12 || (x % 3) == 0) continue;
        }

        if (divisible_by_small_primes(x)) continue;
        if (is_probable_prime(x, rounds, rng)) return x;
    }
}

// ---------- fast mulmod via "wider" type ----------
template<class Int, class Wider>
static inline Int mul_mod(const Int& a, const Int& b, const Int& n) {
    return Int((Wider(a) * Wider(b)) % Wider(n));
}

// ---------- Pollard-Rho (Brent) ----------
template<class Int, class Wider>
static inline Int rho_f(const Int& x, const Int& c, const Int& n) {
    // x^2 + c mod n
    return (mul_mod<Int, Wider>(x, x, n) + c) % n;
}

template<class Int, class Wider>
static Int pollard_rho_brent(const Int& n,
                            std::mt19937_64& rng,
                            int max_restarts,
                            uint32_t m_block,
                            uint64_t timeout_ms,
                            bool& timed_out)
{
    timed_out = false;
    if ((n & 1) == 0) return 2;

    for (uint32_t p : SMALL_PRIMES) {
        if (n % p == 0) return Int(p);
    }

    std::uniform_int_distribution<uint64_t> dist64;
    auto t0 = Clock::now();

    auto time_exceeded = [&]() -> bool {
        if (timeout_ms == 0) return false;
        double ms = sec_since(t0) * 1000.0;
        return ms > double(timeout_ms);
    };

    for (int restart = 0; restart < max_restarts; ++restart) {
        if (time_exceeded()) { timed_out = true; return 0; }

        Int y = Int(dist64(rng)) % (n - 1) + 1;
        Int c = Int(dist64(rng)) % (n - 1) + 1;

        Int m = Int(m_block);
        Int g = 1, r = 1, q = 1;
        Int x = 0, ys = 0;

        while (g == 1) {
            if (time_exceeded()) { timed_out = true; return 0; }

            x = y;
            for (Int i = 0; i < r; ++i) y = rho_f<Int, Wider>(y, c, n);

            Int k = 0;
            while (k < r && g == 1) {
                if (time_exceeded()) { timed_out = true; return 0; }

                ys = y;
                Int rk = r - k;
                Int lim = (m < rk ? m : rk);

                q = 1;
                for (Int i = 0; i < lim; ++i) {
                    y = rho_f<Int, Wider>(y, c, n);
                    Int diff = (x > y) ? (x - y) : (y - x);
                    if (diff != 0) q = mul_mod<Int, Wider>(q, diff, n);
                }
                g = gcd_int(q, n);
                k += lim;
            }
            r <<= 1;
        }

        if (g == n) {
            do {
                if (time_exceeded()) { timed_out = true; return 0; }
                ys = rho_f<Int, Wider>(ys, c, n);
                Int diff = (x > ys) ? (x - ys) : (ys - x);
                g = gcd_int(diff, n);
            } while (g == 1);
        }

        if (g != 1 && g != n) return g;
    }

    return 0;
}

template<class Int, class Wider>
static void factorize_rho(Int n,
                          int mr_rounds,
                          std::mt19937_64& rng,
                          std::vector<Int>& out,
                          int rho_restarts,
                          uint32_t rho_m,
                          uint64_t rho_timeout_ms)
{
    if (n == 1) return;

    // primality fast path
    if (is_probable_prime(n, mr_rounds, rng)) {
        out.push_back(n);
        return;
    }

    // try to split; if it fails, keep retrying with higher restarts/rounds (bounded by timeout)
    for (int attempt = 0; attempt < 6; ++attempt) {
        bool timed_out = false;
        Int d = pollard_rho_brent<Int, Wider>(n, rng, rho_restarts, rho_m, rho_timeout_ms, timed_out);

        if (d != 0 && d != 1 && d != n) {
            factorize_rho<Int, Wider>(d, mr_rounds, rng, out, rho_restarts, rho_m, rho_timeout_ms);
            factorize_rho<Int, Wider>(n / d, mr_rounds, rng, out, rho_restarts, rho_m, rho_timeout_ms);
            return;
        }

        // if timed out: give up gracefully, but do NOT claim success
        if (timed_out) break;

        // maybe we misclassified: increase MR rounds a bit (rare) and re-check
        if (is_probable_prime(n, mr_rounds + 8 + 2*attempt, rng)) {
            out.push_back(n);
            return;
        }

        // increase effort for next attempt
        rho_restarts = std::min(rho_restarts * 2, 2048);
        rho_m = std::min<uint32_t>(rho_m * 2, 2048);
    }

    // Could not factor: return n as "unknown composite/prime-ish" marker.
    out.push_back(n);
}

// ---------- dataset ----------
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

template<class Int>
static std::string join_factors(const std::vector<Int>& v, int max_show = 6) {
    if (v.empty()) return "-";
    if ((int)v.size() <= max_show) {
        std::string s;
        for (size_t i=0;i<v.size();++i) { if (i) s += "*"; s += to_dec(v[i]); }
        return s;
    }
    return to_dec(v.front()) + "*...*" + to_dec(v.back());
}

template<class Int>
struct Sample {
    Kind kind;
    Int n;
    std::vector<Int> true_factors;
    std::vector<Int> found_factors;
    double t_mr_ms = 0.0;
    double t_rho_ms = 0.0;
    bool ok = false;
};

struct Args {
    unsigned bits = 128;
    int count = 5;
    int mr_rounds = 16;
    uint64_t seed = 42;
    bool bamberg = true;
    int rho_restarts = 128;
    uint32_t rho_m = 128;
    uint64_t rho_timeout_ms = 0; // 0 = off
    std::string csv_path;
    std::string json_path;
    bool resample_on_fail = true;
    int max_resample = 50;
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
        else if (k=="--rho-restarts") a.rho_restarts = std::stoi(need(k));
        else if (k=="--rho-m") a.rho_m = (uint32_t)std::stoul(need(k));
        else if (k=="--rho-timeout-ms") a.rho_timeout_ms = (uint64_t)std::stoull(need(k));
        else if (k=="--csv") a.csv_path = need(k);
        else if (k=="--json") a.json_path = need(k);
        else if (k=="--resample-on-fail") a.resample_on_fail = (need(k)=="on");
        else if (k=="--max-resample") a.max_resample = std::stoi(need(k));
        else if (k=="--help" || k=="-h") {
            std::cout <<
                "Usage: ./ring_bench_dataset_rho_fast [options]\n"
                "  --bits N              (<=128 uses uint128_t, <=256 uses uint256_t, else cpp_int)\n"
                "  --count N             (samples per kind)\n"
                "  --mr-rounds R         (MR rounds)\n"
                "  --seed SEED\n"
                "  --bamberg on|off\n"
                "  --rho-restarts N       (default 128)\n"
                "  --rho-m M              (default 128)\n"
                "  --rho-timeout-ms T      (0=off; hard timeout per factorization attempt)\n"
                "  --csv out.csv\n"
                "  --json out.json\n";
            std::exit(0);
        }
    }
    if (a.bits < 64) throw std::runtime_error("--bits should be >= 64");
    if (a.count <= 0) throw std::runtime_error("--count must be > 0");
    if (a.mr_rounds < 1) throw std::runtime_error("--mr-rounds must be >= 1");
    if (a.rho_restarts < 8) a.rho_restarts = 8;
    if (a.rho_m < 8) a.rho_m = 8;
    if (a.max_resample < 0) a.max_resample = 0;
    return a;
}

template<class Int>
static Sample<Int> make_sample(Kind kind, unsigned bits, int mr_rounds, std::mt19937_64& rng) {
    Sample<Int> s;
    s.kind = kind;

    auto make_product = [&](int k, const int* want_mod12) {
        std::vector<Int> facs;
        facs.reserve(k);
        unsigned each = std::max(16u, bits / (unsigned)k);
        for (int i=0;i<k;++i) {
            int want = want_mod12 ? want_mod12[i] : -1;
            facs.push_back(random_probable_prime<Int>(each, mr_rounds, rng, want));
        }
        std::sort(facs.begin(), facs.end());
        Int n = 1;
        for (auto& f : facs) n *= f;
        s.n = n;
        s.true_factors = facs;
    };

    if (kind == Kind::PRIME) {
        s.n = random_probable_prime<Int>(bits, mr_rounds, rng, -1);
        s.true_factors = {s.n};
    } else if (kind == Kind::SEMIPRIME) make_product(2, nullptr);
    else if (kind == Kind::TRIPLE) make_product(3, nullptr);
    else if (kind == Kind::QUADRUPLE) make_product(4, nullptr);
    else if (kind == Kind::BAMBERG_QUAD) { int want[4] = {1,5,7,11}; make_product(4, want); }

    return s;
}

template<class Int, class Wider>
static void run_sample(Sample<Int>& s, int mr_rounds, std::mt19937_64& rng,
                       int rho_restarts, uint32_t rho_m, uint64_t rho_timeout_ms)
{
    auto t0 = Clock::now();
    (void)is_probable_prime(s.n, mr_rounds, rng);
    s.t_mr_ms = sec_since(t0) * 1000.0;

    auto t1 = Clock::now();
    std::vector<Int> facs;
    factorize_rho<Int, Wider>(s.n, mr_rounds, rng, facs, rho_restarts, rho_m, rho_timeout_ms);
    std::sort(facs.begin(), facs.end());
    s.t_rho_ms = sec_since(t1) * 1000.0;
    s.found_factors = facs;
    s.ok = (s.found_factors == s.true_factors);
}

static void print_header() {
    std::cout
      << "i  kind          bits  ok  t_MR(ms)  t_RHO(ms)  #f  n  factors_true  factors_found\n"
      << "-------------------------------------------------------------------------------------------------------------\n";
}

template<class Int>
static void print_row(int idx, unsigned bits, const Sample<Int>& s) {
    std::cout << std::setw(2) << idx << "  "
              << std::setw(12) << kind_name(s.kind) << "  "
              << std::setw(4) << bits << "  "
              << (s.ok ? "Y " : "N ")
              << std::fixed << std::setprecision(3)
              << std::setw(8) << s.t_mr_ms << "  "
              << std::setw(9) << s.t_rho_ms << "  "
              << std::setw(2) << s.true_factors.size() << "  "
              << to_dec(s.n) << "  "
              << join_factors(s.true_factors) << "  "
              << join_factors(s.found_factors)
              << "\n";
}

template<class Int>
static void write_csv(const std::string& path, unsigned bits, const std::vector<Sample<Int>>& samples) {
    std::ofstream f(path);
    f << "kind,bits,n,ok,t_mr_ms,t_rho_ms,true_factors,found_factors\n";
    for (auto& s : samples) {
        f << kind_name(s.kind) << "," << bits << "," << to_dec(s.n) << "," << (s.ok?1:0) << ","
          << s.t_mr_ms << "," << s.t_rho_ms << ",\""
          << join_factors(s.true_factors, 1000) << "\",\""
          << join_factors(s.found_factors, 1000) << "\"\n";
    }
}

template<class Int>
static void write_json(const std::string& path, unsigned bits, const std::vector<Sample<Int>>& samples) {
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
        for (size_t j=0;j<s.true_factors.size();++j) { if (j) f << ", "; f << "\"" << to_dec(s.true_factors[j]) << "\""; }
        f << "],\n";
        f << "      \"found_factors\": [";
        for (size_t j=0;j<s.found_factors.size();++j) { if (j) f << ", "; f << "\"" << to_dec(s.found_factors[j]) << "\""; }
        f << "]\n";
        f << "    }" << (i+1<samples.size() ? "," : "") << "\n";
    }
    f << "  ]\n}\n";
}

template<class Int, class Wider>
static int run_suite(const Args& args) {
    std::mt19937_64 rng(args.seed);

    std::vector<Kind> kinds = {Kind::PRIME, Kind::SEMIPRIME, Kind::TRIPLE, Kind::QUADRUPLE};
    if (args.bamberg) kinds.push_back(Kind::BAMBERG_QUAD);

    std::cout << "\n=== Dataset + Pollard-Rho Benchmark (prefilter DISABLED, fixed-size fast path) ===\n";
    std::cout << "bits=" << args.bits
              << "  count_per_kind=" << args.count
              << "  mr_rounds=" << args.mr_rounds
              << "  seed=" << args.seed
              << "  bamberg_quad=" << (args.bamberg ? "on" : "off")
              << "  rho_restarts=" << args.rho_restarts
              << "  rho_m=" << args.rho_m
              << "  rho_timeout_ms=" << args.rho_timeout_ms
              << "\n\n" << std::flush;

    print_header();

    std::vector<Sample<Int>> all;
    all.reserve((size_t)args.count * kinds.size());

    double sum_mr=0.0, sum_rho=0.0;
    int ok_cnt=0, total=0;

    int idx = 1;
    for (auto kind : kinds) {
        for (int i=0;i<args.count;++i) {
            std::cout << "# progress: generating sample " << idx
                      << " / " << (args.count*(int)kinds.size())
                      << " (" << kind_name(kind) << ")\n" << std::flush;

            Sample<Int> s;
            bool done = false;
            for (int rep = 0; rep <= (args.resample_on_fail ? args.max_resample : 0); ++rep) {
                s = make_sample<Int>(kind, args.bits, args.mr_rounds, rng);
                run_sample<Int, Wider>(s, args.mr_rounds, rng, args.rho_restarts, args.rho_m, args.rho_timeout_ms);
                if (s.ok || !args.resample_on_fail) { done = true; break; }
                std::cout << "# note: factorization failed/timeout -> resample (" << (rep+1) << "/" << args.max_resample << ")\\n" << std::flush;
            }
            (void)done;

            print_row<Int>(idx, args.bits, s);

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

    if (!args.csv_path.empty()) { write_csv<Int>(args.csv_path, args.bits, all); std::cout << "Wrote CSV: " << args.csv_path << "\n"; }
    if (!args.json_path.empty()) { write_json<Int>(args.json_path, args.bits, all); std::cout << "Wrote JSON: " << args.json_path << "\n"; }

    return (ok_cnt == total) ? 0 : 2;
}

int main(int argc, char** argv) {
    try {
        Args args = parse(argc, argv);

        if (args.bits <= 128) {
            using Int = bmp::uint128_t;
            using Wider = bmp::uint256_t;
            return run_suite<Int, Wider>(args);
        } else if (args.bits <= 256) {
            using Int = bmp::uint256_t;
            using Wider = bmp::uint512_t;
            return run_suite<Int, Wider>(args);
        } else {
            using Int = bmp::cpp_int;
            using Wider = bmp::cpp_int;
            return run_suite<Int, Wider>(args);
        }

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\nTry: --help\n";
        return 1;
    }
}
