
// ring_bench_e8_compare_v5b.cpp (see chat for usage)
// NOTE: E8 here is a placeholder "E8-Lite" factor finder so E8 can actually participate in benchmarks.
// Replace e8_find_factor(...) with your true E8 algorithm.

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/multiprecision/integer.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace bmp = boost::multiprecision;
using Clock = std::chrono::steady_clock;

static inline double ms_since(const Clock::time_point& t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

static inline std::string trim(const std::string& s) {
    size_t a=0,b=s.size();
    while (a<b && std::isspace((unsigned char)s[a])) a++;
    while (b>a && std::isspace((unsigned char)s[b-1])) b--;
    return s.substr(a,b-a);
}

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

template<class Int>
static inline Int gcd_int(Int a, Int b) {
    while (b != 0) { Int t = a % b; a = b; b = t; }
    return a;
}
template<class Int>
static inline std::string to_dec(const Int& x) { return x.template convert_to<std::string>(); }

template<class Int>
static inline std::string short_dec(const Int& x, int d, bool show_full) {
    std::string s = to_dec(x);
    if (show_full) return s;
    if ((int)s.size() <= 2*d + 3) return s;
    return s.substr(0, d) + "..." + s.substr((int)s.size() - d);
}

template<class Int>
static inline std::string join_factors_full(const std::vector<Int>& v) {
    if (v.empty()) return "";
    std::ostringstream os;
    for (size_t i=0;i<v.size();++i) { if (i) os << "*"; os << to_dec(v[i]); }
    return os.str();
}
template<class Int>
static inline std::string join_factors_short(const std::vector<Int>& v, int d, bool show_full) {
    if (v.empty()) return "-";
    if (show_full) return join_factors_full(v);
    std::ostringstream os;
    for (size_t i=0;i<v.size();++i) {
        if (i) os << "*";
        os << short_dec(v[i], d, false);
        if (i >= 3 && v.size() > 4) { os << "*..."; break; }
    }
    return os.str();
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
    if (bits > 0) x |= (Int(1) << (bits - 1));
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

template<class Int>
static Int random_probable_prime(unsigned bits, int rounds, std::mt19937_64& rng, int want_mod12 = -1) {
    if (bits < 16) bits = 16;
    for (;;) {
        Int x = rand_odd_candidate<Int>(bits, rng);

        if (want_mod12 != -1) {
            for (int tries=0; tries<8192; ++tries) {
                int r12 = int((x % 12).template convert_to<int>());
                if (r12 < 0) r12 += 12;
                if (r12 == want_mod12 && (x % 3) != 0) break;
                x += 2;
            }
            int r12 = int((x % 12).template convert_to<int>());
            if (r12 < 0) r12 += 12;
            if (r12 != want_mod12 || (x % 3) == 0) continue;
        }

        bool div = false;
        for (uint32_t p : SMALL_PRIMES) { if ((x % p) == 0 && x != p) { div = true; break; } }
        if (div) continue;

        if (is_probable_prime(x, rounds, rng)) return x;
    }
}

struct Deadline {
    Clock::time_point t0;
    uint64_t hard_ms; // 0 unlimited
    Deadline(uint64_t hard_ms_=0): t0(Clock::now()), hard_ms(hard_ms_) {}
    inline bool exceeded() const {
        if (!hard_ms) return false;
        return std::chrono::duration<double, std::milli>(Clock::now() - t0).count() > double(hard_ms);
    }
};

template<class Int, class Wider>
static inline Int mul_mod(const Int& a, const Int& b, const Int& n) {
    return Int((Wider(a) * Wider(b)) % Wider(n));
}
template<class Int, class Wider>
static inline Int pow_mod(Int a, Int e, const Int& n) {
    Int r = 1 % n;
    a %= n;
    while (e > 0) {
        if ((e & 1) != 0) r = mul_mod<Int, Wider>(r, a, n);
        e >>= 1;
        a = mul_mod<Int, Wider>(a, a, n);
    }
    return r;
}
template<class Int, class Wider>
static inline Int rho_f(const Int& x, const Int& c, const Int& n) {
    return (mul_mod<Int, Wider>(x, x, n) + c) % n;
}

template<class Int, class Wider>
static Int pollard_rho_brent(const Int& n, std::mt19937_64& rng, int max_restarts, uint32_t m_block, Deadline& dl, bool& timed_out) {
    timed_out = false;
    if ((n & 1) == 0) return 2;
    for (uint32_t p : SMALL_PRIMES) if (n % p == 0) return Int(p);

    std::uniform_int_distribution<uint64_t> dist64;
    for (int restart=0; restart<max_restarts; ++restart) {
        if (dl.exceeded()) { timed_out = true; return 0; }
        Int y = Int(dist64(rng)) % (n - 1) + 1;
        Int c = Int(dist64(rng)) % (n - 1) + 1;

        Int g = 1, r = 1, q = 1;
        Int x = 0, ys = 0;
        Int m = Int(m_block);

        while (g == 1) {
            if (dl.exceeded()) { timed_out = true; return 0; }

            x = y;
            for (Int i=0;i<r;++i) {
                y = rho_f<Int, Wider>(y, c, n);
                if (dl.exceeded()) { timed_out = true; return 0; }
            }

            Int k = 0;
            while (k < r && g == 1) {
                if (dl.exceeded()) { timed_out = true; return 0; }

                ys = y;
                Int rk = r - k;
                Int lim = (m < rk ? m : rk);

                q = 1;
                for (Int i=0;i<lim;++i) {
                    y = rho_f<Int, Wider>(y, c, n);
                    Int diff = (x > y) ? (x - y) : (y - x);
                    if (diff != 0) q = mul_mod<Int, Wider>(q, diff, n);
                    if (dl.exceeded()) { timed_out = true; return 0; }
                }
                g = gcd_int(q, n);
                k += lim;
            }
            r <<= 1;
        }

        if (g == n) {
            do {
                if (dl.exceeded()) { timed_out = true; return 0; }
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
static bool rho_factorize(Int n, int mr_rounds, std::mt19937_64& rng, std::vector<Int>& out, int rho_restarts, uint32_t rho_m, Deadline& dl, bool& timed_out) {
    timed_out = false;
    if (n == 1) return true;
    if (dl.exceeded()) { timed_out = true; return false; }

    if (is_probable_prime(n, mr_rounds, rng)) { out.push_back(n); return true; }

    for (int attempt=0; attempt<12; ++attempt) {
        if (dl.exceeded()) { timed_out = true; return false; }
        bool to=false;
        Int d = pollard_rho_brent<Int, Wider>(n, rng, rho_restarts, rho_m, dl, to);
        if (to) { timed_out = true; return false; }
        if (d != 0 && d != 1 && d != n) {
            bool to1=false,to2=false;
            if (!rho_factorize<Int, Wider>(d, mr_rounds, rng, out, rho_restarts, rho_m, dl, to1)) { timed_out = to1; return false; }
            if (!rho_factorize<Int, Wider>(n/d, mr_rounds, rng, out, rho_restarts, rho_m, dl, to2)) { timed_out = to2; return false; }
            return true;
        }
        rho_restarts = std::min(rho_restarts*2, 8192);
        rho_m = std::min<uint32_t>(rho_m*2, 8192);
    }
    return false;
}

// ---------- E8 hook (E8-Lite placeholder) ----------
enum class E8Status { DISABLED, HIT, MISS, TIMEOUT, ERROR };
static inline const char* e8_status_name(E8Status s) {
    switch (s) {
        case E8Status::DISABLED: return "DISABLED";
        case E8Status::HIT: return "HIT";
        case E8Status::MISS: return "MISS";
        case E8Status::TIMEOUT: return "TIMEOUT";
        case E8Status::ERROR: return "ERROR";
        default: return "?";
    }
}
template<class Int>
struct E8Run {
    E8Status status = E8Status::DISABLED;
    Int factor = 0;
    double t_ms = 0.0;
};

template<class Int, class Wider>
static E8Run<Int> e8_find_factor(const Int& n, std::mt19937_64& rng, Deadline& dl,
                                 uint32_t trial_bound = 200000,
                                 uint32_t fermat_steps = 20000,
                                 uint32_t p1_B1 = 20000)
{
    (void)rng;
    E8Run<Int> r;
    auto t0 = Clock::now();
    try {
        if (dl.exceeded()) { r.status = E8Status::TIMEOUT; r.t_ms = ms_since(t0); return r; }
        if (n < 4) { r.status = E8Status::MISS; r.t_ms = ms_since(t0); return r; }
        if ((n & 1) == 0) { r.status = E8Status::HIT; r.factor = 2; r.t_ms = ms_since(t0); return r; }

        for (uint32_t p : SMALL_PRIMES) {
            if (dl.exceeded()) { r.status = E8Status::TIMEOUT; r.t_ms = ms_since(t0); return r; }
            if (n % p == 0) { r.status = E8Status::HIT; r.factor = Int(p); r.t_ms = ms_since(t0); return r; }
        }

        for (uint32_t f = 1001; f <= trial_bound; f += 2) {
            if (dl.exceeded()) { r.status = E8Status::TIMEOUT; r.t_ms = ms_since(t0); return r; }
            if (f % 3 == 0 || f % 5 == 0) continue;
            if (n % f == 0) { r.status = E8Status::HIT; r.factor = Int(f); r.t_ms = ms_since(t0); return r; }
        }

        // Fermat window (balanced semiprimes only)
        Int a = bmp::sqrt(n);
        if (a*a < n) a += 1;
        for (uint32_t i=0;i<fermat_steps;++i) {
            if (dl.exceeded()) { r.status = E8Status::TIMEOUT; r.t_ms = ms_since(t0); return r; }
            Int b2 = a*a - n;
            if (b2 >= 0) {
                Int b = bmp::sqrt(b2);
                if (b*b == b2) {
                    Int f1 = a - b;
                    if (f1 > 1 && f1 < n && (n % f1) == 0) { r.status = E8Status::HIT; r.factor = f1; r.t_ms = ms_since(t0); return r; }
                }
            }
            a += 1;
        }

        // Pollard p-1 stage 1
        static std::vector<uint32_t> primes;
        if (primes.empty() || primes.back() < p1_B1) {
            uint32_t N = p1_B1;
            std::vector<bool> is_p(N+1,true);
            is_p[0]=is_p[1]=false;
            for (uint32_t i=2;i*i<=N;++i) if (is_p[i]) for (uint32_t j=i*i;j<=N;j+=i) is_p[j]=false;
            primes.clear();
            for (uint32_t i=2;i<=N;++i) if (is_p[i]) primes.push_back(i);
        }
        Int A = 2 % n;
        for (uint32_t p : primes) {
            if (dl.exceeded()) { r.status = E8Status::TIMEOUT; r.t_ms = ms_since(t0); return r; }
            uint32_t pk = p;
            while (pk <= p1_B1 / p) pk *= p;
            A = pow_mod<Int, Wider>(A, Int(pk), n);
        }
        Int g = gcd_int(A > 1 ? (A-1) : Int(0), n);
        if (g > 1 && g < n) { r.status = E8Status::HIT; r.factor = g; r.t_ms = ms_since(t0); return r; }

        r.status = E8Status::MISS;
        r.t_ms = ms_since(t0);
        return r;
    } catch (...) {
        r.status = E8Status::ERROR;
        r.t_ms = ms_since(t0);
        return r;
    }
}

// -------- dataset --------
enum class Kind { PRIME, SEMI_BALANCED, SEMI_UNBALANCED, TRIPLE, QUADRUPLE, BAMBERG_QUAD };

static inline const char* kind_name(Kind k) {
    switch (k) {
        case Kind::PRIME: return "prime";
        case Kind::SEMI_BALANCED: return "semi_bal";
        case Kind::SEMI_UNBALANCED: return "semi_unbal";
        case Kind::TRIPLE: return "triple";
        case Kind::QUADRUPLE: return "quadruple";
        case Kind::BAMBERG_QUAD: return "bamberg_quad";
        default: return "?";
    }
}
static inline Kind kind_from(const std::string& s) {
    if (s=="prime") return Kind::PRIME;
    if (s=="semi_bal") return Kind::SEMI_BALANCED;
    if (s=="semi_unbal") return Kind::SEMI_UNBALANCED;
    if (s=="triple") return Kind::TRIPLE;
    if (s=="quadruple") return Kind::QUADRUPLE;
    return Kind::BAMBERG_QUAD;
}

template<class Int>
struct Sample { Kind kind; unsigned bits; Int n; std::vector<Int> truth; };

template<class Int>
static Sample<Int> make_sample(Kind kind, unsigned bits, int mr_rounds, std::mt19937_64& rng) {
    Sample<Int> s; s.kind=kind; s.bits=bits;
    auto prod_k = [&](int k, const int* mod12) {
        std::vector<Int> facs; facs.reserve(k);
        unsigned each = std::max(24u, bits/(unsigned)k);
        for (int i=0;i<k;++i) facs.push_back(random_probable_prime<Int>(each, mr_rounds, rng, mod12?mod12[i]:-1));
        std::sort(facs.begin(), facs.end());
        s.n = 1; for (auto& f: facs) s.n *= f;
        s.truth = facs;
    };
    if (kind==Kind::PRIME) { s.n=random_probable_prime<Int>(bits,mr_rounds,rng,-1); s.truth={s.n}; return s; }
    if (kind==Kind::SEMI_BALANCED) {
        unsigned a=std::max(32u,bits/2), b=std::max(32u,bits-a);
        Int p=random_probable_prime<Int>(a,mr_rounds,rng,-1), q=random_probable_prime<Int>(b,mr_rounds,rng,-1);
        if (p>q) std::swap(p,q);
        s.truth={p,q}; s.n=p*q; return s;
    }
    if (kind==Kind::SEMI_UNBALANCED) {
        unsigned small=40u; if (small>=bits-24) small=bits/3;
        unsigned big=std::max(32u,bits-small);
        Int p=random_probable_prime<Int>(small,mr_rounds,rng,-1), q=random_probable_prime<Int>(big,mr_rounds,rng,-1);
        if (p>q) std::swap(p,q);
        s.truth={p,q}; s.n=p*q; return s;
    }
    if (kind==Kind::TRIPLE) { prod_k(3,nullptr); return s; }
    if (kind==Kind::QUADRUPLE) { prod_k(4,nullptr); return s; }
    int want[4]={1,5,7,11}; prod_k(4,want); return s;
}

template<class Int>
static bool save_dataset_csv(const std::string& path, const std::vector<Sample<Int>>& ds) {
    std::ofstream f(path);
    if (!f) return false;
    f << "kind,bits,n,truth\n";
    for (const auto& s: ds) {
        f << kind_name(s.kind) << "," << s.bits << "," << to_dec(s.n) << "," << join_factors_full(s.truth) << "\n";
    }
    return true;
}
template<class Int>
static bool load_dataset_csv(const std::string& path, std::vector<Sample<Int>>& ds_out) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    if (!std::getline(f,line)) return false;
    while (std::getline(f,line)) {
        line=trim(line); if (line.empty()) continue;
        size_t p1=line.find(','); if (p1==std::string::npos) continue;
        size_t p2=line.find(',',p1+1); if (p2==std::string::npos) continue;
        size_t p3=line.find(',',p2+1); if (p3==std::string::npos) continue;
        std::string k=trim(line.substr(0,p1));
        std::string bits_s=trim(line.substr(p1+1,p2-p1-1));
        std::string n_s=trim(line.substr(p2+1,p3-p2-1));
        std::string t_s=trim(line.substr(p3+1));
        Sample<Int> s; s.kind=kind_from(k); s.bits=(unsigned)std::stoul(bits_s);
        bmp::cpp_int nc(n_s); s.n=Int(nc);
        s.truth.clear();
        if (!t_s.empty()) {
            std::stringstream ss(t_s);
            std::string part;
            while (std::getline(ss,part,'*')) {
                part=trim(part); if (part.empty()) continue;
                bmp::cpp_int fc(part); s.truth.push_back(Int(fc));
            }
            std::sort(s.truth.begin(), s.truth.end());
        }
        ds_out.push_back(std::move(s));
    }
    return true;
}

// -------- args --------
struct Args {
    bool generate=true;
    unsigned bits=128;
    int count=30;
    uint64_t seed=42;
    int mr_rounds=16;
    bool include_bamberg=true;

    std::string save_dataset;
    std::string load_dataset;

    bool show_full=false;
    int short_digits=10;

    uint64_t hard_ms=300000;
    uint64_t soft_ms=120000;

    int rho_restarts=1024;
    uint32_t rho_m=512;

    bool e8_on=false;
    std::string out_csv="results.csv";
};

static Args parse(int argc, char** argv) {
    Args a;
    for (int i=1;i<argc;++i) {
        std::string k=argv[i];
        auto need=[&](const std::string& name)->std::string{
            if (i+1>=argc) throw std::runtime_error("missing value for "+name);
            return argv[++i];
        };
        if (k=="--bits") a.bits=(unsigned)std::stoul(need(k));
        else if (k=="--count") a.count=std::stoi(need(k));
        else if (k=="--seed") a.seed=std::stoull(need(k));
        else if (k=="--mr-rounds") a.mr_rounds=std::stoi(need(k));
        else if (k=="--bamberg") a.include_bamberg=(need(k)=="on");
        else if (k=="--save-dataset") a.save_dataset=need(k);
        else if (k=="--load-dataset") { a.load_dataset=need(k); a.generate=false; }
        else if (k=="--show-full") a.show_full=true;
        else if (k=="--short-digits") a.short_digits=std::stoi(need(k));
        else if (k=="--hard-ms") a.hard_ms=(uint64_t)std::stoull(need(k));
        else if (k=="--soft-ms") a.soft_ms=(uint64_t)std::stoull(need(k));
        else if (k=="--rho-restarts") a.rho_restarts=std::stoi(need(k));
        else if (k=="--rho-m") a.rho_m=(uint32_t)std::stoul(need(k));
        else if (k=="--e8") a.e8_on=(need(k)=="on");
        else if (k=="--out") a.out_csv=need(k);
        else if (k=="--help" || k=="-h") {
            std::cout <<
              "Usage: ./ring_bench_e8_compare_v5b [options]\n"
              "Dataset:\n"
              "  --bits N --count N --seed S --bamberg on|off\n"
              "  --save-dataset file.csv  (also runs bench)\n"
              "  --load-dataset file.csv  (bench existing dataset)\n"
              "Bench:\n"
              "  --mr-rounds R --soft-ms T --hard-ms T\n"
              "  --rho-restarts N --rho-m M\n"
              "  --e8 on|off\n"
              "Output:\n"
              "  --out results.csv --short-digits D --show-full\n";
            std::exit(0);
        } else {
            throw std::runtime_error("unknown arg: " + k);
        }
    }
    if (a.count<=0) throw std::runtime_error("--count must be >0");
    if (a.short_digits<1) a.short_digits=1;
    if (a.hard_ms && a.soft_ms > a.hard_ms) a.soft_ms = a.hard_ms;
    return a;
}

// quantile
static inline double quantile(std::vector<double> v, double q) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    double pos = q * (v.size()-1);
    size_t i = (size_t)pos;
    double frac = pos - i;
    if (i+1 < v.size()) return v[i]*(1-frac) + v[i+1]*frac;
    return v.back();
}

template<class Int, class Wider>
static int run(const Args& args) {
    std::mt19937_64 rng(args.seed);
    std::vector<Sample<Int>> ds;

    if (!args.generate) {
        if (!load_dataset_csv<Int>(args.load_dataset, ds)) {
            std::cerr << "ERROR: could not load dataset: " << args.load_dataset << "\n";
            return 2;
        }
    } else {
        std::vector<Kind> kinds = {Kind::PRIME,Kind::SEMI_BALANCED,Kind::SEMI_UNBALANCED,Kind::TRIPLE,Kind::QUADRUPLE};
        if (args.include_bamberg) kinds.push_back(Kind::BAMBERG_QUAD);
        for (auto k : kinds) for (int i=0;i<args.count;++i) ds.push_back(make_sample<Int>(k,args.bits,args.mr_rounds,rng));
        if (!args.save_dataset.empty()) {
            if (save_dataset_csv<Int>(args.save_dataset, ds)) std::cout << "Saved dataset: " << args.save_dataset << " (rows="<<ds.size()<<")\n";
            else std::cerr << "WARNING: could not save dataset\n";
        }
    }

    std::ofstream out(args.out_csv);
    if (!out) { std::cerr << "ERROR: cannot write " << args.out_csv << "\n"; return 3; }
    out << "i,kind,bits,n,truth,MR_ms,MR_is_prime,E8_ms,E8_status,E8_factor,RHO_ms,RHO_status,RHO_factors,ok\n";

    std::cout << "\n=== FULL FACTORIZATION BENCH (B) ===\n"
              << "samples="<<ds.size()<<" bits="<<args.bits<<" mr_rounds="<<args.mr_rounds
              << " soft_ms="<<args.soft_ms<<" hard_ms="<<args.hard_ms
              << " e8="<<(args.e8_on?"on":"off")<<" out="<<args.out_csv<<"\n\n";

    std::cout << "i  kind         MR(ms)  E8(ms)  RHO(ms)  ok  n  truth  found\n"
              << "-----------------------------------------------------------------------------------------------------\n";

    struct Agg{int total=0,ok=0,e8_hit=0,timeouts=0; std::vector<double> mr,e8,rho,tot;};
    std::map<std::string,Agg> agg;

    for (size_t i=0;i<ds.size();++i) {
        const auto& s = ds[i];
        Deadline dl_hard(args.hard_ms);

        auto t_mr0=Clock::now();
        bool mr_prime = is_probable_prime(s.n,args.mr_rounds,rng);
        double mr_ms = ms_since(t_mr0);

        std::vector<Int> truth=s.truth; std::sort(truth.begin(),truth.end());

        if (mr_prime) {
            bool ok = (truth.size()==1 && truth[0]==s.n);
            out << (i+1)<<","<<kind_name(s.kind)<<","<<s.bits<<","<<to_dec(s.n)<<","<<join_factors_full(truth)<<","
                << std::fixed<<std::setprecision(6)<<mr_ms<<",1,"
                << "0.0,DISABLED,,0.0,SKIP,,"
                << (ok?"1":"0") << "\n";

            std::cout<<std::setw(2)<<(i+1)<<"  "<<std::setw(11)<<kind_name(s.kind)<<"  "
                     <<std::fixed<<std::setprecision(3)<<std::setw(6)<<mr_ms<<"  "<<std::setw(6)<<0.0<<"  "<<std::setw(7)<<0.0<<"  "
                     <<(ok?"Y ":"N ")<<short_dec(s.n,args.short_digits,args.show_full)<<"  "
                     <<join_factors_short(truth,args.short_digits,args.show_full)<<"  "
                     <<"PRIME:"<<short_dec(s.n,args.short_digits,args.show_full)<<"\n";

            auto& A=agg[std::string(kind_name(s.kind))];
            A.total++; A.ok += ok?1:0;
            A.mr.push_back(mr_ms); A.e8.push_back(0.0); A.rho.push_back(0.0); A.tot.push_back(mr_ms);
            continue;
        }

        // E8 factor-finder
        double e8_ms=0.0; std::string e8_status="DISABLED"; Int e8_factor=0; bool e8_hit=false;
        if (args.e8_on) {
            auto e8r = e8_find_factor<Int,Wider>(s.n,rng,dl_hard);
            e8_ms = e8r.t_ms;
            e8_status = e8_status_name(e8r.status);
            e8_factor = e8r.factor;
            e8_hit = (e8r.status==E8Status::HIT && e8_factor>1 && e8_factor<s.n && (s.n % e8_factor)==0);
        }

        // Complete with RHO (soft budget)
        Deadline dl_soft(args.soft_ms ? std::min(args.soft_ms,args.hard_ms) : args.hard_ms);
        double rho_ms=0.0; std::string rho_status="SKIP"; std::vector<Int> found; bool rho_to=false;

        auto t_rho0=Clock::now();
        if (e8_hit) {
            std::vector<Int> a,b;
            bool to1=false,to2=false;
            bool ok1 = rho_factorize<Int,Wider>(e8_factor,args.mr_rounds,rng,a,args.rho_restarts,args.rho_m,dl_soft,to1);
            bool ok2 = rho_factorize<Int,Wider>(s.n/e8_factor,args.mr_rounds,rng,b,args.rho_restarts,args.rho_m,dl_soft,to2);
            rho_to = to1||to2;
            if (!rho_to && ok1 && ok2) {
                found.insert(found.end(),a.begin(),a.end());
                found.insert(found.end(),b.begin(),b.end());
                std::sort(found.begin(),found.end());
                rho_status="OK";
            } else rho_status = rho_to ? "TIMEOUT" : "FAIL";
        } else {
            bool okr = rho_factorize<Int,Wider>(s.n,args.mr_rounds,rng,found,args.rho_restarts,args.rho_m,dl_soft,rho_to);
            std::sort(found.begin(),found.end());
            rho_status = rho_to ? "TIMEOUT" : (okr ? "OK" : "FAIL");
        }
        rho_ms = ms_since(t_rho0);

        bool ok = (rho_status=="OK" && found==truth);

        out << (i+1)<<","<<kind_name(s.kind)<<","<<s.bits<<","<<to_dec(s.n)<<","<<join_factors_full(truth)<<","
            << std::fixed<<std::setprecision(6)<<mr_ms<<",0,"
            << e8_ms<<","<<e8_status<<","<<(e8_hit?to_dec(e8_factor):"")<<","
            << rho_ms<<","<<rho_status<<","<<join_factors_full(found)<<","
            << (ok?"1":"0") << "\n";

        std::string found_s = ok ? (std::string(e8_hit?"E8+RHO:":"RHO:")+join_factors_short(found,args.short_digits,args.show_full)) : "-";

        std::cout<<std::setw(2)<<(i+1)<<"  "<<std::setw(11)<<kind_name(s.kind)<<"  "
                 <<std::fixed<<std::setprecision(3)<<std::setw(6)<<mr_ms<<"  "<<std::setw(6)<<e8_ms<<"  "<<std::setw(7)<<rho_ms<<"  "
                 <<(ok?"Y ":"N ")<<short_dec(s.n,args.short_digits,args.show_full)<<"  "
                 <<join_factors_short(truth,args.short_digits,args.show_full)<<"  "<<found_s<<"\n";

        auto& A=agg[std::string(kind_name(s.kind))];
        A.total++; A.ok += ok?1:0; A.e8_hit += e8_hit?1:0;
        if (rho_status=="TIMEOUT") A.timeouts++;
        A.mr.push_back(mr_ms); A.e8.push_back(e8_ms); A.rho.push_back(rho_ms); A.tot.push_back(mr_ms+e8_ms+rho_ms);
    }

    std::cout << "\nSummary (per kind):\n";
    std::cout << "kind         n   ok  e8_hit  timeouts   MR_med  E8_med  RHO_med  TOT_med  RHO_p90  TOT_p90\n";
    std::cout << "-----------------------------------------------------------------------------------------------------------\n";
    for (auto& kv: agg) {
        auto& A=kv.second;
        std::cout<<std::setw(11)<<kv.first<<"  "<<std::setw(3)<<A.total<<"  "<<std::setw(2)<<A.ok<<"  "<<std::setw(6)<<A.e8_hit<<"  "<<std::setw(8)<<A.timeouts<<"  "
                 <<std::fixed<<std::setprecision(2)
                 <<std::setw(6)<<quantile(A.mr,0.5)<<"  "<<std::setw(6)<<quantile(A.e8,0.5)<<"  "<<std::setw(7)<<quantile(A.rho,0.5)<<"  "<<std::setw(7)<<quantile(A.tot,0.5)<<"  "
                 <<std::setw(7)<<quantile(A.rho,0.9)<<"  "<<std::setw(7)<<quantile(A.tot,0.9)<<"\n";
    }

    std::cout << "\nWrote results CSV: " << args.out_csv << "\n";
    return 0;
}

int main(int argc, char** argv) {
    try {
        Args args = parse(argc, argv);
        if (args.bits <= 128) {
            using Int = bmp::uint128_t;
            using Wider = bmp::uint256_t;
            return run<Int,Wider>(args);
        } else if (args.bits <= 256) {
            using Int = bmp::uint256_t;
            using Wider = bmp::uint512_t;
            return run<Int,Wider>(args);
        } else {
            using Int = bmp::cpp_int;
            using Wider = bmp::cpp_int;
            return run<Int,Wider>(args);
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\nTry: --help\n";
        return 1;
    }
}
