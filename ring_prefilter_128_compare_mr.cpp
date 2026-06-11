
/*
 * ring_prefilter_128_compare_mr.cpp
 *
 * 128-bit capable (and beyond) via boost::multiprecision::cpp_int.
 * Correctness-first reference implementation:
 *  - Ring-prefilter factor heuristic in Z[i] and Z[ω] (Eisenstein integers)
 *  - Optional Shadow-Shift for C-like class
 *  - Miller-Rabin primality test for big integers (baseline)
 *  - Speed comparison: MR primality test vs Prefilter attempt on the same dataset
 *
 * Build (macOS/Linux):
 *   g++ -std=c++17 -O3 -march=native -o ring_prefilter_128_compare_mr ring_prefilter_128_compare_mr.cpp
 */

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <optional>
#include <queue>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

using GI = std::pair<cpp_int, cpp_int>; // (a,b) -> a + b i
using EI = std::pair<cpp_int, cpp_int>; // (a,b) -> a + b ω

// ============================================================
// 0) Helpers
// ============================================================

static inline cpp_int abs_cpp(const cpp_int& x) { return x < 0 ? -x : x; }

static inline cpp_int gcd_cpp(cpp_int a, cpp_int b) {
    a = abs_cpp(a); b = abs_cpp(b);
    while (b != 0) {
        cpp_int t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static inline cpp_int isqrt_cpp(const cpp_int& n) {
    if (n <= 0) return 0;
    cpp_int x = n;
    cpp_int y = (x + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
}

static inline cpp_int round_nearest_rational(const cpp_int& num, const cpp_int& den) {
    // Round num/den to nearest integer; ties away from zero.
    if (den <= 0) throw std::runtime_error("den must be positive");
    if (num >= 0) {
        cpp_int q = num / den;
        cpp_int r = num % den;
        cpp_int two_r = r * 2;
        if (two_r > den) return q + 1;
        if (two_r < den) return q;
        return q + 1;
    } else {
        return -round_nearest_rational(-num, den);
    }
}

// ============================================================
// 1) Jacobi + A/B/C/E-like classification
// ============================================================

static inline int jacobi(cpp_int a, cpp_int n) {
    a %= n;
    if (a < 0) a += n;
    if (n <= 0 || (n & 1) == 0) return 0;
    int result = 1;
    while (a != 0) {
        while ((a & 1) == 0) {
            a >>= 1;
            cpp_int r = n & 7; // n % 8
            if (r == 3 || r == 5) result = -result;
        }
        std::swap(a, n);
        if ( (a & 3) == 3 && (n & 3) == 3 ) result = -result; // mod 4 == 3
        a %= n;
        if (a < 0) a += n;
    }
    return (n == 1) ? result : 0;
}

static inline std::pair<int,int> v4_signature(const cpp_int& n) {
    return { jacobi(cpp_int(-1), n), jacobi(cpp_int(-3), n) };
}

static inline char abc_like_class(const cpp_int& n) {
    auto [chi4, chi3] = v4_signature(n);
    if (chi4 == -1 && chi3 == -1) return 'C';
    if (chi4 ==  1 && chi3 ==  1) return 'E';
    if (chi4 ==  1 && chi3 == -1) return 'A';
    if (chi4 == -1 && chi3 ==  1) return 'B';
    return 'U';
}

// ============================================================
// 2) Directions cache (int offsets)
// ============================================================

struct PairHash {
    std::size_t operator()(const std::pair<int,int>& p) const noexcept {
        return (std::uint64_t(std::uint32_t(p.first)) << 32) ^ std::uint32_t(p.second);
    }
};

static std::unordered_map<std::pair<int,int>, std::vector<std::pair<int,int>>, PairHash> DIRS_CACHE;

static std::vector<std::pair<int,int>> surrogate_dirs(int m, std::uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> dist(-200, 200);
    std::vector<std::pair<int,int>> dirs;
    dirs.reserve(m);

    auto gcd_i = [](int a, int b) -> int {
        a = std::abs(a); b = std::abs(b);
        while (b) { int t = a % b; a = b; b = t; }
        return a ? a : 1;
    };

    std::unordered_map<long long, bool> seen;
    seen.reserve(m * 4);

    auto pack = [](int x, int y) -> long long {
        return ( (long long)x << 32 ) ^ (unsigned int)y;
    };

    while ((int)dirs.size() < m) {
        int x = dist(rng), y = dist(rng);
        if (x == 0 && y == 0) continue;
        int g = gcd_i(x, y);
        x /= g; y /= g;
        long long k  = pack(x, y);
        long long kn = pack(-x, -y);
        if (seen.count(k) || seen.count(kn)) continue;
        seen[k] = true;
        dirs.push_back({x, y});
    }
    return dirs;
}

static const std::vector<std::pair<int,int>>& get_surrogate_dirs_cached(int m, int seed) {
    std::pair<int,int> key = {m, seed};
    auto it = DIRS_CACHE.find(key);
    if (it != DIRS_CACHE.end()) return it->second;
    auto v = surrogate_dirs(m, (std::uint32_t)seed);
    auto res = DIRS_CACHE.emplace(key, std::move(v));
    return res.first->second;
}

// ============================================================
// 3) Scoring
// ============================================================

static inline int mod_int(const cpp_int& x, int M) {
    cpp_int r = x % M;
    if (r < 0) r += M;
    return (int)r;
}

static inline double score_phase_gaussian_cls(const GI& g, char cls) {
    const int M = 97;
    int tx = 0, ty = 0;
    if (cls == 'A') tx = ty = M/3;
    else if (cls == 'B') { tx = 2*M/3; ty = M/3; }
    else if (cls == 'C') { tx = M/2; ty = 2*M/3; }
    int x = mod_int(g.first, M);
    int y = mod_int(g.second, M);
    int dx = std::min((x - tx + M) % M, (tx - x + M) % M);
    int dy = std::min((y - ty + M) % M, (ty - y + M) % M);
    return -(dx*dx + dy*dy);
}

static inline double score_phase_eisenstein_cls(const EI& g, char cls) {
    const int M = 97;
    int ta = 0, tb = 0;
    if (cls == 'A') ta = tb = M/3;
    else if (cls == 'B') { ta = 2*M/3; tb = M/3; }
    else if (cls == 'C') { ta = M/2; tb = 2*M/3; }
    int a = mod_int(g.first, M);
    int b = mod_int(g.second, M);
    int da = std::min((a - ta + M) % M, (ta - a + M) % M);
    int db = std::min((b - tb + M) % M, (tb - b + M) % M);
    return -(da*da + db*db);
}

static inline double score_l1_gaussian(const GI& g, const cpp_int& base) {
    return -double(abs_cpp(g.first - base) + abs_cpp(g.second));
}

static inline double score_l1_eisenstein(const EI& g, const cpp_int& base) {
    return -double(abs_cpp(g.first - base) + abs_cpp(g.second));
}

// ============================================================
// 4) Z[i]
// ============================================================

static inline GI gi_mul(const GI& x, const GI& y) {
    return { x.first*y.first - x.second*y.second,
             x.first*y.second + x.second*y.first };
}

static inline GI gi_conj(const GI& x) { return {x.first, -x.second}; }
static inline GI gi_sub(const GI& x, const GI& y) { return {x.first - y.first, x.second - y.second}; }
static inline cpp_int gi_norm(const GI& x) { return x.first*x.first + x.second*x.second; }
static inline bool gi_is_zero(const GI& x) { return x.first == 0 && x.second == 0; }

static inline std::pair<GI,GI> gi_divmod(const GI& alpha, const GI& beta) {
    if (gi_is_zero(beta)) throw std::runtime_error("Z[i] division by zero");
    cpp_int N = gi_norm(beta);
    GI z = gi_mul(alpha, gi_conj(beta));
    cpp_int qa = round_nearest_rational(z.first, N);
    cpp_int qb = round_nearest_rational(z.second, N);
    GI q = {qa, qb};
    GI r = gi_sub(alpha, gi_mul(q, beta));
    return {q, r};
}

static inline GI gi_gcd(GI x, GI y) {
    while (!gi_is_zero(y)) {
        auto qr = gi_divmod(x, y);
        x = y;
        y = qr.second;
    }
    return x;
}

// ============================================================
// 5) Z[ω] with exact rounding (no floats)
// ============================================================

static inline EI ei_mul(const EI& x, const EI& y) {
    return { x.first*y.first - x.second*y.second,
             x.first*y.second + x.second*y.first - x.second*y.second };
}

static inline EI ei_conj(const EI& x) { return { x.first - x.second, -x.second }; }
static inline EI ei_sub(const EI& x, const EI& y) { return { x.first - y.first, x.second - y.second }; }

static inline cpp_int ei_norm(const EI& x) {
    return x.first*x.first - x.first*x.second + x.second*x.second;
}

static inline bool ei_is_zero(const EI& x) { return x.first == 0 && x.second == 0; }

static inline EI nearest_eisenstein_from_t_over_N(const EI& t, const cpp_int& N) {
    const cpp_int& A = t.first;
    const cpp_int& B = t.second;
    cpp_int ru = round_nearest_rational(A, N);
    cpp_int rv = round_nearest_rational(B, N);
    cpp_int C  = -(A + B);
    cpp_int rw = round_nearest_rational(C, N);

    cpp_int du = abs_cpp(ru*N - A);
    cpp_int dv = abs_cpp(rv*N - B);
    cpp_int dw = abs_cpp(rw*N - C);

    if (du > dv && du > dw) ru = -rv - rw;
    else if (dv > dw) rv = -ru - rw;
    else rw = -ru - rv;

    return {ru, rv};
}

static int EI_DIV_CALLS = 0;
static int EI_NEIGHBOR_CALLS = 0;

static inline std::pair<EI,EI> ei_divmod_adaptive(const EI& alpha, const EI& beta,
                                                  bool neighbor_safety = true,
                                                  bool adaptive = true) {
    EI_DIV_CALLS++;
    if (ei_is_zero(beta)) throw std::runtime_error("Z[ω] division by zero");
    cpp_int Nbeta = ei_norm(beta);

    EI t = ei_mul(alpha, ei_conj(beta));
    EI q0 = nearest_eisenstein_from_t_over_N(t, Nbeta);
    EI r0 = ei_sub(alpha, ei_mul(q0, beta));

    if (!neighbor_safety) return {q0, r0};

    if (adaptive) {
        // threshold 0.92 exactly as 92/100
        cpp_int nr0 = ei_norm(r0);
        if (nr0 * 100 <= Nbeta * 92) return {q0, r0};
    }

    EI_NEIGHBOR_CALLS++;
    static const std::pair<int,int> neigh[] = {{0,0},{1,0},{-1,0},{0,1},{0,-1},{1,1},{-1,-1}};

    EI best_q = q0;
    EI best_r = r0;
    cpp_int best_n = ei_norm(r0);

    for (auto [da, db] : neigh) {
        EI q = { q0.first + da, q0.second + db };
        EI r = ei_sub(alpha, ei_mul(q, beta));
        cpp_int nr = ei_norm(r);
        if (nr < best_n) { best_n = nr; best_q = q; best_r = r; }
    }
    return {best_q, best_r};
}

static inline EI ei_gcd(EI x, EI y, bool neighbor_safety = true) {
    while (!ei_is_zero(y)) {
        auto qr = ei_divmod_adaptive(x, y, neighbor_safety, true);
        x = y;
        y = qr.second;
    }
    return x;
}

// ============================================================
// 6) Candidate generation with streaming Top-B heap
// ============================================================

struct CandidateSpec {
    int T = 40;
    int S = 40;
    int base_samples = 8;
    int dir_count = 120;
    std::vector<int> lambdas = {1,2,3,5,8};
    int dir_seed_offset = 12345;
};

template<typename Cand>
struct HeapItem {
    double score;
    std::uint64_t counter;
    Cand cand;
};

template<typename Cand>
struct HeapComp {
    bool operator()(const HeapItem<Cand>& a, const HeapItem<Cand>& b) const {
        if (a.score != b.score) return a.score > b.score;
        return a.counter > b.counter;
    }
};

static inline std::vector<GI> topB_gaussian_stream(uint32_t seed, const cpp_int& base, const CandidateSpec& spec,
                                                   const std::function<double(const GI&)>& score_fn, int B) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> distT(-spec.T, spec.T), distS(-spec.S, spec.S);
    const auto& dirs = get_surrogate_dirs_cached(spec.dir_count, (int)seed + spec.dir_seed_offset);

    std::priority_queue<HeapItem<GI>, std::vector<HeapItem<GI>>, HeapComp<GI>> heap;
    std::uint64_t counter = 0;

    auto push = [&](const GI& cand) {
        double s = score_fn(cand);
        HeapItem<GI> item{s, counter++, cand};
        if ((int)heap.size() < B) heap.push(item);
        else if (s > heap.top().score) { heap.pop(); heap.push(item); }
    };

    for (int i = 0; i < spec.base_samples; ++i) {
        int t = distT(rng), s = distS(rng);
        GI g0 = { base + t, cpp_int(s) };
        push(g0);
        for (auto [dx, dy] : dirs) {
            for (int lam : spec.lambdas) {
                GI g = { g0.first + cpp_int(lam*dx), g0.second + cpp_int(lam*dy) };
                push(g);
            }
        }
    }

    std::vector<GI> out;
    out.reserve(heap.size());
    while (!heap.empty()) { out.push_back(heap.top().cand); heap.pop(); }
    std::reverse(out.begin(), out.end());
    return out;
}

static inline std::vector<EI> topB_eisenstein_stream(uint32_t seed, const cpp_int& base, const CandidateSpec& spec,
                                                     const std::function<double(const EI&)>& score_fn, int B) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> distT(-spec.T, spec.T), distS(-spec.S, spec.S);
    const auto& dirs = get_surrogate_dirs_cached(spec.dir_count, (int)seed + spec.dir_seed_offset + 777);

    std::priority_queue<HeapItem<EI>, std::vector<HeapItem<EI>>, HeapComp<EI>> heap;
    std::uint64_t counter = 0;

    auto push = [&](const EI& cand) {
        double s = score_fn(cand);
        HeapItem<EI> item{s, counter++, cand};
        if ((int)heap.size() < B) heap.push(item);
        else if (s > heap.top().score) { heap.pop(); heap.push(item); }
    };

    for (int i = 0; i < spec.base_samples; ++i) {
        int t = distT(rng), s = distS(rng);
        EI g0 = { base + t, cpp_int(s) };
        push(g0);
        for (auto [da, db] : dirs) {
            for (int lam : spec.lambdas) {
                EI g = { g0.first + cpp_int(lam*da), g0.second + cpp_int(lam*db) };
                push(g);
            }
        }
    }

    std::vector<EI> out;
    out.reserve(heap.size());
    while (!heap.empty()) { out.push_back(heap.top().cand); heap.pop(); }
    std::reverse(out.begin(), out.end());
    return out;
}

// ============================================================
// 7) Prefilter pipeline
// ============================================================

struct PrefilterResult {
    bool success = false;
    std::optional<cpp_int> factor;
    std::optional<int> rank_hit;
    int tested = 0;
    double t_prefilter = 0.0;
    std::string ring_used;
};

static inline std::optional<cpp_int> try_factor_gaussian(const cpp_int& n, const GI& g) {
    GI d = gi_gcd({n, 0}, g);
    cpp_int f = gcd_cpp(gi_norm(d), n);
    if (f > 1 && f < n) return f;
    return std::nullopt;
}

static inline std::optional<cpp_int> try_factor_eisenstein(const cpp_int& n, const EI& g, bool neighbor_safety) {
    EI d = ei_gcd({n, 0}, g, neighbor_safety);
    cpp_int f = gcd_cpp(ei_norm(d), n);
    if (f > 1 && f < n) return f;
    return std::nullopt;
}

static inline PrefilterResult prefilter_gaussian_fast(const cpp_int& n, const CandidateSpec& spec,
                                                      int B, const std::string& score_mode, uint32_t seed) {
    auto t0 = std::chrono::steady_clock::now();
    PrefilterResult out;
    out.ring_used = "gaussian";

    char cls = abc_like_class(n);
    cpp_int base = isqrt_cpp(n);

    std::function<double(const GI&)> score_fn;
    if (score_mode == "phase") score_fn = [cls](const GI& g){ return score_phase_gaussian_cls(g, cls); };
    else score_fn = [base](const GI& g){ return score_l1_gaussian(g, base); };

    auto ordered = topB_gaussian_stream(seed, base, spec, score_fn, B);

    for (std::size_t i = 0; i < ordered.size(); ++i) {
        out.tested++;
        auto f = try_factor_gaussian(n, ordered[i]);
        if (f) { out.success = true; out.factor = *f; out.rank_hit = (int)i + 1; break; }
    }
    out.t_prefilter = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return out;
}

static inline PrefilterResult prefilter_eisenstein_fast(const cpp_int& n, const CandidateSpec& spec,
                                                        int B, const std::string& score_mode, uint32_t seed,
                                                        bool neighbor_safety) {
    auto t0 = std::chrono::steady_clock::now();
    PrefilterResult out;
    out.ring_used = "eisenstein";

    char cls = abc_like_class(n);
    cpp_int base = isqrt_cpp(n);

    std::function<double(const EI&)> score_fn;
    if (score_mode == "phase") score_fn = [cls](const EI& g){ return score_phase_eisenstein_cls(g, cls); };
    else score_fn = [base](const EI& g){ return score_l1_eisenstein(g, base); };

    auto ordered = topB_eisenstein_stream(seed, base, spec, score_fn, B);

    for (std::size_t i = 0; i < ordered.size(); ++i) {
        out.tested++;
        auto f = try_factor_eisenstein(n, ordered[i], neighbor_safety);
        if (f) { out.success = true; out.factor = *f; out.rank_hit = (int)i + 1; break; }
    }
    out.t_prefilter = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return out;
}

static inline PrefilterResult prefilter_pipeline(const cpp_int& n, const CandidateSpec& spec, int B,
                                                 const std::string& score_mode,
                                                 uint32_t seed, bool neighbor_safety,
                                                 bool shadow_shift) {
    char cls = abc_like_class(n);

    if (shadow_shift && cls == 'C') {
        for (int k : {5,7}) {
            cpp_int N = n * k;
            auto res = prefilter_pipeline(N, spec, B, score_mode, seed + 1000 + k, neighbor_safety, false);
            if (res.success && res.factor) {
                cpp_int f = gcd_cpp(*res.factor, n);
                if (f > 1 && f < n) { res.factor = f; res.ring_used = "shadow"; return res; }
            }
        }
        PrefilterResult out;
        out.ring_used = "shadow";
        return out;
    }

    if (cls == 'A') return prefilter_gaussian_fast(n, spec, B, score_mode, seed);
    if (cls == 'B') return prefilter_eisenstein_fast(n, spec, B, score_mode, seed, neighbor_safety);

    if (cls == 'E') {
        int Bg = B/2;
        int Be = B - Bg;
        auto r1 = prefilter_gaussian_fast(n, spec, Bg, score_mode, seed + 1);
        if (r1.success) return r1;
        auto r2 = prefilter_eisenstein_fast(n, spec, Be, score_mode, seed + 2, neighbor_safety);
        if (r2.success) { r2.tested += r1.tested; r2.t_prefilter += r1.t_prefilter; return r2; }
        PrefilterResult out;
        out.ring_used = "E-both";
        out.tested = r1.tested + r2.tested;
        out.t_prefilter = r1.t_prefilter + r2.t_prefilter;
        return out;
    }

    return prefilter_gaussian_fast(n, spec, B, score_mode, seed);
}

// ============================================================
// 8) Miller-Rabin for cpp_int
// ============================================================

static inline cpp_int mod_mul(const cpp_int& a, const cpp_int& b, const cpp_int& mod) {
    return (a * b) % mod;
}

static inline cpp_int mod_pow(cpp_int a, cpp_int d, const cpp_int& mod) {
    cpp_int r = 1 % mod;
    a %= mod;
    while (d > 0) {
        if ((d & 1) != 0) r = mod_mul(r, a, mod);
        a = mod_mul(a, a, mod);
        d >>= 1;
    }
    return r;
}

static inline cpp_int rand_below(const cpp_int& limit, std::mt19937_64& rng) {
    unsigned bits = boost::multiprecision::msb(limit) + 1;
    while (true) {
        cpp_int x = 0;
        for (unsigned i = 0; i < bits; i += 64) { x <<= 64; x += (cpp_int)rng(); }
        if (bits % 64) x &= ((cpp_int(1) << bits) - 1);
        if (x < limit) return x;
    }
}

static inline bool miller_rabin_is_probable_prime(const cpp_int& n, int rounds, std::mt19937_64& rng) {
    if (n < 2) return false;
    static const int small_primes[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int p : small_primes) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }
    if ((n & 1) == 0) return false;

    cpp_int d = n - 1;
    unsigned s = 0;
    while ((d & 1) == 0) { d >>= 1; s++; }

    for (int i = 0; i < rounds; ++i) {
        cpp_int a = 2 + rand_below(n - 3, rng); // [2, n-2]
        cpp_int x = mod_pow(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool witness = true;
        for (unsigned r = 0; r < s - 1; ++r) {
            x = mod_mul(x, x, n);
            if (x == n - 1) { witness = false; break; }
        }
        if (witness) return false;
    }
    return true;
}

// ============================================================
// 9) Dataset: semiprimes
// ============================================================

static inline cpp_int random_bits(unsigned bits, std::mt19937_64& rng) {
    cpp_int x = 0;
    for (unsigned i = 0; i < bits; i += 64) { x <<= 64; x += (cpp_int)rng(); }
    if (bits % 64) x &= ((cpp_int(1) << bits) - 1);
    return x;
}

static inline cpp_int random_odd_with_bits(unsigned bits, std::mt19937_64& rng) {
    cpp_int x = random_bits(bits, rng);
    x |= 1;
    x |= (cpp_int(1) << (bits - 1));
    return x;
}

static inline cpp_int random_probable_prime(unsigned bits, std::mt19937_64& rng, int mr_rounds) {
    while (true) {
        cpp_int x = random_odd_with_bits(bits, rng);
        if (miller_rabin_is_probable_prime(x, mr_rounds, rng)) return x;
    }
}

static inline std::tuple<cpp_int, cpp_int, cpp_int> make_semiprime(unsigned bits, std::mt19937_64& rng, int mr_rounds) {
    cpp_int p = random_probable_prime(bits/2, rng, mr_rounds);
    cpp_int q = random_probable_prime(bits/2, rng, mr_rounds);
    return {p*q, p, q};
}

// ============================================================
// 10) Benchmark
// ============================================================

struct BenchRow {
    int i;
    char cls;
    double t_mr;
    bool mr_probably_prime;
    double t_pref;
    bool pref_ok;
    int tested;
};

static inline void print_rows(const std::vector<BenchRow>& rows, int max_rows = 12) {
    std::cout << "i  cls  t_MR(s)   MR_prime?  t_pref(s)  pref_ok  tested\n";
    std::cout << "------------------------------------------------------\n";
    for (int k = 0; k < (int)rows.size() && k < max_rows; ++k) {
        const auto& r = rows[k];
        std::cout
            << std::setw(2) << r.i << "  "
            << r.cls << "    "
            << std::setw(8) << std::fixed << std::setprecision(5) << r.t_mr << "   "
            << (r.mr_probably_prime ? "Y" : "N") << "         "
            << std::setw(8) << std::fixed << std::setprecision(5) << r.t_pref << "    "
            << (r.pref_ok ? "Y" : "N") << "      "
            << r.tested
            << "\n";
    }
}

static inline unsigned to_uint(const std::string& s, unsigned def) {
    try { return (unsigned)std::stoul(s); } catch(...) { return def; }
}
static inline int to_int(const std::string& s, int def) {
    try { return std::stoi(s); } catch(...) { return def; }
}

int main(int argc, char** argv) {
    int trials = 20;
    unsigned bits = 128;
    int B = 400;
    std::string score = "phase";
    bool neighbor_safety = true;
    bool shadow_shift = true;
    std::uint64_t seed = 1234567;
    int mr_rounds = 12;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--trials" && i+1 < argc) { trials = to_int(argv[++i], trials); continue; }
        if (a == "--bits" && i+1 < argc) { bits = to_uint(argv[++i], bits); continue; }
        if (a == "--B" && i+1 < argc) { B = to_int(argv[++i], B); continue; }
        if (a == "--score" && i+1 < argc) { score = argv[++i]; continue; }
        if (a == "--neighbor-safety") { neighbor_safety = true; continue; }
        if (a == "--no-neighbor-safety") { neighbor_safety = false; continue; }
        if (a == "--shadow") { shadow_shift = true; continue; }
        if (a == "--no-shadow") { shadow_shift = false; continue; }
        if (a == "--seed" && i+1 < argc) { seed = (std::uint64_t)std::stoull(argv[++i]); continue; }
        if (a == "--mr-rounds" && i+1 < argc) { mr_rounds = to_int(argv[++i], mr_rounds); continue; }
    }

    CandidateSpec spec;
    std::mt19937_64 rng(seed);

    std::vector<cpp_int> dataset;
    dataset.reserve(trials);

    for (int i = 0; i < trials; ++i) {
        auto [n, p, q] = make_semiprime(bits, rng, mr_rounds);
        dataset.push_back(n);
    }

    std::vector<BenchRow> rows;
    rows.reserve(trials);

    double sum_mr = 0.0, sum_pref = 0.0;
    int pref_succ = 0;
    long long tested_sum = 0;

    for (int i = 0; i < trials; ++i) {
        const cpp_int& n = dataset[i];
        char cls = abc_like_class(n);

        auto t0 = std::chrono::steady_clock::now();
        bool mr_ok = miller_rabin_is_probable_prime(n, mr_rounds, rng);
        double t_mr = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

        EI_DIV_CALLS = 0;
        EI_NEIGHBOR_CALLS = 0;
        auto t1 = std::chrono::steady_clock::now();
        auto pref = prefilter_pipeline(n, spec, B, score, (uint32_t)rng(), neighbor_safety, shadow_shift);
        double t_pref = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();

        rows.push_back(BenchRow{i+1, cls, t_mr, mr_ok, t_pref, pref.success, pref.tested});
        sum_mr += t_mr;
        sum_pref += t_pref;
        tested_sum += pref.tested;
        if (pref.success) pref_succ++;
    }

    std::cout << "\n=== 128-bit (bigint) benchmark: Prefilter vs Miller-Rabin ===\n";
    std::cout << "trials=" << trials
              << "  bits=" << bits
              << "  B=" << B
              << "  score=" << score
              << "  neighbor_safety=" << (neighbor_safety ? "true" : "false")
              << "  shadow=" << (shadow_shift ? "true" : "false")
              << "  mr_rounds=" << mr_rounds
              << "\n\n";

    print_rows(rows, std::min(12, trials));

    std::cout << "\nSummary:\n";
    std::cout << "  mean_t_MR      = " << std::fixed << std::setprecision(6) << (sum_mr / trials) << " s\n";
    std::cout << "  mean_t_pref    = " << std::fixed << std::setprecision(6) << (sum_pref / trials) << " s\n";
    std::cout << "  prefilter_hit  = " << pref_succ << "/" << trials << "\n";
    std::cout << "  mean_tested    = " << (double)tested_sum / trials << "\n";
    std::cout << "\nNote:\n";
    std::cout << "  - MR is a primality test; Prefilter tries to find a factor (harder problem).\n";
    std::cout << "  - Timings are reported side-by-side as requested.\n\n";

    return 0;
}
