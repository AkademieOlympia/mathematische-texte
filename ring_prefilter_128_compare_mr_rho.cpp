
/*
 * ring_prefilter_128_compare_mr_rho.cpp
 *
 * Erweiterung deines "compare_mr"-Benchmarks:
 *   - nimmt Pollard-Rho (Brent) + MR in die Auswertung auf
 *   - berichtet pro n: t_MR, t_pref, t_rho_direct, t_pref+rho, pref_ok, tested
 *   - Summary: mean_t_MR, mean_t_pref, mean_t_rho_direct, mean_t_pref+rho, prefilter_hit, rho_after_pref, mean_tested
 *
 * Fokus: 128-bit (uint128_t) Benchmark.
 *   - MR und Rho laufen auf boost::multiprecision::uint128_t/uint256_t (schnell)
 *   - Prefilter (Z[i], Z[ω]) bleibt auf cpp_int (Bigint), aber n kommt aus uint128_t.
 *
 * Build:
 *   g++ -std=c++17 -O3 -march=native -o ring_prefilter_128_compare_mr_rho ring_prefilter_128_compare_mr_rho.cpp
 *
 * Run (ähnlich wie vorher):
 *   ./ring_prefilter_128_compare_mr_rho --trials 20 --bits 128 --B 400 --score phase --neighbor-safety --shadow --mr-rounds 12
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
#include <boost/multiprecision/integer.hpp>

using boost::multiprecision::cpp_int;
using boost::multiprecision::uint128_t;
using boost::multiprecision::uint256_t;

// ------------------ small helpers ------------------

static inline cpp_int abs_cpp(const cpp_int& x) { return x < 0 ? -x : x; }
static inline cpp_int gcd_cpp(cpp_int a, cpp_int b) {
    a = abs_cpp(a); b = abs_cpp(b);
    while (b != 0) { cpp_int t = a % b; a = b; b = t; }
    return a;
}
static inline cpp_int isqrt_cpp(const cpp_int& n) {
    if (n <= 0) return 0;
    cpp_int x = n;
    cpp_int y = (x + 1) / 2;
    while (y < x) { x = y; y = (x + n / x) / 2; }
    return x;
}
static inline cpp_int round_nearest_rational(const cpp_int& num, const cpp_int& den) {
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
static inline int mod_int(const cpp_int& x, int M) {
    cpp_int r = x % M;
    if (r < 0) r += M;
    return (int)r;
}

// ------------------ Jacobi + class ------------------

static inline int jacobi(cpp_int a, cpp_int n) {
    a %= n; if (a < 0) a += n;
    if (n <= 0 || (n & 1) == 0) return 0;
    int result = 1;
    while (a != 0) {
        while ((a & 1) == 0) {
            a >>= 1;
            cpp_int r = n & 7;
            if (r == 3 || r == 5) result = -result;
        }
        std::swap(a, n);
        if ((a & 3) == 3 && (n & 3) == 3) result = -result;
        a %= n; if (a < 0) a += n;
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

// ------------------ dirs cache ------------------

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
        return ((long long)x << 32) ^ (unsigned int)y;
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

// ------------------ scoring ------------------

using GI = std::pair<cpp_int, cpp_int>;
using EI = std::pair<cpp_int, cpp_int>;

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
    cpp_int v = abs_cpp(g.first - base) + abs_cpp(g.second);
    return -v.convert_to<double>();
}
static inline double score_l1_eisenstein(const EI& g, const cpp_int& base) {
    cpp_int v = abs_cpp(g.first - base) + abs_cpp(g.second);
    return -v.convert_to<double>();
}

// ------------------ Z[i] ------------------

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

// ------------------ Z[ω] ------------------

static inline EI ei_mul(const EI& x, const EI& y) {
    return { x.first*y.first - x.second*y.second,
             x.first*y.second + x.second*y.first - x.second*y.second };
}
static inline EI ei_conj(const EI& x) { return { x.first - x.second, -x.second }; }
static inline EI ei_sub(const EI& x, const EI& y) { return { x.first - y.first, x.second - y.second }; }
static inline cpp_int ei_norm(const EI& x) { return x.first*x.first - x.first*x.second + x.second*x.second; }
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
                                                  bool neighbor_safety, bool adaptive = true) {
    EI_DIV_CALLS++;
    if (ei_is_zero(beta)) throw std::runtime_error("Z[ω] division by zero");
    cpp_int Nbeta = ei_norm(beta);

    EI t = ei_mul(alpha, ei_conj(beta));
    EI q0 = nearest_eisenstein_from_t_over_N(t, Nbeta);
    EI r0 = ei_sub(alpha, ei_mul(q0, beta));

    if (!neighbor_safety) return {q0, r0};

    if (adaptive) {
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
static inline EI ei_gcd(EI x, EI y, bool neighbor_safety) {
    while (!ei_is_zero(y)) {
        auto qr = ei_divmod_adaptive(x, y, neighbor_safety, true);
        x = y;
        y = qr.second;
    }
    return x;
}

// ------------------ candidates top-B streaming ------------------

struct CandidateSpec {
    int T = 40;
    int S = 40;
    int base_samples = 8;
    int dir_count = 120;
    std::vector<int> lambdas = {1,2,3,5,8};
    int dir_seed_offset = 12345;
};

template<typename Cand>
struct HeapItem { double score; std::uint64_t counter; Cand cand; };

template<typename Cand>
struct HeapComp {
    bool operator()(const HeapItem<Cand>& a, const HeapItem<Cand>& b) const {
        if (a.score != b.score) return a.score > b.score; // min-heap by score
        return a.counter > b.counter;
    }
};

static inline std::vector<GI> topB_gaussian_stream(uint32_t seed, const cpp_int& base, char cls,
                                                   const CandidateSpec& spec, const std::string& score_mode, int B) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> distT(-spec.T, spec.T), distS(-spec.S, spec.S);
    const auto& dirs = get_surrogate_dirs_cached(spec.dir_count, (int)seed + spec.dir_seed_offset);

    std::function<double(const GI&)> score_fn;
    if (score_mode == "phase") score_fn = [cls](const GI& g){ return score_phase_gaussian_cls(g, cls); };
    else score_fn = [base](const GI& g){ return score_l1_gaussian(g, base); };

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

static inline std::vector<EI> topB_eisenstein_stream(uint32_t seed, const cpp_int& base, char cls,
                                                     const CandidateSpec& spec, const std::string& score_mode, int B) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> distT(-spec.T, spec.T), distS(-spec.S, spec.S);
    const auto& dirs = get_surrogate_dirs_cached(spec.dir_count, (int)seed + spec.dir_seed_offset + 777);

    std::function<double(const EI&)> score_fn;
    if (score_mode == "phase") score_fn = [cls](const EI& g){ return score_phase_eisenstein_cls(g, cls); };
    else score_fn = [base](const EI& g){ return score_l1_eisenstein(g, base); };

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

// ------------------ Prefilter pipeline (Shadow zählt tested) ------------------

struct PrefilterResult {
    bool success = false;
    std::optional<cpp_int> factor;
    long long tested = 0;
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

    auto ordered = topB_gaussian_stream(seed, base, cls, spec, score_mode, B);

    for (std::size_t i = 0; i < ordered.size(); ++i) {
        out.tested++;
        auto f = try_factor_gaussian(n, ordered[i]);
        if (f) { out.success = true; out.factor = *f; break; }
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

    auto ordered = topB_eisenstein_stream(seed, base, cls, spec, score_mode, B);

    for (std::size_t i = 0; i < ordered.size(); ++i) {
        out.tested++;
        auto f = try_factor_eisenstein(n, ordered[i], neighbor_safety);
        if (f) { out.success = true; out.factor = *f; break; }
    }
    out.t_prefilter = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return out;
}

static inline PrefilterResult prefilter_pipeline(const cpp_int& n, const CandidateSpec& spec, int B,
                                                 const std::string& score_mode,
                                                 uint32_t seed, bool neighbor_safety,
                                                 bool shadow_shift) {
    PrefilterResult out;
    char cls = abc_like_class(n);

    if (shadow_shift && cls == 'C') {
        out.ring_used = "shadow";
        for (int k : {5,7}) {
            cpp_int N = n * k;
            auto res = prefilter_pipeline(N, spec, B, score_mode, seed + 1000 + (uint32_t)k, neighbor_safety, false);
            out.tested += res.tested;
            out.t_prefilter += res.t_prefilter;
            if (res.success && res.factor) {
                cpp_int f = gcd_cpp(*res.factor, n);
                if (f > 1 && f < n) { out.success = true; out.factor = f; return out; }
            }
        }
        return out;
    }

    if (cls == 'A') return prefilter_gaussian_fast(n, spec, B, score_mode, seed);
    if (cls == 'B') return prefilter_eisenstein_fast(n, spec, B, score_mode, seed, neighbor_safety);
    if (cls == 'E') {
        int Bg = B/2;
        int Be = B - Bg;
        auto r1 = (Bg>0) ? prefilter_gaussian_fast(n, spec, Bg, score_mode, seed+1) : PrefilterResult{};
        if (r1.success) return r1;
        auto r2 = (Be>0) ? prefilter_eisenstein_fast(n, spec, Be, score_mode, seed+2, neighbor_safety) : PrefilterResult{};
        r2.tested += r1.tested;
        r2.t_prefilter += r1.t_prefilter;
        r2.ring_used = "E-both";
        return r2;
    }
    return prefilter_gaussian_fast(n, spec, B, score_mode, seed);
}

// ------------------ MR + Rho (fast u128) ------------------

static inline uint128_t mul_mod_u128(uint128_t a, uint128_t b, uint128_t mod) {
    uint256_t r = uint256_t(a) * uint256_t(b);
    r %= uint256_t(mod);
    return uint128_t(r);
}
static inline uint128_t pow_mod_u128(uint128_t a, uint128_t d, uint128_t mod) {
    uint128_t r = 1 % mod;
    a %= mod;
    while (d > 0) {
        if ((d & 1) != 0) r = mul_mod_u128(r, a, mod);
        a = mul_mod_u128(a, a, mod);
        d >>= 1;
    }
    return r;
}
static inline uint128_t gcd_u128(uint128_t a, uint128_t b) {
    while (b != 0) { uint128_t t = a % b; a = b; b = t; }
    return a;
}

static inline bool mr_is_probable_prime_u128(uint128_t n, int rounds, std::mt19937_64& rng) {
    if (n < 2) return false;
    static const uint64_t small_primes[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (uint64_t p : small_primes) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }
    if ((n & 1) == 0) return false;

    uint128_t d = n - 1;
    unsigned s = 0;
    while ((d & 1) == 0) { d >>= 1; s++; }

    std::uniform_int_distribution<uint64_t> dist64;
    for (int i = 0; i < rounds; ++i) {
        uint128_t a = 2;
        if (n > 4) {
            uint128_t span = n - 3;
            uint128_t x = (uint128_t(dist64(rng)) << 64) | uint128_t(dist64(rng));
            a = 2 + (x % span);
        }
        uint128_t x = pow_mod_u128(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool witness = true;
        for (unsigned r = 0; r < s - 1; ++r) {
            x = mul_mod_u128(x, x, n);
            if (x == n - 1) { witness = false; break; }
        }
        if (witness) return false;
    }
    return true;
}

static inline uint128_t rho_f(uint128_t x, uint128_t c, uint128_t n) {
    return (mul_mod_u128(x, x, n) + c) % n;
}

static inline uint128_t pollard_rho_brent_u128(uint128_t n, std::mt19937_64& rng) {
    if ((n & 1) == 0) return 2;
    if (n % 3 == 0) return 3;
    std::uniform_int_distribution<uint64_t> dist64;

    while (true) {
        uint128_t y = (uint128_t(dist64(rng)) << 64) | uint128_t(dist64(rng));
        y %= (n - 1); y += 1;
        uint128_t c = (uint128_t(dist64(rng)) << 64) | uint128_t(dist64(rng));
        c %= (n - 1); c += 1;

        uint128_t m = 128;
        uint128_t g = 1, r = 1, q = 1;
        uint128_t x = 0, ys = 0;

        while (g == 1) {
            x = y;
            for (uint128_t i = 0; i < r; ++i) y = rho_f(y, c, n);

            uint128_t k = 0;
            while (k < r && g == 1) {
                ys = y;
                uint128_t rk = r - k;
                uint128_t lim = (m < rk ? m : rk);
                for (uint128_t i = 0; i < lim; ++i) {
                    y = rho_f(y, c, n);
                    uint128_t diff = x > y ? x - y : y - x;
                    q = mul_mod_u128(q, diff, n);
                }
                g = gcd_u128(q, n);
                k += lim;
            }
            r <<= 1;
        }

        if (g == n) {
            do {
                ys = rho_f(ys, c, n);
                uint128_t diff = x > ys ? x - ys : ys - x;
                g = gcd_u128(diff, n);
            } while (g == 1);
        }
        if (g != n) return g;
    }
}

static inline void factorize_rho_u128(uint128_t n, int mr_rounds, std::mt19937_64& rng, std::vector<uint128_t>& out) {
    if (n == 1) return;
    if (mr_is_probable_prime_u128(n, mr_rounds, rng)) { out.push_back(n); return; }
    uint128_t d = pollard_rho_brent_u128(n, rng);
    factorize_rho_u128(d, mr_rounds, rng, out);
    factorize_rho_u128(n / d, mr_rounds, rng, out);
}

// ------------------ semiprime dataset (u128) ------------------

static inline uint128_t random_u128_bits(unsigned bits, std::mt19937_64& rng) {
    std::uniform_int_distribution<uint64_t> dist64;
    uint128_t x = (uint128_t(dist64(rng)) << 64) | uint128_t(dist64(rng));
    if (bits < 128) {
        uint128_t mask = (uint128_t(1) << bits) - 1;
        x &= mask;
    }
    return x;
}
static inline uint128_t random_odd_u128(unsigned bits, std::mt19937_64& rng) {
    uint128_t x = random_u128_bits(bits, rng);
    x |= 1;
    x |= (uint128_t(1) << (bits - 1));
    return x;
}
static inline uint128_t random_probable_prime_u128(unsigned bits, std::mt19937_64& rng, int mr_rounds) {
    while (true) {
        uint128_t x = random_odd_u128(bits, rng);
        if (mr_is_probable_prime_u128(x, mr_rounds, rng)) return x;
    }
}
static inline std::tuple<uint128_t,uint128_t,uint128_t> make_semiprime_u128(unsigned bits, std::mt19937_64& rng, int mr_rounds) {
    uint128_t p = random_probable_prime_u128(bits/2, rng, mr_rounds);
    uint128_t q = random_probable_prime_u128(bits/2, rng, mr_rounds);
    return {p*q, p, q};
}

// ------------------ argument helpers ------------------

static inline unsigned to_uint(const std::string& s, unsigned def) { try { return (unsigned)std::stoul(s); } catch(...) { return def; } }
static inline int to_int(const std::string& s, int def) { try { return std::stoi(s); } catch(...) { return def; } }

// ------------------ main ------------------

int main(int argc, char** argv) {
    int trials = 20;
    unsigned bits = 128;
    int B = 400;
    std::string score = "phase";
    bool neighbor_safety = false;
    bool shadow = true;
    std::uint64_t master_seed = 1234567;
    int mr_rounds = 12;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--trials" && i+1 < argc) { trials = to_int(argv[++i], trials); continue; }
        if (a == "--bits" && i+1 < argc) { bits = to_uint(argv[++i], bits); continue; }
        if (a == "--B" && i+1 < argc) { B = to_int(argv[++i], B); continue; }
        if (a == "--score" && i+1 < argc) { score = argv[++i]; continue; }
        if (a == "--neighbor-safety") { neighbor_safety = true; continue; }
        if (a == "--no-neighbor-safety") { neighbor_safety = false; continue; }
        if (a == "--shadow") { shadow = true; continue; }
        if (a == "--no-shadow") { shadow = false; continue; }
        if (a == "--seed" && i+1 < argc) { master_seed = (std::uint64_t)std::stoull(argv[++i]); continue; }
        if (a == "--mr-rounds" && i+1 < argc) { mr_rounds = to_int(argv[++i], mr_rounds); continue; }
    }

    if (bits > 128) {
        std::cerr << "Warning: bits>128 not supported in this benchmark build (u128). For >128 need bigint rho.\n";
        bits = 128;
    }

    CandidateSpec spec;
    std::mt19937_64 rng(master_seed);

    // dataset: semiprimes (generated with MR probable primes)
    std::vector<uint128_t> dataset;
    dataset.reserve(trials);
    for (int i = 0; i < trials; ++i) {
        auto [n, p, q] = make_semiprime_u128(bits, rng, mr_rounds);
        dataset.push_back(n);
    }

    std::cout << "\n=== 128-bit benchmark: Prefilter vs Miller-Rabin (+ Pollard-Rho) ===\n";
    std::cout << "trials=" << trials << "  bits=" << bits << "  B=" << B
              << "  score=" << score
              << "  neighbor_safety=" << (neighbor_safety ? "true" : "false")
              << "  shadow=" << (shadow ? "true" : "false")
              << "  mr_rounds=" << mr_rounds
              << "\n\n";

    std::cout << "i  cls  t_MR(s)   MR_prime?  t_pref(s)  pref_ok  tested   t_rho(s)  t_pref+rho(s)  rho_after_pref\n";
    std::cout << "---------------------------------------------------------------------------------------------------\n";

    double sum_mr = 0.0, sum_pref = 0.0, sum_rho = 0.0, sum_pref_rho = 0.0;
    long long sum_tested = 0;
    int pref_hits = 0;
    int rho_after_pref = 0;

    for (int i = 0; i < trials; ++i) {
        uint128_t n_u = dataset[i];
        cpp_int n = cpp_int(n_u);
        char cls = abc_like_class(n);

        // MR timing (primality test)
        auto t0 = std::chrono::steady_clock::now();
        bool mr_prime = mr_is_probable_prime_u128(n_u, mr_rounds, rng);
        double t_mr = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

        // Prefilter timing (factor attempt)
        auto t1 = std::chrono::steady_clock::now();
        auto pref = prefilter_pipeline(n, spec, B, score, (uint32_t)rng(), neighbor_safety, shadow);
        double t_pref = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();

        // Pollard-Rho direct factorization timing
        auto t2 = std::chrono::steady_clock::now();
        std::vector<uint128_t> facs_rho;
        factorize_rho_u128(n_u, mr_rounds, rng, facs_rho);
        std::sort(facs_rho.begin(), facs_rho.end());
        double t_rho = std::chrono::duration<double>(std::chrono::steady_clock::now() - t2).count();

        // Prefilter + (optional) rho
        auto t3 = std::chrono::steady_clock::now();
        bool used_rho = false;
        std::vector<uint128_t> facs_pref_rho;
        if (pref.success && pref.factor) {
            uint128_t f = pref.factor->convert_to<uint128_t>();
            facs_pref_rho.push_back(f);
            facs_pref_rho.push_back(n_u / f);
        } else {
            used_rho = true;
            factorize_rho_u128(n_u, mr_rounds, rng, facs_pref_rho);
        }
        std::sort(facs_pref_rho.begin(), facs_pref_rho.end());
        double t_pref_rho = std::chrono::duration<double>(std::chrono::steady_clock::now() - t3).count();

        // print row
        std::cout << std::setw(2) << (i+1) << "  "
                  << cls << "   "
                  << std::fixed << std::setprecision(5) << std::setw(8) << t_mr << "   "
                  << (mr_prime ? "Y" : "N") << "        "
                  << std::fixed << std::setprecision(5) << std::setw(8) << t_pref << "   "
                  << (pref.success ? "Y" : "N") << "      "
                  << std::setw(6) << pref.tested << "   "
                  << std::fixed << std::setprecision(5) << std::setw(8) << t_rho << "   "
                  << std::fixed << std::setprecision(5) << std::setw(12) << t_pref_rho << "        "
                  << (used_rho ? "Y" : "N")
                  << "\n";

        sum_mr += t_mr;
        sum_pref += t_pref;
        sum_rho += t_rho;
        sum_pref_rho += t_pref_rho;
        sum_tested += pref.tested;
        if (pref.success) pref_hits++;
        if (used_rho) rho_after_pref++;
    }

    std::cout << "\nSummary:\n";
    std::cout << "  mean_t_MR         = " << std::fixed << std::setprecision(6) << (sum_mr / trials) << " s\n";
    std::cout << "  mean_t_pref       = " << std::fixed << std::setprecision(6) << (sum_pref / trials) << " s\n";
    std::cout << "  mean_t_rho_direct = " << std::fixed << std::setprecision(6) << (sum_rho / trials) << " s\n";
    std::cout << "  mean_t_pref+rho   = " << std::fixed << std::setprecision(6) << (sum_pref_rho / trials) << " s\n";
    std::cout << "  prefilter_hit     = " << pref_hits << "/" << trials << "\n";
    std::cout << "  rho_after_pref    = " << rho_after_pref << "/" << trials << "\n";
    std::cout << "  mean_tested       = " << (double)sum_tested / trials << "\n\n";

    std::cout << "Note:\n";
    std::cout << "  - MR ist nur ein Primalitaetstest; Rho/PREF sind Faktorisierungsroutinen.\n";
    std::cout << "  - Der faire Speed-Check ist: mean_t_pref+rho vs mean_t_rho_direct.\n";
    std::cout << "  - Wenn prefilter_hit ~ 0, ist PREF fast sicher nur Overhead.\n";
    return 0;
}
