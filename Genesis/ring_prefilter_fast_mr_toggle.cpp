/*
 * ring_prefilter_fast_mr_toggle.cpp
 * C++ port of ring_prefilter_fast_mr_toggle.py
 * Supports --bits up to 64 (uint64_t). For 96/128 bits use a bigint library.
 *
 * Build: g++ -std=c++17 -O2 -o ring_prefilter_fast_mr_toggle ring_prefilter_fast_mr_toggle.cpp
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using GI = std::pair<int64_t, int64_t>;
using EI = std::pair<int64_t, int64_t>;

// ============================================================
// 0) Utilities
// ============================================================

uint64_t isqrt(uint64_t n) {
    if (n == 0) return 0;
    uint64_t x = n;
    uint64_t y = (x + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
}

uint64_t gcd(uint64_t a, uint64_t b) {
    while (b) { uint64_t t = b; b = a % b; a = t; }
    return a;
}

int64_t gcd_signed(int64_t a, int64_t b) {
    uint64_t u = static_cast<uint64_t>(std::abs(a));
    uint64_t v = static_cast<uint64_t>(std::abs(b));
    return static_cast<int64_t>(gcd(u, v));
}

int64_t round_nearest_rational(int64_t num, int64_t den) {
    if (den <= 0) abort();
    if (num >= 0) {
        int64_t q = num / den;
        int64_t r = num % den;
        if (2 * r > den) return q + 1;
        if (2 * r < den) return q;
        return q + 1;
    }
    return -round_nearest_rational(-num, den);
}

// ============================================================
// 1) Jacobi + signature class
// ============================================================

int jacobi(int64_t a, int64_t n) {
    a %= n;
    if (a < 0) a += n;
    if (n <= 0 || (n & 1) == 0) return 0;
    int result = 1;
    while (a != 0) {
        while ((a & 1) == 0) {
            a /= 2;
            int r = static_cast<int>(n % 8);
            if (r == 3 || r == 5) result = -result;
        }
        std::swap(a, n);
        if ((a % 4) == 3 && (n % 4) == 3) result = -result;
        a %= n;
        if (a < 0) a += n;
    }
    return (n == 1) ? result : 0;
}

std::pair<int, int> v4_signature(uint64_t n) {
    int64_t n64 = static_cast<int64_t>(n);
    return { jacobi(-1, n64), jacobi(-3, n64) };
}

char abc_like_class(uint64_t n) {
    auto [chi4, chi3] = v4_signature(n);
    if (chi4 == -1 && chi3 == -1) return 'C';
    if (chi4 == 1 && chi3 == 1) return 'E';
    if (chi4 == 1 && chi3 == -1) return 'A';
    if (chi4 == -1 && chi3 == 1) return 'B';
    return 'U';
}

// ============================================================
// 2) Directions + caching
// ============================================================

static std::map<std::pair<int, int>, std::vector<GI>> dirs_cache;
static std::map<std::tuple<std::string, int, int>, std::vector<std::pair<double, double>>> dirs_s1_cache;

std::vector<GI> surrogate_dirs(int m, uint32_t seed = 0) {
    std::mt19937 rng(seed);
    std::vector<GI> dirs;
    std::set<GI> seen;
    std::uniform_int_distribution<int> dist(-200, 200);
    while (static_cast<int>(dirs.size()) < m) {
        int x = dist(rng), y = dist(rng);
        if (x == 0 && y == 0) continue;
        int g = static_cast<int>(gcd_signed(x, y));
        x /= g; y /= g;
        GI k = {x, y};
        GI kn = {-x, -y};
        if (seen.count(k) || seen.count(kn)) continue;
        seen.insert(k);
        dirs.push_back(k);
    }
    return dirs;
}

std::pair<double, double> normalize2(double x, double y) {
    double r = std::hypot(x, y);
    if (r == 0.0) return {0.0, 0.0};
    return {x / r, y / r};
}

constexpr double SQRT3 = 1.7320508075688772;
constexpr double SQRT3_OVER_2 = 0.8660254037844386;

std::pair<double, double> dirvec_gaussian(int dx, int dy) {
    return normalize2(static_cast<double>(dx), static_cast<double>(dy));
}

std::pair<double, double> dirvec_eisenstein(int da, int db) {
    double x = da - 0.5 * db;
    double y = SQRT3_OVER_2 * db;
    return normalize2(x, y);
}

const std::vector<GI>& get_surrogate_dirs_cached(int m, int seed) {
    auto key = std::make_pair(m, seed);
    auto it = dirs_cache.find(key);
    if (it != dirs_cache.end()) return it->second;
    dirs_cache[key] = surrogate_dirs(m, static_cast<uint32_t>(seed));
    return dirs_cache[key];
}

const std::vector<std::pair<double, double>>& get_dirs_s1_cached(const std::string& ring, int m, int seed) {
    auto key = std::make_tuple(ring, m, seed);
    auto it = dirs_s1_cache.find(key);
    if (it != dirs_s1_cache.end()) return it->second;
    const auto& raw = get_surrogate_dirs_cached(m, seed);
    std::vector<std::pair<double, double>> s1;
    for (const auto& [dx, dy] : raw) {
        if (ring == "gaussian")
            s1.push_back(dirvec_gaussian(dx, dy));
        else
            s1.push_back(dirvec_eisenstein(dx, dy));
    }
    dirs_s1_cache[key] = std::move(s1);
    return dirs_s1_cache[key];
}

// ============================================================
// 3) Phase / cheap scoring
// ============================================================

double phi_class_offset(char cls) {
    if (cls == 'E') return 0.0;
    if (cls == 'A') return M_PI / 2;
    if (cls == 'B') return M_PI;
    if (cls == 'C') return 3 * M_PI / 2;
    return 0.0;
}

double phi_of_n(uint64_t n) {
    char cls = abc_like_class(n);
    double base = phi_class_offset(cls);
    const double M = 2147483647.0;
    double x = std::sqrt(static_cast<double>(n)) / M;
    double frac = x - std::floor(x);
    return base + 2 * M_PI * frac;
}

std::pair<double, double> target_vec(uint64_t n) {
    double phi = phi_of_n(n);
    return {std::cos(phi), std::sin(phi)};
}

double alignment_topk(const std::vector<std::pair<double, double>>& dirs_s1,
                     const std::pair<double, double>& t, int topk = 8) {
    std::vector<double> dots;
    for (const auto& [rx, ry] : dirs_s1) {
        if (rx != 0.0 || ry != 0.0)
            dots.push_back(rx * t.first + ry * t.second);
    }
    if (dots.empty()) return -1.0;
    std::sort(dots.begin(), dots.end(), std::greater<double>());
    int k = std::min(topk, static_cast<int>(dots.size()));
    double sum = 0;
    for (int i = 0; i < k; ++i) sum += dots[i];
    return sum / k;
}

double ring_alignment(uint64_t n, int m, int seed, const std::string& ring, int topk = 8) {
    auto t = target_vec(n);
    const auto& dirs_s1 = get_dirs_s1_cached(ring, m, seed);
    return alignment_topk(dirs_s1, t, topk);
}

std::pair<int, int> split_budget(int B, double Ag, double Ae) {
    Ag = std::max(0.0, Ag);
    Ae = std::max(0.0, Ae);
    double s = Ag + Ae;
    if (s <= 1e-12) {
        int Bg = B / 2;
        return {Bg, B - Bg};
    }
    int Bg = static_cast<int>(std::round(B * (Ag / s)));
    Bg = std::max(0, std::min(B, Bg));
    return {Bg, B - Bg};
}

double score_phase_gaussian(const GI& g, uint64_t n) {
    const int M = 97;
    char cls = abc_like_class(n);
    int tx = 0, ty = 0;
    if (cls == 'A') tx = ty = M/3;
    else if (cls == 'B') { tx = 2*M/3; ty = M/3; }
    else if (cls == 'C') { tx = M/2; ty = 2*M/3; }
    int x = (static_cast<int>(g.first % M) + M) % M;
    int y = (static_cast<int>(g.second % M) + M) % M;
    int dx = std::min((x - tx + M) % M, (tx - x + M) % M);
    int dy = std::min((y - ty + M) % M, (ty - y + M) % M);
    return -(dx*dx + dy*dy);
}

double score_phase_eisenstein(const EI& g, uint64_t n) {
    const int M = 97;
    char cls = abc_like_class(n);
    int ta = 0, tb = 0;
    if (cls == 'A') ta = tb = M/3;
    else if (cls == 'B') { ta = 2*M/3; tb = M/3; }
    else if (cls == 'C') { ta = M/2; tb = 2*M/3; }
    int a = (static_cast<int>(g.first % M) + M) % M;
    int b = (static_cast<int>(g.second % M) + M) % M;
    int da = std::min((a - ta + M) % M, (ta - a + M) % M);
    int db = std::min((b - tb + M) % M, (tb - b + M) % M);
    return -(da*da + db*db);
}

double score_l1_gaussian(const GI& g, uint64_t n) {
    int64_t base = static_cast<int64_t>(isqrt(n));
    return -(std::abs(g.first - base) + std::abs(g.second));
}

double score_l1_eisenstein(const EI& g, uint64_t n) {
    int64_t base = static_cast<int64_t>(isqrt(n));
    return -(std::abs(g.first - base) + std::abs(g.second));
}

// ============================================================
// 4) Z[i]
// ============================================================

GI gi_mul(const GI& x, const GI& y) {
    int64_t a = x.first, b = x.second, c = y.first, d = y.second;
    return {a*c - b*d, a*d + b*c};
}

GI gi_conj(const GI& x) { return {x.first, -x.second}; }

GI gi_sub(const GI& x, const GI& y) {
    return {x.first - y.first, x.second - y.second};
}

int64_t gi_norm(const GI& x) {
    int64_t a = x.first, b = x.second;
    return a*a + b*b;
}

bool gi_is_zero(const GI& x) { return x.first == 0 && x.second == 0; }

std::pair<GI, GI> gi_divmod(const GI& alpha, const GI& beta) {
    if (gi_is_zero(beta)) abort();
    int64_t N = gi_norm(beta);
    GI z = gi_mul(alpha, gi_conj(beta));
    int64_t qa = round_nearest_rational(z.first, N);
    int64_t qb = round_nearest_rational(z.second, N);
    GI q = {qa, qb};
    GI r = gi_sub(alpha, gi_mul(q, beta));
    return {q, r};
}

GI gi_gcd(GI x, GI y) {
    while (!gi_is_zero(y)) {
        auto [q, r] = gi_divmod(x, y);
        x = y;
        y = r;
    }
    return x;
}

// ============================================================
// 5) Z[ω]
// ============================================================

EI ei_mul(const EI& x, const EI& y) {
    int64_t a = x.first, b = x.second, c = y.first, d = y.second;
    return {a*c - b*d, a*d + b*c - b*d};
}

EI ei_conj(const EI& x) { return {x.first - x.second, -x.second}; }

EI ei_sub(const EI& x, const EI& y) {
    return {x.first - y.first, x.second - y.second};
}

int64_t ei_norm(const EI& x) {
    int64_t a = x.first, b = x.second;
    return a*a - a*b + b*b;
}

bool ei_is_zero(const EI& x) { return x.first == 0 && x.second == 0; }

std::pair<double, double> ei_to_complex(const EI& x) {
    double a = x.first, b = x.second;
    return {a - b/2.0, SQRT3_OVER_2 * b};
}

EI nearest_eisenstein_from_complex(double zr, double zi) {
    double v = (2.0 / SQRT3) * zi;
    double u = zr + 0.5 * v;
    double w = -u - v;
    int ru = static_cast<int>(std::round(u));
    int rv = static_cast<int>(std::round(v));
    int rw = static_cast<int>(std::round(w));
    double du = std::abs(ru - u), dv = std::abs(rv - v), dw = std::abs(rw - w);
    if (du > dv && du > dw) ru = -rv - rw;
    else if (dv > dw) rv = -ru - rw;
    else rw = -ru - rv;
    return {ru, rv};
}

int ei_div_calls = 0, ei_neighbor_calls = 0;

std::pair<EI, EI> ei_divmod_adaptive(const EI& alpha, const EI& beta,
    bool neighbor_safety = true, bool adaptive = true, double trigger_ratio = 0.92) {
    ei_div_calls++;
    if (ei_is_zero(beta)) abort();
    int64_t Nbeta = ei_norm(beta);
    EI t = ei_mul(alpha, ei_conj(beta));
    auto [tr, ti] = ei_to_complex(t);
    double zr = tr / Nbeta, zi = ti / Nbeta;
    EI q0 = nearest_eisenstein_from_complex(zr, zi);
    EI r0 = ei_sub(alpha, ei_mul(q0, beta));
    if (!neighbor_safety) return {q0, r0};
    if (adaptive && ei_norm(r0) <= static_cast<int64_t>(trigger_ratio * Nbeta)) return {q0, r0};
    ei_neighbor_calls++;
    std::vector<std::pair<int,int>> neighbors = {{0,0},{1,0},{-1,0},{0,1},{0,-1},{1,1},{-1,-1}};
    EI best_q = q0, best_r = r0;
    int64_t best_n = ei_norm(r0);
    for (const auto& [da, db] : neighbors) {
        EI q = {q0.first + da, q0.second + db};
        EI r = ei_sub(alpha, ei_mul(q, beta));
        int64_t nr = ei_norm(r);
        if (nr < best_n) { best_n = nr; best_q = q; best_r = r; }
    }
    return {best_q, best_r};
}

EI ei_gcd(EI x, EI y, bool neighbor_safety = true) {
    while (!ei_is_zero(y)) {
        auto [q, r] = ei_divmod_adaptive(x, y, neighbor_safety, true, 0.92);
        x = y;
        y = r;
    }
    return x;
}

// ============================================================
// 6) CandidateSpec + streams + top-B heap
// ============================================================

struct CandidateSpec {
    int T = 40, S = 40, base_samples = 8, dir_count = 120;
    std::vector<int> lambdas = {1, 2, 3, 5, 8};
    int dir_seed_offset = 12345;
};

struct HeapItem {
    double score;
    int counter;
    GI cand;
    bool operator>(const HeapItem& o) const { return score > o.score; }
};

std::vector<GI> candidates_gaussian_stream(uint64_t n, const CandidateSpec& spec, uint32_t seed) {
    std::vector<GI> out;
    std::mt19937 rng(seed);
    int64_t base = static_cast<int64_t>(isqrt(n));
    std::uniform_int_distribution<int> distT(-spec.T, spec.T), distS(-spec.S, spec.S);
    const auto& dirs = get_surrogate_dirs_cached(spec.dir_count, seed + spec.dir_seed_offset);
    for (int i = 0; i < spec.base_samples; ++i) {
        int t = distT(rng), s = distS(rng);
        GI g0 = {base + t, s};
        out.push_back(g0);
        for (const auto& [dx, dy] : dirs) {
            for (int lam : spec.lambdas)
                out.push_back({g0.first + lam*dx, g0.second + lam*dy});
        }
    }
    return out;
}

std::vector<EI> candidates_eisenstein_stream(uint64_t n, const CandidateSpec& spec, uint32_t seed) {
    std::vector<EI> out;
    std::mt19937 rng(seed);
    int64_t base = static_cast<int64_t>(isqrt(n));
    std::uniform_int_distribution<int> distT(-spec.T, spec.T), distS(-spec.S, spec.S);
    const auto& dirs = get_surrogate_dirs_cached(spec.dir_count, seed + spec.dir_seed_offset + 777);
    for (int i = 0; i < spec.base_samples; ++i) {
        int t = distT(rng), s = distS(rng);
        EI g0 = {base + t, s};
        out.push_back(g0);
        for (const auto& [da, db] : dirs) {
            for (int lam : spec.lambdas)
                out.push_back({g0.first + lam*da, g0.second + lam*db});
        }
    }
    return out;
}

std::vector<GI> top_b_candidates_gaussian(
    const std::vector<GI>& candidates,
    std::function<double(const GI&)> score_fn, int B) {
    if (B <= 0) return {};
    std::priority_queue<HeapItem, std::vector<HeapItem>, std::greater<HeapItem>> heap;
    int counter = 0;
    for (const auto& cand : candidates) {
        double s = score_fn(cand);
        HeapItem item{s, counter++, cand};
        if (static_cast<int>(heap.size()) < B)
            heap.push(item);
        else if (s > heap.top().score) {
            heap.pop();
            heap.push(item);
        }
    }
    std::vector<GI> out;
    while (!heap.empty()) { out.push_back(heap.top().cand); heap.pop(); }
    std::reverse(out.begin(), out.end());
    return out;
}

std::vector<EI> top_b_candidates_eisenstein(
    const std::vector<EI>& candidates,
    std::function<double(const EI&)> score_fn, int B) {
    if (B <= 0) return {};
    struct EIHeapItem { double score; int counter; EI cand; bool operator>(const EIHeapItem& o) const { return score > o.score; } };
    std::priority_queue<EIHeapItem, std::vector<EIHeapItem>, std::greater<EIHeapItem>> heap;
    int counter = 0;
    for (const auto& cand : candidates) {
        double s = score_fn(cand);
        EIHeapItem item{s, counter++, cand};
        if (static_cast<int>(heap.size()) < B)
            heap.push(item);
        else if (s > heap.top().score) {
            heap.pop();
            heap.push(item);
        }
    }
    std::vector<EI> out;
    while (!heap.empty()) { out.push_back(heap.top().cand); heap.pop(); }
    std::reverse(out.begin(), out.end());
    return out;
}

// ============================================================
// 7) Prefilter
// ============================================================

struct PrefilterResult {
    bool success = false;
    std::optional<uint64_t> factor;
    std::optional<int> rank_hit;
    int tested = 0;
    double t_prefilter = 0;
    std::string ring_used;
};

std::optional<uint64_t> try_factor_gaussian(uint64_t n, const GI& g) {
    GI d = gi_gcd({static_cast<int64_t>(n), 0}, g);
    uint64_t nd = static_cast<uint64_t>(std::abs(gi_norm(d)));
    uint64_t f = gcd(nd, n);
    if (f > 1 && f < n) return f;
    return std::nullopt;
}

std::optional<uint64_t> try_factor_eisenstein(uint64_t n, const EI& g, bool neighbor_safety = true) {
    EI d = ei_gcd({static_cast<int64_t>(n), 0}, g, neighbor_safety);
    uint64_t nd = static_cast<uint64_t>(std::abs(ei_norm(d)));
    uint64_t f = gcd(nd, n);
    if (f > 1 && f < n) return f;
    return std::nullopt;
}

PrefilterResult prefilter_gaussian_fast(uint64_t n, const CandidateSpec& spec, int B,
    const std::string& score_mode, uint32_t seed) {
    auto t0 = std::chrono::steady_clock::now();
    std::function<double(const GI&)> score_fn;
    if (score_mode == "phase")
        score_fn = [n](const GI& g) { return score_phase_gaussian(g, n); };
    else
        score_fn = [n](const GI& g) { return score_l1_gaussian(g, n); };
    std::vector<GI> cands = candidates_gaussian_stream(n, spec, seed);
    std::vector<GI> ordered = top_b_candidates_gaussian(cands, score_fn, B);
    PrefilterResult res;
    res.ring_used = "gaussian";
    for (size_t idx = 0; idx < ordered.size(); ++idx) {
        res.tested++;
        auto f = try_factor_gaussian(n, ordered[idx]);
        if (f) {
            res.success = true;
            res.factor = *f;
            res.rank_hit = static_cast<int>(idx + 1);
            res.t_prefilter = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            return res;
        }
    }
    res.t_prefilter = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return res;
}

PrefilterResult prefilter_eisenstein_fast(uint64_t n, const CandidateSpec& spec, int B,
    const std::string& score_mode, uint32_t seed, bool neighbor_safety = true) {
    auto t0 = std::chrono::steady_clock::now();
    std::function<double(const EI&)> score_fn;
    if (score_mode == "phase")
        score_fn = [n](const EI& g) { return score_phase_eisenstein(g, n); };
    else
        score_fn = [n](const EI& g) { return score_l1_eisenstein(g, n); };
    std::vector<EI> cands = candidates_eisenstein_stream(n, spec, seed);
    std::vector<EI> ordered = top_b_candidates_eisenstein(cands, score_fn, B);
    PrefilterResult res;
    res.ring_used = "eisenstein";
    for (size_t idx = 0; idx < ordered.size(); ++idx) {
        res.tested++;
        auto f = try_factor_eisenstein(n, ordered[idx], neighbor_safety);
        if (f) {
            res.success = true;
            res.factor = *f;
            res.rank_hit = static_cast<int>(idx + 1);
            res.t_prefilter = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            return res;
        }
    }
    res.t_prefilter = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return res;
}

PrefilterResult prefilter_pipeline(uint64_t n, const CandidateSpec& spec, int B,
    const std::string& score_mode, uint32_t seed, bool neighbor_safety = true) {
    char cls = abc_like_class(n);
    if (cls == 'C') {
        for (int k : {5, 7}) {
            uint64_t N = n * k;
            auto res = prefilter_pipeline(N, spec, B, score_mode, seed + 1000 + k, neighbor_safety);
            if (res.success && res.factor) {
                uint64_t f = gcd(*res.factor, n);
                if (f > 1 && f < n) {
                    res.factor = f;
                    res.ring_used = "shadow";
                    return res;
                }
            }
        }
        return PrefilterResult{false, std::nullopt, std::nullopt, 0, 0.0, "shadow"};
    }
    if (cls == 'A') return prefilter_gaussian_fast(n, spec, B, score_mode, seed);
    if (cls == 'B') return prefilter_eisenstein_fast(n, spec, B, score_mode, seed, neighbor_safety);
    if (cls == 'E') {
        int m_align = std::min(spec.dir_count, 120);
        double Ag = ring_alignment(n, m_align, seed + 2222, "gaussian", 8);
        double Ae = ring_alignment(n, m_align, seed + 2222, "eisenstein", 8);
        auto [Bg, Be] = split_budget(B, Ag, Ae);
        PrefilterResult r1, r2;
        if (Ag >= Ae) {
            r1 = (Bg > 0) ? prefilter_gaussian_fast(n, spec, Bg, score_mode, seed + 1) : PrefilterResult{};
            if (r1.success) return r1;
            r2 = (Be > 0) ? prefilter_eisenstein_fast(n, spec, Be, score_mode, seed + 2, neighbor_safety) : PrefilterResult{};
        } else {
            r1 = (Be > 0) ? prefilter_eisenstein_fast(n, spec, Be, score_mode, seed + 1, neighbor_safety) : PrefilterResult{};
            if (r1.success) return r1;
            r2 = (Bg > 0) ? prefilter_gaussian_fast(n, spec, Bg, score_mode, seed + 2) : PrefilterResult{};
        }
        if (r2.success) {
            r2.tested += r1.tested;
            r2.t_prefilter += r1.t_prefilter;
            return r2;
        }
        PrefilterResult out;
        out.tested = r1.tested + r2.tested;
        out.t_prefilter = r1.t_prefilter + r2.t_prefilter;
        out.ring_used = "E-both";
        return out;
    }
    return prefilter_gaussian_fast(n, spec, B, score_mode, seed);
}

// ============================================================
// 8) Primality
// ============================================================

bool trial_division_is_prime(uint64_t n, uint64_t limit = 200000) {
    if (n < 2) return false;
    for (uint64_t p : {2,3,5,7,11,13,17,19,23,29,31,37}) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }
    uint64_t r = isqrt(n);
    if (r > limit) return false;
    for (uint64_t f = 41; f <= r; f += 2)
        if (n % f == 0) return false;
    return true;
}

static uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t n) {
#ifdef __SIZEOF_INT128__
    return static_cast<uint64_t>(__uint128_t(a) * b % n);
#else
    uint64_t r = 0;
    a %= n;
    while (b) {
        if (b & 1) r = (r + a) % n;
        a = (a + a) % n;
        b /= 2;
    }
    return r;
#endif
}

uint64_t pow_mod(uint64_t a, uint64_t d, uint64_t n) {
    uint64_t r = 1;
    a %= n;
    while (d) {
        if (d & 1) r = mul_mod(r, a, n);
        a = mul_mod(a, a, n);
        d /= 2;
    }
    return r;
}

bool miller_rabin_is_prime(uint64_t n, int rounds, std::mt19937_64& rng) {
    if (n < 2) return false;
    for (uint64_t p : {2,3,5,7,11,13,17,19,23,29,31,37}) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }
    uint64_t d = n - 1;
    int s = 0;
    while ((d & 1) == 0) { s++; d /= 2; }
    std::uniform_int_distribution<uint64_t> dist(2, n - 2);
    for (int i = 0; i < rounds; ++i) {
        uint64_t a = dist(rng);
        uint64_t x = pow_mod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int r = 0; r < s - 1; ++r) {
            x = mul_mod(x, x, n);
            if (x == n - 1) { composite = false; break; }
        }
        if (composite) return false;
    }
    return true;
}

bool is_prime(uint64_t n, bool use_mr, std::mt19937_64& rng) {
    return use_mr ? miller_rabin_is_prime(n, 12, rng) : trial_division_is_prime(n, 200000);
}

// ============================================================
// 9) Pollard-Rho + factorize
// ============================================================

uint64_t pollard_rho(uint64_t n, std::mt19937_64& rng) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    std::uniform_int_distribution<uint64_t> dist_c(1, n - 2), dist_x(0, n - 2);
    while (true) {
        uint64_t c = dist_c(rng);
        uint64_t x = dist_x(rng);
        uint64_t y = x;
        uint64_t d = 1;
        auto f = [n, c](uint64_t v) { return (mul_mod(v, v, n) + c) % n; };
        while (d == 1) {
            x = f(x);
            y = f(f(y));
            int64_t diff = static_cast<int64_t>(x) - static_cast<int64_t>(y);
            if (diff < 0) diff = -diff;
            d = gcd(static_cast<uint64_t>(diff), n);
        }
        if (d != n) return d;
    }
}

std::vector<uint64_t> factorize(uint64_t n, bool use_mr, std::mt19937_64& rng) {
    if (n == 1) return {};
    if (is_prime(n, use_mr, rng)) return {n};
    uint64_t d = pollard_rho(n, rng);
    auto a = factorize(d, use_mr, rng);
    auto b = factorize(n / d, use_mr, rng);
    a.insert(a.end(), b.begin(), b.end());
    std::sort(a.begin(), a.end());
    return a;
}

// ============================================================
// 10) Semiprime generation
// ============================================================

uint64_t random_probable_prime(int bits, std::mt19937_64& rng, bool use_mr) {
    uint64_t low = (bits >= 64) ? (1ULL << 63) : (1ULL << (bits - 1));
    uint64_t high = (bits >= 64) ? (UINT64_MAX - 1) : ((1ULL << bits) - 1);
    if (high < low) high = low;
    std::uniform_int_distribution<uint64_t> dist(low, high);
    while (true) {
        uint64_t x = dist(rng);
        x |= 1;
        if (x >= 2 && is_prime(x, use_mr, rng)) return x;
    }
}

std::tuple<uint64_t, uint64_t, uint64_t> make_semiprime(int bits, std::mt19937_64& rng, bool use_mr) {
    uint64_t p = random_probable_prime(bits / 2, rng, use_mr);
    uint64_t q = random_probable_prime(bits / 2, rng, use_mr);
    return {p * q, p, q};
}

// ============================================================
// 11) Runner
// ============================================================

struct RunResult {
    uint64_t n;
    char cls;
    bool prefilter_ok;
    std::optional<uint64_t> prefilter_factor;
    std::string prefilter_ring;
    int prefilter_tested;
    double t_prefilter;
    int ei_div_calls, ei_neighbor_calls;
    std::vector<uint64_t> final_factors;
    double t_total;
    bool fallback_used;
};

RunResult run_once(uint64_t n, const CandidateSpec& spec, int B, const std::string& score_mode,
    uint32_t seed, bool neighbor_safety, bool use_fallback, bool use_mr) {
    std::mt19937_64 rng(seed);
    RunResult out;
    out.n = n;
    out.cls = abc_like_class(n);
    ei_div_calls = 0;
    ei_neighbor_calls = 0;
    auto t0 = std::chrono::steady_clock::now();
    std::uniform_int_distribution<uint32_t> dist_seed(1, 1000000000);
    uint32_t pipe_seed = dist_seed(rng);
    auto res = prefilter_pipeline(n, spec, B, score_mode, pipe_seed, neighbor_safety);
    out.t_prefilter = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    out.prefilter_ok = res.success;
    out.prefilter_factor = res.factor;
    out.prefilter_ring = res.ring_used;
    out.prefilter_tested = res.tested;
    out.ei_div_calls = ei_div_calls;
    out.ei_neighbor_calls = ei_neighbor_calls;
    if (res.success) {
        out.final_factors = {*res.factor, n / *res.factor};
        std::sort(out.final_factors.begin(), out.final_factors.end());
        out.t_total = out.t_prefilter;
        out.fallback_used = false;
        return out;
    }
    if (!use_fallback) {
        out.t_total = out.t_prefilter;
        out.fallback_used = false;
        return out;
    }
    auto t1 = std::chrono::steady_clock::now();
    out.final_factors = factorize(n, use_mr, rng);
    double t_fb = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
    out.fallback_used = true;
    out.t_total = out.t_prefilter + t_fb;
    return out;
}

void print_table(const std::vector<std::vector<std::string>>& rows) {
    if (rows.empty()) return;
    size_t cols = rows[0].size();
    std::vector<size_t> widths(cols, 0);
    for (const auto& r : rows)
        for (size_t i = 0; i < cols && i < r.size(); ++i)
            widths[i] = std::max(widths[i], r[i].size());
    for (size_t idx = 0; idx < rows.size(); ++idx) {
        const auto& r = rows[idx];
        for (size_t i = 0; i < cols; ++i)
            std::cout << std::left << std::setw(static_cast<int>(widths[i])) << (i < r.size() ? r[i] : "") << "  ";
        std::cout << "\n";
        if (idx == 0) {
            std::cout << "\n";
            for (size_t i = 0; i < cols; ++i)
                std::cout << std::string(widths[i], '-') << "  ";
        }
        std::cout << "\n";
    }
}

// ============================================================
// main
// ============================================================

int main(int argc, char** argv) {
    int trials = 20, bits = 64, B = 400;
    std::string score = "phase";
    bool neighbor_safety = false, fallback = false, use_mr = true, compare = false;
    uint32_t master_seed = 1234567;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--trials" && i+1 < argc) { trials = std::stoi(argv[++i]); continue; }
        if (a == "--bits" && i+1 < argc) { bits = std::stoi(argv[++i]); continue; }
        if (a == "--B" && i+1 < argc) { B = std::stoi(argv[++i]); continue; }
        if (a == "--score" && i+1 < argc) { score = argv[++i]; continue; }
        if (a == "--neighbor-safety") { neighbor_safety = true; continue; }
        if (a == "--fallback") { fallback = true; continue; }
        if (a == "--mr" && i+1 < argc) { use_mr = (std::string(argv[++i]) == "on"); continue; }
        if (a == "--compare") { compare = true; continue; }
        if (a == "--seed" && i+1 < argc) { master_seed = static_cast<uint32_t>(std::stoul(argv[++i])); continue; }
    }
    if (bits > 64) {
        std::cerr << "Warning: --bits > 64 not fully supported (uint64_t); use 64 or less.\n";
        bits = 64;
    }
    CandidateSpec spec;
    std::mt19937_64 rng(master_seed);
    std::vector<uint64_t> dataset;
    for (int i = 0; i < trials; ++i) {
        auto [n, p, q] = make_semiprime(bits, rng, true);
        dataset.push_back(n);
    }
    auto run_suite = [&](bool use_mr_flag) {
        double total_t = 0;
        int pref_succ = 0, fallback_used = 0;
        std::vector<std::vector<std::string>> rows = {{"i","cls","pref_ok","ring","tested","t_pref","fallback","t_total","factors"}};
        for (size_t i = 0; i < dataset.size(); ++i) {
            auto res = run_once(dataset[i], spec, B, score, master_seed + static_cast<uint32_t>((i+1)*99991),
                neighbor_safety, fallback, use_mr_flag);
            total_t += res.t_total;
            if (res.prefilter_ok) pref_succ++;
            if (res.fallback_used) fallback_used++;
            std::string facs_s = "-";
            if (!res.final_factors.empty()) {
                std::ostringstream os;
                if (res.final_factors.size() <= 4) {
                    for (size_t j = 0; j < res.final_factors.size(); ++j)
                        os << (j ? "*" : "") << res.final_factors[j];
                } else
                    os << res.final_factors[0] << "*...*" << res.final_factors.back();
                facs_s = os.str();
            }
            std::ostringstream os1, os2;
            os1 << std::fixed << std::setprecision(4) << res.t_prefilter;
            os2 << std::fixed << std::setprecision(4) << res.t_total;
            rows.push_back({
                std::to_string(i+1),
                std::string(1, res.cls),
                res.prefilter_ok ? "Y" : "N",
                res.prefilter_ring,
                std::to_string(res.prefilter_tested),
                os1.str(),
                res.fallback_used ? "Y" : "N",
                os2.str(),
                facs_s
            });
        }
        return std::make_tuple(rows, total_t, total_t / dataset.size(), pref_succ, fallback_used);
    };
    if (compare) {
        auto [rows_on, tot_on, mean_on, pref_on, fb_on] = run_suite(true);
        auto [rows_off, tot_off, mean_off, pref_off, fb_off] = run_suite(false);
        std::cout << "\n=== Dataset (generated with MR=ON) ===\n";
        std::cout << "trials=" << trials << "  bits=" << bits << "  B=" << B << "  score=" << score
                  << "  neighbor_safety=" << neighbor_safety << "  fallback=" << fallback << "\n";
        std::cout << "\n=== Run A: MR=ON ===\n";
        print_table(std::vector(rows_on.begin(), rows_on.begin() + std::min(rows_on.size(), 1 + size_t(std::min(12, trials)))));
        std::cout << "\nSummary MR=ON: prefilter_success=" << pref_on << "/" << trials << "  fallback_used=" << fb_on << "/" << trials << "  mean_total_time=" << std::fixed << std::setprecision(4) << mean_on << "s\n";
        std::cout << "\n=== Run B: MR=OFF ===\n";
        print_table(std::vector(rows_off.begin(), rows_off.begin() + std::min(rows_off.size(), 1 + size_t(std::min(12, trials)))));
        std::cout << "\nSummary MR=OFF: prefilter_success=" << pref_off << "/" << trials << "  fallback_used=" << fb_off << "/" << trials << "  mean_total_time=" << mean_off << "s\n";
        std::cout << "\nNote: With MR=OFF, fallback can be very slow.\n";
        return 0;
    }
    std::cout << "\n=== Single run mode ===\n";
    std::cout << "trials=" << trials << "  bits=" << bits << "  B=" << B << "  score=" << score
              << "  neighbor_safety=" << (neighbor_safety ? "true" : "false") << "  fallback=" << (fallback ? "true" : "false") << "  mr=" << (use_mr ? "on" : "off") << "\n";
    auto [rows, tot, mean, pref_succ, fb_used] = run_suite(use_mr);
    print_table(rows);
    std::cout << "\nSummary: prefilter_success=" << pref_succ << "/" << trials << "  fallback_used=" << fb_used << "/" << trials << "  mean_total_time=" << std::fixed << std::setprecision(4) << mean << "s\n";
    return 0;
}
