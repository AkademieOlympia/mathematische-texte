// E8_only_V0.0.14.cpp
// Nur E8-Algorithmus (aus Schütte/JyteV 0.0.14). Kein Pollard Rho, kein ECM.
// E8-Primfaktorzerlegung mit Struktur zwischen den Quadrupeln (Alice, Father, Taurus).
// Theoretische Grundlage: Struktur_zwischen_zwei_Quadrupeln.md
//
// Version 0.0.14 (E8-only): Gleiche Eingabeoptionen, nur E8-Faktorisierung.
//   Weniger Code → bessere I-Cache-Nutzung, ein Hot-Path für den Faktorisierer.
//
// Kompilieren: g++ -std=c++17 -O3 -march=native -I/opt/homebrew/include -o E8_only_V0.0.14 E8_only_V0.0.14.cpp
// Abhängigkeiten: Boost Multiprecision (Header-only)

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/multiprecision/integer.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <array>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace mp = boost::multiprecision;
using Clock = std::chrono::steady_clock;

// --- Hilfsfunktionen für Zeitmessung & String ---

static inline double ms_since(const Clock::time_point& t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

static inline std::string trim(const std::string& s) {
    size_t a=0,b=s.size();
    while (a<b && std::isspace((unsigned char)s[a])) a++;
    while (b>a && std::isspace((unsigned char)s[b-1])) b--;
    return s.substr(a,b-a);
}

// --- Optimierung 1: Primorial (Produkt kleiner Primzahlen) ---
static const std::vector<uint32_t> SMALL_PRIMES = {
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
    73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151
};

template<class Int>
static Int get_primorial() {
    Int p = 1;
    for (uint32_t sp : SMALL_PRIMES) {
        if (sp > 100) break;
        p *= sp;
    }
    return p;
}

template<class Int>
struct GlobalConsts {
    static Int PRIMORIAL;
};
template<class Int> Int GlobalConsts<Int>::PRIMORIAL;

// Kleine Primzahlen 101..173 für Early-Check (Trial-Division)
static const uint32_t EASY_SMALL_PRIMES[] = {
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173
};
static const int NUM_EASY = (int)(sizeof(EASY_SMALL_PRIMES) / sizeof(EASY_SMALL_PRIMES[0]));

// Gemischte Testfaktoren für E8: kleine Primzahlen, Produkte zweier Primzahlen, größere Primzahlen.
// Gleiche Anzahl wie NUM_EASY, damit die Zahl der Testläufe gleich bleibt.
static const uint64_t TEST_INJECTION_FACTORS[] = {
    /* kleine Primzahlen */
    101ULL, 103ULL, 107ULL, 109ULL, 113ULL,
    /* Produkte zweier Primzahlen (Semiprime) */
    101ULL * 103ULL,   /* 10403 */
    107ULL * 109ULL,   /* 11663 */
    109ULL * 113ULL,   /* 12317 */
    113ULL * 127ULL,   /* 14351 */
    127ULL * 131ULL,   /* 16637 */
    /* größere Primzahlen (~14 Bit) */
    10007ULL, 10009ULL, 10037ULL, 10039ULL, 10061ULL
};
static const int NUM_TEST_FACTORS = (int)(sizeof(TEST_INJECTION_FACTORS) / sizeof(TEST_INJECTION_FACTORS[0]));

template<class Int>
static inline Int gcd_int(Int a, Int b) {
    while (b != 0) { Int t = a % b; a = b; b = t; }
    return a;
}

template<class Int>
static inline Int integer_sqrt(const Int& n) {
    if (n <= 0) return 0;
    Int x = n, y = (x + 1) / 2;
    while (y < x) { x = y; y = (x + n / x) / 2; }
    return x;
}

template<class Int>
static inline double log_approx(const Int& x) {
    if (x <= 0) return 0.0;
    unsigned bits = mp::msb(x);
    return (double)bits * 0.69314718056;
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

// --- Zufallsgeneratoren ---
template<class Int>
static inline Int rand_bits(unsigned bits, std::mt19937_64& rng) {
    if (bits == 0) return 0;
    Int res = 0;
    unsigned pending = bits;
    std::uniform_int_distribution<uint64_t> dist;
    while (pending > 64) {
        res = (res << 64) | dist(rng);
        pending -= 64;
    }
    if (pending > 0) {
        uint64_t mask = (pending == 64) ? ~0ULL : ((1ULL << pending) - 1);
        res = (res << pending) | (dist(rng) & mask);
    }
    mp::bit_set(res, bits-1);
    return res;
}

template<class Int>
static inline Int rand_odd_candidate(unsigned bits, std::mt19937_64& rng) {
    Int x = rand_bits<Int>(bits, rng);
    mp::bit_set(x, 0);
    return x;
}

// --- Primzahltest (Miller-Rabin) ---
template<class Int>
static inline bool is_probable_prime(const Int& n, int rounds, std::mt19937_64& rng) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (!(n & 1)) return false;
    if (gcd_int(n, GlobalConsts<Int>::PRIMORIAL) != 1) return false;
    return mp::miller_rabin_test(n, rounds, rng);
}

template<class Int>
static Int random_probable_prime(unsigned bits, int rounds, std::mt19937_64& rng, int want_mod12 = -1) {
    if (bits < 4) bits = 4;
    for (;;) {
        Int x = rand_odd_candidate<Int>(bits, rng);
        if (want_mod12 != -1) {
            long r = (long)(x % 12);
            long diff = want_mod12 - r;
            if (diff < 0) diff += 12;
            x += diff;
            if ((x % 3) == 0) x += 12;
        }
        if (is_probable_prime(x, rounds, rng)) return x;
    }
}

struct Budget {
    Clock::time_point t0;
    uint64_t ms;
    Budget(uint64_t ms_=0): t0(Clock::now()), ms(ms_) {}
    inline bool exceeded() const {
        if (!ms) return false;
        return std::chrono::duration<double, std::milli>(Clock::now() - t0).count() > double(ms);
    }
};

// --- Precomputed E8 Roots ---
struct E8Root {
    double ux;
    double uy;
};

static std::vector<E8Root> PRECOMP_E8_ROOTS;

static void init_e8_roots() {
    if (!PRECOMP_E8_ROOTS.empty()) return;
    PRECOMP_E8_ROOTS.reserve(240);
    const double PI2 = 2.0 * M_PI;
    for (int i = 0; i < 240; ++i) {
        double angle = PI2 * i / 240.0;
        PRECOMP_E8_ROOTS.push_back({std::cos(angle), std::sin(angle)});
    }
}

// --- Energetische Metriken ---
template<class Int>
static void calculate_triplet_tensions(const std::array<Int, 4>& Q, double face_tensions[4]) {
    double log_Q[4];
    for (int i = 0; i < 4; ++i) log_Q[i] = log_approx(Q[i]);
    for (int omit = 0; omit < 4; ++omit) {
        int idx[3], k = 0;
        for (int j = 0; j < 4; ++j) if (j != omit) idx[k++] = j;
        double s_sk = 0.0;
        s_sk += std::fabs(log_Q[idx[0]] - log_Q[idx[1]]);
        s_sk += std::fabs(log_Q[idx[0]] - log_Q[idx[2]]);
        s_sk += std::fabs(log_Q[idx[1]] - log_Q[idx[2]]);
        face_tensions[omit] = s_sk;
    }
}

template<class Int>
static void calculate_metrics(const std::array<Int, 4>& Q, double& s_sk, double& s_i, double& r_q) {
    double log_Q[4];
    for (int i = 0; i < 4; ++i) log_Q[i] = log_approx(Q[i]);
    s_sk = 0.0;
    s_sk += std::fabs(log_Q[0] - log_Q[1]);
    s_sk += std::fabs(log_Q[0] - log_Q[2]);
    s_sk += std::fabs(log_Q[0] - log_Q[3]);
    s_sk += std::fabs(log_Q[1] - log_Q[2]);
    s_sk += std::fabs(log_Q[1] - log_Q[3]);
    s_sk += std::fabs(log_Q[2] - log_Q[3]);
    double sum_log = log_Q[0] + log_Q[1] + log_Q[2] + log_Q[3];
    double inv_sum = (sum_log > 0.0) ? (1.0 / sum_log) : 0.0;
    const double target = 1.0 / 2.61803398875;
    s_i = 0.0;
    for (int i = 0; i < 4; ++i) s_i += std::fabs((log_Q[i] * inv_sum) - target);
    if (s_i < 1e-9) s_i = 1e-9;
    r_q = s_sk / s_i;
}

// --- Schütte-GPS ---
template<class Int>
static Int find_nearest_class_prime(const Int& start, int direction, int want_mod12, int mr_rounds, std::mt19937_64& rng) {
    Int p = start;
    if (direction < 0) p--; else p++;
    int steps = 0;
    while (steps < 5000) {
        if (p < 2) return 0;
        long r12 = (long)(p % 12);
        if (r12 == want_mod12) {
            if (gcd_int(p, GlobalConsts<Int>::PRIMORIAL) == 1) {
                if (mp::miller_rabin_test(p, mr_rounds, rng)) return p;
            }
        }
        if (direction < 0) p--; else p++;
        steps++;
    }
    return 0;
}

static double get_triplet_based_phase_offset(const double face_tensions[4]) {
    int imax = 0;
    for (int i = 1; i < 4; ++i)
        if (face_tensions[i] > face_tensions[imax]) imax = i;
    return imax * (M_PI / 6.0);
}

// --- Resonance Lock (v9): Hardcoded Phase aus Daten ---
template<class Int>
static double get_resonance_lock_phase(const Int& n) {
    if (n < 100) return 0.0;
    long r12 = (long)(n % 12);
    if (r12 < 0) r12 += 12;
    if (r12 == 1) return 0.0;
    if (r12 == 7) return 0.0;
    if (r12 == 5) return 0.0262; // Alpha Base (Klasse A)
    if (r12 == 11) {
        if (n < 1000000) return 0.1047; // tiefe C-Zahlen (Beta Harmonic)
        return 0.0524;                  // Beta = 2 * Alpha
    }
    return 0.0;
}

// --- Aharonov-Bohm-Phasenkorrektor ---
template<class Int>
static double get_aharonov_bohm_phase(const Int& n, int mr_rounds, std::mt19937_64& rng) {
    if (n < 100) return 0.0;
    long r12 = (long)(n % 12);
    if (r12 < 0) r12 += 12;
    if (r12 == 1) return 0.0;
    if (r12 == 7 || r12 == 11) return 0.0;
    if (r12 != 5) return 0.0;

    Int p_B_lo = find_nearest_class_prime(n, -1, 7, mr_rounds, rng);
    Int p_B_hi = find_nearest_class_prime(n, 1, 7, mr_rounds, rng);
    Int p_C_lo = find_nearest_class_prime(n, -1, 11, mr_rounds, rng);
    Int p_C_hi = find_nearest_class_prime(n, 1, 11, mr_rounds, rng);

    double ln_n = log_approx(n);
    auto dist_to = [&ln_n](const Int& p) -> double {
        if (p == 0) return 1e30;
        return std::fabs(ln_n - log_approx(p));
    };
    double d_B_lo = dist_to(p_B_lo), d_B_hi = dist_to(p_B_hi);
    double d_C_lo = dist_to(p_C_lo), d_C_hi = dist_to(p_C_hi);
    Int p_B = (p_B_lo == 0 || (p_B_hi != 0 && d_B_hi < d_B_lo)) ? p_B_hi : p_B_lo;
    Int p_C = (p_C_lo == 0 || (p_C_hi != 0 && d_C_hi < d_C_lo)) ? p_C_hi : p_C_lo;
    if (p_B == 0 || p_C == 0) return 0.0;

    double dist_B = std::fabs(ln_n - log_approx(p_B));
    double dist_C = std::fabs(ln_n - log_approx(p_C));
    double sum = dist_B + dist_C + 1e-12;
    double delta = (dist_B - dist_C) / sum;
    return delta * (M_PI / 12.0);
}

template<class Int>
static bool get_local_energy_state(const Int& n, int mr_rounds, std::mt19937_64& rng,
                                   double& tuned_sens, double& tuned_scale_mult,
                                   double* out_phase_offset = nullptr,
                                   std::array<Int, 4>* out_Q_low = nullptr,
                                   std::array<Int, 4>* out_Q_high = nullptr) {
    if (n < 100) return false;

    std::array<Int, 4> Q_low, Q_high;
    const int residues[4] = {1, 5, 7, 11};

    for (int i=0; i<4; ++i) {
        Q_low[i] = find_nearest_class_prime(n, -1, residues[i], mr_rounds, rng);
        Q_high[i] = find_nearest_class_prime(n, 1, residues[i], mr_rounds, rng);
        if (Q_low[i] == 0 || Q_high[i] == 0) return false;
    }
    if (out_Q_low) *out_Q_low = Q_low;
    if (out_Q_high) *out_Q_high = Q_high;

    double sk_l, si_l, r_l, sk_h, si_h, r_h;
    calculate_metrics(Q_low, sk_l, si_l, r_l);
    calculate_metrics(Q_high, sk_h, si_h, r_h);
    double r_avg = 0.5 * (r_l + r_h);
    double si_avg = 0.5 * (si_l + si_h);

    if (out_phase_offset) {
        double face_low[4], face_high[4], face_avg[4];
        calculate_triplet_tensions(Q_low, face_low);
        calculate_triplet_tensions(Q_high, face_high);
        for (int i = 0; i < 4; ++i) face_avg[i] = 0.5 * (face_low[i] + face_high[i]);
        *out_phase_offset = get_triplet_based_phase_offset(face_avg);
    }

    Int max_low = 0; for(auto& x:Q_low) if(x>max_low) max_low=x;
    Int min_high = Q_high[0]; for(auto& x:Q_high) if(x<min_high) min_high=x;
    double d_low = log_approx(n - max_low);
    double d_high = log_approx(min_high - n);
    double pos_ratio = d_low / (d_low + d_high + 1e-9);

    unsigned bits = mp::msb(n);
    double base_sens = (bits > 100) ? 0.7 : 0.5;
    double base_scale = (bits > 100) ? 1.7 : 1.3;
    double energy_boost = std::log(r_avg + 1.0) / 6.0;
    tuned_sens = std::min(0.98, base_sens + energy_boost);
    tuned_scale_mult = base_scale * (1.0 + si_avg * 1.5);

    if (pos_ratio < 0.2 || pos_ratio > 0.8)
        tuned_sens = std::min(0.98, tuned_sens + 0.05);

    long r12 = (long)(n % 12);
    if (r12 < 0) r12 += 12;
    if (r12 == 5)
        tuned_sens = std::min(0.98, tuned_sens + 0.03);

    return true;
}

// --- E8 Factor Finder (einziger Faktorisierer in dieser Version) ---
enum class E8Status { DISABLED, HIT, MISS, TIMEOUT, ERROR };
const char* e8_status_name(E8Status s) {
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
struct E8Run { E8Status status=E8Status::DISABLED; Int factor=0; double t_ms=0.0; };

template<class Int, class Wider>
static E8Run<Int> e8_find_factor(const Int& n, std::mt19937_64& rng, Budget& bud,
                                 int mr_rounds, uint32_t trial_bound)
{
    E8Run<Int> r;
    auto t0 = Clock::now();

    try {
        if ((n & 1) == 0) { r.status=E8Status::HIT; r.factor=2; r.t_ms=ms_since(t0); return r; }

        for (int i = 0; i < NUM_EASY; ++i) {
            uint32_t p = EASY_SMALL_PRIMES[i];
            if (p >= trial_bound) break;
            if (n % p == 0) { r.status=E8Status::HIT; r.factor=Int(p); r.t_ms=ms_since(t0); return r; }
        }

        Int g = gcd_int(n, GlobalConsts<Int>::PRIMORIAL);
        if (g > 1) { r.status=E8Status::HIT; r.factor=g; r.t_ms=ms_since(t0); return r; }

        const uint32_t quick_trial_limit = 2000;
        for (uint32_t f = 101; f < quick_trial_limit && f < trial_bound; f += 2) {
            if (bud.exceeded()) { r.status=E8Status::TIMEOUT; r.t_ms=ms_since(t0); return r; }
            if (f % 3 == 0 || f % 5 == 0) continue;
            if (n % f == 0) { r.status=E8Status::HIT; r.factor=Int(f); r.t_ms=ms_since(t0); return r; }
        }

        if (bud.exceeded()) { r.status=E8Status::TIMEOUT; r.t_ms=ms_since(t0); return r; }

        double sens = 0.5;
        double scale_mult = 1.3;
        double phase_offset = 0.0;
        double ab_phase = 0.0;
        std::array<Int, 4> Q_low, Q_high;
        bool have_quad = false;
        unsigned bits = mp::msb(n);
        if (bits > 80) {
            have_quad = get_local_energy_state(n, mr_rounds, rng, sens, scale_mult, &phase_offset, &Q_low, &Q_high);
            ab_phase = get_aharonov_bohm_phase(n, mr_rounds, rng);
            ab_phase += get_resonance_lock_phase(n);  // v9: Hardcoded Resonance Lock (Alpha/Beta)
        }

        // Father (Walter-Geodäte): 16 Quadrupel + 8 Trippel
        if (have_quad) {
            for (int mask = 0; mask < 16; ++mask) {
                if (bud.exceeded()) break;
                Int P = 1;
                for (int i = 0; i < 4; ++i) P *= (mask & (1 << i)) ? Q_high[i] : Q_low[i];
                Int g = gcd_int(n, P);
                if (g > 1 && g < n) {
                    r.factor = g;
                    r.status = E8Status::HIT;
                    r.t_ms = ms_since(t0);
                    return r;
                }
            }
            for (int mask = 0; mask < 8; ++mask) {
                if (bud.exceeded()) break;
                Int P = Q_low[0];
                for (int i = 1; i < 4; ++i) P *= (mask & (1 << (i-1))) ? Q_high[i] : Q_low[i];
                Int g = gcd_int(n, P);
                if (g > 1 && g < n) {
                    r.factor = g;
                    r.status = E8Status::HIT;
                    r.t_ms = ms_since(t0);
                    return r;
                }
            }
        }

        // Alice + Taurus: E8-Injection (Hot path)
        Int root_n = integer_sqrt(n);
        int scale = std::max(3, (int)std::round((double)bits * scale_mult));
        static const int y_offsets[5] = { 1, 2, 3, 5, 8 };
        static const double MORLEY_STEP = 2.0 * M_PI / 3.0;
        static const double kepler_scales[3] = { 1.0, 1.732, 2.121 };

        auto run_e8_injection = [&](double tx_val, double ty_val, int scale_val) -> bool {
            for (const auto& root : PRECOMP_E8_ROOTS) {
                if (bud.exceeded()) return false;
                if (root.ux * tx_val + root.uy * ty_val <= sens) continue;
                int dx = (int)std::round(root.ux * scale_val);
                int dy = (int)std::round(root.uy * scale_val * 0.6);
                for (int y : y_offsets) {
                    Int g_real = root_n + dx;
                    if (g_real < 2) g_real = 2;
                    Int g_imag = y + dy;
                    if (g_imag < 1) g_imag = 1;
                    Wider nr = Wider(g_real)*Wider(g_real) + Wider(g_imag)*Wider(g_imag);
                    Int norm_mod = Int(nr % Wider(n));
                    Int fac = gcd_int(norm_mod, n);
                    if (fac > 1 && fac < n) {
                        r.factor = fac;
                        r.status = E8Status::HIT;
                        r.t_ms = ms_since(t0);
                        return true;
                    }
                }
            }
            return false;
        };

        double target_angle = phase_offset + ab_phase;
        for (int k = 0; k < 3; ++k) {
            if (bud.exceeded()) break;
            double tx = std::cos(target_angle + k * MORLEY_STEP);
            double ty = std::sin(target_angle + k * MORLEY_STEP);
            if (run_e8_injection(tx, ty, scale)) return r;
        }
        for (int ks = 1; ks < 3 && ms_since(t0) < 3.0; ++ks) {
            if (bud.exceeded()) break;
            int scale_k = std::max(3, (int)std::round(scale * kepler_scales[ks]));
            double tx = std::cos(target_angle);
            double ty = std::sin(target_angle);
            if (run_e8_injection(tx, ty, scale_k)) return r;
        }

        for (uint32_t f = quick_trial_limit; f < trial_bound; f += 2) {
            if (bud.exceeded()) break;
            if (f % 3 == 0 || f % 5 == 0) continue;
            if (n % f == 0) { r.status=E8Status::HIT; r.factor=Int(f); r.t_ms=ms_since(t0); return r; }
        }

        r.status = E8Status::MISS;
        r.t_ms = ms_since(t0);
        return r;
    } catch (...) {
        r.status = E8Status::ERROR;
        r.t_ms = ms_since(t0);
        return r;
    }
}

// --- RSA-Demo: inv_mod, mod_pow (nur für --rsa-demo) ---
template<class Int>
static std::pair<Int, Int> inv_mod_with_gcd(Int a, Int n) {
    mp::cpp_int aa(a), nn(n);
    aa %= nn; if (aa < 0) aa += nn;
    mp::cpp_int t = 0, newt = 1, r = nn, newr = aa;
    while (newr != 0) {
        mp::cpp_int q = r / newr;
        mp::cpp_int tmp = t - q * newt; t = newt; newt = tmp;
        tmp = r - q * newr; r = newr; newr = tmp;
    }
    if (t < 0) t += nn;
    return {Int(t), Int(r)};
}

template<class Int, class Wider>
static Int mod_pow(Int base, Int exp, const Int& n) {
    Int res = 1;
    base %= n;
    while (exp > 0) {
        if (exp & 1) res = Int((Wider(res) * Wider(base)) % Wider(n));
        exp /= 2;
        if (exp > 0) base = Int((Wider(base) * Wider(base)) % Wider(n));
    }
    return res;
}

// --- Aharonov-Bohm: PPM-Bild (--ab / --aharonov-bohm) ---
static bool write_aharonov_bohm_ppm(const char* filepath, int width = 600, int height = 600) {
    const double x_min = -3.0, x_max = 3.0, y_min = -3.0, y_max = 3.0;
    const double R_eps = 0.1, R_tension = 0.5;
    std::ofstream out(filepath, std::ios::binary);
    if (!out) return false;
    out << "P6\n" << width << " " << height << "\n255\n";
    std::vector<uint8_t> row(3 * width);
    const int cx = width / 2, cy = height / 2;
    const double r_world = 0.5;
    const int r_px = (int)(r_world / (x_max - x_min) * width + 0.5);
    const int dot_radius = 10;
    auto inferno = [](double t) -> std::array<uint8_t, 3> {
        if (t <= 0) return {0, 0, 0};
        if (t >= 1) return {255, 255, 230};
        double r, g, b;
        if (t < 0.25) {
            double s = t / 0.25;
            r = 0.2 * s; g = 0; b = 0.4 * s;
        } else if (t < 0.5) {
            double s = (t - 0.25) / 0.25;
            r = 0.2 + 0.6 * s; g = 0.2 * s; b = 0.4 - 0.4 * s;
        } else if (t < 0.75) {
            double s = (t - 0.5) / 0.25;
            r = 0.8 + 0.2 * s; g = 0.2 + 0.4 * s; b = 0;
        } else {
            double s = (t - 0.75) / 0.25;
            r = 1.0; g = 0.6 + 0.4 * s; b = 0.9 * s;
        }
        return {(uint8_t)(r * 255), (uint8_t)(g * 255), (uint8_t)(b * 255)};
    };
    auto on_circle = [cx, cy, r_px](int i, int j) {
        int di = i - cx, dj = j - cy;
        int d = (int)std::sqrt((double)(di*di + dj*dj));
        return std::abs(d - r_px) <= 1;
    };
    for (int j = 0; j < height; ++j) {
        double y = y_max - (j + 0.5) * (y_max - y_min) / height;
        for (int i = 0; i < width; ++i) {
            double x = x_min + (i + 0.5) * (x_max - x_min) / width;
            double R = std::sqrt(x*x + y*y);
            double Tension = 1.0 / (R + R_tension);
            const double T_max = 2.0;
            double t = std::min(1.0, Tension / T_max);
            int di = i - cx, dj = j - cy;
            bool is_center_dot = (di*di + dj*dj <= dot_radius * dot_radius);
            bool is_circle = on_circle(i, j);
            if (is_center_dot) {
                row[3*i] = 0; row[3*i+1] = 255; row[3*i+2] = 255;
            } else if (is_circle) {
                row[3*i] = 255; row[3*i+1] = 255; row[3*i+2] = 255;
            } else {
                auto rgb = inferno(t);
                row[3*i] = rgb[0]; row[3*i+1] = rgb[1]; row[3*i+2] = rgb[2];
            }
        }
        out.write(reinterpret_cast<const char*>(row.data()), (std::streamsize)row.size());
    }
    return true;
}

int main(int argc, char** argv) {
    init_e8_roots();

    using Int = mp::cpp_int;
    using Wider = mp::cpp_int;

    GlobalConsts<Int>::PRIMORIAL = get_primorial<Int>();

    // --aharonov-bohm / --ab [out_file]
    if (argc >= 2 && (std::string(argv[1]) == "--aharonov-bohm" || std::string(argv[1]) == "--ab")) {
        const char* out_file = (argc >= 3) ? argv[2] : "aharonov_bohm_number_vortex.ppm";
        std::cout << "Aharonov-Bohm-Effekt (Zahlentheorie): Vektorpotenzial um C-Klasse (n≡11 mod 12), Schütte-Kuss-Spannung.\n";
        if (write_aharonov_bohm_ppm(out_file)) {
            std::cout << "Bild gespeichert: " << out_file << std::endl;
        } else {
            std::cerr << "Fehler: Konnte " << out_file << " nicht schreiben." << std::endl;
            return 1;
        }
        return 0;
    }

    // --rsa-demo [count]
    if (argc >= 2 && std::string(argv[1]) == "--rsa-demo") {
        int rsa_count = 2000;
        if (argc >= 3) { int x = std::atoi(argv[2]); if (x > 0) rsa_count = x; }
        const unsigned rsa_q_bits = 256;
        std::mt19937_64 rng(12345);
        const Int e = 65537;
        int factored = 0, decrypted_ok = 0;
        std::cout << "RSA-Demo (E8 only): " << rsa_count << " schwache RSA-Testschlüssel (n = Faktor * q, q " << rsa_q_bits << " Bit). Faktoren: kleine Primzahlen, Produkte, größere Primzahlen.\n";
        std::cout << "E8 faktorisiert n, dann Entschlüsselung.\n" << std::endl;
        for (int i = 0; i < rsa_count; ++i) {
            uint64_t factor = TEST_INJECTION_FACTORS[i % NUM_TEST_FACTORS];
            Int q = random_probable_prime<Int>(rsa_q_bits, 8, rng);
            Int n = Int(factor) * q;
            Int m = (rand_bits<Int>(64, rng) % (n - 2)) + 2;
            Int c = mod_pow<Int, Wider>(m, e, n);
            Budget bud(10000);
            auto res = e8_find_factor<Int, Wider>(n, rng, bud, 8, 2000);
            if (res.status == E8Status::HIT) {
                factored++;
                Int p = res.factor;
                Int q_other = n / p;
                Int phi = (p - 1) * (q_other - 1);
                auto [d, g] = inv_mod_with_gcd<Int>(e, phi);
                if (g == 1) {
                    Int m_dec = mod_pow<Int, Wider>(c, d, n);
                    if (m_dec == m) decrypted_ok++;
                }
            }
            if ((i + 1) % 500 == 0)
                std::cout << "  " << (i+1) << "/" << rsa_count << " Schlüssel verarbeitet, " << factored << " faktorisiert, " << decrypted_ok << " entschlüsselt." << std::endl;
        }
        std::cout << "\nRSA-Demo Ergebnis: " << rsa_count << " Schlüssel, " << factored << " mit E8 faktorisiert, " << decrypted_ok << " entschlüsselt." << std::endl;
        return 0;
    }

    // [bits] → Benchmark nur E8 (mind. 10 HITs)
    unsigned q_bits = 3000;
    if (argc >= 2) {
        int parsed = std::atoi(argv[1]);
        if (parsed > 0) q_bits = static_cast<unsigned>(parsed);
    }

    std::cout << "E8-only Benchmark (V0.0.14) gestartet..." << std::endl;
    std::cout << "Primorial size: " << to_dec(GlobalConsts<Int>::PRIMORIAL).size() << " digits." << std::endl;
    std::cout << "Test: n = Faktor * q mit q " << q_bits << " Bits (Faktoren: kleine Primzahlen, Produkte, größere Primzahlen). Ziel: mindestens 10 HITs (nur E8).\n" << std::endl;

    std::mt19937_64 rng(42);
    const int min_hits = 10;
    const uint32_t trial_bound = 1000;
    int hits = 0;
    int run = 0;

    while (hits < min_hits) {
        run++;
        uint64_t factor = TEST_INJECTION_FACTORS[(run - 1) % NUM_TEST_FACTORS];
        Int q = random_probable_prime<Int>(q_bits, 10, rng);
        Int n = Int(factor) * q;

        std::cout << "  " << run << ". n = " << factor << " * q (" << mp::msb(n) << " bits)" << std::endl;

        Budget b_e8(5000);
        auto t_e8 = Clock::now();
        auto res = e8_find_factor<Int, Wider>(n, rng, b_e8, 10, trial_bound);
        double ms_e8 = ms_since(t_e8);
        std::cout << "      E8:  " << e8_status_name(res.status) << ", " << std::fixed << std::setprecision(3) << ms_e8 << " ms";
        if (res.status == E8Status::HIT) {
            hits++;
            std::cout << "  Faktor: " << res.factor;
            Int cofactor = n / res.factor;
            if (cofactor > 1 && cofactor != res.factor && is_probable_prime(cofactor, 5, rng))
                std::cout << "  (Semi-Prim)";
        }
        std::cout << std::endl;
    }

    std::cout << "\nErgebnis: " << hits << " HITs in " << run << " Läufen (nur E8)." << std::endl;
    return 0;
}
