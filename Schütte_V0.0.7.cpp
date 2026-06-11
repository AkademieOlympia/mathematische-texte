// Schütte_V0.0.7.cpp
// E8-Primfaktorzerlegung mit Struktur zwischen den Quadrupeln (Alice, Father, Taurus).
// Theoretische Grundlage: Struktur_zwischen_zwei_Quadrupeln.md
//
// Version 0.0.7: Precomputed E8 Roots, Primorial-GCD, Schütte-GPS, drei Geodäten:
//   Father (Walter): 16 Quadrupel-Produkte Q_low/Q_high + 8 Trippel-Produkte (A,B,C).
//   Alice (Morley):  E8-Injektion mit 3 Phasen (0, 2π/3, 4π/3).
//   Taurus (Kepler): E8-Injektion mit 3 log-Skalen. Semi-Prim-Erkennung.
//
// Kompilieren: g++ -std=c++17 -O3 -march=native -I/opt/homebrew/include -o Schütte_V0.0.7 Schütte_V0.0.7.cpp
// Abhängigkeiten: Boost Multiprecision (Header-only)

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/multiprecision/integer.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
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
// Statt hunderter Divisionen machen wir einen GCD mit diesem Produkt.
// Enthält alle Primzahlen bis 100 (Produkt passt locker in uint128/256/cpp_int).
// 2*3*5*...*97 ist ca. 10^38, passt also in 128 bit.

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

// Globale Cache-Struktur für Primorials (wird in main init)
template<class Int>
struct GlobalConsts {
    static Int PRIMORIAL;
};
template<class Int> Int GlobalConsts<Int>::PRIMORIAL;

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

// Optimierte Log-Approx (schneller Bit-Scan)
template<class Int>
static inline double log_approx(const Int& x) {
    if (x <= 0) return 0.0;
    // Boost msb liefert den Index des höchsten Bits (0-basiert)
    // log2(x) approx msb(x). ln(x) = log2(x) * ln(2)
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
    // Generiert n Bits Zufall
    // Für Multiprecision müssen wir 64-bit Blöcke zusammensetzen
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
    // Setze höchstes Bit um Bitlänge zu garantieren
    mp::bit_set(res, bits-1);
    return res;
}

template<class Int>
static inline Int rand_odd_candidate(unsigned bits, std::mt19937_64& rng) {
    Int x = rand_bits<Int>(bits, rng);
    mp::bit_set(x, 0); // Mach ungerade
    return x;
}

// --- Primzahltest (Miller-Rabin) ---

template<class Int>
static inline bool is_probable_prime(const Int& n, int rounds, std::mt19937_64& rng) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (!(n & 1)) return false;

    // Schneller GCD-Check vor teurem Miller-Rabin
    if (gcd_int(n, GlobalConsts<Int>::PRIMORIAL) != 1) return false;

    return mp::miller_rabin_test(n, rounds, rng);
}

template<class Int>
static Int random_probable_prime(unsigned bits, int rounds, std::mt19937_64& rng, int want_mod12 = -1) {
    if (bits < 4) bits = 4;
    for (;;) {
        Int x = rand_odd_candidate<Int>(bits, rng);
        
        // Modulo 12 Anpassung
        if (want_mod12 != -1) {
             long r = (long)(x % 12);
             long diff = want_mod12 - r;
             if (diff < 0) diff += 12;
             x += diff;
             // Prüfen ob x durch 3 teilbar wurde (kann passieren bei +diff)
             if ((x % 3) == 0) x += 12; // Zum nächsten gleichen Modulo springen
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

// --- State-of-the-art: Pollard Rho (Brent) ---
// Klassischer Faktorisierungsalgorithmus zum Vergleich mit E8/Schütte.
template<class Int>
static Int rho_f(const Int& x, const Int& c, const Int& n) {
    return (x * x + c) % n;
}

template<class Int>
static Int pollard_rho_brent(const Int& n, std::mt19937_64& rng, Budget& bud,
                             int max_restarts = 256, uint32_t m_block = 256) {
    if ((n & 1) == 0) return 2;
    for (uint32_t p : SMALL_PRIMES) if (n % p == 0) return Int(p);

    std::uniform_int_distribution<uint64_t> dist64;
    for (int restart = 0; restart < max_restarts; ++restart) {
        if (bud.exceeded()) return 0;
        Int y = Int(dist64(rng)) % (n - 1) + 1;
        Int c = Int(dist64(rng)) % (n - 1) + 1;

        Int g = 1, r = 1, q = 1;
        Int x = 0, ys = 0;
        Int m = Int(m_block);

        while (g == 1) {
            if (bud.exceeded()) return 0;
            x = y;
            for (Int i = 0; i < r; ++i) {
                y = rho_f<Int>(y, c, n);
                if (bud.exceeded()) return 0;
            }

            Int k = 0;
            while (k < r && g == 1) {
                if (bud.exceeded()) return 0;
                ys = y;
                Int rk = r - k;
                Int lim = (m < rk ? m : rk);
                q = 1;
                for (Int i = 0; i < lim; ++i) {
                    y = rho_f<Int>(y, c, n);
                    Int diff = (x > y) ? (x - y) : (y - x);
                    if (diff != 0) q = (q * diff) % n;
                    if (bud.exceeded()) return 0;
                }
                g = gcd_int(q, n);
                k += lim;
            }
            r *= 2;
        }

        if (g == n) {
            do {
                if (bud.exceeded()) return 0;
                ys = rho_f<Int>(ys, c, n);
                Int diff = (x > ys) ? (x - ys) : (ys - x);
                g = gcd_int(diff, n);
            } while (g == 1);
        }

        if (g != 1 && g != n) return g;
    }
    return 0;
}

// --- Optimierung 2: Precomputed E8 Roots ---
// Wir speichern die normalisierten Richtungsvektoren
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

// Tetraeder-Spannungsanalyse: Q = [E, A, B, C] (Indizes 0,1,2,3).
// Vier Dreier-Anordnungen (Faces): Face i = Quadrupel ohne Index i → S_SK für dieses Triplett.
// Face 0 (A-B-C) ohne E, Face 1 (E-B-C) ohne A, Face 2 (E-A-C) ohne B, Face 3 (E-A-B) ohne C.
// Divergenz: Wenn Face 3 stark abweicht → C-Wand Warnsignal. Knotenpunkte: Kreuzungen = spannungsfrei.
// Shadow-Shift: Face mit max. Spannung bestimmt Phase (tx, ty) für E8-Interferenz.
template<class Int>
static void calculate_triplet_tensions(const std::array<Int, 4>& Q, double face_tensions[4]) {
    double log_Q[4];
    for (int i = 0; i < 4; ++i) log_Q[i] = log_approx(Q[i]);
    // Für jede Face: 3 Indizes (ohne i), S_SK = Summe |log_a - log_b| über die 3 Paare
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
    // Optimierte Berechnung ohne vector-allocs
    double log_Q[4];
    for (int i = 0; i < 4; ++i) log_Q[i] = log_approx(Q[i]);
    
    s_sk = 0.0;
    // Unrolled loops für 4 Elemente
    s_sk += std::fabs(log_Q[0] - log_Q[1]);
    s_sk += std::fabs(log_Q[0] - log_Q[2]);
    s_sk += std::fabs(log_Q[0] - log_Q[3]);
    s_sk += std::fabs(log_Q[1] - log_Q[2]);
    s_sk += std::fabs(log_Q[1] - log_Q[3]);
    s_sk += std::fabs(log_Q[2] - log_Q[3]);

    double sum_log = log_Q[0] + log_Q[1] + log_Q[2] + log_Q[3];
    double inv_sum = (sum_log > 0.0) ? (1.0 / sum_log) : 0.0;
    
    const double target = 1.0 / 2.61803398875; // 1/phi^2
    s_i = 0.0;
    for (int i = 0; i < 4; ++i) s_i += std::fabs((log_Q[i] * inv_sum) - target);
    
    if (s_i < 1e-9) s_i = 1e-9;
    r_q = s_sk / s_i;
}

// --- Schütte-GPS (Optimierte Suche) ---

template<class Int>
static Int find_nearest_class_prime(const Int& start, int direction, int want_mod12, int mr_rounds, std::mt19937_64& rng) {
    Int p = start;
    if (direction < 0) p--; else p++;
    
    int steps = 0;
    while (steps < 3000) { // Limit steps
        if (p < 2) return 0;
        
        // Fast Modulo check
        long r12 = (long)(p % 12);
        if (r12 == want_mod12) {
            // Erst Primorial check, dann MR
            if (gcd_int(p, GlobalConsts<Int>::PRIMORIAL) == 1) {
                if (mp::miller_rabin_test(p, mr_rounds, rng)) return p;
            }
        }
        
        if (direction < 0) p--; else p++;
        steps++;
    }
    return 0; // Nicht gefunden im Radius
}

// Shadow-Shift aus Face-Spannungen: welche Face hat maximale Spannung?
// Face 0 (A-B-C) ohne E → Shift 0; Face 1 ohne A → π/6; Face 2 ohne B → π/3; Face 3 ohne C → π/2.
// Liefert Winkel in Radiant für (tx, ty) = (cos(angle), sin(angle)).
static double get_triplet_based_phase_offset(const double face_tensions[4]) {
    int imax = 0;
    for (int i = 1; i < 4; ++i)
        if (face_tensions[i] > face_tensions[imax]) imax = i;
    return imax * (M_PI / 6.0); // 0, π/6, π/3, π/2
}

template<class Int>
static bool get_local_energy_state(const Int& n, int mr_rounds, std::mt19937_64& rng, 
                                   double& tuned_sens, double& tuned_scale_mult,
                                   double* out_phase_offset = nullptr,
                                   std::array<Int, 4>* out_Q_low = nullptr,
                                   std::array<Int, 4>* out_Q_high = nullptr) {
    if (n < 100) return false;

    // Wir brauchen Q_low und Q_high für Klassen 1, 5, 7, 11 (E, A, B, C)
    std::array<Int, 4> Q_low, Q_high;
    const int residues[4] = {1, 5, 7, 11};
    
    // Suche Ankerpunkte
    for (int i=0; i<4; ++i) {
        Q_low[i] = find_nearest_class_prime(n, -1, residues[i], mr_rounds, rng);
        Q_high[i] = find_nearest_class_prime(n, 1, residues[i], mr_rounds, rng);
        if (Q_low[i] == 0 || Q_high[i] == 0) return false; // Abbruch wenn GPS Signal fehlt
    }
    if (out_Q_low) *out_Q_low = Q_low;
    if (out_Q_high) *out_Q_high = Q_high;

    double sk_l, si_l, r_l, sk_h, si_h, r_h;
    calculate_metrics(Q_low, sk_l, si_l, r_l);
    calculate_metrics(Q_high, sk_h, si_h, r_h);
    
    double r_avg = 0.5 * (r_l + r_h);
    double si_avg = 0.5 * (si_l + si_h);
    
    // Tetraeder-Face-Spannungen (Dreier-Signale) für Shadow-Shift
    if (out_phase_offset) {
        double face_low[4], face_high[4], face_avg[4];
        calculate_triplet_tensions(Q_low, face_low);
        calculate_triplet_tensions(Q_high, face_high);
        for (int i = 0; i < 4; ++i) face_avg[i] = 0.5 * (face_low[i] + face_high[i]);
        *out_phase_offset = get_triplet_based_phase_offset(face_avg);
    }
    
    // Berechne Position ("Schütte-Position")
    Int max_low = 0; for(auto& x:Q_low) if(x>max_low) max_low=x;
    Int min_high = Q_high[0]; for(auto& x:Q_high) if(x<min_high) min_high=x;
    
    double d_low = log_approx(n - max_low);
    double d_high = log_approx(min_high - n);
    double pos_ratio = d_low / (d_low + d_high + 1e-9);

    // --- Tuning Logik ---
    unsigned bits = mp::msb(n);
    double base_sens = (bits > 100) ? 0.7 : 0.5;
    double base_scale = (bits > 100) ? 1.7 : 1.3;
    
    double energy_boost = std::log(r_avg + 1.0) / 6.0;
    
    tuned_sens = std::min(0.98, base_sens + energy_boost);
    tuned_scale_mult = base_scale * (1.0 + si_avg * 1.5);

    if (pos_ratio < 0.2 || pos_ratio > 0.8) {
        tuned_sens = std::min(0.98, tuned_sens + 0.05);
    }

    return true;
}

// --- E8 Factor Finder ---

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
        // 1. Trivial Checks
        if ((n & 1) == 0) { r.status=E8Status::HIT; r.factor=2; r.t_ms=ms_since(t0); return r; }
        
        // Primorial Check (Super schnell gegen kleine Faktoren)
        Int g = gcd_int(n, GlobalConsts<Int>::PRIMORIAL);
        if (g > 1) { r.status=E8Status::HIT; r.factor=g; r.t_ms=ms_since(t0); return r; }

        // 1b. Schnelle Trial Division VOR E8 (wie Rho/ECM) – kleine Faktoren sofort
        const uint32_t quick_trial_limit = 2000;
        for (uint32_t f = 101; f < quick_trial_limit && f < trial_bound; f += 2) {
            if (bud.exceeded()) { r.status=E8Status::TIMEOUT; r.t_ms=ms_since(t0); return r; }
            if (f % 3 == 0 || f % 5 == 0) continue;
            if (n % f == 0) { r.status=E8Status::HIT; r.factor=Int(f); r.t_ms=ms_since(t0); return r; }
        }

        if (bud.exceeded()) { r.status=E8Status::TIMEOUT; r.t_ms=ms_since(t0); return r; }

        // 2. Energetisches Tuning nur für große n (Bitlänge > 80); sonst feste Parameter
        double sens = 0.5;
        double scale_mult = 1.3;
        double phase_offset = 0.0;
        std::array<Int, 4> Q_low, Q_high;
        bool have_quad = false;
        unsigned bits = mp::msb(n);
        if (bits > 80)
            have_quad = get_local_energy_state(n, mr_rounds, rng, sens, scale_mult, &phase_offset, &Q_low, &Q_high);

        // 2b. Father (Walter-Geodäte): 16 Quadrupel-Produkte aus Q_low/Q_high, GCD mit n (Bamberg-Quadrupel)
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
            // Trippel-Raum (nur A, B, C): E fest (Index 0), 8 Produkte aus (A,B,C) low/high
            for (int mask = 0; mask < 8; ++mask) {
                if (bud.exceeded()) break;
                Int P = Q_low[0]; // E fest
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

        // 3. Alice (Morley-Geodäte) + Taurus (Kepler-Geodäte): E8 Injection mit mehreren Phasen und Skalen
        Int root_n = integer_sqrt(n);
        double base_scale = (double)bits; 
        int scale = std::max(3, (int)std::round(base_scale * scale_mult));
        
        const int y_offsets[5] = { 1, 2, 3, 5, 8 }; // Fibonacci
        
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

        // Alice (Morley-Geodäte): 3 Phasen (0, 2π/3, 4π/3) entlang des Phasenkreises
        const double MORLEY_STEP = 2.0 * M_PI / 3.0;
        for (int k = 0; k < 3; ++k) {
            if (bud.exceeded()) break;
            double tx = std::cos(phase_offset + k * MORLEY_STEP);
            double ty = std::sin(phase_offset + k * MORLEY_STEP);
            if (run_e8_injection(tx, ty, scale)) return r;
        }
        // Taurus (Kepler-Geodäte): 3 log-verteilte Skalen (1.0, ~sqrt(3), ~sqrt(4.5))
        const double kepler_scales[3] = { 1.0, 1.732, 2.121 }; // 1, sqrt(3), sqrt(4.5)
        for (int ks = 1; ks < 3 && ms_since(t0) < 3.0; ++ks) {
            if (bud.exceeded()) break;
            int scale_k = std::max(3, (int)std::round(scale * kepler_scales[ks]));
            double tx = std::cos(phase_offset);
            double ty = std::sin(phase_offset);
            if (run_e8_injection(tx, ty, scale_k)) return r;
        }
        
        // Fallback: Trial Division ab quick_trial_limit bis trial_bound
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

// --- ECM (Elliptic Curve Method) Stage-1 ---
// Lenstra ECM: Vergleichsfigur State-of-the-art für mittlere Faktoren.

template<class Int, class Wider>
static Int mod_add_ecm(const Int& a, const Int& b, const Int& n) { Int r = a + b; r %= n; return r; }
template<class Int, class Wider>
static Int mod_sub_ecm(const Int& a, const Int& b, const Int& n) {
    if (a >= b) return (a - b) % n;
    return (n - (b - a) % n) % n;
}
template<class Int, class Wider>
static Int mod_mul_ecm(const Int& a, const Int& b, const Int& n) { return Int((Wider(a) * Wider(b)) % Wider(n)); }

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
struct XZ { Int X; Int Z; };

template<class Int, class Wider>
static XZ<Int,Wider> xDBL(const XZ<Int,Wider>& P, const Int& A24, const Int& n) {
    Int t1 = mod_add_ecm<Int,Wider>(P.X, P.Z, n);
    Int t2 = mod_sub_ecm<Int,Wider>(P.X, P.Z, n);
    Int t1s = mod_mul_ecm<Int,Wider>(t1, t1, n);
    Int t2s = mod_mul_ecm<Int,Wider>(t2, t2, n);
    Int t3 = mod_sub_ecm<Int,Wider>(t1s, t2s, n);
    Int X2 = mod_mul_ecm<Int,Wider>(t1s, t2s, n);
    Int Z2 = mod_mul_ecm<Int,Wider>(t3, mod_add_ecm<Int,Wider>(t2s, mod_mul_ecm<Int,Wider>(A24, t3, n), n), n);
    return {X2, Z2};
}

template<class Int, class Wider>
static XZ<Int,Wider> xADD(const XZ<Int,Wider>& P, const XZ<Int,Wider>& Q, const XZ<Int,Wider>& D, const Int& n) {
    Int t1 = mod_add_ecm<Int,Wider>(P.X, P.Z, n), t2 = mod_sub_ecm<Int,Wider>(P.X, P.Z, n);
    Int t3 = mod_add_ecm<Int,Wider>(Q.X, Q.Z, n), t4 = mod_sub_ecm<Int,Wider>(Q.X, Q.Z, n);
    Int t5 = mod_mul_ecm<Int,Wider>(t1, t4, n), t6 = mod_mul_ecm<Int,Wider>(t2, t3, n);
    Int t7 = mod_add_ecm<Int,Wider>(t5, t6, n), t8 = mod_sub_ecm<Int,Wider>(t5, t6, n);
    Int X = mod_mul_ecm<Int,Wider>(D.Z, mod_mul_ecm<Int,Wider>(t7, t7, n), n);
    Int Z = mod_mul_ecm<Int,Wider>(D.X, mod_mul_ecm<Int,Wider>(t8, t8, n), n);
    return {X, Z};
}

template<class Int, class Wider>
static XZ<Int,Wider> xMUL(const XZ<Int,Wider>& P, const Int& k, const Int& A24, const Int& n) {
    XZ<Int,Wider> R0 = P, R1 = xDBL<Int,Wider>(P, A24, n), D = P;
    unsigned msb = mp::msb(k);
    for (int i = (int)msb - 1; i >= 0; --i) {
        if (mp::bit_test(k, i)) {
            R0 = xADD<Int,Wider>(R0, R1, D, n);
            R1 = xDBL<Int,Wider>(R1, A24, n);
        } else {
            R1 = xADD<Int,Wider>(R0, R1, D, n);
            R0 = xDBL<Int,Wider>(R0, A24, n);
        }
    }
    return R0;
}

enum class ECMStatus { HIT, MISS, TIMEOUT, ERROR };
static const char* ecm_status_name(ECMStatus s) {
    switch (s) { case ECMStatus::HIT: return "HIT"; case ECMStatus::MISS: return "MISS";
        case ECMStatus::TIMEOUT: return "TIMEOUT"; case ECMStatus::ERROR: return "ERROR"; default: return "?"; }
}

template<class Int>
struct ECMRun { ECMStatus status=ECMStatus::MISS; Int factor=0; double t_ms=0.0; int curves_tried=0; };

static std::vector<uint32_t> primes_up_to_B1(uint32_t N) {
    std::vector<bool> is_p(N+1, true);
    is_p[0]=is_p[1]=false;
    for (uint32_t i=2; i*i<=N; ++i) if (is_p[i]) for (uint32_t j=i*i; j<=N; j+=i) is_p[j]=false;
    std::vector<uint32_t> ps;
    for (uint32_t i=2; i<=N; ++i) if (is_p[i]) ps.push_back(i);
    return ps;
}

static const std::vector<uint32_t>& primes_cached_B1(uint32_t N) {
    static std::map<uint32_t, std::vector<uint32_t>> cache;
    auto it = cache.find(N);
    if (it != cache.end()) return it->second;
    auto v = primes_up_to_B1(N);
    cache[N] = std::move(v);
    return cache[N];
}

template<class Int, class Wider>
static ECMRun<Int> ecm_stage1_find_factor(const Int& n, std::mt19937_64& rng, Budget& bud,
                                         uint32_t B1 = 5000, int curves = 16) {
    ECMRun<Int> R;
    auto t0 = Clock::now();
    try {
        if (bud.exceeded()) { R.status=ECMStatus::TIMEOUT; R.t_ms=ms_since(t0); return R; }
        if (n < 4) { R.status=ECMStatus::MISS; R.t_ms=ms_since(t0); return R; }
        if ((n & 1) == 0) { R.status=ECMStatus::HIT; R.factor=2; R.t_ms=ms_since(t0); return R; }
        for (uint32_t p : SMALL_PRIMES) {
            if (bud.exceeded()) { R.status=ECMStatus::TIMEOUT; R.t_ms=ms_since(t0); return R; }
            if (n % p == 0) { R.status=ECMStatus::HIT; R.factor=Int(p); R.t_ms=ms_since(t0); return R; }
        }
        const auto& ps = primes_cached_B1(B1);
        std::uniform_int_distribution<uint64_t> dist64(6, std::numeric_limits<uint64_t>::max()-1);
        for (int cidx=0; cidx<curves; ++cidx) {
            if (bud.exceeded()) { R.status=ECMStatus::TIMEOUT; R.t_ms=ms_since(t0); return R; }
            R.curves_tried = cidx + 1;
            Int sigma = Int(dist64(rng)) % (n - 7) + 6;
            Int u = (sigma*sigma - 5) % n; if (u < 0) u += n;
            Int v = (Int(4) * sigma) % n;
            Int u2 = mod_mul_ecm<Int,Wider>(u, u, n), u3 = mod_mul_ecm<Int,Wider>(u2, u, n);
            Int v2 = mod_mul_ecm<Int,Wider>(v, v, n), v3 = mod_mul_ecm<Int,Wider>(v2, v, n);
            Int denom = mod_mul_ecm<Int,Wider>(Int(4), mod_mul_ecm<Int,Wider>(u3, v, n), n);
            auto [inv_denom, g] = inv_mod_with_gcd<Int>(denom, n);
            if (g != 1 && g != n) { R.status=ECMStatus::HIT; R.factor=g; R.t_ms=ms_since(t0); return R; }
            if (g == n) continue;
            Int vm_u = mod_sub_ecm<Int,Wider>(v, u, n);
            Int vm_u3 = mod_mul_ecm<Int,Wider>(mod_mul_ecm<Int,Wider>(vm_u, vm_u, n), vm_u, n);
            Int threeu_p_v = mod_add_ecm<Int,Wider>(mod_mul_ecm<Int,Wider>(Int(3), u, n), v, n);
            Int num = mod_mul_ecm<Int,Wider>(vm_u3, threeu_p_v, n);
            Int A = mod_sub_ecm<Int,Wider>(mod_mul_ecm<Int,Wider>(num, inv_denom, n), Int(2) % n, n);
            auto [inv4, g4] = inv_mod_with_gcd<Int>(Int(4), n);
            if (g4 != 1 && g4 != n) { R.status=ECMStatus::HIT; R.factor=g4; R.t_ms=ms_since(t0); return R; }
            Int A24 = mod_mul_ecm<Int,Wider>(mod_add_ecm<Int,Wider>(A, Int(2), n), inv4, n);
            XZ<Int,Wider> P{u3, v3};
            if (P.Z == 0) continue;
            XZ<Int,Wider> Q = P;
            for (uint32_t p : ps) {
                if (bud.exceeded()) { R.status=ECMStatus::TIMEOUT; R.t_ms=ms_since(t0); return R; }
                uint32_t pe = p;
                while (pe <= B1 / p) pe *= p;
                Q = xMUL<Int,Wider>(Q, Int(pe), A24, n);
                Int gz = gcd_int(Q.Z, n);
                if (gz != 1 && gz != n) { R.status=ECMStatus::HIT; R.factor=gz; R.t_ms=ms_since(t0); return R; }
                if (gz == n) break;
            }
        }
        R.status=ECMStatus::MISS;
        R.t_ms=ms_since(t0);
        return R;
    } catch (...) {
        R.status=ECMStatus::ERROR;
        R.t_ms=ms_since(t0);
        return R;
    }
}

// --- Main Loop & Args Parsing (gekürzt) ---

struct Args {
    unsigned bits=64;
    int count=10;
    uint64_t seed=42;
    std::string out_csv="optimized_res.csv";
};

// Kleine Primzahlen im Trial-Bereich [101, trial_bound), damit E8 HIT liefert
static const uint32_t EASY_SMALL_PRIMES[] = {
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173
};
static const int NUM_EASY = (int)(sizeof(EASY_SMALL_PRIMES) / sizeof(EASY_SMALL_PRIMES[0]));

int main(int argc, char** argv) {
    init_e8_roots(); // Precompute init
    
    using Int = mp::cpp_int;
    using Wider = mp::cpp_int;
    
    GlobalConsts<Int>::PRIMORIAL = get_primorial<Int>();
    
    const unsigned q_bits = 1000;  // Bitlänge des großen Faktors q; n hat dann ca. q_bits+7 Bits (cpp_int = beliebig)
    
    std::cout << "Optimierter E8 Benchmark (v9) gestartet..." << std::endl;
    std::cout << "Primorial size: " << to_dec(GlobalConsts<Int>::PRIMORIAL).size() << " digits." << std::endl;
    std::cout << "Test: n = kleiner Faktor * q mit q " << q_bits << " Bits (n ca. " << (q_bits+7) << " Bits)." << std::endl;
    std::cout << "Ziel: mindestens 10 HITs. Vergleich: E8/Schütte (opt. 2 Phasen) | Pollard Rho | ECM Stage-1.\n" << std::endl;
    
    std::mt19937_64 rng(42);
    const int min_hits = 10;
    const uint32_t trial_bound = 1000;
    const uint32_t ecm_B1 = 5000;
    const int ecm_curves = 16;
    int hits = 0;
    int run = 0;
    
    while (hits < min_hits) {
        run++;
        Budget b_e8(5000);
        Budget b_rho(5000);
        Budget b_ecm(5000);
        
        uint32_t small_p = EASY_SMALL_PRIMES[(run - 1) % NUM_EASY];
        Int q = random_probable_prime<Int>(q_bits, 10, rng);
        Int n = Int(small_p) * q;
        
        std::cout << "  " << run << ". n = " << small_p << " * q (" << mp::msb(n) << " bits)" << std::endl;
        
        // E8 / Schütte (optimiert: 2 Phasen)
        auto t_e8 = Clock::now();
        auto res = e8_find_factor<Int, Wider>(n, rng, b_e8, 10, trial_bound);
        double ms_e8 = ms_since(t_e8);
        std::cout << "      E8/Schütte:  " << e8_status_name(res.status) << ", " << std::fixed << std::setprecision(3) << ms_e8 << " ms";
        if (res.status == E8Status::HIT) {
            hits++;
            std::cout << "  Faktor: " << res.factor;
            // Semi-Prim-Erkennung (Struktur-Dokument Abschn. 7): n = f * (n/f) mit n/f prim?
            Int cofactor = n / res.factor;
            if (cofactor > 1 && cofactor != res.factor && is_probable_prime(cofactor, 5, rng))
                std::cout << "  (Semi-Prim)";
        }
        std::cout << std::endl;
        
        // Pollard Rho (Brent)
        auto t_rho = Clock::now();
        Int rho_factor = pollard_rho_brent<Int>(n, rng, b_rho);
        double ms_rho = ms_since(t_rho);
        std::cout << "      Pollard Rho: " << (rho_factor != 0 ? "HIT" : "MISS") << ", " << ms_rho << " ms";
        if (rho_factor != 0) std::cout << "  Faktor: " << rho_factor;
        std::cout << std::endl;
        
        // ECM Stage-1 (State-of-the-art)
        auto t_ecm = Clock::now();
        auto ecm_res = ecm_stage1_find_factor<Int, Wider>(n, rng, b_ecm, ecm_B1, ecm_curves);
        double ms_ecm = ms_since(t_ecm);
        std::cout << "      ECM (B1=" << ecm_B1 << "): " << ecm_status_name(ecm_res.status) << ", " << ms_ecm << " ms, curves=" << ecm_res.curves_tried;
        if (ecm_res.status == ECMStatus::HIT) std::cout << "  Faktor: " << ecm_res.factor;
        std::cout << "  (State-of-the-art)" << std::endl;
    }
    
    std::cout << "\nErgebnis: " << hits << " HITs in " << run << " Läufen (E8/Schütte). Vergleich: Pollard Rho, ECM siehe oben." << std::endl;
    return 0;
}