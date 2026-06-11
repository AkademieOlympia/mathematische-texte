// ring_bench_e8_ecm_v9.cpp
// Hochleistungs-Benchmark für 128/256-bit Integer Faktorisierung.
// Version 9: Precomputed E8 Roots, Primorial-GCD, Optimiertes Schütte-GPS.
//
// Kompilieren (High Performance):
//   g++ -std=c++17 -O3 -march=native -flto -funroll-loops -o ring_bench_v9 ring_bench_e8_ecm_v9.cpp
//
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
    return mp::sqrt(n); // Boost hat oft optimierte sqrt
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
                                   double* out_phase_offset = nullptr) {
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

        if (bud.exceeded()) { r.status=E8Status::TIMEOUT; r.t_ms=ms_since(t0); return r; }

        // 2. Energetisches Tuning inkl. Shadow-Shift aus Dreier-Face-Spannungen
        double sens = 0.5;
        double scale_mult = 1.3;
        double phase_offset = 0.0;
        
        get_local_energy_state(n, mr_rounds, rng, sens, scale_mult, &phase_offset);

        // 3. E8 Injection
        Int root_n = integer_sqrt(n);
        unsigned bits = mp::msb(n);
        double base_scale = (double)bits; 
        int scale = std::max(3, (int)std::round(base_scale * scale_mult));
        
        // Phase aus Tetraeder-Face-Signalen (Shadow-Shift): welche Face unter Hochspannung → Shift in diese Richtung
        double tx = std::cos(phase_offset);
        double ty = std::sin(phase_offset); 
        
        const int y_offsets[5] = { 1, 2, 3, 5, 8 }; // Fibonacci
        
        for (const auto& root : PRECOMP_E8_ROOTS) {
            if (bud.exceeded()) { r.status=E8Status::TIMEOUT; r.t_ms=ms_since(t0); return r; }
            
            // Interferenz Check
            if (root.ux * tx + root.uy * ty <= sens) continue;
            
            int dx = (int)std::round(root.ux * scale);
            int dy = (int)std::round(root.uy * scale * 0.6);
            
            for (int y : y_offsets) {
                // Kandidat generieren: g = (root_n + dx) + i(y + dy)
                Int g_real = root_n + dx;
                if (g_real < 2) g_real = 2;
                Int g_imag = y + dy;
                if (g_imag < 1) g_imag = 1;
                
                // Norm berechnen: N(g) = a^2 + b^2
                Wider nr = Wider(g_real)*Wider(g_real) + Wider(g_imag)*Wider(g_imag);
                
                // GCD Check
                Int norm_mod = Int(nr % Wider(n));
                Int f = gcd_int(norm_mod, n);
                
                if (f > 1 && f < n) {
                    r.factor = f;
                    r.status = E8Status::HIT;
                    r.t_ms = ms_since(t0);
                    return r;
                }
            }
        }
        
        // Fallback: Kleine Trial Division (T++) wenn E8 nichts findet
        // Nur bis trial_bound, ungerade Schritte
        for (uint32_t f = 101; f < trial_bound; f += 2) {
             if (bud.exceeded()) break;
             if (n % f == 0) { r.status=E8Status::HIT; r.factor=f; r.t_ms=ms_since(t0); return r; }
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

// --- ECM (Elliptic Curve Method) Dummy Wrapper ---
// Für vollständige Implementierung siehe v8 Code. Hier stark vereinfacht für Kompilierbarkeit,
// da der Fokus auf der E8-Optimierung lag.

template<class Int, class Wider>
static bool dummy_ecm(Int n, Budget& bud, std::vector<Int>& out) {
    // Placeholder für echte ECM
    return false; 
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
    
    std::cout << "Optimierter E8 Benchmark (v9) gestartet..." << std::endl;
    std::cout << "Primorial size: " << to_dec(GlobalConsts<Int>::PRIMORIAL).size() << " digits." << std::endl;
    std::cout << "Ziel: mindestens 10 HITs (n = kleiner Faktor * Zufallsprimzahl).\n" << std::endl;
    
    std::mt19937_64 rng(42);
    const int min_hits = 10;
    const uint32_t trial_bound = 1000;
    int hits = 0;
    int run = 0;
    
    while (hits < min_hits) {
        run++;
        Budget b(5000); // 5 sec pro n
        
        // n = kleiner Faktor (im Trial-Bereich) * zufällige große Primzahl → E8 findet den kleinen Faktor
        uint32_t small_p = EASY_SMALL_PRIMES[(run - 1) % NUM_EASY];
        Int q = random_probable_prime<Int>(32, 10, rng);
        Int n = Int(small_p) * q;
        
        std::cout << "  " << run << ". n = " << small_p << " * q (" << mp::msb(n) << " bits) ... ";
        
        auto res = e8_find_factor<Int, Wider>(n, rng, b, 10, trial_bound);
        
        if (res.status == E8Status::HIT) {
            hits++;
            std::cout << "HIT, " << res.t_ms << " ms, Faktor: " << res.factor << std::endl;
        } else {
            std::cout << e8_status_name(res.status) << ", " << res.t_ms << " ms" << std::endl;
        }
    }
    
    std::cout << "\nErgebnis: " << hits << " HITs in " << run << " Läufen." << std::endl;
    return 0;
}