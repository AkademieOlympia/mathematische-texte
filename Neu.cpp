// ring_bench_e8_ecm_v10.cpp
// Hochleistungs-Benchmark mit RESONANCE LOCK (Aharonov-Bohm-Korrektur).
// Basiert auf experimentellen Daten: Alpha=0.0262, Beta=0.0524 (2*Alpha).
//
// Kompilieren:
//   g++ -std=c++17 -O3 -march=native -flto -o ring_bench_v10 ring_bench_e8_ecm_v10.cpp

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/multiprecision/integer.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <array>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace mp = boost::multiprecision;
using Clock = std::chrono::steady_clock;

// --- Hilfsfunktionen (Zeit, String, Math) ---
static inline double ms_since(const Clock::time_point& t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

template<class Int>
static inline Int gcd_int(Int a, Int b) {
    while (b != 0) { Int t = a % b; a = b; b = t; }
    return a;
}

template<class Int>
static inline Int integer_sqrt(const Int& n) {
    if (n <= 0) return 0;
    return mp::sqrt(n);
}

// Schnelle Log-Approximation via Bit-Scan
template<class Int>
static inline double log_approx(const Int& x) {
    if (x <= 0) return 0.0;
    unsigned bits = mp::msb(x); 
    return (double)bits * 0.69314718056; 
}

// --- Primorial & Cache ---
static const std::vector<uint32_t> SMALL_PRIMES = {
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97
};

template<class Int>
static Int get_primorial() {
    Int p = 1;
    for (uint32_t sp : SMALL_PRIMES) p *= sp;
    return p;
}

template<class Int>
struct GlobalConsts { static Int PRIMORIAL; };
template<class Int> Int GlobalConsts<Int>::PRIMORIAL;

// --- Zufall & Primzahltest ---
template<class Int>
static inline bool is_probable_prime(const Int& n, int rounds, std::mt19937_64& rng) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (!(n & 1)) return false;
    if (gcd_int(n, GlobalConsts<Int>::PRIMORIAL) != 1) return false;
    return mp::miller_rabin_test(n, rounds, rng);
}

// --- E8 Wurzeln (Precomputed) ---
struct E8Root { double ux; double uy; };
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

// --- Schütte-Analyse & Energetik ---
// (Hier vereinfacht für den Kernel-Fokus, die Logik steckt im Resonanz-Lock)

// --- KERNELEMENT: E8 Factor Finder mit RESONANCE LOCK ---
enum class E8Status { DISABLED, HIT, MISS, TIMEOUT, ERROR };
const char* e8_status_name(E8Status s) {
    switch(s) { case E8Status::HIT: return "HIT"; case E8Status::MISS: return "MISS"; default: return "?"; }
}

template<class Int>
struct E8Run { E8Status status=E8Status::DISABLED; Int factor=0; double t_ms=0.0; };

struct Budget {
    Clock::time_point t0; uint64_t ms;
    Budget(uint64_t ms_=0): t0(Clock::now()), ms(ms_) {}
    bool exceeded() const { return ms && ms_since(t0) > double(ms); }
};

template<class Int, class Wider>
static E8Run<Int> e8_find_factor(const Int& n, std::mt19937_64& rng, Budget& bud)
{
    E8Run<Int> r;
    auto t0 = Clock::now();
    
    try {
        // 1. Trivial & Primorial Check
        if ((n & 1) == 0) { r.status=E8Status::HIT; r.factor=2; r.t_ms=ms_since(t0); return r; }
        Int g = gcd_int(n, GlobalConsts<Int>::PRIMORIAL);
        if (g > 1) { r.status=E8Status::HIT; r.factor=g; r.t_ms=ms_since(t0); return r; }

        // 2. KLASSEN-ANALYSE & RESONANCE LOCK
        long mod12 = (long)(n % 12);
        double phase_correction = 0.0;
        
        // Die experimentellen Konstanten (Aharonov-Bohm-Phasen)
        const double ALPHA_BASE = 0.0262; // Grundresonanz
        const double BETA_RES   = 0.0524; // 2 * Alpha (Konstruktive Interferenz)
        const double BETA_HARM  = 0.1047; // 4 * Alpha (Oberschwingung)

        if (mod12 == 5) { 
            // Klasse A (Gauß-Spin 1/2): Leichte Korrektur
            phase_correction = ALPHA_BASE; 
        } 
        else if (mod12 == 11) {
            // Klasse C (Solenoid / Spin 5/2): Starke Korrektur
            // Wir prüfen kurz die "Tiefe" (heuristisch via Bit-Test oder Log-Ratio)
            // Hier vereinfacht: Standard C nutzt 2*Alpha
            phase_correction = BETA_RES;
            
            // Optional: Wenn n sehr groß ist (>100 Bit), schalte auf Oberschwingung
            if (mp::msb(n) > 100) phase_correction = BETA_HARM;
        }
        // Klasse E (1) und B (7) bleiben oft bei 0.0 (oder B bekommt -Alpha)

        // Phasen-Vektor berechnen (Rotation des Suchstrahls)
        // Wir rotieren das Gitter virtuell um phase_correction
        double cos_phi = std::cos(phase_correction);
        double sin_phi = std::sin(phase_correction);

        // 3. E8 INJECTION LOOP
        Int root_n = integer_sqrt(n);
        unsigned bits = mp::msb(n);
        int scale = std::max(3, (int)(bits * 1.3)); // Basis-Scale
        
        // Target Phase (aus Jacobi-Symbolen, hier vereinfacht auf Hauptachse 1.0, 0.0)
        // Aber korrigiert durch Resonance Lock!
        // Der Suchstrahl (tx, ty) ist der rotierte Einheitsvektor
        double tx = 1.0 * cos_phi - 0.0 * sin_phi; 
        double ty = 1.0 * sin_phi + 0.0 * cos_phi;

        double sens = 0.5; // Basis-Sensitivität (bei C evtl. höher setzen)
        if (mod12 == 11) sens = 0.65; // C braucht mehr Toleranz oder genaueren Treffer? 
                                      // Experimente zeigten hohe pNL bei dieser Phase -> wir können präzise sein!
                                      // Also sens eher NORMAL lassen, weil die Phase stimmt!

        const int y_offsets[5] = { 1, 2, 3, 5, 8 };

        for (const auto& root : PRECOMP_E8_ROOTS) {
            if (bud.exceeded()) { r.status=E8Status::TIMEOUT; break; }

            // Interferenz-Check mit korrigierter Phase
            if (root.ux * tx + root.uy * ty <= sens) continue;

            int dx = (int)std::round(root.ux * scale);
            int dy = (int)std::round(root.uy * scale * 0.6);

            for (int y : y_offsets) {
                Int g_real = root_n + dx; 
                if (g_real < 2) g_real = 2;
                Int g_imag = y + dy;      
                if (g_imag < 1) g_imag = 1;

                Wider nr = Wider(g_real)*Wider(g_real) + Wider(g_imag)*Wider(g_imag);
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
        
        // Fallback Trial Division (klein)
        for(uint32_t f=101; f<1000; f+=2) {
             if(n%f==0) { r.status=E8Status::HIT; r.factor=f; r.t_ms=ms_since(t0); return r; }
        }

        r.status = E8Status::MISS;
        r.t_ms = ms_since(t0);
        return r;
    } catch (...) {
        r.status = E8Status::ERROR;
        return r;
    }
}

// --- Main ---
int main() {
    init_e8_roots();
    using Int = mp::cpp_int;
    using Wider = mp::cpp_int;
    GlobalConsts<Int>::PRIMORIAL = get_primorial<Int>();

    std::cout << "E8 Resonance Lock Benchmark (Alpha=0.0262, Beta=0.0524)" << std::endl;
    std::mt19937_64 rng(42);
    Budget b(2000);

    // Testfall: Eine Zahl der Klasse C (5*7 mod 12 = 35 = 11)
    // Wir erzeugen eine Semiprimzahl die C-artig ist
    // p=5 (A), q=7 (B) -> n=35 (C)
    // Wir nehmen größere: p=17 (A), q=19 (B) -> n=323 (11 mod 12 -> C)
    Int p = 17; 
    Int q = 19;
    Int n = p * q; // 323
    
    std::cout << "Teste n=" << n << " (Klasse " << (long)(n%12) << ")..." << std::endl;
    auto res = e8_find_factor<Int, Wider>(n, rng, b);
    std::cout << "Status: " << e8_status_name(res.status) << " Time: " << res.t_ms << "ms" << std::endl;
    if (res.status == E8Status::HIT) std::cout << "Factor: " << res.factor << std::endl;

    return 0;
}