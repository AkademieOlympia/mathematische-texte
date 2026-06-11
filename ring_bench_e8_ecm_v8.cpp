
// ring_bench_e8_ecm_v8.cpp
// Benchmark harness for 128/256-bit integers (Version 8: Schütte-Nachbarn).
// - Miller–Rabin timing (primality test)
// - E8 factor-finder mit energetischem Vor-Verstärker und Schütte-Einspannung (Q_low, Q_high)
// - ECM (Lenstra) stage-1 factorization (replaces Pollard–Rho)
//
// Schütte-GPS: n wird zwischen unterem (Q_low) und oberem (Q_high) Schütte-Quadrupel
// energetisch eingespannt; Compression/Position steuert Sensitivity/Scale.
//
// Build (macOS/Linux, Boost is header-only for multiprecision):
//   g++ -std=c++17 -O3 -march=native -I/opt/homebrew/include -o ring_bench_e8_ecm_v8 ring_bench_e8_ecm_v8.cpp
//
// Quick start:
//   ./ring_bench_e8_ecm_v8 --bits 128 --count 10 --seed 42 --save-dataset dataset.csv --e8 off --out res_no_e8.csv
//   ./ring_bench_e8_ecm_v8 --load-dataset dataset.csv --e8 on --e8-ms 50 --ecm-ms 3000 --out res_e8.csv
//
// Output:
//   - By default prints shortened numbers (10 digits front/back). Use --show-full to show all digits.
//   - Writes a CSV with full exact n and factors for post-analysis.

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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace mp = boost::multiprecision;
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

// small primes for quick screens
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
static inline Int integer_sqrt(Int n) {
    if (n <= 0) return 0;
    Int x = n, y = (x + 1) / 2;
    while (y < x) { x = y; y = (x + n / x) / 2; }
    return x;
}

// Log-Approximation für energetische Metriken (Heuristik, Double reicht)
template<class Int>
static inline double log_approx(const Int& x) {
    if (x <= 0) return 0.0;
    unsigned bits = 0;
    for (Int t = x; t != 0; t >>= 1) ++bits;
    return (bits - 0.5) * 0.69314718055994531; // ln(2)
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
    return mp::miller_rabin_test(n, rounds, rng);
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
        for (uint32_t p : SMALL_PRIMES) {
            if ((x % p) == 0 && x != p) { div = true; break; }
        }
        if (div) continue;
        if (is_probable_prime(x, rounds, rng)) return x;
    }
}

struct Budget {
    Clock::time_point t0;
    uint64_t ms; // 0 => unlimited
    Budget(uint64_t ms_=0): t0(Clock::now()), ms(ms_) {}
    inline bool exceeded() const {
        if (!ms) return false;
        return std::chrono::duration<double, std::milli>(Clock::now() - t0).count() > double(ms);
    }
};

// ---------------- Energetische Metriken (thermodynamische Ebene) ----------------
// S_SK = Schütte-Kuss-Spannung, S_I = Ikosaeder-Spannung, R(Q) = Impuls

template<class Int>
static void calculate_metrics(const Int Q[4], double& s_sk, double& s_i, double& r_q) {
    double log_Q[4];
    for (int i = 0; i < 4; ++i) log_Q[i] = log_approx(Q[i]);
    s_sk = 0.0;
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j)
            s_sk += std::fabs(log_Q[i] - log_Q[j]);
    double sum_log = log_Q[0] + log_Q[1] + log_Q[2] + log_Q[3];
    double theta[4];
    for (int i = 0; i < 4; ++i) theta[i] = (sum_log > 0.0) ? (log_Q[i] / sum_log) : 0.25;
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    const double target = 1.0 / (phi * phi); // ~0.3819
    s_i = 0.0;
    for (int i = 0; i < 4; ++i) s_i += std::fabs(theta[i] - target);
    if (s_i < 1e-9) s_i = 1e-9;
    r_q = s_sk / s_i;
}

// Anker-Primzahl nahe `near` mit Restklasse want_mod12 (1, 5, 7 oder 11)
template<class Int>
static Int find_anchor_prime(const Int& near, int want_mod12, int mr_rounds, std::mt19937_64& rng,
                             int max_steps = 2000) {
    if (near < 2) return 0;
    Int up = near, down = (near > 2) ? (near - 1) : Int(2);
    for (int step = 0; step < max_steps; ++step) {
        if (up > 0) {
            int r12 = int((up % 12).template convert_to<long>());
            if (r12 < 0) r12 += 12;
            if (r12 == want_mod12 && (up % 2) != 0 && (up % 3) != 0) {
                bool div = false;
                for (uint32_t p : SMALL_PRIMES) {
                    if (p < up && (up % p) == 0) { div = true; break; }
                }
                if (!div && is_probable_prime(up, mr_rounds, rng)) return up;
            }
            ++up;
        }
        if (down >= 2) {
            int r12 = int((down % 12).template convert_to<long>());
            if (r12 < 0) r12 += 12;
            if (r12 == want_mod12 && (down % 2) != 0 && (down % 3) != 0) {
                bool div = false;
                for (uint32_t p : SMALL_PRIMES) {
                    if (p < down && (down % p) == 0) { div = true; break; }
                }
                if (!div && is_probable_prime(down, mr_rounds, rng)) return down;
            }
            if (down <= 2) break;
            --down;
        }
    }
    return 0;
}

// --------------- Schütte-Nachbarn (lokales GPS, Modul 12) ---------------
// Q_low: unteres Schütte-Quadrupel (nächste Primzahlen unter n in 1,5,7,11 mod 12)
// Q_high: oberes Schütte-Quadrupel (nächste Primzahlen über n in 1,5,7,11 mod 12)
// Reihenfolge: E(1), A(5), B(7), C(11)

template<class Int>
static Int find_previous_prime_in_class(const Int& n, int want_mod12, int mr_rounds, std::mt19937_64& rng,
                                        int max_steps = 5000) {
    if (n <= 2) return 0;
    Int p = n - 1;
    for (int step = 0; step < max_steps && p >= 2; ++step, --p) {
        int r12 = int((p % 12).template convert_to<long>());
        if (r12 < 0) r12 += 12;
        if (r12 != want_mod12 || (p % 2) == 0 || (p % 3) == 0) continue;
        bool div = false;
        for (uint32_t q : SMALL_PRIMES) {
            if (q < p && (p % q) == 0) { div = true; break; }
        }
        if (!div && is_probable_prime(p, mr_rounds, rng)) return p;
    }
    return 0;
}

template<class Int>
static Int find_next_prime_in_class(const Int& n, int want_mod12, int mr_rounds, std::mt19937_64& rng,
                                     int max_steps = 5000) {
    Int p = n + 1;
    for (int step = 0; step < max_steps; ++step, ++p) {
        int r12 = int((p % 12).template convert_to<long>());
        if (r12 < 0) r12 += 12;
        if (r12 != want_mod12 || (p % 2) == 0 || (p % 3) == 0) continue;
        bool div = false;
        for (uint32_t q : SMALL_PRIMES) {
            if (q < p && (p % q) == 0) { div = true; break; }
        }
        if (!div && is_probable_prime(p, mr_rounds, rng)) return p;
    }
    return 0;
}

template<class Int>
static bool get_schuette_neighbors(const Int& n, int mr_rounds, std::mt19937_64& rng,
                                    Int Q_low[4], Int Q_high[4]) {
    if (n < 2) return false;
    const int residues[4] = { 1, 5, 7, 11 };
    for (int i = 0; i < 4; ++i) {
        Q_low[i]  = find_previous_prime_in_class(n, residues[i], mr_rounds, rng);
        Q_high[i] = find_next_prime_in_class(n, residues[i], mr_rounds, rng);
        if (Q_low[i] == 0 || Q_high[i] == 0) return false;
    }
    return true;
}

// Position von n zwischen Q_low und Q_high: 0 = nahe Q_low („gepresst“), 0.5 = Mitte („freier Fall“), 1 = nahe Q_high
template<class Int>
static double schuette_position(const Int& n, const Int Q_low[4], const Int Q_high[4]) {
    Int max_low = Q_low[0];
    Int min_high = Q_high[0];
    for (int i = 1; i < 4; ++i) {
        if (Q_low[i] > max_low) max_low = Q_low[i];
        if (Q_high[i] < min_high) min_high = Q_high[i];
    }
    Int d_low = n - max_low;   // Abstand n über dem oberen Rand von Q_low
    Int d_high = min_high - n; // Abstand n unter dem unteren Rand von Q_high
    if (d_low <= 0) return 0.0;
    if (d_high <= 0) return 1.0;
    double denom = double(d_low + d_high);
    if (denom < 1e-30) return 0.5;
    return double(d_high) / denom;
}

template<class Int>
static bool get_local_energy_state(const Int& n, int mr_rounds, std::mt19937_64& rng,
                                   double& s_sk, double& s_i, double& r_q,
                                   double* out_position = nullptr) {
    if (n < 4) return false;
    Int Q_low[4], Q_high[4];
    if (!get_schuette_neighbors(n, mr_rounds, rng, Q_low, Q_high)) return false;
    double s_sk_low, s_i_low, r_low, s_sk_high, s_i_high, r_high;
    calculate_metrics<Int>(Q_low, s_sk_low, s_i_low, r_low);
    calculate_metrics<Int>(Q_high, s_sk_high, s_i_high, r_high);
    s_sk = 0.5 * (s_sk_low + s_sk_high);
    s_i  = 0.5 * (s_i_low + s_i_high);
    if (s_i < 1e-9) s_i = 1e-9;
    r_q = s_sk / s_i;
    if (out_position) *out_position = schuette_position(n, Q_low, Q_high);
    return true;
}

// Sensitivity und Scale-Multiplikator aus energetischem Zustand (inkl. Schütte-Position)
template<class Int>
static void get_energy_tuned_params(const Int& n, int mr_rounds, std::mt19937_64& rng,
                                    double& tuned_sens, double& tuned_scale_mult) {
    unsigned bits = 0;
    for (Int t = n; t != 0; t >>= 1) ++bits;
    double base_sens = 0.5, base_scale_mult = 1.3;
    if (bits > 80)  { base_sens = 0.6;  base_scale_mult = 1.5; }
    if (bits > 100) { base_sens = 0.7;  base_scale_mult = 1.7; }
    tuned_sens = base_sens;
    tuned_scale_mult = base_scale_mult;
    double s_sk, s_i, r_q, position = 0.5;
    if (!get_local_energy_state(n, mr_rounds, rng, s_sk, s_i, r_q, &position)) return;
    double energy_factor = std::log(r_q + 1.0) / 5.0;
    tuned_sens = std::min(0.95, base_sens + energy_factor * 0.1);
    tuned_scale_mult = base_scale_mult * (1.0 + s_i * 2.0);
    // Schütte-Position: „gepresst“ (position < 0.3) → engere Sensitivity; „freier Fall“ (0.4–0.6) → neutral
    if (position < 0.3) tuned_sens = std::min(0.95, tuned_sens + 0.05);
}

// ---------------- E8 hook (mit energetischem Vor-Verstärker, v8: Schütte-GPS) ----------------
// Thermodynamische Ebene: S_SK (Schütte-Kuss-Spannung), S_I (Ikosaeder-Spannung),
// R(Q)=S_SK/S_I (Impuls). Hoher Impuls → engere Sensitivity; hohe S_I → größerer Scale.
// v8: n wird zwischen Q_low und Q_high eingespannt; Position (gepresst / freier Fall)
// kalibriert den Startpunkt für die Interferenz (engere Sensitivity wenn „gepresst“).
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
struct E8Run { E8Status status=E8Status::DISABLED; Int factor=0; double t_ms=0.0; };

// E8EnergyCracker: energetische Heuristik (Sensitivity/Scale) + E8-Injection + Trial
template<class Int, class Wider>
static E8Run<Int> e8_find_factor(const Int& n, std::mt19937_64& rng, Budget& bud,
                                 int mr_rounds, uint32_t trial_bound)
{
    E8Run<Int> r;
    auto t0 = Clock::now();
    try {
        if (bud.exceeded()) { r.status = E8Status::TIMEOUT; r.t_ms = ms_since(t0); return r; }
        if (n < 4) { r.status = E8Status::MISS; r.t_ms = ms_since(t0); return r; }
        if ((n & 1) == 0) { r.status = E8Status::HIT; r.factor = 2; r.t_ms = ms_since(t0); return r; }
        for (uint32_t p : SMALL_PRIMES) {
            if (bud.exceeded()) { r.status = E8Status::TIMEOUT; r.t_ms = ms_since(t0); return r; }
            if (n % p == 0) { r.status = E8Status::HIT; r.factor = Int(p); r.t_ms = ms_since(t0); return r; }
        }

        // Energetisches Tuning: Sensitivity und Scale aus lokalem BM-Quadruppel
        double sens = 0.5, scale_mult = 1.3;
        get_energy_tuned_params(n, mr_rounds, rng, sens, scale_mult);

        // E8-Injection: 240 Wurzeln, gefiltert nach Impuls, dann Norm-GCD
        Int root_n = integer_sqrt(n);
        if (root_n >= 2) {
            unsigned bits = 0;
            for (Int t = n; t != 0; t >>= 1) ++bits;
            double base_scale = (bits > 0) ? (bits - 0.5) : 1.0;
            int scale = std::max(3, (int)std::round(base_scale * scale_mult));
            double tx = 1.0, ty = 0.0; // Platzhalter Phase; v4_signature könnte (tx,ty) setzen
            static const int E8_N = 240;
            static const double PI2 = 2.0 * M_PI;
            const int y_offsets[5] = { 1, 2, 3, 5, 8 };
            for (int i = 0; i < E8_N; ++i) {
                if (bud.exceeded()) { r.status = E8Status::TIMEOUT; r.t_ms = ms_since(t0); return r; }
                double ux = std::cos(PI2 * i / E8_N), uy = std::sin(PI2 * i / E8_N);
                if (ux * tx + uy * ty <= sens) continue;
                int dx = (int)std::round(ux * scale), dy = (int)std::round(uy * scale * 0.6);
                for (int y : y_offsets) {
                    Int g_real = root_n + dx;
                    if (g_real < 2) g_real = 2;
                    Int g_imag = Int(y + dy);
                    if (g_imag < 1) g_imag = 1;
                    Wider nr = Wider(g_real) * Wider(g_real) + Wider(g_imag) * Wider(g_imag);
                    Int norm_mod = Int(nr % Wider(n));
                    Int g = gcd_int(norm_mod, n);
                    if (g > 1 && g < n) {
                        r.factor = g;
                        r.status = E8Status::HIT;
                        r.t_ms = ms_since(t0);
                        return r;
                    }
                }
            }
        }

        // Trial division (T++)
        for (uint32_t f = 1001; f <= trial_bound; f += 2) {
            if (bud.exceeded()) { r.status = E8Status::TIMEOUT; r.t_ms = ms_since(t0); return r; }
            if (f % 3 == 0 || f % 5 == 0) continue;
            if (n % f == 0) { r.status = E8Status::HIT; r.factor = Int(f); r.t_ms = ms_since(t0); return r; }
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

// ---------------- ECM stage-1 ----------------
template<class Int, class Wider>
static inline Int mod_add(const Int& a, const Int& b, const Int& n) { Int r = a + b; r %= n; return r; }
template<class Int, class Wider>
static inline Int mod_sub(const Int& a, const Int& b, const Int& n) {
    Int r;
    if (a >= b) r = a - b;
    else {
        Int t = (b - a) % n;
        r = n - t;
    }
    r %= n;
    return r;
}

template<class Int, class Wider>
static inline Int mod_mul(const Int& a, const Int& b, const Int& n) { return Int((Wider(a) * Wider(b)) % Wider(n)); }

template<class Int>
static std::pair<Int, Int> inv_mod_with_gcd(Int a, Int n) {
    a %= n;
    if (a < 0) a += n;
    Int t = 0, newt = 1;
    Int r = n, newr = a;
    while (newr != 0) {
        Int q = r / newr;
        Int tmp = t - q * newt; t = newt; newt = tmp;
        tmp = r - q * newr; r = newr; newr = tmp;
    }
    if (t < 0) t += n;
    return {t, r};
}

template<class Int, class Wider>
struct XZ { Int X; Int Z; };

template<class Int, class Wider>
static inline XZ<Int,Wider> xDBL(const XZ<Int,Wider>& P, const Int& A24, const Int& n) {
    Int X = P.X, Z = P.Z;
    Int t1 = mod_add<Int,Wider>(X, Z, n);
    Int t2 = mod_sub<Int,Wider>(X, Z, n);
    Int t1s = mod_mul<Int,Wider>(t1, t1, n);
    Int t2s = mod_mul<Int,Wider>(t2, t2, n);
    Int t3 = mod_sub<Int,Wider>(t1s, t2s, n);
    Int X2 = mod_mul<Int,Wider>(t1s, t2s, n);
    Int Z2 = mod_mul<Int,Wider>(t3, mod_add<Int,Wider>(t2s, mod_mul<Int,Wider>(A24, t3, n), n), n);
    return {X2, Z2};
}

template<class Int, class Wider>
static inline XZ<Int,Wider> xADD(const XZ<Int,Wider>& P, const XZ<Int,Wider>& Q, const XZ<Int,Wider>& D, const Int& n) {
    Int t1 = mod_add<Int,Wider>(P.X, P.Z, n);
    Int t2 = mod_sub<Int,Wider>(P.X, P.Z, n);
    Int t3 = mod_add<Int,Wider>(Q.X, Q.Z, n);
    Int t4 = mod_sub<Int,Wider>(Q.X, Q.Z, n);
    Int t5 = mod_mul<Int,Wider>(t1, t4, n);
    Int t6 = mod_mul<Int,Wider>(t2, t3, n);
    Int t7 = mod_add<Int,Wider>(t5, t6, n);
    Int t8 = mod_sub<Int,Wider>(t5, t6, n);
    Int X = mod_mul<Int,Wider>(D.Z, mod_mul<Int,Wider>(t7, t7, n), n);
    Int Z = mod_mul<Int,Wider>(D.X, mod_mul<Int,Wider>(t8, t8, n), n);
    return {X, Z};
}

template<class Int, class Wider>
static XZ<Int,Wider> xMUL(const XZ<Int,Wider>& P, const Int& k, const Int& A24, const Int& n) {
    XZ<Int,Wider> R0 = P;
    XZ<Int,Wider> R1 = xDBL<Int,Wider>(P, A24, n);
    XZ<Int,Wider> D  = P;
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
static inline const char* ecm_status_name(ECMStatus s) {
    switch (s) {
        case ECMStatus::HIT: return "HIT";
        case ECMStatus::MISS: return "MISS";
        case ECMStatus::TIMEOUT: return "TIMEOUT";
        case ECMStatus::ERROR: return "ERROR";
        default: return "?";
    }
}

template<class Int>
struct ECMRun { ECMStatus status=ECMStatus::MISS; Int factor=0; double t_ms=0.0; int curves_tried=0; };

static std::vector<uint32_t> primes_up_to(uint32_t N) {
    std::vector<bool> is_p(N+1, true);
    is_p[0]=is_p[1]=false;
    for (uint32_t i=2;i*i<=N;++i) if (is_p[i]) for (uint32_t j=i*i;j<=N;j+=i) is_p[j]=false;
    std::vector<uint32_t> ps;
    for (uint32_t i=2;i<=N;++i) if (is_p[i]) ps.push_back(i);
    return ps;
}

static const std::vector<uint32_t>& primes_cached(uint32_t N) {
    // simple in-process cache (good enough for benchmarking)
    static std::map<uint32_t, std::vector<uint32_t>> cache;
    auto it = cache.find(N);
    if (it != cache.end()) return it->second;
    auto v = primes_up_to(N);
    auto res = cache.emplace(N, std::move(v));
    return res.first->second;
}

template<class Int, class Wider>
static ECMRun<Int> ecm_stage1_find_factor(const Int& n, std::mt19937_64& rng, Budget& bud,
                                         uint32_t B1 = 50000, int curves = 64)
{
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

        const auto& ps = primes_cached(B1);
        std::uniform_int_distribution<uint64_t> dist64(6, std::numeric_limits<uint64_t>::max()-1);

        for (int cidx=0; cidx<curves; ++cidx) {
            if (bud.exceeded()) { R.status=ECMStatus::TIMEOUT; R.t_ms=ms_since(t0); return R; }
            R.curves_tried = cidx + 1;

            Int sigma = Int(dist64(rng)) % (n - 7) + 6;
            Int u = (sigma*sigma - 5) % n; if (u < 0) u += n;
            Int v = (Int(4) * sigma) % n;

            Int u2 = mod_mul<Int,Wider>(u, u, n);
            Int u3 = mod_mul<Int,Wider>(u2, u, n);
            Int v2 = mod_mul<Int,Wider>(v, v, n);
            Int v3 = mod_mul<Int,Wider>(v2, v, n);

            Int denom = mod_mul<Int,Wider>(Int(4), mod_mul<Int,Wider>(u3, v, n), n);
            auto [inv_denom, g] = inv_mod_with_gcd<Int>(denom, n);
            if (g != 1 && g != n) { R.status=ECMStatus::HIT; R.factor=g; R.t_ms=ms_since(t0); return R; }
            if (g == n) continue;

            Int vm_u = mod_sub<Int,Wider>(v, u, n);
            Int vm_u2 = mod_mul<Int,Wider>(vm_u, vm_u, n);
            Int vm_u3 = mod_mul<Int,Wider>(vm_u2, vm_u, n);
            Int threeu_p_v = mod_add<Int,Wider>(mod_mul<Int,Wider>(Int(3), u, n), v, n);
            Int num = mod_mul<Int,Wider>(vm_u3, threeu_p_v, n);

            Int A = mod_sub<Int,Wider>(mod_mul<Int,Wider>(num, inv_denom, n), Int(2) % n, n);

            auto [inv4, g4] = inv_mod_with_gcd<Int>(Int(4), n);
            if (g4 != 1 && g4 != n) { R.status=ECMStatus::HIT; R.factor=g4; R.t_ms=ms_since(t0); return R; }
            Int A24 = mod_mul<Int,Wider>(mod_add<Int,Wider>(A, Int(2), n), inv4, n);

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

template<class Int, class Wider>
static bool factorize_ecm(Int n, int mr_rounds, std::mt19937_64& rng,
                          uint32_t B1, int curves_per_round,
                          uint32_t trial_bound,
                          Budget& bud, std::vector<Int>& out,
                          bool& timed_out)
{
    timed_out = false;
    if (n == 1) return true;
    if (bud.exceeded()) { timed_out = true; return false; }

    if (is_probable_prime(n, mr_rounds, rng)) { out.push_back(n); return true; }

    for (uint32_t p : SMALL_PRIMES) {
        if (bud.exceeded()) { timed_out = true; return false; }
        if (n % p == 0) {
            out.push_back(Int(p));
            return factorize_ecm<Int,Wider>(n / Int(p), mr_rounds, rng, B1, curves_per_round, trial_bound, bud, out, timed_out);
        }
    }

    // Trial division (T): same logic as e8_find_factor, called inside factorizer
    for (uint32_t f = 1001; f <= trial_bound; f += 2) {
        if (bud.exceeded()) { timed_out = true; return false; }
        if (f % 3 == 0 || f % 5 == 0) continue;
        if (n % f == 0) {
            out.push_back(Int(f));
            return factorize_ecm<Int,Wider>(n / Int(f), mr_rounds, rng, B1, curves_per_round, trial_bound, bud, out, timed_out);
        }
    }

    for (int round=0; round<8; ++round) {
        if (bud.exceeded()) { timed_out = true; return false; }
        uint32_t B1r = B1 * (round < 3 ? 1u : (round < 5 ? 2u : 4u));
        int curvesr = curves_per_round * (round < 3 ? 1 : (round < 5 ? 2 : 4));

        auto ecmr = ecm_stage1_find_factor<Int,Wider>(n, rng, bud, B1r, curvesr);
        if (ecmr.status == ECMStatus::TIMEOUT) { timed_out = true; return false; }
        if (ecmr.status == ECMStatus::HIT && ecmr.factor > 1 && ecmr.factor < n && (n % ecmr.factor) == 0) {
            bool to1=false,to2=false;
            bool ok1 = factorize_ecm<Int,Wider>(ecmr.factor, mr_rounds, rng, B1, curves_per_round, trial_bound, bud, out, to1);
            bool ok2 = factorize_ecm<Int,Wider>(n/ecmr.factor, mr_rounds, rng, B1, curves_per_round, trial_bound, bud, out, to2);
            timed_out = to1 || to2;
            return ok1 && ok2 && !timed_out;
        }
    }
    return false;
}

// ---------------- Dataset ----------------
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
    if (!std::getline(f,line)) return false; // header
    while (std::getline(f,line)) {
        line=trim(line); if (line.empty()) continue;
        size_t p1=line.find(','); if (p1==std::string::npos) continue;
        size_t p2=line.find(',',p1+1); if (p2==std::string::npos) continue;
        size_t p3=line.find(',',p2+1); if (p3==std::string::npos) continue;
        std::string k=trim(line.substr(0,p1));
        std::string bits_s=trim(line.substr(p1+1,p2-p1-1));
        std::string n_s=trim(line.substr(p2+1, p3-p2-1));
        std::string t_s=trim(line.substr(p3+1));

        Sample<Int> s; s.kind=kind_from(k); s.bits=(unsigned)std::stoul(bits_s);
        mp::cpp_int nc(n_s); s.n=Int(nc);

        s.truth.clear();
        if (!t_s.empty()) {
            std::stringstream ss(t_s);
            std::string part;
            while (std::getline(ss,part,'*')) {
                part=trim(part); if (part.empty()) continue;
                mp::cpp_int fc(part); s.truth.push_back(Int(fc));
            }
            std::sort(s.truth.begin(), s.truth.end());
        }
        ds_out.push_back(std::move(s));
    }
    return true;
}

// ---------------- Args ----------------
struct Args {
    bool generate=true;
    unsigned bits=128;
    int count=10;          // per kind
    uint64_t seed=42;
    int mr_rounds=16;
    bool include_bamberg=true;

    std::string save_dataset;
    std::string load_dataset;

    bool show_full=false;
    int short_digits=10;

    bool e8_on=false;
    uint64_t e8_ms=50;
    uint64_t ecm_ms=3000;

    uint32_t ecm_B1=50000;
    int ecm_curves=64;
    uint32_t trial_bound=200000;

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
        else if (k=="--e8") a.e8_on=(need(k)=="on");
        else if (k=="--e8-ms") a.e8_ms=(uint64_t)std::stoull(need(k));
        else if (k=="--ecm-ms") a.ecm_ms=(uint64_t)std::stoull(need(k));
        else if (k=="--ecm-b1") a.ecm_B1=(uint32_t)std::stoul(need(k));
        else if (k=="--ecm-curves") a.ecm_curves=std::stoi(need(k));
        else if (k=="--trial-bound") a.trial_bound=(uint32_t)std::stoul(need(k));
        else if (k=="--out") a.out_csv=need(k);
        else if (k=="--help" || k=="-h") {
            std::cout <<
              "Usage: ./ring_bench_e8_ecm_v8 [options]\n\n"
              "Dataset:\n"
              "  --bits N --count N --seed S --bamberg on|off\n"
              "  --save-dataset file.csv\n"
              "  --load-dataset file.csv\n\n"
              "Bench:\n"
              "  --mr-rounds R\n"
              "  --e8 on|off --e8-ms T\n"
              "  --ecm-ms T --ecm-b1 B1 --ecm-curves C\n"
              "  --trial-bound T (T++ inside factorizer; default 200000)\n\n"
              "Output:\n"
              "  --out results.csv --short-digits D --show-full\n" << std::flush;
            std::exit(0);
        } else {
            throw std::runtime_error("unknown arg: " + k);
        }
    }
    if (a.count<=0) throw std::runtime_error("--count must be >0");
    if (a.short_digits<1) a.short_digits=1;
    return a;
}

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
    std::cout << "ring_bench_e8_ecm_v8: starting...\n" << std::flush;
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
        int total = (int)(kinds.size() * args.count);
        std::cout << "Generating dataset (" << total << " samples)...\n" << std::flush;
        for (auto k : kinds) for (int i=0;i<args.count;++i) ds.push_back(make_sample<Int>(k,args.bits,args.mr_rounds,rng));
        if (!args.save_dataset.empty()) {
            if (save_dataset_csv<Int>(args.save_dataset, ds)) std::cout << "Saved dataset: " << args.save_dataset << " (rows="<<ds.size()<<")\n";
            else std::cerr << "WARNING: could not save dataset\n";
        }
    }

    std::ofstream out(args.out_csv);
    if (!out) { std::cerr << "ERROR: cannot write " << args.out_csv << "\n"; return 3; }

    out << "i,kind,bits,n,truth,MR_ms,MR_is_prime,"
           "E8_ms,E8_status,E8_factor,"
           "ECM_ms,ECM_status,ECM_factors,ok\n";

    std::cout << "\n=== FULL FACTORIZATION BENCH (ECM baseline) ===\n"
              << "samples="<<ds.size()<<" bits="<<args.bits<<" mr_rounds="<<args.mr_rounds
              << " e8="<<(args.e8_on?"on":"off")<<" e8_ms="<<args.e8_ms
              << " ecm_ms="<<args.ecm_ms<<" ecm_B1="<<args.ecm_B1<<" ecm_curves="<<args.ecm_curves
              << " out="<<args.out_csv<<"\n\n" << std::flush;

    std::cout << "i  kind         MR(ms)  E8(ms)  E8_stat   ECM(ms)  ECM_stat  ok  n  truth  found\n"
              << "---------------------------------------------------------------------------------------------------------------------\n" << std::flush;

    struct Agg{int total=0, ok=0, e8_hit=0, ecm_to=0; std::vector<double> mr,e8,ecm,tot;};
    std::map<std::string,Agg> agg;

    for (size_t i=0;i<ds.size();++i) {
        const auto& s = ds[i];
        std::vector<Int> truth=s.truth; std::sort(truth.begin(), truth.end());

        auto t_mr0=Clock::now();
        bool mr_prime = is_probable_prime(s.n,args.mr_rounds,rng);
        double mr_ms = ms_since(t_mr0);

        if (mr_prime) {
            bool ok = (truth.size()==1 && truth[0]==s.n);
            out << (i+1)<<","<<kind_name(s.kind)<<","<<s.bits<<","<<to_dec(s.n)<<","<<join_factors_full(truth)<<","
                << std::fixed<<std::setprecision(6)<<mr_ms<<",1,0.0,DISABLED,,0.0,SKIP,,"<<(ok?"1":"0")<<"\n";

            std::cout<<std::setw(2)<<(i+1)<<"  "<<std::setw(11)<<kind_name(s.kind)<<"  "
                     <<std::fixed<<std::setprecision(3)<<std::setw(6)<<mr_ms<<"  "
                     <<std::setw(6)<<0.0<<"  "<<std::setw(7)<<"-"<<"  "
                     <<std::setw(7)<<0.0<<"  "<<std::setw(8)<<"SKIP"<<"  "
                     <<(ok?"Y ":"N ")<<short_dec(s.n,args.short_digits,args.show_full)<<"  "
                     <<join_factors_short(truth,args.short_digits,args.show_full)<<"  "
                     <<"PRIME:"<<short_dec(s.n,args.short_digits,args.show_full)<<"\n" << std::flush;

            auto& A=agg[std::string(kind_name(s.kind))];
            A.total++; A.ok += ok?1:0;
            A.mr.push_back(mr_ms); A.e8.push_back(0.0); A.ecm.push_back(0.0); A.tot.push_back(mr_ms);
            continue;
        }

        double e8_ms=0.0; std::string e8_stat="DISABLED"; Int e8_factor=0; bool e8_hit=false;
        if (args.e8_on) {
            Budget b_e8(args.e8_ms);
            auto e8r = e8_find_factor<Int,Wider>(s.n, rng, b_e8, args.mr_rounds, args.trial_bound);
            e8_ms = e8r.t_ms;
            e8_stat = e8_status_name(e8r.status);
            e8_factor = e8r.factor;
            e8_hit = (e8r.status==E8Status::HIT && e8_factor>1 && e8_factor<s.n && (s.n % e8_factor)==0);
        }

        Budget b_ecm(args.ecm_ms);
        std::vector<Int> found;
        bool to=false;
        bool okfact=false;
        std::string ecm_stat="MISS";

        auto t_ecm0=Clock::now();
        if (e8_hit) {
            std::vector<Int> a,b;
            bool to1=false,to2=false;
            bool ok1 = factorize_ecm<Int,Wider>(e8_factor, args.mr_rounds, rng, args.ecm_B1, args.ecm_curves, args.trial_bound, b_ecm, a, to1);
            bool ok2 = factorize_ecm<Int,Wider>(s.n/e8_factor, args.mr_rounds, rng, args.ecm_B1, args.ecm_curves, args.trial_bound, b_ecm, b, to2);
            to = to1 || to2;
            if (!to && ok1 && ok2) {
                found.insert(found.end(), a.begin(), a.end());
                found.insert(found.end(), b.begin(), b.end());
                std::sort(found.begin(), found.end());
                okfact=true;
                ecm_stat="OK";
            } else ecm_stat = to ? "TIMEOUT" : "FAIL";
        } else {
            okfact = factorize_ecm<Int,Wider>(s.n, args.mr_rounds, rng, args.ecm_B1, args.ecm_curves, args.trial_bound, b_ecm, found, to);
            std::sort(found.begin(), found.end());
            ecm_stat = to ? "TIMEOUT" : (okfact ? "OK" : "MISS");
        }
        double ecm_ms = ms_since(t_ecm0);

        bool ok = (ecm_stat=="OK" && found==truth);

        out << (i+1)<<","<<kind_name(s.kind)<<","<<s.bits<<","<<to_dec(s.n)<<","<<join_factors_full(truth)<<","
            << std::fixed<<std::setprecision(6)<<mr_ms<<",0,"
            << e8_ms<<","<<e8_stat<<","<<(e8_hit?to_dec(e8_factor):"")<<","
            << ecm_ms<<","<<ecm_stat<<","<<join_factors_full(found)<<","<<(ok?"1":"0")<<"\n";

        std::string found_s = ok ? (std::string(e8_hit?"E8+ECM:":"ECM:")+join_factors_short(found,args.short_digits,args.show_full)) : "-";

        std::cout<<std::setw(2)<<(i+1)<<"  "<<std::setw(11)<<kind_name(s.kind)<<"  "
                 <<std::fixed<<std::setprecision(3)<<std::setw(6)<<mr_ms<<"  "
                 <<std::setw(6)<<e8_ms<<"  "<<std::setw(7)<<e8_stat<<"  "
                 <<std::setw(7)<<ecm_ms<<"  "<<std::setw(8)<<ecm_stat<<"  "
                 <<(ok?"Y ":"N ")<<short_dec(s.n,args.short_digits,args.show_full)<<"  "
                 <<join_factors_short(truth,args.short_digits,args.show_full)<<"  "<<found_s<<"\n" << std::flush;

        auto& A=agg[std::string(kind_name(s.kind))];
        A.total++; A.ok += ok?1:0; A.e8_hit += e8_hit?1:0; if (ecm_stat=="TIMEOUT") A.ecm_to++;
        A.mr.push_back(mr_ms); A.e8.push_back(e8_ms); A.ecm.push_back(ecm_ms); A.tot.push_back(mr_ms+e8_ms+ecm_ms);
    }

    std::cout << "\nSummary (per kind):\n";
    std::cout << "kind         n   ok  e8_hit  ecm_to   MR_med  E8_med  ECM_med  TOT_med  ECM_p90  TOT_p90\n";
    std::cout << "-----------------------------------------------------------------------------------------------------------\n";
    for (auto& kv: agg) {
        auto& A=kv.second;
        std::cout<<std::setw(11)<<kv.first<<"  "<<std::setw(3)<<A.total<<"  "<<std::setw(2)<<A.ok<<"  "
                 <<std::setw(6)<<A.e8_hit<<"  "<<std::setw(6)<<A.ecm_to<<"  "
                 <<std::fixed<<std::setprecision(2)
                 <<std::setw(6)<<quantile(A.mr,0.5)<<"  "<<std::setw(6)<<quantile(A.e8,0.5)<<"  "
                 <<std::setw(7)<<quantile(A.ecm,0.5)<<"  "<<std::setw(7)<<quantile(A.tot,0.5)<<"  "
                 <<std::setw(7)<<quantile(A.ecm,0.9)<<"  "<<std::setw(7)<<quantile(A.tot,0.9)<<"\n";
    }

    std::cout << "\nWrote results CSV: " << args.out_csv << "\n" << std::flush;
    return 0;
}

int main(int argc, char** argv) {
    try {
        Args args = parse(argc, argv);
        if (args.bits <= 128) {
            using Int = mp::uint128_t;
            using Wider = mp::uint256_t;
            return run<Int,Wider>(args);
        } else if (args.bits <= 256) {
            using Int = mp::uint256_t;
            using Wider = mp::uint512_t;
            return run<Int,Wider>(args);
        } else {
            using Int = mp::cpp_int;
            using Wider = mp::cpp_int;
            return run<Int,Wider>(args);
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\nTry: --help\n";
        return 1;
    }
}
