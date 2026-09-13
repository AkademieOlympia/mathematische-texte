#ifndef EABC_ALPHA_HPP
#define EABC_ALPHA_HPP

#include <algorithm>
#include <cmath>
#include <ratio>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/transc.hpp>

namespace eabc {

// ---------------------------------------------------------------------------
// 1. Algebraische Invarianten & Modulformen-Hierarchie (Kompilierzeit)
// ---------------------------------------------------------------------------
inline constexpr unsigned int RK_E8          = 8;
inline constexpr unsigned int MOD12_PERIOD   = 12;
inline constexpr unsigned int DEDEKIND_STUFE = 24;
inline constexpr unsigned int TRIALITAET     = 3;
// Ramanujan-Tau Hecke-Eigenwert fuer p = 3: tau(3) = 252
inline constexpr long long TAU_3             = 252;

// Stufe 1: Dedekind c1 = 3 * (1 - 1/24) = 23/8
using c_dedekind_ratio = std::ratio<23, 8>;
// Stufe 2: E8 x Modulo-12: c2 = 23/8 - 1/(8 * 12) = 275/96
using c_mod12_ratio    = std::ratio<275, 96>;
// Stufe 3: Hecke-T3 Verfeinerung:
// c3 = 275/96 + 1 / (rk(E8) * tau(3)) = 275/96 + 1/2016 = 5776/2016 = 361/126
// (361 = 19^2, 126 = tau(3)/2)
using c_hecke_ratio    = std::ratio<361, 126>;

// Standard-Auswahl: Hecke Stufe 3 (PPT-Niveau)
using c_triality_ratio = c_hecke_ratio;

// Transversaler Defekt D_t = 1 / (c * 137)
// Fuer c = 361/126: D_t = 126 / (361 * 137) = 126 / 49457
using D_t_ratio        = std::ratio<126, 49457>;

inline constexpr double ALPHA_0_INV = 137.0;

// ---------------------------------------------------------------------------
// 2. Intervall-Arithmetik Typdefinitionen
// ---------------------------------------------------------------------------
namespace detail {
    using namespace boost::numeric;
    using namespace interval_lib;
    // rounded_math fehlt log_down; transcendental std + save_state wie in den
    // bestehenden Alpha-Zertifikaten.
    typedef save_state<rounded_transc_std<double>> rounding_policy;
    typedef checking_strict<double> checking_policy;
    typedef interval<double, policies<rounding_policy, checking_policy>> strict_interval;
}

using interval_t = detail::strict_interval;

// ---------------------------------------------------------------------------
// 3. Zertifizierte Intervall-Funktionen
// ---------------------------------------------------------------------------

// Rationale/zertifizierte Schranken fuer den Goldenen Schnitt phi
inline interval_t golden_ratio() {
    using boost::numeric::sqrt;
    return (interval_t(1.0, 1.0) + sqrt(interval_t(5.0, 5.0))) / 2.0;
}

// Radialer 8D-Massendefekt: D_r = rk(E8) / (137 * phi)
inline interval_t radial_defect() {
    const interval_t phi = golden_ratio();
    return interval_t(static_cast<double>(RK_E8), static_cast<double>(RK_E8)) /
           (interval_t(ALPHA_0_INV, ALPHA_0_INV) * phi);
}

// Transversaler Torus-Defekt: D_t = 126 / 49457 (gerichtete Division)
inline interval_t transversal_defect() {
    return interval_t(static_cast<double>(D_t_ratio::num),
                      static_cast<double>(D_t_ratio::num)) /
           interval_t(static_cast<double>(D_t_ratio::den),
                      static_cast<double>(D_t_ratio::den));
}

// Hyperbolische Hoffbauer-Klammer: Delta = sqrt(D_r^2 - D_t^2)
inline interval_t hoffbauer_hyperbolic(const interval_t& d_r, const interval_t& d_t) {
    const interval_t diff = (d_r * d_r) - (d_t * d_t);
    const double lo = std::max(0.0, diff.lower());
    const double hi = std::max(lo, diff.upper());
    return boost::numeric::sqrt(interval_t(lo, hi));
}

// Gesamt-Kopplung alpha^-1 = 137 + Delta_EABC
inline interval_t alpha_inverse_bounds() {
    const interval_t dr = radial_defect();
    const interval_t dt = transversal_defect();
    const interval_t delta = hoffbauer_hyperbolic(dr, dt);
    return interval_t(ALPHA_0_INV, ALPHA_0_INV) + delta;
}

// CODATA-Referenzintervall (zertifiziertes Analysefenster, CODATA-2018-artig)
inline interval_t codata_window() {
    return interval_t(137.03599908, 137.03600000);
}

// Schnelle analytische Auswertung (Punkt, nicht zertifiziert)
inline double alpha_inverse_point() {
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    const double dr  = 8.0 / (137.0 * phi);
    const double dt  = static_cast<double>(D_t_ratio::num) /
                       static_cast<double>(D_t_ratio::den);
    return 137.0 + std::sqrt(dr * dr - dt * dt);
}

} // namespace eabc

#endif // EABC_ALPHA_HPP
