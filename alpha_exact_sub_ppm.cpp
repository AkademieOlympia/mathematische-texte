#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/io.hpp>
#include <boost/numeric/interval/transc.hpp>

using namespace boost::numeric;
using namespace interval_lib;

typedef save_state<rounded_transc_std<double>> rounding_policy;
typedef checking_strict<double> checking_policy;
typedef interval<double, policies<rounding_policy, checking_policy>> strict_interval;

const strict_interval PHI =
    (strict_interval(1.0, 1.0) + boost::numeric::sqrt(strict_interval(5.0, 5.0))) / 2.0;
const strict_interval E8_BASE(137.0, 137.0);
const strict_interval DIM_E8(8.0, 8.0);
const strict_interval CODATA(137.03599908, 137.03600000);
const double CODATA_MID = 0.5 * (137.03599908 + 137.03600000);

double midpoint(const strict_interval& x) {
    return 0.5 * (x.lower() + x.upper());
}

strict_interval hoffbauer_hyperbolic(const strict_interval& d_r, const strict_interval& d_t) {
    strict_interval diff = (d_r * d_r) - (d_t * d_t);
    const double lo = std::max(0.0, diff.lower());
    const double hi = std::max(lo, diff.upper());
    return boost::numeric::sqrt(strict_interval(lo, hi));
}

void check_candidate(const std::string& label, const strict_interval& c_cand) {
    const strict_interval dr = DIM_E8 / (E8_BASE * PHI);
    const strict_interval dt = 1.0 / (c_cand * E8_BASE);
    const strict_interval delta = hoffbauer_hyperbolic(dr, dt);
    const strict_interval alpha_inv = E8_BASE + delta;
    const double mid = midpoint(alpha_inv);
    const double diff = mid - CODATA_MID;
    const bool hit = (alpha_inv.upper() >= CODATA.lower() && alpha_inv.lower() <= CODATA.upper());

    std::cout << std::left << std::setw(38) << label
              << " c=" << std::fixed << std::setprecision(10) << midpoint(c_cand)
              << " : [" << alpha_inv.lower() << ", " << alpha_inv.upper() << "]"
              << " | Res: " << std::showpos << std::scientific << std::setprecision(3) << diff
              << std::noshowpos;
    if (hit) {
        std::cout << "  --> [CODATA HIT EXAKT!]";
    }
    std::cout << std::endl;
}

int main() {
    std::cout << "==========================================================================================" << std::endl;
    std::cout << " EABC SUB-PPM FEINKORREKTUR: MODULO-12 / RAMANUJAN ORDNUNGEN                             " << std::endl;
    std::cout << "==========================================================================================" << std::endl;

    const strict_interval eight(8.0, 8.0);
    const strict_interval twenty_three(23.0, 23.0);
    const strict_interval c_dedekind = twenty_three / eight;
    check_candidate("Dedekind: 23/8", c_dedekind);

    // Ramanujan-E8 Schnittmenge: 23/8 - 1/(8*12) = 275/96
    const strict_interval c_96 = c_dedekind - (1.0 / strict_interval(96.0, 96.0));
    check_candidate("E8 x Mod12: 23/8 - 1/96", c_96);

    // E8 x 24: 23/8 - 1/(8*24) = 23/8 - 1/192
    const strict_interval c_192 = c_dedekind - (1.0 / strict_interval(192.0, 192.0));
    check_candidate("E8 x Ramanujan: 23/8 - 1/192", c_192);

    // Modulo-12 quadratisch: 23/8 - 1/(12^2 - 12) = 23/8 - 1/132
    const strict_interval c_132 = c_dedekind - (1.0 / strict_interval(132.0, 132.0));
    check_candidate("Mod12-Spektrum: 23/8 - 1/132", c_132);

    // Coxeter-Nachkorrektur auf Dedekind: 23/8 - 1/(30 φ^4)
    const strict_interval c_cox_phi =
        c_dedekind - (1.0 / (strict_interval(30.0, 30.0) * PHI * PHI * PHI * PHI));
    check_candidate("Coxeter/phi^4: 23/8 - 1/(30 phi^4)", c_cox_phi);

    // Zweite Ordnung: 3*(1 - 1/24 - 1/(24*137*phi))
    const strict_interval c_2nd =
        strict_interval(3.0, 3.0) *
        (strict_interval(1.0, 1.0) - (1.0 / strict_interval(24.0, 24.0)) -
         (1.0 / (strict_interval(24.0, 24.0) * E8_BASE * PHI)));
    check_candidate("2. Ordnung: 3(1-1/24-1/(24*137*phi))", c_2nd);

    const strict_interval c_exact(2.865078475, 2.865078475);
    check_candidate("c_soll (Referenz)", c_exact);

    std::cout << "==========================================================================================" << std::endl;
    return 0;
}
