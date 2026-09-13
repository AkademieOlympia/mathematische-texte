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

void test_spin8_mod(const std::string& name, const strict_interval& c_tri, bool feedback_alpha) {
    strict_interval alpha_guess = E8_BASE;
    strict_interval dt(0.0, 0.0);
    for (int i = 0; i < 3; ++i) {
        strict_interval denom = feedback_alpha ? alpha_guess : E8_BASE;
        strict_interval dr = DIM_E8 / (denom * PHI);
        dt = 1.0 / (c_tri * denom);
        strict_interval delta = hoffbauer_hyperbolic(dr, dt);
        alpha_guess = E8_BASE + delta;
    }
    const double mid = midpoint(alpha_guess);
    const double diff = mid - CODATA_MID;
    const bool hit = (alpha_guess.upper() >= CODATA.lower() && alpha_guess.lower() <= CODATA.upper());

    std::cout << std::left << std::setw(34) << name
              << std::right << std::fixed << std::setprecision(12)
              << " : [" << alpha_guess.lower() << ", " << alpha_guess.upper() << "]"
              << " | Dt=" << std::scientific << std::setprecision(6) << midpoint(dt)
              << " | Res=" << std::showpos << diff << std::noshowpos;
    if (hit) {
        std::cout << "  --> [CODATA HIT]";
    }
    std::cout << std::endl;
}

int main() {
    const strict_interval dr0 = DIM_E8 / (E8_BASE * PHI);
    const double eps = CODATA_MID - 137.0;
    const double dt_soll = std::sqrt(std::max(0.0, midpoint(dr0) * midpoint(dr0) - eps * eps));
    const double c_soll = 1.0 / (dt_soll * 137.0);

    std::cout << "==========================================================================================" << std::endl;
    std::cout << " SPIN(8) TRIALITAETS-FEINKORREKTUR & INTERVALL-EINSCHLUSS                                " << std::endl;
    std::cout << "==========================================================================================" << std::endl;
    std::cout << std::setprecision(12);
    std::cout << "D_t^soll = sqrt(D_r^2 - eps^2) : " << dt_soll << std::endl;
    std::cout << "c_tri^soll = 1/(D_t^soll * 137): " << c_soll << std::endl;
    std::cout << "CODATA Fenster                 : [" << CODATA.lower() << ", " << CODATA.upper() << "]" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;

    test_spin8_mod("Trialitaet c = 3", strict_interval(3.0, 3.0), false);
    test_spin8_mod("Trialitaet c = 3 (Feedback)", strict_interval(3.0, 3.0), true);

    const strict_interval c_dedekind = strict_interval(23.0, 23.0) / strict_interval(8.0, 8.0);
    test_spin8_mod("Dedekind c = 23/8 = 2.875", c_dedekind, false);
    test_spin8_mod("Dedekind c = 23/8 (Feedback)", c_dedekind, true);

    const strict_interval c_coxeter =
        strict_interval(3.0, 3.0) * (strict_interval(29.0, 29.0) / strict_interval(30.0, 30.0));
    test_spin8_mod("Coxeter c = 3*(29/30) = 2.90", c_coxeter, false);
    test_spin8_mod("Coxeter c = 3*(29/30) (Feedback)", c_coxeter, true);

    const strict_interval c_phi = strict_interval(3.0, 3.0) - (1.0 / (2.0 * PHI * PHI));
    test_spin8_mod("Modul-Fluss c = 3 - 1/(2phi^2)", c_phi, false);
    test_spin8_mod("Modul-Fluss c (Feedback)", c_phi, true);

    const strict_interval c_exact(c_soll, c_soll);
    test_spin8_mod("Exakt c_soll (Referenz)", c_exact, false);
    test_spin8_mod("Exakt c_soll (Feedback)", c_exact, true);

    std::cout << "==========================================================================================" << std::endl;
    return 0;
}
