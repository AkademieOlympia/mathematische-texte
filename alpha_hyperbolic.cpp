#include <iostream>
#include <iomanip>
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
const strict_interval TRIALITY(3.0, 3.0);
const strict_interval CODATA(137.03599908, 137.03600000);
const double CODATA_MID = 0.5 * (137.03599908 + 137.03600000);

double midpoint(const strict_interval& x) {
    return 0.5 * (x.lower() + x.upper());
}

void print_interval(const char* label, const strict_interval& x) {
    std::cout << label << "[" << x.lower() << ", " << x.upper() << "]" << std::endl;
}

// Hyperbolische Hoffbauer-Klammer: sqrt(D_r^2 - D_t^2), Signatur (+,-).
strict_interval hoffbauer_hyperbolic(const strict_interval& d_r, const strict_interval& d_t) {
    strict_interval diff = (d_r * d_r) - (d_t * d_t);
    const double lo = std::max(0.0, diff.lower());
    const double hi = std::max(lo, diff.upper());
    return boost::numeric::sqrt(strict_interval(lo, hi));
}

strict_interval radial_from_base(const strict_interval& alpha_inv) {
    return DIM_E8 / (alpha_inv * PHI);
}

int main() {
    std::cout << std::setprecision(12);
    std::cout << "==========================================================" << std::endl;
    std::cout << " EABC HYPERBOLISCHE METRIK & TRIALITAETS-ANALYSE          " << std::endl;
    std::cout << "==========================================================" << std::endl;

    const strict_interval dr = radial_from_base(E8_BASE);
    const double eps_codata = CODATA_MID - 137.0;
    const double delta_res = midpoint(dr) - eps_codata;
    const double rel_match = 100.0 * eps_codata / midpoint(dr);

    std::cout << std::setprecision(13);
    print_interval("D_r (radial, 8/(137 φ))  : ", dr);
    std::cout << "ε_CODATA                 : " << eps_codata << std::endl;
    std::cout << "Relative Uebereinstimmung: " << std::setprecision(6) << rel_match << " %" << std::endl;
    std::cout << std::setprecision(13);
    std::cout << "Delta_res = D_r - ε      : " << delta_res << std::endl;

    const double dt_required = std::sqrt(std::max(0.0, midpoint(dr) * midpoint(dr) - eps_codata * eps_codata));
    const double dt_approx = std::sqrt(std::max(0.0, 2.0 * midpoint(dr) * delta_res));
    std::cout << "D_t fuer exakten CODATA  : " << dt_required << std::endl;
    std::cout << "sqrt(2 D_r Delta_res)    : " << dt_approx << std::endl;

    // Variante A: Trialitaets-Daempfung D_t = 1/(3*137)
    const strict_interval dt_triality = 1.0 / (TRIALITY * E8_BASE);
    const strict_interval delta_hyp = hoffbauer_hyperbolic(dr, dt_triality);
    const strict_interval alpha_hyp = E8_BASE + delta_hyp;
    const double res_a = midpoint(alpha_hyp) - CODATA_MID;

    std::cout << "----------------------------------------------------------" << std::endl;
    std::cout << " Variante A: hyperbolisch + 1/(3*137)                     " << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;
    print_interval("D_t (Trialitaet 1/(3*137)): ", dt_triality);
    print_interval("Gesamtdefekt (hyperbol.) : ", delta_hyp);
    std::cout << std::setprecision(12);
    print_interval("alpha^-1 Schranken       : ", alpha_hyp);
    print_interval("CODATA Referenz          : ", CODATA);
    std::cout << std::setprecision(7);
    std::cout << "Residuum zu CODATA-Mitte : " << res_a << std::endl;

    // Variante B: Selbstkonsistenz  α = 137 + 8/(α φ)
    // Quadratisch: α^2 - 137 α - 8/φ = 0
    // α = [137 + sqrt(137^2 + 32/φ)] / 2
    const strict_interval disc = E8_BASE * E8_BASE + 32.0 / PHI;
    const strict_interval alpha_b = (E8_BASE + boost::numeric::sqrt(disc)) / 2.0;
    const strict_interval dr_b = radial_from_base(alpha_b);
    const double res_b = midpoint(alpha_b) - CODATA_MID;

    std::cout << std::setprecision(13);
    std::cout << "----------------------------------------------------------" << std::endl;
    std::cout << " Variante B: Nenner-Selbstkonsistenz α = 137 + 8/(α φ)    " << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;
    print_interval("D_r (8/(α φ))            : ", dr_b);
    std::cout << std::setprecision(12);
    print_interval("alpha^-1 Schranken       : ", alpha_b);
    std::cout << std::setprecision(7);
    std::cout << "Residuum zu CODATA-Mitte : " << res_b << std::endl;

    // Kombination A+B: Fixpunkt α = 137 + sqrt( (8/(α φ))^2 - (1/(3*137))^2 )
    strict_interval alpha_ab = alpha_hyp;
    for (int i = 0; i < 20; ++i) {
        strict_interval dr_ab = radial_from_base(alpha_ab);
        alpha_ab = E8_BASE + hoffbauer_hyperbolic(dr_ab, dt_triality);
    }
    const strict_interval dr_ab = radial_from_base(alpha_ab);
    const strict_interval delta_ab = hoffbauer_hyperbolic(dr_ab, dt_triality);
    const double res_ab = midpoint(alpha_ab) - CODATA_MID;

    std::cout << std::setprecision(13);
    std::cout << "----------------------------------------------------------" << std::endl;
    std::cout << " Kombination A+B: hyperbolisch mit rueckgekoppeltem Nenner" << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;
    print_interval("D_r (8/(α φ))            : ", dr_ab);
    print_interval("Gesamtdefekt (hyperbol.) : ", delta_ab);
    std::cout << std::setprecision(12);
    print_interval("alpha^-1 Schranken       : ", alpha_ab);
    std::cout << std::setprecision(7);
    std::cout << "Residuum zu CODATA-Mitte : " << res_ab << std::endl;

    // Volle Selbstkonsistenz: auch D_t = 1/(3 α)
    strict_interval alpha_full = alpha_ab;
    strict_interval dt_full = 1.0 / (TRIALITY * alpha_full);
    for (int i = 0; i < 20; ++i) {
        strict_interval dr_full = radial_from_base(alpha_full);
        dt_full = 1.0 / (TRIALITY * alpha_full);
        alpha_full = E8_BASE + hoffbauer_hyperbolic(dr_full, dt_full);
    }
    const strict_interval dr_full = radial_from_base(alpha_full);
    const strict_interval delta_full = hoffbauer_hyperbolic(dr_full, dt_full);
    const double res_full = midpoint(alpha_full) - CODATA_MID;

    std::cout << std::setprecision(13);
    std::cout << "----------------------------------------------------------" << std::endl;
    std::cout << " Volle Rueckkopplung: D_t = 1/(3 α)                       " << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;
    print_interval("D_t (1/(3 α))            : ", dt_full);
    print_interval("Gesamtdefekt (hyperbol.) : ", delta_full);
    std::cout << std::setprecision(12);
    print_interval("alpha^-1 Schranken       : ", alpha_full);
    std::cout << std::setprecision(7);
    std::cout << "Residuum zu CODATA-Mitte : " << res_full << std::endl;

    std::cout << "==========================================================" << std::endl;
    return 0;
}
