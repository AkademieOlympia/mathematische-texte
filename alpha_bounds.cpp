#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/io.hpp>
#include <boost/numeric/interval/transc.hpp>

using namespace boost::numeric;
using namespace interval_lib;

typedef save_state<rounded_transc_std<double>> rounding_policy;
typedef checking_strict<double> checking_policy;
typedef interval<double, policies<rounding_policy, checking_policy>> strict_interval;

// E8-Gitter: Rang 8, 240 Wurzeln. Komplement der Cartan-Richtungen: 232.
const strict_interval E8_RANK(8.0, 8.0);
const strict_interval E8_ROOTS(240.0, 240.0);
const strict_interval E8_BASE(137.0, 137.0);

const strict_interval GOLDEN_RATIO = (strict_interval(1.0, 1.0) + sqrt(strict_interval(5.0, 5.0))) / 2.0;

// CODATA 2018-Fenster (unteres Ende) bis gerundete 137.036-Kappe.
const strict_interval CODATA(137.03599908, 137.036);

void print_interval(const char* label, const strict_interval& x) {
    std::cout << label << "[" << x.lower() << ", " << x.upper() << "]" << std::endl;
}

strict_interval radial_defect_8d(const strict_interval& base) {
    // 8D-Projektion: D_r = 8 / (137 * φ)
    return E8_RANK / (base * GOLDEN_RATIO);
}

strict_interval transversal_undamped(const strict_interval& base) {
    // Ungedämpfter E8-Torusterm: 232 Wurzelrichtungen * log / 137^2
    // (eine inverse Kernkraft weniger als die physikalische 2-Schleifen-Dämpfung)
    return (E8_ROOTS - E8_RANK) * boost::numeric::log(base) / (base * base);
}

strict_interval damp_inverse_kernel(const strict_interval& dt, const strict_interval& base, int power) {
    strict_interval scale(1.0, 1.0);
    for (int i = 0; i < power; ++i) {
        scale = scale / base;
    }
    return dt * scale;
}

strict_interval hoffbauer_bracket(const strict_interval& d_r, const strict_interval& d_t) {
    return boost::numeric::sqrt(d_r * d_r + d_t * d_t);
}

double midpoint(const strict_interval& x) {
    return 0.5 * (x.lower() + x.upper());
}

void print_stage(const char* title,
                 const strict_interval& dr,
                 const strict_interval& dt,
                 const char* dr_label,
                 const char* dt_label) {
    strict_interval delta = hoffbauer_bracket(dr, dt);
    strict_interval alpha_inv = E8_BASE + delta;
    const double diff_mid = midpoint(alpha_inv) - midpoint(CODATA);

    std::cout << "=================================================" << std::endl;
    std::cout << title << std::endl;
    std::cout << "=================================================" << std::endl;

    std::cout << std::setprecision(13);
    print_interval(dr_label, dr);
    print_interval(dt_label, dt);
    print_interval("EABC Gesamtdefekt        : ", delta);

    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << std::setprecision(12);
    print_interval("alpha^{-1}_{EABC} BOUNDS : ", alpha_inv);
    print_interval("CODATA Referenz          : ", CODATA);

    std::cout << std::setprecision(7);
    std::cout << "Differenz zu CODATA-Mitte: " << diff_mid << std::endl;

    const double adiff = std::fabs(diff_mid);
    if (adiff < 5.0e-4) {
        std::cout << "[DELTA] Naeherung liegt im Bereich von 10^-3 bis 10^-4 zum Soll." << std::endl;
    } else if (adiff < 5.0e-2) {
        std::cout << "[DELTA] Residuum noch O(10^-2); transversale Kerndämpfung erforderlich." << std::endl;
    } else {
        std::cout << "[DELTA] Grobe Abweichung; Modellskala pruefen." << std::endl;
    }
}

int main() {
    const strict_interval dr = radial_defect_8d(E8_BASE);
    const strict_interval dt0 = transversal_undamped(E8_BASE);

    // 1) E8-normiert, Torusterm noch ungedämpft (Diagnose)
    print_stage(
        " EABC E8-Normierte Defektanalyse fuer Alpha^-1   ",
        dr,
        dt0,
        "Radialer Defekt D_r (8D) : ",
        "Transversaler Defekt D_t : "
    );

    // 2) Eine inverse Kernkraft ~ α : D_t / 137  → O(10^{-4})
    const strict_interval dt1 = damp_inverse_kernel(dt0, E8_BASE, 1);
    std::cout << std::endl;
    print_stage(
        " EABC mit 1/137-Kerndaempfung (1-Schleife)       ",
        dr,
        dt1,
        "Radialer Defekt D_r (8D) : ",
        "Transversaler Defekt D_t : "
    );

    // 3) Quadratische Dämpfung ~ α^2 : D_t / 137^2  → O(10^{-6})
    const strict_interval dt2 = damp_inverse_kernel(dt0, E8_BASE, 2);
    std::cout << std::endl;
    print_stage(
        " EABC mit 1/137^2-Kerndaempfung (2-Schleife)     ",
        dr,
        dt2,
        "Radialer Defekt D_r (8D) : ",
        "Transversaler Defekt D_t : "
    );

    std::cout << std::endl;
    std::cout << std::setprecision(13);
    std::cout << "Isolierter Radialrest zu CODATA-eps:" << std::endl;
    const double eps_codata = midpoint(CODATA) - 137.0;
    std::cout << "  8/(137 φ) - ε_CODATA = "
              << (midpoint(dr) - eps_codata) << std::endl;
    std::cout << "  D_t (ungedaempft)    = " << midpoint(dt0) << std::endl;
    std::cout << "  D_t / 137            = " << midpoint(dt1) << std::endl;
    std::cout << "  D_t / 137^2          = " << midpoint(dt2) << std::endl;

    return 0;
}
