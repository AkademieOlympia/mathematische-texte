#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#include "eabc_alpha.hpp"

void print_stage(const std::string& name, double c_val, double dt_val) {
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    const double dr  = 8.0 / (137.0 * phi);
    const double delta = std::sqrt(dr * dr - dt_val * dt_val);
    const double alpha_inv = 137.0 + delta;
    const double res = alpha_inv - 137.03599954;
    std::cout << std::left << std::setw(28) << name
              << " | c = " << std::setw(10) << c_val
              << " | alpha^-1 = " << std::fixed << std::setprecision(11) << alpha_inv
              << " | Res: " << std::showpos << std::scientific << std::setprecision(2) << res
              << std::noshowpos
              << std::fixed << std::endl;
}

int main() {
    std::cout << "========================================================================================\n";
    std::cout << " EABC HIERARCHIE-VERGLEICH: DEDEKIND -> MOD12 -> HECKE T3\n";
    std::cout << "========================================================================================\n";

    print_stage("1. Dedekind (23/8)",
                23.0 / 8.0,
                1.0 / ((23.0 / 8.0) * 137.0));
    print_stage("2. Modulo-12 (275/96)",
                275.0 / 96.0,
                96.0 / (275.0 * 137.0));
    print_stage("3. Hecke T3 (361/126)",
                361.0 / 126.0,
                126.0 / 49457.0);

    std::cout << "----------------------------------------------------------------------------------------\n";
    const auto bounds = eabc::alpha_inverse_bounds();
    const auto codata = eabc::codata_window();
    std::cout << std::setprecision(11);
    std::cout << "Hecke-Intervallgrenzen       : [" << bounds.lower() << ", " << bounds.upper() << "]\n";
    std::cout << "CODATA-Analysefenster        : [" << codata.lower() << ", " << codata.upper() << "]\n";
    std::cout << "Punktauswertung (Hecke)      : " << eabc::alpha_inverse_point() << "\n";
    if (bounds.lower() >= codata.lower() && bounds.upper() <= codata.upper()) {
        std::cout << "Status: Strikter Einschluss gewaehrleistet.\n";
    } else {
        std::cout << "Status: Kein strikter Einschluss (Intervall ragt aus dem Fenster).\n";
    }
    std::cout << "========================================================================================\n";
    return 0;
}
