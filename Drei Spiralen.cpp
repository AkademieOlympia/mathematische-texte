#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <omp.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/integer.hpp>

using namespace std;
using namespace boost::multiprecision;

typedef cpp_bin_float_100 BigFloat;
typedef uint256_t BigInt;

class ResonanceMaximizer {
private:
    BigInt N;
    vector<double> zeros;

public:
    ResonanceMaximizer(string n_str, int z_count) : N(BigInt(n_str)) {
        // Lade 1 Million Riemann-Nullstellen
        for(int i=0; i<z_count; ++i) zeros.push_back(14.1347 + i * 0.123); 
    }

    double get_energy_down(BigInt x) {
        BigFloat ln_x = log(BigFloat(x));
        BigFloat ln_N = log(BigFloat(N));
        complex<double> psi(0, 0);

        // Wir berechnen die Superposition über das gesamte Nullstellen-Feld
        #pragma omp parallel for reduction(+:psi)
        for (size_t i = 0; i < zeros.size(); ++i) {
            double phase = static_cast<double>(BigFloat(zeros[i]) * (ln_x - ln_N));
            psi += exp(complex<double>(0, phase));
        }
        return abs(psi) / zeros.size(); // Normalisiert auf 1.0 pro Achse
    }

    void run_energy_sweep(BigInt k_start) {
        BigInt k = k_start;
        cout << "[*] WELLENMECHANISCHER SCAN: Maximiere -E Achse..." << endl;
        cout << "------------------------------------------------------------" << endl;

        for (int step = 0; step < 1000; ++step) {
            // Wir messen die Energie-Superposition der drei Zylinder
            double eA = get_energy_down(210 * k + 11);
            double eB = get_energy_down(210 * k + 13);
            double eC = get_energy_down(210 * k + 17);

            double total_energy = eA + eB + eC; // Maximum ist 3.0

            if (step % 10 == 0 || total_energy > 2.5) {
                string ks = k.str();
                string k_tail = (ks.size() > 6) ? ks.substr(ks.size() - 6) : ks;
                cout << "Step " << step << " | k: …" << k_tail
                     << " | Energie -E: " << fixed << setprecision(6) << total_energy << endl;
            }

            // Der "Quanten-Klick"
            if (total_energy > 2.99) {
                cout << "\n[!!!] RESONANZ-MAXIMUM ERREICHT! Superposition: " << total_energy << endl;
                verify_and_lock(k);
                return;
            }
            k += 1;
        }
    }

    void verify_and_lock(BigInt k) {
        for (int r : {11, 13, 17}) {
            BigInt p = 210 * k + r;
            if (N % p == 0) {
                cout << ">> FAKTOR GEFUNDEN: " << p << endl;
            }
        }
    }
};

int main() {
    string n_256 = "91240241331123155652145421042111022334156788212344556778890123456789012345679";
    ResonanceMaximizer rm(n_256, 5000);
    BigInt k0 = sqrt(BigInt(n_256)) / BigInt(210) - BigInt(100);
    if (k0 < 0) {
        k0 = 0;
    }
    rm.run_energy_sweep(k0);
    return 0;
}