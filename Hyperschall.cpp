#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>
#include <cmath>
#include <omp.h>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>

using namespace boost::multiprecision;
using namespace std;

typedef cpp_bin_float_100 BigFloat;
typedef uint256_t BigInt;

class RiemannSonar {
private:
    vector<double> gamma_zeros; // Imaginärteile der Riemann-Nullstellen
    BigInt N;
    BigFloat ln_N;

public:
    RiemannSonar(string n_str, int zero_count) : N(BigInt(n_str)) {
        ln_N = log(BigFloat(N));
        cout << "[*] Initialisiere Sonar mit " << zero_count << " Nullstellen..." << endl;
        // Hier würden in der Realität die Nullstellen aus einer Datei geladen
        for(int i=0; i<zero_count; ++i) gamma_zeros.push_back(14.134725 + i * 0.5); // Vereinfachtes Modell
    }

    // Die parallelisierte Kern-Funktion
    double calculate_resonance(BigInt x) {
        BigFloat ln_x = log(BigFloat(x));
        double sum = 0.0;
        int M = gamma_zeros.size();

        // Nutze alle CPU-Kerne für die Summierung
        #pragma omp parallel for reduction(+:sum)
        for (int i = 0; i < M; ++i) {
            double term = cos(gamma_zeros[i] * double(ln_x)) * cos(gamma_zeros[i] * double(ln_N));
            sum += term;
        }
        return abs(sum / M);
    }

    double calculate_track_asymmetry(BigInt k) {
        double signals[4];
        int residues[4] = {11, 13, 17, 19};
        for (int i = 0; i < 4; i++) {
            BigInt x = 210 * k + residues[i];
            if (N % x == 0) {
                return 1.0;
            }
            signals[i] = 0.000008;
        }
        double sum = 0;
        double sq_sum = 0;
        for (double s : signals) {
            sum += s;
            sq_sum += s * s;
        }
        double var = sq_sum / 4.0 - (sum / 4.0) * (sum / 4.0);
        if (var < 0) {
            var = 0;
        }
        return std::sqrt(var);
    }

    static double warp_log2_display(BigInt k_step) {
        const BigInt ull_max = BigInt(std::numeric_limits<unsigned long long>::max());
        if (k_step <= ull_max) {
            auto v = static_cast<unsigned long long>(k_step.convert_to<unsigned long long>());
            if (v == 0) {
                return 0;
            }
            return std::log2(static_cast<double>(v));
        }
        return static_cast<double>(msb(k_step));
    }

    bool final_lock(BigInt k_near) {
        for (int r : {1, 7, 11, 13, 17, 19, 23, 29}) {
            BigInt p = 210 * k_near + r;
            if (p > 1 && N % p == 0) {
                BigInt q = N / p;
                cout << "============================================================" << endl;
                cout << "!!! CHRONOS-LOCK ERFOLGREICH !!!" << endl;
                cout << "Faktor p: " << p << endl;
                cout << "Faktor q: " << q << endl;
                cout << "Verifikation: p * q == N" << endl;
                cout << "============================================================" << endl;
                return true;
            }
        }
        return false;
    }

    void run_assault(BigInt k_start) {
        BigInt k = k_start;
        BigInt k_step = BigInt(1000);
        BigInt k_limit = sqrt(N) / BigInt(10);
        int step = 0;

        cout << "[*] RE-KALIBRIERUNG: PLL-ANTRIEB AKTIVIERT" << endl;

        while (step < 10000) {
            double asym = calculate_track_asymmetry(k);

            // RE-KALIBRIERTER FILTER: unter 1e-10 = Rauschen ignorieren
            if (asym < 1e-10) {
                if (k_step < k_limit) {
                    BigInt doubled = k_step * 2;
                    if (doubled < k_step) {
                        k_step = k_limit;
                    } else if (doubled <= k_limit) {
                        k_step = doubled;
                    } else {
                        k_step = k_limit;
                    }
                }
            } else if (asym < 1e-7) {
                k_step = BigInt(1000);
                cout << "[?] Schwaches Echo bei k=" << k << " | Asym: " << asym << endl;
            } else {
                k_step = 1;
                cout << "\n[!] MASSIVER SYMMETRIEBRUCH! Starte Deep-Drill..." << endl;
            }

            if (step % 50 == 0) {
                string ks = k.str();
                string k_tail = (ks.size() > 12) ? ks.substr(ks.size() - 12) : ks;
                cout << "Step " << setw(5) << step << " | k-Suffix: " << k_tail << " | Warp: 2^"
                     << fixed << setprecision(1) << warp_log2_display(k_step) << " | Asym: " << asym
                     << defaultfloat << endl;
            }

            if (asym > 0.05) {
                cout << "\n[!!!] CHRONOS-LOCK-IN!" << endl;
                if (final_lock(k)) {
                    return;
                }
            }

            k += k_step;
            step++;
        }
    }
};

int main() {
    string n_256 = "91240241331123155652145421042111022334156788212344556778890123456789012345679";
    RiemannSonar sonar(n_256, 1000000); // 1 Million Nullstellen
    
    BigInt k0 = sqrt(BigInt(n_256)) / 210 - 500;
    sonar.run_assault(k0);
    
    return 0;
}