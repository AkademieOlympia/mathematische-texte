#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <complex>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <omp.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/integer.hpp>

#include "ZetaNavigation.h"

using namespace std;
using namespace boost::multiprecision;

typedef cpp_bin_float_100 BigFloat;
typedef uint256_t BigInt;

static string bigint_tail(const BigInt& v, size_t max_chars = 14) {
    string s = static_cast<string>(v.str());
    return s.size() > max_chars ? string("…") + s.substr(s.size() - max_chars) : s;
}

static BigInt big_pow10(int exponent) {
    BigInt r(1);
    for (int e = 0; e < exponent; ++e)
        r *= 10;
    return r;
}

// |drift|·1e40 als Schrittweite (gekappt für uint256_t / double-Rand)
static BigInt drift_scaled_step(double drift) {
    double s = std::fabs(drift * 1e40);
    if (!std::isfinite(s) || s < 1.0)
        return BigInt("1000000");
    constexpr double cap = 1e35;
    if (s > cap)
        s = cap;
    if (s <= static_cast<double>(std::numeric_limits<unsigned long long>::max()))
        return BigInt(static_cast<unsigned long long>(s));
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(0) << s;
    return BigInt(oss.str());
}

class SingularityDrill {
private:
    BigInt N;
    BigFloat ln_N;
    vector<double> real_zeros; // Im(ρ) auf Re(ρ)=1/2, aus zeros6.npy (Projektordner)
    double inv_three_nzeros_{0}; // 1 / (3 * n) für Feldstärke-Normierung

    void load_real_zeros_from_project() {
        namespace fs = std::filesystem;
        const fs::path path = fs::path(__FILE__).parent_path() / "zeros6.npy";
        real_zeros = load_riemann_gammas_from_npy(path.string());
        if (real_zeros.empty())
            throw runtime_error(
                "zeros6.npy nicht lesbar oder kein 1D float64 (<f8) — erwartet in: " + path.string());

        inv_three_nzeros_ = 1.0 / (3.0 * static_cast<double>(real_zeros.size()));
        clog << "[*] " << real_zeros.size() << " Nullstellen aus " << path.string() << " geladen." << endl;
    }

    // Einzelgleis r ∈ {11,13,17,19,…}: |Ψ_r|/n (k ist der arithmetische Index, x = 210k+r)
    double calculate_track_intensity(BigInt k, int r) {
        const size_t nz = real_zeros.size();
        if (nz == 0)
            return 0;

        const BigFloat ln_x = log(BigFloat(210 * k + r));
        const double phase_scale = static_cast<double>(ln_x - ln_N);
        const double* const gz = real_zeros.data();

        double sum_re = 0;
        double sum_im = 0;
        #pragma omp parallel for reduction(+ : sum_re, sum_im) schedule(static) if (omp_get_level() == 0)
        for (size_t j = 0; j < nz; ++j) {
            const double pd = gz[j] * phase_scale;
            sum_re += std::cos(pd);
            sum_im += std::sin(pd);
        }
        return std::hypot(sum_re, sum_im) / static_cast<double>(nz);
    }

public:
    SingularityDrill(string n_str) : N(BigInt(n_str)) {
        ln_N = log(BigFloat(N));
        load_real_zeros_from_project();
    }

    // Die -E Energiebilanz (Superposition der 3 Zylinder)
    double calculate_field_intensity(BigInt k) {
        const size_t nz = real_zeros.size();
        if (nz == 0)
            return 0;

        static constexpr int tracks[3] = {11, 13, 17};
        const double* const gz = real_zeros.data();

        complex<double> total_psi(0, 0);

        for (int ti = 0; ti < 3; ++ti) {
            const int r = tracks[ti];
            const BigFloat ln_x = log(BigFloat(210 * k + r));
            const double phase_scale = static_cast<double>(ln_x - ln_N);

            double sum_re = 0;
            double sum_im = 0;
            // if(omp_get_level()==0): kein verschachteltes Team bei Aufrufen aus parallel for
            #pragma omp parallel for reduction(+ : sum_re, sum_im) schedule(static) if (omp_get_level() == 0)
            for (size_t j = 0; j < nz; ++j) {
                const double phase_diff = gz[j] * phase_scale;
                sum_re += std::cos(phase_diff);
                sum_im += std::sin(phase_diff);
            }
            total_psi += complex<double>(sum_re, sum_im);
        }

        return std::abs(total_psi) * inv_three_nzeros_;
    }

    double calculate_amplified_intensity(BigInt k) {
        double raw_intensity = calculate_field_intensity(k);

        double noise_floor = 0.0000062336;
        double signal = raw_intensity - noise_floor;

        if (signal < 0)
            return 0;

        return pow(signal * 1e6, 2.0);
    }

    // Phasendrift Bahn r=11: Mittelwert sin(γ (ln x − ln N)) — Richtung der Phase
    double calculate_phase_drift(BigInt k) {
        const size_t nz = real_zeros.size();
        if (nz == 0)
            return 0;

        const BigFloat ln_x = log(BigFloat(210 * k + 11));
        const BigFloat ln_N_local = log(BigFloat(N));

        double drift = 0;
        #pragma omp parallel for reduction(+ : drift) schedule(static) if (omp_get_level() == 0)
        for (size_t i = 0; i < nz; ++i) {
            const double theta = static_cast<double>(real_zeros[i] * (ln_x - ln_N_local));
            drift += std::sin(theta);
        }
        return drift / static_cast<double>(nz);
    }

    void execute_drill(BigInt k_start) {
        BigInt low = k_start;
        BigInt high = k_start + BigInt("100000000000000000000"); // 10^20 Suchbereich
        int bisection_step = 0;

        cout << "\n[!] LINEARER STILLSTAND DETEKTIERT. SCHALTE UM AUF BINÄRE ZANGE..." << endl;

        const double t0_bisect = omp_get_wtime();

        while (bisection_step < 256) {
            BigInt mid = (low + high) / 2;

            double intensity_low = calculate_field_intensity(low);
            double intensity_mid = calculate_field_intensity(mid);

            const bool go_forward = intensity_mid > intensity_low;
            if (go_forward) {
                low = mid;
            } else {
                high = mid;
            }

            const double elapsed = omp_get_wtime() - t0_bisect;
            string span_str = static_cast<string>((high - low).str());
            const size_t span_digits = span_str.length();
            const double d_mid_low = intensity_mid - intensity_low;

            cout << fixed << setprecision(6) << "Bisection " << setw(3) << (bisection_step + 1) << "/256 | t="
                 << setprecision(2) << elapsed << "s | " << fixed << setprecision(10) << "I_mid=" << intensity_mid
                 << " I_low=" << intensity_low << " | Δ(mid-low)=" << scientific << setprecision(4) << d_mid_low
                 << " | " << fixed << (go_forward ? "→ vorwärts (low←mid)" : "→ rückwärts (high←mid)")
                 << " | Intervall " << span_digits << " Dez-Ziffern | mid " << bigint_tail(mid, 12) << "\r"
                 << flush;

            if ((high - low) < BigInt(2100)) {
                cout << "\n\n[!!!] ZIEL ANVISIERT. Starte materiellen Scan..." << endl;
                BigInt scan_lo = (low > BigInt(10)) ? low - BigInt(10) : BigInt(0);
                for (BigInt test_k = scan_lo; test_k <= high + BigInt(10); ++test_k) {
                    if (final_lock(test_k))
                        return;
                }
            }
            bisection_step++;
        }
        cout << "\n[*] Bisection: 256 Schritte ohne Abbruch — Intervallende erreicht." << endl;
    }

    void execute_final_drill(BigInt center_k) {
        cout << "\n[*] AKTIVIERE MULTI-POINT-DRILL (Scan-Array)..." << endl;
        cout << "[*] center_k " << bigint_tail(center_k, 16) << " | ±5000 | " << omp_get_max_threads()
             << " Threads" << endl;

        constexpr long k_final_lo = -5000;
        constexpr long k_final_hi = 5000;
        const long final_total = k_final_hi - k_final_lo + 1;
        const double t0_final = omp_get_wtime();
        std::atomic<long> final_done{0};

        #pragma omp parallel for schedule(static)
        for (long i = k_final_lo; i <= k_final_hi; ++i) {
            const long n = ++final_done;
            if (n == 1L || n % 1000L == 0L || n == final_total) {
                #pragma omp critical(final_drill_prog)
                {
                    const double te = omp_get_wtime() - t0_final;
                    const double pct = 100.0 * static_cast<double>(n) / static_cast<double>(final_total);
                    const double rate = static_cast<double>(n) / std::max(te, 1e-9);
                    cout << "\r[Final-Drill] " << fixed << setprecision(1) << pct << "% | " << n << "/"
                         << final_total << " | i=" << i << " | " << setprecision(1) << rate << " Schritte/s | t="
                         << setprecision(2) << te << "s   " << flush;
                }
            }

            BigInt test_k;
            if (i >= 0) {
                test_k = center_k + BigInt(static_cast<unsigned long>(i));
            } else {
                BigInt mag = BigInt(static_cast<unsigned long>(-i));
                if (center_k < mag)
                    continue;
                test_k = center_k - mag;
            }

            double E = calculate_field_intensity(test_k);

            if (E > 0.0001) {
                #pragma omp critical
                {
                    cout << "\n[!!!] ECHTE RESONANZ BEI k = " << test_k << " | E = " << E << endl;
                    check_factors(test_k);
                }
            }
        }
        cout << "\n[*] Final-Drill Scan beendet." << endl;
    }

    void execute_jitter_drill(BigInt& k_target) {
        const double noise_floor = 0.0000062336;

        cout << "\n[*] GEISTER-LOCK AUFGEHOBEN. Starte Jitter-Abtastung..." << endl;
        cout << "[*] k_ref " << bigint_tail(k_target, 16) << " | 10 Orbits × 11 Samples" << endl;

        const double t0_jit = omp_get_wtime();

        for (int orbit = 1; orbit <= 10; ++orbit) {
            BigInt radius = BigInt("1000000000000") * orbit;

            for (int i = -5; i <= 5; ++i) {
                BigInt mag = radius * BigInt(static_cast<unsigned long>(std::abs(i)));
                BigInt test_k;
                if (i >= 0) {
                    test_k = k_target + mag;
                } else {
                    if (k_target < mag)
                        continue;
                    test_k = k_target - mag;
                }

                double E = calculate_field_intensity(test_k);
                double signal = std::fabs(E - noise_floor);

                const double tj = omp_get_wtime() - t0_jit;
                const int sample = (orbit - 1) * 11 + (i + 5) + 1;
                const double r_d = radius.convert_to<double>();
                cout << fixed << "Jitter " << setw(2) << sample << "/110 | Orbit " << orbit << "/10 | i="
                     << setw(2) << i << " | r=" << scientific << setprecision(3) << r_d
                     << " | k " << bigint_tail(test_k, 10) << " | E=" << fixed << setprecision(8) << E
                     << " | |ΔNoise|=" << scientific << setprecision(4) << signal << " | t=" << fixed
                     << setprecision(2) << tj << "s\r" << flush;

                if (signal > 1e-9) {
                    cout << "\n\n[!!!] GRADIENT ERFASST! Neuer Vektor bei k = " << test_k << endl;
                    k_target = test_k;
                    return;
                }
            }
        }
        cout << "\n[*] Jitter-Abtastung beendet (kein Schwellen-Treffer)." << endl;
    }

    void execute_spectral_sweep(BigInt center_k) {
        cout << "\n[*] BISECTION ABGEBROCHEN. STARTE SPEKTRAL-FÄCHER-SCAN..." << endl;
        const long rad = 1000000;
        const long step = 100;
        const long spectral_total = (2 * rad) / step + 1;
        cout << "[*] center_k " << bigint_tail(center_k, 16) << " | ±" << rad << " | Schritt " << step << " | "
             << spectral_total << " Gitterpunkte | " << omp_get_max_threads() << " Threads | dynamic" << endl;

        const double t0_spec = omp_get_wtime();
        std::atomic<long> spec_done{0};
        std::atomic<int> spec_hits{0};

        #pragma omp parallel for schedule(dynamic)
        for (long i = -rad; i <= rad; i += step) {
            const long n = ++spec_done;

            BigInt test_k;
            bool evaluated = false;
            if (i >= 0) {
                test_k = center_k + BigInt(static_cast<unsigned long>(i));
                evaluated = true;
            } else {
                BigInt mag = BigInt(static_cast<unsigned long>(-i));
                if (center_k >= mag) {
                    test_k = center_k - mag;
                    evaluated = true;
                }
            }

            double energy = 0;
            if (evaluated)
                energy = calculate_field_intensity(test_k);

            if (n == 1L || n % 400L == 0L || n == spectral_total) {
                #pragma omp critical(spectral_prog)
                {
                    const double ts = omp_get_wtime() - t0_spec;
                    const double pct = 100.0 * static_cast<double>(n) / static_cast<double>(spectral_total);
                    const double rate = static_cast<double>(n) / std::max(ts, 1e-9);
                    const int hits = spec_hits.load(std::memory_order_relaxed);
                    cout << "\r[Spektral] " << fixed << setprecision(1) << pct << "% | " << n << "/"
                         << spectral_total << " | i=" << i;
                    if (evaluated)
                        cout << " | k " << bigint_tail(test_k, 12) << " | E=" << scientific << setprecision(4)
                             << energy;
                    else
                        cout << " | (übersprungen: |i|>center_k)";
                    cout << " | Treffer>E_th: " << hits << " | " << fixed << setprecision(1) << rate
                         << " Pkt/s | t=" << setprecision(1) << ts << "s | thr=" << omp_get_max_threads()
                         << "   " << flush;
                }
            }

            if (evaluated && energy > 0.005) {
                spec_hits.fetch_add(1, std::memory_order_relaxed);
                #pragma omp critical
                {
                    cout << "\n[!!!] SIGNATURE MATCH BEI k = " << test_k << " | E = " << energy << endl;
                    BigInt micro_lo =
                        (test_k > BigInt(100)) ? test_k - BigInt(100) : BigInt(0);
                    BigInt micro_hi = test_k + BigInt(100);
                    for (BigInt micro_k = micro_lo; micro_k <= micro_hi; ++micro_k)
                        check_factors(micro_k);
                }
            }
        }
        cout << "\n[*] Spektral-Scan fertig. Verarbeitet: " << spec_done.load() << " | Treffer: "
             << spec_hits.load() << " | t=" << fixed << setprecision(2) << (omp_get_wtime() - t0_spec) << "s"
             << endl;
    }

    void execute_grid_scan(BigInt center_k) {
        const double noise_floor = 0.0009055397;

        cout << "\n[*] BISECTION-SACKGASSE VERLASSEN. STARTE SPEKTRAL-FÄCHER-SCAN..." << endl;
        cout << "[*] Modus: 2 Mio. Nullstellen | Suche Peak-Abweichung..." << endl;

        // Exponentielle Ringe 10^6 … 10^22 (Bisection-Fehlerkorrektur)
        for (int exponent = 6; exponent <= 22; ++exponent) {
            const BigInt radius = big_pow10(exponent);
            cout << "\n[>] Scan-Ring: +/- 10^" << exponent << " Einheiten..." << endl;

            #pragma omp parallel for schedule(dynamic)
            for (int i = -100; i <= 100; ++i) {
                if (i == 0)
                    continue;

                const unsigned ai = static_cast<unsigned>(std::abs(i));
                const BigInt mag = (radius * BigInt(ai)) / BigInt(100);

                BigInt test_k;
                if (i > 0) {
                    test_k = center_k + mag;
                } else {
                    if (center_k < mag)
                        continue;
                    test_k = center_k - mag;
                }

                const double E = calculate_field_intensity(test_k);

                if (std::fabs(E - noise_floor) > 1e-11) {
                    #pragma omp critical(grid_hit)
                    {
                        cout << "\n[!!!] GRADIENT ERFASST BEI k = " << test_k << " | 10^" << exponent << " | i=" << i
                             << " | E = " << fixed << setprecision(10) << E << endl;
                        check_factors(test_k);
                    }
                }
            }
        }
    }

    void execute_gradient_drill(BigInt center_k) {
        cout << "\n[*] SCHALTE UM AUF PHASENDRIFT-PEILUNG (Berry-Phase)..." << endl;

        BigInt k = center_k;

        for (int i = 0; i < 500; ++i) {
            const double drift = calculate_phase_drift(k);

            if (std::fabs(drift) < 1e-15) {
                cout << "\n[!!!] PHASEN-NULLPUNKT ERREICHT! Starte Deep-Drill bei k = " << k << endl;
                check_factors(k);
                return;
            }

            BigInt step = drift_scaled_step(drift);
            if (step < BigInt(1))
                step = BigInt("1000000");

            if (drift > 0) {
                if (k >= step)
                    k -= step;
                else
                    k = 0;
            } else
                k += step;

            cout << "Step " << i << " | Drift: " << scientific << drift << " | Warp: " << step << "\r" << flush;
        }
        cout << "\n[*] Phasendrill: 500 Schritte — |drift| blieb ≥ 1e-15." << endl;
    }

    void execute_hoffbauer_drill(BigInt center_k) {
        cout << "\n[*] INITIALISIERE HOFFBAUER-TUPEL-SCAN..." << endl;
        cout << "[*] Fokus: Dualitäts-Resonanz (ABCE <-> CEAB)" << endl;
        cout << "[*] center_k " << bigint_tail(center_k, 16) << " | i ∈ [−500000,500000] | 4 Gleise | "
             << omp_get_max_threads() << " Threads" << endl;

        #pragma omp parallel for schedule(dynamic)
        for (long i = -500000; i <= 500000; ++i) {
            BigInt test_k;
            if (i >= 0) {
                test_k = center_k + BigInt(static_cast<unsigned long>(i));
            } else {
                const BigInt mag = BigInt(static_cast<unsigned long>(-i));
                if (center_k < mag)
                    continue;
                test_k = center_k - mag;
            }

            const double e11 = calculate_track_intensity(test_k, 11);
            const double e13 = calculate_track_intensity(test_k, 13);
            const double e17 = calculate_track_intensity(test_k, 17);
            const double e19 = calculate_track_intensity(test_k, 19);

            const double total_resonance = e11 + e13 + e17 + e19;
            const double symmetry_bias = std::fabs((e11 + e19) - (e13 + e17));

            if (total_resonance > 0.005) {
                #pragma omp critical(hoffbauer_hit)
                {
                    cout << "\n[!!!] TUPEL-RESONANZ GEFUNDEN bei k = " << test_k << " | i=" << i << endl;
                    cout << "[>] E_sum=" << fixed << setprecision(10) << total_resonance
                         << " | Symmetrie-Bias=" << scientific << setprecision(6) << symmetry_bias
                         << " | e11,e13,e17,e19=" << fixed << setprecision(8) << e11 << "," << e13 << "," << e17
                         << "," << e19 << endl;
                    check_factors(test_k);
                }
            }
        }
    }

    bool final_lock(BigInt k) {
        for (int r : {11, 13, 17, 19, 23}) {
            BigInt p = 210 * k + r;
            if (p > 1 && N % p == 0) {
                cout << "\n[!!!] HEUREKA! Faktor identifiziert (k = " << k << ")." << endl;
                cout << ">> FAKTOR P IDENTIFIZIERT: " << p << endl;
                return true;
            }
        }
        return false;
    }

    void check_factors(BigInt k) {
        for(int r : {11, 13, 17, 19, 23}) {
            BigInt p = 210 * k + r;
            if (p > 1 && N % p == 0) cout << ">> FAKTOR P IDENTIFIZIERT: " << p << endl;
        }
    }
};

int main() {
    try {
        string n_256 =
            "91240241331123155652145421042111022334156788212344556778890123456789012345679";
        SingularityDrill drill(n_256);
        BigInt k0 = sqrt(BigInt(n_256)) / BigInt(210) - BigInt(100);
        if (k0 < 0) {
            k0 = 0;
        }
        drill.execute_drill(k0);
    } catch (const exception& e) {
        cerr << "Fehler: " << e.what() << endl;
        return 1;
    }
    return 0;
}