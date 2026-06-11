#include <iostream>
#include <gmpxx.h>
#include <mpfr.h>
#include "ZetaNavigation.h"
#include <vector>
#include <set>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <exception>
#ifdef _OPENMP
#include <omp.h>
#endif

// Konstanten des Hofbi-Bohr-Modells
const double KAPPA_H = 2.3457e-305; // Arithmetisch-Geometrische Gravitationskonstante

class HofbiScanner {
public:
    struct ScanCandidate {
        double offset;
        mpf_class energy;
        mpz_class standard_part_value;
    };

    struct ZoomResult {
        double offset;
        mpf_class energy;
        mpf_class gradient;
        mpf_class curvature;
        mpf_class fourth_derivative;
        mpf_class phase_shift;
        double effective_epsilon;
        mpz_class standard_part_value;
        bool converged;
        bool factor_found;
        bool window_exhausted;
    };

    struct AutoLockResult {
        double start_offset;
        double end_offset;
        int steps_taken;
        mpz_class factor;
        mpz_class last_candidate;
        mpf_class start_energy;
        mpf_class end_energy;
        mpf_class last_gradient;
        mpf_class last_fourth_derivative;
        double final_learning_rate;
        int history_points;
        bool success;
        std::string termination_reason;
    };

    struct SondeHistory {
        int step;
        double offset;
        double energy;
        double delta_e;
        double gradient;
        int consecutive_stalls;
        std::string status;
        double learning_rate;
    };

    struct Sonde {
        int id;
        mpf_class current_offset;
        mpf_class learning_rate;
        mpf_class base_learning_rate = 0;
        mpf_class gradient = 0;
        mpf_class curvature = 0;
        std::vector<SondeHistory> history;
        bool active = true;
        bool hit = false;
        int steps_without_gradient = 0;
        std::string status = "tracking";
        bool radar_announced = false;
        mpf_class damping = 0;
        mpf_class lock_threshold = mpf_class("1e-7");

        void adjust_rate(const mpf_class& current_energy, const mpf_class& current_gradient) {
            mpf_class abs_gradient = (current_gradient < 0) ? -current_gradient : current_gradient;
            if (abs_gradient < mpf_class("1e-8")) {
                learning_rate *= mpf_class("1.1");
            } else if (current_energy < mpf_class("0.1")) {
                learning_rate *= mpf_class("0.95");
            }
        }

        void reset_gradient_history() {
            steps_without_gradient = 0;
        }
    };

    mpz_class N_int;
    mpf_class N;
    unsigned int precision_bits;

    HofbiScanner(std::string n_str, unsigned int bits = 2048) {
        precision_bits = bits;
        mpf_set_default_prec(precision_bits);
        N_int = n_str;
        N = n_str;
    }

    // Simulation der Morley-Spannung (Ptolemäische Bindungsenergie)
    mpf_class calculate_morley_energy(mpf_class x) {
        // Wir simulieren einen Vierling an der Stelle x:
        // x, x+2, x+6, x+8. Die alte Ptolemaeus-Formel war eine Identitaet
        // und lieferte daher numerisch fast immer 0. Hier messen wir statt-
        // dessen die diskrete Krummung der lokalen Log-Abstaende.
        mpf_class l1 = log_mpf(x);
        mpf_class l2 = log_shifted_mpf(x, 2);
        mpf_class l3 = log_shifted_mpf(x, 6);
        mpf_class l4 = log_shifted_mpf(x, 8);

        mpf_class d12 = l2 - l1;
        mpf_class d23 = l3 - l2;
        mpf_class d34 = l4 - l3;

        // Diskrete zweite Differenz der Log-Luecken.
        return abs_mpf(d34 - 2 * d23 + d12);
    }

    // Rydberg-Spektroskopie: Suche nach Simuladenvarianz
    ScanCandidate scan_rydberg_boundary(double start_range, double end_range, double step) {
        std::cout << "Starte Hofbi-Rydberg-Spektroskopie..." << std::endl;
        std::cout << "Suche nach Ionisationspunkten (Faktorisierungs-Resonanz)..." << std::endl;
        std::vector<ScanCandidate> top_candidates;

        for (double offset = start_range; offset <= end_range; offset += step) {
            mpf_class test_x = N * offset; // Virtuelle Schalen-Abtastung
            mpf_class energy = calculate_morley_energy(test_x);
            mpz_class candidate_factor = standard_part(test_x);

            update_top_candidates(top_candidates, offset, energy, candidate_factor, 10);

            // Simulation der Rydberg-Grenze: 
            // Wenn die Energie unter einen kritischen Schwellenwert fällt, 
            // detektieren wir eine "Simulade" (einen potenziellen Faktor)
            if (energy < mpf_class("1e-15")) { // Beispielhafter Schwellenwert fuer Resonanz
                std::cout << "[DETEKTION] Resonanz bei Offset: " << offset << std::endl;
                std::cout << "Spektrale Dichte (Span): " << energy << std::endl;
                
                // Hier würde der "Standard Part" (Shadow Casting) greifen:
                mpz_class potential_factor = candidate_factor;
                
                // Validierung am Kern (N)
                if (potential_factor > 1 &&
                    potential_factor < N_int &&
                    mpz_divisible_p(N_int.get_mpz_t(), potential_factor.get_mpz_t())) {
                    std::cout << ">>> FAKTOR GEFUNDEN: " << potential_factor << std::endl;
                }
            }
        }

        print_top_candidates(top_candidates, "Top-Minima aus dem Grobscan");
        return top_candidates.front();
    }

    // Hochaufloesende Nachsuche im detektierten Bereich.
    void scan_high_res_boundary(double start_range, double end_range, double step) {
        std::cout << "\nStarte hochaufloesende Nachsuche..." << std::endl;
        std::cout << "Bereich: [" << start_range << ", " << end_range << "]"
                  << " mit Schrittweite " << step << std::endl;

        const long long steps =
            static_cast<long long>(((end_range - start_range) / step) + 0.5);
        mpf_class best_energy = 0;
        double best_offset = start_range;
        mpz_class best_factor = 0;
        bool found_singularity = false;

        for (long long i = 0; i <= steps; ++i) {
            double offset = start_range + static_cast<double>(i) * step;
            mpf_class test_x = N * offset;
            mpf_class energy = calculate_high_res_morley(test_x);

            if (i == 0 || energy < best_energy) {
                best_energy = energy;
                best_offset = offset;
                best_factor = standard_part(test_x);
            }

            if (is_singular(energy)) {
                mpz_class factor = standard_part(test_x);
                if (factor > 1 &&
                    factor < N_int &&
                    mpz_divisible_p(N_int.get_mpz_t(), factor.get_mpz_t())) {
                    std::cout << "[HOCHAUFLOESUNG] Singulaerer Punkt bei Offset: "
                              << std::setprecision(6) << std::fixed << offset << std::endl;
                    std::cout << "Energie: " << std::scientific << energy << std::endl;
                    std::cout << ">>> FAKTOR GEFUNDEN: " << factor << std::endl;
                    found_singularity = true;
                    break;
                }
            }
        }

        if (!found_singularity) {
            std::cout << "[HOCHAUFLOESUNG] Kein exakter Faktor gefunden." << std::endl;
            std::cout << "Tiefster Trichter bei Offset: "
                      << std::setprecision(6) << std::fixed << best_offset << std::endl;
            std::cout << "Minimale Energie: " << std::scientific << best_energy << std::endl;
            std::cout << "Standard Part: " << best_factor << std::endl;
        }
    }

    // Gradientengefuehrter Zoom um eine vorgegebene Offset-Monade.
    ZoomResult hofbi_singular_zoom(
        double center_offset,
        double window_radius = 1e-6,
        double epsilon = 1e-12,
        int max_iterations = 128
    ) {
        std::cout << "\nInitialisiere Hofbi-Zoom bei Offset "
                  << std::fixed << std::setprecision(12) << center_offset << "..." << std::endl;
        std::cout << "Fenster: [" << (center_offset - window_radius) << ", "
                  << (center_offset + window_radius) << "]"
                  << " | epsilon=" << std::scientific << epsilon << std::endl;

        double left = center_offset - window_radius;
        double right = center_offset + window_radius;
        double current_offset = center_offset;
        bool converged = false;
        bool window_exhausted = false;

        for (int i = 0; i < max_iterations; ++i) {
            DerivativeSnapshot snapshot = evaluate_derivatives(current_offset, epsilon);
            double gradient_value = snapshot.gradient.get_d();
            double curvature_value = snapshot.curvature.get_d();

            if (abs_mpf(snapshot.gradient) < mpf_class("1e-24")) {
                converged = true;
                break;
            }

            if (gradient_value > 0.0) {
                right = current_offset;
            } else {
                left = current_offset;
            }

            double next_offset = 0.5 * (left + right);

            // Gedämpfter Newton-Schritt innerhalb des aktuellen Suchfensters.
            if (std::isfinite(curvature_value) && std::fabs(curvature_value) > 0.0) {
                double newton_offset = current_offset - (gradient_value / curvature_value);
                if (newton_offset > left && newton_offset < right) {
                    next_offset = newton_offset;
                }
            }

            if (std::fabs(next_offset - current_offset) < epsilon) {
                current_offset = next_offset;
                window_exhausted = true;
                break;
            }

            current_offset = next_offset;
        }

        DerivativeSnapshot final_snapshot = evaluate_derivatives(current_offset, epsilon);
        mpz_class factor = standard_part(N * current_offset);
        bool factor_found =
            factor > 1 &&
            factor < N_int &&
            mpz_divisible_p(N_int.get_mpz_t(), factor.get_mpz_t());

        std::cout << "[ZOOM] Finale Offset-Schaetzung: "
                  << std::fixed << std::setprecision(12) << current_offset << std::endl;
        std::cout << "[ZOOM] Energie: " << std::scientific << final_snapshot.energy << std::endl;
        std::cout << "[ZOOM] Gradient: " << final_snapshot.gradient << std::endl;
        std::cout << "[ZOOM] Kruemmung: " << final_snapshot.curvature << std::endl;
        std::cout << "[ZOOM] 4. Ableitung: " << final_snapshot.fourth_derivative << std::endl;
        std::cout << "[ZOOM] Ptolemaeische Phasenverschiebung: "
                  << final_snapshot.phase_shift << std::endl;
        std::cout << "[ZOOM] Effektives epsilon: "
                  << std::scientific << final_snapshot.effective_epsilon << std::endl;
        std::cout << "[ZOOM] Standard Part: " << factor << std::endl;

        if (factor_found) {
            std::cout << ">>> FAKTOR GEFUNDEN: " << factor << std::endl;
        } else if (converged) {
            std::cout << "[ZOOM] Konvergenz erreicht, aber kein exakter Faktor validiert."
                      << std::endl;
        } else if (window_exhausted) {
            std::cout << "[ZOOM] Kein innerer stationaerer Punkt im Zoom-Fenster;"
                      << " Minimum liegt numerisch am Fensterrand." << std::endl;
        } else {
            std::cout << "[ZOOM] Maximale Iterationszahl erreicht; kein exakter Faktor validiert."
                      << std::endl;
        }

        return {
            current_offset,
            final_snapshot.energy,
            final_snapshot.gradient,
            final_snapshot.curvature,
            final_snapshot.fourth_derivative,
            final_snapshot.phase_shift,
            final_snapshot.effective_epsilon,
            factor,
            converged,
            factor_found,
            window_exhausted
        };
    }

    void scan_gradient_sign_changes(
        double center_offset,
        double window_radius,
        double step,
        double epsilon
    ) {
        std::cout << "\nStarte gerichteten Gradienten-Scan..." << std::endl;
        std::cout << "Bereich: [" << std::fixed << std::setprecision(12)
                  << (center_offset - window_radius) << ", "
                  << (center_offset + window_radius) << "]"
                  << " | Schritt=" << std::scientific << step
                  << " | epsilon=" << epsilon << std::endl;

        const long long steps =
            static_cast<long long>(((2 * window_radius) / step) + 0.5);
        double offset = center_offset - window_radius;
        DerivativeSnapshot previous = evaluate_derivatives(offset, epsilon);
        double best_offset = offset;
        mpf_class best_abs_gradient = abs_mpf(previous.gradient);
        bool sign_change_found = false;

        std::cout << "[SIGN] Linker Randgradient: " << previous.gradient << std::endl;

        for (long long i = 1; i <= steps; ++i) {
            offset = center_offset - window_radius + static_cast<double>(i) * step;
            DerivativeSnapshot current = evaluate_derivatives(offset, epsilon);

            if (abs_mpf(current.gradient) < best_abs_gradient) {
                best_abs_gradient = abs_mpf(current.gradient);
                best_offset = offset;
            }

            int prev_sign = sign_of(previous.gradient);
            int current_sign = sign_of(current.gradient);
            if (prev_sign != 0 && current_sign != 0 && prev_sign != current_sign) {
                std::cout << "[SIGN] Vorzeichenwechsel zwischen Offsets "
                          << std::fixed << std::setprecision(12) << (offset - step)
                          << " und " << offset << std::endl;
                std::cout << "[SIGN] Gradienten: " << previous.gradient
                          << " -> " << current.gradient << std::endl;
                std::cout << "[SIGN] 4. Ableitung: " << previous.fourth_derivative
                          << " -> " << current.fourth_derivative << std::endl;
                sign_change_found = true;
            }

            previous = current;
        }

        std::cout << "[SIGN] Rechter Randgradient: " << previous.gradient << std::endl;
        std::cout << "[SIGN] Kleinster Betrag des Gradienten bei Offset "
                  << std::fixed << std::setprecision(12) << best_offset
                  << " mit |grad|=" << std::scientific << best_abs_gradient
                  << std::endl;

        if (!sign_change_found) {
            std::cout << "[SIGN] Kein Gradienten-Vorzeichenwechsel im Fenster detektiert."
                      << std::endl;
        }
    }

    void rydberg_tile_scan(
        double start_offset,
        double tile_width,
        double gradient_step,
        double epsilon,
        int max_tiles = 64
    ) {
        std::cout << "\nStarte kachelweisen Rechts-Scan..." << std::endl;
        double current_start = start_offset;
        bool factor_found = false;
        bool resonance_found = false;

        for (int tile_index = 0; tile_index < max_tiles && current_start <= 1.0; ++tile_index) {
            double current_end = std::min(1.0, current_start + tile_width);
            std::cout << "[KACHEL-SCAN] Teste Fenster: ["
                      << std::fixed << std::setprecision(12) << current_start
                      << ", " << current_end << "]" << std::endl;

            long long steps = static_cast<long long>(((current_end - current_start) / gradient_step) + 0.5);
            if (steps < 1) {
                steps = 1;
            }

            double offset = current_start;
            DerivativeSnapshot previous = evaluate_derivatives(offset, epsilon);
            mpz_class previous_factor = standard_part(N * offset);
            mpf_class tile_best_energy = previous.energy;
            double tile_best_offset = offset;
            mpf_class tile_best_gradient = abs_mpf(previous.gradient);
            double sign_change_left = 0.0;
            double sign_change_right = 0.0;
            bool gradient_sign_changed = false;

            if (is_valid_factor_candidate(previous_factor)) {
                std::cout << ">>> FAKTOR GEFUNDEN AM KACHELSTART: " << previous_factor << std::endl;
                factor_found = true;
                break;
            }

            for (long long i = 1; i <= steps; ++i) {
                offset = current_start + static_cast<double>(i) * gradient_step;
                if (offset > current_end) {
                    offset = current_end;
                }

                DerivativeSnapshot current = evaluate_derivatives(offset, epsilon);
                mpz_class factor_candidate = standard_part(N * offset);

                if (current.energy < tile_best_energy) {
                    tile_best_energy = current.energy;
                    tile_best_offset = offset;
                }
                if (abs_mpf(current.gradient) < tile_best_gradient) {
                    tile_best_gradient = abs_mpf(current.gradient);
                }

                if (is_valid_factor_candidate(factor_candidate)) {
                    std::cout << ">>> FAKTOR GEFUNDEN IN KACHEL: " << factor_candidate
                              << " bei Offset " << std::fixed << std::setprecision(12)
                              << offset << std::endl;
                    factor_found = true;
                    break;
                }

                int prev_sign = sign_of(previous.gradient);
                int current_sign = sign_of(current.gradient);
                if (prev_sign != 0 && current_sign != 0 && prev_sign != current_sign) {
                    sign_change_left = offset - gradient_step;
                    sign_change_right = offset;
                    gradient_sign_changed = true;
                    std::cout << "[RESONANZ] Lokales Minimum in Kachel gefunden zwischen "
                              << std::fixed << std::setprecision(12) << sign_change_left
                              << " und " << sign_change_right << std::endl;
                    std::cout << "[RESONANZ] Gradienten: " << previous.gradient
                              << " -> " << current.gradient << std::endl;
                    break;
                }

                previous = current;
                previous_factor = factor_candidate;
            }

            std::cout << "[KACHEL-SCAN] Tiefste Energie in Kachel bei Offset "
                      << std::fixed << std::setprecision(12) << tile_best_offset
                      << " | Energie=" << std::scientific << tile_best_energy
                      << " | min |grad|=" << tile_best_gradient << std::endl;

            if (factor_found) {
                break;
            }

            if (gradient_sign_changed) {
                resonance_found = true;
                double zoom_center = 0.5 * (sign_change_left + sign_change_right);
                std::cout << "[RESONANZ] Starte Tiefen-Zoom..." << std::endl;
                hofbi_singular_zoom(zoom_center, 0.5 * tile_width, epsilon, 256);
                break;
            }

            current_start += tile_width;
        }

        if (!factor_found && !resonance_found) {
            std::cout << "[KACHEL-SCAN] Kein Faktor und kein lokales Minimum im Scanbereich gefunden."
                      << std::endl;
        }
    }

    void wideband_aperture_scan(
        double start_offset,
        double end_offset,
        double multiplier,
        double epsilon,
        double anomaly_ratio = 0.9995
    ) {
        std::cout << "\nStarte Breitband-Spektroskopie (Aperture Scan)..." << std::endl;
        std::cout << "Bereich: [" << std::scientific << start_offset << ", " << end_offset
                  << "] | Multiplikator=" << multiplier
                  << " | epsilon=" << epsilon << std::endl;

        if (start_offset <= 0.0 || end_offset <= start_offset || multiplier <= 1.0) {
            std::cout << "[WIDEBAND] Ungueltige Scan-Parameter." << std::endl;
            return;
        }

        double current_offset = start_offset;
        DerivativeSnapshot previous = evaluate_derivatives(current_offset, epsilon);
        mpf_class background_energy = previous.energy;
        mpf_class best_energy = previous.energy;
        mpf_class strongest_anomaly_score = 1;
        double best_offset = current_offset;
        double strongest_anomaly_offset = current_offset;
        bool resonance_found = false;
        bool factor_found = false;

        while (current_offset <= end_offset) {
            DerivativeSnapshot snapshot = evaluate_derivatives(current_offset, epsilon);
            mpz_class factor_candidate = standard_part(N * current_offset);

            if (snapshot.energy < best_energy) {
                best_energy = snapshot.energy;
                best_offset = current_offset;
            }

            bool anomalous = snapshot.energy < background_energy * anomaly_ratio;
            bool curvature_flip =
                sign_of(previous.curvature) != 0 &&
                sign_of(snapshot.curvature) != 0 &&
                sign_of(previous.curvature) != sign_of(snapshot.curvature);

            if (background_energy > 0) {
                mpf_class anomaly_score = snapshot.energy / background_energy;
                if (anomaly_score < strongest_anomaly_score) {
                    strongest_anomaly_score = anomaly_score;
                    strongest_anomaly_offset = current_offset;
                }
            }

            if (is_valid_factor_candidate(factor_candidate)) {
                std::cout << ">>> FAKTOR GEFUNDEN IM BREITBAND-SCAN: "
                          << factor_candidate << " bei Offset "
                          << std::fixed << std::setprecision(12) << current_offset << std::endl;
                factor_found = true;
                break;
            }

            if (curvature_flip) {
                std::cout << "[RESONANZ-KANDIDAT] Offset: "
                          << std::fixed << std::setprecision(12) << current_offset << std::endl;
                std::cout << "[RESONANZ-KANDIDAT] Energie=" << std::scientific << snapshot.energy
                          << " | Hintergrund=" << background_energy
                          << " | Kruemmung=" << snapshot.curvature << std::endl;
                std::cout << "[RESONANZ-KANDIDAT] Vorzeichenwechsel der 2. Ableitung detektiert."
                          << std::endl;

                double zoom_radius = std::max(current_offset * (multiplier - 1.0), start_offset);
                hofbi_singular_zoom(current_offset, zoom_radius, epsilon, 256);
                resonance_found = true;
            }

            background_energy = mpf_class("0.95") * background_energy + mpf_class("0.05") * snapshot.energy;
            previous = snapshot;
            current_offset *= multiplier;
        }

        std::cout << "[WIDEBAND] Tiefste Energie bei Offset "
                  << std::fixed << std::setprecision(12) << best_offset
                  << " | Energie=" << std::scientific << best_energy << std::endl;
        std::cout << "[WIDEBAND] Staerkster Energiedip relativ zum Hintergrund bei Offset "
                  << std::fixed << std::setprecision(12) << strongest_anomaly_offset
                  << " | Verhaeltnis=" << std::scientific << strongest_anomaly_score
                  << std::endl;

        if (!factor_found && !resonance_found) {
            if (strongest_anomaly_score < anomaly_ratio) {
                std::cout << "[WIDEBAND] Kein Kruemmungswechsel gefunden; staerkster Kandidat bleibt ein glatter Energiedip."
                          << std::endl;
            } else {
                std::cout << "[WIDEBAND] Kein Resonanzkandidat und kein Faktor im Breitbandfenster gefunden."
                          << std::endl;
            }
        }
    }

    void target_factor_resonance(
        unsigned long target_factor,
        double scan_width,
        double step,
        const mpf_class& ionization_threshold = mpf_class("1e-100")
    ) {
        mpf_class target_offset = mpf_class(target_factor) / N;
        double target_offset_d = target_offset.get_d();
        double start_offset = target_offset_d - scan_width;
        double end_offset = target_offset_d + scan_width;

        std::cout << "\n[TARGET-LOCK] Fokussiere auf " << target_factor
                  << "er-Monade..." << std::endl;
        std::cout << "[TARGET-LOCK] Ziel-Offset: "
                  << std::scientific << std::setprecision(18) << target_offset_d << std::endl;
        std::cout << "[TARGET-LOCK] Scan-Fenster: [" << start_offset
                  << ", " << end_offset << "] | Schritt=" << step << std::endl;

        const long long steps =
            static_cast<long long>(((end_offset - start_offset) / step) + 0.5);
        mpf_class best_energy = 0;
        double best_offset = start_offset;
        mpz_class best_standard_part = 0;
        bool ionization_detected = false;
        bool factor_detected = false;

        for (long long i = 0; i <= steps; ++i) {
            double current_offset = start_offset + static_cast<double>(i) * step;
            mpf_class x = N * current_offset;
            mpf_class energy = calculate_morley_energy(x);
            mpz_class candidate = standard_part(x);

            if (i == 0 || energy < best_energy) {
                best_energy = energy;
                best_offset = current_offset;
                best_standard_part = candidate;
            }

            if (energy < ionization_threshold) {
                std::cout << ">>> IONISATION DETEKTIERT! <<<" << std::endl;
                std::cout << "Resonanz bei Offset: " << std::fixed << std::setprecision(18)
                          << current_offset << std::endl;
                std::cout << "Energie: " << std::scientific << energy << std::endl;
                std::cout << "Standard Part: " << candidate << std::endl;
                ionization_detected = true;
                break;
            }

            if (candidate == mpz_class(target_factor) && is_valid_factor_candidate(candidate)) {
                std::cout << ">>> ZIELFAKTOR DETEKTIERT! <<<" << std::endl;
                std::cout << "Offset: " << std::fixed << std::setprecision(18)
                          << current_offset << std::endl;
                std::cout << "Energie: " << std::scientific << energy << std::endl;
                std::cout << "Standard Part: " << candidate << std::endl;
                factor_detected = true;
                break;
            }
        }

        if (!ionization_detected && !factor_detected) {
            std::cout << "[TARGET-LOCK] Keine Ionisation im Fenster." << std::endl;
            std::cout << "[TARGET-LOCK] Tiefste Energie bei Offset "
                      << std::fixed << std::setprecision(18) << best_offset
                      << " | Energie=" << std::scientific << best_energy
                      << " | Standard Part=" << best_standard_part << std::endl;
        }
    }

    AutoLockResult auto_lock_factor_resonance(
        Sonde& sonde,
        double epsilon,
        int max_steps = 1000,
        double max_step_change = 5e-6,
        bool verbose = true
    ) {
        double start_offset = sonde.current_offset.get_d();
        mpf_class start_energy = calculate_morley_energy(N * sonde.current_offset);
        mpf_class end_energy = start_energy;
        mpf_class last_gradient = 0;
        mpf_class last_fourth_derivative = 0;
        mpz_class last_candidate = round_standard_part(N * sonde.current_offset);
        mpf_class previous_energy = start_energy;
        double starvation_window_delta = 0.0;
        int starvation_window_steps = 0;
        int consecutive_stalls = 0;
        constexpr int hyper_drift_window = 50;
        const double hyper_drift_threshold = epsilon * 100.0;
        const mpf_class curvature_threshold("1e5");

        if (verbose) {
            std::cout << "\n[AUTO-LOCK] Suche gestartet bei Offset: "
                      << std::scientific << std::setprecision(18) << start_offset << std::endl;
            std::cout << "[AUTO-LOCK] learning_rate=" << sonde.learning_rate
                      << " | epsilon=" << epsilon
                      << " | max_step_change=" << max_step_change << std::endl;
        }

        for (int step_index = 0; step_index < max_steps; ++step_index) {
            mpf_class x = N * sonde.current_offset;
            mpz_class candidate = round_standard_part(x);
            mpf_class energy = calculate_morley_energy(x);
            DerivativeSnapshot derivative_snapshot = evaluate_derivatives(sonde.current_offset.get_d(), epsilon);
            end_energy = energy;
            last_gradient = derivative_snapshot.gradient;
            last_fourth_derivative = derivative_snapshot.fourth_derivative;
            last_candidate = candidate;
            sonde.gradient = derivative_snapshot.gradient;
            sonde.curvature = derivative_snapshot.curvature;

            mpf_class dist_to_q = abs_mpf(sonde.current_offset - mpf_class("0.058826"));
            if (dist_to_q < mpf_class("0.001") && !sonde.radar_announced) {
                std::cout << "[RADAR] Sonde " << sonde.id
                          << " naehrt sich der 61.681-Monade!" << std::endl;
                sonde.radar_announced = true;
            }

            mpf_class abs_gradient = abs_mpf(derivative_snapshot.gradient);
            if (abs_gradient < mpf_class("1e-8")) {
                ++sonde.steps_without_gradient;
            } else {
                sonde.reset_gradient_history();
            }

            if (is_valid_factor_candidate(candidate)) {
                sonde.hit = true;
                sonde.active = false;
                sonde.status = "factor_found";
                if (verbose) {
                    std::cout << ">>> IONISATION ERFOLGREICH! FAKTOR: " << candidate << std::endl;
                    std::cout << "[AUTO-LOCK] Offset: "
                              << std::fixed << std::setprecision(18)
                              << sonde.current_offset.get_d() << std::endl;
                    std::cout << "[AUTO-LOCK] Energie: " << std::scientific << energy << std::endl;
                }
                return {
                    start_offset,
                    sonde.current_offset.get_d(),
                    step_index,
                    candidate,
                    candidate,
                    start_energy,
                    energy,
                    last_gradient,
                    last_fourth_derivative,
                    sonde.learning_rate.get_d(),
                    static_cast<int>(sonde.history.size()),
                    true,
                    "factor_found"
                };
            }

            sonde.adjust_rate(energy, derivative_snapshot.gradient);
            mpf_class abs_curvature = abs_mpf(sonde.curvature);
            mpf_class delta_energy = abs_mpf(energy - previous_energy);
            update_sonde_speed(sonde, N, verbose);
            apply_re_entry_damping(sonde, energy, derivative_snapshot.gradient, epsilon, verbose);
            if (sonde.status == "searching_wave" || sonde.status == "zeta_surfing" ||
                sonde.status == "ZETA_BOOST") {
                if (abs_curvature > curvature_threshold) {
                    sonde.learning_rate = sonde.base_learning_rate;
                    sonde.status = "re_entry_homing";
                    if (verbose) {
                        std::cout << "[LOCK] Sonde " << sonde.id
                                  << " bremst in Monade ein bei Offset: "
                                  << std::fixed << std::setprecision(18)
                                  << sonde.current_offset.get_d() << std::endl;
                    }
                }
            }
            if (sonde.status == "hyper_drift") {
                if (abs_curvature > curvature_threshold && energy < mpf_class("0.1")) {
                    sonde.learning_rate *= mpf_class("0.1");
                    sonde.status = "re_entry";
                    consecutive_stalls = 0;
                    if (verbose) {
                        std::cout << "[LOCK] Sonde " << sonde.id
                                  << " hat Struktur gefasst."
                                  << " | |curvature|=" << abs_curvature
                                  << " | Lernrate=" << sonde.learning_rate << std::endl;
                    }
                } else {
                    sonde.learning_rate *= mpf_class("1.1");
                }
            }
            if (consecutive_stalls > 10 &&
                sonde.status != "relativistic_boost" &&
                sonde.status != "zeta_surfing" &&
                sonde.status != "ZETA_BOOST") {
                sonde.learning_rate *= mpf_class("2.0");
                sonde.status = "relativistic_boost";
                if (verbose) {
                    std::cout << "[TELEMETRIE] Sonde " << sonde.id
                              << " zuendet Boost bei Step " << step_index
                              << " | consecutive_stalls=" << consecutive_stalls
                              << " | Lernrate=" << sonde.learning_rate << std::endl;
                }
            }
            if (sonde.status == "relativistic_boost" && delta_energy > mpf_class("1e-10")) {
                sonde.learning_rate = mpf_class("1e-15");
                sonde.status = "re_entry";
                consecutive_stalls = 0;
                if (verbose) {
                    std::cout << "[BREMSFALLSCHIRM] Sonde " << sonde.id
                              << " geht in re_entry."
                              << " | delta_E=" << delta_energy
                              << std::endl;
                }
            }

            double delta = 0.0;
            if (sonde.status == "target_lock_q") {
                const double gradient_value = derivative_snapshot.gradient.get_d();
                const double damping_value = sonde.damping.get_d();
                const double denominator = std::fabs(gradient_value) + damping_value;
                if (denominator > 0.0) {
                    delta = -sonde.learning_rate.get_d() * (gradient_value / denominator);
                }
            } else {
                delta = -sonde.learning_rate.get_d() * derivative_snapshot.gradient.get_d();
            }

            if (!std::isfinite(delta)) {
                sonde.active = false;
                sonde.status = "invalid_step";
                if (verbose) {
                    std::cout << "[AUTO-LOCK] Numerischer Abbruch: ungueltiger Schritt." << std::endl;
                }
                return {
                    start_offset,
                    sonde.current_offset.get_d(),
                    step_index,
                    0,
                    candidate,
                    start_energy,
                    energy,
                    last_gradient,
                    last_fourth_derivative,
                    sonde.learning_rate.get_d(),
                    static_cast<int>(sonde.history.size()),
                    false,
                    "invalid_step"
                };
            }

            delta = std::clamp(delta, -max_step_change, max_step_change);
            double next_offset = std::clamp(sonde.current_offset.get_d() + delta, epsilon, 1.0);
            double applied_delta = std::fabs(next_offset - sonde.current_offset.get_d());

            starvation_window_delta += applied_delta;
            ++starvation_window_steps;
            if (step_index > hyper_drift_window &&
                starvation_window_steps >= hyper_drift_window &&
                starvation_window_delta < hyper_drift_threshold) {
                ++consecutive_stalls;
                if (!apply_secondary_zeta_boost(sonde, N, consecutive_stalls, verbose)) {
                    sonde.learning_rate *= mpf_class("1.5");
                    sonde.status = "hyper_drift";
                    if (verbose) {
                        std::cout << "[HYPER-DRIFT] Sonde verhungert im flachen Raum."
                                  << " | Delta-Fenster=" << starvation_window_delta
                                  << " | consecutive_stalls=" << consecutive_stalls
                                  << " | neue Lernrate=" << sonde.learning_rate << std::endl;
                    }
                } else {
                    consecutive_stalls = 0;
                }
                starvation_window_delta = 0.0;
                starvation_window_steps = 0;
            } else if (starvation_window_steps >= hyper_drift_window) {
                consecutive_stalls = 0;
                starvation_window_delta = 0.0;
                starvation_window_steps = 0;
            }

            sonde.history.push_back({
                step_index,
                sonde.current_offset.get_d(),
                energy.get_d(),
                delta_energy.get_d(),
                derivative_snapshot.gradient.get_d(),
                consecutive_stalls,
                sonde.status,
                sonde.learning_rate.get_d()
            });

            if (verbose && step_index % 100 == 0) {
                std::cout << "[AUTO-LOCK] Schritt " << step_index
                          << " | Offset: " << std::fixed << std::setprecision(18)
                          << sonde.current_offset.get_d()
                          << " | Energie: " << std::scientific << energy
                          << " | Gradient: " << derivative_snapshot.gradient
                          << " | Kandidat: " << candidate
                          << " | Lernrate: " << sonde.learning_rate
                          << std::endl;
            }

            if (std::fabs(next_offset - sonde.current_offset.get_d()) < epsilon) {
                if (sonde.status != "accelerating" && sonde.status != "relativistic_boost") {
                    sonde.status = "stalled";
                    if (verbose) {
                        std::cout << "[AUTO-LOCK] Schrittweite unter epsilon gefallen." << std::endl;
                        std::cout << "[AUTO-LOCK] Letzter Kandidat: " << candidate
                                  << " | Offset: " << std::fixed << std::setprecision(18)
                                  << sonde.current_offset.get_d() << std::endl;
                        std::cout << "[TELEPORT-IMPULS] Lernrate wird massiv erhoeht." << std::endl;
                    }
                    sonde.learning_rate *= mpf_class("100.0");
                    sonde.status = "accelerating";
                    sonde.reset_gradient_history();
                    continue;
                }

                sonde.active = false;
                sonde.status = "stalled";
                if (verbose) {
                    std::cout << "[AUTO-LOCK] Schrittweite unter epsilon gefallen." << std::endl;
                    std::cout << "[AUTO-LOCK] Letzter Kandidat: " << candidate
                              << " | Offset: " << std::fixed << std::setprecision(18)
                              << sonde.current_offset.get_d() << std::endl;
                }
                return {
                    start_offset,
                    sonde.current_offset.get_d(),
                    step_index,
                    0,
                    candidate,
                    start_energy,
                    energy,
                    last_gradient,
                    last_fourth_derivative,
                    sonde.learning_rate.get_d(),
                    static_cast<int>(sonde.history.size()),
                    false,
                    "stalled"
                };
            }

            if (sonde.status == "re_entry" || sonde.status == "re_entry_homing") {
                sonde.status = "tracking";
            } else if (sonde.status != "hyper_drift" &&
                       sonde.status != "zeta_surfing" &&
                       sonde.status != "ZETA_BOOST" &&
                       sonde.status != "stalling_prevention" &&
                       sonde.status != "searching_wave" &&
                       sonde.status != "target_lock_q" &&
                       sonde.status != "accelerating" &&
                       sonde.status != "relativistic_boost") {
                sonde.status = "tracking";
            }
            sonde.current_offset = next_offset;
            previous_energy = energy;
        }

        sonde.active = false;
        sonde.status = "max_steps";
        if (verbose) {
            std::cout << "[AUTO-LOCK] Maximale Schrittzahl erreicht, ohne Faktor zu validieren."
                      << std::endl;
        }
        return {
            start_offset,
            sonde.current_offset.get_d(),
            max_steps,
            0,
            last_candidate,
            start_energy,
            end_energy,
            last_gradient,
            last_fourth_derivative,
            sonde.learning_rate.get_d(),
            static_cast<int>(sonde.history.size()),
            false,
            "max_steps"
        };
    }

    AutoLockResult auto_lock_factor_resonance(
        double start_offset,
        double learning_rate,
        double epsilon,
        int max_steps = 1000,
        double max_step_change = 5e-6,
        bool verbose = true
    ) {
        Sonde sonde{
            -1,
            mpf_class(start_offset),
            mpf_class(learning_rate),
            mpf_class(learning_rate),
            0,
            0,
            {},
            true,
            false,
            0,
            "tracking",
            false
        };
        return auto_lock_factor_resonance(sonde, epsilon, max_steps, max_step_change, verbose);
    }

    std::vector<double> generate_optimized_hofbi_seeds(
        std::size_t seed_count,
        std::size_t sqrt_cluster_points = 5
    ) {
        std::vector<double> seeds;
        if (seed_count == 0) {
            return seeds;
        }

        mpfr_t n_fr;
        mpfr_t sqrt_n_fr;
        mpfr_t min_offset_fr;
        mpfr_t max_offset_fr;
        mpfr_t log_min_fr;
        mpfr_t log_max_fr;
        mpfr_t current_log_fr;
        mpfr_t seed_fr;

        mpfr_init2(n_fr, precision_bits);
        mpfr_init2(sqrt_n_fr, precision_bits);
        mpfr_init2(min_offset_fr, precision_bits);
        mpfr_init2(max_offset_fr, precision_bits);
        mpfr_init2(log_min_fr, precision_bits);
        mpfr_init2(log_max_fr, precision_bits);
        mpfr_init2(current_log_fr, precision_bits);
        mpfr_init2(seed_fr, precision_bits);

        mpfr_set_f(n_fr, N.get_mpf_t(), MPFR_RNDN);
        mpfr_sqrt(sqrt_n_fr, n_fr, MPFR_RNDN);

        // Goldloeckchen-Grenze: 1 / sqrt(N)
        mpfr_div(max_offset_fr, sqrt_n_fr, n_fr, MPFR_RNDN);

        // Untere Grenze: 3 / N
        mpfr_set_ui(min_offset_fr, 3ul, MPFR_RNDN);
        mpfr_div(min_offset_fr, min_offset_fr, n_fr, MPFR_RNDN);

        if (mpfr_cmp(min_offset_fr, max_offset_fr) >= 0) {
            seeds.push_back(mpf_from_mpfr(min_offset_fr).get_d());
            mpfr_clear(n_fr);
            mpfr_clear(sqrt_n_fr);
            mpfr_clear(min_offset_fr);
            mpfr_clear(max_offset_fr);
            mpfr_clear(log_min_fr);
            mpfr_clear(log_max_fr);
            mpfr_clear(current_log_fr);
            mpfr_clear(seed_fr);
            return seeds;
        }

        if (seed_count == 1) {
            seeds.push_back(mpf_from_mpfr(min_offset_fr).get_d());
            mpfr_clear(n_fr);
            mpfr_clear(sqrt_n_fr);
            mpfr_clear(min_offset_fr);
            mpfr_clear(max_offset_fr);
            mpfr_clear(log_min_fr);
            mpfr_clear(log_max_fr);
            mpfr_clear(current_log_fr);
            mpfr_clear(seed_fr);
            return seeds;
        }

        mpfr_log(log_min_fr, min_offset_fr, MPFR_RNDN);
        mpfr_log(log_max_fr, max_offset_fr, MPFR_RNDN);

        mpf_class log_min = mpf_from_mpfr(log_min_fr);
        mpf_class log_max = mpf_from_mpfr(log_max_fr);
        mpf_class step = (log_max - log_min) / mpf_class(static_cast<unsigned long>(seed_count - 1));

        for (std::size_t i = 0; i < seed_count; ++i) {
            mpf_class current_log = log_min + mpf_class(static_cast<unsigned long>(i)) * step;
            mpfr_set_f(current_log_fr, current_log.get_mpf_t(), MPFR_RNDN);
            mpfr_exp(seed_fr, current_log_fr, MPFR_RNDN);
            seeds.push_back(mpf_from_mpfr(seed_fr).get_d());
        }

        add_sqrt_cluster(seeds, max_offset_fr, sqrt_cluster_points);
        std::sort(seeds.begin(), seeds.end());
        seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());

        mpfr_clear(n_fr);
        mpfr_clear(sqrt_n_fr);
        mpfr_clear(min_offset_fr);
        mpfr_clear(max_offset_fr);
        mpfr_clear(log_min_fr);
        mpfr_clear(log_max_fr);
        mpfr_clear(current_log_fr);
        mpfr_clear(seed_fr);

        return seeds;
    }

    std::set<mpz_class> run_multi_target_array(
        std::size_t seed_count,
        double learning_rate,
        double epsilon,
        int max_steps = 1000,
        double max_step_change = 5e-6,
        std::size_t sqrt_cluster_points = 5
    ) {
        std::vector<double> start_offsets =
            generate_optimized_hofbi_seeds(seed_count, sqrt_cluster_points);
        if (start_offsets.size() >= 14) {
            static const double probe_seed_matrix[] = {
                0.010000,
                0.025000,
                0.040000,
                0.050000,
                0.055000,
                0.057000,
                0.058000,
                0.058500,
                0.058800,
                0.059000,
                0.062000,
                0.065000,
                0.080000,
                0.100000
            };
            for (std::size_t i = 0; i < 14; ++i) {
                start_offsets[i] = probe_seed_matrix[i];
            }
            start_offsets.resize(14);
            std::cout << "[ARRAY] Sonden-Matrix 0-13 aktiviert:"
                      << " 0.010000, 0.025000, 0.040000, 0.050000,"
                      << " 0.055000, 0.057000, 0.058000, 0.058500,"
                      << " 0.058800, 0.059000, 0.062000, 0.065000,"
                      << " 0.080000, 0.100000" << std::endl;
        }
        std::set<mpz_class> found_factors;
        std::vector<AutoLockResult> results(start_offsets.size());
        std::vector<Sonde> sondes(start_offsets.size());

        std::cout << "\n[ARRAY] Starte Mehrsonden-Modus mit "
                  << start_offsets.size() << " Seeds." << std::endl;
        if (!start_offsets.empty()) {
            std::cout << "[ARRAY] Hofbi-Seed-Bereich: ["
                      << std::scientific << std::setprecision(18)
                      << start_offsets.front() << ", " << start_offsets.back() << "]"
                      << std::endl;
        }

#pragma omp parallel for if(start_offsets.size() > 1)
        for (int i = 0; i < static_cast<int>(start_offsets.size()); ++i) {
            sondes[static_cast<std::size_t>(i)] = Sonde{
                i,
                mpf_class(start_offsets[static_cast<std::size_t>(i)]),
                mpf_class(learning_rate),
                mpf_class(learning_rate),
                0,
                0,
                {},
                true,
                false,
                0,
                "tracking",
                false
            };
            apply_deep_space_launch_boost(
                sondes[static_cast<std::size_t>(i)],
                N,
                false
            );
            AutoLockResult result = auto_lock_factor_resonance(
                sondes[static_cast<std::size_t>(i)],
                epsilon,
                max_steps,
                max_step_change,
                false
            );
            results[static_cast<std::size_t>(i)] = result;

            if (result.factor > 1) {
#pragma omp critical
                {
                    if (found_factors.find(result.factor) == found_factors.end()) {
                        std::cout << "[ARRAY-TREFFER] Sonde " << i
                                  << " fand Faktor: " << result.factor
                                  << " | Start-Offset="
                                  << std::scientific << std::setprecision(18)
                                  << start_offsets[static_cast<std::size_t>(i)] << std::endl;
                        found_factors.insert(result.factor);
                    }
                }
            }
        }

        print_probe_table(results);
        export_trajectories_to_csv(sondes);

        if (found_factors.empty()) {
            std::cout << "[ARRAY] Keine Faktoren im Mehrsonden-Modus gefunden." << std::endl;
        }

        return found_factors;
    }

private:
    void export_trajectories_to_csv(
        const std::vector<Sonde>& sondes,
        const std::string& filename = "trajectories.csv"
    ) {
        std::ofstream file(filename);
        if (!file) {
            std::cout << "[EXPORT] Konnte " << filename << " nicht schreiben." << std::endl;
            return;
        }

        file << "ProbeID,Step,Offset,Energy,DeltaE,Gradient,ConsecutiveStalls,Status,LearningRate,TargetLockQ\n";
        file << std::scientific << std::setprecision(18);

        for (const auto& sonde : sondes) {
            for (const auto& history_entry : sonde.history) {
                file << sonde.id << ","
                     << history_entry.step << ","
                     << history_entry.offset << ","
                     << history_entry.energy << ","
                     << history_entry.delta_e << ","
                     << history_entry.gradient << ","
                     << history_entry.consecutive_stalls << ","
                     << history_entry.status << ","
                     << history_entry.learning_rate << ","
                     << ((history_entry.status == "target_lock_q") ? 1 : 0) << "\n";
            }
        }

        file.close();
        std::cout << "[EXPORT] Trajektorien in " << filename << " gespeichert." << std::endl;
    }

    struct DerivativeSnapshot {
        mpf_class energy;
        mpf_class gradient;
        mpf_class curvature;
        mpf_class fourth_derivative;
        mpf_class phase_shift;
        double effective_epsilon;
    };

    // Hilfsfunktionen für mpf/mpfr
    mpf_class mpf_from_mpfr(const mpfr_t& value) {
        mpf_class result(0, precision_bits);
        mpfr_get_f(result.get_mpf_t(), value, MPFR_RNDN);
        return result;
    }
    mpf_class log_mpf(const mpf_class& x) {
        mpfr_t x_fr;
        mpfr_t log_fr;
        mpfr_init2(x_fr, precision_bits);
        mpfr_init2(log_fr, precision_bits);

        mpfr_set_f(x_fr, x.get_mpf_t(), MPFR_RNDN);
        mpfr_log(log_fr, x_fr, MPFR_RNDN);

        mpf_class result = mpf_from_mpfr(log_fr);
        mpfr_clear(x_fr);
        mpfr_clear(log_fr);
        return result;
    }
    mpf_class abs_mpf(mpf_class x) {
        return (x < 0) ? -x : x;
    }
    double calculate_zeta_momentum(double x) {
        const std::vector<double>& gammas = get_riemann_gammas();
        const std::size_t term_count = std::min<std::size_t>(get_zeta_momentum_max_terms(), gammas.size());
        if (x <= 1.0 || term_count == 0) {
            return 0.0;
        }

        double log_x = std::log(x);
        double interference = 0.0;
        for (std::size_t i = 0; i < term_count; ++i) {
            interference += std::cos(gammas[i] * log_x);
        }
        return interference / static_cast<double>(term_count);
    }
    void apply_re_entry_damping(
        Sonde& s,
        const mpf_class& current_energy,
        const mpf_class& current_gradient,
        double epsilon,
        bool verbose = false
    ) {
        static const mpf_class eta_61681("3.9321e-06");
        static const mpf_class selective_lock_threshold("1e-6");

        if (s.status != "zeta_surfing" && s.status != "searching_wave" &&
            s.status != "ZETA_BOOST") {
            return;
        }

        if (s.current_offset > mpf_class("0.04")) {
            s.damping = eta_61681;
            s.lock_threshold = selective_lock_threshold;
        }

        if (current_energy < s.lock_threshold) {
            s.learning_rate = mpf_class(epsilon);
            s.damping = eta_61681;
            s.status = "target_lock_q";

            if (verbose) {
                std::cout << "[LOCK] Sonde " << s.id
                          << " tritt in Ereignishorizont der 61.681 ein."
                          << " | Energie=" << current_energy
                          << " | Lock-Schwelle=" << s.lock_threshold
                          << " | Gradient=" << current_gradient
                          << std::endl;
                std::cout << "[LOCK] Daempfung aktiviert: eta=" << eta_61681
                          << " | Offset="
                          << std::fixed << std::setprecision(18)
                          << s.current_offset.get_d()
                          << std::endl;
            }
        }
    }
    void update_sonde_speed(Sonde& s, const mpf_class& current_N, bool verbose = false) {
        const mpf_class fallback_learning_rate("1e-12");
        const mpf_class baseline = (s.base_learning_rate > 0) ? s.base_learning_rate : fallback_learning_rate;
        const mpf_class min_learning_rate = baseline * mpf_class("0.25");
        const mpf_class max_learning_rate = baseline * mpf_class("12.0");
        auto clamp_learning_rate = [&](const mpf_class& value) {
            if (value < min_learning_rate) {
                return min_learning_rate;
            }
            if (value > max_learning_rate) {
                return max_learning_rate;
            }
            return value;
        };

        mpf_class x_val_mpf = current_N * s.current_offset;
        double x_val = x_val_mpf.get_d();
        if (!std::isfinite(x_val) || x_val <= 1.0) {
            s.learning_rate = clamp_learning_rate(s.learning_rate);
            if (s.status != "factor_found") {
                s.status = "searching_wave";
            }
            if (verbose) {
                std::cout << "[ZETA] Sonde " << s.id
                          << " ueberspringt Surfing, da x=" << x_val
                          << " nicht im log-Bereich liegt." << std::endl;
            }
            return;
        }

        double momentum = calculate_zeta_momentum(x_val);
        if (!std::isfinite(momentum)) {
            momentum = 0.0;
        }

        if (s.id >= 9 && s.status != "factor_found") {
            if (momentum > 0.05) {
                s.learning_rate = clamp_learning_rate(
                    s.learning_rate * mpf_class(1.0 + 10.0 * momentum)
                );
                s.status = "zeta_surfing";
                if (verbose) {
                    std::cout << "[ZETA] Sonde " << s.id
                              << " surft mit momentum=" << momentum
                              << " | Lernrate=" << s.learning_rate << std::endl;
                }
            } else {
                s.learning_rate = clamp_learning_rate(s.learning_rate * mpf_class("1.5"));
                s.status = "searching_wave";
            }
        } else if (s.status == "stalled") {
            if (momentum > 0.05) {
                s.learning_rate = clamp_learning_rate(
                    s.learning_rate * mpf_class(1.0 + 10.0 * momentum)
                );
                s.status = "zeta_surfing";
                if (verbose) {
                    std::cout << "[ZETA] Sonde " << s.id
                              << " surft mit momentum=" << momentum
                              << " | Lernrate=" << s.learning_rate << std::endl;
                }
            } else {
                s.learning_rate = clamp_learning_rate(s.learning_rate * mpf_class("1.5"));
            }
        }
    }
    void apply_deep_space_launch_boost(Sonde& s, const mpf_class& current_N, bool verbose = false) {
        if (s.current_offset <= mpf_class("0.01")) {
            return;
        }

        s.base_learning_rate = mpf_class("1e-7");
        s.learning_rate = s.base_learning_rate;

        const mpf_class x_val_mpf = current_N * s.current_offset;
        const double x_val = x_val_mpf.get_d();
        double momentum = calculate_zeta_momentum(x_val);
        if (!std::isfinite(momentum)) {
            momentum = 0.0;
        }

        if (momentum > 0.02) {
            s.learning_rate *= mpf_class("100.0");
            s.base_learning_rate = s.learning_rate;
            s.status = "ZETA_BOOST";
        } else {
            s.status = "deep_space_cruise";
        }

        if (verbose) {
            std::cout << "[DEEP-SPACE] Sonde " << s.id
                      << " | Offset=" << std::fixed << std::setprecision(6)
                      << s.current_offset.get_d()
                      << " | momentum=" << momentum
                      << " | Lernrate=" << s.learning_rate
                      << " | Status=" << s.status
                      << std::endl;
        }
    }
    bool apply_secondary_zeta_boost(
        Sonde& s,
        const mpf_class& current_N,
        int consecutive_stalls,
        bool verbose = false
    ) {
        if (s.current_offset <= mpf_class("0.04") || consecutive_stalls < 2) {
            return false;
        }

        const mpf_class x_val_mpf = current_N * s.current_offset;
        const double x_val = x_val_mpf.get_d();
        double momentum = calculate_zeta_momentum(x_val);
        if (!std::isfinite(momentum) || momentum <= 0.02) {
            return false;
        }

        s.learning_rate *= mpf_class("100.0");
        if (s.learning_rate > s.base_learning_rate) {
            s.base_learning_rate = s.learning_rate;
        }
        s.status = "ZETA_BOOST";

        if (verbose) {
            std::cout << "[RE-BOOST] Sonde " << s.id
                      << " verlaesst Stagnation erneut per ZETA_BOOST."
                      << " | consecutive_stalls=" << consecutive_stalls
                      << " | momentum=" << momentum
                      << " | Lernrate=" << s.learning_rate
                      << std::endl;
        }

        return true;
    }
    mpf_class log_shifted_mpf(const mpf_class& x, unsigned long shift) {
        mpfr_t relative_shift_fr;
        mpfr_t log1p_fr;
        mpfr_init2(relative_shift_fr, precision_bits);
        mpfr_init2(log1p_fr, precision_bits);

        mpf_class relative_shift = mpf_class(shift) / x;
        mpfr_set_f(relative_shift_fr, relative_shift.get_mpf_t(), MPFR_RNDN);
        mpfr_log1p(log1p_fr, relative_shift_fr, MPFR_RNDN);

        mpf_class result = log_mpf(x) + mpf_from_mpfr(log1p_fr);
        mpfr_clear(relative_shift_fr);
        mpfr_clear(log1p_fr);
        return result;
    }
    mpf_class calculate_high_res_morley(const mpf_class& x) {
        return calculate_morley_energy(x);
    }
    void add_sqrt_cluster(
        std::vector<double>& seeds,
        const mpfr_t& central_offset_fr,
        std::size_t points_in_cluster,
        double spread = 1e-4
    ) {
        if (points_in_cluster == 0) {
            return;
        }

        double central_offset = mpf_from_mpfr(central_offset_fr).get_d();
        std::cout << "[CLUSTER] Setze Root-Cluster um Offset: "
                  << std::scientific << std::setprecision(18) << central_offset << std::endl;

        if (points_in_cluster == 1) {
            if (central_offset > 0.0) {
                seeds.push_back(central_offset);
            }
            return;
        }

        double step = spread / static_cast<double>(points_in_cluster);
        long long half_span = static_cast<long long>(points_in_cluster / 2);

        for (long long i = -half_span; i <= half_span; ++i) {
            double offset = central_offset + static_cast<double>(i) * step;
            if (offset > 0.0 && offset <= 1.0) {
                seeds.push_back(offset);
            }
        }
    }
    std::string energy_trend_label(const mpf_class& start_energy, const mpf_class& end_energy) {
        mpf_class diff = end_energy - start_energy;
        mpf_class tolerance = abs_mpf(start_energy) * mpf_class("1e-9") + mpf_class("1e-18");
        if (diff < -tolerance) {
            return "sinkend";
        }
        if (diff > tolerance) {
            return "steigend";
        }
        return "stabil";
    }
    void print_probe_table(const std::vector<AutoLockResult>& results) {
        std::cout << "\n[ARRAY-TABELLE]" << std::endl;
        std::cout << "Sonde | Start-Offset | End-Offset | Schritte | Delta-Offset | Letzter Kandidat | Faktor | Status | Energie-Tendenz | End-4.Abl. | Lernrate-Ende | Historie" << std::endl;

        for (std::size_t i = 0; i < results.size(); ++i) {
            const AutoLockResult& result = results[i];
            double delta_offset = result.end_offset - result.start_offset;
            std::string factor_text = (result.factor > 1) ? result.factor.get_str() : "-";

            std::cout << std::setw(5) << i
                      << " | " << std::scientific << std::setprecision(6) << result.start_offset
                      << " | " << std::scientific << std::setprecision(6) << result.end_offset
                      << " | " << std::setw(8) << result.steps_taken
                      << " | " << std::scientific << std::setprecision(6) << delta_offset
                      << " | " << result.last_candidate
                      << " | " << factor_text
                      << " | " << result.termination_reason
                      << " | " << energy_trend_label(result.start_energy, result.end_energy)
                      << " | " << std::scientific << std::setprecision(6)
                      << result.last_fourth_derivative
                      << " | " << std::scientific << std::setprecision(6)
                      << result.final_learning_rate
                      << " | " << result.history_points
                      << std::endl;
        }
    }
    bool is_singular(const mpf_class& energy) {
        return energy < mpf_class("1e-18");
    }
    DerivativeSnapshot evaluate_derivatives(double offset, double epsilon) {
        double effective_epsilon = epsilon;

        for (int attempt = 0; attempt < 8; ++attempt) {
            mpf_class x0 = N * offset;
            mpf_class h = N * effective_epsilon;

            mpf_class e_m2 = calculate_high_res_morley(x0 - 2 * h);
            mpf_class e_m1 = calculate_high_res_morley(x0 - h);
            mpf_class e_0 = calculate_high_res_morley(x0);
            mpf_class e_p1 = calculate_high_res_morley(x0 + h);
            mpf_class e_p2 = calculate_high_res_morley(x0 + 2 * h);

            mpf_class gradient = (e_p1 - e_m1) / (2 * effective_epsilon);
            mpf_class curvature =
                (e_p1 - 2 * e_0 + e_m1) / (effective_epsilon * effective_epsilon);
            mpf_class fourth_derivative =
                (e_m2 - 4 * e_m1 + 6 * e_0 - 4 * e_p1 + e_p2) /
                (effective_epsilon * effective_epsilon * effective_epsilon * effective_epsilon);
            mpf_class phase_shift = abs_mpf(gradient) / (1 + abs_mpf(fourth_derivative));

            if (abs_mpf(curvature) > 0 || abs_mpf(gradient) > 0) {
                return {
                    e_0,
                    gradient,
                    curvature,
                    fourth_derivative,
                    phase_shift,
                    effective_epsilon
                };
            }

            effective_epsilon *= 10.0;
        }

        mpf_class x0 = N * offset;
        mpf_class e_0 = calculate_high_res_morley(x0);
        return {e_0, 0, 0, 0, 0, effective_epsilon};
    }
    mpz_class standard_part(const mpf_class& x) {
        mpz_class integer_part;
        mpz_set_f(integer_part.get_mpz_t(), x.get_mpf_t());
        return integer_part;
    }
    mpz_class round_standard_part(const mpf_class& x) {
        mpfr_t x_fr;
        mpfr_init2(x_fr, precision_bits);
        mpfr_set_f(x_fr, x.get_mpf_t(), MPFR_RNDN);

        mpz_class rounded;
        mpfr_get_z(rounded.get_mpz_t(), x_fr, MPFR_RNDN);
        mpfr_clear(x_fr);
        return rounded;
    }
    bool is_valid_factor_candidate(const mpz_class& candidate) {
        return candidate > 1 &&
               candidate < N_int &&
               mpz_divisible_p(N_int.get_mpz_t(), candidate.get_mpz_t());
    }
    int sign_of(const mpf_class& value) {
        if (value > 0) {
            return 1;
        }
        if (value < 0) {
            return -1;
        }
        return 0;
    }
    void update_top_candidates(
        std::vector<ScanCandidate>& candidates,
        double offset,
        const mpf_class& energy,
        const mpz_class& standard_value,
        std::size_t limit
    ) {
        candidates.push_back({offset, energy, standard_value});
        std::sort(
            candidates.begin(),
            candidates.end(),
            [](const ScanCandidate& a, const ScanCandidate& b) {
                return a.energy < b.energy;
            }
        );

        if (candidates.size() > limit) {
            candidates.resize(limit);
        }
    }
    void print_top_candidates(
        const std::vector<ScanCandidate>& candidates,
        const std::string& title
    ) {
        std::cout << "\n" << title << ":" << std::endl;
        for (std::size_t i = 0; i < candidates.size(); ++i) {
            const ScanCandidate& candidate = candidates[i];
            std::cout << "  [" << (i + 1) << "] Offset="
                      << std::fixed << std::setprecision(6) << candidate.offset
                      << " | Energie=" << std::scientific << candidate.energy
                      << " | Standard Part=" << candidate.standard_part_value
                      << std::endl;
        }
    }
};

int main(int argc, char** argv) {
    std::string rsa_n = "1048577";
    int target_factor = 17;
    std::size_t zeta_terms = DEFAULT_ZETA_MOMENTUM_MAX_TERMS;
    std::size_t array_sonden = 14;

    if (argc > 1) {
        const std::string first_arg = argv[1];
        if (first_arg == "--help" || first_arg == "-h") {
            std::cout << "Aufruf: ./rydberg_scanner [N] [ziel_faktor] [zeta_terme] [array_sonden]" << std::endl;
            std::cout << "  N            Standard: 1048577" << std::endl;
            std::cout << "  ziel_faktor  Standard: 17" << std::endl;
            std::cout << "  zeta_terme   Standard: " << DEFAULT_ZETA_MOMENTUM_MAX_TERMS << std::endl;
            std::cout << "  array_sonden Standard: 14" << std::endl;
            return 0;
        }
    }

    try {
        if (argc > 1) {
            rsa_n = argv[1];
        }
        if (argc > 2) {
            target_factor = std::stoi(argv[2]);
        }
        if (argc > 3) {
            zeta_terms = static_cast<std::size_t>(std::stoull(argv[3]));
        }
        if (argc > 4) {
            array_sonden = static_cast<std::size_t>(std::stoull(argv[4]));
        }
    } catch (const std::exception& ex) {
        std::cerr << "Ungueltige Argumente: " << ex.what() << std::endl;
        std::cerr << "Aufruf: ./rydberg_scanner [N] [ziel_faktor] [zeta_terme] [array_sonden]" << std::endl;
        return 1;
    }

    set_zeta_momentum_max_terms(zeta_terms);

    std::cout << "[CLI] N=" << rsa_n
              << " | Ziel-Faktor=" << target_factor
              << " | Zeta-Terme=" << get_zeta_momentum_max_terms()
              << " | Array-Sonden=" << array_sonden
              << std::endl;

    // Beispiel: Eine 1024-Bit RSA Zahl (stark verkürzt für Demo)
    HofbiScanner scanner(rsa_n, 8192);
    
    // Wir scannen das ptolemäische Feld um die Rydberg-Grenze
    HofbiScanner::ScanCandidate best_candidate =
        scanner.scan_rydberg_boundary(0.1, 0.9, 0.001);

    const double half_window = 0.01;
    const double high_res_start = std::max(0.1, best_candidate.offset - half_window);
    const double high_res_end = std::min(0.9, best_candidate.offset + half_window);
    scanner.scan_high_res_boundary(high_res_start, high_res_end, 0.000001);

    std::cout << "\nFolge-Run: Hyper-Gleit-Modus um den Randbefund" << std::endl;
    scanner.hofbi_singular_zoom(0.900001000000, 1e-8, 1e-12, 256);
    scanner.hofbi_singular_zoom(0.900001000000, 1e-8, 3e-13, 256);
    scanner.hofbi_singular_zoom(0.900001000000, 1e-8, 1e-13, 256);
    scanner.scan_gradient_sign_changes(0.900001000000, 1e-8, 1e-10, 1e-12);
    scanner.scan_gradient_sign_changes(0.900001000000, 1e-8, 1e-10, 3e-13);
    scanner.rydberg_tile_scan(0.900001000000, 1e-5, 1e-7, 1e-12, 8);
    scanner.wideband_aperture_scan(1e-6, 0.5, 1.1, 1e-12);
    scanner.target_factor_resonance(target_factor, 1e-10, 1e-14);
    if (target_factor == 17) {
        scanner.auto_lock_factor_resonance(1.0e-5, 1.0e-10, 1.0e-15, 400, 5.0e-6);
    } else {
        std::cout << "[AUTO-LOCK] Ueberspringe legacy-17-Auto-Lock fuer Ziel-Faktor "
                  << target_factor << std::endl;
    }
    scanner.run_multi_target_array(array_sonden, 1.0e-10, 1.0e-15, 5000, 5.0e-6, 5);

    return 0;
}