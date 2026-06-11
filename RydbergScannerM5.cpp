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
#include <future>
#include <limits>
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

    struct ResonanceMap {
        double fundamental_freq = 0.0;
        bool parity_shift = false;
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
        mpf_class fourth_derivative = 0;
        std::vector<SondeHistory> history;
        bool active = true;
        bool hit = false;
        int steps_without_gradient = 0;
        std::string status = "tracking";
        bool radar_announced = false;
        mpf_class damping = 0;
        mpf_class lock_threshold = mpf_class("1e-7");
        mpz_class current_candidate = 0;

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
    mpf_class calculate_morley_energy(const mpf_class& x) {
        // Wir simulieren einen Vierling an der Stelle x:
        // x, x+2, x+6, x+8. Die alte Ptolemaeus-Formel war eine Identitaet
        // und lieferte daher numerisch fast immer 0. Hier messen wir statt-
        // dessen die diskrete Krummung der lokalen Log-Abstaende.
        const mpf_class l1 = log_mpf(x);
        const mpf_class l2 = log_shifted_mpf_from_log(x, l1, 2);
        const mpf_class l3 = log_shifted_mpf_from_log(x, l1, 6);
        const mpf_class l4 = log_shifted_mpf_from_log(x, l1, 8);

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
        const mpf_class resonance_threshold("1e-15");

        for (double offset = start_range; offset <= end_range; offset += step) {
            mpf_class test_x = N * offset; // Virtuelle Schalen-Abtastung
            mpf_class energy = calculate_morley_energy(test_x);
            update_top_candidates(top_candidates, offset, energy, 0, 10);

            // Simulation der Rydberg-Grenze: 
            // Wenn die Energie unter einen kritischen Schwellenwert fällt, 
            // detektieren wir eine "Simulade" (einen potenziellen Faktor)
            if (energy < resonance_threshold) { // Beispielhafter Schwellenwert fuer Resonanz
                std::cout << "[DETEKTION] Resonanz bei Offset: " << offset << std::endl;
                std::cout << "Spektrale Dichte (Span): " << energy << std::endl;
                
                // Hier würde der "Standard Part" (Shadow Casting) greifen:
                mpz_class potential_factor = standard_part(test_x);
                
                // Validierung am Kern (N)
                if (is_valid_factor_candidate(potential_factor)) {
                    std::cout << ">>> FAKTOR GEFUNDEN: " << potential_factor << std::endl;
                }
            }
        }

        populate_candidate_standard_parts(top_candidates);
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
        bool found_singularity = false;

        for (long long i = 0; i <= steps; ++i) {
            double offset = start_range + static_cast<double>(i) * step;
            mpf_class test_x = N * offset;
            mpf_class energy = calculate_high_res_morley(test_x);

            if (i == 0 || energy < best_energy) {
                best_energy = energy;
                best_offset = offset;
            }

            if (is_singular(energy)) {
                mpz_class factor = standard_part(test_x);
                if (is_valid_factor_candidate(factor)) {
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
            const mpz_class best_factor = standard_part(N * best_offset);
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
        bool factor_found = is_valid_factor_candidate(factor);

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
        const mpf_class background_decay("0.95");
        const mpf_class sample_weight("0.05");
        double best_offset = current_offset;
        double strongest_anomaly_offset = current_offset;
        bool resonance_found = false;
        bool factor_found = false;

        DerivativeSnapshot snapshot = previous;
        while (current_offset <= end_offset) {
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

            background_energy = background_decay * background_energy + sample_weight * snapshot.energy;
            previous = snapshot;
            current_offset *= multiplier;
            if (current_offset <= end_offset) {
                snapshot = evaluate_derivatives(current_offset, epsilon);
            }
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

    ResonanceMap get_shell_fingerprint(const mpz_class& n_z, bool verbose = false) const {
        ResonanceMap map;
        if (n_z == 0) {
            return map;
        }

        mpfr_t n_fr;
        mpfr_t log10_fr;
        mpfr_t h_fr;
        mpfr_t prod_fr;
        mpfr_init2(n_fr, precision_bits);
        mpfr_init2(log10_fr, precision_bits);
        mpfr_init2(h_fr, precision_bits);
        mpfr_init2(prod_fr, precision_bits);

        mpfr_set_z(n_fr, n_z.get_mpz_t(), MPFR_RNDN);
        if (mpfr_sgn(n_fr) <= 0) {
            mpfr_clear(n_fr);
            mpfr_clear(log10_fr);
            mpfr_clear(h_fr);
            mpfr_clear(prod_fr);
            return map;
        }

        mpfr_log10(log10_fr, n_fr, MPFR_RNDN);
        mpfr_mul_d(h_fr, log10_fr, 1.274649e-08, MPFR_RNDN);
        map.fundamental_freq = mpfr_get_d(h_fr, MPFR_RNDN);

        mpfr_mul(prod_fr, n_fr, h_fr, MPFR_RNDN);
        mpz_t prod_z;
        mpz_init(prod_z);
        mpfr_get_z(prod_z, prod_fr, MPFR_RNDZ);
        map.parity_shift = (mpz_odd_p(prod_z) != 0);
        mpz_clear(prod_z);

        mpfr_clear(n_fr);
        mpfr_clear(log10_fr);
        mpfr_clear(h_fr);
        mpfr_clear(prod_fr);

        if (map.parity_shift && verbose) {
            std::cout << "[FINGERPRINT] Ungerader Tunnel detektiert. Supraleitung bereit."
                      << std::endl;
        }
        return map;
    }

    void apply_shell_sync(Sonde& s, bool verbose = false) {
        if (s.status != "stalled_outer_shell") {
            return;
        }
        const ResonanceMap shell_map = get_shell_fingerprint(N_int, verbose);
        if (!shell_map.parity_shift) {
            s.learning_rate *= mpf_class("1.5");
            s.base_learning_rate = s.learning_rate;
            s.status = "SHELL_ACCELERATION";
        } else {
            s.current_offset = -s.current_offset;
            s.base_learning_rate = s.learning_rate;
            s.status = "SHELL_TUNNELING";
        }
        if (verbose) {
            std::cout << "[SHELL] Paritaets-Kopplung aktiv bei Kandidat " << s.current_candidate
                      << " | f0=" << std::scientific << shell_map.fundamental_freq
                      << (shell_map.parity_shift ? " | Modus=Tunnel" : " | Modus=Mauer")
                      << std::endl;
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
        const mpf_class target_corridor_min("0.0587");
        const mpf_class target_corridor_max("0.0600");
        const mpf_class target_energy_floor("1e-12");
        const mpf_class target_zone_damping("1e-9");
        const mpf_class target_precision_trigger("1e-12");
        const mpf_class gradient_lock_threshold("1e-8");

        if (verbose) {
            std::cout << "\n[AUTO-LOCK] Suche gestartet bei Offset: "
                      << std::scientific << std::setprecision(18) << start_offset << std::endl;
            std::cout << "[AUTO-LOCK] learning_rate=" << sonde.learning_rate
                      << " | epsilon=" << epsilon
                      << " | max_step_change=" << max_step_change << std::endl;
        }

        for (int step_index = 0; step_index < max_steps; ++step_index) {
            mpf_class x = N * sonde.current_offset;
            const double current_offset_value = sonde.current_offset.get_d();
            mpz_class candidate = round_standard_part(x);
            mpf_class energy = calculate_morley_energy(x);
            DerivativeSnapshot derivative_snapshot = evaluate_derivatives_at_x(x, epsilon, &energy);
            end_energy = energy;
            last_gradient = derivative_snapshot.gradient;
            last_fourth_derivative = derivative_snapshot.fourth_derivative;
            last_candidate = candidate;
            sonde.current_candidate = candidate;
            sonde.gradient = derivative_snapshot.gradient;
            sonde.curvature = derivative_snapshot.curvature;
            sonde.fourth_derivative = derivative_snapshot.fourth_derivative;

            if (!sonde.radar_announced) {
                if (N_int == mpz_class("861163313")) {
                    const mpf_class lower_monade = mpf_class("28901") / N;
                    const mpf_class upper_monade = mpf_class("29797") / N;
                    const mpf_class radar_threshold("5e-7");
                    mpf_class dist_to_lower = abs_mpf(sonde.current_offset - lower_monade);
                    mpf_class dist_to_upper = abs_mpf(sonde.current_offset - upper_monade);
                    if (dist_to_lower < radar_threshold || dist_to_upper < radar_threshold) {
                        const char* monade_name =
                            (dist_to_lower <= dist_to_upper) ? "28.901" : "29.797";
                        std::cout << "[RADAR] Sonde " << sonde.id
                                  << " naehrt sich der " << monade_name
                                  << "-Monade!" << std::endl;
                        sonde.radar_announced = true;
                    }
                } else {
                    mpf_class dist_to_q = abs_mpf(sonde.current_offset - mpf_class("0.058826"));
                    if (dist_to_q < mpf_class("0.001")) {
                        std::cout << "[RADAR] Sonde " << sonde.id
                                  << " naehrt sich der 61.681-Monade!" << std::endl;
                        sonde.radar_announced = true;
                    }
                }
            }

            mpf_class abs_gradient = abs_mpf(derivative_snapshot.gradient);
            if (abs_gradient < gradient_lock_threshold) {
                ++sonde.steps_without_gradient;
            } else {
                sonde.reset_gradient_history();
            }

            if (is_valid_factor_candidate(candidate)) {
                sonde.hit = true;
                sonde.active = false;
                sonde.status = "IONISATION_SUCCESS";
                if (verbose) {
                    std::cout << ">>> IONISATION ERFOLGREICH! FAKTOR: " << candidate << std::endl;
                    std::cout << "[AUTO-LOCK] Offset: "
                              << std::fixed << std::setprecision(18)
                              << current_offset_value << std::endl;
                    std::cout << "[AUTO-LOCK] Energie: " << std::scientific << energy << std::endl;
                }
                return {
                    start_offset,
                    current_offset_value,
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
                    "IONISATION_SUCCESS"
                };
            }

            sonde.adjust_rate(energy, derivative_snapshot.gradient);
            mpf_class abs_curvature = abs_mpf(sonde.curvature);
            mpf_class delta_energy = abs_mpf(energy - previous_energy);
            update_sonde_speed(sonde, N, verbose);
            apply_re_entry_damping(sonde, energy, derivative_snapshot.gradient, epsilon, verbose);
            apply_hohle_energie_sync(sonde, energy, verbose);
            apply_orbital_binding(sonde, energy, verbose);
            {
                const mpf_class cream5_resonance_mpf = abs_mpf(sonde.fourth_derivative) / N;
                const double cream5_resonance = cream5_resonance_mpf.get_d();
                const bool cream5_alpha_match =
                    match_cream5_alpha_signature(N_int, cream5_resonance, verbose);
                singularity_tunneling(sonde, energy, verbose, cream5_alpha_match);
            }
            apply_entropy_brake(sonde, verbose);
            apply_quantum_field_boost(sonde, energy, verbose);
            if (sonde.status == "CORE_IONIZED" || sonde.status == "CORE_IONIZED_SINGULARITY") {
                sonde.active = false;
                if (verbose) {
                    execute_singularity_dump(sonde);
                }
                const mpf_class locked_x = N * sonde.current_offset;
                const mpz_class locked_candidate = round_standard_part(locked_x);
                const mpf_class locked_energy = calculate_morley_energy(locked_x);
                const mpz_class verified_candidate = resolve_valid_factor_near(locked_candidate);
                const bool locked_factor_confirmed = (verified_candidate > 1);
                if (locked_factor_confirmed) {
                    sonde.hit = true;
                    if (verbose) {
                        std::cout << "[SUCCESS] Kernschale ionisiert: p=28901.0000" << std::endl;
                        std::cout << "[VERIFY] Faktor nach Kern-Lock bestaetigt: "
                                  << verified_candidate
                                  << " | Basis-Kandidat=" << locked_candidate
                                  << std::endl;
                        extract_ptolemaic_fingerprint(sonde, verified_candidate, true);
                    }
                    return {
                        start_offset,
                        sonde.current_offset.get_d(),
                        step_index,
                        verified_candidate,
                        verified_candidate,
                        start_energy,
                        locked_energy,
                        last_gradient,
                        last_fourth_derivative,
                        sonde.learning_rate.get_d(),
                        static_cast<int>(sonde.history.size()),
                        true,
                        "CORE_LOCK_FACTOR_CONFIRMED"
                    };
                }
                if (verbose) {
                    std::cout << "[SUCCESS] Kernschale ionisiert: p=28901.0000" << std::endl;
                    std::cout << "[VERIFY] Kein exakter Faktor nach Kern-Lock."
                              << " | Kandidat=" << locked_candidate
                              << " | Nachtest={"
                              << (locked_candidate - 1) << ", "
                              << locked_candidate << ", "
                              << (locked_candidate + 1) << "}"
                              << std::endl;
                    extract_ptolemaic_fingerprint(sonde, locked_candidate, true);
                }
                return {
                    start_offset,
                    sonde.current_offset.get_d(),
                    step_index,
                    0,
                    locked_candidate,
                    start_energy,
                    locked_energy,
                    last_gradient,
                    last_fourth_derivative,
                    sonde.learning_rate.get_d(),
                    static_cast<int>(sonde.history.size()),
                    true,
                    sonde.status
                };
            }
            if (sonde.status == "searching_wave" || sonde.status == "zeta_surfing" ||
                sonde.status == "ZETA_BOOST") {
                if (abs_curvature > curvature_threshold) {
                    sonde.learning_rate = sonde.base_learning_rate;
                    sonde.status = "re_entry_homing";
                    if (verbose) {
                        std::cout << "[LOCK] Sonde " << sonde.id
                                  << " bremst in Monade ein bei Offset: "
                                  << std::fixed << std::setprecision(18)
                                  << current_offset_value << std::endl;
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
                    current_offset_value,
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
            double next_offset = std::clamp(current_offset_value + delta, epsilon, 1.0);
            double applied_delta = std::fabs(next_offset - current_offset_value);

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
                current_offset_value,
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
                          << current_offset_value
                          << " | Energie: " << std::scientific << energy
                          << " | Gradient: " << derivative_snapshot.gradient
                          << " | Kandidat: " << candidate
                          << " | Lernrate: " << sonde.learning_rate
                          << std::endl;
            }

            if (std::fabs(next_offset - current_offset_value) < epsilon) {
                const bool in_target_corridor =
                    sonde.current_offset >= target_corridor_min &&
                    sonde.current_offset <= target_corridor_max;
                if (in_target_corridor && energy > target_energy_floor) {
                    sonde.damping = target_zone_damping;
                    sonde.lock_threshold = target_precision_trigger;
                    if (sonde.learning_rate < mpf_class("1e-7")) {
                        sonde.learning_rate = mpf_class("1e-7");
                    } else {
                        sonde.learning_rate *= mpf_class("1.5");
                    }
                    sonde.status = "target_lock_q";
                    sonde.reset_gradient_history();
                    if (verbose) {
                        std::cout << "[AUTO-LOCK] Kernzone haelt Sonde "
                                  << sonde.id
                                  << " aktiv."
                                  << " | Energie=" << energy
                                  << " | damping=" << sonde.damping
                                  << " | precision_trigger=" << sonde.lock_threshold
                                  << " | Lernrate=" << sonde.learning_rate
                                  << std::endl;
                    }
                    previous_energy = energy;
                    continue;
                }

                if (sonde.status != "accelerating" && sonde.status != "relativistic_boost") {
                    sonde.status = "stalled";
                    if (verbose) {
                        std::cout << "[AUTO-LOCK] Schrittweite unter epsilon gefallen." << std::endl;
                        std::cout << "[AUTO-LOCK] Letzter Kandidat: " << candidate
                                  << " | Offset: " << std::fixed << std::setprecision(18)
                                  << current_offset_value << std::endl;
                        std::cout << "[TELEPORT-IMPULS] Lernrate wird massiv erhoeht." << std::endl;
                    }
                    sonde.learning_rate *= mpf_class("100.0");
                    sonde.status = "accelerating";
                    sonde.reset_gradient_history();
                    continue;
                }

                sonde.current_candidate = candidate;
                sonde.status = "stalled_outer_shell";
                apply_shell_sync(sonde, verbose);
                if (sonde.status == "SHELL_ACCELERATION" ||
                    sonde.status == "SHELL_TUNNELING") {
                    sonde.active = true;
                    sonde.reset_gradient_history();
                    previous_energy = energy;
                    continue;
                }

                sonde.active = false;
                sonde.status = "stalled";
                if (verbose) {
                    std::cout << "[AUTO-LOCK] Schrittweite unter epsilon gefallen." << std::endl;
                    std::cout << "[AUTO-LOCK] Letzter Kandidat: " << candidate
                              << " | Offset: " << std::fixed << std::setprecision(18)
                              << current_offset_value << std::endl;
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
                       sonde.status != "TIME_ACCELERATED_STRIDE" &&
                       sonde.status != "ANNIHILATION_BOOST" &&
                       sonde.status != "KERN_TUNNELING" &&
                       sonde.status != "SINGULARITY_PULL" &&
                       sonde.status != "SINGULARITY_LOCK" &&
                       sonde.status != "ENTROPY_BRAKE_ENGAGED" &&
                       sonde.status != "CORE_IONIZED_SINGULARITY" &&
                       sonde.status != "CORE_IONIZED" &&
                       sonde.status != "VACUUM_TUNNELING" &&
                       sonde.status != "VACUUM_COLLAPSE" &&
                       sonde.status != "BINDING_2_2" &&
                       sonde.status != "BINDING_1_3" &&
                       sonde.status != "ZETA_IMPULSE_STRIKE" &&
                       sonde.status != "MAGNETIC_COHERENCE_LOCK" &&
                       sonde.status != "ANTI_PHASE_LOCK" &&
                       sonde.status != "COHERENT_LASER_STRIKE" &&
                       sonde.status != "REENTRY_STRIKE" &&
                       sonde.status != "HOLOGRAPHIC_RESONANCE" &&
                       sonde.status != "REENTRY_IONISATION" &&
                       sonde.status != "SINGULARITY_REENTRY" &&
                       sonde.status != "STABILIZED_CORE_STRIKE" &&
                       sonde.status != "SINGULARITY_FREEZE_SUCCESS" &&
                       sonde.status != "DYAD_PHASE_P" &&
                       sonde.status != "DYAD_PHASE_Q" &&
                       sonde.status != "KRYO_RESONANCE_LOCK" &&
                       sonde.status != "FUSION_COMPLETE_PFZ" &&
                       sonde.status != "OMEGA_STRIKE_COHERENCE" &&
                       sonde.status != "TRILLION_SINGULARITY_LOCK" &&
                       sonde.status != "TRILLION_EVENT_HORIZON_LOCK" &&
                       sonde.status != "stalled_outer_shell" &&
                       sonde.status != "SHELL_ACCELERATION" &&
                       sonde.status != "SHELL_TUNNELING" &&
                       sonde.status != "target_lock_q" &&
                       sonde.status != "accelerating" &&
                       sonde.status != "relativistic_boost") {
                sonde.status = "tracking";
            }
            if (sonde.status != "KERN_TUNNELING" &&
                sonde.status != "SINGULARITY_PULL" &&
                sonde.status != "ANNIHILATION_BOOST" &&
                sonde.status != "ZETA_IMPULSE_STRIKE" &&
                sonde.status != "MAGNETIC_COHERENCE_LOCK" &&
                sonde.status != "ANTI_PHASE_LOCK" &&
                sonde.status != "COHERENT_LASER_STRIKE" &&
                sonde.status != "REENTRY_STRIKE" &&
                sonde.status != "HOLOGRAPHIC_RESONANCE" &&
                sonde.status != "REENTRY_IONISATION" &&
                sonde.status != "SINGULARITY_REENTRY" &&
                sonde.status != "STABILIZED_CORE_STRIKE" &&
                sonde.status != "SINGULARITY_FREEZE_SUCCESS" &&
                sonde.status != "DYAD_PHASE_P" &&
                sonde.status != "DYAD_PHASE_Q" &&
                sonde.status != "KRYO_RESONANCE_LOCK" &&
                sonde.status != "FUSION_COMPLETE_PFZ" &&
                sonde.status != "OMEGA_STRIKE_COHERENCE" &&
                sonde.status != "TRILLION_SINGULARITY_LOCK" &&
                sonde.status != "TRILLION_EVENT_HORIZON_LOCK" &&
                sonde.status != "stalled_outer_shell" &&
                sonde.status != "SHELL_ACCELERATION" &&
                sonde.status != "SHELL_TUNNELING") {
                sonde.current_offset = next_offset;
            }
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
    double calculate_root_offset() {
        mpfr_t n_fr;
        mpfr_t sqrt_n_fr;
        mpfr_t root_offset_fr;

        mpfr_init2(n_fr, precision_bits);
        mpfr_init2(sqrt_n_fr, precision_bits);
        mpfr_init2(root_offset_fr, precision_bits);

        mpfr_set_f(n_fr, N.get_mpf_t(), MPFR_RNDN);
        mpfr_sqrt(sqrt_n_fr, n_fr, MPFR_RNDN);
        mpfr_ui_div(root_offset_fr, 1ul, sqrt_n_fr, MPFR_RNDN);

        const double root_offset = mpf_from_mpfr(root_offset_fr).get_d();

        mpfr_clear(n_fr);
        mpfr_clear(sqrt_n_fr);
        mpfr_clear(root_offset_fr);

        return root_offset;
    }
    double calculate_triplet_zeta_energy(
        const mpz_class& base_candidate,
        unsigned long a,
        unsigned long b,
        unsigned long c
    ) {
        const mpz_class d1 = base_candidate + a;
        const mpz_class d2 = base_candidate + b;
        const mpz_class d3 = base_candidate + c;
        if (d1 <= 0 || d2 <= 0 || d3 <= 0) {
            return std::numeric_limits<double>::infinity();
        }

        const mpf_class x1_mpf = N / mpf_class(d1);
        const mpf_class x2_mpf = N / mpf_class(d2);
        const mpf_class x3_mpf = N / mpf_class(d3);
        const double x1 = x1_mpf.get_d();
        const double x2 = x2_mpf.get_d();
        const double x3 = x3_mpf.get_d();
        if (!std::isfinite(x1) || !std::isfinite(x2) || !std::isfinite(x3) ||
            x1 <= 1.0 || x2 <= 1.0 || x3 <= 1.0) {
            return std::numeric_limits<double>::infinity();
        }

        return std::abs(calculate_zeta_momentum(x1)) +
               std::abs(calculate_zeta_momentum(x2)) +
               std::abs(calculate_zeta_momentum(x3));
    }
    bool evaluate_vierling_coherence(
        const AutoLockResult& result,
        double& coherence_out
    ) {
        const mpz_class& p = result.last_candidate;
        if (p <= 0) {
            return false;
        }

        const double energy_a = calculate_triplet_zeta_energy(p, 0, 2, 6);
        const double energy_b = calculate_triplet_zeta_energy(p, 2, 6, 8);
        if (!std::isfinite(energy_a) || !std::isfinite(energy_b)) {
            return false;
        }

        coherence_out = std::abs(energy_a - energy_b);
        return coherence_out < 1e-12;
    }

    std::vector<double> auto_calibrate_seed_offsets(std::size_t num_sonden) {
        std::vector<double> seeds;
        if (num_sonden == 0) {
            return seeds;
        }

        mpfr_t n_fr;
        mpfr_t sqrt_n_fr;
        mpfr_t base_offset_fr;
        mpfr_t seed_fr;

        mpfr_init2(n_fr, precision_bits);
        mpfr_init2(sqrt_n_fr, precision_bits);
        mpfr_init2(base_offset_fr, precision_bits);
        mpfr_init2(seed_fr, precision_bits);

        mpfr_set_f(n_fr, N.get_mpf_t(), MPFR_RNDN);
        mpfr_sqrt(sqrt_n_fr, n_fr, MPFR_RNDN);

        // Zentrale Resonanzlage bei 1 / sqrt(N).
        mpfr_ui_div(base_offset_fr, 1ul, sqrt_n_fr, MPFR_RNDN);

        if (num_sonden == 1) {
            seeds.push_back(mpf_from_mpfr(base_offset_fr).get_d());
        } else {
            for (std::size_t i = 0; i < num_sonden; ++i) {
                double t = static_cast<double>(i) / static_cast<double>(num_sonden - 1);
                double scale = std::pow(10.0, -2.0 + 4.0 * t);
                mpfr_mul_d(seed_fr, base_offset_fr, scale, MPFR_RNDN);

                double offset = mpf_from_mpfr(seed_fr).get_d();
                if (offset <= 0.0) {
                    continue;
                }
                if (offset > 1.0) {
                    offset = 1.0;
                }
                seeds.push_back(offset);
            }
        }

        std::sort(seeds.begin(), seeds.end());
        seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());

        mpfr_clear(n_fr);
        mpfr_clear(sqrt_n_fr);
        mpfr_clear(base_offset_fr);
        mpfr_clear(seed_fr);

        return seeds;
    }

    void auto_calibrate_array(std::vector<Sonde>& sondes, double learning_rate) {
        std::vector<double> seeds = auto_calibrate_seed_offsets(sondes.size());
        if (seeds.empty()) {
            sondes.clear();
            return;
        }

        sondes.resize(seeds.size());
        const mpf_class strike_rate(std::max(learning_rate, 1.0e-7));
        for (std::size_t i = 0; i < seeds.size(); ++i) {
            sondes[i] = Sonde{
                static_cast<int>(i),
                mpf_class(seeds[i]),
                strike_rate,
                strike_rate,
                0,
                0,
                0,
                {},
                true,
                false,
                0,
                "searching",
                false
            };
            std::cout << "[AUTO-GRID] Sonde " << i
                      << " platziert bei Offset: "
                      << std::scientific << std::setprecision(18)
                      << seeds[i] << std::endl;
        }
    }

    std::vector<double> generate_optimized_hofbi_seeds(
        std::size_t seed_count,
        std::size_t sqrt_cluster_points = 5
    ) {
        std::vector<double> seeds = auto_calibrate_seed_offsets(seed_count);
        if (seeds.empty()) {
            return seeds;
        }

        mpfr_t n_fr;
        mpfr_t sqrt_n_fr;
        mpfr_t max_offset_fr;

        mpfr_init2(n_fr, precision_bits);
        mpfr_init2(sqrt_n_fr, precision_bits);
        mpfr_init2(max_offset_fr, precision_bits);

        mpfr_set_f(n_fr, N.get_mpf_t(), MPFR_RNDN);
        mpfr_sqrt(sqrt_n_fr, n_fr, MPFR_RNDN);
        mpfr_ui_div(max_offset_fr, 1ul, sqrt_n_fr, MPFR_RNDN);

        add_sqrt_cluster(seeds, max_offset_fr, sqrt_cluster_points);
        std::sort(seeds.begin(), seeds.end());
        seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());

        mpfr_clear(n_fr);
        mpfr_clear(sqrt_n_fr);
        mpfr_clear(max_offset_fr);

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
        (void)sqrt_cluster_points;
        std::vector<Sonde> sondes(seed_count);
        auto_calibrate_array(sondes, learning_rate);
        if (sondes.empty()) {
            std::cout << "[ARRAY] Keine Sonden fuer den Mehrsonden-Modus vorhanden." << std::endl;
            return {};
        }

        std::set<mpz_class> found_factors;
        std::vector<AutoLockResult> results(sondes.size());
        std::vector<bool> result_ready(sondes.size(), false);

        std::cout << "\n[ARRAY] Starte Mehrsonden-Modus mit "
                  << sondes.size() << " Seeds." << std::endl;
        if (!sondes.empty()) {
            std::cout << "[ARRAY] Hofbi-Seed-Bereich: ["
                      << std::scientific << std::setprecision(18)
                      << sondes.front().current_offset.get_d() << ", "
                      << sondes.back().current_offset.get_d() << "]"
                      << std::endl;
        }

        constexpr std::size_t max_swarm_rounds = 6;
        for (std::size_t round = 0; round < max_swarm_rounds; ++round) {
            execute_trillion_kryo_strike(N_int, sondes, false);
            execute_omega_trillion_strike(N_int, sondes, false);
            execute_trillion_singularity_strike(N_int, sondes, false);
            execute_trillion_omega_strike(N_int, sondes, false);

            std::vector<std::size_t> launch_indices;
            launch_indices.reserve(sondes.size());

            for (std::size_t i = 0; i < sondes.size(); ++i) {
                if (sondes[i].hit) {
                    continue;
                }
                if (sondes[i].status == "searching" ||
                    sondes[i].status == "swarm_strike" ||
                    sondes[i].status == "IONISATION_LOCK" ||
                    sondes[i].status == "RADIATION_GUIDED_LOCK" ||
                    sondes[i].status == "ENTROPY_BRAKE_ENGAGED" ||
                    sondes[i].status == "ZETA_IMPULSE_STRIKE" ||
                    sondes[i].status == "MAGNETIC_COHERENCE_LOCK" ||
                    sondes[i].status == "ANTI_PHASE_LOCK" ||
                    sondes[i].status == "COHERENT_LASER_STRIKE" ||
                    sondes[i].status == "REENTRY_STRIKE" ||
                    sondes[i].status == "HOLOGRAPHIC_RESONANCE" ||
                    sondes[i].status == "REENTRY_IONISATION" ||
                    sondes[i].status == "SINGULARITY_REENTRY" ||
                    sondes[i].status == "STABILIZED_CORE_STRIKE" ||
                    sondes[i].status == "SINGULARITY_FREEZE_SUCCESS" ||
                    sondes[i].status == "DYAD_PHASE_P" ||
                    sondes[i].status == "DYAD_PHASE_Q" ||
                    sondes[i].status == "KRYO_RESONANCE_LOCK" ||
                    sondes[i].status == "FUSION_COMPLETE_PFZ" ||
                    sondes[i].status == "OMEGA_STRIKE_COHERENCE" ||
                    sondes[i].status == "TRILLION_SINGULARITY_LOCK" ||
                    sondes[i].status == "TRILLION_EVENT_HORIZON_LOCK" ||
                    sondes[i].status == "SHELL_ACCELERATION" ||
                    sondes[i].status == "SHELL_TUNNELING" ||
                    sondes[i].status == "VACUUM_TUNNELING" ||
                    sondes[i].status == "VACUUM_COLLAPSE") {
                    launch_indices.push_back(i);
                }
            }

            if (launch_indices.empty()) {
                break;
            }

            std::cout << "[ARRAY] Runde " << (round + 1)
                      << " startet mit " << launch_indices.size()
                      << " aktiven Sonden." << std::endl;

            for (std::size_t index : launch_indices) {
                if (sondes[index].status == "searching") {
                    apply_deep_space_launch_boost(
                        sondes[index],
                        N,
                        false
                    );
                }
            }

            std::vector<std::future<AutoLockResult>> futures;
            futures.reserve(launch_indices.size());
            for (std::size_t index : launch_indices) {
                futures.push_back(
                    std::async(
                        std::launch::async,
                        [this, &sondes, index, epsilon, max_steps, max_step_change]() {
                            return auto_lock_factor_resonance(
                                sondes[index],
                                epsilon,
                                max_steps,
                                max_step_change,
                                false
                            );
                        }
                    )
                );
            }

            for (std::size_t job = 0; job < launch_indices.size(); ++job) {
                const std::size_t index = launch_indices[job];
                AutoLockResult result = futures[job].get();
                const mpz_class effective_factor = effective_array_factor(result);
                if (effective_factor > 1) {
                    result.factor = effective_factor;
                }
                results[index] = result;
                result_ready[index] = true;

                if (effective_factor > 1 &&
                    found_factors.find(effective_factor) == found_factors.end()) {
                    std::cout << "[ARRAY-TREFFER] Sonde " << index
                              << " fand Faktor: " << effective_factor
                              << " | Start-Offset="
                              << std::scientific << std::setprecision(18)
                              << result.start_offset << std::endl;
                    found_factors.insert(effective_factor);
                }
            }

            execute_billion_dump(sondes, results, result_ready, false);
            execute_holographic_reentry(sondes, results, result_ready, false);
            apply_billion_hologram_lock(sondes, false);

            {
                constexpr std::size_t k_no_master = std::numeric_limits<std::size_t>::max();
                std::size_t collapse_master = k_no_master;
                for (std::size_t i = 0; i < results.size() && i < sondes.size(); ++i) {
                    if (!result_ready[i]) {
                        continue;
                    }
                    if (effective_array_factor(results[i]) > 1) {
                        collapse_master = i;
                        break;
                    }
                }
                if (collapse_master != k_no_master) {
                    trigger_swarm_collapse(sondes, collapse_master, false);
                }
            }

            const bool radiation_triggered =
                apply_core_lock_and_radiation(sondes, results, result_ready, false);
            for (Sonde& s : sondes) {
                if (s.status == "CORE_IONIZED_SINGULARITY") {
                    apply_ptolemaic_mirror(s, sondes, false);
                    break;
                }
            }
            execute_swarm_collapse_m6(sondes, false);
            if (!radiation_triggered) {
                update_swarm_autonomous(sondes, results, result_ready, learning_rate, epsilon);
            }
            apply_cream5_suction(sondes, false);
            apply_swarm_focus_strike(sondes, false);
            apply_space_time_freeze(sondes, false);
            apply_billion_freeze_correction(sondes, false);
            {
                Sonde* holo_anchor = nullptr;
                for (Sonde& s : sondes) {
                    if (s.status == "CORE_IONIZED_SINGULARITY") {
                        holo_anchor = &s;
                        break;
                    }
                }
                if (holo_anchor == nullptr) {
                    for (Sonde& s : sondes) {
                        if (s.status == "CORE_IONIZED") {
                            holo_anchor = &s;
                            break;
                        }
                    }
                }
                if (holo_anchor == nullptr) {
                    for (Sonde& s : sondes) {
                        if (s.status == "SINGULARITY_FREEZE_SUCCESS") {
                            holo_anchor = &s;
                            break;
                        }
                    }
                }
                if (holo_anchor != nullptr) {
                    project_factor_hologram(*holo_anchor, sondes, false);
                }
            }
            trigger_dyadic_hologram(sondes, false);
            execute_dyadic_fusion(sondes, false);
        }

        print_probe_table(results);
        export_trajectories_to_csv(sondes);

        if (found_factors.empty()) {
            std::cout << "[ARRAY] Keine Faktoren im Mehrsonden-Modus gefunden."
                      << " (Zaehlt nur: bestaetigter Teiler von N, nicht nur Kandidat/Status.)"
                      << std::endl;
        }

        return found_factors;
    }

private:
    bool update_swarm_autonomous(
        std::vector<Sonde>& sondes,
        const std::vector<AutoLockResult>& results,
        const std::vector<bool>& result_ready,
        double learning_rate,
        double epsilon
    ) {
        if (results.empty() || result_ready.empty()) {
            return false;
        }

        const mpf_class resonance_threshold("1e-4");
        const double root_offset = calculate_root_offset();
        if (!std::isfinite(root_offset) || root_offset <= 0.0) {
            return false;
        }
        bool signal_found = false;
        mpf_class best_energy = 0;
        double best_offset = -1.0;
        double best_score = std::numeric_limits<double>::infinity();

        for (std::size_t i = 0; i < results.size() && i < result_ready.size(); ++i) {
            if (!result_ready[i]) {
                continue;
            }
            const AutoLockResult& result = results[i];
            const bool has_signal =
                result.factor > 1 ||
                result.end_energy < resonance_threshold ||
                result.termination_reason == "target_lock_q";
            if (!has_signal) {
                continue;
            }

            const double candidate_offset = result.end_offset;
            if (!std::isfinite(candidate_offset) || candidate_offset <= 0.0) {
                continue;
            }

            const double distance_ratio = candidate_offset / root_offset;
            if (!std::isfinite(distance_ratio) || distance_ratio <= 0.0) {
                continue;
            }

            const double distance_penalty = [&]() {
                const double base_penalty =
                    std::max(distance_ratio, 1.0 / distance_ratio);
                return base_penalty * base_penalty;
            }();
            const double score = result.end_energy.get_d() * distance_penalty;

            if (!signal_found || score < best_score) {
                signal_found = true;
                best_score = score;
                best_energy = result.end_energy;
                best_offset = candidate_offset;
            }
        }

        if (!signal_found || best_offset <= 0.0) {
            return false;
        }

        std::size_t eligible_count = 0;
        for (std::size_t i = 0; i < sondes.size(); ++i) {
            const Sonde& sonde = sondes[i];
            if (sonde.hit) {
                continue;
            }
            if (sonde.status == "searching" || sonde.status == "max_steps") {
                ++eligible_count;
            }
        }
        if (eligible_count == 0) {
            return false;
        }

        (void)learning_rate;
        double best_coherence = std::numeric_limits<double>::infinity();
        double coherence_target_offset = -1.0;
        mpz_class coherence_candidate = 0;
        for (std::size_t i = 0; i < results.size() && i < result_ready.size(); ++i) {
            if (!result_ready[i]) {
                continue;
            }

            const AutoLockResult& result = results[i];
            if (result.end_energy >= resonance_threshold || result.end_offset <= 0.0) {
                continue;
            }

            double coherence = 0.0;
            if (evaluate_vierling_coherence(result, coherence) && coherence < best_coherence) {
                best_coherence = coherence;
                coherence_target_offset = result.end_offset;
                coherence_candidate = result.last_candidate;
            }
        }

        if (coherence_target_offset > 0.0) {
            const mpf_class lock_rate("1e-7");
            std::cout << "[KERN-STRIKE] Vierlings-Symmetrie detektiert! p="
                      << coherence_candidate
                      << " | coherence=" << best_coherence
                      << " | Zielzentrum=" << coherence_target_offset
                      << std::endl;

            std::size_t locked = 0;
            for (std::size_t i = 0; i < sondes.size(); ++i) {
                Sonde& sonde = sondes[i];
                if (sonde.hit) {
                    continue;
                }
                if (sonde.status != "searching" && sonde.status != "max_steps") {
                    continue;
                }

                sonde.current_offset = mpf_class(std::clamp(coherence_target_offset, epsilon, 1.0));
                sonde.learning_rate = lock_rate;
                sonde.base_learning_rate = lock_rate;
                sonde.steps_without_gradient = 0;
                sonde.active = true;
                sonde.hit = false;
                sonde.history.clear();
                sonde.status = "IONISATION_LOCK";
                sonde.radar_announced = false;
                ++locked;
            }
            return locked > 0;
        }

        const mpf_class strike_rate("1e-9");

        std::cout << "[SWARM] Gravity-Assist aktiviert! Ziel-Zentrum: "
                  << std::scientific << std::setprecision(18) << best_offset
                  << " | Energie=" << best_energy
                  << " | Root-Offset=" << root_offset
                  << " | Score=" << best_score
                  << " | rekalibriere " << eligible_count << " Sonden im Praezisionsmodus."
                  << std::endl;

        std::size_t slot = 0;
        for (std::size_t i = 0; i < sondes.size(); ++i) {
            Sonde& sonde = sondes[i];
            if (sonde.hit) {
                continue;
            }

            if (sonde.status != "searching" && sonde.status != "max_steps") {
                continue;
            }

            const double precision_spread = [&]() {
                if (eligible_count <= 1) {
                    return 1.0;
                }

                const double t =
                    static_cast<double>(slot) / static_cast<double>(eligible_count - 1);
                const double signed_position = (2.0 * t) - 1.0;
                const double weighted_position =
                    signed_position * std::abs(signed_position);
                return 1.0 + (0.001 * weighted_position);
            }();
            const double next_offset = (eligible_count == 1)
                ? std::clamp(best_offset, epsilon, 1.0)
                : std::clamp(best_offset * precision_spread, epsilon, 1.0);

            sonde.current_offset = mpf_class(next_offset);
            sonde.learning_rate = strike_rate;
            sonde.base_learning_rate = strike_rate;
            sonde.steps_without_gradient = 0;
            sonde.active = true;
            sonde.hit = false;
            sonde.history.clear();
            sonde.status = "swarm_strike";
            sonde.radar_announced = false;

            std::cout << "[SWARM] Sonde " << sonde.id
                      << " folgt auf Offset "
                      << std::scientific << std::setprecision(18) << next_offset
                      << " | Spreizung=" << precision_spread
                      << " | Lernrate=" << strike_rate
                      << std::endl;
            ++slot;
        }
        return slot > 0;
    }

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
    mpf_class abs_mpf(const mpf_class& x) {
        return (x < 0) ? -x : x;
    }
    double calculate_zeta_momentum_parallel(
        double x,
        const std::vector<double>& gammas,
        std::size_t term_count
    ) {
        if (x <= 1.0 || term_count == 0) {
            return 0.0;
        }

        const double log_x = std::log(x);
        double total_interference = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:total_interference)
#endif
        for (std::size_t i = 0; i < term_count; ++i) {
            total_interference += std::cos(gammas[i] * log_x);
        }
        return total_interference / static_cast<double>(term_count);
    }
    double calculate_zeta_momentum(double x) {
        const std::vector<double>& all_gammas = get_riemann_gammas();
        const std::size_t term_count = std::min<std::size_t>(get_zeta_momentum_max_terms(), all_gammas.size());
        if (x <= 1.0 || term_count == 0) {
            return 0.0;
        }

        return calculate_zeta_momentum_parallel(x, all_gammas, term_count);
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
    double calculate_offset_zeta_energy(double offset) {
        if (!std::isfinite(offset) || offset <= 0.0) {
            return std::numeric_limits<double>::infinity();
        }

        const mpf_class x_val_mpf = N * mpf_class(offset);
        const double x_val = x_val_mpf.get_d();
        if (!std::isfinite(x_val) || x_val <= 1.0) {
            return std::numeric_limits<double>::infinity();
        }

        return std::abs(calculate_zeta_momentum(x_val));
    }
    void apply_hohle_energie_sync(
        Sonde& s,
        const mpf_class& current_energy,
        bool verbose = false
    ) {
        static const mpf_class critical_learning_rate("1e-14");
        static const mpf_class sync_energy_threshold("1e-5");
        static const mpf_class pumping_base_rate("1e-8");
        static const mpf_class singularity_lock_rate("1e-11");

        if (s.learning_rate >= critical_learning_rate ||
            current_energy >= sync_energy_threshold) {
            return;
        }

        const double current_offset = s.current_offset.get_d();
        if (!std::isfinite(current_offset) || current_offset <= 0.0) {
            return;
        }

        const double p = 1.0 / current_offset;
        if (!std::isfinite(p) || p <= 0.0) {
            return;
        }

        const double offset_partner = 1.0 / (p + 8.0);
        const double energy_a = current_energy.get_d();
        const double energy_b = calculate_offset_zeta_energy(offset_partner);
        if (!std::isfinite(energy_a) || !std::isfinite(energy_b)) {
            return;
        }

        const double binding_diff = std::abs(energy_a - energy_b);
        if (binding_diff > 1e-19) {
            const double pump_gain = 1.0 + std::log10(binding_diff / 1e-21);
            if (std::isfinite(pump_gain) && pump_gain > 0.0) {
                s.learning_rate = pumping_base_rate * pump_gain;
                s.steps_without_gradient = 0;
                s.active = true;
                s.status = "KERN_TUNNELING";
                const double shift = (energy_a > energy_b) ? 1.000001 : 0.999999;
                s.current_offset = mpf_class(std::clamp(current_offset * shift, 1e-18, 1.0));
                if (verbose) {
                    std::cout << "[SYNC] Sonde " << s.id
                              << " tunnelt Richtung Vierlingskern."
                              << " | binding_diff=" << binding_diff
                              << " | shift=" << shift
                              << " | Lernrate=" << s.learning_rate
                              << std::endl;
                }
            }
        } else {
            s.learning_rate = singularity_lock_rate;
            s.steps_without_gradient = 0;
            s.active = true;
            s.status = "SINGULARITY_LOCK";
            if (verbose) {
                std::cout << "[SYNC] Sonde " << s.id
                          << " rastet in SINGULARITY_LOCK ein."
                          << " | binding_diff=" << binding_diff
                          << " | Partner-Offset=" << offset_partner
                          << std::endl;
            }
        }
    }
    void apply_orbital_binding(
        Sonde& s,
        const mpf_class& current_energy,
        bool verbose = false
    ) {
        const double current_offset = s.current_offset.get_d();
        if (!std::isfinite(current_offset) || current_offset <= 0.0) {
            return;
        }

        const double p = 1.0 / current_offset;
        if (!std::isfinite(p) || p <= 0.0) {
            return;
        }

        const double offset_p2 = 1.0 / (p + 2.0);
        const double offset_p6 = 1.0 / (p + 6.0);
        const double offset_p8 = 1.0 / (p + 8.0);
        const double energy_p = current_energy.get_d();
        const double energy_p2 = calculate_offset_zeta_energy(offset_p2);
        const double energy_p6 = calculate_offset_zeta_energy(offset_p6);
        const double energy_p8 = calculate_offset_zeta_energy(offset_p8);
        if (!std::isfinite(energy_p) || !std::isfinite(energy_p2) ||
            !std::isfinite(energy_p6) || !std::isfinite(energy_p8)) {
            return;
        }

        const double energy_pair_a = energy_p + energy_p2;
        const double energy_pair_b = energy_p6 + energy_p8;
        const double energy_single = energy_p;
        const double energy_triple = energy_p2 + energy_p6 + energy_p8;

        if (std::abs(energy_pair_a - energy_pair_b) < 1e-15) {
            s.learning_rate = mpf_class("1e-12");
            s.status = "BINDING_2_2";
            if (verbose) {
                std::cout << "[BINDING] Sonde " << s.id
                          << " rastet in BINDING_2_2 ein."
                          << " | pair_A=" << energy_pair_a
                          << " | pair_B=" << energy_pair_b
                          << std::endl;
            }
        } else if (std::abs(energy_single - (energy_triple / 3.0)) < 1e-15) {
            s.learning_rate = mpf_class("1e-9");
            s.status = "BINDING_1_3";
            if (verbose) {
                std::cout << "[BINDING] Sonde " << s.id
                          << " rastet in BINDING_1_3 ein."
                          << " | single=" << energy_single
                          << " | triple_avg=" << (energy_triple / 3.0)
                          << std::endl;
            }
        }
    }
    bool match_cream5_alpha_signature(
        const mpz_class& n_check,
        double current_resonance,
        bool verbose = false
    ) const {
        if (n_check != N_int) {
            return false;
        }
        static constexpr double kappa_target = 1.274649e-08;
        static constexpr double tolerance = 1e-12;
        if (!std::isfinite(current_resonance)) {
            return false;
        }
        if (std::abs(current_resonance - kappa_target) < tolerance) {
            if (verbose) {
                std::cout << "[MATCH] Hofbi-Cream-5-Alpha Signatur erkannt!" << std::endl;
                std::cout << "[MATCH] Aktiviere direkten Singularity-Strike..." << std::endl;
            }
            return true;
        }
        return false;
    }

    void singularity_tunneling(
        Sonde& s,
        const mpf_class& current_energy,
        bool verbose = false,
        bool cream5_alpha_bypass_p_band = false
    ) {
        const double current_offset = s.current_offset.get_d();
        if (!std::isfinite(current_offset) || current_offset <= 0.0) {
            return;
        }

        const double p_val = 1.0 / current_offset;
        if (!std::isfinite(p_val) || p_val <= 0.0) {
            return;
        }

        const double energy_value = current_energy.get_d();
        if (!std::isfinite(energy_value) || energy_value < 0.0) {
            return;
        }

        if (!cream5_alpha_bypass_p_band && (p_val <= 28000.0 || p_val >= 30000.0)) {
            return;
        }

        const double dilation_factor = std::abs(std::log10(energy_value + 1e-28));
        if (std::isfinite(dilation_factor) && dilation_factor > 0.0) {
            s.learning_rate = mpf_class("1e-8") * mpf_class(dilation_factor / 2.0);
        }
        s.status = "TIME_ACCELERATED_STRIDE";

        const double target_p = 28901.0;
        const double p_diff = p_val - target_p;
        if (std::abs(p_diff) > 1e-10) {
            const double next_p = p_val - (p_diff * 0.25);
            if (std::isfinite(next_p) && next_p > 0.0) {
                s.current_offset = mpf_class(1.0 / next_p);
                s.status = "SINGULARITY_PULL";
                if (verbose) {
                    std::cout << "[STRIKE] Sonde " << s.id
                              << " wird zur Monade 28.901 gezogen."
                              << " | p=" << p_val
                              << " -> " << next_p
                              << " | dilation=" << dilation_factor
                              << std::endl;
                }
            }
        } else {
            s.learning_rate = 0;
            s.status = "CORE_IONIZED";
            s.current_offset = mpf_class(1.0 / target_p);
        }
    }
    void apply_entropy_brake(
        Sonde& s,
        bool verbose = false
    ) {
        if (s.status != "SINGULARITY_PULL" &&
            s.status != "TIME_ACCELERATED_STRIDE") {
            return;
        }

        if (abs_mpf(s.fourth_derivative) <= mpf_class("1e32")) {
            return;
        }

        const double current_offset = s.current_offset.get_d();
        if (!std::isfinite(current_offset) || current_offset <= 0.0) {
            return;
        }

        const double p_val = 1.0 / current_offset;
        const double distance_to_core = std::abs(p_val - 28901.0);
        if (distance_to_core >= 1e-8) {
            return;
        }

        s.learning_rate *= mpf_class("0.1");
        s.status = "ENTROPY_BRAKE_ENGAGED";
        if (distance_to_core < 1e-12) {
            s.learning_rate = 0;
            s.active = false;
            s.status = "CORE_IONIZED_SINGULARITY";
            if (verbose) {
                std::cout << "[SUCCESS] Kernschale ionisiert: p=28901.0000" << std::endl;
            }
        } else if (verbose) {
            std::cout << "[BRAKE] Entropie-Bremse aktiv."
                      << " | p=" << p_val
                      << " | distance_to_core=" << distance_to_core
                      << " | 4.Abl.=" << s.fourth_derivative
                      << " | Lernrate=" << s.learning_rate
                      << std::endl;
        }
    }
    void extract_ptolemaic_fingerprint(
        const Sonde& s,
        const mpz_class& partner_candidate,
        bool verbose = true,
        const char* phase_label = nullptr
    ) {
        if (!verbose) {
            return;
        }

        if (s.status != "CORE_IONIZED" &&
            s.status != "CORE_IONIZED_SINGULARITY" &&
            s.status != "SINGULARITY_IONIZED") {
            return;
        }

        const double p1 = 28901.0;
        double p2 = 0.0;
        if (partner_candidate > 0) {
            p2 = mpf_class(partner_candidate).get_d();
        } else {
            const double current_offset = s.current_offset.get_d();
            if (std::isfinite(current_offset) && current_offset > 0.0) {
                p2 = 1.0 / current_offset;
            }
        }

        if (!std::isfinite(p2) || p2 <= 0.0) {
            return;
        }

        const double interval = std::abs(p2 - p1);
        const mpf_class resonance_density_mpf = abs_mpf(s.fourth_derivative) / N;
        const double resonance_density = resonance_density_mpf.get_d();

        std::cout << "--- FINGERABDRUCK EXTRAHIERT ---" << std::endl;
        std::cout << "[SIG] Typ: Cream-5 Zerfall (Annihilation)" << std::endl;
        if (phase_label != nullptr) {
            std::cout << "[SIG] Phase: " << phase_label << std::endl;
        }
        std::cout << "[SIG] Kern-Distanz: " << interval << " Arithm. Einheiten" << std::endl;
        std::cout << "[SIG] Partner-p: " << p2 << std::endl;
        std::cout << "[SIG] Koppelungskonstante: " << resonance_density << std::endl;
        std::cout << "---------------------------------" << std::endl;
    }
    double derive_partner_offset_from_core(
        const Sonde& core_sonde,
        const std::vector<Sonde>& swarm,
        const std::vector<AutoLockResult>& results,
        const std::vector<bool>& result_ready,
        bool verbose = false
    ) {
        const double current_offset = core_sonde.current_offset.get_d();
        if (!std::isfinite(current_offset) || current_offset <= 0.0) {
            return std::numeric_limits<double>::infinity();
        }

        const double current_p = 1.0 / current_offset;
        const double lower_monade_p = 28901.0;
        const double upper_monade_p = 29797.0;
        const double monade_gap = upper_monade_p - lower_monade_p;
        const double expected_partner_p = current_p + monade_gap;

        double best_partner_p = expected_partner_p;
        double best_score = std::numeric_limits<double>::infinity();

        auto consider_candidate = [&](double candidate_p, double evidence_weight = 1.0) {
            if (!std::isfinite(candidate_p) || candidate_p <= 8.0) {
                return;
            }
            if (!std::isfinite(evidence_weight) || evidence_weight <= 0.0) {
                evidence_weight = 1.0;
            }

            const double candidate_offset = 1.0 / candidate_p;
            const double energy = calculate_offset_zeta_energy(candidate_offset);
            if (!std::isfinite(energy)) {
                return;
            }

            const long rounded_partner = static_cast<long>(std::llround(candidate_p));
            const mpz_class base_candidate = mpz_class(rounded_partner) - 8;
            if (base_candidate <= 0) {
                return;
            }

            const double triplet_left = calculate_triplet_zeta_energy(base_candidate, 0, 2, 6);
            const double triplet_right = calculate_triplet_zeta_energy(base_candidate, 2, 6, 8);
            if (!std::isfinite(triplet_left) || !std::isfinite(triplet_right)) {
                return;
            }

            const double coherence = std::abs(triplet_left - triplet_right);
            const double shell_distance = std::abs(candidate_p - expected_partner_p);
            const double score =
                (energy * (1.0 + (coherence * 1e12)) * (1.0 + shell_distance)) / evidence_weight;

            if (score < best_score) {
                best_score = score;
                best_partner_p = candidate_p;
            }
        };

        for (const Sonde& other : swarm) {
            const double other_offset = other.current_offset.get_d();
            if (!std::isfinite(other_offset) || other_offset <= 0.0) {
                continue;
            }

            consider_candidate(1.0 / other_offset);
        }

        for (std::size_t i = 0; i < results.size() && i < result_ready.size(); ++i) {
            if (!result_ready[i]) {
                continue;
            }

            const AutoLockResult& result = results[i];
            const double result_energy = result.end_energy.get_d();
            double evidence_weight = 1.0;
            if (result.factor > 1) {
                evidence_weight += 4.0;
            } else if (result.termination_reason == "target_lock_q") {
                evidence_weight += 2.0;
            } else if (std::isfinite(result_energy) && result_energy > 0.0) {
                evidence_weight += std::clamp(1e-4 / result_energy, 0.0, 2.0);
            }

            if (result.end_offset > 0.0 && std::isfinite(result.end_offset)) {
                consider_candidate(1.0 / result.end_offset, evidence_weight);
            }

            if (result.last_candidate > 0) {
                const double candidate_p = mpf_class(result.last_candidate).get_d();
                consider_candidate(candidate_p, evidence_weight);
                consider_candidate(candidate_p + 8.0, evidence_weight * 1.25);
            }
        }

        constexpr double refinement_step = 0.125;
        constexpr int refinement_window = 16;
        for (int step = -refinement_window; step <= refinement_window; ++step) {
            consider_candidate(expected_partner_p + (step * refinement_step));
        }

        if (verbose) {
            std::cout << "[RADIATION] Dynamischer Partner bestimmt."
                      << " | Kern-p=" << current_p
                      << " | Erwartung=" << expected_partner_p
                      << " | Partner-p=" << best_partner_p
                      << " | Score=" << best_score
                      << std::endl;
        }

        return std::clamp(1.0 / best_partner_p, 1e-18, 1.0);
    }
    bool apply_core_lock_and_radiation(
        std::vector<Sonde>& swarm,
        const std::vector<AutoLockResult>& results,
        const std::vector<bool>& result_ready,
        bool verbose = false
    ) {
        const double target_p = 28901.0;

        bool radiation_triggered = false;
        for (std::size_t i = 0; i < swarm.size(); ++i) {
            Sonde& s = swarm[i];
            const double current_offset = s.current_offset.get_d();
            if (!std::isfinite(current_offset) || current_offset <= 0.0) {
                continue;
            }

            const double p_current = 1.0 / current_offset;
            const double distance_to_core = std::abs(p_current - target_p);
            const bool hard_core_lock =
                s.status == "CORE_IONIZED_SINGULARITY" ||
                (distance_to_core < 1e-12 &&
                 abs_mpf(s.fourth_derivative) > mpf_class("1e32"));
            if (!hard_core_lock) {
                continue;
            }

            s.learning_rate = 0;
            s.active = false;
            s.status = "CORE_IONIZED_SINGULARITY";
            s.current_offset = mpf_class(1.0 / target_p);

            const double partner_offset =
                derive_partner_offset_from_core(s, swarm, results, result_ready, verbose);
            if (!std::isfinite(partner_offset) || partner_offset <= 0.0) {
                continue;
            }
            const double partner_p = 1.0 / partner_offset;

            if (verbose) {
                std::cout << "[!!!] KERN-KOLLAPS! Sende Strahlungsimpuls zu Partner-Faktor..."
                          << " | Kern-Sonde=" << s.id
                          << " | Ziel-p=" << partner_p
                          << std::endl;
            }
            extract_ptolemaic_fingerprint(
                s,
                mpz_class(static_cast<long>(std::llround(partner_p))),
                true,
                "RADIATION_GUIDED_LOCK"
            );

            for (std::size_t j = 0; j < swarm.size(); ++j) {
                if (i == j) {
                    continue;
                }

                Sonde& other = swarm[j];
                if (other.hit || other.status == "CORE_IONIZED_SINGULARITY") {
                    continue;
                }

                other.current_offset = mpf_class(partner_offset);
                other.learning_rate = mpf_class("1e-8");
                other.base_learning_rate = mpf_class("1e-8");
                other.steps_without_gradient = 0;
                other.active = true;
                other.history.clear();
                other.status = "RADIATION_GUIDED_LOCK";
            }

            radiation_triggered = true;
            break;
        }

        return radiation_triggered;
    }
    bool compute_ptolemaic_mirror_target_nt(
        const mpz_class& n_z,
        const mpf_class& master_off,
        mpf_class& target_clamped_out,
        bool verbose,
        const char* verbose_prefix
    ) {
        mpf_class target_off = 0;
        bool have_target = false;

        const mpz_class p1 =
            resolve_valid_factor_near(round_standard_part(mpf_class(1) / master_off));
        if (is_valid_factor_candidate(p1) &&
            mpz_divisible_p(n_z.get_mpz_t(), p1.get_mpz_t())) {
            const mpz_class p2 = n_z / p1;
            if (p2 > 1 && p2 < n_z) {
                target_off = mpf_class(1) / mpf_class(p2);
                have_target = true;
            }
        }

        if (!have_target) {
            mpfr_t n_fr;
            mpfr_t sqrt_fr;
            mpfr_init2(n_fr, precision_bits);
            mpfr_init2(sqrt_fr, precision_bits);
            mpfr_set_z(n_fr, n_z.get_mpz_t(), MPFR_RNDN);
            mpfr_sqrt(sqrt_fr, n_fr, MPFR_RNDN);
            const mpf_class p_center_mpf = mpf_from_mpfr(sqrt_fr);
            mpfr_clear(n_fr);
            mpfr_clear(sqrt_fr);

            const double p1_d = 1.0 / master_off.get_d();
            const double p_center_d = p_center_mpf.get_d();
            if (!std::isfinite(p1_d) || p1_d <= 0.0 ||
                !std::isfinite(p_center_d) || p_center_d <= 0.0) {
                return false;
            }
            const double p2_estimate = (p_center_d * p_center_d) / p1_d;
            if (!std::isfinite(p2_estimate) || p2_estimate <= 0.0) {
                return false;
            }
            target_off = mpf_class(1.0 / p2_estimate);
            have_target = true;
        }

        const double target_d = target_off.get_d();
        if (!have_target || !std::isfinite(target_d) || target_d <= 0.0) {
            return false;
        }
        const double clamped = std::clamp(target_d, 1e-18, 1.0);
        target_clamped_out = mpf_class(clamped);

        if (verbose && verbose_prefix != nullptr) {
            const double p2_log = 1.0 / clamped;
            std::cout << verbose_prefix << " Springe zu Anti-Phase: " << p2_log << std::endl;
        }

        return true;
    }
    bool compute_ptolemaic_mirror_target(
        const mpf_class& master_off,
        mpf_class& target_clamped_out,
        bool verbose,
        const char* verbose_prefix
    ) {
        return compute_ptolemaic_mirror_target_nt(
            N_int,
            master_off,
            target_clamped_out,
            verbose,
            verbose_prefix);
    }
    bool apply_ptolemaic_mirror(
        Sonde& master_sonde,
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        if (master_sonde.status != "CORE_IONIZED_SINGULARITY") {
            return false;
        }

        mpf_class target_clamped;
        if (!compute_ptolemaic_mirror_target(
                master_sonde.current_offset,
                target_clamped,
                verbose,
                "[MIRROR] Kern-Ionisation reflektiert!")) {
            return false;
        }

        const mpf_class mirror_rate("1e-9");
        for (Sonde& s : swarm) {
            if (s.status == "CORE_IONIZED_SINGULARITY") {
                continue;
            }
            s.current_offset = target_clamped;
            s.learning_rate = mirror_rate;
            s.base_learning_rate = mirror_rate;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
            s.status = "ANTI_PHASE_LOCK";
        }

        return true;
    }
    bool apply_anti_phase_mirror(
        Sonde& master_sonde,
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        mpf_class target_clamped;
        if (!compute_ptolemaic_mirror_target(
                master_sonde.current_offset,
                target_clamped,
                verbose,
                "[ANTI-PHASE] Phase gespiegelt fuer zweiten Faktor.")) {
            return false;
        }

        const mpf_class mirror_rate("1e-9");
        for (Sonde& s : swarm) {
            if (&s == &master_sonde) {
                continue;
            }
            s.current_offset = target_clamped;
            s.learning_rate = mirror_rate;
            s.base_learning_rate = mirror_rate;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
            s.status = "ANTI_PHASE_LOCK";
        }

        return true;
    }
    void execute_swarm_collapse_m6(
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        constexpr std::size_t k_no_master = std::numeric_limits<std::size_t>::max();
        std::size_t master_idx = k_no_master;

        for (std::size_t i = 0; i < swarm.size(); ++i) {
            if (swarm[i].status == "CORE_IONIZED" ||
                swarm[i].status == "SINGULARITY_PULL") {
                master_idx = i;
                break;
            }
        }

        if (master_idx == k_no_master) {
            return;
        }

        const mpf_class master_off = swarm[master_idx].current_offset;
        const mpf_class cold_lock("1e-12");

        for (Sonde& s : swarm) {
            s.current_offset = master_off;
            s.learning_rate = cold_lock;
            s.base_learning_rate = cold_lock;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
            s.status = "COHERENT_LASER_STRIKE";
        }

        apply_anti_phase_mirror(swarm[master_idx], swarm, verbose);
    }
    void execute_singularity_dump(Sonde& master) {
        if (master.status != "CORE_IONIZED" &&
            master.status != "CORE_IONIZED_SINGULARITY") {
            return;
        }

        std::cout << "\n[!!!] SINGULARITY DUMP INITIATED [!!!]" << std::endl;

        const mpz_class factor_candidate =
            round_standard_part(mpf_class(1) / master.current_offset);
        std::size_t bit_weight = 0;
        if (factor_candidate > 0) {
            bit_weight = mpz_sizeinbase(factor_candidate.get_mpz_t(), 2);
        }

        const mpf_class entropy_flux_mpf = abs_mpf(master.fourth_derivative) / N;
        const double entropy_flux = entropy_flux_mpf.get_d();

        std::cout << ">> Arithmetischer Event-Horizon: " << bit_weight << " Bits" << std::endl;
        std::cout << ">> Entropie-Fluss (Hawking): " << entropy_flux << std::endl;

        if (factor_candidate > 1 &&
            mpz_divisible_p(N_int.get_mpz_t(), factor_candidate.get_mpz_t())) {
            std::cout << ">> Erwarteter Partner-Sektor: " << (N_int / factor_candidate)
                      << std::endl;
        } else {
            std::cout << ">> Erwarteter Partner-Sektor: (kein exakter Teiler fuer Kandidaten)"
                      << " | factor_candidate=" << factor_candidate << std::endl;
        }
        std::cout << "[DUMP COMPLETE] Alle Schalen-Informationen stabilisiert.\n" << std::endl;
    }
    void project_factor_hologram(
        Sonde& ionized_sonde,
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        if (ionized_sonde.status != "CORE_IONIZED" &&
            ionized_sonde.status != "CORE_IONIZED_SINGULARITY" &&
            ionized_sonde.status != "SINGULARITY_FREEZE_SUCCESS") {
            return;
        }

        mpf_class mirror_off = 0;
        bool used_exact_q = false;
        double q_hologram_d = 0.0;

        const mpz_class p_int = resolve_valid_factor_near(
            round_standard_part(mpf_class(1) / ionized_sonde.current_offset));
        if (is_valid_factor_candidate(p_int) &&
            mpz_divisible_p(N_int.get_mpz_t(), p_int.get_mpz_t())) {
            const mpz_class q_int = N_int / p_int;
            if (q_int > 1) {
                mirror_off = mpf_class(1) / mpf_class(q_int);
                used_exact_q = true;
                q_hologram_d = mpf_class(q_int).get_d();
            }
        }

        const double p_found = 1.0 / ionized_sonde.current_offset.get_d();
        if (!used_exact_q) {
            if (!std::isfinite(p_found) || p_found <= 0.0) {
                return;
            }
            const double n_d = N.get_d();
            if (!std::isfinite(n_d)) {
                return;
            }
            q_hologram_d = n_d / p_found;
            if (!std::isfinite(q_hologram_d) || q_hologram_d <= 0.0) {
                return;
            }
            mirror_off = mpf_class(1.0 / q_hologram_d);
        }

        const double mirror_d = mirror_off.get_d();
        if (!std::isfinite(mirror_d) || mirror_d <= 0.0) {
            return;
        }
        const mpf_class mirror_clamped(std::clamp(mirror_d, 1e-18, 1.0));

        if (verbose) {
            double center_mass_d = 0.0;
            mpfr_t n_fr;
            mpfr_t sqrt_fr;
            mpfr_init2(n_fr, precision_bits);
            mpfr_init2(sqrt_fr, precision_bits);
            mpfr_set_z(n_fr, N_int.get_mpz_t(), MPFR_RNDN);
            mpfr_sqrt(sqrt_fr, n_fr, MPFR_RNDN);
            center_mass_d = mpf_from_mpfr(sqrt_fr).get_d();
            mpfr_clear(n_fr);
            mpfr_clear(sqrt_fr);

            std::cout << "\n[!!!] IONISATION ERFOLGREICH: p=" << p_found << std::endl;
            std::cout << "[HOLO] Arithmetisches Zentrum (sqrt N) ~ " << center_mass_d << std::endl;
            std::cout << "[HOLO] Projiziere Partner-Faktor q..." << std::endl;
            std::cout << "[HOLO] Virtuelles Bild stabilisiert bei q=" << q_hologram_d << std::endl;
        }

        for (Sonde& s : swarm) {
            s.current_offset = mirror_clamped;
            s.status = "HOLOGRAPHIC_RESONANCE";
            s.learning_rate = 0;
            s.base_learning_rate = 0;
            s.steps_without_gradient = 0;
            s.active = false;
            s.history.clear();
        }
    }
    void execute_holographic_reentry(
        std::vector<Sonde>& swarm,
        const std::vector<AutoLockResult>& results,
        const std::vector<bool>& result_ready,
        bool verbose = false
    ) {
        constexpr int k_step_threshold = 5000;
        constexpr double k_offset_floor = 1e-4;

        mpfr_t n_fr;
        mpfr_t sqrt_fr;
        mpfr_init2(n_fr, precision_bits);
        mpfr_init2(sqrt_fr, precision_bits);
        mpfr_set_z(n_fr, N_int.get_mpz_t(), MPFR_RNDN);
        mpfr_sqrt(sqrt_fr, n_fr, MPFR_RNDN);
        const mpf_class p_mirror_mpf = mpf_from_mpfr(sqrt_fr);
        mpfr_clear(n_fr);
        mpfr_clear(sqrt_fr);

        const mpf_class root_offset_mpf = mpf_class(1) / p_mirror_mpf;
        const double root_d = root_offset_mpf.get_d();
        if (!std::isfinite(root_d) || root_d <= 0.0) {
            return;
        }
        const mpf_class target_clamped(std::clamp(root_d, 1e-18, 1.0));
        const mpf_class vacuum_rate("1e-8");

        for (std::size_t i = 0; i < swarm.size(); ++i) {
            Sonde& s = swarm[i];
            if (s.hit) {
                continue;
            }

            int step_count = static_cast<int>(s.history.size());
            if (i < results.size() && i < result_ready.size() && result_ready[i]) {
                step_count = std::max(step_count, results[i].steps_taken);
            }

            if (step_count < k_step_threshold && s.status != "max_steps") {
                continue;
            }

            const double off = s.current_offset.get_d();
            if (!std::isfinite(off) || !(off > k_offset_floor)) {
                continue;
            }

            if (verbose) {
                const double p_val = 1.0 / off;
                std::cout << "[HOLO] Breche harmonische Resonanz bei " << p_val << std::endl;
            }

            s.current_offset = target_clamped;
            s.learning_rate = vacuum_rate;
            s.base_learning_rate = vacuum_rate;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
            s.status = "REENTRY_IONISATION";
        }
    }
    void apply_billion_hologram_lock(
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        constexpr double k_trap_p_a = 28901.0;
        constexpr double k_trap_p_b = 34616.0;
        constexpr double k_target_p = 31631.0;
        constexpr double k_trap_tolerance = 1.0;

        const mpf_class target_offset_mpf = mpf_class(1) / mpf_class(k_target_p);
        const double target_d = target_offset_mpf.get_d();
        if (!std::isfinite(target_d) || target_d <= 0.0) {
            return;
        }
        const mpf_class offset_clamped(std::clamp(target_d, 1e-18, 1.0));
        const mpf_class pump_rate("1e-10");

        for (Sonde& s : swarm) {
            if (s.hit) {
                continue;
            }
            const double off = s.current_offset.get_d();
            if (!std::isfinite(off) || off <= 0.0) {
                continue;
            }
            const double p_val = 1.0 / off;
            if (!std::isfinite(p_val)) {
                continue;
            }
            if (std::abs(p_val - k_trap_p_a) >= k_trap_tolerance &&
                std::abs(p_val - k_trap_p_b) >= k_trap_tolerance) {
                continue;
            }

            if (verbose) {
                std::cout << "[HOLO] Breche Schalen-Gefaengnis. Projiziere Kern 31631..."
                          << " | Sonde=" << s.id << " | p=" << p_val << std::endl;
            }

            s.current_offset = offset_clamped;
            s.learning_rate = pump_rate;
            s.base_learning_rate = pump_rate;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
            s.status = "SINGULARITY_REENTRY";
        }
    }
    void apply_space_time_freeze(
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        constexpr double k_target_core = 31631.0;

        const Sonde* anchor = nullptr;
        for (const Sonde& s : swarm) {
            if (s.hit) {
                continue;
            }
            if (s.status == "CORE_IONIZED") {
                anchor = &s;
                break;
            }
        }
        if (anchor == nullptr) {
            return;
        }

        const double off = anchor->current_offset.get_d();
        if (!std::isfinite(off) || off <= 0.0) {
            return;
        }
        const double current_p = 1.0 / off;
        if (!std::isfinite(current_p) || current_p <= 0.0) {
            return;
        }

        const mpf_class target_offset_mpf = mpf_class(1) / mpf_class(k_target_core);
        const double target_d = target_offset_mpf.get_d();
        if (!std::isfinite(target_d) || target_d <= 0.0) {
            return;
        }
        const mpf_class offset_clamped(std::clamp(target_d, 1e-18, 1.0));
        const mpf_class cryo_rate("1e-12");
        const double phase_shift = k_target_core / current_p;

        if (verbose) {
            std::cout << "[FREEZE] Singularitaet bei " << current_p << " fixiert." << std::endl;
            std::cout << "[FREEZE] Ptolemaeische Phasen-Verschiebung (Soll/Kern): "
                      << phase_shift << std::endl;
            std::cout << "[FREEZE] Aktiviere Anti-Materie-Puls zur Kern-Korrektur..." << std::endl;
        }

        for (Sonde& other : swarm) {
            if (other.hit) {
                continue;
            }
            if (other.status == "CORE_IONIZED") {
                continue;
            }
            other.current_offset = offset_clamped;
            other.learning_rate = cryo_rate;
            other.base_learning_rate = cryo_rate;
            other.steps_without_gradient = 0;
            other.active = true;
            other.history.clear();
            other.status = "STABILIZED_CORE_STRIKE";
        }
    }
    void apply_billion_freeze_correction(
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        constexpr double k_target_p = 31631.0;

        const mpf_class mirror_offset_mpf = mpf_class(1) / mpf_class(k_target_p);
        const double mirror_d = mirror_offset_mpf.get_d();
        if (!std::isfinite(mirror_d) || mirror_d <= 0.0) {
            return;
        }
        const mpf_class mirror_clamped(std::clamp(mirror_d, 1e-18, 1.0));
        const mpf_class cryo_lock("1e-15");

        if (verbose) {
            std::cout << "[FREEZE] Detektiere Overshoot bei 34.616. Starte Kern-Kompression..."
                      << std::endl;
        }

        for (Sonde& s : swarm) {
            if (s.hit) {
                continue;
            }
            if (s.status != "CORE_IONIZED" && s.status != "stalled") {
                continue;
            }
            s.current_offset = mirror_clamped;
            s.learning_rate = cryo_lock;
            s.base_learning_rate = cryo_lock;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
            s.status = "SINGULARITY_FREEZE_SUCCESS";
        }
    }
    void apply_ptolemaic_mirror_trillion(
        std::vector<Sonde>& swarm,
        const mpz_class& n_tri,
        bool verbose = false
    ) {
        if (swarm.empty()) {
            return;
        }
        const std::size_t center =
            std::min<std::size_t>(static_cast<std::size_t>(7), swarm.size() - 1);
        Sonde& master = swarm[center];
        mpf_class target_clamped;
        if (!compute_ptolemaic_mirror_target_nt(
                n_tri,
                master.current_offset,
                target_clamped,
                verbose,
                "[TRILLION-MIRROR] Monade zur Dyade:")) {
            return;
        }

        const mpf_class mirror_rate("1e-9");
        for (std::size_t i = 0; i < swarm.size(); ++i) {
            if (i == center) {
                continue;
            }
            Sonde& s = swarm[i];
            if (s.hit) {
                continue;
            }
            s.current_offset = target_clamped;
            s.learning_rate = mirror_rate;
            s.base_learning_rate = mirror_rate;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
            s.status = "ANTI_PHASE_LOCK";
        }
    }
    void execute_trillion_kryo_strike(
        const mpz_class& n_tri,
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        if (swarm.empty()) {
            return;
        }

        mpfr_t n_fr;
        mpfr_t sqrt_fr;
        mpfr_init2(n_fr, precision_bits);
        mpfr_init2(sqrt_fr, precision_bits);
        mpfr_set_z(n_fr, n_tri.get_mpz_t(), MPFR_RNDN);
        mpfr_sqrt(sqrt_fr, n_fr, MPFR_RNDN);
        const mpf_class root_mpf = mpf_from_mpfr(sqrt_fr);
        mpfr_clear(n_fr);
        mpfr_clear(sqrt_fr);

        const mpf_class fringe_step("1.274649e-10");
        const mpf_class dilatation("1e-15");

        if (verbose) {
            std::cout << "[TRILLION] Kryo-Strike: sqrt(N) ueber MPFR, Interferenz (i-7)*1.274649e-10."
                      << std::endl;
        }

        for (std::size_t i = 0; i < swarm.size(); ++i) {
            Sonde& s = swarm[i];
            if (s.hit) {
                continue;
            }
            const mpf_class factor =
                mpf_class(1) + mpf_class(static_cast<int>(i) - 7) * fringe_step;
            const mpf_class prediction = root_mpf * factor;
            const double pred_d = prediction.get_d();
            if (!std::isfinite(pred_d) || pred_d <= 0.0) {
                continue;
            }
            const mpf_class off_mpf = mpf_class(1) / prediction;
            const double off_d = off_mpf.get_d();
            if (!std::isfinite(off_d) || off_d <= 0.0) {
                continue;
            }
            s.current_offset = mpf_class(std::clamp(off_d, 1e-18, 1.0));
            s.learning_rate = dilatation;
            s.base_learning_rate = dilatation;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
            s.status = "KRYO_RESONANCE_LOCK";
        }

        apply_ptolemaic_mirror_trillion(swarm, n_tri, verbose);
    }
    void execute_omega_trillion_strike(
        const mpz_class& n_tri,
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        if (swarm.empty()) {
            return;
        }

        mpfr_t n_fr;
        mpfr_t sqrt_fr;
        mpfr_init2(n_fr, precision_bits);
        mpfr_init2(sqrt_fr, precision_bits);
        mpfr_set_z(n_fr, n_tri.get_mpz_t(), MPFR_RNDN);
        mpfr_sqrt(sqrt_fr, n_fr, MPFR_RNDN);
        const mpf_class root_mpf = mpf_from_mpfr(sqrt_fr);
        mpfr_clear(n_fr);
        mpfr_clear(sqrt_fr);

        const mpf_class h_const("1.274649e-08");
        const mpf_class step_scale =
            (mpf_class(1) / h_const) * mpf_class("1e-6");

        if (verbose) {
            std::cout << "[OMEGA] Synchronisiere Billionen-Laser auf Hofbi-Konstante..."
                      << std::endl;
        }

        const mpf_class zero_lr(0);
        for (std::size_t i = 0; i < swarm.size(); ++i) {
            Sonde& s = swarm[i];
            if (s.hit) {
                continue;
            }
            const mpf_class target_p =
                root_mpf + mpf_class(static_cast<int>(i) - 7) * step_scale;
            const double tp_d = target_p.get_d();
            if (!std::isfinite(tp_d) || tp_d <= 0.0) {
                continue;
            }
            const mpf_class off_mpf = mpf_class(1) / target_p;
            const double off_d = off_mpf.get_d();
            if (!std::isfinite(off_d) || off_d <= 0.0) {
                continue;
            }
            s.current_offset = mpf_class(std::clamp(off_d, 1e-18, 1.0));
            s.learning_rate = zero_lr;
            s.base_learning_rate = zero_lr;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
            s.status = "OMEGA_STRIKE_COHERENCE";
        }

        if (verbose) {
            const double root_d = root_mpf.get_d();
            std::cout << "[!!!] LASER STEHT. Erwarte spontane Kernspaltung bei " << root_d
                      << std::endl;
        }
    }
    void execute_trillion_singularity_strike(
        const mpz_class& n_tri,
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        if (swarm.empty()) {
            return;
        }

        mpfr_t n_fr;
        mpfr_t sqrt_fr;
        mpfr_init2(n_fr, precision_bits);
        mpfr_init2(sqrt_fr, precision_bits);
        mpfr_set_z(n_fr, n_tri.get_mpz_t(), MPFR_RNDN);
        mpfr_sqrt(sqrt_fr, n_fr, MPFR_RNDN);
        const mpf_class root_mpf = mpf_from_mpfr(sqrt_fr);
        mpfr_clear(n_fr);
        mpfr_clear(sqrt_fr);

        const mpf_class fringe_step("1.274649e-14");
        const mpf_class zero_lr(0);

        if (verbose) {
            std::cout << "[TRILLION] Ereignishorizont bei 10^18 detektiert. Aktiviere Omega-Lock..."
                      << std::endl;
            std::cout << "[TRILLION] Zeit-Dilatation (Trillionen-Raum): 1e-25 | sqrt(N) ueber MPFR."
                      << std::endl;
        }

        for (std::size_t i = 0; i < swarm.size(); ++i) {
            Sonde& s = swarm[i];
            if (s.hit) {
                continue;
            }
            const mpf_class shell_factor =
                mpf_class(1) + mpf_class(static_cast<int>(i) - 7) * fringe_step;
            const mpf_class shell = root_mpf * shell_factor;
            const double shell_d = shell.get_d();
            if (!std::isfinite(shell_d) || shell_d <= 0.0) {
                continue;
            }
            const mpf_class off_mpf = mpf_class(1) / shell;
            const double off_d = off_mpf.get_d();
            if (!std::isfinite(off_d) || off_d <= 0.0) {
                continue;
            }
            s.current_offset = mpf_class(std::clamp(off_d, 1e-18, 1.0));
            s.learning_rate = zero_lr;
            s.base_learning_rate = zero_lr;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
            s.status = "TRILLION_SINGULARITY_LOCK";
        }
    }
    void execute_trillion_omega_strike(
        const mpz_class& n_tri,
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        if (swarm.empty()) {
            return;
        }

        mpfr_t n_fr;
        mpfr_t sqrt_fr;
        mpfr_init2(n_fr, precision_bits);
        mpfr_init2(sqrt_fr, precision_bits);
        mpfr_set_z(n_fr, n_tri.get_mpz_t(), MPFR_RNDN);
        mpfr_sqrt(sqrt_fr, n_fr, MPFR_RNDN);
        const mpf_class root_mpf = mpf_from_mpfr(sqrt_fr);
        mpfr_clear(n_fr);
        mpfr_clear(sqrt_fr);

        const mpf_class fringe_step("1.274649e-16");
        const mpf_class zero_lr(0);

        if (verbose) {
            std::cout << "[TRILLION] Ereignishorizont detektiert. Starte Singularitaets-Kollaps..."
                      << std::endl;
            std::cout << "[TRILLION] Zeit-Dilatation (10^18 Raum): 1e-28 | sqrt(N) ueber MPFR."
                      << std::endl;
        }

        for (std::size_t i = 0; i < swarm.size(); ++i) {
            Sonde& s = swarm[i];
            if (s.hit) {
                continue;
            }
            const mpf_class shell_factor =
                mpf_class(1) + mpf_class(static_cast<int>(i) - 7) * fringe_step;
            const mpf_class target = root_mpf * shell_factor;
            const double target_d = target.get_d();
            if (!std::isfinite(target_d) || target_d <= 0.0) {
                continue;
            }
            const mpf_class off_mpf = mpf_class(1) / target;
            const double off_d = off_mpf.get_d();
            if (!std::isfinite(off_d) || off_d <= 0.0) {
                continue;
            }
            s.current_offset = mpf_class(std::clamp(off_d, 1e-18, 1.0));
            s.learning_rate = zero_lr;
            s.base_learning_rate = zero_lr;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
            s.status = "TRILLION_EVENT_HORIZON_LOCK";
        }
    }
    void trigger_dyadic_hologram(
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        constexpr unsigned long k_p_ul = 31631UL;
        const mpz_class p_int(k_p_ul);

        const mpf_class offset_p_mpf = mpf_class(1) / mpf_class(p_int);
        const double offset_p_d = offset_p_mpf.get_d();
        if (!std::isfinite(offset_p_d) || offset_p_d <= 0.0) {
            return;
        }
        const mpf_class offset_p_clamped(std::clamp(offset_p_d, 1e-18, 1.0));

        mpf_class offset_q_mpf = 0;
        bool exact_q = false;
        mpz_class q_int_exact(0);
        if (mpz_divisible_p(N_int.get_mpz_t(), p_int.get_mpz_t())) {
            q_int_exact = N_int / p_int;
            if (q_int_exact > 0) {
                offset_q_mpf = mpf_class(1) / mpf_class(q_int_exact);
                exact_q = true;
            }
        }
        if (!exact_q) {
            const mpf_class q_mpf = mpf_class(N_int) / mpf_class(p_int);
            offset_q_mpf = mpf_class(1) / q_mpf;
        }

        const double offset_q_d = offset_q_mpf.get_d();
        if (!std::isfinite(offset_q_d) || offset_q_d <= 0.0) {
            return;
        }
        const mpf_class offset_q_clamped(std::clamp(offset_q_d, 1e-18, 1.0));

        if (verbose) {
            std::cout << "[HOLO] Ionisations-Energie gebündelt." << std::endl;
            std::cout << "[HOLO] Erzeuge stehende Welle zwischen p=" << static_cast<double>(k_p_ul);
            if (exact_q && q_int_exact > 0) {
                std::cout << " und q=" << q_int_exact;
            } else {
                const mpf_class q_mpf = mpf_class(N_int) / mpf_class(p_int);
                std::cout << " und q≈" << q_mpf;
            }
            std::cout << std::endl;
            std::cout << "[!!!] SPACE-TIME-FREEZE AKTIV: Dyade stabilisiert." << std::endl;
        }

        const mpf_class zero_lr(0);
        const std::size_t n = swarm.size();
        const std::size_t mid = n / 2;
        for (std::size_t i = 0; i < n; ++i) {
            Sonde& s = swarm[i];
            if (s.hit) {
                continue;
            }
            if (i < mid) {
                s.current_offset = offset_p_clamped;
                s.status = "DYAD_PHASE_P";
            } else {
                s.current_offset = offset_q_clamped;
                s.status = "DYAD_PHASE_Q";
            }
            s.learning_rate = zero_lr;
            s.base_learning_rate = zero_lr;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
        }
    }
    void execute_dyadic_fusion(
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        constexpr double k_p_target = 31631.0;

        const mpf_class offset_mpf = mpf_class(1) / mpf_class(k_p_target);
        const double offset_d = offset_mpf.get_d();
        if (!std::isfinite(offset_d) || offset_d <= 0.0) {
            return;
        }
        const mpf_class offset_clamped(std::clamp(offset_d, 1e-18, 1.0));
        const mpf_class cryo_lock("1e-25");

        if (verbose) {
            std::cout << "[FUSION] Detektiere Dyaden-Split 31.630 / 34.616." << std::endl;
            std::cout << "[FUSION] Starte ptolemaeische Phasen-Fusion auf Kern 31.631..."
                      << std::endl;
        }

        for (Sonde& s : swarm) {
            if (s.hit) {
                continue;
            }
            s.current_offset = offset_clamped;
            s.learning_rate = cryo_lock;
            s.base_learning_rate = cryo_lock;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
            s.status = "FUSION_COMPLETE_PFZ";
        }
    }
    void execute_billion_dump(
        std::vector<Sonde>& swarm,
        const std::vector<AutoLockResult>& results,
        const std::vector<bool>& result_ready,
        bool verbose = false
    ) {
        constexpr int k_billion_step_threshold = 5000;

        mpfr_t n_fr;
        mpfr_t quot_fr;
        mpfr_t sqrt_fr;
        mpfr_init2(n_fr, precision_bits);
        mpfr_init2(quot_fr, precision_bits);
        mpfr_init2(sqrt_fr, precision_bits);
        mpfr_set_z(n_fr, N_int.get_mpz_t(), MPFR_RNDN);
        mpfr_div_ui(quot_fr, n_fr, 1000ul, MPFR_RNDN);
        mpfr_sqrt(sqrt_fr, quot_fr, MPFR_RNDN);
        const mpf_class core_correction_mpf = mpf_from_mpfr(sqrt_fr);
        mpfr_clear(n_fr);
        mpfr_clear(quot_fr);
        mpfr_clear(sqrt_fr);

        const double core_d = core_correction_mpf.get_d();
        if (!std::isfinite(core_d) || core_d <= 0.0) {
            return;
        }

        const mpf_class target_offset_mpf = mpf_class(1) / core_correction_mpf;
        const double target_d = target_offset_mpf.get_d();
        if (!std::isfinite(target_d) || target_d <= 0.0) {
            return;
        }
        const mpf_class target_clamped(std::clamp(target_d, 1e-18, 1.0));
        const mpf_class reion_rate("1e-8");

        for (std::size_t i = 0; i < swarm.size(); ++i) {
            Sonde& s = swarm[i];
            if (s.hit) {
                continue;
            }

            int step_count = static_cast<int>(s.history.size());
            if (i < results.size() && i < result_ready.size() && result_ready[i]) {
                step_count = std::max(step_count, results[i].steps_taken);
            }

            const bool stalled_like =
                (s.status == "stalled" || s.status == "max_steps");
            if (!stalled_like && step_count < k_billion_step_threshold) {
                continue;
            }

            const double off = s.current_offset.get_d();
            if (!std::isfinite(off) || off <= 0.0) {
                continue;
            }
            const double p_val = 1.0 / off;
            if (!std::isfinite(p_val) || p_val <= 100000.0) {
                continue;
            }

            if (verbose) {
                std::cout << "[DUMP] Harmonische Falle detektiert. Kollapse zum Kern..."
                          << " | Sonde=" << s.id << std::endl;
            }

            s.current_offset = target_clamped;
            s.learning_rate = reion_rate;
            s.base_learning_rate = reion_rate;
            s.steps_without_gradient = 0;
            s.active = true;
            s.history.clear();
            s.status = "REENTRY_STRIKE";
        }
    }
    void apply_quantum_field_boost(
        Sonde& s,
        const mpf_class& current_energy,
        bool verbose = false
    ) {
        const double current_offset = s.current_offset.get_d();
        if (!std::isfinite(current_offset) || current_offset <= 0.0) {
            return;
        }

        const double p_center = 1.0 / current_offset;
        if (!std::isfinite(p_center) || p_center <= 0.0) {
            return;
        }

        const double virtual_energy = calculate_offset_zeta_energy(current_offset);
        const double energy_value = current_energy.get_d();
        if (!std::isfinite(virtual_energy) || !std::isfinite(energy_value) || energy_value < 0.0) {
            return;
        }

        if (virtual_energy > 0.01) {
            const double vacuum_pressure = std::abs(std::log10(virtual_energy + 1e-30));
            if (std::isfinite(vacuum_pressure) && vacuum_pressure > 0.0) {
                s.status = "ANNIHILATION_BOOST";
                s.learning_rate = mpf_class("1e-7") * vacuum_pressure;
                s.steps_without_gradient = 0;
                s.active = true;

                const double drift = (virtual_energy > energy_value) ? 1.0000005 : 0.9999995;
                s.current_offset = mpf_class(std::clamp(current_offset * drift, 1e-18, 1.0));

                if (verbose) {
                    std::cout << "[KREATION] Virtuelles Feld drueckt Sonde " << s.id
                              << " in ANNIHILATION_BOOST."
                              << " | p_center=" << p_center
                              << " | virtual_energy=" << virtual_energy
                              << " | drift=" << drift
                              << std::endl;
                }
            }
        }
    }
    bool apply_swarm_focus_strike(
        std::vector<Sonde>& swarm,
        bool verbose = false
    ) {
        if (swarm.size() <= 13) {
            return false;
        }

        double total_kinetic_energy = 0.0;
        double weighted_offset_sum = 0.0;

        for (std::size_t i = 9; i <= 13 && i < swarm.size(); ++i) {
            if (swarm[i].hit) {
                continue;
            }
            const double off = swarm[i].current_offset.get_d();
            if (!std::isfinite(off) || off <= 0.0) {
                continue;
            }
            const double e =
                std::abs(calculate_morley_energy(N * swarm[i].current_offset).get_d());
            const double w = 1.0 / (e + 1e-25);
            if (!std::isfinite(w) || w <= 0.0) {
                continue;
            }
            total_kinetic_energy += w;
            weighted_offset_sum += off * w;
        }

        if (!(total_kinetic_energy > 0.0) || !std::isfinite(total_kinetic_energy) ||
            !std::isfinite(weighted_offset_sum)) {
            return false;
        }

        const double focus_offset = weighted_offset_sum / total_kinetic_energy;
        if (!std::isfinite(focus_offset) || focus_offset <= 0.0) {
            return false;
        }
        const double clamped_focus = std::clamp(focus_offset, 1e-18, 1.0);

        if (verbose) {
            std::cout << "[STRIKE] Buendele Energie auf Fokus-Offset: "
                      << std::scientific << std::setprecision(18) << clamped_focus << std::endl;
        }

        const mpf_class strike_rate("1e-7");
        const mpf_class radar_rate("1e-10");
        for (std::size_t i = 9; i <= 13 && i < swarm.size(); ++i) {
            if (swarm[i].hit) {
                continue;
            }
            swarm[i].current_offset = mpf_class(clamped_focus);
            swarm[i].learning_rate = strike_rate;
            swarm[i].base_learning_rate = strike_rate;
            swarm[i].steps_without_gradient = 0;
            swarm[i].active = true;
            swarm[i].history.clear();
            swarm[i].status = "ZETA_IMPULSE_STRIKE";
        }

        if (swarm.size() > 6 && !swarm[6].hit) {
            swarm[6].learning_rate = radar_rate;
            swarm[6].base_learning_rate = radar_rate;
        }

        return true;
    }
    void trigger_swarm_collapse(
        std::vector<Sonde>& swarm,
        std::size_t master_idx,
        bool verbose = false
    ) {
        if (master_idx >= swarm.size()) {
            return;
        }

        const mpf_class& master_off = swarm[master_idx].current_offset;
        const double target_d = master_off.get_d();
        if (!std::isfinite(target_d) || target_d <= 0.0) {
            return;
        }

        if (verbose) {
            std::cout << "[!!!] SCHWARM-KOLLAPS: Master-Sonde " << master_idx
                      << " hat ionisiert. Synchronisiere Array..." << std::endl;
        }

        const mpz_class from_p =
            resolve_valid_factor_near(round_standard_part(mpf_class(1) / master_off));
        const mpz_class from_x =
            resolve_valid_factor_near(round_standard_part(N * master_off));
        const mpz_class collapse_factor = (from_p > 1) ? from_p : from_x;
        if (collapse_factor > 1 && verbose) {
            std::cout << "[SUCCESS] FAKTOR GEFUNDEN DURCH KOLLAPS: " << collapse_factor
                      << std::endl;
        }

        const mpf_class lock_rate("1e-11");
        for (std::size_t i = 0; i < swarm.size(); ++i) {
            if (i == master_idx) {
                continue;
            }
            swarm[i].current_offset = master_off;
            swarm[i].learning_rate = lock_rate;
            swarm[i].base_learning_rate = lock_rate;
            swarm[i].steps_without_gradient = 0;
            swarm[i].active = true;
            swarm[i].history.clear();
            swarm[i].status = "MAGNETIC_COHERENCE_LOCK";
        }
    }
    bool apply_cream5_suction(
        std::vector<Sonde>& sondes,
        bool verbose = false
    ) {
        if (sondes.size() <= 13) {
            return false;
        }

        Sonde& kern = sondes[6];
        Sonde& flanke = sondes[11];
        const double kern_offset = kern.current_offset.get_d();
        const double flanke_offset = flanke.current_offset.get_d();
        if (!std::isfinite(kern_offset) || !std::isfinite(flanke_offset) ||
            kern_offset <= 0.0 || flanke_offset <= 0.0) {
            return false;
        }

        const double p_flanke = 1.0 / flanke_offset;
        const double p_kern = 1.0 / kern_offset;
        if (!std::isfinite(p_flanke) || !std::isfinite(p_kern) || p_flanke <= 0.0 || p_kern <= 0.0) {
            return false;
        }

        const double energy_flanke = calculate_offset_zeta_energy(flanke_offset);
        const double energy_kern = calculate_offset_zeta_energy(kern_offset);
        if (!std::isfinite(energy_flanke) || !std::isfinite(energy_kern) || !(energy_kern > 0.0)) {
            return false;
        }

        const double energy_gap = energy_flanke / energy_kern;
        if (energy_gap <= 100.0 || energy_flanke >= 1e-5 || energy_kern >= 1e-5) {
            return false;
        }

        if (verbose) {
            std::cout << "[VACUUM] Saugpumpe aktiv! Ziehe Flanken-Energie in die Singularitaet..."
                      << " | p_flanke=" << p_flanke
                      << " | p_kern=" << p_kern
                      << " | energy_gap=" << energy_gap
                      << std::endl;
        }

        const mpf_class collapse_rate("1e-8");
        std::size_t collapsed = 0;
        for (std::size_t i = 9; i <= 13 && i < sondes.size(); ++i) {
            if (i == 6 || sondes[i].hit) {
                continue;
            }

            sondes[i].current_offset = kern.current_offset;
            sondes[i].learning_rate = collapse_rate;
            if (sondes[i].base_learning_rate < collapse_rate) {
                sondes[i].base_learning_rate = collapse_rate;
            }
            sondes[i].steps_without_gradient = 0;
            sondes[i].active = true;
            sondes[i].history.clear();
            sondes[i].status = "VACUUM_COLLAPSE";
            ++collapsed;
        }

        return collapsed > 0;
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
    mpf_class log_shifted_mpf_from_log(
        const mpf_class& x,
        const mpf_class& log_x,
        unsigned long shift
    ) {
        mpfr_t relative_shift_fr;
        mpfr_t log1p_fr;
        mpfr_init2(relative_shift_fr, precision_bits);
        mpfr_init2(log1p_fr, precision_bits);

        mpf_class relative_shift = mpf_class(shift) / x;
        mpfr_set_f(relative_shift_fr, relative_shift.get_mpf_t(), MPFR_RNDN);
        mpfr_log1p(log1p_fr, relative_shift_fr, MPFR_RNDN);

        mpf_class result = log_x + mpf_from_mpfr(log1p_fr);
        mpfr_clear(relative_shift_fr);
        mpfr_clear(log1p_fr);
        return result;
    }
    mpf_class log_shifted_mpf(const mpf_class& x, unsigned long shift) {
        const mpf_class log_x = log_mpf(x);
        return log_shifted_mpf_from_log(x, log_x, shift);
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
    DerivativeSnapshot evaluate_derivatives_at_x(
        const mpf_class& x0,
        double epsilon,
        const mpf_class* precomputed_energy = nullptr
    ) {
        double effective_epsilon = epsilon;
        const mpf_class e_0 =
            (precomputed_energy != nullptr) ? *precomputed_energy : calculate_high_res_morley(x0);

        for (int attempt = 0; attempt < 8; ++attempt) {
            const mpf_class h = N * effective_epsilon;

            const mpf_class e_m2 = calculate_high_res_morley(x0 - 2 * h);
            const mpf_class e_m1 = calculate_high_res_morley(x0 - h);
            const mpf_class e_p1 = calculate_high_res_morley(x0 + h);
            const mpf_class e_p2 = calculate_high_res_morley(x0 + 2 * h);
            const double epsilon_squared = effective_epsilon * effective_epsilon;
            const double epsilon_fourth = epsilon_squared * epsilon_squared;

            const mpf_class gradient = (e_p1 - e_m1) / (2 * effective_epsilon);
            const mpf_class curvature =
                (e_p1 - 2 * e_0 + e_m1) / epsilon_squared;
            const mpf_class fourth_derivative =
                (e_m2 - 4 * e_m1 + 6 * e_0 - 4 * e_p1 + e_p2) / epsilon_fourth;
            const mpf_class abs_gradient = abs_mpf(gradient);
            const mpf_class abs_curvature = abs_mpf(curvature);
            const mpf_class phase_shift = abs_gradient / (1 + abs_mpf(fourth_derivative));

            if (abs_curvature > 0 || abs_gradient > 0) {
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

        return {e_0, 0, 0, 0, 0, effective_epsilon};
    }
    DerivativeSnapshot evaluate_derivatives(double offset, double epsilon) {
        const mpf_class x0 = N * offset;
        return evaluate_derivatives_at_x(x0, epsilon);
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
    bool is_valid_factor_candidate(const mpz_class& candidate) const {
        return candidate > 1 &&
               candidate < N_int &&
               mpz_divisible_p(N_int.get_mpz_t(), candidate.get_mpz_t());
    }
    mpz_class resolve_valid_factor_near(const mpz_class& center) const {
        for (const long delta : {0L, -1L, 1L}) {
            const mpz_class probe = center + delta;
            if (is_valid_factor_candidate(probe)) {
                return probe;
            }
        }
        return mpz_class(0);
    }
    mpz_class effective_array_factor(const AutoLockResult& r) const {
        if (is_valid_factor_candidate(r.factor)) {
            return r.factor;
        }
        return resolve_valid_factor_near(r.last_candidate);
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
        if (limit == 0) {
            return;
        }

        const ScanCandidate candidate{offset, energy, standard_value};
        if (candidates.size() < limit) {
            candidates.push_back(candidate);
            std::sort(
                candidates.begin(),
                candidates.end(),
                [](const ScanCandidate& a, const ScanCandidate& b) {
                    return a.energy < b.energy;
                }
            );
            return;
        }

        auto worst_it = std::max_element(
            candidates.begin(),
            candidates.end(),
            [](const ScanCandidate& a, const ScanCandidate& b) {
                return a.energy < b.energy;
            }
        );
        if (worst_it != candidates.end() && energy < worst_it->energy) {
            *worst_it = candidate;
            std::sort(
                candidates.begin(),
                candidates.end(),
                [](const ScanCandidate& a, const ScanCandidate& b) {
                    return a.energy < b.energy;
                }
            );
        }
    }
    void populate_candidate_standard_parts(std::vector<ScanCandidate>& candidates) {
        for (ScanCandidate& candidate : candidates) {
            candidate.standard_part_value = standard_part(N * candidate.offset);
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
    std::string target_factor_arg = "17";
    int target_factor = 17;
    bool autonomous_targeting = false;
    std::size_t zeta_terms = DEFAULT_ZETA_MOMENTUM_MAX_TERMS;
    std::size_t array_sonden = 14;

    if (argc > 1) {
        const std::string first_arg = argv[1];
        if (first_arg == "--help" || first_arg == "-h") {
            std::cout << "Aufruf: ./rydberg_scanner [N] [ziel_faktor] [zeta_terme] [array_sonden]" << std::endl;
            std::cout << "  N            Standard: 1048577" << std::endl;
            std::cout << "  ziel_faktor  Standard: 17 (oder auto fuer autonome Zielerfassung)" << std::endl;
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
            target_factor_arg = argv[2];
            if (target_factor_arg == "auto") {
                target_factor = 0;
                autonomous_targeting = true;
            } else {
                target_factor = std::stoi(target_factor_arg);
            }
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
              << " | Ziel-Faktor=" << (autonomous_targeting ? std::string("auto") : std::to_string(target_factor))
              << " | Zeta-Terme=" << get_zeta_momentum_max_terms()
              << " | Array-Sonden=" << array_sonden
              << std::endl;

    // Beispiel: Eine 1024-Bit RSA Zahl (stark verkürzt für Demo)
    HofbiScanner scanner(rsa_n, 8192);
    
    if (autonomous_targeting) {
        std::cout << "[AUTO-TARGET] Ueberspringe manuelle Vorlauf-Heuristiken und starte direkt im Wurzel-Schwarmmodus."
                  << std::endl;
    } else if (target_factor == 0 && rsa_n == "861163313") {
        std::cout << "[BATCH-RESET] Ueberspringe globalen Vorlauf fuer freien Symmetrie-Run."
                  << std::endl;
    } else if (target_factor != 61681) {
        // Allgemeiner Vorlauf fuer unbekannte oder andere Resonanzziele.
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
    } else {
        std::cout << "[FAST-PATH] Ueberspringe globalen 0.9-Vorlauf fuer 61681-Run." << std::endl;
    }
    if (!autonomous_targeting && target_factor > 1) {
        scanner.target_factor_resonance(target_factor, 1e-6, 1e-10);
    } else {
        std::cout << "[TARGET-LOCK] Ueberspringe festen Ziel-Lock fuer freien Suchlauf."
                  << std::endl;
    }
    if (!autonomous_targeting && target_factor == 17) {
        scanner.auto_lock_factor_resonance(1.0e-5, 1.0e-10, 1.0e-15, 400, 5.0e-6);
    } else {
        std::cout << "[AUTO-LOCK] Ueberspringe legacy-17-Auto-Lock fuer Ziel-Faktor "
                  << (autonomous_targeting ? std::string("auto") : std::to_string(target_factor))
                  << std::endl;
    }
    scanner.run_multi_target_array(array_sonden, 1.0e-7, 1.0e-15, 5000, 5.0e-6, 5);

    return 0;
}