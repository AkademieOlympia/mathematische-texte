#include <iostream>
#include <vector>
#include <string>
#include <cctype>
#include <cmath>
#include <complex>
#include <numeric>
#include <chrono>
#include <algorithm>
#include <iomanip>

// ==============================================================================
// TEILMODELL KOSMOLOGIE: C++ SCHWESTER-IMPLEMENTIERUNG (GROVER + EABC E8 TOE)
// ==============================================================================

using Complex = std::complex<double>;

struct CosmologyTopic {
    int id;
    std::string topic;
    std::vector<std::string> keywords;
    std::string summary;
    std::vector<double> vector_8d;
};

// --- 1. KOSMOLOGIE WISSENSBASIS ---
const std::vector<CosmologyTopic> COSMOLOGY_KB = {
    {0, "Urknall & Kosmische Inflation", {"urknall", "inflation", "flrw", "inflaton"}, 
     "Expansion des Universums aus einer Inflaton-getriebenen Phase in die FLRW-Metrik.",
     {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}},
    {1, "Dunkle Materie & Halos", {"dunkle materie", "cdm", "wimp", "axion", "halo"}, 
     "Nicht-baryonische kalte Materie zur Erklärung galaktischer Rotationskurven.",
     {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
    {2, "Dunkle Energie & Lambda-CDM", {"dunkle energie", "lambda", "vakuumenergie", "beschleunigung"}, 
     "Kosmologische Konstante Lambda treibt die beschleunigte Expansion der Raumzeit.",
     {1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
    {3, "Kosmischer Mikrowellenhintergrund (CMB)", {"cmb", "planck", "anisotropie", "strahlungszeitalter", "schwarzkörper", "schwarzer körper", "schwarzkoerper", "blackbody"}, 
     "Relikte der Rekombinations-Epoche bei z ≈ 1100, nahezu perfektes Schwarzkörper-Spektrum mit Anisotropien im microKelvin-Bereich.",
     {0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5}},
    {4, "E8-Kosmologie & Spin(8)-Trialität", {"e8 kosmologie", "e8-kosmologie", "e8", "trialitaet", "trialität", "fano", "oktonionen", "spin(8)", "vereinheitlichung"}, 
     "Einbettung der 4 Fundamentalkräfte und der kosmischen Geometrie im E8-Lie-Gitter.",
     {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0}},
    {5, "Penrose-Hawking Singularitäten & Quantengravitation", {"schwarzes loch", "hawking", "holographie", "quantengravitation"}, 
     "Thermodynamik schwarzer Löcher und holographische Informations-Erhaltung.",
     {0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0}},
    {6, "Hubble-Tension & Großräumige Strukturen", {"hubble", "h0", "rotverschiebung", "spannung"}, 
     "Diskrepanz zwischen lokaler H0-Messung (LMC/Supernovae) und CMB-Frühzeitmodellen.",
     {0.5, 0.5, -0.5, -0.5, 0.5, 0.5, -0.5, -0.5}},
    {7, "Baryonische Akustische Oszillationen (BAO)", {"bao", "schallhorizont", "strukturformation"}, 
     "Gefrorene Dichtewellen der baryonisch-photonischen Plasmaphase als kosmischer Maßstab.",
     {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0}}
};

// --- 2. GROVER QUANTUM SEARCH ENGINE ---
class GroverSearchEngine {
private:
    int N;
public:
    GroverSearchEngine(int size) : N(size) {}

    std::vector<double> runGrover(const std::vector<int>& targets, int iterations) {
        std::vector<double> state(N, 1.0 / std::sqrt(N));
        std::vector<char> marked(static_cast<size_t>(N), 0);
        for (int target : targets) {
            if (target >= 0 && target < N) marked[static_cast<size_t>(target)] = 1;
        }

        for (int it = 0; it < iterations; ++it) {
            for (int i = 0; i < N; ++i) {
                if (marked[static_cast<size_t>(i)]) state[static_cast<size_t>(i)] *= -1.0;
            }

            double mean = 0.0;
            for (double val : state) mean += val;
            mean /= N;

            for (int i = 0; i < N; ++i) {
                state[static_cast<size_t>(i)] = 2.0 * mean - state[static_cast<size_t>(i)];
            }
        }

        std::vector<double> probs(N);
        for (int i = 0; i < N; ++i) {
            probs[i] = state[i] * state[i];
        }
        return probs;
    }
};

// --- 3. EABC TOE E8-FANO ENGINE ---
class EABCTOEEngine {
public:
    std::vector<Complex> encodePhaseE8(const std::vector<double>& vec8d) {
        double norm = 0.0;
        for (double v : vec8d) norm += v * v;
        norm = std::sqrt(norm);

        std::vector<Complex> state(8);
        for (size_t i = 0; i < 8; ++i) {
            double angle = M_PI * (norm > 0 ? vec8d[i] / norm : 0.0);
            state[i] = Complex(std::cos(angle / 2.0), std::sin(angle / 2.0));
        }
        
        double st_norm = 0.0;
        for (const auto& c : state) st_norm += std::norm(c);
        st_norm = std::sqrt(st_norm);
        for (size_t i = 0; i < 8; ++i) state[i] /= st_norm;

        return state;
    }

    std::pair<std::vector<Complex>, std::vector<std::vector<Complex>>> applySO8Triality(const std::vector<Complex>& state) {
        double theta = M_PI / 4.0;
        double c = std::cos(theta), s = std::sin(theta);

        std::vector<std::vector<Complex>> U(8, std::vector<Complex>(8, 0.0));
        for (int i = 0; i < 8; ++i) U[i][i] = 1.0;

        U[0][0] = c; U[0][1] = -s;
        U[1][0] = s; U[1][1] = c;
        U[2][2] = c; U[2][3] = -s;
        U[3][2] = s; U[3][3] = c;

        std::vector<Complex> transformed(8, 0.0);
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                transformed[i] += U[i][j] * state[j];
            }
        }
        return {transformed, U};
    }

    double invertAndVerifyFidelity(const std::vector<Complex>& orig, 
                                   const std::vector<Complex>& trans, 
                                   const std::vector<std::vector<Complex>>& U) {
        // Inversion: U^dagger * trans
        std::vector<Complex> inv(8, 0.0);
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                inv[i] += std::conj(U[j][i]) * trans[j];
            }
        }

        Complex inner_prod = 0.0;
        for (int i = 0; i < 8; ++i) {
            inner_prod += std::conj(orig[i]) * inv[i];
        }

        return std::norm(inner_prod) * 100.0;
    }
};

// setw zählt Bytes. Deutsche Umlaute sind zwei Bytes und ein sichtbares Zeichen.
int sichtbare_spalten(const std::string& text) {
    int spalten = 0;
    for (unsigned char zeichen : text) {
        if ((zeichen & 0xC0) != 0x80) ++spalten;
    }
    return spalten;
}

std::string rechts_auffuellen(const std::string& text, int breite) {
    int fehlend = breite - sichtbare_spalten(text);
    if (fehlend <= 0) return text;
    return text + std::string(static_cast<size_t>(fehlend), ' ');
}

struct ThemenTreffer {
    int id = -1;
    int score = 0;
    int span = 0;
    std::vector<std::string> keywords;
};

std::string ascii_klein(std::string text) {
    for (char& zeichen : text) {
        auto wert = static_cast<unsigned char>(zeichen);
        if (wert >= 'A' && wert <= 'Z') zeichen = static_cast<char>(wert - 'A' + 'a');
    }
    return text;
}

std::string beschnitten(const std::string& text) {
    size_t anfang = 0;
    while (anfang < text.size() && std::isspace(static_cast<unsigned char>(text[anfang]))) ++anfang;
    size_t ende = text.size();
    while (ende > anfang && std::isspace(static_cast<unsigned char>(text[ende - 1]))) --ende;
    return text.substr(anfang, ende - anfang);
}

std::vector<std::string> teilfragen(const std::string& query) {
    std::vector<std::string> teile;
    std::string aktuell;
    for (char zeichen : query) {
        if (zeichen == ',' || zeichen == ';') {
            std::string stueck = beschnitten(aktuell);
            if (!stueck.empty()) teile.push_back(stueck);
            aktuell.clear();
        } else {
            aktuell.push_back(zeichen);
        }
    }
    std::string stueck = beschnitten(aktuell);
    if (!stueck.empty()) teile.push_back(stueck);
    if (teile.empty()) teile.push_back(query);
    return teile;
}

std::vector<ThemenTreffer> treffer_fuer(const std::string& text) {
    std::string q = ascii_klein(text);
    std::vector<ThemenTreffer> treffer;
    for (const auto& item : COSMOLOGY_KB) {
        ThemenTreffer hit;
        hit.id = item.id;
        for (const auto& kw : item.keywords) {
            if (q.find(kw) != std::string::npos) {
                hit.score += 2;
                hit.span += static_cast<int>(kw.size());
                hit.keywords.push_back(kw);
            }
        }
        if (hit.score > 0) treffer.push_back(hit);
    }
    std::sort(treffer.begin(), treffer.end(), [](const ThemenTreffer& a, const ThemenTreffer& b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.span != b.span) return a.span > b.span;
        return a.id < b.id;
    });
    return treffer;
}

std::string schluessel_liste(const std::vector<std::string>& keywords) {
    std::string text;
    for (size_t i = 0; i < keywords.size(); ++i) {
        if (i > 0) text += ", ";
        text += keywords[i];
    }
    return text;
}

int grover_iterationen(int markiert, int n) {
    if (markiert <= 0 || markiert >= n) return 0;
    double theta = std::asin(std::sqrt(static_cast<double>(markiert) / static_cast<double>(n)));
    int k = static_cast<int>(std::lround(M_PI / (4.0 * theta) - 0.5));
    return std::max(0, k);
}

void phase_ausgeben(EABCTOEEngine& toe_engine, const std::vector<double>& vec8d, const std::string& einzug) {
    auto q_state = toe_engine.encodePhaseE8(vec8d);
    auto [trans_state, U_matrix] = toe_engine.applySO8Triality(q_state);
    double fidelity = toe_engine.invertAndVerifyFidelity(q_state, trans_state, U_matrix);
    double vec_norm = 0.0;
    for (double v : vec8d) vec_norm += v * v;
    vec_norm = std::sqrt(vec_norm);
    double phase_loss = std::max(0.0, 1.0 - fidelity / 100.0);
    std::cout << einzug << "E8 Root Vector Norm:         " << std::fixed << std::setprecision(3) << vec_norm << "\n";
    std::cout << einzug << "Phase Reconstruction Loss:    " << std::scientific << std::setprecision(2) << phase_loss
              << std::fixed << " (" << std::setprecision(6) << fidelity << " % Fidelity)\n";
}

// --- 4. MAIN INTERFACE ---
int main(int argc, char* argv[]) {
    std::string query = "E8 Kosmologie und Dunkle Energie";
    if (argc > 1) {
        query = argv[1];
    }

    std::cout << "======================================================================\n";
    std::cout << "  TEILMODELL KOSMOLOGIE: OLLAMA + EABC TOE + GROVER ENGINE (C++20)\n";
    std::cout << "======================================================================\n";
    std::cout << "• Such-Anfrage: '" << query << "'\n\n";

    auto start_time = std::chrono::high_resolution_clock::now();

    auto fragen = teilfragen(query);
    auto alle = treffer_fuer(query);
    if (alle.empty()) {
        ThemenTreffer standard;
        standard.id = 4;
        standard.score = 0;
        alle.push_back(standard);
    }

    std::vector<int> targets;
    std::vector<char> ist_ziel(COSMOLOGY_KB.size(), 0);
    for (const auto& hit : alle) {
        targets.push_back(hit.id);
        ist_ziel[static_cast<size_t>(hit.id)] = 1;
    }

    std::cout << "[STUFE 1: OLLAMA & ORAKEL MATCHING]\n";
    std::cout << "  • Teilfragen:                  " << fragen.size() << "\n";
    for (size_t f = 0; f < fragen.size(); ++f) {
        auto lokale = treffer_fuer(fragen[f]);
        std::cout << "    " << (f + 1) << ". " << fragen[f] << "\n";
        if (lokale.empty()) {
            std::cout << "       → kein Thema in der Wissensbasis\n";
        } else {
            for (const auto& hit : lokale) {
                std::cout << "       → [" << hit.id << "] " << COSMOLOGY_KB[static_cast<size_t>(hit.id)].topic
                          << " (" << schluessel_liste(hit.keywords) << ")\n";
            }
        }
    }
    std::cout << "  • Markierte Zustände:          " << targets.size() << "\n";
    for (const auto& hit : alle) {
        const auto& item = COSMOLOGY_KB[static_cast<size_t>(hit.id)];
        std::cout << "    • [" << item.id << "] " << item.topic << "\n";
        std::cout << "      " << item.summary << "\n";
    }
    std::cout << "\n";

    int num_iterations = grover_iterationen(static_cast<int>(targets.size()), 8);
    GroverSearchEngine grover(8);
    auto probs = grover.runGrover(targets, num_iterations);
    double ziel_summe = 0.0;
    for (int id : targets) ziel_summe += probs[static_cast<size_t>(id)];

    std::cout << "[STUFE 2: GROVER AMPLITUDEN-VERSTÄRKUNG (C++ ENGINE)]\n";
    std::cout << "  • Grover Iterationen:          " << num_iterations << "\n";
    std::cout << "  • Ziel-Wahrscheinlichkeit:     " << std::fixed << std::setprecision(2) << ziel_summe * 100.0 << " %\n";
    std::cout << "  • Wahrscheinlichkeitsverteilung über alle 8 Kosmologie-Zustände:\n";

    for (int i = 0; i < 8; ++i) {
        int bar_len = static_cast<int>(probs[static_cast<size_t>(i)] * 30);
        std::string bar(static_cast<size_t>(bar_len), '#');
        std::string marker = ist_ziel[static_cast<size_t>(i)] ? " <--- [ZIEL]" : "";
        std::cout << "    State |" << i << "> " << rechts_auffuellen(COSMOLOGY_KB[static_cast<size_t>(i)].topic, 52)
                  << std::right << std::setw(5) << std::fixed << std::setprecision(1)
                  << probs[static_cast<size_t>(i)] * 100.0 << "% " << bar << marker << "\n";
    }

    EABCTOEEngine toe_engine;
    std::cout << "\n[STUFE 3: EABC TOE E8-FANO & PHASE INVERSION (C++)]\n";
    std::vector<double> summe(8, 0.0);
    for (int id : targets) {
        const auto& item = COSMOLOGY_KB[static_cast<size_t>(id)];
        std::cout << "  • [" << item.id << "] " << item.topic << "\n";
        phase_ausgeben(toe_engine, item.vector_8d, "    ");
        for (size_t k = 0; k < summe.size(); ++k) summe[k] += item.vector_8d[k];
    }
    if (targets.size() > 1) {
        std::cout << "  • Kombinierter Root Vector (Summe der markierten Themen)\n";
        phase_ausgeben(toe_engine, summe, "    ");
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed_us = std::chrono::duration<double, std::micro>(end_time - start_time).count();
    std::cout << "  • Gesamtlaufzeit C++ Pipeline:  " << std::fixed << std::setprecision(2) << elapsed_us << " µs (" << elapsed_us / 1000.0 << " ms)\n";
    std::cout << "======================================================================\n";

    return 0;
}
