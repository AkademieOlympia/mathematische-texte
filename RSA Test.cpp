#include <iostream>
#include <fstream>
#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <array>

// Konfiguration der Präzision: 2048 Bit für die Geometrie-Analyse
// Das entspricht etwa 616 Dezimalstellen Genauigkeit.
const int PRECISION_BITS = 2048;

struct GDA_Report {
    mpf_class n_val;
    mpf_class morley_span;
    mpf_class resonance;
    unsigned long bases_passed;
    unsigned long bases_total;
    std::string label;
};

struct GDA_Result {
    std::string typ;
    int bits;
    mpf_class ms_span;
    double resonanz_prozent;
    mpf_class rel_fehler;
    std::string urteil;
    std::string faktor_p;
    std::string faktor_q;
};

struct SemiprimeData {
    mpz_class p;
    mpz_class q;
    mpz_class n;
};

double log_mpf(const mpf_class& value) {
    mp_exp_t exp2 = 0;
    double mantissa = mpf_get_d_2exp(&exp2, value.get_mpf_t());
    return std::log(mantissa) + static_cast<double>(exp2) * std::log(2.0);
}

mpf_class abs_mpf(const mpf_class& value) {
    return value >= 0 ? value : -value;
}

mpf_class berechne_rel_fehler(const mpf_class& span, const mpz_class& N) {
    mpf_class n_f(N);
    return span / log_mpf(n_f);
}

mpf_class log1p_series(const mpf_class& x) {
    if (x == 0) {
        return 0;
    }

    mpf_class sum = x;
    mpf_class term = x;
    for (int k = 2; k <= 8; ++k) {
        term *= x;
        if (k % 2 == 0) {
            sum -= term / k;
        } else {
            sum += term / k;
        }
    }
    return sum;
}

mpz_class generate_prime_with_bits(gmp_randstate_t state, unsigned long bits) {
    mpz_class candidate;

    while (true) {
        mpz_urandomb(candidate.get_mpz_t(), state, bits);
        mpz_setbit(candidate.get_mpz_t(), bits - 1);
        mpz_setbit(candidate.get_mpz_t(), 0);
        mpz_nextprime(candidate.get_mpz_t(), candidate.get_mpz_t());

        if (mpz_sizeinbase(candidate.get_mpz_t(), 2) == bits) {
            return candidate;
        }
    }
}

SemiprimeData generate_semiprime_with_bits(gmp_randstate_t state, unsigned long bits) {
    const unsigned long factor_bits = bits / 2;

    while (true) {
        mpz_class p = generate_prime_with_bits(state, factor_bits);
        mpz_class q = generate_prime_with_bits(state, factor_bits);
        mpz_class product = p * q;

        if (mpz_sizeinbase(product.get_mpz_t(), 2) == bits) {
            return {p, q, product};
        }
    }
}

bool is_strong_probable_prime_base(const mpz_class& n, unsigned long base_value) {
    if (n < 2) {
        return false;
    }
    if (n == 2 || n == 3) {
        return true;
    }
    if (mpz_even_p(n.get_mpz_t())) {
        return false;
    }

    mpz_class base = base_value % n;
    if (base == 0) {
        return true;
    }

    mpz_class d = n - 1;
    unsigned long s = 0;
    while (mpz_even_p(d.get_mpz_t())) {
        d /= 2;
        ++s;
    }

    mpz_class x;
    mpz_powm(x.get_mpz_t(), base.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());
    if (x == 1 || x == n - 1) {
        return true;
    }

    for (unsigned long r = 1; r < s; ++r) {
        mpz_powm_ui(x.get_mpz_t(), x.get_mpz_t(), 2, n.get_mpz_t());
        if (x == n - 1) {
            return true;
        }
    }

    return false;
}

GDA_Report analyze_geometry(mpz_class N, std::string name) {
    mpf_set_default_prec(PRECISION_BITS);
    
    mpf_class n_f(N);
    mpf_class log_n = log_mpf(n_f); // Logarithmische Skalierung des Raums

    // Definition der ptolemäischen Testpunkte (Vierlings-Struktur)
    // Wir untersuchen die Umgebung von N auf harmonische Einrastung
    unsigned long offsets[] = {0, 2, 6, 8};
    std::vector<mpf_class> points;
    for (unsigned long off : offsets) {
        mpf_class relative_offset = mpf_class(off) / n_f;
        points.push_back(log_n + log1p_series(relative_offset));
    }

    // Wir betrachten die lokalen Log-Abstände und ihre zweite Differenz.
    // Diese Größe ist nicht identisch null und misst eine diskrete Krümmung.
    mpf_class d12 = points[1] - points[0];
    mpf_class d23 = points[2] - points[1];
    mpf_class d34 = points[3] - points[2];
    
    // Die Krümmungsspannung wird aus der zweiten Differenz gebildet.
    mpf_class span = abs_mpf(d34 - 2 * d23 + d12);

    // Arithmetische Resonanz: Anteil kleiner Basen, die den starken
    // Miller-Rabin-Test bestehen. Primzahlen liegen nahe 1, Komposita deutlich darunter.
    const std::array<unsigned long, 12> bases = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    unsigned long passed = 0;

    for (unsigned long base : bases) {
        if (is_strong_probable_prime_base(N, base)) {
            ++passed;
        }
    }

    mpf_class resonance = mpf_class(passed) / bases.size();
    return {n_f, span, resonance, passed, static_cast<unsigned long>(bases.size()), name};
}

std::string bestimme_urteil(const GDA_Report& rep) {
    if (rep.morley_span == 0) {
        return "EXAKT FLACH";
    }
    if (rep.resonance > 9e-1) {
        return "PRIMVERTRAEGLICH";
    }
    if (rep.resonance > 2e-1) {
        return "ARITHMETISCH GRENZFALL";
    }
    return "KOMPOSIT-VERDACHT";
}

GDA_Result to_result(
    const GDA_Report& rep,
    const mpz_class& N,
    const std::string& faktor_p = "",
    const std::string& faktor_q = ""
) {
    return {
        rep.label,
        static_cast<int>(mpz_sizeinbase(N.get_mpz_t(), 2)),
        rep.morley_span,
        rep.resonance.get_d() * 100.0,
        berechne_rel_fehler(rep.morley_span, N),
        bestimme_urteil(rep),
        faktor_p,
        faktor_q
    };
}

void write_csv(const std::vector<GDA_Result>& results) {
    std::ofstream file("gda_results.csv");
    if (!file) {
        std::cerr << ">> Fehler: 'gda_results.csv' konnte nicht geschrieben werden.\n";
        return;
    }

    file << std::scientific << std::setprecision(30);
    file << "Typ;Bits;MS-Span;Resonanz_Prozent;Rel_Fehler;Urteil;Faktor_P;Faktor_Q\n";

    for (const auto& r : results) {
        file << r.typ << ";"
             << r.bits << ";"
             << r.ms_span << ";"
             << r.resonanz_prozent << ";"
             << r.rel_fehler << ";"
             << r.urteil << ";"
             << r.faktor_p << ";"
             << r.faktor_q << "\n";
    }

    std::cout << ">> Ergebnisse in 'gda_results.csv' gespeichert.\n";
}

void print_report(
    const GDA_Report& rep,
    const std::string& faktor_p = "",
    const std::string& faktor_q = ""
) {
    double resonance_percent = rep.resonance.get_d() * 100.0;
    std::string urteil = bestimme_urteil(rep);

    std::cout << "--- ANALYSE: " << rep.label << " ---" << std::endl;
    // Wissenschaftliche Notation macht sehr kleine Unterschiede sichtbar.
    std::cout << std::scientific << std::setprecision(30);
    std::cout << "MS-Span (Spannung): " << rep.morley_span << std::endl;
    std::cout << "Arithm. Resonanz:   " << rep.resonance << std::endl;
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "Bestandene Basen:   " << rep.bases_passed << "/" << rep.bases_total << std::endl;
    std::cout << "Resonanz in %:      " << resonance_percent << "%" << std::endl;
    if (!faktor_p.empty() && !faktor_q.empty()) {
        std::cout << "Faktor p (" << faktor_p.size() << " Ziffern): " << faktor_p << std::endl;
        std::cout << "Faktor q (" << faktor_q.size() << " Ziffern): " << faktor_q << std::endl;
    }
    std::cout << "URTEIL: " << urteil << std::endl;
    std::cout << std::endl;
}

int main() {
    mpf_set_default_prec(PRECISION_BITS);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));
    std::vector<GDA_Result> all_results;

    std::cout << "GDA TESTBED v1.0 - HOFFBAUER GEOMETRIE" << std::endl;
    std::cout << "========================================" << std::endl << std::endl;

    // 1. Erzeugung einer exakten 1024-Bit Semiprimzahl
    SemiprimeData semi_data = generate_semiprime_with_bits(state, 1024);
    mpz_class n_semi = semi_data.n;

    // 2. Erzeugung einer echten 1024-Bit Primzahl
    mpz_class n_prime = generate_prime_with_bits(state, 1024);

    // 3. Zusätzliche Vergleichszahlen für die Einordnung des Basistests
    mpz_class n_carmichael = 561;
    mpz_class n_pseudoprime = 3215031751UL;

    // 4. Durchführung der Geometrischen Analyse
    GDA_Report report_prime = analyze_geometry(n_prime, "ECHTE PRIMZAHL (1024-Bit)");
    GDA_Report report_semi  = analyze_geometry(n_semi,  "SEMIPRIMZAHL (p*q, 1024-Bit)");
    GDA_Report report_carmichael = analyze_geometry(n_carmichael, "CARMICHAEL-ZAHL 561");
    GDA_Report report_pseudoprime = analyze_geometry(n_pseudoprime, "STARKE PSEUDOPRIMZAHL 3215031751");

    // 5. Ergebnisse ausgeben
    print_report(report_prime);
    print_report(report_semi, semi_data.p.get_str(), semi_data.q.get_str());
    print_report(report_carmichael);
    print_report(report_pseudoprime);

    all_results.push_back(to_result(report_prime, n_prime));
    all_results.push_back(to_result(report_semi, n_semi, semi_data.p.get_str(), semi_data.q.get_str()));
    all_results.push_back(to_result(report_carmichael, n_carmichael));
    all_results.push_back(to_result(report_pseudoprime, n_pseudoprime));

    write_csv(all_results);

    return 0;
}