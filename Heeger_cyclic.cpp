#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using std::array;
using std::cout;
using std::endl;
using std::size_t;
using std::vector;

// ============================================================
// Kleine Statistik-Helfer
// ============================================================

struct Stat {
    long long n = 0;
    double sum = 0.0;
    double sum2 = 0.0;

    void add(double x) {
        n++;
        sum += x;
        sum2 += x * x;
    }

    double mean() const {
        return n > 0 ? sum / n : 0.0;
    }

    double var() const {
        if (n <= 1) return 0.0;
        double m = mean();
        return (sum2 - n * m * m) / (n - 1);
    }

    double sd() const {
        double v = var();
        return v > 0.0 ? std::sqrt(v) : 0.0;
    }
};

struct OpenTransition {
    double D;
    int dphi;
};

struct DriftResult {
    double d16;
    double d25;
    double d34;
};

DriftResult compute_drifts_from_labels(
    const vector<double>& Dvals,
    const vector<int>& labels
) {
    array<Stat, 7> st{};
    for (size_t i = 0; i < Dvals.size(); ++i) {
        st[(size_t)labels[i]].add(Dvals[i]);
    }
    return {
        st[1].mean() - st[6].mean(),
        st[2].mean() - st[5].mean(),
        st[3].mean() - st[4].mean()
    };
}

void run_phase_shuffle_test(
    const vector<OpenTransition>& open,
    int B = 1000,
    uint64_t seed = 1234567ULL
) {
    vector<double> Dvals;
    vector<int> labels;
    Dvals.reserve(open.size());
    labels.reserve(open.size());
    for (const auto& x : open) {
        Dvals.push_back(x.D);
        labels.push_back(x.dphi);
    }

    DriftResult obs = compute_drifts_from_labels(Dvals, labels);

    std::mt19937_64 rng(seed);
    int c16 = 0, c25 = 0, c34 = 0;
    for (int b = 0; b < B; ++b) {
        std::shuffle(labels.begin(), labels.end(), rng);
        DriftResult sh = compute_drifts_from_labels(Dvals, labels);
        if (std::fabs(sh.d16) >= std::fabs(obs.d16)) c16++;
        if (std::fabs(sh.d25) >= std::fabs(obs.d25)) c25++;
        if (std::fabs(sh.d34) >= std::fabs(obs.d34)) c34++;
    }
    double p16 = (c16 + 1.0) / (B + 1.0);
    double p25 = (c25 + 1.0) / (B + 1.0);
    double p34 = (c34 + 1.0) / (B + 1.0);

    cout << "\n--- Phase-label shuffle test ---\n";
    cout << "B=" << B << "\n";
    cout << "obs d16=" << obs.d16 << " shuffle_p=" << p16 << "\n";
    cout << "obs d25=" << obs.d25 << " shuffle_p=" << p25 << "\n";
    cout << "obs d34=" << obs.d34 << " shuffle_p=" << p34 << "\n";
}

double normal_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double welch_z_pvalue(const Stat& a, const Stat& b) {
    if (a.n < 2 || b.n < 2) return 1.0;
    double se2 = a.var() / a.n + b.var() / b.n;
    if (se2 <= 0.0) return 1.0;
    double z = (a.mean() - b.mean()) / std::sqrt(se2);
    double p = 2.0 * (1.0 - normal_cdf(std::fabs(z)));
    return p;
}

// ============================================================
// Primzahlsieb
// ============================================================

vector<bool> sieve_prime(size_t N) {
    vector<bool> is_prime(N + 1, true);
    if (N >= 0) is_prime[0] = false;
    if (N >= 1) is_prime[1] = false;

    for (size_t p = 2; p * p <= N; ++p) {
        if (is_prime[p]) {
            for (size_t q = p * p; q <= N; q += p) {
                is_prime[q] = false;
            }
        }
    }
    return is_prime;
}

// ============================================================
// EABC-Tetraeder-Einbettung
// E: 1 mod 12, A: 5 mod 12, B: 7 mod 12, C: 11 mod 12
// ============================================================

struct Vec4 {
    double x0, x1, x2, x3;
};

array<double, 3> family_vec(int pmod12) {
    const double s = 1.0 / std::sqrt(3.0);

    if (pmod12 == 1)  return { s,  s,  s};   // E
    if (pmod12 == 5)  return { s, -s, -s};   // A
    if (pmod12 == 7)  return {-s,  s, -s};   // B
    if (pmod12 == 11) return {-s, -s,  s};   // C

    return {0.0, 0.0, 0.0};
}

Vec4 embed_anchor(long long p) {
    int r12 = (int)(p % 12);
    if (r12 < 0) r12 += 12;
    auto v = family_vec(r12);
    return {std::log((double)p), v[0], v[1], v[2]};
}

double dist4(const Vec4& a, const Vec4& b) {
    double d0 = a.x0 - b.x0;
    double d1 = a.x1 - b.x1;
    double d2 = a.x2 - b.x2;
    double d3 = a.x3 - b.x3;
    return std::sqrt(d0*d0 + d1*d1 + d2*d2 + d3*d3);
}

// Normierter Ptolemäus-Defekt für vier Punkte P0,P1,P2,P3.
// Seiten: a=P0P1, b=P1P2, c=P2P3, d=P3P0.
// Diagonalen: e=P0P2, f=P1P3.
double ptolemy_defect(const Vec4& P0, const Vec4& P1, const Vec4& P2, const Vec4& P3) {
    double a = dist4(P0, P1);
    double b = dist4(P1, P2);
    double c = dist4(P2, P3);
    double d = dist4(P3, P0);
    double e = dist4(P0, P2);
    double f = dist4(P1, P3);

    double lhs = e * f;
    double rhs = a * c + b * d;
    double denom = std::fabs(lhs) + std::fabs(rhs);

    if (denom <= 1e-30) return 0.0;
    return std::fabs(lhs - rhs) / denom;
}

// ============================================================
// Quaternionen und Linkregular-Darstellung
// q = a + b i + c j + d k
// ============================================================

struct Quat {
    double a, b, c, d;

    Quat(double A=0, double B=0, double C=0, double D=0)
        : a(A), b(B), c(C), d(D) {}

    Quat conj() const {
        return Quat(a, -b, -c, -d);
    }

    double reduced_norm() const {
        return a*a + b*b + c*c + d*d;
    }
};

Quat operator+(const Quat& x, const Quat& y) {
    return Quat(x.a+y.a, x.b+y.b, x.c+y.c, x.d+y.d);
}

Quat operator-(const Quat& x, const Quat& y) {
    return Quat(x.a-y.a, x.b-y.b, x.c-y.c, x.d-y.d);
}

Quat operator*(const Quat& x, const Quat& y) {
    return Quat(
        x.a*y.a - x.b*y.b - x.c*y.c - x.d*y.d,
        x.a*y.b + x.b*y.a + x.c*y.d - x.d*y.c,
        x.a*y.c - x.b*y.d + x.c*y.a + x.d*y.b,
        x.a*y.d + x.b*y.c - x.c*y.b + x.d*y.a
    );
}

// 4x4 Linkregular-Matrix für Quaternion q.
// Multiplikation links mit q.
array<array<double,4>,4> left_regular_4(const Quat& q) {
    double a = q.a, b = q.b, c = q.c, d = q.d;
    return {{
        {{ a, -b, -c, -d }},
        {{ b,  a, -d,  c }},
        {{ c,  d,  a, -b }},
        {{ d, -c,  b,  a }}
    }};
}

using Mat16 = array<array<double,16>,16>;
using Mat4Q = array<array<Quat,4>,4>;

Mat16 to_real_16x16(const Mat4Q& M) {
    Mat16 B{};
    for (int bi = 0; bi < 4; ++bi) {
        for (int bj = 0; bj < 4; ++bj) {
            auto block = left_regular_4(M[bi][bj]);
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    B[4*bi + i][4*bj + j] = block[i][j];
                }
            }
        }
    }
    return B;
}

double frobenius_quat(const Mat4Q& M) {
    double s = 0.0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            s += M[i][j].reduced_norm();
    return std::sqrt(s);
}

Mat16 transpose_multiply(const Mat16& B) {
    Mat16 C{};
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            double s = 0.0;
            for (int k = 0; k < 16; ++k) s += B[k][i] * B[k][j];
            C[i][j] = s;
        }
    }
    return C;
}

// Jacobi-Eigenwertverfahren für symmetrische 16x16 Matrix.
// Für diesen Test reicht das völlig.
array<double,16> jacobi_eigenvalues_symmetric(Mat16 A, int max_iter = 2000, double eps = 1e-12) {
    for (int iter = 0; iter < max_iter; ++iter) {
        int p = 0, q = 1;
        double max_off = 0.0;

        for (int i = 0; i < 16; ++i) {
            for (int j = i + 1; j < 16; ++j) {
                double v = std::fabs(A[i][j]);
                if (v > max_off) {
                    max_off = v;
                    p = i;
                    q = j;
                }
            }
        }

        if (max_off < eps) break;

        double app = A[p][p];
        double aqq = A[q][q];
        double apq = A[p][q];

        double tau = (aqq - app) / (2.0 * apq);
        double t = (tau >= 0.0)
            ? 1.0 / (tau + std::sqrt(1.0 + tau*tau))
            : -1.0 / (-tau + std::sqrt(1.0 + tau*tau));
        double c = 1.0 / std::sqrt(1.0 + t*t);
        double s = t * c;

        for (int k = 0; k < 16; ++k) {
            if (k == p || k == q) continue;
            double aik = A[k][p];
            double akq = A[k][q];

            A[k][p] = A[p][k] = c*aik - s*akq;
            A[k][q] = A[q][k] = s*aik + c*akq;
        }

        double new_app = c*c*app - 2.0*s*c*apq + s*s*aqq;
        double new_aqq = s*s*app + 2.0*s*c*apq + c*c*aqq;

        A[p][p] = new_app;
        A[q][q] = new_aqq;
        A[p][q] = A[q][p] = 0.0;
    }

    array<double,16> evals{};
    for (int i = 0; i < 16; ++i) evals[i] = A[i][i];
    return evals;
}

vector<double> singular_values_16(const Mat16& B) {
    Mat16 C = transpose_multiply(B);
    auto evals = jacobi_eigenvalues_symmetric(C);

    vector<double> svals;
    svals.reserve(16);
    for (double x : evals) {
        if (x < 0.0 && x > -1e-8) x = 0.0;
        if (x < 0.0) x = 0.0;
        svals.push_back(std::sqrt(x));
    }

    std::sort(svals.begin(), svals.end(), std::greater<double>());
    return svals;
}

double svd_purity(const vector<double>& svals) {
    double s_sum = 0.0;
    double sq_sum = 0.0;
    for (double s : svals) {
        s_sum += s;
        sq_sum += s*s;
    }
    if (s_sum <= 1e-30) return 0.0;
    return sq_sum / (s_sum * s_sum);
}

// ============================================================
// Heegner/Ramanujan 163-Anker
// C++ double reicht nicht für echte 1e-12-Präzision der Lücke,
// deshalb geben wir hier die bekannte hochpräzise Referenz aus.
// ============================================================

void print_heegner_anchor() {
    cout << "------------------------------------------------------------\n";
    cout << "Heegner-/Ramanujan-Anker\n";
    cout << "e^(pi*sqrt(163)) ist extrem nahe an einer ganzen Zahl.\n";
    cout << "Referenzluecke: Delta_163 ≈ 7.499274028018143e-13\n";
    cout << "Heegner-Kohaerenz H_163 = -log10(Delta_163) ≈ "
         << -std::log10(7.499274028018143e-13) << "\n";
    cout << "Hinweis: Fuer die echte Luecke braucht man MPFR/mpmath/Sage,\n";
    cout << "double reicht hier nicht stabil aus.\n";
    cout << "------------------------------------------------------------\n\n";
}

// ============================================================
// Phase-Quaternion aus Statistiken
// Q_d = S + A i + V j + C k
// S = symmetrischer Defekt
// A = antisymmetrische Drift
// V = mittlere Varianz
// C = log Count Ratio
// ============================================================

Quat phase_quaternion(const array<Stat,7>& phase_stats, int d) {
    int dm = (7 - d) % 7;

    double Dp = phase_stats[d].mean();
    double Dm = phase_stats[dm].mean();

    double Vp = phase_stats[d].var();
    double Vm = phase_stats[dm].var();

    double np = (double)phase_stats[d].n;
    double nm = (double)phase_stats[dm].n;

    double S = 0.5 * (Dp + Dm);
    double A = 0.5 * (Dp - Dm);
    double V = 0.5 * (Vp + Vm);
    double C = std::log((np + 1.0) / (nm + 1.0));

    return Quat(S, A, V, C);
}

Mat4Q build_phase_coherence_matrix(const array<Stat,7>& phase_stats) {
    Quat Q0 = phase_quaternion(phase_stats, 0);
    Quat Q1 = phase_quaternion(phase_stats, 1);
    Quat Q2 = phase_quaternion(phase_stats, 2);
    Quat Q3 = phase_quaternion(phase_stats, 3);

    Mat4Q K{};

    K[0] = { Q0,          Q1,          Q2,          Q3          };
    K[1] = { Q1.conj(),   Q0,          Q3,          Q2          };
    K[2] = { Q2.conj(),   Q3.conj(),   Q0,          Q1          };
    K[3] = { Q3.conj(),   Q2.conj(),   Q1.conj(),   Q0          };

    return K;
}

// ============================================================
// Haupttest
// ============================================================

int main(int argc, char** argv) {
    long long N = 10000000LL;

    // Variable modulo-60 sectors. Default: original 11/41 sector.
    int R1 = 11;
    int R2 = 41;

    if (argc >= 2) {
        N = std::stoll(argv[1]);
    }
    if (argc >= 3) {
        R1 = std::stoi(argv[2]);
    }
    if (argc >= 4) {
        R2 = std::stoi(argv[3]);
    }

    R1 = ((R1 % 60) + 60) % 60;
    R2 = ((R2 % 60) + 60) % 60;

    if (R1 == R2) {
        std::cerr << "Fehler: R1 und R2 muessen verschiedene modulo-60-Sektoren sein.\n";
        return 1;
    }

    int raw_dist = std::abs(R2 - R1);
    int cyclic_dist = std::min(raw_dist, 60 - raw_dist);

    cout << std::fixed << std::setprecision(9);

    print_heegner_anchor();

    cout << "EABC / Onsager / 420-Feinstruktur-Test\n";
    cout << "N = " << N << "\n";
    cout << "Sektor R1/R2 = " << R1 << "/" << R2 << " mod 60\n";
    cout << "Abstand |R2-R1| linear=" << raw_dist
         << " zyklisch_mod_60=" << cyclic_dist << "\n\n";

    if (N < 3) {
        cout << "N muss mindestens 3 sein.\n";
        return 1;
    }
    {
        unsigned long long n2 = (unsigned long long)N + 2uLL;
        if (n2 > (unsigned long long)std::numeric_limits<size_t>::max()) {
            cout << "N zu gross fuer size_t (Sieb).\n";
            return 1;
        }
    }

    auto is_prime = sieve_prime(static_cast<size_t>(N + 2));

    vector<long long> twins;
    twins.reserve((size_t)(N / 100 + 1));

    for (long long p = 7; p + 2 <= N; ++p) {
        if (is_prime[(size_t)p] && is_prime[(size_t)(p + 2)]) {
            twins.push_back(p);
        }
    }

    cout << "Twin-Anker p mit (p,p+2): " << twins.size() << "\n";

    vector<long long> sector;
    sector.reserve(twins.size());

    for (long long p : twins) {
        int r60 = (int)(p % 60);
        if (r60 < 0) r60 += 60;
        if (r60 == R1 || r60 == R2) {
            sector.push_back(p);
        }
    }

    cout << R1 << "/" << R2 << "-Sektor-Anker: " << sector.size() << "\n";

    if (sector.size() < 4) {
        cout << "Zu wenige Daten.\n";
        return 0;
    }

    Stat sR1_R1, sR1_R2, sR2_R1, sR2_R2;
    array<Stat,7> delta_phase_stats{};

    vector<OpenTransition> open;
    open.reserve(sector.size());
    long long open_transitions = 0;

    for (size_t k = 1; k + 2 < sector.size(); ++k) {
        long long p0 = sector[k - 1];
        long long p1 = sector[k];
        long long p2 = sector[k + 1];
        long long p3 = sector[k + 2];

        int from = (int)(p1 % 60);
        int to   = (int)(p2 % 60);
        if (from < 0) from += 60;
        if (to < 0) to += 60;

        Vec4 P0 = embed_anchor(p0);
        Vec4 P1 = embed_anchor(p1);
        Vec4 P2 = embed_anchor(p2);
        Vec4 P3 = embed_anchor(p3);

        double D = ptolemy_defect(P0, P1, P2, P3);

        if (from == R1 && to == R1) sR1_R1.add(D);
        if (from == R1 && to == R2) sR1_R2.add(D);
        if (from == R2 && to == R1) sR2_R1.add(D);
        if (from == R2 && to == R2) sR2_R2.add(D);

        if ((from == R1 && to == R2) || (from == R2 && to == R1)) {
            int phi_from = (int)(p1 % 7);
            int phi_to   = (int)(p2 % 7);
            int dphi = (phi_to - phi_from) % 7;
            if (dphi < 0) dphi += 7;

            delta_phase_stats[dphi].add(D);
            open.push_back({D, dphi});
            open_transitions++;
        }
    }

    auto print_stat = [](const std::string& label, const Stat& s) {
        cout << std::setw(8) << label
             << " n=" << std::setw(8) << s.n
             << " mean=" << std::setw(14) << s.mean()
             << " sd=" << std::setw(14) << s.sd()
             << "\n";
    };

    cout << "\n--- " << R1 << "/" << R2 << "-Ptolemaeus-Matrix ---\n";
    print_stat(std::to_string(R1) + "->" + std::to_string(R1), sR1_R1);
    print_stat(std::to_string(R1) + "->" + std::to_string(R2), sR1_R2);
    print_stat(std::to_string(R2) + "->" + std::to_string(R1), sR2_R1);
    print_stat(std::to_string(R2) + "->" + std::to_string(R2), sR2_R2);

    double delta_ons = sR1_R2.mean() - sR2_R1.mean();
    double scale = 0.5 * (sR1_R2.mean() + sR2_R1.mean());
    double p_ons = welch_z_pvalue(sR1_R2, sR2_R1);
    cout << "\nDelta_Ons = mean("
         << R1 << "->" << R2 << ")-mean("
         << R2 << "->" << R1 << ") = "
         << delta_ons << "\n";
    cout << "Relative Abweichung = " << (scale > 0 ? delta_ons / scale : 0.0) << "\n";
    cout << "Welch-normal-approx p = " << p_ons << "\n";

    cout << "\n--- Offene Uebergaenge nach Delta phi_7 ---\n";
    cout << "Offene " << R1 << "<->" << R2
         << "-Uebergaenge: " << open_transitions << "\n";
    for (int d = 0; d < 7; ++d) {
        cout << "dphi=" << d
             << " n=" << std::setw(8) << delta_phase_stats[d].n
             << " mean=" << std::setw(14) << delta_phase_stats[d].mean()
             << " sd=" << std::setw(14) << delta_phase_stats[d].sd()
             << "\n";
    }

    if (!open.empty()) {
        run_phase_shuffle_test(open, 1000);
    }

    cout << "\n--- Gegenphasen-Tests d gegen -d mod 7 ---\n";
    for (int d : {1,2,3}) {
        int dm = (7 - d) % 7;
        double diff = delta_phase_stats[d].mean() - delta_phase_stats[dm].mean();
        double pval = welch_z_pvalue(delta_phase_stats[d], delta_phase_stats[dm]);

        cout << d << " vs " << dm
             << " diff=" << std::setw(14) << diff
             << " p≈" << pval
             << "\n";
    }

    cout << "\n--- Quaternionischer 16x16-Kohaerenzfilter ---\n";
    Mat4Q K = build_phase_coherence_matrix(delta_phase_stats);
    double qfro = frobenius_quat(K);
    Mat16 B = to_real_16x16(K);
    vector<double> svals = singular_values_16(B);
    double pur = svd_purity(svals);

    cout << "Frobenius-Norm Quaternionmatrix: " << qfro << "\n";
    cout << "Top-8-Singulaerwerte: ";
    for (int i = 0; i < 8 && i < (int)svals.size(); ++i) {
        cout << svals[i] << (i + 1 < 8 ? ", " : "");
    }
    cout << "\n";
    cout << "SVD-Reinheit P = " << pur << "\n";
    double deff = (pur > 1e-30) ? 1.0 / pur : 0.0;
    double purity_over_iso = pur / (1.0 / 16.0);
    cout << "Effektive Kohaerenzdimension d_eff = 1/P = " << deff << "\n";
    cout << "Normiert gegen 16D-Isotropie: P/(1/16) = " << purity_over_iso << "\n";
    cout << "Isotroper Grenzwert 1/16 = " << (1.0/16.0) << "\n";

    double d16 = delta_phase_stats[1].mean() - delta_phase_stats[6].mean();
    double p16 = welch_z_pvalue(delta_phase_stats[1], delta_phase_stats[6]);
    double d25 = delta_phase_stats[2].mean() - delta_phase_stats[5].mean();
    double p25 = welch_z_pvalue(delta_phase_stats[2], delta_phase_stats[5]);
    double d34 = delta_phase_stats[3].mean() - delta_phase_stats[4].mean();
    double p34 = welch_z_pvalue(delta_phase_stats[3], delta_phase_stats[4]);

    std::string dominant_mode = "d16";
    double a16 = std::fabs(d16);
    double a25 = std::fabs(d25);
    double a34 = std::fabs(d34);
    if (a25 >= a16 && a25 >= a34) {
        dominant_mode = "d25";
    } else if (a34 >= a16 && a34 >= a25) {
        dominant_mode = "d34";
    }
    cout << "Dominant mode (|d|): " << dominant_mode << "\n";

    cout << "\nCSV_SUMMARY\n";
    cout << "N,R1,R2,dist,dominant_mode,twin_count,sector_count,delta_ons,p_ons,"
            "d16,p16,d25,p25,d34,p34,purity,d_eff,purity_over_iso\n";
    // defaultfloat+hohe Praezision: winzige p-Werte gehen in CSV nicht als 0.000000000 verloren
    // (oben: fixed setprecision(9) fuer Lesetext).
    cout << std::defaultfloat << std::setprecision(17);
    cout << N << ","
         << R1 << ","
         << R2 << ","
         << cyclic_dist << ","
         << dominant_mode << ","
         << twins.size() << ","
         << sector.size() << ","
         << delta_ons << ","
         << p_ons << ","
         << d16 << ","
         << p16 << ","
         << d25 << ","
         << p25 << ","
         << d34 << ","
         << p34 << ","
         << pur << ","
         << deff << ","
         << purity_over_iso
         << "\n";
    cout << std::fixed << std::setprecision(9);

    cout << "\nInterpretation:\n";
    cout << "- Diagonalwerte " << R1 << "->" << R1
         << " und " << R2 << "->" << R2
         << " nahe 0: ptolemaeisch geschlossen.\n";
    cout << "- Nebendiagonale " << R1 << "<->" << R2
         << " offen, aber aggregiert auf Reziprozitaet testbar.\n";
    cout << "- dphi=1/6 und 2/5 sollten Drift zeigen; dphi=3/4 eher reziprok.\n";
    cout << "- Der 16x16-Kohaerenzfilter misst, ob die Phasenstruktur isotrop oder konzentriert ist.\n";

    return 0;
}