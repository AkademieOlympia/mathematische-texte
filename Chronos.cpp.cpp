#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <complex>
#include <cctype>
#include <cstring>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <filesystem>
#include <omp.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/integer.hpp>

using namespace boost::multiprecision;
using namespace std;
namespace fs = std::filesystem;

// Definition des 256-Bit Typs
typedef uint256_t BigInt;
typedef cpp_bin_float_100 BigFloat;

namespace {

constexpr size_t k_max_sonar_zeros = 100000;
// Berry-Drift / Hyperschall v3: begrenzt für testbare Laufzeit (höher = langsamer)
constexpr int k_berry_gamma_cap = 100000;
constexpr int k_hyperschall_v3_steps = 1000;

bool parse_npy_shape_first_dim(const string& header, size_t& out_n) {
    size_t pos = header.find("'shape':");
    if (pos == string::npos) {
        pos = header.find("\"shape\":");
    }
    if (pos == string::npos) {
        return false;
    }
    size_t lparen = header.find('(', pos);
    if (lparen == string::npos) {
        return false;
    }
    ++lparen;
    while (lparen < header.size() && isspace(static_cast<unsigned char>(header[lparen]))) {
        ++lparen;
    }
    out_n = 0;
    while (lparen < header.size() && isdigit(static_cast<unsigned char>(header[lparen]))) {
        out_n = out_n * 10 + static_cast<size_t>(header[lparen] - '0');
        ++lparen;
    }
    return out_n > 0;
}

bool load_npy_float64_prefix(const string& path, vector<double>& out, size_t max_count) {
    ifstream f(path, ios::binary);
    if (!f) {
        return false;
    }
    char magic[6];
    f.read(magic, 6);
    static const char k_numpy_magic[] = {'\x93', 'N', 'U', 'M', 'P', 'Y'};
    if (!f || memcmp(magic, k_numpy_magic, 6) != 0) {
        return false;
    }
    auto v_major = static_cast<unsigned char>(f.get());
    auto v_minor = static_cast<unsigned char>(f.get());
    (void)v_minor;
    size_t header_len = 0;
    if (v_major == 1) {
        uint16_t hl = 0;
        f.read(reinterpret_cast<char*>(&hl), 2);
        header_len = hl;
    } else if (v_major == 2) {
        uint32_t hl = 0;
        f.read(reinterpret_cast<char*>(&hl), 4);
        header_len = hl;
    } else {
        return false;
    }
    string header(header_len, '\0');
    f.read(&header[0], static_cast<streamsize>(header_len));
    if (!f) {
        return false;
    }
    size_t n = 0;
    if (!parse_npy_shape_first_dim(header, n)) {
        return false;
    }
    size_t take = min(n, max_count);
    out.resize(take);
    f.read(reinterpret_cast<char*>(out.data()), static_cast<streamsize>(take * sizeof(double)));
    return static_cast<size_t>(f.gcount()) == take * sizeof(double);
}

void fill_simulated_zeros(vector<double>& z) {
    constexpr int n = 100000;
    z.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        z[static_cast<size_t>(i)] =
            14.0 + (static_cast<double>(i) / static_cast<double>(n - 1)) * 25000.0;
    }
}

vector<string> zeros_npy_candidates() {
    vector<string> c;
    c.push_back("zeros6.npy");
    try {
        fs::path here(__FILE__);
        c.push_back((here.parent_path() / "zeros6.npy").string());
    } catch (...) {
    }
    return c;
}

bool try_load_zeros(vector<double>& z) {
    for (const string& p : zeros_npy_candidates()) {
        if (load_npy_float64_prefix(p, z, k_max_sonar_zeros)) {
            cerr << "[*] Zeta-Nullstellen geladen: " << z.size() << " γ aus \"" << p << "\"\n";
            return true;
        }
    }
    return false;
}

} // namespace

class ChronosEngine {
private:
    BigInt N;
    vector<double> precomputed_zeros;

public:
    ChronosEngine(string n_str) : N(BigInt(n_str)) {
        if (!try_load_zeros(precomputed_zeros)) {
            fill_simulated_zeros(precomputed_zeros);
            cerr << "[*] zeros6.npy nicht lesbar — simulierte γ (" << precomputed_zeros.size()
                 << " Einträge).\n";
        }
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

    double calculate_vector_force(BigInt k) {
        const BigInt hund = BigInt(100);
        BigInt k_left = (k >= hund) ? (k - hund) : BigInt(0);
        BigInt k_right = k + hund;
        double res_left = calculate_track_asymmetry(k_left);
        double res_right = calculate_track_asymmetry(k_right);
        return res_right - res_left;
    }

    double calculate_real_riemann_sonar(BigInt k) {
        using LogFloat = cpp_bin_float_100;
        LogFloat x_blk = LogFloat(k);
        x_blk *= 210;
        LogFloat n_val = LogFloat(N);
        double ln_x = static_cast<double>(log(x_blk));
        double ln_N = static_cast<double>(log(n_val));

        double total_coherence = 0;
        const int n = static_cast<int>(precomputed_zeros.size());
        for (int i = 0; i < n; ++i) {
            double gamma = precomputed_zeros[i];
            total_coherence += std::cos(gamma * ln_x) * std::cos(gamma * ln_N);
        }
        return std::fabs(total_coherence / static_cast<double>(n));
    }

    static double warp_step_log2(BigInt k_step) {
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

    double calculate_berry_drift(BigInt k) {
        BigFloat ln_x = log(BigFloat(k * BigInt(210)));
        BigFloat ln_N = log(BigFloat(N));
        BigFloat dln = ln_x - ln_N;

        const int M = static_cast<int>(std::min(precomputed_zeros.size(), size_t(k_berry_gamma_cap)));
        if (M == 0) {
            return 0;
        }
        const double* gz = precomputed_zeros.data();

        double phase_drift = 0;
#pragma omp parallel for reduction(+ : phase_drift)
        for (int i = 0; i < M; ++i) {
            double gamma = gz[i];
            BigFloat th = BigFloat(gamma) * dln;
            phase_drift += std::sin(static_cast<double>(th));
        }
        return std::fabs(phase_drift / static_cast<double>(M));
    }

    void run_hyperschall_v3(BigInt k_start) {
        BigInt k = k_start;
        BigInt k_step = BigInt(1000);
        const BigInt k_step_cap = BigInt("1000000000000000000000000000000");
        int step = 0;

        cout << "[*] BERRY-PHASEN-SONAR AKTIVIERT. Suche spektrale Neigung..." << endl;

        while (step < k_hyperschall_v3_steps) {
            double drift = calculate_berry_drift(k);

            if (drift < 1e-8) {
                if (k_step < k_step_cap) {
                    BigInt doubled = k_step * 2;
                    if (doubled < k_step) {
                        k_step = k_step_cap;
                    } else if (doubled <= k_step_cap) {
                        k_step = doubled;
                    } else {
                        k_step = k_step_cap;
                    }
                }
            } else {
                k_step = 1;
                cout << "\n[!] BERRY-DRIFT detektiert bei k=" << k << " | Drift: " << drift << endl;
            }

            if (step % 50 == 0) {
                string ks = k.str();
                string k_tail = (ks.size() > 12) ? ks.substr(ks.size() - 12) : ks;
                cout << "Step " << setw(5) << step << " | k-Suffix: " << k_tail << " | Warp: 2^"
                     << fixed << setprecision(1) << warp_step_log2(k_step) << " | Drift: " << drift
                     << defaultfloat << endl;
            }

            if (drift > 0.001) {
                if (final_lock(k)) {
                    return;
                }
            }

            k += k_step;
            step++;
        }
    }

    void run_assault(BigInt k_start) {
        BigInt k = k_start;
        BigInt warp_factor = BigInt("1000000000000000");
        int step = 0;

        cout << "[*] CHRONOS-ENGINE: FILTER-MODUS AKTIVIERT" << endl;

        while (step < 20000) {
            double asymmetry = calculate_track_asymmetry(k);

            if (asymmetry < 1e-9) {
                if (step % 100 == 0) {
                    warp_factor *= 10;
                }
            } else {
                warp_factor = 1;
                cout << "\n[!] SYMMETRIEBRUCH DETEKTIERT bei k=" << k << endl;
                if (final_lock(k)) {
                    return;
                }
            }

            if (step % 50 == 0) {
                cout << "Step " << setw(5) << step << " | k: " << k.str().substr(0, 10) << "..."
                     << " | Warp: 10^" << (static_cast<int>(warp_factor.str().length()) - 1)
                     << " | Asym: " << asymmetry << endl;
            }

            k += warp_factor;
            step++;
        }
    }

    void run_vector_nav(BigInt k_start) {
        BigInt k = k_start;
        double momentum = 0;
        const double damping = 0.85;
        int step = 0;

        cout << "[*] VEKTOR-NAVIGATION AKTIVIERT. Suche Gradienten..." << endl;
        cout << "------------------------------------------------------------" << endl;

        while (step < 5000) {
            double current_asym = calculate_track_asymmetry(k);
            double force = calculate_vector_force(k);

            momentum = (momentum * damping) + (force * 1000000.0);

            double mag = std::fabs(momentum);
            if (!std::isfinite(mag) || mag < 1.0) {
                mag = 1.0;
            }
            const double k_step_cap = 1e12;
            if (mag > k_step_cap) {
                mag = k_step_cap;
            }
            unsigned long long step_u = static_cast<unsigned long long>(mag);
            if (step_u < 1) {
                step_u = 1;
            }
            BigInt move_dist = BigInt(step_u);
            if (current_asym > 1e-6) {
                move_dist = BigInt(1);
            }

            if (momentum > 0) {
                k += move_dist;
            } else if (k >= move_dist) {
                k -= move_dist;
            } else {
                k = 0;
            }

            string ks = k.str();
            string k_tail = (ks.size() > 10) ? ks.substr(ks.size() - 10) : ks;

            if (step % 100 == 0 || current_asym > 1e-7) {
                cout << "Step " << setw(5) << step << " | k: …" << k_tail
                     << " | Force: " << scientific << setprecision(2) << force
                     << " | Asym: " << current_asym << defaultfloat << endl;
            }

            if (current_asym > 0.1) {
                if (final_lock(k)) {
                    return;
                }
            }

            step++;
        }
    }

    bool final_lock(BigInt k_near) {
        // Suche in der unmittelbaren Umgebung des Resonanz-Peaks
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
};

int main() {
    // 1. Die 256-Bit Zielzahl N
    string n_256 = "91240241331123155652145421042111022334156788212344556778890123456789012345679";
    BigInt N(n_256);

    // 2. Ereignishorizont: sqrt(N) / Gitter-Basis 210
    BigInt k_target = sqrt(N) / BigInt(210);

    // 3. Sicherheits-Puffer (1 Mio. Blöcke vor dem Ziel-Sektor)
    BigInt k0 = k_target - BigInt(1000000);
    if (k0 < 0) {
        k0 = 0;
    }

    ChronosEngine engine(n_256);

    cout << "[*] WARP-ANTRIEB AKTIVIERT." << endl;
    cout << "[*] Ziel-Sektor k: " << k_target << endl;

    engine.run_hyperschall_v3(k0);

    return 0;
}