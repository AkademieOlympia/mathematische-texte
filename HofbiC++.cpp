#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

namespace {

constexpr std::array<u64, 9> kLocalFilterPrimes{
    3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL
};

constexpr std::array<u64, 12> kTrialDivisionPrimes{
    2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL
};

constexpr std::array<u64, 7> kMillerRabinBases{
    2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL
};

struct ClusterMetrics {
    double resonance;
    double walter_ratio;
    double morley_diff;
};

const char* classify_cluster(const ClusterMetrics& metrics, bool local_filter_ok) {
    if (
        local_filter_ok
        && std::abs(1.0 - metrics.resonance) < 1e-4
        && std::abs(1.0 - metrics.walter_ratio) < 1e-2
        && metrics.morley_diff < 1e-5
    ) {
        return "GEOMETRISCH REIF (Prim-Struktur)";
    }
    return "STRUKTURBRUCH (Zusammengesetzt/Instabil)";
}

u64 gcd_u64(u64 a, u64 b) {
    while (b != 0) {
        const u64 rest = a % b;
        a = b;
        b = rest;
    }
    return a;
}

}  // namespace

std::array<double, 4> log_cluster_embedding(u64 N, int k = 0) {
    const u64 m_k = static_cast<u64>(6 + 12 * k);
    return {
        std::log(static_cast<double>(N)),
        std::log(static_cast<double>(N + 2)),
        std::log(static_cast<double>(N + m_k)),
        std::log(static_cast<double>(N + m_k + 2))
    };
}

ClusterMetrics compute_cluster_metrics(u64 N, int k = 0) {
    const auto x = log_cluster_embedding(N, k);
    const double delta_e = std::abs(-2.0 * (x[2] - x[1]) * (x[3] - x[0]));
    const double theory_e = (64.0 + 288.0 * k + 288.0 * k * k)
        / (static_cast<double>(N) * static_cast<double>(N));

    const double walter_val = 0.05 * (x[0] + x[2]) * (x[3] - x[1]);
    const double expected = (static_cast<double>(6 + 12 * k) * x[0]) / (10.0 * static_cast<double>(N));
    const double left_edge = x[1] - x[0];
    const double right_edge = x[3] - x[2];

    return {
        theory_e > 0.0 ? delta_e / theory_e : 0.0,
        expected > 0.0 ? walter_val / expected : 0.0,
        std::abs(left_edge - right_edge)
    };
}

double get_hoeffbauer_resonance(u64 N, int k = 0) {
    const auto x = log_cluster_embedding(N, k);
    const double delta_e = std::abs(-2.0 * (x[2] - x[1]) * (x[3] - x[0]));
    const double theory_e = (64.0 + 288.0 * k + 288.0 * k * k)
        / (static_cast<double>(N) * static_cast<double>(N));
    return theory_e > 0.0 ? delta_e / theory_e : 0.0;
}

double calculate_walter_split(u64 N, int k = 0) {
    const auto x = log_cluster_embedding(N, k);
    return 0.05 * (x[0] + x[2]) * (x[3] - x[1]);
}

double expected_walter(u64 N, int k = 0) {
    const double m_k = static_cast<double>(6 + 12 * k);
    return (m_k * std::log(static_cast<double>(N))) / (10.0 * static_cast<double>(N));
}

double calculate_morley_diff(u64 N, int k = 0) {
    const auto x = log_cluster_embedding(N, k);
    const double left_edge = x[1] - x[0];
    const double right_edge = x[3] - x[2];
    return std::abs(left_edge - right_edge);
}

bool passes_local_prime_filter(u64 N) {
    if (N < 2) {
        return false;
    }
    if (N == 2) {
        return true;
    }
    if (N % 2 == 0) {
        return false;
    }
    for (u64 p : kLocalFilterPrimes) {
        if (N == p) {
            return true;
        }
        if (N % p == 0) {
            return false;
        }
    }
    return true;
}

std::string h_factor_v2(u64 N, int k = 0) {
    const bool local_filter_ok = passes_local_prime_filter(N);
    const auto metrics = compute_cluster_metrics(N, k);
    return classify_cluster(metrics, local_filter_ok);
}

u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % static_cast<u128>(mod));
}

u64 pow_mod(u64 base, u64 exp, u64 mod) {
    u64 result = 1;
    while (exp > 0) {
        if (exp & 1) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1;
    }
    return result;
}

bool miller_rabin_check(u64 n, u64 a, u64 d, int s) {
    u64 x = pow_mod(a % n, d, n);
    if (x == 1 || x == n - 1) {
        return true;
    }
    for (int r = 1; r < s; ++r) {
        x = mul_mod(x, x, n);
        if (x == n - 1) {
            return true;
        }
    }
    return false;
}

bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }
    for (u64 p : kTrialDivisionPrimes) {
        if (n == p) {
            return true;
        }
        if (n % p == 0) {
            return false;
        }
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }

    for (u64 a : kMillerRabinBases) {
        if (!miller_rabin_check(n, a, d, s)) {
            return false;
        }
    }
    return true;
}

u64 pollard_rho(u64 n) {
    if (n % 2 == 0) {
        return 2;
    }
    if (n % 3 == 0) {
        return 3;
    }

    static std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<u64> dist(2, n - 2);

    while (true) {
        const u64 c = dist(rng);
        u64 x = dist(rng);
        u64 y = x;
        u64 d = 1;

        auto f = [n, c](u64 v) {
            return (mul_mod(v, v, n) + c) % n;
        };

        while (d == 1) {
            x = f(x);
            y = f(f(y));
            const u64 diff = x > y ? x - y : y - x;
            d = gcd_u64(diff, n);
        }

        if (d != n) {
            return d;
        }
    }
}

void run_benchmark(const std::vector<u64>& test_numbers) {
    std::cout
        << std::setw(15) << "Zahl N" << " | "
        << std::setw(6) << "Prim?" << " | "
        << std::setw(10) << "Resonanz" << " | "
        << std::setw(42) << "GDA-Urteil" << " | "
        << "Pollard-Rho Time" << '\n';
    std::cout << std::string(110, '-') << '\n';

    for (u64 N : test_numbers) {
        const bool local_filter_ok = passes_local_prime_filter(N);
        const auto metrics = compute_cluster_metrics(N, 0);
        const char* gda_judgment = classify_cluster(metrics, local_filter_ok);

        const auto start_pollard = std::chrono::steady_clock::now();
        const bool actual_prime = is_prime(N);
        if (!actual_prime) {
            volatile u64 factor = pollard_rho(N);
            (void)factor;
        }
        const auto end_pollard = std::chrono::steady_clock::now();
        const std::chrono::duration<double> pollard_time = end_pollard - start_pollard;

        std::cout
            << std::setw(15) << N << " | "
            << std::setw(6) << (actual_prime ? "true" : "false") << " | "
            << std::setw(10) << std::fixed << std::setprecision(6) << metrics.resonance << " | "
            << std::setw(42) << gda_judgment << " | "
            << std::fixed << std::setprecision(6) << pollard_time.count() << "s\n";
    }
}

int main() {
    const std::vector<u64> test_set = {
        1000003ULL,
        1000005ULL,
        999999999999989ULL,
        999999999999991ULL,
        12345678901234567ULL
    };

    run_benchmark(test_set);
    return EXIT_SUCCESS;
}
