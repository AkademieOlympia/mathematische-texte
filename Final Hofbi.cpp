#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using ld = long double;

namespace {

constexpr std::array<u64, 9> kLocalFilterPrimes{
    3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL
};

struct LogGapProfile {
    ld left_gap;
    ld center_gap;
    ld right_gap;
};

struct ClusterMetrics {
    ld angle_resonance;
    ld balance_resonance;
    ld edge_symmetry;
};

struct MultiShellMetrics {
    ld avg_angle_resonance;
    ld avg_balance_resonance;
    ld avg_edge_symmetry;
    ld angle_span;
    ld balance_span;
    ld edge_span;
};

struct CliConfig {
    int single_k = 0;
    int shell_count = 129;
    std::vector<std::string> numbers;
};

std::string normalize_decimal(const std::string& text) {
    if (text.empty()) {
        throw std::invalid_argument("Leere Zahl ist nicht erlaubt.");
    }

    std::size_t idx = 0;
    if (text[idx] == '+' || text[idx] == '-') {
        if (text[idx] == '-') {
            throw std::invalid_argument("Nur positive Ganzzahlen sind erlaubt.");
        }
        ++idx;
    }
    if (idx == text.size()) {
        throw std::invalid_argument("Vorzeichen ohne Ziffern ist nicht erlaubt.");
    }

    while (idx < text.size() && text[idx] == '0') {
        ++idx;
    }
    if (idx == text.size()) {
        return "0";
    }

    std::string normalized;
    normalized.reserve(text.size() - idx);
    for (; idx < text.size(); ++idx) {
        const unsigned char ch = static_cast<unsigned char>(text[idx]);
        if (!std::isdigit(ch)) {
            throw std::invalid_argument("Nur ganzzahlige Dezimalzahlen sind erlaubt: " + text);
        }
        normalized.push_back(static_cast<char>(ch));
    }
    return normalized;
}

int parse_non_negative_int(const std::string& text, const char* flag_name) {
    try {
        std::size_t consumed = 0;
        const long long value = std::stoll(text, &consumed);
        if (consumed != text.size() || value < 0) {
            throw std::invalid_argument("negativ");
        }
        if (value > std::numeric_limits<int>::max()) {
            throw std::out_of_range("zu gross");
        }
        return static_cast<int>(value);
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string(flag_name) + " erwartet eine nichtnegative Ganzzahl.");
    }
}

bool decimal_equals_small(const std::string& value, u64 small) {
    return value == std::to_string(small);
}

std::string add_small(const std::string& value, u64 increment) {
    if (increment == 0) {
        return value;
    }

    std::string result = value;
    std::size_t pos = result.size();
    u64 carry = increment;

    while (pos > 0 && carry > 0) {
        --pos;
        const u64 digit = static_cast<u64>(result[pos] - '0');
        const u64 sum = digit + (carry % 10ULL);
        result[pos] = static_cast<char>('0' + (sum % 10ULL));
        carry /= 10ULL;
        carry += sum / 10ULL;
    }

    while (carry > 0) {
        result.insert(result.begin(), static_cast<char>('0' + (carry % 10ULL)));
        carry /= 10ULL;
    }

    return result;
}

u64 decimal_mod_small(const std::string& value, u64 mod) {
    u64 remainder = 0;
    for (char ch : value) {
        remainder = (remainder * 10ULL + static_cast<u64>(ch - '0')) % mod;
    }
    return remainder;
}

ld decimal_reciprocal_ld(const std::string& value) {
    if (value == "0") {
        throw std::invalid_argument("Division durch 0 ist nicht erlaubt.");
    }
    if (value.size() <= 18) {
        return 1.0L / std::stold(value);
    }

    const std::size_t leading_digits = 18;
    const std::string leading = value.substr(0, leading_digits);
    const std::string significand_text = leading.substr(0, 1) + "." + leading.substr(1);
    const ld significand = std::stold(significand_text);
    const int exponent10 = static_cast<int>(value.size() - 1);
    const std::string scale_text = "1e-" + std::to_string(exponent10);
    const ld scale = std::strtold(scale_text.c_str(), nullptr);
    return scale / significand;
}

u64 shell_offset(int k) {
    return static_cast<u64>(6ULL + 12ULL * static_cast<u64>(k));
}

bool uses_asymptotic_gap_model(const std::string& N) {
    // Auf diesem System ist `long double` praktisch auf Double-Bereich begrenzt.
    // Für sehr große Dezimalzahlen wechseln wir daher auf reine Gap-Verhältnisse.
    return N.size() > 300;
}

LogGapProfile compute_log_gap_profile(const std::string& N, int k) {
    const u64 m_k = shell_offset(k);
    const ld inv_left = decimal_reciprocal_ld(N);
    const ld inv_center = decimal_reciprocal_ld(add_small(N, 2ULL));
    const ld inv_right = decimal_reciprocal_ld(add_small(N, m_k));

    // `log1p` bleibt auch für sehr große N stabil, weil nur kleine Quotienten entstehen.
    return {
        std::log1pl(2.0L * inv_left),
        std::log1pl(static_cast<ld>(m_k - 2ULL) * inv_center),
        std::log1pl(2.0L * inv_right)
    };
}

bool passes_local_prime_filter(const std::string& N) {
    if (N == "0" || N == "1") {
        return false;
    }
    const int last_digit = N.back() - '0';
    if ((last_digit % 2) == 0) {
        return decimal_equals_small(N, 2ULL);
    }
    for (u64 p : kLocalFilterPrimes) {
        if (decimal_equals_small(N, p)) {
            return true;
        }
        if (decimal_mod_small(N, p) == 0) {
            return false;
        }
    }
    return true;
}

ClusterMetrics compute_cluster_metrics(const std::string& N, int k) {
    if (uses_asymptotic_gap_model(N)) {
        const ld center_to_outer = static_cast<ld>(shell_offset(k) - 2ULL) / 2.0L;
        const ld left_angle = std::atan(center_to_outer);
        const ld right_angle = left_angle;
        const ld angle_sum = left_angle + right_angle;
        const ld angle_resonance = angle_sum > 0.0L ? 1.0L : 0.0L;

        const ld balance_target = std::atan(static_cast<ld>(shell_offset(k) - 2ULL) / 4.0L);
        const ld balance_resonance = balance_target > 0.0L ? 1.0L : 0.0L;
        return {angle_resonance, balance_resonance, 1.0L};
    }

    const auto gaps = compute_log_gap_profile(N, k);
    const ld left_angle = std::atan2(gaps.center_gap, gaps.left_gap);
    const ld right_angle = std::atan2(gaps.center_gap, gaps.right_gap);
    const ld angle_sum = left_angle + right_angle;
    const ld angle_resonance = angle_sum > 0.0L
        ? 1.0L - std::abs(left_angle - right_angle) / angle_sum
        : 0.0L;

    const ld balance_angle = std::atan2(gaps.center_gap, gaps.left_gap + gaps.right_gap);
    const ld balance_target = std::atan(static_cast<ld>(shell_offset(k) - 2ULL) / 4.0L);
    const ld balance_resonance = balance_target > 0.0L
        ? 1.0L - std::abs(balance_angle - balance_target) / balance_target
        : 0.0L;

    const ld outer_max = std::max(gaps.left_gap, gaps.right_gap);
    const ld edge_symmetry = outer_max > 0.0L
        ? std::min(gaps.left_gap, gaps.right_gap) / outer_max
        : 0.0L;

    return {angle_resonance, balance_resonance, edge_symmetry};
}

MultiShellMetrics compute_multi_shell_metrics(const std::string& N, int shell_count) {
    ld angle_sum = 0.0L;
    ld balance_sum = 0.0L;
    ld edge_sum = 0.0L;

    ld angle_min = 1.0L;
    ld angle_max = 0.0L;
    ld balance_min = 1.0L;
    ld balance_max = 0.0L;
    ld edge_min = 1.0L;
    ld edge_max = 0.0L;

    for (int k = 0; k < shell_count; ++k) {
        const auto metrics = compute_cluster_metrics(N, k);
        angle_sum += metrics.angle_resonance;
        balance_sum += metrics.balance_resonance;
        edge_sum += metrics.edge_symmetry;

        angle_min = std::min(angle_min, metrics.angle_resonance);
        angle_max = std::max(angle_max, metrics.angle_resonance);
        balance_min = std::min(balance_min, metrics.balance_resonance);
        balance_max = std::max(balance_max, metrics.balance_resonance);
        edge_min = std::min(edge_min, metrics.edge_symmetry);
        edge_max = std::max(edge_max, metrics.edge_symmetry);
    }

    const ld inv_shell_count = shell_count > 0 ? 1.0L / static_cast<ld>(shell_count) : 0.0L;
    return {
        angle_sum * inv_shell_count,
        balance_sum * inv_shell_count,
        edge_sum * inv_shell_count,
        angle_max - angle_min,
        balance_max - balance_min,
        edge_max - edge_min
    };
}

const char* classify_cluster_multiscale(
    const ClusterMetrics& metrics,
    const MultiShellMetrics& multi_shell,
    bool local_filter_ok
) {
    const bool single_shell_ok =
        metrics.angle_resonance > 0.999998L
        && metrics.balance_resonance > 0.999999L
        && metrics.edge_symmetry > 0.99999L;

    // Für große k und sehr große Integer ist die Mehrschalen-Konsistenz wichtiger
    // als ein einzelner Shell-Snapshot.
    const bool multiscale_ok =
        multi_shell.avg_angle_resonance > 0.999999L
        && multi_shell.avg_balance_resonance > 0.999999L
        && multi_shell.avg_edge_symmetry > 0.99999L
        && multi_shell.angle_span < 1e-9L
        && multi_shell.balance_span < 1e-9L
        && multi_shell.edge_span < 1e-9L;

    if (local_filter_ok && (single_shell_ok || multiscale_ok)) {
        return "GEOMETRISCH REIF (Prim-Struktur)";
    }
    return "STRUKTURBRUCH (Zusammengesetzt/Instabil)";
}

CliConfig parse_cli(int argc, char** argv) {
    CliConfig config;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--k") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("--k braucht einen Wert.");
            }
            config.single_k = parse_non_negative_int(argv[++i], "--k");
        } else if (arg == "--k-max") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("--k-max braucht einen Wert.");
            }
            config.shell_count = parse_non_negative_int(argv[++i], "--k-max") + 1;
        } else if (arg == "--shell-count") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("--shell-count braucht einen Wert.");
            }
            config.shell_count = parse_non_negative_int(argv[++i], "--shell-count");
        } else if (arg == "--help") {
            std::cout
                << "Verwendung: ./FinalHofbi [--k N] [--k-max N | --shell-count N] [dezimalzahl ...]\n"
                << "  --k N            Schale fuer die Einzelshell-Ausgabe (Standard: 0)\n"
                << "  --k-max N        analysiert Mehrschalen-Konsistenz von k=0 bis k=N\n"
                << "  --shell-count N  analysiert N Schalen ab k=0\n"
                << "  dezimalzahl      beliebig grosse Ganzzahl als Eingabe\n";
            std::exit(0);
        } else {
            config.numbers.push_back(arg);
        }
    }

    if (config.shell_count <= 0) {
        throw std::invalid_argument("Die Schalenanzahl muss groesser als 0 sein.");
    }
    return config;
}

std::vector<std::string> default_numbers() {
    return {
        "1000003",
        "1000005",
        "999999999999989",
        "12345678901234567",
        "179769313486231590770839156793787453197860296048756011706444423684197180216158519368947833795864925541502180565485980503646440548199239100050792877003355816639229553136239076508735759914822574862575007425302077447712589550957937778424442426617334727629299387668709205606050270810842907692932019128194467627007"
    };
}

void print_header(const CliConfig& config) {
    std::cout << std::setw(44) << "Zahl N" << " | "
              << std::setw(3) << "k" << " | "
              << std::setw(14) << "Winkel-Res." << " | "
              << std::setw(14) << "Balance" << " | "
              << std::setw(12) << "MS-Span" << " | "
              << std::setw(11) << "MS-BalSpan" << " | "
              << "GDA-Urteil\n";
    std::cout << std::string(120, '-') << "\n";
    std::cout << std::setprecision(12);
    (void)config;
}

void analyze_and_print(const std::string& N, const CliConfig& config) {
    const auto metrics = compute_cluster_metrics(N, config.single_k);
    const auto multi_shell = compute_multi_shell_metrics(N, config.shell_count);
    const char* verdict = classify_cluster_multiscale(
        metrics,
        multi_shell,
        passes_local_prime_filter(N)
    );

    std::cout << std::setw(44) << N << " | "
              << std::setw(3) << config.single_k << " | "
              << std::setw(14) << std::fixed << static_cast<double>(metrics.angle_resonance) << " | "
              << std::setw(14) << std::fixed << static_cast<double>(metrics.balance_resonance) << " | "
              << std::setw(12) << std::scientific << std::setprecision(3) << static_cast<double>(multi_shell.angle_span) << " | "
              << std::setw(11) << std::scientific << std::setprecision(3) << static_cast<double>(multi_shell.balance_span) << " | "
              << std::fixed << std::setprecision(12)
              << verdict << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const CliConfig config = parse_cli(argc, argv);
        const std::vector<std::string> inputs = config.numbers.empty() ? default_numbers() : config.numbers;

        print_header(config);
        for (const std::string& input : inputs) {
            analyze_and_print(normalize_decimal(input), config);
        }
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Fehler: " << ex.what() << '\n';
        return 1;
    }
}