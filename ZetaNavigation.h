#ifndef ZETA_NAVIGATION_H
#define ZETA_NAVIGATION_H

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

/**
 * Standardzahl der aus der Datei geladenen Riemann-Nullstellen, die
 * fuer die Zeta-Interferenzsumme verwendet werden.
 */
constexpr std::size_t DEFAULT_ZETA_MOMENTUM_MAX_TERMS = 5000;

inline std::size_t& configured_zeta_momentum_max_terms() {
    static std::size_t value = DEFAULT_ZETA_MOMENTUM_MAX_TERMS;
    return value;
}

inline void set_zeta_momentum_max_terms(std::size_t value) {
    configured_zeta_momentum_max_terms() = std::max<std::size_t>(1, value);
}

inline std::size_t get_zeta_momentum_max_terms() {
    return configured_zeta_momentum_max_terms();
}

inline const std::vector<double>& fallback_riemann_gammas() {
    static const std::vector<double> fallback = {
        14.1347251417, 21.0220396388, 25.0108575801, 30.4248761259, 32.9350615877,
        37.5861781588, 40.9187190121, 43.3270732809, 48.0051508812, 49.7738324777,
        52.9703214777, 56.4462476971, 59.3470484566, 60.8317785246, 65.1125440481,
        67.0798105291, 69.5464017112, 72.0671576744, 75.7046906991, 77.1448400689,
        79.3373750202, 82.9103808541, 84.7354929805, 87.4252746131, 88.8091112076,
        92.4918992706, 94.6513440405, 95.8706342282, 98.8311942182, 101.317851006,
        103.725538040, 105.446623052, 107.168611184, 111.029535543, 111.874659177,
        114.320220915, 116.226680321, 118.790782866, 121.370125002, 122.946829294,
        124.256818554, 127.516483120, 129.578704200, 131.087688856, 133.497737203,
        134.756509753, 138.116042055, 139.736208952, 141.119793613, 143.111845808
    };
    return fallback;
}

inline std::size_t parse_npy_entry_count(const std::string& header) {
    const std::size_t shape_start = header.find('(');
    const std::size_t shape_end = header.find(')', shape_start);
    if (shape_start == std::string::npos || shape_end == std::string::npos || shape_end <= shape_start + 1) {
        return 0;
    }

    std::string digits;
    for (std::size_t i = shape_start + 1; i < shape_end; ++i) {
        const char ch = header[i];
        if (ch >= '0' && ch <= '9') {
            digits.push_back(ch);
        } else if (!digits.empty()) {
            break;
        }
    }

    return digits.empty() ? 0 : static_cast<std::size_t>(std::stoull(digits));
}

inline std::vector<double> load_riemann_gammas_from_npy(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        return {};
    }

    char magic[6] = {};
    input.read(magic, sizeof(magic));
    if (!input || std::string(magic, sizeof(magic)) != std::string("\x93NUMPY", 6)) {
        return {};
    }

    unsigned char major = 0;
    unsigned char minor = 0;
    input.read(reinterpret_cast<char*>(&major), 1);
    input.read(reinterpret_cast<char*>(&minor), 1);
    if (!input) {
        return {};
    }

    std::uint32_t header_length = 0;
    if (major == 1) {
        std::uint16_t header16 = 0;
        input.read(reinterpret_cast<char*>(&header16), sizeof(header16));
        header_length = header16;
    } else if (major == 2 || major == 3) {
        input.read(reinterpret_cast<char*>(&header_length), sizeof(header_length));
    } else {
        return {};
    }

    if (!input || header_length == 0) {
        return {};
    }

    std::string header(header_length, '\0');
    input.read(&header[0], static_cast<std::streamsize>(header.size()));
    if (!input) {
        return {};
    }

    if (header.find("'descr': '<f8'") == std::string::npos &&
        header.find("\"descr\": \"<f8\"") == std::string::npos) {
        return {};
    }

    if (header.find("True") != std::string::npos) {
        return {};
    }

    const std::size_t count = parse_npy_entry_count(header);
    if (count == 0) {
        return {};
    }

    std::vector<double> gammas(count);
    input.read(reinterpret_cast<char*>(gammas.data()), static_cast<std::streamsize>(count * sizeof(double)));
    if (!input) {
        return {};
    }

    return gammas;
}

inline const std::vector<double>& get_riemann_gammas() {
    static const std::vector<double> gammas = [] {
        std::vector<double> loaded = load_riemann_gammas_from_npy("zeros6.npy");
        if (!loaded.empty()) {
            return loaded;
        }
        return fallback_riemann_gammas();
    }();
    return gammas;
}

#endif
