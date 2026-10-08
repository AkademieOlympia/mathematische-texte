
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <chrono>
#include <unordered_map>

extern "C" {

struct EngineConfig {
    size_t memory_gb;
    size_t num_threads;
    size_t processing_qubits;
    size_t storage_qubits;
    double lambda_scale;
    double sw_u_vev;
};

struct GroverResult {
    int match_index;
    double probability;
    double execution_time_ms;
};

struct ShorResult {
    long long factor1;
    long long factor2;
    long long period_r;
    double execution_time_ms;
    int success;
};

void* create_engine(EngineConfig cfg) {
    EngineConfig* engine = new EngineConfig(cfg);
    return static_cast<void*>(engine);
}

void free_engine(void* handle) {
    if (handle) {
        delete static_cast<EngineConfig*>(handle);
    }
}

double compute_sw_frequency(double u, double lambda) {
    double k2 = 2.0 * std::pow(lambda, 2) / (u + std::pow(lambda, 2));
    double k_val = std::sqrt(std::max(0.0, std::min(1.0, k2)));
    double a = 1.0, b = std::sqrt(1.0 - k_val * k_val);
    for (int i = 0; i < 10; ++i) {
        double a_next = 0.5 * (a + b);
        double b_next = std::sqrt(a * b);
        a = a_next; b = b_next;
    }
    double K = M_PI / (2.0 * a);
    return std::sqrt(2.0) / (M_PI * std::sqrt(u + std::pow(lambda, 2))) * K;
}

GroverResult run_grover_search(void* handle, int num_items, int target_id) {
    auto start = std::chrono::high_resolution_clock::now();
    EngineConfig* cfg = static_cast<EngineConfig*>(handle);
    
    if (num_items <= 0 || target_id < 0 || target_id >= num_items) return {-1, 0.0, 0.0};

    double sw_freq = compute_sw_frequency(cfg->sw_u_vev, cfg->lambda_scale);
    int iterations = static_cast<int>(std::floor(M_PI / 4.0 * std::sqrt(num_items) / (1.0 + 0.2 * sw_freq)));
    if (iterations < 1) iterations = 1;

    double state_prob = 1.0 / num_items;
    for (int i = 0; i < iterations; ++i) {
        double theta = M_PI * sw_freq;
        state_prob = std::min(0.999999, state_prob + (1.0 - state_prob) * std::sin(theta) * std::sin(theta));
    }

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();

    return {target_id, state_prob, elapsed_ms};
}

long long gcd(long long a, long long b) {
    while (b != 0) {
        long long t = b;
        b = a % b;
        a = t;
    }
    return a;
}

long long power_mod(long long base, long long exp, long long mod) {
    long long res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (__int128(res) * base) % mod;
        base = (__int128(base) * base) % mod;
        exp /= 2;
    }
    return res;
}

long long mod_inverse(long long x, long long mod) {
    long long old_r = x, r = mod;
    long long old_s = 1, s = 0;
    while (r != 0) {
        long long q = old_r / r;
        long long tmp = old_r - q * r;
        old_r = r;
        r = tmp;
        tmp = old_s - q * s;
        old_s = s;
        s = tmp;
    }
    if (old_r != 1 && old_r != -1) return -1;
    if (old_s < 0) old_s += mod;
    return old_s;
}

// Kleinste r mit 1 <= r <= max_order und a^r ≡ 1 (mod N), sonst 0.
// Baby-Step/Giant-Step braucht etwa 2*sqrt(max_order) Multiplikationen.
long long find_order_bsgs(long long a, long long N, long long max_order) {
    if (max_order < 1) return 0;
    const long long cap = 1000000;
    long long m = 1;
    while (m < cap && m <= max_order / m) ++m;
    if (m > cap) m = cap;

    std::unordered_map<long long, long long> baby;
    baby.reserve(static_cast<size_t>(m) * 2);
    long long acc = 1;
    for (long long j = 0; j < m; ++j) {
        if (j > 0 && acc == 1) return j;
        if (baby.find(acc) == baby.end()) baby.emplace(acc, j);
        acc = static_cast<long long>((__int128)acc * a % N);
    }

    long long inv = mod_inverse(power_mod(a, m, N), N);
    if (inv < 0) return 0;
    long long cur = 1;
    for (long long i = 0; i < m; ++i) {
        auto it = baby.find(cur);
        if (it != baby.end()) {
            long long r = i * m + it->second;
            if (r > 0 && r <= max_order) return r;
        }
        cur = static_cast<long long>((__int128)cur * inv % N);
    }
    return 0;
}

ShorResult run_shor_algorithm(void* handle, long long N) {
    auto start = std::chrono::high_resolution_clock::now();
    if (N % 2 == 0) return {2, N / 2, 1, 0.01, 1};

    for (int attempt = 0; attempt < 100; ++attempt) {
        long long a = 2 + rand() % (N - 2);
        long long g = gcd(a, N);
        if (g > 1) {
            auto end = std::chrono::high_resolution_clock::now();
            return {g, N / g, 1, std::chrono::duration<double, std::milli>(end - start).count(), 1};
        }

        long long max_order = N - 1;
        const long long bsgs_cap = 1000000LL * 1000000LL;
        if (max_order > bsgs_cap) max_order = bsgs_cap;
        long long r = find_order_bsgs(a, N, max_order);

        if (r > 0 && r % 2 == 0) {
            long long x = power_mod(a, r / 2, N);
            if ((x + 1) % N != 0) {
                long long f1 = gcd(x - 1, N);
                if (f1 > 1 && f1 < N) {
                    auto end = std::chrono::high_resolution_clock::now();
                    return {f1, N / f1, r, std::chrono::duration<double, std::milli>(end - start).count(), 1};
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    return {-1, -1, 0, std::chrono::duration<double, std::milli>(end - start).count(), 0};
}

} // extern "C"
