#ifndef E8_QRAM_SHOR_AVX512_KERNEL_HPP
#define E8_QRAM_SHOR_AVX512_KERNEL_HPP

#include <iostream>
#include <vector>
#include <array>
#include <complex>
#include <cmath>
#include <chrono>
#include <concepts>
#include <memory>
#include <numeric>
#include <algorithm>

// Cross-Platform Architecture Switch: x86_64 (AVX-512) vs. ARM64 (Apple Silicon NEON)
#if defined(__x86_64__) || defined(_M_X64)
    #include <immintrin.h>
#elif defined(__ARM_NEON) || defined(__aarch64__)
    #include <arm_neon.h>
#endif

namespace d9_shor {

template<typename T>
concept SIMDFloat = std::same_as<T, double> || std::same_as<T, float>;

template<SIMDFloat T = double>
class E8QRAMShorAVX512Kernel {
private:
    uint64_t target_N_;
    uint64_t generator_a_;
    size_t qram_capacity_;
    
    alignas(64) std::array<double, 240 * 8> e8_roots_flat_;

public:
    explicit E8QRAMShorAVX512Kernel(uint64_t N, uint64_t a, size_t capacity = 65536)
        : target_N_(N), generator_a_(a), qram_capacity_(capacity) {
        init_e8_root_lattice();
    }

private:
    void init_e8_root_lattice() noexcept {
        size_t idx = 0;
        for (int i = 0; i < 8 && idx < 240; ++i) {
            for (int j = i + 1; j < 8 && idx < 240; ++j) {
                for (int s1 : {-1, 1}) {
                    for (int s2 : {-1, 1}) {
                        if (idx < 240) {
                            for (int k = 0; k < 8; ++k) e8_roots_flat_[idx * 8 + k] = 0.0;
                            e8_roots_flat_[idx * 8 + i] = static_cast<double>(s1);
                            e8_roots_flat_[idx * 8 + j] = static_cast<double>(s2);
                            idx++;
                        }
                    }
                }
            }
        }
        for (int b = 0; b < 256 && idx < 240; ++b) {
            int minus_count = 0;
            for (int k = 0; k < 8; ++k) if ((b >> k) & 1) minus_count++;
            if (minus_count % 2 == 0) {
                for (int k = 0; k < 8; ++k) {
                    e8_roots_flat_[idx * 8 + k] = ((b >> k) & 1) ? -0.5 : 0.5;
                }
                idx++;
            }
        }
    }

public:
    inline size_t quantize_e8_avx512(const double* vec8) const noexcept {
#if defined(__AVX512F__)
        __m512d v_input0 = _mm512_loadu_pd(vec8);
        __m512d v_input1 = _mm512_loadu_pd(vec8 + 4);

        double min_dist_sq = 1e30;
        size_t best_idx = 0;

        for (size_t r = 0; r < 240; ++r) {
            const double* root_ptr = &e8_roots_flat_[r * 8];
            __m512d v_r0 = _mm512_loadu_pd(root_ptr);
            __m512d v_r1 = _mm512_loadu_pd(root_ptr + 4);

            __m512d diff0 = _mm512_sub_pd(v_input0, v_r0);
            __m512d diff1 = _mm512_sub_pd(v_input1, v_r1);

            __m512d sq0 = _mm512_mul_pd(diff0, diff0);
            __m512d sq1 = _mm512_mul_pd(diff1, diff1);

            double dist_sq = _mm512_reduce_add_pd(sq0) + _mm512_reduce_add_pd(sq1);

            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                best_idx = r;
            }
        }
        return best_idx;
#elif defined(__ARM_NEON) || defined(__aarch64__)
        float64x2_t v_in0 = vld1q_f64(vec8);
        float64x2_t v_in1 = vld1q_f64(vec8 + 2);
        float64x2_t v_in2 = vld1q_f64(vec8 + 4);
        float64x2_t v_in3 = vld1q_f64(vec8 + 6);

        double min_dist_sq = 1e30;
        size_t best_idx = 0;

        for (size_t r = 0; r < 240; ++r) {
            const double* root_ptr = &e8_roots_flat_[r * 8];
            
            float64x2_t d0 = vsubq_f64(v_in0, vld1q_f64(root_ptr));
            float64x2_t d1 = vsubq_f64(v_in1, vld1q_f64(root_ptr + 2));
            float64x2_t d2 = vsubq_f64(v_in2, vld1q_f64(root_ptr + 4));
            float64x2_t d3 = vsubq_f64(v_in3, vld1q_f64(root_ptr + 6));

            float64x2_t sq0 = vmulq_f64(d0, d0);
            float64x2_t sq1 = vmulq_f64(d1, d1);
            float64x2_t sq2 = vmulq_f64(d2, d2);
            float64x2_t sq3 = vmulq_f64(d3, d3);

            double dist_sq = vaddvq_f64(sq0) + vaddvq_f64(sq1) + vaddvq_f64(sq2) + vaddvq_f64(sq3);

            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                best_idx = r;
            }
        }
        return best_idx;
#else
        double min_dist_sq = 1e30;
        size_t best_idx = 0;
        for (size_t r = 0; r < 240; ++r) {
            double dist_sq = 0.0;
            for (int k = 0; k < 8; ++k) {
                double diff = vec8[k] - e8_roots_flat_[r * 8 + k];
                dist_sq += diff * diff;
            }
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                best_idx = r;
            }
        }
        return best_idx;
#endif
    }

    inline bool check_topos_1_omega_norm(double z_real, double z_imag) const noexcept {
        double phi = (1.0 + std::sqrt(5.0)) * 0.5;
        double norm = (z_real * z_real) - phi * (z_imag * z_imag);
        return std::abs(std::abs(norm) - 1.0) < 0.25;
    }

    uint64_t execute_shor_period_finding_avx512() const noexcept {
        uint64_t val = 1;
        uint64_t period = 0;

        alignas(64) double state_vec[8] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        for (uint64_t step = 1; step <= 65536; ++step) {
            val = (val * generator_a_) % target_N_;

            state_vec[0] = std::cos(2.0 * M_PI * val / target_N_);
            state_vec[1] = std::sin(2.0 * M_PI * val / target_N_);
            state_vec[2] = static_cast<double>(val) / target_N_;

            size_t root_idx = quantize_e8_avx512(state_vec);
            bool is_valid = check_topos_1_omega_norm(state_vec[0], state_vec[1]);

            if (is_valid && val == 1) {
                period = step;
                break;
            }
        }

        return period;
    }

    std::pair<uint64_t, uint64_t> run_factorization() const {
        if (target_N_ % 2 == 0) return {2, target_N_ / 2};

        uint64_t r = execute_shor_period_finding_avx512();
        if (r == 0 || r % 2 != 0) return {0, 0};

        uint64_t half_power = 1;
        for (uint64_t i = 0; i < r / 2; ++i) {
            half_power = (half_power * generator_a_) % target_N_;
        }

        if ((half_power + 1) % target_N_ == 0) return {0, 0};

        uint64_t p = std::gcd(half_power - 1, target_N_);
        uint64_t q = std::gcd(half_power + 1, target_N_);

        if (p * q == target_N_ && p > 1 && q > 1) {
            return {p, q};
        }
        return {0, 0};
    }
};

} // namespace d9_shor

#endif // E8_QRAM_SHOR_AVX512_KERNEL_HPP
