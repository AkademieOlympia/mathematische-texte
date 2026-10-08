#include <gmpxx.h>

#include <iostream>
#include <random>
#include <string>

// Arbitrary-precision: the next 10 primes starting at the given number (inclusive).
// Input: one non-negative integer in decimal, unlimited bit width.
// Output: ten primes, one per line.
//
// Build:
//   clang++ -std=c++17 -O2 -Wall -Wextra -Wno-deprecated-literal-operator \
//     -I/opt/homebrew/include -L/opt/homebrew/lib \
//     -o Zehn_ab_X "Zehn ab X.cpp" \
//     -lgmpxx -lgmp

namespace {

constexpr int kCount = 10;
constexpr int kProbabilisticRounds = 40;

const unsigned int kSmallPrimes[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
    157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
    239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
    331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
    421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
    509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
    709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811,
    821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911,
    919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997
};

const unsigned long kDeterministic64Bases[] = {
    2UL,
    325UL,
    9375UL,
    28178UL,
    450775UL,
    9780504UL,
    1795265022UL
};

gmp_randclass& rng()
{
    static gmp_randclass instance(gmp_randinit_mt);
    static const bool seeded = [] {
        std::random_device rd;
        unsigned long seed = rd();
        if (seed == 0) {
            seed = 1;
        }
        instance.seed(seed);
        return true;
    }();
    (void)seeded;
    return instance;
}

bool miller_rabin_round(const mpz_class& n,
                        const mpz_class& a,
                        const mpz_class& d,
                        unsigned int s)
{
    mpz_class x;
    mpz_powm(x.get_mpz_t(), a.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());

    if (x == 1 || x == n - 1) {
        return true;
    }

    for (unsigned int r = 1; r < s; ++r) {
        mpz_powm_ui(x.get_mpz_t(), x.get_mpz_t(), 2, n.get_mpz_t());
        if (x == n - 1) {
            return true;
        }
    }

    return false;
}

void split_odd_times_power_of_two(const mpz_class& n_minus_one,
                                  mpz_class& d,
                                  unsigned int& s)
{
    d = n_minus_one;
    s = static_cast<unsigned int>(mpz_scan1(d.get_mpz_t(), 0));
    mpz_fdiv_q_2exp(d.get_mpz_t(), d.get_mpz_t(), s);
}

bool is_prime(const mpz_class& n)
{
    if (n < 2) {
        return false;
    }

    for (unsigned int p : kSmallPrimes) {
        if (n == p) {
            return true;
        }
        if (mpz_divisible_ui_p(n.get_mpz_t(), p) != 0) {
            return false;
        }
    }

    mpz_class d;
    unsigned int s = 0;
    split_odd_times_power_of_two(n - 1, d, s);

    static const mpz_class kTwo64 = mpz_class(1) << 64;

    if (n < kTwo64) {
        for (unsigned long base : kDeterministic64Bases) {
            mpz_class a(base);
            if (a % n == 0) {
                continue;
            }
            if (!miller_rabin_round(n, a, d, s)) {
                return false;
            }
        }
        return true;
    }

    for (int i = 0; i < kProbabilisticRounds; ++i) {
        mpz_class a = rng().get_z_range(n - 3) + 2;
        if (!miller_rabin_round(n, a, d, s)) {
            return false;
        }
    }

    return true;
}

mpz_class next_prime(const mpz_class& n)
{
    if (n < 2) {
        return 2;
    }

    mpz_class candidate = n + 1;
    if (mpz_even_p(candidate.get_mpz_t()) != 0) {
        ++candidate;
    }

    while (true) {
        if (is_prime(candidate)) {
            return candidate;
        }
        candidate += 2;
    }
}

mpz_class first_prime_from(const mpz_class& n)
{
    if (n <= 2) {
        return 2;
    }
    if (is_prime(n)) {
        return n;
    }
    return next_prime(n);
}

}  // namespace

int main()
{
    std::string token;
    if (!(std::cin >> token)) {
        return 1;
    }

    mpz_class n;
    if (mpz_set_str(n.get_mpz_t(), token.c_str(), 10) != 0 || n < 0) {
        return 1;
    }

    mpz_class p = first_prime_from(n);
    for (int i = 0; i < kCount; ++i) {
        std::cout << p << '\n';
        if (i + 1 < kCount) {
            p = next_prime(p);
        }
    }

    return 0;
}
