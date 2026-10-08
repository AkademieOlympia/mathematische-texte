#include <cstdint>
#include <iostream>
#include <limits>

// Multiplication modulo mod without 64-bit overflow.
// Requires a compiler with unsigned __int128 support (GCC / Clang).
static std::uint64_t mul_mod(std::uint64_t a,
                             std::uint64_t b,
                             std::uint64_t mod)
{
    return static_cast<std::uint64_t>(
        (static_cast<unsigned __int128>(a) * b) % mod
    );
}

static std::uint64_t pow_mod(std::uint64_t base,
                             std::uint64_t exp,
                             std::uint64_t mod)
{
    std::uint64_t result = 1;

    while (exp > 0) {
        if (exp & 1ULL) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1ULL;
    }

    return result;
}

// Deterministic Miller-Rabin primality test for all 64-bit unsigned integers.
static bool is_prime(std::uint64_t n)
{
    if (n < 2) return false;

    // Remove small primes immediately.
    static constexpr std::uint64_t small_primes[] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
    };

    for (std::uint64_t p : small_primes) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    // Write n - 1 = d * 2^s, with d odd.
    std::uint64_t d = n - 1;
    unsigned int s = 0;

    while ((d & 1ULL) == 0) {
        d >>= 1ULL;
        ++s;
    }

    // Deterministic bases for n < 2^64.
    // This particular set is sufficient for the entire uint64_t range.
    static constexpr std::uint64_t bases[] = {
        2ULL,
        325ULL,
        9375ULL,
        28178ULL,
        450775ULL,
        9780504ULL,
        1795265022ULL
    };

    for (std::uint64_t a : bases) {
        if (a % n == 0) {
            continue;
        }

        std::uint64_t x = pow_mod(a % n, d, n);

        if (x == 1 || x == n - 1) {
            continue;
        }

        bool witness_passed = false;

        for (unsigned int r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);

            if (x == n - 1) {
                witness_passed = true;
                break;
            }
        }

        if (!witness_passed) {
            return false;
        }
    }

    return true;
}

// Returns the smallest prime strictly greater than n.
// Returns 0 if the uint64_t range is exhausted.
static std::uint64_t next_prime(std::uint64_t n)
{
    if (n < 2) {
        return 2;
    }

    if (n >= std::numeric_limits<std::uint64_t>::max() - 2) {
        return 0;
    }

    std::uint64_t candidate = n + 1;

    // Except for 2, primes are odd.
    if ((candidate & 1ULL) == 0) {
        ++candidate;
    }

    while (candidate > n) {
        if (is_prime(candidate)) {
            // Final Miller-Rabin verification before returning.
            if (is_prime(candidate)) {
                return candidate;
            }
        }

        if (candidate >
            std::numeric_limits<std::uint64_t>::max() - 2) {
            return 0;
        }

        candidate += 2;
    }

    return 0;
}

int main()
{
    std::uint64_t n;

    // Exactly one natural-number input.
    if (!(std::cin >> n)) {
        return 1;
    }

    const std::uint64_t p = next_prime(n);

    if (p == 0) {
        return 2;
    }

    // Exactly one output: the next prime.
    std::cout << p << '\n';

    return 0;
}
