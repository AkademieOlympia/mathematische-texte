// Commercial benchmark: random L-bit seed -> 10 following L-bit primes.
//
// Protocol
//   seed  = SHA-256-CTR("eabc-commercial-v1" || bits || trial), top bit set
//   p_1   = nextprime(seed)     // strictly greater than seed
//   p_k   = nextprime(p_{k-1})  // k = 2..10
//   engines, in isolation on the same process/hardware:
//     gmp_stock      mpz_nextprime (GMP 6.3.0 as shipped)
//     openssl_stock  odd walk + BN_check_prime (OpenSSL 3 as shipped)
//     prime_peano    EABC/210-wheel + small-prime sieve + counted Miller-Rabin
//   Independent check: every output prime is retested with BOTH libraries
//   (GMP mpz_probab_prime_p reps=64 and OpenSSL BN_check_prime). Engines must
//   also agree on the same 10 primes for a given seed.
//
// Metrics: T_prime, T_10, thread CPU ns, Miller-Rabin rounds (Prime-Peano),
// estimated heap, ru_maxrss.
//
// This is a next-prime chain, not RSA keygen of independent random primes.
// Prime-Peano implements SP(n) from prime_peano_rekursion_de.tex with a
// commercial candidate filter; it is not a closed-form prime formula.

#include <gmpxx.h>
#include <openssl/bn.h>
#include <openssl/crypto.h>
#include <openssl/sha.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#include <vector>

namespace {

constexpr int kPrimeCount = 10;
constexpr int kPeanoMrRounds = 25;  // match GMP mpz_nextprime default reps
constexpr int kVerifyReps = 64;
constexpr const char* kMaster = "eabc-commercial-v1";
constexpr int kBitListDefault[] = {512, 1024, 1536, 2048, 3072, 4096};

const unsigned int kSievePrimes[] = {
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
    919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013,
    1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091,
    1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181,
    1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1277, 1279,
    1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367,
    1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453,
    1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543,
    1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613,
    1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709,
    1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801,
    1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901,
    1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999,
    2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087,
    2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179,
    2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281,
    2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371,
    2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447,
    2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557,
    2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671,
    2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731,
    2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833,
    2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927,
    2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037,
    3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163,
    3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253,
    3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343,
    3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449, 3457,
    3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539,
    3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623,
    3631, 3637, 3643, 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719,
    3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823,
    3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919,
    3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019,
    4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127,
    4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229,
    4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327,
    4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441,
    4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523,
    4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643,
    4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733,
    4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861,
    4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951, 4967,
    4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051,
    5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167,
    5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
    5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399,
    5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471, 5477, 5479,
    5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569, 5573,
    5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683,
    5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
    5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, 5861, 5867,
    5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987, 6007,
    6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101,
    6113, 6121, 6131, 6133, 6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211,
    6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
    6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379,
    6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481, 6491, 6521, 6529,
    6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637,
    6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733,
    6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833,
    6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949,
    6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027,
    7039, 7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151, 7159,
    7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253,
    7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411,
    7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517,
    7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 7573, 7577, 7583, 7589,
    7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687, 7691,
    7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817,
    7823, 7829, 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919
};

int gcd_int(int a, int b) {
    while (b) {
        int t = a % b;
        a = b;
        b = t;
    }
    return a;
}

std::vector<int> wheel210() {
    std::vector<int> w;
    w.reserve(48);
    for (int r = 1; r < 210; ++r) {
        if (gcd_int(r, 210) == 1) {
            w.push_back(r);
        }
    }
    return w;
}

bool residue_ok_210(unsigned long r, const std::vector<int>& wheel) {
    r %= 210;
    return std::binary_search(wheel.begin(), wheel.end(), static_cast<int>(r));
}

// FIPS 186-5 / OpenSSL BN_prime_checks_for_size approximation.
int fips_mr_rounds(int bits) {
    if (bits >= 3747) return 3;
    if (bits >= 1345) return 4;
    if (bits >= 476) return 5;
    if (bits >= 400) return 6;
    if (bits >= 347) return 7;
    if (bits >= 308) return 8;
    return 27;
}

long long timespec_ns(const timespec& t) {
    return static_cast<long long>(t.tv_sec) * 1000000000LL + t.tv_nsec;
}

long long thread_cpu_ns() {
    timespec t{};
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t) != 0) {
        return 0;
    }
    return timespec_ns(t);
}

long long wall_ns() {
    timespec t{};
    clock_gettime(CLOCK_MONOTONIC_RAW, &t);
    return timespec_ns(t);
}

long long rss_kb() {
    rusage ru{};
    getrusage(RUSAGE_SELF, &ru);
    return static_cast<long long>(ru.ru_maxrss);
}

std::string hex_of(const mpz_class& n) {
    char* s = mpz_get_str(nullptr, 16, n.get_mpz_t());
    std::string out(s);
    free(s);
    return out;
}

std::string sha256_hex(const std::string& data) {
    unsigned char out[SHA256_DIGEST_LENGTH];
    SHA256(reinterpret_cast<const unsigned char*>(data.data()), data.size(), out);
    static const char* kHex = "0123456789abcdef";
    std::string hex;
    hex.resize(64);
    for (int i = 0; i < SHA256_DIGEST_LENGTH; ++i) {
        hex[2 * i] = kHex[(out[i] >> 4) & 0xf];
        hex[2 * i + 1] = kHex[out[i] & 0xf];
    }
    return hex;
}

mpz_class make_seed(int bits, int trial) {
    std::vector<unsigned char> buf((static_cast<size_t>(bits) + 7) / 8);
    uint32_t counter = 0;
    size_t filled = 0;
    while (filled < buf.size()) {
        unsigned char block[SHA256_DIGEST_LENGTH];
        SHA256_CTX ctx;
        SHA256_Init(&ctx);
        SHA256_Update(&ctx, kMaster, std::strlen(kMaster));
        unsigned char hdr[12];
        hdr[0] = static_cast<unsigned char>((bits >> 24) & 0xff);
        hdr[1] = static_cast<unsigned char>((bits >> 16) & 0xff);
        hdr[2] = static_cast<unsigned char>((bits >> 8) & 0xff);
        hdr[3] = static_cast<unsigned char>(bits & 0xff);
        hdr[4] = static_cast<unsigned char>((trial >> 24) & 0xff);
        hdr[5] = static_cast<unsigned char>((trial >> 16) & 0xff);
        hdr[6] = static_cast<unsigned char>((trial >> 8) & 0xff);
        hdr[7] = static_cast<unsigned char>(trial & 0xff);
        hdr[8] = static_cast<unsigned char>((counter >> 24) & 0xff);
        hdr[9] = static_cast<unsigned char>((counter >> 16) & 0xff);
        hdr[10] = static_cast<unsigned char>((counter >> 8) & 0xff);
        hdr[11] = static_cast<unsigned char>(counter & 0xff);
        SHA256_Update(&ctx, hdr, sizeof(hdr));
        SHA256_Final(block, &ctx);
        size_t take = std::min(sizeof(block), buf.size() - filled);
        std::memcpy(buf.data() + filled, block, take);
        filled += take;
        ++counter;
    }
    mpz_class n;
    mpz_import(n.get_mpz_t(), buf.size(), 1, 1, 0, 0, buf.data());
    mpz_setbit(n.get_mpz_t(), bits - 1);
    mpz_clrbit(n.get_mpz_t(), bits);
    if (mpz_sizeinbase(n.get_mpz_t(), 2) > static_cast<size_t>(bits)) {
        mpz_clrbit(n.get_mpz_t(), mpz_sizeinbase(n.get_mpz_t(), 2) - 1);
        mpz_setbit(n.get_mpz_t(), bits - 1);
    }
    return n;
}

struct Stats {
    long long t_prime_ns = 0;
    long long t_10_ns = 0;
    long long cpu_ns = 0;
    long long rss_kb_after = 0;
    long long mr_rounds = 0;
    long long candidates_tested = 0;
    long long sieve_survivors = 0;
    int bitlen_ok = 1;
    int verified = 0;
    std::vector<mpz_class> primes;
};

bool miller_rabin_round(const mpz_class& n,
                        const mpz_class& a,
                        const mpz_class& d,
                        unsigned int s) {
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

bool peano_is_prime(const mpz_class& n, int rounds, Stats& st) {
    if (n < 2) {
        return false;
    }
    for (unsigned int p : kSievePrimes) {
        if (n == p) {
            return true;
        }
        if (mpz_divisible_ui_p(n.get_mpz_t(), p) != 0) {
            return false;
        }
    }
    mpz_class d = n - 1;
    unsigned int s = static_cast<unsigned int>(mpz_scan1(d.get_mpz_t(), 0));
    mpz_fdiv_q_2exp(d.get_mpz_t(), d.get_mpz_t(), s);
    // Deterministic bases mixed with SHA-derived bases from n itself.
    static const unsigned long kBases[] = {2UL, 3UL, 5UL, 7UL, 11UL, 13UL, 17UL, 19UL, 23UL};
    int used = 0;
    for (unsigned long b : kBases) {
        if (used >= rounds) {
            break;
        }
        mpz_class a(b);
        if (a >= n) {
            continue;
        }
        ++st.mr_rounds;
        ++used;
        if (!miller_rabin_round(n, a, d, s)) {
            return false;
        }
    }
    mpz_class n_copy = n;
    while (used < rounds) {
        unsigned char digest[SHA256_DIGEST_LENGTH];
        std::string hex = hex_of(n_copy) + ":" + std::to_string(used);
        SHA256(reinterpret_cast<const unsigned char*>(hex.data()), hex.size(), digest);
        mpz_class a;
        mpz_import(a.get_mpz_t(), SHA256_DIGEST_LENGTH, 1, 1, 0, 0, digest);
        a %= (n - 3);
        a += 2;
        ++st.mr_rounds;
        ++used;
        if (!miller_rabin_round(n, a, d, s)) {
            return false;
        }
    }
    return true;
}

void check_bitlen(Stats& st, int bits) {
    for (const auto& p : st.primes) {
        if (mpz_sizeinbase(p.get_mpz_t(), 2) != static_cast<size_t>(bits)) {
            st.bitlen_ok = 0;
        }
    }
}

Stats run_gmp_stock(const mpz_class& seed, int bits) {
    Stats st;
    long long cpu0 = thread_cpu_ns();
    long long t0 = wall_ns();
    mpz_class p;
    mpz_nextprime(p.get_mpz_t(), seed.get_mpz_t());
    st.t_prime_ns = wall_ns() - t0;
    st.primes.push_back(p);
    for (int i = 1; i < kPrimeCount; ++i) {
        mpz_class nxt;
        mpz_nextprime(nxt.get_mpz_t(), p.get_mpz_t());
        p = nxt;
        st.primes.push_back(p);
    }
    long long t1 = wall_ns();
    st.t_10_ns = t1 - t0;
    st.cpu_ns = thread_cpu_ns() - cpu0;
    st.rss_kb_after = rss_kb();
    check_bitlen(st, bits);
    return st;
}

BIGNUM* mpz_to_bn(const mpz_class& n) {
    size_t count = 0;
    void* buf = mpz_export(nullptr, &count, 1, 1, 0, 0, n.get_mpz_t());
    BIGNUM* bn = BN_bin2bn(static_cast<unsigned char*>(buf), static_cast<int>(count), nullptr);
    free(buf);
    return bn;
}

mpz_class bn_to_mpz(const BIGNUM* bn) {
    int len = BN_num_bytes(bn);
    std::vector<unsigned char> buf(static_cast<size_t>(len));
    BN_bn2bin(bn, buf.data());
    mpz_class n;
    mpz_import(n.get_mpz_t(), buf.size(), 1, 1, 0, 0, buf.data());
    return n;
}

Stats run_openssl_stock(const mpz_class& seed, int bits) {
    Stats st;
    long long cpu0 = thread_cpu_ns();
    long long t0 = wall_ns();
    BIGNUM* cand = mpz_to_bn(seed);
    BN_add_word(cand, 1);
    if (!BN_is_odd(cand)) {
        BN_add_word(cand, 1);
    }
    BN_CTX* ctx = BN_CTX_new();
    auto next_prime = [&](BIGNUM* n) {
        while (true) {
            int rc = BN_check_prime(n, ctx, nullptr);
            ++st.candidates_tested;
            if (rc == 1) {
                return;
            }
            BN_add_word(n, 2);
        }
    };
    next_prime(cand);
    st.t_prime_ns = wall_ns() - t0;
    st.primes.push_back(bn_to_mpz(cand));
    for (int i = 1; i < kPrimeCount; ++i) {
        BN_add_word(cand, 2);
        next_prime(cand);
        st.primes.push_back(bn_to_mpz(cand));
    }
    st.t_10_ns = wall_ns() - t0;
    st.cpu_ns = thread_cpu_ns() - cpu0;
    st.rss_kb_after = rss_kb();
    BN_free(cand);
    BN_CTX_free(ctx);
    check_bitlen(st, bits);
    return st;
}

Stats run_prime_peano(const mpz_class& seed, int bits) {
    Stats st;
    const auto wheel = wheel210();
    const int rounds = kPeanoMrRounds;
    long long cpu0 = thread_cpu_ns();
    long long t0 = wall_ns();

    mpz_class start = seed + 1;
    size_t window = std::max(static_cast<size_t>(16384), static_cast<size_t>(bits) * 96);

    while (st.primes.size() < static_cast<size_t>(kPrimeCount)) {
        std::vector<uint8_t> mark(window, 1);
        for (unsigned int p : kSievePrimes) {
            mpz_class rem = start % p;
            unsigned long r = mpz_get_ui(rem.get_mpz_t());
            unsigned long offset = (r == 0) ? 0 : (p - r);
            for (unsigned long i = offset; i < window; i += p) {
                mark[i] = 0;
            }
        }
        for (size_t i = 0; i < window && st.primes.size() < static_cast<size_t>(kPrimeCount); ++i) {
            if (!mark[i]) {
                continue;
            }
            mpz_class cand = start + static_cast<unsigned long>(i);
            unsigned long r210 = mpz_fdiv_ui(cand.get_mpz_t(), 210);
            if (!residue_ok_210(r210, wheel)) {
                continue;
            }
            ++st.sieve_survivors;
            ++st.candidates_tested;
            if (peano_is_prime(cand, rounds, st)) {
                if (st.primes.empty()) {
                    st.t_prime_ns = wall_ns() - t0;
                }
                st.primes.push_back(cand);
            }
        }
        start += window;
        window *= 2;
    }

    st.t_10_ns = wall_ns() - t0;
    st.cpu_ns = thread_cpu_ns() - cpu0;
    st.rss_kb_after = rss_kb();
    check_bitlen(st, bits);
    return st;
}

bool verify_prime_both(const mpz_class& n) {
    if (mpz_probab_prime_p(n.get_mpz_t(), kVerifyReps) == 0) {
        return false;
    }
    BN_CTX* ctx = BN_CTX_new();
    BIGNUM* bn = mpz_to_bn(n);
    int rc = BN_check_prime(bn, ctx, nullptr);
    BN_free(bn);
    BN_CTX_free(ctx);
    return rc == 1;
}

void verify_stats(Stats& st) {
    if (st.primes.size() != static_cast<size_t>(kPrimeCount)) {
        st.verified = 0;
        return;
    }
    for (size_t i = 0; i < st.primes.size(); ++i) {
        if (i > 0 && st.primes[i] <= st.primes[i - 1]) {
            st.verified = 0;
            return;
        }
        if (!verify_prime_both(st.primes[i])) {
            st.verified = 0;
            return;
        }
    }
    st.verified = 1;
}

bool same_primes(const Stats& a, const Stats& b) {
    if (a.primes.size() != b.primes.size()) {
        return false;
    }
    for (size_t i = 0; i < a.primes.size(); ++i) {
        if (a.primes[i] != b.primes[i]) {
            return false;
        }
    }
    return true;
}

void emit_engine_json(std::ostream& os, const char* name, const Stats& st) {
    std::string concat;
    for (const auto& p : st.primes) {
        concat += hex_of(p);
        concat += ",";
    }
    os << "    \"" << name << "\": {\n";
    os << "      \"t_prime_ns\": " << st.t_prime_ns << ",\n";
    os << "      \"t_10_ns\": " << st.t_10_ns << ",\n";
    os << "      \"cpu_ns\": " << st.cpu_ns << ",\n";
    os << "      \"rss_kb_after\": " << st.rss_kb_after << ",\n";
    os << "      \"mr_rounds\": " << st.mr_rounds << ",\n";
    os << "      \"candidates_tested\": " << st.candidates_tested << ",\n";
    os << "      \"sieve_survivors\": " << st.sieve_survivors << ",\n";
    os << "      \"bitlen_ok\": " << st.bitlen_ok << ",\n";
    os << "      \"verified\": " << st.verified << ",\n";
    os << "      \"primes_sha256\": \"" << sha256_hex(concat) << "\"\n";
    os << "    }";
}

void write_prime_list(const std::string& dir,
                      int bits,
                      int trial,
                      const char* engine,
                      const Stats& st) {
    std::ostringstream path;
    path << dir << "/" << bits << "_" << trial << "_" << engine << ".txt";
    std::ofstream out(path.str());
    for (const auto& p : st.primes) {
        out << p << "\n";
    }
}

void usage() {
    std::cerr
        << "Usage: bench_next10 [--trials N] [--bits 512,1024,...] [--out DIR]\n";
}

std::vector<int> parse_bits(const std::string& s) {
    std::vector<int> bits;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) {
            bits.push_back(std::atoi(item.c_str()));
        }
    }
    return bits;
}

int parse_args(int argc, char** argv, int& trials, std::vector<int>& bits, std::string& outdir) {
    trials = 5;
    bits.assign(std::begin(kBitListDefault), std::end(kBitListDefault));
    outdir = "out";
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--trials" && i + 1 < argc) {
            trials = std::atoi(argv[++i]);
        } else if (a == "--bits" && i + 1 < argc) {
            bits = parse_bits(argv[++i]);
        } else if (a == "--out" && i + 1 < argc) {
            outdir = argv[++i];
        } else if (a == "--help") {
            usage();
            return 1;
        } else {
            usage();
            return 2;
        }
    }
    if (trials < 1) {
        trials = 1;
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    int trials = 5;
    std::vector<int> bits;
    std::string outdir;
    int pa = parse_args(argc, argv, trials, bits, outdir);
    if (pa != 0) {
        return pa == 1 ? 0 : 1;
    }

    std::string cmd = "mkdir -p '" + outdir + "/primes'";
    if (std::system(cmd.c_str()) != 0) {
        std::cerr << "cannot create " << outdir << "\n";
        return 1;
    }

    std::ofstream json(outdir + "/results.json");
    json << "{\n";
    json << "  \"protocol\": \"random L-bit seed -> 10 nextprimes, same seeds, same hardware\",\n";
    json << "  \"master_seed\": \"" << kMaster << "\",\n";
    json << "  \"hardware\": {\n";
    json << "    \"cpu\": \"Apple M5\",\n";
    json << "    \"ncpu\": 10,\n";
    json << "    \"ram_bytes\": 34359738368,\n";
    json << "    \"os\": \"macOS 26.7 arm64\"\n";
    json << "  },\n";
    json << "  \"libraries\": {\n";
    json << "    \"gmp\": \"6.3.0\",\n";
    json << "    \"openssl\": \"" << OpenSSL_version(OPENSSL_VERSION) << "\",\n";
    json << "    \"compiler\": \"Apple clang 21.0.0 -O3 -march=native\"\n";
    json << "  },\n";
    json << "  \"peano_mr_rounds\": " << kPeanoMrRounds << ",\n";
    json << "  \"verify_reps\": " << kVerifyReps << ",\n";
    json << "  \"prime_count\": " << kPrimeCount << ",\n";
    json << "  \"trials\": " << trials << ",\n";
    json << "  \"runs\": [\n";

    bool first_run = true;
    int failures = 0;

    for (int L : bits) {
        for (int trial = 0; trial < trials; ++trial) {
            mpz_class seed = make_seed(L, trial);
            std::cerr << "[bench] bits=" << L << " trial=" << trial
                      << " seed_bits=" << mpz_sizeinbase(seed.get_mpz_t(), 2) << "\n";

            Stats gmp = run_gmp_stock(seed, L);
            Stats ossl = run_openssl_stock(seed, L);
            Stats peano = run_prime_peano(seed, L);

            verify_stats(gmp);
            verify_stats(ossl);
            verify_stats(peano);

            bool agree = same_primes(gmp, ossl) && same_primes(gmp, peano);
            if (!gmp.verified || !ossl.verified || !peano.verified || !agree) {
                ++failures;
                std::cerr << "  FAIL verify/agree bits=" << L << " trial=" << trial
                          << " gmp=" << gmp.verified << " openssl=" << ossl.verified
                          << " peano=" << peano.verified << " agree=" << agree << "\n";
            } else {
                std::cerr << "  T10_ms gmp=" << (gmp.t_10_ns / 1e6)
                          << " openssl=" << (ossl.t_10_ns / 1e6)
                          << " peano=" << (peano.t_10_ns / 1e6)
                          << " mr_peano=" << peano.mr_rounds << "\n";
            }

            write_prime_list(outdir + "/primes", L, trial, "gmp", gmp);
            write_prime_list(outdir + "/primes", L, trial, "openssl", ossl);
            write_prime_list(outdir + "/primes", L, trial, "peano", peano);

            if (!first_run) {
                json << ",\n";
            }
            first_run = false;
            json << "    {\n";
            json << "      \"bits\": " << L << ",\n";
            json << "      \"trial\": " << trial << ",\n";
            json << "      \"fips_mr_rounds\": " << fips_mr_rounds(L) << ",\n";
            json << "      \"seed_hex\": \"" << hex_of(seed) << "\",\n";
            json << "      \"engines_agree\": " << (agree ? 1 : 0) << ",\n";
            json << "      \"engines\": {\n";
            emit_engine_json(json, "gmp_stock", gmp);
            json << ",\n";
            emit_engine_json(json, "openssl_stock", ossl);
            json << ",\n";
            emit_engine_json(json, "prime_peano", peano);
            json << "\n      }\n";
            json << "    }";
            json.flush();
        }
    }

    json << "\n  ],\n";
    json << "  \"failures\": " << failures << "\n";
    json << "}\n";
    json.close();

    std::cerr << "[bench] wrote " << outdir << "/results.json failures=" << failures << "\n";
    return failures == 0 ? 0 : 3;
}
