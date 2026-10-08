// Primes_EABC.cpp
//
// EABC next-prime search SP(n) = min{q prime : q > n} with a commercial
// candidate filter. This is not a closed-form prime formula and not a
// primality proof. Outputs are probable primes under a documented test T.
//
// EABC filter: for q > 3 every prime lies in {1, 5, 7, 11} (mod 12),
// the families E, A, B, C. Small-prime sieve (p <= 7919) then Miller-Rabin
// or OpenSSL BN_check_prime on survivors.
//
// Benchmark protocol (same seeds, child-process isolation, fairness matrix):
//   seed = SHA-256 counter-hash ("eabc-v2" || bits || trial || counter)
//          then n := n mod 2^L, bit L-1 set. Not a NIST CTR-DRBG.
//   p_1  = nextprime(seed)   strictly greater than seed, p < 2^L
//   p_k  = nextprime(p_{k-1})
//   engines:
//     gmp_stock        mpz_nextprime as shipped
//     openssl_stock    odd walk + BN_check_prime as shipped
//     openssl_sieved   same EABC sieve as eabc_mr, BN_check_prime
//     eabc_mr          same EABC sieve, counted Miller-Rabin
//   Independent check: GMP mpz_probab_prime_p(reps=64) and OpenSSL
//   BN_check_prime. bitlen_ok and range_exhausted participate in failure.
//
// Build:
//   make Primes_EABC
//
// Usage:
//   ./Primes_EABC --next 100
//   ./Primes_EABC --bench --trials 5 --bits 512,1024 --out out_eabc

#include <gmpxx.h>
#include <openssl/bn.h>
#include <openssl/crypto.h>
#include <openssl/sha.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#include <vector>

#ifdef __APPLE__
#include <sys/sysctl.h>
#endif

namespace {

constexpr int kPrimeCountDefault = 10;
constexpr int kMrRoundsStock = 25;  // GMP mpz_probab_prime_p default reps
constexpr int kVerifyReps = 64;
constexpr unsigned kSieveLimit = 7919;
constexpr size_t kWindowCap = 4u * 1024u * 1024u;
constexpr const char* kMaster = "eabc-v2";
constexpr int kBitListDefault[] = {512, 1024, 1536, 2048, 3072, 4096};
constexpr const char* kEngineNames[] = {
    "gmp_stock", "openssl_stock", "openssl_sieved", "eabc_mr"};
constexpr int kEngineCount = 4;

enum class Tester { EabcMr, Openssl };

struct Stats {
    long long t_prime_ns = 0;
    long long t_10_ns = 0;
    long long cpu_ns = 0;
    long long rss_bytes = 0;
    long long mr_rounds_executed = 0;
    long long candidates_tested = 0;
    long long sieve_survivors = 0;
    int bitlen_ok = 1;
    int range_exhausted = 0;
    int verified = 0;
    int error = 0;
    std::string error_msg;
    std::vector<mpz_class> primes;
};

char eabc_family(const mpz_class& n) {
    if (n == 2 || n == 3) {
        return 'S';  // 2, 3 lie outside the four families
    }
    switch (mpz_fdiv_ui(n.get_mpz_t(), 12)) {
        case 1:
            return 'E';
        case 5:
            return 'A';
        case 7:
            return 'B';
        case 11:
            return 'C';
        default:
            return '?';
    }
}

bool eabc_residue_ok(unsigned long r12) {
    r12 %= 12;
    return r12 == 1 || r12 == 5 || r12 == 7 || r12 == 11;
}

void snap_eabc(mpz_class& n) {
    static const int kDelta[12] = {1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0};
    unsigned long r = mpz_fdiv_ui(n.get_mpz_t(), 12);
    n += kDelta[r];
}

std::vector<unsigned> make_sieve_primes(unsigned limit) {
    std::vector<uint8_t> is(static_cast<size_t>(limit) + 1, 1);
    is[0] = 0;
    is[1] = 0;
    for (unsigned p = 2; static_cast<unsigned long long>(p) * p <= limit; ++p) {
        if (!is[p]) {
            continue;
        }
        for (unsigned long long q = static_cast<unsigned long long>(p) * p;
             q <= limit; q += p) {
            is[static_cast<size_t>(q)] = 0;
        }
    }
    std::vector<unsigned> ps;
    ps.reserve(limit / 8);
    for (unsigned i = 2; i <= limit; ++i) {
        if (is[i]) {
            ps.push_back(i);
        }
    }
    return ps;
}

const std::vector<unsigned>& sieve_primes() {
    static const std::vector<unsigned> kPs = make_sieve_primes(kSieveLimit);
    return kPs;
}

int fips_mr_rounds(int bits) {
    if (bits >= 3747) {
        return 3;
    }
    if (bits >= 1345) {
        return 4;
    }
    if (bits >= 476) {
        return 5;
    }
    if (bits >= 400) {
        return 6;
    }
    if (bits >= 347) {
        return 7;
    }
    if (bits >= 308) {
        return 8;
    }
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

long long rusage_rss_bytes(const rusage& ru) {
#ifdef __linux__
    return static_cast<long long>(ru.ru_maxrss) * 1024LL;
#else
    return static_cast<long long>(ru.ru_maxrss);
#endif
}

long long self_rss_bytes() {
    rusage ru{};
    getrusage(RUSAGE_SELF, &ru);
    return rusage_rss_bytes(ru);
}

std::string hex_of(const mpz_class& n) {
    char* s = mpz_get_str(nullptr, 16, n.get_mpz_t());
    std::string out(s);
    free(s);
    return out;
}

std::string json_escape(const std::string& s) {
    std::string o;
    o.reserve(s.size() + 8);
    for (unsigned char c : s) {
        switch (c) {
            case '"':
                o += "\\\"";
                break;
            case '\\':
                o += "\\\\";
                break;
            case '\n':
                o += "\\n";
                break;
            case '\r':
                o += "\\r";
                break;
            case '\t':
                o += "\\t";
                break;
            default:
                if (c < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    o += buf;
                } else {
                    o.push_back(static_cast<char>(c));
                }
        }
    }
    return o;
}

void append_be32(std::vector<unsigned char>& m, uint32_t v) {
    m.push_back(static_cast<unsigned char>((v >> 24) & 0xff));
    m.push_back(static_cast<unsigned char>((v >> 16) & 0xff));
    m.push_back(static_cast<unsigned char>((v >> 8) & 0xff));
    m.push_back(static_cast<unsigned char>(v & 0xff));
}

void sha256_bytes(const unsigned char* data, size_t n, unsigned char out[SHA256_DIGEST_LENGTH]) {
    SHA256(data, n, out);
}

std::string sha256_hex(const std::string& data) {
    unsigned char out[SHA256_DIGEST_LENGTH];
    sha256_bytes(reinterpret_cast<const unsigned char*>(data.data()), data.size(), out);
    static const char* kHex = "0123456789abcdef";
    std::string hex(64, '0');
    for (int i = 0; i < SHA256_DIGEST_LENGTH; ++i) {
        hex[static_cast<size_t>(2 * i)] = kHex[(out[i] >> 4) & 0xf];
        hex[static_cast<size_t>(2 * i) + 1] = kHex[out[i] & 0xf];
    }
    return hex;
}

mpz_class make_seed(int bits, int trial) {
    std::vector<unsigned char> buf(static_cast<size_t>((bits + 7) / 8));
    uint32_t counter = 0;
    size_t filled = 0;
    const unsigned char* master =
        reinterpret_cast<const unsigned char*>(kMaster);
    const size_t master_len = std::strlen(kMaster);
    while (filled < buf.size()) {
        std::vector<unsigned char> msg;
        msg.insert(msg.end(), master, master + master_len);
        append_be32(msg, static_cast<uint32_t>(bits));
        append_be32(msg, static_cast<uint32_t>(trial));
        append_be32(msg, counter);
        unsigned char block[SHA256_DIGEST_LENGTH];
        sha256_bytes(msg.data(), msg.size(), block);
        size_t take = std::min(sizeof(block), buf.size() - filled);
        std::memcpy(buf.data() + filled, block, take);
        filled += take;
        ++counter;
    }
    mpz_class n;
    mpz_import(n.get_mpz_t(), buf.size(), 1, 1, 0, 0, buf.data());
    mpz_fdiv_r_2exp(n.get_mpz_t(), n.get_mpz_t(), static_cast<mp_bitcnt_t>(bits));
    mpz_setbit(n.get_mpz_t(), static_cast<mp_bitcnt_t>(bits - 1));
    return n;
}

mpz_class bit_limit(int bits) {
    if (bits <= 0) {
        return mpz_class(0);
    }
    mpz_class lim(1);
    mpz_mul_2exp(lim.get_mpz_t(), lim.get_mpz_t(), static_cast<mp_bitcnt_t>(bits));
    return lim;
}

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

bool miller_rabin_probable(const mpz_class& n, int rounds, Stats& st) {
    if (n < 2) {
        return false;
    }
    mpz_class d = n - 1;
    unsigned int s = static_cast<unsigned int>(mpz_scan1(d.get_mpz_t(), 0));
    mpz_fdiv_q_2exp(d.get_mpz_t(), d.get_mpz_t(), s);
    const auto& bases = sieve_primes();
    int used = 0;
    for (unsigned p : bases) {
        if (used >= rounds) {
            break;
        }
        mpz_class a(p);
        if (a >= n) {
            continue;
        }
        ++st.mr_rounds_executed;
        ++used;
        if (!miller_rabin_round(n, a, d, s)) {
            return false;
        }
    }
    if (used == 0) {
        return n == 2 || n == 3;
    }
    return true;
}

BIGNUM* mpz_to_bn(const mpz_class& n) {
    if (n == 0) {
        return BN_new();
    }
    size_t count = 0;
    void* buf = mpz_export(nullptr, &count, 1, 1, 0, 0, n.get_mpz_t());
    if (!buf || count == 0) {
        return nullptr;
    }
    BIGNUM* bn = BN_bin2bn(static_cast<unsigned char*>(buf), static_cast<int>(count), nullptr);
    free(buf);
    return bn;
}

mpz_class bn_to_mpz(const BIGNUM* bn) {
    int len = BN_num_bytes(bn);
    if (len <= 0) {
        return mpz_class(0);
    }
    std::vector<unsigned char> buf(static_cast<size_t>(len));
    BN_bn2bin(bn, buf.data());
    mpz_class n;
    mpz_import(n.get_mpz_t(), buf.size(), 1, 1, 0, 0, buf.data());
    return n;
}

int openssl_is_prime(const mpz_class& n, Stats& st) {
    BN_CTX* ctx = BN_CTX_new();
    BIGNUM* bn = mpz_to_bn(n);
    if (!ctx || !bn) {
        st.error = 1;
        st.error_msg = "BN_CTX_new/mpz_to_bn failed";
        if (bn) {
            BN_free(bn);
        }
        if (ctx) {
            BN_CTX_free(ctx);
        }
        return -1;
    }
    int rc = BN_check_prime(bn, ctx, nullptr);
    BN_free(bn);
    BN_CTX_free(ctx);
    if (rc < 0) {
        st.error = 1;
        st.error_msg = "BN_check_prime returned -1";
    }
    return rc;
}

void mark_multiples(std::vector<uint8_t>& mark,
                    const mpz_class& start,
                    unsigned p) {
    const size_t window = mark.size();
    mpz_class rem = start % p;
    unsigned long r = mpz_get_ui(rem.get_mpz_t());
    unsigned long offset = (r == 0) ? 0 : (p - r);
    mpz_class p2 = mpz_class(p) * p;
    mpz_class first = start + offset;
    if (first < p2) {
        mpz_class skip = p2 - first;
        mpz_class addz = (skip + (p - 1)) / p;
        offset += mpz_get_ui(addz.get_mpz_t()) * p;
    }
    for (unsigned long i = offset; i < window; i += p) {
        mark[i] = 0;
    }
}

void apply_eabc_and_sieve(std::vector<uint8_t>& mark, const mpz_class& start) {
    const size_t window = mark.size();
    unsigned long r0 = mpz_fdiv_ui(start.get_mpz_t(), 12);
    for (size_t i = 0; i < window; ++i) {
        if (!eabc_residue_ok(r0 + i)) {
            mark[i] = 0;
        }
    }
    for (unsigned p : sieve_primes()) {
        if (p == 2 || p == 3) {
            continue;  // already implied by the EABC residue filter
        }
        mark_multiples(mark, start, p);
    }
}

bool test_survivor(const mpz_class& cand, Tester tester, int mr_rounds, Stats& st) {
    ++st.sieve_survivors;
    ++st.candidates_tested;
    if (tester == Tester::EabcMr) {
        return miller_rabin_probable(cand, mr_rounds, st);
    }
    int rc = openssl_is_prime(cand, st);
    if (st.error) {
        return false;
    }
    return rc == 1;
}

void check_bitlen(Stats& st, int bits) {
    if (bits <= 0) {
        return;
    }
    for (const auto& p : st.primes) {
        if (mpz_sizeinbase(p.get_mpz_t(), 2) != static_cast<size_t>(bits)) {
            st.bitlen_ok = 0;
        }
    }
}

void finish_timing(Stats& st, long long t0, long long cpu0) {
    st.t_10_ns = wall_ns() - t0;
    st.cpu_ns = thread_cpu_ns() - cpu0;
    st.rss_bytes = self_rss_bytes();
}

Stats run_sieved(const mpz_class& seed, int bits, int count, int mr_rounds, Tester tester) {
    Stats st;
    long long cpu0 = thread_cpu_ns();
    long long t0 = wall_ns();
    const mpz_class lim = bit_limit(bits);

    mpz_class start = seed + 1;
    if (start < 2) {
        start = 2;
    }
    if (start == 2) {
        st.t_prime_ns = wall_ns() - t0;
        st.primes.emplace_back(2);
        start = 3;
    }
    if (st.primes.size() < static_cast<size_t>(count) && start <= 3) {
        if (st.primes.empty()) {
            st.t_prime_ns = wall_ns() - t0;
        }
        st.primes.emplace_back(3);
        start = 4;
    }
    snap_eabc(start);

    size_t window = std::max(static_cast<size_t>(16384), static_cast<size_t>(std::max(bits, 1)) * 96);
    window = std::min(window, kWindowCap);

    while (st.primes.size() < static_cast<size_t>(count) && !st.error) {
        if (bits > 0 && start >= lim) {
            st.range_exhausted = 1;
            st.bitlen_ok = 0;
            break;
        }
        size_t w = window;
        if (bits > 0) {
            mpz_class room = lim - start;
            if (room < 1) {
                st.range_exhausted = 1;
                st.bitlen_ok = 0;
                break;
            }
            if (mpz_fits_ulong_p(room.get_mpz_t()) != 0) {
                w = std::min(w, static_cast<size_t>(mpz_get_ui(room.get_mpz_t())));
            }
        }
        std::vector<uint8_t> mark(w, 1);
        apply_eabc_and_sieve(mark, start);
        for (size_t i = 0; i < w && st.primes.size() < static_cast<size_t>(count) && !st.error;
             ++i) {
            if (!mark[i]) {
                continue;
            }
            mpz_class cand = start + static_cast<unsigned long>(i);
            if (bits > 0 && cand >= lim) {
                st.range_exhausted = 1;
                st.bitlen_ok = 0;
                break;
            }
            if (test_survivor(cand, tester, mr_rounds, st)) {
                if (st.primes.empty()) {
                    st.t_prime_ns = wall_ns() - t0;
                }
                st.primes.push_back(cand);
            }
        }
        start += w;
        if (window < kWindowCap) {
            window = std::min(kWindowCap, window * 2);
        }
    }

    finish_timing(st, t0, cpu0);
    check_bitlen(st, bits);
    return st;
}

Stats run_gmp_stock(const mpz_class& seed, int bits, int count, int /*mr_rounds*/) {
    Stats st;
    long long cpu0 = thread_cpu_ns();
    long long t0 = wall_ns();
    const mpz_class lim = bit_limit(bits);
    mpz_class p = seed;
    for (int i = 0; i < count; ++i) {
        mpz_class nxt;
        mpz_nextprime(nxt.get_mpz_t(), p.get_mpz_t());
        if (bits > 0 && nxt >= lim) {
            st.range_exhausted = 1;
            st.bitlen_ok = 0;
            break;
        }
        if (st.primes.empty()) {
            st.t_prime_ns = wall_ns() - t0;
        }
        st.primes.push_back(nxt);
        p = nxt;
    }
    finish_timing(st, t0, cpu0);
    check_bitlen(st, bits);
    return st;
}

Stats run_openssl_stock(const mpz_class& seed, int bits, int count, int /*mr_rounds*/) {
    Stats st;
    long long cpu0 = thread_cpu_ns();
    long long t0 = wall_ns();
    const mpz_class lim = bit_limit(bits);
    BIGNUM* cand = mpz_to_bn(seed);
    BN_CTX* ctx = BN_CTX_new();
    if (!cand || !ctx) {
        st.error = 1;
        st.error_msg = "OpenSSL alloc failed";
        finish_timing(st, t0, cpu0);
        if (cand) {
            BN_free(cand);
        }
        if (ctx) {
            BN_CTX_free(ctx);
        }
        return st;
    }
    if (BN_add_word(cand, 1) != 1) {
        st.error = 1;
        st.error_msg = "BN_add_word failed";
        BN_free(cand);
        BN_CTX_free(ctx);
        finish_timing(st, t0, cpu0);
        return st;
    }
    if (!BN_is_odd(cand) && BN_add_word(cand, 1) != 1) {
        st.error = 1;
        st.error_msg = "BN_add_word failed";
        BN_free(cand);
        BN_CTX_free(ctx);
        finish_timing(st, t0, cpu0);
        return st;
    }

    auto over_limit = [&](const BIGNUM* n) {
        if (bits <= 0) {
            return false;
        }
        return BN_num_bits(n) > bits;
    };

    while (st.primes.size() < static_cast<size_t>(count) && !st.error) {
        if (over_limit(cand)) {
            st.range_exhausted = 1;
            st.bitlen_ok = 0;
            break;
        }
        int rc = BN_check_prime(cand, ctx, nullptr);
        ++st.candidates_tested;
        if (rc < 0) {
            st.error = 1;
            st.error_msg = "BN_check_prime returned -1";
            break;
        }
        if (rc == 1) {
            if (st.primes.empty()) {
                st.t_prime_ns = wall_ns() - t0;
            }
            st.primes.push_back(bn_to_mpz(cand));
        }
        if (BN_add_word(cand, 2) != 1) {
            st.error = 1;
            st.error_msg = "BN_add_word failed";
            break;
        }
    }
    BN_free(cand);
    BN_CTX_free(ctx);
    finish_timing(st, t0, cpu0);
    check_bitlen(st, bits);
    return st;
}

Stats run_engine(const std::string& name,
                 const mpz_class& seed,
                 int bits,
                 int count,
                 int mr_rounds) {
    if (name == "gmp_stock") {
        return run_gmp_stock(seed, bits, count, mr_rounds);
    }
    if (name == "openssl_stock") {
        return run_openssl_stock(seed, bits, count, mr_rounds);
    }
    if (name == "openssl_sieved") {
        return run_sieved(seed, bits, count, mr_rounds, Tester::Openssl);
    }
    if (name == "eabc_mr") {
        return run_sieved(seed, bits, count, mr_rounds, Tester::EabcMr);
    }
    Stats st;
    st.error = 1;
    st.error_msg = "unknown engine " + name;
    return st;
}

bool verify_prime_both(const mpz_class& n) {
    if (mpz_probab_prime_p(n.get_mpz_t(), kVerifyReps) == 0) {
        return false;
    }
    Stats tmp;
    int rc = openssl_is_prime(n, tmp);
    return rc == 1 && !tmp.error;
}

void verify_stats(Stats& st, int bits, int count) {
    if (st.error || st.range_exhausted || !st.bitlen_ok) {
        st.verified = 0;
        return;
    }
    if (static_cast<int>(st.primes.size()) != count) {
        st.verified = 0;
        return;
    }
    for (size_t i = 0; i < st.primes.size(); ++i) {
        if (i > 0 && st.primes[i] <= st.primes[i - 1]) {
            st.verified = 0;
            return;
        }
        if (bits > 0 && mpz_sizeinbase(st.primes[i].get_mpz_t(), 2) != static_cast<size_t>(bits)) {
            st.bitlen_ok = 0;
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
    os << "      \"rss_bytes\": " << st.rss_bytes << ",\n";
    os << "      \"mr_rounds_executed\": " << st.mr_rounds_executed << ",\n";
    os << "      \"candidates_tested\": " << st.candidates_tested << ",\n";
    os << "      \"sieve_survivors\": " << st.sieve_survivors << ",\n";
    os << "      \"bitlen_ok\": " << st.bitlen_ok << ",\n";
    os << "      \"range_exhausted\": " << st.range_exhausted << ",\n";
    os << "      \"error\": " << st.error << ",\n";
    os << "      \"error_msg\": \"" << json_escape(st.error_msg) << "\",\n";
    os << "      \"probable_prime_verified\": " << st.verified << ",\n";
    os << "      \"primes_sha256\": \"" << sha256_hex(concat) << "\"\n";
    os << "    }";
}

bool write_prime_list(const std::string& dir,
                      int bits,
                      int trial,
                      const char* engine,
                      const Stats& st) {
    std::ostringstream path;
    path << dir << "/" << bits << "_" << trial << "_" << engine << ".txt";
    std::ofstream out(path.str());
    if (!out) {
        return false;
    }
    for (const auto& p : st.primes) {
        out << p << " " << eabc_family(p) << "\n";
    }
    out.flush();
    return static_cast<bool>(out);
}

void write_worker_file(const std::string& path, const Stats& st) {
    std::ofstream out(path);
    out << "t_prime_ns=" << st.t_prime_ns << "\n";
    out << "t_10_ns=" << st.t_10_ns << "\n";
    out << "cpu_ns=" << st.cpu_ns << "\n";
    out << "rss_bytes=" << st.rss_bytes << "\n";
    out << "mr_rounds_executed=" << st.mr_rounds_executed << "\n";
    out << "candidates_tested=" << st.candidates_tested << "\n";
    out << "sieve_survivors=" << st.sieve_survivors << "\n";
    out << "bitlen_ok=" << st.bitlen_ok << "\n";
    out << "range_exhausted=" << st.range_exhausted << "\n";
    out << "error=" << st.error << "\n";
    out << "error_msg=" << st.error_msg << "\n";
    for (const auto& p : st.primes) {
        out << "prime=" << hex_of(p) << "\n";
    }
}

Stats read_worker_file(const std::string& path) {
    Stats st;
    std::ifstream in(path);
    if (!in) {
        st.error = 1;
        st.error_msg = "cannot read worker file";
        return st;
    }
    std::string line;
    while (std::getline(in, line)) {
        auto eq = line.find('=');
        if (eq == std::string::npos) {
            continue;
        }
        std::string k = line.substr(0, eq);
        std::string v = line.substr(eq + 1);
        if (k == "t_prime_ns") {
            st.t_prime_ns = std::strtoll(v.c_str(), nullptr, 10);
        } else if (k == "t_10_ns") {
            st.t_10_ns = std::strtoll(v.c_str(), nullptr, 10);
        } else if (k == "cpu_ns") {
            st.cpu_ns = std::strtoll(v.c_str(), nullptr, 10);
        } else if (k == "rss_bytes") {
            st.rss_bytes = std::strtoll(v.c_str(), nullptr, 10);
        } else if (k == "mr_rounds_executed") {
            st.mr_rounds_executed = std::strtoll(v.c_str(), nullptr, 10);
        } else if (k == "candidates_tested") {
            st.candidates_tested = std::strtoll(v.c_str(), nullptr, 10);
        } else if (k == "sieve_survivors") {
            st.sieve_survivors = std::strtoll(v.c_str(), nullptr, 10);
        } else if (k == "bitlen_ok") {
            st.bitlen_ok = std::atoi(v.c_str());
        } else if (k == "range_exhausted") {
            st.range_exhausted = std::atoi(v.c_str());
        } else if (k == "error") {
            st.error = std::atoi(v.c_str());
        } else if (k == "error_msg") {
            st.error_msg = v;
        } else if (k == "prime") {
            mpz_class p;
            if (mpz_set_str(p.get_mpz_t(), v.c_str(), 16) != 0) {
                st.error = 1;
                st.error_msg = "bad prime hex";
            } else {
                st.primes.push_back(p);
            }
        }
    }
    return st;
}

std::string sysctl_str(const char* name) {
#ifdef __APPLE__
    size_t n = 0;
    if (sysctlbyname(name, nullptr, &n, nullptr, 0) != 0 || n == 0) {
        return "";
    }
    std::string buf(n, '\0');
    if (sysctlbyname(name, buf.data(), &n, nullptr, 0) != 0) {
        return "";
    }
    while (!buf.empty() && buf.back() == '\0') {
        buf.pop_back();
    }
    return buf;
#else
    (void)name;
    return "";
#endif
}

long long sysctl_i64(const char* name) {
#ifdef __APPLE__
    long long v = 0;
    size_t n = sizeof(v);
    if (sysctlbyname(name, &v, &n, nullptr, 0) == 0) {
        return v;
    }
#else
    (void)name;
#endif
    return 0;
}

std::string compiler_id() {
#if defined(__clang__)
    return std::string("clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("gcc ") + __VERSION__;
#else
    return "unknown";
#endif
}

std::string os_id() {
    utsname u{};
    if (uname(&u) != 0) {
        return "unknown";
    }
    std::ostringstream os;
    os << u.sysname << " " << u.release << " " << u.machine;
    return os.str();
}

void emit_hardware_json(std::ostream& os) {
    std::string cpu = sysctl_str("machdep.cpu.brand_string");
    if (cpu.empty()) {
        cpu = sysctl_str("hw.model");
    }
    if (cpu.empty()) {
        cpu = "unknown";
    }
    long long ncpu = sysctl_i64("hw.ncpu");
    if (ncpu <= 0) {
        ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    }
    long long ram = sysctl_i64("hw.memsize");
    os << "  \"hardware\": {\n";
    os << "    \"cpu\": \"" << json_escape(cpu) << "\",\n";
    os << "    \"ncpu\": " << ncpu << ",\n";
    os << "    \"ram_bytes\": " << ram << ",\n";
    os << "    \"os\": \"" << json_escape(os_id()) << "\",\n";
    os << "    \"rss_unit\": \"bytes\"\n";
    os << "  },\n";
    os << "  \"libraries\": {\n";
    os << "    \"gmp\": \"" << json_escape(gmp_version) << "\",\n";
    os << "    \"openssl\": \"" << json_escape(OpenSSL_version(OPENSSL_VERSION)) << "\",\n";
    os << "    \"compiler\": \"" << json_escape(compiler_id()) << "\"\n";
    os << "  },\n";
}

bool parse_positive_int(const char* s, int& out, int lo, int hi) {
    if (!s || !*s) {
        return false;
    }
    char* end = nullptr;
    long v = std::strtol(s, &end, 10);
    if (end == s || *end != '\0' || v < lo || v > hi) {
        return false;
    }
    out = static_cast<int>(v);
    return true;
}

std::vector<int> parse_bits(const std::string& s) {
    std::vector<int> bits;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) {
            continue;
        }
        int v = 0;
        if (!parse_positive_int(item.c_str(), v, 8, 16384)) {
            return {};
        }
        bits.push_back(v);
    }
    return bits;
}

Stats run_isolated(const char* self,
                   const std::string& engine,
                   const mpz_class& seed,
                   int bits,
                   int count,
                   int mr_rounds,
                   const std::string& workdir) {
    std::ostringstream path;
    path << workdir << "/worker_" << engine << "_" << getpid() << "_" << bits << ".txt";
    const std::string outpath = path.str();
    const std::string seed_hex = hex_of(seed);
    const std::string bits_s = std::to_string(bits);
    const std::string count_s = std::to_string(count);
    const std::string mr_s = std::to_string(mr_rounds);

    pid_t pid = fork();
    if (pid < 0) {
        Stats st;
        st.error = 1;
        st.error_msg = "fork failed; running in-process";
        Stats inproc = run_engine(engine, seed, bits, count, mr_rounds);
        inproc.error_msg = st.error_msg + "; " + inproc.error_msg;
        return inproc;
    }
    if (pid == 0) {
        execl(self, self, "--worker", engine.c_str(), "--bits", bits_s.c_str(),
              "--count", count_s.c_str(), "--mr-rounds", mr_s.c_str(),
              "--seed-hex", seed_hex.c_str(), "--worker-out", outpath.c_str(),
              static_cast<char*>(nullptr));
        _exit(127);
    }
    int status = 0;
    rusage ru{};
    if (wait4(pid, &status, 0, &ru) < 0) {
        Stats st;
        st.error = 1;
        st.error_msg = "wait4 failed";
        return st;
    }
    Stats st = read_worker_file(outpath);
    st.rss_bytes = rusage_rss_bytes(ru);
    std::error_code ec;
    std::filesystem::remove(outpath, ec);
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        st.error = 1;
        if (st.error_msg.empty()) {
            st.error_msg = "worker exit " + std::to_string(status);
        }
    }
    return st;
}

void usage() {
    std::cerr
        << "Usage:\n"
        << "  Primes_EABC --next N [--count K] [--mr-rounds R]\n"
        << "  Primes_EABC --bench [--trials N] [--bits 512,1024,...] [--out DIR]\n"
        << "             [--count K] [--mr-rounds R] [--mr-policy stock|fips]\n"
        << "             [--no-isolate]\n";
}

struct Args {
    bool bench = false;
    bool next_mode = false;
    bool worker = false;
    bool isolate = true;
    std::string engine;
    std::string seed_hex;
    std::string worker_out;
    std::string next_n;
    std::string outdir = "out_eabc";
    std::string mr_policy = "stock";
    int trials = 5;
    int bits = 0;
    int count = kPrimeCountDefault;
    int mr_rounds = kMrRoundsStock;
    std::vector<int> bit_list;
};

int parse_args(int argc, char** argv, Args& a) {
    a.bit_list.assign(std::begin(kBitListDefault), std::end(kBitListDefault));
    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        auto need = [&](int& dst, int lo, int hi) {
            if (i + 1 >= argc || !parse_positive_int(argv[i + 1], dst, lo, hi)) {
                return false;
            }
            ++i;
            return true;
        };
        if (s == "--help") {
            usage();
            return 1;
        } else if (s == "--bench") {
            a.bench = true;
        } else if (s == "--next" && i + 1 < argc) {
            a.next_mode = true;
            a.next_n = argv[++i];
        } else if (s == "--worker" && i + 1 < argc) {
            a.worker = true;
            a.engine = argv[++i];
        } else if (s == "--no-isolate") {
            a.isolate = false;
        } else if (s == "--trials") {
            if (!need(a.trials, 1, 10000)) {
                return 2;
            }
        } else if (s == "--count") {
            if (!need(a.count, 1, 100)) {
                return 2;
            }
        } else if (s == "--mr-rounds") {
            if (!need(a.mr_rounds, 1, 64)) {
                return 2;
            }
        } else if (s == "--bits" && i + 1 < argc) {
            std::string raw = argv[++i];
            if (a.worker || a.next_mode) {
                if (!parse_positive_int(raw.c_str(), a.bits, 0, 16384)) {
                    return 2;
                }
            } else {
                a.bit_list = parse_bits(raw);
                if (a.bit_list.empty()) {
                    return 2;
                }
            }
        } else if (s == "--out" && i + 1 < argc) {
            a.outdir = argv[++i];
        } else if (s == "--seed-hex" && i + 1 < argc) {
            a.seed_hex = argv[++i];
        } else if (s == "--worker-out" && i + 1 < argc) {
            a.worker_out = argv[++i];
        } else if (s == "--mr-policy" && i + 1 < argc) {
            a.mr_policy = argv[++i];
            if (a.mr_policy != "stock" && a.mr_policy != "fips") {
                return 2;
            }
        } else {
            usage();
            return 2;
        }
    }
    if (!a.bench && !a.next_mode && !a.worker) {
        usage();
        return 2;
    }
    return 0;
}

int run_worker(const Args& a) {
    mpz_class seed;
    if (a.seed_hex.empty() || mpz_set_str(seed.get_mpz_t(), a.seed_hex.c_str(), 16) != 0) {
        std::cerr << "worker: bad --seed-hex\n";
        return 1;
    }
    Stats st = run_engine(a.engine, seed, a.bits, a.count, a.mr_rounds);
    if (a.worker_out.empty()) {
        write_worker_file("/dev/stdout", st);
    } else {
        write_worker_file(a.worker_out, st);
    }
    return st.error ? 4 : 0;
}

int run_next(const Args& a) {
    mpz_class n;
    if (mpz_set_str(n.get_mpz_t(), a.next_n.c_str(), 10) != 0) {
        std::cerr << "not a decimal integer: " << a.next_n << "\n";
        return 1;
    }
    Stats st = run_sieved(n, 0, a.count, a.mr_rounds, Tester::EabcMr);
    verify_stats(st, 0, a.count);
    for (const auto& p : st.primes) {
        std::cout << p << " " << eabc_family(p) << "\n";
    }
    if (!st.verified) {
        std::cerr << "probable-prime check failed or incomplete\n";
        return 3;
    }
    return 0;
}

int mr_rounds_for(const Args& a, int bits) {
    if (a.mr_policy == "fips") {
        return fips_mr_rounds(bits);
    }
    return a.mr_rounds;
}

int run_bench(char* argv0, const Args& a) {
    std::error_code ec;
    const std::string primeddir = a.outdir + "/primes";
    std::filesystem::create_directories(primeddir, ec);
    if (ec) {
        std::cerr << "cannot create " << primeddir << ": " << ec.message() << "\n";
        return 1;
    }
    const std::string workdir = a.outdir + "/workers";
    std::filesystem::create_directories(workdir, ec);
    if (ec) {
        std::cerr << "cannot create " << workdir << "\n";
        return 1;
    }

    std::ofstream json(a.outdir + "/results.json");
    if (!json) {
        std::cerr << "cannot write " << a.outdir << "/results.json\n";
        return 1;
    }

    json << "{\n";
    json << "  \"protocol\": \"random L-bit seed -> " << a.count
         << " next probable primes, child isolation, fairness matrix\",\n";
    json << "  \"claim_level\": \"probable_prime\",\n";
    json << "  \"master_seed\": \"" << kMaster << "\",\n";
    json << "  \"seed_construction\": \"SHA-256 counter hash, not NIST CTR-DRBG\",\n";
    json << "  \"eabc_filter\": \"mod 12 in {1,5,7,11} after handling 2 and 3\",\n";
    json << "  \"sieve_limit\": " << kSieveLimit << ",\n";
    emit_hardware_json(json);
    json << "  \"mr_policy\": \"" << json_escape(a.mr_policy) << "\",\n";
    json << "  \"mr_rounds_stock\": " << a.mr_rounds << ",\n";
    json << "  \"verify_reps\": " << kVerifyReps << ",\n";
    json << "  \"prime_count\": " << a.count << ",\n";
    json << "  \"trials\": " << a.trials << ",\n";
    json << "  \"isolate\": " << (a.isolate ? 1 : 0) << ",\n";
    json << "  \"engines\": [\"gmp_stock\",\"openssl_stock\",\"openssl_sieved\",\"eabc_mr\"],\n";
    json << "  \"runs\": [\n";

    bool first_run = true;
    int failures = 0;
    const char* self = argv0 ? argv0 : "./Primes_EABC";

    for (int L : a.bit_list) {
        const int rounds = mr_rounds_for(a, L);
        for (int trial = 0; trial < a.trials; ++trial) {
            mpz_class seed = make_seed(L, trial);
            std::cerr << "[bench] bits=" << L << " trial=" << trial
                      << " seed_bits=" << mpz_sizeinbase(seed.get_mpz_t(), 2)
                      << " isolate=" << a.isolate << " mr=" << rounds << "\n";

            std::vector<int> order(kEngineCount);
            for (int i = 0; i < kEngineCount; ++i) {
                order[i] = i;
            }
            std::mt19937 rng(0xEABCu + static_cast<unsigned>(L) * 10007u +
                             static_cast<unsigned>(trial));
            std::shuffle(order.begin(), order.end(), rng);

            Stats results[kEngineCount];
            for (int idx : order) {
                const char* name = kEngineNames[idx];
                if (a.isolate) {
                    results[idx] = run_isolated(self, name, seed, L, a.count, rounds, workdir);
                } else {
                    results[idx] = run_engine(name, seed, L, a.count, rounds);
                }
                verify_stats(results[idx], L, a.count);
            }

            bool agree = true;
            for (int i = 1; i < kEngineCount; ++i) {
                if (!same_primes(results[0], results[i])) {
                    agree = false;
                }
            }
            bool ok = agree;
            for (int i = 0; i < kEngineCount; ++i) {
                if (!results[i].verified || !results[i].bitlen_ok || results[i].error ||
                    results[i].range_exhausted) {
                    ok = false;
                }
            }
            if (!ok) {
                ++failures;
                std::cerr << "  FAIL bits=" << L << " trial=" << trial
                          << " agree=" << agree;
                for (int i = 0; i < kEngineCount; ++i) {
                    std::cerr << " " << kEngineNames[i] << "(v=" << results[i].verified
                              << ",b=" << results[i].bitlen_ok << ",e=" << results[i].error
                              << ",x=" << results[i].range_exhausted << ")";
                }
                std::cerr << "\n";
            } else {
                std::cerr << "  T10_ms";
                for (int i = 0; i < kEngineCount; ++i) {
                    std::cerr << " " << kEngineNames[i] << "="
                              << (results[i].t_10_ns / 1e6);
                }
                std::cerr << " mr_eabc=" << results[3].mr_rounds_executed << "\n";
            }

            for (int i = 0; i < kEngineCount; ++i) {
                if (!write_prime_list(primeddir, L, trial, kEngineNames[i], results[i])) {
                    ++failures;
                    std::cerr << "  FAIL write primes " << kEngineNames[i] << "\n";
                }
            }

            if (!first_run) {
                json << ",\n";
            }
            first_run = false;
            json << "    {\n";
            json << "      \"bits\": " << L << ",\n";
            json << "      \"trial\": " << trial << ",\n";
            json << "      \"mr_rounds_configured\": " << rounds << ",\n";
            json << "      \"fips_mr_rounds\": " << fips_mr_rounds(L) << ",\n";
            json << "      \"seed_hex\": \"" << hex_of(seed) << "\",\n";
            json << "      \"engine_order\": [";
            for (size_t i = 0; i < order.size(); ++i) {
                if (i) {
                    json << ",";
                }
                json << "\"" << kEngineNames[order[i]] << "\"";
            }
            json << "],\n";
            json << "      \"engines_agree\": " << (agree ? 1 : 0) << ",\n";
            json << "      \"engines\": {\n";
            for (int i = 0; i < kEngineCount; ++i) {
                emit_engine_json(json, kEngineNames[i], results[i]);
                if (i + 1 < kEngineCount) {
                    json << ",\n";
                }
            }
            json << "\n      }\n";
            json << "    }";
            json.flush();
            if (!json) {
                std::cerr << "JSON write failed\n";
                return 1;
            }
        }
    }

    json << "\n  ],\n";
    json << "  \"failures\": " << failures << "\n";
    json << "}\n";
    json.close();
    if (!json) {
        std::cerr << "JSON close failed\n";
        return 1;
    }
    std::cerr << "[bench] wrote " << a.outdir << "/results.json failures=" << failures << "\n";
    return failures == 0 ? 0 : 3;
}

}  // namespace

int main(int argc, char** argv) {
    Args a;
    int pa = parse_args(argc, argv, a);
    if (pa != 0) {
        return pa == 1 ? 0 : 1;
    }
    if (a.worker) {
        return run_worker(a);
    }
    if (a.next_mode) {
        return run_next(a);
    }
    return run_bench(argv[0], a);
}
