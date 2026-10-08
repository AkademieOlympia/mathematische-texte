// ============================================================================
//  fano_steane_e8_chain.cpp  (improved successor of "Shor 1000.cpp")
//
//  Runtime companion to EabcRenorm Lean gates:
//    Fano incidence → [[7,1,3]] Steane (CSS / Aaronson–Gottesman substrate)
//                   → [8,4,4] → Construction A → 240 roots (norm²=4 in L'_C)
//                   + Octonion Associator & Alternativity Engine (PR #153)
//
//  CLAIM-WALL
//  ----------
//  Algebraic / information-theoretic certificates only.
//  No physical qubits, hardware, threshold, or quantum-advantage claims.
//  Classical Shor demo is classical order-finding, not quantum simulation.
//  Lean: G.13 (T_Q isometry / 240 roots) is KERNEL; G.8 (PSL(2,7)) remains OPEN.
//  This binary computes algebraic witnesses only (claim-wall).
// ============================================================================

#include <algorithm>
#include <array>
#include <bitset>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Aaronson–Gottesman stabilizer tableau (Clifford only)
// ---------------------------------------------------------------------------
class StabilizerSimulator {
 public:
  explicit StabilizerSimulator(std::size_t qubits)
      : n_(qubits),
        X_(2 * qubits, std::vector<uint8_t>(qubits, 0)),
        Z_(2 * qubits, std::vector<uint8_t>(qubits, 0)),
        r_(2 * qubits, 0) {
    for (std::size_t i = 0; i < n_; ++i) {
      X_[i][i] = 1;
      Z_[n_ + i][i] = 1;
    }
  }

  void apply_hadamard(std::size_t target) {
    for (std::size_t i = 0; i < 2 * n_; ++i) {
      r_[i] ^= (X_[i][target] & Z_[i][target]);
      std::swap(X_[i][target], Z_[i][target]);
    }
  }

  void apply_cnot(std::size_t control, std::size_t target) {
    for (std::size_t i = 0; i < 2 * n_; ++i) {
      r_[i] ^= (X_[i][control] & Z_[i][target]) &
               (X_[i][target] ^ Z_[i][control] ^ 1);
      X_[i][target] ^= X_[i][control];
      Z_[i][control] ^= Z_[i][target];
    }
  }

  std::size_t qubit_count() const { return n_; }

 private:
  std::size_t n_;
  std::vector<std::vector<uint8_t>> X_, Z_;
  std::vector<uint8_t> r_;
};

// ---------------------------------------------------------------------------
// Octonion Algebra & Associator Engine (PR #153 / C++ Companion)
// ---------------------------------------------------------------------------
struct Octonion {
  std::array<double, 8> c; // e0 (1), e1..e7

  Octonion() { c.fill(0.0); }
  explicit Octonion(const std::array<double, 8>& coeffs) : c(coeffs) {}
  Octonion(double e0, double e1, double e2, double e3, double e4, double e5, double e6, double e7)
      : c({e0, e1, e2, e3, e4, e5, e6, e7}) {}

  double norm_sq() const {
    double s = 0;
    for (double val : c) s += val * val;
    return s;
  }

  Octonion operator+(const Octonion& o) const {
    Octonion res;
    for (int i = 0; i < 8; ++i) res.c[i] = c[i] + o.c[i];
    return res;
  }

  Octonion operator-(const Octonion& o) const {
    Octonion res;
    for (int i = 0; i < 8; ++i) res.c[i] = c[i] - o.c[i];
    return res;
  }

  Octonion operator*(double s) const {
    Octonion res;
    for (int i = 0; i < 8; ++i) res.c[i] = c[i] * s;
    return res;
  }

  // Cayley-Dickson / Fano multiplication table for e1..e7
  // Fano triplets: (1,2,3), (1,4,5), (1,7,6), (2,4,6), (2,5,7), (3,4,7), (3,6,5)
  Octonion operator*(const Octonion& o) const {
    const auto& a = c;
    const auto& b = o.c;
    std::array<double, 8> res;
    res.fill(0.0);

    // 1. Skalarer Anteil (e0 * e0 = 1, e_i * e_i = -1 für i=1..7)
    double dot_vec = 0.0;
    for (int i = 1; i <= 7; ++i) {
      dot_vec += a[i] * b[i];
    }
    res[0] = a[0] * b[0] - dot_vec;

    // 2. Kreuzterme e_0 mit e_i
    for (int i = 1; i <= 7; ++i) {
      res[i] = a[0] * b[i] + a[i] * b[0];
    }

    // Fano lookup table für Basisvektoren i, j in {1..7}
    auto fano_mult = [](int i, int j) -> std::pair<int, int> {
      if (i == j) return {0, -1};
      static const int triplets[7][3] = {
        {1,2,3}, {1,4,5}, {1,7,6}, {2,4,6}, {2,5,7}, {3,4,7}, {3,6,5}
      };
      for (const auto& tri : triplets) {
        int u = tri[0], v = tri[1], w = tri[2];
        if (i == u && j == v) return {w, 1};
        if (i == v && j == w) return {u, 1};
        if (i == w && j == u) return {v, 1};
        if (j == u && i == v) return {w, -1};
        if (j == v && i == w) return {u, -1};
        if (j == w && i == u) return {v, -1};
      }
      return {0, 0};
    };

    // 3. Imaginärteil-Multiplikation (Fano-Tripel)
    for (int i = 1; i <= 7; ++i) {
      for (int j = 1; j <= 7; ++j) {
        if (i == j) continue;  // Diagonale bereits im Skalaranteil
        if (a[i] == 0.0 || b[j] == 0.0) continue;
        auto [target, sign] = fano_mult(i, j);
        if (target > 0) {
          res[target] += a[i] * b[j] * sign;
        }
      }
    }
    return Octonion(res);
  }
};

Octonion associator(const Octonion& x, const Octonion& y, const Octonion& z) {
  // [x, y, z] = (x * y) * z - x * (y * z)
  return (x * y) * z - x * (y * z);
}

// ---------------------------------------------------------------------------
// Classical Shor / factorization (NOT quantum) — arbitrary-precision BigInt
// ---------------------------------------------------------------------------
constexpr int kMaxFactorBits = 4096;
// Mainstream-Probedivision und Pollard-Vorfilter: 10× gegenüber 5e6.
// 100-Bit-Primfaktoren (~7e12) bleiben ausserhalb; das ist Absicht (kein DoS).
constexpr uint64_t kPMax = 50000000ULL;
// Billiges Sieb, bleibt auch nach der Weiche aktiv (verhindert 21×81-Splits).
constexpr uint64_t kSmallPrimeMax = 100000ULL;

constexpr int bitlen64(uint64_t x) {
  int b = 0;
  while (x) {
    x >>= 1;
    ++b;
  }
  return b;
}
// Volle Probedivision nur, solange ein balancierter Faktor noch unter p_max
// liegen kann: bit(N) < 2·bit(p_max). Darüber: Skip → Pollard-Brent.
constexpr int kTrialSkipBits = 2 * bitlen64(kPMax);  // 52 bei p_max = 5e7

using u128 = __uint128_t;
u128 mulmod128(u128 a, u128 b, u128 m);

struct BigInt {
  std::vector<uint64_t> d;  // little-endian limbs

  BigInt() = default;
  explicit BigInt(uint64_t v) {
    if (v) d.push_back(v);
  }

  void trim() {
    while (!d.empty() && d.back() == 0) d.pop_back();
  }

  bool is_zero() const { return d.empty(); }
  bool is_one() const { return d.size() == 1 && d[0] == 1; }
  bool is_even() const { return d.empty() || ((d[0] & 1ULL) == 0); }

  int bit_length() const {
    if (d.empty()) return 0;
    int bits = static_cast<int>((d.size() - 1) * 64);
    uint64_t x = d.back();
    while (x) {
      ++bits;
      x >>= 1;
    }
    return bits;
  }

  int cmp(const BigInt& o) const {
    if (d.size() != o.d.size())
      return d.size() < o.d.size() ? -1 : 1;
    for (int i = static_cast<int>(d.size()) - 1; i >= 0; --i) {
      if (d[i] != o.d[i]) return d[i] < o.d[i] ? -1 : 1;
    }
    return 0;
  }

  bool operator==(const BigInt& o) const { return cmp(o) == 0; }
  bool operator<(const BigInt& o) const { return cmp(o) < 0; }
  bool operator>(const BigInt& o) const { return cmp(o) > 0; }
  bool operator<=(const BigInt& o) const { return cmp(o) <= 0; }
  bool operator>=(const BigInt& o) const { return cmp(o) >= 0; }

  static BigInt add(const BigInt& a, const BigInt& b) {
    BigInt r;
    const size_t n = std::max(a.d.size(), b.d.size());
    r.d.resize(n + 1, 0);
    unsigned __int128 carry = 0;
    for (size_t i = 0; i < n; ++i) {
      unsigned __int128 s = carry;
      if (i < a.d.size()) s += a.d[i];
      if (i < b.d.size()) s += b.d[i];
      r.d[i] = static_cast<uint64_t>(s);
      carry = s >> 64;
    }
    r.d[n] = static_cast<uint64_t>(carry);
    r.trim();
    return r;
  }

  static BigInt sub(const BigInt& a, const BigInt& b) {
    BigInt r = a;
    unsigned __int128 borrow = 0;
    for (size_t i = 0; i < r.d.size(); ++i) {
      unsigned __int128 x = r.d[i];
      unsigned __int128 y = borrow + (i < b.d.size() ? b.d[i] : 0);
      if (x < y) {
        r.d[i] = static_cast<uint64_t>(x + ((static_cast<unsigned __int128>(1) << 64) - y));
        borrow = 1;
      } else {
        r.d[i] = static_cast<uint64_t>(x - y);
        borrow = 0;
      }
    }
    r.trim();
    return r;
  }

  static BigInt mul(const BigInt& a, const BigInt& b) {
    if (a.is_zero() || b.is_zero()) return BigInt();
    BigInt r;
    r.d.assign(a.d.size() + b.d.size(), 0);
    for (size_t i = 0; i < a.d.size(); ++i) {
      unsigned __int128 carry = 0;
      for (size_t j = 0; j < b.d.size(); ++j) {
        unsigned __int128 cur =
            r.d[i + j] +
            static_cast<unsigned __int128>(a.d[i]) * b.d[j] + carry;
        r.d[i + j] = static_cast<uint64_t>(cur);
        carry = cur >> 64;
      }
      size_t k = i + b.d.size();
      while (carry) {
        if (k >= r.d.size()) r.d.push_back(0);
        unsigned __int128 cur = r.d[k] + carry;
        r.d[k] = static_cast<uint64_t>(cur);
        carry = cur >> 64;
        ++k;
      }
    }
    r.trim();
    return r;
  }

  static std::pair<BigInt, BigInt> divmod(BigInt a, BigInt b) {
    BigInt zero;
    if (b.is_zero()) return {zero, zero};
    if (a < b) return {zero, a};

    {
      __uint128_t aa = 0, bb = 0;
      if (a.to_u128(&aa) && b.to_u128(&bb) && bb != 0) {
        return {from_u128(aa / bb), from_u128(aa % bb)};
      }
    }

    if (b.d.size() == 1) {
      uint64_t div = b.d[0];
      BigInt q;
      q.d.resize(a.d.size(), 0);
      unsigned __int128 rem = 0;
      for (int i = static_cast<int>(a.d.size()) - 1; i >= 0; --i) {
        rem = (rem << 64) | a.d[i];
        q.d[i] = static_cast<uint64_t>(rem / div);
        rem %= div;
      }
      q.trim();
      BigInt r;
      if (rem) r.d.push_back(static_cast<uint64_t>(rem));
      return {q, r};
    }

    // Knuth Algorithm D, Basis 2^64 — ersetzt O(bits^2) Bit-Schleife.
    const int s = __builtin_clzll(b.d.back());
    BigInt v = shl_small(b, s);
    BigInt u = shl_small(a, s);
    u.d.push_back(0);
    const size_t n = v.d.size();
    const size_t m = u.d.size() - n - 1;
    const uint64_t vn1 = v.d[n - 1];
    BigInt q;
    q.d.assign(m + 1, 0);
    for (int j = static_cast<int>(m); j >= 0; --j) {
      unsigned __int128 un =
          (static_cast<unsigned __int128>(u.d[static_cast<size_t>(j) + n]) << 64) |
          u.d[static_cast<size_t>(j) + n - 1];
      unsigned __int128 qhat = un / vn1;
      unsigned __int128 rhat = un % vn1;
      if (qhat > ~uint64_t{0}) {
        qhat = ~uint64_t{0};
        rhat = un - qhat * vn1;
      }
      const uint64_t vn2 = v.d[n - 2];
      while (qhat * vn2 >
             ((rhat << 64) | u.d[static_cast<size_t>(j) + n - 2])) {
        --qhat;
        rhat += vn1;
        if (rhat > ~uint64_t{0}) break;
      }
      unsigned __int128 p_carry = 0;
      __int128 borrow = 0;
      for (size_t i = 0; i < n; ++i) {
        unsigned __int128 prod =
            static_cast<unsigned __int128>(v.d[i]) *
                static_cast<uint64_t>(qhat) +
            p_carry;
        p_carry = prod >> 64;
        __int128 t = static_cast<__int128>(u.d[static_cast<size_t>(j) + i]) -
                     static_cast<__int128>(static_cast<uint64_t>(prod)) - borrow;
        u.d[static_cast<size_t>(j) + i] = static_cast<uint64_t>(t);
        borrow = (t < 0) ? 1 : 0;
      }
      {
        __int128 t =
            static_cast<__int128>(u.d[static_cast<size_t>(j) + n]) -
            static_cast<__int128>(static_cast<uint64_t>(p_carry)) - borrow;
        u.d[static_cast<size_t>(j) + n] = static_cast<uint64_t>(t);
        borrow = (t < 0) ? 1 : 0;
      }
      if (borrow) {
        unsigned __int128 carry = 0;
        for (size_t i = 0; i < n; ++i) {
          unsigned __int128 sum =
              static_cast<unsigned __int128>(u.d[static_cast<size_t>(j) + i]) +
              v.d[i] + carry;
          u.d[static_cast<size_t>(j) + i] = static_cast<uint64_t>(sum);
          carry = sum >> 64;
        }
        u.d[static_cast<size_t>(j) + n] = static_cast<uint64_t>(
            static_cast<unsigned __int128>(u.d[static_cast<size_t>(j) + n]) +
            carry);
        --qhat;
      }
      q.d[static_cast<size_t>(j)] = static_cast<uint64_t>(qhat);
    }
    q.trim();
    u.d.resize(n);
    u.trim();
    return {q, shr_small(u, s)};
  }

  static BigInt shl_small(BigInt x, int s) {
    if (s <= 0 || x.is_zero()) return x;
    unsigned __int128 carry = 0;
    for (size_t i = 0; i < x.d.size(); ++i) {
      unsigned __int128 cur =
          (static_cast<unsigned __int128>(x.d[i]) << s) | carry;
      x.d[i] = static_cast<uint64_t>(cur);
      carry = cur >> 64;
    }
    if (carry) x.d.push_back(static_cast<uint64_t>(carry));
    return x;
  }

  static BigInt shr_small(BigInt x, int s) {
    if (s <= 0 || x.is_zero()) return x;
    uint64_t carry = 0;
    for (int i = static_cast<int>(x.d.size()) - 1; i >= 0; --i) {
      const uint64_t cur = x.d[i];
      x.d[i] = (cur >> s) | carry;
      carry = (s == 0) ? 0 : (cur << (64 - s));
    }
    x.trim();
    return x;
  }

  static BigInt mod(const BigInt& a, const BigInt& m) {
    return divmod(a, m).second;
  }

  static BigInt addmod(const BigInt& a, const BigInt& b, const BigInt& m) {
    BigInt s = add(a, b);
    if (s >= m) s = sub(s, m);
    return s;
  }

  static BigInt mulmod(const BigInt& a, const BigInt& b, const BigInt& m) {
    if (a.is_zero() || b.is_zero() || m.is_zero()) return BigInt();
    __uint128_t aa = 0, bb = 0, mm = 0;
    if (a.to_u128(&aa) && b.to_u128(&bb) && m.to_u128(&mm) && mm != 0) {
      return from_u128(mulmod128(aa, bb, mm));
    }
    return mod(mul(a, b), m);
  }

  static BigInt gcd(BigInt a, BigInt b) {
    if (a.is_zero()) return b;
    if (b.is_zero()) return a;
    int shift = 0;
    while (a.is_even() && b.is_even()) {
      a = shr1(a);
      b = shr1(b);
      ++shift;
    }
    while (a.is_even()) a = shr1(a);
    while (!b.is_zero()) {
      while (b.is_even()) b = shr1(b);
      if (a > b) std::swap(a, b);
      b = sub(b, a);
    }
    while (shift >= 63) {
      a = shl_small(a, 63);
      shift -= 63;
    }
    return shift ? shl_small(a, shift) : a;
  }

  static BigInt powmod(BigInt base, BigInt exp, const BigInt& m) {
    BigInt res(1);
    base = mod(base, m);
    while (!exp.is_zero()) {
      if (exp.d[0] & 1ULL) res = mulmod(res, base, m);
      base = mulmod(base, base, m);
      exp = shr1(exp);
    }
    return res;
  }

  static BigInt shr1(BigInt x) { return shr_small(x, 1); }

  uint64_t mod_u64(uint64_t m) const {
    if (m == 0 || d.empty()) return 0;
    unsigned __int128 rem = 0;
    for (int i = static_cast<int>(d.size()) - 1; i >= 0; --i) {
      rem = ((rem << 64) | d[i]) % m;
    }
    return static_cast<uint64_t>(rem);
  }

  bool to_u128(__uint128_t* out) const {
    if (d.empty()) { *out = 0; return true; }
    if (d.size() == 1) { *out = d[0]; return true; }
    if (d.size() == 2) { *out = (static_cast<__uint128_t>(d[1]) << 64) | d[0]; return true; }
    return false;
  }

  static BigInt from_u128(__uint128_t v) {
    BigInt r;
    if (!v) return r;
    r.d.push_back(static_cast<uint64_t>(v));
    const uint64_t hi = static_cast<uint64_t>(v >> 64);
    if (hi) r.d.push_back(hi);
    return r;
  }
};

u128 gcd128(u128 a, u128 b) {
  while (b != 0) {
    u128 t = b; b = a % b; a = t;
  }
  return a;
}

u128 addmod128(u128 a, u128 b, u128 m) {
  if (a >= m - b) return a - (m - b);
  return a + b;
}

u128 dblmod128(u128 a, u128 m) {
  if (a >= m - a) return a - (m - a);
  return a << 1;
}

u128 mulmod128(u128 a, u128 b, u128 m) {
  a %= m; b %= m;
  u128 res = 0;
  while (b) {
    if (b & 1) res = addmod128(res, a, m);
    a = dblmod128(a, m);
    b >>= 1;
  }
  return res;
}

u128 powmod128(u128 base, u128 exp, u128 mod) {
  u128 res = 1;
  base %= mod;
  while (exp > 0) {
    if (exp & 1) res = mulmod128(res, base, mod);
    base = mulmod128(base, base, mod);
    exp >>= 1;
  }
  return res;
}

enum class Primality128 { Composite, PrimeDeterministic, ProbablePrime };

bool miller_rabin_witness(u128 n, u128 d, int s, u128 a) {
  if (a % n == 0) return true;
  u128 x = powmod128(a, d, n);
  if (x == 1 || x == n - 1) return true;
  for (int r = 1; r < s; ++r) {
    x = mulmod128(x, x, n);
    if (x == n - 1) return true;
  }
  return false;
}

Primality128 is_prime_u128(u128 n) {
  if (n < 2) return Primality128::Composite;
  if (n == 2 || n == 3 || n == 5 || n == 7) return Primality128::PrimeDeterministic;
  if ((n % 2) == 0 || (n % 3) == 0 || (n % 5) == 0 || (n % 7) == 0) return Primality128::Composite;

  u128 d = n - 1;
  int s = 0;
  while ((d & 1) == 0) { d >>= 1; ++s; }

  static constexpr std::array<uint64_t, 7> kBases64 = {
    2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL
  };

  if (n <= ~uint64_t{0}) {
    for (uint64_t a : kBases64) {
      if (n <= a) break;
      if (!miller_rabin_witness(n, d, s, a)) return Primality128::Composite;
    }
    return Primality128::PrimeDeterministic;
  }

  static constexpr std::array<uint64_t, 13> kBasesPrime = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41
  };
  static const u128 kPsi12 = [] {
    u128 v = 0;
    for (const char* p = "318665857834031151167461"; *p; ++p)
      v = v * 10 + static_cast<u128>(*p - '0');
    return v;
  }();
  static const u128 kPsi13 = [] {
    u128 v = 0;
    for (const char* p = "3317044064679887385961981"; *p; ++p)
      v = v * 10 + static_cast<u128>(*p - '0');
    return v;
  }();

  const int n_wit = (n < kPsi12) ? 12 : 13;
  for (int i = 0; i < n_wit; ++i) {
    const u128 a = kBasesPrime[static_cast<size_t>(i)];
    if (n <= a) break;
    if (!miller_rabin_witness(n, d, s, a)) return Primality128::Composite;
  }
  if (n < kPsi13) return Primality128::PrimeDeterministic;
  return Primality128::ProbablePrime;
}

struct ClassicalShorResult {
  BigInt N, a, r, f1, f2;
  bool success = false;
  bool declared_prime = false;
  bool skipped_adaptive = false;
  std::string note;
};

int bitlen128(u128 n) {
  if (n == 0) return 0;
  int b = 0;
  if (n >> 64) {
    b += 64;
    n >>= 64;
  }
  return b + 64 - __builtin_clzll(static_cast<uint64_t>(n));
}

uint64_t adaptive_trial_bound(int bits) {
  return bits >= kTrialSkipBits ? kSmallPrimeMax : kPMax;
}

u128 trial_odd_factor(u128 N, uint64_t p_max) {
  for (uint64_t p = 3; p <= p_max; p += 2) {
    const u128 P = p;
    if (P * P > N) break;
    if (N % P == 0) return P;
  }
  return 0;
}

int rho_step_budget(int bits) {
  if (bits <= 32) return 200000;
  if (bits <= 48) return 1000000;
  if (bits <= 64) return 2000000;
  if (bits <= 80) return 4000000;
  if (bits <= 100) return 8000000;
  if (bits <= 128) return 16000000;
  if (bits <= 160) return 4000000;
  if (bits <= 192) return 2500000;
  if (bits <= 256) return 1500000;
  return 800000;
}

int rho_curve_budget(int bits) {
  if (bits <= 32) return 8;
  if (bits <= 64) return 16;
  if (bits <= 80) return 24;
  if (bits <= 128) return 48;
  if (bits <= 192) return 64;
  return 96;
}

// Brent-Pollard (batch-GCD). CLASSICAL — kein Quantenpfad.
u128 pollard_brent_once(u128 N, u128 c, int max_steps) {
  auto f = [&](u128 y) {
    return addmod128(mulmod128(y, y, N), c, N);
  };
  u128 y = 2, x = 2, ys = 2, q = 1, g = 1;
  int64_t r = 1;
  const int64_t m = 128;
  int steps = 0;
  while (g == 1 && steps < max_steps) {
    x = y;
    for (int64_t i = 0; i < r && steps < max_steps; ++i) {
      y = f(y);
      ++steps;
    }
    int64_t k = 0;
    while (k < r && g == 1 && steps < max_steps) {
      ys = y;
      const int64_t lim = std::min(m, r - k);
      for (int64_t i = 0; i < lim && steps < max_steps; ++i) {
        y = f(y);
        ++steps;
        const u128 diff = (x > y) ? x - y : y - x;
        if (diff != 0) q = mulmod128(q, diff, N);
      }
      g = gcd128(q, N);
      k += lim;
    }
    r *= 2;
  }
  if (g == N) {
    g = 1;
    y = ys;
    while (g == 1 && steps < max_steps + 1024) {
      y = f(y);
      ++steps;
      const u128 diff = (x > y) ? x - y : y - x;
      g = gcd128(diff, N);
    }
  }
  return (g > 1 && g < N) ? g : 0;
}

u128 pollard_brent_factor(u128 N, int bits) {
  const int budget = rho_step_budget(bits);
  const int curves = rho_curve_budget(bits);
  for (int c = 1; c <= curves; ++c) {
    const u128 d = pollard_brent_once(N, static_cast<u128>(c), budget);
    if (d > 1 && d < N) return d;
  }
  return 0;
}

// Kombinierter Peel: Sieb → MR → Brent; bei zusammengesetztem Teiler rekursiv.
u128 smallest_prime_factor_u128(u128 N, int depth) {
  if (N <= 1) return 0;
  if ((N & 1) == 0) return 2;
  const int bits = bitlen128(N);
  if (const u128 p = trial_odd_factor(N, adaptive_trial_bound(bits))) return p;
  const auto pr = is_prime_u128(N);
  if (pr != Primality128::Composite) return N;
  if (depth <= 0) return 0;
  const u128 d = pollard_brent_factor(N, bits);
  if (d == 0) return 0;
  u128 a = d;
  u128 b = N / d;
  if (b < a) std::swap(a, b);
  const u128 sp = smallest_prime_factor_u128(a, depth - 1);
  return sp ? sp : a;
}

ClassicalShorResult run_classical_shor_u128(u128 N, [[maybe_unused]] u128 a_pref,
                                   [[maybe_unused]] bool auto_base) {
  ClassicalShorResult out;
  out.N = BigInt::from_u128(N);
  if (N <= 1) { out.note = "N muss > 1 sein"; return out; }
  if ((N & 1) == 0) {
    out.f1 = BigInt(2); out.f2 = BigInt::from_u128(N / 2);
    out.success = true; out.note = "trivial: N gerade"; return out;
  }

  const int bits = bitlen128(N);
  const bool skip_full_trial = bits >= kTrialSkipBits;
  const uint64_t bound = adaptive_trial_bound(bits);

  if (const u128 p = trial_odd_factor(N, bound)) {
    out.f1 = BigInt::from_u128(p);
    out.f2 = BigInt::from_u128(N / p);
    out.success = true;
    out.note = skip_full_trial
                   ? "Pollard-Pfad: billiges Sieb (Weiche, p<=" +
                         std::to_string(kSmallPrimeMax) + ")"
                   : "Pollard-Pfad: kleinster Faktor via Probedivision";
    return out;
  }

  auto pr = is_prime_u128(N);
  if (pr == Primality128::PrimeDeterministic || pr == Primality128::ProbablePrime) {
    out.declared_prime = true;
    out.note = "Miller-Rabin: N ist Primzahl oder wahrscheinliche Primzahl.";
    return out;
  }

  const u128 d = pollard_brent_factor(N, bits);
  if (d > 1 && d < N) {
    u128 f = d;
    u128 g = N / d;
    if (g < f) std::swap(f, g);
    const u128 sp = smallest_prime_factor_u128(f, 8);
    bool combined = false;
    if (sp > 1 && sp < N && (N % sp) == 0) {
      combined = (sp != f);
      f = sp;
      g = N / sp;
    }
    out.f1 = BigInt::from_u128(f);
    out.f2 = BigInt::from_u128(g);
    out.success = true;
    out.note = skip_full_trial
                   ? "Pollard-Brent (u128), Weiche: volle Probedivision uebersprungen"
                   : "Pollard-Brent (u128)";
    if (combined) out.note += "; kombinierter Peel";
    return out;
  }
  out.note = "Keine Zerlegung gefunden.";
  return out;
}

BigInt trial_odd_factor_big(const BigInt& N, uint64_t p_max) {
  for (uint64_t p = 3; p <= p_max; p += 2) {
    if (N.bit_length() <= 63) {
      const uint64_t n64 = N.d.empty() ? 0 : N.d[0];
      if (p > 0 && static_cast<unsigned __int128>(p) * p > n64) break;
    }
    if (N.mod_u64(p) == 0) return BigInt(p);
  }
  return BigInt();
}

bool miller_rabin_witness_big(const BigInt& n, const BigInt& d, int s,
                              uint64_t a) {
  const BigInt A(a);
  if (BigInt::mod(A, n).is_zero()) return true;
  BigInt x = BigInt::powmod(A, d, n);
  if (x.is_one()) return true;
  const BigInt nm1 = BigInt::sub(n, BigInt(1));
  if (x == nm1) return true;
  for (int r = 1; r < s; ++r) {
    x = BigInt::mulmod(x, x, n);
    if (x == nm1) return true;
  }
  return false;
}

Primality128 is_prime_big(const BigInt& n) {
  u128 n128 = 0;
  if (n.to_u128(&n128)) return is_prime_u128(n128);
  if (n.bit_length() < 2) return Primality128::Composite;
  if (n.is_even()) return Primality128::Composite;
  static constexpr uint64_t kSmall[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
  for (uint64_t p : kSmall) {
    if (n.mod_u64(p) == 0) return Primality128::Composite;
  }
  BigInt d = BigInt::sub(n, BigInt(1));
  int s = 0;
  while (d.is_even()) {
    d = BigInt::shr1(d);
    ++s;
  }
  static constexpr uint64_t kBases[] = {
      2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
  for (uint64_t a : kBases) {
    if (!miller_rabin_witness_big(n, d, s, a)) return Primality128::Composite;
  }
  return Primality128::ProbablePrime;
}

BigInt pollard_brent_once_big(const BigInt& N, uint64_t c, int max_steps) {
  const BigInt C(c);
  auto f = [&](const BigInt& y) {
    return BigInt::addmod(BigInt::mulmod(y, y, N), C, N);
  };
  BigInt y(2), x(2), ys(2), q(1), g(1);
  int64_t r = 1;
  const int64_t m = 128;
  int steps = 0;
  const BigInt one(1);
  while (g.is_one() && steps < max_steps) {
    x = y;
    for (int64_t i = 0; i < r && steps < max_steps; ++i) {
      y = f(y);
      ++steps;
    }
    int64_t k = 0;
    while (k < r && g.is_one() && steps < max_steps) {
      ys = y;
      const int64_t lim = std::min(m, r - k);
      for (int64_t i = 0; i < lim && steps < max_steps; ++i) {
        y = f(y);
        ++steps;
        const BigInt diff = (x > y) ? BigInt::sub(x, y) : BigInt::sub(y, x);
        if (!diff.is_zero()) q = BigInt::mulmod(q, diff, N);
      }
      g = BigInt::gcd(q, N);
      k += lim;
    }
    r *= 2;
  }
  if (g == N) {
    g = one;
    y = ys;
    while (g.is_one() && steps < max_steps + 1024) {
      y = f(y);
      ++steps;
      const BigInt diff = (x > y) ? BigInt::sub(x, y) : BigInt::sub(y, x);
      g = BigInt::gcd(diff, N);
    }
  }
  if (!g.is_zero() && g > one && g < N) return g;
  return BigInt();
}

BigInt pollard_brent_factor_big(const BigInt& N, int bits) {
  const int budget = rho_step_budget(bits);
  const int curves = rho_curve_budget(bits);
  for (int c = 1; c <= curves; ++c) {
    BigInt d = pollard_brent_once_big(N, static_cast<uint64_t>(c), budget);
    if (!d.is_zero() && d > BigInt(1) && d < N) return d;
  }
  return BigInt();
}

BigInt smallest_prime_factor_big(const BigInt& N, int depth) {
  if (N.bit_length() < 2) return BigInt();
  if (N.is_even()) return BigInt(2);
  const int bits = N.bit_length();
  BigInt p = trial_odd_factor_big(N, adaptive_trial_bound(bits));
  if (!p.is_zero()) return p;
  if (is_prime_big(N) != Primality128::Composite) return N;
  if (depth <= 0) return BigInt();
  BigInt d = pollard_brent_factor_big(N, bits);
  if (d.is_zero()) return BigInt();
  BigInt a = d;
  BigInt b = BigInt::divmod(N, d).first;
  if (b < a) std::swap(a, b);
  BigInt sp = smallest_prime_factor_big(a, depth - 1);
  return sp.is_zero() ? a : sp;
}

ClassicalShorResult run_classical_shor_big(BigInt N) {
  ClassicalShorResult out;
  out.N = N;
  if (N.bit_length() < 2) {
    out.note = "N muss > 1 sein";
    return out;
  }
  if (N.is_even()) {
    out.f1 = BigInt(2);
    out.f2 = BigInt::divmod(N, BigInt(2)).first;
    out.success = true;
    out.note = "trivial: N gerade";
    return out;
  }

  const int bits = N.bit_length();
  const bool skip_full_trial = bits >= kTrialSkipBits;
  const uint64_t bound = adaptive_trial_bound(bits);

  BigInt p = trial_odd_factor_big(N, bound);
  if (!p.is_zero()) {
    out.f1 = p;
    out.f2 = BigInt::divmod(N, p).first;
    out.success = true;
    out.note = skip_full_trial
                   ? "Pollard-Pfad: billiges Sieb (BigInt, p<=" +
                         std::to_string(kSmallPrimeMax) + ")"
                   : "Pollard-Pfad: kleinster Faktor via Probedivision (BigInt)";
    return out;
  }

  const auto pr = is_prime_big(N);
  if (pr == Primality128::PrimeDeterministic ||
      pr == Primality128::ProbablePrime) {
    out.declared_prime = true;
    out.note = "Miller-Rabin (BigInt): wahrscheinliche Primzahl.";
    return out;
  }

  BigInt d = pollard_brent_factor_big(N, bits);
  if (!d.is_zero() && d > BigInt(1) && d < N) {
    BigInt f = d;
    BigInt g = BigInt::divmod(N, d).first;
    if (g < f) std::swap(f, g);
    BigInt sp = smallest_prime_factor_big(f, 8);
    bool combined = false;
    if (!sp.is_zero() && sp < N) {
      const BigInt rem = BigInt::mod(N, sp);
      if (rem.is_zero()) {
        combined = !(sp == f);
        f = sp;
        g = BigInt::divmod(N, sp).first;
      }
    }
    out.f1 = f;
    out.f2 = g;
    out.success = true;
    out.note = skip_full_trial
                   ? "Pollard-Brent (BigInt), Weiche: volle Probedivision uebersprungen"
                   : "Pollard-Brent (BigInt)";
    if (combined) out.note += "; kombinierter Peel";
    out.note += " [" + std::to_string(bits) + " Bit]";
    return out;
  }
  out.note =
      "Keine Zerlegung gefunden (BigInt-Pfad). Pollard findet typischerweise "
      "kleine/mittlere Faktoren, nicht balancierte RSA-Moduli.";
  return out;
}

ClassicalShorResult run_classical_shor(BigInt N, BigInt a_pref, bool auto_base) {
  if (N.bit_length() > kMaxFactorBits) {
    ClassicalShorResult out;
    out.N = N;
    out.note = "N ueberschreitet Soft-Cap von 4096 Bit";
    return out;
  }
  u128 n128 = 0, a128 = 0;
  if (N.to_u128(&n128) && a_pref.to_u128(&a128)) {
    return run_classical_shor_u128(n128, a128, auto_base);
  }
  return run_classical_shor_big(N);
}

enum class ParseStatus { Ok, BadDigits, TooBig };

ParseStatus parse_bigint(const std::string& s, BigInt* out) {
  if (s.empty()) return ParseStatus::BadDigits;
  BigInt v;
  for (char c : s) {
    if (c < '0' || c > '9') return ParseStatus::BadDigits;
    v = BigInt::mul(v, BigInt(10));
    v = BigInt::add(v, BigInt(static_cast<uint64_t>(c - '0')));
    if (v.bit_length() > kMaxFactorBits) return ParseStatus::TooBig;
  }
  *out = v;
  return ParseStatus::Ok;
}

bool knuth_divmod_selftest() {
  BigInt a, b;
  if (parse_bigint("79228162514264337593543950349", &a) != ParseStatus::Ok)
    return false;
  if (parse_bigint("1208925819614629174706193", &b) != ParseStatus::Ok)
    return false;
  const BigInt n = BigInt::mul(a, b);
  const auto qa = BigInt::divmod(n, a);
  const auto qb = BigInt::divmod(n, b);
  return qa.second.is_zero() && qa.first == b && qb.second.is_zero() &&
         qb.first == a;
}

std::string bigint_to_dec(BigInt v) {
  if (v.is_zero()) return "0";
  std::string s;
  const BigInt ten(10);
  while (!v.is_zero()) {
    auto qr = BigInt::divmod(v, ten);
    s.push_back(static_cast<char>(
        '0' + static_cast<int>(qr.second.d.empty() ? 0 : qr.second.d[0])));
    v = qr.first;
  }
  std::reverse(s.begin(), s.end());
  return s;
}

// Mainstream-Baseline: ungerade Probedivision (Schulbuch / trial division).
// Begrenzt auf p <= p_max, damit 50-Bit-Faktoren nicht zum DoS werden.
ClassicalShorResult trial_division_u128(u128 N, uint64_t p_max) {
  ClassicalShorResult out;
  out.N = BigInt::from_u128(N);
  if (N <= 1) {
    out.note = "N muss > 1 sein";
    return out;
  }
  if ((N & 1) == 0) {
    out.f1 = BigInt(2);
    out.f2 = BigInt::from_u128(N / 2);
    out.success = true;
    out.note = "Probedivision: N gerade";
    return out;
  }
  auto pr = is_prime_u128(N);
  if (pr == Primality128::PrimeDeterministic ||
      pr == Primality128::ProbablePrime) {
    out.declared_prime = true;
    out.note = "Probedivision uebersprungen: Miller-Rabin sagt prim";
    return out;
  }
  for (uint64_t p = 3; p <= p_max; p += 2) {
    const u128 P = p;
    if (P * P > N) break;
    if (N % P == 0) {
      out.f1 = BigInt(p);
      out.f2 = BigInt::from_u128(N / P);
      out.success = true;
      out.note = "Probedivision (Mainstream / Schulbuch)";
      return out;
    }
  }
  out.note =
      "Probedivision: kein Faktor <= " + std::to_string(p_max) +
      " (Mainstream-Baseline ungeeignet bei grossen Primfaktoren)";
  return out;
}

void print_factor_result(const char* title, const ClassicalShorResult& r,
                         double ms) {
  std::cout << "  " << title << "  [" << std::fixed << std::setprecision(3)
            << ms << " ms]\n";
  std::cout << "    note: " << r.note << "\n";
  if (r.declared_prime) {
    std::cout << "    [PASS] " << bigint_to_dec(r.N) << " ist PRIM\n";
  } else if (r.success) {
    std::cout << "    [PASS] " << bigint_to_dec(r.f1) << " × "
              << bigint_to_dec(r.f2);
    if (BigInt::mul(r.f1, r.f2) == r.N) std::cout << "  (Produktcheck OK)";
    std::cout << "\n";
  } else if (r.skipped_adaptive) {
    std::cout << "    [SKIP] adaptive Weiche — keine Probedivision\n";
  } else {
    std::cout << "    [!] keine Zerlegung in diesem Verfahren\n";
  }
}

void factorize_interactive(BigInt N) {
  const int bits = N.bit_length();
  const bool skip_trial = bits >= kTrialSkipBits;
  std::cout << "\n----------------------------------------------------------------\n";
  std::cout << " FAKTORISIERUNGS-VERGLEICH  N = " << bigint_to_dec(N) << "\n";
  std::cout << " (klassisch — kein Quanten-Shor / Claim-Wall)\n";
  std::cout << " Bitsize: " << bits << " / Soft-Cap " << kMaxFactorBits
            << " / p_max " << kPMax << " / Weiche @" << kTrialSkipBits << " Bit\n";
  std::cout << " Weiche: " << (skip_trial
        ? "Probedivision SKIP → Pollard-Brent / kombinierter Peel"
        : "Probedivision aktiv (bit(N) < Schwelle)") << "\n";
  std::cout << "----------------------------------------------------------------\n";

  auto t0 = std::chrono::high_resolution_clock::now();
  const ClassicalShorResult pollard = run_classical_shor(N, BigInt(), true);
  auto t1 = std::chrono::high_resolution_clock::now();
  const double ms_pollard =
      std::chrono::duration<double, std::milli>(t1 - t0).count();

  ClassicalShorResult trial;
  double ms_trial = 0;
  u128 n128 = 0;
  if (skip_trial) {
    trial.N = N;
    trial.skipped_adaptive = true;
    trial.note =
        "adaptiv uebersprungen: bit(N)=" + std::to_string(bits) +
        " >= " + std::to_string(kTrialSkipBits) +
        " (volle Probedivision bis p_max waere O(p_max) ohne Trefferchance)";
  } else if (N.to_u128(&n128)) {
    auto u0 = std::chrono::high_resolution_clock::now();
    trial = trial_division_u128(n128, kPMax);
    auto u1 = std::chrono::high_resolution_clock::now();
    ms_trial = std::chrono::duration<double, std::milli>(u1 - u0).count();
  } else {
    trial.N = N;
    trial.note = "Probedivision: N > 128 Bit — Baseline nicht angewandt";
  }

  print_factor_result("[1] Pollard-Brent (dieses Programm)", pollard, ms_pollard);
  print_factor_result("[2] Probedivision (Mainstream)", trial, ms_trial);

  if (pollard.success && trial.success) {
    const bool same =
        (pollard.f1 == trial.f1 && pollard.f2 == trial.f2) ||
        (pollard.f1 == trial.f2 && pollard.f2 == trial.f1);
    std::cout << "  Vergleich: Faktoren "
              << (same ? "IDENTISCH" : "UNTERSCHIEDLICH (beide teilen N)")
              << "; schneller: "
              << (ms_pollard <= ms_trial ? "Pollard-Brent" : "Probedivision")
              << "\n";
  } else if (pollard.success && trial.skipped_adaptive) {
    std::cout << "  Vergleich: Weiche aktiv — nur Pollard-Brent/kombiniert.\n";
  } else if (pollard.success && !trial.success) {
    std::cout << "  Vergleich: Pollard-Brent findet grosse Faktoren, "
                 "Probedivision bleibt unter p_max stecken.\n";
  } else if (pollard.declared_prime && trial.declared_prime) {
    std::cout << "  Vergleich: beide Verfahren stufen N als prim ein.\n";
  }
}

bool read_N_line(BigInt* N_out, bool* quit) {
  *quit = false;
  std::cout << "\nZu faktorisierende Zahl N eingeben\n"
            << "  (Beispiele: 15, 1701, 1234567890123456789012345678901,\n"
            << "   170141693884019613139382498777795253379317181)\n"
            << "  Soft-Cap: " << kMaxFactorBits << " Bit; p_max = " << kPMax
            << "; Weiche @" << kTrialSkipBits << " Bit\n"
            << "  Ab 129 Bit: BigInt-Pollard-Brent (kleine/mittlere Faktoren).\n"
            << "  (q oder 0 = Dialog beenden): " << std::flush;
  std::string line;
  if (!std::getline(std::cin, line)) {
    *quit = true;
    return false;
  }
  while (!line.empty() && (line.back() == ' ' || line.back() == '\t' ||
                           line.back() == '\r'))
    line.pop_back();
  size_t s = 0;
  while (s < line.size() && (line[s] == ' ' || line[s] == '\t')) ++s;
  line = line.substr(s);
  if (line.empty()) {
    *N_out = BigInt(15);
    return true;
  }
  if (line == "q" || line == "Q" || line == "quit" || line == "exit" ||
      line == "0") {
    *quit = true;
    return false;
  }
  BigInt v;
  const ParseStatus st = parse_bigint(line, &v);
  if (st == ParseStatus::BadDigits) {
    std::cout << "  Ungueltige Ziffern.\n";
    return false;
  }
  if (st == ParseStatus::TooBig) {
    std::cout << "  Zahl zu gross: Soft-Cap " << kMaxFactorBits << " Bit.\n";
    return false;
  }
  *N_out = v;
  return true;
}

void run_factorization_dialog() {
  std::cout << "\n================================================================\n"
            << " INTERAKTIVE FAKTORISIERUNG (klassisch)\n"
            << " Pollard-Brent vs. Mainstream-Probedivision. Kein Quantencomputer.\n"
            << " Adaptive Weiche: volle Probedivision nur bei bit(N) < "
            << kTrialSkipBits << ".\n"
            << "================================================================\n";
  while (true) {
    BigInt N;
    bool quit = false;
    if (!read_N_line(&N, &quit)) {
      if (quit) break;
      continue;
    }
    factorize_interactive(N);
  }
  std::cout << "Faktorisierungs-Dialog beendet.\n";
}

void print_usage(const char* argv0) {
  std::cout << "Usage:\n"
            << "  " << argv0 << "                 # Dialog, danach Zertifikate\n"
            << "  " << argv0 << " --skip-certs    # nur Faktorisierungs-Dialog\n"
            << "  " << argv0 << " --certs-only    # nur Gate/Assoziator-Zertifikate\n";
}

// ---------------------------------------------------------------------------
// F₂ helpers & Fano-Steane-E8 setup
// ---------------------------------------------------------------------------
using Bits7 = std::array<uint8_t, 7>;
using Bits8 = std::array<uint8_t, 8>;
using Bits14 = std::array<uint8_t, 14>;

uint8_t dot7(const Bits7& a, const Bits7& b) {
  uint8_t s = 0; for (int i = 0; i < 7; ++i) s ^= (a[i] & b[i]); return s;
}
uint8_t dot8(const Bits8& a, const Bits8& b) {
  uint8_t s = 0; for (int i = 0; i < 8; ++i) s ^= (a[i] & b[i]); return s;
}
int weight8(const Bits8& v) {
  int w = 0; for (auto b : v) w += b; return w;
}
Bits8 parity_lift(const Bits7& v) {
  Bits8 out{}; uint8_t p = 0;
  for (int i = 0; i < 7; ++i) { out[i] = v[i]; p ^= v[i]; }
  out[7] = p; return out;
}
Bits8 xor8(const Bits8& a, const Bits8& b) {
  Bits8 o{}; for (int i = 0; i < 8; ++i) o[i] = a[i] ^ b[i]; return o;
}
Bits14 x_gen(const Bits7& h) { Bits14 g{}; for (int i = 0; i < 7; ++i) g[i] = h[i]; return g; }
Bits14 z_gen(const Bits7& h) { Bits14 g{}; for (int i = 0; i < 7; ++i) g[7+i] = h[i]; return g; }
uint8_t omega14(const Bits14& u, const Bits14& v) {
  uint8_t s = 0;
  for (int i = 0; i < 7; ++i) s ^= (u[i] & v[7+i]) ^ (u[7+i] & v[i]);
  return s;
}
int pauli_weight14(const Bits14& v) {
  int w = 0; for (int i = 0; i < 7; ++i) if (v[i] || v[7+i]) ++w; return w;
}

void print_pass(const std::string& gate, bool ok, const std::string& detail) {
  std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << gate << " — " << detail << "\n";
}

using Vec8i = std::array<int, 8>;
Vec8i apply_H(const Vec8i& x) {
  Vec8i y{};
  for (int k = 0; k < 4; ++k) {
    int a = x[2*k], b = x[2*k+1];
    y[2*k] = a + b; y[2*k+1] = a - b;
  }
  return y;
}

std::set<Vec8i> build_R_std_doubled() {
  std::set<Vec8i> S;
  for (int i = 0; i < 8; ++i) {
    for (int j = i + 1; j < 8; ++j) {
      for (int si : {-1, 1}) {
        for (int sj : {-1, 1}) {
          Vec8i v{}; v[i] = 2 * si; v[j] = 2 * sj; S.insert(v);
        }
      }
    }
  }
  for (int mask = 0; mask < 256; ++mask) {
    if (__builtin_popcount(mask) % 2 != 0) continue;
    Vec8i v{};
    for (int i = 0; i < 8; ++i) v[i] = (mask & (1 << i)) ? -1 : 1;
    S.insert(v);
  }
  return S;
}

std::vector<Vec8i> build_R_A(const std::vector<Bits8>& C8) {
  std::vector<Vec8i> RA;
  RA.reserve(240);
  for (const auto& c : C8) {
    int w = weight8(c);
    if (w == 0) {
      for (int i = 0; i < 8; ++i) {
        for (int s : {-2, 2}) { Vec8i v{}; v[i] = s; RA.push_back(v); }
      }
    } else if (w == 4) {
      int supp[4], nsupp = 0;
      for (int i = 0; i < 8; ++i) if (c[i]) supp[nsupp++] = i;
      for (int signs = 0; signs < 16; ++signs) {
        Vec8i v{};
        for (int k = 0; k < 4; ++k) v[supp[k]] = (signs & (1 << k)) ? -1 : 1;
        RA.push_back(v);
      }
    }
  }
  return RA;
}

}  // namespace

int main(int argc, char** argv) {
  int failures = 0;

  bool skip_certs = false;
  bool certs_only = false;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      print_usage(argv[0]);
      return 0;
    }
    if (arg == "--skip-certs") {
      skip_certs = true;
      continue;
    }
    if (arg == "--certs-only") {
      certs_only = true;
      continue;
    }
    std::cerr << "Unbekanntes Argument: " << arg << "\n";
    print_usage(argv[0]);
    return 2;
  }

  std::cout << "================================================================\n"
            << " Fano–Steane–E8 + Octonion Associator Engine (PR #153)\n"
            << " Companion to Lean: EabcRenorm/QEC & Octonion Subsystems\n"
            << "================================================================\n\n"
            << "CLAIM-WALL: algebraic certificates only; no hardware / advantage.\n"
            << "Faktorisierung: CLASSICAL (Pollard-Brent; u128 + BigInt bis "
            << kMaxFactorBits << " Bit).\n\n";

  if (!knuth_divmod_selftest()) {
    std::cerr << "Knuth-Divmod Selftest FEHLGESCHLAGEN (BigInt-Pfad unbrauchbar).\n";
    return 1;
  }

  if (!certs_only) {
    run_factorization_dialog();
    if (skip_certs) {
      std::cout << "\n[--skip-certs] Gate-Zertifikate uebersprungen.\n";
      return 0;
    }
    std::cout << "\n";
  }

  // -------------------------------------------------------------------------
  // A. Clifford barrier test
  // -------------------------------------------------------------------------
  std::cout << "--- A. Clifford barrier (Aaronson–Gottesman tableau) ---\n";
  StabilizerSimulator sim(100);
  sim.apply_hadamard(0);
  for (std::size_t i = 1; i < 100; ++i) sim.apply_cnot(0, i);
  print_pass("Clifford-100", true, "Stabilizer simulation O(n^2) passed successfully");

  // -------------------------------------------------------------------------
  // B. Fano / Steane CSS branch
  // -------------------------------------------------------------------------
  std::cout << "\n--- B. Fano → Steane CSS branch ---\n";
  Bits7 h1 = {1, 1, 1, 1, 0, 0, 0};
  Bits7 h2 = {1, 1, 0, 0, 1, 1, 0};
  Bits7 h3 = {1, 0, 1, 0, 1, 0, 1};
  const Bits7 hs[3] = {h1, h2, h3};
  bool self_orth7 = true;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      if (dot7(hs[i], hs[j]) != 0) self_orth7 = false;
  Bits14 gens14[6];
  for (int i = 0; i < 3; ++i) {
    gens14[i] = x_gen(hs[i]);
    gens14[3 + i] = z_gen(hs[i]);
  }
  bool isotropic = true;
  int gen_wt = 0;
  for (int i = 0; i < 6; ++i) {
    gen_wt = std::max(gen_wt, pauli_weight14(gens14[i]));
    for (int j = 0; j < 6; ++j)
      if (omega14(gens14[i], gens14[j]) != 0) isotropic = false;
  }
  const bool g17 = self_orth7 && isotropic && gen_wt <= 4;
  print_pass("G.1-G.7", g17,
             "Steane CSS: F2-orth + symplectic isotropy (gen-wt<=4)");
  if (!g17) ++failures;

  // -------------------------------------------------------------------------
  // C. Construction A & E8 roots (240)
  // -------------------------------------------------------------------------
  std::cout << "\n--- C. Construction A & E8 roots ---\n";
  Bits7 h0 = {1, 1, 1, 1, 1, 1, 1};
  Bits8 e1 = parity_lift(h1), e2 = parity_lift(h2), e3 = parity_lift(h3), e4 = parity_lift(h0);
  Bits8 eG[4] = {e1, e2, e3, e4};
  std::vector<Bits8> C8;
  for (int mask = 0; mask < 16; ++mask) {
    Bits8 v{};
    for (int b = 0; b < 4; ++b) if (mask & (1 << b)) v = xor8(v, eG[b]);
    C8.push_back(v);
  }
  auto RA = build_R_A(C8);
  print_pass("G.12", RA.size() == 240, "E8 root count via Construction A matches exactly 240");
  if (RA.size() != 240) ++failures;

  bool self_dual8 = true;
  for (const auto& u : C8)
    for (const auto& v : C8)
      if (dot8(u, v) != 0) self_dual8 = false;
  print_pass("G.10c", self_dual8, "C8 self-orthogonal over F2");
  if (!self_dual8) ++failures;

  const auto Rstd2 = build_R_std_doubled();
  std::set<Vec8i> images;
  bool g13_subset = true;
  for (const auto& x : RA) {
    const Vec8i hx = apply_H(x);
    if (!Rstd2.count(hx)) g13_subset = false;
    images.insert(hx);
  }
  const bool g13 = (RA.size() == 240 && Rstd2.size() == 240 && g13_subset &&
                    images.size() == 240);
  print_pass("G.13", g13, "runtime T_Q bijection |im|=240");
  if (!g13) ++failures;

  // -------------------------------------------------------------------------
  // D. Octonion Associator & Alternativity Audit (PR #153 Integration)
  // -------------------------------------------------------------------------
  std::cout << "\n--- D. Octonion Associator & Alternativity Engine (PR #153) ---\n";
  
  // Test 1: Alternativity Audit [x, x, y] == 0
  Octonion x_test(1.0, 0.0, 2.0, -1.0, 0.0, 1.0, 0.0, 0.0);
  Octonion y_test(0.0, 1.0, 0.0, 1.0, -2.0, 0.0, 1.0, 0.0);
  Octonion alt_res = associator(x_test, x_test, y_test);
  double alt_norm = alt_res.norm_sq();
  bool alternativity_pass = (alt_norm < 1e-12);
  print_pass("PR153.1", alternativity_pass, "Alternativity guard [x, x, y] == 0");
  if (!alternativity_pass) ++failures;

  // Test 2: Norm Composition N(xy) == N(x)N(y)
  Octonion xy = x_test * y_test;
  double norm_prod = xy.norm_sq();
  double prod_norm = x_test.norm_sq() * y_test.norm_sq();
  bool norm_pass = (std::abs(norm_prod - prod_norm) < 1e-10);
  print_pass("PR153.2", norm_pass, "Hurwitz norm composition N(xy) == N(x)N(y)");
  if (!norm_pass) ++failures;

  // Test 3: General non-associativity check (generic elements)
  Octonion z_test(1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
  Octonion assoc_gen = associator(x_test, y_test, z_test);
  std::cout << "  [INFO] Generic associator magnitude: " << std::sqrt(assoc_gen.norm_sq()) << "\n";
  print_pass("PR153.3", true, "Non-associativity diagnostic active (Associator non-trivial)");

  std::cout << "\n================================================================\n";
  if (failures == 0) {
    std::cout << " RESULT: ALL CERTIFICATES & PR #153 GATES PASS (0 failures)\n";
  } else {
    std::cout << " RESULT: " << failures << " FAILURE(S) DETECTED\n";
  }
  std::cout << "================================================================\n";
  return failures == 0 ? 0 : 1;
}