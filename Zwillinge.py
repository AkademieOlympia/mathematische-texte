import numpy as np
import math
import time

# ----------------------------
# Config
# ----------------------------
x = 10_000_000
h_list = [2, 4, 6, 10, 12, 18, 30, 210]  # frei anpassen
H = max(h_list)
N = x + H

# Twin-prime constant C2 (Hardy–Littlewood)
C2 = 0.6601618158468696  # ausreichend genau für numerische Tests

# Mod-12 Familien (deine Notation)
# E:1, A:5, B:7, C:11 (mod 12)
classes = {
    "E": 1,
    "A": 5,
    "B": 7,
    "C": 11,
}

# ----------------------------
# Sieve primes up to N
# ----------------------------
t0 = time.time()
is_prime = np.ones(N + 1, dtype=bool)
is_prime[:2] = False
limit = int(math.isqrt(N))
for p in range(2, limit + 1):
    if is_prime[p]:
        is_prime[p*p:N+1:p] = False
primes = np.flatnonzero(is_prime)
t1 = time.time()
print(f"[sieve] N={N:,} primes={len(primes):,} time={t1-t0:.2f}s")

# ----------------------------
# Build von Mangoldt arrays:
# Lambda(n) = log p if n = p^k (k>=1), else 0
# plus channel Lambdas by base prime class
# ----------------------------
Lam = np.zeros(N + 1, dtype=np.float32)
Lam_cls = {k: np.zeros(N + 1, dtype=np.float32) for k in classes}

# fill primes and prime powers
for p in primes:
    lp = math.log(p)
    # determine class by p mod 12 (only meaningful for odd primes > 3)
    pmod = p % 12
    cls_key = None
    for k, r in classes.items():
        if pmod == r:
            cls_key = k
            break

    # mark p^k
    pp = p
    while pp <= N:
        Lam[pp] = lp
        if cls_key is not None:
            Lam_cls[cls_key][pp] = lp
        # next power; guard overflow
        if pp > N // p:
            break
        pp *= p

t2 = time.time()
print(f"[Lambda] built time={t2-t1:.2f}s")

# ----------------------------
# Hardy–Littlewood singular series S(h)
# For even h:
#   S(h) = 2*C2 * Π_{p|h, p>2} (p-1)/(p-2)
# For odd h: S(h)=0
# ----------------------------
def prime_factors_unique(n: int):
    fac = []
    d = 2
    while d*d <= n:
        if n % d == 0:
            fac.append(d)
            while n % d == 0:
                n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        fac.append(n)
    return fac

def singular_series(h: int):
    if h % 2 == 1:
        return 0.0
    fac = prime_factors_unique(abs(h))
    prod = 1.0
    for p in fac:
        if p == 2:
            continue
        prod *= (p - 1.0) / (p - 2.0)
    return 2.0 * C2 * prod

# ----------------------------
# Vectorized correlations
# C_h(x) = sum_{n<=x} Lam[n] Lam[n+h]
# And channel versions: O_h^{XY}(x) = sum Lam_X[n] Lam_Y[n+h]
# ----------------------------
lam0 = Lam[1:x+1].astype(np.float64)  # float64 for stable summation

# marginals for baseline (channel densities)
A_cls = {}
for k in classes:
    A_cls[k] = float(np.sum(Lam_cls[k][1:x+1], dtype=np.float64))

A_all = float(np.sum(Lam[1:x+1], dtype=np.float64))
print(f"[marginals] sum Lambda ~ {A_all:.3e}")
for k in classes:
    print(f"  sum Lambda_{k} ~ {A_cls[k]:.3e}")

def corr(L1, L2, h):
    # sum_{n<=x} L1[n] * L2[n+h]
    return float(np.sum(L1 * L2[1+h:x+1+h], dtype=np.float64))

t3 = time.time()

print("\n=== Global C_h(x) ===")
for h in h_list:
    ch = corr(lam0, Lam.astype(np.float64), h)
    Sh = singular_series(h)
    # normalized by x
    ch_over_x = ch / x
    # fixpoint ratio vs HL prediction: should be ~ S(h) (roughly), so (ch/x)/S(h) ~ 1
    ratio = (ch_over_x / Sh) if Sh > 0 else float("nan")
    print(f"h={h:>4}  C_h/x={ch_over_x:.6f}  S(h)={Sh:.6f}  (C_h/x)/S(h)={ratio:.4f}")

print("\n=== Channel Cooper operators O_h^{XY}(x) ===")
pairs = [("E","A"), ("E","B"), ("E","C"), ("A","B"), ("A","C"), ("B","C")]
for h in h_list:
    print(f"\n-- h={h} --")
    for X,Y in pairs:
        LX = Lam_cls[X].astype(np.float64)
        LY = Lam_cls[Y].astype(np.float64)
        val = corr(LX[1:x+1], LY, h) / x

        # independence baseline from marginals:
        # E[ sum L_X(n) L_Y(n+h) ] ~ (sum L_X)(sum L_Y)/x
        base = (A_cls[X] * A_cls[Y]) / (x * x)

        # "excess" factor
        exc = (val / base) if base > 0 else float("nan")
        print(f"{X}{Y}: O_h/x={val:.6e}  baseline={base:.6e}  factor={exc:.4f}")

t4 = time.time()
print(f"\n[done] correlations time={t4-t3:.2f}s total={t4-t0:.2f}s")