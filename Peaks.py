import numpy as np
from sympy import primerange

def family_mod12(p):
    r = p % 12
    if r == 1:  return "E"
    if r == 5:  return "A"
    if r == 7:  return "B"
    if r == 11: return "C"
    return None

def build_mask(primes, start_value=137):
    # ABC-only ab p>=137
    m = np.zeros(len(primes), dtype=float)
    for i,p in enumerate(primes):
        if p < start_value:
            continue
        fam = family_mod12(p)
        if fam in ("A","B","C"):
            m[i] = 1.0
    return m

def periodogram(x):
    x = np.asarray(x, float)
    x = x - x.mean()
    n = len(x)
    X = np.fft.rfft(x)
    P = (np.abs(X)**2) / n
    freqs = np.fft.rfftfreq(n, d=1.0)
    return freqs, P

def bandpower(freqs, P, L0, rel_width=0.01):
    f0 = 1.0 / L0
    lo, hi = f0*(1-rel_width), f0*(1+rel_width)
    m = (freqs >= lo) & (freqs <= hi)
    return float(P[m].sum())

def block_shuffle(x, block, rng):
    n = len(x)
    idx = np.arange(n)
    blocks = [idx[i:i+block] for i in range(0, n, block)]
    rng.shuffle(blocks)
    perm = np.concatenate(blocks)
    return x[perm]

def pvalue_bandpower_mask(mask, L0, n_perm=300, block=200, rel_width=0.01, seed=0):
    freqs, P = periodogram(mask)
    obs = bandpower(freqs, P, L0, rel_width=rel_width)

    rng = np.random.default_rng(seed)
    sims = np.empty(n_perm, float)
    for k in range(n_perm):
        m2 = block_shuffle(mask, block, rng)
        f2, P2 = periodogram(m2)
        sims[k] = bandpower(f2, P2, L0, rel_width=rel_width)

    p = (np.sum(sims >= obs) + 1) / (n_perm + 1)
    return obs, float(p)

if __name__ == "__main__":
    limit = 100000
    primes = list(primerange(1, limit))
    mask = build_mask(primes, start_value=137)

    print(f"mask density (mean): {mask.mean():.4f}  | N={len(mask)}")

    for L0 in (120, 137, 240, 274):
        obs, p = pvalue_bandpower_mask(mask, L0, n_perm=300, block=200, rel_width=0.01, seed=L0)
        print(f"L={L0:>3}: bandpower={obs:.4e}  p={p:.4f}")
