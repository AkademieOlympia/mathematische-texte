import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange

# -----------------------
# 1) Mod12-Familien & Maske
# -----------------------
def family_mod12(p):
    r = p % 12
    if r == 1:  return "E"
    if r == 5:  return "A"
    if r == 7:  return "B"
    if r == 11: return "C"
    return None

def build_ab_mask_abc_after_137(primes, start_value=137):
    mask = np.zeros(len(primes), dtype=float)
    for i, p in enumerate(primes):
        if p < start_value:
            continue
        fam = family_mod12(p)
        if fam in ("A","B","C"):
            mask[i] = 1.0
    return mask

# -----------------------
# 2) Periodogramm
# -----------------------
def periodogram(x):
    """
    Einfache FFT-Power. Entfernt Mittelwert (DC) und gibt Frequenzen/Power zurück.
    Frequenzen in Zyklen pro Sample (0..0.5).
    """
    x = np.asarray(x, float)
    x = x - x.mean()
    n = len(x)
    X = np.fft.rfft(x)
    P = (np.abs(X)**2) / n
    freqs = np.fft.rfftfreq(n, d=1.0)
    return freqs, P

def top_peaks(freqs, P, k=12, fmin=1e-5):
    # ignoriere DC und ganz niedrige f
    idx = np.where(freqs > fmin)[0]
    j = idx[np.argsort(P[idx])[::-1][:k]]
    peaks = []
    for t in j:
        f = freqs[t]
        L = 1.0 / f if f > 0 else np.inf
        peaks.append((P[t], f, L))
    return peaks

# -----------------------
# 3) Zieltest: Energie nahe Periode L0
# -----------------------
def band_power_near_L(freqs, P, L0, rel_width=0.01):
    """
    Macht einen "Bandpower"-Test um f0=1/L0 herum.
    rel_width=0.01 => ±1% Band um f0.
    """
    f0 = 1.0 / L0
    lo, hi = f0*(1-rel_width), f0*(1+rel_width)
    m = (freqs >= lo) & (freqs <= hi)
    return float(P[m].sum())

def shuffle_mask(mask, rng, block=None):
    """
    Permutationstest:
    - block=None: komplette Permutation (maximal zerstörend)
    - block=int: Block-Permutation (erhält lokale Cluster)
    """
    mask = np.asarray(mask, float)
    n = len(mask)
    if block is None or block <= 1:
        perm = rng.permutation(n)
        return mask[perm]
    idx = np.arange(n)
    blocks = [idx[i:i+block] for i in range(0, n, block)]
    rng.shuffle(blocks)
    perm = np.concatenate(blocks)
    return mask[perm]

def bandpower_pvalue(mask, L0, n_perm=200, rel_width=0.01, seed=0, block=None):
    freqs, P = periodogram(mask)
    obs = band_power_near_L(freqs, P, L0, rel_width=rel_width)

    rng = np.random.default_rng(seed)
    sims = np.empty(n_perm, float)
    for k in range(n_perm):
        m2 = shuffle_mask(mask, rng, block=block)
        f2, P2 = periodogram(m2)
        sims[k] = band_power_near_L(f2, P2, L0, rel_width=rel_width)

    p = (np.sum(sims >= obs) + 1) / (n_perm + 1)
    return obs, float(p), sims

# -----------------------
# 4) Run
# -----------------------
limit = 100000
primes = list(primerange(1, limit))
mask = build_ab_mask_abc_after_137(primes, start_value=137)

freqs, P = periodogram(mask)
peaks = top_peaks(freqs, P, k=15)

print("Top-Peaks der MASKE (Power, freq, Periode~1/f):")
for power, f, L in peaks:
    print(f"  power={power:.4e}  f={f:.6e}  L~{L:.2f}")

# explizit 137 und 274 testen
for L0 in (120,137,240,274):
    obs, p, sims = bandpower_pvalue(mask, L0, n_perm=200, rel_width=0.01, seed=L0, block=200)
    print(f"\nBandpower-Test Maske @ L={L0}: obs={obs:.4e}, p={p:.4f}  (Block=200)")

# Plot Periodogramm (optional)
plt.figure(figsize=(12,5))
plt.plot(freqs[1:], P[1:], linewidth=1.0)
plt.title("Periodogramm der AB-Maske (ABC-only, ab p≥137)")
plt.xlabel("Frequenz (Zyklen pro Schritt)")
plt.ylabel("Power")
plt.grid(True, alpha=0.25)
plt.tight_layout()
plt.show()
