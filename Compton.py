from sage.all import *
import numpy as np
import mpmath as mp
import math

# -----------------------------
# FIXED HARD SETTINGS
# -----------------------------
P_MAX  = 850000
alpha  = 1.0/137.0
T_total   = 40000       # <- guter Start
burn_frac = 0.30
M_per_family = 800      # <- schneller als 2000
R_runs = 10             # Kugeln

np.random.seed(1)

# -----------------------------
# Riemann target (mod 2pi)
# -----------------------------
mp.mp.dps = 50
gamma_850k = float(mp.im(mp.zetazero(850000)))
gamma_mod  = gamma_850k % (2.0*math.pi)

def wrap_pm_pi(x):
    return (x + math.pi)%(2.0*math.pi) - math.pi

# -----------------------------
# Families e,a,b,c = mod12 1,5,7,11
# -----------------------------
fam_list = ['e','a','b','c']
residue  = {'e':1,'a':5,'b':7,'c':11}

# Get primes and split by family
pr = list(primes(P_MAX))
fam_primes = {f: [] for f in fam_list}

for p in pr:
    if p in (2,3):
        continue
    r = p % 12
    if r == 1:   fam_primes['e'].append(p)
    elif r == 5: fam_primes['a'].append(p)
    elif r == 7: fam_primes['b'].append(p)
    elif r == 11:fam_primes['c'].append(p)

Nfam = {f: len(fam_primes[f]) for f in fam_list}

# Precompute omega distributions log(p) for sampling
omega_list = {f: np.log(np.array(fam_primes[f], dtype=np.float64)) for f in fam_list}

# Coupling K_ij(alpha): +1 same, -alpha different
K = np.ones((4,4), dtype=np.float64)
for i in range(4):
    for j in range(4):
        if i != j:
            K[i,j] = -alpha

# Map family index
idx = {'e':0,'a':1,'b':2,'c':3}

def one_run(seed):
    rng = np.random.default_rng(seed)

    # sample phases + omegas per family
    phi = {}
    omg = {}
    for f in fam_list:
        phi[f] = 2.0*math.pi*rng.random(M_per_family)
        w = omega_list[f]
        omg[f] = w[rng.integers(0, len(w), size=M_per_family)]

    u = np.zeros(T_total, dtype=np.float64)

    for t in range(T_total):
        # order parameters
        Z = {}
        for f in fam_list:
            Z[f] = np.mean(np.exp(1j*phi[f]))

        # u(t) ~ sum N_f Re(Z_f)
        u[t] = sum(Nfam[f]*np.real(Z[f]) for f in fam_list)

        # update (vectorized)
        for f in fam_list:
            i = idx[f]
            # precompute sin/cos once
            s = np.sin(phi[f])
            c = np.cos(phi[f])

            # force = sum_g N_g K_ig Im(Z_g e^{-i phi})
            # Im(Z e^{-i phi}) = Re(Z)*sin(phi) - Im(Z)*cos(phi)
            force = 0.0
            for g in fam_list:
                j = idx[g]
                ReZ = np.real(Z[g])
                ImZ = np.imag(Z[g])
                force += Nfam[g] * K[i,j] * (ReZ*s - ImZ*c)

            phi[f] = (phi[f] + omg[f] + alpha*force) % (2.0*math.pi)

    # FFT(u)
    burn = int(burn_frac*T_total)
    u2 = u[burn:] - np.mean(u[burn:])

    # mild window to reduce leakage
    u2 = u2 * np.hanning(len(u2))

    fft_vals = np.fft.rfft(u2)
    power = np.abs(fft_vals)**2
    freqs = np.fft.rfftfreq(len(u2), d=1.0)

    power[0] = 0.0
    kmax = int(np.argmax(power))

    f_peak = float(freqs[kmax])
    Omega_peak = 2.0*math.pi*f_peak
    delta_mod = wrap_pm_pi(Omega_peak - gamma_mod)

    return f_peak, Omega_peak, delta_mod

# -----------------------------
# Run many Kugeln
# -----------------------------
peaks = []
deltas = []

for r in range(R_runs):
    fpk, Om, d = one_run(seed=100+r)
    peaks.append(fpk)
    deltas.append(abs(d))

peaks = np.array(peaks)
deltas = np.array(deltas)

print("=== 850k Kugel Mean-Field FFT(u) Riemann check ===")
print("alpha =", alpha)
print("gamma_850000 =", gamma_850k)
print("gamma_mod_2pi =", gamma_mod)
print()
print("Runs:", R_runs, "T_total:", T_total, "M_per_family:", M_per_family)
print("Mean f_peak:", peaks.mean(), "Std:", peaks.std())
print("Mean |delta_mod|:", deltas.mean())
print("Min  |delta_mod|:", deltas.min())
print("Hits |delta_mod|<0.2:", int(np.sum(deltas < 0.2)))
print("Hits |delta_mod|<0.1:", int(np.sum(deltas < 0.1)))
print("Hits |delta_mod|<0.05:", int(np.sum(deltas < 0.05)))
