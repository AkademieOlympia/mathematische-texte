import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# -----------------------------
# Parameter
# -----------------------------

N = 1000                # Anzahl Primzahlen
ZEROS_FILE = "zeros6.npy"

# -----------------------------
# Primzahlen
# -----------------------------

print("generating primes...")

primes = list(sp.primerange(2, 50000))[:N]  # 50000 liefert >4000 Primzahlen
if len(primes) < N:
    N = len(primes)
    print(f"Hinweis: Nur {N} Primzahlen verfügbar.")
primes = primes[:N]
p = np.array(primes, dtype=float)

logp = np.log(p)

# -----------------------------
# BM Matrix
# -----------------------------

print("building BM operator...")

Q = np.zeros((N, N))

for i in range(N):
    for j in range(N):
        if i != j:
            Q[i,j] = np.log(p[i]*p[j]) / np.sqrt(p[i]*p[j])

D = np.block([
    [np.zeros((N,N)), Q],
    [Q.T, np.zeros((N,N))]
])

# -----------------------------
# Eigenwerte
# -----------------------------

print("computing spectrum...")

eigvals = np.linalg.eigvals(D)
spec = np.sort(np.abs(eigvals))
spec = spec[:N]

# -----------------------------
# Riemann Nullstellen
# -----------------------------

print("loading zeros...")

zeros = np.load(ZEROS_FILE)
gammas = zeros[:N]

# -----------------------------
# Spektrumvergleich
# -----------------------------

corr = np.corrcoef(spec, gammas)[0,1]

plt.figure()
plt.scatter(gammas, spec, s=10)
plt.xlabel("Riemann zero ordinate γ")
plt.ylabel("BM eigenvalue |λ|")
plt.title("BM spectrum vs Riemann zeros")
plt.savefig("bm_spectrum_vs_riemann.png", dpi=300)

print("correlation =", corr)

# -----------------------------
# Spacing Statistik
# -----------------------------

s = np.diff(spec)
s_norm = s / np.mean(s)

plt.figure()
plt.hist(s_norm, bins=40, density=True)
plt.xlabel("normalized spacing")
plt.ylabel("density")
plt.title("BM spacing distribution")
plt.savefig("bm_spacing.png", dpi=300)

# -----------------------------
# Unfolded Spectrum
# -----------------------------

x = np.arange(len(spec))
coeff = np.polyfit(x, spec, 3)
smooth = np.polyval(coeff, x)

unfolded = spec / np.mean(np.diff(smooth))

s_un = np.diff(unfolded)
s_un = s_un / np.mean(s_un)

plt.figure()
plt.hist(s_un, bins=40, density=True)
plt.xlabel("unfolded spacing")
plt.ylabel("density")
plt.title("BM unfolded spectrum")
plt.savefig("bm_unfolded.png", dpi=300)

# -----------------------------
# Trace Test
# -----------------------------

tvals = np.linspace(0.01, 2, 150)

trace_vals = []
prime_vals = []

for t in tvals:

    trace_vals.append(
        np.sum(np.exp(-t * spec**2))
    )

    prime_vals.append(
        np.sum(np.log(p) * np.exp(-t*(np.log(p))**2))
    )

trace_vals = np.array(trace_vals)
prime_vals = np.array(prime_vals)

trace_vals /= np.max(trace_vals)
prime_vals /= np.max(prime_vals)

plt.figure()
plt.plot(tvals, trace_vals, label="spectral trace")
plt.plot(tvals, prime_vals, label="prime sum")
plt.xlabel("t")
plt.ylabel("normalized value")
plt.title("Trace comparison")
plt.legend()
plt.savefig("bm_trace.png", dpi=300)

print("done.")