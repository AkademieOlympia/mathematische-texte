import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Load zeros (fast)
# ----------------------------
zeros = np.load("zeros6.npy")  # oder "/mnt/data/zeros6.npy" im Sandbox-Kontext
t = zeros.astype(np.float64)
# t sind die Imaginärteile der Nullstellen auf der kritischen Geraden (sortiert)

# Optional: nur ein Bereich (macht Tests schneller/sauberer)
# z.B. t in [Tmin, Tmax]
Tmin, Tmax = 1e4, 2e5
mask = (t >= Tmin) & (t <= Tmax)
t = t[mask]
print("N zeros in window:", len(t), "range:", (t[0], t[-1]))

# ----------------------------
# Unfolding via mean counting function N(T)
# Riemann–von Mangoldt main term:
# N(T) ≈ (T/2π) log(T/2π) - (T/2π) + 7/8
# (ignores S(T) fluctuations; good enough for spacing/statistics)
# ----------------------------
twopi = 2.0 * np.pi

def N_smooth(T):
    T = np.asarray(T, dtype=np.float64)
    x = T / twopi
    return x * np.log(x) - x + 7.0/8.0

u = N_smooth(t)  # unfolded "levels"
# spacings on unfolded scale:
s = np.diff(u)
# normalize to mean 1 (should already be ~1, but do it for robustness)
s = s / np.mean(s)
print("mean spacing:", np.mean(s), "std:", np.std(s))

# ----------------------------
# GUE Wigner surmise PDF for spacing:
# p(s) = (32/pi^2) s^2 exp(-4 s^2 / pi)
# We'll compare histogram + a quick distance metric.
# ----------------------------
def p_gue(s):
    s = np.asarray(s, dtype=np.float64)
    return (32.0 / (np.pi**2)) * (s**2) * np.exp(-(4.0/np.pi) * s**2)

# Histogram comparison
bins = np.linspace(0, 4, 81)
hist, edges = np.histogram(s, bins=bins, density=True)
centers = 0.5 * (edges[:-1] + edges[1:])
gue_pdf = p_gue(centers)

# L1 distance over bins (quick diagnostic)
dx = edges[1] - edges[0]
L1 = np.sum(np.abs(hist - gue_pdf)) * dx
print("L1 distance to GUE (lower is better):", L1)

plt.figure()
plt.plot(centers, hist, label="data (NN spacing)")
plt.plot(centers, gue_pdf, label="GUE Wigner surmise")
plt.xlabel("s (unfolded spacing, mean=1)")
plt.ylabel("density")
plt.title("Nearest-neighbor spacing vs GUE")
plt.legend()
plt.show()

# ----------------------------
# Simple spectral form factor on unfolded levels
# Define:
# S(tau) = sum_n exp(2π i tau * u_n)
# K(tau) = |S(tau)|^2 / N
# We also window u to reduce end effects by subtracting mean
# ----------------------------
u0 = u - np.mean(u)
N = len(u0)

taus = np.linspace(0.0, 5.0, 400)  # adjust range/resolution
K = np.empty_like(taus)

# To speed up: compute in blocks if N is huge
for j, tau in enumerate(taus):
    ph = np.exp(2j * np.pi * tau * u0)
    S_tau = ph.sum()
    K[j] = (np.abs(S_tau)**2) / N

plt.figure()
plt.plot(taus, K)
plt.xlabel("tau")
plt.ylabel("K(tau)")
plt.title("Spectral form factor (unfolded levels)")
plt.show()

# ----------------------------
# Hooks for your resonance tests (C1/C2) if you have F(t)
# Provide F(t) on some grid t_grid, then:
# - C1: distance peaks -> nearest zero
# - C2: hit mass in windows around zeros
# ----------------------------
def nearest_zero(t_query):
    i = np.searchsorted(zeros, t_query, side="left")
    if i <= 0:
        return zeros[0]
    if i >= len(zeros):
        return zeros[-1]
    left, right = zeros[i-1], zeros[i]
    return left if (t_query-left) <= (right-t_query) else right

def c1_peak_to_zero(peaks_t):
    peaks_t = np.asarray(peaks_t, dtype=np.float64)
    dists = np.empty_like(peaks_t)
    nz = np.empty_like(peaks_t)
    for k, tp in enumerate(peaks_t):
        z = nearest_zero(float(tp))
        nz[k] = z
        dists[k] = abs(tp - z)
    return dists, nz

def c2_hit_mass(t_grid, F, delta):
    t_grid = np.asarray(t_grid, dtype=np.float64)
    F = np.asarray(F, dtype=np.float64)
    total = np.trapz(F, t_grid)
    if total <= 0:
        return 0.0
    inside = np.zeros(t_grid.size, dtype=bool)
    for i, tt in enumerate(t_grid):
        lo = np.searchsorted(zeros, tt - delta, side="left")
        hi = np.searchsorted(zeros, tt + delta, side="right")
        inside[i] = (hi > lo)
    hit = np.trapz(F * inside.astype(np.float64), t_grid)
    return float(hit / total)

print("Done.")