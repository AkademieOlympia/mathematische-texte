import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Load zeros (fast)
# -----------------------------
zeros = np.load("zeros6.npy")  # ggf. "/mnt/data/zeros6.npy"
t_all = zeros.astype(np.float64)

# -----------------------------
# Optional: restrict to a T-range (recommended for speed & stability)
# -----------------------------
Tmin, Tmax = 1e4, 1e6
mask = (t_all >= Tmin) & (t_all <= Tmax)
t = t_all[mask]
print("N zeros:", len(t), "range:", (t[0], t[-1]))

# -----------------------------
# Smooth counting function Nbar(T) (Riemann-von-Mangoldt main term)
# -----------------------------
twopi = 2.0 * np.pi

def Nbar(T):
    T = np.asarray(T, dtype=np.float64)
    x = T / twopi
    return x * np.log(x) - x + 7.0/8.0

u = Nbar(t)  # unfolded levels

# spacings (unfolded)
s = np.diff(u)
s = s / np.mean(s)  # robust normalization
print("mean spacing:", np.mean(s), "std:", np.std(s))

# -----------------------------
# Windowed RG / "arithmetical HR diagram"
# -----------------------------
m = 50000          # window size in number of zeros (try 5k, 20k, 50k)
step = 10000       # shift in zeros between windows (overlap allowed)
qs = [0.1, 0.5, 0.9]

# optional form factor taus
taus = [0.2, 0.5, 1.0]  # a few points; keep small for speed

Xs = []
VarS = []
Quant = {q: [] for q in qs}
Ktau = {tau: [] for tau in taus}

# We use spacings s_n inside the window: s[idx:idx+m-1]
# and the corresponding t-window: t[idx:idx+m]
for idx in range(0, len(t) - m, step):
    t_win = t[idx:idx+m]
    u_win = u[idx:idx+m]
    s_win = s[idx:idx+m-1]

    Tj = np.median(t_win)
    Xj = np.log(Tj / twopi)

    Xs.append(Xj)
    VarS.append(np.var(s_win))

    for q in qs:
        Quant[q].append(np.quantile(s_win, q))

    # Simple form factor on unfolded levels (window-centered)
    u0 = u_win - np.mean(u_win)
    Nw = len(u0)
    for tau in taus:
        S_tau = np.exp(2j * np.pi * tau * u0).sum()
        Ktau[tau].append((np.abs(S_tau)**2) / Nw)

Xs = np.array(Xs)
VarS = np.array(VarS)
for q in qs:
    Quant[q] = np.array(Quant[q])
for tau in taus:
    Ktau[tau] = np.array(Ktau[tau])

# -----------------------------
# Plot: "Main sequence" = Var(spacings) vs log(T/2pi)
# -----------------------------
plt.figure()
plt.plot(Xs, VarS, marker=".", linestyle="none")
plt.xlabel(r"$X=\log(T/2\pi)$")
plt.ylabel(r"$Y=\mathrm{Var}(s)$ (unfolded)")
plt.title("Arithmetisches HR/RG-Diagramm: Fluktuationsstärke vs Skala")
plt.show()

# Plot: quantile bands
plt.figure()
for q in qs:
    plt.plot(Xs, Quant[q], marker=".", linestyle="none", label=f"Q{int(q*100)}")
plt.xlabel(r"$X=\log(T/2\pi)$")
plt.ylabel("Quantile der unfolded spacings s")
plt.title("Quantile-Bänder der Spacings (Fenster-RG)")
plt.legend()
plt.show()

# Plot: form factor points
plt.figure()
for tau in taus:
    plt.plot(Xs, Ktau[tau], marker=".", linestyle="none", label=f"tau={tau}")
plt.xlabel(r"$X=\log(T/2\pi)$")
plt.ylabel(r"$K(\tau)$")
plt.title("Fenster-Formfaktor (unfolded levels)")
plt.legend()
plt.show()