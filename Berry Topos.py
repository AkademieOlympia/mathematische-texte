from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SCRIPT_DIR = Path(__file__).resolve().parent
_CSV_NAME = "Vier_bis_4700_mit_Familien.csv"
_csv_candidates = [
    Path.home() / "Downloads" / _CSV_NAME,
    _SCRIPT_DIR / _CSV_NAME,
]
_csv_path = next((p for p in _csv_candidates if p.is_file()), None)
if _csv_path is None:
    raise FileNotFoundError(
        f"{_CSV_NAME} nicht gefunden (erwartet in Downloads oder {_SCRIPT_DIR})."
    )

_zeros_path = _SCRIPT_DIR / "zeros6.npy"
if not _zeros_path.is_file():
    raise FileNotFoundError(
        f"zeros6.npy fehlt: {_zeros_path} — Datei ins Projektverzeichnis legen."
    )

# --- Daten ---
quad_df = pd.read_csv(_csv_path)
zeros = np.load(_zeros_path)

quadruples = quad_df.iloc[:60, :4].values.astype(int)

# --- Knoten ---
primes = sorted(set(quadruples.flatten()))
n = len(primes)


def state(p):
    """Zustandsvektor im 5D-Hilbertraum (Features aus ln p und √p)."""
    x = np.log(p)
    return np.array(
        [
            np.cos(x),
            np.sin(x),
            np.cos(2 * x),
            np.sin(2 * x),
            np.cos(np.sqrt(p)),
        ]
    )


states = np.array([state(p) for p in primes])
states = states / np.linalg.norm(states, axis=1, keepdims=True)

# --- Quantum-Metrik M (reell, symmetrisch) ---
M = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i != j:
            overlap = np.dot(states[i], states[j])
            d = abs(np.log(primes[i]) - np.log(primes[j]))
            if d < 1.0:
                w = 1.0
            else:
                w = 0.0
            M[i, j] = w * (1.0 - overlap**2)

M += 1e-6 * np.eye(n)

# --- Berry-Term F (schiefsymmetrisch) ---
F = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i != j:
            overlap = np.dot(states[i], states[j])
            d = abs(np.log(primes[i]) - np.log(primes[j]))
            w = np.exp(-d)
            F[i, j] = w * np.sin(overlap)

# antisymmetrisch machen
F = F - F.T

alpha = 0.2
# H ist hermitesch: M symmetrisch, F schief → (iαF)† = iαF
H = M + 1j * alpha * F

# --- Eigenwerte (komplex → real extrahieren)
eigvals = np.linalg.eigvals(H)
eigvals = np.real(eigvals)

# sortieren
eigvals = np.sort(eigvals)

# Filter (robust!)
vals = eigvals[np.abs(eigvals) > 1e-6]

print("Eigenwerte:", len(vals))

# --- Spacings ---
sp = np.diff(vals)

if len(sp) > 10:
    sp = sp / np.mean(sp)
    print("Var:", np.var(sp))
else:
    print("ZU WENIG DATEN!")

# --- Zeta ---
zsp = np.diff(zeros[:2000])
zsp = zsp / np.mean(zsp)

# --- GUE ---
s = np.linspace(0, 3, 200)
wigner = (32 / (np.pi**2)) * s**2 * np.exp(-4 * s**2 / np.pi)

plt.figure(figsize=(8, 5))
if len(sp) > 10:
    plt.hist(
        sp,
        bins=20,
        density=True,
        alpha=0.5,
        label=rf"$\lambda(H)$, $M_{{ij}}=1-\langle i|j\rangle^2$, $F_{{ij}}=\langle i|j\rangle$, $\alpha={alpha}$",
    )
else:
    plt.text(
        0.5,
        0.5,
        "λ(H): zu wenig Spacings für Histogramm",
        ha="center",
        va="center",
        transform=plt.gca().transAxes,
    )
plt.hist(zsp, bins=20, density=True, alpha=0.5, label="Zeta")
plt.plot(s, wigner, label="GUE", linewidth=2)

plt.legend()
plt.xlabel("Spacing")
plt.ylabel("Dichte")
plt.title(
    r"Berry Topos — $|\psi_p\rangle$ in $\mathbb{R}^5$; "
    r"$H=M+\mathrm{i}\alpha F$ hermitesch; Spacings $\lambda(H)$"
)

if vals.size:
    print(
        f"λ(H) gefiltert: min={vals.min():.4g}, max={vals.max():.4g}, "
        f"Spacings={len(sp)}"
    )
print("Zeta Var:", np.var(zsp))

plt.show()
