from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

# --- Daten laden ---
quad_df = pd.read_csv("Vier_bis_4700_mit_Familien.csv")
zeros = np.load("zeros6.npy")

# --- kleine Teilmenge (wichtig!) ---
quadruples = quad_df.iloc[:40, :4].values.astype(int)

# --- Graph ---
primes = sorted(set(quadruples.flatten()))
index = {p:i for i,p in enumerate(primes)}

edges = []
for quad in quadruples:
    for i in range(4):
        for j in range(i+1,4):
            edges.append((index[quad[i]], index[quad[j]]))

n = len(primes)
m = len(edges)

# --- Inzidenzmatrix ---
B = np.zeros((n,m))
for k,(u,v) in enumerate(edges):
    B[u,k] = 1
    B[v,k] = -1

# --- BB^T Trick ---
M = B @ B.T
eigvals = np.linalg.eigvalsh(M)

# --- Dirac-Eigenwerte ---
vals = np.sqrt(np.maximum(eigvals, 0))
vals = np.sort(vals[vals > 1e-8])

# --- Spacings ---
sp = np.diff(vals)
sp = sp / np.mean(sp)

# --- Zeta ---
zsp = np.diff(zeros[:2000])
zsp = zsp / np.mean(zsp)

# --- Plot ---
s = np.linspace(0,3,200)
wigner = (32/(np.pi**2))*s**2*np.exp(-4*s**2/np.pi)

plt.figure(figsize=(8,5))
plt.hist(sp, bins=20, density=True, alpha=0.5, label="Dirac (Vierlinge)")
plt.hist(zsp, bins=20, density=True, alpha=0.5, label="Zeta")
plt.plot(s, wigner, label="GUE", linewidth=2)

plt.legend()
plt.xlabel("Spacing")
plt.ylabel("Dichte")
plt.title("Erster globaler Dirac-Test")

plt.show()

# --- Varianzen ---
print("Dirac Var:", np.var(sp))
print("Zeta Var :", np.var(zsp))