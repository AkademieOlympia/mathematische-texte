import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from sympy import factorint, primerange


def build_H(m: np.ndarray, K: int = 200, scale: float = 1.0) -> np.ndarray:
    """
    Symmetrische Matrix H mit H_ij = cos(y / scale), y = |m_i - m_j| / 2
    (halber Abstand der Mittelpunkte), auf den ersten ``min(K, len(m))`` Einträgen von m.

    ``scale > 0`` skaliert das Argument von ``numpy.cos`` (gleiche Einheit wie ``y``).
    """
    m_arr = np.asarray(m, dtype=float).ravel()[: int(K)]
    kk = m_arr.size
    if kk < 1:
        raise ValueError("build_H: nach Abschneiden keine Mittelpunkte übrig")
    scl = float(scale)
    if not np.isfinite(scl) or scl <= 0.0:
        raise ValueError("build_H: scale muss endlich und positiv sein")
    H = np.zeros((kk, kk))
    for i in range(kk):
        for j in range(kk):
            y = abs(m_arr[i] - m_arr[j]) / 2.0
            H[i, j] = float(np.cos(y / scl))
    return H


def build_W(m: np.ndarray, K: int, k_neighbors: int = 10) -> np.ndarray:
    """
    Pro Zeile ``i``: die ``k_neighbors`` kleinsten Werte unter
    ``y = |m_i - m_j| // 2`` (ganzzahliger halber Mittelpunktsabstand),
    bezogen auf alle ``j != i``. Setze dort ``W_ij = 1/(1+y)``, ``W_ji = W_ij``.

    Wenn ``k_neighbors > K - 1`` ist, gilt automatisch höchstens ``K - 1`` Nachbarn.
    """
    ms = np.asarray(m, dtype=float).ravel()[: int(K)]
    kk = ms.size
    if kk < 2:
        raise ValueError("build_W: mindestens zwei Mittelpunkte nötig")
    mi = ms.astype(np.int64, copy=False)
    kp = max(0, int(k_neighbors))

    Wmat = np.zeros((kk, kk), dtype=float)
    limit = max(0, min(kp, kk - 1))

    for i in range(kk):
        distances = [
            (j, abs(int(mi[i]) - int(mi[j])) // 2) for j in range(kk) if j != i
        ]
        distances.sort(key=lambda t: t[1])
        for j, y in distances[:limit]:
            w = 1.0 / (1.0 + float(y))
            Wmat[i, j] = w
            Wmat[j, i] = w

    return Wmat


def graph_laplacian(W: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Diagonalmatrix ``D`` mit Zeilensummen ``D_ii = sum_j W_ij`` und Graph-Laplacian
    ``L = D - W`` (bei symmetrischen ``W`` auch Spaltensummen).
    """
    W_arr = np.asarray(W, dtype=float)
    sums = np.sum(W_arr, axis=1)
    Ddiag = np.diag(sums)
    lap = Ddiag - W_arr
    return Ddiag, lap


def unfold_spectrum(E: np.ndarray, *, spline_s: float | None = None) -> np.ndarray:
    """
    Glatte empirische Zählfunktion Ñ(E): Spline durch ``(E_i, i)`` bei sortiertem Spektrum.
    """
    e = np.sort(np.asarray(E, dtype=float).ravel())
    n_pts = e.size
    if n_pts < 3:
        raise ValueError("unfold_spectrum: mindestens drei Eigenwerte nötig")
    for ii in range(1, n_pts):
        if e[ii] <= e[ii - 1]:
            e[ii] = np.nextafter(e[ii - 1], np.inf)
    idx = np.arange(n_pts, dtype=float)
    if spline_s is None:
        spline_s = n_pts * 0.1
    spl = UnivariateSpline(e, idx, s=float(spline_s))
    return np.asarray(spl(e), dtype=float)


def level_spacings_unfolded(E: np.ndarray, *, spline_s: float | None = None) -> np.ndarray:
    """Abstände ΔÑ der entfalteten Zählfunktion, auf ⟨s⟩ = 1 normiert."""
    xi = unfold_spectrum(E, spline_s=spline_s)
    gaps = np.diff(xi)
    mu = float(np.mean(gaps))
    if not np.isfinite(mu) or mu <= 0:
        raise ValueError(f"level_spacings_unfolded: ungültiges ⟨ΔÑ⟩: {mu!r}")
    return gaps / mu


def normalize_spacings(sorted_eigenvalues: np.ndarray) -> np.ndarray:
    """
    Rohe NN-Abstände in der Energie, normiert auf ⟨s⟩ = 1 (kein Unfold über Ñ(E)).
    """
    e = np.asarray(sorted_eigenvalues, dtype=float).ravel()
    if e.size < 2:
        raise ValueError("normalize_spacings: mindestens zwei Eigenwerte nötig")
    s = np.diff(e)
    mu = float(np.mean(s))
    if not np.isfinite(mu) or mu <= 0:
        raise ValueError(f"normalize_spacings: ungültiger mittlerer Abstand {mu!r}")
    return s / mu


# --------------------------------------------------
# 1. Primzahlen bis N
# --------------------------------------------------
N = 10_000_000
primes = list(primerange(2, N))

prime_set = set(primes)

# --------------------------------------------------
# 2. Vierlinge finden
# --------------------------------------------------
quadruplets = []

for p in primes:
    if p + 8 >= N:
        break
    if (p+2 in prime_set and 
        p+6 in prime_set and 
        p+8 in prime_set):
        quadruplets.append(p)

print(f"Gefundene Vierlinge: {len(quadruplets)}")

# --------------------------------------------------
# 3. Mittelpunkte
# --------------------------------------------------
m = np.array([p + 4 for p in quadruplets])

# --------------------------------------------------
# 4. y-Werte (nur Nachbarn oder alle Paare)
# --------------------------------------------------

# Variante A: nur aufeinanderfolgende Vierlinge
y_values = []
for i in range(len(m)-1):
    y = abs(m[i+1] - m[i]) // 2
    y_values.append(y)

y_values = np.array(y_values)

print("Erste y-Werte:", y_values[:10])

for val in y_values[:20]:
    print(val, factorint(val // 3))

# --------------------------------------------------
# 5. Histogramm
# --------------------------------------------------
plt.figure()
plt.hist(np.log(y_values.astype(float)), bins=50)
plt.title("Histogramm der y-Werte")
plt.xlabel(r"$\ln(y)$")
plt.ylabel("Häufigkeit")
plt.show()

# --------------------------------------------------
# 6. Operator W — k nächste Nachbarn unter y (pro Zeile symmetrisiert)
# --------------------------------------------------

K = min(200, len(m))
k = 10  # oder z. B. 20
W = build_W(m, K=K, k_neighbors=k)
Dmat, Lap = graph_laplacian(W)

# --------------------------------------------------
# 7. Eigenwerte von L = D − W  (Diagonalmatrix D_ii = sum_j W_ij)
# --------------------------------------------------
eigvals = np.linalg.eigvalsh(Lap)
eigvals = np.sort(eigvals)

uniq, counts = np.unique(np.round(eigvals, 10), return_counts=True)
print("Gerundete Eigenwerte — Anzahl verschiedene Werte:", len(uniq))
print("uniq (Auszug bis 40):", uniq[:40])
print("counts (Auszug bis 40):", counts[:40])
print("Eigenwerte mit Vielfachheit > 1 (Anzahl):", int(np.sum(counts > 1)))

# --------------------------------------------------
# 8. Eigenwert-Histogramm
# --------------------------------------------------
plt.figure()
plt.hist(eigvals, bins=50)
plt.title("Eigenwerte des Graph-Laplace L = D - W (k-NN Gewichte)")
plt.xlabel("Eigenwert")
plt.ylabel("Häufigkeit")
plt.show()

# --------------------------------------------------
# 9. Spacing-Analyse: entfaltet (Spline) und Roh-Energie
# --------------------------------------------------

SPLINE_S: float | None = None  # None → s = 0.1 * N im UnivariateSpline

s_unfold = level_spacings_unfolded(eigvals, spline_s=SPLINE_S)
plt.figure()
plt.hist(s_unfold, bins=50, density=True)
plt.title(r"$L=D{-}$W$: normierte NN-Abstände ($\Delta\tilde N$, Spline-Unfold)")
plt.xlabel(r"$s / \langle s\rangle$")
plt.ylabel("Dichte")

s_raw = normalize_spacings(eigvals)
plt.figure()
plt.hist(s_raw, bins=50, density=True)
plt.title(r"$L=D{-}$W$: normierte NN-Abstände in $E$ (ohne Unfold)")
plt.xlabel(r"$s / \langle s\rangle$")
plt.ylabel("Dichte")
plt.show()

print(
    "Unfold ⟨s⟩:",
    float(np.mean(s_unfold)),
    "  Var:",
    float(np.var(s_unfold)),
)
print("Roh-E   ⟨s⟩:", float(np.mean(s_raw)), "  Var:", float(np.var(s_raw)))