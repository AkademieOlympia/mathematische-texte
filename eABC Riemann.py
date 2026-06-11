import numpy as np
import matplotlib.pyplot as plt
from math import isqrt

# optional für schönere 3D-Projektion
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

def safe_show():
    try:
        plt.show()
    except KeyboardInterrupt:
        print("\nPlot durch Benutzer abgebrochen.")

# ============================================================
# EINSTELLUNGEN
# ============================================================
N = 10_000_000
ZERO_FILE = "zeros6.npy"

NUM_ZEROS = 300
TOP_MODES = 12

SAVE_PREFIX = "vierlinge_riemann_e8"

# ============================================================
# PRIMZAHL-SIEB
# ============================================================
def prime_sieve(n: int) -> np.ndarray:
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for k in range(2, isqrt(n) + 1):
        if sieve[k]:
            sieve[k * k:n + 1:k] = False
    return sieve

print(f"Siebe Primzahlen bis {N:,} ...")
sieve = prime_sieve(N)

# ============================================================
# RIEMANN-NULLSTELLEN LADEN
# ============================================================
zeros = np.load(ZERO_FILE).astype(float)
print("Geladene Nullstellen:", len(zeros))
gammas = zeros[:NUM_ZEROS]

# ============================================================
# VIERLINGE BERECHNEN
# Standardform: (p, p+2, p+6, p+8)
# ABCE: p % 12 == 5
# CEAB: p % 12 == 11
# ============================================================
abce_rows = []
ceab_rows = []
all_rows = []

print("Berechne Primzahlvierlinge ...")

for p in range(5, N - 8 + 1, 2):
    if p % 12 not in (5, 11):
        continue
    if not (sieve[p] and sieve[p + 2] and sieve[p + 6] and sieve[p + 8]):
        continue

    p1, p2, p3, p4 = p, p + 2, p + 6, p + 8

    if p % 12 == 5:
        # ABCE: (a,b,c,e) = (p, p+2, p+6, p+8)
        a, b, c, e = p1, p2, p3, p4
        x = a
        delta = e - (a * b * c) ** (1.0 / 3.0) - 16.0 / 3.0
        trend = (
            28.0 / (9.0 * x)
            - 832.0 / (81.0 * x * x)
            + 10288.0 / (243.0 * x ** 3)
        )
        delta_detr = delta - trend
        fam_codes = [1, 2, 3, 0]  # A B C E
        row = [p1, p2, p3, p4, x, delta, delta_detr, 0] + fam_codes
        abce_rows.append(row)
        all_rows.append(row)

    else:
        # CEAB: (c,e,a,b) = (p, p+2, p+6, p+8)
        c, e, a, b = p1, p2, p3, p4
        x = c
        delta = e - (a * b * c) ** (1.0 / 3.0) + 4.0
        trend = (
            52.0 / (9.0 * x)
            - 1624.0 / (81.0 * x * x)
            + 23008.0 / (243.0 * x ** 3)
        )
        delta_detr = delta - trend
        fam_codes = [3, 0, 1, 2]  # C E A B
        row = [p1, p2, p3, p4, x, delta, delta_detr, 1] + fam_codes
        ceab_rows.append(row)
        all_rows.append(row)

abce = np.array(abce_rows, dtype=float)
ceab = np.array(ceab_rows, dtype=float)
allQ = np.array(all_rows, dtype=float)

print("ABCE-Vierlinge:", len(abce))
print("CEAB-Vierlinge:", len(ceab))
print("Gesamt:", len(allQ))

# ============================================================
# HILFSFUNKTIONEN
# ============================================================
def analyze_riemann_modes(name, arr, gammas, top_modes=12):
    """
    arr columns:
    0..3: p1..p4
    4: x (Basis)
    5: delta
    6: delta_detr
    7: type flag
    8..11: family codes
    """
    x = arr[:, 4].astype(float)
    u = np.log(x)
    y = arr[:, 6].astype(float)
    y = y - np.mean(y)

    # Design-Matrizen
    C = np.column_stack([np.cos(g * u) for g in gammas])
    S = np.column_stack([np.sin(g * u) for g in gammas])

    # Spalten normieren
    Cn = C / np.maximum(np.linalg.norm(C, axis=0), 1e-15)
    Sn = S / np.maximum(np.linalg.norm(S, axis=0), 1e-15)

    proj_cos = Cn.T @ y
    proj_sin = Sn.T @ y
    amps = np.sqrt(proj_cos ** 2 + proj_sin ** 2)

    idx = np.argsort(amps)[-top_modes:][::-1]

    print(f"\nTop-{top_modes} Riemann-Moden für {name}:")
    for j in idx:
        print(
            f"gamma={gammas[j]:.12f}   "
            f"amp={amps[j]:.12e}   "
            f"cos={proj_cos[j]:.12e}   "
            f"sin={proj_sin[j]:.12e}"
        )

    # Rekonstruktion
    y_fit = np.zeros_like(y)
    for j in idx:
        y_fit += proj_cos[j] * Cn[:, j] + proj_sin[j] * Sn[:, j]

    rss = np.sum((y - y_fit) ** 2)
    tss = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - rss / tss if tss > 0 else np.nan
    print(f"R^2 ({name}) = {r2:.12f}")

    # Gram-Matrizen
    G_cos = Cn.T @ Cn
    G_sin = Sn.T @ Sn
    G_mix = Cn.T @ Sn

    # Sortierung für Plot
    order = np.argsort(u)

    result = {
        "x": x,
        "u": u,
        "y": y,
        "y_fit": y_fit,
        "proj_cos": proj_cos,
        "proj_sin": proj_sin,
        "amps": amps,
        "top_idx": idx,
        "r2": r2,
        "G_cos": G_cos,
        "G_sin": G_sin,
        "G_mix": G_mix,
        "order": order,
    }
    return result


def phi_from_code(code):
    # 0=E, 1=A, 2=B, 3=C
    return {
        0: 0.0,
        1: 0.5 * np.pi,
        2: np.pi,
        3: 1.5 * np.pi
    }[int(code)]


def build_8d_embedding(arr):
    p = arr[:, :4].astype(float)
    f = arr[:, 8:12].astype(int)

    phi = np.vectorize(phi_from_code)(f)

    # 8D-Einbettung
    V = np.column_stack([
        np.log(p[:, 0]),
        np.log(p[:, 1]),
        np.log(p[:, 2]),
        np.log(p[:, 3]),
        np.cos(phi[:, 0]),
        np.cos(phi[:, 1]),
        np.sin(phi[:, 2]),
        np.sin(phi[:, 3]),
    ])

    V = V - np.mean(V, axis=0, keepdims=True)
    return V


def project_8d(V):
    dirs = []
    labels = []
    dim = V.shape[1]

    for i in range(dim):
        for j in range(i + 1, dim):
            d1 = np.zeros(dim)
            d1[i] = 1.0
            d1[j] = -1.0
            d1 /= np.linalg.norm(d1)
            dirs.append(d1)
            labels.append(f"e{i+1}-e{j+1}")

            d2 = np.zeros(dim)
            d2[i] = 1.0
            d2[j] = 1.0
            d2 /= np.linalg.norm(d2)
            dirs.append(d2)
            labels.append(f"e{i+1}+e{j+1}")

    dirs = np.array(dirs)
    proj = V @ dirs.T
    proj_var = np.var(proj, axis=0)
    idx = np.argsort(proj_var)[::-1]

    print("\nTop-20 8D-Projektionsrichtungen:")
    for k in idx[:20]:
        print(f"{labels[k]:>10s}   var={proj_var[k]:.12f}")

    return {
        "proj": proj,
        "dirs": dirs,
        "labels": labels,
        "proj_var": proj_var,
        "idx": idx
    }

# ============================================================
# RIEMANN-ANALYSEN
# ============================================================
res_abce = analyze_riemann_modes("ABCE", abce, gammas, top_modes=TOP_MODES)
res_ceab = analyze_riemann_modes("CEAB", ceab, gammas, top_modes=TOP_MODES)
res_all  = analyze_riemann_modes("BOTH", allQ, gammas, top_modes=TOP_MODES)

# ============================================================
# 8D-PROJEKTION
# ============================================================
V_abce = build_8d_embedding(abce)
V_ceab = build_8d_embedding(ceab)
V_all  = build_8d_embedding(allQ)

proj_abce = project_8d(V_abce)
proj_ceab = project_8d(V_ceab)
proj_all  = project_8d(V_all)

# ============================================================
# SPEICHERN
# ============================================================
np.save(f"{SAVE_PREFIX}_abce.npy", abce)
np.save(f"{SAVE_PREFIX}_ceab.npy", ceab)
np.save(f"{SAVE_PREFIX}_all.npy", allQ)

np.save(f"{SAVE_PREFIX}_abce_amps.npy", res_abce["amps"])
np.save(f"{SAVE_PREFIX}_ceab_amps.npy", res_ceab["amps"])
np.save(f"{SAVE_PREFIX}_all_amps.npy", res_all["amps"])

np.save(f"{SAVE_PREFIX}_abce_proj8d.npy", proj_abce["proj"])
np.save(f"{SAVE_PREFIX}_ceab_proj8d.npy", proj_ceab["proj"])
np.save(f"{SAVE_PREFIX}_all_proj8d.npy", proj_all["proj"])

print("\nDateien gespeichert.")

# ============================================================
# PLOTS: RIEMANN
# ============================================================
def plot_riemann_results(name, res, gammas):
    order = res["order"]
    u = res["u"][order]
    y = res["y"][order]
    fit = res["y_fit"][order]
    amps = res["amps"]

    plt.figure(figsize=(10, 5))
    plt.plot(u, y, ".", ms=3, label="detrendete Daten")
    plt.plot(u, fit, "-", lw=2, label="Top-Moden-Fit")
    plt.xlabel("u = log(x)")
    plt.ylabel("detrendeter Defekt")
    plt.title(f"Riemann-Rekonstruktion: {name}")
    plt.legend()
    plt.tight_layout()
    safe_show()

    plt.figure(figsize=(10, 5))
    plt.plot(gammas, amps, "-", lw=1)
    plt.xlabel("gamma")
    plt.ylabel("Amplitude")
    plt.title(f"Riemann-Projektionen: {name}")
    plt.tight_layout()
    safe_show()

    plt.figure(figsize=(7, 6))
    plt.imshow(res["G_cos"], origin="lower", aspect="auto")
    plt.colorbar(label="cos-Überlappung")
    plt.title(f"Gram-Matrix cos: {name}")
    plt.tight_layout()
    safe_show()

    plt.figure(figsize=(7, 6))
    plt.imshow(res["G_sin"], origin="lower", aspect="auto")
    plt.colorbar(label="sin-Überlappung")
    plt.title(f"Gram-Matrix sin: {name}")
    plt.tight_layout()
    safe_show()

plot_riemann_results("ABCE", res_abce, gammas)
plot_riemann_results("CEAB", res_ceab, gammas)
plot_riemann_results("BOTH", res_all, gammas)

# ============================================================
# PLOTS: 8D / "E8-nahe" Projektionen
# ============================================================
def plot_8d_results(name, proj_res):
    idx = proj_res["idx"]
    proj = proj_res["proj"]
    labels = proj_res["labels"]
    proj_var = proj_res["proj_var"]

    i1, i2 = idx[:2]
    x = proj[:, i1]
    y = proj[:, i2]

    plt.figure(figsize=(7, 7))
    plt.plot(x, y, ".", ms=4)
    plt.xlabel(labels[i1])
    plt.ylabel(labels[i2])
    plt.title(f"8D-Projektion: {name}")
    plt.tight_layout()
    safe_show()

    plt.figure(figsize=(10, 5))
    plt.hist(proj[:, i1], bins=40)
    plt.xlabel(labels[i1])
    plt.ylabel("Häufigkeit")
    plt.title(f"Histogramm stärkste 8D-Projektion: {name}")
    plt.tight_layout()
    safe_show()

    plt.figure(figsize=(10, 5))
    plt.plot(np.sort(proj_var)[::-1], ".")
    plt.xlabel("Rang")
    plt.ylabel("Varianz")
    plt.title(f"Varianzspektrum der 8D-Projektionen: {name}")
    plt.tight_layout()
    safe_show()

plot_8d_results("ABCE", proj_abce)
plot_8d_results("CEAB", proj_ceab)
plot_8d_results("BOTH", proj_all)

print("\nFertig.")