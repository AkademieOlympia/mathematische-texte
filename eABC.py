import numpy as np
import matplotlib.pyplot as plt
from math import isqrt
import sys

N = 2_000_000

# --------------------------------------------------
# Primzahlsieb
# --------------------------------------------------
def prime_sieve(n: int):
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    limit = isqrt(n)
    for k in range(2, limit + 1):
        if sieve[k]:
            sieve[k * k : n + 1 : k] = False
    return sieve

print(f"Siebe Primzahlen bis {N}...")
sieve = prime_sieve(N)

def is_prime_vec(p):
    # Vectorized is_prime check
    # p must be within bounds 0..N
    mask = (p >= 0) & (p <= N)
    result = np.zeros_like(p, dtype=bool)
    result[mask] = sieve[p[mask]]
    return result

# --------------------------------------------------
# Datenaufbereitung
# --------------------------------------------------
# Wir nutzen numpy arrays für schnellen Zugriff
primes_all = np.nonzero(sieve)[0]
# Filter für Familien
A_primes = primes_all[primes_all % 12 == 5]
C_primes = primes_all[primes_all % 12 == 11]

# Samples für Heatmap (Performance)
step = 80
A_sample = A_primes[::step]
C_sample = C_primes[::step]

print(f"Sample-Größen: A={len(A_sample)}, C={len(C_sample)}")

# --------------------------------------------------
# Resonanzdefekte (Vektorisiert)
# --------------------------------------------------
def delta_abce_vec(a, b, c, e):
    # e - (abc)^(1/3) - 16/3
    return e - (a * b * c) ** (1.0 / 3.0) - 16.0 / 3.0

def delta_ceab_vec(a, b, c, e):
    # e - (abc)^(1/3) + 4
    return e - (a * b * c) ** (1.0 / 3.0) + 4.0

# --------------------------------------------------
# Analyse A-B-C-E
# --------------------------------------------------
print("Berechne A-B-C-E Heatmap...")

# Gitterdefinition
db_vals = np.array([2, 14, 26, 38])
dc_vals = np.array([6, 18, 30, 42])
de_vals = np.array([8, 20, 32, 44])

grid_abce = []
for db in db_vals:
    for dc in dc_vals:
        for de in de_vals:
            if db < dc < de:
                grid_abce.append([db, dc, de])
grid_abce = np.array(grid_abce)  # Shape (M, 3)

# Broadcasting: (M, 1) + (1, K) -> (M, K)
# Zeilen: Gitterpunkte (Offsets), Spalten: A-Samples
M = len(grid_abce)
K = len(A_sample)

DB = grid_abce[:, 0:1]
DC = grid_abce[:, 1:2]
DE = grid_abce[:, 2:3]

A_mat = A_sample.reshape(1, K)
B_mat = A_mat + DB
C_mat = A_mat + DC
E_mat = A_mat + DE

# Filterung: Bereich und Residuen
valid_mask = (E_mat <= N) & \
             (B_mat % 12 == 7) & \
             (C_mat % 12 == 11) & \
             (E_mat % 12 == 1)

# Defekte berechnen
deltas = delta_abce_vec(A_mat, B_mat, C_mat, E_mat)
heat_abce = np.abs(deltas)
# Ungültige Einträge ausblenden (NaN)
heat_abce[~valid_mask] = np.nan

# Echte Treffer finden (alle sind prim)
hits_mask = valid_mask & \
            is_prime_vec(A_mat * np.ones((M, 1), dtype=int)) & \
            is_prime_vec(B_mat) & \
            is_prime_vec(C_mat) & \
            is_prime_vec(E_mat)

# Koordinaten der Treffer für Scatterplot
hit_y, hit_x = np.where(hits_mask)
# Werte sammeln
real_hits_abce = []
for y, x in zip(hit_y, hit_x):
    a = A_sample[x]
    db, dc, de = grid_abce[y]
    d = deltas[y, x]
    real_hits_abce.append((a, a+db, a+dc, a+de, d))


# --------------------------------------------------
# Analyse C-E-A-B
# --------------------------------------------------
print("Berechne C-E-A-B Heatmap...")

de_vals_c = np.array([2, 14, 26, 38])
da_vals_c = np.array([6, 18, 30, 42])
db_vals_c = np.array([8, 20, 32, 44])

grid_ceab = []
for de in de_vals_c:
    for da in da_vals_c:
        for db in db_vals_c:
            if de < da < db:
                grid_ceab.append([de, da, db])
grid_ceab = np.array(grid_ceab)

M_c = len(grid_ceab)
K_c = len(C_sample)

DE_c = grid_ceab[:, 0:1]
DA_c = grid_ceab[:, 1:2]
DB_c = grid_ceab[:, 2:3]

C_mat = C_sample.reshape(1, K_c)
E_mat = C_mat + DE_c
A_mat = C_mat + DA_c
B_mat = C_mat + DB_c

valid_mask_c = (B_mat <= N) & \
               (E_mat % 12 == 1) & \
               (A_mat % 12 == 5) & \
               (B_mat % 12 == 7)

deltas_c = delta_ceab_vec(A_mat, B_mat, C_mat, E_mat)
heat_ceab = np.abs(deltas_c)
heat_ceab[~valid_mask_c] = np.nan

hits_mask_c = valid_mask_c & \
              is_prime_vec(C_mat * np.ones((M_c, 1), dtype=int)) & \
              is_prime_vec(E_mat) & \
              is_prime_vec(A_mat) & \
              is_prime_vec(B_mat)

hit_y_c, hit_x_c = np.where(hits_mask_c)
real_hits_ceab = []
for y, x in zip(hit_y_c, hit_x_c):
    c = C_sample[x]
    de, da, db = grid_ceab[y]
    d = deltas_c[y, x]
    # Tupel Reihenfolge: (c, e, a, b)
    e = c + de
    a = c + da
    b = c + db
    real_hits_ceab.append((c, e, a, b, d))


# --------------------------------------------------
# Plots
# --------------------------------------------------
try:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 16))

    # Plot 1: A-B-C-E
    im1 = ax1.imshow(
        heat_abce,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap="viridis_r" # Reverse colormap -> dunkel = kleiner Defekt = gut
    )
    fig.colorbar(im1, ax=ax1, label="|zentrierter Resonanzdefekt|")
    
    if len(hit_x) > 0:
        ax1.scatter(hit_x, hit_y, c="red", marker="o", s=20, label="Primzahl-Vierlinge")
        ax1.legend()

    ax1.set_title("Heatmap A-B-C-E: |Delta| (dunkel = gut)")
    ax1.set_xlabel(f"Index des gesampelten A-Primzahl-Basiswerts (Step={step})")
    ax1.set_ylabel("Offset-Kombination Index")

    # Plot 2: C-E-A-B
    im2 = ax2.imshow(
        heat_ceab,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap="viridis_r"
    )
    fig.colorbar(im2, ax=ax2, label="|zentrierter Resonanzdefekt|")

    if len(hit_x_c) > 0:
        ax2.scatter(hit_x_c, hit_y_c, c="red", marker="o", s=20, label="Primzahl-Vierlinge")
        ax2.legend()

    ax2.set_title("Heatmap C-E-A-B: |Delta| (dunkel = gut)")
    ax2.set_xlabel(f"Index des gesampelten C-Primzahl-Basiswerts (Step={step})")
    ax2.set_ylabel("Offset-Kombination Index")

    plt.tight_layout()
    print("Zeige Plots... (Fenster schließen zum Beenden)")
    plt.show()

except KeyboardInterrupt:
    print("\nPlot durch Benutzer abgebrochen.")

# --------------------------------------------------
# Textausgabe
# --------------------------------------------------
print("\n" + "="*50)
print(f"Echte Treffer im A-B-C-E-Sample ({len(real_hits_abce)}):")
for r in real_hits_abce[:15]:
    print(f"  (a,b,c,e)={r[:4]}, delta={r[4]:.6f}")

print("\n" + "="*50)
print(f"Echte Treffer im C-E-A-B-Sample ({len(real_hits_ceab)}):")
for r in real_hits_ceab[:15]:
    print(f"  (c,e,a,b)={r[:4]}, delta={r[4]:.6f}")
