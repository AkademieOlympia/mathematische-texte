import numpy as np
import matplotlib.pyplot as plt
from math import isqrt
import sys

N = 2_000_000

print(f"Siebe Primzahlen bis {N}...")
def prime_sieve(n: int):
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    limit = isqrt(n)
    for k in range(2, limit + 1):
        if sieve[k]:
            sieve[k * k : n + 1 : k] = False
    return sieve

sieve = prime_sieve(N)

def is_prime_vec(p):
    mask = (p >= 0) & (p <= N)
    result = np.zeros_like(p, dtype=bool)
    result[mask] = sieve[p[mask]]
    return result

primes = np.nonzero(sieve)[0]
A_primes = primes[primes % 12 == 5]

def delta_abce_vec(a, b, c, e):
    return e - (a * b * c) ** (1.0 / 3.0) - 16.0 / 3.0

# Gitter
db_vals = np.arange(2, 22, 2)
dc_vals = np.arange(4, 24, 2)
de_vals = np.arange(6, 26, 2)

offset_grid = []
for db in db_vals:
    for dc in dc_vals:
        for de in de_vals:
            if db < dc < de:
                offset_grid.append([db, dc, de])
offset_grid = np.array(offset_grid) # (M, 3)

A_sample = A_primes[::200]
print(f"Sample-Größe A: {len(A_sample)}")
print(f"Gitter-Punkte: {len(offset_grid)}")

# Broadcasting
M = len(offset_grid)
K = len(A_sample)

DB = offset_grid[:, 0:1]
DC = offset_grid[:, 1:2]
DE = offset_grid[:, 2:3]

A_mat = A_sample.reshape(1, K)
B_mat = A_mat + DB
C_mat = A_mat + DC
E_mat = A_mat + DE

# Filterung
valid_mask = (E_mat <= N) & \
             (B_mat % 12 == 7) & \
             (C_mat % 12 == 11) & \
             (E_mat % 12 == 1)

# Deltas
deltas = delta_abce_vec(A_mat, B_mat, C_mat, E_mat)
abs_deltas = np.abs(deltas)

# Primzahltreffer
hits_mask = valid_mask & \
            is_prime_vec(B_mat) & \
            is_prime_vec(C_mat) & \
            is_prime_vec(E_mat)

# Aggregation pro Gitterpunkt
# Wir mitteln nur über gültige Einträge
# np.nanmean ignoriert NaNs
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    VAL = np.nanmean(abs_deltas_masked, axis=1)
REAL_HITS = np.sum(hits_mask, axis=1)

# Ergebnisse aufbereiten
X = offset_grid[:, 0]
Y = offset_grid[:, 1]
Z = offset_grid[:, 2]

# Plot
try:
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    sizes = 30 + 60 * (REAL_HITS > 0) + 10 * REAL_HITS
    # Filtern: Nur Punkte anzeigen, die gültige Deltas haben (nicht NaN)
    valid_idx = ~np.isnan(VAL)
    
    sc = ax.scatter(X[valid_idx], Y[valid_idx], Z[valid_idx], c=VAL[valid_idx], s=sizes[valid_idx])

    ax.set_xlabel("db")
    ax.set_ylabel("dc")
    ax.set_zlabel("de")
    ax.set_title("3D-Resonanzfläche A-B-C-E: Mittelwert von |Delta|")
    fig.colorbar(sc, ax=ax, label="mittlerer |zentrierter Resonanzdefekt|")
    plt.tight_layout()
    print("Zeige Plot... (Fenster schließen zum Beenden)")
    plt.show()

except KeyboardInterrupt:
    print("\nPlot durch Benutzer abgebrochen.")

# Beste Ergebnisse
if np.any(~np.isnan(VAL)):
    order = np.argsort(VAL)
    # NaNs stehen bei argsort am Ende
    
    print("\nBeste Offset-Tripel nach mittlerem |Delta|:")
    count = 0
    for k in order:
        if np.isnan(VAL[k]):
            continue
        print(
            f"(db,dc,de)=({X[k]},{Y[k]},{Z[k]})   "
            f"mean|Delta|={VAL[k]:.6f}   "
            f"echte Treffer im Sample={REAL_HITS[k]}"
        )
        count += 1
        if count >= 20:
            break
else:
    print("\nKeine gültigen Datenpunkte gefunden.")