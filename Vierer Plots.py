import numpy as np
import json
import matplotlib.pyplot as plt
import os

def safe_show():
    try:
        plt.show()
    except KeyboardInterrupt:
        print("\nPlot durch Benutzer abgebrochen.")

# --------------------------------------------------
# Konfiguration & Hilfsfunktionen
# --------------------------------------------------
N_LIMIT = 2_000_000

def family(p):
    r = p % 12
    if r == 1: return 0  # E
    if r == 5: return 1  # A
    if r == 7: return 2  # B
    if r == 11: return 3 # C
    return -1

def generate_quadruples(limit):
    print(f"Generiere Primzahlvierlinge bis {limit}...")
    sieve = np.ones(limit + 10, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i*i:limit+1:i] = False
    
    primes = np.nonzero(sieve)[0]
    quads = []
    for p in primes:
        if p + 8 > limit:
            break
        if sieve[p+2] and sieve[p+6] and sieve[p+8]:
            quads.append([
                p, p+2, p+6, p+8,
                family(p), family(p+2), family(p+6), family(p+8)
            ])
    return np.array(quads, dtype=np.int64), primes

# --------------------------------------------------
# Daten laden oder generieren
# --------------------------------------------------

# 1. Vierlinge (Q)
if os.path.exists("prime_quadruples_EABC_2M.npy"):
    print("Lade prime_quadruples_EABC_2M.npy...")
    Q = np.load("prime_quadruples_EABC_2M.npy")
else:
    Q, _ = generate_quadruples(N_LIMIT)
    # Optional: Speichern für nächstes Mal
    # np.save("prime_quadruples_EABC_2M.npy", Q)

# 2. Andere Dateien (optional)
def load_if_exists(fname):
    if os.path.exists(fname):
        print(f"Lade {fname}...")
        return np.load(fname)
    return None

P = load_if_exists("eabc_primes_2M.npy")
C = load_if_exists("quadruple_EABC_coords_2M.npy")
CF = load_if_exists("quadruple_EABC_coords_flat_2M.npy")
T = load_if_exists("quadruple_torus_embedding_2M.npy")
Z = load_if_exists("riemann_zeros_approx_10000.npy")

RJ = None
if os.path.exists("quadruple_riemann_compare_2M.json"):
    with open("quadruple_riemann_compare_2M.json", "r", encoding="utf-8") as f:
        RJ = json.load(f)

# --------------------------------------------------
# Familienkodierung
# --------------------------------------------------
code_to_family = {0: "E", 1: "A", 2: "B", 3: "C"}
family_to_residue = {"E": 1, "A": 5, "B": 7, "C": 11}

def fams_of_row(row):
    return tuple(code_to_family[int(x)] for x in row[4:])

# --------------------------------------------------
# Grundausgabe
# --------------------------------------------------
if P is not None:
    print("Anzahl EABC-Primzahlen:", len(P))
print("Anzahl Primzahlvierlinge:", len(Q))
print()

patterns = {}
for row in Q:
    pat = fams_of_row(row)
    patterns[pat] = patterns.get(pat, 0) + 1

print("Familienmuster:")
for k, v in sorted(patterns.items(), key=lambda kv: (-kv[1], kv[0])):
    print(f"{k}: {v}")

print()
print("Erste 10 Vierlinge:")
for i, row in enumerate(Q[:10], start=1):
    p1, p2, p3, p4 = row[:4]
    fam = fams_of_row(row)
    print(f"{i:2d}: ({p1}, {p2}, {p3}, {p4})  ->  {fam}")

if RJ is not None:
    print()
    print("Riemann-Vergleich:")
    print(json.dumps(RJ, indent=2, ensure_ascii=False))

# --------------------------------------------------
# Plot 1: Startwerte p der Vierlinge
# --------------------------------------------------
p_start = Q[:, 0]

plt.figure(figsize=(10, 4))
plt.plot(np.arange(1, len(p_start) + 1), p_start, marker="o", linestyle="", markersize=2)
plt.xlabel("Index des Vierlings")
    plt.ylabel("Startprime p")
    plt.title("Startwerte der Primzahlvierlinge bis 2 Mio")
    plt.tight_layout()
    safe_show()

# --------------------------------------------------
# Plot 2: Familienmuster als Balkendiagramm
# --------------------------------------------------
labels = ["-".join(k) for k in patterns.keys()]
values = list(patterns.values())

plt.figure(figsize=(8, 4))
plt.bar(labels, values)
plt.xlabel("Familienmuster")
plt.ylabel("Anzahl")
plt.title("Häufigkeit der EABC-Muster")
plt.xticks(rotation=45)
plt.tight_layout()
safe_show()

# --------------------------------------------------
# Plot 3: 3D-Toruspunkte des ersten Vierlings (falls vorhanden)
# --------------------------------------------------
if T is not None:
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")

    first = T[0]  # 4 Punkte des ersten Vierlings
    ax.scatter(first[:, 0], first[:, 1], first[:, 2], s=60)

    for i, point in enumerate(first):
        ax.text(point[0], point[1], point[2], f"p{i+1}")

    ax.set_title("Torus-Einbettung des ersten Vierlings")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.tight_layout()
    safe_show()
else:
    print("\n[Info] Torus-Einbettungsdaten (T) nicht gefunden. Plot 3 wird übersprungen.")
