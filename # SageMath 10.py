# SageMath 10.5 Code für das erweiterte Bamberg-Modell
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def get_patterned_intersection(n):
    """
    Berechnet die Schnittmenge aus Teilern und Ziffern (ohne 0).
    Entspricht der Definition aus Campbell (2026).
    """
    # Ziffern holen
    n_digits = set([int(d) for d in str(n)])
    # Teiler holen (nur die bis 9 sind laut Paper relevant für 'Patterned',
    # aber für das Modell nutzen wir alle einstelligen Teiler)
    n_divisors = set([d for d in divisors(n) if d <= 9])
    
    # Schnittmenge
    intersection = n_digits.intersection(n_divisors)
    
    # Wir entfernen die 0, falls sie als Teiler auftauchen könnte (unmöglich),
    # aber Ziffer 0 ist im Paper explizit ausgeschlossen für die Pattern-Regel.
    if 0 in intersection:
        intersection.remove(0)
        
    return len(intersection)

# Initialisierung des Walks
# Position im 4D Raum: [Real, i, j, k]
pos = np.array([0, 0, 0, 0])
path = [pos.copy()]
prime_indices = []

max_n = 5000

print(f"Berechne Trajektorie für n = 1 bis {max_n}...")

for n in range(1, max_n + 1):
    # 1. Bestimme die "Pattern-Stärke"
    pattern_count = get_patterned_intersection(n)
    
    # 2. Mappe auf Quaternionen-Basis (Modulo 4 Zyklus)
    # Dies ist der "Turn Operator" erweitert auf 4D
    step = np.zeros(4)
    direction = pattern_count % 4
    
    if direction == 0:   # Richtung 1 (Real)
        step[0] = 1
    elif direction == 1: # Richtung i
        step[1] = 1
    elif direction == 2: # Richtung j
        step[2] = 1
    elif direction == 3: # Richtung k
        step[3] = 1
        
    # 3. Gehe den Schritt
    pos = pos + step
    path.append(pos.copy())
    
    # 4. Prüfe, ob Primzahl (für Visualisierung)
    if is_prime(n):
        prime_indices.append(len(path) - 1)

# Daten für Plot vorbereiten
path = np.array(path)
xs = path[:, 1] # i-Achse
ys = path[:, 2] # j-Achse
zs = path[:, 3] # k-Achse
# Der Realteil (path[:, 0]) könnte als Farbe oder Zeitachse genutzt werden

# Primzahlen extrahieren
primes_x = xs[prime_indices]
primes_y = ys[prime_indices]
primes_z = zs[prime_indices]

print("Erstelle 3D-Visualisierung...")

# 3D Plot
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Den Pfad plotten (Blau = Composite/Patterned Mix)
ax.plot(xs, ys, zs, alpha=0.5, linewidth=0.8, label='Quaternion Walk (Composite)', color='blue')

# Die Primzahlen hervorheben (Rot)
ax.scatter(primes_x, primes_y, primes_z, s=5, c='red', depthshade=False, label='Primzahlen (P-Nodes)')

# Start und Ende
ax.text(xs[0], ys[0], zs[0], "START", color='green', fontweight='bold')
ax.text(xs[-1], ys[-1], zs[-1], f"n={max_n}", color='black', fontweight='bold')

ax.set_xlabel('i - Achse (Pattern mod 4 = 1)')
ax.set_ylabel('j - Achse (Pattern mod 4 = 2)')
ax.set_zlabel('k - Achse (Pattern mod 4 = 3)')
ax.set_title(f'Bamberg-Modell: Patterned Number Walk (n={max_n})\nRot = Primzahlen im Strukturfluss')

plt.legend()
plt.show()

# Zusatzanalyse: Wo landen die Primzahlen meistens?
print("Verteilung der Schritte (Modulo 4):")
counts = {0:0, 1:0, 2:0, 3:0}
for n in range(1, max_n + 1):
    c = get_patterned_intersection(n) % 4
    counts[c] += 1
print(f"Richtung 1 (Real): {counts[0]}")
print(f"Richtung i:       {counts[1]}")
print(f"Richtung j:       {counts[2]}")
print(f"Richtung k:       {counts[3]}")