import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import det

# ---------------------------------------------------------
# 1. Hilfsfunktionen (Geometrie & Atiyah-Logik)
# ---------------------------------------------------------

def is_prime(n):
    """Einfacher Primzahltest."""
    if n <= 1: return False
    if n <= 3: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def generate_primes(limit):
    """Erzeugt Liste aller Primzahlen bis limit."""
    return [p for p in range(2, limit) if is_prime(p)]

def stereographic_projection(vec):
    """
    Projiziert R3 -> C (Atiyahs Methode).
    Vermeidet Division durch Null am Nordpol.
    """
    x, y, z = vec
    if np.isclose(z, 1.0):
        return complex(0, 0) # Fallback für Nordpol (eigentlich unendlich)
    return complex(x, y) / (1.0 - z)

def poly_from_roots(roots):
    """Erstellt Koeffizienten aus Nullstellen (Atiyahs Polynome)."""
    coeffs = [1.0]
    for r in roots:
        new_coeffs = [0.0] * (len(coeffs) + 1)
        for i in range(len(coeffs)):
            new_coeffs[i] += coeffs[i]
            new_coeffs[i+1] -= r * coeffs[i]
        coeffs = new_coeffs
    return coeffs

def prim_zu_R3(p):
    """
    Mapping: Primzahl -> R3 (Imaginärteil des Quaternions).
    Hier: Helix-Mapping, um die 'Drehung' der Primzahlen zu simulieren.
    """
    # Skalierung, damit die Werte nicht explodieren (wichtig für numerische Stabilität)
    scale = 0.1 
    p_s = p * scale
    return np.array([
        np.cos(p_s), 
        np.sin(p_s), 
        np.tanh(p_s) # Begrenzte Höhe, projiziert auf Einheitssphäre
    ], dtype=float)

def calculate_atiyah_energy(primes):
    """
    Berechnet die Energie E = -ln(|det(M)|) für eine Menge von Primzahlen.
    M ist die Matrix der Koeffizienten der Atiyah-Polynome.
    """
    n = len(primes)
    points = [prim_zu_R3(p) for p in primes]
    
    # Normalisierung der Punkte auf die Einheitssphäre (wichtig für Atiyah!)
    points = [v / np.linalg.norm(v) for v in points]
    
    matrix_rows = []
    for i in range(n):
        roots = []
        for j in range(n):
            if i == j: continue
            diff = points[j] - points[i]
            dist = np.linalg.norm(diff)
            if dist == 0: continue
            u_vec = diff / dist
            roots.append(stereographic_projection(u_vec))
        
        coeffs = poly_from_roots(roots)
        # Padding, falls Grad variiert (sollte bei n distinkten Punkten nicht passieren)
        if len(coeffs) < n:
            coeffs = [0.0]*(n - len(coeffs)) + coeffs
        matrix_rows.append(coeffs)
    
    mat = np.array(matrix_rows)
    
    # Berechnung der Determinante
    # Wir nutzen log-determinante via SVD um Overflow/Underflow zu vermeiden
    # |det| = prod(singular_values) -> log(|det|) = sum(log(sv))
    _, sv, _ = np.linalg.svd(mat)
    
    # Filter sehr kleine Singulärwerte (numerisches Rauschen)
    sv = sv[sv > 1e-15]
    
    log_det = np.sum(np.log(sv))
    
    # Energie ist negativ Log-Det (Hohe Determinante = Tiefe Energie)
    energy = -log_det
    return energy

# ---------------------------------------------------------
# 2. Datenerzeugung: Bamberg vs. Zufall
# ---------------------------------------------------------

print("Generiere Primzahlen und Konstellationen...")
max_val = 2000
all_primes = generate_primes(max_val)

# Sortierung nach Modulo 12 (Bamberg-Modell)
mods = {1: [], 5: [], 7: [], 11: []}
for p in all_primes:
    m = p % 12
    if m in mods:
        mods[m].append(p)

# Erzeuge Bamberg-Quadrupel (E, A, B, C)
bamberg_tuples = []
min_len = min(len(v) for v in mods.values())
for i in range(min_len):
    # Tupel: (p_1, p_5, p_7, p_11)
    t = [mods[1][i], mods[5][i], mods[7][i], mods[11][i]]
    bamberg_tuples.append(t)

# Erzeuge Zufalls-Quadrupel zum Vergleich
# Wir nehmen einfach zufällige 4er-Gruppen aus der Gesamtliste
import random
random_tuples = []
for _ in range(len(bamberg_tuples)):
    # Wähle 4 zufällige Primzahlen aus dem gleichen Bereich
    rt = random.sample(all_primes, 4)
    random_tuples.append(rt)

# ---------------------------------------------------------
# 3. Energie-Berechnung
# ---------------------------------------------------------

print(f"Berechne Atiyah-Energie für {len(bamberg_tuples)} Quadrupel...")

bamberg_energies = []
random_energies = []
x_axis = [] # Indizes oder Summe der Primzahlen als x-Achse

for i in range(len(bamberg_tuples)):
    # Bamberg
    b_tup = bamberg_tuples[i]
    e_b = calculate_atiyah_energy(b_tup)
    bamberg_energies.append(e_b)
    
    # Random
    r_tup = random_tuples[i]
    e_r = calculate_atiyah_energy(r_tup)
    random_energies.append(e_r)
    
    x_axis.append(i+1)

# ---------------------------------------------------------
# 4. Visualisierung
# ---------------------------------------------------------

plt.figure(figsize=(10, 6))

# Plot Bamberg (Blau)
plt.plot(x_axis, bamberg_energies, label='Bamberg-Quadrupel (Mod 12)', 
         color='blue', alpha=0.7, linewidth=2)

# Plot Random (Orange/Rot)
plt.plot(x_axis, random_energies, label='Zufalls-Quadrupel', 
         color='red', alpha=0.4, linestyle='dashed')

# Trendlinien (Gleitender Durchschnitt)
def moving_average(a, n=10):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

if len(x_axis) > 10:
    ma_bamberg = moving_average(bamberg_energies)
    ma_random = moving_average(random_energies)
    plt.plot(x_axis[9:], ma_bamberg, color='darkblue', linewidth=3, label='Trend Bamberg')
    plt.plot(x_axis[9:], ma_random, color='darkred', linewidth=3, label='Trend Zufall')

plt.title('Atiyah-Energie: Bamberg-Modell vs. Zufall\n(Niedrigere Energie = Stabilere Geometrie)', fontsize=14)
plt.xlabel('Index des Quadrupels (wachsende Primzahlen)', fontsize=12)
plt.ylabel('Energie E = -ln(|Det|)', fontsize=12)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()

# ---------------------------------------------------------
# 5. Statistische Auswertung
# ---------------------------------------------------------
avg_b = np.mean(bamberg_energies)
avg_r = np.mean(random_energies)

print("\n=== Statistische Zusammenfassung ===")
print(f"Durchschnitts-Energie Bamberg: {avg_b:.4f}")
print(f"Durchschnitts-Energie Zufall:  {avg_r:.4f}")

if avg_b < avg_r:
    print("\nERGEBNIS: Das Bamberg-Modell ist 'energetisch günstiger'.")
    print("Die Primzahlen verteilen sich geordneter als der Zufall.")
else:
    print("\nERGEBNIS: Kein signifikanter energetischer Vorteil gegenüber Zufall gefunden.")