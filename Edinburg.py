"""
Edinburg: Atiyah-Polynom-Matrix aus Primzahlen (S² → C, beta_i(z) = prod(z - u_ij)).

Was man mit den hier gezeigten Matrizen und Determinanten berechnen kann:
  - Determinante (bzw. Betrag): Test der Atiyah-Vermutung (Polynome linear unabhängig? det ≠ 0).
  - Rang / Singulärwerte (SVD): lineare Unabhängigkeit, Kondition, Stabilität der Konstellation.
  - Konditionszahl (σ_max/σ_min): wie „entartet“ die Richtungen u_ij sind; groß = nahe an Abhängigkeit.
  - Invarianten: Spur, Frobenius-Norm der Matrix; evtl. Eigenwerte von M*M für spektrale Deutung.

Bezug zum Bamberger Modell:
  - Bamberg-Kugel = S² mit Walter/Rest-Kanälen und Morley-Phase; hier ist S² die Riemann-Kugel
    (stereographische Projektion). Die Richtungen u_ij liegen auf S² und können mit Morley-Großkreisen
    (z.B. 2π/3-Phasen) verglichen werden.
  - Quadrupel: Wählt man Primzahlen als ein Bamberg-Quadrupel (je eine aus E, A, B, C mod 12),
    so ist die Atiyah-Matrix die Matrix zu dieser Quadrupel-Konstellation; det ≠ 0 bedeutet
    lineare Unabhängigkeit der zugehörigen Polynome für dieses Quadrupel.
  - Schütte-Spannung: Man kann für verschiedene Primzahlmengen (z.B. mit hoher vs. niedriger
    Spannung T auf E) die Determinante oder die Konditionszahl vergleichen; mögliche Korrelation
    von |det| oder σ_min mit T.
  - E8 / Alice–Taurus: Die Zeilen der Matrix entsprechen „Primzahl-Positionen“; die Struktur
    (E,A,B,C) und die Skalen (Kepler) lassen sich in einer Bamberg-konsistenten Einbettung
    prim_zu_R3(p) unterbringen (z.B. mod-12-Klasse → Richtung auf der Kugel).
"""
import numpy as np
from scipy.linalg import det, svd

# ---------------------------------------------------------
# 1. Konfiguration & Hilfsfunktionen
# ---------------------------------------------------------

def stereographic_projection(vec):
    """
    Projiziert einen Einheitsvektor (x, y, z) auf die komplexe Ebene (u).
    Atiyah nutzt die Nordpol-Projektion. N = (0,0,1).
    Formel: u = (x + iy) / (1 - z)
    """
    x, y, z = vec
    # Vermeidung von Division durch Null am Nordpol
    if np.isclose(z, 1.0):
        return complex(np.inf) 
    return complex(x, y) / (1.0 - z)

def poly_from_roots(roots):
    """
    Erstellt die Koeffizienten eines Polynoms aus seinen Wurzeln.
    p(z) = (z - r1)(z - r2)...
    Rückgabe: Liste der Koeffizienten [an, ..., a1, a0]
    """
    # Startet mit Polynom "1" (Grad 0)
    coeffs = [1.0]
    for r in roots:
        # Multipliziere aktuelles Polynom mit (z - r)
        # (a_n z^n + ... + a_0) * (z - r)
        # = a_n z^{n+1} + (a_{n-1} - r*a_n) z^n + ... - r*a_0
        new_coeffs = [0.0] * (len(coeffs) + 1)
        for i in range(len(coeffs)):
            new_coeffs[i] += coeffs[i]       # Term mal z
            new_coeffs[i+1] -= r * coeffs[i] # Term mal -r
        coeffs = new_coeffs
    return coeffs

# ---------------------------------------------------------
# 2. Das Bamberg-Modell: Mapping der Primzahlen
# ---------------------------------------------------------

def prim_zu_R3(p):
    """
    PLATZHALTER: Hier muss Ihre Quaternionen-Logik hin.
    Momentan: Eine 3D-Spirale ("Prime Helix"), um Kollisionen zu vermeiden.
    x = p * cos(p)
    y = p * sin(p)
    z = p
    (Dies simuliert den Imaginärteil eines Quaternions q = 0 + xi + yj + zk)
    """
    return np.array([
        p * np.cos(p), 
        p * np.sin(p), 
        p
    ], dtype=float)

# ---------------------------------------------------------
# 3. Atiyahs Algorithmus (Lecture 1)
# ---------------------------------------------------------

def atiyah_polynomials(primes, silent=False):
    n = len(primes)
    points = [prim_zu_R3(p) for p in primes]
    
    # Speichert die Polynom-Koeffizienten für jeden Punkt i
    # Jedes Polynom hat Grad n-1 -> n Koeffizienten
    matrix_rows = []
    
    if not silent:
        print(f"Analysiere Konstellation für {n} Primzahlen...")
    
    for i in range(n):
        roots = []
        for j in range(n):
            if i == j:
                continue
            
            # Vektor von i nach j
            diff = points[j] - points[i]
            
            # Gleichung (1.1) in Atiyahs Skript:
            # u_ij = (x_j - x_i) / ||x_j - x_i||
            dist = np.linalg.norm(diff)
            if dist == 0:
                continue # Sollte bei distinkten Primzahlen nicht passieren
                
            u_vec = diff / dist
            
            # Projektion auf Riemannsche Zahlenkugel (S^2 -> C)
            u_complex = stereographic_projection(u_vec)
            roots.append(u_complex)
        
        # Gleichung (1.2): beta_i(z) = prod(z - u_ij) 
        coeffs = poly_from_roots(roots)
        
        # Wir fügen die Koeffizienten als Zeile zur Matrix hinzu
        # coeffs ist [c_n-1, ..., c_0]
        matrix_rows.append(coeffs)

    return np.array(matrix_rows)


def primes_up_to(N):
    """Eratosthenes: Liste aller Primzahlen <= N."""
    if N < 2:
        return []
    is_prime = np.ones(N + 1, dtype=bool)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(np.sqrt(N)) + 1):
        if is_prime[i]:
            is_prime[i * i : N + 1 : i] = False
    return [p for p in range(2, N + 1) if is_prime[p]]


def primes_by_class(r, count, primes_full=None):
    """
    Erste `count` Primzahlen p > 3 mit p ≡ r (mod 12).
    r in {1, 5, 7, 11} (E, A, B, C).
    """
    if primes_full is None:
        # Überschätzung: k-te Primzahl in einer Klasse ~ 4*k-te Primzahl
        primes_full = primes_up_to(max(5000, 50 * 20))
    out = [p for p in primes_full if p > 3 and p % 12 == r]
    return out[:count]


def first_n_quadrupel(n=50):
    """
    Die ersten n Bamberg-Quadrupel: je (p_E, p_A, p_B, p_C) mit
    p_E ≡ 1, p_A ≡ 5, p_B ≡ 7, p_C ≡ 11 (mod 12).
    """
    plist = primes_up_to(max(5000, n * 25))
    E = primes_by_class(1, n, plist)
    A = primes_by_class(5, n, plist)
    B = primes_by_class(7, n, plist)
    C = primes_by_class(11, n, plist)
    if len(E) < n or len(A) < n or len(B) < n or len(C) < n:
        raise ValueError("Nicht genug Primzahlen pro Klasse; n oder Sieb-Grenze erhöhen.")
    return list(zip(E, A, B, C))


def atiyah_det_and_cond(primes, silent=True):
    """Atiyah-Matrix zu `primes`; liefert (|det|, Konditionszahl)."""
    M = atiyah_polynomials(primes, silent=silent)
    d = det(M)
    s = np.linalg.svd(M, compute_uv=False)
    cond = np.max(s) / (np.min(s) + 1e-20)
    return float(np.abs(d)), float(cond)


# ---------------------------------------------------------
# 4. Ausführung
# ---------------------------------------------------------
# Die definierten Matrizen (Atiyah-Polynom-Koeffizienten) passen zu unserer
# Faktorisierung: dieselben Restklassen (mod 12) und die geometrische
# Struktur (E, A, B, C; Quadrupel, E8-Injection) tauchen in beiden auf.

# Die ersten n Primzahlen (für Bamberg-Bezug: wähle z.B. ein Quadrupel pro Klasse E,A,B,C mod 12)
primes = [2, 3, 5, 7, 11, 13, 17]  # 2,3; 5,13; 7; 11,17 → E=13,17 A=5,13 B=7 C=11
# (Erhöhen Sie n, um die Stabilität zu testen. Atiyah sagt, n > 5 ist unbewiesen!)

# Berechne die Matrix der Polynom-Koeffizienten
coeff_matrix = atiyah_polynomials(primes)

# Berechne Determinante (Normalisierte Determinante D_R im Limit)
# Da wir komplexe Zahlen haben, schauen wir auf den Betrag.
determinant = det(coeff_matrix)
norm_det = np.abs(determinant)

# Singulärwertzerlegung für numerische Stabilität (Rang-Check) und Kondition
_, singular_values, _ = svd(coeff_matrix)
rank = np.sum(singular_values > 1e-10)
sigma_max = np.max(singular_values)
sigma_min = np.min(singular_values)
cond = sigma_max / (sigma_min + 1e-20)  # Konditionszahl

print("-" * 30)
print(f"Atiyah-Matrix ({len(primes)}x{len(primes)}) konstruiert.")
print(f"Betrag der Determinante: {norm_det:.8e}")
print(f"Numerischer Rang: {rank} (Erwartet: {len(primes)})")
print(f"Konditionszahl (σ_max/σ_min): {cond:.4e}")

if rank == len(primes):
    print("ERGEBNIS: Die Polynome sind linear UNABHÄNGIG.")
    print("-> Die Atiyah-Vermutung hält für diese Primzahl-Konfiguration.")
else:
    print("ERGEBNIS: Lineare Abhängigkeit gefunden (oder numerischer Fehler).")

# Optional: Zeige die Koeffizienten des ersten Polynoms (für p=2)
print("-" * 30)
print(f"Beispiel-Polynom beta_0 (für p={primes[0]}):")
print(np.round(coeff_matrix[0], 2))

# ---------------------------------------------------------
# 5. Erste 50 E-A-B-C-Quadrupel (Bamberger Modell)
# ---------------------------------------------------------
print("\n" + "=" * 60)
print("Erste 50 Bamberg-Quadrupel (E, A, B, C mod 12): Atiyah |det| und Kondition")
print("=" * 60)

quadrupel = first_n_quadrupel(50)
results = []
for (pE, pA, pB, pC) in quadrupel:
    qprimes = [pE, pA, pB, pC]
    abs_det, cond = atiyah_det_and_cond(qprimes, silent=True)
    results.append((pE, pA, pB, pC, abs_det, cond))

# Tabelle: erste 10 + Übersicht
print(f"{'#':>3}  {'p_E (E)':>8}  {'p_A (A)':>8}  {'p_B (B)':>8}  {'p_C (C)':>8}  {'|det|':>12}  {'cond':>10}")
print("-" * 70)
for i, (pE, pA, pB, pC, abs_det, cond) in enumerate(results[:10], 1):
    print(f"{i:>3}  {pE:>8}  {pA:>8}  {pB:>8}  {pC:>8}  {abs_det:>12.4e}  {cond:>10.2e}")
print("  ...")
for i, (pE, pA, pB, pC, abs_det, cond) in enumerate(results[45:50], 46):
    print(f"{i:>3}  {pE:>8}  {pA:>8}  {pB:>8}  {pC:>8}  {abs_det:>12.4e}  {cond:>10.2e}")

abs_dets = [r[4] for r in results]
conds = [r[5] for r in results]
nonzero = sum(1 for d in abs_dets if d > 1e-10)
print("-" * 70)
print(f"Über alle 50 Quadrupel:")
print(f"  |det|:  Min = {min(abs_dets):.4e},  Max = {max(abs_dets):.4e},  Mittel = {np.mean(abs_dets):.4e}")
print(f"  Kondition:  Min = {min(conds):.2e},  Max = {max(conds):.2e}")
print(f"  Atiyah (det ≠ 0): {nonzero}/50 Quadrupel mit linear unabhängigen Polynomen.")