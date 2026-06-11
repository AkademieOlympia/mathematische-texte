import math
from collections import defaultdict
import numpy as np

def is_prime_basic(n):
    """Klassischer Primzahltest für die Normen."""
    if n < 2:
        return False
    for i in range(2, int(math.isqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def generate_hurwitz_primes_pure(norm_bound):
    """
    Generiert Hurwitz-Quaternionen in reinem Python bis zu einer maximalen Norm.
    Ein Hurwitz-Quaternion hat entweder nur ganzzahlige oder nur halbganzzahlige Koeffizienten.
    """
    hurwitz_primes = []
    # Da wir mit Halben arbeiten, verdoppeln wir die Schleifengrenzen zur Ganzzahlarithmetik
    bound = int(2 * math.sqrt(norm_bound))
    
    for a0 in range(-bound, bound + 1):
        for a1 in range(-bound, bound + 1):
            for a2 in range(-bound, bound + 1):
                for a3 in range(-bound, bound + 1):
                    # Hurwitz-Bedingung: Alle Koordinaten müssen die gleiche Parität haben
                    if (a0 % 2 == a1 % 2) and (a1 % 2 == a2 % 2) and (a2 % 2 == a3 % 2):
                        # Echte Norm berechnen: (a0/2)^2 + (a1/2)^2 + ... = (a0^2 + ...)/4
                        norm_quad = a0**2 + a1**2 + a2**2 + a3**2
                        if norm_quad % 4 == 0:
                            norm = norm_quad // 4
                        else:
                            norm = norm_quad / 4
                            
                        # Wir betrachten ganzzahlige Normen für die Primzahl-Gitterstruktur
                        if 0 < norm <= norm_bound and isinstance(norm, int) or (isinstance(norm, float) and norm.is_integer()):
                            norm_int = int(norm)
                            if is_prime_basic(norm_int):
                                # Speichern als (Scalar, i, j, k, Norm)
                                hurwitz_primes.append((a0/2, a1/2, a2/2, a3/2, norm_int))
                                
    return hurwitz_primes

def analyze_eabc_monopoles(primes):
    """
    Analysiert die topologischen Defekte (Monopol-Kandidaten).
    Sucht nach reinen Imaginär-Primzahlen (Skalarteil = 0), die eine
    lokale Symmetriebrechung im harmonischen Fluss der Eichleitung erzeugen.
    """
    defects = defaultdict(int)
    
    for p in primes:
        a0, a1, a2, a3, norm = p
        
        # Monopol-Bedingung im Modell: Verschwinden des zeitartigen/skalaren Anteils
        # Dies isoliert den reinen Raumfluss (transversal/longitudinal)
        if abs(a0) < 1e-9:
            # Diskrete Gitterkoordinate im 3D-Raum (i, j, k)
            coord = (round(a1, 1), round(a2, 1), round(a3, 1))
            
            # Die "Ladungsstärke" korreliert im Modell mit der Dichte der Primzahl-Norm
            defects[coord] += 1
            
    return defects

# --- Hauptprogramm ---
if __name__ == "__main__":
    # Definition der oberen Schranke für die numerische Untersuchung
    MAX_NORM = 150
    
    print(f"--- Starte arithmetische Monopol-Analyse im EABC-Modell ---")
    print(f"Berechne Hurwitz-Gitter bis Norm: {MAX_NORM}...\n")
    
    # 1. Gitter-Primes generieren
    primes = generate_hurwitz_primes_pure(MAX_NORM)
    print(f"Gefundene Hurwitz-Primelemente im Gitter: {len(primes)}")
    
    # 2. Topologische Defekte isolieren
    monopoles = analyze_eabc_monopoles(primes)
    print(f"Identifizierte topologische Defekt-Knotenpunkte (Monopole): {len(monopoles)}\n")
    
    # 3. Strukturierte Ausgabe der signifikantesten Defekte
    print(f"{'3D-Gitterkoordinate (i, j, k)':<35} | {'Akkumulierte Feld-Divergenz (Ladung)':<15}")
    print("-" * 65)
    
    # Sortiert nach der Stärke des Defekts
    sorted_monopoles = sorted(monopoles.items(), key=lambda x: x[1], reverse=True)
    
    for coord, charge in sorted_monopoles[:12]:
        print(f"{str(coord):<35} | {charge:<15}")
        
    # 4. Numerische Verknüpfung zur Eichleitung (Matrix-Repräsentation für Folgeberechnungen)
    # Wir überführen die Defekte in ein Numpy-Array für physikalische Tensorenberechnungen
    if sorted_monopoles:
        matrix_data = np.array([[c[0], c[1], c[2], ch] for c, ch in sorted_monopoles])
        print(f"\n[Info] Datenstruktur erfolgreich für Tensor-Analyse vorbereitet. Form: {matrix_data.shape}")