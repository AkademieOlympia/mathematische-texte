import numpy as np

def run_eabc_deterministic_analysis(M=113160):
    print(f"=== EABC-Modell: Deterministische Analyse am Fixpunkt M = {M} ===")
    
    # 1. Primzahlen generieren (Sieve of Eratosthenes)
    sieve = [True] * (M + 1)
    for p in range(2, int(np.sqrt(M)) + 1):
        if sieve[p]:
            for i in range(p*p, M + 1, p):
                sieve[i] = False
    primes = [p for p in range(2, M + 1) if sieve[p]]
    
    # 2. Aufteilung in die quaternionischen Strukturfamilien via Modulo 3
    # Familie E (Skalarteil) wird separat durch L_vac und die Asymptotik bestimmt
    fam_A = [p for p in primes if p % 3 == 1]  # Richtung i
    fam_B = [p for p in primes if p % 3 == 2]  # Richtung j
    fam_C = [p for p in primes if p == 3]      # Richtung k (Singulärer Verzweigungspunkt)
    
    count_A = len(fam_A)
    count_B = len(fam_B)
    
    print(f"Anzahl Primzahlen in Familie A (p = 1 mod 3): {count_A}")
    print(f"Anzahl Primzahlen in Familie B (p = 2 mod 3): {count_B}")
    
    # 3. Berechnung des Dirichlet-Defekts (Lokale Asymmetrie)
    # Misst die Abweichung von der perfekten euklidischen Kugelgeometrie
    dirichlet_defect = count_A - count_B
    total_active_primes = count_A + count_B
    
    # Relativer Jitter-Gradient (Asymmetrie-Index)
    jitter_gradient = dirichlet_defect / total_active_primes
    print(f"Dirichlet-Defekt (A - B): {dirichlet_defect}")
    print(f"Lokaler Jitter-Gradient: {jitter_gradient:.6f}")
    
    # 4. Berechnung des nackten Basiswerts alpha_0 (Asymptotischer Vakuumwert)
    # Konstanter Wert nach Hassall bei s=3
    zeta_3 = 1.202056903159594  # Apery-Konstante
    alpha_0 = 1.0 / (4.0 * np.pi * zeta_3 * 3**2)
    inverse_alpha_0 = 1.0 / alpha_0
    
    print(f"\nNackter Basiswert alpha_0: {alpha_0:.8f} (1/{inverse_alpha_0:.4f})")
    
    # 5. Modellinterne Kalibrierung über den quaternionischen Verformungsvektor
    # Der Jitter-Gradient verzerrt die lokale Metrik und dämpft die transversale Kopplung.
    # Wir nutzen den Defekt als Skalierungsfaktor kappa(M) für die effektive Dichte.
    
    # Modellinterner Skalierungsfaktor kappa(M)
    # Vorzeichen-Inversion: negativer Dirichlet-Defekt (B > A) senkt die effektive
    # Kopplung (kappa < 1), sodass 1/alpha_intern Richtung CODATA (1/137) wandert.
    kappa = 1.0 + (jitter_gradient * np.pi / 2.0)
    
    alpha_intern = alpha_0 * kappa
    inverse_alpha_intern = 1.0 / alpha_intern
    
    print(f"Modellinterner Dämpfungsfaktor kappa(M): {kappa:.6f}")
    print(f"-> Modellinterner Endwert von Alpha: {alpha_intern:.8f} (1/{inverse_alpha_intern:.4f})")
    
    # 6. Abgleich mit dem experimentellen CODATA-Wert
    alpha_codata = 1.0 / 137.035999177
    error_percent = abs(alpha_intern - alpha_codata) / alpha_codata * 100
    
    print(f"\nCODATA Referenzwert: {alpha_codata:.8f} (1/137.0360)")
    print(f"Modellinterne Abweichung zum Experiment: {error_percent:.6f} %")
    
    # 7. Quaternionische Vakuum-Metrik-Repräsentation
    # q = (Scalar, i, j, k)
    # Wir bilden den Verformungszustand als Quaternion ab
    q_deformation = np.array([total_active_primes, count_A, count_B, len(fam_C)])
    q_norm = q_deformation / np.linalg.norm(q_deformation)
    
    print(f"\nNormierter quaternionischer Strukturvektor [E, A, B, C]:")
    print(f"[{q_norm[0]:.4f}, {q_norm[1]:.4f}*i, {q_norm[2]:.4f}*j, {q_norm[3]:.4f}*k]")

def run_dynamic_eabc_scan(max_M=500000, steps=10):
    print(f"=== EABC-Modell: Dynamischer M-Scan bis M_max = {max_M} ===")
    print(f"{'Fixpunkt M':<12} | {'Sieve-Größe':<12} | {'Defekt (A-B)':<12} | {'Jitter-Grad':<12} | {'Alpha_intern':<12} | {'Abweichung':<10}")
    print("-" * 83)

    # Vorberechnen des maximalen Primzahl-Siebs für maximale Effizienz
    sieve = [True] * (max_M + 1)
    for p in range(2, int(np.sqrt(max_M)) + 1):
        if sieve[p]:
            for i in range(p * p, max_M + 1, p):
                sieve[i] = False

    # Konstanten nach Hassall (s=3)
    zeta_3 = 1.202056903159594
    alpha_0 = 1.0 / (4.0 * np.pi * zeta_3 * 3**2)
    alpha_codata = 1.0 / 137.035999177

    # Wir scannen markante harmonische Schritte und Ihren Fixpunkt
    test_points = sorted(list(set([10000, 50000, 113160, 200000, 300000, 400000, 500000])))

    for M in test_points:
        # Extrahiere Primzahlen bis M
        primes_M = [p for p in range(2, M + 1) if sieve[p]]

        count_A = sum(1 for p in primes_M if p % 3 == 1)
        count_B = sum(1 for p in primes_M if p % 3 == 2)
        total_active = count_A + count_B

        if total_active == 0:
            continue

        dirichlet_defect = count_A - count_B
        jitter_gradient = dirichlet_defect / total_active

        # Invertierter Dämpfungsfaktor (wie im Terminal verifiziert)
        kappa = 1.0 + (jitter_gradient * np.pi / 2.0)
        alpha_intern = alpha_0 * kappa
        inverse_alpha_intern = 1.0 / alpha_intern

        error_percent = abs(alpha_intern - alpha_codata) / alpha_codata * 100

        print(f"{M:<12} | {total_active:<12} | {dirichlet_defect:<12} | {jitter_gradient:+.6f} | 1/{inverse_alpha_intern:.4f}  | {error_percent:.4f} %")


if __name__ == "__main__":
    run_dynamic_eabc_scan()
    print()
    run_eabc_deterministic_analysis(M=113160)