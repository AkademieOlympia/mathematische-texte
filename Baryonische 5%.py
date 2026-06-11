import numpy as np

def check_matter_ratio_at_timestamp(filepath='zeros6.npy'):
    print("--- #Energiedoku: Analyse des Materieverhältnisses ---")
    zeros = np.load(filepath)
    norms = np.round(zeros).astype(int)
    
    # Primzahl-Sieb für die Identifikation der baryonischen Periodizität
    max_norm = int(np.max(norms)) + 10
    is_p = np.ones(max_norm, dtype=bool)
    is_p[0] = is_p[1] = False
    for p in range(2, int(np.sqrt(max_norm)) + 1):
        if is_p[p]:
            is_p[p*p::p] = False
            
    # Extraktion der Primzahlvierlinge
    quartet_shells = set()
    for p in range(11, max_norm - 8):
        if is_p[p] and is_p[p+2] and is_p[p+6] and is_p[p+8]:
            quartet_shells.update([p, p+2, p+6, p+8])

    # Zählung der besetzten Zustände im Gitter
    omega_c_nodes = 0  # Dunkle Materie (Skelett/Halbzahlen)
    omega_b_nodes = 0  # Baryonische Materie (Periodische Vierlinge)
    
    for N in norms:
        # Im Hurwitz-Gitter korrespondieren die Normen der dichten Packung
        # mit den Filtereigenschaften des Selberg-Siebs
        if N in quartet_shells:
            omega_b_nodes += 1
        else:
            omega_c_nodes += 1

    # Berechnung des asymptotischen Quoten
    ratio = omega_c_nodes / omega_b_nodes
    
    print(f"Zustände im Selberg-Sieb (Dunkle Materie Omega_c): {omega_c_nodes}")
    print(f"Zustände in Periodischen Vierlingen (Baryonisch Omega_b): {omega_b_nodes}")
    print(f"--> Numerisch ermitteltes Verhältnis: {ratio:.4f}")
    
    # Weicht die Simulation vom Planck-Wert ab, zeigt dies den Grad 
    # der unvollständigen Kristallisation (Szenario C) an.
    return ratio

if __name__ == "__main__":
    check_matter_ratio_at_timestamp()