import numpy as np

def correlate_eigenvalues_to_riemann(eigenvalues, gamma_list):
    """
    Findet die Kopplung zwischen eabc-Eigenwerten und Riemann-Nullstellen.
    """
    # 1. Normalisierung der Eigenwerte (Energieskalierung)
    # Wir setzen den stärksten Eigenwert in Relation zur ersten Nullstelle
    norm_evals = eigenvalues / np.max(eigenvalues)
    
    # 2. Die 'theoretische Energie' der Nullstellen
    # In der Quanten-Zahlentheorie ist E ~ log(gamma) oder gamma selbst
    riemann_energies = np.array([np.log(g) for g in gamma_list[:len(eigenvalues)]])
    norm_riemann = riemann_energies / np.max(riemann_energies)
    
    # 3. Mapping-Analyse
    mapping = []
    for i, e_val in enumerate(norm_evals):
        # Wir suchen die Nullstelle, deren relative Energie am nächsten liegt
        diffs = np.abs(norm_riemann - e_val)
        best_match_idx = np.argmin(diffs)
        
        mapping.append({
            'Dimension': i,
            'Eigenwert': eigenvalues[i],
            'Zugeordnete_Nullstelle': gamma_list[best_match_idx],
            'Match_Qualität': 1 - diffs[best_match_idx]
        })
        
    return mapping

# --- Anwendung ---
# gamma_top_8 = [14.1347, 21.0220, 25.0108, 30.4248, 32.9350, 37.5861, 40.9187, 43.3270]
# results = correlate_eigenvalues_to_riemann(evals, gamma_top_8)

print(f"{'Dim':<5} | {'Eigenwert':<12} | {'Gamma':<10} | {'Konfidenz'}")
print("-" * 45)
for res in results:
    print(f"{res['Dimension']:<5} | {res['Eigenwert']:<12.4e} | {res['Zugeordnete_Nullstelle']:<10.2f} | {res['Match_Qualität']:.2%}")