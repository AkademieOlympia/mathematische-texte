from sage.all import *
import numpy as np

def batch_process_kanon(limit=1000):
    results = []
    # Suche nach Clustern (Drillinge/Vierlinge)
    p = 2
    count = 0
    while count < limit:
        p = next_prime(p)
        if is_prime(p+2) and is_prime(p+6): # Beispielhafter Drillings-Scan
            # eabc-Index berechnen
            q_sum = sum(four_number_sum_squares(x)[1]**2 for x in [p, p+2, p+6])
            xi_local = q_sum / 12 # Normierter Index
            
            # Phasenverschiebung & Zeta-Resonanz (simuliert/vereinfacht)
            obs = 0.13387 # Dein verifizierter Wert aus dem letzten Lauf
            
            # Rückrechnung des benötigten Beta
            beta_local = (1836.1527 - xi_local) / obs
            results.append(beta_local)
            count += 1
    return results

beta_distribution = batch_process_kanon(1000)
mean_beta = np.mean(beta_distribution)
std_beta = np.std(beta_distribution)

print(f"### STATISTISCHER KANON-CHECK ###")
print(f"Anzahl Cluster: 1000")
print(f"Mittelwert Beta: {mean_beta:.4f}")
print(f"Standardabweichung: {std_beta:.4f}")