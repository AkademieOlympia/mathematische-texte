# SageMath/Python Skript zur Verifikation des Divergenz-Kollapses
import numpy as np

def verify_divergence_collapse(n=1000):
    # Konstruktion von zwei diagonalen Operatoren (Dichtematrizen rho_T, rho_0)
    # in einem diskreten Hilbertraum H_E
    rng = np.random.default_rng(42)
    rho_T = rng.dirichlet(np.ones(n)) # Empirischer Transport
    rho_0 = rng.dirichlet(np.ones(n)) # Vakuum-Referenz
    
    # 1. Quantenrelative Entropie (Araki): Tr(rho_T(ln(rho_T) - ln(rho_0)))
    # Da rho diagonal, reduziert sich dies auf die diskrete KL-Divergenz
    araki_entropy = np.sum(rho_T * (np.log(rho_T) - np.log(rho_0)))
    
    # 2. Kullback-Leibler (Klassisch)
    kl_divergence = np.sum(rho_T * np.log(rho_T / rho_0))
    
    return araki_entropy, kl_divergence

# Ausführung
ara, kl = verify_divergence_collapse()
print(f"Araki-Entropie: {ara}")
print(f"Kullback-Leibler: {kl}")
print(f"Abweichung: {abs(ara - kl)}")