from sage.all import *
import numpy as np

def is_mersenne_prime_ll(p):
    """
    Lucas-Lehmer-Test für M_p = 2^p - 1.
    [cite: 900, 1002]
    """
    if p == 2: return True
    if not is_prime(p): return False
    
    # Der Lucas-Lehmer-Algorithmus
    s = 4
    m = 2**p - 1
    for _ in range(p - 2):
        s = (s*s - 2) % m
    return s == 0

def bamberger_candidate_search(iterations, k, p0):
    """
    Kombiniert Bhattacharyas Bayes-Logik mit dem quaternionischen Filter.
   
    """
    current_t = 1.0 # Startwert im TMCMC-Raum
    found_primes = []
    mersenne_kandidaten = []

    for i in range(iterations):
        # 1. TMCMC-Schritt (Transformation-based MCMC) [cite: 17, 91]
        current_t = tmcmc_step_bm(current_t, k, p0)
        
        # 2. Primzahlprüfung (Miller-Rabin) [cite: 889]
        candidate_p = int(current_t + p0)
        
        # 3. Energiedoku-Check: Nur p aus den quaternionischen Klassen A,B,C,D
        if candidate_p % 12 in [1, 5, 7, 11]:
            if is_prime(candidate_p):
                if candidate_p not in found_primes:
                    found_primes.append(candidate_p)
                    
                    # 4. Mersenne-Validierung [cite: 900, 927]
                    print(f"BM-Resonanz gefunden: p = {candidate_p}")
                    if is_mersenne_prime_ll(candidate_p):
                        print(f"!!! NEUE MERSENNE-PRIMZAHL ENTDECKT: 2^{candidate_p}-1 !!!")
                        mersenne_kandidaten.append(candidate_p)
                        
    return mersenne_kandidaten

# Parameter basierend auf Bhattacharyas Mersenne-Kandidaten [cite: 945]
p0_start = 140000000
k_eff = p0_start / log(p0_start) 

# Suche starten
# bamberger_candidate_search(1000, k_eff, p0_start)