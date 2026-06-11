# SageMath: RH-Stabilität für Modulformen (Hecke-Eigenwerte)
# Prüft die Ramanujan-Petersson-Bedingung |λ| ≤ 2√p

def certify_rh_stability(p_limit):
    M = ModularForms(1, 12).new_subspace()
    for p in primes(p_limit):
        evs = M.T(p).matrix().eigenvalues()
        if any(abs(ev) > 2 * sqrt(p) for ev in evs):
            return False
    return True
