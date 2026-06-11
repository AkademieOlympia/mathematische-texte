import numpy as np
import time
from sympy import isprime, pollard_rho

def _log_cluster_embedding(N, k=0):
    """Liefert die logarithmische Einbettung des Viererclusters Q_k(N)."""
    m_k = 6 + 12 * k
    cluster = np.array([N, N + 2, N + m_k, N + m_k + 2], dtype=np.float64)
    return np.log(cluster)

def get_hoeffbauer_resonance(N, k=0):
    """Berechnet die ptolemäische Resonanz eines Clusters Q_k(N)."""
    x = _log_cluster_embedding(N, k)
    
    # Realer metrischer Kreuzdefekt (Ptolemäus) [cite: 72, 90]
    delta_e = abs(-2 * (x[2] - x[1]) * (x[3] - x[0]))
    
    # Theoretisches Skalengesetz nach Satz 7 [cite: 87, 95]
    theory_e = (64 + 288*k + 288*k**2) / (N**2)
    
    return delta_e / theory_e if theory_e > 0 else 0

def calculate_walter_split(N, k=0):
    """Berechnet den Walter-Split in der logarithmischen EABC-Einbettung."""
    x1, x2, x3, x4 = _log_cluster_embedding(N, k)
    return (1 / 20) * (x1 + x3) * (x4 - x2)

def expected_walter(N, k=0):
    """Erwartungswert des Walter-Splits nach dem Skalengesetz."""
    m_k = 6 + 12 * k
    return (m_k * np.log(N)) / (10 * N)

def calculate_morley_diff(N, k=0):
    """Misst die Symmetrie der beiden Randkanten als Morley-Güte."""
    x1, x2, x3, x4 = _log_cluster_embedding(N, k)
    left_edge = x2 - x1
    right_edge = x4 - x3
    return abs(left_edge - right_edge)

def passes_local_prime_filter(N, small_primes=(3, 5, 7, 11, 13, 17, 19, 23, 29)):
    """
    Verwirft offensichtliche Strukturbrüche durch kleine modulare Obstruktionen.
    Das ist kein vollständiger Primtest, aber ein wirksamer lokaler Sperrfilter.
    """
    if N < 2:
        return False
    if N == 2:
        return True
    if N % 2 == 0:
        return False
    for p in small_primes:
        if N == p:
            return True
        if N % p == 0:
            return False
    return True

def h_factor_v2(N, k=0):
    # Ebene 1: Ptolemäische Resonanz (Globaler Rahmen)
    res_p = get_hoeffbauer_resonance(N, k)
    
    # Ebene 2: Walter-Stabilität (Lastverteilung)
    # Berechne die Flächendifferenz in der EABC-Einbettung
    walter_val = calculate_walter_split(N, k)
    res_w = walter_val / expected_walter(N, k)
    
    # Ebene 3: Morley-Güte (Gestaltreife)
    morley_val = calculate_morley_diff(N, k)
    local_filter_ok = passes_local_prime_filter(N)
    
    # Gesamt-Urteil: Nur wenn alle drei Resonanzen harmonieren
    if (
        abs(1 - res_p) < 1e-4
        and abs(1 - res_w) < 1e-2
        and morley_val < 1e-5
        and local_filter_ok
    ):
        return "GEOMETRISCH REIF (Prim-Struktur)"
    else:
        return "STRUKTURBRUCH (Zusammengesetzt/Instabil)"

def run_benchmark(test_numbers):
    print(f"{'Zahl N':>15} | {'Prim?':>6} | {'Resonanz':>10} | {'GDA-Urteil':>42} | {'Pollard-Rho Time'}")
    print("-" * 110)
    
    for N in test_numbers:
        # 1. Geometrische Analyse (GDA)
        start_gda = time.time()
        res = get_hoeffbauer_resonance(N, k=0)
        gda_judgment = h_factor_v2(N, k=0)
        gda_time = time.time() - start_gda
        
        # 2. Mainstream-Referenz (Pollard-Rho)
        start_pollard = time.time()
        actual_prime = isprime(N)
        if not actual_prime:
            factor = pollard_rho(N)
        pollard_time = time.time() - start_pollard
        
        print(f"{N:15d} | {str(actual_prime):6} | {res:10.6f} | {gda_judgment:42} | {pollard_time:.6f}s")

# Test-Set: Mischung aus Primzahlen und zusammengesetzten Zahlen
test_set = [
    1000003,            # Primzahl
    1000005,            # Zusammengesetzt (5 * ...)
    999999999999989,    # Große Primzahl
    999999999999991,    # Zusammengesetzt
    12345678901234567   # Zusammengesetzt
]

if __name__ == "__main__":
    run_benchmark(test_set)