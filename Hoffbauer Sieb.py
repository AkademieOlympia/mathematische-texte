import numpy as np
import pandas as pd
import time

def fast_prime_sieve(limit):
    """Ein hocheffizientes Sieb des Eratosthenes mit Numpy."""
    sieve = np.ones(limit // 3 + (limit % 6 == 2), dtype=bool)
    for i in range(1, int(limit**0.5) // 3 + 1):
        if sieve[i]:
            k = 3 * i + 1 | 1
            sieve[k * k // 3 :: 2 * k] = False
            sieve[k * (k - 2 * (i & 1) + 4) // 3 :: 2 * k] = False
    return np.r_[2, 3, ((3 * np.nonzero(sieve)[0][1:] + 1) | 1)]

def analyze_hoffbauer_sieve(limit, k_range=[0, 1, 2]):
    """
    Geometrisches Sieb zur Klassifikation von Viererordnungen Q_k.
    Implementiert Ptolemäus- und Walter-Metriken[cite: 72, 122].
    """
    print(f"Generiere Primzahlen bis {limit}...")
    start_time = time.time()
    primes = fast_prime_sieve(limit)
    prime_set = set(primes)
    print(f"Gefunden: {len(primes)} Primzahlen in {time.time() - start_time:.2f}s")

    results = []

    for k in k_range:
        m_k = 6 + 12 * k  # Dehnungsabstand [cite: 43]
        print(f"Analysiere Dehnungsstufe Q_{k} (m={m_k})...")
        
        # Suche nach Q_k(p) = (p, p+2, p+m_k, p+m_k+2) [cite: 44]
        for p in primes:
            if p + m_k + 2 > limit:
                break
            
            if (p+2 in prime_set) and (p+m_k in prime_set) and (p+m_k+2 in prime_set):
                # Geometrische Punkte (logarithmisch) [cite: 56, 57, 58]
                x1, x2, x3, x4 = np.log([p, p+2, p+m_k, p+m_k+2])
                
                # 1. Metrischer Kreuzdefekt (Ptolemäus) 
                delta_e = abs(-2 * (x3 - x2) * (x4 - x1))
                theory_e = (64 + 288*k + 288*k**2) / (p**2) # 
                
                # 2. Walter-Split AC [cite: 122]
                delta_w = (1/20) * (x1 + x3) * (x4 - x2)
                theory_w = (m_k * np.log(p)) / (10 * p) # [cite: 132]
                
                results.append({
                    'k': k,
                    'p': p,
                    'delta_E': delta_e,
                    'E_Asymptotik': theory_e,
                    'delta_W': delta_w,
                    'W_Asymptotik': theory_w,
                    'E_Resonanz': delta_e / theory_e if theory_e > 0 else 0
                })

    return pd.DataFrame(results)

# --- Testlauf ---
if __name__ == "__main__":
    # Testbereich bis 1 Million, Stufen Q0 bis Q2
    df = analyze_hoffbauer_sieve(1000000, k_range=[0, 1, 2])
    
    print("\nErgebnis der geometrischen Klassifikation (Auszug):")
    print(df.groupby('k')[['p', 'E_Resonanz']].head(3))
    
    # Speichern der Ergebnisse
    df.to_csv("hoffbauer_geometrie_test.csv", index=False)
    print("\nKlassifikationsvektoren wurden in 'hoffbauer_geometrie_test.csv' gespeichert.")