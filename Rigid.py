import numpy as np
from scipy.stats import kruskal

def generate_primes(start, stop):
    """Erzeugt Primzahlen im Bereich [start, stop) ohne externe Abhängigkeiten."""
    if stop <= 2:
        return []

    sieve = [True] * stop
    sieve[0] = False
    sieve[1] = False

    limit = int(stop ** 0.5) + 1
    for candidate in range(2, limit):
        if sieve[candidate]:
            first_multiple = candidate * candidate
            step = candidate
            sieve[first_multiple:stop:step] = [False] * len(range(first_multiple, stop, step))

    return [n for n in range(max(2, start), stop) if sieve[n]]

def calculate_f_B(p):
    """
    Berechnet f_B(p): Die exakte Anzahl der geordneten Stammbruch-Tripel
    unter der strikten Halbordnung x <= y <= z.
    """
    p = int(p)
    solutions_count = 0
    start_x = p // 4 + 1
    end_x = p // 3 + 1
    
    for x in range(start_x, end_x):
        num = 4 * x - p
        den = p * x
        start_y = den // num + 1
        end_y = 2 * den // num + 1
        
        for y in range(start_y, end_y):
            rem = num * y - den
            if rem > 0 and (den * y) % rem == 0:
                solutions_count += 1
                
    return solutions_count

def run_asymptotic_diophantine_test(max_x):
    """
    Führt den unbestechlichen Skalierungstest über die EABC-Klassen durch.
    Sortiert die echten Primzahlen in die vier Stichproben S_1, S_5, S_7, S_11
    und berechnet die Mittelwerte, Mediane und den Kruskal-Wallis p-Wert.
    """
    print(f"=== STARTE NUMERISCHEN TESTING-PROZESS BIS x = {max_x} ===")
    
    # Primzahlgenerierung ohne SageMath (ohne die 2 und 3)
    primes = generate_primes(5, max_x)
    
    # Stichproben-Container für die vier Modulo-12-Klassen
    samples = {1: [], 5: [], 7: [], 11: []}
    
    # Daten-Kollation
    for p in primes:
        rem12 = p % 12
        if rem12 in samples:
            f_B = calculate_f_B(p)
            samples[rem12].append(f_B)
            
    print(f"-> {len(primes)} Primzahlen erfolgreich verarbeitet und klassifiziert.\n")
    
    # Statistische Auswertung
    print(f"{'Klasse p mod 12':<15} | {'Anzahl (N)':<10} | {'Mittelwert (mu)':<15} | {'Median (~f)':<10}")
    print("-" * 62)
    
    for r in [1, 5, 7, 11]:
        data = samples[r]
        if data:
            mu = float(np.mean(data))
            median = float(np.median(data))
            print(f"{r:<15} | {len(data):<10} | {mu:<15.4f} | {median:<10.1f}")
            
    # Kruskal-Wallis-Test zur Absicherung gegen den Small-Number-Bias
    # Prüft, ob die Unterschiede zwischen den Gruppen statistisch signifikant sind
    try:
        stat, p_value = kruskal(samples[1], samples[5], samples[7], samples[11])
        print("\n===== STATISTISCHE ABSICHERUNG =====")
        print(f"Kruskal-Wallis H-Statistik: {stat:.4f}")
        print(f"Asymptotischer p-Wert:      {p_value:.6e}")
        
        if p_value < 0.05:
            print("-> Ergebnis: H0 wird VERWORFEN. Die Modulo-12-Anisotropie ist statistisch signifikant!")
        else:
            print("-> Ergebnis: H0 kann NICHT verworfen werden. Die Unterschiede sind statistisch nicht signifikant (Rauschen).")
    except Exception as e:
        print(f"\n[WARNUNG] Statistischer Test fehlgeschlagen: {e}")

# Dieser Befehl stößt den echten Lauf lokal an.
# Starte im Labor mit 2000 (zur Verifikation deiner Daten), dann 10^4, 10^5...
run_asymptotic_diophantine_test(max_x=2000)