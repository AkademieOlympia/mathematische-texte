import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def analyze_zeros_file(file_path):
    # 1. VERSUCH: Datei als binary float64 lesen (Standard für solche Datensätze)
    try:
        # Odlyzko-Dateien sind oft einfach rohe Doubles
        zeros = np.fromfile(file_path, dtype=np.float64)
        
        # Plausibilitäts-Check: Die erste Nullstelle muss ca. 14.1347 sein
        # Wenn die Werte riesig sind, ist es ein Fenster bei hoher Energie.
        # Wenn sie klein sind, fängt es bei n=1 an.
        print(f"Anzahl geladener Nullstellen: {len(zeros)}")
        print(f"Erste 5 Werte: {zeros[:5]}")
        
        is_start_of_zeta = False
        if 14.0 < zeros[0] < 14.2:
            print("Verifiziert: Datensatz beginnt bei der ersten Nullstelle.")
            is_start_of_zeta = True
        else:
            print(f"Datensatz scheint ein Ausschnitt zu sein (Startwert: {zeros[0]}).")
            # Wir fahren trotzdem fort, die Resonanz ist global.

    except Exception as e:
        return f"Fehler beim Lesen der Datei: {e}"

    # Wir nutzen eine Untermenge für den Scan, um Rechenzeit zu sparen,
    # aber genug für Statistik (z.B. 100.000 Werte).
    # Wenn wir zu viele nehmen, wird das Rauschen bei falschem Alpha zu stark.
    sample_size = min(len(zeros), 100000) 
    zeros_subset = zeros[:sample_size]

    # 2. DEFINITION DER KOSTENFUNKTION (Bamberg-Resonanz)
    # Ziel: Finde Alpha, sodass (E_n * ln(1/alpha)) mod 2pi konstant ist (oder minimal streut)
    
    def get_variance_for_alpha(alpha_inv, zero_data):
        # Der Faktor X = ln(alpha_inv)
        factor = np.log(alpha_inv)
        # Phase = E * factor
        phases = zero_data * factor
        # Modulo 2pi
        phases_mod = (phases % (2 * np.pi))
        
        # Wir wollen, dass sich diese Phasen bei einem bestimmten Wert (Delta) häufen.
        # Eine gute Metrik ist 1 minus die Länge des Vektormittelwerts (Kuramoto-Order-Parameter).
        # Wenn alle Phasen gleich sind, ist R = 1. Wenn sie random sind, ist R = 0.
        # Wir wollen R maximieren (bzw. 1-R minimieren).
        
        # Umwandlung in komplexe Zeiger auf dem Einheitskreis
        complex_phases = np.exp(1j * phases)
        order_parameter = np.abs(np.mean(complex_phases))
        
        # Wir geben (1 - Ordnung) zurück, da Optimizer minimieren wollen.
        return 1.0 - order_parameter

    # 3. DER SCAN (Grid Search)
    # Wir scannen um den physikalischen Wert 137.036
    # Bereich: 136.0 bis 138.0
    print("Starte Scan um Alpha^-1 = 137...")
    
    scan_range = np.linspace(136.0, 138.0, 2000) # Hohe Auflösung
    costs = []
    
    for val in scan_range:
        c = get_variance_for_alpha(val, zeros_subset)
        costs.append(c)
        
    costs = np.array(costs)
    
    # Finde das Minimum
    min_idx = np.argmin(costs)
    best_alpha_inv = scan_range[min_idx]
    best_score = 1.0 - costs[min_idx] # Das ist der "Resonanz-Faktor" (0 bis 1)

    # 4. BESTIMMUNG VON DELTA (PHASENVERSATZ) für den besten Alpha-Wert
    factor_best = np.log(best_alpha_inv)
    phases_best = (zeros_subset * factor_best) % (2 * np.pi)
    
    # Histogramm der Phasen
    hist, bin_edges = np.histogram(phases_best, bins=100, density=True)
    
    # Visualisierung vorbereiten
    plt.figure(figsize=(10, 5))
    plt.plot(scan_range, 1.0 - costs) # Plotten der "Ordnung" (höher ist besser)
    plt.axvline(137.035999, color='r', linestyle='--', label='CODATA 1/alpha')
    plt.axvline(best_alpha_inv, color='g', linestyle='-', label=f'Best Fit: {best_alpha_inv:.4f}')
    plt.title("Bamberg-Resonanz-Scan: Ordnungsparameter vs. 1/Alpha")
    plt.xlabel("Inverse Feinstrukturkonstante 1/alpha")
    plt.ylabel("Resonanz-Stärke (0=Chaos, 1=Perfekt)")
    plt.legend()
    plt.grid(True)
    
    return {
        "best_alpha_inv": best_alpha_inv,
        "resonance_score": best_score,
        "n_samples": sample_size,
        "first_zero": zeros[0],
        "plot_obj": plt
    }

# Ausführung mit der hochgeladenen Datei
results = analyze_zeros_file("zeros6.gz/zeros6")
print(results)