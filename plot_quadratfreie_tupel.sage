#!/usr/bin/env sage
"""
Visualisierung der gefundenen Zahlen
====================================
Liest die CSV-Datei und erstellt ein Diagramm mit:
- Primzahlfunktion π(n) angewendet auf die gefundenen Zahlen
- Logarithmus von n
"""

from sage.all import *
import csv

# Versuche matplotlib zu importieren
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Hinweis: matplotlib nicht verfügbar. Verwende SageMath's plot-Funktionen.")

def plot_quadratfreie_tupel(csv_file="quadratfreie_tupel_mod12.csv", output_file="quadratfreie_tupel_plot.png"):
    """
    Liest die CSV-Datei und erstellt ein Diagramm.
    
    Args:
        csv_file: Name der CSV-Datei
        output_file: Name der Ausgabedatei für das Diagramm
    """
    print(f"=== Visualisierung der gefundenen Zahlen ===")
    print(f"Lese Daten aus '{csv_file}'...\n")
    
    # Lese CSV-Datei
    numbers = []
    try:
        with open(csv_file, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                n_val = int(row["n"])
                numbers.append(n_val)
    except FileNotFoundError:
        print(f"✗ Fehler: Datei '{csv_file}' nicht gefunden!")
        return
    except Exception as e:
        print(f"✗ Fehler beim Lesen der CSV-Datei: {e}")
        return
    
    if not numbers:
        print("⚠️  Warnung: Keine Daten gefunden!")
        return
    
    print(f"Gefunden: {len(numbers)} Zahlen")
    print(f"Bereich: {min(numbers)} bis {max(numbers)}\n")
    
    # Sortiere Zahlen
    numbers.sort()
    
    # Berechne π(n) für jede Zahl (Anzahl der Primzahlen ≤ n)
    print("Berechne Primzahlfunktion π(n) für jede Zahl...")
    pi_values = []
    for n in numbers:
        pi_n = prime_pi(n)  # SageMath's prime_pi Funktion
        pi_values.append(pi_n)
    
    # Berechne Logarithmus von n
    print("Berechne Logarithmus von n...")
    log_values = [log(n) for n in numbers]
    
    # Erstelle Diagramm
    print(f"\nErstelle Diagramm...")
    
    if HAS_MATPLOTLIB:
        # Verwende matplotlib
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Diagramm 1: π(n) vs. n
        ax1.plot(numbers, pi_values, 'b-', linewidth=1.5, label='π(n)')
        ax1.set_xlabel('n (gefundene Zahl)', fontsize=12)
        ax1.set_ylabel('π(n) (Anzahl Primzahlen ≤ n)', fontsize=12)
        ax1.set_title('Primzahlfunktion π(n) angewendet auf gefundene Zahlen', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Diagramm 2: log(n) vs. n
        ax2.plot(numbers, log_values, 'r-', linewidth=1.5, label='log(n)')
        ax2.set_xlabel('n (gefundene Zahl)', fontsize=12)
        ax2.set_ylabel('log(n) (natürlicher Logarithmus)', fontsize=12)
        ax2.set_title('Logarithmus der gefundenen Zahlen', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Diagramm gespeichert nach '{output_file}'")
        
        try:
            plt.show()
        except:
            pass
    else:
        # Verwende SageMath's plot-Funktionen
        # Erstelle Liste von Punkten für π(n)
        points_pi = list(zip(numbers, pi_values))
        plot_pi = list_plot(points_pi, plotjoined=True, color='blue', 
                           axes_labels=['n', 'π(n)'],
                           title='Primzahlfunktion π(n) angewendet auf gefundene Zahlen',
                           figsize=(10, 5))
        
        # Erstelle Liste von Punkten für log(n)
        points_log = list(zip(numbers, log_values))
        plot_log = list_plot(points_log, plotjoined=True, color='red',
                            axes_labels=['n', 'log(n)'],
                            title='Logarithmus der gefundenen Zahlen',
                            figsize=(10, 5))
        
        # Kombiniere beide Plots
        combined_plot = graphics_array([plot_pi, plot_log], nrows=2, ncols=1)
        combined_plot.save(output_file, figsize=(12, 10), dpi=300)
        print(f"✓ Diagramm gespeichert nach '{output_file}'")
        
        # Zeige Plot
        show(combined_plot)
    
    # Zeige Statistiken
    print(f"\nStatistiken:")
    print(f"  Anzahl Zahlen: {len(numbers)}")
    print(f"  Minimum: {min(numbers)}")
    print(f"  Maximum: {max(numbers)}")
    print(f"  π(min): {prime_pi(min(numbers))}")
    print(f"  π(max): {prime_pi(max(numbers))}")
    print(f"  log(min): {log(min(numbers)):.4f}")
    print(f"  log(max): {log(max(numbers)):.4f}")

# Hauptausführung
if __name__ == "__main__":
    import sys
    
    csv_file = "quadratfreie_tupel_mod12.csv"
    output_file = "quadratfreie_tupel_plot.png"
    
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    
    plot_quadratfreie_tupel(csv_file, output_file)
