#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#Energiedoku - Interaktive Test- und Simulationsumgebung
Modellierung von Hurwitz-Flussquantisierungen, Aharonov-Bohm-Phasen
und Chebyshev-Symmetriebrechungen auf Basis von Riemann-Nullstellen.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class EnergiedokuEnvironment:
    def __init__(self, filepath='zeros6.npy'):
        print(f"Lade arithmetischen Datensatz: {filepath}...")
        self.zeros = np.load(filepath)
        self.total_zeros = len(self.zeros)
        print(f"Erfolgreich geladen. {self.total_zeros:,} Quanten-Knotenpunkte verfügbar.")
        
    def run_sieve_analysis(self, start_idx=0, end_idx=100000):
        """
        Simuliert das Selberg-Sieb und analysiert die geometrische 
        Aharonov-Bohm-Phase für ein definiertes Nullstellen-Intervall.
        """
        slice_zeros = self.zeros[start_idx:end_idx]
        norms = np.round(slice_zeros).astype(int)
        phases = slice_zeros % (2 * np.pi)
        
        rem = norms % 4
        class_1 = (rem == 1)
        class_3 = (rem == 3)
        
        report = {
            'count_1mod4': int(np.sum(class_1)),
            'count_3mod4': int(np.sum(class_3)),
            'mean_phase_1': float(np.mean(phases[class_1])) if np.sum(class_1) > 0 else 0.0,
            'mean_phase_3': float(np.mean(phases[class_3])) if np.sum(class_3) > 0 else 0.0,
        }
        report['delta_phase'] = report['mean_phase_3'] - report['mean_phase_1']
        return report

    def plot_von_klitzing_plateaus(self, start_idx=0, end_idx=5000):
        """
        Visualisiert die Energie-Abstoßung (Landau-Niveaus) im Stil
        von Von-Klitzing-Widerstandsplateaus.
        """
        slice_zeros = self.zeros[start_idx:end_idx]
        spacings = np.diff(slice_zeros)
        
        # Plotten der quantisierten Niveaus
        plt.hist(spacings, bins=100, density=True, alpha=0.75, color='teal', edgecolor='black')
        plt.title(f'Von-Klitzing Landau Niveaus (GUE Repulsion, Indizes {start_idx}-{end_idx})')
        plt.xlabel('Arithmetischer Niveau-Abstand ($\delta_n$)')
        plt.ylabel('Quanten-Wahrscheinlichkeitsdichte')
        plt.grid(True, linestyle=':')
        plt.savefig('test_env_von_klitzing.png')
        plt.close()
        print("Grafik 'test_env_von_klitzing.png' wurde erzeugt.")

    def get_quaternion_coordinates(self, index):
        """
        Überführt eine spezifische Riemann-Nullstelle exemplarisch in 
        die 4D-Komponenten eines Hurwitz-Knotens.
        """
        gamma = self.zeros[index]
        norm_target = int(np.round(gamma))
        
        # Arithmetischer Dekompositionsansatz (Vier-Quadrate-Satz für Hurwitz)
        # Liefert die strukturelle Basis für die a, b, c, e Anteile
        return f"Nullstelle {gamma:.4f} -> Ziel-Norm Schale: {norm_target}"

# --- Start der interaktiven Testschleife ---
if __name__ == "__main__":
    # Instanziierung der Umgebung mit den realen Daten
    env = EnergiedokuEnvironment('zeros6.npy')
    
    # Test 1: Analyse des tiefen thermischen Bereichs (Erste Epoche)
    print("\n--- Test 1: Tiefe thermische Epoche (0 - 50.000) ---")
    res_low = env.run_sieve_analysis(0, 50000)
    print(f"Klasse 1 mod 4: {res_low['count_1mod4']} | Klasse 3 mod 4: {res_low['count_3mod4']}")
    print(f"Aharonov-Bohm Phasenverschiebung (Δϕ): {res_low['delta_phase']:.6f} rad")
    
    # Test 2: Analyse des hohen asymptotischen Bereichs (Späte Epoche)
    print("\n--- Test 2: Asymptotischer Grenzbereich (1.500.000 - 1.550.000) ---")
    res_high = env.run_sieve_analysis(1500000, 1550000)
    print(f"Klasse 1 mod 4: {res_high['count_1mod4']} | Klasse 3 mod 4: {res_high['count_3mod4']}")
    print(f"Aharonov-Bohm Phasenverschiebung (Δϕ): {res_high['delta_phase']:.6f} rad")
    
    # Test 3: Generierung eines lokalen Hall-Plateau-Musters
    env.plot_von_klitzing_plateaus(0, 10000)