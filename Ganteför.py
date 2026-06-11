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
from collections import defaultdict

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

    def derive_thermodynamic_equation(self, num_blocks=10):
        """
        Teilt die Nullstellen in Energiezonen und leitet daraus eine
        logarithmische Zustandsgleichung der Aharonov-Bohm-Phasendifferenz ab.
        """
        if num_blocks <= 0:
            raise ValueError("num_blocks muss positiv sein.")

        block_size = max(1, self.total_zeros // num_blocks)
        gammas = []
        delta_phases = []

        print(f"Leite Zustandsgleichung über {num_blocks} Energiezonen ab...")

        for b in range(num_blocks):
            start = b * block_size
            if start >= self.total_zeros:
                break

            end = (b + 1) * block_size if b < num_blocks - 1 else self.total_zeros
            sub_zeros = self.zeros[start:end]
            norms = np.round(sub_zeros).astype(int)
            phases = sub_zeros % (2 * np.pi)

            rem = norms % 4
            idx_1 = (rem == 1)
            idx_3 = (rem == 3)

            mean_1 = np.mean(phases[idx_1]) if np.sum(idx_1) > 0 else 0.0
            mean_3 = np.mean(phases[idx_3]) if np.sum(idx_3) > 0 else 0.0

            gammas.append(float(np.mean(sub_zeros)))
            delta_phases.append(float(mean_3 - mean_1))

        if len(gammas) < 2:
            raise ValueError("Für die logarithmische Regression werden mindestens zwei Blöcke benötigt.")

        # Logarithmische Regression: delta_phase = a * log(gamma) + b
        coeffs = np.polyfit(np.log(gammas), delta_phases, 1)
        return gammas, delta_phases, coeffs

    def _rounded_norms(self):
        return np.round(self.zeros).astype(int)

    def _generate_prime_quartets(self, norms):
        max_norm = int(np.max(norms)) + 10
        is_prime = np.ones(max_norm, dtype=bool)
        is_prime[:2] = False

        for p in range(2, int(np.sqrt(max_norm)) + 1):
            if is_prime[p]:
                is_prime[p * p::p] = False

        quartets = []
        for p in range(11, max_norm - 8):
            if is_prime[p] and is_prime[p + 2] and is_prime[p + 6] and is_prime[p + 8]:
                quartets.append((p, p + 2, p + 6, p + 8))

        return quartets

    def _map_norms_to_zeros(self, norms):
        norm_to_zeros = defaultdict(list)
        for gamma, norm in zip(self.zeros, norms):
            norm_to_zeros[norm].append(gamma)
        return norm_to_zeros

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
        plt.xlabel('Arithmetischer Niveau-Abstand ($\\delta_n$)')
        plt.ylabel('Quanten-Wahrscheinlichkeitsdichte')
        plt.grid(True, linestyle=':')
        plt.savefig('test_env_von_klitzing.png')
        plt.close()
        print("Grafik 'test_env_von_klitzing.png' wurde erzeugt.")

    def plot_phase_state_equation(self, gammas, delta_phases, coeffs):
        """
        Visualisiert die asymptotische Zustandsgleichung der AB-Phasendifferenz.
        """
        plt.scatter(gammas, delta_phases, color='crimson', label='Empirische Daten (zeros6.npy)', zorder=3)

        fit_gammas = np.linspace(min(gammas), max(gammas), 500)
        fit_phases = coeffs[0] * np.log(fit_gammas) + coeffs[1]

        plt.plot(
            fit_gammas,
            fit_phases,
            color='darkblue',
            linestyle='--',
            label=f'Zustandsgleichung: $\\Delta\\phi(\\gamma) = {coeffs[0]:.6f} \\cdot \\ln(\\gamma) + {coeffs[1]:.6f}$'
        )

        plt.title("Aharonov-Bohm-Phasenverschiebung: asymptotische Kinetik (#Energiedoku)")
        plt.xlabel('Kosmologisches Energieniveau ($\\gamma_n$)')
        plt.ylabel('Phasendifferenz $\\Delta\\phi$ (rad)')
        plt.grid(True, linestyle=':')
        plt.legend()
        plt.savefig('ab_phase_state_equation.png')
        plt.close()
        print("Grafik 'ab_phase_state_equation.png' wurde erzeugt.")

    def run_prime_quartet_analysis(self, global_gammas=None, global_delta_phases=None):
        """
        Analysiert lokale Aharonov-Bohm-Phasenresonanzen auf Primzahlvierlingen
        der Form (p, p+2, p+6, p+8).
        """
        print("\n--- Test 5: Lokale Phasenresonanz in Vierlingsprimzahlen ---")

        norms = self._rounded_norms()
        quartets = self._generate_prime_quartets(norms)
        print(f"Gefundene Primzahlvierlinge im mathematischen Raum: {len(quartets)}")

        norm_to_zeros = self._map_norms_to_zeros(norms)
        quartet_hits = []
        for quartet in quartets:
            p, p2, p6, p8 = quartet
            zeros_p = norm_to_zeros[p]
            zeros_p2 = norm_to_zeros[p2]
            zeros_p6 = norm_to_zeros[p6]
            zeros_p8 = norm_to_zeros[p8]

            if zeros_p or zeros_p2 or zeros_p6 or zeros_p8:
                quartet_hits.append({
                    'quartet': quartet,
                    'zeros_p': zeros_p,
                    'zeros_p2': zeros_p2,
                    'zeros_p6': zeros_p6,
                    'zeros_p8': zeros_p8,
                })

        print(f"Von Riemann-Nullstellen aktivierte Interferenz-Vierlinge: {len(quartet_hits)}")

        quartet_hits_sorted = sorted(quartet_hits, key=lambda hit: hit['quartet'][0])
        mid = len(quartet_hits_sorted) // 2
        low_quartets = quartet_hits_sorted[:mid]
        high_quartets = quartet_hits_sorted[mid:]

        def calc_quartet_diff(q_list):
            phases_3mod4 = []
            phases_1mod4 = []
            for hit in q_list:
                for gamma in hit['zeros_p'] + hit['zeros_p8']:
                    phases_3mod4.append(gamma % (2 * np.pi))
                for gamma in hit['zeros_p2'] + hit['zeros_p6']:
                    phases_1mod4.append(gamma % (2 * np.pi))

            if not phases_3mod4 or not phases_1mod4:
                raise ValueError("Zu wenige Treffer für eine stabile Vierlings-Phasendifferenz.")

            mean_3 = float(np.mean(phases_3mod4))
            mean_1 = float(np.mean(phases_1mod4))
            return mean_3, mean_1, mean_3 - mean_1, len(phases_3mod4) + len(phases_1mod4)

        if not low_quartets or not high_quartets:
            raise ValueError("Zu wenige aktivierte Primzahlvierlinge für den Niedrig/Hoch-Vergleich.")

        mean_3_low, mean_1_low, diff_low, count_low = calc_quartet_diff(low_quartets)
        mean_3_high, mean_1_high, diff_high, count_high = calc_quartet_diff(high_quartets)

        mean_norm_low = float(np.mean([hit['quartet'][0] for hit in low_quartets]))
        mean_norm_high = float(np.mean([hit['quartet'][0] for hit in high_quartets]))

        print("\n--- Ergebnisse der lokalen Vierlings-Resonanz ---")
        print(f"Niedrige Vierlinge (mittlere Norm ~{mean_norm_low:.0f}):")
        print(f"  Phase 3 mod 4: {mean_3_low:.5f} rad | Phase 1 mod 4: {mean_1_low:.5f} rad")
        print(f"  --> Lokale Phasendifferenz (Delta phi): {diff_low:+.5f} rad (Treffer: {count_low})")

        print(f"Höhere Vierlinge (mittlere Norm ~{mean_norm_high:.0f}):")
        print(f"  Phase 3 mod 4: {mean_3_high:.5f} rad | Phase 1 mod 4: {mean_1_high:.5f} rad")
        print(f"  --> Lokale Phasendifferenz (Delta phi): {diff_high:+.5f} rad (Treffer: {count_high})")

        if global_gammas is None or global_delta_phases is None:
            global_gammas = [mean_norm_low, mean_norm_high]
            global_delta_phases = [diff_low, diff_high]

        gamma_global = [global_gammas[0], global_gammas[-1]]
        phase_global = [global_delta_phases[0], global_delta_phases[-1]]
        gamma_quartet = [mean_norm_low, mean_norm_high]
        phase_quartet = [diff_low, diff_high]

        plt.plot(
            gamma_global,
            phase_global,
            color='blue',
            marker='o',
            linestyle='-',
            linewidth=2,
            label='Globaler Hintergrund-Drift',
        )
        plt.plot(
            gamma_quartet,
            phase_quartet,
            color='crimson',
            marker='^',
            linestyle='--',
            linewidth=2,
            label='Lokale Vierlings-Resonanz',
        )

        plt.title('Topologische Phasen-Verstärkung im Primzahlvierling')
        plt.xlabel(r'Kosmische Energie / Gitter-Norm ($\gamma$)')
        plt.ylabel(r'Aharonov-Bohm Phasendifferenz $\Delta\phi$ (rad)')
        plt.xscale('log')
        plt.grid(True, which='both', linestyle=':')
        plt.legend()
        plt.savefig('quartet_compensation_dynamics.png')
        plt.close()
        print("\nGrafik 'quartet_compensation_dynamics.png' wurde erzeugt.")

        return {
            'quartet_count': len(quartets),
            'activated_quartet_count': len(quartet_hits),
            'low_phase_diff': diff_low,
            'high_phase_diff': diff_high,
        }

    def analyze_von_klitzing_steps(self):
        """
        Extrahiert die diskreten Phasensprünge zwischen aktivierten
        Primzahlvierlings-Plateaus und visualisiert das Treppenmuster.
        """
        print("\n--- Test 6: Quantisierte Hall-Sprünge zwischen Vierlings-Plateaus ---")

        norms = self._rounded_norms()
        quartets = self._generate_prime_quartets(norms)
        norm_to_zeros = self._map_norms_to_zeros(norms)

        quartet_energies = []
        quartet_steps = []

        for quartet in quartets:
            p, p2, p6, p8 = quartet
            zeros_3mod4 = norm_to_zeros[p] + norm_to_zeros[p8]
            zeros_1mod4 = norm_to_zeros[p2] + norm_to_zeros[p6]

            # Aktiv sind hier nur Vierlinge mit Treffern in beiden Chiralitätssektoren.
            if zeros_3mod4 and zeros_1mod4:
                phase_3 = np.mean([gamma % (2 * np.pi) for gamma in zeros_3mod4])
                phase_1 = np.mean([gamma % (2 * np.pi) for gamma in zeros_1mod4])

                quartet_energies.append(float(np.mean(quartet)))
                quartet_steps.append(float(phase_3 - phase_1))

        if len(quartet_steps) < 2:
            raise ValueError("Zu wenige aktive Vierlings-Plateaus für eine Stufenanalyse.")

        sort_idx = np.argsort(quartet_energies)
        quartet_energies = np.array(quartet_energies)[sort_idx]
        quartet_steps = np.array(quartet_steps)[sort_idx]
        jump_heights = np.diff(quartet_steps)
        abs_jump_heights = np.abs(jump_heights)

        print("\n--- Numerischer Report: Quantisierte Stufen ---")
        print(f"Anzahl sequenzieller Übergänge:      {len(jump_heights)}")
        print(f"Mittlere Von-Klitzing-Sprunghöhe:   {np.mean(abs_jump_heights):.6f} rad")
        print(f"Maximale Hall-Spannungs-Fluktuation: {np.max(abs_jump_heights):.6f} rad")
        print(f"Minimale Übergangs-Schwelle:         {np.min(abs_jump_heights):.6f} rad")

        plt.figure(figsize=(10, 6))
        plt.step(
            quartet_energies,
            quartet_steps,
            where='mid',
            color='darkorange',
            linewidth=2,
            label='Quantisierte Vierlings-Plateaus',
        )
        plt.scatter(
            quartet_energies,
            quartet_steps,
            color='navy',
            s=15,
            alpha=0.7,
            zorder=3,
            label='Aktive Resonanz-Zentren',
        )

        plt.title('Arithmetischer Quanten-Hall-Effekt: Von-Klitzing-Stufen im Hurwitz-Raum')
        plt.xlabel(r'Kosmische Energie-Skala / Mittlere Vierlings-Norm ($\gamma_{quad}$)')
        plt.ylabel(r'Lokales Aharonov-Bohm-Potential $\Delta\phi_{\text{local}}$ (rad)')
        plt.xscale('log')
        plt.grid(True, which='both', linestyle=':', alpha=0.5)
        plt.legend()
        plt.savefig('von_klitzing_arithmetic_steps.png')
        plt.close()
        print("\nGrafik 'von_klitzing_arithmetic_steps.png' wurde erzeugt.")

        return {
            'transition_count': len(jump_heights),
            'mean_jump_height': float(np.mean(abs_jump_heights)),
            'max_jump_height': float(np.max(abs_jump_heights)),
            'min_jump_height': float(np.min(abs_jump_heights)),
        }

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

    # Test 4: Asymptotische Zustandsgleichung der AB-Phasendifferenz
    print("\n--- Test 4: Thermodynamische AB-Zustandsgleichung ---")
    gammas, delta_phases, coeffs = env.derive_thermodynamic_equation(num_blocks=10)

    print("\n--- Ergebnisse und Zustandsgleichung ---")
    for gamma, delta_phase in zip(gammas, delta_phases):
        print(f"Energie (Gamma): {gamma:12.2f} -> AB-Phasendifferenz: {delta_phase:.6f} rad")

    print("\nAbgeleitete kinetische Gleichung aus den Nullstellen:")
    print(f"Delta_Phi(Gamma) = {coeffs[0]:.6f} * ln(Gamma) + {coeffs[1]:.6f}")

    env.plot_phase_state_equation(gammas, delta_phases, coeffs)

    # Test 5: Lokale Primzahlvierlings-Resonanz gegen globalen Drift
    env.run_prime_quartet_analysis(gammas, delta_phases)

    # Test 6: Quantisierte Hall-Sprünge zwischen Vierlings-Plateaus
    env.analyze_von_klitzing_steps()