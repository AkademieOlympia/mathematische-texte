#!/usr/bin/env python3
"""
E8 Shell Spacing and Fano-Projection Audit for Prime Quadruplets
Author: EABC / Kepler-Hurwitz Framework
"""

import numpy as np
import matplotlib.pyplot as plt

def is_prime(n):
    if n < 2: return False
    if n in (2, 3): return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_prime_quadruplets(limit_k):
    """
    Finds prime quadruplets Q(p) = (p, p+2, p+6, p+8) 
    for p of the form 30k + 11 up to limit_k.
    """
    quadruplets = []
    for k in range(limit_k):
        p = 30 * k + 11
        if is_prime(p) and is_prime(p+2) and is_prime(p+6) and is_prime(p+8):
            quadruplets.append((k, p))
    return quadruplets

def run_audit():
    print("==========================================================")
    print(" E8 SHELL & FANO ENERGY LEVEL SPACING AUDIT (LOCAL M-MAC)")
    print("==========================================================")
    
    limit_k = 50000
    print(f"[Suche] Scanne Primvierlings-Anker p = 30k + 11 bis k = {limit_k:,} ...")
    quads = find_prime_quadruplets(limit_k)
    print(f"[Erfolg] Gefundene Primvierlinge: {len(quads):,}")
    
    if len(quads) < 10:
        print("[Hinweis] Zu wenige Punkte für statistische Analyse.")
        return

    energies = []
    for k, p in quads:
        n_p = 2 * p**2 + 16 * p + 52
        e_total = 6 * n_p
        energies.append(e_total)

    energies = np.array(energies, dtype=np.float64)
    
    # Nächste-Nachbarn-Abstände
    spacings = np.diff(energies)
    mean_spacing = np.mean(spacings)
    normalized_spacings = spacings / mean_spacing
    
    spacing_std = np.std(normalized_spacings)
    spacing_median = np.median(normalized_spacings)

    print("\n" + "-"*50)
    print(" STATISTISCHE ENERGIE-SPEKTRUM AUSWERTUNG")
    print("-" * 50)
    print(f" Anzahl analysierte Energieniveaus    : {len(energies):,}")
    print(f" Mittlerer Niveau-Abstand <sə>        : {mean_spacing:.2f}")
    print(f" Standardabweichung der norm. Abstände  : {spacing_std:.4f}")
    print(f" Median der norm. Abstände            : {spacing_median:.4f}")
    print("-" * 50)

    # Matplotlib Visualisierung
    print("[Plot] Generiere Level-Spacing-Histogramm...")
    plt.figure(figsize=(8, 5))
    plt.hist(normalized_spacings, bins=50, density=True, alpha=0.7, color='indigo', edgecolor='black', label='Empirische Spacing-Verteilung')
    
    # Zum Vergleich: Poisson-Verteilung (ungeordnete Unabhängigkeit) vs Wigner-Dyson (GUE)
    s = np.linspace(0, 3, 100)
    poisson = np.exp(-s)
    wigner_gue = (32 / np.pi**2) * s**2 * np.exp(-4 * s**2 / np.pi)
    
    plt.plot(s, poisson, 'r--', lw=2, label='Poisson-Verteilung (Zufall)')
    plt.plot(s, wigner_gue, 'g-', lw=2, label='Wigner-Dyson GUE (Quantenchaos)')
    
    plt.title(r'$E_8$ Fano-Energy Level Spacing Distribution ($N = ' + str(len(quads)) + '$)')
    plt.xlabel('Normierter Niveauabstand $s = \\Delta E / \\langle \\Delta E \\rangle$')
    plt.ylabel('Dichte $P(s)$')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()
    
    output_path = 'scripts/e8_level_spacing_plot.png'
    plt.savefig(output_path, dpi=300)
    print(f"[Erfolg] Plot gespeichert unter: {output_path}")
    plt.show()

if __name__ == "__main__":
    run_audit()