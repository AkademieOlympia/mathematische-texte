import numpy as np
from scipy.stats import pearsonr

def calculate_smoking_gun_correlation(diff_map, y_edges, top_modes):
    """
    Berechnet die Korrelation zwischen der Vakuumpolarisation 
    und den theoretischen Riemann-Oszillationen.
    
    diff_map: Die 2D-Differenz-Matrix (ABCE - CEAB)
    y_edges: Die vertikalen Koordinaten der Heatmap
    top_modes: Liste von (gamma, amplitude) aus deinem Log
    """
    # 1. Extraktion des vertikalen Profils (Signal)
    # Wir summieren horizontal, um das Hauptsignal der rechten Kante zu isolieren
    signal = np.sum(diff_map, axis=1) 
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    
    # 2. Erzeugung des theoretischen Riemann-Signals
    # Wir nutzen die Superposition der stärksten Moden
    theory = np.zeros_like(y_centers)
    for gamma, amp, cos_val, sin_val in top_modes:
        # Die Phase wird aus deinen cos/sin Komponenten im Log bestimmt
        theory += amp * (cos_val * np.cos(gamma * y_centers) + 
                         sin_val * np.sin(gamma * y_centers))
    
    # 3. Normalisierung
    signal = (signal - np.mean(signal)) / np.std(signal)
    theory = (theory - np.mean(theory)) / np.std(theory)
    
    # 4. Pearson-Korrelation
    r_coeff, p_value = pearsonr(signal, theory)
    
    return r_coeff, p_value, y_centers, signal, theory

# Daten aus deinem Log einsetzen:
my_top_modes = [
    (14.1347, 1.1785e-02, 8.63e-03, 8.02e-03), # gamma1
    (204.1896, 1.1721e-02, 3.87e-03, -1.10e-02), # Beispiel weitere Mode
    # ... hier weitere Moden aus deinem Log ergänzen
]

# r, p, y, s, t = calculate_smoking_gun_correlation(diff_map, y_edges, my_top_modes)
# print(f"Smoking Gun Korrelation r = {r:.4f} (p-Wert: {p_value:.3e})")