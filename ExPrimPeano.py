# SageMath 10.5 Skript für die #Energiedoku
# Korrelation von N=3333 Nullstellen mit G60 Facetten
# Läuft mit SageMath ODER mit matplotlib (pip install matplotlib)

import numpy as np

# SageMath oder Fallback (nur matplotlib, keine mpmath nötig)
try:
    from sage.all import zeta_zeros, bar_chart, list_plot, show
    USE_SAGE = True
except ModuleNotFoundError:
    import matplotlib.pyplot as plt
    USE_SAGE = False

# Erste 100 Riemann-Nullstellen (exakt, für Fallback ohne Sage)
_ZEROS_100 = [
    14.1347, 21.0220, 25.0109, 30.4249, 32.9351, 37.5862, 40.9187, 43.3271,
    48.0052, 49.7738, 52.9703, 56.4462, 59.3470, 60.8318, 65.1125, 67.0798,
    69.5464, 72.0672, 75.7047, 77.1448, 79.3374, 82.9104, 84.7355, 87.4253,
    88.8091, 92.4919, 94.6513, 95.8706, 98.8312, 101.3179, 103.7255, 105.4466,
    107.1686, 111.0295, 111.8747, 114.3202, 116.2270, 118.7908, 121.3701,
    122.9468, 124.2568, 127.5167, 129.5787, 131.0877, 133.4977, 134.7565,
    138.1165, 139.7362, 141.1238, 143.1118, 146.0009, 147.4228, 150.0535,
    150.9253, 153.0247, 156.0129, 157.5976, 158.8493, 161.1889, 163.0307,
    165.5371, 167.1844, 169.0945, 169.9119, 173.4115, 174.7542, 176.4414,
    178.3772, 179.9165, 182.2071, 184.8745, 185.5982, 187.2289, 189.4162,
    192.0267, 193.0797, 195.2654, 196.8765, 198.0153, 201.2647, 202.4936,
    204.1896, 205.3947, 207.9063, 209.5765, 211.6909, 213.3479, 214.5470,
    216.1695, 219.0676, 220.7145, 221.4307, 224.0070, 224.9833, 227.4214,
    229.3374, 231.2502, 231.9872, 233.6934, 236.5242, 237.7697, 239.5555,
    241.0492, 242.8233, 244.0719, 245.8835, 247.1369, 249.0522, 251.1517,
    253.0697, 253.6983, 256.4235, 257.9825, 259.7414, 261.4475, 263.0003,
]

def _approx_zero(n):
    """Asymptotische Näherung (Riemann-von Mangoldt): γ_n ≈ 2πn / ln(n/(2π))."""
    return 2 * np.pi * n / np.log(n / (2 * np.pi))

def _get_zeros(N_max):
    """Riemann-Nullstellen (Imaginärteile): Sage oder Fallback (exakt + Näherung)."""
    if USE_SAGE:
        zz = zeta_zeros()
        return [float(zz[i]) for i in range(N_max)]
    # Fallback: erste 100 exakt, Rest asymptotische Näherung (schnell)
    result = []
    for i in range(N_max):
        n = i + 1
        if n <= len(_ZEROS_100):
            result.append(_ZEROS_100[n - 1])
        else:
            result.append(_approx_zero(n))
    return result

def energiedoku_simulation(N_max=3333):
    G60 = 60
    alpha_inv = 137.03599
    
    # 1. Generierung der Riemann-Nullstellen (Imaginärteil gamma)
    print(f"Initialisiere Filter für N={N_max} Nullstellen...")
    gammas = _get_zeros(N_max)
    
    # 2. Arithmetische Peano-Transformation
    # Wir projizieren jede Nullstelle auf die 60 Facetten
    facette_counts = np.zeros(G60)
    resonances = []
    
    for i, gamma in enumerate(gammas):
        # Der Phasenwinkel im 60-Flächen-Raum
        # Die 137er Resonanz wirkt als Modulo-Filter
        phase = (float(gamma) * alpha_inv) % G60
        facette_index = min(int(phase), G60 - 1)  # Sicherheit gegen Rundungsfehler
        
        facette_counts[facette_index] += 1
        
        # Speichern der Resonanz-Punkte für die Grafik
        if i < 500: # Fokus auf die ersten 500 für die Übersicht
            resonances.append((i, phase))
            
    # 3. Visualisierung
    if USE_SAGE:
        p1 = bar_chart(list(facette_counts), title=f"Besetzung der 60 Facetten (N={N_max})",
                      axes_labels=['Facette Index', 'Anzahl Resonanzen'], rgbcolor=(0, 0, 0.55))
        p2 = list_plot(resonances, plotjoined=True, title="Der arithmetische Peano-Faden",
                      axes_labels=['Nullstelle n', 'Phase mod 60'], color='red')
        return p1, p2, facette_counts

    # Fallback: matplotlib
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    ax1.bar(range(G60), facette_counts, color=(0, 0, 0.55))
    ax1.set_title(f"Besetzung der 60 Facetten (N={N_max})")
    ax1.set_xlabel("Facette Index")
    ax1.set_ylabel("Anzahl Resonanzen")
    xs, ys = zip(*resonances)
    ax2.plot(xs, ys, 'r-')
    ax2.set_title("Der arithmetische Peano-Faden")
    ax2.set_xlabel("Nullstelle n")
    ax2.set_ylabel("Phase mod 60")
    plt.tight_layout()
    plt.show()
    return None, None, facette_counts

# Ausführung
if __name__ == "__main__":
    plot_facettes, plot_peano, data = energiedoku_simulation(3333)
    if USE_SAGE:
        show(plot_facettes)
        show(plot_peano)