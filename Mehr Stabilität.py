# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Spektrale Analyse der Primzahlen
#
# Ziel dieser Analyse ist es, die Verteilung der Primzahlen im Kontext der Riemannschen Zeta-Funktion zu untersuchen.

# %%
import math
import sys
from pathlib import Path
import numpy as np
import matplotlib
if 'ipykernel' not in sys.modules:
    matplotlib.use('Agg')  # Nicht-interaktiv, wenn als Skript ausgeführt
import matplotlib.pyplot as plt
from numpy.fft import fft as np_fft

# %%
def spektrale_analyse_primzahlen(n_max=1000):
    """
    Spektrale Analyse der Riemann-Nullstellen im Kontext des 8. Hilbert-Problems.
    Berechnet Abstände zwischen Nullstellen, Skalierungs-Resonanz und erstellt
    Visualisierungen (riemann_rigidity.png, Power Spectrum).
    """
    # %% [markdown]
    # ## 1. Konstanten definieren

    c_stern = math.sqrt(math.pi * math.e / 2)
    
    # %% [markdown]
    # ## 2. Riemann-Nullstellen (Imaginärteile der ersten Nullstellen)

    try:
        script_dir = Path(__file__).resolve().parent
    except NameError:
        script_dir = Path('.')
    zeros_arr = None
    for name in ('zeros6.npz', 'zeros6.npy'):
        p = script_dir / name
        if p.exists():
            try:
                data = np.load(p)
                zeros_arr = data['zeros'] if name.endswith('.npz') else data
                break
            except Exception:
                pass
    if zeros_arr is not None:
        zeros = zeros_arr[:n_max].tolist()
    else:
        zeros = [14.1347, 21.0220, 25.0109, 30.4249, 32.9351, 37.5862, 40.9187,
                 43.3271, 48.0052, 49.7738, 52.9703, 56.4462, 59.3470, 60.8318,
                 65.1125, 67.0798, 69.5464, 72.0672, 75.7047, 77.1448, 79.3374,
                 82.9104, 84.7355, 87.4253, 88.8091, 92.4919, 94.6513, 95.8706,
                 98.8312, 101.3179, 103.7255, 105.4466, 107.1686, 111.0295,
                 111.8747, 114.3202, 116.2270, 118.7908, 121.3701, 122.9468,
                 124.2568, 127.5167, 129.5787, 131.0887, 133.4977, 134.7570,
                 138.1160, 139.7362, 141.1237, 143.1118][:n_max]
    
    # %% [markdown]
    # ## 3. Normalisierung der Abstände (Spektrale Rigidität)
    diffs = [zeros[i+1] - zeros[i] for i in range(len(zeros)-1)]
    avg_diff = sum(diffs) / len(diffs)
    
    # %% [markdown]
    # ## 4. Abgleich mit dem Modell
    skalierte_dichte = avg_diff / c_stern
    
    print(f"--- Spektrale Analyse des 8. Hilbert-Problems ---")
    print(f"C_* Skala:                {c_stern}")
    print(f"Mittlerer Nullstellenabstand: {avg_diff}")
    print(f"Skalierungs-Resonanz:      {skalierte_dichte}")
    
    # %% [markdown]
    # ## 5. Visualisierung der Korrelation
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
    x_vals = list(range(len(diffs)))
    ax1.step(x_vals, diffs, where='mid', color='blue', label='Nullstellen-Abstände')
    ax1.axhline(avg_diff, color='red', linestyle='--', label='Mittelwert')
    ax1.set_xlabel('Index')
    ax1.set_ylabel('Abstand')
    ax1.legend()
    ax1.set_title('Riemann-Nullstellen-Abstände')
    
    # %% [markdown]
    # ## 5a. Power Spectrum der Abstände
    fft_data = np_fft(diffs)
    power_spectrum = [abs(x)**2 for x in fft_data]
    ax2.plot(power_spectrum, color='green', label='Power Spectrum')
    ax2.set_xlabel('Frequenz')
    ax2.set_ylabel('Power')
    ax2.legend()

    # %% [markdown]
    # ## 6. Speichern der Visualisierung
    fig.tight_layout()
    fig.savefig('riemann_rigidity.png')
    plt.close(fig)
    print("\nGrafik 'riemann_rigidity.png' wurde erstellt.")
    return power_spectrum

# %%
power_spectrum = spektrale_analyse_primzahlen()

# %% [markdown]
# Dies ist das Ende der spektralen Analyse. Die Grafik 'riemann_rigidity.png' visualisiert die Abstände zwischen den Nullstellen und ihren Mittelwert.

# %%
if __name__ == "__main__":
    pass  # Grafik bereits in spektrale_analyse_primzahlen() gespeichert