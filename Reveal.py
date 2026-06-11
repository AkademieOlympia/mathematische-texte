import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.ndimage import gaussian_filter

def safe_show():
    try:
        plt.show()
    except KeyboardInterrupt:
        print("\nPlot durch Benutzer abgebrochen.")

def plot_lamb_shift_vector(ax, data_abce, data_ceab, dim_x=0, dim_y=1):
    """
    ax: Das Matplotlib Axis-Objekt deiner 8D-Projektion
    dim_x, dim_y: Die Indizes der projizierten Dimensionen (z.B. e1 und e2)
    """
    # Zentren berechnen
    c_abce = np.mean(data_abce, axis=0)
    c_ceab = np.mean(data_ceab, axis=0)
    
    # Startpunkt (Zentrum ABCE) und Vektor (zu CEAB)
    start_x, start_y = c_abce[dim_x], c_abce[dim_y]
    dx = c_ceab[dim_x] - c_abce[dim_x]
    dy = c_ceab[dim_y] - c_abce[dim_y]
    
    # Den Vektor zeichnen (stark vergrößert, damit er sichtbar ist)
    scale_factor = 10  # Vergrößerung für die Sichtbarkeit des winzigen Shifts
    ax.arrow(start_x, start_y, dx*scale_factor, dy*scale_factor, 
             head_width=0.05, head_length=0.1, fc='red', ec='red', 
             label=f'Lamb-Shift (x{scale_factor})')
    
    # Zentren markieren
    ax.scatter(c_abce[dim_x], c_abce[dim_y], color='yellow', s=100, edgecolors='black', label='Fokus ABCE')
    ax.scatter(c_ceab[dim_x], c_ceab[dim_y], color='blue', s=100, edgecolors='black', label='Fokus CEAB')

def plot_eabc_polarization_map(ax, data_abce, data_ceab, bins=100, sigma=2):
    """
    Erstellt eine Differenz-Heatmap (Vakuumpolarisation).
    Blau = CEAB-Überschuss, Rot = ABCE-Überschuss.
    """
    # 1. Gemeinsamen Bereich festlegen
    x_min = min(data_abce[:, 0].min(), data_ceab[:, 0].min())
    x_max = max(data_abce[:, 0].max(), data_ceab[:, 0].max())
    y_min = min(data_abce[:, 1].min(), data_ceab[:, 1].min())
    y_max = max(data_abce[:, 1].max(), data_ceab[:, 1].max())

    # 2. 2D-Histogramme berechnen
    hist_abce, x_edges, y_edges = np.histogram2d(
        data_abce[:, 0], data_abce[:, 1], bins=bins, 
        range=[[x_min, x_max], [y_min, y_max]]
    )
    hist_ceab, _, _ = np.histogram2d(
        data_ceab[:, 0], data_ceab[:, 1], bins=bins, 
        range=[[x_min, x_max], [y_min, y_max]]
    )

    # 3. Glätten (simuliert die Unschärfe der Wellenfunktion)
    map_abce = gaussian_filter(hist_abce, sigma=sigma)
    map_ceab = gaussian_filter(hist_ceab, sigma=sigma)

    # 4. Differenz bilden (Die Polarisation)
    # Wir normalisieren, damit die Amplituden vergleichbar sind
    diff_map = (map_abce / np.sum(map_abce)) - (map_ceab / np.sum(map_ceab))

    # 5. Plotten
    extent = [x_min, x_max, y_min, y_max]
    im = ax.imshow(diff_map.T, extent=extent, origin='lower', 
                   cmap='RdBu_r', aspect='auto', interpolation='bilinear')
    
    plt.colorbar(im, ax=ax, label='Polarisations-Dichte (ABCE - CEAB)')
    ax.set_title("Arithmetische Vakuumpolarisation (eabc-Feld)")

def overlay_riemann_grid(ax, gamma, y_range, amplitude):
    """
    Legt ein Raster der Riemann-Frequenz über die Heatmap.
    """
    # Wellenlänge der ersten Nullstelle im eabc-Raum (skaliert)
    # Da y in deinen Daten bis -17.5 geht, müssen wir die Phase anpassen
    y_vals = np.linspace(y_range[0], y_range[1], 500)
    
    # Wir simulieren die harmonische Mode: cos(gamma * log(y)) oder einfach linear
    # Je nachdem wie deine 8D-Projektion skaliert ist. 
    # Testen wir zuerst die lineare Resonanz:
    wave = np.cos(gamma * y_vals * 0.5) # Der Faktor 0.5 ist ein Skalierungs-Fit
    
    for i in range(len(y_vals)-1):
        if wave[i] * wave[i+1] < 0: # Nulldurchgang
            ax.axhline(y_vals[i], color='white', linestyle='--', alpha=0.3)

def get_arithmetic_density(n):
    """
    Berechnet die 'Kern-Interaktion' basierend auf Primfaktoren.
    Hohe Werte = starke Interaktion mit dem arithmetischen Vakuum.
    """
    factors = [2, 3, 5, 7, 11] # Die 'inneren Schalen' des Zahlen-Atoms
    density = sum([1.0/p for p in factors if n % p == 0])
    return density

def eabc_refined_lamb_shift(n_max, alpha=0.137):
    n_values = np.arange(2, n_max)
    
    # 1. Dirac-Basis (Primzahlsatz / Li(x) Näherung)
    E_dirac = n_values / np.log(n_values)
    
    # 2. Refinierter S-Orbital-Faktor (Arithmetische Kopplung)
    # Wir berechnen für jede Zahl, wie stark sie im 'Kern' schwingt
    psi_0_sq = np.array([get_arithmetic_density(i) for i in n_values])
    
    # 3. Der Lamb-Shift mit logarithmischer Renormierung
    # Die Verschiebung ist stärker für Zahlen mit hoher Prim-Resonanz
    shift = alpha * (psi_0_sq / np.log(n_values))
    
    # 4. Resultierende Eigenwerte (Das eabc-Spektrum)
    E_lamb = E_dirac + shift
    
    return n_values, E_dirac, E_lamb, shift

# Berechnung
n_range = 200
n, e_d, e_l, s = eabc_refined_lamb_shift(n_range)

# Plotting
plt.figure(figsize=(12, 7))
plt.scatter(n, e_l, c=s, cmap='viridis', s=20, label='eabc-Spektrum (mit Lamb-Shift)')
plt.plot(n, e_d, color='red', linestyle='--', alpha=0.6, label='Dirac-Linie (Ideal)')

# Markierung von Primzahlen (wo der Shift minimal ist)
primes = [i for i in range(2, n_range) if all(i % d != 0 for d in range(2, int(i**0.5)+1))]
plt.scatter(primes, [e_l[i-2] for i in primes], facecolors='none', edgecolors='r', s=80, label='Primzahl-Zustände')

plt.title("Refiniertes eabc-Modell: Aufhebung der Entartung durch Prim-Resonanz")
plt.xlabel("Zahl / Quantenzustand (n)")
plt.ylabel("Eigenwert (Energie)")
plt.legend()
safe_show()

# ------------------------------------------------------------
# Neuer Teil: Lamb-Shift Vektor aus 8D-Projektionen
# ------------------------------------------------------------
abce_file = "vierlinge_riemann_e8_abce_proj8d.npy"
ceab_file = "vierlinge_riemann_e8_ceab_proj8d.npy"

if os.path.exists(abce_file) and os.path.exists(ceab_file):
    print("\nLade 8D-Projektionen für Lamb-Shift-Analyse...")
    abce_proj = np.load(abce_file)
    ceab_proj = np.load(ceab_file)
    
    # Wir nehmen die Dimensionen 0 und 1 (stärkste Varianz)
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Scatter der Punktwolken
    ax.scatter(abce_proj[:,0], abce_proj[:,1], alpha=0.3, c='gold', label='ABCE-Vierlinge', s=15)
    ax.scatter(ceab_proj[:,0], ceab_proj[:,1], alpha=0.3, c='royalblue', label='CEAB-Vierlinge', s=15)
    
    # Lamb-Shift Vektor einzeichnen
    plot_lamb_shift_vector(ax, abce_proj, ceab_proj, dim_x=0, dim_y=1)
    
    ax.set_title("Struktureller Lamb-Shift zwischen ABCE und CEAB (8D-Projektion)")
    ax.set_xlabel("Hauptkomponente 1")
    ax.set_ylabel("Hauptkomponente 2")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    safe_show()

    # Neuer Plot: Vakuumpolarisation
    fig, ax = plt.subplots(figsize=(10, 8))
    plot_eabc_polarization_map(ax, abce_proj, ceab_proj)
    
    # Overlay der ersten Riemann-Nullstelle
    # Gamma1 ~ 14.1347, Skalierung heuristisch
    overlay_riemann_grid(ax, 14.1347, [-17.5, 2.5], 1.17e-2)
    
    safe_show()
else:
    print(f"\n[Info] Projektionsdateien nicht gefunden ({abce_file}, {ceab_file}). Überspringe Vektor-Plot.")