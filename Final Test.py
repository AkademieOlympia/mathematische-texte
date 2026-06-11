import matplotlib.pyplot as plt

def plot_smoking_gun_validation(y_centers, signal, theory, r_coeff):
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Achse 1: Das extrahierte eabc-Signal (Empirie)
    color = 'tab:red'
    ax1.set_xlabel('Projektions-Koordinate (y)')
    ax1.set_ylabel('eabc Polarisations-Signal', color=color)
    ax1.plot(y_centers, signal, color=color, label='eabc Signal (Extrahiert)', linewidth=2, alpha=0.8)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, which='both', linestyle='--', alpha=0.5)

    # Achse 2: Die Riemann-Theorie (Vakuum-Moden)
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('Theoretische Riemann-Mode (Σ amp_k)', color=color)
    ax2.plot(y_centers, theory, color=color, label='Riemann Theorie', linewidth=2, linestyle=':')
    ax2.tick_params(axis='y', labelcolor=color)

    # Titel und Statistik
    plt.title(f"Smoking Gun: Korrelation zwischen eabc-Shift und Riemann-Vakuum\nPearson r = {r_coeff:.4f}")
    
    # Legende zusammenführen
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper left')

    plt.tight_layout()
    plt.show()

# Aufruf nach der Korrelationsberechnung:
# plot_smoking_gun_validation(y_centers, signal, theory, r_coeff)