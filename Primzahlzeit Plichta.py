"""
Chronobiologische Resonanzanalyse (24h, Monate, Jahr) – Plichta-Primzahlstrahlen.
Ohne Sage – NumPy + Matplotlib.
"""

import numpy as np
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

MONATE = ["Jan", "Feb", "Mär", "Apr", "Mai", "Jun", "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"]


def analyze_time_resonance(jahr: int | None = None):
    """
    Simuliert die Resonanzpunkte der 8 Primzahlstrahlen
    in 24h-Zyklus, 12 Monaten und Jahr.
    """
    p_stunden = [1, 5, 7, 11, 13, 17, 19, 23]
    p_monate = sorted(set(h % 12 or 12 for h in p_stunden))  # 1, 5, 7, 11

    if jahr is None:
        jahr = datetime.now().year

    print("### Chronobiologische Resonanzanalyse (24h, Monate, Jahr) ###\n")
    print(f"{'Stunde':>8} | {'Winkel':>8} | {'Bedeutung im Kreuz'}")
    print("-" * 45)
    for h in p_stunden:
        grad = h * 15
        print(f"{h:8} | {grad:8}° | Primzahlstrahl (Resonanz)")

    print(f"\n{'Monat':>8} | {'Name':>8} | {'Resonanz'}")
    print("-" * 35)
    for m in p_monate:
        print(f"{m:8} | {MONATE[m-1]:>8} | Primzahl-Monat")
    print(f"\nJahr: {jahr}")

    # Visualisierung: 3 Subplots
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharey=False)

    # 1. 24h-Zyklus
    t = np.linspace(0, 24, 500)
    welle = np.sin(2 * np.pi * t / 24)
    axes[0].plot(t, welle, "b-", lw=2, label="24h-Energiewelle")
    for h in p_stunden:
        y_h = np.sin(2 * np.pi * h / 24)
        axes[0].scatter([h], [y_h], c="red", s=80, zorder=5)
    axes[0].axhline(0, color="gray", linestyle="--", alpha=0.5)
    axes[0].set_xlim(0, 24)
    axes[0].set_ylim(-1.2, 1.2)
    axes[0].set_xlabel("Stunde (24h-Zyklus)")
    axes[0].set_ylabel("Amplitude")
    axes[0].set_title("24h-Energiewelle und Primzahl-Resonanz (Plichta)")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # 2. 12 Monate
    t_m = np.linspace(1, 12, 500)
    welle_m = np.sin(2 * np.pi * (t_m - 1) / 12)
    axes[1].plot(t_m, welle_m, "g-", lw=2, label="Jahres-Energiewelle")
    for m in p_monate:
        y_m = np.sin(2 * np.pi * (m - 1) / 12)
        axes[1].scatter([m], [y_m], c="red", s=80, zorder=5)
    axes[1].axhline(0, color="gray", linestyle="--", alpha=0.5)
    axes[1].set_xlim(0.5, 12.5)
    axes[1].set_ylim(-1.2, 1.2)
    axes[1].set_xticks(range(1, 13))
    axes[1].set_xticklabels(MONATE)
    axes[1].set_xlabel("Monat")
    axes[1].set_ylabel("Amplitude")
    axes[1].set_title("12 Monate – Primzahl-Resonanz (Jan, Mai, Jul, Nov)")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    # 3. Jahr
    axes[2].set_xlim(0, 2)
    axes[2].set_ylim(0, 1)
    axes[2].add_patch(Rectangle((0.3, 0.2), 1.4, 0.6, facecolor="steelblue", alpha=0.3, edgecolor="navy"))
    axes[2].text(1, 0.5, str(jahr), fontsize=28, ha="center", va="center", fontweight="bold")
    axes[2].set_title(f"Jahr {jahr}")
    axes[2].axis("off")

    plt.tight_layout()
    return fig, axes


if __name__ == "__main__":
    fig, axes = analyze_time_resonance()
    plt.savefig("Primzahlzeit_Plichta.png", dpi=150)
    print("\nPlot gespeichert: Primzahlzeit_Plichta.png")
