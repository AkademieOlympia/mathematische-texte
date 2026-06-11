"""
Das Primzahlkreuz (24er-System) – #Energiedoku.
Ohne Sage – NumPy + Matplotlib.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def prime_sieve(limit: int) -> np.ndarray:
    """Sieb des Eratosthenes – schnell für viele Zahlen."""
    is_p = np.ones(limit + 1, dtype=bool)
    is_p[:2] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_p[i]:
            is_p[i * i : limit + 1 : i] = False
    return is_p


def prime_squares(limit: int, is_prime_arr: np.ndarray) -> set:
    """Quadrate von Primzahlen bis sqrt(limit)."""
    squares = set()
    for p in range(2, min(int(limit**0.5) + 1, len(is_prime_arr))):
        if is_prime_arr[p]:
            squares.add(p * p)
    return squares


def plot_primzahlkreuz(schalen: int = 10):
    """
    Erzeugt eine grafische Darstellung des Primzahlkreuzes (24er-System).
    """
    p_color = "red"
    q_color = "blue"
    bg_color = "lightgray"

    fig, ax = plt.subplots(figsize=(20, 20))  # Verzehnfachte Auflösung
    ax.set_aspect("equal")

    # 1. Acht Strahlen (koprim zu 24)
    strahlen = [1, 5, 7, 11, 13, 17, 19, 23]
    for s in strahlen:
        angle = np.radians(90 - s * 15)
        ax.plot([0, schalen * np.cos(angle)], [0, schalen * np.sin(angle)],
                color=bg_color, linestyle="--", zorder=1)

    # 2. Quadrate von Primzahlen (für Markierung)
    max_n = 24 * schalen
    is_prime_arr = prime_sieve(max_n)
    squares = prime_squares(max_n, is_prime_arr)

    # 3. Zahlen auf den Schalen verteilen
    primes_x, primes_y, primes_labels = [], [], []
    quads_x, quads_y, quads_labels = [], [], []
    other_x, other_y = [], []

    for n in range(1, max_n + 1):
        schale = int(np.ceil(n / 24))
        pos_in_24 = n % 24
        if pos_in_24 == 0:
            pos_in_24 = 24

        angle = np.radians(90 - pos_in_24 * 15)
        x = schale * np.cos(angle)
        y = schale * np.sin(angle)

        if n < len(is_prime_arr) and is_prime_arr[n]:
            primes_x.append(x)
            primes_y.append(y)
            primes_labels.append(str(n))
        elif n in squares:
            quads_x.append(x)
            quads_y.append(y)
            quads_labels.append(str(n))
        else:
            other_x.append(x)
            other_y.append(y)

    ax.scatter(other_x, other_y, c="gray", s=10, alpha=0.3, zorder=2)
    ax.scatter(primes_x, primes_y, c=p_color, s=30, zorder=3)
    ax.scatter(quads_x, quads_y, c=q_color, s=80, marker="s", zorder=4)

    for x, y, lbl in zip(primes_x, primes_y, primes_labels):
        ax.text(x * 1.1, y * 1.1, lbl, fontsize=8, color=p_color, ha="center", va="center")
    for x, y, lbl in zip(quads_x, quads_y, quads_labels):
        ax.text(x * 1.1, y * 1.1, lbl, fontsize=9, color=q_color, fontweight="bold", ha="center", va="center")

    ax.set_xlim(-schalen - 2, schalen + 2)
    ax.set_ylim(-schalen - 2, schalen + 2)
    ax.axis("off")
    ax.set_title("Das Primzahlkreuz (24er-System) – #Energiedoku")

    return fig, ax


if __name__ == "__main__":
    fig, ax = plot_primzahlkreuz(80)  # Verzehnfacht (8 → 80 Schalen)
    plt.savefig("Primzahlkreuz.png", dpi=150)
    print("Plot gespeichert: Primzahlkreuz.png")
