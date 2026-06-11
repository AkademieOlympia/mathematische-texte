# Prinzip-Skizze für die Abweichungsanalyse im Modell
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _prime_pi_values(bound):
    """Gibt pi(x) fuer alle x < bound zurueck."""
    if bound <= 2:
        return [0] * max(bound, 0)

    is_prime = [True] * bound
    is_prime[0] = False
    is_prime[1] = False

    for candidate in range(2, int(math.sqrt(bound - 1)) + 1):
        if is_prime[candidate]:
            start = candidate * candidate
            is_prime[start:bound:candidate] = [False] * len(is_prime[start:bound:candidate])

    values = [0] * bound
    count = 0
    for x in range(bound):
        if is_prime[x]:
            count += 1
        values[x] = count

    return values


def _li_values(bound):
    """Approximation der offset-logarithmischen Integralfunktion Li(x)."""
    values = [0.0] * max(bound, 0)
    if bound <= 2:
        return values

    total = 0.0
    previous = 1.0 / math.log(2)
    for x in range(3, bound):
        current = 1.0 / math.log(x)
        total += 0.5 * (previous + current)
        values[x] = total
        previous = current

    return values

def prime_deviation_energy(bound):
    """
    Berechnet die energetische Fluktuation (pi(x) - Li(x)) bis zu einer höheren Schranke.
    """
    prime_counts = _prime_pi_values(bound)
    li_approximations = _li_values(bound)

    results = []
    for x in range(2, bound):
        pi_x = prime_counts[x]
        li_x = li_approximations[x]
        deviation = pi_x - li_x
        # In unserem Modell korreliert dies mit der Spur der quaternionischen Übergangsmatrix
        results.append((x, deviation))
    return results


def print_summary(bound, daten, prime_counts, li_approximations):
    """Gibt die letzten Werte zur Kontrolle im Terminal aus."""
    print(f"--- Analyse bis x = {bound} ---")
    print("x     | pi(x)  | Li(x)     | Abweichung (Energie-Fluktuation)")
    print("-" * 58)

    for x, dev in daten[-5:]:
        pi_x = prime_counts[x]
        li_x = li_approximations[x]
        print(f"{x:<5} | {pi_x:<6} | {li_x:<9.2f} | {dev:.4f}")


def save_plot(daten, bound, filename):
    """Speichert die arithmetische Interferenz als PNG."""
    x_values = [x for x, _ in daten]
    deviations = [dev for _, dev in daten]

    riemann_obere = [
        (1 / (8 * math.pi)) * math.sqrt(x) * math.log(x)
        for x in x_values
    ]
    riemann_untere = [-value for value in riemann_obere]

    plt.figure(figsize=(12, 7))
    plt.plot(x_values, deviations, color="blue", label="pi(x) - Li(x)")
    plt.plot(
        x_values,
        riemann_obere,
        color="red",
        linestyle="--",
        label="Rigiditäts-Schranke (von Koch)",
    )
    plt.plot(x_values, riemann_untere, color="red", linestyle="--")
    plt.title(f"Arithmetische Interferenz und Rigidität bis x = {bound}")
    plt.xlabel("x")
    plt.ylabel("Energie-Fluktuation")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=160)
    plt.close()


if __name__ == "__main__":
    schranke = 100000000
    dateiname = "bamberger_modell_10k.png"

    daten = prime_deviation_energy(schranke)
    prime_counts = _prime_pi_values(schranke)
    li_approximations = _li_values(schranke)

    print_summary(schranke, daten, prime_counts, li_approximations)
    save_plot(daten, schranke, dateiname)
    print(f"\nErweiterte Grafik wurde als '{dateiname}' im aktuellen Verzeichnis gespeichert.")