"""
Schütte-Prognose-Test: Spannungsvorhersage für E-Primzahlen (mod 12).
Behebt Indexfehler durch saubere Zuordnung tension_history[i] <-> history_primes[i].
"""
import numpy as np
import matplotlib.pyplot as plt
try:
    from sympy import primerange
except ImportError:
    def primerange(a, b):
        """Einfacher Eratosthenes-Sieb als Fallback wenn sympy fehlt."""
        if a < 2:
            a = 2
        n = max(b - 1, 2)
        is_prime = np.ones(n + 1, dtype=bool)
        is_prime[0] = is_prime[1] = False
        for i in range(2, int(n**0.5) + 1):
            if is_prime[i]:
                is_prime[i*i:n+1:i] = False
        return (p for p in range(a, n + 1) if is_prime[p])


def analyze_schuette_prediction(limit=100000, lookahead=5000):
    primes = list(primerange(1, limit + lookahead))
    history_primes = [p for p in primes if p <= limit]

    c1, c5, c7, c11 = 0, 0, 0, 0
    # tension_history[j] = Spannung *nach* history_primes[j]; für j=0,1 (2,3) setzen wir 0
    tension_history = []

    for i, p in enumerate(history_primes):
        if p <= 3:
            tension_history.append(0.0)  # Platzhalter, damit Index = i
            continue
        r = p % 12
        if r == 1:
            c1 += 1
        elif r == 5:
            c5 += 1
        elif r == 7:
            c7 += 1
        elif r == 11:
            c11 += 1

        avg_non_res = (c5 + c7 + c11) / 3.0
        tension = avg_non_res - c1
        tension_history.append(tension)

    avg_tension = np.mean(tension_history)
    std_tension = np.std(tension_history)
    print(f"Global Stats (N={limit}):")
    print(f"Mean Tension on E: {avg_tension:.2f}")
    print(f"Max Tension on E:  {max(tension_history):.2f}")

    high_tension_threshold = avg_tension + 1.5 * std_tension
    start_idx = max(2, int(len(history_primes) * 0.9))  # >=2 wegen tension_history[i-1] bei Bedarf

    print(f"\nAnalyzing prediction accuracy in range [{history_primes[start_idx]}, {limit}]...")

    hits = 0
    trials = 0
    for i in range(start_idx, len(history_primes)):
        curr_t = tension_history[i]
        if curr_t > high_tension_threshold:
            trials += 1
            current_p = history_primes[i]
            dist_to_next_e = 0
            found = False
            for j in range(i + 1, len(primes)):
                if primes[j] % 12 == 1:
                    dist_to_next_e = primes[j] - current_p
                    found = True
                    break
            if found:
                avg_gap_e = 4 * np.log(current_p)
                if dist_to_next_e < avg_gap_e:
                    hits += 1

    if trials > 0:
        print(f"Number of High Tension situations: {trials}")
        print(f"Next E-Prime appeared 'quickly' (< avg gap): {hits}")
        print(f"Prediction Accuracy: {hits/trials*100:.1f}%")
    else:
        print("No high-tension events in the analyzed range.")

    return history_primes, tension_history, high_tension_threshold


if __name__ == "__main__":
    h_primes, h_tension, threshold = analyze_schuette_prediction(100000)

    # Plot: letzte 1000 Punkte
    subset = 1000
    n = len(h_primes)
    start = max(0, n - subset)
    x_sub = h_primes[start:]
    t_sub = h_tension[start:]

    plt.figure(figsize=(12, 6))
    plt.plot(x_sub, t_sub, label=r'Spannung auf E ($T_{vac}$)', color='purple')
    plt.axhline(threshold, color='red', linestyle='--', label='Kritische Schwelle (Vorhersage-Trigger)')
    plt.fill_between(x_sub, t_sub, threshold, where=(np.array(t_sub) > threshold),
                    color='red', alpha=0.3, label='Hohe Wahrscheinlichkeit für E')

    plt.title('Schütte-Mikroskop in Aktion: Spannungsvorhersage bei N=100.000')
    plt.xlabel('n')
    plt.ylabel('Vakuum-Spannung auf Klasse E')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('schuette_prediction.png')
    print("Grafik gespeichert: schuette_prediction.png")
