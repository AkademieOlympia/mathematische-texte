"""
Prüft, ob der Übergang von Quadrupel Q1 zu Q2 energetisch möglich ist.
Nach Barandes' HMM-Modell (2026). Ohne Sage – reines Python.
Calibriert gegen echte Primzahl-Vierlinge (Quadranten mod 12).
"""

import os

def is_prime(n):
    """Einfacher Primzahltest."""
    if n < 2:
        return False
    for d in range(2, int(n**0.5) + 1):
        if n % d == 0:
            return False
    return True

def quat_mul(q1, q2):
    """Quaternionen-Multiplikation."""
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    return (
        a1*a2 - b1*b2 - c1*c2 - d1*d2,
        a1*b2 + b1*a2 + c1*d2 - d1*c2,
        a1*c2 - b1*d2 + c1*a2 + d1*b2,
        a1*d2 + b1*c2 - c1*b2 + d1*a2,
    )

def quat_conj(q):
    """Konjugiertes Quaternion."""
    return (q[0], -q[1], -q[2], -q[3])

def quat_norm_sq(q):
    """|q|^2."""
    return sum(x*x for x in q)

def quat_inverse(q):
    """q^{-1} = conj(q) / |q|^2."""
    n2 = quat_norm_sq(q)
    c = quat_conj(q)
    return tuple(x / n2 for x in c)

def energy_factor(Q1, Q2):
    """Energie-Faktor für Übergang Q1 -> Q2 (ohne Ausgabe)."""
    q1, q2 = tuple(Q1), tuple(Q2)
    q1_inv = quat_inverse(q1)
    delta = quat_mul(q1_inv, q2)
    return quat_norm_sq(delta) ** 0.5

def check_quadruple_transition(Q1, Q2, verbose=True):
    """
    Prüft, ob der Übergang von Quadrupel Q1 zu Q2 energetisch möglich ist.
    """
    q1 = tuple(Q1)
    q2 = tuple(Q2)

    required_energy = energy_factor(q1, q2)
    n1 = quat_norm_sq(q1) ** 0.5
    n2 = quat_norm_sq(q2) ** 0.5

    if verbose:
        print(f"Start-Quadrupel (E={n1:.2f}): {Q1}")
        print(f"Ziel-Quadrupel  (E={n2:.2f}): {Q2}")
        print(f"Benötigter Energie-Faktor: {required_energy:.4f}")

    is_resonant = is_prime(int(round(required_energy)))
    return is_resonant

def quadrant(p):
    """Quadrant: p mod 12 (Familien 1, 5, 7, 11)."""
    return p % 12

def calibrate_against_quadruplets(csv_path="prime_quadruplets.csv", max_pairs=500):
    """
    Kalibriert das Modell gegen echte Primzahl-Vierlinge.
    Gruppiert nach Quadrant (p mod 12) des Start-Vierlings.
    """
    if not os.path.exists(csv_path):
        print(f"Datei nicht gefunden: {csv_path}")
        print("Führe zuerst 'Erzeuge Vierlinge.py' aus.")
        return

    quads = []
    with open(csv_path, "r", encoding="utf-8") as f:
        next(f)  # header
        for line in f:
            p, p2, p6, p8 = map(int, line.strip().split(","))
            quads.append((p, p2, p6, p8))

    print(f"\n--- Kalibrierung gegen {len(quads)} Primzahl-Vierlinge ---")
    print(f"Quadranten: p mod 12 ∈ {{1, 5, 7, 11}}\n")

    n_resonant = 0
    n_total = 0
    by_quadrant = {}  # quadrant -> (resonant, total)

    for i in range(min(len(quads) - 1, max_pairs)):
        Q1, Q2 = quads[i], quads[i + 1]
        e = energy_factor(Q1, Q2)
        r = int(round(e))
        resonant = is_prime(r) if r >= 2 else False

        n_resonant += int(resonant)
        n_total += 1

        q_start = quadrant(Q1[0])
        if q_start not in by_quadrant:
            by_quadrant[q_start] = [0, 0]
        by_quadrant[q_start][1] += 1
        if resonant:
            by_quadrant[q_start][0] += 1

    print("Gesamt:")
    print(f"  Übergänge: {n_total}")
    print(f"  Resonant (round(E) prim): {n_resonant} ({100*n_resonant/n_total:.1f}%)")
    print(f"  Zufallserwartung: ~{100*0.38:.0f}% (Dichte kleiner Primzahlen)")

    print("\nNach Quadrant (Start-Vierling p mod 12):")
    for q in sorted(by_quadrant.keys()):
        r, t = by_quadrant[q]
        fam = {1: "E", 5: "A", 7: "B", 11: "C"}.get(q, str(q))
        print(f"  Quadrant {q} ({fam}): {r}/{t} resonant ({100*r/t:.1f}%)")

# Beispiel: Kann [2, 3, 5, 7] zu [11, 13, 17, 19] springen?
if __name__ == "__main__":
    Q_start = [2, 3, 5, 7]
    Q_ziel = [11, 13, 17, 19]

    resonant = check_quadruple_transition(Q_start, Q_ziel)
    print(f"Übergang resonant: {resonant}")

    calibrate_against_quadruplets(csv_path="prime_quadruplets.csv", max_pairs=500)
