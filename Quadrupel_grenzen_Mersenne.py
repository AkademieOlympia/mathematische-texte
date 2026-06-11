#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zwei begrenzende Quadrupel für die größte bekannte Primzahl M = 2^82589933 - 1.

Definition (Familien E, A, B, C nach Schütte/Gitter):
- E: Primzahlen ≡ 1 (mod 12)
- A: Primzahlen ≡ 5 (mod 12)
- B: Primzahlen ≡ 7 (mod 12)
- C: Primzahlen ≡ 11 (mod 12)

Quadrupel = Produkt von genau einer Primzahl aus jeder Familie: P = p_E * p_A * p_B * p_C.

Die größte bekannte Primzahl M = 2^82589933 - 1 wird durch zwei Quadrupel eingeschlossen:
- Unteres Quadrupel: N_low = p_E * p_A * p_B * p_C maximal mit N_low ≤ M
- Oberes Quadrupel: N_high = p_E' * p_A' * p_B' * p_C' minimal mit N_high ≥ M

Dann: N_low ≤ M ≤ N_high  (mit „relativ wenigen Schritten“, da nur je vier Primzahlen).

Für M = 2^82589933 - 1 können die konkreten Primzahlen nicht angegeben werden (Größenordnung ~10^6 Ziffern pro Faktor).
Dieses Skript demonstriert den Prozess an einem handhabbaren M (z.B. 2^50 - 1) und gibt die exakten Definitionen für M.
"""

from math import isqrt

# Familien: Rest mod 12
FAM_E, FAM_A, FAM_B, FAM_C = 1, 5, 7, 11
RESIDUES = (FAM_E, FAM_A, FAM_B, FAM_C)
FAM_NAMES = {1: "E", 5: "A", 7: "B", 11: "C"}


def is_prime_small(n):
    """Einfacher Primzahltest für kleine n."""
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    d = 5
    while d * d <= n:
        if n % d == 0 or n % (d + 2) == 0:
            return False
        d += 6
    return True


def next_prime_in_class(start, residue_mod12, direction=+1):
    """Nächste Primzahl mit p ≡ residue_mod12 (mod 12). direction +1 = aufwärts, -1 = abwärts."""
    p = start
    if direction > 0:
        # Gehe zum nächsten Kongruenten
        r = p % 12
        diff = (residue_mod12 - r) % 12
        if diff == 0 and p > 0:
            p += 12
        else:
            p += diff
        while True:
            if p < 2:
                return None
            if is_prime_small(p):
                return p
            p += 12
    else:
        r = p % 12
        diff = (r - residue_mod12) % 12
        if diff == 0:
            p -= 12
        else:
            p -= diff
        while True:
            if p < 2:
                return None
            if is_prime_small(p):
                return p
            p -= 12


def largest_prime_in_class_up_to(n_bound, residue_mod12):
    """Größte Primzahl p ≡ residue_mod12 (mod 12) mit p ≤ n_bound."""
    p = (n_bound // 12) * 12 + residue_mod12
    if p > n_bound:
        p -= 12
    while p >= 2:
        if is_prime_small(p):
            return p
        p -= 12
    return None


def smallest_prime_in_class_at_least(n_bound, residue_mod12):
    """Kleinste Primzahl p ≡ residue_mod12 (mod 12) mit p ≥ n_bound."""
    r = n_bound % 12
    diff = (residue_mod12 - r) % 12
    p = n_bound + diff
    if p < 2:
        p = max(2, residue_mod12 if residue_mod12 != 1 else 13)
    while True:
        if is_prime_small(p):
            return p
        p += 12


def find_lower_quadrupel(M):
    """
    Größtes Quadrupel N_low = p_E * p_A * p_B * p_C mit N_low ≤ M.
    Gibt (N_low, (p_E, p_A, p_B, p_C)) zurück.
    """
    t = int(M ** 0.25)
    if t < 12:
        t = 12
    # Start: größte Primzahl pro Familie ≤ t
    primes = []
    for r in RESIDUES:
        p = largest_prime_in_class_up_to(t, r)
        if p is None:
            return None, None
        primes.append(p)
    N = primes[0] * primes[1] * primes[2] * primes[3]
    while N > M and t > 2:
        t -= 1
        primes = []
        for r in RESIDUES:
            p = largest_prime_in_class_up_to(t, r)
            if p is None:
                return None, None
            primes.append(p)
        N = primes[0] * primes[1] * primes[2] * primes[3]
    if N > M:
        return None, None
    # Optimiere: versuche eine Komponente zu vergrößern (bleibe ≤ M)
    improved = True
    while improved:
        improved = False
        for i in range(4):
            # Ersetze primes[i] durch nächstgrößere Primzahl derselben Familie
            p_next = next_prime_in_class(primes[i] + 1, RESIDUES[i], +1)
            if p_next is None:
                continue
            new_primes = primes[:]
            new_primes[i] = p_next
            new_N = new_primes[0] * new_primes[1] * new_primes[2] * new_primes[3]
            if new_N <= M and new_N > N:
                N = new_N
                primes = new_primes
                improved = True
                break
    return N, tuple(primes)


def find_upper_quadrupel(M):
    """
    Kleinstes Quadrupel N_high = p_E * p_A * p_B * p_C mit N_high ≥ M.
    Gibt (N_high, (p_E, p_A, p_B, p_C)) zurück.
    """
    t = int(M ** 0.25) + 1
    primes = []
    for r in RESIDUES:
        p = smallest_prime_in_class_at_least(t, r)
        if p is None:
            return None, None
        primes.append(p)
    N = primes[0] * primes[1] * primes[2] * primes[3]
    while N < M:
        t += 1
        primes = []
        for r in RESIDUES:
            p = smallest_prime_in_class_at_least(t, r)
            if p is None:
                return None, None
            primes.append(p)
        N = primes[0] * primes[1] * primes[2] * primes[3]
    # Optimiere: versuche eine Komponente zu verkleinern (bleibe ≥ M)
    improved = True
    while improved:
        improved = False
        for i in range(4):
            # Ersetze primes[i] durch nächstkleinere Primzahl derselben Familie (aber noch so dass Produkt ≥ M)
            p_prev = next_prime_in_class(primes[i] - 1, RESIDUES[i], -1)
            if p_prev is None:
                continue
            new_primes = primes[:]
            new_primes[i] = p_prev
            new_N = new_primes[0] * new_primes[1] * new_primes[2] * new_primes[3]
            if new_N >= M and new_N < N:
                N = new_N
                primes = new_primes
                improved = True
                break
    return N, tuple(primes)


def main():
    # Demo mit M = 2^80 - 1 (Achtzig)
    M_demo = (1 << 80) - 1
    print("=== Demo: Zwei begrenzende Quadrupel (M = 2^80 - 1) ===\n")
    print("Zielzahl M = 2^80 - 1 =", M_demo)
    print("(Analog für die größte bekannte Primzahl M = 2^82589933 - 1.)\n")

    N_low, quad_low = find_lower_quadrupel(M_demo)
    N_high, quad_high = find_upper_quadrupel(M_demo)

    if N_low is None or quad_low is None:
        print("Unteres Quadrupel konnte nicht bestimmt werden.")
    else:
        print("--- Unteres Quadrupel (größtes Quadrupel ≤ M) ---")
        print("  p_E (≡1 mod 12) =", quad_low[0])
        print("  p_A (≡5 mod 12) =", quad_low[1])
        print("  p_B (≡7 mod 12) =", quad_low[2])
        print("  p_C (≡11 mod 12) =", quad_low[3])
        print("  N_low = p_E * p_A * p_B * p_C =", N_low)
        print("  Prüfung: N_low ≤ M?", N_low <= M_demo)

    print()
    if N_high is None or quad_high is None:
        print("Oberes Quadrupel konnte nicht bestimmt werden.")
    else:
        print("--- Oberes Quadrupel (kleinstes Quadrupel ≥ M) ---")
        print("  p_E' (≡1 mod 12) =", quad_high[0])
        print("  p_A' (≡5 mod 12) =", quad_high[1])
        print("  p_B' (≡7 mod 12) =", quad_high[2])
        print("  p_C' (≡11 mod 12) =", quad_high[3])
        print("  N_high = p_E' * p_A' * p_B' * p_C' =", N_high)
        print("  Prüfung: N_high ≥ M?", N_high >= M_demo)

    if N_low is not None and N_high is not None:
        print()
        print("=== Einschließung ===")
        print("  N_low ≤ M ≤ N_high:", N_low, "≤", M_demo, "≤", N_high)

    print()
    print("=== Exakte Definition für die größte bekannte Primzahl ===")
    print("  M = 2^82589933 - 1  (größte bekannte Primzahl, Mersenne-Primzahl)")
    print()
    print("  Unteres Quadrupel N_low:")
    print("    N_low = p_E * p_A * p_B * p_C")
    print("    mit p_E ≡ 1, p_A ≡ 5, p_B ≡ 7, p_C ≡ 11 (mod 12), alle prim,")
    print("    und N_low ist das Maximum aller solchen Produkte mit N_low ≤ M.")
    print()
    print("  Oberes Quadrupel N_high:")
    print("    N_high = p_E' * p_A' * p_B' * p_C'")
    print("    mit p_E' ≡ 1, p_A' ≡ 5, p_B' ≡ 7, p_C' ≡ 11 (mod 12), alle prim,")
    print("    und N_high ist das Minimum aller solchen Produkte mit N_high ≥ M.")
    print()
    print("  Dann gilt:  N_low ≤ M ≤ N_high .")
    print("  Die konkreten Primzahlen für M = 2^82589933 - 1 haben je etwa 6,2 Mio. Dezimalstellen")
    print("  und können nicht explizit angegeben werden; die beiden Quadrupel sind durch")
    print("  obige Maximal-/Minimal-Eigenschaft eindeutig bestimmt.")


if __name__ == "__main__":
    main()
