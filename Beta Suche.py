from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isqrt
from typing import Dict, List, Set, Tuple


# =========================================================
# Primzahlen und Zwillingszentren
# =========================================================

def prime_sieve(limit: int) -> List[bool]:
    sieve = [True] * (limit + 1)
    if limit >= 0:
        sieve[0] = False
    if limit >= 1:
        sieve[1] = False

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:limit + 1:p] = [False] * (((limit - start) // p) + 1)

    return sieve


def is_twin_center(m: int, prime_flags: List[bool]) -> bool:
    return m - 1 >= 2 and m + 1 < len(prime_flags) and prime_flags[m - 1] and prime_flags[m + 1]


def stable_twin_centers(bound: int, prime_flags: List[bool]) -> Set[int]:
    return {m for m in range(2, bound) if is_twin_center(m, prime_flags)}


# =========================================================
# Resonanzklassen
# =========================================================

def resonance_classes(T: Set[int]) -> Tuple[Set[int], Set[int], Set[int], Set[int]]:
    RL = {m for m in T if (m + 6) in T}
    RR = {m for m in T if (m - 6) in T}
    R = RL | RR
    Rbi = RL & RR
    return R, RL, RR, Rbi


# =========================================================
# Beta-Kanal
# =========================================================

@dataclass
class BetaStats:
    bound: int
    T_size: int
    R_size: int
    RL_size: int
    RR_size: int
    beta_plus_in_T: int
    beta_minus_in_T: int
    beta_plus_L_to_R: int
    beta_minus_R_to_L: int
    beta_plus_into_R: int
    beta_minus_into_R: int
    beta_active_total: int
    beta_active_fraction_in_T: float
    sample_plus_transitions: List[Tuple[int, int]]
    sample_minus_transitions: List[Tuple[int, int]]
    mod30_beta_active: Dict[int, int]
    mod210_beta_active: Dict[int, int]


def histogram_mod(values: Set[int], modulus: int) -> Dict[int, int]:
    h = Counter(v % modulus for v in values)
    return dict(sorted(h.items()))


def analyze_beta_channel(bound: int = 100_000) -> BetaStats:
    prime_flags = prime_sieve(bound + 20)

    T = stable_twin_centers(bound, prime_flags)
    R, RL, RR, Rbi = resonance_classes(T)

    beta_plus_in_T = 0
    beta_minus_in_T = 0
    beta_plus_L_to_R = 0
    beta_minus_R_to_L = 0
    beta_plus_into_R = 0
    beta_minus_into_R = 0

    plus_transitions: List[Tuple[int, int]] = []
    minus_transitions: List[Tuple[int, int]] = []

    beta_active: Set[int] = set()

    for m in sorted(T):
        mp = m + 6
        mm = m - 6

        if mp in T:
            beta_plus_in_T += 1
            plus_transitions.append((m, mp))
            beta_active.add(m)
            beta_active.add(mp)

            if m in RL and mp in RR:
                beta_plus_L_to_R += 1

            if mp in R:
                beta_plus_into_R += 1

        if mm in T:
            beta_minus_in_T += 1
            minus_transitions.append((m, mm))
            beta_active.add(m)
            beta_active.add(mm)

            if m in RR and mm in RL:
                beta_minus_R_to_L += 1

            if mm in R:
                beta_minus_into_R += 1

    return BetaStats(
        bound=bound,
        T_size=len(T),
        R_size=len(R),
        RL_size=len(RL),
        RR_size=len(RR),
        beta_plus_in_T=beta_plus_in_T,
        beta_minus_in_T=beta_minus_in_T,
        beta_plus_L_to_R=beta_plus_L_to_R,
        beta_minus_R_to_L=beta_minus_R_to_L,
        beta_plus_into_R=beta_plus_into_R,
        beta_minus_into_R=beta_minus_into_R,
        beta_active_total=len(beta_active),
        beta_active_fraction_in_T=(len(beta_active) / len(T) if T else 0.0),
        sample_plus_transitions=plus_transitions[:30],
        sample_minus_transitions=minus_transitions[:30],
        mod30_beta_active=histogram_mod(beta_active, 30),
        mod210_beta_active=histogram_mod(beta_active, 210),
    )


def print_beta_stats(stats: BetaStats) -> None:
    print(f"Schranke: {stats.bound}")
    print()
    print("Grundmengen:")
    print(f"|T|      = {stats.T_size}")
    print(f"|R|      = {stats.R_size}")
    print(f"|R_L|    = {stats.RL_size}")
    print(f"|R_R|    = {stats.RR_size}")
    print()

    print("Beta-Kanal:")
    print(f"B_+(m)=m+6 bleibt in T      : {stats.beta_plus_in_T}")
    print(f"B_-(m)=m-6 bleibt in T      : {stats.beta_minus_in_T}")
    print(f"B_+ : R_L -> R_R            : {stats.beta_plus_L_to_R}")
    print(f"B_- : R_R -> R_L            : {stats.beta_minus_R_to_L}")
    print(f"B_+ landet in R             : {stats.beta_plus_into_R}")
    print(f"B_- landet in R             : {stats.beta_minus_into_R}")
    print()

    print(f"Beta-aktive Zentren insgesamt: {stats.beta_active_total}")
    print(f"Anteil beta-aktiver Zentren in T: {stats.beta_active_fraction_in_T:.6f}")
    print()

    print("Erste B_+-Transitionen:")
    print(stats.sample_plus_transitions)
    print()

    print("Erste B_--Transitionen:")
    print(stats.sample_minus_transitions)
    print()

    print("Modulo-30-Verteilung beta-aktiver Zentren:")
    print(stats.mod30_beta_active)
    print()

    print("Modulo-210-Verteilung beta-aktiver Zentren:")
    nonzero = {k: v for k, v in stats.mod210_beta_active.items() if v > 0}
    print(nonzero)
    print()


if __name__ == "__main__":
    stats = analyze_beta_channel(1_000_000)
    print_beta_stats(stats)