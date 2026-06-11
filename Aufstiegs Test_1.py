from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isqrt
from typing import Dict, List, Set, Tuple


@dataclass
class ResonanceStats:
    bound: int
    T_size: int
    R_size: int
    RL_size: int
    RR_size: int
    Rbi_size: int
    density_R_in_T: float
    mod30_T: Dict[int, int]
    mod30_R: Dict[int, int]
    mod210_R: Dict[int, int]
    sample_T: List[int]
    sample_R: List[int]
    sample_Rbi: List[int]
    fourling_centers: List[int]


def prime_sieve(limit: int) -> List[bool]:
    sieve = [True] * (limit + 1)
    if limit >= 0:
        sieve[0] = False
    if limit >= 1:
        sieve[1] = False
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            step = p
            sieve[start:limit + 1:step] = [False] * (((limit - start) // step) + 1)
    return sieve


def stable_twin_centers(bound: int, prime_flags: List[bool]) -> Set[int]:
    """
    Stabile Zwillingszentren m mit (m-1,m+1) prim und m+1 <= bound.
    """
    T = set()
    for m in range(2, bound):
        if m - 1 >= 2 and m + 1 <= bound and prime_flags[m - 1] and prime_flags[m + 1]:
            T.add(m)
    return T


def resonance_classes(T: Set[int]) -> Tuple[Set[int], Set[int], Set[int], Set[int]]:
    RL = {m for m in T if (m + 6) in T}
    RR = {m for m in T if (m - 6) in T}
    R = RL | RR
    Rbi = RL & RR
    return R, RL, RR, Rbi


def standard_fourling_centers(bound: int, prime_flags: List[bool]) -> List[int]:
    """
    Zentren c = p+4 von Standardvierlingen (p,p+2,p+6,p+8).
    """
    centers = []
    for p in range(2, bound - 8 + 1):
        if prime_flags[p] and prime_flags[p + 2] and prime_flags[p + 6] and prime_flags[p + 8]:
            centers.append(p + 4)
    return centers


def histogram_mod(values: Set[int] | List[int], modulus: int) -> Dict[int, int]:
    hist = Counter(v % modulus for v in values)
    return dict(sorted(hist.items()))


def analyze_resonance(bound: int = 1_000_000) -> ResonanceStats:
    prime_flags = prime_sieve(bound + 10)

    T = stable_twin_centers(bound, prime_flags)
    R, RL, RR, Rbi = resonance_classes(T)
    four_centers = standard_fourling_centers(bound, prime_flags)

    return ResonanceStats(
        bound=bound,
        T_size=len(T),
        R_size=len(R),
        RL_size=len(RL),
        RR_size=len(RR),
        Rbi_size=len(Rbi),
        density_R_in_T=(len(R) / len(T) if T else 0.0),
        mod30_T=histogram_mod(T, 30),
        mod30_R=histogram_mod(R, 30),
        mod210_R=histogram_mod(R, 210),
        sample_T=sorted(T)[:30],
        sample_R=sorted(R)[:30],
        sample_Rbi=sorted(Rbi)[:30],
        fourling_centers=four_centers[:50],
    )


def print_stats(stats: ResonanceStats) -> None:
    print(f"Schranke: {stats.bound}")
    print()
    print("Größen:")
    print(f"|T|      = {stats.T_size}")
    print(f"|R|      = {stats.R_size}")
    print(f"|R_L|    = {stats.RL_size}")
    print(f"|R_R|    = {stats.RR_size}")
    print(f"|R_bi|   = {stats.Rbi_size}")
    print(f"|R|/|T|  = {stats.density_R_in_T:.6f}")
    print()

    print("Modulo-30-Verteilung von T:")
    print(stats.mod30_T)
    print()

    print("Modulo-30-Verteilung von R:")
    print(stats.mod30_R)
    print()

    print("Modulo-210-Verteilung von R (nichtleere Klassen):")
    nonzero_210 = {k: v for k, v in stats.mod210_R.items() if v > 0}
    print(nonzero_210)
    print()

    print("Erste Elemente von T:")
    print(stats.sample_T)
    print()

    print("Erste Elemente von R:")
    print(stats.sample_R)
    print()

    print("Erste Elemente von R_bi:")
    print(stats.sample_Rbi)
    print()

    print("Erste Vierlingszentren c = p+4:")
    print(stats.fourling_centers)
    print()


if __name__ == "__main__":
    stats = analyze_resonance(10**6)
    print_stats(stats)