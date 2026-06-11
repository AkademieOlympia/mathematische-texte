from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isqrt
from statistics import mean
from typing import List, Tuple, Optional, Set, Dict


# =========================================================
# Primzahlen, Zwillinge, Resonanz
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


def is_prime(n: int, prime_flags: List[bool]) -> bool:
    return 0 <= n < len(prime_flags) and prime_flags[n]


def is_twin_center(m: int, prime_flags: List[bool]) -> bool:
    return m - 1 >= 2 and m + 1 < len(prime_flags) and prime_flags[m - 1] and prime_flags[m + 1]


def stable_twin_centers(bound: int, prime_flags: List[bool]) -> Set[int]:
    return {m for m in range(2, bound) if is_twin_center(m, prime_flags)}


def resonance_classes(T: Set[int]) -> Tuple[Set[int], Set[int], Set[int], Set[int]]:
    RL = {m for m in T if (m + 6) in T}
    RR = {m for m in T if (m - 6) in T}
    R = RL | RR
    Rbi = RL & RR
    return R, RL, RR, Rbi


# =========================================================
# Vierlingszentren
# =========================================================

def standard_fourling_centers(limit: int, prime_flags: List[bool]) -> List[int]:
    centers = []
    for p in range(2, limit - 8 + 1):
        if prime_flags[p] and prime_flags[p + 2] and prime_flags[p + 6] and prime_flags[p + 8]:
            centers.append(p + 4)
    return centers


def generate_states_from_fourling_centers(
    four_centers: List[int],
    max_states: Optional[int] = None
) -> List[Tuple[int, int, int, int]]:
    states = []
    for c in four_centers:
        states.append((c - 4, c - 2, c + 2, c + 4))
        if max_states is not None and len(states) >= max_states:
            break
    return states


def generate_shifted_states_from_fourling_centers(
    four_centers: List[int],
    shifts: List[int]
) -> List[Tuple[int, int, int, int]]:
    states = []
    for c in four_centers:
        for s in shifts:
            states.append((c - 4 + s, c - 2 + s, c + 2 + s, c + 4 + s))
    return states


# =========================================================
# Bertrand-Anker
# =========================================================

def all_primes_in_open_interval(lo: int, hi: int, prime_flags: List[bool]) -> List[int]:
    return [p for p in range(lo + 1, hi) if is_prime(p, prime_flags)]


def choose_bertrand_prime(
    lo: int,
    hi: int,
    prime_flags: List[bool],
    mode: str = "first"
) -> Optional[int]:
    ps = all_primes_in_open_interval(lo, hi, prime_flags)
    if not ps:
        return None

    if mode == "first":
        return ps[0]
    elif mode == "last":
        return ps[-1]
    elif mode == "middle":
        target = (lo + hi) / 2
        return min(ps, key=lambda p: abs(p - target))
    else:
        raise ValueError(f"Unbekannter Bertrand-Modus: {mode}")


# =========================================================
# Triaden-Relaxation
# =========================================================

@dataclass
class TriadCheck:
    success: bool
    kind: Optional[str]          # prime / twin_center / resonance / none
    witness: Optional[int]
    residues: Tuple[int, int, int]


def triad_relaxation_check(
    triad: Tuple[int, int, int],
    pB: int,
    prime_flags: List[bool],
    T: Set[int],
    R: Set[int],
    criterion: str = "prime_or_twin_or_resonance"
) -> TriadCheck:
    residues = tuple(abs(x - pB) for x in triad)

    def classify(r: int) -> Optional[str]:
        if criterion == "prime_only":
            return "prime" if is_prime(r, prime_flags) else None

        if criterion == "twin_only":
            return "twin_center" if r in T else None

        if criterion == "resonance_only":
            return "resonance" if r in R else None

        if criterion == "prime_or_twin":
            if is_prime(r, prime_flags):
                return "prime"
            if r in T:
                return "twin_center"
            return None

        if criterion == "prime_or_twin_or_resonance":
            if is_prime(r, prime_flags):
                return "prime"
            if r in T:
                return "twin_center"
            if r in R:
                return "resonance"
            return None

        raise ValueError(f"Unbekanntes Kriterium: {criterion}")

    for r in residues:
        kind = classify(r)
        if kind is not None:
            return TriadCheck(True, kind, r, residues)

    return TriadCheck(False, None, None, residues)


# =========================================================
# Hierarchischer Niveautest
# =========================================================

@dataclass
class LevelHit:
    state: Tuple[int, int, int, int]
    sorted_state: Tuple[int, int, int, int]
    level: Optional[int]
    absorber_value: Optional[int]
    radius: Optional[int]
    bertrand_prime: Optional[int]
    triad: Optional[Tuple[int, int, int]]
    success_kind: Optional[str]
    witness: Optional[int]
    residues: Optional[Tuple[int, int, int]]


def hierarchical_absorption_test(
    states: List[Tuple[int, int, int, int]],
    prime_flags: List[bool],
    T: Set[int],
    R: Set[int],
    bertrand_mode: str = "first",
    criterion: str = "prime_or_twin_or_resonance"
) -> List[LevelHit]:
    hits: List[LevelHit] = []

    for state in states:
        ms = tuple(sorted(state))  # m1 <= m2 <= m3 <= m4
        chosen: Optional[LevelHit] = None

        for idx in range(4):
            mk = ms[idx]
            triad = tuple(ms[j] for j in range(4) if j != idx)

            pB = choose_bertrand_prime(mk, 2 * mk, prime_flags, mode=bertrand_mode)
            if pB is None:
                continue

            chk = triad_relaxation_check(
                triad, pB, prime_flags, T, R, criterion=criterion
            )

            if chk.success:
                chosen = LevelHit(
                    state=state,
                    sorted_state=ms,
                    level=idx + 1,
                    absorber_value=mk,
                    radius=2 * mk,
                    bertrand_prime=pB,
                    triad=triad,
                    success_kind=chk.kind,
                    witness=chk.witness,
                    residues=chk.residues,
                )
                break

        if chosen is None:
            chosen = LevelHit(
                state=state,
                sorted_state=ms,
                level=None,
                absorber_value=None,
                radius=None,
                bertrand_prime=None,
                triad=None,
                success_kind=None,
                witness=None,
                residues=None,
            )

        hits.append(chosen)

    return hits


# =========================================================
# Auswertung
# =========================================================

def summarize_hits(hits: List[LevelHit]) -> None:
    level_counter = Counter(h.level for h in hits)
    kind_counter = Counter(h.success_kind for h in hits if h.success_kind is not None)
    witness_counter = Counter(h.witness for h in hits if h.witness is not None)

    n = len(hits)
    success = sum(1 for h in hits if h.level is not None)

    print(f"Anzahl Zustände               : {n}")
    print(f"Erfolgreiche Niveautreffer    : {success}")
    print(f"Erfolgsquote                  : {success / n:.6f}" if n else "Erfolgsquote: n/a")
    print()

    print("Niveaustatistik:")
    for k in [1, 2, 3, 4, None]:
        label = f"Level {k}" if k is not None else "kein Treffer"
        print(f"  {label:>12}: {level_counter.get(k, 0)}")
    print()

    print("Erfolgsarten:")
    for kind in ["prime", "twin_center", "resonance"]:
        print(f"  {kind:>12}: {kind_counter.get(kind, 0)}")
    print()

    print("Häufigste Witness-Werte:")
    for val, cnt in witness_counter.most_common(20):
        print(f"  {str(val):>8}: {cnt}")
    print()

    print("Erste erfolgreiche Treffer:")
    shown = 0
    for h in hits:
        if h.level is not None:
            print(
                f"  Zustand={h.state} | sortiert={h.sorted_state} | "
                f"Level={h.level} | m_k={h.absorber_value} | R={h.radius} | "
                f"p_B={h.bertrand_prime} | Triade={h.triad} | "
                f"Residuen={h.residues} | Art={h.success_kind} | Witness={h.witness}"
            )
            shown += 1
            if shown >= 20:
                break


# =========================================================
# Hauptprogramm
# =========================================================

def main() -> None:
    limit = 1_000_000
    prime_flags = prime_sieve(limit + 100)

    T = stable_twin_centers(limit, prime_flags)
    R, RL, RR, Rbi = resonance_classes(T)

    four_centers = standard_fourling_centers(limit, prime_flags)

    print(f"Schranke: {limit}")
    print(f"Anzahl Vierlingszentren: {len(four_centers)}")
    print(f"|T|={len(T)}, |R|={len(R)}")
    print()

    states_A = generate_states_from_fourling_centers(four_centers[:200])
    states_B = generate_shifted_states_from_fourling_centers(four_centers[:100], shifts=[-3, -1, 1, 3])

    print("=== Test A: reine Vierlingszustände ===")
    hits_A = hierarchical_absorption_test(
        states_A,
        prime_flags,
        T,
        R,
        bertrand_mode="first",
        criterion="prime_or_twin_or_resonance"
    )
    summarize_hits(hits_A)

    print("\n=== Test B: verschobene Zustände ===")
    hits_B = hierarchical_absorption_test(
        states_B,
        prime_flags,
        T,
        R,
        bertrand_mode="first",
        criterion="prime_or_twin_or_resonance"
    )
    summarize_hits(hits_B)


if __name__ == "__main__":
    main()