from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isqrt
from typing import List, Tuple, Optional, Set


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
# Bertrand-Anker
# =========================================================

def all_primes_in_open_interval(lo: int, hi: int, prime_flags: List[bool]) -> List[int]:
    return [p for p in range(lo + 1, hi) if is_prime(p, prime_flags)]


def choose_bertrand_prime(lo: int, hi: int, prime_flags: List[bool], mode: str = "first") -> Optional[int]:
    ps = all_primes_in_open_interval(lo, hi, prime_flags)
    if not ps:
        return None
    if mode == "first":
        return ps[0]
    if mode == "last":
        return ps[-1]
    if mode == "middle":
        target = (lo + hi) / 2
        return min(ps, key=lambda p: abs(p - target))
    raise ValueError(f"Unbekannter Bertrand-Modus: {mode}")


# =========================================================
# Hoch-Beta-Kandidat
# =========================================================

@dataclass
class HighBetaHit:
    triplet: Tuple[int, int, int]
    sorted_triplet: Tuple[int, int, int]
    m: int
    radius: int
    pB: int
    residues: Tuple[int, int, int]
    norm_residues: Tuple[int, int, int]
    resonance_hits: Tuple[int, ...]


def normalize_residues(residues: Tuple[int, int, int]) -> Tuple[int, int, int]:
    return tuple(sorted(residues))


def classify_high_beta_triplet(
    triplet: Tuple[int, int, int],
    prime_flags: List[bool],
    T: Set[int],
    R: Set[int],
    m_mode: str = "min",
    bertrand_mode: str = "first",
    forbid_zero: bool = True,
    require_m_ge: int = 10,
) -> Optional[HighBetaHit]:
    xs = tuple(sorted(triplet))

    if m_mode == "min":
        m = xs[0]
    elif m_mode == "middle":
        m = xs[1]
    elif m_mode == "max":
        m = xs[2]
    else:
        raise ValueError(f"Unbekannter m_mode: {m_mode}")

    if m < require_m_ge:
        return None

    pB = choose_bertrand_prime(m, 2 * m, prime_flags, mode=bertrand_mode)
    if pB is None:
        return None

    residues = tuple(abs(x - pB) for x in xs)

    if forbid_zero and any(r == 0 for r in residues):
        return None

    resonance_hits = tuple(sorted(r for r in residues if r in R))
    prime_hits = tuple(r for r in residues if is_prime(r, prime_flags))
    twin_hits = tuple(r for r in residues if r in T)

    # Hoch-beta:
    # mindestens ein Resonanztreffer,
    # keine Primresiduen,
    # und keine nichtresonanten Zwillingszentren
    nonresonant_twin_hits = tuple(r for r in twin_hits if r not in R)

    if resonance_hits and not prime_hits and not nonresonant_twin_hits:
        return HighBetaHit(
            triplet=triplet,
            sorted_triplet=xs,
            m=m,
            radius=2 * m,
            pB=pB,
            residues=residues,
            norm_residues=normalize_residues(residues),
            resonance_hits=resonance_hits,
        )

    return None


# =========================================================
# Konstruktive Sucher
# =========================================================

def search_constructive_high_beta(
    prime_flags: List[bool],
    T: Set[int],
    R: Set[int],
    m_values: List[int],
    resonance_values: List[int],
    a_offsets: List[int],
    require_m_ge: int = 10,
    bertrand_mode: str = "first",
    limit_hits: int = 200,
) -> List[HighBetaHit]:
    hits: List[HighBetaHit] = []
    seen = set()

    for m in m_values:
        if m < require_m_ge:
            continue

        pB = choose_bertrand_prime(m, 2 * m, prime_flags, mode=bertrand_mode)
        if pB is None:
            continue

        for r in resonance_values:
            # Familie 1: (m, pB+a, pB+r)
            for a in a_offsets:
                trip = tuple(sorted((m, pB + a, pB + r)))
                if trip[0] < 2:
                    continue
                hit = classify_high_beta_triplet(
                    trip, prime_flags, T, R,
                    require_m_ge=require_m_ge,
                    bertrand_mode=bertrand_mode,
                )
                if hit and hit.triplet not in seen:
                    hits.append(hit)
                    seen.add(hit.triplet)
                    if len(hits) >= limit_hits:
                        return hits

            # Familie 2: (m, pB-r, pB+r)
            trip = tuple(sorted((m, pB - r, pB + r)))
            if trip[0] >= 2:
                hit = classify_high_beta_triplet(
                    trip, prime_flags, T, R,
                    require_m_ge=require_m_ge,
                    bertrand_mode=bertrand_mode,
                )
                if hit and hit.triplet not in seen:
                    hits.append(hit)
                    seen.add(hit.triplet)
                    if len(hits) >= limit_hits:
                        return hits

            # Familie 3: (m, pB+r1, pB+r2)
            for r2 in resonance_values:
                trip = tuple(sorted((m, pB + r, pB + r2)))
                if trip[0] < 2:
                    continue
                hit = classify_high_beta_triplet(
                    trip, prime_flags, T, R,
                    require_m_ge=require_m_ge,
                    bertrand_mode=bertrand_mode,
                )
                if hit and hit.triplet not in seen:
                    hits.append(hit)
                    seen.add(hit.triplet)
                    if len(hits) >= limit_hits:
                        return hits

    return hits


# =========================================================
# Ausgabe
# =========================================================

def summarize_high_beta_hits(hits: List[HighBetaHit]) -> None:
    print(f"Anzahl gefundener Hoch-Beta-Tripel: {len(hits)}")
    print()

    norm_counter = Counter(h.norm_residues for h in hits)
    resonance_counter = Counter(h.resonance_hits for h in hits)
    pB_counter = Counter(h.pB for h in hits)
    m_counter = Counter(h.m for h in hits)

    print("Häufigste normalisierte Residuen:")
    for res, cnt in norm_counter.most_common(20):
        print(f"  {str(res):>18}: {cnt}")
    print()

    print("Häufigste Resonanz-Hits:")
    for res, cnt in resonance_counter.most_common(20):
        print(f"  {str(res):>18}: {cnt}")
    print()

    print("Häufigste Bertrand-Anker:")
    for p, cnt in pB_counter.most_common(20):
        print(f"  {str(p):>8}: {cnt}")
    print()

    print("Häufigste Referenzwerte m:")
    for m, cnt in m_counter.most_common(20):
        print(f"  {str(m):>8}: {cnt}")
    print()

    print("Erste Treffer:")
    for h in hits[:30]:
        print(
            f"  T={h.triplet} | sort={h.sorted_triplet} | m={h.m} | "
            f"R={h.radius} | pB={h.pB} | residues={h.residues} | "
            f"norm={h.norm_residues} | R-hit={h.resonance_hits}"
        )


# =========================================================
# Hauptprogramm
# =========================================================

def main() -> None:
    limit = 1_000_000
    prime_flags = prime_sieve(limit + 100)

    T = stable_twin_centers(limit, prime_flags)
    R, RL, RR, Rbi = resonance_classes(T)

    print(f"Schranke: {limit}")
    print(f"|T|={len(T)}, |R|={len(R)}")
    print()

    # Nur höhere Resonanzwerte, um triviale 6er-Randfälle zu vermeiden
    resonance_values = sorted([r for r in R if r >= 12 and r <= 300])[:40]

    # Referenzwerte m in höheren Schalen
    m_values = list(range(10, 301))

    # kleine Offsets relativ zu pB
    a_offsets = [-18, -12, -10, -8, -6, -4, -2, 2, 4, 6, 8, 10, 12, 18]

    hits = search_constructive_high_beta(
        prime_flags=prime_flags,
        T=T,
        R=R,
        m_values=m_values,
        resonance_values=resonance_values,
        a_offsets=a_offsets,
        require_m_ge=10,
        bertrand_mode="first",
        limit_hits=200,
    )

    summarize_high_beta_hits(hits)


if __name__ == "__main__":
    main()