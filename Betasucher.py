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


def choose_bertrand_prime(
    lo: int,
    hi: int,
    prime_flags: List[bool],
    mode: str = "first",
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
# Beta-Kandidat
# =========================================================

@dataclass
class BetaTripletHit:
    triplet: Tuple[int, int, int]
    sorted_triplet: Tuple[int, int, int]
    reference_m: int
    radius: int
    bertrand_prime: int
    residues: Tuple[int, int, int]
    normalized_residues: Tuple[int, int, int]
    resonance_residues: Tuple[int, ...]


def normalize_residues(residues: Tuple[int, int, int]) -> Tuple[int, int, int]:
    return tuple(sorted(residues))


def classify_triplet_for_beta(
    triplet: Tuple[int, int, int],
    prime_flags: List[bool],
    T: Set[int],
    R: Set[int],
    m_mode: str = "min",
    bertrand_mode: str = "first",
) -> Optional[BetaTripletHit]:
    xs = tuple(sorted(triplet))

    if m_mode == "min":
        m = xs[0]
    elif m_mode == "middle":
        m = xs[1]
    elif m_mode == "max":
        m = xs[2]
    else:
        raise ValueError(f"Unbekannter m_mode: {m_mode}")

    if m < 2:
        m = 2

    pB = choose_bertrand_prime(m, 2 * m, prime_flags, mode=bertrand_mode)
    if pB is None:
        return None

    residues = tuple(abs(x - pB) for x in xs)

    # echte Beta-Kandidaten:
    # - mindestens ein Residuum in R
    # - aber kein Residuum prim
    resonance_hits = tuple(r for r in residues if r in R)
    prime_hits = tuple(r for r in residues if is_prime(r, prime_flags))

    if resonance_hits and not prime_hits:
        return BetaTripletHit(
            triplet=triplet,
            sorted_triplet=xs,
            reference_m=m,
            radius=2 * m,
            bertrand_prime=pB,
            residues=residues,
            normalized_residues=normalize_residues(residues),
            resonance_residues=tuple(sorted(resonance_hits)),
        )

    return None


# =========================================================
# Suchraum
# =========================================================

def search_beta_triplets_box(
    max_value: int,
    prime_flags: List[bool],
    T: Set[int],
    R: Set[int],
    m_mode: str = "min",
    bertrand_mode: str = "first",
    limit_hits: int = 200
) -> List[BetaTripletHit]:
    hits: List[BetaTripletHit] = []

    for x in range(2, max_value + 1):
        for y in range(x, max_value + 1):
            for z in range(y, max_value + 1):
                hit = classify_triplet_for_beta(
                    (x, y, z),
                    prime_flags,
                    T,
                    R,
                    m_mode=m_mode,
                    bertrand_mode=bertrand_mode,
                )
                if hit is not None:
                    hits.append(hit)
                    if len(hits) >= limit_hits:
                        return hits
    return hits


def search_beta_triplets_targeted(
    resonance_values: List[int],
    prime_flags: List[bool],
    T: Set[int],
    R: Set[int],
    offsets: List[int],
    m_mode: str = "min",
    bertrand_mode: str = "first",
    limit_hits: int = 200
) -> List[BetaTripletHit]:
    """
    Baut gezielt Tripel um Resonanzwerte herum.
    """
    hits: List[BetaTripletHit] = []

    for r in resonance_values:
        for a in offsets:
            for b in offsets:
                for c in offsets:
                    triplet = tuple(sorted((r + a, r + b, r + c)))
                    if triplet[0] < 2:
                        continue

                    hit = classify_triplet_for_beta(
                        triplet,
                        prime_flags,
                        T,
                        R,
                        m_mode=m_mode,
                        bertrand_mode=bertrand_mode,
                    )
                    if hit is not None:
                        hits.append(hit)
                        if len(hits) >= limit_hits:
                            return hits
    return hits


# =========================================================
# Ausgabe
# =========================================================

def summarize_beta_hits(hits: List[BetaTripletHit]) -> None:
    print(f"Anzahl gefundener Beta-aktiver Tripel: {len(hits)}")
    print()

    residue_counter = Counter(h.normalized_residues for h in hits)
    resonance_counter = Counter(h.resonance_residues for h in hits)
    pB_counter = Counter(h.bertrand_prime for h in hits)

    print("Häufigste normalisierte Residuen:")
    for res, cnt in residue_counter.most_common(20):
        print(f"  {str(res):>16}: {cnt}")
    print()

    print("Häufigste Resonanz-Residuen:")
    for res, cnt in resonance_counter.most_common(20):
        print(f"  {str(res):>16}: {cnt}")
    print()

    print("Häufigste Bertrand-Anker:")
    for p, cnt in pB_counter.most_common(20):
        print(f"  {str(p):>8}: {cnt}")
    print()

    print("Erste Treffer:")
    for h in hits[:30]:
        print(
            f"  T={h.triplet} | sort={h.sorted_triplet} | m={h.reference_m} | "
            f"R={h.radius} | pB={h.bertrand_prime} | residues={h.residues} | "
            f"norm={h.normalized_residues} | R-hit={h.resonance_residues}"
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

    resonance_seed = sorted([6, 12, 18, 102, 108, 192, 198, 822, 828, 1482, 1488])
    offsets = [-12, -6, -4, -2, 0, 2, 4, 6, 12]

    print("=== Gezielsuche um Resonanzwerte ===")
    hits_targeted = search_beta_triplets_targeted(
        resonance_values=resonance_seed,
        prime_flags=prime_flags,
        T=T,
        R=R,
        offsets=offsets,
        m_mode="min",
        bertrand_mode="first",
        limit_hits=200,
    )
    summarize_beta_hits(hits_targeted)

    print("\n=== Box-Suche im kleinen Bereich ===")
    hits_box = search_beta_triplets_box(
        max_value=120,
        prime_flags=prime_flags,
        T=T,
        R=R,
        m_mode="min",
        bertrand_mode="first",
        limit_hits=100,
    )
    summarize_beta_hits(hits_box)


if __name__ == "__main__":
    main()