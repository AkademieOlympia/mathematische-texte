from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isqrt
from statistics import mean
from typing import List, Tuple, Optional, Dict, Set


# =========================================================
# Primzahlwerkzeuge
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
# Linienfamilien
# =========================================================

def normalize_residues(residues: Tuple[int, int, int]) -> Tuple[int, int, int]:
    return tuple(sorted(residues))


def classify_line_family(norm_res: Tuple[int, int, int]) -> str:
    """
    Erste Heuristik für Linienfamilien.
    """
    if norm_res == (0, 4, 6):
        return "B0_bound"

    if norm_res == (1, 3, 5):
        return "P1_prime"

    if norm_res == (1, 5, 7):
        return "P2_prime"

    # allgemeine Heuristiken
    a, b, c = norm_res

    if a == 0 and b % 2 == 0 and c % 2 == 0:
        return "bound_like"

    if a == 1 and is_odd_triple(norm_res):
        return "prime_like"

    return "unclassified"


def is_odd_triple(triple: Tuple[int, int, int]) -> bool:
    return all(x % 2 == 1 for x in triple)


# =========================================================
# Einzelklassifikation eines Tripels
# =========================================================

@dataclass
class TripletClassification:
    triplet: Tuple[int, int, int]
    sorted_triplet: Tuple[int, int, int]
    reference_m: int
    radius: int
    bertrand_prime: Optional[int]
    residues: Optional[Tuple[int, int, int]]
    normalized_residues: Optional[Tuple[int, int, int]]
    prime_hits: List[int]
    twin_hits: List[int]
    resonance_hits: List[int]
    dominant_channel: str
    line_family: str


def classify_triplet(
    triplet: Tuple[int, int, int],
    prime_flags: List[bool],
    T: Set[int],
    R: Set[int],
    m_mode: str = "min",
    bertrand_mode: str = "first",
) -> TripletClassification:
    xs = tuple(sorted(triplet))

    if m_mode == "min":
        m = xs[0]
    elif m_mode == "middle":
        m = xs[1]
    elif m_mode == "max":
        m = xs[2]
    elif m_mode == "mean_floor":
        m = sum(xs) // 3
    else:
        raise ValueError(f"Unbekannter m_mode: {m_mode}")

    if m < 2:
        m = 2

    pB = choose_bertrand_prime(m, 2 * m, prime_flags, mode=bertrand_mode)

    if pB is None:
        return TripletClassification(
            triplet=triplet,
            sorted_triplet=xs,
            reference_m=m,
            radius=2 * m,
            bertrand_prime=None,
            residues=None,
            normalized_residues=None,
            prime_hits=[],
            twin_hits=[],
            resonance_hits=[],
            dominant_channel="no_bertrand_prime",
            line_family="none",
        )

    residues = tuple(abs(x - pB) for x in xs)
    norm_res = normalize_residues(residues)

    prime_hits = [r for r in residues if is_prime(r, prime_flags)]
    twin_hits = [r for r in residues if r in T]
    resonance_hits = [r for r in residues if r in R]

    # dominanter Kanal
    if twin_hits:
        dominant = "alpha_twin"
    elif resonance_hits:
        dominant = "beta_resonance"
    elif prime_hits:
        dominant = "gamma_prime"
    else:
        dominant = "none"

    family = classify_line_family(norm_res)

    return TripletClassification(
        triplet=triplet,
        sorted_triplet=xs,
        reference_m=m,
        radius=2 * m,
        bertrand_prime=pB,
        residues=residues,
        normalized_residues=norm_res,
        prime_hits=prime_hits,
        twin_hits=twin_hits,
        resonance_hits=resonance_hits,
        dominant_channel=dominant,
        line_family=family,
    )


# =========================================================
# Batch-Auswertung
# =========================================================

def summarize_triplet_results(results: List[TripletClassification]) -> None:
    channel_counter = Counter(r.dominant_channel for r in results)
    family_counter = Counter(r.line_family for r in results)
    residue_counter = Counter(r.normalized_residues for r in results if r.normalized_residues is not None)
    pB_counter = Counter(r.bertrand_prime for r in results if r.bertrand_prime is not None)

    print(f"Anzahl Tripel: {len(results)}")
    print()

    print("Kanäle:")
    for k in ["alpha_twin", "beta_resonance", "gamma_prime", "none", "no_bertrand_prime"]:
        print(f"  {k:>16}: {channel_counter.get(k, 0)}")
    print()

    print("Linienfamilien:")
    for fam, cnt in family_counter.most_common(15):
        print(f"  {fam:>16}: {cnt}")
    print()

    print("Häufigste normalisierte Residuen:")
    for res, cnt in residue_counter.most_common(20):
        print(f"  {str(res):>16}: {cnt}")
    print()

    print("Häufigste Bertrand-Anker:")
    for p, cnt in pB_counter.most_common(20):
        print(f"  {str(p):>8}: {cnt}")
    print()

    print("Erste Beispiele:")
    shown = 0
    for r in results:
        print(
            f"  T={r.triplet} | sort={r.sorted_triplet} | m={r.reference_m} | "
            f"R={r.radius} | pB={r.bertrand_prime} | residues={r.residues} | "
            f"norm={r.normalized_residues} | channel={r.dominant_channel} | family={r.line_family}"
        )
        shown += 1
        if shown >= 20:
            break


# =========================================================
# Testtripel erzeugen
# =========================================================

def generate_standard_triplets_from_fourling_centers(
    four_centers: List[int],
    max_count: Optional[int] = None
) -> List[Tuple[int, int, int]]:
    """
    Aus (p,p+2,p+6,p+8) entstehen vier natürliche Triaden:
    (p,p+2,p+6), (p,p+2,p+8), (p,p+6,p+8), (p+2,p+6,p+8)
    """
    triplets = []
    for c in four_centers:
        p = c - 4
        quad = (p, p + 2, p + 6, p + 8)
        triplets.extend([
            (quad[0], quad[1], quad[2]),
            (quad[0], quad[1], quad[3]),
            (quad[0], quad[2], quad[3]),
            (quad[1], quad[2], quad[3]),
        ])
        if max_count is not None and len(triplets) >= max_count:
            return triplets[:max_count]
    return triplets


def generate_shifted_triplets(
    base_triplets: List[Tuple[int, int, int]],
    shifts: List[int]
) -> List[Tuple[int, int, int]]:
    out = []
    for t in base_triplets:
        for s in shifts:
            out.append((t[0] + s, t[1] + s, t[2] + s))
    return out


# =========================================================
# Vierlingszentren
# =========================================================

def standard_fourling_centers(limit: int, prime_flags: List[bool]) -> List[int]:
    centers = []
    for p in range(2, limit - 8 + 1):
        if prime_flags[p] and prime_flags[p + 2] and prime_flags[p + 6] and prime_flags[p + 8]:
            centers.append(p + 4)
    return centers


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
    print(f"Vierlingszentren: {len(four_centers)}")
    print(f"|T|={len(T)}, |R|={len(R)}")
    print()

    triplets_A = generate_standard_triplets_from_fourling_centers(four_centers[:100], max_count=400)
    triplets_B = generate_shifted_triplets(triplets_A[:100], shifts=[-3, -1, 1, 3])

    print("=== Test A: natürliche Tripel aus Vierlingen ===")
    results_A = [
        classify_triplet(t, prime_flags, T, R, m_mode="min", bertrand_mode="first")
        for t in triplets_A
    ]
    summarize_triplet_results(results_A)

    print("\n=== Test B: verschobene Tripel ===")
    results_B = [
        classify_triplet(t, prime_flags, T, R, m_mode="min", bertrand_mode="first")
        for t in triplets_B
    ]
    summarize_triplet_results(results_B)


if __name__ == "__main__":
    main()