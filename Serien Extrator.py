from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
from math import isqrt
from typing import List, Tuple, Set, Dict

import matplotlib.pyplot as plt


# =========================================================
# Grundwerkzeuge
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


def normalize_residues(residues: Tuple[int, int, int]) -> Tuple[int, int, int]:
    return tuple(sorted(residues))


# =========================================================
# Klassifikation
# =========================================================

@dataclass
class SeriesInfo:
    normalized: Tuple[int, int, int]
    series: str
    subtype: str
    confidence: str


def classify_series(
    residues: Tuple[int, int, int],
    prime_flags: List[bool],
    T: Set[int],
    R: Set[int],
) -> SeriesInfo:
    nr = normalize_residues(residues)
    a, b, c = nr

    alpha_exact = {
        (0, 2, 4): "alpha_bound_024",
        (0, 2, 6): "alpha_bound_026",
        (2, 4, 6): "alpha_bound_246",
        (0, 4, 6): "alpha_bound_046",
    }

    beta_exact = {
        (1, 6, 12): "beta_res_1612",
        (1, 12, 18): "beta_res_11218",
        (1, 12, 102): "beta_res_112102",
        (1, 12, 108): "beta_res_112108",
        (1, 12, 192): "beta_res_112192",
        (1, 12, 198): "beta_res_112198",
        (1, 18, 102): "beta_res_118102",
        (1, 18, 108): "beta_res_118108",
        (1, 18, 192): "beta_res_118192",
        (1, 18, 198): "beta_res_118198",
        (1, 102, 108): "beta_res_1102108",
        (1, 102, 192): "beta_res_1102192",
        (1, 102, 198): "beta_res_1102198",
        (1, 108, 192): "beta_res_1108192",
        (1, 108, 198): "beta_res_1108198",
    }

    gamma_exact = {
        (1, 3, 5): "gamma_prime_135",
        (1, 5, 7): "gamma_prime_157",
        (1, 1, 5): "gamma_prime_115",
        (1, 1, 7): "gamma_prime_117",
        (1, 3, 3): "gamma_prime_133",
        (3, 3, 5): "gamma_prime_335",
    }

    if nr in alpha_exact:
        return SeriesInfo(nr, "alpha", alpha_exact[nr], "high")
    if nr in beta_exact:
        return SeriesInfo(nr, "beta", beta_exact[nr], "high")
    if nr in gamma_exact:
        return SeriesInfo(nr, "gamma", gamma_exact[nr], "high")

    prime_hits = [x for x in nr if is_prime(x, prime_flags)]
    twin_hits = [x for x in nr if x in T]
    resonance_hits = [x for x in nr if x in R]
    even_count = sum(x % 2 == 0 for x in nr)
    odd_count = 3 - even_count

    if (a == 0 and even_count >= 2) or (even_count == 3 and any(x in T for x in nr)):
        return SeriesInfo(nr, "alpha", "alpha_like", "medium")

    if a == 1 and len(resonance_hits) >= 1:
        return SeriesInfo(nr, "beta", "beta_like", "medium")

    if a == 1 and odd_count >= 2 and len(prime_hits) >= 1 and len(resonance_hits) == 0:
        return SeriesInfo(nr, "gamma", "gamma_like", "medium")

    return SeriesInfo(nr, "unknown", "unclassified", "low")


# =========================================================
# Quantenzahlen
# =========================================================

@dataclass
class QuantumNumbers:
    residues: Tuple[int, int, int]
    normalized: Tuple[int, int, int]
    series: str
    subtype: str
    confidence: str
    n: int
    ell: int
    rho: int
    mu: int


def extract_quantum_numbers(
    residues: Tuple[int, int, int],
    prime_flags: List[bool],
    T: Set[int],
    R: Set[int],
    force_alpha_rho_zero: bool = False,
) -> QuantumNumbers:
    info = classify_series(residues, prime_flags, T, R)
    nr = info.normalized
    a, b, c = nr

    if info.series == "alpha":
        n = 1
    elif info.series == "beta":
        n = 2
    elif info.series == "gamma":
        n = 3
    else:
        n = 0

    resonance_hits = sorted([x for x in nr if x in R])

    if info.series == "alpha" and force_alpha_rho_zero:
        rho = 0
        mu = 0
    else:
        rho = max(resonance_hits) if resonance_hits else 0
        mu = len(resonance_hits)

    ell = b

    return QuantumNumbers(
        residues=residues,
        normalized=nr,
        series=info.series,
        subtype=info.subtype,
        confidence=info.confidence,
        n=n,
        ell=ell,
        rho=rho,
        mu=mu,
    )


# =========================================================
# Scanner nach festem rho
# =========================================================

def scan_rho_families(
    patterns: List[Tuple[int, int, int]],
    prime_flags: List[bool],
    T: Set[int],
    R: Set[int],
    force_alpha_rho_zero: bool = False,
) -> Dict[Tuple[str, int], List[QuantumNumbers]]:
    buckets: Dict[Tuple[str, int], List[QuantumNumbers]] = defaultdict(list)

    for pat in patterns:
        q = extract_quantum_numbers(
            pat, prime_flags, T, R,
            force_alpha_rho_zero=force_alpha_rho_zero,
        )
        key = (q.series, q.rho)
        buckets[key].append(q)

    for key in buckets:
        buckets[key].sort(key=lambda x: (x.ell, x.normalized))

    return buckets


def summarize_rho_families(
    buckets: Dict[Tuple[str, int], List[QuantumNumbers]],
    only_series: str | None = None,
    min_size: int = 1,
) -> None:
    keys = sorted(buckets.keys(), key=lambda k: (k[0], k[1]))

    if only_series is not None:
        keys = [k for k in keys if k[0] == only_series]

    print("Serien nach festem rho")
    print("-" * 100)

    for series, rho in keys:
        items = buckets[(series, rho)]
        if len(items) < min_size:
            continue

        ell_values = [q.ell for q in items]
        ell_unique = sorted(set(ell_values))
        mu_counter = Counter(q.mu for q in items)
        subtype_counter = Counter(q.subtype for q in items)
        patterns = [q.normalized for q in items]

        print(f"Serie={series:>6} | rho={rho:>4} | Anzahl={len(items):>3}")
        print(f"  ell-Werte      : {ell_unique}")
        print(f"  mu-Verteilung  : {dict(mu_counter)}")
        print(f"  Untertypen     : {dict(subtype_counter)}")
        print(f"  Beispiele      : {patterns[:10]}")
        print()


# =========================================================
# Demo
# =========================================================

def main() -> None:
    limit = 1_000_000
    prime_flags = prime_sieve(limit + 100)
    T = stable_twin_centers(limit, prime_flags)
    R, RL, RR, Rbi = resonance_classes(T)

    # Demo-Menge: harte + weiche Beispiele
    patterns = [
        # alpha
        (0, 2, 4),
        (0, 2, 6),
        (2, 4, 6),
        (2, 8, 12),
        (4, 6, 12),

        # gamma
        (1, 3, 5),
        (1, 5, 7),
        (1, 1, 5),
        (1, 3, 3),

        # beta
        (1, 6, 12),
        (1, 8, 12),
        (1, 10, 12),
        (1, 12, 12),
        (1, 12, 18),
        (1, 6, 18),
        (1, 8, 18),
        (1, 10, 18),
        (1, 18, 18),
        (1, 12, 102),
        (1, 18, 102),
        (1, 102, 102),
        (1, 18, 108),
        (1, 102, 108),
        (1, 102, 192),
        (1, 108, 192),
        (1, 12, 198),
        (1, 18, 198),
        (1, 102, 198),
        (1, 6, 8),
        (1, 6, 14),
    ]

    buckets = scan_rho_families(
        patterns,
        prime_flags,
        T,
        R,
        force_alpha_rho_zero=False,
    )

    summarize_rho_families(buckets, only_series="beta", min_size=1)

    print()
    print_beta_ladder(beta_ladder)
    plot_beta_ladder(beta_ladder, savepath="beta_schalenleiter.png")


if __name__ == "__main__":
    main()