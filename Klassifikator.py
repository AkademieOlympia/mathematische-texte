from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isqrt
from typing import List, Tuple, Optional, Set


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
# Serienklassifikation
# =========================================================

@dataclass
class SeriesClassification:
    residues: Tuple[int, int, int]
    normalized: Tuple[int, int, int]
    series: str                 # alpha / beta / gamma / mixed / unknown
    subtype: str
    confidence: str             # high / medium / low
    reason: str


def classify_series(
    residues: Tuple[int, int, int],
    prime_flags: List[bool],
    T: Set[int],
    R: Set[int],
) -> SeriesClassification:
    nr = normalize_residues(residues)
    a, b, c = nr

    # -----------------------------------------------------
    # 1. Harte exakte Muster
    # -----------------------------------------------------
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
        return SeriesClassification(
            residues=residues,
            normalized=nr,
            series="alpha",
            subtype=alpha_exact[nr],
            confidence="high",
            reason="exaktes gebundenes Linienmuster",
        )

    if nr in beta_exact:
        return SeriesClassification(
            residues=residues,
            normalized=nr,
            series="beta",
            subtype=beta_exact[nr],
            confidence="high",
            reason="exaktes resonanzvermitteltes Linienmuster",
        )

    if nr in gamma_exact:
        return SeriesClassification(
            residues=residues,
            normalized=nr,
            series="gamma",
            subtype=gamma_exact[nr],
            confidence="high",
            reason="exaktes primares Linienmuster",
        )

    # -----------------------------------------------------
    # 2. Strukturmerkmale
    # -----------------------------------------------------
    prime_hits = [x for x in nr if is_prime(x, prime_flags)]
    twin_hits = [x for x in nr if x in T]
    resonance_hits = [x for x in nr if x in R]

    even_count = sum(x % 2 == 0 for x in nr)
    odd_count = 3 - even_count

    # Alpha-artig:
    # viel gerade Struktur, oft 0 oder 2, keine führende 1 nötig
    if (a == 0 and even_count >= 2) or (even_count == 3 and any(x in T for x in nr)):
        return SeriesClassification(
            residues=residues,
            normalized=nr,
            series="alpha",
            subtype="alpha_like",
            confidence="medium",
            reason="gebundene gerade Reststruktur mit Zwillingsnähe",
        )

    # Beta-artig:
    # führende 1, mindestens ein Resonanzwert, nicht primdominiert
    if a == 1 and len(resonance_hits) >= 1:
        return SeriesClassification(
            residues=residues,
            normalized=nr,
            series="beta",
            subtype="beta_like",
            confidence="medium",
            reason="ankernahe 1 mit Resonanzrest(en)",
        )

    # Gamma-artig:
    # kleine ungerade Primstruktur, keine Resonanzdominanz
    if a == 1 and odd_count >= 2 and len(prime_hits) >= 1 and len(resonance_hits) == 0:
        return SeriesClassification(
            residues=residues,
            normalized=nr,
            series="gamma",
            subtype="gamma_like",
            confidence="medium",
            reason="primartige ungerade Reststruktur ohne Resonanzdominanz",
        )

    # Gemischte Fälle:
    if resonance_hits and prime_hits:
        return SeriesClassification(
            residues=residues,
            normalized=nr,
            series="mixed",
            subtype="beta_gamma_mix",
            confidence="low",
            reason="gleichzeitig Prim- und Resonanzanteile",
        )

    if twin_hits and resonance_hits:
        return SeriesClassification(
            residues=residues,
            normalized=nr,
            series="mixed",
            subtype="alpha_beta_mix",
            confidence="low",
            reason="gleichzeitig Zwillings- und Resonanzanteile",
        )

    return SeriesClassification(
        residues=residues,
        normalized=nr,
        series="unknown",
        subtype="unclassified",
        confidence="low",
        reason="kein bekanntes exaktes oder strukturelles Muster",
    )


# =========================================================
# Batch-Auswertung
# =========================================================

def summarize_series_classifications(results: List[SeriesClassification]) -> None:
    series_counter = Counter(r.series for r in results)
    subtype_counter = Counter(r.subtype for r in results)
    norm_counter = Counter(r.normalized for r in results)

    print(f"Anzahl klassifizierter Muster: {len(results)}")
    print()

    print("Serienklassen:")
    for key in ["alpha", "beta", "gamma", "mixed", "unknown"]:
        print(f"  {key:>8}: {series_counter.get(key, 0)}")
    print()

    print("Häufigste Untertypen:")
    for st, cnt in subtype_counter.most_common(20):
        print(f"  {st:>20}: {cnt}")
    print()

    print("Häufigste Muster:")
    for nr, cnt in norm_counter.most_common(20):
        print(f"  {str(nr):>20}: {cnt}")
    print()

    print("Erste Beispiele:")
    for r in results[:30]:
        print(
            f"  residues={r.residues} | norm={r.normalized} | "
            f"series={r.series} | subtype={r.subtype} | "
            f"confidence={r.confidence} | reason={r.reason}"
        )


# =========================================================
# Demo
# =========================================================

def main() -> None:
    limit = 1_000_000
    prime_flags = prime_sieve(limit + 100)
    T = stable_twin_centers(limit, prime_flags)
    R, RL, RR, Rbi = resonance_classes(T)

    demo_patterns = [
        (0, 2, 4),
        (0, 2, 6),
        (2, 4, 6),
        (1, 3, 5),
        (1, 5, 7),
        (1, 1, 5),
        (1, 3, 3),
        (1, 6, 12),
        (1, 12, 18),
        (1, 12, 102),
        (1, 18, 108),
        (1, 102, 192),
        (1, 6, 8),
        (1, 6, 14),
        (2, 8, 12),
        (4, 6, 12),
    ]

    results = [classify_series(p, prime_flags, T, R) for p in demo_patterns]
    summarize_series_classifications(results)


if __name__ == "__main__":
    main()