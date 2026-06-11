from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import isqrt
from statistics import mean
from typing import Dict, List, Tuple


@dataclass
class CenterStats:
    M: int
    twin_count: int
    quadruple_count_all: int
    quadruple_count_small_delta: int
    exact_symmetric_count: int
    near_symmetric_count_6: int
    near_symmetric_count_12: int
    weighted_cupola_all: float
    weighted_cupola_small_delta: float
    delta_histogram_small: Dict[int, int]
    eta_histogram_small: Dict[int, int]


def prime_sieve(limit: int) -> List[bool]:
    """
    Bool-Sieb für Primzahlen bis einschließlich 'limit'.
    """
    if limit < 2:
        return [False] * (limit + 1)

    sieve = [True] * (limit + 1)
    sieve[0] = False
    sieve[1] = False

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            step = p
            sieve[start : limit + 1 : step] = [False] * (((limit - start) // step) + 1)

    return sieve


def collect_twin_primes(prime_flags: List[bool]) -> List[Tuple[int, int]]:
    """
    Alle Primzahlzwillinge (p, p+2) innerhalb des Siebbereichs.
    """
    limit = len(prime_flags) - 1
    twins = []
    for p in range(2, limit - 1):
        if prime_flags[p] and prime_flags[p + 2]:
            twins.append((p, p + 2))
    return twins


def twins_by_center(twins: List[Tuple[int, int]]) -> Dict[int, List[Tuple[int, int]]]:
    """
    Ordnet jedem Zwillingszentrum z = p+1 die Zwillinge (p, p+2) zu.
    """
    result: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
    for p, q in twins:
        result[p + 1].append((p, q))
    return result


def analyze_center(
    M: int,
    radius: int,
    center_to_twins: Dict[int, List[Tuple[int, int]]],
    delta_max: int,
) -> CenterStats:
    """
    Analysiert die Schale [M-radius, M+radius] um ein Zentrum M.

    Gemessen werden:
    - alle Quadrupel links/rechts von M
    - kleine Delta-Werte bis delta_max
    - exakte Symmetrie
    - Nahsymmetrie eta <= 6, eta <= 12
    """
    local_centers = [z for z in center_to_twins.keys() if abs(z - M) <= radius]
    local_centers.sort()

    local_twins = []
    for z in local_centers:
        local_twins.extend(center_to_twins[z])

    left_centers = [z for z in local_centers if z < M]
    right_centers = [z for z in local_centers if z > M]

    quadruple_count_all = 0
    quadruple_count_small_delta = 0
    exact_symmetric_count = 0
    near_symmetric_count_6 = 0
    near_symmetric_count_12 = 0

    weighted_cupola_all = 0.0
    weighted_cupola_small_delta = 0.0

    delta_hist_small = Counter()
    eta_hist_small = Counter()

    for z1 in left_centers:
        for z2 in right_centers:
            delta = z2 - z1
            if delta <= 0:
                continue
            if delta % 6 != 0:
                continue

            quadruple_count_all += 1
            weighted_cupola_all += 1.0 / delta

            # Nahsymmetrie
            eta = abs((M - z1) - (z2 - M))

            if eta == 0:
                exact_symmetric_count += 1
            if eta <= 6:
                near_symmetric_count_6 += 1
            if eta <= 12:
                near_symmetric_count_12 += 1

            # nur kleine delta separat zählen
            if delta <= delta_max:
                quadruple_count_small_delta += 1
                weighted_cupola_small_delta += 1.0 / delta
                delta_hist_small[delta] += 1
                eta_hist_small[eta] += 1

    return CenterStats(
        M=M,
        twin_count=len(local_twins),
        quadruple_count_all=quadruple_count_all,
        quadruple_count_small_delta=quadruple_count_small_delta,
        exact_symmetric_count=exact_symmetric_count,
        near_symmetric_count_6=near_symmetric_count_6,
        near_symmetric_count_12=near_symmetric_count_12,
        weighted_cupola_all=weighted_cupola_all,
        weighted_cupola_small_delta=weighted_cupola_small_delta,
        delta_histogram_small=dict(sorted(delta_hist_small.items())),
        eta_histogram_small=dict(sorted(eta_hist_small.items())),
    )


def analyze_centers(
    centers: List[int],
    radius: int,
    center_to_twins: Dict[int, List[Tuple[int, int]]],
    delta_max: int,
) -> List[CenterStats]:
    return [
        analyze_center(
            M=M,
            radius=radius,
            center_to_twins=center_to_twins,
            delta_max=delta_max,
        )
        for M in centers
    ]


def summarize_family(name: str, stats: List[CenterStats]) -> None:
    """
    Saubere Zusammenfassung einer Zentrumsfamilie.
    """
    if not stats:
        print(f"\n=== {name} ===")
        print("Keine Daten.")
        return

    avg_twins = mean(s.twin_count for s in stats)
    avg_quad_all = mean(s.quadruple_count_all for s in stats)
    avg_quad_small = mean(s.quadruple_count_small_delta for s in stats)
    avg_exact = mean(s.exact_symmetric_count for s in stats)
    avg_near6 = mean(s.near_symmetric_count_6 for s in stats)
    avg_near12 = mean(s.near_symmetric_count_12 for s in stats)
    avg_kw_all = mean(s.weighted_cupola_all for s in stats)
    avg_kw_small = mean(s.weighted_cupola_small_delta for s in stats)

    combined_delta = Counter()
    combined_eta = Counter()
    for s in stats:
        combined_delta.update(s.delta_histogram_small)
        combined_eta.update(s.eta_histogram_small)

    print(f"\n=== Zusammenfassung: {name} ===")
    print(f"Anzahl Zentren:                         {len(stats)}")
    print(f"Mittlere lokale Zwillingszahl:         {avg_twins:.6f}")
    print(f"Mittlere Quadrupelzahl (alle):         {avg_quad_all:.6f}")
    print(f"Mittlere Quadrupelzahl (kleines δ):    {avg_quad_small:.6f}")
    print(f"Mittlere exakte Symmetrien:            {avg_exact:.6f}")
    print(f"Mittlere Nahsymmetrien η<=6:           {avg_near6:.6f}")
    print(f"Mittlere Nahsymmetrien η<=12:          {avg_near12:.6f}")
    print(f"Mittlere gewichtete Kuppel (alle):     {avg_kw_all:.6f}")
    print(f"Mittlere gewichtete Kuppel (kleines δ):{avg_kw_small:.6f}")

    print("\nHäufigste kleine δ-Werte:")
    for delta, count in combined_delta.most_common(12):
        print(f"  δ = {delta:>3} : {count}")

    print("\nHäufigste kleine η-Werte:")
    for eta, count in combined_eta.most_common(12):
        print(f"  η = {eta:>3} : {count}")


def print_top_centers(
    name: str,
    stats: List[CenterStats],
    top_k: int = 10,
    use_small_delta: bool = True,
) -> None:
    """
    Zeigt die auffälligsten Zentren.
    """
    if use_small_delta:
        ranked = sorted(
            stats,
            key=lambda s: (
                s.weighted_cupola_small_delta,
                s.quadruple_count_small_delta,
                s.exact_symmetric_count,
                s.near_symmetric_count_6,
            ),
            reverse=True,
        )
        mode_text = "kleines δ"
    else:
        ranked = sorted(
            stats,
            key=lambda s: (
                s.weighted_cupola_all,
                s.quadruple_count_all,
                s.exact_symmetric_count,
                s.near_symmetric_count_6,
            ),
            reverse=True,
        )
        mode_text = "alle δ"

    print(f"\nTop-{top_k} Zentren [{name}] nach gewichteter Kuppelfunktion ({mode_text}):")
    for s in ranked[:top_k]:
        print(
            f"M={s.M:>8} | Zwillinge={s.twin_count:>3} | "
            f"Q_all={s.quadruple_count_all:>4} | "
            f"Q_small={s.quadruple_count_small_delta:>4} | "
            f"symm={s.exact_symmetric_count:>3} | "
            f"η<=6={s.near_symmetric_count_6:>3} | "
            f"Ksmall={s.weighted_cupola_small_delta:.6f}"
        )


def build_centers_by_divisibility(
    min_center: int,
    max_center: int,
    divisibility: int,
) -> List[int]:
    """
    Zentren M im Intervall [min_center, max_center] mit divisibility | M.
    """
    start = ((min_center + divisibility - 1) // divisibility) * divisibility
    return list(range(start, max_center + 1, divisibility))


def build_control_centers(
    min_center: int,
    max_center: int,
    residue: int,
    modulus: int = 30,
) -> List[int]:
    """
    Kontrollzentren M ≡ residue mod modulus im Intervall [min_center, max_center].
    """
    centers = []
    first = min_center
    while first % modulus != residue % modulus:
        first += 1
    for M in range(first, max_center + 1, modulus):
        centers.append(M)
    return centers


def main() -> None:
    # Parameter
    min_center = 10_000
    max_center = 300_000
    radius = 300
    delta_max = 120

    # Siebgrenze
    sieve_limit = max_center + radius + 10
    prime_flags = prime_sieve(sieve_limit)
    twins = collect_twin_primes(prime_flags)
    center_map = twins_by_center(twins)

    # Familien
    centers_30 = build_centers_by_divisibility(min_center, max_center, 30)
    centers_210 = build_centers_by_divisibility(min_center, max_center, 210)
    centers_2310 = build_centers_by_divisibility(min_center, max_center, 2310)

    # Kontrollen
    controls_1 = build_control_centers(min_center, max_center, residue=1, modulus=30)
    controls_7 = build_control_centers(min_center, max_center, residue=7, modulus=30)
    controls_11 = build_control_centers(min_center, max_center, residue=11, modulus=30)
    controls_13 = build_control_centers(min_center, max_center, residue=13, modulus=30)

    # Analyse
    stats_30 = analyze_centers(centers_30, radius, center_map, delta_max)
    stats_210 = analyze_centers(centers_210, radius, center_map, delta_max)
    stats_2310 = analyze_centers(centers_2310, radius, center_map, delta_max)

    stats_c1 = analyze_centers(controls_1, radius, center_map, delta_max)
    stats_c7 = analyze_centers(controls_7, radius, center_map, delta_max)
    stats_c11 = analyze_centers(controls_11, radius, center_map, delta_max)
    stats_c13 = analyze_centers(controls_13, radius, center_map, delta_max)

    # Ausgabe
    print(f"Parameter: min_center={min_center}, max_center={max_center}, radius={radius}, delta_max={delta_max}")

    summarize_family("30 | M", stats_30)
    summarize_family("210 | M", stats_210)
    summarize_family("2310 | M", stats_2310)

    summarize_family("Kontrolle M ≡ 1 mod 30", stats_c1)
    summarize_family("Kontrolle M ≡ 7 mod 30", stats_c7)
    summarize_family("Kontrolle M ≡ 11 mod 30", stats_c11)
    summarize_family("Kontrolle M ≡ 13 mod 30", stats_c13)

    print_top_centers("30 | M", stats_30, top_k=12, use_small_delta=True)
    print_top_centers("210 | M", stats_210, top_k=12, use_small_delta=True)
    print_top_centers("2310 | M", stats_2310, top_k=12, use_small_delta=True)


if __name__ == "__main__":
    main()