from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from math import isqrt, sqrt
from pathlib import Path
from statistics import mean, pstdev
from typing import Dict, Iterable, List, Tuple


# =========================
# Datenklassen
# =========================

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


@dataclass
class FamilySummary:
    name: str
    n: int
    mean_twin_count: float
    mean_quad_all: float
    mean_quad_small: float
    mean_exact_sym: float
    mean_near6: float
    mean_near12: float
    mean_kw_all: float
    mean_kw_small: float

    std_twin_count: float
    std_quad_all: float
    std_quad_small: float
    std_exact_sym: float
    std_near6: float
    std_near12: float
    std_kw_all: float
    std_kw_small: float

    se_twin_count: float
    se_quad_all: float
    se_quad_small: float
    se_exact_sym: float
    se_near6: float
    se_near12: float
    se_kw_all: float
    se_kw_small: float


# =========================
# Primzahlen / Zwillinge
# =========================

def prime_sieve(limit: int) -> List[bool]:
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
    limit = len(prime_flags) - 1
    twins = []
    for p in range(2, limit - 1):
        if prime_flags[p] and prime_flags[p + 2]:
            twins.append((p, p + 2))
    return twins


def twins_by_center(twins: List[Tuple[int, int]]) -> Dict[int, List[Tuple[int, int]]]:
    result: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
    for p, q in twins:
        result[p + 1].append((p, q))
    return result


# =========================
# Analyse einzelner Zentren
# =========================

def analyze_center(
    M: int,
    radius: int,
    center_to_twins: Dict[int, List[Tuple[int, int]]],
    delta_max: int,
) -> CenterStats:
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

            eta = abs((M - z1) - (z2 - M))

            if eta == 0:
                exact_symmetric_count += 1
            if eta <= 6:
                near_symmetric_count_6 += 1
            if eta <= 12:
                near_symmetric_count_12 += 1

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


# =========================
# Zentren erzeugen
# =========================

def build_centers_by_divisibility(
    min_center: int,
    max_center: int,
    divisibility: int,
) -> List[int]:
    start = ((min_center + divisibility - 1) // divisibility) * divisibility
    return list(range(start, max_center + 1, divisibility))


def build_control_centers(
    min_center: int,
    max_center: int,
    residue: int,
    modulus: int = 30,
) -> List[int]:
    centers = []
    first = min_center
    while first % modulus != residue % modulus:
        first += 1
    for M in range(first, max_center + 1, modulus):
        centers.append(M)
    return centers


# =========================
# Statistik-Helfer
# =========================

def safe_pstdev(values: List[float]) -> float:
    if len(values) <= 1:
        return 0.0
    return pstdev(values)


def safe_se(std_value: float, n: int) -> float:
    if n <= 0:
        return 0.0
    return std_value / sqrt(n)


def summarize_family(name: str, stats: List[CenterStats]) -> FamilySummary:
    if not stats:
        return FamilySummary(
            name=name,
            n=0,
            mean_twin_count=0.0,
            mean_quad_all=0.0,
            mean_quad_small=0.0,
            mean_exact_sym=0.0,
            mean_near6=0.0,
            mean_near12=0.0,
            mean_kw_all=0.0,
            mean_kw_small=0.0,
            std_twin_count=0.0,
            std_quad_all=0.0,
            std_quad_small=0.0,
            std_exact_sym=0.0,
            std_near6=0.0,
            std_near12=0.0,
            std_kw_all=0.0,
            std_kw_small=0.0,
            se_twin_count=0.0,
            se_quad_all=0.0,
            se_quad_small=0.0,
            se_exact_sym=0.0,
            se_near6=0.0,
            se_near12=0.0,
            se_kw_all=0.0,
            se_kw_small=0.0,
        )

    twin_vals = [s.twin_count for s in stats]
    quad_all_vals = [s.quadruple_count_all for s in stats]
    quad_small_vals = [s.quadruple_count_small_delta for s in stats]
    exact_vals = [s.exact_symmetric_count for s in stats]
    near6_vals = [s.near_symmetric_count_6 for s in stats]
    near12_vals = [s.near_symmetric_count_12 for s in stats]
    kw_all_vals = [s.weighted_cupola_all for s in stats]
    kw_small_vals = [s.weighted_cupola_small_delta for s in stats]

    n = len(stats)

    std_twin = safe_pstdev(twin_vals)
    std_quad_all = safe_pstdev(quad_all_vals)
    std_quad_small = safe_pstdev(quad_small_vals)
    std_exact = safe_pstdev(exact_vals)
    std_near6 = safe_pstdev(near6_vals)
    std_near12 = safe_pstdev(near12_vals)
    std_kw_all = safe_pstdev(kw_all_vals)
    std_kw_small = safe_pstdev(kw_small_vals)

    return FamilySummary(
        name=name,
        n=n,
        mean_twin_count=mean(twin_vals),
        mean_quad_all=mean(quad_all_vals),
        mean_quad_small=mean(quad_small_vals),
        mean_exact_sym=mean(exact_vals),
        mean_near6=mean(near6_vals),
        mean_near12=mean(near12_vals),
        mean_kw_all=mean(kw_all_vals),
        mean_kw_small=mean(kw_small_vals),
        std_twin_count=std_twin,
        std_quad_all=std_quad_all,
        std_quad_small=std_quad_small,
        std_exact_sym=std_exact,
        std_near6=std_near6,
        std_near12=std_near12,
        std_kw_all=std_kw_all,
        std_kw_small=std_kw_small,
        se_twin_count=safe_se(std_twin, n),
        se_quad_all=safe_se(std_quad_all, n),
        se_quad_small=safe_se(std_quad_small, n),
        se_exact_sym=safe_se(std_exact, n),
        se_near6=safe_se(std_near6, n),
        se_near12=safe_se(std_near12, n),
        se_kw_all=safe_se(std_kw_all, n),
        se_kw_small=safe_se(std_kw_small, n),
    )


def z_like(diff: float, se_a: float, se_b: float) -> float:
    """
    Einfacher Vergleichswert:
    diff / sqrt(se_a^2 + se_b^2)

    Kein exakter Hypothesentest, aber ein guter Schnellindikator.
    """
    denom = sqrt(se_a * se_a + se_b * se_b)
    if denom == 0:
        return 0.0
    return diff / denom


def print_family_summary(summary: FamilySummary) -> None:
    print(f"\n=== Zusammenfassung: {summary.name} ===")
    print(f"n = {summary.n}")
    print(f"twin_count        : mean={summary.mean_twin_count:.6f}, se={summary.se_twin_count:.6f}")
    print(f"quad_all          : mean={summary.mean_quad_all:.6f}, se={summary.se_quad_all:.6f}")
    print(f"quad_small        : mean={summary.mean_quad_small:.6f}, se={summary.se_quad_small:.6f}")
    print(f"exact_sym         : mean={summary.mean_exact_sym:.6f}, se={summary.se_exact_sym:.6f}")
    print(f"near6             : mean={summary.mean_near6:.6f}, se={summary.se_near6:.6f}")
    print(f"near12            : mean={summary.mean_near12:.6f}, se={summary.se_near12:.6f}")
    print(f"kw_all            : mean={summary.mean_kw_all:.6f}, se={summary.se_kw_all:.6f}")
    print(f"kw_small          : mean={summary.mean_kw_small:.6f}, se={summary.se_kw_small:.6f}")


def print_histograms(name: str, stats: List[CenterStats], top_k: int = 12) -> None:
    delta_hist = Counter()
    eta_hist = Counter()
    for s in stats:
        delta_hist.update(s.delta_histogram_small)
        eta_hist.update(s.eta_histogram_small)

    print(f"\n--- Häufigste kleine δ-Werte [{name}] ---")
    for delta, count in delta_hist.most_common(top_k):
        print(f"δ = {delta:>3} : {count}")

    print(f"\n--- Häufigste kleine η-Werte [{name}] ---")
    for eta, count in eta_hist.most_common(top_k):
        print(f"η = {eta:>3} : {count}")


def print_top_centers(
    name: str,
    stats: List[CenterStats],
    top_k: int = 10,
) -> None:
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

    print(f"\nTop-{top_k} Zentren [{name}] nach kw_small:")
    for s in ranked[:top_k]:
        print(
            f"M={s.M:>8} | twins={s.twin_count:>3} | "
            f"Qall={s.quadruple_count_all:>3} | "
            f"Qsmall={s.quadruple_count_small_delta:>3} | "
            f"exact={s.exact_symmetric_count:>3} | "
            f"near6={s.near_symmetric_count_6:>3} | "
            f"near12={s.near_symmetric_count_12:>3} | "
            f"kw_small={s.weighted_cupola_small_delta:.6f}"
        )


# =========================
# Vergleichstabellen
# =========================

METRICS = [
    ("twin_count", "mean_twin_count", "se_twin_count"),
    ("quad_all", "mean_quad_all", "se_quad_all"),
    ("quad_small", "mean_quad_small", "se_quad_small"),
    ("exact_sym", "mean_exact_sym", "se_exact_sym"),
    ("near6", "mean_near6", "se_near6"),
    ("near12", "mean_near12", "se_near12"),
    ("kw_all", "mean_kw_all", "se_kw_all"),
    ("kw_small", "mean_kw_small", "se_kw_small"),
]


def print_difference_table(
    base: FamilySummary,
    other: FamilySummary,
) -> None:
    print(f"\n=== Differenztabelle: {base.name} minus {other.name} ===")
    print(f"{'Metrik':<12} {'Diff':>12} {'z-like':>12}")
    print("-" * 38)
    for label, mean_attr, se_attr in METRICS:
        diff = getattr(base, mean_attr) - getattr(other, mean_attr)
        zval = z_like(diff, getattr(base, se_attr), getattr(other, se_attr))
        print(f"{label:<12} {diff:>12.6f} {zval:>12.3f}")


# =========================
# CSV-Export
# =========================

def export_center_stats_csv(path: Path, stats: List[CenterStats]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "M",
            "twin_count",
            "quadruple_count_all",
            "quadruple_count_small_delta",
            "exact_symmetric_count",
            "near_symmetric_count_6",
            "near_symmetric_count_12",
            "weighted_cupola_all",
            "weighted_cupola_small_delta",
            "delta_histogram_small",
            "eta_histogram_small",
        ])
        for s in stats:
            writer.writerow([
                s.M,
                s.twin_count,
                s.quadruple_count_all,
                s.quadruple_count_small_delta,
                s.exact_symmetric_count,
                s.near_symmetric_count_6,
                s.near_symmetric_count_12,
                f"{s.weighted_cupola_all:.12f}",
                f"{s.weighted_cupola_small_delta:.12f}",
                repr(s.delta_histogram_small),
                repr(s.eta_histogram_small),
            ])


def export_family_summaries_csv(path: Path, summaries: List[FamilySummary]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(asdict(summaries[0]).keys()) if summaries else []
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for s in summaries:
            writer.writerow(asdict(s))


# =========================
# Hauptprogramm
# =========================

def main() -> None:
    # -------- Parameter --------
    min_center = 10_000
    max_center = 1_000_000
    radius = 300
    delta_max = 120

    export_dir = Path("zwillings_export_v3")
    do_csv_export = True

    # -------- Primzahlen / Zwillinge --------
    sieve_limit = max_center + radius + 10
    prime_flags = prime_sieve(sieve_limit)
    twins = collect_twin_primes(prime_flags)
    center_map = twins_by_center(twins)

    # -------- Familien --------
    centers_30 = build_centers_by_divisibility(min_center, max_center, 30)
    centers_210 = build_centers_by_divisibility(min_center, max_center, 210)
    centers_2310 = build_centers_by_divisibility(min_center, max_center, 2310)

    controls_1 = build_control_centers(min_center, max_center, residue=1, modulus=30)
    controls_7 = build_control_centers(min_center, max_center, residue=7, modulus=30)
    controls_11 = build_control_centers(min_center, max_center, residue=11, modulus=30)
    controls_13 = build_control_centers(min_center, max_center, residue=13, modulus=30)

    print(f"Parameter: min_center={min_center}, max_center={max_center}, radius={radius}, delta_max={delta_max}")
    print("Starte Analyse ...")

    # -------- Analyse --------
    stats_30 = analyze_centers(centers_30, radius, center_map, delta_max)
    stats_210 = analyze_centers(centers_210, radius, center_map, delta_max)
    stats_2310 = analyze_centers(centers_2310, radius, center_map, delta_max)

    stats_c1 = analyze_centers(controls_1, radius, center_map, delta_max)
    stats_c7 = analyze_centers(controls_7, radius, center_map, delta_max)
    stats_c11 = analyze_centers(controls_11, radius, center_map, delta_max)
    stats_c13 = analyze_centers(controls_13, radius, center_map, delta_max)

    # -------- Summaries --------
    sum_30 = summarize_family("30 | M", stats_30)
    sum_210 = summarize_family("210 | M", stats_210)
    sum_2310 = summarize_family("2310 | M", stats_2310)
    sum_c1 = summarize_family("Kontrolle M ≡ 1 mod 30", stats_c1)
    sum_c7 = summarize_family("Kontrolle M ≡ 7 mod 30", stats_c7)
    sum_c11 = summarize_family("Kontrolle M ≡ 11 mod 30", stats_c11)
    sum_c13 = summarize_family("Kontrolle M ≡ 13 mod 30", stats_c13)

    summaries = [sum_30, sum_210, sum_2310, sum_c1, sum_c7, sum_c11, sum_c13]

    # -------- Ausgabe --------
    for s in summaries:
        print_family_summary(s)

    print_histograms("30 | M", stats_30)
    print_histograms("210 | M", stats_210)
    print_histograms("2310 | M", stats_2310)

    print_top_centers("30 | M", stats_30, top_k=12)
    print_top_centers("210 | M", stats_210, top_k=12)
    print_top_centers("2310 | M", stats_2310, top_k=12)

    # Differenztabellen
    print_difference_table(sum_210, sum_30)
    print_difference_table(sum_2310, sum_30)
    print_difference_table(sum_30, sum_c1)
    print_difference_table(sum_30, sum_c7)
    print_difference_table(sum_30, sum_c11)
    print_difference_table(sum_30, sum_c13)
    print_difference_table(sum_210, sum_c1)
    print_difference_table(sum_210, sum_c7)
    print_difference_table(sum_210, sum_c11)
    print_difference_table(sum_210, sum_c13)

    # CSV-Export
    if do_csv_export:
        export_center_stats_csv(export_dir / "centers_30.csv", stats_30)
        export_center_stats_csv(export_dir / "centers_210.csv", stats_210)
        export_center_stats_csv(export_dir / "centers_2310.csv", stats_2310)
        export_center_stats_csv(export_dir / "controls_r1.csv", stats_c1)
        export_center_stats_csv(export_dir / "controls_r7.csv", stats_c7)
        export_center_stats_csv(export_dir / "controls_r11.csv", stats_c11)
        export_center_stats_csv(export_dir / "controls_r13.csv", stats_c13)

        export_family_summaries_csv(export_dir / "family_summaries.csv", summaries)
        print(f"\nCSV-Export geschrieben nach: {export_dir.resolve()}")


if __name__ == "__main__":
    main()