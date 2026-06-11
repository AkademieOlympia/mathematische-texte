from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from statistics import mean
from typing import Dict, List, Sequence, Tuple


FINE_CLASSES = {1: "E", 5: "A", 7: "B", 11: "C"}
COARSE_MAP = {"E": "EA", "A": "EA", "B": "BC", "C": "BC"}


class PrimeSieve:
    def __init__(self, limit: int) -> None:
        self.limit = limit
        self.is_prime = self._build(limit)
        self.primes = [i for i, flag in enumerate(self.is_prime) if flag]

    @staticmethod
    def _build(limit: int) -> List[bool]:
        sieve = [False, False] + [True] * (limit - 1)
        for p in range(2, int(limit**0.5) + 1):
            if sieve[p]:
                start = p * p
                sieve[start : limit + 1 : p] = [False] * (((limit - start) // p) + 1)
        return sieve

    def __contains__(self, n: int) -> bool:
        return 0 <= n <= self.limit and self.is_prime[n]


def dense_quadruplets(limit: int, sieve: PrimeSieve) -> List[Tuple[int, int, int, int]]:
    out: List[Tuple[int, int, int, int]] = []
    for p in sieve.primes:
        if p <= 3:
            continue
        if p + 8 > limit:
            break
        quad = (p, p + 2, p + 6, p + 8)
        if all(q in sieve for q in quad):
            out.append(quad)
    return out


def coarse_class(p: int) -> str:
    r = p % 12
    if r not in FINE_CLASSES:
        raise ValueError(f"Primzahl {p} liegt nicht in 1,5,7,11 mod 12")
    return COARSE_MAP[FINE_CLASSES[r]]


def quadruplet_features(quad: Sequence[int], pivot: int, sieve: PrimeSieve) -> Dict[str, object]:
    n_minus_ea = n_minus_bc = n_plus_ea = n_plus_bc = 0
    coarse_pattern: List[str] = []

    for q in quad:
        cc = coarse_class(q)
        coarse_pattern.append(cc)

        if q < pivot:
            if cc == "EA":
                n_minus_ea += 1
            else:
                n_minus_bc += 1
        elif q > pivot:
            if cc == "EA":
                n_plus_ea += 1
            else:
                n_plus_bc += 1

    delta = (n_minus_ea + n_minus_bc) - (n_plus_ea + n_plus_bc)
    sigma = (n_minus_ea - n_minus_bc) - (n_plus_ea - n_plus_bc)
    B_log = sum(math.log(q) for q in quad if q < pivot) - sum(math.log(q) for q in quad if q > pivot)
    beta5 = B_log / math.log(5)

    return {
        "pivot": pivot,
        "start": quad[0],
        "quadruplet": quad,
        "coarse_pattern": "-".join(coarse_pattern),
        "delta": delta,
        "sigma": sigma,
        "B_log": B_log,
        "beta5": beta5,
        "extendable": (quad[0] + 12) in sieve,
    }


def quantile_groups(rows: Sequence[Dict[str, object]], num_groups: int) -> List[List[Dict[str, object]]]:
    sorted_rows = sorted(rows, key=lambda r: float(r["B_log"]))
    n = len(sorted_rows)
    groups: List[List[Dict[str, object]]] = []
    for i in range(num_groups):
        a = (i * n) // num_groups
        b = ((i + 1) * n) // num_groups
        groups.append(sorted_rows[a:b])
    return groups


def summarize_quantile(rows: Sequence[Dict[str, object]], pivot: int, label: str) -> Dict[str, object]:
    if not rows:
        return {
            "pivot": pivot,
            "quantile": label,
            "N": 0,
            "B_min": None,
            "B_max": None,
            "B_mean": None,
            "beta5_mean": None,
            "extendable_rate": None,
            "delta_mean": None,
            "sigma_mean": None,
        }

    return {
        "pivot": pivot,
        "quantile": label,
        "N": len(rows),
        "B_min": min(float(r["B_log"]) for r in rows),
        "B_max": max(float(r["B_log"]) for r in rows),
        "B_mean": mean(float(r["B_log"]) for r in rows),
        "beta5_mean": mean(float(r["beta5"]) for r in rows),
        "extendable_rate": sum(1 for r in rows if bool(r["extendable"])) / len(rows),
        "delta_mean": mean(float(r["delta"]) for r in rows),
        "sigma_mean": mean(float(r["sigma"]) for r in rows),
    }


def summarize_threshold(rows: Sequence[Dict[str, object]], pivot: int, threshold: float) -> Dict[str, object]:
    sub = [r for r in rows if float(r["B_log"]) >= threshold]
    if not sub:
        return {
            "pivot": pivot,
            "threshold": threshold,
            "N": 0,
            "extendable_rate": None,
            "B_mean": None,
            "beta5_mean": None,
        }

    return {
        "pivot": pivot,
        "threshold": threshold,
        "N": len(sub),
        "extendable_rate": sum(1 for r in sub if bool(r["extendable"])) / len(sub),
        "B_mean": mean(float(r["B_log"]) for r in sub),
        "beta5_mean": mean(float(r["beta5"]) for r in sub),
    }


def write_csv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_pivots(s: str) -> List[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def print_quantile_summary(rows: Sequence[Dict[str, object]], quantiles: int) -> None:
    print("=" * 148)
    print("MEHRPIVOT-QUANTILTEST NACH B_pi(Q4)")
    print("=" * 148)
    header = (
        f"{'π':>5} {'Q':>4} {'N':>6} {'B_min':>12} {'B_max':>12} "
        f"{'<B>':>12} {'<β5>':>12} {'Erg.-Rate':>12} {'<Δ>':>8} {'<Σ>':>8}"
    )
    print(header)
    print("-" * len(header))
    for row in rows:
        if row["N"] == 0:
            continue
        print(
            f"{int(row['pivot']):>5} {str(row['quantile']):>4} {int(row['N']):>6} "
            f"{float(row['B_min']):>12.6f} {float(row['B_max']):>12.6f} "
            f"{float(row['B_mean']):>12.6f} {float(row['beta5_mean']):>12.6f} "
            f"{float(row['extendable_rate'])*100:>11.2f}% "
            f"{float(row['delta_mean']):>8.3f} {float(row['sigma_mean']):>8.3f}"
        )


def print_top_quantile_overview(rows: Sequence[Dict[str, object]], quantiles: int) -> None:
    print("\n" + "=" * 88)
    print("OBERSTES QUANTIL JE PIVOT")
    print("=" * 88)
    header = f"{'π':>5} {'Q_top':>8} {'N':>6} {'<B>':>12} {'Erg.-Rate':>12} {'Grundrate':>12}"
    print(header)
    print("-" * len(header))

    pivots = sorted(set(int(r["pivot"]) for r in rows))
    for pivot in pivots:
        sub = [r for r in rows if int(r["pivot"]) == pivot]
        top = [r for r in sub if str(r["quantile"]) == f"Q{quantiles}"][0]
        total_n = sum(int(r["N"]) for r in sub)
        total_hits = sum(float(r["extendable_rate"]) * int(r["N"]) for r in sub)
        base_rate = total_hits / total_n if total_n else 0.0

        print(
            f"{pivot:>5} {str(top['quantile']):>8} {int(top['N']):>6} "
            f"{float(top['B_mean']):>12.6f} {float(top['extendable_rate'])*100:>11.2f}% "
            f"{base_rate*100:>11.2f}%"
        )


def print_threshold_summary(rows: Sequence[Dict[str, object]]) -> None:
    print("\n" + "=" * 96)
    print("MEHRPIVOT-SCHWELLENTEST")
    print("=" * 96)
    header = f"{'π':>5} {'c':>12} {'N':>8} {'Erg.-Rate':>12} {'<B>':>12} {'<β5>':>12}"
    print(header)
    print("-" * len(header))
    for row in rows:
        if row["N"] == 0:
            print(f"{int(row['pivot']):>5} {float(row['threshold']):>12.6f} {0:>8} {'-':>12} {'-':>12} {'-':>12}")
        else:
            print(
                f"{int(row['pivot']):>5} {float(row['threshold']):>12.6f} {int(row['N']):>8} "
                f"{float(row['extendable_rate'])*100:>11.2f}% "
                f"{float(row['B_mean']):>12.6f} {float(row['beta5_mean']):>12.6f}"
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Vergleicht Quantil- und Schwellentests für B_pi(Q4) über mehrere Pivot-Primzahlen."
    )
    parser.add_argument("--limit", type=int, default=1_000_000, help="obere Grenze")
    parser.add_argument("--pivots", type=str, default="11,13,17,19,29,37", help="kommagetrennte Pivot-Primzahlen")
    parser.add_argument("--quantiles", type=int, default=5, help="Anzahl der Quantile")
    parser.add_argument(
        "--threshold-mode",
        type=str,
        default="auto",
        choices=["auto", "none"],
        help="automatische Schwellen aus Quantilgrenzen oder keine",
    )
    parser.add_argument("--quantile-csv", type=Path, default=None, help="CSV für Quantilergebnisse")
    parser.add_argument("--threshold-csv", type=Path, default=None, help="CSV für Schwellenergebnisse")
    args = parser.parse_args()

    pivots = parse_pivots(args.pivots)
    sieve = PrimeSieve(args.limit + 12)
    quads = dense_quadruplets(args.limit, sieve)

    quantile_rows: List[Dict[str, object]] = []
    threshold_rows: List[Dict[str, object]] = []

    for pivot in pivots:
        rows = [quadruplet_features(q, pivot, sieve) for q in quads]
        groups = quantile_groups(rows, args.quantiles)

        for i, g in enumerate(groups):
            quantile_rows.append(summarize_quantile(g, pivot=pivot, label=f"Q{i+1}"))

        if args.threshold_mode == "auto":
            Bs = sorted(float(r["B_log"]) for r in rows)
            thresholds = [
                Bs[int(0.50 * len(Bs))],
                Bs[int(0.70 * len(Bs))],
                Bs[int(0.85 * len(Bs))],
                Bs[int(0.95 * len(Bs))],
            ]
            # doppelte Schwellen entfernen, Reihenfolge bewahren
            seen = set()
            uniq_thresholds = []
            for t in thresholds:
                key = round(t, 12)
                if key not in seen:
                    seen.add(key)
                    uniq_thresholds.append(t)
            for t in uniq_thresholds:
                threshold_rows.append(summarize_threshold(rows, pivot=pivot, threshold=t))

    print_quantile_summary(quantile_rows, args.quantiles)
    print_top_quantile_overview(quantile_rows, args.quantiles)
    if threshold_rows:
        print_threshold_summary(threshold_rows)

    if args.quantile_csv is not None:
        write_csv(args.quantile_csv, quantile_rows)
        print(f"\nQuantil-CSV gespeichert: {args.quantile_csv}")
    if args.threshold_csv is not None and threshold_rows:
        write_csv(args.threshold_csv, threshold_rows)
        print(f"Schwellen-CSV gespeichert: {args.threshold_csv}")


if __name__ == "__main__":
    main()