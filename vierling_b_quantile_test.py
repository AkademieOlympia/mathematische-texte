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


def fine_class(p: int) -> str:
    r = p % 12
    if r not in FINE_CLASSES:
        raise ValueError(f"Primzahl {p} liegt nicht in 1,5,7,11 mod 12")
    return FINE_CLASSES[r]


def coarse_class(p: int) -> str:
    return COARSE_MAP[fine_class(p)]


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


def quadruplet_features(quad: Sequence[int], pivot: int, sieve: PrimeSieve) -> Dict[str, object]:
    n_minus_ea = n_minus_bc = n_plus_ea = n_plus_bc = 0
    fine_pattern: List[str] = []
    coarse_pattern: List[str] = []

    for q in quad:
        fc = fine_class(q)
        cc = coarse_class(q)
        fine_pattern.append(fc)
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
        "fine_pattern": "-".join(fine_pattern),
        "coarse_pattern": "-".join(coarse_pattern),
        "delta": delta,
        "sigma": sigma,
        "B_log": B_log,
        "beta5": beta5,
        "extendable": (quad[0] + 12) in sieve,
    }


def quantile_groups(rows: Sequence[Dict[str, object]], num_groups: int) -> List[List[Dict[str, object]]]:
    if num_groups < 1:
        raise ValueError("num_groups muss >= 1 sein")
    sorted_rows = sorted(rows, key=lambda r: float(r["B_log"]))
    n = len(sorted_rows)
    groups: List[List[Dict[str, object]]] = []
    for i in range(num_groups):
        a = (i * n) // num_groups
        b = ((i + 1) * n) // num_groups
        groups.append(sorted_rows[a:b])
    return groups


def summarize_group(rows: Sequence[Dict[str, object]], label: str) -> Dict[str, object]:
    if not rows:
        return {
            "label": label,
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
        "label": label,
        "N": len(rows),
        "B_min": min(float(r["B_log"]) for r in rows),
        "B_max": max(float(r["B_log"]) for r in rows),
        "B_mean": mean(float(r["B_log"]) for r in rows),
        "beta5_mean": mean(float(r["beta5"]) for r in rows),
        "extendable_rate": sum(1 for r in rows if bool(r["extendable"])) / len(rows),
        "delta_mean": mean(float(r["delta"]) for r in rows),
        "sigma_mean": mean(float(r["sigma"]) for r in rows),
    }


def threshold_summary(rows: Sequence[Dict[str, object]], thresholds: Sequence[float]) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    for c in thresholds:
        sub = [r for r in rows if float(r["B_log"]) >= c]
        if sub:
            out.append(
                {
                    "threshold": c,
                    "N": len(sub),
                    "extendable_rate": sum(1 for r in sub if bool(r["extendable"])) / len(sub),
                    "B_mean": mean(float(r["B_log"]) for r in sub),
                    "beta5_mean": mean(float(r["beta5"]) for r in sub),
                }
            )
        else:
            out.append(
                {
                    "threshold": c,
                    "N": 0,
                    "extendable_rate": None,
                    "B_mean": None,
                    "beta5_mean": None,
                }
            )
    return out


def write_csv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def print_quantiles(rows: Sequence[Dict[str, object]]) -> None:
    print("=" * 120)
    print("QUANTILTEST NACH B_pi(Q4)")
    print("=" * 120)
    header = f"{'Quantil':>8} {'N':>6} {'B_min':>12} {'B_max':>12} {'<B>':>12} {'<β5>':>12} {'Erg.-Rate':>12} {'<Δ>':>8} {'<Σ>':>8}"
    print(header)
    print("-" * len(header))
    for row in rows:
        if row["N"] == 0:
            print(f"{row['label']:>8} {0:>6}")
        else:
            print(
                f"{row['label']:>8} {int(row['N']):>6} "
                f"{float(row['B_min']):>12.6f} {float(row['B_max']):>12.6f} "
                f"{float(row['B_mean']):>12.6f} {float(row['beta5_mean']):>12.6f} "
                f"{float(row['extendable_rate'])*100:>11.2f}% "
                f"{float(row['delta_mean']):>8.3f} {float(row['sigma_mean']):>8.3f}"
            )


def print_thresholds(rows: Sequence[Dict[str, object]]) -> None:
    print("\n" + "=" * 96)
    print("SCHWELLENTEST: B_pi(Q4) >= c")
    print("=" * 96)
    header = f"{'c':>12} {'N':>8} {'Erg.-Rate':>12} {'<B>':>12} {'<β5>':>12}"
    print(header)
    print("-" * len(header))
    for row in rows:
        if row["N"] == 0:
            print(f"{float(row['threshold']):>12.6f} {0:>8} {'-':>12} {'-':>12} {'-':>12}")
        else:
            print(
                f"{float(row['threshold']):>12.6f} {int(row['N']):>8} "
                f"{float(row['extendable_rate'])*100:>11.2f}% "
                f"{float(row['B_mean']):>12.6f} {float(row['beta5_mean']):>12.6f}"
            )


def print_examples(rows: Sequence[Dict[str, object]], title: str, limit: int = 20) -> None:
    print("\n" + "=" * 132)
    print(title)
    print("=" * 132)
    if not rows:
        print("Keine Daten.")
        return
    header = f"{'Start':>8} {'Q4 grob':>15} {'B':>12} {'β5':>12} {'Δ':>4} {'Σ':>4} {'ergänzbar':>10}"
    print(header)
    print("-" * len(header))
    for r in rows[:limit]:
        print(
            f"{int(r['start']):>8} {str(r['coarse_pattern']):>15} "
            f"{float(r['B_log']):>12.6f} {float(r['beta5']):>12.6f} "
            f"{int(r['delta']):>4} {int(r['sigma']):>4} {str(r['extendable']):>10}"
        )
    if len(rows) > limit:
        print(f"... weitere {len(rows) - limit} nicht angezeigt.")


def parse_thresholds(s: str) -> List[float]:
    if not s.strip():
        return []
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Testet, ob B_pi(Q4) die Ergänzbarkeit dichter Vierlinge zu Fünflingen vorhersagt."
    )
    parser.add_argument("--limit", type=int, default=1_000_000, help="obere Grenze")
    parser.add_argument("--pivot", type=int, default=11, help="feste Pivot-Primzahl")
    parser.add_argument("--quantiles", type=int, default=5, help="Anzahl Quantile")
    parser.add_argument(
        "--thresholds",
        type=str,
        default="",
        help="kommagetrennte Schwellen c für Test B_pi(Q4) >= c",
    )
    parser.add_argument("--show", type=int, default=20, help="wie viele Beispiele gezeigt werden")
    parser.add_argument("--csv", type=Path, default=None, help="Detaildaten als CSV")
    parser.add_argument("--quantile-csv", type=Path, default=None, help="Quantil-Zusammenfassung als CSV")
    parser.add_argument("--threshold-csv", type=Path, default=None, help="Schwellen-Zusammenfassung als CSV")
    args = parser.parse_args()

    sieve = PrimeSieve(args.limit + 12)
    quads = dense_quadruplets(args.limit, sieve)
    rows = [quadruplet_features(q, args.pivot, sieve) for q in quads]

    q_groups = quantile_groups(rows, args.quantiles)
    q_summary = [summarize_group(g, label=f"Q{i+1}") for i, g in enumerate(q_groups)]

    thresholds = parse_thresholds(args.thresholds)
    if not thresholds:
        Bs = sorted(float(r["B_log"]) for r in rows)
        # automatische sinnvolle Schwellen aus dem oberen Bereich
        thresholds = [
            Bs[int(0.50 * len(Bs))],
            Bs[int(0.70 * len(Bs))],
            Bs[int(0.85 * len(Bs))],
            Bs[int(0.95 * len(Bs))],
        ]
    t_summary = threshold_summary(rows, thresholds)

    print(f"Pivot π = {args.pivot}")
    print(f"Anzahl dichter Vierlinge: {len(rows)}")
    print(f"Gesamte Ergänzungsrate: {100 * sum(1 for r in rows if r['extendable']) / len(rows):.2f}%\n")

    print_quantiles(q_summary)
    print_thresholds(t_summary)

    rows_sorted_high_B = sorted(rows, key=lambda r: float(r["B_log"]), reverse=True)
    print_examples(rows_sorted_high_B, "VIERLINGE MIT HOHEM B_pi(Q4)", limit=args.show)

    if args.csv is not None:
        write_csv(args.csv, rows)
        print(f"\nDetail-CSV gespeichert: {args.csv}")
    if args.quantile_csv is not None:
        write_csv(args.quantile_csv, q_summary)
        print(f"Quantil-CSV gespeichert: {args.quantile_csv}")
    if args.threshold_csv is not None:
        write_csv(args.threshold_csv, t_summary)
        print(f"Schwellen-CSV gespeichert: {args.threshold_csv}")


if __name__ == "__main__":
    main()