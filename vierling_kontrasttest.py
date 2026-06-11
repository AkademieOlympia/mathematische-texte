from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
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


def sign(x: float, eps: float = 1e-12) -> int:
    if x > eps:
        return 1
    if x < -eps:
        return -1
    return 0


def quadruplet_features(quad: Sequence[int], pivot: int) -> Dict[str, object]:
    n_minus_ea = n_minus_bc = n_plus_ea = n_plus_bc = 0
    fine_pattern: List[str] = []
    coarse_pattern: List[str] = []

    for q in quad:
        fc = fine_class(q)
        cc = COARSE_MAP[fc]
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
        "n_minus_ea": n_minus_ea,
        "n_minus_bc": n_minus_bc,
        "n_plus_ea": n_plus_ea,
        "n_plus_bc": n_plus_bc,
        "delta": delta,
        "sigma": sigma,
        "B_log": B_log,
        "beta5": beta5,
        "sgn_B": sign(B_log),
    }


def summarize_group(rows: Sequence[Dict[str, object]]) -> Dict[str, float]:
    if not rows:
        return {
            "N": 0,
            "mean_delta": 0.0,
            "mean_sigma": 0.0,
            "mean_B": 0.0,
            "mean_beta5": 0.0,
            "std_B": 0.0,
            "frac_delta_eq_0": 0.0,
            "frac_sigma_eq_0": 0.0,
            "frac_B_pos": 0.0,
            "frac_B_neg": 0.0,
        }

    Bs = [float(r["B_log"]) for r in rows]
    mean_B = mean(Bs)
    var_B = mean([(x - mean_B) ** 2 for x in Bs])
    std_B = math.sqrt(var_B)

    return {
        "N": len(rows),
        "mean_delta": mean(float(r["delta"]) for r in rows),
        "mean_sigma": mean(float(r["sigma"]) for r in rows),
        "mean_B": mean_B,
        "mean_beta5": mean(float(r["beta5"]) for r in rows),
        "std_B": std_B,
        "frac_delta_eq_0": sum(1 for r in rows if int(r["delta"]) == 0) / len(rows),
        "frac_sigma_eq_0": sum(1 for r in rows if int(r["sigma"]) == 0) / len(rows),
        "frac_B_pos": sum(1 for r in rows if int(r["sgn_B"]) == 1) / len(rows),
        "frac_B_neg": sum(1 for r in rows if int(r["sgn_B"]) == -1) / len(rows),
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


def parse_pivots(pivot_str: str) -> List[int]:
    return [int(x.strip()) for x in pivot_str.split(",") if x.strip()]


def print_summary(summary_rows: Sequence[Dict[str, object]]) -> None:
    print("=" * 154)
    print("KONTRASTTEST: ERGÄNZBARE VS. NICHT ERGÄNZBARE VIERLINGE")
    print("=" * 154)
    header = (
        f"{'π':>5} {'Gruppe':>14} {'N':>6} {'<Δ>':>8} {'<Σ>':>8} "
        f"{'<B>':>12} {'<β5>':>10} {'std(B)':>10} {'Δ=0':>8} {'Σ=0':>8} {'B>0':>8} {'B<0':>8}"
    )
    print(header)
    print("-" * len(header))
    for row in summary_rows:
        print(
            f"{row['pivot']:>5} {row['group']:>14} {int(row['N']):>6} "
            f"{row['mean_delta']:>8.3f} {row['mean_sigma']:>8.3f} "
            f"{row['mean_B']:>12.6f} {row['mean_beta5']:>10.6f} {row['std_B']:>10.6f} "
            f"{row['frac_delta_eq_0']*100:>7.2f}% {row['frac_sigma_eq_0']*100:>7.2f}% "
            f"{row['frac_B_pos']*100:>7.2f}% {row['frac_B_neg']*100:>7.2f}%"
        )


def print_examples(rows: Sequence[Dict[str, object]], title: str, limit: int = 20) -> None:
    print("\n" + "=" * 154)
    print(title)
    print("=" * 154)
    if not rows:
        print("Keine Daten.")
        return
    header = f"{'π':>5} {'Start':>8} {'Q4 grob':>15} {'Δ':>4} {'Σ':>4} {'B':>12} {'β5':>10} {'ergänzbar':>10}"
    print(header)
    print("-" * len(header))
    for r in rows[:limit]:
        print(
            f"{r['pivot']:>5} {r['start']:>8} {r['coarse_pattern']:>15} "
            f"{int(r['delta']):>4} {int(r['sigma']):>4} {float(r['B_log']):>12.6f} "
            f"{float(r['beta5']):>10.6f} {str(r['extendable']):>10}"
        )
    if len(rows) > limit:
        print(f"... weitere {len(rows) - limit} nicht angezeigt.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Vergleicht ergänzbare und nicht ergänzbare dichte Primzahlvierlinge relativ zu Pivot-Primzahlen."
    )
    parser.add_argument("--limit", type=int, default=1_000_000, help="obere Grenze")
    parser.add_argument("--pivots", type=str, default="11,13,17,19,29,37", help="kommagetrennte Pivot-Primzahlen")
    parser.add_argument("--show", type=int, default=20, help="wie viele Beispiele angezeigt werden")
    parser.add_argument("--csv", type=Path, default=None, help="Detaildaten als CSV")
    parser.add_argument("--summary-csv", type=Path, default=None, help="Zusammenfassung als CSV")
    args = parser.parse_args()

    pivots = parse_pivots(args.pivots)
    sieve = PrimeSieve(args.limit + 12)
    quads = dense_quadruplets(args.limit, sieve)

    detail_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []

    for pivot in pivots:
        ext_rows: List[Dict[str, object]] = []
        non_rows: List[Dict[str, object]] = []

        for quad in quads:
            row = quadruplet_features(quad, pivot)
            row["extendable"] = (quad[0] + 12) in sieve
            if row["extendable"]:
                ext_rows.append(row)
            else:
                non_rows.append(row)
            detail_rows.append(row)

        ext_sum = summarize_group(ext_rows)
        ext_sum["pivot"] = pivot
        ext_sum["group"] = "ergänzbar"
        summary_rows.append(ext_sum)

        non_sum = summarize_group(non_rows)
        non_sum["pivot"] = pivot
        non_sum["group"] = "nicht ergänzbar"
        summary_rows.append(non_sum)

    print_summary(summary_rows)

    # Beispiele: je Pivot zuerst ergänzbar, dann nicht ergänzbar
    ext_examples = [r for r in detail_rows if r["extendable"]]
    non_examples = [r for r in detail_rows if not r["extendable"]]
    print_examples(ext_examples, "BEISPIELE: ERGÄNZBARE VIERLINGE", limit=args.show)
    print_examples(non_examples, "BEISPIELE: NICHT ERGÄNZBARE VIERLINGE", limit=args.show)

    if args.csv is not None:
        write_csv(args.csv, detail_rows)
        print(f"\nDetail-CSV gespeichert: {args.csv}")
    if args.summary_csv is not None:
        write_csv(args.summary_csv, summary_rows)
        print(f"Zusammenfassungs-CSV gespeichert: {args.summary_csv}")


if __name__ == "__main__":
    main()