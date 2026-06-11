from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

FINE_CLASSES = {1: "E", 5: "A", 7: "B", 11: "C"}
COARSE_MAP = {"E": "EA", "A": "EA", "B": "BC", "C": "BC"}


@dataclass
class TupleBalance:
    pivot: int
    tuple_type: str
    start: int
    values: Tuple[int, ...]
    fine_pattern: Tuple[str, ...]
    coarse_pattern: Tuple[str, ...]
    n_minus_ea: int
    n_minus_bc: int
    n_plus_ea: int
    n_plus_bc: int
    delta: int
    sigma: int
    B_log: float
    beta5: float


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
                step = p
                start = p * p
                sieve[start : limit + 1 : step] = [False] * (((limit - start) // step) + 1)
        return sieve

    def __contains__(self, n: int) -> bool:
        return 0 <= n <= self.limit and self.is_prime[n]


def fine_class(p: int) -> str:
    r = p % 12
    if r not in FINE_CLASSES:
        raise ValueError(f"Primzahl {p} liegt nicht in einer erlaubten Klasse mod 12")
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


def is_dense_quintuplet(quad: Tuple[int, int, int, int], sieve: PrimeSieve) -> bool:
    return quad[0] + 12 in sieve


def tuple_balance(values: Sequence[int], pivot: int, tuple_type: str) -> TupleBalance:
    n_minus_ea = n_minus_bc = n_plus_ea = n_plus_bc = 0
    fine_pattern: List[str] = []
    coarse_pattern: List[str] = []

    for q in values:
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
    B_log = sum(math.log(q) for q in values if q < pivot) - sum(math.log(q) for q in values if q > pivot)
    beta5 = B_log / math.log(5)

    return TupleBalance(
        pivot=pivot,
        tuple_type=tuple_type,
        start=values[0],
        values=tuple(values),
        fine_pattern=tuple(fine_pattern),
        coarse_pattern=tuple(coarse_pattern),
        n_minus_ea=n_minus_ea,
        n_minus_bc=n_minus_bc,
        n_plus_ea=n_plus_ea,
        n_plus_bc=n_plus_bc,
        delta=delta,
        sigma=sigma,
        B_log=B_log,
        beta5=beta5,
    )


def sign(x: float, eps: float = 1e-12) -> int:
    if x > eps:
        return 1
    if x < -eps:
        return -1
    return 0


def compare_quad_quint(quad: Tuple[int, int, int, int], pivot: int) -> Dict[str, object]:
    quint = quad + (quad[0] + 12,)
    qb = tuple_balance(quad, pivot, "Q4")
    fb = tuple_balance(quint, pivot, "Q5")

    return {
        "pivot": pivot,
        "start": quad[0],
        "quadruplet": quad,
        "quintuplet": quint,
        "quad_fine": "-".join(qb.fine_pattern),
        "quad_coarse": "-".join(qb.coarse_pattern),
        "quint_fine": "-".join(fb.fine_pattern),
        "quint_coarse": "-".join(fb.coarse_pattern),
        "delta_q4": qb.delta,
        "delta_q5": fb.delta,
        "delta_jump": fb.delta - qb.delta,
        "sigma_q4": qb.sigma,
        "sigma_q5": fb.sigma,
        "sigma_jump": fb.sigma - qb.sigma,
        "B_q4": qb.B_log,
        "B_q5": fb.B_log,
        "B_jump": fb.B_log - qb.B_log,
        "beta5_q4": qb.beta5,
        "beta5_q5": fb.beta5,
        "beta5_jump": fb.beta5 - qb.beta5,
        "sgn_B_q4": sign(qb.B_log),
        "sgn_B_q5": sign(fb.B_log),
        "B_sign_flip": sign(qb.B_log) != sign(fb.B_log),
        "q5_less_than_pivot": quint[-1] < pivot,
        "q5_greater_than_pivot": quint[-1] > pivot,
    }


def summarize(rows: Sequence[Dict[str, object]], all_quads: int, completed_quads: int, pivot: int) -> Dict[str, object]:
    n = len(rows)
    if n == 0:
        return {
            "pivot": pivot,
            "all_quadruplets": all_quads,
            "completed_quintuplets": completed_quads,
            "completion_rate": completed_quads / all_quads if all_quads else 0.0,
            "tested_extensions": 0,
            "B_sign_flip_rate": 0.0,
            "mean_delta_jump": 0.0,
            "mean_sigma_jump": 0.0,
            "mean_B_jump": 0.0,
            "mean_beta5_jump": 0.0,
            "q5_above_pivot_rate": 0.0,
            "q5_below_pivot_rate": 0.0,
        }

    return {
        "pivot": pivot,
        "all_quadruplets": all_quads,
        "completed_quintuplets": completed_quads,
        "completion_rate": completed_quads / all_quads if all_quads else 0.0,
        "tested_extensions": n,
        "B_sign_flip_rate": sum(1 for r in rows if r["B_sign_flip"]) / n,
        "mean_delta_jump": sum(float(r["delta_jump"]) for r in rows) / n,
        "mean_sigma_jump": sum(float(r["sigma_jump"]) for r in rows) / n,
        "mean_B_jump": sum(float(r["B_jump"]) for r in rows) / n,
        "mean_beta5_jump": sum(float(r["beta5_jump"]) for r in rows) / n,
        "q5_above_pivot_rate": sum(1 for r in rows if r["q5_greater_than_pivot"]) / n,
        "q5_below_pivot_rate": sum(1 for r in rows if r["q5_less_than_pivot"]) / n,
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


def print_summary_table(summary_rows: Sequence[Dict[str, object]]) -> None:
    print("=" * 140)
    print("VIERLING -> FÜNFLING: BALANCE- UND 5-SPRUNG-TEST")
    print("=" * 140)
    header = (
        f"{'Pivot':>7} {'Q4':>6} {'Q5':>6} {'Rate':>8} {'Tests':>7} {'B-Flip':>9} "
        f"{'<Δ>':>8} {'<Σ>':>8} {'<ΔB>':>12} {'<Δβ5>':>10} {'Q5>π':>8}"
    )
    print(header)
    print("-" * len(header))
    for row in summary_rows:
        print(
            f"{row['pivot']:>7} {row['all_quadruplets']:>6} {row['completed_quintuplets']:>6} "
            f"{row['completion_rate']*100:>7.2f}% {row['tested_extensions']:>7} "
            f"{row['B_sign_flip_rate']*100:>8.2f}% {row['mean_delta_jump']:>8.3f} "
            f"{row['mean_sigma_jump']:>8.3f} {row['mean_B_jump']:>12.6f} "
            f"{row['mean_beta5_jump']:>10.6f} {row['q5_above_pivot_rate']*100:>7.2f}%"
        )


def print_examples(rows: Sequence[Dict[str, object]], limit: int = 20) -> None:
    print("\n" + "=" * 140)
    print("AUSGEWÄHLTE ERWEITERUNGEN Q4 -> Q5")
    print("=" * 140)
    if not rows:
        print("Keine Erweiterungen gefunden.")
        return
    header = (
        f"{'π':>5} {'Start':>7} {'Q4 grob':>15} {'Q5 grob':>18} {'Δ':>4} {'Σ':>4} {'ΔB':>11} {'Flip':>6}"
    )
    print(header)
    print("-" * len(header))
    for row in rows[:limit]:
        print(
            f"{row['pivot']:>5} {row['start']:>7} {row['quad_coarse']:>15} {row['quint_coarse']:>18} "
            f"{int(row['delta_jump']):>4} {int(row['sigma_jump']):>4} {float(row['B_jump']):>11.6f} {str(row['B_sign_flip']):>6}"
        )
    if len(rows) > limit:
        print(f"... weitere {len(rows) - limit} nicht angezeigt.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Testet, ob Vierling->Fünfling-Erweiterungen mit Balancewechseln relativ zu Pivot-Primzahlen korrelieren."
    )
    parser.add_argument("--limit", type=int, default=1_000_000, help="obere Grenze für Start-p des Musters")
    parser.add_argument("--pivots", type=str, default="11,13,17,19,29,37", help="kommagetrennte Pivot-Primzahlen")
    parser.add_argument("--show", type=int, default=20, help="wie viele Beispielzeilen gezeigt werden")
    parser.add_argument("--csv", type=Path, default=None, help="CSV-Datei für Detaildaten")
    parser.add_argument("--summary-csv", type=Path, default=None, help="CSV-Datei für Zusammenfassung")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    pivots = parse_pivots(args.pivots)

    sieve = PrimeSieve(args.limit + 12)
    quads = dense_quadruplets(args.limit, sieve)
    completed = [q for q in quads if is_dense_quintuplet(q, sieve)]

    detail_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []

    for pivot in pivots:
        rows = [compare_quad_quint(q, pivot) for q in completed]
        for row in rows:
            row["pivot"] = pivot
        detail_rows.extend(rows)
        summary_rows.append(summarize(rows, all_quads=len(quads), completed_quads=len(completed), pivot=pivot))

    print_summary_table(summary_rows)
    print_examples(detail_rows, limit=args.show)

    if args.csv is not None:
        write_csv(args.csv, detail_rows)
        print(f"\nDetail-CSV gespeichert: {args.csv}")
    if args.summary_csv is not None:
        write_csv(args.summary_csv, summary_rows)
        print(f"Zusammenfassungs-CSV gespeichert: {args.summary_csv}")


if __name__ == "__main__":
    main()
