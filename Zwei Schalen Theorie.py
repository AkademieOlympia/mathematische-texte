from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


FINE_CLASSES = {1: "E", 5: "A", 7: "B", 11: "C"}
COARSE_MAP = {"E": "EA", "A": "EA", "B": "BC", "C": "BC"}
FINE_ORDER = ["E", "A", "B", "C"]
COARSE_ORDER = ["EA", "BC"]


class SPFTable:
    def __init__(self, limit: int) -> None:
        if limit < 2:
            raise ValueError("limit muss mindestens 2 sein")
        self.limit = limit
        self.spf = self._build(limit)

    @staticmethod
    def _build(limit: int) -> List[int]:
        spf = list(range(limit + 1))
        spf[0] = 0
        spf[1] = 1
        for p in range(2, int(limit**0.5) + 1):
            if spf[p] == p:
                for m in range(p * p, limit + 1, p):
                    if spf[m] == m:
                        spf[m] = p
        return spf

    def factor(self, n: int) -> Dict[int, int]:
        fac: Dict[int, int] = {}
        while n > 1:
            p = self.spf[n]
            fac[p] = fac.get(p, 0) + 1
            n //= p
        return fac


def rough_core_from_factors(factors: Dict[int, int]) -> int:
    m = 1
    for p, e in factors.items():
        if p not in (2, 3, 5):
            m *= p**e
    return m


def fine_class_of_rough(m: int) -> Optional[str]:
    if m == 1:
        return None
    r = m % 12
    return FINE_CLASSES.get(r)


def coarse_class_of_fine(fine: Optional[str]) -> Optional[str]:
    if fine is None:
        return None
    return COARSE_MAP[fine]


def entropy_from_counter(counter: Counter[str]) -> float:
    total = sum(counter.values())
    if total == 0:
        return 0.0
    H = 0.0
    for c in counter.values():
        p = c / total
        H -= p * math.log(p)
    return H


def primes_up_to(limit: int) -> List[int]:
    if limit < 2:
        return []
    sieve = [True] * (limit + 1)
    sieve[0:2] = [False, False]
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            for m in range(p * p, limit + 1, p):
                sieve[m] = False
    return [p for p in range(2, limit + 1) if sieve[p]]


def find_prime_quadruplets(limit: int) -> List[Tuple[int, int, int, int]]:
    ps = set(primes_up_to(limit + 8))
    out: List[Tuple[int, int, int, int]] = []
    for p in sorted(ps):
        if p <= 3:
            continue
        quad = (p, p + 2, p + 6, p + 8)
        if all(x in ps for x in quad):
            out.append(quad)
    return out


def analyze(limit: int) -> Dict[str, object]:
    spf = SPFTable(limit)

    fine_counter: Counter[str] = Counter()
    coarse_counter: Counter[str] = Counter()

    fine_transition: Dict[Tuple[str, str], int] = defaultdict(int)
    coarse_transition: Dict[Tuple[str, str], int] = defaultdict(int)

    rows: List[Dict[str, object]] = []

    prev_fine: Optional[str] = None
    prev_coarse: Optional[str] = None

    for n in range(1, limit + 1):
        factors = spf.factor(n)
        m = rough_core_from_factors(factors)

        fine = fine_class_of_rough(m)
        coarse = coarse_class_of_fine(fine)

        if fine is not None:
            fine_counter[fine] += 1
        if coarse is not None:
            coarse_counter[coarse] += 1

        if prev_fine is not None and fine is not None:
            fine_transition[(prev_fine, fine)] += 1
        if prev_coarse is not None and coarse is not None:
            coarse_transition[(prev_coarse, coarse)] += 1

        if fine is not None:
            prev_fine = fine
        if coarse is not None:
            prev_coarse = coarse

        rows.append(
            {
                "n": n,
                "rough_core": m,
                "fine_class": fine,
                "coarse_class": coarse,
            }
        )

    H4 = entropy_from_counter(fine_counter)
    H2 = entropy_from_counter(coarse_counter)

    quadruplets = find_prime_quadruplets(limit)
    quad_rows: List[Dict[str, object]] = []
    correct_projection_count = 0

    for quad in quadruplets:
        fines = [fine_class_of_rough(q) for q in quad]
        coarses = [coarse_class_of_fine(f) for f in fines]

        fine_pattern = tuple(fines)
        coarse_pattern = tuple(coarses)

        # Erwartete Grobprojektion eines Vierlings:
        # feine Struktur zyklisch aus (A,B,C,E), grob dann (EA,BC,BC,EA)
        # bis auf zyklische Verschiebung
        target = ("EA", "BC", "BC", "EA")
        rotations = [
            target,
            (target[1], target[2], target[3], target[0]),
            (target[2], target[3], target[0], target[1]),
            (target[3], target[0], target[1], target[2]),
        ]
        projection_ok = coarse_pattern in rotations
        if projection_ok:
            correct_projection_count += 1

        quad_rows.append(
            {
                "quadruplet": quad,
                "fine_pattern": fine_pattern,
                "coarse_pattern": coarse_pattern,
                "projection_ok": projection_ok,
            }
        )

    return {
        "rows": rows,
        "fine_counter": fine_counter,
        "coarse_counter": coarse_counter,
        "fine_transition": fine_transition,
        "coarse_transition": coarse_transition,
        "H4": H4,
        "H2": H2,
        "delta_H": H4 - H2,
        "quadruplet_rows": quad_rows,
        "quadruplet_count": len(quadruplets),
        "quadruplet_projection_ok_count": correct_projection_count,
    }


def print_counter(counter: Counter[str], order: List[str], title: str) -> None:
    total = sum(counter.values())
    print(title)
    for key in order:
        value = counter.get(key, 0)
        pct = 100.0 * value / total if total else 0.0
        print(f"  {key:>3}: {value:>10}   ({pct:6.2f}%)")
    print(f"  Summe: {total}")


def print_transition_matrix(
    transitions: Dict[Tuple[str, str], int],
    order: List[str],
    title: str,
) -> None:
    print(title)
    header = "       " + "".join(f"{c:>10}" for c in order)
    print(header)
    print("-" * len(header))
    for a in order:
        row = f"{a:>5}  "
        for b in order:
            row += f"{transitions.get((a, b), 0):>10}"
        print(row)


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Kopplungstest zwischen feiner EABC-Viererschale und grober EA/BC-Zweierschale."
    )
    parser.add_argument("--limit", type=int, default=10**6, help="oberes n")
    parser.add_argument("--csv", type=Path, default=None, help="CSV-Datei für n -> rough_core -> Klassen")
    parser.add_argument("--quad-csv", type=Path, default=None, help="CSV-Datei für Primzahlvierlinge")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    result = analyze(args.limit)

    print("=" * 100)
    print(f"KOPPLUNGSTEST BIS n = {args.limit}")
    print("=" * 100)

    print_counter(result["fine_counter"], FINE_ORDER, "\nFeine Klassenverteilung (E,A,B,C):")
    print()
    print_counter(result["coarse_counter"], COARSE_ORDER, "Grobe Klassenverteilung (EA,BC):")

    print()
    print_transition_matrix(result["fine_transition"], FINE_ORDER, "Feine Übergangsmatrix T4:")
    print()
    print_transition_matrix(result["coarse_transition"], COARSE_ORDER, "Grobe Übergangsmatrix T2:")

    print("\nEntropien:")
    print(f"  H4      = {result['H4']:.12f}")
    print(f"  H2      = {result['H2']:.12f}")
    print(f"  Delta H = {result['delta_H']:.12f}")

    print("\nPrimzahlvierlinge:")
    print(f"  Anzahl gefundener Vierlinge bis {args.limit}: {result['quadruplet_count']}")
    print(
        f"  Projektion ok: {result['quadruplet_projection_ok_count']} / {result['quadruplet_count']}"
    )

    if args.csv is not None:
        write_csv(args.csv, result["rows"])
        print(f"\nCSV gespeichert: {args.csv}")

    if args.quad_csv is not None:
        write_csv(args.quad_csv, result["quadruplet_rows"])
        print(f"Vierlings-CSV gespeichert: {args.quad_csv}")


if __name__ == "__main__":
    main()