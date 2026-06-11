from __future__ import annotations

import argparse
import math
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


def parse_pivots(s: str) -> List[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


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
    coarse_pattern = "-".join(coarse_class(q) for q in quad)

    n_minus_ea = n_minus_bc = n_plus_ea = n_plus_bc = 0
    for q in quad:
        cc = coarse_class(q)
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
        "start": quad[0],
        "coarse_pattern": coarse_pattern,
        "delta": delta,
        "sigma": sigma,
        "B_log": B_log,
        "beta5": beta5,
        "extendable": int((quad[0] + 12) in sieve),
    }


def summarize(rows: Sequence[Dict[str, object]]) -> Dict[str, object]:
    if not rows:
        return {
            "N": 0,
            "rate": None,
            "B_mean": None,
            "B_min": None,
            "B_max": None,
            "delta_mean": None,
            "sigma_mean": None,
        }
    return {
        "N": len(rows),
        "rate": mean(int(r["extendable"]) for r in rows),
        "B_mean": mean(float(r["B_log"]) for r in rows),
        "B_min": min(float(r["B_log"]) for r in rows),
        "B_max": max(float(r["B_log"]) for r in rows),
        "delta_mean": mean(float(r["delta"]) for r in rows),
        "sigma_mean": mean(float(r["sigma"]) for r in rows),
    }


def bucket_by_threshold(rows: Sequence[Dict[str, object]], threshold: float) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    low = [r for r in rows if float(r["B_log"]) < threshold]
    high = [r for r in rows if float(r["B_log"]) >= threshold]
    return low, high


def choose_threshold(rows: Sequence[Dict[str, object]], mode: str) -> float:
    Bs = sorted(float(r["B_log"]) for r in rows)
    if mode == "median":
        return Bs[len(Bs) // 2]
    if mode == "mean":
        return mean(Bs)
    if mode == "q80":
        return Bs[int(0.80 * len(Bs))]
    if mode == "q90":
        return Bs[int(0.90 * len(Bs))]
    raise ValueError("unbekannter threshold-mode")


def print_table(results: Sequence[Dict[str, object]]) -> None:
    print("=" * 158)
    print("MEHRPIVOT-2D-TEST: MUSTERTYP × B-HOCH/NIEDRIG")
    print("=" * 158)
    header = (
        f"{'π':>5} {'Schwelle':>12} {'Typ':>20} {'B-Bereich':>12} {'N':>6} "
        f"{'Rate':>10} {'<B>':>12} {'B_min':>12} {'B_max':>12} {'<Δ>':>8} {'<Σ>':>8}"
    )
    print(header)
    print("-" * len(header))
    for row in results:
        if row["N"] == 0:
            print(
                f"{int(row['pivot']):>5} {float(row['threshold']):>12.6f} {row['pattern']:>20} "
                f"{row['bucket']:>12} {0:>6}"
            )
        else:
            print(
                f"{int(row['pivot']):>5} {float(row['threshold']):>12.6f} {row['pattern']:>20} "
                f"{row['bucket']:>12} {int(row['N']):>6} {float(row['rate'])*100:>9.2f}% "
                f"{float(row['B_mean']):>12.6f} {float(row['B_min']):>12.6f} {float(row['B_max']):>12.6f} "
                f"{float(row['delta_mean']):>8.3f} {float(row['sigma_mean']):>8.3f}"
            )


def print_compact(results: Sequence[Dict[str, object]]) -> None:
    print("\n" + "=" * 108)
    print("KOMPAKT: HOCH minus NIEDRIG je Typ und Pivot")
    print("=" * 108)
    header = f"{'π':>5} {'Typ':>20} {'Rate hoch':>12} {'Rate niedrig':>14} {'Differenz':>12}"
    print(header)
    print("-" * len(header))

    pivots = sorted(set(int(r["pivot"]) for r in results))
    patterns = sorted(set(str(r["pattern"]) for r in results))

    for pivot in pivots:
        for pattern in patterns:
            low = [r for r in results if int(r["pivot"]) == pivot and str(r["pattern"]) == pattern and str(r["bucket"]) == "niedrig"][0]
            high = [r for r in results if int(r["pivot"]) == pivot and str(r["pattern"]) == pattern and str(r["bucket"]) == "hoch"][0]
            rl = float(low["rate"]) if low["rate"] is not None else float("nan")
            rh = float(high["rate"]) if high["rate"] is not None else float("nan")
            diff = rh - rl
            print(f"{pivot:>5} {pattern:>20} {rh*100:>11.2f}% {rl*100:>13.2f}% {diff*100:>11.2f}%")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mehrpivot-Variante des 2D-Tests: Mustertyp × hoher/niedriger B_pi-Wert."
    )
    parser.add_argument("--limit", type=int, default=1_000_000, help="obere Grenze")
    parser.add_argument("--pivots", type=str, default="11,13,17,19,29,37", help="kommagetrennte Pivot-Primzahlen")
    parser.add_argument(
        "--threshold-mode",
        type=str,
        default="median",
        choices=["median", "mean", "q80", "q90"],
        help="wie die B-Schwelle je Pivot gewählt wird",
    )
    args = parser.parse_args()

    pivots = parse_pivots(args.pivots)
    sieve = PrimeSieve(args.limit + 12)
    quads = dense_quadruplets(args.limit, sieve)

    results: List[Dict[str, object]] = []

    for pivot in pivots:
        rows = [quadruplet_features(q, pivot, sieve) for q in quads]
        threshold = choose_threshold(rows, args.threshold_mode)
        patterns = sorted(set(str(r["coarse_pattern"]) for r in rows))

        for pattern in patterns:
            sub = [r for r in rows if str(r["coarse_pattern"]) == pattern]
            low, high = bucket_by_threshold(sub, threshold)

            low_sum = summarize(low)
            low_sum["pivot"] = pivot
            low_sum["threshold"] = threshold
            low_sum["pattern"] = pattern
            low_sum["bucket"] = "niedrig"
            results.append(low_sum)

            high_sum = summarize(high)
            high_sum["pivot"] = pivot
            high_sum["threshold"] = threshold
            high_sum["pattern"] = pattern
            high_sum["bucket"] = "hoch"
            results.append(high_sum)

    print_table(results)
    print_compact(results)


if __name__ == "__main__":
    main()