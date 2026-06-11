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
        "quadruplet": quad,
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


def print_table(results: List[Dict[str, object]], threshold: float, threshold_label: str) -> None:
    print("=" * 132)
    print("2D-TEST: B_pi-HOCH/NIEDRIG × VIERLINGSMUSTER-TYP")
    print("=" * 132)
    print(f"Schwelle ({threshold_label}) = {threshold:.6f}\n")
    header = f"{'Typ':>20} {'B-Bereich':>12} {'N':>6} {'Rate':>10} {'<B>':>12} {'B_min':>12} {'B_max':>12} {'<Δ>':>8} {'<Σ>':>8}"
    print(header)
    print("-" * len(header))
    for row in results:
        if row["N"] == 0:
            print(f"{row['pattern']:>20} {row['bucket']:>12} {0:>6}")
        else:
            print(
                f"{row['pattern']:>20} {row['bucket']:>12} {int(row['N']):>6} "
                f"{float(row['rate'])*100:>9.2f}% {float(row['B_mean']):>12.6f} "
                f"{float(row['B_min']):>12.6f} {float(row['B_max']):>12.6f} "
                f"{float(row['delta_mean']):>8.3f} {float(row['sigma_mean']):>8.3f}"
            )


def print_examples(rows: Sequence[Dict[str, object]], title: str, limit: int = 12) -> None:
    print("\n" + "=" * 110)
    print(title)
    print("=" * 110)
    if not rows:
        print("Keine Daten.")
        return
    header = f"{'Start':>8} {'Muster':>18} {'B':>12} {'Δ':>4} {'Σ':>4} {'ergänzbar':>10}"
    print(header)
    print("-" * len(header))
    for r in rows[:limit]:
        print(
            f"{int(r['start']):>8} {str(r['coarse_pattern']):>18} "
            f"{float(r['B_log']):>12.6f} {int(r['delta']):>4} {int(r['sigma']):>4} {int(r['extendable']):>10}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Vergleicht Ergänzungsraten in einem 2D-Test: hoher/niedriger B_pi-Wert × Vierlingsmuster-Typ."
    )
    parser.add_argument("--limit", type=int, default=1_000_000, help="obere Grenze")
    parser.add_argument("--pivot", type=int, default=11, help="Pivot-Primzahl π")
    parser.add_argument(
        "--threshold-mode",
        type=str,
        default="median",
        choices=["median", "mean", "q80", "q90"],
        help="wie die B-Schwelle gewählt wird",
    )
    parser.add_argument("--show", type=int, default=12, help="wie viele Beispiele pro Feld gezeigt werden")
    args = parser.parse_args()

    sieve = PrimeSieve(args.limit + 12)
    quads = dense_quadruplets(args.limit, sieve)
    rows = [quadruplet_features(q, args.pivot, sieve) for q in quads]

    Bs = sorted(float(r["B_log"]) for r in rows)
    if args.threshold_mode == "median":
        threshold = Bs[len(Bs) // 2]
    elif args.threshold_mode == "mean":
        threshold = mean(Bs)
    elif args.threshold_mode == "q80":
        threshold = Bs[int(0.80 * len(Bs))]
    elif args.threshold_mode == "q90":
        threshold = Bs[int(0.90 * len(Bs))]
    else:
        raise ValueError("unbekannter threshold-mode")

    patterns = sorted(set(str(r["coarse_pattern"]) for r in rows))
    results: List[Dict[str, object]] = []

    for pattern in patterns:
        sub = [r for r in rows if str(r["coarse_pattern"]) == pattern]
        low, high = bucket_by_threshold(sub, threshold)

        low_sum = summarize(low)
        low_sum["pattern"] = pattern
        low_sum["bucket"] = "niedrig"
        results.append(low_sum)

        high_sum = summarize(high)
        high_sum["pattern"] = pattern
        high_sum["bucket"] = "hoch"
        results.append(high_sum)

    print(f"Pivot π = {args.pivot}")
    print(f"Anzahl dichter Vierlinge: {len(rows)}")
    print(f"Gesamt-Ergänzungsrate: {100 * mean(int(r['extendable']) for r in rows):.2f}%")
    print_table(results, threshold, args.threshold_mode)

    for pattern in patterns:
        sub = [r for r in rows if str(r["coarse_pattern"]) == pattern]
        low, high = bucket_by_threshold(sub, threshold)
        low_sorted = sorted(low, key=lambda r: float(r["B_log"]), reverse=True)
        high_sorted = sorted(high, key=lambda r: float(r["B_log"]), reverse=True)

        print_examples(high_sorted, f"BEISPIELE: {pattern} mit hohem B", limit=args.show)
        print_examples(low_sorted, f"BEISPIELE: {pattern} mit niedrigem B", limit=args.show)


if __name__ == "__main__":
    main()