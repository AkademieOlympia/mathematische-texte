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
        "delta": delta,
        "sigma": sigma,
        "B_log": B_log,
        "beta5": beta5,
        "extendable": int((quad[0] + 12) in sieve),
    }


def rankdata_average(values: Sequence[float]) -> List[float]:
    indexed = sorted(enumerate(values), key=lambda t: t[1])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i
        while j + 1 < len(indexed) and indexed[j + 1][1] == indexed[i][1]:
            j += 1
        avg_rank = (i + j + 2) / 2.0
        for k in range(i, j + 1):
            ranks[indexed[k][0]] = avg_rank
        i = j + 1
    return ranks


def pearson_corr(x: Sequence[float], y: Sequence[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        return float("nan")
    mx = mean(x)
    my = mean(y)
    sx = sum((xi - mx) ** 2 for xi in x)
    sy = sum((yi - my) ** 2 for yi in y)
    if sx == 0 or sy == 0:
        return float("nan")
    sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    return sxy / math.sqrt(sx * sy)


def spearman_corr(x: Sequence[float], y: Sequence[float]) -> float:
    rx = rankdata_average(x)
    ry = rankdata_average(y)
    return pearson_corr(rx, ry)


def point_biserial_corr(x: Sequence[float], y01: Sequence[int]) -> float:
    return pearson_corr(x, [float(v) for v in y01])


def summarize_bucket(rows: Sequence[Dict[str, object]], label: str) -> Dict[str, object]:
    if not rows:
        return {
            "label": label,
            "N": 0,
            "extendable_rate": None,
            "B_mean": None,
            "B_min": None,
            "B_max": None,
            "delta_mean": None,
            "sigma_mean": None,
        }
    return {
        "label": label,
        "N": len(rows),
        "extendable_rate": mean(int(r["extendable"]) for r in rows),
        "B_mean": mean(float(r["B_log"]) for r in rows),
        "B_min": min(float(r["B_log"]) for r in rows),
        "B_max": max(float(r["B_log"]) for r in rows),
        "delta_mean": mean(float(r["delta"]) for r in rows),
        "sigma_mean": mean(float(r["sigma"]) for r in rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Top-/Bottom-Test und Rangkorrelation für B_pi(Q4) gegen Fünflingsergänzbarkeit."
    )
    parser.add_argument("--limit", type=int, default=1_000_000, help="obere Grenze")
    parser.add_argument("--pivot", type=int, default=11, help="Pivot-Primzahl π")
    parser.add_argument("--tail", type=float, default=0.10, help="Anteil für Top/Bottom-Vergleich, z.B. 0.10")
    parser.add_argument("--show", type=int, default=15, help="wie viele Beispiele aus Top/Bottom gezeigt werden")
    args = parser.parse_args()

    if not (0 < args.tail < 0.5):
        raise ValueError("--tail muss zwischen 0 und 0.5 liegen")

    sieve = PrimeSieve(args.limit + 12)
    quads = dense_quadruplets(args.limit, sieve)
    rows = [quadruplet_features(q, args.pivot, sieve) for q in quads]

    rows_sorted = sorted(rows, key=lambda r: float(r["B_log"]))
    n = len(rows_sorted)
    k = max(1, int(round(args.tail * n)))

    bottom = rows_sorted[:k]
    top = rows_sorted[-k:]

    all_B = [float(r["B_log"]) for r in rows]
    all_Y = [int(r["extendable"]) for r in rows]

    pb = point_biserial_corr(all_B, all_Y)
    sp = spearman_corr(all_B, all_Y)

    overall_rate = mean(all_Y)

    bottom_sum = summarize_bucket(bottom, f"Bottom {args.tail:.0%}")
    top_sum = summarize_bucket(top, f"Top {args.tail:.0%}")

    print("=" * 110)
    print("TOP-/BOTTOM-TEST UND RANGKORRELATION FÜR B_pi(Q4)")
    print("=" * 110)
    print(f"Pivot π                 : {args.pivot}")
    print(f"Anzahl Vierlinge        : {n}")
    print(f"Gesamt-Ergänzungsrate   : {overall_rate*100:.2f}%")
    print(f"Tail-Anteil             : {args.tail:.0%}")
    print()

    print("Korrelationen")
    print(f"  Point-biserial corr(B, Ergänzbarkeit) : {pb:.6f}")
    print(f"  Spearman corr(B, Ergänzbarkeit)       : {sp:.6f}")
    print()

    print("Top-/Bottom-Vergleich")
    header = f"{'Gruppe':>12} {'N':>6} {'Rate':>10} {'<B>':>12} {'B_min':>12} {'B_max':>12} {'<Δ>':>8} {'<Σ>':>8}"
    print(header)
    print("-" * len(header))
    for row in [bottom_sum, top_sum]:
        print(
            f"{row['label']:>12} {int(row['N']):>6} {float(row['extendable_rate'])*100:>9.2f}% "
            f"{float(row['B_mean']):>12.6f} {float(row['B_min']):>12.6f} {float(row['B_max']):>12.6f} "
            f"{float(row['delta_mean']):>8.3f} {float(row['sigma_mean']):>8.3f}"
        )

    print("\n" + "=" * 110)
    print("BEISPIELE AUS DEM BOTTOM-BEREICH")
    print("=" * 110)
    header2 = f"{'Start':>8} {'B':>12} {'Δ':>4} {'Σ':>4} {'ergänzbar':>10}"
    print(header2)
    print("-" * len(header2))
    for r in bottom[: args.show]:
        print(
            f"{int(r['start']):>8} {float(r['B_log']):>12.6f} "
            f"{int(r['delta']):>4} {int(r['sigma']):>4} {bool(r['extendable']):>10}"
        )

    print("\n" + "=" * 110)
    print("BEISPIELE AUS DEM TOP-BEREICH")
    print("=" * 110)
    print(header2)
    print("-" * len(header2))
    for r in reversed(top[-args.show:]):
        print(
            f"{int(r['start']):>8} {float(r['B_log']):>12.6f} "
            f"{int(r['delta']):>4} {int(r['sigma']):>4} {bool(r['extendable']):>10}"
        )


if __name__ == "__main__":
    main()