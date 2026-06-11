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


def parse_int_list(s: str) -> List[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def parse_float_list(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


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

    return {
        "start": quad[0],
        "coarse_pattern": coarse_pattern,
        "delta": delta,
        "sigma": sigma,
        "B_log": B_log,
        "active_channel": int(coarse_pattern == "EA-BC-BC-EA"),
        "extendable": int((quad[0] + 12) in sieve),
    }


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


def spearman_corr(x: Sequence[float], y: Sequence[float]) -> float:
    return pearson_corr(rankdata_average(x), rankdata_average(y))


def top_bottom_stats(rows: Sequence[Dict[str, object]], key: str, tail_frac: float) -> Dict[str, float]:
    rows_sorted = sorted(rows, key=lambda r: float(r[key]))
    n = len(rows_sorted)
    k = max(1, int(round(tail_frac * n)))
    bottom = rows_sorted[:k]
    top = rows_sorted[-k:]
    return {
        "bottom_rate": mean(int(r["extendable"]) for r in bottom),
        "top_rate": mean(int(r["extendable"]) for r in top),
        "diff": mean(int(r["extendable"]) for r in top) - mean(int(r["extendable"]) for r in bottom),
        "top_active_frac": mean(int(r["active_channel"]) for r in top),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mehr-Pivot-Test für kalibrierte additive Trigger T_c = B_pi + c*1_{EA-BC-BC-EA}."
    )
    parser.add_argument("--limit", type=int, default=1_000_000, help="obere Grenze")
    parser.add_argument("--pivots", type=str, default="11,13,17,19,29,37", help="kommagetrennte Pivots")
    parser.add_argument("--cs", type=str, default="6,9", help="kommagetrennte c-Kandidaten")
    parser.add_argument("--tail", type=float, default=0.10, help="Top/Bottom-Anteil")
    args = parser.parse_args()

    pivots = parse_int_list(args.pivots)
    c_values = parse_float_list(args.cs)

    sieve = PrimeSieve(args.limit + 12)
    quads = dense_quadruplets(args.limit, sieve)

    print("=" * 144)
    print("MEHR-PIVOT-TEST FÜR KALIBRIERTE EINZELTRIGGER")
    print("=" * 144)
    header = (
        f"{'π':>5} {'c':>6} {'Pearson':>12} {'Spearman':>12} "
        f"{'Grundrate':>10} {'Top10%':>10} {'Bottom10%':>12} {'Diff.':>10} {'Top aktiv':>11}"
    )
    print(header)
    print("-" * len(header))

    rows_out: List[Dict[str, object]] = []

    for pivot in pivots:
        base_rows = [quadruplet_features(q, pivot, sieve) for q in quads]
        y_vals = [float(r["extendable"]) for r in base_rows]
        base_rate = mean(y_vals)

        for c in c_values:
            augmented = []
            for r in base_rows:
                rr = dict(r)
                rr["T_c"] = float(rr["B_log"]) + c * int(rr["active_channel"])
                augmented.append(rr)

            x_vals = [float(r["T_c"]) for r in augmented]
            pear = pearson_corr(x_vals, y_vals)
            spear = spearman_corr(x_vals, y_vals)
            tb = top_bottom_stats(augmented, "T_c", args.tail)

            row = {
                "pivot": pivot,
                "c": c,
                "pearson": pear,
                "spearman": spear,
                "base_rate": base_rate,
                "top_rate": tb["top_rate"],
                "bottom_rate": tb["bottom_rate"],
                "diff": tb["diff"],
                "top_active_frac": tb["top_active_frac"],
            }
            rows_out.append(row)

            print(
                f"{pivot:>5} {c:>6.2f} {pear:>12.6f} {spear:>12.6f} "
                f"{base_rate*100:>9.2f}% {tb['top_rate']*100:>9.2f}% {tb['bottom_rate']*100:>11.2f}% "
                f"{tb['diff']*100:>9.2f}% {tb['top_active_frac']*100:>10.2f}%"
            )

    print("\n" + "=" * 80)
    print("KOMPAKT: MITTELWERTE ÜBER ALLE PIVOTS")
    print("=" * 80)
    for c in c_values:
        sub = [r for r in rows_out if abs(r["c"] - c) < 1e-12]
        print(
            f"c = {c:.2f} | "
            f"Pearson = {mean(r['pearson'] for r in sub):.6f}, "
            f"Spearman = {mean(r['spearman'] for r in sub):.6f}, "
            f"Top10% = {mean(r['top_rate'] for r in sub)*100:.2f}%, "
            f"Bottom10% = {mean(r['bottom_rate'] for r in sub)*100:.2f}%, "
            f"Diff. = {mean(r['diff'] for r in sub)*100:.2f}%"
        )


if __name__ == "__main__":
    main()