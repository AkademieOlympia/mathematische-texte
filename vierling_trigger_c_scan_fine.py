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


def frange(start: float, stop: float, step: float) -> List[float]:
    vals: List[float] = []
    x = start
    while x <= stop + 1e-12:
        vals.append(round(x, 10))
        x += step
    return vals


def print_table(rows: Sequence[Dict[str, object]]) -> None:
    print("=" * 128)
    print("FEINSCAN DES KANALBONUS c")
    print("=" * 128)
    header = (
        f"{'c':>8} {'Pearson':>12} {'Spearman':>12} "
        f"{'Top10%':>10} {'Bottom10%':>12} {'Diff.':>10} {'Top aktiv':>11}"
    )
    print(header)
    print("-" * len(header))
    for r in rows:
        print(
            f"{r['c']:>8.2f} {r['pearson']:>12.6f} {r['spearman']:>12.6f} "
            f"{r['top_rate']*100:>9.2f}% {r['bottom_rate']*100:>11.2f}% "
            f"{r['diff']*100:>9.2f}% {r['top_active_frac']*100:>10.2f}%"
        )


def print_best(rows: Sequence[Dict[str, object]]) -> None:
    best_p = max(rows, key=lambda r: r["pearson"])
    best_s = max(rows, key=lambda r: r["spearman"])
    best_d = max(rows, key=lambda r: r["diff"])

    print("\n" + "=" * 80)
    print("BESTE WERTE IM FEINSCAN")
    print("=" * 80)
    print(f"Bestes Pearson        : c = {best_p['c']:.2f}, Wert = {best_p['pearson']:.6f}")
    print(f"Bestes Spearman       : c = {best_s['c']:.2f}, Wert = {best_s['spearman']:.6f}")
    print(f"Beste Top-Bottom-Diff.: c = {best_d['c']:.2f}, Wert = {best_d['diff']*100:.2f}%")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Feinscan des kanal-korrigierten Triggers T_c = B_pi + c*1_{EA-BC-BC-EA}."
    )
    parser.add_argument("--limit", type=int, default=1_000_000, help="obere Grenze")
    parser.add_argument("--pivot", type=int, default=11, help="Pivot-Primzahl π")
    parser.add_argument("--c-min", type=float, default=5.0, help="minimales c")
    parser.add_argument("--c-max", type=float, default=9.0, help="maximales c")
    parser.add_argument("--c-step", type=float, default=0.5, help="Schrittweite")
    parser.add_argument("--tail", type=float, default=0.10, help="Top/Bottom-Anteil")
    args = parser.parse_args()

    sieve = PrimeSieve(args.limit + 12)
    quads = dense_quadruplets(args.limit, sieve)
    rows = [quadruplet_features(q, args.pivot, sieve) for q in quads]
    y_vals = [float(r["extendable"]) for r in rows]

    c_values = frange(args.c_min, args.c_max, args.c_step)
    out_rows: List[Dict[str, object]] = []

    for c in c_values:
        augmented = []
        for r in rows:
            rr = dict(r)
            rr["T_c"] = float(rr["B_log"]) + c * int(rr["active_channel"])
            augmented.append(rr)

        x_vals = [float(r["T_c"]) for r in augmented]
        pear = pearson_corr(x_vals, y_vals)
        spear = spearman_corr(x_vals, y_vals)
        tb = top_bottom_stats(augmented, "T_c", args.tail)

        out_rows.append(
            {
                "c": c,
                "pearson": pear,
                "spearman": spear,
                "top_rate": tb["top_rate"],
                "bottom_rate": tb["bottom_rate"],
                "diff": tb["diff"],
                "top_active_frac": tb["top_active_frac"],
            }
        )

    print(f"Pivot π = {args.pivot}")
    print(f"Anzahl Vierlinge = {len(rows)}")
    print(f"Gesamt-Ergänzungsrate = {mean(y_vals)*100:.2f}%\n")

    print_table(out_rows)
    print_best(out_rows)


if __name__ == "__main__":
    main()