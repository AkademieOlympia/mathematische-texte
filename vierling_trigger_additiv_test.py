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
        "quadruplet": quad,
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


def quantile_groups(rows: Sequence[Dict[str, object]], num_groups: int, key: str) -> List[List[Dict[str, object]]]:
    rows_sorted = sorted(rows, key=lambda r: float(r[key]))
    n = len(rows_sorted)
    groups: List[List[Dict[str, object]]] = []
    for i in range(num_groups):
        a = (i * n) // num_groups
        b = ((i + 1) * n) // num_groups
        groups.append(rows_sorted[a:b])
    return groups


def summarize_group(rows: Sequence[Dict[str, object]], label: str, key: str) -> Dict[str, object]:
    return {
        "label": label,
        "N": len(rows),
        "x_min": min(float(r[key]) for r in rows),
        "x_max": max(float(r[key]) for r in rows),
        "x_mean": mean(float(r[key]) for r in rows),
        "rate": mean(int(r["extendable"]) for r in rows),
        "active_frac": mean(int(r["active_channel"]) for r in rows),
    }


def parse_cs(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def print_summary(rows: Sequence[Dict[str, object]]) -> None:
    print("=" * 116)
    print("ADDITIVER EINZELTRIGGERTEST  T_c(Q4)=B_pi(Q4)+c*1_{EA-BC-BC-EA}")
    print("=" * 116)
    header = f"{'c':>8} {'Pearson':>12} {'Spearman':>12} {'Top10%':>10} {'Bottom10%':>12} {'Diff.':>10}"
    print(header)
    print("-" * len(header))
    for r in rows:
        print(
            f"{r['c']:>8.2f} {r['pearson']:>12.6f} {r['spearman']:>12.6f} "
            f"{r['top_rate']*100:>9.2f}% {r['bottom_rate']*100:>11.2f}% {r['diff']*100:>9.2f}%"
        )


def print_quantiles(rows: Sequence[Dict[str, object]], c_value: float) -> None:
    print("\n" + "=" * 110)
    print(f"QUANTILE FÜR c = {c_value}")
    print("=" * 110)
    header = f"{'Q':>4} {'N':>6} {'min':>12} {'max':>12} {'<T_c>':>12} {'Rate':>10} {'aktiver Kanal':>14}"
    print(header)
    print("-" * len(header))
    for r in rows:
        print(
            f"{r['label']:>4} {r['N']:>6} {r['x_min']:>12.6f} {r['x_max']:>12.6f} "
            f"{r['x_mean']:>12.6f} {r['rate']*100:>9.2f}% {r['active_frac']*100:>13.2f}%"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Testet additive Einzeltrigger T_c = B_pi + c*1_{EA-BC-BC-EA}."
    )
    parser.add_argument("--limit", type=int, default=1_000_000, help="obere Grenze")
    parser.add_argument("--pivot", type=int, default=11, help="Pivot-Primzahl π")
    parser.add_argument("--cs", type=str, default="0,2,5,10,15", help="kommagetrennte c-Werte")
    parser.add_argument("--quantiles", type=int, default=5, help="Anzahl Quantile")
    parser.add_argument("--show-best", type=float, default=0.0, help="für diesen c-Wert zusätzlich Quantile ausgeben")
    args = parser.parse_args()

    c_values = parse_cs(args.cs)

    sieve = PrimeSieve(args.limit + 12)
    quads = dense_quadruplets(args.limit, sieve)
    base_rows = [quadruplet_features(q, args.pivot, sieve) for q in quads]

    y_vals = [float(r["extendable"]) for r in base_rows]
    n = len(base_rows)
    k = max(1, int(round(0.10 * n)))

    summary_rows: List[Dict[str, object]] = []
    quantile_dump: List[Dict[str, object]] | None = None

    for c in c_values:
        rows = []
        for r in base_rows:
            rr = dict(r)
            rr["T_c"] = float(rr["B_log"]) + c * int(rr["active_channel"])
            rows.append(rr)

        x_vals = [float(r["T_c"]) for r in rows]
        pear = pearson_corr(x_vals, y_vals)
        spear = spearman_corr(x_vals, y_vals)

        rows_sorted = sorted(rows, key=lambda r: float(r["T_c"]))
        bottom = rows_sorted[:k]
        top = rows_sorted[-k:]

        bottom_rate = mean(int(r["extendable"]) for r in bottom)
        top_rate = mean(int(r["extendable"]) for r in top)

        summary_rows.append(
            {
                "c": c,
                "pearson": pear,
                "spearman": spear,
                "top_rate": top_rate,
                "bottom_rate": bottom_rate,
                "diff": top_rate - bottom_rate,
            }
        )

        if abs(c - args.show_best) < 1e-12:
            groups = quantile_groups(rows, args.quantiles, "T_c")
            quantile_dump = [summarize_group(g, f"Q{i+1}", "T_c") for i, g in enumerate(groups)]

    print(f"Pivot π = {args.pivot}")
    print(f"Anzahl Vierlinge = {n}")
    print(f"Gesamt-Ergänzungsrate = {mean(y_vals)*100:.2f}%\n")
    print_summary(summary_rows)

    if quantile_dump is not None:
        print_quantiles(quantile_dump, args.show_best)


if __name__ == "__main__":
    main()