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

    active_channel = (coarse_pattern == "EA-BC-BC-EA")
    trigger = B_log if active_channel else 0.0

    return {
        "start": quad[0],
        "quadruplet": quad,
        "coarse_pattern": coarse_pattern,
        "delta": delta,
        "sigma": sigma,
        "B_log": B_log,
        "beta5": beta5,
        "active_channel": int(active_channel),
        "trigger_T": trigger,
        "extendable": int((quad[0] + 12) in sieve),
    }


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
    if not rows:
        return {
            "label": label,
            "N": 0,
            "x_min": None,
            "x_max": None,
            "x_mean": None,
            "rate": None,
            "delta_mean": None,
            "sigma_mean": None,
            "active_frac": None,
        }

    return {
        "label": label,
        "N": len(rows),
        "x_min": min(float(r[key]) for r in rows),
        "x_max": max(float(r[key]) for r in rows),
        "x_mean": mean(float(r[key]) for r in rows),
        "rate": mean(int(r["extendable"]) for r in rows),
        "delta_mean": mean(float(r["delta"]) for r in rows),
        "sigma_mean": mean(float(r["sigma"]) for r in rows),
        "active_frac": mean(int(r["active_channel"]) for r in rows),
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


def print_quantiles(rows: Sequence[Dict[str, object]], title: str) -> None:
    print(title)
    header = f"{'Gruppe':>8} {'N':>6} {'min':>12} {'max':>12} {'<x>':>12} {'Rate':>10} {'aktiver Kanal':>14} {'<Δ>':>8} {'<Σ>':>8}"
    print(header)
    print("-" * len(header))
    for row in rows:
        print(
            f"{row['label']:>8} {int(row['N']):>6} "
            f"{float(row['x_min']):>12.6f} {float(row['x_max']):>12.6f} {float(row['x_mean']):>12.6f} "
            f"{float(row['rate'])*100:>9.2f}% {float(row['active_frac'])*100:>13.2f}% "
            f"{float(row['delta_mean']):>8.3f} {float(row['sigma_mean']):>8.3f}"
        )


def print_examples(rows: Sequence[Dict[str, object]], title: str, key: str, limit: int = 15, reverse: bool = True) -> None:
    rows_sorted = sorted(rows, key=lambda r: float(r[key]), reverse=reverse)
    print("\n" + "=" * 118)
    print(title)
    print("=" * 118)
    header = f"{'Start':>8} {'Typ':>20} {'B':>12} {'T':>12} {'Δ':>4} {'Σ':>4} {'ergänzbar':>10}"
    print(header)
    print("-" * len(header))
    for r in rows_sorted[:limit]:
        print(
            f"{int(r['start']):>8} {str(r['coarse_pattern']):>20} "
            f"{float(r['B_log']):>12.6f} {float(r['trigger_T']):>12.6f} "
            f"{int(r['delta']):>4} {int(r['sigma']):>4} {int(r['extendable']):>10}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Testet den einzelnen Trigger T(Q4)=1_{EA-BC-BC-EA} * B_pi(Q4)."
    )
    parser.add_argument("--limit", type=int, default=1_000_000, help="obere Grenze")
    parser.add_argument("--pivot", type=int, default=11, help="Pivot-Primzahl π")
    parser.add_argument("--quantiles", type=int, default=5, help="Anzahl Quantile")
    parser.add_argument("--show", type=int, default=15, help="Anzahl Beispiele")
    args = parser.parse_args()

    sieve = PrimeSieve(args.limit + 12)
    quads = dense_quadruplets(args.limit, sieve)
    rows = [quadruplet_features(q, args.pivot, sieve) for q in quads]

    overall_rate = mean(int(r["extendable"]) for r in rows)

    B_vals = [float(r["B_log"]) for r in rows]
    T_vals = [float(r["trigger_T"]) for r in rows]
    Y_vals = [float(r["extendable"]) for r in rows]

    corr_B = pearson_corr(B_vals, Y_vals)
    corr_T = pearson_corr(T_vals, Y_vals)
    spear_B = spearman_corr(B_vals, Y_vals)
    spear_T = spearman_corr(T_vals, Y_vals)

    qB = quantile_groups(rows, args.quantiles, "B_log")
    qT = quantile_groups(rows, args.quantiles, "trigger_T")

    qB_sum = [summarize_group(g, f"Q{i+1}", "B_log") for i, g in enumerate(qB)]
    qT_sum = [summarize_group(g, f"Q{i+1}", "trigger_T") for i, g in enumerate(qT)]

    print("=" * 118)
    print("TEST DES EINZELTRIGGERS  T(Q4)=1_{EA-BC-BC-EA} * B_pi(Q4)")
    print("=" * 118)
    print(f"Pivot π               : {args.pivot}")
    print(f"Anzahl Vierlinge      : {len(rows)}")
    print(f"Gesamt-Ergänzungsrate : {overall_rate*100:.2f}%")
    print()
    print("Korrelationen")
    print(f"  Pearson  corr(B, Ergänzbarkeit) : {corr_B:.6f}")
    print(f"  Pearson  corr(T, Ergänzbarkeit) : {corr_T:.6f}")
    print(f"  Spearman corr(B, Ergänzbarkeit) : {spear_B:.6f}")
    print(f"  Spearman corr(T, Ergänzbarkeit) : {spear_T:.6f}")
    print()

    print_quantiles(qB_sum, "Quantile nach B_pi(Q4)")
    print()
    print_quantiles(qT_sum, "Quantile nach Trigger T(Q4)")

    print_examples(rows, "GRÖSSTE T-WERTE", "trigger_T", limit=args.show, reverse=True)
    print_examples(rows, "KLEINSTE T-WERTE", "trigger_T", limit=args.show, reverse=False)


if __name__ == "__main__":
    main()