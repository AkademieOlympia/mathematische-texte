from __future__ import annotations

import argparse
import math
from statistics import mean, median
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
        raise ValueError(f"Primzahl {p} hat unpassende Restklasse {r} mod 12")
    return FINE_CLASSES[r]


def coarse_class(p: int) -> str:
    return COARSE_MAP[fine_class(p)]


def quadruplet_pattern(quad: Sequence[int]) -> str:
    return "-".join(coarse_class(q) for q in quad)


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


def twin_primes_after(n: int, sieve: PrimeSieve) -> Tuple[int, int]:
    p = n + 1
    if p % 2 == 0:
        p += 1
    while p + 2 <= sieve.limit:
        if p in sieve and (p + 2) in sieve:
            return (p, p + 2)
        p += 2
    raise RuntimeError("Kein weiterer Primzahlzwilling gefunden")


def B_pi_of_quad(quad: Sequence[int], pivot: int) -> float:
    return sum(math.log(q) for q in quad if q < pivot) - sum(math.log(q) for q in quad if q > pivot)


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


def quantile_groups(rows: Sequence[Dict[str, object]], key: str, q: int) -> List[List[Dict[str, object]]]:
    rows_sorted = sorted(rows, key=lambda r: float(r[key]))
    n = len(rows_sorted)
    groups: List[List[Dict[str, object]]] = []
    for i in range(q):
        a = (i * n) // q
        b = ((i + 1) * n) // q
        groups.append(rows_sorted[a:b])
    return groups


def analyze_quad(quad: Tuple[int, int, int, int], sieve: PrimeSieve, pivot: int) -> Dict[str, object]:
    pattern = quadruplet_pattern(quad)
    nxt = twin_primes_after(quad[-1] + 1, sieve)
    gap = nxt[0] - quad[-1]
    B = B_pi_of_quad(quad, pivot)
    active = int(pattern == "EA-BC-BC-EA")
    T_star = B + 6.0 * active
    return {
        "start": quad[0],
        "pattern": pattern,
        "next_twin": nxt,
        "gap": gap,
        "log_gap": math.log(gap),
        "B_pi": B,
        "T_star": T_star,
        "active_channel": active,
    }


def summarize_group(rows: Sequence[Dict[str, object]], label: str, key: str) -> Dict[str, object]:
    return {
        "label": label,
        "N": len(rows),
        "x_min": min(float(r[key]) for r in rows),
        "x_max": max(float(r[key]) for r in rows),
        "x_mean": mean(float(r[key]) for r in rows),
        "gap_mean": mean(float(r["gap"]) for r in rows),
        "gap_median": median(float(r["gap"]) for r in rows),
        "active_frac": mean(int(r["active_channel"]) for r in rows),
    }


def print_quantiles(title: str, rows: Sequence[Dict[str, object]]) -> None:
    print("\n" + "=" * 112)
    print(title)
    print("=" * 112)
    header = f"{'Q':>4} {'N':>6} {'min':>12} {'max':>12} {'<x>':>12} {'<Gap>':>12} {'MedianGap':>12} {'aktiv':>10}"
    print(header)
    print("-" * len(header))
    for r in rows:
        print(
            f"{r['label']:>4} {r['N']:>6} {r['x_min']:>12.6f} {r['x_max']:>12.6f} "
            f"{r['x_mean']:>12.6f} {r['gap_mean']:>12.4f} {r['gap_median']:>12.4f} "
            f"{r['active_frac']*100:>9.2f}%"
        )


def print_examples(title: str, rows: Sequence[Dict[str, object]], key: str, reverse: bool, limit: int = 15) -> None:
    rows_sorted = sorted(rows, key=lambda r: float(r[key]), reverse=reverse)
    print("\n" + "=" * 120)
    print(title)
    print("=" * 120)
    header = f"{'Start':>8} {'Muster':>20} {'Gap':>8} {'B_pi':>12} {'T*':>12} {'Zwilling':>14}"
    print(header)
    print("-" * len(header))
    for r in rows_sorted[:limit]:
        print(
            f"{int(r['start']):>8} {str(r['pattern']):>20} {int(r['gap']):>8} "
            f"{float(r['B_pi']):>12.6f} {float(r['T_star']):>12.6f} {str(r['next_twin']):>14}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Testet die Distanz zum nächsten Primzahlzwilling gegen B_pi, T* und Vierlingsmuster.")
    parser.add_argument("--limit", type=int, default=2_000_000)
    parser.add_argument("--pivot", type=int, default=11)
    parser.add_argument("--quantiles", type=int, default=5)
    parser.add_argument("--show", type=int, default=15)
    args = parser.parse_args()

    sieve = PrimeSieve(args.limit + 2000)
    quads = dense_quadruplets(args.limit, sieve)
    rows = [analyze_quad(q, sieve, args.pivot) for q in quads]

    gaps = [float(r["gap"]) for r in rows]
    log_gaps = [float(r["log_gap"]) for r in rows]
    B = [float(r["B_pi"]) for r in rows]
    T = [float(r["T_star"]) for r in rows]
    A = [float(r["active_channel"]) for r in rows]

    print("=" * 120)
    print(f"GAP-TEST BIS n = {args.limit}")
    print("=" * 120)
    print(f"Anzahl dichter Vierlinge: {len(rows)}")
    print(f"<Gap>      = {mean(gaps):.6f}")
    print(f"Median Gap = {median(gaps):.6f}")
    print()

    print("Korrelationen mit Gap")
    print(f"corr(B_pi, Gap)      = {pearson_corr(B, gaps):.6f}")
    print(f"corr(T*, Gap)        = {pearson_corr(T, gaps):.6f}")
    print(f"corr(aktiv, Gap)     = {pearson_corr(A, gaps):.6f}")
    print(f"spearman(B_pi, Gap)  = {spearman_corr(B, gaps):.6f}")
    print(f"spearman(T*, Gap)    = {spearman_corr(T, gaps):.6f}")
    print(f"spearman(aktiv, Gap) = {spearman_corr(A, gaps):.6f}")
    print()
    print("Korrelationen mit log(Gap)")
    print(f"corr(B_pi, logGap)      = {pearson_corr(B, log_gaps):.6f}")
    print(f"corr(T*, logGap)        = {pearson_corr(T, log_gaps):.6f}")
    print(f"corr(aktiv, logGap)     = {pearson_corr(A, log_gaps):.6f}")
    print(f"spearman(B_pi, logGap)  = {spearman_corr(B, log_gaps):.6f}")
    print(f"spearman(T*, logGap)    = {spearman_corr(T, log_gaps):.6f}")
    print(f"spearman(aktiv, logGap) = {spearman_corr(A, log_gaps):.6f}")

    by_pattern: Dict[str, List[Dict[str, object]]] = {}
    for p in sorted({str(r["pattern"]) for r in rows}):
        sub = [r for r in rows if str(r["pattern"]) == p]
        by_pattern[p] = sub

    print("\n" + "=" * 80)
    print("NACH VIERLINGSMUSTER")
    print("=" * 80)
    print(f"{'Muster':>20} {'N':>8} {'<Gap>':>12} {'MedianGap':>12} {'<B_pi>':>12} {'<T*>':>12}")
    print("-" * 76)
    for p, sub in by_pattern.items():
        print(
            f"{p:>20} {len(sub):>8} {mean(float(r['gap']) for r in sub):>12.4f} "
            f"{median(float(r['gap']) for r in sub):>12.4f} "
            f"{mean(float(r['B_pi']) for r in sub):>12.6f} {mean(float(r['T_star']) for r in sub):>12.6f}"
        )

    qB = [summarize_group(g, f"Q{i+1}", "B_pi") for i, g in enumerate(quantile_groups(rows, "B_pi", args.quantiles))]
    qT = [summarize_group(g, f"Q{i+1}", "T_star") for i, g in enumerate(quantile_groups(rows, "T_star", args.quantiles))]

    print_quantiles("QUANTILE NACH B_pi", qB)
    print_quantiles("QUANTILE NACH T*", qT)

    print_examples("KLEINSTE GAPS", rows, "gap", reverse=False, limit=args.show)
    print_examples("GRÖSSTE GAPS", rows, "gap", reverse=True, limit=args.show)


if __name__ == "__main__":
    main()