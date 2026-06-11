from __future__ import annotations

import argparse
import math
from collections import Counter, defaultdict
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


def fine_class(p: int) -> str:
    r = p % 12
    if r not in FINE_CLASSES:
        raise ValueError(f"Primzahl {p} hat unpassende Restklasse {r} mod 12")
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


def twin_primes_after(n: int, sieve: PrimeSieve) -> Tuple[int, int]:
    p = n + 1
    if p % 2 == 0:
        p += 1
    while p + 2 <= sieve.limit:
        if p in sieve and (p + 2) in sieve:
            return (p, p + 2)
        p += 2
    raise RuntimeError("Kein weiterer Primzahlzwilling im Siebbereich gefunden")


def twin_type(twin: Tuple[int, int]) -> str:
    a, b = fine_class(twin[0]), fine_class(twin[1])
    if (a, b) == ("A", "B"):
        return "AB"
    if (a, b) == ("C", "E"):
        return "CE"
    raise ValueError(f"Unerwarteter Zwillingstyp: {twin} -> {(a, b)}")


def quadruplet_pattern(quad: Sequence[int]) -> str:
    return "-".join(coarse_class(q) for q in quad)


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


def analyze_quad_to_next_twin(
    quad: Tuple[int, int, int, int],
    sieve: PrimeSieve,
    pivot: int,
) -> Dict[str, object]:
    q_end = quad[-1]
    next_twin = twin_primes_after(q_end + 1, sieve)
    ttype = twin_type(next_twin)
    pattern = quadruplet_pattern(quad)

    if ttype == "AB":
        target_families = {"A", "B"}
        foreign_families = {"C", "E"}
    elif ttype == "CE":
        target_families = {"C", "E"}
        foreign_families = {"A", "B"}
    else:
        raise ValueError("unbekannter Zwillingstyp")

    between_primes = [p for p in sieve.primes if q_end < p < next_twin[0]]
    between_classes = [fine_class(p) for p in between_primes]

    family_counts = Counter(between_classes)
    foreign_count = sum(family_counts[f] for f in foreign_families)
    target_count = sum(family_counts[f] for f in target_families)

    B = B_pi_of_quad(quad, pivot)
    active = int(pattern == "EA-BC-BC-EA")
    T_star = B + 6.0 * active

    return {
        "start": quad[0],
        "quad": quad,
        "pattern": pattern,
        "next_twin_start": next_twin[0],
        "next_twin": next_twin,
        "next_twin_type": ttype,
        "gap_to_twin": next_twin[0] - q_end,
        "between_count": len(between_primes),
        "family_counts": dict(family_counts),
        "foreign_count": foreign_count,
        "target_count": target_count,
        "foreign_ratio": (foreign_count / len(between_primes)) if between_primes else 0.0,
        "B_pi": B,
        "T_star": T_star,
        "active_channel": active,
    }


def fmt_counts(d: Dict[str, int]) -> str:
    return f"E:{d.get('E',0)} A:{d.get('A',0)} B:{d.get('B',0)} C:{d.get('C',0)}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Testet Zwischenfamilien zwischen Primzahlvierling und nächstem Primzahlzwilling."
    )
    parser.add_argument("--limit", type=int, default=2_000_000, help="obere Grenze")
    parser.add_argument("--pivot", type=int, default=11, help="Pivot π für B_pi und T*")
    parser.add_argument("--show", type=int, default=20, help="Anzahl Beispielzeilen")
    args = parser.parse_args()

    sieve = PrimeSieve(args.limit + 1000)
    quads = dense_quadruplets(args.limit, sieve)
    rows = [analyze_quad_to_next_twin(q, sieve, args.pivot) for q in quads]

    X_vals = [float(r["foreign_count"]) for r in rows]
    B_vals = [float(r["B_pi"]) for r in rows]
    T_vals = [float(r["T_star"]) for r in rows]
    gap_vals = [float(r["gap_to_twin"]) for r in rows]

    by_type: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    by_pattern: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for r in rows:
        by_type[str(r["next_twin_type"])].append(r)
        by_pattern[str(r["pattern"])].append(r)

    X_dist = Counter(int(r["foreign_count"]) for r in rows)

    print("=" * 132)
    print(f"ZWISCHENVEKTOR-TEST BIS n = {args.limit}")
    print("=" * 132)
    print(f"Anzahl dichter Vierlinge: {len(rows)}")
    print(f"Mittelwert X = #nicht-zielige Zwischenprimzahlen: {mean(X_vals):.4f}")
    print(f"Mittelwert Gap zum nächsten Zwilling: {mean(gap_vals):.4f}")
    print()

    print("Verteilung von X")
    print(f"{'X':>6} {'Anzahl':>10} {'Anteil':>10}")
    print("-" * 30)
    for x in sorted(X_dist):
        print(f"{x:>6} {X_dist[x]:>10} {100*X_dist[x]/len(rows):>9.2f}%")

    print("\n" + "=" * 80)
    print("NACH ZWILLINGSTYP")
    print("=" * 80)
    print(f"{'Typ':>8} {'N':>8} {'<X>':>10} {'X=0':>10} {'<Gap>':>10}")
    print("-" * 50)
    for t in sorted(by_type):
        sub = by_type[t]
        x0 = sum(1 for r in sub if int(r["foreign_count"]) == 0)
        print(
            f"{t:>8} {len(sub):>8} {mean(float(r['foreign_count']) for r in sub):>10.4f} "
            f"{100*x0/len(sub):>9.2f}% {mean(float(r['gap_to_twin']) for r in sub):>10.4f}"
        )

    print("\n" + "=" * 80)
    print("NACH VIERLINGSMUSTER")
    print("=" * 80)
    print(f"{'Muster':>20} {'N':>8} {'<X>':>10} {'X=0':>10} {'<Gap>':>10}")
    print("-" * 60)
    for p in sorted(by_pattern):
        sub = by_pattern[p]
        x0 = sum(1 for r in sub if int(r["foreign_count"]) == 0)
        print(
            f"{p:>20} {len(sub):>8} {mean(float(r['foreign_count']) for r in sub):>10.4f} "
            f"{100*x0/len(sub):>9.2f}% {mean(float(r['gap_to_twin']) for r in sub):>10.4f}"
        )

    print("\n" + "=" * 80)
    print("KORRELATIONEN")
    print("=" * 80)
    print(f"corr(X, B_pi)   = {pearson_corr(X_vals, B_vals):.6f}")
    print(f"corr(X, T_star) = {pearson_corr(X_vals, T_vals):.6f}")
    print(f"corr(X, Gap)    = {pearson_corr(X_vals, gap_vals):.6f}")

    clean = [r for r in rows if int(r["foreign_count"]) == 0]
    messy = [r for r in rows if int(r["foreign_count"]) > 0]

    print("\n" + "=" * 80)
    print("SAUBERE VS. NICHT-SAUBERE FÄLLE")
    print("=" * 80)
    print(f"{'Gruppe':>14} {'N':>8} {'<B_pi>':>12} {'<T*>':>12} {'<Gap>':>12}")
    print("-" * 62)
    for name, sub in [("X=0", clean), ("X>0", messy)]:
        if not sub:
            continue
        print(
            f"{name:>14} {len(sub):>8} {mean(float(r['B_pi']) for r in sub):>12.6f} "
            f"{mean(float(r['T_star']) for r in sub):>12.6f} {mean(float(r['gap_to_twin']) for r in sub):>12.6f}"
        )

    print("\n" + "=" * 132)
    print("BEISPIELE MIT X = 0")
    print("=" * 132)
    print(f"{'Start':>8} {'Muster':>20} {'Zwilling':>14} {'Typ':>6} {'Gap':>6} {'X':>4} {'Counts':>22} {'B_pi':>12} {'T*':>12}")
    print("-" * 118)
    for r in clean[: args.show]:
        print(
            f"{int(r['start']):>8} {str(r['pattern']):>20} "
            f"{str(r['next_twin']):>14} {str(r['next_twin_type']):>6} "
            f"{int(r['gap_to_twin']):>6} {int(r['foreign_count']):>4} "
            f"{fmt_counts(r['family_counts']):>22} {float(r['B_pi']):>12.6f} {float(r['T_star']):>12.6f}"
        )

    print("\n" + "=" * 132)
    print("BEISPIELE MIT GROSSEM X")
    print("=" * 132)
    print(f"{'Start':>8} {'Muster':>20} {'Zwilling':>14} {'Typ':>6} {'Gap':>6} {'X':>4} {'Counts':>22} {'B_pi':>12} {'T*':>12}")
    print("-" * 118)
    bigX = sorted(rows, key=lambda r: (int(r["foreign_count"]), int(r["gap_to_twin"])), reverse=True)
    for r in bigX[: args.show]:
        print(
            f"{int(r['start']):>8} {str(r['pattern']):>20} "
            f"{str(r['next_twin']):>14} {str(r['next_twin_type']):>6} "
            f"{int(r['gap_to_twin']):>6} {int(r['foreign_count']):>4} "
            f"{fmt_counts(r['family_counts']):>22} {float(r['B_pi']):>12.6f} {float(r['T_star']):>12.6f}"
        )


if __name__ == "__main__":
    main()