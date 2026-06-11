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


def analyze_quad(
    quad: Tuple[int, int, int, int],
    sieve: PrimeSieve,
    pivot: int,
) -> Dict[str, object]:
    pattern = quadruplet_pattern(quad)
    q_end = quad[-1]
    nxt = twin_primes_after(q_end + 1, sieve)
    ttype = twin_type(nxt)
    gap = nxt[0] - q_end

    B = B_pi_of_quad(quad, pivot)
    active = int(pattern == "EA-BC-BC-EA")
    T_star = B + 6.0 * active

    y_ab = 1 if ttype == "AB" else 0
    y_ce = 1 if ttype == "CE" else 0

    return {
        "start": quad[0],
        "quad": quad,
        "pattern": pattern,
        "next_twin": nxt,
        "next_twin_type": ttype,
        "gap": gap,
        "B_pi": B,
        "T_star": T_star,
        "active_channel": active,
        "is_AB": y_ab,
        "is_CE": y_ce,
    }


def print_examples(rows: Sequence[Dict[str, object]], title: str, limit: int = 15) -> None:
    print("\n" + "=" * 118)
    print(title)
    print("=" * 118)
    header = f"{'Start':>8} {'Muster':>20} {'Zwilling':>14} {'Typ':>6} {'Gap':>6} {'B_pi':>12} {'T*':>12}"
    print(header)
    print("-" * len(header))
    for r in rows[:limit]:
        print(
            f"{int(r['start']):>8} {str(r['pattern']):>20} {str(r['next_twin']):>14} "
            f"{str(r['next_twin_type']):>6} {int(r['gap']):>6} "
            f"{float(r['B_pi']):>12.6f} {float(r['T_star']):>12.6f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Testet, ob Vierlingsstruktur den Typ des nächsten Primzahlzwillings (AB oder CE) vorhersagt."
    )
    parser.add_argument("--limit", type=int, default=2_000_000, help="obere Grenze")
    parser.add_argument("--pivot", type=int, default=11, help="Pivot π")
    parser.add_argument("--show", type=int, default=15, help="Beispielanzahl")
    args = parser.parse_args()

    sieve = PrimeSieve(args.limit + 1000)
    quads = dense_quadruplets(args.limit, sieve)
    rows = [analyze_quad(q, sieve, args.pivot) for q in quads]

    by_type: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    by_pattern: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for r in rows:
        by_type[str(r["next_twin_type"])].append(r)
        by_pattern[str(r["pattern"])].append(r)

    B = [float(r["B_pi"]) for r in rows]
    T = [float(r["T_star"]) for r in rows]
    G = [float(r["gap"]) for r in rows]
    Yab = [float(r["is_AB"]) for r in rows]
    Yce = [float(r["is_CE"]) for r in rows]

    print("=" * 132)
    print(f"TYPTEST DES NÄCHSTEN PRIMZAHLZWILLINGS BIS n = {args.limit}")
    print("=" * 132)
    print(f"Anzahl dichter Vierlinge: {len(rows)}")
    print()

    type_counts = Counter(str(r["next_twin_type"]) for r in rows)
    print("Verteilung des nächsten Zwillingstyps")
    print(f"{'Typ':>8} {'Anzahl':>10} {'Anteil':>10}")
    print("-" * 30)
    for t in sorted(type_counts):
        print(f"{t:>8} {type_counts[t]:>10} {100*type_counts[t]/len(rows):>9.2f}%")

    print("\n" + "=" * 84)
    print("KONTINGENZTAFEL: VIERLINGSMUSTER × NÄCHSTER ZWILLINGSTYP")
    print("=" * 84)
    patterns = sorted(by_pattern)
    print(f"{'Muster':>20} {'AB':>8} {'CE':>8} {'AB-Anteil':>12}")
    print("-" * 52)
    for p in patterns:
        sub = by_pattern[p]
        ab = sum(1 for r in sub if str(r['next_twin_type']) == "AB")
        ce = sum(1 for r in sub if str(r['next_twin_type']) == "CE")
        print(f"{p:>20} {ab:>8} {ce:>8} {100*ab/len(sub):>11.2f}%")

    print("\n" + "=" * 84)
    print("MITTELWERTE NACH ZWILLINGSTYP")
    print("=" * 84)
    print(f"{'Typ':>8} {'N':>8} {'<Gap>':>10} {'<B_pi>':>12} {'<T*>':>12} {'aktiver Kanal':>14}")
    print("-" * 68)
    for t in sorted(by_type):
        sub = by_type[t]
        print(
            f"{t:>8} {len(sub):>8} {mean(float(r['gap']) for r in sub):>10.4f} "
            f"{mean(float(r['B_pi']) for r in sub):>12.6f} "
            f"{mean(float(r['T_star']) for r in sub):>12.6f} "
            f"{100*mean(int(r['active_channel']) for r in sub):>13.2f}%"
        )

    print("\n" + "=" * 84)
    print("KORRELATIONEN MIT ZWILLINGSTYP AB")
    print("=" * 84)
    print(f"corr(B_pi, AB)  = {pearson_corr(B, Yab):.6f}")
    print(f"corr(T*, AB)    = {pearson_corr(T, Yab):.6f}")
    print(f"corr(Gap, AB)   = {pearson_corr(G, Yab):.6f}")
    print(f"spearman(B, AB) = {spearman_corr(B, Yab):.6f}")
    print(f"spearman(T, AB) = {spearman_corr(T, Yab):.6f}")

    print("\n" + "=" * 84)
    print("KORRELATIONEN MIT ZWILLINGSTYP CE")
    print("=" * 84)
    print(f"corr(B_pi, CE)  = {pearson_corr(B, Yce):.6f}")
    print(f"corr(T*, CE)    = {pearson_corr(T, Yce):.6f}")
    print(f"corr(Gap, CE)   = {pearson_corr(G, Yce):.6f}")
    print(f"spearman(B, CE) = {spearman_corr(B, Yce):.6f}")
    print(f"spearman(T, CE) = {spearman_corr(T, Yce):.6f}")

    ab_rows = [r for r in rows if str(r["next_twin_type"]) == "AB"]
    ce_rows = [r for r in rows if str(r["next_twin_type"]) == "CE"]

    ab_sorted = sorted(ab_rows, key=lambda r: float(r["T_star"]), reverse=True)
    ce_sorted = sorted(ce_rows, key=lambda r: float(r["T_star"]), reverse=True)

    print_examples(ab_sorted, "BEISPIELE: NÄCHSTER ZWILLING VOM TYP AB", limit=args.show)
    print_examples(ce_sorted, "BEISPIELE: NÄCHSTER ZWILLING VOM TYP CE", limit=args.show)


if __name__ == "__main__":
    main()