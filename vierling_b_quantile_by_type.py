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


def quantile_groups(rows: Sequence[Dict[str, object]], num_groups: int) -> List[List[Dict[str, object]]]:
    rows_sorted = sorted(rows, key=lambda r: float(r["B_log"]))
    n = len(rows_sorted)
    groups: List[List[Dict[str, object]]] = []
    for i in range(num_groups):
        a = (i * n) // num_groups
        b = ((i + 1) * n) // num_groups
        groups.append(rows_sorted[a:b])
    return groups


def summarize_quantile(rows: Sequence[Dict[str, object]], pattern: str, label: str) -> Dict[str, object]:
    if not rows:
        return {
            "pattern": pattern,
            "quantile": label,
            "N": 0,
            "B_min": None,
            "B_max": None,
            "B_mean": None,
            "beta5_mean": None,
            "extendable_rate": None,
            "delta_mean": None,
            "sigma_mean": None,
        }

    return {
        "pattern": pattern,
        "quantile": label,
        "N": len(rows),
        "B_min": min(float(r["B_log"]) for r in rows),
        "B_max": max(float(r["B_log"]) for r in rows),
        "B_mean": mean(float(r["B_log"]) for r in rows),
        "beta5_mean": mean(float(r["beta5"]) for r in rows),
        "extendable_rate": mean(int(r["extendable"]) for r in rows),
        "delta_mean": mean(float(r["delta"]) for r in rows),
        "sigma_mean": mean(float(r["sigma"]) for r in rows),
    }


def print_quantile_table(rows: Sequence[Dict[str, object]]) -> None:
    print("=" * 144)
    print("GETRENNTER QUANTILTEST JE MUSTERTYP")
    print("=" * 144)
    header = (
        f"{'Typ':>20} {'Q':>4} {'N':>6} {'B_min':>12} {'B_max':>12} "
        f"{'<B>':>12} {'<β5>':>12} {'Erg.-Rate':>12} {'<Δ>':>8} {'<Σ>':>8}"
    )
    print(header)
    print("-" * len(header))
    for row in rows:
        if row["N"] == 0:
            continue
        print(
            f"{row['pattern']:>20} {row['quantile']:>4} {int(row['N']):>6} "
            f"{float(row['B_min']):>12.6f} {float(row['B_max']):>12.6f} "
            f"{float(row['B_mean']):>12.6f} {float(row['beta5_mean']):>12.6f} "
            f"{float(row['extendable_rate'])*100:>11.2f}% "
            f"{float(row['delta_mean']):>8.3f} {float(row['sigma_mean']):>8.3f}"
        )


def print_compact(rows: Sequence[Dict[str, object]], quantiles: int) -> None:
    print("\n" + "=" * 104)
    print("KOMPAKT: Q1 -> Q5 JE MUSTERTYP")
    print("=" * 104)
    header = f"{'Typ':>20} {'Q1-Rate':>10} {'Q5-Rate':>10} {'Differenz':>12} {'Monotonie grob':>18}"
    print(header)
    print("-" * len(header))

    patterns = sorted(set(str(r["pattern"]) for r in rows))
    for pattern in patterns:
        sub = [r for r in rows if str(r["pattern"]) == pattern]
        q1 = [r for r in sub if str(r["quantile"]) == "Q1"][0]
        q5 = [r for r in sub if str(r["quantile"]) == f"Q{quantiles}"][0]

        rates = [float(r["extendable_rate"]) for r in sorted(sub, key=lambda x: int(str(x["quantile"])[1:]))]
        monotone_non_decreasing = all(rates[i] <= rates[i + 1] + 1e-15 for i in range(len(rates) - 1))
        mono_text = "ja" if monotone_non_decreasing else "nein"

        diff = float(q5["extendable_rate"]) - float(q1["extendable_rate"])
        print(
            f"{pattern:>20} {float(q1['extendable_rate'])*100:>9.2f}% "
            f"{float(q5['extendable_rate'])*100:>9.2f}% {diff*100:>11.2f}% {mono_text:>18}"
        )


def print_examples(rows: Sequence[Dict[str, object]], title: str, limit: int = 12) -> None:
    print("\n" + "=" * 116)
    print(title)
    print("=" * 116)
    if not rows:
        print("Keine Daten.")
        return
    header = f"{'Start':>8} {'Typ':>20} {'B':>12} {'Δ':>4} {'Σ':>4} {'ergänzbar':>10}"
    print(header)
    print("-" * len(header))
    for r in rows[:limit]:
        print(
            f"{int(r['start']):>8} {str(r['coarse_pattern']):>20} "
            f"{float(r['B_log']):>12.6f} {int(r['delta']):>4} {int(r['sigma']):>4} {int(r['extendable']):>10}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Getrennter Quantiltest nach B_pi(Q4) für die beiden groben Vierlingsmuster."
    )
    parser.add_argument("--limit", type=int, default=1_000_000, help="obere Grenze")
    parser.add_argument("--pivot", type=int, default=11, help="Pivot-Primzahl π")
    parser.add_argument("--quantiles", type=int, default=5, help="Anzahl Quantile je Typ")
    parser.add_argument("--show", type=int, default=12, help="Beispielanzahl")
    args = parser.parse_args()

    sieve = PrimeSieve(args.limit + 12)
    quads = dense_quadruplets(args.limit, sieve)
    rows = [quadruplet_features(q, args.pivot, sieve) for q in quads]

    patterns = sorted(set(str(r["coarse_pattern"]) for r in rows))
    results: List[Dict[str, object]] = []

    print(f"Pivot π = {args.pivot}")
    print(f"Anzahl dichter Vierlinge: {len(rows)}")
    print(f"Gesamt-Ergänzungsrate: {100 * mean(int(r['extendable']) for r in rows):.2f}%\n")

    for pattern in patterns:
        sub = [r for r in rows if str(r["coarse_pattern"]) == pattern]
        groups = quantile_groups(sub, args.quantiles)
        for i, g in enumerate(groups):
            results.append(summarize_quantile(g, pattern=pattern, label=f"Q{i+1}"))

    print_quantile_table(results)
    print_compact(results, args.quantiles)

    for pattern in patterns:
        sub = [r for r in rows if str(r["coarse_pattern"]) == pattern]
        sub_sorted = sorted(sub, key=lambda r: float(r["B_log"]), reverse=True)
        print_examples(sub_sorted, f"BEISPIELE: {pattern} mit hohem B", limit=args.show)
        print_examples(list(reversed(sub_sorted[-args.show:])), f"BEISPIELE: {pattern} mit niedrigem B", limit=args.show)


if __name__ == "__main__":
    main()