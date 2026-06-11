from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

MOD12_CLASSES = {1: "E", 5: "A", 7: "B", 11: "C"}


class SPFTable:
    def __init__(self, limit: int) -> None:
        if limit < 2:
            raise ValueError("limit muss mindestens 2 sein")
        self.limit = limit
        self.spf = self._build(limit)

    @staticmethod
    def _build(limit: int) -> List[int]:
        spf = list(range(limit + 1))
        spf[0] = 0
        spf[1] = 1
        for p in range(2, int(limit**0.5) + 1):
            if spf[p] == p:
                for m in range(p * p, limit + 1, p):
                    if spf[m] == m:
                        spf[m] = p
        return spf

    def factor(self, n: int) -> Dict[int, int]:
        if not (1 <= n <= self.limit):
            raise ValueError(f"n muss zwischen 1 und {self.limit} liegen")
        fac: Dict[int, int] = {}
        while n > 1:
            p = self.spf[n]
            fac[p] = fac.get(p, 0) + 1
            n //= p
        return fac


def primes_up_to(limit: int) -> List[int]:
    if limit < 2:
        return []
    sieve = [True] * (limit + 1)
    sieve[0:2] = [False, False]
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            for m in range(p * p, limit + 1, p):
                sieve[m] = False
    return [p for p in range(2, limit + 1) if sieve[p]]


def first_primes(count: int) -> List[int]:
    if count < 1:
        raise ValueError("count muss positiv sein")
    out: List[int] = []
    x = 2
    while len(out) < count:
        is_prime = True
        root = int(x**0.5)
        for p in out:
            if p > root:
                break
            if x % p == 0:
                is_prime = False
                break
        if is_prime:
            out.append(x)
        x += 1 if x == 2 else 2
    return out


def class_of_prime(p: int) -> str:
    if p <= 3:
        return "shell"
    residue = p % 12
    if residue not in MOD12_CLASSES:
        raise ValueError(f"unerwartete Restklasse für Primzahl {p}")
    return MOD12_CLASSES[residue]


def families_for_pair(p: int, q: int) -> List[str]:
    cp = class_of_prime(p)
    cq = class_of_prime(q)
    if cp == "E" and cq == "E":
        return ["E*E"]
    if cp in ("A", "B", "C") and cq in ("A", "B", "C"):
        return ["ABC*ABC"]
    if (cp in ("A", "B", "C") and cq == "E") or (cp == "E" and cq in ("A", "B", "C")):
        return ["ABC*E"]
    return []


def expected_mode_for_family(family: str) -> str:
    if family == "E*E":
        return "E"
    if family == "ABC*ABC":
        return "ABC"
    if family == "ABC*E":
        return "R"
    raise ValueError(f"unbekannte Familie: {family}")


def shell_part_from_factors(factors: Dict[int, int], shell_primes: Sequence[int]) -> int:
    part = 1
    for p in shell_primes:
        exp = factors.get(p, 0)
        if exp:
            part *= p**exp
    return part


def split_rest_classes(rest_factors: Dict[int, int]) -> Tuple[int, int, int, int]:
    E = A = B = C = 1
    for p, exp in rest_factors.items():
        residue = p % 12
        value = p**exp
        label = MOD12_CLASSES.get(residue)
        if label == "E":
            E *= value
        elif label == "A":
            A *= value
        elif label == "B":
            B *= value
        elif label == "C":
            C *= value
        else:
            raise ValueError(
                f"Primzahl {p} im Rest hat Restklasse {residue} mod 12; "
                "das sollte nach Schalenabzug nicht passieren."
            )
    return E, A, B, C


def safe_log(x: int) -> float:
    return math.log(x) if x > 1 else 0.0


def sign_of(x: float, eps: float = 1e-15) -> int:
    if x > eps:
        return 1
    if x < -eps:
        return -1
    return 0


def mode_value(mode: str, shell_part: int, rest: int, E: int, A: int, B: int, C: int) -> float:
    if mode == "E":
        return safe_log(E) - safe_log(shell_part)
    if mode == "ABC":
        return safe_log(A * B * C) - safe_log(shell_part)
    if mode == "R":
        return safe_log(rest) - safe_log(shell_part)
    raise ValueError(f"unbekannter Modus: {mode}")


def compute_mode_data(
    n: int,
    factors: Dict[int, int],
    shell_primes: Sequence[int],
    mode: str,
) -> Dict[str, Optional[float]]:
    values: List[float] = []
    signs: List[int] = []

    best_nt_abs = float("inf")
    best_nt_prime: Optional[int] = None
    best_nt_k: Optional[int] = None
    true_switch_prime: Optional[int] = None
    true_switch_k: Optional[int] = None

    prev_nonzero: Optional[int] = None

    for k in range(1, len(shell_primes) + 1):
        current_shell = tuple(shell_primes[:k])
        shell_prime = current_shell[-1]

        shell_part = shell_part_from_factors(factors, current_shell)
        rest = n // shell_part
        rest_factors = {p: e for p, e in factors.items() if p not in current_shell}
        E, A, B, C = split_rest_classes(rest_factors)

        L = mode_value(mode, shell_part, rest, E, A, B, C)
        s = sign_of(L)
        values.append(L)
        signs.append(s)

        if s in (-1, 1):
            if prev_nonzero is not None and true_switch_prime is None and s != prev_nonzero:
                true_switch_prime = shell_prime
                true_switch_k = k
            prev_nonzero = s

        if mode == "E":
            nontrivial = shell_part > 1 and E > 1
        elif mode == "ABC":
            nontrivial = shell_part > 1 and (A * B * C) > 1
        else:
            nontrivial = shell_part > 1 and rest > 1

        if nontrivial and abs(L) < best_nt_abs:
            best_nt_abs = abs(L)
            best_nt_prime = shell_prime
            best_nt_k = k

    return {
        "true_switch_prime": true_switch_prime,
        "true_switch_k": true_switch_k,
        "nontrivial_balance_prime": best_nt_prime,
        "nontrivial_balance_k": best_nt_k,
        "nontrivial_balance_L_abs": (best_nt_abs if best_nt_prime is not None else None),
        "initial_sign": signs[0] if signs else None,
        "final_sign": signs[-1] if signs else None,
    }


def generate_cases(prime_limit: int, max_exp: int) -> List[Tuple[str, int, int, int, int]]:
    ps = [p for p in primes_up_to(prime_limit) if p > 3]
    cases: List[Tuple[str, int, int, int, int]] = []
    for i, p in enumerate(ps):
        for q in ps[i + 1 :]:
            fams = families_for_pair(p, q)
            if not fams:
                continue
            for a in range(1, max_exp + 1):
                for b in range(1, max_exp + 1):
                    for fam in fams:
                        cases.append((fam, p, a, q, b))
    return cases


def evaluate_counterexamples(prime_limit: int, max_exp: int, max_k: int) -> List[Dict[str, object]]:
    cases = generate_cases(prime_limit, max_exp)
    if not cases:
        return []

    shell_primes = first_primes(max_k)

    rows: List[Dict[str, object]] = []

    for family, p, a, q, b in cases:
        n = (p**a) * (q**b)
        factors = {p: a, q: b}
        mode = expected_mode_for_family(family)
        data = compute_mode_data(n, factors, shell_primes, mode)

        true_switch_prime = data["true_switch_prime"]
        if true_switch_prime is None:
            continue
        if true_switch_prime == q:
            continue

        rows.append(
            {
                "family": family,
                "mode": mode,
                "n": n,
                "p": p,
                "a": a,
                "q": q,
                "b": b,
                "factor_classes": f"{class_of_prime(p)},{class_of_prime(q)}",
                "true_switch_prime": true_switch_prime,
                "true_switch_k": data["true_switch_k"],
                "nontrivial_balance_prime": data["nontrivial_balance_prime"],
                "nontrivial_balance_k": data["nontrivial_balance_k"],
                "nontrivial_balance_L_abs": data["nontrivial_balance_L_abs"],
                "initial_sign": data["initial_sign"],
                "final_sign": data["final_sign"],
            }
        )
    return rows


def summarize_by_exponents(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[str, int, int], List[Dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["family"]), int(row["a"]), int(row["b"]))].append(row)

    out: List[Dict[str, object]] = []
    for (family, a, b), block in sorted(grouped.items()):
        out.append(
            {
                "family": family,
                "a": a,
                "b": b,
                "count": len(block),
            }
        )
    return out


def print_counterexamples(rows: Sequence[Dict[str, object]], limit: int = 50) -> None:
    print("=" * 130)
    print("GEGENBEISPIELE ZU: Wechsel = größerer Primfaktor")
    print("=" * 130)
    if not rows:
        print("Keine Gegenbeispiele gefunden.")
        return

    header = (
        f"{'Fam.':<10} {'n':>8} {'p':>4} {'a':>2} {'q':>4} {'b':>2} "
        f"{'Wechsel':>8} {'Balance':>8} {'|L|_min':>12} {'init':>5} {'final':>5}"
    )
    print(header)
    print("-" * len(header))
    for row in rows[:limit]:
        print(
            f"{row['family']:<10} {row['n']:>8} {row['p']:>4} {row['a']:>2} "
            f"{row['q']:>4} {row['b']:>2} "
            f"{str(row['true_switch_prime']):>8} "
            f"{str(row['nontrivial_balance_prime']):>8} "
            f"{str(round(float(row['nontrivial_balance_L_abs']), 6)) if row['nontrivial_balance_L_abs'] is not None else 'None':>12} "
            f"{str(row['initial_sign']):>5} {str(row['final_sign']):>5}"
        )
    if len(rows) > limit:
        print(f"... weitere {len(rows) - limit} Gegenbeispiele nicht angezeigt.")


def print_exponent_summary(rows: Sequence[Dict[str, object]]) -> None:
    summary = summarize_by_exponents(rows)
    print("\n" + "=" * 80)
    print("GRUPPIERUNG NACH EXPONENTENPAAREN (a,b)")
    print("=" * 80)
    if not summary:
        print("Keine Gegenbeispiele.")
        return

    header = f"{'Familie':<12} {'a':>3} {'b':>3} {'Anzahl Gegenbeispiele':>22}"
    print(header)
    print("-" * len(header))
    for row in summary:
        print(f"{row['family']:<12} {row['a']:>3} {row['b']:>3} {row['count']:>22}")


def write_csv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Filtert Gegenbeispiele zur Regel 'Wechsel = größerer Primfaktor' für p^a q^b."
    )
    parser.add_argument("--prime-limit", type=int, default=100, help="obere Primgrenze")
    parser.add_argument("--max-exp", type=int, default=3, help="maximaler Exponent")
    parser.add_argument("--max-k", type=int, default=25, help="Anzahl der Schalenprimzahlen")
    parser.add_argument("--show", type=int, default=50, help="wie viele Gegenbeispiele anzeigen")
    parser.add_argument("--csv", type=Path, default=None, help="CSV-Datei für Gegenbeispiele")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    rows = evaluate_counterexamples(args.prime_limit, args.max_exp, args.max_k)
    print_counterexamples(rows, limit=args.show)
    print_exponent_summary(rows)

    if args.csv is not None:
        write_csv(args.csv, rows)
        print(f"\nCSV gespeichert: {args.csv}")


if __name__ == "__main__":
    main()