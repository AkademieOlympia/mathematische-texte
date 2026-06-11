from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

MOD12_CLASSES = {1: "E", 5: "A", 7: "B", 11: "C"}


@dataclass
class ModeSummary:
    raw_switch_prime: Optional[int]
    raw_switch_k: Optional[int]
    true_switch_prime: Optional[int]
    true_switch_k: Optional[int]
    balance_prime: Optional[int]
    balance_k: Optional[int]
    balance_L_abs: Optional[float]
    nontrivial_balance_prime: Optional[int]
    nontrivial_balance_k: Optional[int]
    nontrivial_balance_L_abs: Optional[float]
    initial_sign: int
    final_sign: int


@dataclass
class NumberAnalysis:
    n: int
    factor_data: Tuple[Tuple[int, int], ...]  # ((p,a),(q,b))
    factor_classes: Tuple[str, ...]
    mode_E: ModeSummary
    mode_ABC: ModeSummary
    mode_R: ModeSummary


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


def class_of_prime(p: int) -> str:
    if p <= 3:
        return "shell"
    residue = p % 12
    if residue not in MOD12_CLASSES:
        raise ValueError(f"Primzahl {p} hat unerwartete Restklasse {residue} mod 12")
    return MOD12_CLASSES[residue]


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


def compute_mode_summary(
    shell_primes: Sequence[int],
    shell_parts: Sequence[int],
    rests: Sequence[int],
    Es: Sequence[int],
    As: Sequence[int],
    Bs: Sequence[int],
    Cs: Sequence[int],
    mode: str,
) -> ModeSummary:
    values: List[float] = []
    signs: List[int] = []

    for shell_part, rest, E, A, B, C in zip(shell_parts, rests, Es, As, Bs, Cs):
        L = mode_value(mode, shell_part, rest, E, A, B, C)
        values.append(L)
        signs.append(sign_of(L))

    raw_switch_prime: Optional[int] = None
    raw_switch_k: Optional[int] = None
    true_switch_prime: Optional[int] = None
    true_switch_k: Optional[int] = None

    prev_nonzero: Optional[int] = None
    for i in range(1, len(signs)):
        prev_sign = signs[i - 1]
        current_sign = signs[i]

        if raw_switch_prime is None and current_sign != prev_sign and not (current_sign == 0 and prev_sign == 0):
            raw_switch_prime = shell_primes[i]
            raw_switch_k = i + 1

        if current_sign in (-1, 1):
            if prev_nonzero is not None and true_switch_prime is None and current_sign != prev_nonzero:
                true_switch_prime = shell_primes[i]
                true_switch_k = i + 1
            prev_nonzero = current_sign

    best_abs = float("inf")
    best_prime: Optional[int] = None
    best_k: Optional[int] = None

    best_nt_abs = float("inf")
    best_nt_prime: Optional[int] = None
    best_nt_k: Optional[int] = None

    for i, (L, shell_part, rest, E, A, B, C) in enumerate(zip(values, shell_parts, rests, Es, As, Bs, Cs)):
        if abs(L) < best_abs:
            best_abs = abs(L)
            best_prime = shell_primes[i]
            best_k = i + 1

        if mode == "E":
            nontrivial = shell_part > 1 and E > 1
        elif mode == "ABC":
            nontrivial = shell_part > 1 and (A * B * C) > 1
        elif mode == "R":
            nontrivial = shell_part > 1 and rest > 1
        else:
            raise ValueError(f"unbekannter Modus: {mode}")

        if nontrivial and abs(L) < best_nt_abs:
            best_nt_abs = abs(L)
            best_nt_prime = shell_primes[i]
            best_nt_k = i + 1

    return ModeSummary(
        raw_switch_prime=raw_switch_prime,
        raw_switch_k=raw_switch_k,
        true_switch_prime=true_switch_prime,
        true_switch_k=true_switch_k,
        balance_prime=best_prime,
        balance_k=best_k,
        balance_L_abs=best_abs if best_prime is not None else None,
        nontrivial_balance_prime=best_nt_prime,
        nontrivial_balance_k=best_nt_k,
        nontrivial_balance_L_abs=best_nt_abs if best_nt_prime is not None else None,
        initial_sign=signs[0],
        final_sign=signs[-1],
    )


def analyse_from_factors(
    n: int,
    factors: Dict[int, int],
    shell_primes: Sequence[int],
) -> NumberAnalysis:
    factor_data = tuple(sorted(factors.items()))
    factor_classes = tuple(class_of_prime(p) for p, _ in factor_data)

    shell_parts: List[int] = []
    rests: List[int] = []
    Es: List[int] = []
    As: List[int] = []
    Bs: List[int] = []
    Cs: List[int] = []

    for k in range(1, len(shell_primes) + 1):
        current_shell = tuple(shell_primes[:k])
        shell_part = shell_part_from_factors(factors, current_shell)
        rest = n // shell_part
        rest_factors = {p: e for p, e in factors.items() if p not in current_shell}
        E, A, B, C = split_rest_classes(rest_factors)

        shell_parts.append(shell_part)
        rests.append(rest)
        Es.append(E)
        As.append(A)
        Bs.append(B)
        Cs.append(C)

    return NumberAnalysis(
        n=n,
        factor_data=factor_data,
        factor_classes=factor_classes,
        mode_E=compute_mode_summary(shell_primes, shell_parts, rests, Es, As, Bs, Cs, "E"),
        mode_ABC=compute_mode_summary(shell_primes, shell_parts, rests, Es, As, Bs, Cs, "ABC"),
        mode_R=compute_mode_summary(shell_primes, shell_parts, rests, Es, As, Bs, Cs, "R"),
    )


def analyse_number(n: int, spf: SPFTable, shell_primes: Sequence[int]) -> NumberAnalysis:
    return analyse_from_factors(n, spf.factor(n), shell_primes)


def mode_for_family(family: str) -> str:
    if family == "E*E":
        return "E"
    if family == "ABC*ABC":
        return "ABC"
    if family == "ABC*E":
        return "R"
    raise ValueError(f"unbekannte Familie: {family}")


def expected_small_large(factor_data: Tuple[Tuple[int, int], ...]) -> Tuple[int, int]:
    primes = sorted(p for p, _ in factor_data)
    return primes[0], primes[-1]


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


def generate_paqb_cases(prime_limit: int, max_exp: int) -> List[Tuple[str, int, int, int, int]]:
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


def get_mode(analysis: NumberAnalysis, mode_name: str) -> ModeSummary:
    if mode_name == "E":
        return analysis.mode_E
    if mode_name == "ABC":
        return analysis.mode_ABC
    if mode_name == "R":
        return analysis.mode_R
    raise ValueError(f"unbekannter Modus: {mode_name}")


def evaluate_cases(prime_limit: int, max_exp: int, max_k: int) -> List[Dict[str, object]]:
    cases = generate_paqb_cases(prime_limit, max_exp)
    if not cases:
        return []

    shell_primes = first_primes(max_k)

    rows: List[Dict[str, object]] = []
    for family, p, a, q, b in cases:
        n = (p**a) * (q**b)
        analysis = analyse_from_factors(n, {p: a, q: b}, shell_primes)
        expected_balance, expected_switch = expected_small_large(analysis.factor_data)

        for mode_name in ("E", "ABC", "R"):
            mode = get_mode(analysis, mode_name)
            rows.append(
                {
                    "family": family,
                    "mode": mode_name,
                    "n": n,
                    "p": p,
                    "a": a,
                    "q": q,
                    "b": b,
                    "factor_classes": ",".join(analysis.factor_classes),
                    "expected_mode": mode_for_family(family),
                    "true_switch_prime": mode.true_switch_prime,
                    "nontrivial_balance_prime": mode.nontrivial_balance_prime,
                    "has_true_switch": mode.true_switch_prime is not None,
                    "has_nontrivial_balance": mode.nontrivial_balance_prime is not None,
                    "balance_is_smallest": mode.nontrivial_balance_prime == expected_balance,
                    "switch_is_largest": mode.true_switch_prime == expected_switch,
                    "both_pattern": (
                        mode.nontrivial_balance_prime == expected_balance
                        and mode.true_switch_prime == expected_switch
                    ),
                }
            )
    return rows


def summarize_rows(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[str, str], List[Dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["family"]), str(row["mode"]))].append(row)

    out: List[Dict[str, object]] = []
    for (family, mode), block in sorted(grouped.items()):
        total = len(block)
        true_switch_count = sum(1 for r in block if r["has_true_switch"])
        nt_balance_count = sum(1 for r in block if r["has_nontrivial_balance"])
        balance_smallest_count = sum(1 for r in block if r["balance_is_smallest"])
        switch_largest_count = sum(1 for r in block if r["switch_is_largest"])
        both_count = sum(1 for r in block if r["both_pattern"])

        out.append(
            {
                "family": family,
                "mode": mode,
                "total": total,
                "true_switch_rate": true_switch_count / total if total else 0.0,
                "nt_balance_rate": nt_balance_count / total if total else 0.0,
                "balance_smallest_rate": balance_smallest_count / total if total else 0.0,
                "switch_largest_rate": switch_largest_count / total if total else 0.0,
                "both_rate": both_count / total if total else 0.0,
            }
        )
    return out


def print_summary_table(summary: Sequence[Dict[str, object]]) -> None:
    print("=" * 120)
    print("P^A Q^B-TEST: STABILITÄT VON 'BALANCE UNTEN, WECHSEL OBEN'")
    print("=" * 120)
    header = (
        f"{'Familie':<12} {'Modus':<6} {'N':>6} {'echter W.':>10} {'nt. Bal.':>10} "
        f"{'Bal.=klein':>11} {'Wech.=groß':>12} {'beides':>10}"
    )
    print(header)
    print("-" * len(header))
    for row in summary:
        print(
            f"{row['family']:<12} {row['mode']:<6} "
            f"{row['total']:>6} "
            f"{row['true_switch_rate']*100:>9.1f}% "
            f"{row['nt_balance_rate']*100:>9.1f}% "
            f"{row['balance_smallest_rate']*100:>10.1f}% "
            f"{row['switch_largest_rate']*100:>11.1f}% "
            f"{row['both_rate']*100:>9.1f}%"
        )


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
        description="Systematische Tests für Schalenwechsel auf Zahlen der Form p^a q^b."
    )
    parser.add_argument("--prime-limit", type=int, default=100, help="obere Primgrenze für p und q")
    parser.add_argument("--max-exp", type=int, default=3, help="maximaler Exponent a,b")
    parser.add_argument("--max-k", type=int, default=25, help="Anzahl der ersten Primzahlen in der Schalenfolge")
    parser.add_argument("--csv", type=Path, default=None, help="CSV-Datei für Rohdaten")
    parser.add_argument("--summary-csv", type=Path, default=None, help="CSV-Datei für Zusammenfassung")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    rows = evaluate_cases(args.prime_limit, args.max_exp, args.max_k)
    summary = summarize_rows(rows)
    print_summary_table(summary)

    if args.csv is not None:
        write_csv(args.csv, rows)
        print(f"\nRohdaten gespeichert: {args.csv}")

    if args.summary_csv is not None:
        write_csv(args.summary_csv, summary)
        print(f"Zusammenfassung gespeichert: {args.summary_csv}")


if __name__ == "__main__":
    main()
