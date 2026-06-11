from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

MOD12_CLASSES = {1: "E", 5: "A", 7: "B", 11: "C"}


@dataclass
class StepData:
    k: int
    shell_prime: int
    shell_primes: Tuple[int, ...]
    shell_part: int
    rest: int
    rest_two: int
    rest_three: int
    E: int
    A: int
    B: int
    C: int
    L_E: float
    sign_E: int
    L_ABC: float
    sign_ABC: int
    L_R: float
    sign_R: int


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
    n_mod_12: int
    mode_E: ModeSummary
    mode_ABC: ModeSummary
    mode_R: ModeSummary
    steps: List[StepData]


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


def shell_part_from_factors(factors: Dict[int, int], shell_primes: Sequence[int]) -> int:
    part = 1
    for p in shell_primes:
        exp = factors.get(p, 0)
        if exp:
            part *= p**exp
    return part


def split_rest_classes(rest_factors: Dict[int, int]) -> Tuple[int, int, int, int, int, int]:
    rest_two = rest_three = 1
    E = A = B = C = 1
    for p, exp in rest_factors.items():
        if p == 2:
            rest_two *= p**exp
            continue
        if p == 3:
            rest_three *= p**exp
            continue
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
    return rest_two, rest_three, E, A, B, C


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


def restricted_iter(start: int, stop: int, coprime_to_6_only: bool) -> Iterator[int]:
    if coprime_to_6_only:
        for n in range(start, stop + 1):
            if math.gcd(n, 6) == 1:
                yield n
    else:
        yield from range(start, stop + 1)


def compute_mode_summary(steps: Sequence[StepData], mode: str) -> ModeSummary:
    if not steps:
        raise RuntimeError("keine Schrittdaten vorhanden")

    values: List[float] = []
    signs: List[int] = []

    for step in steps:
        if mode == "E":
            values.append(step.L_E)
            signs.append(step.sign_E)
        elif mode == "ABC":
            values.append(step.L_ABC)
            signs.append(step.sign_ABC)
        elif mode == "R":
            values.append(step.L_R)
            signs.append(step.sign_R)
        else:
            raise ValueError(f"unbekannter Modus: {mode}")

    raw_switch_prime: Optional[int] = None
    raw_switch_k: Optional[int] = None
    true_switch_prime: Optional[int] = None
    true_switch_k: Optional[int] = None

    prev_nonzero: Optional[int] = signs[0] if signs[0] in (-1, 1) else None
    for i in range(1, len(signs)):
        prev_sign = signs[i - 1]
        current_sign = signs[i]

        if raw_switch_prime is None and current_sign != prev_sign and not (current_sign == 0 and prev_sign == 0):
            raw_switch_prime = steps[i].shell_prime
            raw_switch_k = steps[i].k

        if current_sign in (-1, 1):
            if prev_nonzero is not None and true_switch_prime is None and current_sign != prev_nonzero:
                true_switch_prime = steps[i].shell_prime
                true_switch_k = steps[i].k
            prev_nonzero = current_sign

    best_abs = float("inf")
    best_prime: Optional[int] = None
    best_k: Optional[int] = None

    best_nt_abs = float("inf")
    best_nt_prime: Optional[int] = None
    best_nt_k: Optional[int] = None

    for step, value in zip(steps, values):
        if abs(value) < best_abs:
            best_abs = abs(value)
            best_prime = step.shell_prime
            best_k = step.k

        if mode == "E":
            nontrivial = step.shell_part > 1 and step.E > 1
        elif mode == "ABC":
            nontrivial = step.shell_part > 1 and (step.A * step.B * step.C) > 1
        else:
            nontrivial = step.shell_part > 1 and step.rest > 1

        if nontrivial and abs(value) < best_nt_abs:
            best_nt_abs = abs(value)
            best_nt_prime = step.shell_prime
            best_nt_k = step.k

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
        nontrivial_balance_L_abs=(best_nt_abs if best_nt_prime is not None else None),
        initial_sign=signs[0],
        final_sign=signs[-1],
    )


def analyse_number(n: int, spf: SPFTable, shell_primes: Sequence[int]) -> NumberAnalysis:
    factors = spf.factor(n)
    steps: List[StepData] = []

    for k in range(1, len(shell_primes) + 1):
        current_shell = tuple(shell_primes[:k])
        shell_prime = current_shell[-1]

        shell_part = shell_part_from_factors(factors, current_shell)
        rest = n // shell_part
        rest_factors = {p: e for p, e in factors.items() if p not in current_shell}

        rest_two, rest_three, E, A, B, C = split_rest_classes(rest_factors)
        L_E = mode_value("E", shell_part, rest, E, A, B, C)
        L_ABC = mode_value("ABC", shell_part, rest, E, A, B, C)
        L_R = mode_value("R", shell_part, rest, E, A, B, C)

        steps.append(
            StepData(
                k=k,
                shell_prime=shell_prime,
                shell_primes=current_shell,
                shell_part=shell_part,
                rest=rest,
                rest_two=rest_two,
                rest_three=rest_three,
                E=E,
                A=A,
                B=B,
                C=C,
                L_E=L_E,
                sign_E=sign_of(L_E),
                L_ABC=L_ABC,
                sign_ABC=sign_of(L_ABC),
                L_R=L_R,
                sign_R=sign_of(L_R),
            )
        )

    return NumberAnalysis(
        n=n,
        n_mod_12=n % 12,
        mode_E=compute_mode_summary(steps, "E"),
        mode_ABC=compute_mode_summary(steps, "ABC"),
        mode_R=compute_mode_summary(steps, "R"),
        steps=steps,
    )


def analyse_range(start: int, stop: int, max_k: int, *, coprime_to_6_only: bool = True) -> List[NumberAnalysis]:
    if start < 1 or stop < start:
        raise ValueError("ungültiger Bereich")
    spf = SPFTable(stop)
    primes = first_primes(max_k)
    return [analyse_number(n, spf, primes) for n in restricted_iter(start, stop, coprime_to_6_only)]


def _print_mode_block(label: str, mode: ModeSummary) -> None:
    print(f"  [{label}]")
    print(f"    rohe Wechselprimzahl      : {mode.raw_switch_prime}")
    print(f"    rohe Wechsel-k           : {mode.raw_switch_k}")
    print(f"    echte Wechselprimzahl    : {mode.true_switch_prime}")
    print(f"    echte Wechsel-k          : {mode.true_switch_k}")
    print(f"    Balanceprimzahl          : {mode.balance_prime}")
    print(f"    Balance-k                : {mode.balance_k}")
    print(f"    |L|_min                  : {mode.balance_L_abs}")
    print(f"    nt. Balanceprimzahl      : {mode.nontrivial_balance_prime}")
    print(f"    nt. Balance-k            : {mode.nontrivial_balance_k}")
    print(f"    nt. |L|_min              : {mode.nontrivial_balance_L_abs}")
    print(f"    Anfangszeichen           : {mode.initial_sign}")
    print(f"    Endzeichen               : {mode.final_sign}")


def print_analysis(analysis: NumberAnalysis, show_steps: bool = False) -> None:
    print(f"n = {analysis.n}")
    print(f"  n mod 12        : {analysis.n_mod_12}")
    _print_mode_block("E", analysis.mode_E)
    _print_mode_block("ABC", analysis.mode_ABC)
    _print_mode_block("R", analysis.mode_R)
    if show_steps:
        print("  Schritte:")
        for s in analysis.steps:
            shell_str = "{" + ",".join(map(str, s.shell_primes)) + "}"
            print(
                f"    k={s.k:2d}, p_k={s.shell_prime:2d}, shell={shell_str:<20} "
                f"S={s.shell_part:<8d} R={s.rest:<8d} R2={s.rest_two:<6d} R3={s.rest_three:<6d} "
                f"E={s.E:<8d} A={s.A:<8d} B={s.B:<8d} C={s.C:<8d} "
                f"L_E={s.L_E:+.6f} L_ABC={s.L_ABC:+.6f} L_R={s.L_R:+.6f}"
            )


def write_summary_csv(path: Path, analyses: Sequence[NumberAnalysis]) -> None:
    fieldnames = [
        "n",
        "n_mod_12",
        "E_raw_switch_prime",
        "E_true_switch_prime",
        "E_balance_prime",
        "E_nontrivial_balance_prime",
        "ABC_raw_switch_prime",
        "ABC_true_switch_prime",
        "ABC_balance_prime",
        "ABC_nontrivial_balance_prime",
        "R_raw_switch_prime",
        "R_true_switch_prime",
        "R_balance_prime",
        "R_nontrivial_balance_prime",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for a in analyses:
            writer.writerow(
                {
                    "n": a.n,
                    "n_mod_12": a.n_mod_12,
                    "E_raw_switch_prime": a.mode_E.raw_switch_prime,
                    "E_true_switch_prime": a.mode_E.true_switch_prime,
                    "E_balance_prime": a.mode_E.balance_prime,
                    "E_nontrivial_balance_prime": a.mode_E.nontrivial_balance_prime,
                    "ABC_raw_switch_prime": a.mode_ABC.raw_switch_prime,
                    "ABC_true_switch_prime": a.mode_ABC.true_switch_prime,
                    "ABC_balance_prime": a.mode_ABC.balance_prime,
                    "ABC_nontrivial_balance_prime": a.mode_ABC.nontrivial_balance_prime,
                    "R_raw_switch_prime": a.mode_R.raw_switch_prime,
                    "R_true_switch_prime": a.mode_R.true_switch_prime,
                    "R_balance_prime": a.mode_R.balance_prime,
                    "R_nontrivial_balance_prime": a.mode_R.nontrivial_balance_prime,
                }
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Berechnet Wechsel- und Balanceprimzahlen in den Modi E, ABC und R."
    )
    parser.add_argument("--start", type=int, default=1, help="Startwert")
    parser.add_argument("--stop", type=int, default=100, help="Endwert")
    parser.add_argument("--max-k", type=int, default=10, help="Anzahl der ersten Primzahlen in der Schalenfolge")
    parser.add_argument(
        "--all-numbers",
        action="store_true",
        help="analysiert alle Zahlen im Bereich; Restanteile von 2 und 3 werden separat mitgefuehrt",
    )
    parser.add_argument("--show-steps", action="store_true", help="Gibt die gesamte Schalenfolge pro Zahl aus")
    parser.add_argument("--csv", type=Path, default=None, help="CSV-Datei für die Zusammenfassung")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    analyses = analyse_range(
        args.start,
        args.stop,
        args.max_k,
        coprime_to_6_only=not args.all_numbers,
    )

    for a in analyses:
        print_analysis(a, show_steps=args.show_steps)
        print()

    if args.csv is not None:
        write_summary_csv(args.csv, analyses)
        print(f"\nCSV gespeichert: {args.csv}")


if __name__ == "__main__":
    main()
