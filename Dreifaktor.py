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
    true_switch_prime: Optional[int]
    true_switch_k: Optional[int]
    nontrivial_balance_prime: Optional[int]
    nontrivial_balance_k: Optional[int]
    nontrivial_balance_L_abs: Optional[float]
    initial_sign: int
    final_sign: int


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


def classify_triple(p: int, q: int, r: int) -> str:
    classes = sorted(class_of_prime(x) for x in (p, q, r))
    if classes == ["E", "E", "E"]:
        return "E*E*E"
    if all(c in ("A", "B", "C") for c in classes):
        return "ABC*ABC*ABC"
    if classes.count("E") == 1:
        return "ABC*ABC*E"
    if classes.count("E") == 2:
        return "ABC*E*E"
    raise ValueError(f"unbekannte Klassenkombination: {classes}")


def expected_mode_for_family(family: str) -> str:
    if family == "E*E*E":
        return "E"
    if family == "ABC*ABC*ABC":
        return "ABC"
    if family in ("ABC*ABC*E", "ABC*E*E"):
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


def compute_mode_summary(
    n: int,
    factors: Dict[int, int],
    shell_primes: Sequence[int],
    mode: str,
) -> ModeSummary:
    signs: List[int] = []
    true_switch_prime: Optional[int] = None
    true_switch_k: Optional[int] = None
    prev_nonzero: Optional[int] = None

    best_nt_abs = float("inf")
    best_nt_prime: Optional[int] = None
    best_nt_k: Optional[int] = None

    for k in range(1, len(shell_primes) + 1):
        current_shell = tuple(shell_primes[:k])
        shell_prime = current_shell[-1]
        shell_part = shell_part_from_factors(factors, current_shell)
        rest = n // shell_part
        rest_factors = {p: e for p, e in factors.items() if p not in current_shell}
        E, A, B, C = split_rest_classes(rest_factors)

        L = mode_value(mode, shell_part, rest, E, A, B, C)
        s = sign_of(L)
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
        elif mode == "R":
            nontrivial = shell_part > 1 and rest > 1
        else:
            raise ValueError(f"unbekannter Modus: {mode}")

        if nontrivial and abs(L) < best_nt_abs:
            best_nt_abs = abs(L)
            best_nt_prime = shell_prime
            best_nt_k = k

    return ModeSummary(
        true_switch_prime=true_switch_prime,
        true_switch_k=true_switch_k,
        nontrivial_balance_prime=best_nt_prime,
        nontrivial_balance_k=best_nt_k,
        nontrivial_balance_L_abs=(best_nt_abs if best_nt_prime is not None else None),
        initial_sign=signs[0] if signs else 0,
        final_sign=signs[-1] if signs else 0,
    )


def generate_triple_cases(prime_limit: int, pattern: str) -> List[Tuple[str, int, int, int, int, int, int, int, int]]:
    ps = [p for p in primes_up_to(prime_limit) if p > 3]
    cases: List[Tuple[str, int, int, int, int, int, int, int, int]] = []
    patterns = {
        "pqr": [(1, 1, 1)],
        "p2qr": [(2, 1, 1)],
        "pq2r": [(1, 2, 1)],
        "pqr2": [(1, 1, 2)],
        "all_small": [(1, 1, 1), (2, 1, 1), (1, 2, 1), (1, 1, 2)],
    }
    if pattern not in patterns:
        raise ValueError(f"unbekanntes Muster: {pattern}")

    for i, p in enumerate(ps):
        for j, q in enumerate(ps[i + 1 :], start=i + 1):
            for r in ps[j + 1 :]:
                family = classify_triple(p, q, r)
                for a, b, c in patterns[pattern]:
                    cases.append((family, p, a, q, b, r, c, p**a * q**b * r**c, 0))
    return cases


def evaluate_cases(prime_limit: int, pattern: str, max_k: int) -> List[Dict[str, object]]:
    cases = generate_triple_cases(prime_limit, pattern)
    if not cases:
        return []

    max_n = max(n for *_, n, _ in cases)
    spf = SPFTable(max_n)
    shell_primes = first_primes(max_k)
    rows: List[Dict[str, object]] = []

    for family, p, a, q, b, r, c, n, _ in cases:
        factors = spf.factor(n)
        mode = expected_mode_for_family(family)
        mode_data = compute_mode_summary(n, factors, shell_primes, mode)

        balance_expected = p
        switch_expected_largest_prime = r

        mass_p = a * math.log(p)
        mass_q = b * math.log(q)
        mass_r = c * math.log(r)
        mass_map = {p: mass_p, q: mass_q, r: mass_r}
        switch_expected_logmass = max(mass_map.items(), key=lambda kv: (kv[1], kv[0]))[0]

        rows.append(
            {
                "family": family,
                "mode": mode,
                "p": p,
                "a": a,
                "q": q,
                "b": b,
                "r": r,
                "c": c,
                "n": n,
                "true_switch_prime": mode_data.true_switch_prime,
                "nontrivial_balance_prime": mode_data.nontrivial_balance_prime,
                "nt_balance_L_abs": mode_data.nontrivial_balance_L_abs,
                "balance_is_smallest_prime": mode_data.nontrivial_balance_prime == balance_expected,
                "switch_is_largest_prime": mode_data.true_switch_prime == switch_expected_largest_prime,
                "switch_is_max_logmass": mode_data.true_switch_prime == switch_expected_logmass,
                "has_true_switch": mode_data.true_switch_prime is not None,
                "has_nontrivial_balance": mode_data.nontrivial_balance_prime is not None,
                "mass_p": mass_p,
                "mass_q": mass_q,
                "mass_r": mass_r,
                "argmax_logmass_prime": switch_expected_logmass,
                "initial_sign": mode_data.initial_sign,
                "final_sign": mode_data.final_sign,
            }
        )
    return rows


def summarize(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[str, str], List[Dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["family"]), str(row["mode"]))].append(row)

    out: List[Dict[str, object]] = []
    for (family, mode), block in sorted(grouped.items()):
        total = len(block)
        has_switch = sum(1 for r in block if r["has_true_switch"])
        has_balance = sum(1 for r in block if r["has_nontrivial_balance"])
        bal_small = sum(1 for r in block if r["balance_is_smallest_prime"])
        sw_largest = sum(1 for r in block if r["switch_is_largest_prime"])
        sw_logmass = sum(1 for r in block if r["switch_is_max_logmass"])
        out.append(
            {
                "family": family,
                "mode": mode,
                "total": total,
                "true_switch_rate": has_switch / total if total else 0.0,
                "nt_balance_rate": has_balance / total if total else 0.0,
                "balance_is_smallest_rate": bal_small / total if total else 0.0,
                "switch_is_largest_rate": sw_largest / total if total else 0.0,
                "switch_is_logmass_rate": sw_logmass / total if total else 0.0,
            }
        )
    return out


def print_summary(summary: Sequence[Dict[str, object]], pattern: str) -> None:
    print("=" * 150)
    print(f"DREIFAKTOR-TEST ({pattern})")
    print("=" * 150)
    header = (
        f"{'Familie':<16} {'Modus':<6} {'N':>7} {'echter W.':>10} {'nt. Bal.':>10} "
        f"{'Bal.=p':>9} {'Wech.=r':>10} {'Wech.=maxLog':>14}"
    )
    print(header)
    print("-" * len(header))
    for row in summary:
        print(
            f"{row['family']:<16} {row['mode']:<6} {row['total']:>7} "
            f"{row['true_switch_rate']*100:>9.1f}% "
            f"{row['nt_balance_rate']*100:>9.1f}% "
            f"{row['balance_is_smallest_rate']*100:>8.1f}% "
            f"{row['switch_is_largest_rate']*100:>9.1f}% "
            f"{row['switch_is_logmass_rate']*100:>13.1f}%"
        )


def print_sample(rows: Sequence[Dict[str, object]], limit: int = 30) -> None:
    interesting = [
        r for r in rows
        if r["has_true_switch"] and (
            not r["switch_is_largest_prime"] or not r["switch_is_max_logmass"] or not r["balance_is_smallest_prime"]
        )
    ]
    print("\n" + "=" * 150)
    print("AUSGEWÄHLTE INTERESSANTE FÄLLE")
    print("=" * 150)
    if not interesting:
        print("Keine abweichenden Fälle in der Auswahl gefunden.")
        return

    header = (
        f"{'Fam.':<16} {'n':>10} {'p^a':>8} {'q^b':>8} {'r^c':>8} {'Bal.':>6} {'Wech.':>6} "
        f"{'maxLog':>7} {'|L|_min':>12}"
    )
    print(header)
    print("-" * len(header))
    for row in interesting[:limit]:
        print(
            f"{row['family']:<16} {row['n']:>10} "
            f"{str(row['p']) + '^' + str(row['a']):>8} "
            f"{str(row['q']) + '^' + str(row['b']):>8} "
            f"{str(row['r']) + '^' + str(row['c']):>8} "
            f"{str(row['nontrivial_balance_prime']):>6} "
            f"{str(row['true_switch_prime']):>6} "
            f"{str(row['argmax_logmass_prime']):>7} "
            f"{str(round(float(row['nt_balance_L_abs']), 6)) if row['nt_balance_L_abs'] is not None else 'None':>12}"
        )
    if len(interesting) > limit:
        print(f"... weitere {len(interesting) - limit} nicht angezeigt.")


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
        description="Erste systematische Tests für den Dreifaktorfall p^a q^b r^c im aktiven Modus."
    )
    parser.add_argument("--prime-limit", type=int, default=40, help="obere Primgrenze")
    parser.add_argument(
        "--pattern",
        type=str,
        default="all_small",
        choices=["pqr", "p2qr", "pq2r", "pqr2", "all_small"],
        help="welche Exponentenmuster getestet werden sollen",
    )
    parser.add_argument("--max-k", type=int, default=25, help="Anzahl der Schalenprimzahlen")
    parser.add_argument("--show", type=int, default=30, help="wie viele interessante Fälle gezeigt werden sollen")
    parser.add_argument("--csv", type=Path, default=None, help="CSV-Datei für Rohdaten")
    parser.add_argument("--summary-csv", type=Path, default=None, help="CSV-Datei für Zusammenfassung")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    rows = evaluate_cases(args.prime_limit, args.pattern, args.max_k)
    summary = summarize(rows)
    print_summary(summary, pattern=args.pattern)
    print_sample(rows, limit=args.show)

    if args.csv is not None:
        write_csv(args.csv, rows)
        print(f"\nRohdaten gespeichert: {args.csv}")
    if args.summary_csv is not None:
        write_csv(args.summary_csv, summary)
        print(f"Zusammenfassung gespeichert: {args.summary_csv}")


if __name__ == "__main__":
    main()
