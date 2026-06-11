from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple


MOD12_CLASSES = {
    1: "E",
    5: "A",
    7: "B",
    11: "C",
}


@dataclass(frozen=True)
class ShellResult:
    k: int
    shell_primes: Tuple[int, ...]
    count: int
    mean_log_E: float
    mean_log_A: float
    mean_log_B: float
    mean_log_C: float
    weight_E: float
    weight_A: float
    weight_B: float
    weight_C: float
    diff_B_minus_C: float
    diff_C_minus_A: float
    diff_A_minus_E: float
    ordering: Tuple[str, ...]


@dataclass(frozen=True)
class SampleRecord:
    n: int
    shell_part: int
    rest: int
    E: int
    A: int
    B: int
    C: int


class SPFTable:
    """Smallest-prime-factor table for fast factorizations up to a fixed bound."""

    def __init__(self, limit: int) -> None:
        if limit < 2:
            raise ValueError("limit must be at least 2")
        self.limit = limit
        self.spf = self._build(limit)

    @staticmethod
    def _build(limit: int) -> List[int]:
        spf = list(range(limit + 1))
        spf[0] = 0
        spf[1] = 1
        for p in range(2, int(limit**0.5) + 1):
            if spf[p] == p:
                for multiple in range(p * p, limit + 1, p):
                    if spf[multiple] == multiple:
                        spf[multiple] = p
        return spf

    def factor(self, n: int) -> Dict[int, int]:
        if n < 1 or n > self.limit:
            raise ValueError(f"n must satisfy 1 <= n <= {self.limit}")
        factors: Dict[int, int] = {}
        while n > 1:
            p = self.spf[n]
            factors[p] = factors.get(p, 0) + 1
            n //= p
        return factors


def first_primes(count: int) -> List[int]:
    if count < 1:
        raise ValueError("count must be positive")
    primes: List[int] = []
    candidate = 2
    while len(primes) < count:
        is_prime = True
        root = int(candidate**0.5)
        for p in primes:
            if p > root:
                break
            if candidate % p == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(candidate)
        candidate += 1 if candidate == 2 else 2
    return primes


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
                f"prime {p} in rest has residue {residue} mod 12; "
                "this indicates that the shell did not remove all shell primes."
            )
    return E, A, B, C


def restricted_iter(limit: int, coprime_to_6_only: bool) -> Iterator[int]:
    if coprime_to_6_only:
        for n in range(1, limit + 1):
            if math.gcd(n, 6) == 1:
                yield n
    else:
        yield from range(1, limit + 1)


def _safe_log(x: int) -> float:
    return math.log(x) if x > 1 else 0.0


def _ordering(weights: Dict[str, float]) -> Tuple[str, ...]:
    return tuple(name for name, _ in sorted(weights.items(), key=lambda item: item[1], reverse=True))


def analyse_shell(
    limit: int,
    shell_primes: Sequence[int],
    spf: SPFTable,
    *,
    coprime_to_6_only: bool = True,
    collect_samples: int = 0,
) -> Tuple[ShellResult, List[SampleRecord]]:
    if not shell_primes:
        raise ValueError("shell_primes must not be empty")

    shell_set = set(shell_primes)
    sample_rows: List[SampleRecord] = []

    sum_log_E = 0.0
    sum_log_A = 0.0
    sum_log_B = 0.0
    sum_log_C = 0.0
    count = 0

    for n in restricted_iter(limit, coprime_to_6_only):
        factors = spf.factor(n)
        shell_part = shell_part_from_factors(factors, shell_primes)

        rest_factors = {p: exp for p, exp in factors.items() if p not in shell_set}
        rest = n // shell_part
        E, A, B, C = split_rest_classes(rest_factors)

        sum_log_E += _safe_log(E)
        sum_log_A += _safe_log(A)
        sum_log_B += _safe_log(B)
        sum_log_C += _safe_log(C)
        count += 1

        if len(sample_rows) < collect_samples:
            sample_rows.append(
                SampleRecord(
                    n=n,
                    shell_part=shell_part,
                    rest=rest,
                    E=E,
                    A=A,
                    B=B,
                    C=C,
                )
            )

    if count == 0:
        raise RuntimeError("no numbers were analysed; check the input range")

    mean_log_E = sum_log_E / count
    mean_log_A = sum_log_A / count
    mean_log_B = sum_log_B / count
    mean_log_C = sum_log_C / count
    total = mean_log_E + mean_log_A + mean_log_B + mean_log_C
    if total == 0.0:
        raise RuntimeError("all logarithmic means vanished; this should not happen for a nontrivial range")

    weight_E = mean_log_E / total
    weight_A = mean_log_A / total
    weight_B = mean_log_B / total
    weight_C = mean_log_C / total

    weights = {"E": weight_E, "A": weight_A, "B": weight_B, "C": weight_C}

    result = ShellResult(
        k=len(shell_primes),
        shell_primes=tuple(shell_primes),
        count=count,
        mean_log_E=mean_log_E,
        mean_log_A=mean_log_A,
        mean_log_B=mean_log_B,
        mean_log_C=mean_log_C,
        weight_E=weight_E,
        weight_A=weight_A,
        weight_B=weight_B,
        weight_C=weight_C,
        diff_B_minus_C=mean_log_B - mean_log_C,
        diff_C_minus_A=mean_log_C - mean_log_A,
        diff_A_minus_E=mean_log_A - mean_log_E,
        ordering=_ordering(weights),
    )
    return result, sample_rows


def analyse_shell_family(
    limit: int,
    max_k: int,
    *,
    coprime_to_6_only: bool = True,
    collect_samples: int = 0,
) -> Tuple[List[ShellResult], Dict[int, List[SampleRecord]]]:
    if max_k < 1:
        raise ValueError("max_k must be at least 1")
    spf = SPFTable(limit)
    primes = first_primes(max_k)

    results: List[ShellResult] = []
    samples: Dict[int, List[SampleRecord]] = {}
    for k in range(1, max_k + 1):
        shell = primes[:k]
        result, rows = analyse_shell(
            limit,
            shell,
            spf,
            coprime_to_6_only=coprime_to_6_only,
            collect_samples=collect_samples,
        )
        results.append(result)
        if rows:
            samples[k] = rows
    return results, samples


def write_results_csv(path: Path, results: Sequence[ShellResult]) -> None:
    fieldnames = [
        "k",
        "shell_primes",
        "count",
        "mean_log_E",
        "mean_log_A",
        "mean_log_B",
        "mean_log_C",
        "weight_E",
        "weight_A",
        "weight_B",
        "weight_C",
        "diff_B_minus_C",
        "diff_C_minus_A",
        "diff_A_minus_E",
        "ordering",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            row = asdict(result)
            row["shell_primes"] = ",".join(str(p) for p in result.shell_primes)
            row["ordering"] = ">".join(result.ordering)
            writer.writerow(row)


def write_samples_json(path: Path, samples: Dict[int, List[SampleRecord]]) -> None:
    payload = {
        str(k): [asdict(row) for row in rows]
        for k, rows in samples.items()
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def render_summary(results: Sequence[ShellResult]) -> str:
    lines = []
    for result in results:
        shell_str = "{" + ",".join(str(p) for p in result.shell_primes) + "}"
        lines.append(
            f"k={result.k:>2} shell={shell_str:<30} "
            f"count={result.count:<10} "
            f"W=(E={result.weight_E:.6f}, A={result.weight_A:.6f}, "
            f"B={result.weight_B:.6f}, C={result.weight_C:.6f}) "
            f"order={' > '.join(result.ordering)}"
        )
    return "\n".join(lines)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Testet die numerische Hypothese zu normierten logarithmischen Restklassenbeiträgen "
            "nach Schalenreduktion."
        )
    )
    parser.add_argument("--limit", type=int, default=10**6, help="obere Grenze N")
    parser.add_argument("--max-k", type=int, default=10, help="Anzahl der ersten Primzahlen in der größten Schale")
    parser.add_argument(
        "--all-numbers",
        action="store_true",
        help="benutzt alle Zahlen 1..N statt nur der Zahlen mit gcd(n,6)=1",
    )
    parser.add_argument(
        "--collect-samples",
        type=int,
        default=0,
        help="speichert pro Schale die ersten Beispielzeilen",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Pfad für CSV-Export der Hauptergebnisse",
    )
    parser.add_argument(
        "--samples-json",
        type=Path,
        default=None,
        help="Pfad für JSON-Export der Beispielzeilen",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    results, samples = analyse_shell_family(
        limit=args.limit,
        max_k=args.max_k,
        coprime_to_6_only=not args.all_numbers,
        collect_samples=args.collect_samples,
    )

    print(render_summary(results))

    if args.csv is not None:
        write_results_csv(args.csv, results)
        print(f"\nCSV gespeichert: {args.csv}")

    if args.samples_json is not None:
        write_samples_json(args.samples_json, samples)
        print(f"Beispieldaten gespeichert: {args.samples_json}")


if __name__ == "__main__":
    main()
