"""Command-line interface for basel-with-a."""

from __future__ import annotations

import argparse

from .core import analyze_first_primes


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("count muss mindestens 1 sein")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Berechnet die Quadratsumme der ersten n Primzahlen ohne SageMath."
    )
    parser.add_argument(
        "--count",
        type=positive_int,
        default=100,
        help="Anzahl der Primzahlen, standardmaessig 100.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    summary = analyze_first_primes(args.count)

    print(f"Anzahl Primzahlen: {summary.prime_count}")
    print(f"{summary.prime_count}. Primzahl: {summary.nth_prime}")
    print(f"Summe p^2: {summary.prime_square_sum}")
    print(f"zeta(2): {summary.zeta_2:.12f}")
    print(f"Differenz: {summary.difference:.12f}")
    print(f"Indizwert a: {summary.index_value_a:.12f}")
