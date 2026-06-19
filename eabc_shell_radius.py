#!/usr/bin/env python3
"""
EABC-Schalenradius — numerischer Stub (§4.8.20–§4.8.21 in collatz_eabc_zirkulationshypothese.md).

Primverteilung auf Schalen: eabc_shell_prime_distribution.py (§4.8.21).

Berechnet pythagoräische Schalenradien
    R_{n,m} = sqrt(2^n + 3^m)
für kleine Exponenten n, m.

Verknüpfungen: glatte (2,3)-Schale, Radikalschichten, EABC mod 12.
"""

from __future__ import annotations

import argparse
import math


def shell_radius(n: int, m: int) -> float:
    """Schalenradius R_{n,m} = sqrt(2^n + 3^m)."""
    if n < 0 or m < 0:
        raise ValueError("n und m müssen nichtnegativ sein")
    return math.sqrt(2**n + 3**m)


def shell_radius_squared(n: int, m: int) -> int:
    """Pythagoras-Untergrund R_{n,m}^2 = 2^n + 3^m."""
    if n < 0 or m < 0:
        raise ValueError("n und m müssen nichtnegativ sein")
    return 2**n + 3**m


def main() -> None:
    parser = argparse.ArgumentParser(
        description="EABC-Schalenradien R_{n,m} = sqrt(2^n + 3^m)"
    )
    parser.add_argument(
        "--n-max",
        type=int,
        default=6,
        help="maximaler 2-Exponent (default: 6)",
    )
    parser.add_argument(
        "--m-max",
        type=int,
        default=6,
        help="maximaler 3-Exponent (default: 6)",
    )
    args = parser.parse_args()

    print("n\tm\t2^n\t3^m\tR_{n,m}^2\tR_{n,m}")
    for n in range(args.n_max + 1):
        for m in range(args.m_max + 1):
            r2 = shell_radius_squared(n, m)
            r = shell_radius(n, m)
            print(f"{n}\t{m}\t{2**n}\t{3**m}\t{r2}\t{r:.6f}")


if __name__ == "__main__":
    main()
