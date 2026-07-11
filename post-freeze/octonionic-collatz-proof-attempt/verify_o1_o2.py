#!/usr/bin/env python3
"""Reproducible checks for post-freeze Collatz modules O1 and O2.

O1: exact odd-core log cocycle (multiplicative form, exact rationals where possible).
O2: valuation-one ladder at n0 = 2^(L+1) - 1.
"""

from __future__ import annotations

import argparse
import math
import random
import sys
from fractions import Fraction


def nu2(n: int) -> int:
    if n <= 0:
        raise ValueError(f"nu2 requires positive n, got {n}")
    v = 0
    while n % 2 == 0:
        v += 1
        n //= 2
    return v


def odd_core_step(n: int) -> int:
    if n % 2 == 0:
        raise ValueError("odd_core_step expects odd n")
    y = 3 * n + 1
    return y >> nu2(y)


def orbit(n0: int, k: int) -> list[int]:
    path = [n0]
    n = n0
    for _ in range(k):
        n = odd_core_step(n)
        path.append(n)
    return path


def check_o1_multiplicative(n0: int, k: int) -> None:
    path = orbit(n0, k)
    lhs = Fraction(path[k], path[0])
    rhs = Fraction(1, 1)
    for j in range(k):
        nj = path[j]
        v = nu2(3 * nj + 1)
        rhs *= Fraction(3, 2**v) * Fraction(3 * nj + 1, 3 * nj)
    if lhs != rhs:
        raise AssertionError(f"O1 failed: n0={n0}, k={k}, lhs={lhs}, rhs={rhs}")


def check_o1_log(n0: int, k: int, tol: float = 1e-12) -> None:
    path = orbit(n0, k)
    log_lhs = math.log(path[k] / path[0])
    log_rhs = k * math.log(3.0)
    for j in range(k):
        nj = path[j]
        log_rhs -= math.log(2.0) * nu2(3 * nj + 1)
        log_rhs += math.log1p(1.0 / (3 * nj))
    if abs(log_lhs - log_rhs) > tol:
        raise AssertionError(
            f"O1 log form failed: n0={n0}, k={k}, delta={log_lhs - log_rhs}"
        )


def check_o2(L: int) -> None:
    n0 = (1 << (L + 1)) - 1
    n = n0
    for j in range(L + 1):
        expected = (1 << (L + 1 - j)) * (3**j) - 1
        if n != expected:
            raise AssertionError(f"O2 trajectory: L={L}, j={j}, n={n}, expected={expected}")
        if j < L and nu2(3 * n + 1) != 1:
            raise AssertionError(f"O2 valuation: L={L}, j={j}, nu2={nu2(3 * n + 1)}")
        if j < L:
            n = odd_core_step(n)


def run_o1(*, samples: int, k_max: int, seed: int) -> None:
    rng = random.Random(seed)
    # fixed small cases
    for n0 in [1, 3, 5, 7, 27, 255, 6171]:
        for k in range(1, min(k_max, 12) + 1):
            check_o1_multiplicative(n0, k)
            check_o1_log(n0, k)
    # random odd starts
    for _ in range(samples):
        n0 = rng.randrange(1, 10**6) | 1
        k = rng.randint(1, k_max)
        check_o1_multiplicative(n0, k)
        check_o1_log(n0, k)
    print(f"O1 OK: {samples} random + fixed cases, k_max={k_max}")


def run_o2(*, L_max: int) -> None:
    for L in range(1, L_max + 1):
        check_o2(L)
    print(f"O2 OK: L=1..{L_max} (valuation-one ladder + closed form)")


def main() -> int:
    p = argparse.ArgumentParser(description="Verify post-freeze Collatz modules O1/O2")
    p.add_argument(
        "--module",
        choices=("O1", "O2", "all"),
        default="all",
        help="Which module to run",
    )
    p.add_argument("--samples", type=int, default=200, help="Random O1 samples")
    p.add_argument("--k-max", type=int, default=24, help="Max orbit length for O1")
    p.add_argument("--L-max", type=int, default=32, help="Max L for O2 ladder")
    p.add_argument("--seed", type=int, default=20260711)
    args = p.parse_args()

    if args.module in ("O1", "all"):
        run_o1(samples=args.samples, k_max=args.k_max, seed=args.seed)
    if args.module in ("O2", "all"):
        run_o2(L_max=args.L_max)
    return 0


if __name__ == "__main__":
    sys.exit(main())
