#!/usr/bin/env python3
"""
Klassische Identität ζ(2n) vs. Bernoulli-Zahlen
=================================================

Verifiziert für n = 1 … N_MAX:

    ζ(2n) = (-1)^(n+1) * B_{2n} * (2π)^{2n} / (2 * (2n)!)

Referenz: collatz_eabc_zeta_exponential_gedankenexperiment.md §6.2, §8.1
Kein EABC-Anspruch — reine Zahlentheorie-Diagnose.
"""

from __future__ import annotations

import sys
from math import factorial, pi

from mpmath import bernoulli, mp, zeta

N_MAX = 10
TOL = mp.mpf("1e-12")


def zeta_even_from_bernoulli(n: int) -> mp.mpf:
    """n ≥ 1: closed form for ζ(2n)."""
    b = bernoulli(2 * n)
    return (-1) ** (n + 1) * b * (2 * pi) ** (2 * n) / (2 * factorial(2 * n))


def main() -> int:
    mp.dps = 50
    ok = True
    print(f"{'n':>3}  {'ζ(2n)':>22}  {'Bernoulli':>22}  {'|Δ|':>12}")
    print("-" * 65)
    for n in range(1, N_MAX + 1):
        s = 2 * n
        direct = zeta(s)
        via_b = zeta_even_from_bernoulli(n)
        delta = abs(direct - via_b)
        flag = "ok" if delta < TOL else "FAIL"
        if delta >= TOL:
            ok = False
        print(
            f"{n:3d}  {float(direct):22.12f}  {float(via_b):22.12f}  "
            f"{float(delta):12.2e}  {flag}"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
