#!/usr/bin/env python3
"""
Fibonacci: Binet vs. Matrixpotenz
=================================

Verifiziert für n = 0 … N_MAX:

    F_n = (φ^n - ψ^n) / √5,   ψ = -φ^{-1}

gegen

    (M^n)_{1,0}  mit  M = [[1,1],[1,0]]

Referenz: collatz_eabc_zeta_exponential_gedankenexperiment.md §9.2
Kein EABC-Anspruch — reine Zahlentheorie-Diagnose.
"""

from __future__ import annotations

import sys

N_MAX = 20

PHI = (1 + 5**0.5) / 2
PSI = -1 / PHI


def binet(n: int) -> int:
    return round((PHI**n - PSI**n) / 5**0.5)


def fib_matrix(n: int) -> int:
    """(M^n)_{0,1} = F_n for M = [[1,1],[1,0]], F_0 = 0."""
    if n == 0:
        return 0.0

    def mul(a: list, b: list) -> list:
        return [
            [a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]],
            [a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]],
        ]

    m = [[1, 1], [1, 0]]
    r = [[1, 0], [0, 1]]
    exp = n
    while exp:
        if exp & 1:
            r = mul(r, m)
        m = mul(m, m)
        exp >>= 1
    return r[0][1]


def main() -> int:
    ok = True
    print(f"{'n':>3}  {'Binet':>18}  {'Matrix':>18}  {'|Δ|':>12}")
    print("-" * 58)
    for n in range(N_MAX + 1):
        via_binet = binet(n)
        via_mat = fib_matrix(n)
        delta = abs(via_binet - via_mat)
        flag = "ok" if delta == 0 else "FAIL"
        if delta != 0:
            ok = False
        print(f"{n:3d}  {via_binet:18.6f}  {via_mat:18.6f}  {delta:12.2e}  {flag}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
