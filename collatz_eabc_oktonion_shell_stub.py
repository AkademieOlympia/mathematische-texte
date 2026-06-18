#!/usr/bin/env python3
"""
Oktanionische Normschalen — Minimal-Stub (nur r_8(n), kein μ_n).

Kanonsiche Theorie: collatz_eabc_oktonion_singularitaet.md

Zählt Darstellungen n = x_1^2 + ... + x_8^2 mit x_i ∈ Z (Λ_O^(8) = Z^8 Stub).
Kein Hurwitz-Oktanionen-Gitter, keine glatt-EABC-Γ-Signatur, kein Schalenmaß μ_n.

Theorem-Referenz: r_8(n) = Koeffizient von q^n in θ_3(q)^8 (Jacobi / Acht-Quadrate).

Ausführung:
    python3 collatz_eabc_oktonion_shell_stub.py
    python3 collatz_eabc_oktonion_shell_stub.py --max-n 10
"""

from __future__ import annotations

import argparse
import json
from math import isqrt
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_oktonion_shell_stub.json"


def r8_brute(n: int) -> int:
    """Zähle (x_1,...,x_8) ∈ Z^8 mit Σ x_i^2 = n."""
    if n < 0:
        return 0
    if n == 0:
        return 1
    count = 0
    bound = isqrt(n)
    for x1 in range(-bound, bound + 1):
        s1 = x1 * x1
        if s1 > n:
            continue
        for x2 in range(-bound, bound + 1):
            s2 = s1 + x2 * x2
            if s2 > n:
                continue
            for x3 in range(-bound, bound + 1):
                s3 = s2 + x3 * x3
                if s3 > n:
                    continue
                for x4 in range(-bound, bound + 1):
                    s4 = s3 + x4 * x4
                    if s4 > n:
                        continue
                    for x5 in range(-bound, bound + 1):
                        s5 = s4 + x5 * x5
                        if s5 > n:
                            continue
                        for x6 in range(-bound, bound + 1):
                            s6 = s5 + x6 * x6
                            if s6 > n:
                                continue
                            for x7 in range(-bound, bound + 1):
                                s7 = s6 + x7 * x7
                                if s7 > n:
                                    continue
                                rem = n - s7
                                d = isqrt(rem)
                                if d * d != rem:
                                    continue
                                if d == 0:
                                    count += 1
                                else:
                                    count += 2
    return count


def shell_table(max_n: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for n in range(1, max_n + 1):
        brute = r8_brute(n)
        rows.append(
            {
                "n": n,
                "r_8(n)": brute,
                "is_prime": _is_prime(n),
                "omega": _omega(n),
                "growth_note": f"r_8({n})={brute}; OEIS A000118; asymptotisch ~ C·n^(3/2)",
            }
        )
    return rows


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def _omega(n: int) -> int:
    if n < 2:
        return 0
    count = 0
    d = 2
    m = n
    while d * d <= m:
        if m % d == 0:
            count += 1
            while m % d == 0:
                m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        count += 1
    return count


def main() -> None:
    parser = argparse.ArgumentParser(description="Oktanionische Schalen-Stub: r_8(n) nur")
    parser.add_argument("--max-n", type=int, default=10, help="Oberes n (Default 10)")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    rows = shell_table(args.max_n)

    payload = {
        "meta": {
            "hypothesis_doc": "collatz_eabc_oktonion_singularitaet.md",
            "lattice": "Z^8 stub (nicht Hurwitz O_H)",
            "computed": ["r_8(n) shell size |Sigma_n^(8)|"],
            "not_computed": [
                "mu_n (glatt-EABC auf 8 Beinen)",
                "H_n, chi_n, K_n, D(n)",
                "U_O orbit partition (|U_O|=240)",
                "Hopf S^7 -> S^4 Faserstatistik",
            ],
            "max_n": args.max_n,
            "formula": "Koeffizient von q^n in θ_3(q)^8 (Jacobi; vgl. OEIS A000118)",
            "theorem_ref": "Jacobi theta_3^8 coefficient / eight squares",
            "quaternion_crossref": "collatz_eabc_shell_defekt.json (D(n) prime anomaly NOT at n<=200, rolling)",
        },
        "rows": rows,
        "verification": {
            "method": "brute_enumeration_Z8",
            "oeis_crosscheck": "A000118 first terms: 1,16,112,448,1136,...",
            "epistemic_note": "Nur Schalengröße; keine EABC-Maßtheorie in 8D",
        },
    }

    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    print(f"Wrote {args.output} ({len(rows)} shells)")


if __name__ == "__main__":
    main()
