#!/usr/bin/env python3
"""
Kommutator-Norm-Stub auf Hurwitz-Normschalen Σ_n in H.

Kanonsiche Theorie: collatz_eabc_kommutator_assoziator.md

Zeigt für n ≤ max_n: generische Paare (x,y) ∈ Σ_n × Σ_n haben
  α_com(x,y) = N([x,y]) = N(xy - yx) > 0
(Quaternionen sind assoziativ → Assoziator ≡ 0, nicht getestet).

Ausführung:
    python3 collatz_eabc_quaternion_commutator_stub.py
    python3 collatz_eabc_quaternion_commutator_stub.py --max-n 20
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

from collatz_eabc_hurwitz_orbit_test import (
    QuatFrac,
    coords_to_frac,
    hurwitz_norm,
    hurwitz_shell_elements,
    q_mul,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_quaternion_commutator_stub.json"


def q_sub(x: QuatFrac, y: QuatFrac) -> QuatFrac:
    return tuple(a - b for a, b in zip(x, y))


def q_norm_sq(x: QuatFrac) -> Fraction:
    return sum(c * c for c in x)


def commutator(x: QuatFrac, y: QuatFrac) -> QuatFrac:
    """[x,y] = xy - yx."""
    return q_sub(q_mul(x, y), q_mul(y, x))


def commutator_norm_int(x: QuatFrac, y: QuatFrac) -> int:
    """N([x,y]) als ganzzahlige Norm (Hurwitz-Koordinaten)."""
    c = commutator(x, y)
    num = sum(int(co.numerator) ** 2 for co in c)
    den = c[0].denominator
    assert all(co.denominator == den for co in c)
    return num // (den * den)


def shell_commutator_report(n: int, max_pairs: int = 8) -> dict[str, Any]:
    elems = hurwitz_shell_elements(n)
    if len(elems) < 2:
        return {"n": n, "shell_size": len(elems), "skipped": True}
    fracs = [coords_to_frac(*e) for e in elems]
    nonzero = 0
    total = 0
    samples: list[dict[str, Any]] = []
    limit = min(len(fracs), max_pairs)
    for i in range(limit):
        for j in range(i + 1, limit):
            x, y = fracs[i], fracs[j]
            nc = commutator_norm_int(x, y)
            total += 1
            if nc > 0:
                nonzero += 1
            if len(samples) < 4 and nc > 0:
                samples.append(
                    {
                        "x": elems[i],
                        "y": elems[j],
                        "commutator_norm": nc,
                    }
                )
    frac = nonzero / total if total else 0.0
    return {
        "n": n,
        "shell_size": len(elems),
        "pairs_tested": total,
        "pairs_nonzero_commutator": nonzero,
        "fraction_nonzero": frac,
        "exists_nonzero_commutator": nonzero > 0,
        "samples": samples,
    }


def run(max_n: int = 20, output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    per_n = [shell_commutator_report(n) for n in range(2, max_n + 1)]
    tested = [r for r in per_n if not r.get("skipped")]
    all_have_nonzero = all(r["exists_nonzero_commutator"] for r in tested)
    min_frac = min(r["fraction_nonzero"] for r in tested) if tested else 0.0
    report: dict[str, Any] = {
        "meta": {
            "max_n": max_n,
            "theory": "collatz_eabc_kommutator_assoziator.md",
            "note": "H is associative; associator not tested (identically zero).",
        },
        "summary": {
            "levels_tested": len(tested),
            "all_shells_have_nonzero_pair": all_have_nonzero,
            "min_fraction_nonzero_pairs": min_frac,
        },
        "per_n": per_n,
        "output_path": str(output),
    }
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> None:
    p = argparse.ArgumentParser(description="Quaternion commutator norm stub on Σ_n")
    p.add_argument("--max-n", type=int, default=20)
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = p.parse_args()
    report = run(max_n=args.max_n, output=args.output)
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
