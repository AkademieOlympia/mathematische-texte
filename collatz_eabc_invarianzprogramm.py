#!/usr/bin/env python3
"""
EABC-Invarianzprogramm: strikte Zählvektoren, Simplex-Anteile und Quadrupel-Signaturen.

Kanonsiche Definitionen: collatz_eabc_invarianzprogramm.md
Philosophischer Kontext: collatz_eabc_bernoulli_uebersetzung.md §22 (Querverweis).

Ausführung:
    python3 collatz_eabc_invarianzprogramm.py
    python3 collatz_eabc_invarianzprogramm.py --max-x 50000 --output collatz_eabc_invarianzprogramm.json
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from math import isqrt
from pathlib import Path
from typing import Any

from eabc_from_lean import EClass, class_of, is_prime_quadruplet, q

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_invarianzprogramm.json"

EABC_ORDER = (EClass.E, EClass.A, EClass.B, EClass.C)
SIGMA4_ALPHABET = "EABC"


def _sieve_primes(limit: int) -> list[int]:
    if limit < 2:
        return []
    flags = [True] * (limit + 1)
    flags[0] = flags[1] = False
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            step = p
            start = p * p
            flags[start : limit + 1 : step] = [False] * (((limit - start) // step) + 1)
    return [i for i, ok in enumerate(flags) if ok]


def kappa(p: int) -> EClass | None:
    """κ: P_{>3} → {E,A,B,C}; für p ≤ 3 oder nicht-prim: None."""
    if p <= 3:
        return None
    return class_of(p)


@dataclass(frozen=True, slots=True)
class CountVector:
    """V(x) = (E(x), A(x), B(x), C(x)) — Zählvektor über p ≤ x in P_{>3}."""

    e: int
    a: int
    b: int
    c: int

    def as_tuple(self) -> tuple[int, int, int, int]:
        return (self.e, self.a, self.b, self.c)

    @property
    def total(self) -> int:
        return self.e + self.a + self.b + self.c

    def to_dict(self) -> dict[str, int]:
        return {"E": self.e, "A": self.a, "B": self.b, "C": self.c}


@dataclass(frozen=True, slots=True)
class SimplexPoint:
    """S(x) = (1/π(x))·V(x) ∈ Δ_3; π(x) hier = |{p ≤ x : p ∈ P_{>3}}|."""

    e: float
    a: float
    b: float
    c: float

    def as_tuple(self) -> tuple[float, float, float, float]:
        return (self.e, self.a, self.b, self.c)

    def to_dict(self) -> dict[str, float]:
        return {"E": self.e, "A": self.a, "B": self.b, "C": self.c}


def v_at_x(x: int, primes: list[int] | None = None) -> CountVector:
    """V(x): zählt Primzahlen p ≤ x, p > 3, nach EABC-Klasse."""
    if x < 5:
        return CountVector(0, 0, 0, 0)
    if primes is None:
        primes = _sieve_primes(x)
    counts = {cls: 0 for cls in EABC_ORDER}
    for p in primes:
        if p > x:
            break
        cls = kappa(p)
        if cls is not None:
            counts[cls] += 1
    return CountVector(
        e=counts[EClass.E],
        a=counts[EClass.A],
        b=counts[EClass.B],
        c=counts[EClass.C],
    )


def pi_eabc(x: int, v: CountVector | None = None) -> int:
    """π_{>3}(x) = E(x)+A(x)+B(x)+C(x); Nenner für S(x) und χ(x)."""
    if v is None:
        v = v_at_x(x)
    return v.total


def s_at_x(x: int, v: CountVector | None = None) -> SimplexPoint:
    """S(x) = (1/π(x))·V(x) auf dem Standard-3-Simplex (Summe 1)."""
    if v is None:
        v = v_at_x(x)
    denom = v.total
    if denom == 0:
        return SimplexPoint(0.0, 0.0, 0.0, 0.0)
    return SimplexPoint(
        e=v.e / denom,
        a=v.a / denom,
        b=v.b / denom,
        c=v.c / denom,
    )


def chi_at_x(x: int, v: CountVector | None = None) -> float:
    """χ(x) = ((E+C)−(A+B)) / π(x) — Beispiel-EABC-Invariante (Def. 4)."""
    if v is None:
        v = v_at_x(x)
    denom = v.total
    if denom == 0:
        return 0.0
    return ((v.e + v.c) - (v.a + v.b)) / denom


def sigma_quadruplet(p: int) -> str | None:
    """σ(Q) = (κ(p), κ(p+2), κ(p+6), κ(p+8)) als String über Σ_4."""
    if not is_prime_quadruplet(p):
        return None
    parts: list[str] = []
    for val in q(p):
        cls = kappa(val)
        if cls is None:
            return None
        parts.append(cls.value)
    return "".join(parts)


def quadruplet_signature_frequencies(
    max_p: int, primes: list[int] | None = None
) -> dict[str, Any]:
    """Schätzt μ auf Σ_4 durch Häufigkeiten der σ(Q) für Prim-Vierlinge mit p ≤ max_p."""
    if primes is None:
        primes = _sieve_primes(max_p + 8)
    counts: dict[str, int] = {}
    total = 0
    for p in primes:
        if p > max_p:
            break
        sig = sigma_quadruplet(p)
        if sig is None:
            continue
        counts[sig] = counts.get(sig, 0) + 1
        total += 1
    frequencies = {sig: c / total for sig, c in sorted(counts.items())} if total else {}
    return {
        "max_p": max_p,
        "quadruplet_count": total,
        "signature_counts": dict(sorted(counts.items())),
        "signature_frequencies": frequencies,
        "sigma4_alphabet": SIGMA4_ALPHABET,
        "sigma4_cardinality": len(SIGMA4_ALPHABET) ** 4,
    }


@dataclass(frozen=True, slots=True)
class Snapshot:
    x: int
    pi: int
    v: CountVector
    s: SimplexPoint
    chi: float


def snapshot_at_x(x: int, primes: list[int] | None = None) -> Snapshot:
    v = v_at_x(x, primes)
    return Snapshot(
        x=x,
        pi=v.total,
        v=v,
        s=s_at_x(x, v),
        chi=chi_at_x(x, v),
    )


def run_program(
    max_x: int,
    sample_points: list[int] | None = None,
) -> dict[str, Any]:
    primes = _sieve_primes(max_x + 8)
    if sample_points is None:
        sample_points = sorted(
            {
                20,
                100,
                500,
                1000,
                5000,
                max_x,
            }
            & {x for x in range(5, max_x + 1)}
        )
        if not sample_points:
            sample_points = [max_x] if max_x >= 5 else []

    samples = []
    chi_series: list[dict[str, float | int]] = []
    for x in range(5, max_x + 1):
        v = v_at_x(x, primes)
        chi_val = chi_at_x(x, v)
        chi_series.append({"x": x, "chi": chi_val, "pi": v.total})
        if x in sample_points:
            snap = snapshot_at_x(x, primes)
            samples.append(
                {
                    "x": snap.x,
                    "pi_eabc": snap.pi,
                    "V": snap.v.to_dict(),
                    "S": snap.s.to_dict(),
                    "chi": snap.chi,
                }
            )

    quad_stats = quadruplet_signature_frequencies(max_x, primes)

    chi_values = [row["chi"] for row in chi_series]
    chi_trend = {
        "min": min(chi_values) if chi_values else 0.0,
        "max": max(chi_values) if chi_values else 0.0,
        "mean": sum(chi_values) / len(chi_values) if chi_values else 0.0,
        "last": chi_values[-1] if chi_values else 0.0,
        "first": chi_values[0] if chi_values else 0.0,
    }

    return {
        "program": "EABC-Invarianzprogramm",
        "canonical_doc": "collatz_eabc_invarianzprogramm.md",
        "max_x": max_x,
        "definitions": {
            "kappa": "P_{>3} → {E,A,B,C}, E≡1, A≡5, B≡7, C≡11 (mod 12)",
            "V": "Zählvektor (E,A,B,C) für p ≤ x, p > 3",
            "S": "V / π_{>3}(x) ∈ Δ_3",
            "chi": "((E+C)-(A+B)) / π_{>3}(x)",
            "sigma_quadruplet": "σ(Q)=(κ(p),κ(p+2),κ(p+6),κ(p+8)) für Q=(p,p+2,p+6,p+8)",
        },
        "samples": samples,
        "chi_trend": chi_trend,
        "chi_series_tail": chi_series[-20:] if len(chi_series) > 20 else chi_series,
        "quadruplet_signatures": quad_stats,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC-Invarianzprogramm numerisch")
    parser.add_argument("--max-x", type=int, default=10_000, help="Obergrenze für x")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="JSON-Ausgabedatei",
    )
    args = parser.parse_args()
    result = run_program(args.max_x)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(f"Wrote {args.output} (max_x={args.max_x})")
    for row in result["samples"]:
        print(
            f"x={row['x']:>6}  S={tuple(round(row['S'][k], 4) for k in 'EABC')}  "
            f"χ={row['chi']:+.4f}"
        )


if __name__ == "__main__":
    main()
