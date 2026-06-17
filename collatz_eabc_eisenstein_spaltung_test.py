#!/usr/bin/env python3
"""
Eisenstein–EABC-Spaltungstest mit glatt-EABC-Zerlegung.

Kanonsiche Hypothese: collatz_eabc_eisenstein_spaltung.md

Für split-Primzahlen p ≡ 1 (mod 3), p ≠ 3, p = a² − ab + b² (kanonisch 0 < a ≤ b):
  a = 2^α_a · 3^β_a · a',  gcd(a', 6) = 1
  b = 2^α_b · 3^β_b · b',  gcd(b', 6) = 1
  Γ_E(p) = (κ(a'), κ(b')) ∈ {E,A,B,C}²

μ_X(γ) = #{p ≤ X split : Γ_E(p)=γ} / #{p ≤ X split}

Ausführung:
    python3 collatz_eabc_eisenstein_spaltung_test.py
    python3 collatz_eabc_eisenstein_spaltung_test.py --max-p 100000 1000000
"""

from __future__ import annotations

import argparse
import json
import math
import random
import statistics
from collections import Counter
from dataclasses import dataclass
from math import isqrt
from pathlib import Path
from typing import Any

from collatz_eabc_gauss_spaltung_test import (
    CHI2_CRIT_15DF_005,
    EABC_LABELS,
    GAMMA_PAIRS,
    NULL_TRIALS,
    UNIFORM_MU,
    _chi_square_uniform,
    _mu_table,
    _shuffle_null_chi2,
    _sieve_primes,
    _top_deviations,
    _verdict,
    gamma_label,
    kappa_glatt,
    strip_smooth,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_eisenstein_spaltung.json"

EB_CLASSES = frozenset({"E", "B"})
AC_CLASSES = frozenset({"A", "C"})


def eisenstein_norm(a: int, b: int) -> int:
    """N(a + bω) = a² − ab + b² für ω = e^(2πi/3)."""
    return a * a - a * b + b * b


def eisenstein_factor_pair(p: int) -> tuple[int, int] | None:
    """
    Kanonisch 0 < a ≤ b mit p = a² − ab + b²; nur für p ≡ 1 (mod 3), p ≠ 3.

    Nutzt 4p = (2a − b)² + 3b² → O(√p) pro Primzahl.
    """
    if p <= 0 or p % 3 != 1 or p == 3:
        return None
    best: tuple[int, int] | None = None
    limit = isqrt((4 * p) // 3) + 1
    for b in range(1, limit + 1):
        t = 4 * p - 3 * b * b
        if t < 0:
            continue
        u = isqrt(t)
        if u * u != t:
            continue
        for uu in (u, -u):
            if (uu + b) % 2 != 0:
                continue
            a = (uu + b) // 2
            if a <= 0 or eisenstein_norm(a, b) != p:
                continue
            pair = (a, b) if a <= b else (b, a)
            if best is None or pair < best:
                best = pair
    return best


def eisenstein_split_class(p: int) -> str | None:
    """split (p≡1 mod 3), inert (p≡2 mod 3), ramified (p=3)."""
    if p == 3:
        return "ramified"
    if p < 2:
        return None
    r = p % 3
    if r == 1:
        return "split"
    if r == 2:
        return "inert"
    return None


def eabc_mod3_bucket(label: str) -> str:
    """Grobe mod-3-Seite: EB (≡1 mod 3) vs AC (≡2 mod 3)."""
    return "EB" if label in EB_CLASSES else "AC"


@dataclass(frozen=True, slots=True)
class EisensteinSplitRow:
    p: int
    a: int
    b: int
    alpha_a: int
    beta_a: int
    a_prime: int
    kappa_a: str
    alpha_b: int
    beta_b: int
    b_prime: int
    kappa_b: str
    gamma: tuple[str, str]

    def to_dict(self) -> dict[str, Any]:
        return {
            "p": self.p,
            "a": self.a,
            "b": self.b,
            "alpha_a": self.alpha_a,
            "beta_a": self.beta_a,
            "a_prime": self.a_prime,
            "kappa_a": self.kappa_a,
            "alpha_b": self.alpha_b,
            "beta_b": self.beta_b,
            "b_prime": self.b_prime,
            "kappa_b": self.kappa_b,
            "gamma": list(self.gamma),
        }


def classify_eisenstein_split_prime(p: int) -> EisensteinSplitRow | None:
    if p <= 3 or p % 3 != 1:
        return None
    pair = eisenstein_factor_pair(p)
    if pair is None:
        return None
    a, b = pair
    alpha_a, beta_a, a_prime, ka = kappa_glatt(a)
    alpha_b, beta_b, b_prime, kb = kappa_glatt(b)
    return EisensteinSplitRow(
        p=p,
        a=a,
        b=b,
        alpha_a=alpha_a,
        beta_a=beta_a,
        a_prime=a_prime,
        kappa_a=ka,
        alpha_b=alpha_b,
        beta_b=beta_b,
        b_prime=b_prime,
        kappa_b=kb,
        gamma=(ka, kb),
    )


def coarse_mod3_defekt_report(max_p: int) -> dict[str, Any]:
    """Analog zu collatz_eabc_gauss_defekt_test: split/inert vs EB/AC."""
    rows: list[dict[str, Any]] = []
    mismatches = 0
    total = 0
    for p in _sieve_primes(max_p):
        if p <= 3:
            continue
        split_cls = eisenstein_split_class(p)
        if split_cls is None or split_cls == "ramified":
            continue
        from eabc_from_lean import class_of

        ec = class_of(p)
        if ec is None:
            continue
        label = ec.value
        coarse_match = (split_cls == "split" and label in EB_CLASSES) or (
            split_cls == "inert" and label in AC_CLASSES
        )
        if not coarse_match:
            mismatches += 1
        total += 1
        if len(rows) < 10:
            rows.append(
                {
                    "p": p,
                    "mod3": p % 3,
                    "mod12": p % 12,
                    "eisenstein": split_cls,
                    "kappa": label,
                    "mod3_bucket": eabc_mod3_bucket(label),
                    "coarse_match": coarse_match,
                }
            )
    return {
        "max_p": max_p,
        "prime_count": total,
        "mismatches": mismatches,
        "exact_coarse_bipartite": mismatches == 0,
        "bipartition": "split ↔ E∪B (≡1 mod 3), inert ↔ A∪C (≡2 mod 3)",
        "sample_rows": rows,
    }


def spaltung_report(max_p: int) -> dict[str, Any]:
    rows: list[EisensteinSplitRow] = []
    gamma_counts: Counter[tuple[str, str]] = Counter()
    marginal_a: Counter[str] = Counter()
    marginal_b: Counter[str] = Counter()

    for p in _sieve_primes(max_p):
        row = classify_eisenstein_split_prime(p)
        if row is None:
            continue
        rows.append(row)
        gamma_counts[row.gamma] += 1
        marginal_a[row.kappa_a] += 1
        marginal_b[row.kappa_b] += 1

    n = len(rows)
    mu = _mu_table(gamma_counts, n)
    chi2 = _chi_square_uniform(gamma_counts, n)
    top_dev = _top_deviations(mu)
    max_dev = top_dev[0]["abs_delta"] if top_dev else 0.0

    pairs = [r.gamma for r in rows]
    null = _shuffle_null_chi2(pairs)
    z_shuffle = 0.0
    if null["std"] > 0:
        z_shuffle = (chi2 - null["mean"]) / null["std"]

    return {
        "max_p": max_p,
        "split_count": n,
        "mu_X": mu,
        "counts": {gamma_label(*g): gamma_counts.get(g, 0) for g in GAMMA_PAIRS},
        "marginal_kappa_a_prime": dict(marginal_a),
        "marginal_kappa_b_prime": dict(marginal_b),
        "chi2_uniform_16": chi2,
        "chi2_critical_15df_005": CHI2_CRIT_15DF_005,
        "max_deviation_from_uniform": max_dev,
        "top_deviations": top_dev,
        "shuffle_null_chi2": null,
        "z_score_chi2_vs_shuffle": z_shuffle,
        "verdict": _verdict(chi2, z_shuffle, max_dev),
        "coarse_mod3_defekt": coarse_mod3_defekt_report(max_p),
        "sample_rows": [r.to_dict() for r in rows[:20]],
    }


def multi_scale_report(max_ps: list[int]) -> dict[str, Any]:
    scales = {str(x): spaltung_report(x) for x in sorted(set(max_ps))}
    return {
        "meta": {
            "hypothesis_doc": "collatz_eabc_eisenstein_spaltung.md",
            "gaussian_reference": "collatz_eabc_gauss_spaltung_test.py",
            "ring": "Z[ω], ω=e^(2πi/3), N(a+bω)=a²−ab+b²",
            "split_condition": "p ≡ 1 (mod 3), p ≠ 3",
            "glatt_note": (
                "Gleiche strip_smooth wie Gauß: 2^α 3^β abziehen; "
                "3 ramifiziert in Z[ω], aber κ bleibt mod-12-basiert"
            ),
            "uniform_null": UNIFORM_MU,
            "gamma_space_size": len(GAMMA_PAIRS),
            "scales": sorted(set(max_ps)),
        },
        "scales": scales,
    }


def format_table(report: dict[str, Any]) -> str:
    lines = [
        "Eisenstein–EABC-Spaltung (glatt-EABC)",
        "=" * 55,
    ]
    for key, scale in report["scales"].items():
        coarse = scale["coarse_mod3_defekt"]
        lines.extend(
            [
                f"\nX = {key}  (split p>3: {scale['split_count']})",
                f"mod-3 bipartit exakt: {coarse['exact_coarse_bipartite']} "
                f"({coarse['mismatches']} Mismatches / {coarse['prime_count']})",
                f"χ² vs 1/16: {scale['chi2_uniform_16']:.3f}  "
                f"(krit. {scale['chi2_critical_15df_005']:.1f})",
                f"z vs Shuffle-Null: {scale['z_score_chi2_vs_shuffle']:.2f}",
                f"max |μ−1/16|: {scale['max_deviation_from_uniform']:.5f}",
                "Top-Abweichungen:",
            ]
        )
        for row in scale["top_deviations"]:
            lines.append(
                f"  {row['gamma']}: μ={row['mu_X']:.5f}  Δ={row['delta']:+.5f}"
            )
        lines.append(f"verdict: {scale['verdict']}")
    return "\n".join(lines)


def run(max_ps: list[int] | None = None, output: Path | None = None) -> dict[str, Any]:
    scales = max_ps or [10_000, 100_000, 1_000_000]
    report = multi_scale_report(scales)
    out = output or DEFAULT_OUTPUT
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report["output_path"] = str(out)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Eisenstein–EABC-Spaltungstest (glatt-EABC)")
    parser.add_argument(
        "--max-p",
        type=int,
        nargs="+",
        default=[10_000, 100_000, 1_000_000],
        help="eine oder mehrere Obergrenzen X",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="JSON-Ausgabedatei",
    )
    args = parser.parse_args()
    report = run(max_ps=args.max_p, output=args.output)
    print(format_table(report))
    print(f"\nJSON: {report['output_path']}")


if __name__ == "__main__":
    main()
