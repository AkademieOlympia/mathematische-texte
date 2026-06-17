#!/usr/bin/env python3
"""
Gauß–EABC-Spaltungstest mit glatt-EABC-Zerlegung.

Kanonsiche Hypothese: collatz_eabc_gauss_spaltung_hypothese.md

Für split-Primzahlen p ≡ 1 (mod 4), p = a² + b² (kanonisch 0 < a ≤ b):
  a = 2^α_a · 3^β_a · a',  gcd(a', 6) = 1
  b = 2^α_b · 3^β_b · b',  gcd(b', 6) = 1
  Γ(p) = (κ(a'), κ(b')) ∈ {E,A,B,C}²

μ_X(γ) = #{p ≤ X split : Γ(p)=γ} / #{p ≤ X split}

Kritischer Fix gegenüber collatz_eabc_gauss_faktor_eabc_test.py:
  glatte 2^α 3^β-Anteile abziehen → beide Legs EABC-sichtbar (16 Klassen).

Ausführung:
    python3 collatz_eabc_gauss_spaltung_test.py
    python3 collatz_eabc_gauss_spaltung_test.py --max-p 100000 1000000
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

from eabc_from_lean import class_of

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_gauss_spaltung.json"

EABC_LABELS = ("E", "A", "B", "C")
GAMMA_PAIRS: tuple[tuple[str, str], ...] = tuple(
    (ka, kb) for ka in EABC_LABELS for kb in EABC_LABELS
)
UNIFORM_MU = 1.0 / len(GAMMA_PAIRS)
CHI2_CRIT_15DF_005 = 25.0
NULL_TRIALS = 5000


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


def gaussian_factor_pair(p: int) -> tuple[int, int] | None:
    """Kanonisch 0 < a ≤ b mit p = a² + b²; nur für p ≡ 1 (mod 4)."""
    if p <= 0 or p % 4 != 1:
        return None
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        b = isqrt(b2)
        if b > 0 and b * b == b2:
            return (a, b) if a <= b else (b, a)
    return None


def strip_smooth(n: int) -> tuple[int, int, int]:
    """
    Zerlege n = 2^α · 3^β · n' mit gcd(n', 6) = 1.
    Rückgabe: (α, β, n').
    """
    if n < 1:
        raise ValueError("n muss positiv sein")
    alpha = beta = 0
    x = n
    while x % 2 == 0:
        alpha += 1
        x //= 2
    while x % 3 == 0:
        beta += 1
        x //= 3
    return alpha, beta, x


def kappa_glatt(n: int) -> tuple[int, int, int, str]:
    """(α, β, n', κ(n'))."""
    alpha, beta, core = strip_smooth(n)
    ec = class_of(core)
    if ec is None:
        raise ValueError(f"κ undefiniert für glatter Kern {core} von n={n}")
    return alpha, beta, core, ec.value


def gamma_label(ka: str, kb: str) -> str:
    return f"({ka},{kb})"


@dataclass(frozen=True, slots=True)
class SplitRow:
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


def classify_split_prime(p: int) -> SplitRow | None:
    if p <= 3 or p % 4 != 1:
        return None
    pair = gaussian_factor_pair(p)
    if pair is None:
        return None
    a, b = pair
    alpha_a, beta_a, a_prime, ka = kappa_glatt(a)
    alpha_b, beta_b, b_prime, kb = kappa_glatt(b)
    return SplitRow(
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


def _chi_square_uniform(counts: Counter[tuple[str, str]], n: int) -> float:
    if n == 0:
        return 0.0
    expected = n * UNIFORM_MU
    return sum((counts.get(g, 0) - expected) ** 2 / expected for g in GAMMA_PAIRS)


def _mu_table(counts: Counter[tuple[str, str]], n: int) -> dict[str, float]:
    if n == 0:
        return {gamma_label(*g): 0.0 for g in GAMMA_PAIRS}
    return {gamma_label(*g): counts.get(g, 0) / n for g in GAMMA_PAIRS}


def _top_deviations(mu: dict[str, float], k: int = 5) -> list[dict[str, Any]]:
    ranked = sorted(
        (
            {
                "gamma": g,
                "mu_X": mu[g],
                "delta": mu[g] - UNIFORM_MU,
                "abs_delta": abs(mu[g] - UNIFORM_MU),
            }
            for g in mu
        ),
        key=lambda row: (-row["abs_delta"], row["gamma"]),
    )
    return ranked[:k]


def _shuffle_null_chi2(
    pairs: list[tuple[str, str]],
    trials: int = NULL_TRIALS,
    seed: int = 0,
) -> dict[str, float]:
    """Marginal-erhaltende Null: κ_a und κ_b unabhängig permutiert."""
    if not pairs:
        return {"mean": 0.0, "std": 0.0, "trials": 0}
    rng = random.Random(seed)
    n = len(pairs)
    legs_a = [ka for ka, _ in pairs]
    legs_b = [kb for _, kb in pairs]
    chi2_vals: list[float] = []
    for _ in range(trials):
        sh_a = legs_a[:]
        sh_b = legs_b[:]
        rng.shuffle(sh_a)
        rng.shuffle(sh_b)
        c: Counter[tuple[str, str]] = Counter(zip(sh_a, sh_b, strict=True))
        chi2_vals.append(_chi_square_uniform(c, n))
    return {
        "mean": statistics.mean(chi2_vals),
        "std": statistics.stdev(chi2_vals) if len(chi2_vals) > 1 else 0.0,
        "trials": trials,
    }


def _verdict(chi2: float, z_shuffle: float, max_dev: float) -> str:
    if chi2 > CHI2_CRIT_15DF_005 and abs(z_shuffle) > 2.0:
        return (
            "bias_detected: Γ(p) weicht signifikant von Uniform 1/16 ab "
            "(χ² + Shuffle-Null) → Anker für Z[i]→EABC-Brücke"
        )
    if chi2 > CHI2_CRIT_15DF_005:
        return (
            "marginal_chi2_signal: χ² über Schwellwert, aber Shuffle-Null nicht extrem; "
            "vorsichtig interpretieren"
        )
    if max_dev < 0.01:
        return (
            "approximately_uniform: μ_X ≈ 1/16 über alle 16 Γ-Klassen "
            "→ schwächt Spaltungs-Orientierungskonjektur"
        )
    return (
        "no_significant_bias: Fluktuationen konsistent mit Uniform 1/16 "
        "und Shuffle-Null → Konjektur nicht gestützt"
    )


def spaltung_report(max_p: int) -> dict[str, Any]:
    rows: list[SplitRow] = []
    gamma_counts: Counter[tuple[str, str]] = Counter()
    marginal_a: Counter[str] = Counter()
    marginal_b: Counter[str] = Counter()

    for p in _sieve_primes(max_p):
        row = classify_split_prime(p)
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
        "sample_rows": [r.to_dict() for r in rows[:20]],
    }


def multi_scale_report(max_ps: list[int]) -> dict[str, Any]:
    scales = {str(x): spaltung_report(x) for x in sorted(set(max_ps))}
    return {
        "meta": {
            "hypothesis_doc": "collatz_eabc_gauss_spaltung_hypothese.md",
            "supersedes": "collatz_eabc_gauss_faktor_eabc_test.py",
            "fix_note": (
                "Glatt-EABC: 2^α 3^β abziehen vor κ; beide Faktorlegs EABC-sichtbar; "
                "16 Γ-Klassen statt höchstens einer sichtbaren Leg"
            ),
            "uniform_null": UNIFORM_MU,
            "gamma_space_size": len(GAMMA_PAIRS),
            "scales": sorted(set(max_ps)),
        },
        "scales": scales,
    }


def format_table(report: dict[str, Any]) -> str:
    lines = [
        "Gauß–EABC-Spaltung (glatt-EABC)",
        "=" * 55,
    ]
    for key, scale in report["scales"].items():
        lines.extend(
            [
                f"\nX = {key}  (split p>3: {scale['split_count']})",
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
    parser = argparse.ArgumentParser(description="Gauß–EABC-Spaltungstest (glatt-EABC)")
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
