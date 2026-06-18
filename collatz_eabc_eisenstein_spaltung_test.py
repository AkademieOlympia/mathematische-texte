#!/usr/bin/env python3
"""
Eisenstein–EABC-Spaltungstest mit glatt-EABC-Zerlegung (volle Γ-Signatur).

Kanonsiche Hypothese: collatz_eabc_eisenstein_spaltung.md

Für split-Primzahlen p ≡ 1 (mod 3), p ≠ 3, p = a² − ab + b² (kanonisch 0 < a ≤ b):
  a = 2^{α_a} 3^{β_a} a',  gcd(a', 6) = 1
  b = 2^{α_b} 3^{β_b} b',  gcd(b', 6) = 1
  Γ_E(p) = ((α_a, β_a, κ(a')), (α_b, β_b, κ(b')))

Ausführung:
    python3 collatz_eabc_eisenstein_spaltung_test.py
    python3 collatz_eabc_eisenstein_spaltung_test.py --max-p 100000 1000000
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from dataclasses import dataclass
from math import isqrt
from pathlib import Path
from typing import Any

from collatz_eabc_gauss_spaltung_test import (
    CHI2_CRIT_15DF_005,
    CHI2_CRIT_3DF_005,
    GAMMA_PAIRS,
    FullGamma,
    UNIFORM_MARGINAL,
    UNIFORM_MU,
    _chi_square_uniform_pairs,
    _eabc_marginal_bias,
    _marginal_eabc_conditional,
    _mu_table,
    _shuffle_null_chi2,
    _shuffle_null_conditional_chi2,
    _sieve_primes,
    _top_deviations,
    full_gamma_label,
    gamma_label,
    kappa_glatt,
    smooth_pattern_key,
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
    full_gamma: FullGamma

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
            "full_gamma": list(self.full_gamma),
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
    full = (alpha_a, beta_a, ka, alpha_b, beta_b, kb)
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
        full_gamma=full,
    )


def _parity_smooth_summary_eisenstein(rows: list[EisensteinSplitRow]) -> dict[str, Any]:
    a_even = sum(1 for r in rows if r.alpha_a >= 1 and r.alpha_b == 0)
    b_even = sum(1 for r in rows if r.alpha_a == 0 and r.alpha_b >= 1)
    both_even = sum(1 for r in rows if r.alpha_a >= 1 and r.alpha_b >= 1)
    both_odd = sum(1 for r in rows if r.alpha_a == 0 and r.alpha_b == 0)
    one_even = a_even + b_even
    return {
        "exactly_one_even": one_even,
        "both_even": both_even,
        "both_odd": both_odd,
        "a_even_b_odd": a_even,
        "a_odd_b_even": b_even,
        "note": (
            "Eisenstein: p≡1 mod 3 erlaubt beide ungerade Legs "
            "(anders als Gauß p≡1 mod 4)"
        ),
    }


def _smooth_pattern_counts(rows: list[EisensteinSplitRow]) -> dict[str, Any]:
    counts: Counter[tuple[int, int, int, int]] = Counter(
        smooth_pattern_key(r) for r in rows
    )
    n = len(rows)
    top = sorted(
        (
            {
                "pattern": f"(α_a={p[0]},β_a={p[1]},α_b={p[2]},β_b={p[3]})",
                "count": c,
                "fraction": c / n if n else 0.0,
            }
            for p, c in counts.items()
        ),
        key=lambda row: (-row["count"], row["pattern"]),
    )
    return {
        "distinct_patterns": len(counts),
        "top_patterns": top[:8],
        "counts": {str(p): c for p, c in counts.most_common()},
    }


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


def _verdict_eisenstein(
    chi2_joint: float,
    z_shuffle: float,
    z_conditional: float,
    eabc_bias: dict[str, Any],
    max_dev_joint: float,
) -> str:
    marginals_uniform = (
        eabc_bias["chi2_marginal_a_vs_quarter"] < CHI2_CRIT_3DF_005
        and eabc_bias["chi2_marginal_b_vs_quarter"] < CHI2_CRIT_3DF_005
    )
    joint_sig = chi2_joint > CHI2_CRIT_15DF_005 and abs(z_shuffle) > 2.0
    conditional_sig = abs(z_conditional) > 2.0

    if joint_sig and conditional_sig:
        return (
            "eabc_anchor: stabile (κ_a',κ_b')-Abweichung nach (α,β)-Stratifizierung "
            "→ reale Z[ω]→EABC-Orientierung"
        )
    if marginals_uniform and not conditional_sig:
        return (
            "falsification_weak: κ(a'), κ(b') nach glatt-strip marginal ~ uniform 1/4; "
            "keine robuste bedingte Kopplung → Konjektur nicht gestützt"
        )
    if chi2_joint > CHI2_CRIT_15DF_005 and not conditional_sig:
        return (
            "marginal_driven: erhöhte χ² vs 1/16 durch Randverteilungen, "
            "nicht durch bedingte EABC-Kopplung"
        )
    if max_dev_joint < 0.01:
        return (
            "approximately_uniform: μ_X ≈ 1/16 → schwächt Eisenstein-Spaltungskonjektur"
        )
    return "inconclusive: mehr Skala oder stärkere Signale nötig"


def spaltung_report(max_p: int) -> dict[str, Any]:
    rows: list[EisensteinSplitRow] = []
    gamma_counts: Counter[tuple[str, str]] = Counter()
    full_counts: Counter[FullGamma] = Counter()

    for p in _sieve_primes(max_p):
        row = classify_eisenstein_split_prime(p)
        if row is None:
            continue
        rows.append(row)
        gamma_counts[row.gamma] += 1
        full_counts[row.full_gamma] += 1

    n = len(rows)
    mu_joint = _mu_table(gamma_counts, n)
    chi2_joint = _chi_square_uniform_pairs(gamma_counts, n)
    top_dev_joint = _top_deviations(mu_joint, UNIFORM_MU)
    mu_full = (
        {full_gamma_label(g): full_counts[g] / n for g in full_counts} if n else {}
    )
    top_dev_full = _top_deviations(mu_full, 1.0 / n if n else 0.0, k=5)
    max_dev_joint = top_dev_joint[0]["abs_delta"] if top_dev_joint else 0.0

    pairs = [r.gamma for r in rows]
    null = _shuffle_null_chi2(pairs)
    z_shuffle = 0.0
    if null["std"] > 0:
        z_shuffle = (chi2_joint - null["mean"]) / null["std"]

    null_cond = _shuffle_null_conditional_chi2(rows)  # type: ignore[arg-type]
    z_conditional = 0.0
    if null_cond["std"] > 0:
        z_conditional = (chi2_joint - null_cond["mean"]) / null_cond["std"]

    eabc_bias = _eabc_marginal_bias(rows)  # type: ignore[arg-type]

    return {
        "max_p": max_p,
        "split_count": n,
        "full_signature_space_size": len(full_counts),
        "smooth_pattern_space_size": _smooth_pattern_counts(rows)["distinct_patterns"],
        "parity_smooth": _parity_smooth_summary_eisenstein(rows),
        "smooth_patterns": _smooth_pattern_counts(rows),
        "mu_X_joint_kappa": mu_joint,
        "mu_X_full_gamma": mu_full,
        "counts_joint": {gamma_label(*g): gamma_counts.get(g, 0) for g in GAMMA_PAIRS},
        "marginal_eabc_conditional_on_smooth": _marginal_eabc_conditional(rows),  # type: ignore[arg-type]
        "eabc_marginal_bias": eabc_bias,
        "chi2_joint_kappa_vs_16": chi2_joint,
        "chi2_critical_15df_005": CHI2_CRIT_15DF_005,
        "max_deviation_joint_from_16": max_dev_joint,
        "top_deviations_joint": top_dev_joint,
        "top_deviations_full_gamma": top_dev_full,
        "shuffle_null_chi2": null,
        "shuffle_null_conditional_chi2": null_cond,
        "z_score_chi2_vs_shuffle": z_shuffle,
        "z_score_chi2_vs_conditional_shuffle": z_conditional,
        "verdict": _verdict_eisenstein(
            chi2_joint, z_shuffle, z_conditional, eabc_bias, max_dev_joint
        ),
        "coarse_mod3_defekt": coarse_mod3_defekt_report(max_p),
        "sample_rows": [r.to_dict() for r in rows[:20]],
    }


def multi_scale_report(max_ps: list[int]) -> dict[str, Any]:
    scales = {str(x): spaltung_report(x) for x in sorted(set(max_ps))}
    largest = scales[str(max(sorted(set(max_ps))))]
    return {
        "meta": {
            "hypothesis_doc": "collatz_eabc_eisenstein_spaltung.md",
            "gaussian_reference": "collatz_eabc_gauss_spaltung_test.py",
            "ring": "Z[ω], ω=e^(2πi/3), N(a+bω)=a²−ab+b²",
            "split_condition": "p ≡ 1 (mod 3), p ≠ 3",
            "gamma_definition": (
                "Γ_E(p)=((α_a,β_a,κ(a')),(α_b,β_b,κ(b'))) "
                "kompakt (ν_2(a),ν_3(a),κ(a'),ν_2(b),ν_3(b),κ(b'))"
            ),
            "uniform_null_joint": UNIFORM_MU,
            "uniform_null_marginal": UNIFORM_MARGINAL,
            "joint_kappa_space_size": len(GAMMA_PAIRS),
            "scales": sorted(set(max_ps)),
            "largest_scale_full_signature_space_size": largest[
                "full_signature_space_size"
            ],
        },
        "scales": scales,
    }


def format_table(report: dict[str, Any]) -> str:
    lines = [
        "Eisenstein–EABC-Spaltung (volle Γ-Signatur, glatt-EABC)",
        "=" * 60,
    ]
    for key, scale in report["scales"].items():
        coarse = scale["coarse_mod3_defekt"]
        parity = scale["parity_smooth"]
        eabc = scale["eabc_marginal_bias"]
        lines.extend(
            [
                f"\nX = {key}  (split p>3: {scale['split_count']})",
                f"mod-3 bipartit exakt: {coarse['exact_coarse_bipartite']} "
                f"({coarse['mismatches']} Mismatches / {coarse['prime_count']})",
                f"|Γ|_beobachtet: {scale['full_signature_space_size']}  "
                f"Parität 1 gerade: {parity['exactly_one_even']}/{scale['split_count']}",
                f"χ² κ-Paar vs 1/16: {scale['chi2_joint_kappa_vs_16']:.3f}  "
                f"z_cond={scale['z_score_chi2_vs_conditional_shuffle']:.2f}",
                f"χ² κ(a') vs 1/4: {eabc['chi2_marginal_a_vs_quarter']:.3f}  "
                f"κ(b') vs 1/4: {eabc['chi2_marginal_b_vs_quarter']:.3f}",
                "Top-5 κ-Paar-Abweichungen:",
            ]
        )
        for row in scale["top_deviations_joint"][:5]:
            lines.append(
                f"  {row['gamma']}: μ={row['mu_X']:.5f}  Δ={row['delta']:+.5f}"
            )
        lines.append(f"verdict: {scale['verdict']}")
    return "\n".join(lines)


def run(max_ps: list[int] | None = None, output: Path | None = None) -> dict[str, Any]:
    scales = max_ps or [10_000, 100_000, 1_000_000]
    report = multi_scale_report(scales)
    out = output or DEFAULT_OUTPUT
    out.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    report["output_path"] = str(out)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Eisenstein–EABC-Spaltungstest (volle Γ)"
    )
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
