#!/usr/bin/env python3
"""
Gauß–EABC-Spaltungstest mit glatt-EABC-Zerlegung (volle Γ-Signatur).

Kanonsiche Hypothese: collatz_eabc_gauss_spaltung_hypothese.md

Für split-Primzahlen p ≡ 1 (mod 4), p = a² + b² (kanonisch 0 < a ≤ b):
  a = 2^{α_a} 3^{β_a} a',  gcd(a', 6) = 1
  b = 2^{α_b} 3^{β_b} b',  gcd(b', 6) = 1
  Γ(p) = ((α_a, β_a, κ(a')), (α_b, β_b, κ(b')))
  kompakt: (ν_2(a), ν_3(a), κ(a'), ν_2(b), ν_3(b), κ(b'))

Parität: genau eine der Legs ist gerade (α=ν_2 ≥ 1 auf genau einer Seite).

μ_X(γ) = #{p ≤ X split : Γ(p)=γ} / #{p ≤ X split}

Ausführung:
    python3 collatz_eabc_gauss_spaltung_test.py
    python3 collatz_eabc_gauss_spaltung_test.py --max-p 100000 1000000
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
from collections import Counter, defaultdict
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
UNIFORM_MARGINAL = 0.25
CHI2_CRIT_15DF_005 = 25.0
CHI2_CRIT_3DF_005 = 7.81
NULL_TRIALS = 5000
MIN_SMOOTH_BUCKET = 8

FullGamma = tuple[int, int, str, int, int, str]
SmoothPattern = tuple[int, int, int, int]


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


def full_gamma_label(sig: FullGamma) -> str:
    return f"({sig[0]},{sig[1]},{sig[2]},{sig[3]},{sig[4]},{sig[5]})"


def smooth_pattern_key(row: Any) -> SmoothPattern:
    return row.alpha_a, row.beta_a, row.alpha_b, row.beta_b


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


def classify_split_prime(p: int) -> SplitRow | None:
    if p <= 3 or p % 4 != 1:
        return None
    pair = gaussian_factor_pair(p)
    if pair is None:
        return None
    a, b = pair
    alpha_a, beta_a, a_prime, ka = kappa_glatt(a)
    alpha_b, beta_b, b_prime, kb = kappa_glatt(b)
    full = (alpha_a, beta_a, ka, alpha_b, beta_b, kb)
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
        full_gamma=full,
    )


def _chi_square_uniform(
    counts: Counter[Any],
    labels: tuple[Any, ...],
    n: int,
    expected_frac: float,
) -> float:
    if n == 0:
        return 0.0
    expected = n * expected_frac
    return sum((counts.get(g, 0) - expected) ** 2 / expected for g in labels)


def _chi_square_uniform_pairs(counts: Counter[tuple[str, str]], n: int) -> float:
    return _chi_square_uniform(counts, GAMMA_PAIRS, n, UNIFORM_MU)


def _mu_table(counts: Counter[tuple[str, str]], n: int) -> dict[str, float]:
    if n == 0:
        return {gamma_label(*g): 0.0 for g in GAMMA_PAIRS}
    return {gamma_label(*g): counts.get(g, 0) / n for g in GAMMA_PAIRS}


def _mu_full_table(counts: Counter[FullGamma], n: int) -> dict[str, float]:
    if n == 0:
        return {}
    return {full_gamma_label(g): counts.get(g, 0) / n for g in counts}


def _top_deviations(
    mu: dict[str, float],
    uniform: float,
    k: int = 5,
) -> list[dict[str, Any]]:
    ranked = sorted(
        (
            {
                "gamma": g,
                "mu_X": mu[g],
                "delta": mu[g] - uniform,
                "abs_delta": abs(mu[g] - uniform),
            }
            for g in mu
        ),
        key=lambda row: (-row["abs_delta"], row["gamma"]),
    )
    return ranked[:k]


def _parity_smooth_summary(rows: list[SplitRow]) -> dict[str, Any]:
    n = len(rows)
    if n == 0:
        return {
            "exactly_one_even": 0,
            "both_even": 0,
            "both_odd": 0,
            "a_even_b_odd": 0,
            "a_odd_b_even": 0,
            "parity_constraint_holds": True,
            "note": "p ungerade ⇒ genau eine Leg gerade (ν_2≥1)",
        }
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
        "parity_constraint_holds": one_even == n and both_even == 0 and both_odd == 0,
        "even_leg_alpha_distribution": dict(
            Counter(r.alpha_a if r.alpha_a >= 1 else r.alpha_b for r in rows)
        ),
        "note": "p ungerade ⇒ genau eine Leg gerade; ν_2≥1 = glatte 2-Komponente",
    }


def _smooth_pattern_counts(rows: list[SplitRow]) -> dict[str, Any]:
    counts: Counter[SmoothPattern] = Counter(smooth_pattern_key(r) for r in rows)
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


def _marginal_eabc_conditional(
    rows: list[SplitRow],
    min_bucket: int = MIN_SMOOTH_BUCKET,
) -> dict[str, Any]:
    buckets: dict[SmoothPattern, list[tuple[str, str]]] = defaultdict(list)
    for r in rows:
        buckets[smooth_pattern_key(r)].append(r.gamma)

    bucket_reports: list[dict[str, Any]] = []
    pooled_chi2 = 0.0
    pooled_df = 0

    for pattern, pairs in sorted(buckets.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        n = len(pairs)
        if n < min_bucket:
            continue
        counts: Counter[tuple[str, str]] = Counter(pairs)
        mu = _mu_table(counts, n)
        chi2 = _chi_square_uniform_pairs(counts, n)
        top = _top_deviations(mu, UNIFORM_MU, k=3)
        bucket_reports.append(
            {
                "smooth_pattern": {
                    "alpha_a": pattern[0],
                    "beta_a": pattern[1],
                    "alpha_b": pattern[2],
                    "beta_b": pattern[3],
                },
                "count": n,
                "chi2_vs_uniform_16": chi2,
                "top_deviations": top,
                "mu_X": mu,
            }
        )
        pooled_chi2 += chi2
        pooled_df += 15

    return {
        "min_bucket": min_bucket,
        "buckets_reported": len(bucket_reports),
        "pooled_chi2_vs_uniform_16": pooled_chi2,
        "pooled_df": pooled_df,
        "buckets": bucket_reports[:12],
    }


def _eabc_marginal_bias(rows: list[SplitRow]) -> dict[str, Any]:
    n = len(rows)
    marginal_a: Counter[str] = Counter(r.kappa_a for r in rows)
    marginal_b: Counter[str] = Counter(r.kappa_b for r in rows)
    joint: Counter[tuple[str, str]] = Counter(r.gamma for r in rows)

    chi2_a = _chi_square_uniform(marginal_a, EABC_LABELS, n, UNIFORM_MARGINAL)
    chi2_b = _chi_square_uniform(marginal_b, EABC_LABELS, n, UNIFORM_MARGINAL)
    chi2_joint = _chi_square_uniform_pairs(joint, n)

    return {
        "marginal_kappa_a_prime": dict(marginal_a),
        "marginal_kappa_b_prime": dict(marginal_b),
        "chi2_marginal_a_vs_quarter": chi2_a,
        "chi2_marginal_b_vs_quarter": chi2_b,
        "chi2_joint_vs_sixteenth": chi2_joint,
        "chi2_critical_3df_005": CHI2_CRIT_3DF_005,
        "chi2_critical_15df_005": CHI2_CRIT_15DF_005,
        "max_marginal_deviation_a": (
            max(abs(marginal_a.get(k, 0) / n - UNIFORM_MARGINAL) for k in EABC_LABELS)
            if n
            else 0.0
        ),
        "max_marginal_deviation_b": (
            max(abs(marginal_b.get(k, 0) / n - UNIFORM_MARGINAL) for k in EABC_LABELS)
            if n
            else 0.0
        ),
    }


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
        chi2_vals.append(_chi_square_uniform_pairs(c, n))
    return {
        "mean": statistics.mean(chi2_vals),
        "std": statistics.stdev(chi2_vals) if len(chi2_vals) > 1 else 0.0,
        "trials": trials,
    }


def _shuffle_null_conditional_chi2(
    rows: list[SplitRow],
    trials: int = NULL_TRIALS,
    seed: int = 0,
) -> dict[str, float]:
    """Null: innerhalb jedes (α,β)-Musters κ_a, κ_b unabhängig permutiert."""
    if not rows:
        return {"mean": 0.0, "std": 0.0, "trials": 0}
    rng = random.Random(seed)
    buckets: dict[SmoothPattern, list[tuple[str, str]]] = defaultdict(list)
    for r in rows:
        buckets[smooth_pattern_key(r)].append(r.gamma)

    n = len(rows)
    chi2_vals: list[float] = []
    for _ in range(trials):
        shuffled: list[tuple[str, str]] = []
        for pairs in buckets.values():
            if len(pairs) < 2:
                shuffled.extend(pairs)
                continue
            legs_a = [ka for ka, _ in pairs]
            legs_b = [kb for _, kb in pairs]
            rng.shuffle(legs_a)
            rng.shuffle(legs_b)
            shuffled.extend(zip(legs_a, legs_b, strict=True))
        c = Counter(shuffled)
        chi2_vals.append(_chi_square_uniform_pairs(c, n))
    return {
        "mean": statistics.mean(chi2_vals),
        "std": statistics.stdev(chi2_vals) if len(chi2_vals) > 1 else 0.0,
        "trials": trials,
    }


def _verdict(
    chi2_joint: float,
    z_shuffle: float,
    z_conditional: float,
    eabc_bias: dict[str, Any],
    max_dev_joint: float,
    parity_ok: bool,
) -> str:
    marginals_uniform = (
        eabc_bias["chi2_marginal_a_vs_quarter"] < CHI2_CRIT_3DF_005
        and eabc_bias["chi2_marginal_b_vs_quarter"] < CHI2_CRIT_3DF_005
    )
    joint_sig = chi2_joint > CHI2_CRIT_15DF_005 and abs(z_shuffle) > 2.0
    conditional_sig = abs(z_conditional) > 2.0

    if not parity_ok:
        return "parity_violation: genau-eine-gerade-Regel verletzt — Datenfehler prüfen"

    if joint_sig and conditional_sig:
        return (
            "eabc_anchor: stabile (κ_a',κ_b')-Abweichung auch nach (α,β)-Stratifizierung "
            "→ reale Z[i]→EABC-Orientierung jenseits glatter Parität"
        )
    if marginals_uniform and not conditional_sig:
        return (
            "falsification_weak: κ(a'), κ(b') nach glatt-strip marginal ~ uniform 1/4; "
            "keine robuste Kopplung über (α,β)-Muster → Konjektur nicht gestützt"
        )
    if chi2_joint > CHI2_CRIT_15DF_005 and not conditional_sig:
        return (
            "marginal_driven: erhöhte χ² vs 1/16 erklärt durch Randverteilungen, "
            "nicht durch bedingte EABC-Kopplung — vorsichtig interpretieren"
        )
    if max_dev_joint < 0.01:
        return (
            "approximately_uniform: μ_X ≈ 1/16 auf κ-Paar und volle Γ-Signatur dünn "
            "→ schwächt Spaltungs-Orientierungskonjektur"
        )
    return (
        "inconclusive: Fluktuationen nicht eindeutig von Uniform-Null oder "
        "Shuffle-Null trennbar — mehr Skala nötig"
    )


def spaltung_report(max_p: int) -> dict[str, Any]:
    rows: list[SplitRow] = []
    gamma_counts: Counter[tuple[str, str]] = Counter()
    full_counts: Counter[FullGamma] = Counter()

    for p in _sieve_primes(max_p):
        row = classify_split_prime(p)
        if row is None:
            continue
        rows.append(row)
        gamma_counts[row.gamma] += 1
        full_counts[row.full_gamma] += 1

    n = len(rows)
    mu_joint = _mu_table(gamma_counts, n)
    mu_full = _mu_full_table(full_counts, n)
    chi2_joint = _chi_square_uniform_pairs(gamma_counts, n)
    top_dev_joint = _top_deviations(mu_joint, UNIFORM_MU)
    top_dev_full = _top_deviations(mu_full, 1.0 / n if n else 0.0, k=5)
    max_dev_joint = top_dev_joint[0]["abs_delta"] if top_dev_joint else 0.0

    pairs = [r.gamma for r in rows]
    null = _shuffle_null_chi2(pairs)
    z_shuffle = 0.0
    if null["std"] > 0:
        z_shuffle = (chi2_joint - null["mean"]) / null["std"]

    null_cond = _shuffle_null_conditional_chi2(rows)
    z_conditional = 0.0
    if null_cond["std"] > 0:
        z_conditional = (chi2_joint - null_cond["mean"]) / null_cond["std"]

    parity = _parity_smooth_summary(rows)
    smooth = _smooth_pattern_counts(rows)
    conditional = _marginal_eabc_conditional(rows)
    eabc_bias = _eabc_marginal_bias(rows)

    return {
        "max_p": max_p,
        "split_count": n,
        "full_signature_space_size": len(full_counts),
        "smooth_pattern_space_size": smooth["distinct_patterns"],
        "parity_smooth": parity,
        "smooth_patterns": smooth,
        "mu_X_joint_kappa": mu_joint,
        "mu_X_full_gamma": mu_full,
        "counts_joint": {gamma_label(*g): gamma_counts.get(g, 0) for g in GAMMA_PAIRS},
        "counts_full_top20": dict(
            sorted(
                ((full_gamma_label(g), c) for g, c in full_counts.items()),
                key=lambda kv: (-kv[1], kv[0]),
            )[:20]
        ),
        "marginal_eabc_conditional_on_smooth": conditional,
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
        "verdict": _verdict(
            chi2_joint,
            z_shuffle,
            z_conditional,
            eabc_bias,
            max_dev_joint,
            parity["parity_constraint_holds"],
        ),
        "sample_rows": [r.to_dict() for r in rows[:20]],
    }


def multi_scale_report(max_ps: list[int]) -> dict[str, Any]:
    scales = {str(x): spaltung_report(x) for x in sorted(set(max_ps))}
    largest = scales[str(max(sorted(set(max_ps))))]
    return {
        "meta": {
            "hypothesis_doc": "collatz_eabc_gauss_spaltung_hypothese.md",
            "supersedes": "collatz_eabc_gauss_faktor_eabc_test.py",
            "fix_note": (
                "Volle Γ-Signatur (ν_2,ν_3,κ) pro Leg; Parität: genau eine Leg gerade; "
                "EABC-Marginaltest nach glatt-strip"
            ),
            "gamma_definition": (
                "Γ(p)=((α_a,β_a,κ(a')),(α_b,β_b,κ(b'))) "
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
        "Gauß–EABC-Spaltung (volle Γ-Signatur, glatt-EABC)",
        "=" * 60,
    ]
    for key, scale in report["scales"].items():
        parity = scale["parity_smooth"]
        eabc = scale["eabc_marginal_bias"]
        lines.extend(
            [
                f"\nX = {key}  (split p>3: {scale['split_count']})",
                f"|Γ|_beobachtet: {scale['full_signature_space_size']}  "
                f"|(α,β)-Muster|: {scale['smooth_pattern_space_size']}",
                f"Parität (genau eine gerade): {parity['exactly_one_even']}/{scale['split_count']}  "
                f"OK={parity['parity_constraint_holds']}",
                f"χ² κ-Paar vs 1/16: {scale['chi2_joint_kappa_vs_16']:.3f}  "
                f"z_shuffle={scale['z_score_chi2_vs_shuffle']:.2f}  "
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
        lines.append("Top-5 volle Γ-Signaturen:")
        for row in scale["top_deviations_full_gamma"][:5]:
            lines.append(f"  {row['gamma']}: μ={row['mu_X']:.5f}")
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
        description="Gauß–EABC-Spaltungstest (volle Γ-Signatur)"
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
