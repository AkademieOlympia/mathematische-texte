#!/usr/bin/env python3
"""
Hurwitz–EABC-Orbit-Experiment (volle Γ-Signatur auf O_p).

Kanonsiche Theorie: collatz_eabc_hurwitz_spaltung.md

Für jede Primzahl p ≤ X:
  O_p = { q ∈ H_H : N(q) = p }  (Ganzzahl- und Halbganzzahl-Hurwitz-Koordinaten)
  q = (a, b, c, e)  — vierte Komponente e (User-Nomenklatur)
  strip_smooth auf |Koordinate| (bei Halbganzzahlen: ungerader Zähler)
  Γ(q) = (κ(a'), κ(b'), κ(c'), κ(e'))

Pro Orbit: Größe, Γ-Verteilung, mittlere Chiralität χ = #(E∪C) − #(A∪B),
Kanal-Korrelationen (a',b') vs (c',e'), Unabhängigkeitstest vs. Randprodukt.

Ausführung:
    python3 collatz_eabc_hurwitz_orbit_test.py
    python3 collatz_eabc_hurwitz_orbit_test.py --max-p 10000
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from math import isqrt
from pathlib import Path
from typing import Any

from collatz_eabc_gauss_spaltung_test import (
    CHI2_CRIT_15DF_005,
    EABC_LABELS,
    UNIFORM_MARGINAL,
    _chi_square_uniform,
    _chi_square_uniform_pairs,
    _sieve_primes,
    gamma_label,
    kappa_glatt,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_hurwitz_orbit.json"

Gamma4 = tuple[str, str, str, str]
ZERO_LEG = "0"
CHI2_CRIT_255DF_005 = 310.0  # 255 df for 4^4−1 nontrivial cells with zeros folded
CHI2_CRIT_15DF_PAIR = CHI2_CRIT_15DF_005


def hurwitz_norm(a: int, b: int, c: int, e: int, denom: int = 1) -> int:
    """N(a,b,c,e) mit Koordinaten in Z/denom."""
    if denom == 1:
        return a * a + b * b + c * c + e * e
    return (a * a + b * b + c * c + e * e) // (denom * denom)


def _is_hurwitz_coords(a: int, b: int, c: int, e: int, denom: int) -> bool:
    if denom == 1:
        return True
    # Halbganzzahlig: (m_1+...+m_4)/2 ∈ Z ⇔ Summe der Zähler gerade
    return (a + b + c + e) % 2 == 0


def enum_integer_solutions(p: int) -> set[tuple[int, int, int, int]]:
    sols: set[tuple[int, int, int, int]] = set()
    if p < 1:
        return sols
    bound = isqrt(p)
    for a in range(-bound, bound + 1):
        a2 = a * a
        if a2 > p:
            continue
        for b in range(-bound, bound + 1):
            ab2 = a2 + b * b
            if ab2 > p:
                continue
            for c in range(-bound, bound + 1):
                abc2 = ab2 + c * c
                if abc2 > p:
                    continue
                rem = p - abc2
                d = isqrt(rem)
                if d * d != rem:
                    continue
                if d == 0:
                    sols.add((a, b, c, 0))
                else:
                    for e in (d, -d):
                        sols.add((a, b, c, e))
    return sols


def enum_half_integer_solutions(p: int) -> set[tuple[int, int, int, int]]:
    """Koordinaten als ungerade Zähler (m_a,m_b,m_c,m_e) mit ∑ m_i² = 4p."""
    sols: set[tuple[int, int, int, int]] = set()
    target = 4 * p
    bound = isqrt(target)
    odd_vals = [m for m in range(-bound, bound + 1) if m % 2 != 0]
    for ma in odd_vals:
        s1 = ma * ma
        if s1 > target:
            continue
        for mb in odd_vals:
            s2 = s1 + mb * mb
            if s2 > target:
                continue
            for mc in odd_vals:
                s3 = s2 + mc * mc
                if s3 > target:
                    continue
                rem = target - s3
                me = isqrt(rem)
                if me * me != rem or me % 2 == 0:
                    continue
                for sign in (1, -1) if me != 0 else (1,):
                    me_signed = me * sign
                    if not _is_hurwitz_coords(ma, mb, mc, me_signed, 2):
                        continue
                    sols.add((ma, mb, mc, me_signed))
    return sols


@lru_cache(maxsize=None)
def hurwitz_orbit_elements(p: int) -> tuple[tuple[int, int, int, int, int], ...]:
    """
    Alle q ∈ H_H mit N(q)=p.

    Rückgabe: (a, b, c, e, denom) mit denom ∈ {1, 2}.
    Ganzzahlig: denom=1. Halbganzzahlig: Koordinaten = Zähler/2, denom=2.
    """
    out: list[tuple[int, int, int, int, int]] = []
    for a, b, c, e in sorted(enum_integer_solutions(p)):
        out.append((a, b, c, e, 1))
    for ma, mb, mc, me in sorted(enum_half_integer_solutions(p)):
        out.append((ma, mb, mc, me, 2))
    return tuple(out)


def kappa_leg(coord: int) -> str:
    """κ auf einer Koordinate; 0 → ZERO_LEG (kein EABC-Kern)."""
    if coord == 0:
        return ZERO_LEG
    _alpha, _beta, core, ka = kappa_glatt(abs(coord))
    return ka


def gamma4_label(g: Gamma4) -> str:
    return f"({g[0]},{g[1]},{g[2]},{g[3]})"


def chirality_score(gamma: Gamma4) -> int:
    """χ = #(E∪C) − #(A∪B) über definierte (nicht-Null) Legs."""
    ec = ac = 0
    for k in gamma:
        if k == ZERO_LEG:
            continue
        if k in ("E", "C"):
            ec += 1
        elif k in ("A", "B"):
            ac += 1
    return ec - ac


def mu_p_distribution(gamma_counts: Counter[Gamma4], n: int) -> dict[str, float]:
    """μ_p(γ) = Anteil der Orbit-Punkte mit Signatur γ."""
    if n == 0:
        return {}
    return {gamma4_label(g): c / n for g, c in gamma_counts.items()}


def shannon_entropy(mu: dict[str, float]) -> float:
    """H(p) = −∑_γ μ_p(γ) log μ_p(γ) (natürlicher Logarithmus)."""
    return -sum(p * math.log(p) for p in mu.values() if p > 0)


def entropy_from_counts(gamma_counts: Counter[Gamma4], n: int) -> float:
    if n == 0:
        return 0.0
    return shannon_entropy(mu_p_distribution(gamma_counts, n))


def _H_by_p_mod(mu_rows: list[dict[str, Any]], modulus: int = 12) -> dict[str, Any]:
    """Mittelwert und Streuung von H(p) nach p mod modulus."""
    buckets: dict[int, list[float]] = defaultdict(list)
    for row in mu_rows:
        buckets[row["p"] % modulus].append(row["H_p"])
    by_residue = {
        str(r): {
            "count": len(vals),
            "mean_H": statistics.mean(vals) if vals else 0.0,
            "std_H": statistics.pstdev(vals) if len(vals) > 1 else 0.0,
        }
        for r, vals in sorted(buckets.items())
    }
    means = [v["mean_H"] for v in by_residue.values() if v["count"] > 0]
    spread = max(means) - min(means) if len(means) > 1 else 0.0
    return {
        "modulus": modulus,
        "by_residue": by_residue,
        "mean_H_spread_across_residues": spread,
        "residue_bias_detected": spread > 0.15,
    }


def _H_distribution_stats(H_values: list[float]) -> dict[str, float]:
    if not H_values:
        return {"mean": 0.0, "std": 0.0, "min": 0.0, "max": 0.0, "median": 0.0}
    return {
        "mean": statistics.mean(H_values),
        "std": statistics.pstdev(H_values) if len(H_values) > 1 else 0.0,
        "min": min(H_values),
        "max": max(H_values),
        "median": statistics.median(H_values),
    }


def ab_ce_channels(gamma: Gamma4) -> tuple[tuple[str, str], tuple[str, str]]:
    return (gamma[0], gamma[1]), (gamma[2], gamma[3])


@dataclass(frozen=True, slots=True)
class OrbitPoint:
    a: int
    b: int
    c: int
    e: int
    denom: int
    gamma: Gamma4
    chi: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "coords": [self.a, self.b, self.c, self.e],
            "denom": self.denom,
            "gamma": list(self.gamma),
            "chi": self.chi,
        }


def classify_orbit_point(a: int, b: int, c: int, e: int, denom: int) -> OrbitPoint:
    gamma: Gamma4 = tuple(kappa_leg(x) for x in (a, b, c, e))  # type: ignore[assignment]
    return OrbitPoint(
        a=a, b=b, c=c, e=e, denom=denom, gamma=gamma, chi=chirality_score(gamma)
    )


def _channel_correlation(points: list[OrbitPoint]) -> dict[str, Any]:
    """Korrelation der Kanalpaare (κ_a,κ_b) vs (κ_c,κ_e) auf O_p."""
    n = len(points)
    if n == 0:
        return {"n": 0}
    ab_counts: Counter[tuple[str, str]] = Counter()
    ce_counts: Counter[tuple[str, str]] = Counter()
    joint_ab_ce: Counter[tuple[tuple[str, str], tuple[str, str]]] = Counter()
    for pt in points:
        ab, ce = ab_ce_channels(pt.gamma)
        ab_counts[ab] += 1
        ce_counts[ce] += 1
        joint_ab_ce[(ab, ce)] += 1

    chi2_indep = 0.0
    for (ab, ce), c in joint_ab_ce.items():
        expected = ab_counts[ab] * ce_counts[ce] / n
        if expected > 0:
            chi2_indep += (c - expected) ** 2 / expected

    return {
        "n": n,
        "chi2_ab_ce_independence": chi2_indep,
        "top_ab_pairs": [
            {"pair": list(k), "count": v, "frac": v / n}
            for k, v in ab_counts.most_common(5)
        ],
        "top_ce_pairs": [
            {"pair": list(k), "count": v, "frac": v / n}
            for k, v in ce_counts.most_common(5)
        ],
    }


def _independence_vs_marginals(gamma_counts: Counter[Gamma4], n: int) -> dict[str, Any]:
    """Vergleich μ(Γ) mit Produkt der vier Bein-Marginalen."""
    if n == 0:
        return {"chi2_vs_product_marginals": 0.0, "n": 0}
    marginals: list[Counter[str]] = [Counter() for _ in range(4)]
    for g, c in gamma_counts.items():
        for i, leg in enumerate(g):
            marginals[i][leg] += c

    chi2 = 0.0
    cells = 0
    for g, obs in gamma_counts.items():
        expected = n
        for i, leg in enumerate(g):
            expected *= marginals[i][leg] / n
        if expected > 1e-12:
            chi2 += (obs - expected) ** 2 / expected
            cells += 1

    return {
        "n": n,
        "distinct_gamma": len(gamma_counts),
        "chi2_vs_product_marginals": chi2,
        "cells_used": cells,
        "marginal_legs": [
            {leg: marginals[i][leg] / n for leg in sorted(marginals[i])}
            for i in range(4)
        ],
    }


def _report_from_points(p: int, points: list[OrbitPoint]) -> dict[str, Any]:
    n = len(points)
    gamma_counts: Counter[Gamma4] = Counter(pt.gamma for pt in points)
    chi_values = [pt.chi for pt in points]
    mean_chi = sum(chi_values) / n if n else 0.0
    mu_p = mu_p_distribution(gamma_counts, n)
    H_p = entropy_from_counts(gamma_counts, n)

    int_count = sum(1 for pt in points if pt.denom == 1)
    half_count = n - int_count

    return {
        "p": p,
        "p_mod_12": p % 12,
        "orbit_size": n,
        "integer_coords": int_count,
        "half_integer_coords": half_count,
        "distinct_gamma": len(gamma_counts),
        "H_p": round(H_p, 6),
        "mean_chi": mean_chi,
        "mu_p_top": [
            {"gamma": g, "prob": round(prob, 6)}
            for g, prob in sorted(mu_p.items(), key=lambda kv: -kv[1])[:6]
        ],
        "chi_histogram": dict(Counter(chi_values)),
        "gamma_counts_top": [
            {"gamma": gamma4_label(g), "count": c, "frac": c / n}
            for g, c in gamma_counts.most_common(8)
        ],
        "channel_correlation": _channel_correlation(points),
        "independence": _independence_vs_marginals(gamma_counts, n),
        "sample_points": [pt.to_dict() for pt in points[:4]],
    }


def prime_orbit_report(p: int) -> dict[str, Any]:
    points = [classify_orbit_point(*coords) for coords in hurwitz_orbit_elements(p)]
    return _report_from_points(p, points)


def _sample_prime_brief_rows(
    rows: list[dict[str, Any]], max_rows: int = 40
) -> list[dict[str, Any]]:
    if len(rows) <= max_rows:
        return rows
    head = rows[:15]
    tail = rows[-5:]
    step = max(1, (len(rows) - 20) // (max_rows - 20))
    mid = [rows[i] for i in range(15, len(rows) - 5, step)]
    return head + mid[: max_rows - len(head) - len(tail)] + tail


def aggregate_orbit_report(max_p: int, sample_prime_detail: int = 12) -> dict[str, Any]:
    primes = [p for p in _sieve_primes(max_p) if p >= 2]
    all_gamma: Counter[Gamma4] = Counter()
    all_ab: Counter[tuple[str, str]] = Counter()
    all_ce: Counter[tuple[str, str]] = Counter()
    orbit_sizes: list[int] = []
    mean_chi_by_prime: list[float] = []
    H_by_prime: list[float] = []
    mu_rows: list[dict[str, Any]] = []
    chi2_indep_sum = 0.0
    chi2_indep_count = 0
    total_points = 0

    per_prime_brief: list[dict[str, Any]] = []
    per_prime_detail: list[dict[str, Any]] = []

    for p in primes:
        raw = hurwitz_orbit_elements(p)
        points = [classify_orbit_point(*coords) for coords in raw]
        n = len(points)
        orbit_sizes.append(n)
        total_points += n
        if n == 0:
            continue

        gamma_counts: Counter[Gamma4] = Counter(pt.gamma for pt in points)
        all_gamma.update(gamma_counts)
        for pt in points:
            ab, ce = ab_ce_channels(pt.gamma)
            all_ab[ab] += 1
            all_ce[ce] += 1

        chis = [pt.chi for pt in points]
        mean_chi = sum(chis) / n
        mean_chi_by_prime.append(mean_chi)
        H_p = entropy_from_counts(gamma_counts, n)
        H_by_prime.append(H_p)

        indep = _independence_vs_marginals(gamma_counts, n)
        chi2_indep_sum += indep["chi2_vs_product_marginals"]
        chi2_indep_count += 1

        brief = {
            "p": p,
            "p_mod_12": p % 12,
            "orbit_size": n,
            "distinct_gamma": len(gamma_counts),
            "H_p": round(H_p, 4),
            "mean_chi": round(mean_chi, 4),
        }
        per_prime_brief.append(brief)
        mu_rows.append({"p": p, "H_p": H_p, "mean_chi": mean_chi, "orbit_size": n})
        if len(per_prime_detail) < sample_prime_detail:
            per_prime_detail.append(_report_from_points(p, points))

    n_pts = total_points
    mu_gamma = (
        {gamma4_label(g): c / n_pts for g, c in all_gamma.items()} if n_pts else {}
    )
    mu_ab = {gamma_label(*k): v / n_pts for k, v in all_ab.items()} if n_pts else {}
    mu_ce = {gamma_label(*k): v / n_pts for k, v in all_ce.items()} if n_pts else {}

    marginals: list[Counter[str]] = [Counter() for _ in range(4)]
    for g, c in all_gamma.items():
        for i, leg in enumerate(g):
            marginals[i][leg] += c
    marginal_tables = [
        {leg: marginals[i][leg] / n_pts for leg in sorted(marginals[i])}
        for i in range(4)
    ]

    agg_indep = _independence_vs_marginals(all_gamma, n_pts)
    chi2_ab = _chi_square_uniform_pairs(all_ab, n_pts)
    chi2_ce = _chi_square_uniform_pairs(all_ce, n_pts)

    eabc_marginals = [Counter() for _ in range(4)]
    for g, c in all_gamma.items():
        for i, leg in enumerate(g):
            if leg in EABC_LABELS:
                eabc_marginals[i][leg] += c
    chi2_legs = []
    for i in range(4):
        n_leg = sum(eabc_marginals[i].values())
        if n_leg:
            chi2_legs.append(
                _chi_square_uniform(
                    eabc_marginals[i], EABC_LABELS, n_leg, UNIFORM_MARGINAL
                )
            )
        else:
            chi2_legs.append(0.0)

    orbit_size_stats = {
        "min": min(orbit_sizes) if orbit_sizes else 0,
        "max": max(orbit_sizes) if orbit_sizes else 0,
        "mean": sum(orbit_sizes) / len(orbit_sizes) if orbit_sizes else 0.0,
        "median": sorted(orbit_sizes)[len(orbit_sizes) // 2] if orbit_sizes else 0,
        "sample_sizes": dict(Counter(orbit_sizes).most_common(10)),
    }

    global_mean_chi = (
        sum(mean_chi_by_prime) / len(mean_chi_by_prime) if mean_chi_by_prime else 0.0
    )
    chi_bias = abs(global_mean_chi) > 0.05
    H_stats = _H_distribution_stats(H_by_prime)
    H_mod12 = _H_by_p_mod(mu_rows, modulus=12)

    return {
        "max_p": max_p,
        "prime_count": len(primes),
        "total_orbit_points": n_pts,
        "orbit_size_stats": orbit_size_stats,
        "H_stats": H_stats,
        "H_by_p_mod_12": H_mod12,
        "per_prime_brief": _sample_prime_brief_rows(per_prime_brief),
        "per_prime_brief_count": len(per_prime_brief),
        "per_prime_detail": per_prime_detail,
        "aggregate": {
            "distinct_gamma": len(all_gamma),
            "mu_gamma": dict(sorted(mu_gamma.items(), key=lambda kv: -kv[1])[:32]),
            "mu_ab_channel": mu_ab,
            "mu_ce_channel": mu_ce,
            "marginal_by_leg": marginal_tables,
            "independence_test": agg_indep,
            "chi2_ab_channel_vs_16": chi2_ab,
            "chi2_ce_channel_vs_16": chi2_ce,
            "chi2_marginal_per_leg_vs_quarter": chi2_legs,
            "mean_chi_per_prime_avg": global_mean_chi,
            "orbit_level_chi_bias": chi_bias,
            "mean_orbit_chi2_independence": (
                chi2_indep_sum / chi2_indep_count if chi2_indep_count else 0.0
            ),
        },
        "verdict": _verdict(agg_indep, chi2_ab, chi2_ce, global_mean_chi, chi_bias),
    }


def _verdict(
    indep: dict[str, Any],
    chi2_ab: float,
    chi2_ce: float,
    mean_chi: float,
    chi_bias: bool,
) -> str:
    indep_sig = indep.get("chi2_vs_product_marginals", 0.0) > CHI2_CRIT_255DF_005
    channel_sig = chi2_ab > CHI2_CRIT_15DF_PAIR or chi2_ce > CHI2_CRIT_15DF_PAIR
    if chi_bias and (indep_sig or channel_sig):
        return (
            "orbit_chi_bias: mittlere Chiralität ≠ 0 und/oder Kanal-Kopplung "
            "→ nichttriviale Γ-Struktur auf Hurwitz-Orbits"
        )
    if not chi_bias and not indep_sig:
        return (
            "approximately_uniform: kein robuster orbit-weiter χ-Bias; "
            "Γ nahe Produkt der Bein-Marginalen → Einzigartigkeitsverlust wie erwartet"
        )
    if channel_sig:
        return "channel_coupling: (a',b') vs (c',e') zeigt Abhängigkeit jenseits Marginalen"
    return (
        "inconclusive: schwache Signale — größere Skala oder Einheiten-Quotient nötig"
    )


def multi_scale_report(max_ps: list[int]) -> dict[str, Any]:
    scales = {str(x): aggregate_orbit_report(x) for x in sorted(set(max_ps))}
    largest = scales[str(max(sorted(set(max_ps))))]
    return {
        "meta": {
            "hypothesis_doc": "collatz_eabc_hurwitz_spaltung.md",
            "gaussian_reference": "collatz_eabc_gauss_spaltung_test.py",
            "eisenstein_reference": "collatz_eabc_eisenstein_spaltung_test.py",
            "normabstieg": "collatz_eabc_normabstieg_hypothese.md",
            "euklidische_hebung": "collatz_eabc_euklidische_hebung.md",
            "invarianzprogramm": "collatz_eabc_invarianzprogramm.md",
            "morley_parallel": "collatz_morley_tm_numerik.py (G_M, W_M Konfigurationssensoren)",
            "ring": "H_H (Hurwitz-Maximalordnung), N(q)=a²+b²+c²+e²",
            "research_chain": "p ↦ O_p ↦ μ_p ↦ I(μ_p); I-Kandidaten: H(p), χ̄_p, Kanal-χ²",
            "enumeration": (
                "Ganzzahlige und halbganzzahlige Hurwitz-Darstellungen; "
                "keine kanonische Einzelwahl (vs. Z[i], Z[ω])"
            ),
            "gamma_definition": "Γ(q)=(κ(a'),κ(b'),κ(c'),κ(e')) nach strip_smooth; 0-Leg → '0'",
            "chirality": "χ = #(E∪C) − #(A∪B) über nicht-Null-Legs",
            "scales": sorted(set(max_ps)),
            "largest_scale_total_points": largest["total_orbit_points"],
        },
        "scales": scales,
    }


def format_table(report: dict[str, Any]) -> str:
    lines = [
        "Hurwitz–EABC-Orbit (Γ auf O_p, volle Vierer-Signatur)",
        "=" * 62,
    ]
    for key, scale in report["scales"].items():
        agg = scale["aggregate"]
        oss = scale["orbit_size_stats"]
        lines.extend(
            [
                f"\nX = {key}  (Primzahlen: {scale['prime_count']}, "
                f"Punkte: {scale['total_orbit_points']})",
                f"Orbit-Größe: min={oss['min']} max={oss['max']} "
                f"median={oss['median']} mean={oss['mean']:.1f}",
                f"H(p): mean={scale['H_stats']['mean']:.4f} std={scale['H_stats']['std']:.4f} "
                f"min={scale['H_stats']['min']:.4f} max={scale['H_stats']['max']:.4f}",
                f"H(p) mod-12 spread: {scale['H_by_p_mod_12']['mean_H_spread_across_residues']:.4f}  "
                f"bias: {scale['H_by_p_mod_12']['residue_bias_detected']}",
                f"mittl. χ pro Prim-Orbit (Mittel über p): {agg['mean_chi_per_prime_avg']:.4f}  "
                f"Bias: {agg['orbit_level_chi_bias']}",
                f"χ² AB-Kanal vs 1/16: {agg['chi2_ab_channel_vs_16']:.2f}  "
                f"CE-Kanal: {agg['chi2_ce_channel_vs_16']:.2f}",
                f"χ² Γ vs Produkt-Marginalen: {agg['independence_test']['chi2_vs_product_marginals']:.2f}",
                f"verdict: {scale['verdict']}",
            ]
        )
        lines.append("Beispiel Orbit-Größen (p, |O_p|, H(p)):")
        for row in scale["per_prime_brief"][:8]:
            lines.append(
                f"  p={row['p']}: |O_p|={row['orbit_size']}  "
                f"H={row['H_p']}  mean_χ={row['mean_chi']}"
            )
    return "\n".join(lines)


def run(max_ps: list[int] | None = None, output: Path | None = None) -> dict[str, Any]:
    scales = max_ps or [2_000, 10_000]
    report = multi_scale_report(scales)
    out = output or DEFAULT_OUTPUT
    out.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    report["output_path"] = str(out)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Hurwitz–EABC-Orbit-Experiment")
    parser.add_argument(
        "--max-p",
        type=int,
        nargs="+",
        default=[2_000, 10_000],
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
