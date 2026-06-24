#!/usr/bin/env python3
"""
Hurwitz–EABC-Schalen- und Orbit-Experiment (Γ-Signatur auf Σ_p und O_i).

Kanonsiche Theorie: collatz_eabc_hurwitz_spaltung.md

Drei Ebenen:
  p  ↦  Σ_p = {q ∈ H_H : N(q)=p}   (Normschale, nicht „Orbit“)
       ↦  ⊔_i O_i  mit  O_i = {u q v : u,v ∈ U_H}, |U_H|=24
       ↦  μ_p (Schale) und optional μ^{(i)} (pro U_H-Orbit)

Maßtheorie:
  M_p(γ) = #{q ∈ Σ_p : Γ(q)=γ}
  μ_p(γ) = M_p(γ) / Σ_η M_p(η)

Theorem (ganzzahlige Vier-Quadrat-Darstellungen, Primzahl p):
  r_4(p) = #{a²+b²+c²+e²=p} = 8(p+1)

Invarianten auf μ_p: χ_p, H_p, K_p = Cov_μ(κ(a'),κ(b'),κ(c'),κ(e')).

Ausführung:
    python3 collatz_eabc_hurwitz_orbit_test.py
    python3 collatz_eabc_hurwitz_orbit_test.py --max-p 10000
"""

from __future__ import annotations

import argparse
import json
import math
import random
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
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
LegSig = tuple[int, int, str]
FullGamma4 = tuple[LegSig, LegSig, LegSig, LegSig]
QuatFrac = tuple[Fraction, Fraction, Fraction, Fraction]
ZERO_LEG = "0"
LEG_NAMES = ("a", "b", "c", "e")
CHI2_CRIT_255DF_005 = 310.0
CHI2_CRIT_15DF_PAIR = CHI2_CRIT_15DF_005
KAPPA_NUM = {"E": 0.0, "A": 1.0, "B": 2.0, "C": 3.0, ZERO_LEG: -1.0}


def hurwitz_norm(a: int, b: int, c: int, e: int, denom: int = 1) -> int:
    """N(a,b,c,e) mit Koordinaten in Z/denom."""
    if denom == 1:
        return a * a + b * b + c * c + e * e
    return (a * a + b * b + c * c + e * e) // (denom * denom)


def r4_theorem(p: int) -> int:
    """r_4(p) = 8(p+1) für Primzahl p (Theorem, ganzzahlige Vier-Quadrate)."""
    return 8 * (p + 1)


def _is_hurwitz_coords(a: int, b: int, c: int, e: int, denom: int) -> bool:
    if denom == 1:
        return True
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
def hurwitz_shell_elements(p: int) -> tuple[tuple[int, int, int, int, int], ...]:
    """
    Σ_p = {q ∈ H_H : N(q)=p}.

    Rückgabe: (a, b, c, e, denom) mit denom ∈ {1, 2}.
  """
    out: list[tuple[int, int, int, int, int]] = []
    for a, b, c, e in sorted(enum_integer_solutions(p)):
        out.append((a, b, c, e, 1))
    for ma, mb, mc, me in sorted(enum_half_integer_solutions(p)):
        out.append((ma, mb, mc, me, 2))
    return tuple(out)


# Rückwärtskompatibler Alias
hurwitz_orbit_elements = hurwitz_shell_elements


def coords_to_frac(a: int, b: int, c: int, e: int, denom: int) -> QuatFrac:
    return (
        Fraction(a, denom),
        Fraction(b, denom),
        Fraction(c, denom),
        Fraction(e, denom),
    )


def frac_to_coords(q: QuatFrac) -> tuple[int, int, int, int, int]:
    a, b, c, e = q
    if all(x.denominator == 1 for x in q):
        return (a.numerator, b.numerator, c.numerator, e.numerator, 1)
    return (int(2 * a), int(2 * b), int(2 * c), int(2 * e), 2)


def q_mul(x: QuatFrac, y: QuatFrac) -> QuatFrac:
    a, b, c, e = x
    f, g, h, k = y
    return (
        a * f - b * g - c * h - e * k,
        a * g + b * f + c * k - e * h,
        a * h - b * k + c * f + e * g,
        a * k + b * h - c * g + e * f,
    )


@lru_cache(maxsize=1)
def hurwitz_units_frac() -> tuple[QuatFrac, ...]:
    """U_H: 24 Hurwitz-Einheiten (Norm 1)."""
    units: set[QuatFrac] = set()
    basic = [
        (1, 0, 0, 0),
        (-1, 0, 0, 0),
        (0, 1, 0, 0),
        (0, -1, 0, 0),
        (0, 0, 1, 0),
        (0, 0, -1, 0),
        (0, 0, 0, 1),
        (0, 0, 0, -1),
    ]
    for t in basic:
        units.add(coords_to_frac(*t, 1))
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            for s3 in (-1, 1):
                for s4 in (-1, 1):
                    units.add(
                        (
                            Fraction(s1, 2),
                            Fraction(s2, 2),
                            Fraction(s3, 2),
                            Fraction(s4, 2),
                        )
                    )
    if len(units) != 24:
        raise ValueError(f"Erwartet 24 Einheiten, erhalten {len(units)}")
    return tuple(sorted(units))


def double_orbit_of(q: QuatFrac, shell_set: set[QuatFrac]) -> frozenset[QuatFrac]:
    """O(q) = {u q v : u,v ∈ U_H} ∩ Σ_p."""
    units = hurwitz_units_frac()
    orb: set[QuatFrac] = set()
    for u in units:
        uq = q_mul(u, q)
        for v in units:
            w = q_mul(uq, v)
            if w in shell_set:
                orb.add(w)
    return frozenset(orb)


@lru_cache(maxsize=None)
def uh_orbit_partition(p: int) -> tuple[frozenset[QuatFrac], ...]:
    """Zerlegung Σ_p in U_H-Doppelbahnen O_i = {uqv}."""
    shell = [coords_to_frac(*c) for c in hurwitz_shell_elements(p)]
    shell_set = set(shell)
    units = hurwitz_units_frac()
    seen: set[QuatFrac] = set()
    orbits: list[frozenset[QuatFrac]] = []
    for q in shell:
        if q in seen:
            continue
        orb: set[QuatFrac] = set()
        for u in units:
            uq = q_mul(u, q)
            for v in units:
                w = q_mul(uq, v)
                if w in shell_set:
                    orb.add(w)
        fro = frozenset(orb)
        orbits.append(fro)
        seen |= orb
    return tuple(orbits)


def uh_orbit_count(p: int) -> int:
    return len(uh_orbit_partition(p))


def leg_signature(coord: int) -> LegSig:
    """Γ_leg(x) = (α, β, κ(x')) oder (0, 0, '0') bei x=0."""
    if coord == 0:
        return (0, 0, ZERO_LEG)
    alpha, beta, _core, ka = kappa_glatt(abs(coord))
    return (alpha, beta, ka)


def kappa_leg(coord: int) -> str:
    return leg_signature(coord)[2]


def full_gamma4(a: int, b: int, c: int, e: int) -> FullGamma4:
    return (leg_signature(a), leg_signature(b), leg_signature(c), leg_signature(e))


def full_gamma4_label(sig: FullGamma4) -> str:
    parts: list[str] = []
    for name, leg in zip(LEG_NAMES, sig):
        parts.append(f"{name}:({leg[0]},{leg[1]},{leg[2]})")
    return "{" + ",".join(parts) + "}"


def full_gamma4_compact(sig: FullGamma4) -> tuple[int, int, str, int, int, str, int, int, str, int, int, str]:
    return (
        sig[0][0], sig[0][1], sig[0][2],
        sig[1][0], sig[1][1], sig[1][2],
        sig[2][0], sig[2][1], sig[2][2],
        sig[3][0], sig[3][1], sig[3][2],
    )


def smooth_pattern4(sig: FullGamma4) -> tuple[int, int, int, int, int, int, int, int]:
    """(α_a,β_a,α_b,β_b,α_c,β_c,α_e,β_e) — stratifiziert volle Γ."""
    return tuple(x for leg in sig for x in leg[:2])  # type: ignore[return-value]


def gamma4_label(g: Gamma4) -> str:
    return f"({g[0]},{g[1]},{g[2]},{g[3]})"


def chirality_score(gamma: Gamma4) -> int:
    ec = ac = 0
    for k in gamma:
        if k == ZERO_LEG:
            continue
        if k in ("E", "C"):
            ec += 1
        elif k in ("A", "B"):
            ac += 1
    return ec - ac


def M_p_counts(gamma_counts: Counter[Gamma4]) -> dict[str, int]:
    return {gamma4_label(g): c for g, c in gamma_counts.items()}


def mu_p_distribution(gamma_counts: Counter[Gamma4], n: int) -> dict[str, float]:
    if n == 0:
        return {}
    return {gamma4_label(g): c / n for g, c in gamma_counts.items()}


def shannon_entropy(mu: dict[str, float]) -> float:
    return -sum(p * math.log(p) for p in mu.values() if p > 0)


def entropy_from_counts(gamma_counts: Counter[Gamma4], n: int) -> float:
    if n == 0:
        return 0.0
    return shannon_entropy(mu_p_distribution(gamma_counts, n))


def chi_p_from_counts(gamma_counts: Counter[Gamma4], n: int) -> float:
    """χ_p = Σ_γ χ(γ) μ_p(γ)."""
    if n == 0:
        return 0.0
    return sum(chirality_score(g) * c / n for g, c in gamma_counts.items())


@dataclass(frozen=True, slots=True)
class OrbitPoint:
    a: int
    b: int
    c: int
    e: int
    denom: int
    gamma: Gamma4
    full_gamma: FullGamma4
    chi: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "coords": [self.a, self.b, self.c, self.e],
            "denom": self.denom,
            "gamma": list(self.gamma),
            "full_gamma": [list(leg) for leg in self.full_gamma],
            "full_gamma_compact": list(full_gamma4_compact(self.full_gamma)),
            "chi": self.chi,
        }


def K_p_covariance(points: list[OrbitPoint]) -> dict[str, Any]:
    """K_p = Cov_{μ_p}(κ(a'),κ(b'),κ(c'),κ(e'))."""
    n = len(points)
    if n < 2:
        return {"matrix": [[0.0] * 4 for _ in range(4)], "trace": 0.0, "frobenius": 0.0}
    rows = [[KAPPA_NUM[leg] for leg in pt.gamma] for pt in points]
    means = [sum(r[i] for r in rows) / n for i in range(4)]
    cov = [[0.0] * 4 for _ in range(4)]
    for row in rows:
        for i in range(4):
            for j in range(4):
                cov[i][j] += (row[i] - means[i]) * (row[j] - means[j])
    for i in range(4):
        for j in range(4):
            cov[i][j] /= n
    trace = sum(cov[i][i] for i in range(4))
    frob = math.sqrt(sum(cov[i][j] ** 2 for i in range(4) for j in range(4)))
    return {
        "matrix": [[round(cov[i][j], 6) for j in range(4)] for i in range(4)],
        "trace": round(trace, 6),
        "frobenius": round(frob, 6),
    }


def shuffle_null_shell(
    points: list[OrbitPoint], trials: int = 20, seed: int = 0
) -> dict[str, Any]:
    """Null-Kontrolle: Γ-Labels auf Σ_p zufällig permutieren."""
    rng = random.Random(seed)
    gammas = [pt.gamma for pt in points]
    n = len(gammas)
    if n == 0:
        return {"trials": trials, "H_p_mean": 0.0, "chi_p_mean": 0.0, "K_trace_mean": 0.0}

    H_vals: list[float] = []
    chi_vals: list[float] = []
    K_traces: list[float] = []
    for _ in range(trials):
        perm = gammas[:]
        rng.shuffle(perm)
        counts: Counter[Gamma4] = Counter(perm)
        H_vals.append(entropy_from_counts(counts, n))
        chi_vals.append(chi_p_from_counts(counts, n))
        null_pts = [
            OrbitPoint(pt.a, pt.b, pt.c, pt.e, pt.denom, g, pt.full_gamma, chirality_score(g))
            for pt, g in zip(points, perm)
        ]
        K_traces.append(K_p_covariance(null_pts)["trace"])

    return {
        "trials": trials,
        "H_p_mean": round(statistics.mean(H_vals), 6),
        "H_p_std": round(statistics.pstdev(H_vals) if len(H_vals) > 1 else 0.0, 6),
        "chi_p_mean": round(statistics.mean(chi_vals), 6),
        "chi_p_std": round(statistics.pstdev(chi_vals) if len(chi_vals) > 1 else 0.0, 6),
        "K_trace_mean": round(statistics.mean(K_traces), 6),
        "K_trace_std": round(statistics.pstdev(K_traces) if len(K_traces) > 1 else 0.0, 6),
    }


def _H_by_p_mod(mu_rows: list[dict[str, Any]], modulus: int = 12) -> dict[str, Any]:
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


def _H_trend(H_by_p: list[tuple[int, float]]) -> dict[str, Any]:
    """Linearer Trend H_p über p (für Limit-Frage)."""
    if len(H_by_p) < 3:
        return {"slope": 0.0, "intercept": 0.0, "n": len(H_by_p)}
    ps = [x[0] for x in H_by_p]
    Hs = [x[1] for x in H_by_p]
    n = len(ps)
    mp = sum(ps) / n
    mH = sum(Hs) / n
    num = sum((ps[i] - mp) * (Hs[i] - mH) for i in range(n))
    den = sum((ps[i] - mp) ** 2 for i in range(n))
    slope = num / den if den else 0.0
    intercept = mH - slope * mp
    return {"slope": round(slope, 8), "intercept": round(intercept, 6), "n": n}


def ab_ce_channels(gamma: Gamma4) -> tuple[tuple[str, str], tuple[str, str]]:
    return (gamma[0], gamma[1]), (gamma[2], gamma[3])


def classify_shell_point(a: int, b: int, c: int, e: int, denom: int) -> OrbitPoint:
    fg = full_gamma4(a, b, c, e)
    gamma: Gamma4 = tuple(leg[2] for leg in fg)  # type: ignore[assignment]
    return OrbitPoint(
        a=a, b=b, c=c, e=e, denom=denom, gamma=gamma, full_gamma=fg, chi=chirality_score(gamma)
    )


classify_orbit_point = classify_shell_point


def _channel_correlation(points: list[OrbitPoint]) -> dict[str, Any]:
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


def verify_r4_theorem(p: int) -> dict[str, Any]:
    """Prüft r_4(p)=8(p+1) auf ganzzahligen Darstellungen und |Σ_p|."""
    int_count = len(enum_integer_solutions(p))
    shell_size = len(hurwitz_shell_elements(p))
    r4 = r4_theorem(p)
    return {
        "p": p,
        "r4_theorem": r4,
        "integer_reps": int_count,
        "r4_matches": int_count == r4,
        "shell_size": shell_size,
        "shell_over_r4": round(shell_size / r4, 6) if r4 else 0.0,
    }


def _asymptotic_p_buckets(mu_rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Trend von H_p, χ_p in drei p-Größen-Buckets (klein/mittel/groß)."""
    if not mu_rows:
        return {"buckets": {}, "mu_infinity_hint": "insufficient_data"}
    sorted_rows = sorted(mu_rows, key=lambda r: r["p"])
    n = len(sorted_rows)
    cuts = [n // 3, 2 * n // 3]
    labels = ["small_p", "medium_p", "large_p"]
    slices: list[list[dict[str, Any]]] = [
        sorted_rows[: cuts[0]] if cuts[0] else [],
        sorted_rows[cuts[0] : cuts[1]],
        sorted_rows[cuts[1] :],
    ]
    buckets: dict[str, Any] = {}
    H_means: list[float] = []
    chi_means: list[float] = []
    for label, chunk in zip(labels, slices):
        if not chunk:
            buckets[label] = {"count": 0, "p_range": None, "mean_H_p": 0.0, "mean_chi_p": 0.0}
            continue
        H_vals = [r["H_p"] for r in chunk]
        chi_vals = [r["chi_p"] for r in chunk]
        mH = statistics.mean(H_vals)
        mchi = statistics.mean(chi_vals)
        H_means.append(mH)
        chi_means.append(mchi)
        buckets[label] = {
            "count": len(chunk),
            "p_range": [chunk[0]["p"], chunk[-1]["p"]],
            "mean_H_p": round(mH, 6),
            "mean_chi_p": round(mchi, 6),
            "std_H_p": round(statistics.pstdev(H_vals) if len(H_vals) > 1 else 0.0, 6),
        }
    hint = "stable"
    if len(H_means) >= 2:
        dH = H_means[-1] - H_means[0]
        if abs(dH) > 0.05:
            hint = "H_p_drift_large_to_small" if dH < 0 else "H_p_rises_with_p"
        elif abs(chi_means[-1] - chi_means[0]) > 0.03:
            hint = "chi_p_drift_across_buckets"
        else:
            hint = "consistent_with_mu_infinity"
    return {
        "buckets": buckets,
        "mu_infinity_hint": hint,
        "H_delta_large_minus_small": round(H_means[-1] - H_means[0], 6) if len(H_means) >= 2 else 0.0,
    }


def _orbit_level_reports(p: int, points: list[OrbitPoint]) -> dict[str, Any]:
    orbits = uh_orbit_partition(p)
    orbit_sizes = sorted(len(o) for o in orbits)
    per_orbit: list[dict[str, Any]] = []
    for idx, orb in enumerate(orbits[:6]):
        orb_pts = []
        for q in orb:
            a, b, c, e, denom = frac_to_coords(q)
            orb_pts.append(classify_shell_point(a, b, c, e, denom))
        gc: Counter[Gamma4] = Counter(pt.gamma for pt in orb_pts)
        gfc: Counter[FullGamma4] = Counter(pt.full_gamma for pt in orb_pts)
        n = len(orb_pts)
        per_orbit.append(
            {
                "orbit_index": idx,
                "size": n,
                "H_p_orbit": round(entropy_from_counts(gc, n), 6),
                "chi_p_orbit": round(chi_p_from_counts(gc, n), 6),
                "mu_orbit_top": [
                    {"gamma": gamma4_label(g), "prob": round(c / n, 6)}
                    for g, c in gc.most_common(3)
                ],
                "mu_orbit_full_top": [
                    {"full_gamma": full_gamma4_label(g), "prob": round(c / n, 6)}
                    for g, c in gfc.most_common(2)
                ],
            }
        )
    return {
        "orbit_count": len(orbits),
        "orbit_sizes": orbit_sizes,
        "per_orbit_sample": per_orbit,
    }


def _report_from_points(p: int, points: list[OrbitPoint], include_orbits: bool = True) -> dict[str, Any]:
    n = len(points)
    gamma_counts: Counter[Gamma4] = Counter(pt.gamma for pt in points)
    full_gamma_counts: Counter[FullGamma4] = Counter(pt.full_gamma for pt in points)
    smooth_patterns: Counter[tuple[int, ...]] = Counter(
        smooth_pattern4(pt.full_gamma) for pt in points
    )
    mu_p = mu_p_distribution(gamma_counts, n)
    H_p = entropy_from_counts(gamma_counts, n)
    chi_p = chi_p_from_counts(gamma_counts, n)
    K_p = K_p_covariance(points)
    r4_check = verify_r4_theorem(p)

    int_count = sum(1 for pt in points if pt.denom == 1)
    half_count = n - int_count

    report: dict[str, Any] = {
        "p": p,
        "p_mod_12": p % 12,
        "shell_size": n,
        "orbit_size": n,  # Rückwärtskompatibilität
        "integer_coords": int_count,
        "half_integer_coords": half_count,
        "r4_check": r4_check,
        "distinct_gamma": len(gamma_counts),
        "distinct_full_gamma": len(full_gamma_counts),
        "distinct_smooth_patterns": len(smooth_patterns),
        "M_p": M_p_counts(gamma_counts),
        "M_p_full_gamma_top": [
            {"full_gamma": full_gamma4_label(g), "count": c}
            for g, c in full_gamma_counts.most_common(6)
        ],
        "H_p": round(H_p, 6),
        "chi_p": round(chi_p, 6),
        "mean_chi": chi_p,
        "K_p": K_p,
        "mu_p_top": [
            {"gamma": g, "prob": round(prob, 6)}
            for g, prob in sorted(mu_p.items(), key=lambda kv: -kv[1])[:6]
        ],
        "chi_histogram": dict(Counter(pt.chi for pt in points)),
        "gamma_counts_top": [
            {"gamma": gamma4_label(g), "count": c, "frac": c / n}
            for g, c in gamma_counts.most_common(8)
        ],
        "channel_correlation": _channel_correlation(points),
        "independence": _independence_vs_marginals(gamma_counts, n),
        "shuffle_null": shuffle_null_shell(points, trials=20, seed=p),
        "sample_points": [pt.to_dict() for pt in points[:4]],
    }
    if include_orbits:
        report["uh_orbits"] = _orbit_level_reports(p, points)
    return report


def prime_shell_report(p: int) -> dict[str, Any]:
    points = [classify_shell_point(*coords) for coords in hurwitz_shell_elements(p)]
    return _report_from_points(p, points)


prime_orbit_report = prime_shell_report


def _sample_prime_brief_rows(rows: list[dict[str, Any]], max_rows: int = 40) -> list[dict[str, Any]]:
    if len(rows) <= max_rows:
        return rows
    head = rows[:15]
    tail = rows[-5:]
    step = max(1, (len(rows) - 20) // (max_rows - 20))
    mid = [rows[i] for i in range(15, len(rows) - 5, step)]
    return head + mid[: max_rows - len(head) - len(tail)] + tail


def aggregate_shell_report(max_p: int, sample_prime_detail: int = 12) -> dict[str, Any]:
    primes = [p for p in _sieve_primes(max_p) if p >= 2]
    all_gamma: Counter[Gamma4] = Counter()
    all_ab: Counter[tuple[str, str]] = Counter()
    all_ce: Counter[tuple[str, str]] = Counter()
    shell_sizes: list[int] = []
    chi_p_by_prime: list[float] = []
    H_by_prime: list[float] = []
    H_by_p_pairs: list[tuple[int, float]] = []
    mu_rows: list[dict[str, Any]] = []
    r4_checks: list[dict[str, Any]] = []
    orbit_counts: list[int] = []
    chi2_indep_sum = 0.0
    chi2_indep_count = 0
    total_points = 0

    per_prime_brief: list[dict[str, Any]] = []
    per_prime_detail: list[dict[str, Any]] = []

    for p in primes:
        raw = hurwitz_shell_elements(p)
        points = [classify_shell_point(*coords) for coords in raw]
        n = len(points)
        shell_sizes.append(n)
        total_points += n
        r4_checks.append(verify_r4_theorem(p))
        if n == 0:
            continue

        orbits = uh_orbit_partition(p)
        orbit_counts.append(len(orbits))

        gamma_counts: Counter[Gamma4] = Counter(pt.gamma for pt in points)
        all_gamma.update(gamma_counts)
        for pt in points:
            ab, ce = ab_ce_channels(pt.gamma)
            all_ab[ab] += 1
            all_ce[ce] += 1

        chi_p = chi_p_from_counts(gamma_counts, n)
        chi_p_by_prime.append(chi_p)
        H_p = entropy_from_counts(gamma_counts, n)
        H_by_prime.append(H_p)
        H_by_p_pairs.append((p, H_p))

        indep = _independence_vs_marginals(gamma_counts, n)
        chi2_indep_sum += indep["chi2_vs_product_marginals"]
        chi2_indep_count += 1

        brief = {
            "p": p,
            "p_mod_12": p % 12,
            "shell_size": n,
            "orbit_size": n,
            "orbit_count": len(orbits),
            "r4_matches": r4_checks[-1]["r4_matches"],
            "distinct_gamma": len(gamma_counts),
            "H_p": round(H_p, 4),
            "chi_p": round(chi_p, 4),
            "mean_chi": round(chi_p, 4),
        }
        per_prime_brief.append(brief)
        mu_rows.append(
            {
                "p": p,
                "H_p": H_p,
                "chi_p": chi_p,
                "mean_chi": chi_p,
                "shell_size": n,
                "orbit_count": len(orbits),
            }
        )
        if len(per_prime_detail) < sample_prime_detail:
            per_prime_detail.append(_report_from_points(p, points))

    n_pts = total_points
    mu_gamma = {gamma4_label(g): c / n_pts for g, c in all_gamma.items()} if n_pts else {}
    mu_ab = {gamma_label(*k): v / n_pts for k, v in all_ab.items()} if n_pts else {}
    mu_ce = {gamma_label(*k): v / n_pts for k, v in all_ce.items()} if n_pts else {}

    marginals: list[Counter[str]] = [Counter() for _ in range(4)]
    for g, c in all_gamma.items():
        for i, leg in enumerate(g):
            marginals[i][leg] += c
    marginal_tables = [
        {leg: marginals[i][leg] / n_pts for leg in sorted(marginals[i])} for i in range(4)
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
                _chi_square_uniform(eabc_marginals[i], EABC_LABELS, n_leg, UNIFORM_MARGINAL)
            )
        else:
            chi2_legs.append(0.0)

    shell_size_stats = {
        "min": min(shell_sizes) if shell_sizes else 0,
        "max": max(shell_sizes) if shell_sizes else 0,
        "mean": sum(shell_sizes) / len(shell_sizes) if shell_sizes else 0.0,
        "median": sorted(shell_sizes)[len(shell_sizes) // 2] if shell_sizes else 0,
        "sample_sizes": dict(Counter(shell_sizes).most_common(10)),
    }

    global_chi_p = sum(chi_p_by_prime) / len(chi_p_by_prime) if chi_p_by_prime else 0.0
    chi_bias = abs(global_chi_p) > 0.05
    H_stats = _H_distribution_stats(H_by_prime)
    H_mod12 = _H_by_p_mod(mu_rows, modulus=12)
    H_trend = _H_trend(H_by_p_pairs)
    asymptotic = _asymptotic_p_buckets(mu_rows)

    r4_all_match = all(c["r4_matches"] for c in r4_checks)

    return {
        "max_p": max_p,
        "prime_count": len(primes),
        "total_shell_points": n_pts,
        "total_orbit_points": n_pts,
        "shell_size_stats": shell_size_stats,
        "orbit_size_stats": shell_size_stats,
        "r4_verification": {
            "theorem": "r_4(p)=8(p+1) für ganzzahlige Vier-Quadrate",
            "all_primes_match": r4_all_match,
            "sample_checks": r4_checks[:12] + r4_checks[-3:] if len(r4_checks) > 15 else r4_checks,
        },
        "orbit_count_stats": {
            "min": min(orbit_counts) if orbit_counts else 0,
            "max": max(orbit_counts) if orbit_counts else 0,
            "mean": sum(orbit_counts) / len(orbit_counts) if orbit_counts else 0.0,
        },
        "H_stats": H_stats,
        "H_trend": H_trend,
        "H_by_p_mod_12": H_mod12,
        "asymptotic_buckets": asymptotic,
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
            "chi_p_per_prime_avg": global_chi_p,
            "mean_chi_per_prime_avg": global_chi_p,
            "shell_level_chi_bias": chi_bias,
            "orbit_level_chi_bias": chi_bias,
            "mean_shell_chi2_independence": chi2_indep_sum / chi2_indep_count if chi2_indep_count else 0.0,
        },
        "verdict": _verdict(agg_indep, chi2_ab, chi2_ce, global_chi_p, chi_bias),
    }


aggregate_orbit_report = aggregate_shell_report


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
            "shell_chi_bias: χ_p ≠ 0 und/oder Kanal-Kopplung "
            "→ nichttriviale Γ-Struktur auf Hurwitz-Schale Σ_p"
        )
    if not chi_bias and not indep_sig:
        return (
            "approximately_uniform: kein robuster Schalen-χ-Bias; "
            "μ_p nahe Produkt der Bein-Marginalen → Vielfalt wie erwartet"
        )
    if channel_sig:
        return "channel_coupling: (a',b') vs (c',e') zeigt Abhängigkeit jenseits Marginalen"
    return "inconclusive: schwache Signale — größere Skala oder feinere Orbit-Statistik nötig"


def multi_scale_report(max_ps: list[int]) -> dict[str, Any]:
    scales = {str(x): aggregate_shell_report(x) for x in sorted(set(max_ps))}
    largest = scales[str(max(sorted(set(max_ps))))]
    return {
        "meta": {
            "hypothesis_doc": "collatz_eabc_hurwitz_spaltung.md",
            "extended_mass_hypothesis": "collatz_eabc_quaternion_mass_hypothese.md",
            "gaussian_reference": "collatz_eabc_gauss_spaltung_test.py",
            "eisenstein_reference": "collatz_eabc_eisenstein_spaltung_test.py",
            "normabstieg": "collatz_eabc_normabstieg_hypothese.md",
            "euklidische_hebung": "collatz_eabc_euklidische_hebung.md",
            "invarianzprogramm": "collatz_eabc_invarianzprogramm.md",
            "morley_parallel": "collatz_morley_tm_numerik.py (G_M, W_M — nur heuristische Analogie)",
            "ring": "H_H (Hurwitz-Maximalordnung), N(q)=a²+b²+c²+e²",
            "research_chain": "p ↦ Σ_p ↦ μ_p ↦ (χ_p, H_p, K_p, …); Σ_p = ⊔_i O_i unter U_H",
            "shell_vs_orbit": "Σ_p = Normschale; O(q) = {uqv : u,v∈U_H, |U_H|=24}",
            "measure": "M_p(γ)=#{q∈Σ_p:Γ(q)=γ}; μ_p(γ)=M_p(γ)/Σ_η M_p(η)",
            "r4_theorem": "r_4(p)=8(p+1) für Primzahl p (ganzzahlige Vier-Quadrate)",
            "enumeration": (
                "Ganzzahlige und halbganzzahlige Hurwitz-Darstellungen; "
                "keine kanonische Einzelwahl (vs. Z[i], Z[ω])"
            ),
            "gamma_definition": (
                "Γ(q)=((α_a,β_a,κ(a')),(α_b,β_b,κ(b')),(α_c,β_c,κ(c')),(α_e,β_e,κ(e'))) "
                "nach strip_smooth; 0-Leg → (0,0,'0')"
            ),
            "chirality": "χ_p = Σ_γ χ(γ) μ_p(γ); χ(q)=#(E∪C)−#(A∪B)",
            "scales": sorted(set(max_ps)),
            "largest_scale_total_points": largest["total_shell_points"],
        },
        "scales": scales,
    }


def format_table(report: dict[str, Any]) -> str:
    lines = [
        "Hurwitz–EABC (Γ auf Σ_p, U_H-Orbit-Zerlegung)",
        "=" * 62,
    ]
    for key, scale in report["scales"].items():
        agg = scale["aggregate"]
        sss = scale["shell_size_stats"]
        r4 = scale["r4_verification"]
        lines.extend(
            [
                f"\nX = {key}  (Primzahlen: {scale['prime_count']}, "
                f"|Σ|-Punkte: {scale['total_shell_points']})",
                f"r₄-Theorem: alle Primzahlen r₄=8(p+1)? {r4['all_primes_match']}",
                f"|Σ_p|: min={sss['min']} max={sss['max']} "
                f"median={sss['median']} mean={sss['mean']:.1f}",
                f"U_H-Orbits/ p: mean={scale['orbit_count_stats']['mean']:.1f}",
                f"H_p: mean={scale['H_stats']['mean']:.4f} std={scale['H_stats']['std']:.4f} "
                f"trend slope={scale['H_trend']['slope']:.6f}  "
                f"μ_∞ hint={scale['asymptotic_buckets']['mu_infinity_hint']}",
                f"χ_p (Mittel über p): {agg['chi_p_per_prime_avg']:.4f}  Bias: {agg['shell_level_chi_bias']}",
                f"verdict: {scale['verdict']}",
            ]
        )
        lines.append("Beispiel (p, |Σ_p|, r₄-ok, #Orbits, H_p, χ_p):")
        for row in scale["per_prime_brief"][:8]:
            lines.append(
                f"  p={row['p']}: |Σ_p|={row['shell_size']}  r₄={row['r4_matches']}  "
                f"#O={row['orbit_count']}  H={row['H_p']}  χ={row['chi_p']}"
            )
    return "\n".join(lines)


def run(max_ps: list[int] | None = None, output: Path | None = None) -> dict[str, Any]:
    scales = max_ps or [2_000, 10_000]
    report = multi_scale_report(scales)
    out = output or DEFAULT_OUTPUT
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report["output_path"] = str(out)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Hurwitz–EABC-Schalen- und Orbit-Experiment")
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
