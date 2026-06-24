#!/usr/bin/env python3
"""
EABC-Zerlegungsregimen: Z_fact, Z_EABC, Z_regime und ΔZ-Sprünge an Primzahlen.

Kanonsiche Theorie: collatz_eabc_zerlegungsregimen.md
Infrastruktur: collatz_eabc_hurwitz_orbit_test.py (Σ_n, μ_n), collatz_eabc_shell_defekt_test.py (I(μ_n))

Operational definitions (Hurwitz-Schalen, n ≤ max_n):
  Z_fact(n)   — Anzahl ungeordneter Faktorisierungen n=ab mit a,b>1 (klassische Baseline)
  Z_EABC(n)   — Anzahl verschiedener EABC-Zerlegungskanäle (μ_a, μ_b)-Signaturen über alle Faktorisierungen
  Z_regime(n) — 1 wenn μ_n nicht als μ_a⊗μ_b aus irgendeiner Faktorisierung erklärbar; 0 sonst; Prim: 1
  Z(n)        — primär Z_EABC(n) (Decompositionszahl); ΔZ(n)=Z(n)-Z(n-1)

Ausführung:
    python3 collatz_eabc_Z_decomposition_test.py
    python3 collatz_eabc_Z_decomposition_test.py --max-n 200 --output collatz_eabc_Z_decomposition.json
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from collections import Counter
from dataclasses import dataclass
from math import isqrt
from pathlib import Path
from typing import Any

from collatz_eabc_hurwitz_orbit_test import (
    Gamma4,
    K_p_covariance,
    chi_p_from_counts,
    classify_shell_point,
    entropy_from_counts,
    hurwitz_shell_elements,
)
from collatz_eabc_shell_defekt_test import is_prime, omega_distinct, tau_divisor_count

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_Z_decomposition.json"

# Heuristik-Schwelle: L2-Abstand I(μ_n) vs. additiv vorhergesagtes I aus Faktorisierung
REGIME_TOL = 0.75


def unordered_factorizations(n: int, min_factor: int = 2) -> list[tuple[int, int]]:
    """Ungeordnete Paare (a,b) mit a≤b, a*b=n, a,b≥min_factor."""
    if n < min_factor * min_factor:
        return []
    out: list[tuple[int, int]] = []
    for a in range(min_factor, isqrt(n) + 1):
        if n % a == 0:
            b = n // a
            if b >= min_factor:
                out.append((a, b))
    return out


def Z_fact(n: int) -> int:
    return len(unordered_factorizations(n, min_factor=2))


@dataclass(frozen=True, slots=True)
class ShellMeasure:
    n: int
    shell_size: int
    H: float
    chi: float
    K_trace: float
    distinct_gamma: int
    gamma_counts: Counter[Gamma4]

    def I_vector(self) -> tuple[float, float, float]:
        return (self.H, self.chi, self.K_trace)

    def mu_signature(self) -> tuple[float, float, float, int]:
        """Kompakte Signatur für Zerlegungskanäle."""
        return (
            round(self.H, 4),
            round(self.chi, 4),
            round(self.K_trace, 4),
            self.distinct_gamma,
        )

    def marginals(self) -> tuple[tuple[str, ...], ...]:
        """Randverteilungen pro Bein (κ-Label)."""
        legs: list[Counter[str]] = [Counter() for _ in range(4)]
        for g, c in self.gamma_counts.items():
            for i, leg in enumerate(g):
                legs[i][leg] += c
        return tuple(
            tuple(leg for leg, _ in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0])))
            for cnt in legs
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "n": self.n,
            "shell_size": self.shell_size,
            "H": round(self.H, 6),
            "chi": round(self.chi, 6),
            "K_trace": round(self.K_trace, 6),
            "distinct_gamma": self.distinct_gamma,
            "mu_signature": list(self.mu_signature()),
        }


def compute_shell_measure(n: int) -> ShellMeasure | None:
    raw = hurwitz_shell_elements(n)
    if not raw:
        return None
    points = [classify_shell_point(*coords) for coords in raw]
    gamma_counts: Counter[Gamma4] = Counter(pt.gamma for pt in points)
    size = len(points)
    return ShellMeasure(
        n=n,
        shell_size=size,
        H=entropy_from_counts(gamma_counts, size),
        chi=chi_p_from_counts(gamma_counts, size),
        K_trace=K_p_covariance(points)["trace"],
        distinct_gamma=len(gamma_counts),
        gamma_counts=gamma_counts,
    )


def channel_signature(
    sig_a: tuple[float, float, float, int], sig_b: tuple[float, float, float, int]
) -> tuple:
    """Ungeordnetes Kanalpaar."""
    return tuple(sorted((sig_a, sig_b)))


def additive_I_predict(
    mu_a: ShellMeasure, mu_b: ShellMeasure
) -> tuple[float, float, float]:
    """Schematische V_4-Vorhersage: I(μ_ab) ≈ I(μ_a) + I(μ_b) (Heuristik, dokumentiert)."""
    Ia, Ib = mu_a.I_vector(), mu_b.I_vector()
    return (Ia[0] + Ib[0], Ia[1] + Ib[1], Ia[2] + Ib[2])


def I_distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def marginal_tensor_distance(
    mu_n: ShellMeasure, mu_a: ShellMeasure, mu_b: ShellMeasure
) -> float:
    """
    Abstand der Bein-Rangfolgen: für jedes Bein vergleiche sortierte κ-Labels von μ_n
    mit vereinigter Rangfolge aus μ_a, μ_b (schematisches ⊗ auf κ-Ebene).
    """
    ma, mb, mn = mu_a.marginals(), mu_b.marginals(), mu_n.marginals()
    dist = 0.0
    for la, _lb, ln in zip(ma, mb, mn):
        # Jaccard-Abstand der Top-Label-Mengen
        sa, sb = set(la[:3]), set(ln[:3])
        union = sa | sb
        dist += 1.0 - len(sa & sb) / len(union) if union else 0.0
    return dist / 4.0


def is_explainable_by_factorization(
    mu_n: ShellMeasure, mu_a: ShellMeasure, mu_b: ShellMeasure, tol: float = REGIME_TOL
) -> bool:
    pred = additive_I_predict(mu_a, mu_b)
    if I_distance(mu_n.I_vector(), pred) <= tol:
        return True
    # Schwächere Bedingung: marginals ähnlich
    if marginal_tensor_distance(mu_n, mu_a, mu_b) <= 0.35:
        return True
    return False


def Z_EABC_count(n: int, cache: dict[int, ShellMeasure]) -> int:
    facs = unordered_factorizations(n, min_factor=2)
    if not facs:
        return 0
    channels: set[tuple] = set()
    for a, b in facs:
        ma, mb = cache.get(a), cache.get(b)
        if ma is None or mb is None:
            continue
        channels.add(channel_signature(ma.mu_signature(), mb.mu_signature()))
    return len(channels)


def Z_regime_count(
    n: int, cache: dict[int, ShellMeasure], tol: float = REGIME_TOL
) -> int:
    mu_n = cache.get(n)
    if mu_n is None:
        return 0
    if is_prime(n):
        return 1
    facs = unordered_factorizations(n, min_factor=2)
    for a, b in facs:
        ma, mb = cache.get(a), cache.get(b)
        if ma is None or mb is None:
            continue
        if is_explainable_by_factorization(mu_n, ma, mb, tol=tol):
            return 0
    return 1


def primary_Z(n: int, cache: dict[int, ShellMeasure]) -> int:
    """Primäre Decompositionszahl Z(n) := Z_EABC(n)."""
    return Z_EABC_count(n, cache)


def delta_Z(n: int, cache: dict[int, ShellMeasure]) -> int:
    return primary_Z(n, cache) - primary_Z(n - 1, cache)


def build_cache(max_n: int) -> dict[int, ShellMeasure]:
    cache: dict[int, ShellMeasure] = {}
    for n in range(2, max_n + 1):
        m = compute_shell_measure(n)
        if m is not None:
            cache[n] = m
    return cache


def z_row(n: int, cache: dict[int, ShellMeasure]) -> dict[str, Any]:
    mu = cache.get(n)
    zf = Z_fact(n)
    ze = Z_EABC_count(n, cache)
    zr = Z_regime_count(n, cache)
    dz = delta_Z(n, cache)
    facs = unordered_factorizations(n, min_factor=2)
    return {
        "n": n,
        "is_prime": is_prime(n),
        "omega": omega_distinct(n),
        "tau": tau_divisor_count(n),
        "Z_fact": zf,
        "Z_EABC": ze,
        "Z_regime": zr,
        "Z": ze,
        "delta_Z": dz,
        "abs_delta_Z": abs(dz),
        "factorization_count": len(facs),
        "shell_size": mu.shell_size if mu else 0,
        "H_n": round(mu.H, 6) if mu else None,
        "chi_n": round(mu.chi, 6) if mu else None,
    }


def prime_vs_composite_jump_stats(
    rows: list[dict[str, Any]], key: str = "abs_delta_Z"
) -> dict[str, Any]:
    primes = [r for r in rows if r["is_prime"]]
    composites = [r for r in rows if not r["is_prime"]]
    mean_p = statistics.mean(r[key] for r in primes) if primes else 0.0
    mean_c = statistics.mean(r[key] for r in composites) if composites else 0.0
    ratio = mean_p / mean_c if mean_c > 1e-12 else float("inf") if mean_p > 0 else 0.0
    top10 = sorted(rows, key=lambda r: -r[key])[:10]
    primes_in_top10 = sum(1 for t in top10 if t["is_prime"])
    return {
        "metric": key,
        "prime_count": len(primes),
        "composite_count": len(composites),
        "mean_prime": round(mean_p, 6),
        "mean_composite": round(mean_c, 6),
        "ratio_mean_prime_over_composite": round(ratio, 6),
        "primes_larger_on_mean": mean_p > mean_c,
        "primes_in_top10": primes_in_top10,
        "top10_jumps": [
            {
                "n": t["n"],
                "is_prime": t["is_prime"],
                "omega": t["omega"],
                key: t[key],
                "Z": t["Z"],
                "delta_Z": t["delta_Z"],
            }
            for t in top10
        ],
    }


def omega_correlation(rows: list[dict[str, Any]], key: str = "Z") -> dict[str, Any]:
    """Prüft ob Z bzw. ΔZ mit ω(n) korreliert (Falsifikation: nur Faktorisierungsarithmetik?)."""
    if len(rows) < 3:
        return {"pearson_r": 0.0, "note": "insufficient_data"}
    xs = [r["omega"] for r in rows]
    ys = [r[key] for r in rows]
    mx, my = statistics.mean(xs), statistics.mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mx) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - my) ** 2 for y in ys))
    r = num / (den_x * den_y) if den_x and den_y else 0.0
    return {
        "pearson_r": round(r, 6),
        "key": key,
        "interpretation": (
            "strong_omega_repackaging"
            if abs(r) > 0.7
            else "partial_omega_signal" if abs(r) > 0.4 else "weak_omega_link"
        ),
    }


def z_fact_vs_z_eabc_correlation(rows: list[dict[str, Any]]) -> dict[str, Any]:
    composites = [r for r in rows if not r["is_prime"] and r["Z_fact"] > 0]
    if len(composites) < 3:
        return {"pearson_r": 0.0, "note": "insufficient_composites"}
    xs = [r["Z_fact"] for r in composites]
    ys = [r["Z_EABC"] for r in composites]
    mx, my = statistics.mean(xs), statistics.mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mx) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - my) ** 2 for y in ys))
    r = num / (den_x * den_y) if den_x and den_y else 0.0
    exact_match = sum(1 for r in composites if r["Z_fact"] == r["Z_EABC"])
    return {
        "pearson_r": round(r, 6),
        "composite_count": len(composites),
        "exact_Z_fact_eq_Z_EABC": exact_match,
        "exact_match_fraction": round(exact_match / len(composites), 6),
    }


def verdict(report: dict[str, Any]) -> str:
    dz = report["prime_vs_composite"]["delta_Z_abs"]
    zr = report["prime_vs_composite"]["Z_regime"]
    omega_z = report["omega_correlation"]["Z"]
    fact_corr = report["z_fact_vs_z_eabc"]
    parts: list[str] = []
    ratio = dz["ratio_mean_prime_over_composite"]
    if dz["primes_larger_on_mean"] and ratio > 1.1:
        parts.append(f"|ΔZ|-Primüberschuss: Ratio {ratio:.3f}")
    else:
        parts.append(f"kein |ΔZ|-Primüberschuss (Ratio {ratio:.3f})")
    if omega_z["interpretation"] == "strong_omega_repackaging":
        parts.append(
            f"Z korreliert stark mit ω (r={omega_z['pearson_r']}) — eher Arithmetik-Repackaging"
        )
    if fact_corr.get("pearson_r", 0) > 0.85:
        parts.append(
            f"Z_EABC ≈ Z_fact (r={fact_corr['pearson_r']}) — kritische Falsifikation"
        )
    if zr["mean_prime"] > zr["mean_composite"]:
        parts.append(
            f"Z_regime: Prim {zr['mean_prime']:.2f} > comp {zr['mean_composite']:.2f}"
        )
    return "; ".join(parts)


def z_decomposition_report(max_n: int = 100) -> dict[str, Any]:
    cache = build_cache(max_n)
    rows = [z_row(n, cache) for n in sorted(cache)]
    pvc_dz = prime_vs_composite_jump_stats(rows, "abs_delta_Z")
    pvc_z = prime_vs_composite_jump_stats(rows, "Z")
    pvc_zr = prime_vs_composite_jump_stats(rows, "Z_regime")
    return {
        "meta": {
            "hypothesis_doc": "collatz_eabc_zerlegungsregimen.md",
            "max_n": max_n,
            "shell_count": len(rows),
            "Z_definitions": {
                "Z_fact": "unordered n=ab, a,b>1",
                "Z_EABC": "distinct (μ_a,μ_b) channel signatures over factorizations",
                "Z_regime": "1 if irreducible regime (prime or μ_n not explainable by any μ_a⊗μ_b)",
                "Z": "primary = Z_EABC",
                "delta_Z": "Z(n)-Z(n-1)",
            },
            "regime_tolerance": REGIME_TOL,
            "links": [
                "collatz_eabc_quaternion_mass_hypothese.md",
                "collatz_eabc_oktonion_singularitaet.md",
                "collatz_eabc_shell_defekt_test.py",
            ],
        },
        "prime_vs_composite": {
            "delta_Z_abs": pvc_dz,
            "Z": pvc_z,
            "Z_regime": pvc_zr,
        },
        "omega_correlation": {
            "Z": omega_correlation(rows, "Z"),
            "abs_delta_Z": omega_correlation(rows, "abs_delta_Z"),
            "Z_regime": omega_correlation(rows, "Z_regime"),
        },
        "z_fact_vs_z_eabc": z_fact_vs_z_eabc_correlation(rows),
        "verdict": "",
        "epistemic_status": (
            "Experiment — testet ob ΔZ an Primzahlen springt (Zerlegungsregimen-Konjektur). "
            "Wenn Z_EABC ≈ Z_fact und mit ω korreliert, ist die Konjektur epistemisch hohl "
            "(klassische Faktorisierung, kein neues EABC-Phänomen)."
        ),
        "rows": rows,
    }


def format_table(report: dict[str, Any]) -> str:
    dz = report["prime_vs_composite"]["delta_Z_abs"]
    lines = [
        "EABC-Zerlegungsregimen Z(n), ΔZ(n)",
        "=" * 52,
        f"n ∈ [2,{report['meta']['max_n']}], Schalen: {report['meta']['shell_count']}",
        f"Verdict: {report['verdict']}",
        "",
        f"|ΔZ|: Prim {dz['mean_prime']:.4f}  comp {dz['mean_composite']:.4f}  "
        f"ratio={dz['ratio_mean_prime_over_composite']:.3f}",
        f"Top-10 |ΔZ|: {dz['primes_in_top10']}/10 prim",
        "",
        "Top Sprünge:",
    ]
    for t in dz["top10_jumps"][:6]:
        tag = "prim" if t["is_prime"] else f"ω={t['omega']}"
        lines.append(f"  n={t['n']} ({tag}): |ΔZ|={t['abs_delta_Z']}  Z={t['Z']}")
    return "\n".join(lines)


def run(max_n: int = 100, output: Path | None = None) -> dict[str, Any]:
    report = z_decomposition_report(max_n=max_n)
    report["verdict"] = verdict(report)
    out = output or DEFAULT_OUTPUT
    out.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    report["output_path"] = str(out)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC-Zerlegungsregimen Z(n), ΔZ(n)")
    parser.add_argument("--max-n", type=int, default=100, help="Oberes n für Σ_n")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(max_n=args.max_n, output=args.output)
    print(format_table(report))
    print(f"\nJSON: {report['output_path']}")


if __name__ == "__main__":
    main()
