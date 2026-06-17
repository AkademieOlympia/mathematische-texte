#!/usr/bin/env python3
"""
Σ→p-Schalen-Defekt-Experiment: I(μ_n), Baseline I_avg(n), Defekt D(n).

Kanonsiche Theorie: collatz_eabc_quaternion_mass_hypothese.md §11

Perspektivwechsel: nicht p ↦ Σ_p, sondern n ↦ Σ_n ↦ μ_n; Primzahlen als
Singularitäten, an denen |D(n)| = |I(μ_n) − I_avg(n)| außergewöhnlich ist.

Ausführung:
    python3 collatz_eabc_shell_defekt_test.py
    python3 collatz_eabc_shell_defekt_test.py --max-n 60 --output collatz_eabc_shell_defekt.json
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from math import isqrt
from pathlib import Path
from typing import Any

from collatz_eabc_hurwitz_orbit_test import (
    K_p_covariance,
    chi_p_from_counts,
    classify_shell_point,
    entropy_from_counts,
    hurwitz_shell_elements,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_shell_defekt.json"


def _sieve_primes(limit: int) -> set[int]:
    if limit < 2:
        return set()
    flags = [True] * (limit + 1)
    flags[0] = flags[1] = False
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            step = p
            start = p * p
            flags[start : limit + 1 : step] = [False] * (((limit - start) // step) + 1)
    return {i for i, ok in enumerate(flags) if ok}


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = isqrt(n)
    for k in range(3, d + 1, 2):
        if n % k == 0:
            return False
    return True


def omega_distinct(n: int) -> int:
    """ω(n): Anzahl verschiedener Primfaktoren."""
    if n < 2:
        return 0
    count = 0
    x = n
    p = 2
    while p * p <= x:
        if x % p == 0:
            count += 1
            while x % p == 0:
                x //= p
        p += 1 if p == 2 else 2
    if x > 1:
        count += 1
    return count


@dataclass(frozen=True, slots=True)
class ShellInvariants:
    n: int
    shell_size: int
    H_n: float
    chi_n: float
    K_trace: float
    distinct_gamma: int
    is_prime: bool
    omega: int

    def I_vector(self) -> tuple[float, float, float]:
        return (self.H_n, self.chi_n, self.K_trace)

    def to_dict(self) -> dict[str, Any]:
        return {
            "n": self.n,
            "shell_size": self.shell_size,
            "H_n": round(self.H_n, 6),
            "chi_n": round(self.chi_n, 6),
            "K_trace": round(self.K_trace, 6),
            "distinct_gamma": self.distinct_gamma,
            "is_prime": self.is_prime,
            "omega": self.omega,
        }


def compute_shell_invariants(n: int) -> ShellInvariants | None:
    raw = hurwitz_shell_elements(n)
    if not raw:
        return None
    points = [classify_shell_point(*coords) for coords in raw]
    gamma_counts: Counter = Counter(pt.gamma for pt in points)
    size = len(points)
    H_n = entropy_from_counts(gamma_counts, size)
    chi_n = chi_p_from_counts(gamma_counts, size)
    K_trace = K_p_covariance(points)["trace"]
    return ShellInvariants(
        n=n,
        shell_size=size,
        H_n=H_n,
        chi_n=chi_n,
        K_trace=K_trace,
        distinct_gamma=len(gamma_counts),
        is_prime=is_prime(n),
        omega=omega_distinct(n),
    )


def rolling_baseline(
    rows: list[ShellInvariants], idx: int, window: int = 5
) -> tuple[float, float, float]:
    """Mittel von I über Nachbarn n±window (ohne n selbst)."""
    n = rows[idx].n
    lo, hi = max(2, n - window), n + window
    peers = [r for r in rows if lo <= r.n <= hi and r.n != n]
    if not peers:
        return rows[idx].I_vector()
    H = statistics.mean(r.H_n for r in peers)
    chi = statistics.mean(r.chi_n for r in peers)
    K = statistics.mean(r.K_trace for r in peers)
    return H, chi, K


def omega_baseline(rows: list[ShellInvariants], target: ShellInvariants) -> tuple[float, float, float]:
    """Mittel von I über alle m mit gleichem ω(m), m ≠ n."""
    peers = [r for r in rows if r.omega == target.omega and r.n != target.n]
    if not peers:
        return statistics.mean(r.H_n for r in rows), statistics.mean(r.chi_n for r in rows), statistics.mean(
            r.K_trace for r in rows
        )
    return (
        statistics.mean(r.H_n for r in peers),
        statistics.mean(r.chi_n for r in peers),
        statistics.mean(r.K_trace for r in peers),
    )


def defect_magnitude(dH: float, dchi: float, dK: float, K_scale: float = 1.0) -> float:
    """Skalierte euklidische Norm des Defektvektors."""
    return math.sqrt(dH * dH + dchi * dchi + (dK / K_scale) ** 2)


def shell_defekt_row(
    inv: ShellInvariants,
    I_roll: tuple[float, float, float],
    I_omega: tuple[float, float, float],
    I_global: tuple[float, float, float],
    K_scale: float,
) -> dict[str, Any]:
    dH_r = inv.H_n - I_roll[0]
    dchi_r = inv.chi_n - I_roll[1]
    dK_r = inv.K_trace - I_roll[2]
    dH_o = inv.H_n - I_omega[0]
    dchi_o = inv.chi_n - I_omega[1]
    dK_o = inv.K_trace - I_omega[2]
    dH_g = inv.H_n - I_global[0]
    dchi_g = inv.chi_n - I_global[1]
    dK_g = inv.K_trace - I_global[2]

    return {
        **inv.to_dict(),
        "I_avg_rolling": [round(x, 6) for x in I_roll],
        "I_avg_omega": [round(x, 6) for x in I_omega],
        "I_avg_global": [round(x, 6) for x in I_global],
        "D_rolling": {
            "H": round(dH_r, 6),
            "chi": round(dchi_r, 6),
            "K_trace": round(dK_r, 6),
            "magnitude": round(defect_magnitude(dH_r, dchi_r, dK_r, K_scale), 6),
        },
        "D_omega": {
            "H": round(dH_o, 6),
            "chi": round(dchi_o, 6),
            "K_trace": round(dK_o, 6),
            "magnitude": round(defect_magnitude(dH_o, dchi_o, dK_o, K_scale), 6),
        },
        "D_global": {
            "H": round(dH_g, 6),
            "chi": round(dchi_g, 6),
            "K_trace": round(dK_g, 6),
            "magnitude": round(defect_magnitude(dH_g, dchi_g, dK_g, K_scale), 6),
        },
    }


def prime_vs_composite_comparison(rows: list[dict[str, Any]], key: str = "D_rolling") -> dict[str, Any]:
    primes = [r[key]["magnitude"] for r in rows if r["is_prime"]]
    composites = [r[key]["magnitude"] for r in rows if not r["is_prime"]]
    if not primes or not composites:
        return {"key": key, "insufficient_data": True}

    mean_p = statistics.mean(primes)
    mean_c = statistics.mean(composites)
    med_p = statistics.median(primes)
    med_c = statistics.median(composites)
    ratio_mean = mean_p / mean_c if mean_c > 1e-12 else float("inf")

    ranked = sorted(rows, key=lambda r: r[key]["magnitude"], reverse=True)
    top10 = [
        {"n": r["n"], "is_prime": r["is_prime"], "magnitude": r[key]["magnitude"], "omega": r["omega"]}
        for r in ranked[:10]
    ]
    prime_in_top10 = sum(1 for t in top10 if t["is_prime"])

    return {
        "key": key,
        "prime_count": len(primes),
        "composite_count": len(composites),
        "mean_abs_D_prime": round(mean_p, 6),
        "mean_abs_D_composite": round(mean_c, 6),
        "median_abs_D_prime": round(med_p, 6),
        "median_abs_D_composite": round(med_c, 6),
        "ratio_mean_prime_over_composite": round(ratio_mean, 6),
        "primes_larger_on_mean": mean_p > mean_c,
        "top10_outliers": top10,
        "primes_in_top10": prime_in_top10,
        "epistemic_note": (
            "Exploratives Muster — kein Beweis, dass Primzahlen systematisch größere |D(n)| tragen. "
            "Kleine N, Schalengröße korreliert mit n; Baseline-Wahl beeinflusst Verhältnis."
        ),
    }


def bernoulli_discretization_hint(rows: list[ShellInvariants]) -> dict[str, Any]:
    """Heuristik: μ_n − μ_∞ als Diskretisierungskorrektur (Bernoulli-Analogie)."""
    if len(rows) < 4:
        return {"hint": "insufficient_data"}
    H_vals = [r.H_n for r in rows]
    chi_vals = [r.chi_n for r in rows]
    H_inf = statistics.mean(H_vals)
    chi_inf = statistics.mean(chi_vals)
    prime_rows = [r for r in rows if r.is_prime]
    comp_rows = [r for r in rows if not r.is_prime]
    return {
        "mu_infinity_proxy": {"H": round(H_inf, 6), "chi": round(chi_inf, 6)},
        "prime_mean_delta_H": round(statistics.mean(r.H_n - H_inf for r in prime_rows), 6) if prime_rows else 0.0,
        "composite_mean_delta_H": round(statistics.mean(r.H_n - H_inf for r in comp_rows), 6) if comp_rows else 0.0,
        "heuristic": "μ_n − μ_∞ als endliche-Schalen-Korrektur (Bernoulli-Heuristik, nicht bewiesen)",
    }


def shell_defekt_report(max_n: int = 50, rolling_window: int = 5) -> dict[str, Any]:
    invariants: list[ShellInvariants] = []
    for n in range(2, max_n + 1):
        inv = compute_shell_invariants(n)
        if inv is not None:
            invariants.append(inv)

    if not invariants:
        return {"max_n": max_n, "rows": [], "verdict": "no_shells"}

    K_scale = max(statistics.mean(abs(r.K_trace) for r in invariants), 0.01)
    I_global = (
        statistics.mean(r.H_n for r in invariants),
        statistics.mean(r.chi_n for r in invariants),
        statistics.mean(r.K_trace for r in invariants),
    )

    rows: list[dict[str, Any]] = []
    for i, inv in enumerate(invariants):
        I_roll = rolling_baseline(invariants, i, window=rolling_window)
        I_omega = omega_baseline(invariants, inv)
        rows.append(shell_defekt_row(inv, I_roll, I_omega, I_global, K_scale))

    cmp_roll = prime_vs_composite_comparison(rows, "D_rolling")
    cmp_omega = prime_vs_composite_comparison(rows, "D_omega")
    cmp_global = prime_vs_composite_comparison(rows, "D_global")
    bernoulli = bernoulli_discretization_hint(invariants)

    verdict_parts: list[str] = []
    if cmp_roll.get("primes_larger_on_mean"):
        verdict_parts.append(
            f"rolling: Prim-Mittel |D|={cmp_roll['mean_abs_D_prime']:.4f} > "
            f"zusammengesetzt {cmp_roll['mean_abs_D_composite']:.4f} "
            f"(Ratio {cmp_roll['ratio_mean_prime_over_composite']:.3f})"
        )
    else:
        verdict_parts.append(
            f"rolling: kein Prim-Überhang (Prim {cmp_roll['mean_abs_D_prime']:.4f} "
            f"vs zusammengesetzt {cmp_roll['mean_abs_D_composite']:.4f})"
        )
    verdict_parts.append(f"Top-10-|D|-Ausreißer: {cmp_roll['primes_in_top10']}/10 sind prim")

    return {
        "meta": {
            "hypothesis_doc": "collatz_eabc_quaternion_mass_hypothese.md",
            "perspective": "Σ→p: n ↦ Σ_n ↦ μ_n; Primzahlen als Defekt-Singularitäten (Conjecture)",
            "invariants": "I(μ_n) = (H_n, χ_n, K_trace)",
            "defect": "D(n) = I(μ_n) − I_avg(n); Baselines: rolling, ω(n)-Stratum, global",
            "max_n": max_n,
            "rolling_window": rolling_window,
            "shell_count": len(rows),
            "K_scale_for_magnitude": round(K_scale, 6),
        },
        "bernoulli_heuristic": bernoulli,
        "prime_vs_composite": {
            "rolling": cmp_roll,
            "omega_stratum": cmp_omega,
            "global": cmp_global,
        },
        "rows": rows,
        "verdict": "; ".join(verdict_parts),
        "epistemic_status": (
            "Experiment — testet heuristische Σ→p-Vermutung auf kleinem n; "
            "negative oder schwache Prim-Separation falsifiziert die naive Singularitätslesart nicht vollständig, "
            "deutet aber auf nötige Verfeinerung (Schalengröße, ω-Stratum, asymptotisches μ_∞)."
        ),
    }


def format_table(report: dict[str, Any]) -> str:
    lines = [
        "Σ→p-Schalen-Defekt D(n) = I(μ_n) − I_avg(n)",
        "=" * 58,
        f"n ∈ [2,{report['meta']['max_n']}], Schalen: {report['meta']['shell_count']}",
        f"Verdict: {report['verdict']}",
        "",
    ]
    pvc = report["prime_vs_composite"]["rolling"]
    lines.extend(
        [
            f"Prim |D| Mittel: {pvc['mean_abs_D_prime']:.4f}  "
            f"zusammengesetzt: {pvc['mean_abs_D_composite']:.4f}  "
            f"Ratio: {pvc['ratio_mean_prime_over_composite']:.3f}",
            "",
            "Top Ausreißer (rolling |D|):",
        ]
    )
    for t in pvc["top10_outliers"][:6]:
        tag = "prim" if t["is_prime"] else f"ω={t['omega']}"
        lines.append(f"  n={t['n']} ({tag}): |D|={t['magnitude']:.4f}")
    return "\n".join(lines)


def run(max_n: int = 50, rolling_window: int = 5, output: Path | None = None) -> dict[str, Any]:
    report = shell_defekt_report(max_n=max_n, rolling_window=rolling_window)
    out = output or DEFAULT_OUTPUT
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report["output_path"] = str(out)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Σ→p-Schalen-Defekt-Experiment")
    parser.add_argument("--max-n", type=int, default=50, help="Oberes n für Σ_n (klein halten)")
    parser.add_argument("--rolling-window", type=int, default=5, help="Fenster für rollendes I_avg")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(max_n=args.max_n, rolling_window=args.rolling_window, output=args.output)
    print(format_table(report))
    print(f"\nJSON: {report['output_path']}")


if __name__ == "__main__":
    main()
