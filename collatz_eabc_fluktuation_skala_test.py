#!/usr/bin/env python3
"""
Forschungsproblem A: Skalierung des EABC-Fluktuationsfelds bei großem x.

Testet empirisch, ob H(x)/π(x), χ(x)/√π(x), H(x)/√π(x) und Modenkoeffizienten
c_i(x) Grenzwerte haben oder vorhersagbar skalieren (1/√x, 1/x, konstant, …).

Kanonsiche Definitionen: collatz_eabc_invarianzprogramm.md §8, Forschungsproblem A.

Ausführung:
    python3 collatz_eabc_fluktuation_skala_test.py
    python3 collatz_eabc_fluktuation_skala_test.py --max-x 1000000
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from collatz_eabc_invarianzprogramm import (
    CountVector,
    EABC_ORDER,
    _sieve_primes,
    delta_at_x,
    h_at_x,
    kappa,
    mode_coefficients,
    snapshot_at_x,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_fluktuation_skala.json"

DEFAULT_GRID = [100, 500, 1_000, 5_000, 10_000, 50_000, 100_000, 500_000, 1_000_000]


@dataclass(frozen=True, slots=True)
class ScaleRow:
    x: int
    pi: int
    h: float
    chi: float
    chi_fluct: float
    h_over_pi: float
    h_over_sqrt_pi: float
    chi_over_sqrt_pi: float
    c1: float
    c2: float
    c3: float
    c1_over_sqrt_pi: float
    c2_over_sqrt_pi: float
    c3_over_sqrt_pi: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "x": self.x,
            "pi_eabc": self.pi,
            "H": self.h,
            "chi": self.chi,
            "chi_fluct": self.chi_fluct,
            "H_over_pi": self.h_over_pi,
            "H_over_sqrt_pi": self.h_over_sqrt_pi,
            "chi_over_sqrt_pi": self.chi_over_sqrt_pi,
            "mode_c": {"c1": self.c1, "c2": self.c2, "c3": self.c3},
            "mode_c_over_sqrt_pi": {
                "c1": self.c1_over_sqrt_pi,
                "c2": self.c2_over_sqrt_pi,
                "c3": self.c3_over_sqrt_pi,
            },
        }


def snapshots_at_grid(
    grid: list[int],
    primes: list[int] | None = None,
) -> list[ScaleRow]:
    """Einzelnes Sieb; inkrementelle Zählung an Gitterpunkten (O(π(max x)))."""
    if not grid:
        return []
    grid_sorted = sorted(set(x for x in grid if x >= 5))
    max_x = grid_sorted[-1]
    if primes is None:
        primes = _sieve_primes(max_x)

    counts = {cls: 0 for cls in EABC_ORDER}
    rows: list[ScaleRow] = []
    gi = 0

    def emit(x_target: int) -> None:
        v = CountVector(
            e=counts[EABC_ORDER[0]],
            a=counts[EABC_ORDER[1]],
            b=counts[EABC_ORDER[2]],
            c=counts[EABC_ORDER[3]],
        )
        pi = v.total
        delta = delta_at_x(x_target, v)
        h_val = h_at_x(x_target, v, delta)
        chi_fluct = (delta.e + delta.c) - (delta.a + delta.b)
        chi_val = chi_fluct / pi if pi else 0.0
        sqrt_pi = math.sqrt(pi) if pi else 1.0
        c1, c2, c3 = mode_coefficients(delta)
        rows.append(
            ScaleRow(
                x=x_target,
                pi=pi,
                h=h_val,
                chi=chi_val,
                chi_fluct=chi_fluct,
                h_over_pi=h_val / pi if pi else 0.0,
                h_over_sqrt_pi=h_val / sqrt_pi if pi else 0.0,
                chi_over_sqrt_pi=chi_val / sqrt_pi if pi else 0.0,
                c1=c1,
                c2=c2,
                c3=c3,
                c1_over_sqrt_pi=c1 / sqrt_pi if pi else 0.0,
                c2_over_sqrt_pi=c2 / sqrt_pi if pi else 0.0,
                c3_over_sqrt_pi=c3 / sqrt_pi if pi else 0.0,
            )
        )

    for p in primes:
        if p > max_x:
            break
        while gi < len(grid_sorted) and p > grid_sorted[gi]:
            emit(grid_sorted[gi])
            gi += 1
        cls = kappa(p)
        if cls is not None:
            counts[cls] += 1

    while gi < len(grid_sorted):
        emit(grid_sorted[gi])
        gi += 1

    return rows


def _ols_loglog(x_vals: list[float], y_vals: list[float]) -> dict[str, float | None]:
    """log(y) = intercept + slope * log(x); nur positive Paare."""
    pairs = [
        (math.log(x), math.log(y))
        for x, y in zip(x_vals, y_vals)
        if x > 0 and y > 0 and math.isfinite(x) and math.isfinite(y)
    ]
    n = len(pairs)
    if n < 2:
        return {"slope": None, "intercept": None, "r_squared": None, "n": n}
    sx = sum(p[0] for p in pairs)
    sy = sum(p[1] for p in pairs)
    sxx = sum(p[0] ** 2 for p in pairs)
    sxy = sum(p[0] * p[1] for p in pairs)
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-30:
        return {"slope": None, "intercept": None, "r_squared": None, "n": n}
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    y_mean = sy / n
    ss_tot = sum((p[1] - y_mean) ** 2 for p in pairs)
    ss_res = sum((p[1] - (intercept + slope * p[0])) ** 2 for p in pairs)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return {"slope": slope, "intercept": intercept, "r_squared": r2, "n": n}


def _ols_linear(x_vals: list[float], y_vals: list[float]) -> dict[str, float | None]:
    """y = intercept + slope * x (für signed χ/√π vs 1/√π etc.)."""
    pairs = [
        (x, y)
        for x, y in zip(x_vals, y_vals)
        if math.isfinite(x) and math.isfinite(y)
    ]
    n = len(pairs)
    if n < 2:
        return {"slope": None, "intercept": None, "r_squared": None, "n": n}
    sx = sum(p[0] for p in pairs)
    sy = sum(p[1] for p in pairs)
    sxx = sum(p[0] ** 2 for p in pairs)
    sxy = sum(p[0] * p[1] for p in pairs)
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-30:
        return {"slope": None, "intercept": None, "r_squared": None, "n": n}
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    y_mean = sy / n
    ss_tot = sum((p[1] - y_mean) ** 2 for p in pairs)
    ss_res = sum((p[1] - (intercept + slope * p[0])) ** 2 for p in pairs)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return {"slope": slope, "intercept": intercept, "r_squared": r2, "n": n}


def _interpret_h_slope(slope: float | None) -> str:
    """Steigung von log(H) vs log(π): 0→H~const, 1→H~π, 0.5→H~√π."""
    if slope is None:
        return "unbestimmt"
    if abs(slope) < 0.08:
        return "H ~ konstant (H/π ~ 1/π)"
    if abs(slope - 1.0) < 0.12:
        return "H ~ π (H/π ~ konstant)"
    if abs(slope - 0.5) < 0.12:
        return "H ~ √π (H/√π ~ konstant)"
    return f"H ~ π^{slope:.2f}"


def fit_scaling_hypotheses(rows: list[ScaleRow]) -> dict[str, Any]:
    """Log-log- und lineare Fits für Forschungsproblem A."""
    pi_vals = [float(r.pi) for r in rows]
    x_vals = [float(r.x) for r in rows]
    inv_sqrt_pi = [1.0 / math.sqrt(r.pi) if r.pi else 0.0 for r in rows]
    inv_pi = [1.0 / r.pi if r.pi else 0.0 for r in rows]
    sqrt_pi = [math.sqrt(r.pi) for r in rows]

    h_vals = [r.h for r in rows]
    h_over_pi = [r.h_over_pi for r in rows]
    h_over_sqrt_pi = [r.h_over_sqrt_pi for r in rows]
    chi_sqrt = [r.chi_over_sqrt_pi for r in rows]
    chi_fluct = [r.chi_fluct for r in rows]
    c1s = [r.c1 for r in rows]
    c2s = [r.c2 for r in rows]
    c3s = [r.c3 for r in rows]

    fits: dict[str, Any] = {
        "log_H_vs_log_pi": _ols_loglog(pi_vals, h_vals),
        "log_H_over_pi_vs_log_pi": _ols_loglog(pi_vals, h_over_pi),
        "log_H_over_sqrt_pi_vs_log_pi": _ols_loglog(pi_vals, h_over_sqrt_pi),
        "log_abs_chi_fluct_vs_log_pi": _ols_loglog(
            pi_vals, [abs(v) for v in chi_fluct]
        ),
        "log_abs_c1_vs_log_pi": _ols_loglog(pi_vals, [abs(v) for v in c1s]),
        "log_abs_c2_vs_log_pi": _ols_loglog(pi_vals, [abs(v) for v in c2s]),
        "log_abs_c3_vs_log_pi": _ols_loglog(pi_vals, [abs(v) for v in c3s]),
        "chi_over_sqrt_pi_vs_inv_sqrt_pi": _ols_linear(inv_sqrt_pi, chi_sqrt),
        "chi_fluct_vs_sqrt_pi": _ols_linear(sqrt_pi, chi_fluct),
        "H_over_pi_vs_inv_pi": _ols_linear(inv_pi, h_over_pi),
        "H_vs_sqrt_pi": _ols_linear(sqrt_pi, h_vals),
    }

    h_slope = fits["log_H_vs_log_pi"]["slope"]
    h_pi_slope = fits["log_H_over_pi_vs_log_pi"]["slope"]
    chi_lin = fits["chi_over_sqrt_pi_vs_inv_sqrt_pi"]

    hypotheses = [
        {
            "name": "H/π konstant",
            "predicted_log_slope_H_over_pi_vs_pi": 0.0,
            "observed_log_slope": h_pi_slope,
            "residual": abs(h_pi_slope - 0.0) if h_pi_slope is not None else None,
        },
        {
            "name": "H/√π konstant (H ~ √π)",
            "predicted_log_slope_H_vs_pi": 0.5,
            "observed_log_slope_H_vs_pi": h_slope,
            "residual": abs(h_slope - 0.5) if h_slope is not None else None,
        },
        {
            "name": "H konstant (H/π ~ 1/π)",
            "predicted_log_slope_H_vs_pi": 0.0,
            "observed_log_slope_H_vs_pi": h_slope,
            "residual": abs(h_slope - 0.0) if h_slope is not None else None,
        },
        {
            "name": "χ/√π konstant",
            "predicted_slope_vs_1_over_sqrt_pi": 0.0,
            "observed_slope": chi_lin["slope"],
            "residual": abs(chi_lin["slope"]) if chi_lin["slope"] is not None else None,
        },
        {
            "name": "χ_fluct ~ √π",
            "predicted_slope_vs_sqrt_pi": 1.0,
            "observed_slope": fits["chi_fluct_vs_sqrt_pi"]["slope"],
            "residual": (
                abs(fits["chi_fluct_vs_sqrt_pi"]["slope"] - 1.0)
                if fits["chi_fluct_vs_sqrt_pi"]["slope"] is not None
                else None
            ),
        },
    ]

    valid = [h for h in hypotheses if h["residual"] is not None]
    best = min(valid, key=lambda h: h["residual"]) if valid else None

    return {
        "fits": fits,
        "hypothesis_scores": hypotheses,
        "best_fit": best,
        "H_scaling_interpretation": _interpret_h_slope(
            h_slope if isinstance(h_slope, float) else None
        ),
        "chi_over_sqrt_pi_range": {
            "min": min(chi_sqrt),
            "max": max(chi_sqrt),
            "first": chi_sqrt[0],
            "last": chi_sqrt[-1],
            "spread": max(chi_sqrt) - min(chi_sqrt),
        },
        "H_over_pi_range": {
            "min": min(h_over_pi),
            "max": max(h_over_pi),
            "first": h_over_pi[0],
            "last": h_over_pi[-1],
        },
    }


def run_scaling_experiment(
    grid: list[int] | None = None,
) -> dict[str, Any]:
    if grid is None:
        grid = DEFAULT_GRID
    grid = sorted(set(x for x in grid if x >= 5))
    t0 = time.perf_counter()
    rows = snapshots_at_grid(grid)
    elapsed = time.perf_counter() - t0

    # Kreuzcheck gegen snapshot_at_x für kleinste Punkte
    primes = _sieve_primes(grid[-1])
    for r in rows[:3]:
        ref = snapshot_at_x(r.x, primes)
        assert abs(ref.h - r.h) < 1e-9
        assert abs(ref.chi - r.chi) < 1e-12

    scaling = fit_scaling_hypotheses(rows)
    return {
        "program": "EABC-Fluktuation-Skalierung (Forschungsproblem A)",
        "canonical_doc": "collatz_eabc_invarianzprogramm.md",
        "label": "Experiment",
        "grid_x": grid,
        "elapsed_seconds": round(elapsed, 3),
        "rows": [r.to_dict() for r in rows],
        "scaling_analysis": scaling,
        "epistemic_note": (
            "Exploratives Skalierungsexperiment; keine Behauptung von Grenzwerten "
            "oder Beweisen. Log-log-Steigungen sind heuristisch über endliche Gitterpunkte."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Forschungsproblem A: EABC-Fluktuationsfeld-Skalierung"
    )
    parser.add_argument(
        "--max-x",
        type=int,
        default=None,
        help="Obergrenze; Gitter wird auf x ≤ max-x beschnitten",
    )
    parser.add_argument(
        "--grid",
        type=int,
        nargs="*",
        default=None,
        help="Explizite x-Gitterpunkte (Standard: 100..1M)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="JSON-Ausgabedatei",
    )
    args = parser.parse_args()

    grid = list(args.grid) if args.grid else list(DEFAULT_GRID)
    if args.max_x is not None:
        grid = [x for x in grid if x <= args.max_x]
        if not grid:
            grid = [max(5, args.max_x)]

    result = run_scaling_experiment(grid)
    args.output.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    )
    print(f"Wrote {args.output} ({len(grid)} Gitterpunkte, {result['elapsed_seconds']}s)")
    print(
        f"{'x':>10} {'π':>8} {'H/π':>10} {'H/√π':>10} "
        f"{'χ/√π':>12} {'c1':>8} {'c2':>8} {'c3':>8}"
    )
    for row in result["rows"]:
        mc = row["mode_c"]
        print(
            f"{row['x']:10d} {row['pi_eabc']:8d} "
            f"{row['H_over_pi']:10.5f} {row['H_over_sqrt_pi']:10.3f} "
            f"{row['chi_over_sqrt_pi']:+12.6f} "
            f"{mc['c1']:+8.3f} {mc['c2']:+8.3f} {mc['c3']:+8.3f}"
        )
    best = result["scaling_analysis"]["best_fit"]
    if best:
        print(f"\nBest-fit Hypothese (kleinster Residual): {best['name']}")
    print(
        f"H-Skalierung: {result['scaling_analysis']['H_scaling_interpretation']}"
    )


if __name__ == "__main__":
    main()
