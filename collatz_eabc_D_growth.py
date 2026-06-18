#!/usr/bin/env python3
"""
EABC: Wachstumsdiagnostik D_E(X) und Dirichlet-Charakter-Stub mod 12.

Theorie: collatz_eabc_evolution_analytik.md

  D_E(X) = N_plus(X) - N_minus(X)  auf Gitter X in [10^3, 10^6]
  Heuristische Klassifikation: O(1), O(log X), O(sqrt X), Potenzgesetz
  Stub: D_E ~ sum_chi a_chi sum_{p<=X} chi(p)

Ausführung:
    python3 collatz_eabc_D_growth.py
    python3 collatz_eabc_D_growth.py --max-x 1000000
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from collatz_eabc_holonomie_fehlerterm import holonomy_counts
from collatz_eabc_transition_graph import prime_eabc_sequence

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_D_growth.json"
THEORY_EVOLUTION = "collatz_eabc_evolution_analytik.md"
THEORY_ZIRKULATION = "collatz_eabc_zirkulationshypothese.md"
THEORY_HOLONOMY_STAGES = "collatz_eabc_holonomie_stufen.md"

DEFAULT_GRID = (
    1_000,
    2_000,
    5_000,
    10_000,
    20_000,
    50_000,
    100_000,
    200_000,
    500_000,
    1_000_000,
)

GROWTH_SCENARIOS = ("O(1)", "O(log X)", "O(sqrt X)", "power_law")
HOLONOMY_FALLS = ("A", "B", "C")


def default_x_grid(max_x: int = 1_000_000) -> tuple[int, ...]:
    return tuple(x for x in DEFAULT_GRID if x <= max_x)


def d_e_at_x(x: int) -> dict[str, Any]:
    """D_E und Normalisierung an Prim-Obergrenze X."""
    row = holonomy_counts(x)
    n_total = row["N_plus"] + row["N_minus"]
    return {
        "X": x,
        "N_plus": row["N_plus"],
        "N_minus": row["N_minus"],
        "N_total": n_total,
        "D_E": row["D_E"],
        "S_E": row["S_E"],
        "D_tilde_E": row["D_tilde_E"],
        "C_E": row["C_E"],
    }


def d_e_series(grid: tuple[int, ...] | list[int]) -> list[dict[str, Any]]:
    return [d_e_at_x(x) for x in grid]


def _rss(y: list[float], yhat: list[float]) -> float:
    return sum((a - b) ** 2 for a, b in zip(y, yhat, strict=True))


def _fit_affine(
    xs: list[float], ys: list[float], transform
) -> tuple[float, float, float]:
    """y ~ a + b * transform(x); returns (a, b, rss)."""
    z = [transform(x) for x in xs]
    n = len(xs)
    if n < 2:
        mean_y = ys[0] if ys else 0.0
        return mean_y, 0.0, _rss(ys, [mean_y] * n)
    z_mean = sum(z) / n
    y_mean = sum(ys) / n
    var_z = sum((t - z_mean) ** 2 for t in z)
    if var_z < 1e-15:
        return y_mean, 0.0, _rss(ys, [y_mean] * n)
    cov = sum((t - z_mean) * (y - y_mean) for t, y in zip(z, ys, strict=True))
    b = cov / var_z
    a = y_mean - b * z_mean
    yhat = [a + b * t for t in z]
    return a, b, _rss(ys, yhat)


def _fit_power_law(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    """y ~ c * x^alpha on points with y > 0; returns (c, alpha, rss)."""
    pairs = [(x, y) for x, y in zip(xs, ys, strict=True) if y > 0]
    if len(pairs) < 2:
        c = ys[-1] if ys else 0.0
        return c, 0.0, _rss(ys, [c] * len(ys))
    lx = [math.log(x) for x, _ in pairs]
    ly = [math.log(y) for _, y in pairs]
    _, alpha, rss_log = _fit_affine(lx, ly, lambda t: t)
    alpha = max(0.0, alpha)
    c = math.exp(sum(ly) / len(ly) - alpha * sum(lx) / len(lx))
    yhat_all = [c * (x**alpha) if y > 0 else 0.0 for x, y in zip(xs, ys, strict=True)]
    return c, alpha, _rss(ys, yhat_all)


def classify_growth(series: list[dict[str, Any]]) -> dict[str, Any]:
    """
    Heuristische Wachstumsklassifikation über Modell-RSS.

    Szenarien: O(1), O(log X), O(sqrt X), power_law c*X^alpha.
    """
    xs = [float(r["X"]) for r in series]
    ys = [float(r["D_E"]) for r in series]
    n = len(series)

    const = ys[-1] if ys else 0.0
    rss_const = _rss(ys, [const] * n)

    _, b_log, rss_log = _fit_affine(xs, ys, math.log)
    _, b_sqrt, rss_sqrt = _fit_affine(xs, ys, math.sqrt)
    c_pow, alpha, rss_pow = _fit_power_law(xs, ys)

    models = {
        "O(1)": {"rss": rss_const, "params": {"c": const}},
        "O(log X)": {"rss": rss_log, "params": {"b": b_log}},
        "O(sqrt X)": {"rss": rss_sqrt, "params": {"b": b_sqrt}},
        "power_law": {"rss": rss_pow, "params": {"c": c_pow, "alpha": alpha}},
    }
    ranked = sorted(models.items(), key=lambda kv: kv[1]["rss"])
    best = ranked[0][0]

    last = series[-1] if series else {}
    x_max = float(last.get("X", 1))
    d_max = float(last.get("D_E", 0))
    n_total = float(last.get("N_total", 0))

    diagnostics = {
        "D_E_at_max": d_max,
        "log_X_at_max": math.log(x_max) if x_max > 1 else 0.0,
        "sqrt_X_at_max": math.sqrt(x_max),
        "D_over_logX": d_max / math.log(x_max) if x_max > 1 else 0.0,
        "D_over_sqrtX": d_max / math.sqrt(x_max) if x_max > 0 else 0.0,
        "D_over_sqrtN": d_max / math.sqrt(n_total) if n_total > 0 else 0.0,
    }

    # Zusatz-Heuristik: harte Ausschlüsse für Szenarien A und C bei großem X
    reject_O1 = d_max > 10 and (series[-2]["D_E"] if len(series) > 1 else 0) < d_max
    reject_sqrt = d_max < 0.15 * math.sqrt(x_max) if x_max > 0 else True
    reject_log = d_max > 3.0 * math.log(x_max) if x_max > math.e else False

    notes: list[str] = []
    if reject_O1:
        notes.append("Szenario A (O(1)) empirisch ausgeschlossen: D_E wächst mit X.")
    if reject_sqrt:
        notes.append(
            "Szenario C (O(sqrt X)) empirisch ausgeschlossen: |D_E| << sqrt(X)."
        )
    if reject_log:
        notes.append(
            "Szenario B (O(log X)) schwach: |D_E| deutlich größer als c·log X."
        )

    scenario_map = {
        "O(1)": "A",
        "O(log X)": "B",
        "O(sqrt X)": "C",
        "power_law": "D",
    }
    preferred = best
    if reject_O1 and preferred == "O(1)":
        preferred = ranked[1][0]
    if reject_sqrt and preferred == "O(sqrt X)":
        preferred = next((name for name, _ in ranked if name != "O(sqrt X)"), best)

    return {
        "models": models,
        "ranked_by_rss": [name for name, _ in ranked],
        "best_fit_model": best,
        "preferred_scenario": preferred,
        "scenario_letter": scenario_map.get(preferred, "?"),
        "diagnostics_at_max_X": diagnostics,
        "heuristic_notes": notes,
        "epistemic": "Experiment (heuristische Modellwahl, kein Theorem)",
    }


def classify_holonomy_growth(series: list[dict[str, Any]]) -> dict[str, Any]:
    """
    Holonomie-Wachstum in N = N_+ + N_- (Fall A/B/C).

    Fall A: D_E = O(1)           -> Hol_E = 0
    Fall B: D_E = O(sqrt N)      -> Hol_E -> 0  (Nullhypothese)
    Fall C: D_E ~ alpha N        -> Hol_E = alpha != 0
    """
    ns = [float(r["N_total"]) for r in series]
    ys = [float(r["D_E"]) for r in series]
    n = len(series)

    mean_y = sum(ys) / n if n else 0.0
    rss_o1 = _rss(ys, [mean_y] * n)

    _, b_sqrt, rss_sqrt = _fit_affine(ns, ys, math.sqrt)
    _, alpha_lin, rss_linear = _fit_affine(ns, ys, lambda t: t)

    models = {
        "A": {
            "label": "O(1)",
            "rss": rss_o1,
            "params": {"mean_D_E": mean_y},
            "Hol_E_limit": 0.0,
        },
        "B": {
            "label": "O(sqrt N)",
            "rss": rss_sqrt,
            "params": {"b": b_sqrt},
            "Hol_E_limit": 0.0,
        },
        "C": {
            "label": "alpha N",
            "rss": rss_linear,
            "params": {"alpha": alpha_lin},
            "Hol_E_limit": alpha_lin,
        },
    }
    ranked = sorted(models.items(), key=lambda kv: kv[1]["rss"])
    best = ranked[0][0]

    last = series[-1] if series else {}
    n_max = float(last.get("N_total", 0))
    d_max = float(last.get("D_E", 0))
    s_max = float(last.get("S_E", 0))
    d_tilde_max = float(last.get("D_tilde_E", 0))

    diagnostics = {
        "N_at_max": n_max,
        "D_E_at_max": d_max,
        "S_E_at_max": s_max,
        "D_tilde_E_at_max": d_tilde_max,
        "D_over_sqrtN": d_max / math.sqrt(n_max) if n_max > 0 else 0.0,
        "alpha_hat": d_max / n_max if n_max > 0 else 0.0,
    }

    reject_a = d_max > 5 and (series[-2]["D_E"] if len(series) > 1 else 0) < d_max
    reject_b = s_max > 0.08 and d_tilde_max > 1.5
    reject_c = s_max < 0.02 and n_max > 50

    notes: list[str] = []
    if reject_a:
        notes.append("Fall A (D_E=O(1)) ausgeschlossen: |D_E| wächst mit N.")
    if reject_b:
        notes.append(
            "Fall B schwach: S_E bleibt deutlich > 0 und D̃_E nicht O(1)-stabil."
        )
    if reject_c:
        notes.append("Fall C ausgeschlossen: S_E nahe 0 bei großem N.")

    preferred = best
    if reject_a and preferred == "A":
        preferred = ranked[1][0]
    if reject_b and preferred == "B":
        preferred = next((fall for fall, _ in ranked if fall != "B"), best)
    if reject_c and preferred == "C":
        preferred = next((fall for fall, _ in ranked if fall != "C"), best)

    hol_e_reading = {
        "A": "Hol_E = 0 (absoluter Effekt stirbt)",
        "B": "Hol_E -> 0 (Nullhypothese, Random Walk)",
        "C": f"Hol_E = alpha ~ {alpha_lin:.4f} (stabile Chiralitaet)",
    }

    return {
        "theory": THEORY_HOLONOMY_STAGES,
        "models": models,
        "ranked_by_rss": [fall for fall, _ in ranked],
        "best_fit_fall": best,
        "preferred_fall": preferred,
        "Hol_E_reading": hol_e_reading.get(preferred, "?"),
        "diagnostics_at_max_N": diagnostics,
        "heuristic_notes": notes,
        "epistemic": "Experiment (heuristische Modellwahl in N, kein Theorem)",
    }


def _chi_mod12_characters() -> dict[str, Any]:
    """
      Nichttriviale Charaktere auf (Z/12Z)^x = {1,5,7,11} ≅ C2 x C2.

      chi_4: Legendre-Symbol (-1/p) für ungerade p (entspricht p mod 4).
      chi_3: quadratischer Charakter mod 3 (p mod 3 == 1 -> 1, == 2 -> -1).
    chi_12: Produkt chi_4 * chi_3 (primitiver Charakter mod 12).
    """

    def chi_trivial(p: int) -> int:
        return 1

    def chi_4(p: int) -> int:
        if p == 2:
            return 0
        return -1 if p % 4 == 3 else 1

    def chi_3(p: int) -> int:
        if p == 3:
            return 0
        r = p % 3
        if r == 1:
            return 1
        if r == 2:
            return -1
        return 0

    def chi_12(p: int) -> int:
        return chi_4(p) * chi_3(p)

    return {
        "chi_0": {"order": 1, "fn": chi_trivial},
        "chi_4": {"order": 2, "fn": chi_4},
        "chi_3": {"order": 2, "fn": chi_3},
        "chi_12": {"order": 2, "fn": chi_12},
    }


def character_sums(max_x: int) -> dict[str, int]:
    """sum_{p <= X} chi(p) für Primzahlen p > 3."""
    chars = _chi_mod12_characters()
    sums = {name: 0 for name in chars}
    for row in prime_eabc_sequence(max_x):
        p = row["p"]
        for name, spec in chars.items():
            sums[name] += spec["fn"](p)
    return sums


def dirichlet_decomposition_stub(max_x: int) -> dict[str, Any]:
    """
    Stub: Schätze a_chi in D_E(X) ≈ sum_chi a_chi S_chi(X).

    Lineare Projektion auf nichttriviale Charakter-Summen (experimentell).
    """
    row = holonomy_counts(max_x)
    d_e = row["D_E"]
    sums = character_sums(max_x)

    # Design: D_E ~ a4*S4 + a3*S3 + a12*S12 (ohne trivialem chi_0)
    design = [
        ("chi_4", sums["chi_4"]),
        ("chi_3", sums["chi_3"]),
        ("chi_12", sums["chi_12"]),
    ]
    # Einfache proportionale Schätzung: a_chi = D_E * S_chi / sum_j S_j^2
    denom = sum(s * s for _, s in design)
    coeffs: dict[str, float] = {}
    if denom > 0:
        for name, s in design:
            coeffs[name] = d_e * s / denom
    else:
        for name, _ in design:
            coeffs[name] = 0.0

    reconstructed = sum(coeffs[name] * s for name, s in design)
    return {
        "X": max_x,
        "D_E": d_e,
        "character_sums": sums,
        "dirichlet_coefficients": coeffs,
        "reconstructed_D_E": reconstructed,
        "residual": d_e - reconstructed,
        "formula": "D_E(X) ≈ sum_{chi nontrivial} a_chi sum_{p<=X} chi(p)",
        "status": "stub (lineare Projektion, nicht bewiesen)",
        "theory": THEORY_EVOLUTION,
    }


def c4_laplacian_spectrum() -> dict[str, Any]:
    """Kombinatorisches Spec(L_{C4}^{sym}) = {0, 2, 2, 4}."""
    return {
        "eigenvalues_symmetrized": [0.0, 2.0, 2.0, 4.0],
        "spectral_gap": 2.0,
        "note": "Reines C4-Gerüst; prim-gewichtetes L_E(X) variabel (graph_laplacian.py)",
    }


def growth_report(
    max_x: int = 1_000_000,
    grid: tuple[int, ...] | list[int] | None = None,
) -> dict[str, Any]:
    if grid is None:
        grid = default_x_grid(max_x)
    series = d_e_series(grid)
    classification = classify_growth(series)
    last_x = series[-1]["X"] if series else max_x
    holonomy = classify_holonomy_growth(series)
    return {
        "meta": {
            "module": "collatz_eabc_D_growth.py",
            "theory_evolution": THEORY_EVOLUTION,
            "theory_zirkulation": THEORY_ZIRKULATION,
            "theory_holonomy_stages": THEORY_HOLONOMY_STAGES,
            "max_x": max_x,
            "grid": list(grid),
        },
        "D_E_series": series,
        "growth_classification": classification,
        "holonomy_growth": holonomy,
        "dirichlet_stub": dirichlet_decomposition_stub(last_x),
        "c4_spectrum": c4_laplacian_spectrum(),
        "boxed_conclusion": {
            "evolution": "Bell -> Sagnac -> C_E -> Spec(L_E)",
            "strongest_next_step": "D_E als mod-12-Chebyshev-Race via Dirichlet-Charaktere",
            "growth_at_max_X": {
                "X": last_x,
                "D_E": series[-1]["D_E"] if series else 0,
                "scenario": classification["preferred_scenario"],
                "scenario_letter": classification["scenario_letter"],
            },
            "holonomy_at_max_N": {
                "N": series[-1]["N_total"] if series else 0,
                "D_E": series[-1]["D_E"] if series else 0,
                "S_E": series[-1]["S_E"] if series else 0,
                "fall": holonomy["preferred_fall"],
                "Hol_E_reading": holonomy["Hol_E_reading"],
            },
        },
    }


def run(
    max_x: int = 1_000_000,
    output: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    report = growth_report(max_x=max_x)
    output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC D_E Wachstumsdiagnostik")
    parser.add_argument("--max-x", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(max_x=args.max_x, output=args.output)
    gc = report["growth_classification"]
    last = report["D_E_series"][-1]
    print("=== EABC D_E Wachstumsdiagnostik ===")
    print(f"Theorie: {THEORY_EVOLUTION}")
    print()
    for row in report["D_E_series"]:
        print(
            f"X={row['X']:>7}: D_E={row['D_E']:+4d}  "
            f"N={row['N_total']:4d}  D̃_E={row['D_tilde_E']:+.3f}  S_E={row['S_E']:+.4f}"
        )
    print()
    print(f"Best fit (RSS): {gc['best_fit_model']}")
    print(f"Bevorzugtes Szenario: {gc['scenario_letter']} — {gc['preferred_scenario']}")
    diag = gc["diagnostics_at_max_X"]
    print(
        f"@ X={last['X']}: D_E={diag['D_E_at_max']:.0f}, "
        f"D/log X={diag['D_over_logX']:.2f}, D/sqrt X={diag['D_over_sqrtX']:.4f}"
    )
    for note in gc["heuristic_notes"]:
        print(f"  • {note}")
    hg = report["holonomy_growth"]
    print()
    print(f"Holonomie-Fall (in N): {hg['preferred_fall']} — {hg['Hol_E_reading']}")
    diag_n = hg["diagnostics_at_max_N"]
    print(
        f"@ N={diag_n['N_at_max']:.0f}: S_E={diag_n['S_E_at_max']:.4f}, "
        f"alpha_hat={diag_n['alpha_hat']:.4f}, D/sqrt N={diag_n['D_over_sqrtN']:.3f}"
    )
    for note in hg["heuristic_notes"]:
        print(f"  • {note}")
    stub = report["dirichlet_stub"]
    print()
    print("Dirichlet-Stub a_chi:")
    for name, val in stub["dirichlet_coefficients"].items():
        print(f"  {name}: {val:+.4f}")
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
