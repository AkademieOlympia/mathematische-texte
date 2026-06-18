#!/usr/bin/env python3
"""
EABC Brachistochrone / Fermat-Prinzip: fünf Potenzial-Kandidaten, Wegzeit T = Σ Δs/v.

Theorie: collatz_eabc_brachistochrone.md (Modell/Hypothese — kein Physikanspruch)
Chirale Erweiterung: collatz_eabc_chirale_polarisation.md (T_R, T_L, Birefringenz)
Epistemik: collatz_eabc_epistemik_physik.md (Wegfunktion, Nicht-SRT)
Verknüpfung: collatz_eabc_kritische_abbildung.py (ABCEA-Trajektorie, s_v)

  T = ∫ ds / v(x),   v(x) = f(V(x))
  Vergleich: gerader Polygonzug vs. gestörter Pfad (Gradient / Halbkreis)

Ausführung:
    python3 collatz_eabc_brachistochrone.py
    python3 collatz_eabc_brachistochrone.py --max-p 100000
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

from collatz_eabc_holonomie_fehlerterm import holonomy_counts
from collatz_eabc_kritische_abbildung import (
    CANONICAL_GAP_PATTERN,
    GAMMA_1_APPROX,
    SOURCE_X,
    holonomy_sensor_trajectory,
    semicircle_arc_length,
    zeta_imaginary_parts,
)
from collatz_eabc_transition_graph import (
    ABCEA_WORD,
    classes_from_sequence,
    prime_eabc_sequence,
    transition_counts,
    transition_probabilities,
)
from eabc_from_lean import EClass

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_brachistochrone.json"
THEORY = "collatz_eabc_brachistochrone.md"
THEORY_KRITISCHE = "collatz_eabc_kritische_abbildung.md"
THEORY_HOLONOMIE_STUFEN = "collatz_eabc_holonomie_stufen.md"
THEORY_EPISTEMIK = "collatz_eabc_epistemik_physik.md"
THEORY_CHIRAL = "collatz_eabc_chirale_polarisation.md"
THEORY_GENERALANGRIFF = "collatz_generalangriff_2026.md"

V_MIN = 1e-6
N_ZETA_ZEROS = 20
PERTURBATION_SCALE = 0.15

Point2D = tuple[float, float]
VelocityFunc = Callable[[float, float], float]


def velocity_from_potential(
    V: float, model: str, *, v0: float = 1.0, alpha: float = 0.01
) -> float:
    """
    v = f(V) für die fünf Potenzialfamilien (Modellwahl).

    model:
      log           — v ∝ ln(x) bei V = ln(x)
      inverse_log   — v ∝ 1/ln(x)
      zeta          — v = v0 + V (ζ-Summenpotential)
      chirality     — v = v0 + α·V mit V = D_E
      curvature     — v = v0 / (1 + |V|)  (Krümmungshemmung)
      information   — v = exp(-V) = P  (Informationspotential)
    """
    if model == "log":
        return max(V, V_MIN)
    if model == "inverse_log":
        return max(1.0 / max(V, V_MIN), V_MIN)
    if model == "zeta":
        return max(v0 + V, V_MIN)
    if model == "chirality":
        return max(v0 + alpha * V, V_MIN)
    if model == "curvature":
        return max(v0 / (1.0 + abs(V)), V_MIN)
    if model == "information":
        return max(math.exp(-V), V_MIN)
    raise ValueError(f"unknown velocity model: {model}")


def segment_length(p0: Point2D, p1: Point2D) -> float:
    return math.hypot(p1[0] - p0[0], p1[1] - p0[1])


def birefringent_velocity_pair(
    v_base: float,
    *,
    d_e_proxy: float = 0.0,
    v0: float = 1.0,
    alpha: float = 0.01,
) -> tuple[float, float]:
    """Zwei Polarisationskanäle: v_R = v0 + α·V_E, v_L = v0 - α·V_E."""
    v_r = max(v0 + alpha * d_e_proxy, V_MIN)
    v_l = max(v0 - alpha * d_e_proxy, V_MIN)
    if v_base > 0:
        scale = v_base / max(v0, V_MIN)
        v_r = max(v_r * scale, V_MIN)
        v_l = max(v_l * scale, V_MIN)
    return v_r, v_l


def travel_time_birefringent(
    path_points: list[Point2D] | tuple[Point2D, ...],
    v_func: VelocityFunc,
    *,
    d_e_global: float = 0.0,
    alpha: float = 0.01,
    v0: float = 1.0,
) -> dict[str, float]:
    """T_R = Σ Δs/v_R, T_L = Σ Δs/v_L — Birefringenz-Analogie."""
    pts = list(path_points)
    if len(pts) < 2:
        return {
            "T_R": 0.0,
            "T_L": 0.0,
            "delta_T": 0.0,
            "ratio_T_R_over_T_L": float("nan"),
        }
    x_end = pts[-1][0]
    x_start = pts[0][0]
    span = max(x_end - x_start, 1e-9)
    t_r = 0.0
    t_l = 0.0
    for i in range(len(pts) - 1):
        ds = segment_length(pts[i], pts[i + 1])
        xm = 0.5 * (pts[i][0] + pts[i + 1][0])
        gm = 0.5 * (pts[i][1] + pts[i + 1][1])
        frac = (xm - x_start) / span
        d_e_local = d_e_global * frac
        v_base = max(v_func(xm, gm), V_MIN)
        v_r, v_l = birefringent_velocity_pair(
            v_base, d_e_proxy=d_e_local, v0=v0, alpha=alpha
        )
        t_r += ds / v_r
        t_l += ds / v_l
    delta = t_r - t_l
    ratio = t_r / t_l if t_l > 0 else float("nan")
    return {"T_R": t_r, "T_L": t_l, "delta_T": delta, "ratio_T_R_over_T_L": ratio}


def travel_time_integral(
    path_points: list[Point2D] | tuple[Point2D, ...],
    v_func: VelocityFunc,
) -> float:
    """Diskrete Wegzeit T = Σ_j Δs_j / v(x̄_j, γ̄_j) entlang Polygonzug."""
    pts = list(path_points)
    if len(pts) < 2:
        return 0.0
    total = 0.0
    for i in range(len(pts) - 1):
        ds = segment_length(pts[i], pts[i + 1])
        xm = 0.5 * (pts[i][0] + pts[i + 1][0])
        gm = 0.5 * (pts[i][1] + pts[i + 1][1])
        v_mid = max(v_func(xm, gm), V_MIN)
        total += ds / v_mid
    return total


def abcea_trajectory_xy(
    gamma_ref: float = GAMMA_1_APPROX,
    gaps: tuple[int, ...] | None = None,
) -> list[Point2D]:
    """ABCEA-Holonomie-Sensor als (x, γ)-Punkte."""
    traj = holonomy_sensor_trajectory("ABCEA", gaps=gaps, gamma_ref=gamma_ref)
    return [(seg["x"], seg["gamma_cumulative"]) for seg in traj["segments"]]


def straight_path(points: list[Point2D]) -> list[Point2D]:
    """Gerader Polygonzug durch die Trajektorienknoten."""
    return list(points)


def perturbed_path_gradient(
    points: list[Point2D],
    grad_func: Callable[[float, float], Point2D],
    scale: float = PERTURBATION_SCALE,
) -> list[Point2D]:
    """
    Gestörter Pfad: pro Kante Mittelpunkt in Richtung ∇V verschoben (Fermat-Heuristik).
    """
    if len(points) < 2:
        return list(points)
    out: list[Point2D] = [points[0]]
    for i in range(len(points) - 1):
        p0, p1 = points[i], points[i + 1]
        xm = 0.5 * (p0[0] + p1[0])
        gm = 0.5 * (p0[1] + p1[1])
        gx, gy = grad_func(xm, gm)
        norm = math.hypot(gx, gy)
        if norm > 0:
            ds = segment_length(p0, p1)
            bump = scale * ds
            out.append((xm + bump * gx / norm, gm + bump * gy / norm))
        out.append(p1)
    return out


def perturbed_path_semicircle_complex(
    points: list[Point2D],
) -> float:
    """
    Halbkreis-Kette in der komplexen Ebene (s_v-Punkte) — liefert Gesamtwegzeit bei v≡1.
    Für variable v: diskretisiere Bogen in N Segmente.
    """
    z_pts = [complex(x, g) for x, g in points]
    if len(z_pts) < 2:
        return 0.0
    total = 0.0
    for i in range(len(z_pts) - 1):
        total += semicircle_arc_length(z_pts[i], z_pts[i + 1])
    total += semicircle_arc_length(z_pts[-1], z_pts[0])
    return total


def semicircle_path_discretized(
    points: list[Point2D],
    n_segments_per_arc: int = 8,
) -> list[Point2D]:
    """Halbkreis-Bögen als Polygonzug in (x, γ) (obere Halbebene)."""
    z_pts = [complex(x, g) for x, g in points]
    if len(z_pts) < 2:
        return list(points)
    poly: list[Point2D] = []
    pairs = [(z_pts[i], z_pts[i + 1]) for i in range(len(z_pts) - 1)]
    pairs.append((z_pts[-1], z_pts[0]))
    for z1, z2 in pairs:
        poly.append((z1.real, z1.imag))
        for k in range(1, n_segments_per_arc):
            t = k / n_segments_per_arc
            angle1 = math.atan2(z1.imag, z1.real)
            angle2 = math.atan2(z2.imag, z2.real)
            if angle2 < angle1:
                angle2 += 2 * math.pi
            mid_angle = 0.5 * (angle1 + angle2)
            r1 = abs(z1)
            r2 = abs(z2)
            r = r1 + t * (r2 - r1)
            ang = angle1 + t * (angle2 - angle1)
            z = r * complex(math.cos(ang), math.sin(ang))
            poly.append((z.real, z.imag))
    poly.append((z_pts[0].real, z_pts[0].imag))
    return poly


def potential_log(x: float, _gamma: float) -> float:
    """V(x) = ln(x) — Primdichte-Heuristik."""
    return math.log(max(x, SOURCE_X + 1e-9))


def potential_log_log(x: float, _gamma: float) -> float:
    """V(x) = ln ln(x)."""
    lx = math.log(max(x, SOURCE_X + 1e-9))
    return math.log(max(lx, 1e-9))


def potential_zeta_sum(
    gamma: float, gammas: list[float], epsilon: float = 0.5
) -> float:
    """V(γ) = Σ_n 1/((γ - γ_n)² + ε) auf der kritischen Linie."""
    return sum(1.0 / ((gamma - gn) ** 2 + epsilon) for gn in gammas)


def build_d_e_lookup(max_p: int) -> dict[str, Any]:
    """Globales D_E und Schritt-Differenzen entlang der Primfolge."""
    row = holonomy_counts(max_p)
    seq = prime_eabc_sequence(max_p)
    n = len(seq)
    return {
        "D_E_global": row["D_E"],
        "N_plus": row["N_plus"],
        "N_minus": row["N_minus"],
        "prime_count": n,
        "max_p": max_p,
    }


def d_e_at_fraction(frac: float, d_e_global: float) -> float:
    """Lokales D_E-Proxy: skaliert mit Trajektorienfortschritt (Modell)."""
    return d_e_global * frac


def build_pattern_info(max_p: int) -> dict[str, Any]:
    """Übergangswahrscheinlichkeiten und lokale Muster-Entropie."""
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    counts = transition_counts(classes)
    probs = transition_probabilities(counts)
    class_to_idx = {c.value: i for i, c in enumerate(EClass)}
    return {
        "transition_probabilities": probs,
        "class_to_idx": class_to_idx,
        "total_transitions": sum(sum(row) for row in counts),
    }


def _trajectory_fraction(x: float) -> float:
    """Fortschritt entlang ABCEA-Sensor in [0, 1]."""
    x_end = SOURCE_X + sum(CANONICAL_GAP_PATTERN)
    return min(max((x - SOURCE_X) / max(x_end - SOURCE_X, 1e-9), 0.0), 1.0)


def pattern_probability_at_vertex(
    vertex: str,
    pattern_info: dict[str, Any],
) -> float:
    """P — Zeilensumme der Übergangswkt. ausgehend von vertex (Proxy)."""
    idx = pattern_info["class_to_idx"][vertex]
    row = pattern_info["transition_probabilities"][idx]
    return sum(row) / max(len(row), 1)


def make_velocity_func(
    potential_name: str,
    *,
    gammas: list[float] | None = None,
    d_e_global: float = 0.0,
    pattern_info: dict[str, Any] | None = None,
    velocity_model: str | None = None,
) -> tuple[VelocityFunc, Callable[[float, float], Point2D]]:
    """Potenzial + zugehörige Geschwindigkeit und Gradient (numerisch)."""
    if gammas is None:
        gammas = zeta_imaginary_parts(N_ZETA_ZEROS)

    def V_at(x: float, gamma: float) -> float:
        if potential_name == "log":
            return potential_log(x, gamma)
        if potential_name == "log_log":
            return potential_log_log(x, gamma)
        if potential_name == "zeta":
            return potential_zeta_sum(gamma, gammas)
        if potential_name == "chirality":
            return d_e_at_fraction(_trajectory_fraction(x), d_e_global)
        if potential_name == "curvature":
            return abs(gamma) / max(gammas[0], 1.0)
        if potential_name == "information":
            vertices = list(ABCEA_WORD)
            frac = _trajectory_fraction(x)
            idx = min(int(frac * (len(vertices) - 1)), len(vertices) - 1)
            if pattern_info is None:
                return 1.0
            p = pattern_probability_at_vertex(vertices[idx], pattern_info)
            return -math.log(max(p, 1e-12))
        raise ValueError(f"unknown potential: {potential_name}")

    v_model = (
        velocity_model
        or {
            "log": "log",
            "log_log": "inverse_log",
            "zeta": "zeta",
            "chirality": "chirality",
            "curvature": "curvature",
            "information": "information",
        }[potential_name]
    )

    def v_func(x: float, gamma: float) -> float:
        return velocity_from_potential(V_at(x, gamma), v_model)

    def grad_func(x: float, gamma: float, h: float = 1e-4) -> Point2D:
        vx = (V_at(x + h, gamma) - V_at(x - h, gamma)) / (2 * h)
        vg = (V_at(x, gamma + h) - V_at(x, gamma - h)) / (2 * h)
        return (vx, vg)

    return v_func, grad_func


def compare_paths_for_potential(
    potential_name: str,
    points: list[Point2D],
    *,
    gammas: list[float] | None = None,
    d_e_global: float = 0.0,
    pattern_info: dict[str, Any] | None = None,
    perturbation_scale: float = PERTURBATION_SCALE,
) -> dict[str, Any]:
    """Gerade vs. Gradient- und Halbkreis-Pfade für ein Potenzial."""
    v_func, grad_func = make_velocity_func(
        potential_name,
        gammas=gammas,
        d_e_global=d_e_global,
        pattern_info=pattern_info,
    )
    straight = straight_path(points)
    bent = perturbed_path_gradient(points, grad_func, scale=perturbation_scale)
    semi_poly = semicircle_path_discretized(points)

    t_straight = travel_time_integral(straight, v_func)
    t_bent = travel_time_integral(bent, v_func)
    t_semi = travel_time_integral(semi_poly, v_func)
    biref_straight = travel_time_birefringent(straight, v_func, d_e_global=d_e_global)
    biref_bent = travel_time_birefringent(bent, v_func, d_e_global=d_e_global)

    rel_bent = (t_bent - t_straight) / t_straight if t_straight > 0 else float("nan")
    rel_semi = (t_semi - t_straight) / t_straight if t_straight > 0 else float("nan")
    bent_faster = t_bent < t_straight * (1.0 - 1e-6)
    semi_faster = t_semi < t_straight * (1.0 - 1e-6)
    significant_bend = abs(rel_bent) > 0.01 or abs(rel_semi) > 0.01

    return {
        "potential": potential_name,
        "T_straight": t_straight,
        "T_gradient_perturbed": t_bent,
        "T_semicircle": t_semi,
        "rel_gradient_vs_straight": rel_bent,
        "rel_semicircle_vs_straight": rel_semi,
        "gradient_faster_than_straight": bent_faster,
        "semicircle_faster_than_straight": semi_faster,
        "path_bends_significantly": significant_bend,
        "birefringence_straight": biref_straight,
        "birefringence_gradient": biref_bent,
        "n_straight_segments": len(straight) - 1,
        "n_gradient_points": len(bent),
    }


def run_comparison(
    max_p: int = 100_000,
    gamma_ref: float = GAMMA_1_APPROX,
    perturbation_scale: float = PERTURBATION_SCALE,
) -> dict[str, Any]:
    """Vollständiger Vergleich aller fünf Potenziale auf ABCEA-Trajektorie."""
    points = abcea_trajectory_xy(gamma_ref=gamma_ref)
    gammas = zeta_imaginary_parts(N_ZETA_ZEROS)
    d_e_info = build_d_e_lookup(max_p)
    pattern_info = build_pattern_info(max_p)

    potentials = ("log", "zeta", "chirality", "curvature", "information")
    comparisons = [
        compare_paths_for_potential(
            name,
            points,
            gammas=gammas,
            d_e_global=float(d_e_info["D_E_global"]),
            pattern_info=pattern_info,
            perturbation_scale=perturbation_scale,
        )
        for name in potentials
    ]

    any_significant = any(c["path_bends_significantly"] for c in comparisons)
    any_gradient_wins = any(c["gradient_faster_than_straight"] for c in comparisons)

    return {
        "meta": {
            "module": "collatz_eabc_brachistochrone.py",
            "theory": THEORY,
            "theory_kritische_abbildung": THEORY_KRITISCHE,
            "theory_holonomie_stufen": THEORY_HOLONOMIE_STUFEN,
            "theory_epistemik": THEORY_EPISTEMIK,
            "theory_chiral": THEORY_CHIRAL,
            "theory_generalangriff": THEORY_GENERALANGRIFF,
            "epistemic": "Variationsprinzip / Modell — kein Physikanspruch, kein SRT",
            "max_p": max_p,
            "gamma_ref": gamma_ref,
            "n_zeta_zeros": N_ZETA_ZEROS,
            "perturbation_scale": perturbation_scale,
        },
        "fermat_setup": {
            "formula": "T = sum_j Delta_s_j / v(x_j, gamma_j)",
            "velocity": "v = f(V) — siehe velocity_from_potential",
            "birefringence": "T_R = sum ds/v_R, T_L = sum ds/v_L — collatz_eabc_chirale_polarisation.md",
            "potentials": {
                "log": "V(x) = ln(x), v ~ ln(x) oder 1/ln(x)",
                "zeta": "V(gamma) = sum_n 1/((gamma-gamma_n)^2 + eps)",
                "chirality": "V_E = D_E, v = v0 + alpha D_E",
                "curvature": "V ~ |gamma|/gamma_1 (Winkeldefekt-Proxy)",
                "information": "V = -log P(x) aus Übergangszählern",
            },
        },
        "abcea_trajectory": {
            "points_xy": [{"x": x, "gamma": g} for x, g in points],
            "gaps": list(CANONICAL_GAP_PATTERN),
            "x_endpoint": points[-1][0] if points else SOURCE_X,
        },
        "d_e_context": d_e_info,
        "pattern_context": {
            "total_transitions": pattern_info["total_transitions"],
            "ABCEA_word": ABCEA_WORD,
        },
        "comparisons": comparisons,
        "summary": {
            "any_path_bends_significantly": any_significant,
            "any_gradient_perturbation_faster": any_gradient_wins,
            "verdict": (
                "Mindestens ein Potenzial zeigt merkliche Abweichung gerade vs. gestört."
                if any_significant
                else "Auf ABCEA-Sensor bleiben Wegzeiten bei allen V näherungsweise geradlinig."
            ),
            "research_question": (
                "Bleibt der optimale Pfad zwischen 1/2 und X gerade "
                "oder biegt er zu hoher Primdichte / Holonomie?"
            ),
        },
        "boxed": {
            "principle": "Fermat: T = integral ds/v(x) — Modellabbildung",
            "hierarchy": "1. ln(x), 2. zeta-Summe, 3. D_E",
            "not_physics": "Kein SRT-Anspruch — collatz_eabc_epistemik_physik.md",
        },
    }


def run(
    max_p: int = 100_000,
    output: Path = DEFAULT_OUTPUT,
    gamma_ref: float = GAMMA_1_APPROX,
) -> dict[str, Any]:
    report = run_comparison(max_p=max_p, gamma_ref=gamma_ref)
    output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="EABC Brachistochrone / Fermat-Prinzip"
    )
    parser.add_argument("--max-p", type=int, default=100_000)
    parser.add_argument("--gamma-ref", type=float, default=GAMMA_1_APPROX)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(max_p=args.max_p, output=args.output, gamma_ref=args.gamma_ref)
    print("=== EABC Brachistochrone / Fermat ===")
    print(f"ABCEA-Trajektorie: {len(report['abcea_trajectory']['points_xy'])} Knoten")
    print(f"D_E @ max_p={args.max_p}: {report['d_e_context']['D_E_global']:+d}")
    print()
    for row in report["comparisons"]:
        print(
            f"  {row['potential']:12s}: T_str={row['T_straight']:.4f}  "
            f"T_grad={row['T_gradient_perturbed']:.4f} ({row['rel_gradient_vs_straight']:+.2%})  "
            f"bends={row['path_bends_significantly']}"
        )
    print()
    print(report["summary"]["verdict"])
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
