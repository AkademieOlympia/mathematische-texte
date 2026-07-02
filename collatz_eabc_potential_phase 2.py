#!/usr/bin/env python3
"""
EABC reine Potentialverbindungen: Bohm-, AB-, Berry-Stubs (diskret).

Theorie: collatz_eabc_potential_geometrie.md
Brachistochrone: collatz_eabc_brachistochrone.md (v = f(V), T = ∫ ds/v)
Chiraler Transport: collatz_eabc_chirale_transport.py (φ_R, φ_L, U_E)
Holonomie-Stufen: collatz_eabc_holonomie_stufen.md §2–3
Epistemik: collatz_eabc_epistemik_physik.md (v ≠ c, kein SRT)
Wigner: collatz_eabc_wigner_analog.md, collatz_eabc_wigner_field.py

  bohm_like_velocity(∇V)     — Führungsgeschwindigkeit aus Potentialgradient (Stub)
  aharonov_bohm_phase(γ, A)  — diskretes ∮ A auf C₄-Kanten
  berry_phase_difference     — φ_R − φ_L (Berry-/Holonomie-Observable)

Ausführung:
    python3 collatz_eabc_potential_phase.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

from collatz_eabc_chirale_transport import (
    DEFAULT_THETA_EDGE,
    holonomy_phase_difference,
    holonomy_unitary_phases,
)
from collatz_eabc_kritische_abbildung import ABCEA_EDGE_KEYS
from collatz_eabc_transition_graph import ABCEA_WORD, CEABC_WORD

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_potential_phase.json"
THEORY = "collatz_eabc_potential_geometrie.md"
THEORY_BRACHISTOCHRONE = "collatz_eabc_brachistochrone.md"
THEORY_CHIRAL = "collatz_eabc_chirale_polarisation.md"
THEORY_HOLONOMIE_STUFEN = "collatz_eabc_holonomie_stufen.md"
THEORY_EPISTEMIK = "collatz_eabc_epistemik_physik.md"
THEORY_WIGNER = "collatz_eabc_wigner_analog.md"

EdgeKey = str
AField = Mapping[EdgeKey, float]
LoopOrientation = str

CANONICAL_LOOPS: tuple[LoopOrientation, ...] = (ABCEA_WORD, CEABC_WORD)


def default_ab_edge_field(
    *,
    theta_edge: float = DEFAULT_THETA_EDGE,
) -> dict[str, float]:
    """Standard-A-Feld: gleiche Kantenphase θ_edge auf ABCEA-Kanten."""
    return {key: theta_edge for key in ABCEA_EDGE_KEYS}


def _edges_for_loop(loop: LoopOrientation) -> tuple[str, ...]:
    if loop == ABCEA_WORD:
        return ABCEA_EDGE_KEYS
    if loop == CEABC_WORD:
        return tuple(reversed(ABCEA_EDGE_KEYS))
    raise ValueError(f"unknown loop orientation: {loop!r}")


def bohm_like_velocity(
    potential_gradient: float,
    *,
    v0: float = 1.0,
    beta: float = 0.01,
    v_min: float = 1e-6,
) -> float:
    """
    Bohm-Analog: Führungsgeschwindigkeit aus Potentialgradient (Stub).

    v = v0 + β · ∇V  (Modellparameter — nicht Lichtgeschwindigkeit c).
    """
    return max(v0 + beta * potential_gradient, v_min)


def aharonov_bohm_phase(
    loop: LoopOrientation,
    a_field: AField,
    *,
    sign_from_orientation: bool = True,
) -> float:
    """
    Diskrete AB-Phase: Σ_{e∈γ} A(e) auf EABC-C₄-Kanten.

    Bei sign_from_orientation=True wird das Vorzeichen der CEABC-Schleife
    invertiert (ω(γ)=−1), analog zu holonomy_sign.
    """
    edges = _edges_for_loop(loop)
    total = sum(float(a_field.get(e, 0.0)) for e in edges)
    if sign_from_orientation and loop == CEABC_WORD:
        total = -total
    return total


def berry_phase_difference(
    phi_r: float,
    phi_l: float,
) -> float:
    """Berry-/Holonomie-Observable: φ_R − φ_L."""
    return phi_r - phi_l


def berry_phase_from_loop(
    loop: LoopOrientation,
    *,
    theta_edge: float = DEFAULT_THETA_EDGE,
) -> dict[str, float | int]:
    """
    Verknüpft Berry-Stub mit chiralem Transport: φ_R, φ_L, Δφ aus Schleife.
    """
    unitary = holonomy_unitary_phases(loop, theta_edge=theta_edge)
    phi_r = float(unitary["phi_R"])
    phi_l = float(unitary["phi_L"])
    return {
        "loop": loop,
        "phi_R": phi_r,
        "phi_L": phi_l,
        "berry_phase_difference": berry_phase_difference(phi_r, phi_l),
        "holonomy_phase_difference": holonomy_phase_difference(
            loop, theta_edge=theta_edge
        ),
        "omega": int(unitary["omega"]),
    }


def potential_phase_report(
    *,
    theta_edge: float = DEFAULT_THETA_EDGE,
) -> dict[str, Any]:
    """Kompakter Bericht: Bohm-Stub, AB-Phase, Berry-Differenz für beide Schleifen."""
    a_field = default_ab_edge_field(theta_edge=theta_edge)
    grad_v = bohm_like_velocity(100.0) - bohm_like_velocity(0.0)

    loops: dict[str, Any] = {}
    for loop in CANONICAL_LOOPS:
        loops[loop] = {
            "aharonov_bohm_phase": aharonov_bohm_phase(loop, a_field),
            "berry": berry_phase_from_loop(loop, theta_edge=theta_edge),
        }

    return {
        "meta": {
            "module": "collatz_eabc_potential_phase.py",
            "theory": THEORY,
            "theory_brachistochrone": THEORY_BRACHISTOCHRONE,
            "theory_chiral": THEORY_CHIRAL,
            "theory_holonomie_stufen": THEORY_HOLONOMIE_STUFEN,
            "theory_epistemik": THEORY_EPISTEMIK,
            "theory_wigner": THEORY_WIGNER,
            "epistemic": "Bohm/AB/Berry = Analogie — v, θ_edge Modellparameter, nicht c",
            "theta_edge": theta_edge,
        },
        "bohm_stub": {
            "potential_gradient_delta_v": grad_v,
            "v_at_zero_grad": bohm_like_velocity(0.0),
            "v_at_unit_grad": bohm_like_velocity(1.0),
        },
        "a_field": dict(a_field),
        "loops": loops,
        "boxed": {
            "not_c": "Effektive Geschwindigkeiten und Kantenphasen sind Modellparameter",
            "ab": "Σ A(e) auf C₄-Schleife — diskrete AB-Analogie",
            "berry": "φ_R − φ_L = Observable auf chiralem Faserbündel (Stufe 2)",
            "cylinder_dual_path": "geplant — Vierling-Doppelpfad (collatz_eabc_potential_geometrie.md §5)",
        },
    }


def run(
    output: Path = DEFAULT_OUTPUT,
    *,
    theta_edge: float = DEFAULT_THETA_EDGE,
) -> dict[str, Any]:
    report = potential_phase_report(theta_edge=theta_edge)
    output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="EABC Bohm/AB/Berry Potential-Phase Stubs"
    )
    parser.add_argument("--theta-edge", type=float, default=DEFAULT_THETA_EDGE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(output=args.output, theta_edge=args.theta_edge)
    abcea = report["loops"][ABCEA_WORD]
    print("=== EABC Potential-Phase (Bohm/AB/Berry) ===")
    print(f"θ_edge = {report['meta']['theta_edge']:.6f}")
    print(f"AB ABCEA: φ_AB = {abcea['aharonov_bohm_phase']:.4f}")
    print(
        f"Berry Δφ = {abcea['berry']['berry_phase_difference']:.4f}"
    )
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
