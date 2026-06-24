#!/usr/bin/env python3
"""
EABC chiraler Transport: Helizität λ=±1, Phasenkanäle φ_R/φ_L, Holonomie-U_E.

Theorie: collatz_eabc_chirale_polarisation.md
Stufe-2-Upgrade: collatz_eabc_holonomie_stufen.md §2
Zirkulation: collatz_eabc_zirkulationshypothese.md (C_E = N_R - N_L = D_E)

  ABCEA ↔ λ=+1 (R),  CEABC ↔ λ=-1 (L)
  φ += ω(γ) · θ_edge  entlang diskretem Pfad
  U_E = diag(e^{iφ_R}, e^{iφ_L});  Observable: φ_R - φ_L

Ausführung:
    python3 collatz_eabc_chirale_transport.py
    python3 collatz_eabc_chirale_transport.py --max-p 100000
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Literal

from collatz_eabc_holonomie_fehlerterm import holonomy_counts
from collatz_eabc_kritische_abbildung import (
    ABCEA_EDGE_KEYS,
    CANONICAL_GAP_PATTERN,
    holonomy_sign,
)
from collatz_eabc_transition_graph import (
    ABCEA_WORD,
    CEABC_WORD,
    classes_from_sequence,
    omega_hol,
    prime_eabc_sequence,
    sliding_windows,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_chirale_transport.json"
THEORY = "collatz_eabc_chirale_polarisation.md"
THEORY_HOLONOMIE_STUFEN = "collatz_eabc_holonomie_stufen.md"
THEORY_ZIRKULATION = "collatz_eabc_zirkulationshypothese.md"
THEORY_EPISTEMIK = "collatz_eabc_epistemik_physik.md"

HelicityChannel = Literal["R", "L"]
Orientation = Literal["ABCEA", "CEABC"]

CHANNEL_FOR_ORIENTATION: dict[Orientation, HelicityChannel] = {
    "ABCEA": "R",
    "CEABC": "L",
}
ORIENTATION_FOR_CHANNEL: dict[HelicityChannel, Orientation] = {
    "R": "ABCEA",
    "L": "CEABC",
}

DEFAULT_THETA_EDGE = math.pi / 4


def helicity_channel(word: str) -> HelicityChannel | None:
    """Mappe erkanntes 5-Wort auf Helizitätskanal R/L."""
    if word == ABCEA_WORD:
        return "R"
    if word == CEABC_WORD:
        return "L"
    return None


def orientation_from_word(word: str) -> Orientation | None:
    if word == ABCEA_WORD:
        return "ABCEA"
    if word == CEABC_WORD:
        return "CEABC"
    return None


def theta_edge_from_gaps(
    gaps: tuple[int, ...] = CANONICAL_GAP_PATTERN,
    *,
    normalize: bool = True,
) -> float:
    """Kantenphase θ_edge aus Lückenlängen ℓ_j (Holonomie-Sensor-Proxy)."""
    if not gaps:
        return DEFAULT_THETA_EDGE
    total = float(sum(gaps))
    mean_ell = total / len(gaps)
    return (math.pi / 2.0) * mean_ell / total if normalize else mean_ell


def edge_phases_for_loop(
    orientation: Orientation,
    *,
    theta_edge: float = DEFAULT_THETA_EDGE,
    n_edges: int = 4,
) -> dict[str, float]:
    """Phasenbeiträge pro Kante entlang geschlossener Orientierung."""
    sign = holonomy_sign(orientation)
    keys = (
        list(ABCEA_EDGE_KEYS)
        if orientation == "ABCEA"
        else list(reversed(ABCEA_EDGE_KEYS))
    )
    return {key: sign * theta_edge for key in keys[:n_edges]}


def accumulate_phases_along_windows(
    windows: list[dict[str, Any]],
    *,
    theta_edge: float = DEFAULT_THETA_EDGE,
) -> dict[str, float]:
    """Akkumuliere φ_R, φ_L entlang diskretem Pfad (Gleitfenster)."""
    phi_r = 0.0
    phi_l = 0.0
    n_r = 0
    n_l = 0
    for w in windows:
        word = w["word"]
        om = omega_hol(word)
        if om == 1:
            phi_r += om * theta_edge
            n_r += 1
        elif om == -1:
            phi_l += (-om) * theta_edge
            n_l += 1
    return {
        "phi_R": phi_r,
        "phi_L": phi_l,
        "N_R": n_r,
        "N_L": n_l,
        "theta_edge": theta_edge,
    }


def holonomy_phase_difference(
    gamma_loop: str,
    *,
    theta_edge: float = DEFAULT_THETA_EDGE,
) -> float:
    """Phasendifferenz φ_R - φ_L nach Schleife γ_loop (ABCEA oder CEABC)."""
    om = omega_hol(gamma_loop)
    if om == 0:
        return 0.0
    n_edges = 4
    phi_r = max(om, 0) * n_edges * theta_edge
    phi_l = max(-om, 0) * n_edges * theta_edge
    return phi_r - phi_l


def holonomy_unitary_phases(
    gamma_loop: str,
    *,
    theta_edge: float = DEFAULT_THETA_EDGE,
) -> dict[str, complex | float | int]:
    """U_E = diag(e^{iφ_R}, e^{iφ_L}) für eine orientierte Schleife."""
    om = omega_hol(gamma_loop)
    n_edges = 4
    phi_r = max(om, 0) * n_edges * theta_edge
    phi_l = max(-om, 0) * n_edges * theta_edge
    return {
        "phi_R": phi_r,
        "phi_L": phi_l,
        "phase_difference": phi_r - phi_l,
        "U_R": complex(math.cos(phi_r), math.sin(phi_r)),
        "U_L": complex(math.cos(phi_l), math.sin(phi_l)),
        "omega": om,
    }


def chirality_flux_from_counts(n_r: int, n_l: int) -> dict[str, int | float]:
    """Diskreter Chiralitätsfluss C_E = N_R - N_L = D_E."""
    d_e = n_r - n_l
    total = n_r + n_l
    s_e = d_e / total if total > 0 else 0.0
    return {
        "N_R": n_r,
        "N_L": n_l,
        "C_E": d_e,
        "D_E": d_e,
        "S_E": s_e,
        "chirality_flux": d_e,
    }


def link_phase_difference_to_D_E(
    phi_diff: float,
    d_e: int,
    *,
    theta_edge: float = DEFAULT_THETA_EDGE,
) -> dict[str, float | int | bool]:
    """Konsistenzcheck: φ_R - φ_L vs. D_E · θ_edge (Modell)."""
    expected_from_d_e = d_e * theta_edge
    ratio = phi_diff / expected_from_d_e if expected_from_d_e != 0 else float("nan")
    return {
        "phi_difference": phi_diff,
        "D_E": d_e,
        "theta_edge": theta_edge,
        "expected_D_E_times_theta": expected_from_d_e,
        "ratio_phi_to_D_E_theta": ratio,
        "consistent_sign": (phi_diff == 0 and d_e == 0)
        or (phi_diff > 0 and d_e > 0)
        or (phi_diff < 0 and d_e < 0),
    }


def chiral_transport_report(
    max_p: int = 100_000,
    *,
    theta_edge: float | None = None,
) -> dict[str, Any]:
    """Vollständiger Bericht: Phasenakkumulation, Holonomie, D_E-Verknüpfung."""
    if theta_edge is None:
        theta_edge = theta_edge_from_gaps()

    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    windows = sliding_windows(classes, width=5)
    phases = accumulate_phases_along_windows(windows, theta_edge=theta_edge)
    counts = holonomy_counts(max_p)
    flux = chirality_flux_from_counts(phases["N_R"], phases["N_L"])

    phi_diff_abcea = holonomy_phase_difference(ABCEA_WORD, theta_edge=theta_edge)
    phi_diff_ceabc = holonomy_phase_difference(CEABC_WORD, theta_edge=theta_edge)
    u_abcea = holonomy_unitary_phases(ABCEA_WORD, theta_edge=theta_edge)
    u_ceabc = holonomy_unitary_phases(CEABC_WORD, theta_edge=theta_edge)

    accumulated_phi_diff = phases["phi_R"] - phases["phi_L"]
    d_e = int(counts["D_E"])
    link = link_phase_difference_to_D_E(
        accumulated_phi_diff, d_e, theta_edge=theta_edge
    )

    return {
        "meta": {
            "module": "collatz_eabc_chirale_transport.py",
            "theory": THEORY,
            "theory_holonomie_stufen": THEORY_HOLONOMIE_STUFEN,
            "theory_zirkulation": THEORY_ZIRKULATION,
            "theory_epistemik": THEORY_EPISTEMIK,
            "epistemic": "Stufe-2-Modell — Helizität/Phasenkanal, kein SRT-Anspruch",
            "max_p": max_p,
            "theta_edge": theta_edge,
        },
        "helicity_map": {
            "ABCEA": {"channel": "R", "lambda": +1, "omega": +1},
            "CEABC": {"channel": "L", "lambda": -1, "omega": -1},
        },
        "single_loop": {
            "ABCEA": {
                "holonomy_phase_difference": phi_diff_abcea,
                "unitary": {
                    "phi_R": u_abcea["phi_R"],
                    "phi_L": u_abcea["phi_L"],
                    "phase_difference": u_abcea["phase_difference"],
                },
            },
            "CEABC": {
                "holonomy_phase_difference": phi_diff_ceabc,
                "unitary": {
                    "phi_R": u_ceabc["phi_R"],
                    "phi_L": u_ceabc["phi_L"],
                    "phase_difference": u_ceabc["phase_difference"],
                },
            },
        },
        "path_accumulation": phases,
        "chirality_flux": flux,
        "holonomy_counts": {
            "N_plus": counts["N_plus"],
            "N_minus": counts["N_minus"],
            "D_E": counts["D_E"],
            "S_E": counts["S_E"],
        },
        "D_E_link": link,
        "state_vector": {
            "description": "Ψ = (R, L)^T auf Faserbündel; U_E = diag(e^{iφ_R}, e^{iφ_L})",
            "phi_R": phases["phi_R"],
            "phi_L": phases["phi_L"],
            "observable_phase_difference": accumulated_phi_diff,
        },
        "boxed": {
            "helicity": "ω(γ)=±1 als diskrete Helizität, nicht nur Umlaufzahl",
            "flux": "C_E = Σω(γ) = N_R - N_L = D_E",
            "not_proper_time": "Chiraler Polarisationsraum — kein relativistisches Eigenzeit-Modell",
        },
    }


def run(
    max_p: int = 100_000,
    output: Path = DEFAULT_OUTPUT,
    *,
    theta_edge: float | None = None,
) -> dict[str, Any]:
    report = chiral_transport_report(max_p=max_p, theta_edge=theta_edge)
    output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC chiraler Transport / Helizität")
    parser.add_argument("--max-p", type=int, default=100_000)
    parser.add_argument("--theta-edge", type=float, default=None)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(max_p=args.max_p, output=args.output, theta_edge=args.theta_edge)
    phases = report["path_accumulation"]
    flux = report["chirality_flux"]
    print("=== EABC Chiraler Transport ===")
    print(f"θ_edge = {report['meta']['theta_edge']:.6f}")
    print(f"φ_R = {phases['phi_R']:.4f},  φ_L = {phases['phi_L']:.4f}")
    print(f"φ_R - φ_L = {report['state_vector']['observable_phase_difference']:.4f}")
    print(f"N_R = {flux['N_R']},  N_L = {flux['N_L']},  D_E = {flux['D_E']:+d}")
    print(
        f"ABCEA Δφ = {report['single_loop']['ABCEA']['holonomy_phase_difference']:.4f}"
    )
    print(
        f"CEABC Δφ = {report['single_loop']['CEABC']['holonomy_phase_difference']:.4f}"
    )
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
