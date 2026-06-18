#!/usr/bin/env python3
"""
EABC Übergangsraum-Geometrie: C4=S^1, Hodge-Zerlegung, magnetischer Laplace, Flussdichte.

Theorie: collatz_eabc_uebergangsraum.md, collatz_eabc_signierte_massstruktur.md

  Kantenfundament: EA, AB, BC, CE auf Zyklus E→A→B→C→E
  h ∈ H^1(S^1) kanonischer harmonischer Generator
  ⟨ω_E, h⟩ — priminduzierter Fluss entlang harmonischer Klasse
  L = D - W (reell); L_mag = D - U, U_ij = A_ij e^{iθ_ij}
  Φ_E = C_E = N_+ - N_- (diskreter AB-Fluss)
  flux_density = C_E / N_cycles = S_E

Ausführung:
    python3 collatz_eabc_hodge_eabc.py
    python3 collatz_eabc_hodge_eabc.py --max-p 1000000
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from collatz_eabc_holonomie_fehlerterm import holonomy_counts
from collatz_eabc_transition_graph import (
    ABCEA_WORD,
    CEABC_WORD,
    IDX,
    LABELS,
    classes_from_sequence,
    prime_eabc_sequence,
    sliding_windows,
    transition_counts,
)
from collatz_eabc_wigner_field import (
    build_w_transition_matrix,
    w_e_counts,
    w_e_edge_pair_field,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_hodge_eabc.json"
THEORY_UEBERGANGSRAUM = "collatz_eabc_uebergangsraum.md"
THEORY_SIGNIERTE_MASS = "collatz_eabc_signierte_massstruktur.md"
THEORY_ZIRKULATION = "collatz_eabc_zirkulationshypothese.md"

NEAR_ZERO_TOL = 1e-6

# Kanonischer C4-Zyklus E→A→B→C→E; fundamentale Kanten (nicht Knoten)
C4_EDGE_LABELS: tuple[str, ...] = ("EA", "AB", "BC", "CE")
C4_EDGE_ARROWS: tuple[tuple[str, str], ...] = (
    ("E", "A"),
    ("A", "B"),
    ("B", "C"),
    ("C", "E"),
)
C4_EDGE_INDEX = {label: i for i, label in enumerate(C4_EDGE_LABELS)}


def edge_incidence_matrix_c4() -> np.ndarray:
    """
    Orientierte Inzidenz B ∈ R^{4×4}: Zeilen = Knoten {E,A,B,C}, Spalten = Kanten EA,AB,BC,CE.

    B[v,e] = +1 wenn Kante e in v endet, −1 wenn e in v startet (E→A: B[E,EA]=−1, B[A,EA]=+1).
    """
    b = np.zeros((4, 4), dtype=float)
    for j, (src, dst) in enumerate(C4_EDGE_ARROWS):
        b[IDX[src], j] = -1.0
        b[IDX[dst], j] = 1.0
    return b


def harmonic_form_c4(*, normalized: bool = True) -> np.ndarray:
    """Kanonischer Generator h ∈ H^1(S^1) auf C4: alle Vorwärtskanten +1."""
    h = np.ones(4, dtype=float)
    if normalized:
        h /= np.linalg.norm(h)
    return h


def laplacian_from_W(
    w: np.ndarray,
    *,
    symmetrize: bool = True,
    tol: float = NEAR_ZERO_TOL,
) -> dict[str, Any]:
    """
    Graph-Laplace L = D − W auf Knotenraum {E,A,B,C}.

    D_ii = Σ_j W_ij. Optional Symmetrisierung (W+W^T)/2 vor L.
    """
    mat = np.array(w, dtype=float)
    if symmetrize:
        mat = 0.5 * (mat + mat.T)
    deg = mat.sum(axis=1)
    lap = np.diag(deg) - mat
    eigvals = np.linalg.eigvalsh(lap)
    eigvals = np.sort(eigvals)
    eigvecs = np.linalg.eigh(lap)[1]
    near_idx = [i for i, v in enumerate(eigvals) if abs(v) <= tol]
    near_modes: list[dict[str, Any]] = []
    for i in near_idx:
        vec = eigvecs[:, i]
        near_modes.append(
            {
                "index": int(i),
                "eigenvalue": float(eigvals[i]),
                "eigenvector": {LABELS[k]: float(vec[k]) for k in range(4)},
            }
        )
    return {
        "laplacian": lap.tolist(),
        "eigenvalues": [float(v) for v in eigvals],
        "smallest_lambda": float(eigvals[0]),
        "lambda_2": float(eigvals[1]) if len(eigvals) > 1 else None,
        "near_zero_tolerance": tol,
        "near_zero_eigenmodes": near_modes,
        "near_zero_count": len(near_modes),
        "formula": "L = D - W",
    }


def omega_edge_from_holonomy(classes: list[str]) -> dict[str, Any]:
    """
    Priminduzierte Kanten-1-Form ω_E auf {EA,AB,BC,CE}.

    Jede erkannte 5-Zyklus-Orientierung ABCEA (+1) bzw. CEABC (−1) trägt ±1 auf
    die vier Zykluskanten des Wortes.
    """
    omega = np.zeros(4, dtype=float)
    n_plus = n_minus = 0
    windows = sliding_windows(classes, width=5)
    for w in windows:
        sign = int(w["omega"])
        if sign == 0:
            continue
        if sign == 1:
            n_plus += 1
        else:
            n_minus += 1
        word = w["word"]
        for k in range(len(word) - 1):
            src, dst = word[k], word[k + 1]
            pair = (src, dst)
            if pair in C4_EDGE_ARROWS:
                omega[C4_EDGE_INDEX[_pair_label(pair)]] += sign
            elif (dst, src) in C4_EDGE_ARROWS:
                rev_label = _pair_label((dst, src))
                omega[C4_EDGE_INDEX[rev_label]] -= sign
    return {
        "omega_E": {C4_EDGE_LABELS[i]: float(omega[i]) for i in range(4)},
        "omega_vector": omega.tolist(),
        "N_plus_cycles": n_plus,
        "N_minus_cycles": n_minus,
        "C_E": n_plus - n_minus,
        "carrier": "holonomy_5_windows_on_C4_edges",
    }


def _pair_label(pair: tuple[str, str]) -> str:
    return pair[0] + pair[1]


def inner_product_omega_h(omega: np.ndarray, h: np.ndarray | None = None) -> float:
    """⟨ω_E, h⟩ mit kanonischem harmonischem Generator h."""
    if h is None:
        h = harmonic_form_c4(normalized=True)
    return float(np.dot(omega, h))


def discrete_hodge_decomposition(omega: np.ndarray) -> dict[str, Any]:
    """
    Diskrete Hodge-Zerlegung ω = ω_grad + ω_harm auf dem C4-Komplex (Stub).

    ω_grad = B^T φ mit φ = (BB^T)^{-1} B ω (L²-Projektion auf Bild von d).
    ω_harm = ω − ω_grad (Rest; auf C4 ≈ harmonischer Anteil).
    """
    b = edge_incidence_matrix_c4()
    omega = np.asarray(omega, dtype=float)
    bbt = b @ b.T
    phi = np.linalg.solve(bbt + 1e-12 * np.eye(4), b @ omega)
    omega_grad = b.T @ phi
    omega_harm = omega - omega_grad
    h = harmonic_form_c4(normalized=True)
    harm_norm = float(np.linalg.norm(omega_harm))
    return {
        "omega_gradient": omega_grad.tolist(),
        "omega_harmonic": omega_harm.tolist(),
        "harmonic_norm": harm_norm,
        "gradient_norm": float(np.linalg.norm(omega_grad)),
        "omega_total_norm": float(np.linalg.norm(omega)),
        "inner_product_harmonic_generator": inner_product_omega_h(omega_harm, h),
        "formula": "omega = d phi + h (coexact stub omitted on C4)",
    }


def harmonic_holonomy_component(classes: list[str]) -> dict[str, Any]:
    """Harmonischer Anteil von ω_E und sein Bezug zu C_E."""
    edge_data = omega_edge_from_holonomy(classes)
    omega = np.array(edge_data["omega_vector"], dtype=float)
    hodge = discrete_hodge_decomposition(omega)
    h = harmonic_form_c4()
    c_e = edge_data["C_E"]
    inner_full = inner_product_omega_h(omega, h)
    return {
        "C_E": c_e,
        "omega_E": edge_data["omega_E"],
        "inner_product_omega_h": inner_full,
        "harmonic_component": hodge,
        "harmonic_norm": hodge["harmonic_norm"],
        "ratio_harmonic_norm_to_C_E": (
            hodge["harmonic_norm"] / abs(c_e) if c_e != 0 else None
        ),
        "interpretation": "⟨ω_E,h⟩ measures prime-induced flux along harmonic class",
    }


def magnetic_phase_matrix(classes: list[str]) -> np.ndarray:
    """
  Orientierte Phasen θ_ij aus chiralem Kantenfeld W_E(i,j;N).

  θ_ij = (π/2) · W_E(i,j) für besetzte Kanten, sonst 0.
    """
    edge_field = w_e_edge_pair_field(classes)["W_E_edge_field"]
    phases = np.zeros((4, 4), dtype=float)
    for ia, a in enumerate(LABELS):
        for ib, b in enumerate(LABELS):
            val = edge_field[a][b]
            if val is not None:
                phases[ia, ib] = 0.5 * math.pi * float(val)
    return phases


def magnetic_laplacian(
    adj: np.ndarray,
    phases: np.ndarray,
    *,
    hermitian: bool = True,
) -> dict[str, Any]:
    """
    Magnetischer Laplace L_mag = D − U, U_ij = A_ij exp(i θ_ij).

    Optional Hermitisierung: U_sym = (U + U†)/2.
    """
    a = np.array(adj, dtype=float)
    theta = np.array(phases, dtype=float)
    u = a * np.exp(1j * theta)
    if hermitian:
        u = 0.5 * (u + u.conj().T)
        a_eff = np.abs(u)
    else:
        a_eff = a
    deg = a_eff.sum(axis=1)
    l_mag = np.diag(deg) - u
    if hermitian:
        eigvals = np.linalg.eigvalsh(l_mag)
    else:
        eigvals = np.linalg.eigvals(l_mag)
        eigvals = np.sort(np.real(eigvals))
    eigvals = np.sort(eigvals)
    return {
        "laplacian_real": np.real(l_mag).tolist(),
        "laplacian_imag": np.imag(l_mag).tolist(),
        "eigenvalues": [float(v) for v in eigvals],
        "smallest_lambda": float(eigvals[0]),
        "hermitian": hermitian,
        "formula": "L_mag = D - U, U_ij = A_ij exp(i theta_ij)",
    }


def signed_measure_graph(w_matrix: np.ndarray) -> dict[str, Any]:
    """Gerichteter gewichteter Graph G_E=(V,E,w) aus 4×4-W-Matrix."""
    w = np.array(w_matrix, dtype=float)
    edges: list[dict[str, Any]] = []
    for i, src in enumerate(LABELS):
        for j, dst in enumerate(LABELS):
            weight = float(w[i, j])
            if weight != 0.0:
                edges.append(
                    {
                        "source": src,
                        "target": dst,
                        "weight": weight,
                        "abs_weight": abs(weight),
                        "sign": "+" if weight > 0 else "-",
                    }
                )
    return {
        "vertices": list(LABELS),
        "edges": edges,
        "edge_count": len(edges),
        "total_abs_weight": float(np.abs(w).sum()),
        "carrier": "signed_W_matrix",
    }


def flux_density_limit(max_p: int) -> dict[str, Any]:
    """
    Flussdichte Φ_E / N_cycles = C_E(X) / #{erkannte Zyklen} = S_E(X).

    Zentral für Orientierungsklassen-Vermutung: lim ≠ 0 ⇔ bevorzugte Orientierung.
    """
    hol = holonomy_counts(max_p)
    n_cycles = hol["N_plus"] + hol["N_minus"]
    c_e = hol["D_E"]
    density = c_e / n_cycles if n_cycles > 0 else 0.0
    return {
        "max_p": max_p,
        "C_E": c_e,
        "N_cycles": n_cycles,
        "N_plus": hol["N_plus"],
        "N_minus": hol["N_minus"],
        "flux_density": density,
        "S_E": hol["S_E"],
        "formula": "flux_density = C_E / (N_+ + N_-) = S_E",
        "conjecture": "lim_{X→∞} flux_density ≠ 0 ⟺ arithmetische Orientierungsklasse",
    }


def flux_density_series(limits: list[int] | None = None) -> dict[str, Any]:
    """Flussdichte an mehreren Prim-Obergrenzen (10^3 … 10^6)."""
    if limits is None:
        limits = [1_000, 10_000, 100_000, 1_000_000]
    series = [flux_density_limit(x) for x in limits]
    return {
        "series": series,
        "limits": limits,
        "interpretation": "tanh(Θ_E) analog: normalized oriented flux on S^1",
    }


def orientation_information_test(max_p: int) -> dict[str, Any]:
    """
    Kann (N_+, N_−) aus (N_++N_−, S_E) rekonstruiert werden?

    Algebraisch ja; ohne S_E (nur Gesamtzahl) nein.
    """
    hol = holonomy_counts(max_p)
    w = w_e_counts(max_p)
    n_total_5 = hol["N_plus"] + hol["N_minus"]
    s_e = hol["S_E"]
    n_plus_rec = n_total_5 * (1.0 + s_e) / 2.0
    n_minus_rec = n_total_5 * (1.0 - s_e) / 2.0
    exact_5 = (
        abs(n_plus_rec - hol["N_plus"]) < 1e-9
        and abs(n_minus_rec - hol["N_minus"]) < 1e-9
    )
    n_total_4 = w["N_ABCE"] + w["N_CEAB"]
    s_w = w["S_W"]
    n_abce_rec = n_total_4 * (1.0 + s_w) / 2.0
    n_ceab_rec = n_total_4 * (1.0 - s_w) / 2.0
    exact_4 = (
        abs(n_abce_rec - w["N_ABCE"]) < 1e-9
        and abs(n_ceab_rec - w["N_CEAB"]) < 1e-9
    )
    return {
        "max_p": max_p,
        "five_block": {
            "N_plus": hol["N_plus"],
            "N_minus": hol["N_minus"],
            "N_total": n_total_5,
            "S_E": s_e,
            "reconstructed_N_plus": n_plus_rec,
            "reconstructed_N_minus": n_minus_rec,
            "recoverable_from_total_and_S": exact_5,
            "recoverable_from_total_only": False,
        },
        "four_block": {
            "N_ABCE": w["N_ABCE"],
            "N_CEAB": w["N_CEAB"],
            "N_total": n_total_4,
            "S_W": s_w,
            "reconstructed_N_ABCE": n_abce_rec,
            "reconstructed_N_CEAB": n_ceab_rec,
            "recoverable_from_total_and_S": exact_4,
        },
        "hypothesis": "arithmetische Wigner-Negativität: per-edge sign structure not in marginals",
        "epistemic_label": "Definition + Experiment",
    }


def hodge_report(max_p: int = 100_000) -> dict[str, Any]:
    """Vollständiger Übergangsraum-/Hodge-/Fluss-Bericht."""
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    w_trans = build_w_transition_matrix(classes)
    w_mat = np.array(w_trans["matrix"], dtype=float)
    lap_w = laplacian_from_W(w_mat, symmetrize=True)
    phases = magnetic_phase_matrix(classes)
    mag = magnetic_laplacian(w_mat, phases, hermitian=True)
    harm = harmonic_holonomy_component(classes)
    flux = flux_density_limit(max_p)
    orient = orientation_information_test(max_p)
    h = harmonic_form_c4()
    b = edge_incidence_matrix_c4()

    return {
        "meta": {
            "module": "collatz_eabc_hodge_eabc.py",
            "theory_uebergangsraum": THEORY_UEBERGANGSRAUM,
            "theory_signierte_mass": THEORY_SIGNIERTE_MASS,
            "theory_zirkulation": THEORY_ZIRKULATION,
            "max_p": max_p,
            "epistemic": "Geometrie/Analogie/Hypothese — keine Quantenphysik-Behauptung",
        },
        "topology": {
            "cycle": "E→A→B→C→E ≅ S^1",
            "H1": "Z",
            "fundamental_edges": list(C4_EDGE_LABELS),
            "incidence_matrix": b.tolist(),
        },
        "harmonic_form": {
            "h": {C4_EDGE_LABELS[i]: float(h[i]) for i in range(4)},
            "normalized": True,
        },
        "omega_E": harm["omega_E"],
        "inner_product_omega_h": harm["inner_product_omega_h"],
        "harmonic_holonomy": harm,
        "laplacian_from_W": lap_w,
        "magnetic_laplacian": mag,
        "signed_measure_graph": signed_measure_graph(w_mat),
        "flux_density": flux,
        "flux_density_series": flux_density_series(
            [x for x in (1_000, 10_000, 100_000, max_p) if x <= max_p]
            + ([max_p] if max_p not in (1_000, 10_000, 100_000) else [])
        ),
        "orientation_information_test": orient,
        "boxed": {
            "central_conjecture": (
                "lim_{X→∞} C_E(X)/N_cycles(X) ≠ 0 ⟺ nontrivial arithmetic orientation class"
            ),
            "flux_not_geometry": "C_E = discrete AB flux Φ_E, not local edge geometry",
            "W_E_reinterpretation": "S_W = normalized flux density ≈ tanh(Θ_E), Wilson loop not quasi-probability",
        },
    }


def run(max_p: int = 100_000, output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    report = hodge_report(max_p=max_p)
    output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC Übergangsraum Hodge/Fluss")
    parser.add_argument("--max-p", type=int, default=100_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(max_p=args.max_p, output=args.output)
    lap = report["laplacian_from_W"]
    flux = report["flux_density"]
    print("=== EABC Übergangsraum (Hodge/Fluss) ===")
    print(f"Spec(L) = {[f'{v:.4f}' for v in lap['eigenvalues']]}")
    print(f"⟨ω_E, h⟩ = {report['inner_product_omega_h']:+.4f}")
    print(f"C_E = {flux['C_E']:+d}, flux_density = {flux['flux_density']:+.8f}")
    print(f"Harmonic norm = {report['harmonic_holonomy']['harmonic_norm']:.4f}")
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
