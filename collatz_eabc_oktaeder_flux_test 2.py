#!/usr/bin/env python3
"""
EABC Oktaeder-Umgebung: numerische Validierung von Φ_E, ⟨ω_E,h⟩, ABCEA/CEABC.

Theorie: collatz_eabc_oktaeder_test.md, collatz_eabc_diskrete_geometrie.md

Abbildung:
  - C4-Äquator des regulären Oktaeders O_6 (E,A,B,C auf ±e₁, ±e₂)
  - Polare Lift-Achse ±e₃ (ABCEA→P⁺, CEABC→P⁻)
  - Oktonion-Schalen-Gewicht r_8(p) via Jacobi r₄∗r₄

Ausführung:
    python3 collatz_eabc_oktaeder_flux_test.py
    python3 collatz_eabc_oktaeder_flux_test.py --max-p 1000000
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from collatz_eabc_hodge_eabc import (
    C4_EDGE_ARROWS,
    C4_EDGE_INDEX,
    C4_EDGE_LABELS,
    C4_EDGE_NEGATIVE,
    C4_EDGE_POSITIVE,
    Phi_E,
    harmonic_form_c4,
    harmonic_holonomy_component,
    inner_product_omega_h,
    magnetic_laplacian,
    magnetic_phase_matrix,
    omega_edge_from_holonomy,
)
from collatz_eabc_holonomie_fehlerterm import holonomy_counts
from collatz_eabc_transition_graph import (
    LABELS,
    classes_from_sequence,
    prime_eabc_sequence,
    sliding_windows,
    transition_counts,
)
from collatz_eabc_wigner_field import build_w_transition_matrix

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_oktaeder_flux_test.json"
THEORY_OKTAEDER = "collatz_eabc_oktaeder_test.md"
THEORY_DISKRETE = "collatz_eabc_diskrete_geometrie.md"

EQUATORIAL_LABELS = C4_EDGE_LABELS
POLE_PLUS = "P+"
POLE_MINUS = "P-"
OCTA_VERTICES = (*LABELS, POLE_PLUS, POLE_MINUS)
VERTEX_INDEX = {v: i for i, v in enumerate(OCTA_VERTICES)}


@dataclass(frozen=True)
class OctahedronGraph:
    """Regulärer Oktaeder: 6 Knoten, 12 Kanten, C4-Äquator + polare Diagonalen."""

    vertices: tuple[str, ...]
    edges: tuple[tuple[str, str, str], ...]  # (src, dst, label)
    equatorial_edge_indices: tuple[int, ...]
    polar_plus_indices: tuple[int, ...]
    polar_minus_indices: tuple[int, ...]


def _pair_label(pair: tuple[str, str]) -> str:
    return pair[0] + pair[1]


def build_octahedron_graph() -> OctahedronGraph:
    """12 Kanten: 4 äquatoriale (C4) + 4 zu P⁺ + 4 zu P⁻."""
    edges: list[tuple[str, str, str]] = []
    eq_idx: list[int] = []
    plus_idx: list[int] = []
    minus_idx: list[int] = []

    for src, dst in C4_EDGE_ARROWS:
        label = _pair_label((src, dst))
        eq_idx.append(len(edges))
        edges.append((src, dst, label))

    for v in LABELS:
        plus_idx.append(len(edges))
        edges.append((v, POLE_PLUS, f"{v}P+"))
        minus_idx.append(len(edges))
        edges.append((v, POLE_MINUS, f"{v}P-"))

    return OctahedronGraph(
        vertices=OCTA_VERTICES,
        edges=tuple(edges),
        equatorial_edge_indices=tuple(eq_idx),
        polar_plus_indices=tuple(plus_idx),
        polar_minus_indices=tuple(minus_idx),
    )


def r4_jacobi(n: int) -> int:
    """r₄(n): Koeffizient von q^n in θ₃(q)^4 (Hardy–Wright)."""
    if n < 0:
        return 0
    if n == 0:
        return 1
    total = 0
    d = 1
    while d * d <= n:
        if n % d == 0:
            if d % 4 != 0:
                total += d
            d2 = n // d
            if d2 != d and d2 % 4 != 0:
                total += d2
        d += 1
    return 8 * total


def precompute_r8(max_n: int) -> np.ndarray:
    """r₈(n) = (r₄ ∗ r₄)(n) — Oktonionische Schalengröße |Σ_n^(8)| (Z^8-Stub)."""
    r4 = np.array([r4_jacobi(n) for n in range(max_n + 1)], dtype=float)
    conv = np.convolve(r4, r4)
    return conv[: max_n + 1]


def verify_r8_oeis() -> dict[str, Any]:
    """Sanity: erste Terme von A000118."""
    r8 = precompute_r8(5)
    expected = {0: 1.0, 1: 16.0, 2: 112.0, 3: 448.0, 4: 1136.0}
    ok = all(abs(r8[n] - expected[n]) < 1e-9 for n in expected)
    return {"expected": expected, "computed": {n: float(r8[n]) for n in expected}, "ok": ok}


def octahedron_incidence_matrix(graph: OctahedronGraph) -> np.ndarray:
    """Orientierte Inzidenz B ∈ R^{6×12}."""
    n_v = len(graph.vertices)
    n_e = len(graph.edges)
    b = np.zeros((n_v, n_e), dtype=float)
    for j, (src, dst, _label) in enumerate(graph.edges):
        b[VERTEX_INDEX[src], j] = -1.0
        b[VERTEX_INDEX[dst], j] = 1.0
    return b


def equatorial_harmonic_on_octahedron(graph: OctahedronGraph) -> np.ndarray:
    """h auf 12 Kanten: +1/2 auf äquatorialen Vorwärtskanten, 0 auf Pol-Kanten, normiert."""
    h = np.zeros(len(graph.edges), dtype=float)
    scale = 1.0 / math.sqrt(len(EQUATORIAL_LABELS))
    for j in graph.equatorial_edge_indices:
        h[j] = scale
    return h


def omega_octahedron_from_holonomy(
    classes: list[str],
    primes: list[int],
    graph: OctahedronGraph,
    r8: np.ndarray,
    *,
    shell_weighted: bool = False,
) -> dict[str, Any]:
    """
    Priminduzierte 1-Form auf O_6.

    Äquatorial: identisch zu ω_E auf C4.
    Polar: ABCEA (+1) → P⁺-Diagonalen, CEABC (−1) → P⁻-Diagonalen.
    Optional: Gewicht r₈(p) am Fenster-End-Prim p_{i+4}.
    """
    omega = np.zeros(len(graph.edges), dtype=float)
    n_plus = n_minus = 0
    shell_num = shell_den = 0.0
    polar_plus_flux = polar_minus_flux = 0.0

    windows = sliding_windows(classes, width=5)
    for w in windows:
        sign = int(w["omega"])
        if sign == 0:
            continue
        idx = w["index"]
        end_prime = primes[idx + 4] if idx + 4 < len(primes) else primes[-1]
        weight = float(r8[end_prime]) if shell_weighted and end_prime < len(r8) else 1.0

        if sign == 1:
            n_plus += 1
            for j in graph.polar_plus_indices:
                omega[j] += weight
            polar_plus_flux += weight
        else:
            n_minus += 1
            for j in graph.polar_minus_indices:
                omega[j] += weight
            polar_minus_flux += weight

        word = w["word"]
        for k in range(len(word) - 1):
            src, dst = word[k], word[k + 1]
            pair = (src, dst)
            delta = sign * weight
            if pair in C4_EDGE_ARROWS:
                omega[graph.equatorial_edge_indices[C4_EDGE_INDEX[_pair_label(pair)]]] += delta
            elif (dst, src) in C4_EDGE_ARROWS:
                rev = _pair_label((dst, src))
                omega[graph.equatorial_edge_indices[C4_EDGE_INDEX[rev]]] -= delta

        if shell_weighted:
            shell_num += sign * weight
            shell_den += weight

    eq_omega = np.array(
        [omega[j] for j in graph.equatorial_edge_indices], dtype=float
    )
    c_e = n_plus - n_minus
    s_e = n_plus + n_minus
    phi_eq = c_e / s_e if s_e > 0 else 0.0
    phi_shell = shell_num / shell_den if shell_den > 0 else 0.0

    h_eq = harmonic_form_c4()
    inner_eq = inner_product_omega_h(eq_omega, h_eq)
    h_oct = equatorial_harmonic_on_octahedron(graph)
    inner_oct = float(np.dot(omega, h_oct))

    return {
        "omega_equatorial": {EQUATORIAL_LABELS[i]: float(eq_omega[i]) for i in range(4)},
        "omega_vector": omega.tolist(),
        "N_plus": n_plus,
        "N_minus": n_minus,
        "C_E": c_e,
        "S_E": s_e,
        "Phi_equatorial": phi_eq,
        "Phi_shell_weighted": phi_shell if shell_weighted else None,
        "shell_weighted": shell_weighted,
        "polar_plus_flux": polar_plus_flux,
        "polar_minus_flux": polar_minus_flux,
        "polar_preference": (
            (polar_plus_flux - polar_minus_flux) / (polar_plus_flux + polar_minus_flux)
            if (polar_plus_flux + polar_minus_flux) > 0
            else 0.0
        ),
        "inner_product_equatorial": inner_eq,
        "inner_product_octahedron": inner_oct,
    }


def octahedron_magnetic_laplacian(
    classes: list[str],
    graph: OctahedronGraph,
) -> dict[str, Any]:
    """Magnetischer Laplace auf 6 Oktaeder-Knoten via W-Matrix + chirale Phasen."""
    w_trans = build_w_transition_matrix(classes)
    w4 = np.array(w_trans["matrix"], dtype=float)
    phases4 = magnetic_phase_matrix(classes)

    n = len(graph.vertices)
    adj = np.zeros((n, n), dtype=float)
    phases = np.zeros((n, n), dtype=float)

    for i, a in enumerate(LABELS):
        for j, b in enumerate(LABELS):
            adj[i, j] = w4[i, j]
            phases[i, j] = phases4[i, j]

    mag = magnetic_laplacian(adj, phases, hermitian=True)
    near_tol = 1e-6
    near = [v for v in mag["eigenvalues"] if abs(v) <= near_tol]
    return {
        "vertices": list(graph.vertices),
        "eigenvalues": mag["eigenvalues"],
        "smallest_lambda": mag["smallest_lambda"],
        "near_zero_count": len(near),
        "near_zero_eigenvalues": near,
        "formula": "L_mag on equatorial EABC subgraph embedded in O_6",
    }


def compare_c4_vs_octahedron(max_p: int) -> dict[str, Any]:
    """Vergleich 1D-C4 mit Oktaeder-Umgebung."""
    graph = build_octahedron_graph()
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    primes = [row["p"] for row in seq]
    max_shell = max(primes) if primes else max_p
    r8 = precompute_r8(max_shell)

    c4_phi = Phi_E(max_p)
    c4_harm = harmonic_holonomy_component(classes)
    c4_omega = omega_edge_from_holonomy(classes)

    oct_unweighted = omega_octahedron_from_holonomy(
        classes, primes, graph, r8, shell_weighted=False
    )
    oct_weighted = omega_octahedron_from_holonomy(
        classes, primes, graph, r8, shell_weighted=True
    )
    mag = octahedron_magnetic_laplacian(classes, graph)
    hol = holonomy_counts(max_p)

    phi_c4 = c4_phi["Phi_E"]
    phi_oct_eq = oct_unweighted["Phi_equatorial"]
    phi_oct_shell = oct_weighted["Phi_shell_weighted"]

    supports_phi = (
        abs(phi_c4) > 1e-12
        and abs(phi_oct_eq) > 1e-12
        and math.copysign(1.0, phi_c4) == math.copysign(1.0, phi_oct_eq)
    )
    shell_same_sign = (
        phi_oct_shell is not None
        and abs(phi_oct_shell) > 1e-12
        and math.copysign(1.0, phi_c4) == math.copysign(1.0, phi_oct_shell)
    )

    inner_c4 = c4_harm["inner_product_omega_h"]
    inner_oct = oct_unweighted["inner_product_octahedron"]
    supports_harmonic = inner_c4 != 0.0 and inner_oct != 0.0

    return {
        "meta": {
            "module": "collatz_eabc_oktaeder_flux_test.py",
            "theory_oktaeder": THEORY_OKTAEDER,
            "theory_diskrete_geometrie": THEORY_DISKRETE,
            "max_p": max_p,
            "epistemic_label": "Modellabbildung / Experiment — kein Beweis",
            "mapping": {
                "carrier": "regular octahedron O_6 in R^3",
                "equator": "C4 cycle E→A→B→C→E",
                "poles": "P+ (ABCEA), P- (CEABC)",
                "shell_weight": "r_8(p) via Jacobi r4*r4 convolution",
                "8D_projection": "(e1,e2)→E, (e3,e4)→A, (e5,e6)→B, (e7,e8)→C",
            },
        },
        "graph": {
            "vertices": list(graph.vertices),
            "edge_count": len(graph.edges),
            "equatorial_edges": list(EQUATORIAL_LABELS),
            "E_plus": sorted(C4_EDGE_POSITIVE),
            "E_minus": sorted(C4_EDGE_NEGATIVE),
        },
        "r8_verification": verify_r8_oeis(),
        "c4_baseline": {
            "Phi_E": phi_c4,
            "C_E": c4_phi["C_E"],
            "S_E": c4_phi["S_E"],
            "N_plus": c4_phi["N_plus"],
            "N_minus": c4_phi["N_minus"],
            "inner_product_omega_h": inner_c4,
            "omega_E": c4_omega["omega_E"],
            "holonomy_S_E": hol["S_E"],
        },
        "octahedron_unweighted": oct_unweighted,
        "octahedron_shell_weighted": oct_weighted,
        "magnetic_laplacian": mag,
        "comparison": {
            "Phi_equatorial_matches_c4": abs(phi_oct_eq - phi_c4) < 1e-12,
            "Phi_equatorial_delta": phi_oct_eq - phi_c4,
            "C_E_matches": oct_unweighted["C_E"] == c4_phi["C_E"],
            "inner_product_equatorial_matches_c4": abs(
                oct_unweighted["inner_product_equatorial"] - inner_c4
            )
            < 1e-9,
            "polar_orientation_preference": oct_unweighted["polar_preference"],
            "ABCEA_vs_CEABC_pole_ratio": {
                "P_plus_flux": oct_unweighted["polar_plus_flux"],
                "P_minus_flux": oct_unweighted["polar_minus_flux"],
            },
        },
        "verdict": {
            "Phi_E_nonzero_at_scale": abs(phi_c4) > 1e-12,
            "Phi_oct_equatorial_nonzero": abs(phi_oct_eq) > 1e-12,
            "Phi_oct_shell_nonzero": (
                abs(phi_oct_shell) > 1e-12 if phi_oct_shell is not None else None
            ),
            "Phi_oct_shell_same_sign_as_c4": shell_same_sign,
            "supports_Phi_E_conjecture": supports_phi and supports_harmonic,
            "supports_harmonic_pairing": supports_harmonic,
            "interpretation": (
                "Oktaeder-Äquator + harmonische Paarung unterstützen Φ_E ≠ 0"
                + (
                    "; r₈-Schalen-Gewicht verstärkt gleiches Vorzeichen"
                    if shell_same_sign
                    else "; r₈-Schalen-Gewicht kann bei endlichem X Vorzeichen drehen"
                )
                if supports_phi and supports_harmonic
                else (
                    "Oktaeder-Äquator konsistent, harmonische Paarung schwach"
                    if abs(phi_oct_eq - phi_c4) < 1e-12
                    else "Oktaeder-Test widerspricht C4-Baseline — Implementierung prüfen"
                )
            ),
        },
        "flux_density_series": _flux_series(max_p),
    }


def _flux_series(max_p: int) -> list[dict[str, Any]]:
    limits = [x for x in (1_000, 10_000, 100_000, max_p) if x <= max_p]
    if max_p not in limits:
        limits.append(max_p)
    graph = build_octahedron_graph()
    rows: list[dict[str, Any]] = []
    for lim in limits:
        seq = prime_eabc_sequence(lim)
        classes = classes_from_sequence(seq)
        primes = [row["p"] for row in seq]
        r8 = precompute_r8(max(primes) if primes else lim)
        oct_u = omega_octahedron_from_holonomy(
            classes, primes, graph, r8, shell_weighted=False
        )
        oct_w = omega_octahedron_from_holonomy(
            classes, primes, graph, r8, shell_weighted=True
        )
        c4 = Phi_E(lim)
        rows.append(
            {
                "max_p": lim,
                "Phi_E_c4": c4["Phi_E"],
                "Phi_oct_equatorial": oct_u["Phi_equatorial"],
                "Phi_oct_shell": oct_w["Phi_shell_weighted"],
                "C_E": c4["C_E"],
                "S_E": c4["S_E"],
            }
        )
    return rows


def run(max_p: int = 100_000, output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    report = compare_c4_vs_octahedron(max_p)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC Oktaeder-Flux-Test")
    parser.add_argument("--max-p", type=int, default=100_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(max_p=args.max_p, output=args.output)
    v = report["verdict"]
    c4 = report["c4_baseline"]
    oct_u = report["octahedron_unweighted"]
    print("=== EABC Oktaeder-Flux-Test ===")
    print(f"max_p = {args.max_p}")
    print(f"Φ_E (C4)           = {c4['Phi_E']:+.8f}")
    print(f"Φ_oct (äquatorial) = {oct_u['Phi_equatorial']:+.8f}")
    shell = report["octahedron_shell_weighted"].get("Phi_shell_weighted")
    if shell is not None:
        print(f"Φ_oct (r₈-gewichtet)= {shell:+.8f}")
    print(f"⟨ω,h⟩ C4           = {c4['inner_product_omega_h']:+.4f}")
    print(f"⟨ω,h⟩ Oktaeder     = {oct_u['inner_product_octahedron']:+.4f}")
    print(f"Pol-Präferenz      = {oct_u['polar_preference']:+.8f}")
    print(f"Urteil: {v['interpretation']}")
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
