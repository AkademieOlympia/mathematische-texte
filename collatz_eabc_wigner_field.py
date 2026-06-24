#!/usr/bin/env python3
"""
EABC Wigner-Feld-Analogie: signiertes Informationsfeld auf {E,A,B,C}.

Theorie: collatz_eabc_wigner_analog.md
Potential/Berry: collatz_eabc_potential_geometrie.md

  W_E(N) = #ABCE - #CEAB          (4-Pfad-Orientierung, Pfad ≠ Holonomie)
  D_E(N) = #ABCEA - #CEABC        (5-Zyklus, Zirkulationshypothese)

  W(a,b) = Σ_n χ_a(n) χ_b(n) Q(n) — signierte 4×4-Korrelationsmatrix
  Q(n) = Ω_Pfad(P_n^(4)) ∈ {+1,-1,0};  χ_c = Klassenindikator im Fenster

Ausführung:
    python3 collatz_eabc_wigner_field.py
    python3 collatz_eabc_wigner_field.py --max-p 100000
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from collatz_eabc_graph_laplacian import (
    adjacency_matrix,
    eigenvalues_sorted,
    laplacian_symmetrized,
)
from collatz_eabc_holonomie_fehlerterm import holonomy_counts
from collatz_eabc_transition_graph import (
    ABCE_WORD,
    CEAB_WORD,
    IDX,
    LABELS,
    chi_pfad_sliding,
    classes_from_sequence,
    omega_pfad,
    prime_eabc_sequence,
    sliding_windows,
    transition_counts,
)
from eabc_from_lean import Chirality, class_of, is_prime_quadruplet, q

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_wigner_field.json"
THEORY = "collatz_eabc_wigner_analog.md"
THEORY_CHIRAL = "collatz_eabc_chirale_polarisation.md"
THEORY_ZIRKULATION = "collatz_eabc_zirkulationshypothese.md"
THEORY_HOLONOMIE_STUFEN = "collatz_eabc_holonomie_stufen.md"
THEORY_POTENTIAL_GEOMETRIE = "collatz_eabc_potential_geometrie.md"
THEORY_GENERALANGRIFF = "collatz_generalangriff_2026.md"

NEAR_ZERO_TOL = 1e-6


def enumerate_quadruplets(limit: int) -> list[dict[str, Any]]:
    """Arithmetische Prim-Vierlinge (leichtgewichtig, ohne holonomie_test-Import)."""
    rows: list[dict[str, Any]] = []
    for p in range(5, limit - 7, 2):
        if p % 12 not in (5, 11):
            continue
        if not is_prime_quadruplet(p):
            continue
        chi = Chirality.ABCE if p % 12 == 5 else Chirality.CEAB
        classes = [c.value for c in (class_of(n) for n in q(p)) if c is not None]
        rows.append(
            {
                "p": p,
                "word": "".join(classes),
                "chirality": chi.value,
                "omega": 1 if chi is Chirality.ABCE else -1,
                "signature": classes,
            }
        )
    return rows


def class_indicator_vector(window_classes: list[str]) -> list[int]:
    """χ_c(n): Anwesenheit der Klasse c im 4-Fenster (0/1 je Position)."""
    vec = [0, 0, 0, 0]
    for c in window_classes:
        vec[IDX[c]] += 1
    return vec


def quadruplet_indicator(word: str) -> int:
    """Q(n) = Ω_Pfad: +1 ABCE, -1 CEAB, 0 sonst."""
    return int(omega_pfad(word))


def wigner_correlation_entry(
    chi_a: list[int],
    chi_b: list[int],
    q: int,
) -> int:
    """Einzelbeitrag χ_a · χ_b · Q — diskrete Korrelation."""
    return q * sum(x * y for x, y in zip(chi_a, chi_b))


def build_w_transition_matrix(classes: list[str]) -> dict[str, Any]:
    """
    Übergangsgeometrie W_ij = Σ_n χ_i(n) χ_j(n+1).

    Korrelations-/Übergangsgeometrie — nicht literale Wigner-Matrix.
    """
    mat = np.zeros((4, 4), dtype=float)
    for n in range(len(classes) - 1):
        i = IDX[classes[n]]
        j = IDX[classes[n + 1]]
        mat[i, j] += 1.0
    return {
        "matrix": mat.tolist(),
        "matrix_labeled": {
            LABELS[i]: {LABELS[j]: float(mat[i, j]) for j in range(4)} for i in range(4)
        },
        "transition_count": int(mat.sum()),
        "carrier": "consecutive_class_transitions",
        "formula": "W_ij = sum_n chi_i(n) chi_j(n+1)",
    }


def w_e_edge_pair_field(classes: list[str]) -> dict[str, Any]:
    """
    W_E(i,j;N) = (N_ij^(+) - N_ij^(-)) / (N_ij^(+) + N_ij^(-))

    Chirale Quasi-Wahrscheinlichkeit pro Kantenpaar i→j aus 4-Fenstern ABCE/CEAB.
    """
    n_plus: dict[tuple[str, str], int] = {(a, b): 0 for a in LABELS for b in LABELS}
    n_minus: dict[tuple[str, str], int] = {(a, b): 0 for a in LABELS for b in LABELS}
    windows = sliding_windows(classes, width=4)
    for w in windows:
        omega = int(w["omega"])
        if omega == 0:
            continue
        word = w["word"]
        bucket = n_plus if omega == 1 else n_minus
        for k in range(len(word) - 1):
            bucket[(word[k], word[k + 1])] += 1

    field: dict[str, dict[str, float | None]] = {}
    counts_labeled: dict[str, dict[str, dict[str, int]]] = {}
    for a in LABELS:
        field[a] = {}
        counts_labeled[a] = {}
        for b in LABELS:
            np_ = n_plus[(a, b)]
            nm_ = n_minus[(a, b)]
            total = np_ + nm_
            counts_labeled[a][b] = {"N_plus": np_, "N_minus": nm_, "total": total}
            field[a][b] = (np_ - nm_) / total if total > 0 else None

    return {
        "W_E_edge_field": field,
        "edge_counts": counts_labeled,
        "formula": "W_E(i,j;N) = (N_ij^(+) - N_ij^(-)) / (N_ij^(+) + N_ij^(-))",
        "oriented_words": {"positive": ABCE_WORD, "negative": CEAB_WORD},
    }


def marginal_reconstruction_w_e_edge(
    edge_field: dict[str, Any],
    transition: list[list[int]],
    global_s_w: float,
) -> dict[str, Any]:
    """
    Rekonstruiere W_E(i,j) nur aus Marginalen: globales S_W + T_ij.

    Naive Hypothese: alle Kanten tragen dieselbe globale Chiralität S_W.
    """
    total_trans = sum(sum(row) for row in transition)
    reconstructed: dict[str, dict[str, float | None]] = {}
    for ia, a in enumerate(LABELS):
        reconstructed[a] = {}
        for ib, b in enumerate(LABELS):
            t_ij = transition[ia][ib]
            reconstructed[a][b] = global_s_w if t_ij > 0 else None

    actual = edge_field["W_E_edge_field"]
    sq_err = 0.0
    compared = 0
    max_abs_err = 0.0
    for a in LABELS:
        for b in LABELS:
            act = actual[a][b]
            rec = reconstructed[a][b]
            if act is None or rec is None:
                continue
            err = act - rec
            sq_err += err * err
            max_abs_err = max(max_abs_err, abs(err))
            compared += 1

    rmse = (sq_err / compared) ** 0.5 if compared > 0 else 0.0
    recoverable = compared > 0 and rmse < 1e-12 and max_abs_err < 1e-12

    return {
        "reconstructed_from_marginals": reconstructed,
        "marginal_inputs": {
            "global_S_W": global_s_w,
            "total_transitions": total_trans,
            "note": "prediction: W_E(i,j) ≈ S_W wherever T_ij > 0",
        },
        "compared_edge_pairs": compared,
        "rmse_vs_marginals": rmse,
        "max_abs_error": max_abs_err,
        "recoverable_from_marginals_only": recoverable,
        "information_excess": not recoverable,
        "interpretation": (
            "arithmetische Wigner-Negativität: sign structure not in marginals alone"
            if not recoverable
            else "degenerate: all edges share global chirality"
        ),
    }


def information_excess_test(classes: list[str], global_s_w: float) -> dict[str, Any]:
    """Test: kann W_E(i,j) allein aus Marginalen rekonstruiert werden?"""
    edge = w_e_edge_pair_field(classes)
    transition = transition_counts(classes)
    marginal = marginal_reconstruction_w_e_edge(edge, transition, global_s_w)
    return {
        "edge_field": edge,
        "marginal_reconstruction": marginal,
        "hypothesis": "arithmetische Wigner-Negativität",
        "epistemic_label": "Hypothese + Experiment",
    }


def build_w_matrix_from_windows(
    windows: list[dict[str, Any]],
    *,
    carrier: str = "sliding_4_paths",
) -> dict[str, Any]:
    """4×4-Matrix W über {E,A,B,C} aus Gleitfenstern / Vierlingspfaden."""
    mat = np.zeros((4, 4), dtype=float)
    oriented = 0
    for w in windows:
        q = int(w["omega"])
        if q == 0:
            continue
        oriented += 1
        word = w["word"]
        chi = {c: 0 for c in LABELS}
        for ch in word:
            chi[ch] += 1
        for ia, a in enumerate(LABELS):
            for ib, b in enumerate(LABELS):
                mat[ia, ib] += q * chi[a] * chi[b]
    return {
        "matrix": mat.tolist(),
        "matrix_labeled": {
            LABELS[i]: {LABELS[j]: float(mat[i, j]) for j in range(4)} for i in range(4)
        },
        "oriented_window_count": oriented,
        "carrier": carrier,
        "formula": "W(a,b) = sum_n chi_a(n) chi_b(n) Q(n)",
    }


def build_w_matrix_from_quadruplets(max_p: int) -> dict[str, Any]:
    """W-Matrix auf arithmetischen Prim-Vierlingen (χ_E^quad-Träger)."""
    mat = np.zeros((4, 4), dtype=float)
    quads = enumerate_quadruplets(max_p)
    for qd in quads:
        q = int(qd["omega"])
        classes = qd["signature"]
        chi = {c: 0 for c in LABELS}
        for ch in classes:
            chi[ch] += 1
        for ia, a in enumerate(LABELS):
            for ib, b in enumerate(LABELS):
                mat[ia, ib] += q * chi[a] * chi[b]
    return {
        "matrix": mat.tolist(),
        "matrix_labeled": {
            LABELS[i]: {LABELS[j]: float(mat[i, j]) for j in range(4)} for i in range(4)
        },
        "quadruplet_count": len(quads),
        "carrier": "prime_quadruplets",
        "formula": "W(a,b) = sum_{Q in Vierlingen} chi_a chi_b omega(Q)",
    }


def w_e_counts(max_p: int) -> dict[str, Any]:
    """W_E(N) = #ABCE - #CEAB auf 4-Pfad-Gleitfenstern."""
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    pfad = chi_pfad_sliding(classes)
    n_abce = pfad["abce_windows"]
    n_ceab = pfad["ceab_windows"]
    w_e = n_abce - n_ceab
    total = n_abce + n_ceab
    return {
        "max_p": max_p,
        "N_ABCE": n_abce,
        "N_CEAB": n_ceab,
        "W_E": w_e,
        "S_W": w_e / total if total > 0 else 0.0,
        "block_width": 4,
        "words": {"positive": ABCE_WORD, "negative": CEAB_WORD},
    }


def w_e_profile(classes: list[str]) -> list[dict[str, Any]]:
    """Kumulatives W_E(n) entlang der Primfolge (laufende 4-Pfad-Bilanz)."""
    windows = sliding_windows(classes, width=4)
    cum_abce = 0
    cum_ceab = 0
    profile: list[dict[str, Any]] = []
    for w in windows:
        if w["omega"] == 1:
            cum_abce += 1
        elif w["omega"] == -1:
            cum_ceab += 1
        profile.append(
            {
                "index": w["index"],
                "word": w["word"],
                "omega": w["omega"],
                "cum_ABCE": cum_abce,
                "cum_CEAB": cum_ceab,
                "W_E_cumulative": cum_abce - cum_ceab,
            }
        )
    return profile


def sign_domain_analysis(w_e_series: list[int], *, tol: float = 0.0) -> dict[str, Any]:
    """
    Vorzeichendomänen (+), (-), (0) für eine W_E-Zeitreihe.

    Region = maximales zusammenhängendes Intervall konstanten Vorzeichens.
    """
    if not w_e_series:
        return {
            "positive_regions": 0,
            "negative_regions": 0,
            "null_regions": 0,
            "positive_steps": 0,
            "negative_steps": 0,
            "null_steps": 0,
            "sign_changes": 0,
            "dominant_sign": "null",
        }

    def sign(v: int) -> str:
        if v > tol:
            return "positive"
        if v < -tol:
            return "negative"
        return "null"

    signs = [sign(v) for v in w_e_series]
    pos_r = neg_r = null_r = 0
    pos_s = signs.count("positive")
    neg_s = signs.count("negative")
    null_s = signs.count("null")
    changes = 0
    for i in range(1, len(signs)):
        if signs[i] != signs[i - 1]:
            changes += 1
    prev = signs[0]
    for s in signs[1:]:
        if s != prev:
            if prev == "positive":
                pos_r += 1
            elif prev == "negative":
                neg_r += 1
            else:
                null_r += 1
            prev = s
    if prev == "positive":
        pos_r += 1
    elif prev == "negative":
        neg_r += 1
    else:
        null_r += 1

    totals = {"positive": pos_s, "negative": neg_s, "null": null_s}
    dominant = max(totals, key=totals.get)  # type: ignore[arg-type]

    return {
        "positive_regions": pos_r,
        "negative_regions": neg_r,
        "null_regions": null_r,
        "positive_steps": pos_s,
        "negative_steps": neg_s,
        "null_steps": null_s,
        "sign_changes": changes,
        "dominant_sign": dominant,
        "info_in_transitions": changes > 0,
    }


def near_zero_eigenmode_stub(max_p: int, *, tol: float = NEAR_ZERO_TOL) -> dict[str, Any]:
    """
    Near-zero-Moden von L_E^sym — Wigner-Interferenzregion-Analog.

    Stub: nutzt collatz_eabc_graph_laplacian; kein Dirac-Operator im Repo.
    """
    adj, meta = adjacency_matrix(max_p)
    l_sym = laplacian_symmetrized(adj)
    eigvals = eigenvalues_sorted(l_sym)
    eigvecs = np.linalg.eigh(l_sym)[1]
    near_idx = [i for i, v in enumerate(eigvals) if abs(v) <= tol]
    modes: list[dict[str, Any]] = []
    for i in near_idx:
        vec = eigvecs[:, i]
        modes.append(
            {
                "index": i,
                "eigenvalue": float(eigvals[i]),
                "eigenvector": {LABELS[k]: float(vec[k]) for k in range(4)},
                "norm": float(np.linalg.norm(vec)),
            }
        )
    return {
        "operator": "L_E^sym (graph Laplacian stub; no Dirac D in repo)",
        "reference": "collatz_eabc_graph_laplacian.py",
        "dirac_stub": "D psi = lambda psi — not implemented; use L_E near-zero modes",
        "eigenvalues": [float(v) for v in eigvals],
        "near_zero_tolerance": tol,
        "near_zero_mode_count": len(modes),
        "near_zero_modes": modes,
        "pop_mass_analog": "sum of |psi_k|^2 for |lambda_k| < tol",
        "graph_meta": {"max_p": meta["max_p"], "total_edges": meta["total_edges"]},
    }


def wigner_field_report(max_p: int = 100_000) -> dict[str, Any]:
    """Vollständiger Wigner-Feld-Bericht."""
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    windows4 = sliding_windows(classes, width=4)
    w_sliding = build_w_matrix_from_windows(windows4, carrier="prime_sequence_sliding_4")
    w_transition = build_w_transition_matrix(classes)
    w_quad = build_w_matrix_from_quadruplets(max_p)
    w_e = w_e_counts(max_p)
    hol = holonomy_counts(max_p)
    edge_field = w_e_edge_pair_field(classes)
    info_excess = information_excess_test(classes, w_e["S_W"])
    profile = w_e_profile(classes)
    cum_series = [p["W_E_cumulative"] for p in profile]
    domains = sign_domain_analysis(cum_series)
    spectrum_stub = near_zero_eigenmode_stub(max_p)

    w_mat = np.array(w_sliding["matrix"])
    sign_structure = {
        LABELS[i]: {
            LABELS[j]: ("+" if w_mat[i, j] > 0 else "-" if w_mat[i, j] < 0 else "0")
            for j in range(4)
        }
        for i in range(4)
    }

    return {
        "meta": {
            "module": "collatz_eabc_wigner_field.py",
            "theory": THEORY,
            "theory_chirale_polarisation": THEORY_CHIRAL,
            "theory_zirkulation": THEORY_ZIRKULATION,
            "theory_holonomie_stufen": THEORY_HOLONOMIE_STUFEN,
            "theory_potential_geometrie": THEORY_POTENTIAL_GEOMETRIE,
            "theory_generalangriff": THEORY_GENERALANGRIFF,
            "epistemic": "Analogie/Hypothese/Modell — keine Quantenphysik-Behauptung",
            "max_p": max_p,
        },
        "definitions": {
            "W_E_4block": "count(ABCE) - count(CEAB) on 4-path windows",
            "D_E_5block": "count(ABCEA) - count(CEABC) on 5-cycle holonomy",
            "W_transition_ij": "W_ij = sum_n chi_i(n) chi_j(n+1)",
            "W_matrix_signed": "W(a,b) = sum_n chi_a(n) chi_b(n) Q(n)",
            "W_E_edge": "W_E(i,j;N) = (N_ij^(+) - N_ij^(-)) / (N_ij^(+) + N_ij^(-))",
            "Q": "Omega_Pfad in {+1,-1,0}",
        },
        "W_E": w_e,
        "D_E": {
            "N_plus": hol["N_plus"],
            "N_minus": hol["N_minus"],
            "D_E": hol["D_E"],
            "S_E": hol["S_E"],
            "block_width": 5,
        },
        "four_vs_five_block": {
            "W_E": w_e["W_E"],
            "D_E": hol["D_E"],
            "same_sign": (w_e["W_E"] >= 0) == (hol["D_E"] >= 0),
            "difference": w_e["W_E"] - hol["D_E"],
            "note": "4-block Pfadorientierung (CEAB) vs. 5-block Holonomie (CEABC)",
        },
        "W_transition_matrix": w_transition,
        "W_matrix_sliding": w_sliding,
        "W_matrix_quadruplets": w_quad,
        "W_E_edge_field": edge_field,
        "information_excess_test": info_excess["marginal_reconstruction"],
        "information_excess_hypothesis": info_excess["hypothesis"],
        "sign_structure": sign_structure,
        "sign_domains": domains,
        "W_E_profile_tail": profile[-5:] if len(profile) >= 5 else profile,
        "W_E_profile_length": len(profile),
        "near_zero_modes": spectrum_stub,
        "boxed": {
            "information": "signed field on EABC — not point counts",
            "wigner_analog": "quasi-probability: negative domains allowed",
            "state_space": "transition structure E<->A<->B<->C, not single primes",
        },
    }


def run(max_p: int = 100_000, output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    report = wigner_field_report(max_p=max_p)
    output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC Wigner-Feld-Analogie")
    parser.add_argument("--max-p", type=int, default=100_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(max_p=args.max_p, output=args.output)
    w_e = report["W_E"]
    d_e = report["D_E"]
    w_trans = report["W_transition_matrix"]["matrix_labeled"]
    w_edge = report["W_E_edge_field"]["W_E_edge_field"]
    excess = report["information_excess_test"]
    print("=== EABC Wigner-Feld-Analogie ===")
    print(f"W_E (4-block) = {w_e['W_E']:+d}  S_W = {w_e['S_W']:+.6f}")
    print(f"D_E (5-block) = {d_e['D_E']:+d}  S_E = {d_e['S_E']:+.6f}")
    print(f"ABCE={w_e['N_ABCE']}, CEAB={w_e['N_CEAB']}")
    print("W_ij (transition geometry):")
    for row_label in LABELS:
        row = w_trans[row_label]
        vals = " ".join(f"{row[c]:8.0f}" for c in LABELS)
        print(f"  {row_label}: {vals}")
    w_ab = w_edge["A"]["B"]
    print(f"W_E(A,B;N) = {w_ab:+.6f}" if w_ab is not None else "W_E(A,B;N) = n/a")
    print(
        f"Information excess (marginals only): {excess['information_excess']} "
        f"(rmse={excess['rmse_vs_marginals']:.6f})"
    )
    nz = report["near_zero_modes"]["near_zero_mode_count"]
    print(f"Near-zero modes (L_E^sym): {nz}")
    print(f"Sign changes (W_E cumulative): {report['sign_domains']['sign_changes']}")
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
