#!/usr/bin/env python3
"""
EABC-Sagnac-Zirkulation C_E(X): diskrete 1-Form, Kantenorientierung ω, Zyklus-Ω(C).

Theorie: collatz_eabc_zirkulationshypothese.md (kanonisch), collatz_eabc_sagnac.md (Intuition)

  G_E = (V, E), V = {E, A, B, C}
  ω(e) ∈ {+1, -1, 0} auf kanonischen Zykluskanten (Lückenmuster (2,4,2,4))
  A(i→j) diskrete 1-Form auf dem 4-Zyklus A→B→C→E→A
  C_E(X) = Σ_{γ erkannt} ω(γ)  über 5-Fenster ABCEA (+1) / CEABC (-1)

  C_E(X) = N_+(X) - N_-(X) = D_E(X) = Δ_E(X)
  S_E(X) = C_E(X) / (N_+(X) + N_-(X))

Ausführung:
    python3 collatz_eabc_sagnac_circulation.py
    python3 collatz_eabc_sagnac_circulation.py --max-p 1000000
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from collatz_eabc_transition_graph import (
    ABCEA_WORD,
    CEABC_WORD,
    LABELS,
    classes_from_sequence,
    omega_hol,
    prime_eabc_sequence,
    sliding_windows,
)
from eabc_from_lean import EClass

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_sagnac_circulation.json"
THEORY_ZIRKULATION = "collatz_eabc_zirkulationshypothese.md"
THEORY_SAGNAC = "collatz_eabc_sagnac.md"

# Kanonischer 4-Zyklus A→B→C→E→A; Lücken mod 12 = (2,4,2,4)
CANONICAL_FORWARD_EDGES: tuple[tuple[str, str], ...] = (
    ("A", "B"),
    ("B", "C"),
    ("C", "E"),
    ("E", "A"),
)
CANONICAL_REVERSE_EDGES: tuple[tuple[str, str], ...] = tuple(
    (dst, src) for src, dst in CANONICAL_FORWARD_EDGES
)
CANONICAL_GAP_PATTERN = (2, 4, 2, 4)
EABC_RESIDUES = {EClass.E: 1, EClass.A: 5, EClass.B: 7, EClass.C: 11}


def edge_omega(src: str, dst: str) -> int:
    """
    Kantenorientierung ω(e) auf dem EABC-Transportgraphen.

    +1: Kante aligned mit A→B→C→E→A (kanonische Vorwärtsrichtung)
    -1: entgegengesetzte kanonische Kante
     0: keine Zykluskante
    """
    pair = (src, dst)
    if pair in CANONICAL_FORWARD_EDGES:
        return 1
    if pair in CANONICAL_REVERSE_EDGES:
        return -1
    return 0


def gap_from_classes(src: str, dst: str) -> int:
    """Primlücke (r_dst - r_src) mod 12 zwischen EABC-Klassen."""
    return (EABC_RESIDUES[EClass(dst)] - EABC_RESIDUES[EClass(src)]) % 12


def edge_omega_from_gap(src: str, dst: str) -> int:
    """
    Alternative ω-Lesart: Vorwärtskante iff Lücke in (2,4,2,4)-Zyklus passt.

    Für kanonische Zykluskanten stimmt dies mit edge_omega überein.
    """
    gap = gap_from_classes(src, dst)
    if gap not in CANONICAL_GAP_PATTERN:
        return 0
    return edge_omega(src, dst)


def discrete_one_form() -> dict[str, float]:
    """
    Diskrete 1-Form A auf Zykluskanten: A(e) = ω(e)/4.

    Dann ∮_{ABCEA} A = +1 und ∮_{CEABC} A = -∮_{ABCEA} A (Vorzeichen via ω(γ)).
    """
    form: dict[str, float] = {}
    for src, dst in CANONICAL_FORWARD_EDGES:
        key = f"{src}->{dst}"
        form[key] = edge_omega(src, dst) / 4.0
    for src, dst in CANONICAL_REVERSE_EDGES:
        key = f"{src}->{dst}"
        form[key] = edge_omega(src, dst) / 4.0
    return form


def cycle_omega_graph(word: str) -> int:
    """Ω(C) = ∏ ω(v_i, v_{i+1}) entlang eines 5-Worts (0 wenn nicht Zykluskante)."""
    if len(word) < 2:
        return 0
    prod = 1
    for i in range(len(word) - 1):
        w = edge_omega(word[i], word[i + 1])
        if w == 0:
            return 0
        prod *= w
    return prod


def line_integral_one_form(word: str) -> float:
    """
    ∮_γ A = ω(γ) · Σ |A(v_i, v_{i+1})| auf erkannten Zyklusorientierungen.

    Für ABCEA/CEABC mit nur Vorwärtskanten liefert Σ A(e)=+1; ω(γ)=±1 gibt das Vorzeichen.
    """
    if len(word) < 2:
        return 0.0
    form = discrete_one_form()
    magnitude = 0.0
    for i in range(len(word) - 1):
        key = f"{word[i]}->{word[i + 1]}"
        magnitude += abs(form.get(key, 0.0))
    orient = omega_cycle(word)
    if orient == 0:
        return magnitude if magnitude else 0.0
    return orient * magnitude


def omega_cycle(word: str) -> int:
    """ω(γ): ABCEA=+1, CEABC=-1, sonst 0 (Sagnac-Zyklusorientierung)."""
    return omega_hol(word)


def circulation_C_E(max_p: int) -> dict[str, Any]:
    """
    C_E(X) = Σ ω(γ) über erkannte 5-Zyklen bis Prim-Obergrenze X.

    Identisch mit D_E = N_+ - N_-.
    """
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    windows = sliding_windows(classes, width=5)

    c_e = 0
    n_plus = 0
    n_minus = 0
    omega_graph_sum = 0
    line_integral_plus = 0.0
    line_integral_minus = 0.0

    for w in windows:
        word = w["word"]
        om = omega_cycle(word)
        if om == 1:
            n_plus += 1
            c_e += 1
            omega_graph_sum += cycle_omega_graph(word)
            line_integral_plus += line_integral_one_form(word)
        elif om == -1:
            n_minus += 1
            c_e -= 1
            omega_graph_sum += cycle_omega_graph(word)
            line_integral_minus += line_integral_one_form(word)

    total = n_plus + n_minus
    delta_e = n_plus - n_minus
    s_e = delta_e / total if total > 0 else 0.0
    d_tilde = delta_e / math.sqrt(total) if total > 0 else 0.0

    return {
        "X": max_p,
        "max_p": max_p,
        "C_E": c_e,
        "N_plus": n_plus,
        "N_minus": n_minus,
        "Delta_E": delta_e,
        "D_E": delta_e,
        "S_E": s_e,
        "D_tilde_E": d_tilde,
        "omega_graph_sum": omega_graph_sum,
        "line_integral_ABCEA_sum": line_integral_plus,
        "line_integral_CEABC_sum": line_integral_minus,
        "C_E_equals_D_E": c_e == delta_e,
        "detected_cycles": total,
    }


def edge_orientation_table() -> dict[str, Any]:
    """Vollständige ω-Tabelle auf V×V und diskrete 1-Form."""
    edges: dict[str, int] = {}
    for src in LABELS:
        for dst in LABELS:
            if src == dst:
                continue
            key = f"{src}->{dst}"
            edges[key] = edge_omega(src, dst)
    return {
        "vertices": list(LABELS),
        "canonical_forward_edges": [f"{a}->{b}" for a, b in CANONICAL_FORWARD_EDGES],
        "canonical_gap_pattern_mod12": list(CANONICAL_GAP_PATTERN),
        "edge_omega": edges,
        "discrete_one_form_A": discrete_one_form(),
        "ABCEA_Omega": cycle_omega_graph(ABCEA_WORD),
        "CEABC_Omega_graph": cycle_omega_graph(CEABC_WORD),
        "ABCEA_omega_cycle": omega_cycle(ABCEA_WORD),
        "CEABC_omega_cycle": omega_cycle(CEABC_WORD),
        "ABCEA_line_integral_A": line_integral_one_form(ABCEA_WORD),
        "CEABC_line_integral_A": line_integral_one_form(CEABC_WORD),
    }


def circulation_report(max_p: int) -> dict[str, Any]:
    """Vollständiger Sagnac-Zirkulationsbericht."""
    circ = circulation_C_E(max_p)
    orient = edge_orientation_table()
    return {
        "theory": THEORY_ZIRKULATION,
        "theory_intuition": THEORY_SAGNAC,
        "X": max_p,
        "circulation": circ,
        "graph_orientation": orient,
        "relations": {
            "C_E_equals_D_E": circ["C_E_equals_D_E"],
            "C_E_formula": "C_E(X) = sum_{gamma detected} omega(gamma)",
            "D_E_formula": "D_E(X) = N_+(X) - N_-(X)",
            "S_E_formula": "S_E(X) = C_E(X) / (N_+ + N_-)",
            "one_form": "oint_ABCEA A = +1, oint_CEABC A = -oint_ABCEA A (via omega(gamma))",
        },
        "epistemic": {
            "C_E": "Definition",
            "edge_omega": "Definition",
            "C_E_bias_structure": "Hypothese (Chebyshev/mod-q/Dirichlet-L)",
        },
    }


def run(max_p: int = 1_000_000, output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    report = circulation_report(max_p)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC-Sagnac-Zirkulation C_E(X)")
    parser.add_argument("--max-p", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(max_p=args.max_p, output=args.output)
    circ = report["circulation"]
    orient = report["graph_orientation"]
    print("=== EABC Sagnac Zirkulation C_E ===")
    print(f"X={circ['X']}: C_E={circ['C_E']:+d}  D_E={circ['D_E']:+d}  S_E={circ['S_E']:+.4f}")
    print(f"N_plus={circ['N_plus']}  N_minus={circ['N_minus']}  D̃_E={circ['D_tilde_E']:+.3f}")
    print(f"C_E = D_E: {circ['C_E_equals_D_E']}")
    print()
    print("Kantenorientierung ω (kanonischer 4-Zyklus):")
    for e in orient["canonical_forward_edges"]:
        print(f"  ω({e}) = {orient['edge_omega'][e]:+d}")
    print(f"Ω_graph(ABCEA)={orient['ABCEA_Omega']}  ω_cycle(ABCEA)={orient['ABCEA_omega_cycle']}")
    print(f"Ω_graph(CEABC)={orient['CEABC_Omega_graph']}  ω_cycle(CEABC)={orient['CEABC_omega_cycle']}")
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
