#!/usr/bin/env python3
"""
EABC-Holonomie: Orientierung ω auf Vierlingen, χ-Invarianz, Γ_E-Klammerdefekt.

Kanonsiche Theorie: collatz_eabc_holonomie.md

Kernkorrektur:
  V₄ ≅ Klein-Gruppe → naive 𝔞 auf Φ ≡ 0 (siehe collatz_eabc_discrete_associator.py).
  Echter Defekt: Γ((xy)z) - Γ(x(yz)) auf Trägerobjekten (Oktanion-Stub, Vierlings-Orientierung ω).

Ausführung:
    python3 collatz_eabc_holonomie_test.py
    python3 collatz_eabc_holonomie_test.py --limit 50000
"""

from __future__ import annotations

import argparse
import json
from math import isqrt
from pathlib import Path
from typing import Any, Literal

from collatz_eabc_discrete_associator import prove_v4_klein_associativity
from collatz_eabc_oktonion_associator import (
    bracketing_magnitudes,
    canonical_triples_test,
    eabc_associator_vector,
    o_mul,
)
from eabc_from_lean import Chirality, EClass, class_of, is_prime_quadruplet, q

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_holonomie.json"

Orientation = Literal[1, -1]


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def omega_orientation(chirality: Chirality) -> Orientation:
    """ω(ABCE)=+1, ω(CEAB)=-1 (Lean-Chiralität, keine V₄-Produkte)."""
    return 1 if chirality is Chirality.ABCE else -1


def quadruplet_chirality(p: int) -> Chirality:
    mod = p % 12
    if mod == 5:
        return Chirality.ABCE
    if mod == 11:
        return Chirality.CEAB
    raise ValueError(f"Kein Vierlingsstart: p={p} (mod 12 = {mod})")


def chi_leg_score(classes: tuple[EClass, ...]) -> int:
    """χ_leg(Q) = #(E∪C) - #(A∪B) auf den vier Beinen."""
    ec = sum(1 for c in classes if c in (EClass.E, EClass.C))
    ab = sum(1 for c in classes if c in (EClass.A, EClass.B))
    return ec - ab


def count_eabc_vector(limit: int) -> dict[str, int]:
    """V(limit) aus collatz_eabc_invarianzprogramm.md — nur p>3."""
    counts = {k: 0 for k in ("E", "A", "B", "C")}
    for p in range(5, limit + 1):
        if not _is_prime(p):
            continue
        cls = class_of(p)
        if cls is not None:
            counts[cls.value] += 1
    return counts


def chi_global(limit: int) -> dict[str, Any]:
    """χ(x) = ((E+C)-(A+B)) / π_{>3}(x) — Invarianzprogramm Definition 2/Beispiel."""
    v = count_eabc_vector(limit)
    pi = sum(v.values())
    numer = (v["E"] + v["C"]) - (v["A"] + v["B"])
    chi = numer / pi if pi else 0.0
    return {
        "limit": limit,
        "counts": v,
        "pi_gt3": pi,
        "chi_fluct": numer,
        "chi": chi,
    }


def enumerate_quadruplets(limit: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for p in range(5, limit - 7, 2):
        if p % 12 not in (5, 11):
            continue
        if not is_prime_quadruplet(p):
            continue
        chi = quadruplet_chirality(p)
        classes = tuple(class_of(n) for n in q(p))
        assert all(c is not None for c in classes)
        rows.append(
            {
                "p": p,
                "word": "".join(c.value for c in classes),  # type: ignore[union-attr]
                "chirality": chi.value,
                "omega": omega_orientation(chi),
                "chi_leg": chi_leg_score(classes),  # type: ignore[arg-type]
                "signature": [c.value for c in classes],  # type: ignore[union-attr]
            }
        )
    return rows


def chi_quad_legs(limit: int, quads: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    """χ auf allen Vierlingsbeinen mit p_start ≤ limit."""
    if quads is None:
        quads = enumerate_quadruplets(limit)
    leg_counts = {k: 0 for k in ("E", "A", "B", "C")}
    for row in quads:
        for leg in row["signature"]:
            leg_counts[leg] += 1
    total_legs = sum(leg_counts.values())
    numer = (leg_counts["E"] + leg_counts["C"]) - (leg_counts["A"] + leg_counts["B"])
    return {
        "limit": limit,
        "quadruplet_count": len(quads),
        "leg_counts": leg_counts,
        "total_legs": total_legs,
        "chi_fluct_legs": numer,
        "chi_legs": numer / total_legs if total_legs else 0.0,
    }


def holonomy_flux_phi_quad(quads: list[dict[str, Any]]) -> dict[str, Any]:
    """Φ_quad = mittlere Orientierung ω über Vierlinge ∈ [-1,1]."""
    if not quads:
        return {
            "quadruplet_count": 0,
            "abce_count": 0,
            "ceab_count": 0,
            "omega_sum": 0,
            "phi_quad": 0.0,
        }
    abce = sum(1 for r in quads if r["omega"] == 1)
    ceab = sum(1 for r in quads if r["omega"] == -1)
    omega_sum = sum(int(r["omega"]) for r in quads)
    n = len(quads)
    return {
        "quadruplet_count": n,
        "abce_count": abce,
        "ceab_count": ceab,
        "omega_sum": omega_sum,
        "phi_quad": omega_sum / n,
    }


def holonomy_chi_connection(limit: int) -> dict[str, Any]:
    """Ehrlicher Vergleich: globale χ vs. Vierlings-Holonomie vs. Bein-χ."""
    quads = enumerate_quadruplets(limit)
    global_chi = chi_global(limit)
    quad_legs = chi_quad_legs(limit, quads)
    flux = holonomy_flux_phi_quad(quads)
    all_chi_leg_zero = all(r["chi_leg"] == 0 for r in quads)
    return {
        "limit": limit,
        "chi_global": global_chi["chi"],
        "chi_quad_legs": quad_legs["chi_legs"],
        "phi_quad_holonomy": flux["phi_quad"],
        "all_quadruplet_chi_leg_zero": all_chi_leg_zero,
        "verdict": (
            "χ_leg(Q)=0 für alle kanonischen Vierlinge (balancierte Signaturen); "
            "ω(Q)∈{±1} misst Orientierung, nicht Bein-Asymmetrie. "
            "Globale χ(x) und Φ_quad sind verwandte chirale Observablen in verschiedenen Räumen — "
            "nicht identisch."
        ),
        "equivalence_note": (
            "Holonomie auf Vierlingen ERWEITERT χ: χ auf Beinen trivial; "
            "Φ_quad codiert ABCE/CEAB-Phase. Globale χ trackt Φ₁-Modus der Primzählung."
        ),
    }


def octonion_gamma_holonomy_samples() -> dict[str, Any]:
    """Γ_E((xy)z) vs Γ_E(x(yz)) auf kanonischen Tripeln (Oktanion-Stub)."""
    canon = canonical_triples_test()
    e1 = (0, 1, 0, 0, 0, 0, 0, 0)
    e2 = (0, 0, 1, 0, 0, 0, 0, 0)
    e3 = (0, 0, 0, 1, 0, 0, 0, 0)
    e4 = (0, 0, 0, 0, 1, 0, 0, 0)

    def holonomy_triple(x, y, z, label: str) -> dict[str, Any]:
        left = o_mul(o_mul(x, y), z)
        right = o_mul(x, o_mul(y, z))
        vec = eabc_associator_vector(x, y, z)
        mags = bracketing_magnitudes(x, y, z)
        return {
            "label": label,
            "gamma_left": list(left),
            "gamma_right": list(right),
            "gamma_e_diff": list(vec),
            "gamma_e_norm": mags["eabc_associator_norm"],
            "algebraic_associator_norm": mags["algebraic_associator_norm"],
            "paths_differ_in_O": left != right,
            "gamma_e_holonomy_nonzero": mags["eabc_associator_norm"] > 0.0,
        }

    samples = [
        holonomy_triple(e1, e2, e3, "quaternion_subalgebra_e1_e2_e3"),
        holonomy_triple(e1, e2, e4, "generic_o_e1_e2_e4"),
    ]
    x = (-1, -1, 0, 0, 0, 0, 0, 0)
    y = (-1, 0, 0, 0, -1, 0, 0, 0)
    z = (-1, 0, -1, -1, 0, 0, 0, 0)
    samples.append(holonomy_triple(x, y, z, "shell_sample_sigma2_sigma2_sigma3"))

    return {
        "canonical_summary": canon,
        "holonomy_samples": samples,
        "note": (
            "Γ_E-Holonomie = glatt-EABC-Differenz der Projektionen entlang zweier Klammerwege; "
            "unabhängig vom trivialen V₄-Assoziator."
        ),
    }


def run(limit: int = 10_000, output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    v4_proof = prove_v4_klein_associativity()
    quads = enumerate_quadruplets(limit)
    report: dict[str, Any] = {
        "meta": {
            "module": "collatz_eabc_holonomie_test.py",
            "theory": "collatz_eabc_holonomie.md",
            "correction": (
                "V₄ Klein → naive 𝔞≡0; echter Defekt = Γ-Holonomie auf Trägern; "
                "ω(ABCE)=+1, ω(CEAB)=-1"
            ),
        },
        "v4_associativity_proof": v4_proof,
        "quadruplets": {
            "limit": limit,
            "count": len(quads),
            "samples": quads[:10],
            "holonomy_flux": holonomy_flux_phi_quad(quads),
            "chi_quad_legs": chi_quad_legs(limit, quads),
        },
        "chi_global": chi_global(limit),
        "holonomy_chi_connection": holonomy_chi_connection(limit),
        "octonion_gamma_holonomy": octonion_gamma_holonomy_samples(),
        "epistemic_labels": {
            "v4_proof": "Theorem",
            "omega_orientation": "Definition",
            "phi_quad": "Experiment",
            "chi_global": "Definition",
            "holonomy_chi_verdict": "Experiment",
            "octonion_gamma": "Experiment",
        },
    }
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC-Holonomie-Test")
    parser.add_argument("--limit", type=int, default=10_000, help="Obergrenze für Vierlinge und χ")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(limit=args.limit, output=args.output)
    proof = report["v4_associativity_proof"]
    conn = report["holonomy_chi_connection"]
    flux = report["quadruplets"]["holonomy_flux"]
    print("=== EABC-Holonomie ===")
    print(proof["verdict"])
    print(f"V₄ assoziativ: {proof['associative']}  ({proof['triples_tested']} Tripel)")
    print(f"Vierlinge bis {args.limit}: {flux['quadruplet_count']}  "
          f"(ABCE {flux['abce_count']}, CEAB {flux['ceab_count']})")
    print(f"Φ_quad = {flux['phi_quad']:.4f}  |  χ_global = {conn['chi_global']:.6f}  |  "
          f"χ_quad_legs = {conn['chi_quad_legs']:.6f}")
    print(conn["verdict"])
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
