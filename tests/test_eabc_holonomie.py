"""Tests für EABC-Holonomie (collatz_eabc_holonomie_test.py)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_discrete_associator import prove_v4_klein_associativity
from collatz_eabc_holonomie_test import (
    chi_global,
    chi_leg_score,
    chi_quad_legs,
    enumerate_quadruplets,
    holonomy_chi_connection,
    holonomy_flux_phi_quad,
    octonion_gamma_holonomy_samples,
    omega_orientation,
    quadruplet_chirality,
    run,
)
from eabc_from_lean import Chirality, EClass, class_of, is_prime_quadruplet, q


def test_v4_klein_proof_associative():
    proof = prove_v4_klein_associativity()
    assert proof["associative"]
    assert proof["closure_mod12"]
    assert proof["identity_E"]
    assert proof["self_inverse_nonE"]
    assert proof["isomorphism_to_Z2_squared"]
    assert proof["naive_associator_all_zero"]
    assert proof["triples_tested"] == 64
    assert proof["counterexample_count"] == 0


def test_omega_orientation_signs():
    assert omega_orientation(Chirality.ABCE) == 1
    assert omega_orientation(Chirality.CEAB) == -1


def test_quadruplet_q5_abce_omega():
    assert is_prime_quadruplet(5)
    assert quadruplet_chirality(5) is Chirality.ABCE
    classes = tuple(class_of(n) for n in q(5))
    assert chi_leg_score(classes) == 0  # type: ignore[arg-type]
    quads = enumerate_quadruplets(100)
    q5 = next(r for r in quads if r["p"] == 5)
    assert q5["omega"] == 1
    assert q5["word"] == "ABCE"


def test_all_quadruplet_chi_leg_zero_up_to_5000():
    quads = enumerate_quadruplets(5000)
    assert quads
    assert all(r["chi_leg"] == 0 for r in quads)


def test_chi_global_matches_invarianzprogramm_formula():
    report = chi_global(200)
    v = report["counts"]
    pi = report["pi_gt3"]
    expected = ((v["E"] + v["C"]) - (v["A"] + v["B"])) / pi
    assert abs(report["chi"] - expected) < 1e-12


def test_chi_quad_legs_zero():
    legs = chi_quad_legs(5000)
    assert legs["quadruplet_count"] > 0
    assert legs["chi_legs"] == 0.0
    assert legs["chi_fluct_legs"] == 0


def test_holonomy_flux_bounded():
    quads = enumerate_quadruplets(5000)
    flux = holonomy_flux_phi_quad(quads)
    assert -1.0 <= flux["phi_quad"] <= 1.0
    assert flux["abce_count"] + flux["ceab_count"] == flux["quadruplet_count"]


def test_holonomy_chi_connection_honest():
    conn = holonomy_chi_connection(5000)
    assert conn["all_quadruplet_chi_leg_zero"]
    assert "nicht identisch" in conn["verdict"]


def test_octonion_gamma_holonomy_generic_nonzero():
    oct_report = octonion_gamma_holonomy_samples()
    generic = next(
        s for s in oct_report["holonomy_samples"] if s["label"] == "generic_o_e1_e2_e4"
    )
    assert generic["algebraic_associator_norm"] > 0.0
    assert oct_report["canonical_summary"]["generic_associator_nonzero"]


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "hol.json"
    report = run(limit=500, output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["v4_associativity_proof"]["associative"]
    assert loaded["holonomy_chi_connection"]["all_quadruplet_chi_leg_zero"]
    assert report["output_path"] == str(out)
