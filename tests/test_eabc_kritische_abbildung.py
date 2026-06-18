"""Tests für EABC kritische Abbildung — Geschwindigkeitsmodell s_v(x)."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_kritische_abbildung import (
    ABCEA_EDGES,
    CANONICAL_GAP_PATTERN,
    CEABC_EDGES,
    GAMMA_1_APPROX,
    SOURCE_X,
    dual_circuit_report,
    eabc_circuit_report,
    example_gamma_1,
    gamma_v,
    holonomy_sign,
    run,
    s_v,
    x_from_gamma,
    x_n_v,
    zeta_imaginary_parts,
)


def test_s_v_source_point():
    assert s_v(0.5, 1.0) == complex(0.5, 0.0)
    assert s_v(0.5, 10.0) == complex(0.5, 0.0)


def test_s_v_and_inverse():
    for v in (1.0, 2.0, 10.0):
        x = 14.634725
        s = s_v(x, v)
        assert math.isclose(s.real, 0.5)
        assert math.isclose(s.imag, gamma_v(x, v))
        assert math.isclose(x_from_gamma(s.imag, v), x, rel_tol=0, abs_tol=1e-9)


def test_gamma_1_examples():
    ex = example_gamma_1()
    assert math.isclose(ex["v=1"], 0.5 + GAMMA_1_APPROX, rel_tol=0, abs_tol=1e-4)
    assert math.isclose(ex["v=2"], 0.5 + GAMMA_1_APPROX / 2, rel_tol=0, abs_tol=1e-4)
    assert math.isclose(ex["v=10"], 0.5 + GAMMA_1_APPROX / 10, rel_tol=0, abs_tol=1e-4)


def test_x_n_v_first_zero():
    assert math.isclose(x_n_v(1, 1.0), 0.5 + GAMMA_1_APPROX, rel_tol=0, abs_tol=1e-6)
    assert math.isclose(x_n_v(1, 2.0), 0.5 + GAMMA_1_APPROX / 2, rel_tol=0, abs_tol=1e-6)


def test_zeta_imaginary_parts_fallback():
    gammas = zeta_imaginary_parts(3)
    assert len(gammas) == 3
    assert gammas[0] > 14.0


def test_abcea_gap_pattern():
    lengths = tuple(ell for _, _, ell in ABCEA_EDGES)
    assert lengths == CANONICAL_GAP_PATTERN


def test_ceabc_gap_pattern_cyclic_shift():
    lengths = tuple(ell for _, _, ell in CEABC_EDGES)
    assert lengths == (2, 4, 2, 4)
    assert CEABC_EDGES[0][0] == "C"


def test_abcea_circuit_v1_canonical():
    circ = eabc_circuit_report(v=1.0, orientation="ABCEA", length_model="canonical")
    assert circ["holonomy_sign"] == 1
    assert circ["total_length"] == 12
    endpoint = circ["endpoint"]
    assert endpoint is not None
    assert math.isclose(endpoint["x"], 12.5)
    assert math.isclose(endpoint["gamma_cumulative"], 12.0)
    assert math.isclose(endpoint["s_v"]["im"], 12.0)

    velocities = [seg["velocity"] for seg in circ["segments"] if "velocity" in seg]
    assert velocities == [1.0, 1.0, 1.0, 1.0]


def test_ceabc_holonomy_minus_one():
    circ = eabc_circuit_report(v=1.0, orientation="CEABC")
    assert circ["holonomy_sign"] == -1
    assert circ["total_length"] == 12


def test_dual_circuit_same_span():
    dual = dual_circuit_report(v=1.0)
    ab = dual["ABCEA"]
    ce = dual["CEABC"]
    assert ab["total_length"] == ce["total_length"]
    assert holonomy_sign("ABCEA") == 1
    assert holonomy_sign("CEABC") == -1


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "kritische.json"
    report = run(n_zeros=3, v_circuit=1.0, output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["meta"]["theory"] == "collatz_eabc_kritische_abbildung.md"
    assert loaded["meta"]["theory_zirkulation"] == "collatz_eabc_zirkulationshypothese.md"
    assert "eabc_circuits" in loaded
    assert loaded["eabc_circuits"]["ABCEA"]["holonomy_sign"] == 1
    assert loaded["eabc_circuits"]["CEABC"]["holonomy_sign"] == -1
    assert report["output_path"] == str(out)
