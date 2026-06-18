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
    compare_holonomy_sensor_trajectories,
    compare_path_times,
    dual_circuit_report,
    eabc_circuit_report,
    edge_velocities_from_gaps,
    example_gamma_1,
    gaps_abcea_to_eaabc,
    gaps_eaabc_to_abcea,
    holonomy_sensor_trajectory,
    gamma_v,
    holonomy_sign,
    linear_round_trip_time,
    path_time_T,
    prime_window_gap_samples,
    run,
    semicircle_chain_time,
    sensor_trajectory_points,
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
    assert "holonomy_sensor" in loaded
    assert loaded["holonomy_sensor"]["edge_velocities"]["v_EA"] > 0
    assert report["output_path"] == str(out)


def test_gap_order_eaabc_abcea_roundtrip():
    eaabc = (4, 2, 4, 2)
    abcea = gaps_eaabc_to_abcea(eaabc)
    assert abcea == CANONICAL_GAP_PATTERN
    assert gaps_abcea_to_eaabc(abcea) == eaabc


def test_edge_velocities_from_gaps_gamma_1():
    ev = edge_velocities_from_gaps(CANONICAL_GAP_PATTERN, GAMMA_1_APPROX)
    assert math.isclose(ev["v_EA"], GAMMA_1_APPROX / 4, rel_tol=0, abs_tol=1e-9)
    assert math.isclose(ev["v_AB"], GAMMA_1_APPROX / 2, rel_tol=0, abs_tol=1e-9)
    assert math.isclose(ev["v_BC"], GAMMA_1_APPROX / 4, rel_tol=0, abs_tol=1e-9)
    assert math.isclose(ev["v_CE"], GAMMA_1_APPROX / 2, rel_tol=0, abs_tol=1e-9)
    ev_ea = edge_velocities_from_gaps((4, 2, 4, 2), GAMMA_1_APPROX, gap_order="EAABC")
    assert ev_ea["v_EA"] == ev["v_EA"]


def test_holonomy_sensor_constant_delta_gamma():
    traj = holonomy_sensor_trajectory("ABCEA", gamma_ref=GAMMA_1_APPROX)
    deltas = [seg["delta_gamma"] for seg in traj["segments"] if "delta_gamma" in seg]
    assert len(deltas) == 4
    assert all(math.isclose(d, GAMMA_1_APPROX, rel_tol=0, abs_tol=1e-9) for d in deltas)
    assert math.isclose(traj["total_gamma"], 4 * GAMMA_1_APPROX, rel_tol=0, abs_tol=1e-9)


def test_compare_holonomy_sensor_abcea_ceabc():
    cmp = compare_holonomy_sensor_trajectories(gamma_ref=GAMMA_1_APPROX)
    assert cmp["holonomy_contrast"]["sign_ABCEA"] == 1
    assert cmp["holonomy_contrast"]["sign_CEABC"] == -1
    assert cmp["holonomy_contrast"]["same_total_length"]
    assert cmp["holonomy_contrast"]["same_total_gamma"]
    assert cmp["holonomy_contrast"]["same_edge_velocities"]


def test_prime_window_gap_samples():
    samples = prime_window_gap_samples(max_p=5000, limit=1)
    assert len(samples) >= 1
    for row in samples:
        assert row["gaps_abcea"] == [2, 4, 2, 4]
        assert "edge_velocities" in row


def test_linear_round_trip_vertical_chain():
    points = [complex(0.5, 0.0), complex(0.5, 1.0), complex(0.5, 3.0)]
    assert math.isclose(linear_round_trip_time(points), 1.0 + 2.0 + 3.0, rel_tol=0, abs_tol=1e-12)


def test_semicircle_chain_pi_over_two_factor():
    points = [complex(0.5, 0.0), complex(0.5, 2.0)]
    lin = linear_round_trip_time(points)
    semi = semicircle_chain_time(points)
    assert math.isclose(semi / lin, math.pi / 2, rel_tol=0, abs_tol=1e-12)


def test_path_time_T_canonical_gaps():
    """T = Σ ℓ_j² / γ_ref für (2,4,2,4)."""
    t = path_time_T(CANONICAL_GAP_PATTERN, GAMMA_1_APPROX)
    expected = (4 + 16 + 4 + 16) / GAMMA_1_APPROX
    assert math.isclose(t, expected, rel_tol=0, abs_tol=1e-12)


def test_compare_path_times_abcea_gamma_1():
    cmp = compare_path_times("ABCEA", gamma_ref=GAMMA_1_APPROX)
    expected_linear = 8 * GAMMA_1_APPROX
    expected_semi = 4 * math.pi * GAMMA_1_APPROX
    assert math.isclose(cmp["T_linear"], expected_linear, rel_tol=0, abs_tol=1e-9)
    assert math.isclose(cmp["T_semicircle"], expected_semi, rel_tol=0, abs_tol=1e-9)
    assert math.isclose(cmp["ratio_semi_over_linear"], math.pi / 2, rel_tol=0, abs_tol=1e-12)
    assert cmp["same_linear_time_both_orientations"]
    assert cmp["same_semicircle_time_both_orientations"]


def test_sensor_trajectory_points_start_at_p():
    points = sensor_trajectory_points("ABCEA", gamma_ref=GAMMA_1_APPROX)
    assert len(points) == 5
    assert math.isclose(points[0].real, SOURCE_X, rel_tol=0, abs_tol=1e-12)
    assert math.isclose(points[0].imag, 0.0, rel_tol=0, abs_tol=1e-12)
