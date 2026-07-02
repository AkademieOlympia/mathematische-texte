"""Tests für EABC Brachistochrone / Fermat-Prinzip."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_brachistochrone import (
    abcea_trajectory_xy,
    compare_paths_for_potential,
    make_velocity_func,
    perturbed_path_gradient,
    run,
    straight_path,
    travel_time_integral,
    velocity_from_potential,
)
from collatz_eabc_kritische_abbildung import GAMMA_1_APPROX, SOURCE_X


def test_velocity_from_potential_log():
    assert velocity_from_potential(math.e, "log") == math.e
    assert math.isclose(velocity_from_potential(2.0, "inverse_log"), 0.5, rel_tol=0, abs_tol=1e-9)


def test_velocity_from_potential_information():
    v = velocity_from_potential(0.0, "information")
    assert math.isclose(v, 1.0, rel_tol=0, abs_tol=1e-9)


def test_travel_time_integral_uniform_speed():
    pts = [(0.0, 0.0), (3.0, 4.0)]
    t = travel_time_integral(pts, lambda _x, _g: 1.0)
    assert math.isclose(t, 5.0, rel_tol=0, abs_tol=1e-12)


def test_abcea_trajectory_starts_at_source():
    pts = abcea_trajectory_xy(gamma_ref=GAMMA_1_APPROX)
    assert len(pts) == 5
    assert math.isclose(pts[0][0], SOURCE_X, rel_tol=0, abs_tol=1e-12)
    assert math.isclose(pts[0][1], 0.0, rel_tol=0, abs_tol=1e-12)


def test_straight_path_unchanged():
    pts = [(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)]
    assert straight_path(pts) == pts


def test_perturbed_path_has_more_points():
    pts = abcea_trajectory_xy(gamma_ref=GAMMA_1_APPROX)
    _, grad = make_velocity_func("log")
    bent = perturbed_path_gradient(pts, grad, scale=0.1)
    assert len(bent) > len(pts)


def test_compare_paths_for_log_potential():
    pts = abcea_trajectory_xy(gamma_ref=GAMMA_1_APPROX)
    row = compare_paths_for_potential("log", pts, d_e_global=10.0)
    assert row["T_straight"] > 0
    assert "rel_gradient_vs_straight" in row
    assert isinstance(row["path_bends_significantly"], bool)


def test_compare_paths_all_five_potentials():
    pts = abcea_trajectory_xy(gamma_ref=GAMMA_1_APPROX)
    for name in ("log", "zeta", "chirality", "curvature", "information"):
        row = compare_paths_for_potential(name, pts, d_e_global=5.0)
        assert row["potential"] == name
        assert row["T_straight"] > 0


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "brachistochrone.json"
    report = run(max_p=5000, output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["meta"]["theory"] == "collatz_eabc_brachistochrone.md"
    assert loaded["meta"]["theory_epistemik"] == "collatz_eabc_epistemik_physik.md"
    assert len(loaded["comparisons"]) == 5
    assert "summary" in loaded
    assert report["output_path"] == str(out)
