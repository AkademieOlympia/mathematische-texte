"""Tests für EABC Oktaeder-Umgebung (Φ_E, ⟨ω,h⟩, Schalen-Gewicht)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_oktaeder_flux_test import (
    build_octahedron_graph,
    compare_c4_vs_octahedron,
    octahedron_incidence_matrix,
    r4_jacobi,
    run,
    verify_r8_oeis,
)


def test_r8_oeis_crosscheck():
    assert verify_r8_oeis()["ok"]


def test_r4_jacobi_small_values():
    assert r4_jacobi(0) == 1
    assert r4_jacobi(1) == 8
    assert r4_jacobi(2) == 24


def test_octahedron_graph_structure():
    g = build_octahedron_graph()
    assert len(g.vertices) == 6
    assert len(g.edges) == 12
    assert len(g.equatorial_edge_indices) == 4
    assert len(g.polar_plus_indices) == 4
    assert len(g.polar_minus_indices) == 4


def test_incidence_row_sums_zero():
    g = build_octahedron_graph()
    b = octahedron_incidence_matrix(g)
    assert b.shape == (6, 12)
    np.testing.assert_allclose(b.sum(axis=0), 0.0, atol=1e-12)


def test_phi_oct_equatorial_matches_c4():
    max_p = 30_000
    report = compare_c4_vs_octahedron(max_p)
    cmp = report["comparison"]
    assert cmp["Phi_equatorial_matches_c4"]
    assert cmp["C_E_matches"]
    assert cmp["inner_product_equatorial_matches_c4"]
    assert abs(report["c4_baseline"]["Phi_E"] - report["octahedron_unweighted"]["Phi_equatorial"]) < 1e-12


def test_harmonic_pairing_nonzero_at_scale():
    max_p = 50_000
    report = compare_c4_vs_octahedron(max_p)
    c4 = report["c4_baseline"]
    oct_u = report["octahedron_unweighted"]
    assert c4["inner_product_omega_h"] != 0.0
    assert oct_u["inner_product_octahedron"] != 0.0
    assert report["verdict"]["supports_harmonic_pairing"]


def test_shell_weighted_phi_computed():
    max_p = 100_000
    report = compare_c4_vs_octahedron(max_p)
    phi_shell = report["octahedron_shell_weighted"]["Phi_shell_weighted"]
    assert phi_shell is not None
    assert abs(phi_shell) > 1e-12


def test_shell_weighted_same_sign_at_100k():
    max_p = 100_000
    report = compare_c4_vs_octahedron(max_p)
    phi_c4 = report["c4_baseline"]["Phi_E"]
    phi_shell = report["octahedron_shell_weighted"]["Phi_shell_weighted"]
    assert np.sign(phi_c4) == np.sign(phi_shell)


def test_polar_preference_matches_orientation():
    max_p = 20_000
    report = compare_c4_vs_octahedron(max_p)
    oct_u = report["octahedron_unweighted"]
    assert oct_u["polar_plus_flux"] == float(oct_u["N_plus"])
    assert oct_u["polar_minus_flux"] == float(oct_u["N_minus"])


def test_run_writes_json(tmp_path):
    out = tmp_path / "octa_test.json"
    report = run(max_p=8_000, output=out)
    assert out.exists()
    loaded = json.loads(out.read_text())
    assert loaded["meta"]["max_p"] == 8_000
    assert "verdict" in loaded
    assert report["output_path"] == str(out)


def test_flux_series_monotone_s_e():
    report = compare_c4_vs_octahedron(25_000)
    series = report["flux_density_series"]
    s_values = [row["S_E"] for row in series]
    assert s_values == sorted(s_values)
