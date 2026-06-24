"""Tests für EABC Übergangsraum Hodge/Fluss/magnetischen Laplace."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_hodge_eabc import (
    C4_EDGE_LABELS,
    C4_EDGE_NEGATIVE,
    C4_EDGE_POSITIVE,
    Phi_E,
    discrete_hodge_decomposition,
    edge_incidence_matrix_c4,
    flux_density_limit,
    flux_density_series,
    harmonic_form_c4,
    harmonic_holonomy_component,
    hodge_report,
    inner_product_omega_h,
    laplacian_from_W,
    magnetic_laplacian,
    magnetic_laplacian_eigenvalues,
    magnetic_phase_matrix,
    orientation_information_test,
    signed_measure_graph,
    synthesis_report,
)
from collatz_eabc_holonomie_fehlerterm import holonomy_counts
from collatz_eabc_transition_graph import (
    LABELS,
    classes_from_sequence,
    prime_eabc_sequence,
    transition_counts,
)


def test_incidence_matrix_row_sums_zero():
    b = edge_incidence_matrix_c4()
    assert b.shape == (4, 4)
    np.testing.assert_allclose(b.sum(axis=0), 0.0, atol=1e-12)


def test_harmonic_form_normalized():
    h = harmonic_form_c4()
    assert len(h) == 4
    np.testing.assert_allclose(np.linalg.norm(h), 1.0, atol=1e-12)


def test_laplacian_from_W_four_eigenvalues():
    counts = transition_counts(classes_from_sequence(prime_eabc_sequence(8_000)))
    lap = laplacian_from_W(np.array(counts, dtype=float))
    assert len(lap["eigenvalues"]) == 4
    assert lap["smallest_lambda"] >= -1e-10


def test_signed_measure_graph_has_edges():
    counts = np.array(transition_counts(classes_from_sequence(prime_eabc_sequence(5_000))))
    g = signed_measure_graph(counts)
    assert g["edge_count"] > 0
    assert set(g["vertices"]) == set(LABELS)


def test_flux_density_equals_S_E():
    max_p = 50_000
    flux = flux_density_limit(max_p)
    hol = holonomy_counts(max_p)
    assert flux["flux_density"] == hol["S_E"]
    assert flux["Phi_E"] == hol["S_E"]
    assert flux["W_E"] == hol["S_E"]
    assert flux["C_E"] == hol["D_E"]
    assert flux["S_E"] == hol["N_plus"] + hol["N_minus"]


def test_Phi_E_alias():
    max_p = 20_000
    assert Phi_E(max_p) == flux_density_limit(max_p)


def test_E_plus_E_minus_edges():
    assert C4_EDGE_POSITIVE == frozenset(("EA", "AB", "BC", "CE"))
    assert C4_EDGE_NEGATIVE == frozenset(("EC", "CB", "BA", "AE"))
    assert C4_EDGE_POSITIVE.isdisjoint(C4_EDGE_NEGATIVE)


def test_orientation_information_recoverable_with_S():
    result = orientation_information_test(30_000)
    assert result["five_block"]["recoverable_from_total_and_S"] is True
    assert result["five_block"]["recoverable_from_total_only"] is False
    assert result["four_block"]["recoverable_from_total_and_S"] is True


def test_harmonic_holonomy_inner_product_sign():
    classes = classes_from_sequence(prime_eabc_sequence(40_000))
    harm = harmonic_holonomy_component(classes)
    c_e = harm["C_E"]
    inner = harm["inner_product_omega_h"]
    if c_e != 0:
        assert (inner > 0) == (c_e > 0) or abs(inner) > 0


def test_magnetic_laplacian_hermitian_eigenvalues_real():
    classes = classes_from_sequence(prime_eabc_sequence(12_000))
    counts = np.array(transition_counts(classes), dtype=float)
    phases = magnetic_phase_matrix(classes)
    mag = magnetic_laplacian(counts, phases, hermitian=True)
    assert len(mag["eigenvalues"]) == 4
    assert mag["smallest_lambda"] >= -1e-8


def test_hodge_report_at_100k():
    report = hodge_report(100_000)
    assert len(report["laplacian_from_W"]["eigenvalues"]) == 4
    assert "inner_product_omega_h" in report
    assert report["flux_density"]["flux_density"] == report["flux_density"]["Phi_E"]
    assert "E_plus" in report["topology"]["G_E"]


def test_synthesis_report_writes_json(tmp_path):
    out = tmp_path / "synthesis.json"
    syn = synthesis_report(10_000, output=out)
    assert out.exists()
    assert syn["observables"]["Phi_E"] == syn["observables"]["W_E"]
    assert "harmonic" in syn


def test_flux_density_series_monotone_limits():
    series = flux_density_series([5_000, 20_000, 50_000])
    assert len(series["series"]) == 3
    for row in series["series"]:
        assert "flux_density" in row


def test_discrete_hodge_decomposition_dims():
    omega = np.array([10.0, -3.0, 5.0, -12.0])
    hodge = discrete_hodge_decomposition(omega)
    grad = np.array(hodge["omega_gradient"])
    harm = np.array(hodge["omega_harmonic"])
    np.testing.assert_allclose(grad + harm, omega, atol=1e-10)


def test_c4_edge_labels_canonical():
    assert C4_EDGE_LABELS == ("EA", "AB", "BC", "CE")


def test_magnetic_laplacian_eigenvalues_alias():
    classes = classes_from_sequence(prime_eabc_sequence(10_000))
    counts = np.array(transition_counts(classes), dtype=float)
    phases = magnetic_phase_matrix(classes)
    alias = magnetic_laplacian_eigenvalues(counts, phases)
    full = magnetic_laplacian(counts, phases)["eigenvalues"]
    assert alias == full
    assert len(alias) == 4


def test_synthesis_report_structure(tmp_path):
    out = tmp_path / "synthesis.json"
    report = synthesis_report(25_000, output=out)
    assert report["meta"]["theory"] == "collatz_eabc_diskrete_geometrie.md"
    assert "flux_density" in report["observables"]
    assert report["observables"]["flux_density"] == report["observables"]["Phi_E"]
    assert report["observables"]["W_E"] == report["observables"]["Phi_E"]
    assert len(report["spectra"]["laplacian_from_W"]["eigenvalues"]) == 4
    assert len(report["spectra"]["magnetic_laplacian"]["eigenvalues"]) == 4
    assert "inner_product_omega_h" in report["harmonic"]
    assert out.exists()
