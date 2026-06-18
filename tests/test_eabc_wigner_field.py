"""Tests für EABC Wigner-Feld-Analogie."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_holonomie_fehlerterm import holonomy_counts
from collatz_eabc_transition_graph import (
    ABCE_WORD,
    CEAB_WORD,
    LABELS,
    chi_pfad_sliding,
    classes_from_sequence,
    prime_eabc_sequence,
    sliding_windows,
)
from collatz_eabc_wigner_field import (
    build_w_matrix_from_windows,
    quadruplet_indicator,
    sign_domain_analysis,
    w_e_counts,
    w_e_profile,
    wigner_correlation_entry,
    wigner_field_report,
)


def test_quadruplet_indicator():
    assert quadruplet_indicator(ABCE_WORD) == 1
    assert quadruplet_indicator(CEAB_WORD) == -1
    assert quadruplet_indicator("EABC") == 0


def test_wigner_correlation_entry():
    chi_abce = [0, 1, 1, 1]  # A,B,C present once each in ABCE window
    assert wigner_correlation_entry(chi_abce, chi_abce, 1) == 3


def test_w_e_counts_consistency():
    max_p = 30_000
    w = w_e_counts(max_p)
    seq = prime_eabc_sequence(max_p)
    pfad = chi_pfad_sliding(classes_from_sequence(seq))
    assert w["N_ABCE"] == pfad["abce_windows"]
    assert w["N_CEAB"] == pfad["ceab_windows"]
    assert w["W_E"] == w["N_ABCE"] - w["N_CEAB"]


def test_w_e_profile_cumulative():
    max_p = 10_000
    classes = classes_from_sequence(prime_eabc_sequence(max_p))
    profile = w_e_profile(classes)
    if profile:
        last = profile[-1]
        assert last["W_E_cumulative"] == last["cum_ABCE"] - last["cum_CEAB"]


def test_build_w_matrix_shape_and_labels():
    classes = classes_from_sequence(prime_eabc_sequence(5_000))
    windows = sliding_windows(classes, width=4)
    w = build_w_matrix_from_windows(windows)
    mat = np.array(w["matrix"])
    assert mat.shape == (4, 4)
    assert set(w["matrix_labeled"].keys()) == set(LABELS)


def test_sign_domain_analysis():
    domains = sign_domain_analysis([1, 1, -1, -1, 0, 2])
    assert domains["sign_changes"] >= 2
    assert domains["positive_steps"] >= 1
    assert domains["negative_steps"] >= 1


def test_four_vs_five_block_in_report():
    report = wigner_field_report(20_000)
    assert report["W_E"]["block_width"] == 4
    assert report["D_E"]["block_width"] == 5
    hol = holonomy_counts(20_000)
    assert report["D_E"]["D_E"] == hol["D_E"]
    assert "four_vs_five_block" in report


def test_near_zero_modes_stub():
    report = wigner_field_report(15_000)
    stub = report["near_zero_modes"]
    assert stub["near_zero_mode_count"] >= 1
    assert len(stub["eigenvalues"]) == 4
    assert "dirac_stub" in stub


def test_wigner_report_at_100k():
    report = wigner_field_report(100_000)
    w_mat = report["W_matrix_sliding"]["matrix_labeled"]
    assert all(row_label in w_mat for row_label in LABELS)
    assert report["sign_structure"]["E"]["A"] in ("+", "-", "0")
