"""Tests für EABC chiraler Transport / Helizität und Brachistochrone-Doppelkanal."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_brachistochrone import (
    abcea_trajectory_xy,
    birefringent_velocity_pair,
    make_velocity_func,
    travel_time_birefringent,
    travel_time_integral,
)
from collatz_eabc_chirale_transport import (
    ABCEA_WORD,
    CEABC_WORD,
    accumulate_phases_along_windows,
    chiral_transport_report,
    chirality_flux_from_counts,
    helicity_channel,
    holonomy_phase_difference,
    holonomy_unitary_phases,
    link_phase_difference_to_D_E,
    theta_edge_from_gaps,
)
from collatz_eabc_holonomie_fehlerterm import holonomy_counts
from collatz_eabc_kritische_abbildung import CANONICAL_GAP_PATTERN
from collatz_eabc_transition_graph import (
    classes_from_sequence,
    prime_eabc_sequence,
    sliding_windows,
)


def test_helicity_channel_mapping():
    assert helicity_channel(ABCEA_WORD) == "R"
    assert helicity_channel(CEABC_WORD) == "L"
    assert helicity_channel("ABCDE") is None


def test_holonomy_phase_difference_signs():
    theta = math.pi / 4
    phi_abcea = holonomy_phase_difference(ABCEA_WORD, theta_edge=theta)
    phi_ceabc = holonomy_phase_difference(CEABC_WORD, theta_edge=theta)
    assert phi_abcea == 4 * theta
    assert phi_ceabc == -4 * theta
    assert phi_abcea == -phi_ceabc


def test_holonomy_unitary_phases():
    theta = math.pi / 4
    u = holonomy_unitary_phases(ABCEA_WORD, theta_edge=theta)
    assert abs(u["U_R"] - complex(math.cos(4 * theta), math.sin(4 * theta))) < 1e-12
    assert u["U_L"] == 1 + 0j
    assert u["phase_difference"] == 4 * theta


def test_accumulate_phases_matches_counts():
    max_p = 50_000
    theta = theta_edge_from_gaps()
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    windows = sliding_windows(classes, width=5)
    phases = accumulate_phases_along_windows(windows, theta_edge=theta)
    counts = holonomy_counts(max_p)
    assert phases["N_R"] == counts["N_plus"]
    assert phases["N_L"] == counts["N_minus"]
    assert phases["phi_R"] == pytest.approx(counts["N_plus"] * theta)
    assert phases["phi_L"] == pytest.approx(counts["N_minus"] * theta)


def test_chirality_flux_equals_D_E():
    flux = chirality_flux_from_counts(100, 80)
    assert flux["D_E"] == 20
    assert flux["C_E"] == 20
    assert flux["N_R"] - flux["N_L"] == flux["D_E"]
    assert abs(flux["S_E"] - 20 / 180) < 1e-12


def test_link_phase_difference_to_D_E():
    theta = 0.5
    link = link_phase_difference_to_D_E(10.0, 20, theta_edge=theta)
    assert link["expected_D_E_times_theta"] == 10.0
    assert link["consistent_sign"] is True


def test_chiral_transport_report_structure():
    report = chiral_transport_report(30_000)
    assert report["helicity_map"]["ABCEA"]["channel"] == "R"
    assert report["helicity_map"]["CEABC"]["channel"] == "L"
    assert report["chirality_flux"]["D_E"] == report["holonomy_counts"]["D_E"]
    phi_diff = report["state_vector"]["observable_phase_difference"]
    d_e = report["chirality_flux"]["D_E"]
    theta = report["meta"]["theta_edge"]
    assert abs(phi_diff - d_e * theta) < 1e-9


def test_birefringent_velocity_pair():
    v_r, v_l = birefringent_velocity_pair(1.0, d_e_proxy=100.0, v0=1.0, alpha=0.01)
    assert v_r > v_l
    assert v_r == 2.0
    assert v_l > 0


def test_travel_time_birefringent():
    points = abcea_trajectory_xy()
    v_func, _ = make_velocity_func("chirality", d_e_global=50.0)
    result = travel_time_birefringent(points, v_func, d_e_global=50.0)
    assert result["T_R"] > 0
    assert result["T_L"] > 0


def test_brachistochrone_chirality_potential_has_birefringence():
    points = abcea_trajectory_xy()
    v_func, _ = make_velocity_func("chirality", d_e_global=100.0)
    t_single = travel_time_integral(points, v_func)
    biref = travel_time_birefringent(points, v_func, d_e_global=100.0)
    assert t_single > 0
    assert biref["T_R"] > 0
    assert biref["T_L"] > 0


def test_theta_edge_from_canonical_gaps():
    theta = theta_edge_from_gaps(CANONICAL_GAP_PATTERN)
    assert theta > 0
    assert theta == (math.pi / 2) * 3.0 / 12.0
