"""Tests für EABC Bohm-/AB-/Berry-Potential-Phase-Stubs."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_chirale_transport import holonomy_phase_difference
from collatz_eabc_potential_phase import (
    aharonov_bohm_phase,
    berry_phase_difference,
    berry_phase_from_loop,
    bohm_like_velocity,
    default_ab_edge_field,
    potential_phase_report,
)
from collatz_eabc_transition_graph import ABCEA_WORD, CEABC_WORD


def test_bohm_like_velocity_increases_with_gradient():
    v0 = bohm_like_velocity(0.0)
    v1 = bohm_like_velocity(10.0)
    assert v1 > v0
    assert v0 == 1.0
    assert v1 == pytest.approx(1.1)


def test_default_ab_edge_field_has_four_edges():
    a_field = default_ab_edge_field(theta_edge=0.5)
    assert set(a_field.keys()) == {"AB", "BC", "CE", "EA"}
    assert all(v == 0.5 for v in a_field.values())


def test_aharonov_bohm_phase_signs():
    theta = math.pi / 4
    a_field = default_ab_edge_field(theta_edge=theta)
    phi_abcea = aharonov_bohm_phase(ABCEA_WORD, a_field)
    phi_ceabc = aharonov_bohm_phase(CEABC_WORD, a_field)
    assert phi_abcea == 4 * theta
    assert phi_ceabc == -4 * theta
    assert phi_abcea == -phi_ceabc


def test_berry_phase_difference_matches_holonomy():
    theta = math.pi / 4
    for loop in (ABCEA_WORD, CEABC_WORD):
        berry = berry_phase_from_loop(loop, theta_edge=theta)
        assert berry["berry_phase_difference"] == berry_phase_difference(
            berry["phi_R"], berry["phi_L"]
        )
        assert (
            berry["holonomy_phase_difference"]
            == holonomy_phase_difference(loop, theta_edge=theta)
        )


def test_berry_phase_difference_direct():
    assert berry_phase_difference(3.0, 1.0) == 2.0
    assert berry_phase_difference(0.0, 0.0) == 0.0


def test_potential_phase_report_structure():
    report = potential_phase_report(theta_edge=0.25)
    assert report["meta"]["theory"] == "collatz_eabc_potential_geometrie.md"
    assert ABCEA_WORD in report["loops"]
    assert CEABC_WORD in report["loops"]
    assert report["boxed"]["not_c"].startswith("Effektive")
