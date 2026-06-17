from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_phi_pref_test import (
    Z0,
    phi_pref,
    quadruplet_word,
    run_test,
    unit_phase,
)


def test_chirality_unit_phases():
    assert unit_phase("ABCE") == 1j
    assert unit_phase("CEAB") == -1j


def test_quadruplet_words_match_chirality():
    assert quadruplet_word(5) == "ABCE"
    assert quadruplet_word(11) == "CEAB"


def test_phi_pref_abce_ceab_separated():
    z_abce, t_abce = phi_pref("ABCE")
    z_ceab, t_ceab = phi_pref("CEAB")
    assert t_abce == 4.0
    assert t_ceab == 4.0
    assert z_abce != z_ceab
    assert z_abce.imag > z_ceab.imag


def test_run_test_two_tubes_at_1e6():
    _, summary = run_test(10**6, 2000, seed=42)
    assert summary.quadruplet_count == 166
    assert summary.abce_count == 84
    assert summary.ceab_count == 82
    assert summary.real_clusters_in_two_tubes
    assert summary.inter_tube_distance > 0.1
    assert summary.quadruplet_phase_shell_hits["+i"] == summary.abce_count
    assert summary.quadruplet_phase_shell_hits["-i"] == summary.ceab_count


def test_z0_anchor():
    assert Z0 == complex(-1.0, -1.0 / 12.0)
