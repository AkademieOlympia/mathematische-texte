"""Smoke-Tests für A-Absorption / A-Transfer (§4.8.26)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from eabc_a_absorption import (
    absorb_exact_search,
    absorb_greedy_log,
    absorption_defect,
    apply_a_comp,
    compression_grade,
    is_a_class_prime,
    multilayer_r_squared,
    run_a_absorption_analysis,
    simple_shell_r_squared,
)


def test_is_a_class_prime():
    assert is_a_class_prime(5)
    assert is_a_class_prime(17)
    assert not is_a_class_prime(7)
    assert not is_a_class_prime(3)


def test_multilayer_r_squared_basic():
    assert multilayer_r_squared(2, 2, [(5, 1)]) == 4 + 9 + 5


def test_absorption_defect_zero_exact():
    r = simple_shell_r_squared(3, 2)
    n, m, rt, d = absorb_exact_search(r, 6, 6)
    assert d == 0
    assert rt == r


def test_apply_a_comp_5_power():
    out = apply_a_comp(2, 2, [(5, 1)], exp_max=8)
    assert out["k_multi"] == 3
    assert out["a_eabc_classes"] == ["A"]
    assert out["a_comp"]["absorption_defect"] >= 0
    assert out["a_comp"]["compression_grade"] == 1.5


def test_greedy_log_runs():
    n, m, rt, d, method = absorb_greedy_log(2, 2, [(5, 1)])
    assert method == "greedy_log"
    assert rt == simple_shell_r_squared(n, m)
    assert d == absorption_defect(multilayer_r_squared(2, 2, [(5, 1)]), rt)


def test_compression_grade():
    assert compression_grade(3) == 1.5
    assert compression_grade(4) == 2.0


def test_run_smoke():
    report = run_a_absorption_analysis(
        n_max=2,
        m_max=2,
        l_max=1,
        exp_max=8,
        prime_limit=2_000,
        epsilon=2.0,
        seed=0,
    )
    assert report["section"] == "4.8.26"
    assert report["n_shells"] > 0
    assert report["summary"]["all_a_terms_in_A_class"]
    for sh in report["shells"]:
        assert "a_comp" in sh
        assert sh["a_comp"]["compression_grade"] >= 1.0
