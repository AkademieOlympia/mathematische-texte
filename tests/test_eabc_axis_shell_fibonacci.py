"""Smoke-Tests für Achsenschalen κ_F und C_F (§4.8.25)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from eabc_axis_shell_fibonacci import (
    axis_shell_r_squared,
    enrich_shell_fibonacci,
    fibonacci_table,
    generate_axis_shells,
    kappa_F,
    min_fibonacci_summands,
    run_axis_shell_analysis,
)


def test_axis_shell_r_squared_classic():
    assert axis_shell_r_squared([(2, 3), (3, 2)]) == 8 + 9


def test_fibonacci_exact_small():
    fibs = fibonacci_table(100)
    # 8 = F_6 = 8 (F_1=1,F_2=1,F_3=2,F_4=3,F_5=5,F_6=8)
    s, total = min_fibonacci_summands(8, fibs, eps_rel=0.0)
    assert s == 1
    assert total == 8


def test_kappa_F_compression():
    assert kappa_F(3, 1) == 3.0
    assert kappa_F(4, 2) == 2.0
    assert kappa_F(3, 3) == 1.0


def test_generate_k2_shells_nonempty():
    shells = generate_axis_shells(2, exp_max=2)
    assert len(shells) == 2 * 2  # min_exp=1, exp_max=2
    assert all(s["k"] == 2 for s in shells)


def test_enrich_marks_compressible():
    fibs = fibonacci_table(500)
    sh = {"k": 3, "r_squared": 8, "R": 2.0, "terms": [(2, 3), (3, 1), (5, 1)]}
    out = enrich_shell_fibonacci(sh, fibs, eps_rel=0.0)
    assert out["s_F"] >= 1
    assert out["kappa_F"] >= 1.0


def test_run_axis_shell_smoke():
    report = run_axis_shell_analysis(
        k_values=(2, 3),
        exp_max=2,
        prime_limit=2_000,
        epsilon=2.0,
        n_perm=20,
        seed=0,
    )
    assert report["section"] == "4.8.25"
    assert report["n_shells_total"] > 0
    assert report["n_shells_eligible"] >= 0
    corr = report["correlation"]["kappa_F_vs_asymmetry_L1"]
    assert "pearson" in corr
    assert len(report["shells"]) == report["n_shells_total"]
    for row in report["shells"]:
        assert "C_F" in row
        assert "kappa_F" in row
        assert row["kappa_F"] >= 1.0 - 1e-9
