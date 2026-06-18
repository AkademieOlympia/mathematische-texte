"""Tests für Dirichlet-Erzeuger D̂(s) (collatz_eabc_dirichlet_D.py)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_dirichlet_D import (
    bernoulli_B,
    bernoulli_comparison,
    collect_D_terms,
    dirichlet_D_report,
    dirichlet_partial,
    repackaging_assessment,
    run,
    zeta_bernoulli_value,
    zeta_partial,
)


def test_bernoulli_B_known_values():
    assert bernoulli_B(2) is not None
    assert float(bernoulli_B(2)) == 1 / 6
    assert float(bernoulli_B(4)) == -1 / 30


def test_zeta_bernoulli_bridge():
    assert abs(zeta_bernoulli_value(1) - (-1 / 12)) < 1e-12
    assert abs(zeta_bernoulli_value(2) - (1 / 120)) < 1e-12


def test_collect_D_terms_rolling():
    terms = collect_D_terms(25, "rolling", rolling_window=3)
    assert len(terms) >= 20
    assert all(t.D >= 0 for t in terms)
    assert terms[0].n >= 2


def test_dirichlet_partial_s_zero():
    terms = collect_D_terms(20, "mu_infinity", rolling_window=3)
    s0 = dirichlet_partial(terms, 0.0)
    assert s0 == sum(t.D for t in terms)


def test_zeta_partial_matches_hand():
    assert abs(zeta_partial(5, 2.0, start_n=1) - (1 + 1 / 4 + 1 / 9 + 1 / 16 + 1 / 25)) < 1e-12


def test_bernoulli_comparison_keys():
    terms = collect_D_terms(30, "rolling", rolling_window=3)
    cmp = bernoulli_comparison(terms)
    assert "-2" in cmp
    assert "-4" in cmp
    assert "D_hat_N" in cmp["-2"]
    assert "B_2m" in cmp["-2"]
    assert cmp["-2"]["B_2m"]["2m"] == 2


def test_dirichlet_D_report_structure():
    report = dirichlet_D_report(max_n=30, rolling_window=3)
    assert report["meta"]["hypothesis_doc"].endswith("§13")
    assert "rolling" in report["by_I_ref"]
    assert "mu_infinity" in report["by_I_ref"]
    rolling = report["by_I_ref"]["rolling"]
    assert "growth_by_s" in rolling
    assert "bernoulli_exploration" in rolling
    assert "repackaging" in rolling
    assert len(rolling["first_terms"]) >= 10
    s_vals = {row["s"] for row in rolling["growth_by_s"]}
    assert 2.0 in s_vals
    assert -2.0 in s_vals
    assert -4.0 in s_vals


def test_repackaging_assessment_has_verdict():
    terms = collect_D_terms(40, "rolling", rolling_window=3)
    growth = [
        {"s": 2.0, "D_hat_over_zeta_partial": 0.5},
        {"s": 1.0, "D_hat_over_zeta_partial": 0.51},
        {"s": 0.0, "D_hat_N": sum(t.D for t in terms)},
    ]
    assess = repackaging_assessment(terms, growth)
    assert "verdict" in assess
    assert "pearson_D_vs_omega" in assess


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "dirichlet_D.json"
    report = run(max_n=25, output=out, rolling_window=3)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["meta"]["max_n"] == 25
    assert report["output_path"] == str(out)
