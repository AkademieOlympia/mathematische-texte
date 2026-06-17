"""Tests für Z[i]-Faktor → EABC-Verteilung (collatz_eabc_gauss_faktor_eabc_test.py)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_gauss_faktor_eabc_test import (
    EABC_RESIDUES,
    classify_split_prime,
    eabc_label,
    factor_distribution_report,
    gaussian_factor_pair,
    run,
)
from eabc_from_lean import class_of


def test_gaussian_factor_pair_known_values():
    assert gaussian_factor_pair(5) == (1, 2)
    assert gaussian_factor_pair(13) == (2, 3)
    assert gaussian_factor_pair(17) == (1, 4)
    assert gaussian_factor_pair(3) is None
    assert gaussian_factor_pair(7) is None


def test_eabc_label_only_on_quadrant_residues():
    assert eabc_label(1) == "E"
    assert eabc_label(5) == "A"
    assert eabc_label(2) is None
    assert eabc_label(4) is None


def test_split_prime_row_structure():
    row = classify_split_prime(17)
    assert row is not None
    assert row.p == 17
    assert row.a == 1 and row.b == 4
    assert row.kappa == "A"
    assert row.a_eabc == "E"
    assert row.b_eabc is None


def test_no_both_eabc_visible_for_split_primes_up_to_5000():
    report = factor_distribution_report(5000)
    assert report["counts"]["both_eabc_visible"] == 0
    assert report["counts"]["split_primes_p_gt_3"] > 0


def test_kappa_odd_leg_null_not_extreme_at_10k():
    report = factor_distribution_report(10_000)
    z = abs(report["asymmetry"]["z_score_kappa_vs_odd_leg"])
    assert z < 5.0, "κ vs odd leg should not deviate wildly from shuffle null"


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "gauss_faktor.json"
    report = run(max_p=500, output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["counts"]["split_primes_p_gt_3"] == report["counts"]["split_primes_p_gt_3"]


def test_eabc_residues_match_lean():
    for r in EABC_RESIDUES:
        assert class_of(r) is not None
