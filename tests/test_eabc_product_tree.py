"""Tests für EABC-Produktbaum-Stub (collatz_eabc_product_tree_stub.py)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_product_tree_stub import (
    catalan,
    product_tree_report,
    run,
    tree_counts_for_n,
    unordered_compositions,
)


def test_catalan():
    assert catalan(0) == 1
    assert catalan(1) == 1
    assert catalan(2) == 2
    assert catalan(3) == 5


def test_unordered_compositions_12():
    comps = unordered_compositions(12)
    assert (2, 6) in comps
    assert (3, 4) in comps
    assert (2, 2, 3) in comps


def test_tree_counts_binary_catalan_one():
    tc = tree_counts_for_n(12)
    assert tc["Z_fact"] == 2
    assert tc["Z_tree_eff_H"] == 3
    assert all(c["catalan_trees"] == 1 for c in tc["compositions"] if c["k"] == 2)


def test_product_tree_report_structure():
    report = product_tree_report(max_n=20)
    assert "hurwitz_quaternion" in report
    assert "octonion_theoretical" in report
    assert "contrast" in report
    assert report["hurwitz_quaternion"]["algebra"] == "H"
    assert report["hurwitz_quaternion"]["Z_EABC_equals_Z_fact_fraction"] >= 0.0


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "pt.json"
    report = run(max_n=15, output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["meta"]["max_n"] == 15
    assert report["output_path"] == str(out)
