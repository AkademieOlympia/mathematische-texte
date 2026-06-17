"""Tests für Eisenstein–EABC-Spaltung mit glatt-EABC (collatz_eabc_eisenstein_spaltung_test.py)."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_eisenstein_spaltung_test import (
    GAMMA_PAIRS,
    classify_eisenstein_split_prime,
    coarse_mod3_defekt_report,
    eisenstein_factor_pair,
    eisenstein_norm,
    eisenstein_split_class,
    run,
    spaltung_report,
)
from eabc_from_lean import class_of


def test_eisenstein_norm_basic():
    assert eisenstein_norm(1, 3) == 7
    assert eisenstein_norm(1, 4) == 13
    assert eisenstein_norm(1, 6) == 31


def test_eisenstein_factor_pair_known():
    assert eisenstein_factor_pair(7) == (1, 3)  # lex. min. unter a ≤ b
    assert eisenstein_factor_pair(13) == (1, 4)
    assert eisenstein_factor_pair(5) is None  # 5 ≡ 2 mod 3, inert
    assert eisenstein_factor_pair(3) is None  # ramified


def test_eisenstein_split_class():
    assert eisenstein_split_class(7) == "split"
    assert eisenstein_split_class(5) == "inert"
    assert eisenstein_split_class(3) == "ramified"


def test_classify_eisenstein_gamma_both_legs():
    row = classify_eisenstein_split_prime(7)
    assert row is not None
    assert row.a == 1 and row.b == 3
    assert row.gamma[0] in {"E", "A", "B", "C"}
    assert row.gamma[1] in {"E", "A", "B", "C"}


def test_mod3_coarse_bipartite_exact_at_10k():
    report = coarse_mod3_defekt_report(10_000)
    assert report["exact_coarse_bipartite"] is True
    assert report["mismatches"] == 0


def test_split_primes_have_eb_kappa():
    for p in [7, 13, 31, 37, 43]:
        assert eisenstein_split_class(p) == "split"
        ec = class_of(p)
        assert ec is not None
        assert ec.value in {"E", "B"}


def test_inert_primes_have_ac_kappa():
    for p in [5, 11, 17, 23, 29]:
        assert eisenstein_split_class(p) == "inert"
        ec = class_of(p)
        assert ec is not None
        assert ec.value in {"A", "C"}


def test_all_gamma_classes_reachable_by_10k():
    report = spaltung_report(10_000)
    observed = sum(1 for g in GAMMA_PAIRS if report["counts"][f"({g[0]},{g[1]})"] > 0)
    assert observed >= 8


def test_mu_sums_to_one():
    report = spaltung_report(10_000)
    total = sum(report["mu_X"].values())
    assert math.isclose(total, 1.0, rel_tol=1e-9)


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "eisenstein_spaltung.json"
    report = run(max_ps=[500], output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert "scales" in loaded
    assert loaded["scales"]["500"]["split_count"] == report["scales"]["500"]["split_count"]
