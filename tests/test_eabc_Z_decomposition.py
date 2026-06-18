"""Tests für EABC-Zerlegungsregimen-Experiment (collatz_eabc_Z_decomposition_test.py)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_Z_decomposition_test import (
    Z_EABC_count,
    Z_fact,
    Z_regime_count,
    build_cache,
    compute_shell_measure,
    delta_Z,
    run,
    unordered_factorizations,
    z_decomposition_report,
    z_row,
)


def test_Z_fact_basic():
    assert Z_fact(7) == 0
    assert Z_fact(6) == 1  # 2*3
    assert Z_fact(12) == 2  # 2*6, 3*4
    assert Z_fact(36) == 4


def test_unordered_factorizations():
    assert unordered_factorizations(7) == []
    assert (2, 3) in unordered_factorizations(6)
    assert len(unordered_factorizations(12)) == 2


def test_compute_shell_measure():
    mu = compute_shell_measure(7)
    assert mu is not None
    assert mu.shell_size == 192
    assert mu.H > 0
    assert mu.distinct_gamma > 0


def test_Z_EABC_prime_zero():
    cache = build_cache(20)
    assert Z_EABC_count(7, cache) == 0
    assert Z_EABC_count(6, cache) >= 1


def test_Z_regime_prime_one():
    cache = build_cache(20)
    assert Z_regime_count(7, cache) == 1
    assert Z_regime_count(6, cache) in (0, 1)


def test_delta_Z_defined():
    cache = build_cache(15)
    assert isinstance(delta_Z(7, cache), int)


def test_z_row_structure():
    cache = build_cache(20)
    row = z_row(6, cache)
    assert row["n"] == 6
    assert not row["is_prime"]
    assert row["Z_fact"] == 1
    assert "delta_Z" in row
    assert "Z_regime" in row


def test_z_decomposition_report_structure():
    report = z_decomposition_report(max_n=30)
    assert report["meta"]["shell_count"] >= 20
    assert "prime_vs_composite" in report
    assert "omega_correlation" in report
    assert "z_fact_vs_z_eabc" in report
    assert "rows" in report
    pvc = report["prime_vs_composite"]["delta_Z_abs"]
    assert "ratio_mean_prime_over_composite" in pvc
    assert "mean_prime" in pvc


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "z_decomp.json"
    report = run(max_n=25, output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["meta"]["max_n"] == 25
    assert report["output_path"] == str(out)
    assert report["verdict"]
