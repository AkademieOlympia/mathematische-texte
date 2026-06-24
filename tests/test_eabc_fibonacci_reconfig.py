"""Tests für collatz_eabc_fibonacci_reconfig_test.py (§20.5 Experiment)."""

from __future__ import annotations

import json
from pathlib import Path

from collatz_eabc_fibonacci_reconfig_test import fibonacci_up_to, run_scan


def test_fibonacci_up_to():
    assert fibonacci_up_to(0) == []
    assert fibonacci_up_to(1) == [1, 1]
    assert fibonacci_up_to(10) == [1, 1, 2, 3, 5, 8]
    assert fibonacci_up_to(20) == [1, 1, 2, 3, 5, 8, 13]


def test_run_scan_structure():
    report = run_scan(max_n=50, window=1, seed=0)
    assert report["max_n"] == 50
    assert "comparison" in report
    assert "fibonacci_samples" in report
    assert report["comparison"]["fibonacci_window"]["count"] >= 1


def test_run_scan_json_serializable(tmp_path: Path):
    report = run_scan(max_n=30, window=0, seed=1)
    out = tmp_path / "fib.json"
    out.write_text(json.dumps(report), encoding="utf-8")
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["label"].startswith("Experiment")
