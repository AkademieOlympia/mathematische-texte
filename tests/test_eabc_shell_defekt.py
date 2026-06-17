"""Tests für Σ→p-Schalen-Defekt-Experiment (collatz_eabc_shell_defekt_test.py)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_shell_defekt_test import (
    compute_shell_invariants,
    defect_magnitude,
    is_prime,
    multi_scale_summary,
    omega_distinct,
    pick_best_iref,
    prime_vs_composite_comparison,
    run,
    shell_defekt_report,
    tau_divisor_count,
)


def test_is_prime_basic():
    assert is_prime(2)
    assert is_prime(7)
    assert not is_prime(1)
    assert not is_prime(9)
    assert not is_prime(15)


def test_omega_distinct():
    assert omega_distinct(12) == 2
    assert omega_distinct(7) == 1
    assert omega_distinct(30) == 3


def test_tau_divisor_count():
    assert tau_divisor_count(1) == 1
    assert tau_divisor_count(6) == 4
    assert tau_divisor_count(7) == 2


def test_compute_shell_invariants_prime():
    inv = compute_shell_invariants(7)
    assert inv is not None
    assert inv.is_prime
    assert inv.shell_size == 192
    assert inv.H_n > 0
    assert inv.omega == 1


def test_compute_shell_invariants_composite():
    inv = compute_shell_invariants(6)
    assert inv is not None
    assert not inv.is_prime
    assert inv.shell_size == 96
    assert inv.omega == 2


def test_defect_magnitude_zero():
    assert defect_magnitude(0.0, 0.0, 0.0) == 0.0


def test_shell_defekt_report_structure():
    report = shell_defekt_report(max_n=30, rolling_window=3)
    assert report["meta"]["shell_count"] >= 20
    assert "rows" in report
    assert "prime_vs_composite" in report
    assert "best_I_ref" in report
    assert "bernoulli_D_correlation" in report
    row = report["rows"][0]
    assert "D_rolling" in row
    assert "D_omega" in row
    assert "D_tau" in row
    assert "D_mu_infinity" in row
    assert "V_n_bernoulli_proxy" in row
    assert "magnitude" in row["D_rolling"]
    pvc = report["prime_vs_composite"]["rolling"]
    assert "mean_abs_D_prime" in pvc
    assert "mean_abs_D_composite" in pvc
    assert "epistemic_note" in pvc
    assert "I_ref_variants" in report["meta"]


def test_prime_vs_composite_comparison():
    rows = [
        {"is_prime": True, "D_rolling": {"magnitude": 0.5}, "n": 7, "omega": 1},
        {"is_prime": False, "D_rolling": {"magnitude": 0.1}, "n": 6, "omega": 2},
        {"is_prime": True, "D_rolling": {"magnitude": 0.3}, "n": 5, "omega": 1},
        {"is_prime": False, "D_rolling": {"magnitude": 0.2}, "n": 8, "omega": 1},
    ]
    cmp = prime_vs_composite_comparison(rows)
    assert cmp["prime_count"] == 2
    assert cmp["composite_count"] == 2
    assert cmp["primes_larger_on_mean"]


def test_pick_best_iref_ranking():
    comparisons = {
        "rolling": {"ratio_mean_prime_over_composite": 0.9, "primes_larger_on_mean": False},
        "omega": {"ratio_mean_prime_over_composite": 1.2, "primes_larger_on_mean": True},
    }
    best = pick_best_iref(comparisons)
    assert best["best_ratio_name"] == "omega"
    assert best["candidate"] == "omega"
    assert best["stable_prime_anomaly"]


def test_multi_scale_summary():
    summary = multi_scale_summary([25, 30], rolling_window=3)
    assert "25" in summary["scales"]
    assert "30" in summary["scales"]


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "shell_defekt.json"
    report = run(max_n=25, output=out, include_multi_scale=False)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["meta"]["max_n"] == 25
    assert report["output_path"] == str(out)
