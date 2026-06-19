"""Sanity-Tests für Fibonacci-Zeta-Zeugen (kleine Schranken)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from eabc_zeta_fibonacci_witnesses import (
    ABCE_PATTERN,
    CEAB_PATTERN,
    LOG_PHI,
    fibonacci_below,
    fourgram_pattern,
    run_witnesses,
    witness_fibonacci_windows,
    witness_golden_resonance,
    witness_prime_phases,
    zeta_fibonacci,
)


def test_fibonacci_below_monotone():
    fib = fibonacci_below(1000)
    assert fib[0] == 0 and fib[1] == 1
    assert all(fib[i] <= fib[i + 1] for i in range(len(fib) - 1))
    assert all(fib[i] < fib[i + 1] for i in range(2, len(fib) - 1))


def test_fourgram_patterns():
    # künstliche Primliste mit bekannten mod-12-Resten
    # 5→A, 7→B, 11→C, 13→E
    primes = [5, 7, 11, 13, 17]
    assert fourgram_pattern(primes, 0) == "ABCE"
    # 11→C, 13→E, 17→A (17%12=5), 19→B
    primes2 = [11, 13, 17, 19]
    assert fourgram_pattern(primes2, 0) == "CEAB"


def test_zeta_fibonacci_finite():
    s = 2.0 + 0.0j
    val = zeta_fibonacci(s, n_terms=10)
    assert val.real > 0
    assert abs(val.imag) < 1e-10


def test_witness_prime_phases_small():
    w = witness_prime_phases(prime_bound=500, n_bins=10)
    assert w["n_primes"] > 10
    assert len(w["histogram_counts"]) == 10
    assert sum(w["histogram_counts"]) == w["n_primes"]
    assert 0.4 < w["mean_theta"] < 0.6


def test_witness_fibonacci_windows_small():
    w = witness_fibonacci_windows(prime_bound=2000)
    assert w["n_windows"] >= 3
    assert w["aggregate"]["total_ABCE"] + w["aggregate"]["total_CEAB"] >= 0


def test_witness_golden_resonance_structure():
    w = witness_golden_resonance(k_max=3, n_random=20, t_max=30.0, seed=0, zeta_f_terms=20)
    assert len(w["resonance_points"]) == 3
    t1 = 2 * np.pi / LOG_PHI
    assert abs(w["resonance_points"][0]["t_k"] - t1) < 1e-12
    if w["backend_zeta"] == "mpmath":
        assert w["zeta_summary"]["at_resonance"]["count"] == 3


def test_run_witnesses_json_serializable(tmp_path):
    out = tmp_path / "witnesses.json"
    data = run_witnesses(
        k_max=4,
        n_random=10,
        t_max=20.0,
        prime_bound=1000,
        n_bins=5,
        seed=1,
        zeta_f_terms=15,
    )
    out.write_text(json.dumps(data))
    loaded = json.loads(out.read_text())
    assert "witness_1_golden_resonance" in loaded
    assert loaded["epistemic_label"].startswith("Experiment")
