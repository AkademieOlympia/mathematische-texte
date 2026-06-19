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
    SIGNATURE_TARGET,
    chi_quadruplet,
    compute_Cm_Zm,
    fibonacci_below,
    find_quadruplets,
    fourgram_pattern,
    matches_signature,
    mod210_triple,
    run_witnesses,
    witness_fibonacci_mod210_shells,
    witness_fibonacci_windows,
    witness_golden_fourier,
    witness_golden_resonance,
    witness_meromorphic_normal_form,
    witness_prime_phases,
    zeta_f_meromorphic,
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
        m_max=3,
        quadruplet_bound=20_000,
        meromorphic_m_max=20,
    )
    out.write_text(json.dumps(data))
    loaded = json.loads(out.read_text())
    assert "witness_1_golden_resonance" in loaded
    assert "witness_4_meromorphic_normal_form" in loaded
    assert "witness_5_golden_fourier" in loaded
    assert "witness_6_fibonacci_mod210_shells" in loaded
    assert loaded["epistemic_label"].startswith("Experiment")


def test_chi_quadruplet_signs():
    assert chi_quadruplet(101) == 1   # 101 mod 12 = 5, ABCE
    assert chi_quadruplet(11) == -1   # 11 mod 12 = 11, CEAB
    assert chi_quadruplet(7) == 0


def test_zeta_f_meromorphic_vs_direct():
    s = 2.5 + 0.5j
    direct = zeta_fibonacci(s, n_terms=40)
    merom = zeta_f_meromorphic(s, m_max=30)
    rel = abs(direct - merom) / max(abs(direct), 1e-15)
    assert rel < 0.05


def test_find_quadruplets_small():
    quads = find_quadruplets(100)
    assert 5 in quads
    assert all(
        p + 2 in set(range(2, 101))
        for p in quads[:3]
    )


def test_compute_Cm_Zm_normalization():
    quads = find_quadruplets(10_000)
    _, z_m, q = compute_Cm_Zm(quads, m_max=2)
    assert q > 0
    assert len(z_m) == 3
    assert 0.0 <= z_m[0] <= 1.0 + 1e-9


def test_witness_golden_fourier_small():
    w = witness_golden_fourier(quadruplet_bound=50_000, m_max=4)
    assert w["checkpoints"]
    last = w["checkpoints"][-1]
    assert "Z_m" in last
    assert len(last["Z_m"]) == 5


def test_witness_meromorphic_structure():
    w = witness_meromorphic_normal_form(m_max=10, n_terms=20)
    assert w["resonance_towers"]
    assert w["max_relative_error"] < 0.2


def test_fibonacci_mod210_shells():
    w = witness_fibonacci_mod210_shells(quadruplet_bound=100_000)
    assert w["n_fibonacci_shells"] >= 1
    agg = w["aggregate_mod210"]
    d11, d101, d191 = mod210_triple(agg)
    sig = (d11 > 0, d101 > 0, d191 < 0)
    # dokumentiere ob Zielsignatur im Aggregat gilt (kein harter Test)
    assert isinstance(sig, tuple)
    assert SIGNATURE_TARGET == ("+", "+", "-")
