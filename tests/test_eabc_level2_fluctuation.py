"""Tests für EABC Level-2-Fluktuationsgeometrie (Δ_F auf Λ²(ℝ⁴))."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from eabc_level2_fluctuation import (
    CHECKPOINTS,
    analyze_checkpoint,
    build_eabc_word,
    collect_hl_null_vectors,
    collect_markov_null_vectors,
    compute_a_vector,
    delta_F,
    empirical_covariance,
    run_fluctuation_test,
)


def test_compute_a_vector_bounds():
    word = build_eabc_word(500)
    a = compute_a_vector(word[:200])
    assert a.shape == (6,)
    assert np.all(np.abs(a) <= 1.0 + 1e-12)
    assert not np.any(np.isnan(a))


def test_delta_F_zero_for_identical():
    rng = np.random.default_rng(0)
    vecs = rng.normal(size=(30, 6))
    _, sigma = empirical_covariance(vecs)
    assert delta_F(sigma, sigma) == 0.0


def test_delta_F_positive_for_different():
    rng = np.random.default_rng(1)
    _, s1 = empirical_covariance(rng.normal(scale=1.0, size=(40, 6)))
    _, s2 = empirical_covariance(rng.normal(scale=2.0, size=(40, 6)))
    assert delta_F(s1, s2) > 0.05


def test_analyze_checkpoint_small():
    word = build_eabc_word(5_000)
    rng = np.random.default_rng(42)
    row = analyze_checkpoint(word, m=100, B_rand=5, rng=rng)
    assert row["m"] == 100
    assert row["K_prime"] >= 40
    assert len(row["Sigma_A_prime"]) == 6
    assert 0.0 <= row["delta_F_perm"] < 1.0
    assert 0.0 <= row["delta_F_markov"] < 1.0
    assert row["Delta_F"] == row["delta_F_perm"]
    assert row["delta_F_hl"] is None
    assert row["mu_A_prime_norm"] < 0.5


def test_markov_null_differs_from_perm():
    word = build_eabc_word(3_000)
    rng = np.random.default_rng(99)
    perm = collect_markov_null_vectors(word, m=200, B=10, rng=rng)
    # Markov-Resampling erzeugt gültige 6-Vektoren
    assert perm.shape[1] == 6
    assert perm.shape[0] >= 100


def test_hl_null_stub():
    word = build_eabc_word(500)
    rng = np.random.default_rng(0)
    try:
        collect_hl_null_vectors(word, m=100, B=2, rng=rng)
        assert False, "HL-Stub sollte NotImplementedError werfen"
    except NotImplementedError as exc:
        assert "Stufe 3" in str(exc)


def test_run_fluctuation_test_structure():
    report = run_fluctuation_test(
        n_primes=10_000,
        checkpoints=[200, 500],
        B_rand=3,
        seed=7,
    )
    assert report["n_primes"] == 10_000
    assert len(report["results"]) == 2
    for row in report["results"]:
        assert "Delta_F" in row
        assert "delta_F_perm" in row
        assert "delta_F_markov" in row
        assert row["delta_F_hl"] is None
        assert len(row["spec_prime"]) == 6


def test_checkpoints_default():
    assert CHECKPOINTS == [1_000, 2_000, 5_000, 10_000, 20_000]
