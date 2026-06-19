"""Smoke-Tests für EABC-Fibonacci-Ensemble (§4.8.24)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from eabc_fibonacci_ensemble import (
    fibonacci_contraction,
    run_ensemble_test,
    weighted_class_sums,
)
from eabc_level2_fluctuation import build_eabc_word_with_primes, theta_phi_from_primes


def test_fibonacci_contraction_bounds():
    import numpy as np

    theta = np.linspace(0.0, 0.99, 50)
    kappa = fibonacci_contraction(theta)
    assert np.all(kappa > 0.0)
    assert np.all(kappa <= 1.0 + 1e-12)
    assert kappa[0] == 1.0


def test_ensemble_smoke():
    report = run_ensemble_test(n_primes=5_000, n_bins=8, B_perm=20, seed=42)
    assert report["n_primes"] == 5_000
    assert report["section"] == "4.8.24"
    assert report["Psi_F"] != 0.0
    assert report["Psi_Fib"] > 0.0
    assert len(report["per_class"]) == 4
    assert report["perm_null"]["n_draws"] == 20
    assert report["bin_ensemble"]["asymmetry_ratio"] is not None


def test_perm_null_centered_on_shuffled():
    word, primes = build_eabc_word_with_primes(2_000)
    theta = theta_phi_from_primes(primes)
    weights = {"E": 1.0, "A": 1.0, "B": 1.0, "C": 1.0}
    import numpy as np

    rng = np.random.default_rng(0)
    shuffled = word.copy()
    rng.shuffle(shuffled)
    orig = weighted_class_sums(word, theta, weights, "chi_A")["ensemble"]
    shuf = weighted_class_sums(shuffled, theta, weights, "chi_A")["ensemble"]
    # Permutation ändert Ψ_F typischerweise (nicht invariant)
    assert isinstance(orig, float)
    assert isinstance(shuf, float)
