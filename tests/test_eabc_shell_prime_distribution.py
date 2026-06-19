"""Smoke-Tests für EABC-Primverteilung auf Schalen (§4.8.21–4.8.22)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from eabc_shell_prime_distribution import (
    analyze_defect_distribution,
    analyze_distribution,
    assign_shell_defect_nearest,
    assign_shell_factor,
    assign_shell_log_grid,
    charge_observables,
    class_of,
    defect_delta,
    format_summary,
    smooth_shell_exponents,
)
import numpy as np


def test_prime_shell_factor_is_origin():
    assert assign_shell_factor(5) == (0, 0)
    assert assign_shell_factor(7) == (0, 0)
    assert assign_shell_factor(11) == (0, 0)


def test_smooth_shell_exponents():
    assert smooth_shell_exponents(72) == (3, 2)
    assert smooth_shell_exponents(12) == (2, 1)


def test_class_of_primes():
    assert class_of(5) == "A"
    assert class_of(7) == "B"
    assert class_of(11) == "C"
    assert class_of(13) == "E"


def test_log_grid_assigns_valid_pair():
    n, m = assign_shell_log_grid(101, n_max=6, m_max=6)
    assert 0 <= n <= 6
    assert 0 <= m <= 6


def test_defect_delta_nonnegative():
    delta, n, m, r = defect_delta(97, n_max=6, m_max=6)
    assert delta >= 0
    assert 0 <= n <= 6
    assert 0 <= m <= 6
    assert r > 0


def test_defect_nearest_assigns_shell():
    n, m, r, delta = assign_shell_defect_nearest(97, n_max=6, m_max=6)
    assert delta == defect_delta(97, n_max=6, m_max=6)[0]
    assert r > 0


def test_charge_observables_uniform():
    counts = np.array([25, 25, 25, 25], dtype=np.int64)
    obs = charge_observables(counts)
    assert obs["D_AB"] == 0.0
    assert obs["D_AC"] == 0.0
    assert abs(obs["D_E"]) < 1e-9


def test_charge_observables_D_E():
    counts = np.array([10, 9, 9, 9], dtype=np.int64)
    obs = charge_observables(counts)
    assert obs["D_E"] == 10 - (9 + 9 + 9) / 3


def test_analyze_distribution_factor_smoke():
    result = analyze_distribution(limit=2_000, shell_mode="factor", n_perm=49, seed=0)
    assert result["n_primes"] > 100
    assert result["n_buckets"] == 1
    assert result["buckets"][0]["n"] == 0
    assert result["buckets"][0]["m"] == 0
    g = result["global"]
    assert abs(sum(g["fractions"].values()) - 1.0) < 1e-9
    assert 0.0 < g["perm_p_chi2"] <= 1.0
    summary = format_summary(result)
    assert "mode=factor" in summary


def test_analyze_distribution_log_grid_smoke():
    result = analyze_distribution(
        limit=3_000, shell_mode="log-grid", n_max=5, m_max=5, n_perm=49, seed=1
    )
    assert result["n_buckets"] >= 2
    assert all(b["n_primes"] >= 1 for b in result["buckets"])


def test_analyze_defect_distribution_smoke():
    result = analyze_defect_distribution(
        limit=5_000, n_max=6, m_max=6, n_perm=49, seed=2, band_mode="nearest"
    )
    assert result["section"] == "4.8.22"
    assert result["shell_mode"] == "defect"
    assert result["n_primes"] > 100
    assert result["n_buckets"] >= 2
    g = result["global"]
    assert "charge" in g
    assert "D_AB" in g["charge"]
    assert "perm_p_D_AB" in g
    assert result["delta_stats"]["mean"] is not None


def test_analyze_defect_epsilon_band():
    result = analyze_defect_distribution(
        limit=3_000,
        n_max=5,
        m_max=5,
        n_perm=29,
        seed=3,
        band_mode="epsilon",
        epsilon=2.0,
    )
    assert result["band_mode"] == "epsilon"
    assert result["epsilon"] == 2.0
    assert result["n_primes"] >= 1


def test_analyze_defect_interval_band():
    result = analyze_defect_distribution(
        limit=4_000, n_max=4, m_max=4, n_perm=29, seed=4, band_mode="interval"
    )
    assert result["band_mode"] == "interval"
    assert result["n_buckets"] >= 1
