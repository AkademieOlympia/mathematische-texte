"""Tests für Gauß–EABC-Spaltung mit glatt-EABC (collatz_eabc_gauss_spaltung_test.py)."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_gauss_spaltung_test import (
    GAMMA_PAIRS,
    UNIFORM_MU,
    classify_split_prime,
    gaussian_factor_pair,
    kappa_glatt,
    run,
    spaltung_report,
    strip_smooth,
)
from eabc_from_lean import class_of


def test_strip_smooth_basic():
    assert strip_smooth(1) == (0, 0, 1)
    assert strip_smooth(12) == (2, 1, 1)
    assert strip_smooth(18) == (1, 2, 1)
    assert strip_smooth(5) == (0, 0, 5)


def test_kappa_glatt_coprime_to_6():
    alpha, beta, core, k = kappa_glatt(12)
    assert (alpha, beta, core) == (2, 1, 1)
    assert k == "E"
    alpha, beta, core, k = kappa_glatt(5)
    assert core == 5 and k == "A"


def test_gaussian_factor_pair_known():
    assert gaussian_factor_pair(5) == (1, 2)
    assert gaussian_factor_pair(13) == (2, 3)
    assert gaussian_factor_pair(7) is None


def test_classify_split_prime_gamma_both_legs():
    row = classify_split_prime(17)
    assert row is not None
    assert row.a == 1 and row.b == 4
    assert row.a_prime == 1 and row.b_prime == 1
    assert row.gamma == ("E", "E")


def test_all_gamma_classes_reachable_by_10k():
    report = spaltung_report(10_000)
    observed = sum(1 for g in GAMMA_PAIRS if report["counts"][f"({g[0]},{g[1]})"] > 0)
    assert observed >= 8, "mindestens halbe Γ-Klassen sollten bis 10k vorkommen"


def test_mu_sums_to_one():
    report = spaltung_report(10_000)
    total = sum(report["mu_X"].values())
    assert math.isclose(total, 1.0, rel_tol=1e-9)


def test_chi2_vs_shuffle_null_not_extreme_at_100k():
    """Hohe χ² vs 1/16 kann von Randverteilungen stammen; Shuffle-Null prüft Kopplung."""
    report = spaltung_report(100_000)
    z = abs(report["z_score_chi2_vs_shuffle"])
    assert z < 3.0


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "spaltung.json"
    report = run(max_ps=[500], output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert "scales" in loaded
    assert loaded["scales"]["500"]["split_count"] == report["scales"]["500"]["split_count"]


def test_glatt_core_always_has_kappa():
    for p in [5, 13, 17, 29, 37]:
        row = classify_split_prime(p)
        assert row is not None
        assert class_of(row.a_prime) is not None
        assert class_of(row.b_prime) is not None
