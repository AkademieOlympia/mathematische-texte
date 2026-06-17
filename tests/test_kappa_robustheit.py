from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_kappa_robustheit import (
    KappaId,
    bootstrap_kappa1,
    compare_variants,
    entropy_estimates,
    kappa1_letter_at,
    kappa2_letter_at,
    realized_words_variant,
    run_suite,
)


def test_kappa_letters_defined_on_e_class():
    assert kappa1_letter_at(1) == "E"
    assert kappa2_letter_at(1) is not None


def test_kappa1_matches_legacy_buckets():
    legacy = realized_words_variant(4, 5000, KappaId.KAPPA1)
    from collatz_l_arith_test import realized_words_by_length

    ref = realized_words_by_length(4, 5000)
    assert legacy[2] == ref[2]


def test_be_forbidden_kappa1_small():
    buckets = realized_words_variant(4, 20_000, KappaId.KAPPA1)
    assert "BE" not in buckets[2]


def test_compare_variants_smoke():
    result = compare_variants(list(KappaId), limit=5000, fn_lengths=[2, 3], ratio_ks=(4,), max_enumerate=50_000)
    assert len(result["variants"]) == 3
    assert "comparison" in result


def test_entropy_positive():
    est = entropy_estimates({2: 1, 3: 6, 4: 38})
    assert est["h_F_limsup_estimate"] is not None
    assert est["h_F_limsup_estimate"] > 0


def test_bootstrap_smoke():
    rows = bootstrap_kappa1([2000, 5000], fn_lengths=[2], ratio_ks=(4,), max_enumerate=10_000)
    assert len(rows) == 2
    assert rows[0].be_forbidden


def test_run_suite_smoke():
    result = run_suite(limit=3000, bootstrap_limits=[3000], max_enumerate=20_000)
    assert result["stage"] == "2B"
    assert "test1_kappa_comparison" in result
    assert len(result["test1_kappa_comparison"]["variants"]) == 3
