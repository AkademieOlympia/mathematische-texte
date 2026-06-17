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


def test_be_realized_kappa2_key_finding():
    """Stufe-2B-Kernbefund: κ₂ realisiert BE — kein kodierungsunabhängiger Zeuge."""
    buckets = realized_words_variant(4, 50_000, KappaId.KAPPA2)
    assert "BE" in buckets[2]


def test_kappa1_kappa3_align_on_be():
    limit = 30_000
    b1 = realized_words_variant(4, limit, KappaId.KAPPA1)
    b3 = realized_words_variant(4, limit, KappaId.KAPPA3)
    assert ("BE" not in b1[2]) == ("BE" not in b3[2])
    assert b1[2] == b3[2]


def test_compare_variants_be_not_robust():
    result = compare_variants(
        list(KappaId),
        limit=100_000,
        fn_lengths=[2, 3, 4, 5],
        ratio_ks=(4, 10),
        max_enumerate=200_000,
    )
    cmp_ = result["comparison"]
    assert cmp_["be_forbidden_all_variants"] is False
    assert cmp_["minimal_counterexample"]["kappa1_naive_mod12"]["word"] == "BE"
    k2_min = cmp_["minimal_counterexample"]["kappa2_nu2_rotate"]
    assert k2_min["word"] != "BE"
    assert k2_min["length"] >= 3
    by_var = {v["variant"]: v for v in result["variants"]}
    assert by_var["kappa1_naive_mod12"]["be_forbidden"] is True
    assert by_var["kappa2_nu2_rotate"]["be_forbidden"] is False


def test_kappa2_be_in_grammar_length2():
    """Länge 2: BE ∈ L und unter κ₂ realisiert (|F_2|=0 im Vollauf n≤10⁶)."""
    from collatz_l_arith_test import generate_grammar_words

    buckets = realized_words_variant(2, 50_000, KappaId.KAPPA2)
    assert set(generate_grammar_words(2)) <= buckets[2]
    assert "BE" in buckets[2]


def test_entropy_positive():
    est = entropy_estimates({2: 1, 3: 6, 4: 38})
    assert est["h_F_limsup_estimate"] is not None
    assert est["h_F_limsup_estimate"] > 0


def test_bootstrap_kappa1_monotone_limits():
    rows = bootstrap_kappa1([5000, 15_000], fn_lengths=[2, 3], ratio_ks=(4,), max_enumerate=50_000)
    assert len(rows) == 2
    assert all(r.be_forbidden for r in rows)
    assert rows[0].minimal_word == "BE"


def test_run_suite_structure():
    result = run_suite(limit=5000, bootstrap_limits=[5000], max_enumerate=30_000)
    assert result["stage"] == "2B"
    variants = result["test1_kappa_comparison"]["variants"]
    by_id = {v["variant"]: v for v in variants}
    assert by_id["kappa1_naive_mod12"]["be_forbidden"] is True
    assert by_id["kappa2_nu2_rotate"]["be_forbidden"] is False
    assert "κ-sensitiv" in result["test1_kappa_comparison"]["comparison"]["interpretation"]
