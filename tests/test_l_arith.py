from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_l_arith_test import (
    GrammarRules,
    count_grammar_words,
    find_minimal_counterexamples,
    generate_grammar_words,
    is_grammar_valid_word,
    realized_words_by_length,
    run_suite,
)


def test_no_bb_in_grammar():
    rules = GrammarRules()
    for w in generate_grammar_words(6, rules):
        assert "BB" not in w


def test_be_grammar_valid():
    assert is_grammar_valid_word("BE")


def test_count_matches_generate_small():
    k = 5
    assert count_grammar_words(k) == len(generate_grammar_words(k))


def test_minimal_counterexample_be():
    buckets = realized_words_by_length(4, limit=10_000)
    hits = find_minimal_counterexamples(buckets)
    assert hits and hits[0]["word"] == "BE" and hits[0]["length"] == 2


def test_run_suite_smoke():
    result = run_suite([5], limit=2000, max_enumerate=10_000)
    assert result["lengths"][0]["grammar_enumerated"]
    assert result["lengths"][0]["grammar_count"] > 0


def test_full_enum_k4_ratio():
    result = run_suite([4], limit=50_000, max_enumerate=10_000)
    row = result["lengths"][0]
    assert row["grammar_count"] == 89
    assert row["grammar_enumerated"]
    assert 0 < row["ratio"] < 1
    assert result["minimal_counterexamples"][0]["word"] == "BE"


def test_sampling_mode_k10():
    result = run_suite([10], limit=5000, max_enumerate=100)
    row = result["lengths"][0]
    assert not row["grammar_enumerated"]
    assert row["notes"]
