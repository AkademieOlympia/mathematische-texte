from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_forbidden_words import compute_forbidden_at_length, run_forbidden_catalog
from collatz_l_arith_test import realized_words_by_length


def test_f2_contains_be():
    buckets = realized_words_by_length(4, limit=10_000)
    rep = compute_forbidden_at_length(2, buckets[2])
    assert "BE" in rep.forbidden_words
    assert rep.hero == "BE"
    assert rep.forbidden_count >= 1


def test_forbidden_subset_of_grammar():
    buckets = realized_words_by_length(5, limit=5000)
    rep = compute_forbidden_at_length(3, buckets[3])
    assert rep.forbidden_count == rep.grammar_count - rep.arith_count


def test_run_catalog_smoke():
    result = run_forbidden_catalog([2, 3], limit=2000)
    assert len(result["F_n"]) == 2
    f2 = next(r for r in result["F_n"] if r["n"] == 2)
    assert f2["hero"] == "BE"
    assert result["pipeline"]
