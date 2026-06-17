#!/usr/bin/env python3
"""Katalog verbotener Wörter F_n = L(n) \\ L_arith(n) (Stufe 2, Tao-Stil).

Für jedes n: alle grammatisch zulässigen Wörter der Länge n, die auf der
Suchtiefe ungerade n ≤ limit nicht als κ-Präfix realisiert werden.
Nullstellenkatalog-Stil (vgl. IEANTN-Tabellen); kein Collatz-Beweis.

Referenz: collatz_generalangriff_2026.md, collatz_l_arith_test.py.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

from collatz_l_arith_test import (
    FULL_ENUM_MAX_K,
    GrammarRules,
    analyze_length,
    count_grammar_words,
    generate_grammar_words,
    realized_words_by_length,
)

DEFAULT_LENGTHS = list(range(2, 9))  # n = 2, …, 8
DEFAULT_LIMIT = 1_000_000
MAX_LIST_WORDS = 500  # pro Länge; Rest nur Zähler


@dataclass
class ForbiddenLengthReport:
    n: int
    grammar_count: int
    arith_count: int
    forbidden_count: int
    forbidden_words: list[str] = field(default_factory=list)
    forbidden_truncated: bool = False
    enumeration_complete: bool = True
    hero: str | None = None  # erstes/kürzestes bekanntes Zeugenwort (BE bei n=2)


def compute_forbidden_at_length(
    n: int,
    arith: set[str],
    rules: GrammarRules | None = None,
    max_list: int = MAX_LIST_WORDS,
) -> ForbiddenLengthReport:
    rules = rules or GrammarRules()
    grammar_count = count_grammar_words(n, rules)
    grammar_set = set(generate_grammar_words(n, rules))
    enumeration_complete = len(grammar_set) == grammar_count
    forbidden = sorted(grammar_set - arith)
    truncated = len(forbidden) > max_list
    hero = "BE" if n == 2 and "BE" in forbidden else (forbidden[0] if forbidden else None)
    return ForbiddenLengthReport(
        n=n,
        grammar_count=grammar_count,
        arith_count=len(arith & grammar_set),
        forbidden_count=len(forbidden),
        forbidden_words=forbidden[:max_list],
        forbidden_truncated=truncated,
        enumeration_complete=enumeration_complete,
        hero=hero,
    )


def compute_ratios(
    lengths: list[int],
    buckets: list[set[str]],
    limit: int,
    max_enumerate: int = 3_000_000,
) -> list[dict]:
    """R(k)=|L_arith∩L|/|L| — gleiche Methode wie collatz_l_arith_test.analyze_length."""
    rows = []
    for k in lengths:
        rep = analyze_length(k, buckets[k], max_enumerate=max_enumerate)
        rep.limit = limit
        rows.append(
            {
                "k": k,
                "grammar_count": rep.grammar_count,
                "arith_in_L": rep.arith_in_L,
                "arith_total": rep.arith_total,
                "ratio": rep.ratio,
                "grammar_enumerated": rep.grammar_enumerated,
            }
        )
    return rows


def run_forbidden_catalog(
    lengths: list[int],
    limit: int,
    max_list: int = MAX_LIST_WORDS,
    ratio_lengths: list[int] | None = None,
) -> dict:
    rules = GrammarRules()
    max_n = max(lengths)
    bucket_k = max(max_n, max(ratio_lengths or (4, 8, 10)))
    buckets = realized_words_by_length(bucket_k, limit)
    reports = [compute_forbidden_at_length(n, buckets[n], rules, max_list) for n in lengths]

    ratio_lengths = sorted({k for k in (ratio_lengths or (4, 8, 10)) if k < len(buckets)})
    ratios = compute_ratios(ratio_lengths, buckets, limit) if ratio_lengths else []

    return {
        "grammar_rules": rules.as_dict(),
        "limit": limit,
        "lengths_requested": lengths,
        "F_n": [asdict(r) for r in reports],
        "ratios": ratios,
        "hero_word": "BE",
        "hero_length": 2,
        "honest_limits": [
            f"F_n = L(n) \\ L_arith(n) auf Suchtiefe ungerade n ≤ {limit}",
            "n ≡ 3,9 (mod 12): κ=none, tragen nicht zu L_arith bei",
            f"Vollliste L(n) nur für n ≤ {FULL_ENUM_MAX_K} in collatz_l_arith_test.py",
            "F_n ist Zeugen-/Experimentebene, kein Theorem über alle n",
        ],
        "pipeline": "Arithmetik → κ → L_arith → Verbotene Muster → Dynamische Konsequenzen",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="F_n = L(n) \\ L_arith(n)")
    parser.add_argument("--lengths", type=int, nargs="+", default=DEFAULT_LENGTHS)
    parser.add_argument("--limit", type=int, default=DEFAULT_LIMIT)
    parser.add_argument("--max-list", type=int, default=MAX_LIST_WORDS)
    parser.add_argument("--output", type=Path, default=Path("collatz_forbidden_words.json"))
    args = parser.parse_args()

    result = run_forbidden_catalog(args.lengths, args.limit, args.max_list)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")

    for row in result["F_n"]:
        n = row["n"]
        print(f"F_{n}: {row['forbidden_count']} Wörter (|L|={row['grammar_count']}, hero={row['hero']})")
        if row["forbidden_words"]:
            preview = row["forbidden_words"][:8]
            suffix = " …" if row["forbidden_truncated"] or len(row["forbidden_words"]) > 8 else ""
            print(f"  {preview}{suffix}")
    print(f"→ {args.output}")


if __name__ == "__main__":
    main()
