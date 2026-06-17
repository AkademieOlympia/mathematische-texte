#!/usr/bin/env python3
"""Stufe 2 Generalangriff: Grammatik L(k) und arithmetische Realisierbarkeit L_arith(k).

Regeln (vgl. collatz_equivalenz_e_infty.tex, collatz_schlussartikel_arxiv.tex §C-Ketten,
CollatzEabc.Density):
  1. BB ∉ L (B→B verboten)
  2. Endliche C-Ketten: beim Eintritt in eine C-Kette Kapazität cap ≥ 1 (ν₂-Budget)
  3. EA-Zwang nach *maximaler* C-Kette (cap C-Buchstaben → E, A)

Kein Collatz-Beweis.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass, field
from functools import lru_cache
from pathlib import Path

from collatz_kappa_test import iterate_u
from eabc_from_lean import class_of

LETTERS = ("E", "A", "B", "C")
LAST_NONE = 4
LAST = {"E": 0, "A": 1, "B": 2, "C": 3, "": LAST_NONE}

# mode: 0=free, 1=in_c, 2=need_e, 3=need_a

# k ≤ FULL_ENUM_MAX_K: Vollliste L(k); darüber Stichprobe (Grammatik explodiert).
FULL_ENUM_MAX_K = 8
DEFAULT_LENGTHS = [10, 20, 30]
SAMPLE_SIZE = 5000


@dataclass(frozen=True, slots=True)
class GrammarRules:
    forbid_bb: bool = True
    finite_c_chains: bool = True
    ea_after_max_c: bool = True
    max_c_cap: int | None = None

    def as_dict(self) -> dict:
        return {
            "forbid_bb": self.forbid_bb,
            "finite_c_chains": self.finite_c_chains,
            "ea_after_max_c": self.ea_after_max_c,
            "max_c_cap": self.max_c_cap,
            "references": [
                "collatz_equivalenz_e_infty.tex",
                "collatz_schlussartikel_arxiv.tex §C-Ketten",
                "CollatzEabc.Density",
            ],
        }


def _cap_upper(remaining: int, rules: GrammarRules) -> int:
    hi = remaining
    if rules.max_c_cap is not None:
        hi = min(hi, rules.max_c_cap)
    return max(0, hi)


@lru_cache(maxsize=None)
def _count_dp(k: int, pos: int, mode: int, c_run: int, c_cap: int, last: int, rules_key: tuple) -> int:
    rules = GrammarRules(*rules_key)
    if pos == k:
        return 1 if mode == 0 else 0
    remaining = k - pos

    if mode == 2:
        return _count_dp(k, pos + 1, 3, 0, 0, LAST["E"], rules_key)
    if mode == 3:
        return _count_dp(k, pos + 1, 0, 0, 0, LAST["A"], rules_key)

    if mode == 1:
        total = 0
        if c_run < c_cap:
            total += _count_dp(k, pos + 1, 1, c_run + 1, c_cap, LAST["C"], rules_key)
            for letter in ("E", "A", "B"):
                total += _count_dp(k, pos + 1, 0, 0, 0, LAST[letter], rules_key)
        elif remaining >= 2:
            total += _count_dp(k, pos + 2, 0, 0, 0, LAST["A"], rules_key)
        return total

    # free
    total = 0
    for letter in LETTERS:
        if rules.forbid_bb and letter == "B" and last == LAST["B"]:
            continue
        if letter == "C" and rules.finite_c_chains:
            for cap in range(1, _cap_upper(remaining, rules) + 1):
                if rules.ea_after_max_c and cap + 2 > remaining:
                    continue
                total += _count_dp(k, pos + 1, 1, 1, cap, LAST["C"], rules_key)
        else:
            total += _count_dp(k, pos + 1, 0, 0, 0, LAST[letter], rules_key)
    return total


def count_grammar_words(k: int, rules: GrammarRules | None = None) -> int:
    rules = rules or GrammarRules()
    key = (rules.forbid_bb, rules.finite_c_chains, rules.ea_after_max_c, rules.max_c_cap)
    _count_dp.cache_clear()
    return _count_dp(k, 0, 0, 0, 0, LAST_NONE, key)


def _collect_words(
    k: int,
    pos: int,
    mode: int,
    c_run: int,
    c_cap: int,
    last: int,
    prefix: list[str],
    rules: GrammarRules,
    out: list[str],
    max_words: int | None,
) -> None:
    if max_words is not None and len(out) >= max_words:
        return
    if pos == k:
        if mode == 0:
            out.append("".join(prefix))
        return
    remaining = k - pos

    if mode == 2:
        _collect_words(k, pos + 1, 3, 0, 0, LAST["E"], prefix + ["E"], rules, out, max_words)
        return
    if mode == 3:
        _collect_words(k, pos + 1, 0, 0, 0, LAST["A"], prefix + ["A"], rules, out, max_words)
        return

    if mode == 1:
        if c_run < c_cap:
            _collect_words(
                k, pos + 1, 1, c_run + 1, c_cap, LAST["C"], prefix + ["C"], rules, out, max_words
            )
            for letter in ("E", "A", "B"):
                _collect_words(
                    k, pos + 1, 0, 0, 0, LAST[letter], prefix + [letter], rules, out, max_words
                )
        elif remaining >= 2:
            _collect_words(k, pos + 2, 0, 0, 0, LAST["A"], prefix + ["E", "A"], rules, out, max_words)
        return

    for letter in LETTERS:
        if rules.forbid_bb and letter == "B" and last == LAST["B"]:
            continue
        if letter == "C" and rules.finite_c_chains:
            for cap in range(1, _cap_upper(remaining, rules) + 1):
                if rules.ea_after_max_c and cap + 2 > remaining:
                    continue
                _collect_words(
                    k, pos + 1, 1, 1, cap, LAST["C"], prefix + ["C"], rules, out, max_words
                )
        else:
            _collect_words(
                k, pos + 1, 0, 0, 0, LAST[letter], prefix + [letter], rules, out, max_words
            )


def generate_grammar_words(
    k: int,
    rules: GrammarRules | None = None,
    max_words: int | None = None,
) -> list[str]:
    rules = rules or GrammarRules()
    out: list[str] = []
    _collect_words(k, 0, 0, 0, 0, LAST_NONE, [], rules, out, max_words)
    return out


def realized_words_by_length(max_k: int, limit: int) -> list[set[str]]:
    """Für jedes k≤max_k die Menge der vollständigen κ-Präfixe Länge k (ein Lauf pro n)."""
    buckets = [set() for _ in range(max_k + 1)]
    for n in range(1, limit + 1, 2):
        letters: list[str] = []
        for i in range(max_k):
            cls = class_of(iterate_u(n, i))
            if cls is None:
                break
            letters.append(cls.value)
            buckets[len(letters)].add("".join(letters))
    return buckets


def is_grammar_valid_word(w: str, rules: GrammarRules | None = None) -> bool:
    rules = rules or GrammarRules()
    k = len(w)

    def bt(pos: int, mode: int, c_run: int, c_cap: int) -> bool:
        if pos == k:
            return mode == 0
        if mode == 2:
            return w[pos] == "E" and bt(pos + 1, 3, 0, 0)
        if mode == 3:
            return w[pos] == "A" and bt(pos + 1, 0, 0, 0)
        if mode == 1:
            if w[pos] == "C" and c_run < c_cap:
                return bt(pos + 1, 1, c_run + 1, c_cap)
            if c_run < c_cap and w[pos] in "EAB":
                return bt(pos + 1, 0, 0, 0)
            if c_run == c_cap and pos + 1 < k and w[pos : pos + 2] == "EA":
                return bt(pos + 2, 0, 0, 0)
            return False
        # free
        if w[pos] == "B" and pos > 0 and w[pos - 1] == "B":
            return False
        if w[pos] == "C" and rules.finite_c_chains:
            remaining = k - pos
            for cap in range(1, _cap_upper(remaining, rules) + 1):
                if rules.ea_after_max_c and cap + 2 > remaining:
                    continue
                if bt(pos + 1, 1, 1, cap):
                    return True
            return False
        return bt(pos + 1, 0, 0, 0)

    return bt(0, 0, 0, 0)


@dataclass
class LengthReport:
    k: int
    grammar_count: int
    arith_in_L: int
    arith_total: int
    ratio: float | None
    limit: int
    grammar_enumerated: bool
    minimal_non_realizable: list[str] = field(default_factory=list)
    sample_non_realizable: list[str] = field(default_factory=list)
    notes: str = ""


def _should_full_enumerate(k: int, grammar_count: int, max_enumerate: int) -> bool:
    return k <= FULL_ENUM_MAX_K and grammar_count <= max_enumerate


def analyze_length(
    k: int,
    arith: set[str],
    rules: GrammarRules | None = None,
    max_enumerate: int = 500_000,
) -> LengthReport:
    rules = rules or GrammarRules()
    grammar_count = count_grammar_words(k, rules)

    if _should_full_enumerate(k, grammar_count, max_enumerate):
        grammar_set = set(generate_grammar_words(k, rules))
        non_real = sorted(grammar_set - arith)
        arith_in_L = len(arith & grammar_set)
        ratio = arith_in_L / grammar_count if grammar_count else 0.0
        return LengthReport(
            k=k,
            grammar_count=grammar_count,
            arith_in_L=arith_in_L,
            arith_total=len(arith),
            ratio=ratio,
            limit=0,
            grammar_enumerated=True,
            minimal_non_realizable=non_real[:20],
            sample_non_realizable=non_real[:10],
        )

    sample_words = generate_grammar_words(k, rules, max_words=SAMPLE_SIZE)
    non_real = [w for w in sample_words if w not in arith]
    arith_in_sample = sum(1 for w in sample_words if w in arith)
    est_ratio = arith_in_sample / len(sample_words) if sample_words else None
    return LengthReport(
        k=k,
        grammar_count=grammar_count,
        arith_in_L=arith_in_sample,
        arith_total=len(arith),
        ratio=est_ratio,
        limit=0,
        grammar_enumerated=False,
        sample_non_realizable=non_real[:10],
        notes=(
            f"|L({k})|={grammar_count} — Stichprobe n={len(sample_words)} "
            f"(Vollliste ab k>{FULL_ENUM_MAX_K} oder |L|>{max_enumerate})"
        ),
    )


def find_minimal_counterexamples(
    buckets: list[set[str]],
    rules: GrammarRules | None = None,
    max_k: int | None = None,
    max_per_length: int = 5,
) -> list[dict]:
    """Kürzeste Länge mit w ∈ L(k) \\ L_arith(k); optional weitere Beispiele gleicher Länge."""
    rules = rules or GrammarRules()
    upper = min(len(buckets) - 1, max_k or len(buckets) - 1)
    for k in range(1, upper + 1):
        hits: list[str] = []
        for w in generate_grammar_words(k, rules):
            if w not in buckets[k]:
                hits.append(w)
                if len(hits) >= max_per_length:
                    break
        if hits:
            return [{"word": w, "length": k} for w in hits]
    return []


def run_suite(
    lengths: list[int],
    limit: int,
    max_enumerate: int = 500_000,
) -> dict:
    rules = GrammarRules()
    max_k = max(lengths)
    buckets = realized_words_by_length(max_k, limit)
    reports = []
    for k in lengths:
        rep = analyze_length(k, buckets[k], rules, max_enumerate)
        rep.limit = limit
        reports.append(rep)

    minimal = find_minimal_counterexamples(buckets, rules, max_k=FULL_ENUM_MAX_K)

    return {
        "grammar_rules": rules.as_dict(),
        "full_enum_max_k": FULL_ENUM_MAX_K,
        "limit": limit,
        "lengths": [asdict(r) for r in reports],
        "minimal_counterexamples": minimal,
        "pr38_context": {
            "naive_kappa_injective": False,
            "naive_kappa_dynamics": True,
            "note": (
                "PR #38: naive κ dynamiktreu, nicht injektiv — κ_naiv scheidet als finale Brücke aus; "
                "treue κ (FaithfulKappa) bleibt offen und soll verborgene Verbotsregeln von L_arith sichtbar machen"
            ),
        },
        "honest_limits": [
            f"Realisierbarkeit nur für ungerade n ≤ {limit} mit vollständigem κ-Präfix (kein none)",
            "n ≡ 3,9 (mod 12) liefern κ=none und tragen nicht zu L_arith bei",
            f"k ≤ {FULL_ENUM_MAX_K}: Vollliste L(k); k > {FULL_ENUM_MAX_K}: Stichprobe (|L(k)| explodiert)",
            "Ratio bei Stichprobe ist Schätzung, kein exakter |L_arith|/|L|-Wert",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="L_arith vs L")
    parser.add_argument("--lengths", type=int, nargs="+", default=DEFAULT_LENGTHS)
    parser.add_argument("--limit", type=int, default=1_000_000)
    parser.add_argument("--max-enumerate", type=int, default=500_000)
    parser.add_argument("--output", type=Path, default=Path("collatz_l_arith_test.json"))
    args = parser.parse_args()
    result = run_suite(args.lengths, args.limit, args.max_enumerate)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result["lengths"], indent=2))
    if result["minimal_counterexamples"]:
        print("Minimale Gegenbeispiele:", result["minimal_counterexamples"])
    print(f"→ {args.output}")


if __name__ == "__main__":
    main()
