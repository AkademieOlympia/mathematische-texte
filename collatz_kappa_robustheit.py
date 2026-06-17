#!/usr/bin/env python3
"""Stufe 2B Generalangriff: κ-Robustheit und Entropie von L_arith.

Kritische Schwachstelle: Realisierbarkeit hängt an der Wahl von κ.
Fall A: κ trifft Struktur → echte Grammatik-Lücke.
Fall B: κ erzeugt künstliche Ausschlüsse → R(k) misst nur κ-Artefakte.

Tests:
  1. Drei Kodierungen κ₁, κ₂, κ₃ — BE, F_n, R(k) vergleichen
  2. Bootstrap n ≤ 10^6, 10^7, (10^8 wenn machbar)
  3. Entropie h_F aus |F_n|-Folge (symbolische Dynamik)

Kein Collatz-Beweis.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Callable

from collatz_forbidden_words import compute_forbidden_at_length, compute_ratios
from collatz_kappa_test import iterate_u, nu2
from collatz_l_arith_test import (
    FULL_ENUM_MAX_K,
    GrammarRules,
    analyze_length,
    find_minimal_counterexamples,
    generate_grammar_words,
    realized_words_by_length,
)
from eabc_from_lean import EClass, class_of, t

LETTERS = ("E", "A", "B", "C")
RATIO_KS = (4, 8, 10)
FN_LENGTHS = list(range(2, 9))
FN_SIZES_REFERENCE = {2: 1, 3: 6, 4: 38, 5: 183, 6: 807, 7: 3402, 8: 13924}
DEFAULT_BOOTSTRAP = (1_000_000, 10_000_000)


class KappaId(str, Enum):
    KAPPA1 = "kappa1_naive_mod12"
    KAPPA2 = "kappa2_nu2_rotate"
    KAPPA3 = "kappa3_successor_block"


@dataclass(frozen=True, slots=True)
class KappaSpec:
    id: KappaId
    label: str
    description: str
    dynamics_note: str


KAPPA_SPECS: dict[KappaId, KappaSpec] = {
    KappaId.KAPPA1: KappaSpec(
        id=KappaId.KAPPA1,
        label="κ₁",
        description=(
            "Naive mod-12-Klasse entlang odd-to-odd-Bahn: "
            "Buchstabe_i = classOf(iterateU(n,i)); Abbruch bei classOf=none "
            "(n≡3,9 mod 12). Identisch mit collatz_l_arith_test.realized_words_by_length."
        ),
        dynamics_note="Dynamiktreu (Lean naiveKappa_shift); nicht injektiv (PR #38).",
    ),
    KappaId.KAPPA2: KappaSpec(
        id=KappaId.KAPPA2,
        label="κ₂",
        description=(
            "ν₂-Blocktyp: Buchstabe_i = t^{ν₂(3n_i+1) mod 4}(classOf(n_i)) mit n_i=iterateU(n,i), "
            "t die EABC-Rotation E→A→B→C→E. Nutzt die 2-adische Tiefe des geraden Blocks "
            "vor dem nächsten ungeraden Wert — gleiche mod-12-Basis, andere Feinstruktur."
        ),
        dynamics_note="Experimentell; testet ob Verbotsmuster von ν₂-Information stammen.",
    ),
    KappaId.KAPPA3: KappaSpec(
        id=KappaId.KAPPA3,
        label="κ₃",
        description=(
            "Successor-/Block-Kodierung: Buchstabe_i = classOf(iterateU(n,i+1)) für i<k-1, "
            "letzter Buchstabe = classOf(iterateU(n,k-1)). Kodiert die Klasse des *nächsten* "
            "odd-to-odd-Zustands (Nachfolger im EA-Block-Sinne), nicht den aktuellen."
        ),
        dynamics_note="Shift um einen odd-to-odd-Schritt; EA-Marker implizit über Nachfolger.",
    ),
}


def _rotate_class(cls: EClass, steps: int) -> EClass:
    cur = cls
    for _ in range(steps % 4):
        cur = t(cur)
    return cur


def kappa1_letter_at(n: int) -> str | None:
    cls = class_of(n)
    return cls.value if cls else None


def kappa2_letter_at(n: int) -> str | None:
    cls = class_of(n)
    if cls is None:
        return None
    v = nu2(3 * n + 1)
    return _rotate_class(cls, v).value


def kappa3_letter_at(n: int, next_n: int | None) -> str | None:
    target = next_n if next_n is not None else n
    cls = class_of(target)
    return cls.value if cls else None


LetterFn = Callable[[int, int], str | None]


def _letter_fn(variant: KappaId) -> LetterFn:
    if variant == KappaId.KAPPA1:

        def f(n: int, _i: int) -> str | None:
            return kappa1_letter_at(n)

        return f
    if variant == KappaId.KAPPA2:

        def f(n: int, _i: int) -> str | None:
            return kappa2_letter_at(n)

        return f

    def f(n: int, i: int) -> str | None:
        nxt = iterate_u(n, 1)
        return kappa3_letter_at(n, nxt)

    return f


def realized_words_variant(max_k: int, limit: int, variant: KappaId) -> list[set[str]]:
    """Für jedes k≤max_k die Menge der vollständigen κ-Präfixe Länge k."""
    letter_at = _letter_fn(variant)
    buckets = [set() for _ in range(max_k + 1)]
    for n in range(1, limit + 1, 2):
        letters: list[str] = []
        for i in range(max_k):
            cur = iterate_u(n, i)
            letter = letter_at(cur, i)
            if letter is None:
                break
            letters.append(letter)
            buckets[len(letters)].add("".join(letters))
    return buckets


def forbidden_sizes(buckets: list[set[str]], lengths: list[int]) -> dict[int, int]:
    rules = GrammarRules()
    return {
        n: compute_forbidden_at_length(n, buckets[n], rules, max_list=10_000).forbidden_count
        for n in lengths
        if n < len(buckets)
    }


def entropy_estimates(fn_sizes: dict[int, int]) -> dict:
    """h_F ≈ (1/n) log |F_n| pro n; limsup-Schätzer = max über n."""
    rows = []
    for n in sorted(fn_sizes):
        size = fn_sizes[n]
        if size <= 0:
            continue
        h_n = math.log(size) / n
        rows.append({"n": n, "F_n_size": size, "h_n": h_n})
    h_limsup = max((r["h_n"] for r in rows), default=None)
    return {
        "per_length": rows,
        "h_F_limsup_estimate": h_limsup,
        "reference_sizes_kappa1_n1e6": FN_SIZES_REFERENCE,
        "note": (
            "h_F = limsup_{n→∞} (1/n) log |F_n|; hier nur endliche n≤8, "
            "kein asymptotischer Beweis"
        ),
    }


@dataclass
class VariantReport:
    variant: str
    label: str
    description: str
    minimal_counterexamples: list[dict]
    be_forbidden: bool
    F_n_sizes: dict[int, int]
    ratios: list[dict]
    entropy: dict


def analyze_variant(
    variant: KappaId,
    limit: int,
    fn_lengths: list[int],
    ratio_ks: tuple[int, ...],
    max_enumerate: int,
) -> VariantReport:
    max_k = max(max(fn_lengths), max(ratio_ks))
    buckets = realized_words_variant(max_k, limit, variant)
    rules = GrammarRules()

    minimal = find_minimal_counterexamples(buckets, rules, max_k=FULL_ENUM_MAX_K)
    be_forbidden = "BE" not in buckets[2] if len(buckets) > 2 else True

    fn_sizes = forbidden_sizes(buckets, fn_lengths)
    ratios = compute_ratios(list(ratio_ks), buckets, limit, max_enumerate)

    return VariantReport(
        variant=variant.value,
        label=KAPPA_SPECS[variant].label,
        description=KAPPA_SPECS[variant].description,
        minimal_counterexamples=minimal[:5],
        be_forbidden=be_forbidden,
        F_n_sizes=fn_sizes,
        ratios=ratios,
        entropy=entropy_estimates(fn_sizes),
    )


@dataclass
class BootstrapRow:
    limit: int
    be_forbidden: bool
    minimal_word: str | None
    minimal_length: int | None
    F_n_sizes: dict[int, int]
    ratios: list[dict]
    new_minimal_vs_prev: list[str] = field(default_factory=list)


def bootstrap_kappa1(
    limits: list[int],
    fn_lengths: list[int],
    ratio_ks: tuple[int, ...],
    max_enumerate: int,
) -> list[BootstrapRow]:
    """Bootstrap nur für κ₁ (Referenz); neue minimale Wörter vs. vorherige Tiefe."""
    rows: list[BootstrapRow] = []
    prev_minimal: set[str] = set()
    for limit in sorted(limits):
        max_k = max(max(fn_lengths), max(ratio_ks))
        buckets = realized_words_by_length(max_k, limit)
        minimal_hits = find_minimal_counterexamples(buckets, max_k=FULL_ENUM_MAX_K)
        min_word = minimal_hits[0]["word"] if minimal_hits else None
        min_len = minimal_hits[0]["length"] if minimal_hits else None
        fn_sizes = forbidden_sizes(buckets, fn_lengths)
        ratios = compute_ratios(list(ratio_ks), buckets, limit, max_enumerate)
        new_min = [h["word"] for h in minimal_hits if h["word"] not in prev_minimal]
        prev_minimal = {h["word"] for h in minimal_hits}
        rows.append(
            BootstrapRow(
                limit=limit,
                be_forbidden="BE" not in buckets[2],
                minimal_word=min_word,
                minimal_length=min_len,
                F_n_sizes=fn_sizes,
                ratios=ratios,
                new_minimal_vs_prev=new_min,
            )
        )
    return rows


def compare_variants(
    variants: list[KappaId],
    limit: int,
    fn_lengths: list[int],
    ratio_ks: tuple[int, ...],
    max_enumerate: int,
) -> dict:
    reports = [
        analyze_variant(v, limit, fn_lengths, ratio_ks, max_enumerate) for v in variants
    ]
    be_stable = all(r.be_forbidden for r in reports)
    minimal_words = {r.variant: (r.minimal_counterexamples[0] if r.minimal_counterexamples else None) for r in reports}
    ratio_matrix = {}
    for k in ratio_ks:
        ratio_matrix[str(k)] = {r.variant: next((x["ratio"] for x in r.ratios if x["k"] == k), None) for r in reports}

    fn_matrix = {}
    for n in fn_lengths:
        fn_matrix[str(n)] = {r.variant: r.F_n_sizes.get(n) for r in reports}

    return {
        "limit": limit,
        "variants": [asdict(r) for r in reports],
        "comparison": {
            "be_forbidden_all_variants": be_stable,
            "minimal_counterexample": minimal_words,
            "R_k_by_variant": ratio_matrix,
            "F_n_sizes_by_variant": fn_matrix,
            "interpretation": (
                "Robustheit: BE und minimales Gegenbeispiel stabil über κ₁–κ₃ "
                "→ eher Fall A (Struktur). Abweichung → Fall B (κ-Artefakt)."
                if be_stable
                else "Mindestens eine κ-Variante realisiert BE — κ-abhängig, Fall B möglich."
            ),
        },
        "kappa_specs": {k.value: asdict(v) for k, v in KAPPA_SPECS.items()},
    }


def run_suite(
    limit: int,
    bootstrap_limits: list[int],
    max_enumerate: int = 3_000_000,
    include_1e8: bool = False,
) -> dict:
    variants = list(KappaId)
    boots = list(bootstrap_limits)
    if include_1e8 and 100_000_000 not in boots:
        boots.append(100_000_000)

    test1 = compare_variants(variants, limit, FN_LENGTHS, RATIO_KS, max_enumerate)
    test2_rows = bootstrap_kappa1(boots, FN_LENGTHS, RATIO_KS, max_enumerate)

    # Entropie aus κ₁-Referenz bei Haupt-limit
    ref_buckets = realized_words_by_length(max(FN_LENGTHS), limit)
    ref_fn = forbidden_sizes(ref_buckets, FN_LENGTHS)
    test3 = entropy_estimates(ref_fn)

    return {
        "stage": "2B",
        "title": "κ-Robustheit und Entropie von L_arith",
        "pr39_closed": (
            "PR #39 beantwortet: L_arith ⊊ L experimentell (BE, Ratios, F_n-Katalog); "
            "keine weitere Erweiterung — Stufe 2B testet κ-Abhängigkeit."
        ),
        "boxed_open_question": (
            "Welche arithmetischen Regeln erzeugen die verbotenen Wörter? "
            "BE ist Symptom, nicht Erklärung."
        ),
        "test1_kappa_comparison": test1,
        "test2_bootstrap_kappa1": {
            "variant": KappaId.KAPPA1.value,
            "limits": boots,
            "rows": [asdict(r) for r in test2_rows],
        },
        "test3_entropy": test3,
        "honest_limits": [
            "Drei κ-Varianten sind dokumentierte Experimente, keine treue FaithfulKappa",
            f"Bootstrap und Vergleich auf ungeraden n ≤ jeweilige Grenze",
            "n ≡ 3,9 (mod 12): classOf=none — kein Beitrag zu L_arith",
            f"F_n-Vollliste nur n ≤ {FULL_ENUM_MAX_K}; |F_8| aus Referenzlauf",
            "10^8-Bootstrap optional (--include-1e8); Laufzeit speicherintensiv",
            "h_F-Schätzung aus endlichen n — keine asymptotische Aussage",
            "Kein Collatz-Beweis; κ-Robustheit ≠ arithmetische Erklärung von BE",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Stufe 2B: κ-Robustheit")
    parser.add_argument("--limit", type=int, default=1_000_000, help="Hauptlimit für Test 1")
    parser.add_argument(
        "--bootstrap",
        type=int,
        nargs="+",
        default=list(DEFAULT_BOOTSTRAP),
        help="Bootstrap-Limits für κ₁",
    )
    parser.add_argument("--include-1e8", action="store_true", help="10^8 Bootstrap versuchen")
    parser.add_argument("--max-enumerate", type=int, default=3_000_000)
    parser.add_argument("--output", type=Path, default=Path("collatz_kappa_robustheit.json"))
    args = parser.parse_args()

    result = run_suite(
        args.limit,
        args.bootstrap,
        args.max_enumerate,
        args.include_1e8,
    )
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")

    cmp_ = result["test1_kappa_comparison"]["comparison"]
    print(f"BE forbidden (all κ): {cmp_['be_forbidden_all_variants']}")
    print("Minimal counterexamples:", cmp_["minimal_counterexample"])
    print("R(k):", json.dumps(cmp_["R_k_by_variant"], indent=2))
    for row in result["test2_bootstrap_kappa1"]["rows"]:
        print(
            f"Bootstrap n≤{row['limit']}: BE forbidden={row['be_forbidden']}, "
            f"min={row['minimal_word']}"
        )
    h = result["test3_entropy"]["h_F_limsup_estimate"]
    print(f"h_F limsup estimate (n≤8): {h:.4f}" if h else "h_F: n/a")
    print(f"→ {args.output}")


if __name__ == "__main__":
    main()
