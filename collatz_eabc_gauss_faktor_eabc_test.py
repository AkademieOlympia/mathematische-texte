#!/usr/bin/env python3
"""
Z[i]-Faktor → EABC-Test: Gauß-Faktorpaar (a,b) mit p = a² + b² für p ≡ 1 (mod 4).

Kanonsiche Hypothese: collatz_eabc_normabstieg_hypothese.md §11

Für split-Primzahlen p > 3: kanonische Darstellung 0 < a ≤ b, a² + b² = p.
Projektion (a mod 12, b mod 12) auf EABC-Klassen via eabc_from_lean.class_of.

Epistemik:
  - Der bipartite Split↔E∪A-Test (collatz_eabc_gauss_defekt_test.py) ist arithmetisch
    trivial (mod 4 × mod 12).
  - Dieser Test ist die erste nicht-triviale Brücke: feine Faktorgeometrie in Z[i].

EABC-Sichtbarkeit:
  class_of(n) ist nur für n ≡ 1,5,7,11 (mod 12) definiert (E,A,B,C).
  Für p ≡ 1 (mod 4) ist genau eine Faktorleg gerade → höchstens eine Leg EABC-sichtbar
  pro Paar; beide gleichzeitig structurally impossible.

Ausführung:
    python3 collatz_eabc_gauss_faktor_eabc_test.py
    python3 collatz_eabc_gauss_faktor_eabc_test.py --max-p 10000
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
from collections import Counter
from dataclasses import dataclass
from math import isqrt
from pathlib import Path
from typing import Any

from eabc_from_lean import class_of

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_gauss_faktor_eabc.json"

EABC_RESIDUES = frozenset({1, 5, 7, 11})
EABC_LABELS = ("E", "A", "B", "C")
NULL_TRIALS = 5000


def _sieve_primes(limit: int) -> list[int]:
    if limit < 2:
        return []
    flags = [True] * (limit + 1)
    flags[0] = flags[1] = False
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            step = p
            start = p * p
            flags[start : limit + 1 : step] = [False] * (((limit - start) // step) + 1)
    return [i for i, ok in enumerate(flags) if ok]


def gaussian_factor_pair(p: int) -> tuple[int, int] | None:
    """
    Kanonische Darstellung p = a² + b² mit 0 < a ≤ b.
    Nur für p ≡ 1 (mod 4); sonst None.
    """
    if p <= 0 or p % 4 != 1:
        return None
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        b = isqrt(b2)
        if b > 0 and b * b == b2:
            return (a, b) if a <= b else (b, a)
    return None


def eabc_label(n: int) -> str | None:
    """EABC-Klasse von n mod 12, oder None außerhalb {1,5,7,11}."""
    ec = class_of(n)
    return ec.value if ec is not None else None


def residue_mod12(n: int) -> int:
    return n % 12


def odd_leg(pair: tuple[int, int]) -> tuple[str, int]:
    a, b = pair
    if a % 2 == 1:
        return "a", a
    return "b", b


def even_leg(pair: tuple[int, int]) -> tuple[str, int]:
    a, b = pair
    if a % 2 == 0:
        return "a", a
    return "b", b


@dataclass(frozen=True, slots=True)
class FactorRow:
    p: int
    a: int
    b: int
    kappa: str
    a_mod12: int
    b_mod12: int
    a_eabc: str | None
    b_eabc: str | None
    odd_leg: str
    even_mod12: int
    kappa_matches_a: bool
    kappa_matches_b: bool
    kappa_matches_odd: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "p": self.p,
            "a": self.a,
            "b": self.b,
            "kappa": self.kappa,
            "a_mod12": self.a_mod12,
            "b_mod12": self.b_mod12,
            "a_eabc": self.a_eabc,
            "b_eabc": self.b_eabc,
            "odd_leg": self.odd_leg,
            "even_mod12": self.even_mod12,
            "kappa_matches_a": self.kappa_matches_a,
            "kappa_matches_b": self.kappa_matches_b,
            "kappa_matches_odd": self.kappa_matches_odd,
        }


def classify_split_prime(p: int) -> FactorRow | None:
    if p <= 3 or p % 4 != 1:
        return None
    pair = gaussian_factor_pair(p)
    if pair is None:
        return None
    a, b = pair
    kappa = class_of(p)
    if kappa is None:
        return None
    a_eabc = eabc_label(a)
    b_eabc = eabc_label(b)
    odd_name, odd_val = odd_leg(pair)
    _, even_val = even_leg(pair)
    return FactorRow(
        p=p,
        a=a,
        b=b,
        kappa=kappa.value,
        a_mod12=residue_mod12(a),
        b_mod12=residue_mod12(b),
        a_eabc=a_eabc,
        b_eabc=b_eabc,
        odd_leg=odd_name,
        even_mod12=residue_mod12(even_val),
        kappa_matches_a=a_eabc == kappa.value,
        kappa_matches_b=b_eabc == kappa.value,
        kappa_matches_odd=eabc_label(odd_val) == kappa.value,
    )


def _counter_to_nested_dict(c: Counter[tuple[Any, ...]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, count in sorted(c.items(), key=lambda kv: (-kv[1], kv[0])):
        if len(key) == 1:
            out[str(key[0])] = count
        else:
            out[str(key)] = count
    return out


def _chi_square_uniform(observed: dict[str, int], categories: tuple[str, ...]) -> float:
    total = sum(observed.get(c, 0) for c in categories)
    if total == 0:
        return 0.0
    expected = total / len(categories)
    return sum((observed.get(c, 0) - expected) ** 2 / expected for c in categories)


def _null_match_rate(
    pairs: list[tuple[str, str]],
    trials: int = NULL_TRIALS,
    seed: int = 0,
) -> dict[str, float]:
    """Shuffle-Null: κ bleibt, ungerade Leg-Klasse permutiert."""
    if not pairs:
        return {"mean": 0.0, "std": 0.0, "trials": 0}
    rng = random.Random(seed)
    rates: list[float] = []
    legs = [leg for _, leg in pairs]
    for _ in range(trials):
        shuffled = legs[:]
        rng.shuffle(shuffled)
        matches = sum(1 for (k, _), leg in zip(pairs, shuffled) if k == leg)
        rates.append(matches / len(pairs))
    return {
        "mean": statistics.mean(rates),
        "std": statistics.stdev(rates) if len(rates) > 1 else 0.0,
        "trials": trials,
    }


def _independent_uniform_null(
    n: int,
    trials: int = NULL_TRIALS,
    seed: int = 0,
) -> dict[str, float]:
    """Unabhängig uniform auf {E,A,B,C}² — Referenz für 2-leg-Paare (hypothetisch)."""
    rng = random.Random(seed + 1)
    rates: list[float] = []
    for _ in range(trials):
        matches = sum(
            1
            for _ in range(n)
            if rng.choice(EABC_LABELS) == rng.choice(EABC_LABELS)
        )
        rates.append(matches / n if n else 0.0)
    return {
        "mean": statistics.mean(rates),
        "std": statistics.stdev(rates) if len(rates) > 1 else 0.0,
        "trials": trials,
    }


def factor_distribution_report(max_p: int) -> dict[str, Any]:
    rows: list[FactorRow] = []
    raw_mod12_pairs: Counter[tuple[int, int]] = Counter()
    eabc_pair_joint: Counter[tuple[str, str]] = Counter()
    kappa_x_odd_leg: Counter[tuple[str, str]] = Counter()
    kappa_x_leg_site: Counter[tuple[str, str, str]] = Counter()
    even_mod12_dist: Counter[int] = Counter()
    odd_on_smaller_leg = 0
    both_eabc_visible = 0
    exactly_one_eabc = 0
    kappa_match_a = kappa_match_b = kappa_match_odd = 0

    for p in _sieve_primes(max_p):
        row = classify_split_prime(p)
        if row is None:
            continue
        rows.append(row)
        raw_mod12_pairs[(row.a_mod12, row.b_mod12)] += 1
        even_mod12_dist[row.even_mod12] += 1
        if row.a_eabc and row.b_eabc:
            both_eabc_visible += 1
            eabc_pair_joint[(row.a_eabc, row.b_eabc)] += 1
        elif row.a_eabc or row.b_eabc:
            exactly_one_eabc += 1
        if row.odd_leg == "a":
            odd_on_smaller_leg += 1
        if row.kappa_matches_a:
            kappa_match_a += 1
        if row.kappa_matches_b:
            kappa_match_b += 1
        if row.kappa_matches_odd:
            kappa_match_odd += 1
        odd_val = row.a if row.odd_leg == "a" else row.b
        odd_class = eabc_label(odd_val)
        if odd_class:
            kappa_x_odd_leg[(row.kappa, odd_class)] += 1
        for site, label in (("a", row.a_eabc), ("b", row.b_eabc)):
            if label:
                kappa_x_leg_site[(row.kappa, label, site)] += 1

    n = len(rows)
    odd_pairs = [(r.kappa, eabc_label(r.a if r.odd_leg == "a" else r.b)) for r in rows]
    odd_pairs = [(k, c) for k, c in odd_pairs if c is not None]
    observed_odd_match = (
        sum(1 for k, c in odd_pairs if k == c) / len(odd_pairs) if odd_pairs else 0.0
    )
    null_odd = _null_match_rate(odd_pairs)

    # Marginal EABC on visible legs
    leg_marginal: Counter[str] = Counter()
    for r in rows:
        if r.a_eabc:
            leg_marginal[r.a_eabc] += 1
        if r.b_eabc:
            leg_marginal[r.b_eabc] += 1

    even_chi2 = _chi_square_uniform(
        {str(k): v for k, v in even_mod12_dist.items()},
        tuple(str(x) for x in sorted(even_mod12_dist.keys())),
    )

    z_score_odd = 0.0
    if null_odd["std"] > 0:
        z_score_odd = (observed_odd_match - null_odd["mean"]) / null_odd["std"]

    if both_eabc_visible > 0:
        joint_verdict = "unexpected_both_eabc_visible: investigate representation convention"
    elif abs(z_score_odd) > 2.0:
        joint_verdict = (
            f"kappa_odd_leg_asymmetry: |z|={abs(z_score_odd):.2f} vs shuffle null"
        )
    elif even_chi2 > 15.0:
        joint_verdict = (
            "even_leg_mod12_nonuniform: structured even-component residues "
            "(not EABC-native; document as geometric prior)"
        )
    else:
        joint_verdict = (
            "no_significant_eabc_chirality: κ vs odd leg matches shuffle null; "
            "joint EABC pair empty by p≡1(mod4) parity (one leg always even); "
            "bipartite split test remains trivial"
        )

    return {
        "meta": {
            "max_p": max_p,
            "hypothesis_doc": "collatz_eabc_normabstieg_hypothese.md",
            "eabc_residues": sorted(EABC_RESIDUES),
            "eabc_visibility_note": (
                "class_of only for mod12 in {1,5,7,11}; split primes have one even leg "
                "→ at most one EABC-visible factor per pair"
            ),
            "epistemics": {
                "bipartite_split_ea_test": "trivial (mod 4 × mod 12)",
                "this_test": "non-trivial Z[i] factor geometry → EABC fine resolution",
                "both_legs_eabc": "structurally impossible for p≡1(mod4)",
            },
        },
        "counts": {
            "split_primes_p_gt_3": n,
            "both_eabc_visible": both_eabc_visible,
            "exactly_one_eabc_visible": exactly_one_eabc,
            "odd_on_smaller_leg_a": odd_on_smaller_leg,
            "kappa_matches_a": kappa_match_a,
            "kappa_matches_b": kappa_match_b,
            "kappa_matches_odd_eabc_leg": kappa_match_odd,
        },
        "distributions": {
            "raw_mod12_pairs_top20": _counter_to_nested_dict(
                Counter(dict(raw_mod12_pairs.most_common(20)))
            ),
            "eabc_joint_pairs": dict(eabc_pair_joint),
            "kappa_x_odd_leg_eabc": _counter_to_nested_dict(kappa_x_odd_leg),
            "kappa_x_leg_site": _counter_to_nested_dict(kappa_x_leg_site),
            "even_leg_mod12": dict(sorted(even_mod12_dist.items())),
            "visible_leg_marginal_eabc": dict(leg_marginal),
        },
        "asymmetry": {
            "kappa_odd_leg_match_rate": observed_odd_match,
            "null_shuffle_odd_leg": null_odd,
            "null_independent_uniform": _independent_uniform_null(len(odd_pairs)),
            "z_score_kappa_vs_odd_leg": z_score_odd,
            "even_mod12_chi2_vs_uniform": even_chi2,
            "verdict": joint_verdict,
        },
        "sample_rows": [r.to_dict() for r in rows[:30]],
    }


def format_distribution_table(report: dict[str, Any]) -> str:
    d = report["distributions"]
    a = report["asymmetry"]
    c = report["counts"]
    lines = [
        "Z[i]-Faktor → EABC (split-Primzahlen)",
        "=" * 55,
        f"split p>3: {c['split_primes_p_gt_3']}",
        f"beide EABC-sichtbar: {c['both_eabc_visible']} (erwartet: 0)",
        f"genau eine EABC-Leg: {c['exactly_one_eabc_visible']}",
        "",
        "κ × ungerade Leg (EABC):",
    ]
    for key, val in d["kappa_x_odd_leg_eabc"].items():
        lines.append(f"  {key}: {val}")
    lines.extend(
        [
            "",
            "Gerade Leg mod 12:",
            "  " + ", ".join(f"{k}:{v}" for k, v in d["even_leg_mod12"].items()),
            "",
            f"κ = ungerade Leg: {a['kappa_odd_leg_match_rate']:.4f}  "
            f"(Null μ={a['null_shuffle_odd_leg']['mean']:.4f}, "
            f"z={a['z_score_kappa_vs_odd_leg']:.2f})",
            f"verdict: {a['verdict']}",
        ]
    )
    return "\n".join(lines)


def run(max_p: int = 10_000, output: Path | None = None) -> dict[str, Any]:
    report = factor_distribution_report(max_p)
    out = output or DEFAULT_OUTPUT
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report["output_path"] = str(out)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Z[i]-Faktor → EABC-Verteilungstest")
    parser.add_argument("--max-p", type=int, default=10_000, help="obere Schranke für p")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="JSON-Ausgabedatei",
    )
    args = parser.parse_args()
    report = run(max_p=args.max_p, output=args.output)
    print(format_distribution_table(report))
    print(f"\nJSON: {report['output_path']}")


if __name__ == "__main__":
    main()
