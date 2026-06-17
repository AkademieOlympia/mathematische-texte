#!/usr/bin/env python3
"""
EABC-Produktbäume: Klammerung vs. Assoziativität (H-Testbed, O-Theorie).

Kanonsiche Theorie: collatz_eabc_plattenuebergang.md §2.5
Querverweise: collatz_eabc_oktonion_singularitaet.md §3.6,
              collatz_eabc_zerlegungsregimen.md §5

Für n = f_1 · … · f_k (k ≥ 2, f_i ≥ 2) gibt es C_{k-1} binäre Produktbäume (Catalan).
  - H (assoziativ): alle Klammerungen liefern dieselbe Multiplikationsabbildung → 1 effektive Klasse/Faktorisierung
  - O (nicht assoziativ): Klammerung kann geometrisch verschiedene Kanäle liefern → bis zu C_{k-1} Klassen

Ausführung:
    python3 collatz_eabc_product_tree_stub.py
    python3 collatz_eabc_product_tree_stub.py --max-n 50
"""

from __future__ import annotations

import argparse
import json
from math import isqrt
from pathlib import Path
from typing import Any

from collatz_eabc_Z_decomposition_test import Z_EABC_count, Z_fact, build_cache

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_product_tree_stub.json"


def catalan(n: int) -> int:
    """C_n = binomial(2n,n)/(n+1); hier n = k-1 für k Faktoren."""
    if n < 0:
        return 0
    if n <= 1:
        return 1
    c = 1
    for i in range(1, n + 1):
        c = c * (n + i) // i
    return c // (n + 1)


def _factor_tuples(n: int, min_factor: int = 2) -> list[tuple[int, ...]]:
    """Alle geordneten Faktortupel (f_1,...,f_k), k≥2, f_i≥min_factor, Produkt n."""
    if n < min_factor * min_factor:
        return []
    out: list[tuple[int, ...]] = []

    def rec(remaining: int, start: int, acc: list[int]) -> None:
        if remaining == 1:
            if len(acc) >= 2:
                out.append(tuple(acc))
            return
        for f in range(start, remaining + 1):
            if remaining % f != 0:
                continue
            acc.append(f)
            rec(remaining // f, f, acc)
            acc.pop()

    rec(n, min_factor, [])
    return out


def unordered_compositions(n: int, min_factor: int = 2) -> list[tuple[int, ...]]:
    """Ungeordnete Multifaktor-Zerlegungen n = f_1·…·f_k (sortiertes Tupel)."""
    seen: set[tuple[int, ...]] = set()
    for tup in _factor_tuples(n, min_factor=min_factor):
        key = tuple(sorted(tup))
        seen.add(key)
    return sorted(seen)


def tree_counts_for_n(n: int) -> dict[str, Any]:
    comps = unordered_compositions(n)
    if not comps:
        return {
            "n": n,
            "Z_fact": 0,
            "composition_count": 0,
            "Z_tree_raw": 0,
            "Z_tree_eff_H": 0,
            "max_catalan_per_composition": 0,
            "compositions": [],
        }
    per_comp: list[dict[str, Any]] = []
    z_raw = 0
    for comp in comps:
        k = len(comp)
        c_trees = catalan(k - 1)
        z_raw += c_trees
        per_comp.append(
            {
                "factors": list(comp),
                "k": k,
                "catalan_trees": c_trees,
                "effective_H": 1,
            }
        )
    return {
        "n": n,
        "Z_fact": Z_fact(n),
        "composition_count": len(comps),
        "Z_tree_raw": z_raw,
        "Z_tree_eff_H": len(comps),
        "max_catalan_per_composition": max(p["catalan_trees"] for p in per_comp),
        "compositions": per_comp,
    }


def octonion_theoretical_note(max_n: int) -> dict[str, Any]:
    """Keine Schalen-Enumeration in O — nur Catalan-Summen über Kompositionen."""
    rows = [tree_counts_for_n(n) for n in range(2, max_n + 1)]
    composites = [r for r in rows if r["composition_count"] > 0]
    multi_k3 = [r for r in composites if any(c["k"] >= 3 for c in r["compositions"])]
    return {
        "algebra": "O",
        "enumeration": "theoretical_only",
        "note": (
            "Assoziativität fehlt: pro Komposition n=f_1…f_k bis zu C_{k-1} geometrisch "
            "verschiedene Produktbäume denkbar; volle Z_tree^O erfordert μ_n auf Σ_n^(8)."
        ),
        "n_with_k_ge_3_compositions": len(multi_k3),
        "example_k_ge_3": multi_k3[:5],
        "max_Z_tree_raw_up_to_n": max((r["Z_tree_raw"] for r in composites), default=0),
    }


def hurwitz_quaternion_report(max_n: int) -> dict[str, Any]:
    cache = build_cache(max_n)
    rows: list[dict[str, Any]] = []
    for n in range(2, max_n + 1):
        tc = tree_counts_for_n(n)
        ze = Z_EABC_count(n, cache) if n in cache else None
        tc["Z_EABC"] = ze
        tc["tree_adds_beyond_Z_fact"] = tc["Z_tree_eff_H"] > tc["Z_fact"]
        tc["Z_EABC_eq_Z_fact"] = ze == tc["Z_fact"] if ze is not None else None
        rows.append(tc)

    composites = [r for r in rows if r["composition_count"] > 0]
    z_fact_eq_ze = sum(
        1 for r in composites if r["Z_EABC"] is not None and r["Z_EABC"] == r["Z_fact"]
    )
    eff_gt_fact = sum(1 for r in composites if r["tree_adds_beyond_Z_fact"])

    return {
        "algebra": "H",
        "max_n": max_n,
        "composite_count": len(composites),
        "Z_EABC_equals_Z_fact_fraction": round(
            z_fact_eq_ze / len(composites), 6
        )
        if composites
        else 0.0,
        "composition_count_gt_Z_fact_count": eff_gt_fact,
        "verdict": (
            "Auf H kollabiert jede Klammerung zu 1 effektiver Klasse pro Komposition; "
            "Catalan-Bäume erhöhen Z_EABC nicht (weiterhin Z_EABC≈Z_fact). "
            "Mehrfach-Kompositionen (k≥3) zählen mehr Zerlegungsmodi als binäres Z_fact, "
            "ändern aber nicht die bekannte EABC-Kanal-Zählung über Paare (σ_a,σ_b)."
        ),
        "rows": rows,
    }


def product_tree_report(max_n: int = 50) -> dict[str, Any]:
    h = hurwitz_quaternion_report(max_n)
    o = octonion_theoretical_note(max_n)
    return {
        "meta": {
            "hypothesis_doc": "collatz_eabc_plattenuebergang.md §2.5",
            "max_n": max_n,
            "catalan_formula": "C_{k-1} binary trees for k factors",
            "H_associativity": "Z_tree_eff = 1 per composition",
            "O_non_associativity": "up to C_{k-1} classes per composition (theoretical)",
        },
        "hurwitz_quaternion": h,
        "octonion_theoretical": o,
        "contrast": (
            "H: Norm-Multiplikativität und Assoziativität — Produktabbildung Σ_a×Σ_b→Σ_n "
            "unabhängig von Klammerung; binäre Catalan-Zahl ist stets 1. "
            "O: Σ_a×Σ_b→Σ_ab für Norm OK, aber (xy)z≠x(yz); Klammerung kann geometrisch "
            "neue EABC-Rekonstruktionsmodi liefern — das ist die oktanion-spezifische Novelty."
        ),
    }


def format_summary(report: dict[str, Any]) -> str:
    h = report["hurwitz_quaternion"]
    lines = [
        "EABC-Produktbäume (H vs. O)",
        "=" * 44,
        f"n ≤ {report['meta']['max_n']}",
        f"H: Z_EABC=Z_fact bei {100 * h['Z_EABC_equals_Z_fact_fraction']:.1f}% der zusammengesetzten n",
        f"H: Kompositionen > Z_fact bei {h['composition_count_gt_Z_fact_count']} n "
        "(k≥3-Faktoren; Klammerung trotzdem trivial)",
        f"O: {report['octonion_theoretical']['n_with_k_ge_3_compositions']} n mit k≥3-Komposition "
        f"(theoretisch bis max Z_tree_raw={report['octonion_theoretical']['max_Z_tree_raw_up_to_n']})",
        "",
        report["contrast"],
    ]
    return "\n".join(lines)


def run(max_n: int = 50, output: Path | None = None) -> dict[str, Any]:
    report = product_tree_report(max_n=max_n)
    out = output or DEFAULT_OUTPUT
    out.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    report["output_path"] = str(out)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC-Produktbäume H/O")
    parser.add_argument("--max-n", type=int, default=50)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(max_n=args.max_n, output=args.output)
    print(format_summary(report))
    print(f"\nJSON: {report['output_path']}")


if __name__ == "__main__":
    main()
