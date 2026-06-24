#!/usr/bin/env python3
"""
Gauß–EABC-Defekt-Test: Korrelation split/inert in Z[i] mit EABC-Klassen.

Kanonsiche Hypothese: collatz_eabc_normabstieg_hypothese.md §8

Für Primzahlen p > 3:
  - p ≡ 1 (mod 4)  split in Z[i]  ↔  κ(p) ∈ {E, A}  (4n+1)
  - p ≡ 3 (mod 4)  inert in Z[i]  ↔  κ(p) ∈ {B, C}  (4n+3)

Split/inert ist Theorem; die interpretative Gauß–EABC-Brücke ist Conjecture/Heuristik.
Die grobe bipartite Zuordnung ist für p > 3 arithmetisch exakt (mod-12-Restklassen).

Ausführung:
    python3 collatz_eabc_gauss_defekt_test.py
    python3 collatz_eabc_gauss_defekt_test.py --max-p 10000 --output collatz_eabc_gauss_defekt.json
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from math import isqrt
from pathlib import Path
from typing import Any

from eabc_from_lean import EClass, class_of

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_gauss_defekt.json"

EA_CLASSES = frozenset({EClass.E, EClass.A})
BC_CLASSES = frozenset({EClass.B, EClass.C})


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


def gauss_split_class(p: int) -> str | None:
    """
    Zerlegungsverhalten rationaler Primzahl p in Z[i].
    'split' (p≡1 mod 4), 'inert' (p≡3 mod 4), 'ramified' (p=2), None (nicht prim).
    """
    if p == 2:
        return "ramified"
    if p < 2:
        return None
    r = p % 4
    if r == 1:
        return "split"
    if r == 3:
        return "inert"
    return None


def eabc_coarse_bucket(eclass: EClass) -> str:
    """Grobe EABC-Seite: EA (4n+1-Rest) vs BC (4n+3-Rest)."""
    return "EA" if eclass in EA_CLASSES else "BC"


def gauss_coarse_bucket(split_class: str) -> str:
    return "EA" if split_class == "split" else "BC"


@dataclass(frozen=True, slots=True)
class PrimeRow:
    p: int
    mod4: int
    mod12: int
    gauss: str
    eabc: str
    coarse_match: bool
    fine_note: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "p": self.p,
            "mod4": self.mod4,
            "mod12": self.mod12,
            "gauss": self.gauss,
            "eabc": self.eabc,
            "coarse_match": self.coarse_match,
            "fine_note": self.fine_note,
        }


def classify_prime(p: int) -> PrimeRow | None:
    gauss = gauss_split_class(p)
    if gauss is None or gauss == "ramified":
        return None
    eabc = class_of(p)
    if eabc is None:
        return None
    coarse_match = gauss_coarse_bucket(gauss) == eabc_coarse_bucket(eabc)
    fine_note = (
        "split→E/A or inert→B/C (coarse only; E vs A not determined by split/inert)"
    )
    return PrimeRow(
        p=p,
        mod4=p % 4,
        mod12=p % 12,
        gauss=gauss,
        eabc=eabc.value,
        coarse_match=coarse_match,
        fine_note=fine_note,
    )


def correlation_report(max_p: int) -> dict[str, Any]:
    primes = _sieve_primes(max_p)
    rows: list[PrimeRow] = []
    counts = {
        "split_E": 0,
        "split_A": 0,
        "inert_B": 0,
        "inert_C": 0,
        "coarse_match": 0,
        "coarse_mismatch": 0,
        "total_p_gt_3": 0,
    }
    cross: dict[str, dict[str, int]] = {
        "split": {"E": 0, "A": 0, "B": 0, "C": 0},
        "inert": {"E": 0, "A": 0, "B": 0, "C": 0},
    }

    for p in primes:
        row = classify_prime(p)
        if row is None:
            continue
        rows.append(row)
        counts["total_p_gt_3"] += 1
        cross[row.gauss][row.eabc] += 1
        key = f"{row.gauss}_{row.eabc}"
        if key in counts:
            counts[key] += 1
        if row.coarse_match:
            counts["coarse_match"] += 1
        else:
            counts["coarse_mismatch"] += 1

    total = counts["total_p_gt_3"]
    coarse_rate = counts["coarse_match"] / total if total else 0.0

    if coarse_rate == 1.0 and total > 0:
        mapping_verdict = (
            "exact_coarse_bipartition: split↔E∪A and inert↔B∪C for all p>3 "
            "(arithmetically forced by mod 4 and mod 12); "
            "interpretive EABC-bridge remains Conjecture/Heuristik"
        )
    elif coarse_rate >= 0.99:
        mapping_verdict = "approximate_coarse: minor exceptions — investigate"
    else:
        mapping_verdict = "coarse_mapping_fails: unexpected mismatches"

    return {
        "meta": {
            "max_p": max_p,
            "hypothesis_doc": "collatz_eabc_normabstieg_hypothese.md",
            "epistemics": {
                "gauss_split_inert": "Theorem",
                "eabc_gross_mapping": "Heuristik/Conjecture (coarse partition exact for p>3)",
                "eabc_fine_E_vs_A": "not determined by split/inert alone",
            },
        },
        "counts": counts,
        "cross_table": cross,
        "coarse_match_rate": coarse_rate,
        "mapping_verdict": mapping_verdict,
        "sample_primes": [r.to_dict() for r in rows[:40]],
        "mismatches": [r.to_dict() for r in rows if not r.coarse_match],
    }


def format_correlation_table(report: dict[str, Any]) -> str:
    cross = report["cross_table"]
    lines = [
        "Gauss–EABC Korrelation (grobe Zuordnung)",
        "=" * 50,
        f"{'Gauss':<8} {'E':>6} {'A':>6} {'B':>6} {'C':>6}",
    ]
    for gauss in ("split", "inert"):
        row = cross[gauss]
        lines.append(
            f"{gauss:<8} {row['E']:>6} {row['A']:>6} {row['B']:>6} {row['C']:>6}"
        )
    lines.append("")
    c = report["counts"]
    lines.append(
        f"p>3 counted: {c['total_p_gt_3']}  "
        f"coarse match: {c['coarse_match']}  "
        f"rate: {report['coarse_match_rate']:.6f}"
    )
    lines.append(f"verdict: {report['mapping_verdict']}")
    return "\n".join(lines)


def run(max_p: int = 5000, output: Path | None = None) -> dict[str, Any]:
    report = correlation_report(max_p)
    out = output or DEFAULT_OUTPUT
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report["output_path"] = str(out)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Gauß–EABC-Defekt-Korrelationstest")
    parser.add_argument("--max-p", type=int, default=5000, help="obere Schranke für Primzahlen")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="JSON-Ausgabedatei",
    )
    args = parser.parse_args()
    report = run(max_p=args.max_p, output=args.output)
    print(format_correlation_table(report))
    print(f"\nJSON: {report['output_path']}")


if __name__ == "__main__":
    main()
