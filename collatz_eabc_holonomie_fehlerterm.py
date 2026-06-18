#!/usr/bin/env python3
"""
EABC-Holonomie: Fehlerterm D_E / Delta_E, Sagnac-Observable S_E, normalisiertes D̃_E.

Theorie: collatz_eabc_zirkulationshypothese.md (kanonisch), collatz_eabc_fehlerterm_hypothese.md (Teilhypothese)

  N_plus(X), N_minus(X)  — Zählung geschlossener 5-Zyklen (ABCEA / CEABC)
  Delta_E(X) = D_E(X) = N_plus - N_minus   (EABC-Sagnac-Fehlerterm)
  S_E(X) = (N_plus - N_minus) / (N_plus + N_minus)  (analog ΔT/T)
  D̃_E(X) = Delta_E / sqrt(N_plus + N_minus)

Ausführung:
    python3 collatz_eabc_holonomie_fehlerterm.py
    python3 collatz_eabc_holonomie_fehlerterm.py --max-p 1000000
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from collatz_eabc_transition_graph import (
    ABCEA_WORD,
    CEABC_WORD,
    chi_hol_sliding,
    classes_from_sequence,
    prime_eabc_sequence,
)
from eabc_from_lean import EClass, residue

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_holonomie_fehlerterm.json"
THEORY_ZIRKULATION = "collatz_eabc_zirkulationshypothese.md"
THEORY_SAGNAC = "collatz_eabc_sagnac.md"
THEORY_ENDFORM = "collatz_eabc_fehlerterm_hypothese.md"
THEORY_BWEISVERSUCH = "collatz_eabc_holonomie_beweisversuch.md"

CANONICAL_GAP_PATTERN = (2, 4, 2, 4)
EABC_RESIDUES = {EClass.E: 1, EClass.A: 5, EClass.B: 7, EClass.C: 11}
WORD_RESIDUES = {
    ABCEA_WORD: tuple(EABC_RESIDUES[EClass(c)] for c in ABCEA_WORD),
    CEABC_WORD: tuple(EABC_RESIDUES[EClass(c)] for c in CEABC_WORD),
}


def gap_pattern_mod12(residues: tuple[int, ...]) -> tuple[int, ...]:
    """Lücken (r_{i+1}-r_i) mod 12 entlang einer geschlossenen Restklassenfolge."""
    if len(residues) < 2:
        return ()
    gaps: list[int] = []
    for i in range(len(residues) - 1):
        gaps.append((residues[i + 1] - residues[i]) % 12)
    return tuple(gaps)


def verify_gap_patterns() -> dict[str, Any]:
    """Lemma-Skizze §3: ABCEA und CEABC tragen (2,4,2,4)."""
    rows: dict[str, Any] = {}
    for word in (ABCEA_WORD, CEABC_WORD):
        res = WORD_RESIDUES[word]
        gaps = gap_pattern_mod12(res)
        rows[word] = {
            "residues_mod12": list(res),
            "gap_pattern": list(gaps),
            "matches_canonical": gaps == CANONICAL_GAP_PATTERN,
            "start_class_mod12": res[0],
        }
    return {
        "canonical_gap_pattern": list(CANONICAL_GAP_PATTERN),
        "words": rows,
        "only_start_differs": (
            rows[ABCEA_WORD]["start_class_mod12"] == 5
            and rows[CEABC_WORD]["start_class_mod12"] == 11
        ),
        "verdict": (
            "Beide Wörter tragen Lückenmuster (2,4,2,4); "
            f"Startklassen 5 (ABCEA) vs. 11 (CEABC)."
        ),
    }


def holonomy_counts(max_p: int) -> dict[str, Any]:
    """N_plus, N_minus, S_E, Delta_E, D_E, D̃_E für Primfolge bis max_p."""
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    hol = chi_hol_sliding(classes)
    n_plus = hol["abcea_windows"]
    n_minus = hol["ceabc_windows"]
    total = n_plus + n_minus
    delta_e = n_plus - n_minus
    s_e = delta_e / total if total > 0 else 0.0
    d_tilde = delta_e / math.sqrt(total) if total > 0 else 0.0
    return {
        "max_p": max_p,
        "X": max_p,
        "prime_count": len(seq),
        "N_plus": n_plus,
        "N_minus": n_minus,
        "N_ABCEA": n_plus,
        "N_CEABC": n_minus,
        "Delta_E": delta_e,
        "D_E": delta_e,
        "C_E": delta_e,
        "S_E": s_e,
        "chi_Hol": hol["chi_hol"],
        "D_tilde_E": d_tilde,
        "nonzero_windows": hol["nonzero_windows"],
        "omega_sum": hol["omega_sum"],
    }


def sagnac_report(max_p: int) -> dict[str, Any]:
    """Sagnac-Observable am Prim-Obergrenze X=max_p (collatz_eabc_sagnac.md)."""
    row = holonomy_counts(max_p)
    total = row["N_plus"] + row["N_minus"]
    return {
        "X": max_p,
        "gamma_plus": ABCEA_WORD,
        "gamma_minus": CEABC_WORD,
        "N_plus": row["N_plus"],
        "N_minus": row["N_minus"],
        "Delta_E": row["Delta_E"],
        "D_E": row["D_E"],
        "S_E": row["S_E"],
        "D_tilde_E": row["D_tilde_E"],
        "N_total": total,
        "rotationsfrei": abs(row["S_E"]) < 1e-15 if total > 0 else True,
        "theory": THEORY_SAGNAC,
        "epistemic": {
            "Delta_E": "Definition (EABC-Sagnac-Fehlerterm)",
            "S_E": "Definition (normalisierte Sagnac-Observable)",
            "S_E_to_zero": "Vermutung (Hauptterm Hol_E=0)",
        },
    }


def d_tilde_e_series(limits: list[int]) -> list[dict[str, Any]]:
    """Zeitreihe von D̃_E(X) über die angegebenen Prim-Obergrenzen."""
    return [
        {
            "X": row["X"],
            "N_plus": row["N_plus"],
            "N_minus": row["N_minus"],
            "Delta_E": row["Delta_E"],
            "D_E": row["D_E"],
            "S_E": row["S_E"],
            "D_tilde_E": row["D_tilde_E"],
            "chi_Hol": row["chi_Hol"],
        }
        for row in (holonomy_counts(lim) for lim in limits)
    ]


def chebyshev_bias_comparison(max_p: int) -> dict[str, Any]:
    """
    Qualitativer Vergleich: D_E-Oszillation vs. klassischer Chebyshev-Bias mod 4.

    Chebyshev: # {p ≤ X : p ≡ 3 (4)} - # {p ≤ X : p ≡ 1 (4)} oszilliert vorzeichenbehaftet.
    Analog: D_E(X) oszilliert bei endlichen X, während χ_Hol(X) klein bleiben kann.
    """
    limits = sorted({x for x in (1000, 5000, 10_000, 50_000, 100_000, 500_000, max_p) if x <= max_p})
    hol_series: list[dict[str, Any]] = []
    cheb_series: list[dict[str, Any]] = []

    for lim in limits:
        row = holonomy_counts(lim)
        hol_series.append(
            {
                "X": lim,
                "N_plus": row["N_plus"],
                "N_minus": row["N_minus"],
                "D_E": row["D_E"],
                "D_tilde_E": row["D_tilde_E"],
                "chi_Hol": row["chi_Hol"],
                "N_total": row["N_plus"] + row["N_minus"],
            }
        )

        primes = [r["p"] for r in prime_eabc_sequence(lim)]
        mod1 = sum(1 for p in primes if p % 4 == 1)
        mod3 = sum(1 for p in primes if p % 4 == 3)
        cheb_series.append(
            {
                "X": lim,
                "chebyshev_D": mod3 - mod1,
                "chebyshev_D_tilde": (mod3 - mod1) / math.sqrt(mod1 + mod3) if (mod1 + mod3) else 0.0,
                "count_mod1": mod1,
                "count_mod3": mod3,
            }
        )

    d_signs = [1 if r["D_E"] > 0 else (-1 if r["D_E"] < 0 else 0) for r in hol_series]
    sign_changes = sum(
        1 for i in range(1, len(d_signs)) if d_signs[i] != 0 and d_signs[i - 1] != 0 and d_signs[i] != d_signs[i - 1]
    )

    last = hol_series[-1] if hol_series else {}
    chi_values = [r["chi_Hol"] for r in hol_series if r["N_total"] > 0]
    chi_range = max(chi_values) - min(chi_values) if chi_values else 0.0
    chi_bounded = all(abs(c) < 0.25 for c in chi_values)

    return {
        "holonomy_series": hol_series,
        "D_tilde_E_series": [{"X": r["X"], "D_tilde_E": r["D_tilde_E"]} for r in hol_series],
        "chebyshev_mod4_series": cheb_series,
        "qualitative": {
            "D_E_positive_at_max": last.get("D_E", 0) > 0,
            "chi_Hol_at_max": last.get("chi_Hol", 0.0),
            "chi_Hol_range_over_limits": chi_range,
            "D_E_sign_changes": sign_changes,
            "chebyshev_D_at_max": cheb_series[-1]["chebyshev_D"] if cheb_series else 0,
            "analogy": (
                "D_E bleibt vorzeichenbehaftet und wächst mit X (wie Chebyshev-Differenz), "
                "während χ_Hol = D_E/(N_plus+N_minus) durch wachsenden Nenner gedämpft wird — "
                "konsistent mit Hol_E=0 als Hauptterm und strukturiertem Fehlerterm."
            ),
            "supports_Hol_E_to_zero": (
                chi_bounded and abs(last.get("chi_Hol", 1.0)) < 0.25
                if last
                else False
            ),
            "supports_oscillating_D_E": last.get("D_E", 0) != 0 if last else False,
        },
    }


def run_series(
    limits: list[int] | None = None,
    max_p: int = 1_000_000,
) -> dict[str, Any]:
    if limits is None:
        limits = sorted({x for x in (1_000, 10_000, 100_000, 500_000, max_p) if x <= max_p})
    series = [holonomy_counts(lim) for lim in limits]
    gap = verify_gap_patterns()
    cheb = chebyshev_bias_comparison(max_p)
    de_bell = None
    circulation = None
    try:
        from collatz_eabc_bell_inequality_test import de_bell_combined_report

        de_bell = de_bell_combined_report(max_p)
    except ImportError:
        pass
    try:
        from collatz_eabc_sagnac_circulation import circulation_report

        circulation = circulation_report(max_p)
    except ImportError:
        pass
    return {
        "meta": {
            "module": "collatz_eabc_holonomie_fehlerterm.py",
            "theory_zirkulation": THEORY_ZIRKULATION,
            "theory_sagnac": THEORY_SAGNAC,
            "theory_endform": THEORY_ENDFORM,
            "theory": THEORY_BWEISVERSUCH,
            "theory_bell_bridge": "collatz_eabc_bell_holonomie.md §12 (sekundäre Analogie)",
            "max_p": max_p,
            "limits": limits,
        },
        "gap_patterns": gap,
        "sagnac_report": sagnac_report(max_p),
        "holonomy_series": series,
        "D_tilde_E_series": d_tilde_e_series(limits),
        "chebyshev_comparison": cheb,
        "de_bell_combined": de_bell,
        "circulation_report": circulation,
        "boxed_conclusions": {
            "zirkulation_hypothesis": "EABC als Zirkulationstheorie (collatz_eabc_zirkulationshypothese.md)",
            "Hol_E_main_term": "S_E -> 0 unter mod-12-Symmetrie + HL-Äquidistribution (Hauptvermutung)",
            "interesting_question": "D_E(X) als Prime Race; Chebyshev-Analogie / L-Funktionen mod 12",
        },
        "epistemic_labels": {
            "gap_pattern": "Lemma-Skizze",
            "Hol_E_zero": "Vermutung",
            "Delta_E": "Definition",
            "S_E": "Definition",
            "D_E": "Definition (Alias Delta_E)",
            "D_tilde_E": "Definition",
            "fehlerterm_hypothese": "Hypothese",
            "numerics": "Experiment",
            "chebyshev_analogy": "Heuristik",
            "bell_chsh": "Analogie (sekundär)",
        },
    }


def run(
    max_p: int = 1_000_000,
    output: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    report = run_series(max_p=max_p)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC-Holonomie Fehlerterm D_E")
    parser.add_argument("--max-p", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(max_p=args.max_p, output=args.output)
    gap = report["gap_patterns"]
    print("=== EABC Holonomie Fehlerterm ===")
    print(gap["verdict"])
    for word, info in gap["words"].items():
        print(f"  {word}: gaps={info['gap_pattern']}, start={info['start_class_mod12']}")
    print()
    for row in report["holonomy_series"]:
        print(
            f"X={row['max_p']:>7}: N_plus={row['N_plus']:4}, N_minus={row['N_minus']:4}, "
            f"S_E={row['S_E']:+.4f}, Delta_E={row['Delta_E']:+4}, D̃_E={row['D_tilde_E']:+.3f}"
        )
    print()
    sag = report["sagnac_report"]
    print(
        f"Sagnac @ X={sag['X']}: gamma+={sag['gamma_plus']}  gamma-={sag['gamma_minus']}  "
        f"S_E={sag['S_E']:+.4f}  Delta_E={sag['Delta_E']:+d}"
    )
    print()
    print("D̃_E Zeitreihe:")
    for pt in report["D_tilde_E_series"]:
        print(f"  X={pt['X']:>7}: D̃_E={pt['D_tilde_E']:+.3f}")
    qual = report["chebyshev_comparison"]["qualitative"]
    print()
    print(qual["analogy"])
    if report.get("de_bell_combined"):
        db = report["de_bell_combined"]
        print(
            f"Bell bridge @ X={db['X']}: D_E={db['D_E']:+d}  "
            f"D̃_E={db['D_tilde_E']:+.3f}  S_EABC={db['S_EABC']:.4f}  "
            f"|S|={db['abs_S_EABC']:.4f}"
        )
    print(f"supports Hol_E→0 (χ bounded <0.25): {qual['supports_Hol_E_to_zero']}")
    print(f"supports oscillating D_E≠0: {qual['supports_oscillating_D_E']}")
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
