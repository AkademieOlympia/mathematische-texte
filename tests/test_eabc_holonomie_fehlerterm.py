"""Tests für EABC-Holonomie Fehlerterm D_E und Lückenmuster."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_holonomie_fehlerterm import (
    CANONICAL_GAP_PATTERN,
    chebyshev_bias_comparison,
    gap_pattern_mod12,
    holonomy_counts,
    run,
    verify_gap_patterns,
)
from collatz_eabc_transition_graph import ABCEA_WORD, CEABC_WORD, chi_hol_sliding
from collatz_eabc_transition_graph import classes_from_sequence, prime_eabc_sequence


def test_gap_pattern_abcea_ceabc():
    gap = verify_gap_patterns()
    assert gap["words"][ABCEA_WORD]["matches_canonical"]
    assert gap["words"][CEABC_WORD]["matches_canonical"]
    assert gap["words"][ABCEA_WORD]["gap_pattern"] == [2, 4, 2, 4]
    assert gap["words"][CEABC_WORD]["gap_pattern"] == [2, 4, 2, 4]
    assert gap["only_start_differs"]


def test_gap_pattern_closed_loop():
    residues = (5, 7, 11, 1, 5)
    assert gap_pattern_mod12(residues) == CANONICAL_GAP_PATTERN


def test_holonomy_counts_formula():
    row = holonomy_counts(50_000)
    n_ab = row["N_ABCEA"]
    n_ce = row["N_CEABC"]
    total = n_ab + n_ce
    assert row["D_E"] == n_ab - n_ce
    if total > 0:
        assert abs(row["chi_Hol"] - row["D_E"] / total) < 1e-12
        assert abs(row["D_tilde_E"] - row["D_E"] / (total**0.5)) < 1e-12


def test_holonomy_counts_matches_chi_hol_sliding():
    max_p = 20_000
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    hol = chi_hol_sliding(classes)
    row = holonomy_counts(max_p)
    assert row["N_ABCEA"] == hol["abcea_windows"]
    assert row["N_CEABC"] == hol["ceabc_windows"]
    assert row["chi_Hol"] == hol["chi_hol"]


def test_chebyshev_comparison_structure():
    cmp = chebyshev_bias_comparison(30_000)
    assert "holonomy_series" in cmp
    assert "chebyshev_mod4_series" in cmp
    assert "qualitative" in cmp
    assert len(cmp["holonomy_series"]) >= 3
    assert cmp["qualitative"]["supports_oscillating_D_E"] is True


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "fehlerterm.json"
    report = run(max_p=10_000, output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["meta"]["theory"] == "collatz_eabc_holonomie_beweisversuch.md"
    assert loaded["gap_patterns"]["words"][ABCEA_WORD]["matches_canonical"]
    assert "boxed_conclusions" in loaded
    assert report["output_path"] == str(out)
