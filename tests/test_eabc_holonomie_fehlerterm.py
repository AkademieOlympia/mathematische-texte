"""Tests für EABC-Holonomie Fehlerterm D_E und Lückenmuster."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_holonomie_fehlerterm import (
    CANONICAL_GAP_PATTERN,
    THEORY_ENDFORM,
    chebyshev_bias_comparison,
    d_tilde_e_series,
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
    n_plus = row["N_plus"]
    n_minus = row["N_minus"]
    total = n_plus + n_minus
    assert row["N_ABCEA"] == n_plus
    assert row["N_CEABC"] == n_minus
    assert row["D_E"] == n_plus - n_minus
    if total > 0:
        assert abs(row["chi_Hol"] - row["D_E"] / total) < 1e-12
        assert abs(row["D_tilde_E"] - row["D_E"] / (total**0.5)) < 1e-12


def test_holonomy_counts_matches_chi_hol_sliding():
    max_p = 20_000
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    hol = chi_hol_sliding(classes)
    row = holonomy_counts(max_p)
    assert row["N_plus"] == hol["abcea_windows"]
    assert row["N_minus"] == hol["ceabc_windows"]
    assert row["chi_Hol"] == hol["chi_hol"]


def test_d_tilde_e_series():
    limits = [1000, 10_000, 50_000]
    series = d_tilde_e_series(limits)
    assert len(series) == 3
    for pt, lim in zip(series, limits):
        assert pt["X"] == lim
        row = holonomy_counts(lim)
        assert pt["D_tilde_E"] == row["D_tilde_E"]
        assert pt["N_plus"] == row["N_plus"]


def test_chebyshev_comparison_structure():
    cmp = chebyshev_bias_comparison(30_000)
    assert "holonomy_series" in cmp
    assert "D_tilde_E_series" in cmp
    assert "chebyshev_mod4_series" in cmp
    assert "qualitative" in cmp
    assert len(cmp["holonomy_series"]) >= 3
    assert len(cmp["D_tilde_E_series"]) == len(cmp["holonomy_series"])
    assert cmp["qualitative"]["supports_oscillating_D_E"] is True


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "fehlerterm.json"
    report = run(max_p=10_000, output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["meta"]["theory_endform"] == THEORY_ENDFORM
    assert loaded["gap_patterns"]["words"][ABCEA_WORD]["matches_canonical"]
    assert "D_tilde_E_series" in loaded
    assert "boxed_conclusions" in loaded
    assert "de_bell_combined" in loaded
    if loaded["de_bell_combined"]:
        db = loaded["de_bell_combined"]
        assert db["D_E"] == db["N_plus"] - db["N_minus"]
        assert db["S_EABC"] is not None
    assert report["output_path"] == str(out)
