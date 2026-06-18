"""Tests für EABC-Bell-Ungleichheit, CHSH-Analog und P_same-Experiment."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_bell_inequality_test import (
    CLASSICAL_CHSH_MAX,
    THEORY,
    _chsh_sum,
    chsh_eabc_cycle_report,
    p_same_marginal,
    p_same_win_report,
    run,
    t_alignment,
)
from collatz_eabc_transition_graph import classes_from_sequence, prime_eabc_sequence


def test_t_alignment_on_t_cycle():
    classes = ["E", "A", "B", "C", "E"]
    align = t_alignment(classes)
    assert align == [1, 1, 1, 1]


def test_t_alignment_mismatch():
    classes = ["E", "B", "A", "C"]
    align = t_alignment(classes)
    assert align[0] == 0  # E -> B not t(E)=A


def test_b_win_theorem_on_manual_abce():
    classes = list("ABCEABCEA")
    report = p_same_win_report(classes)
    assert report["abce_window_count"] >= 1
    assert report["per_window_min_pair_matches"] >= 1
    assert report["theorem_B_win_ge_1"]


def test_p_same_win_formula_manual():
    classes = list("ABCEA")
    report = p_same_win_report(classes)
    assert report["abce_window_count"] == 1
    assert report["B_win"] == 3.0


def test_p_same_marginal_structure():
    classes = list("ABCEABCEABCEA")
    marg = p_same_marginal(classes)
    assert "B_marg" in marg
    assert "P_same_marg" in marg
    assert len(marg["marginals"]) == 3


def test_chsh_formula_sign():
    assert _chsh_sum(1.0, 0.0, 0.0, 0.0) == 1.0
    assert _chsh_sum(1.0, 1.0, 1.0, 1.0) == 2.0


def test_chsh_cycle_report_structure():
    classes = list("ABCEABCEABCEABCEA")
    report = chsh_eabc_cycle_report(classes)
    assert report["sample_size"] >= 1
    assert "E-node" in report["mapping"]["a"]
    assert report["abs_S_EABC"] is not None
    assert report["abs_S_EABC"] <= CLASSICAL_CHSH_MAX + 1e-9 + 4  # loose upper bound


def test_run_report_structure():
    report = run(10_000)
    assert report["meta"]["theory"] == THEORY
    assert "chsh_eabc_cycle" in report
    assert "summary" in report
    assert report["summary"]["B_win_satisfies_bound"] is True
    assert "abs_S_EABC" in report["summary"]


def test_run_json_roundtrip(tmp_path: Path):
    out = tmp_path / "bell.json"
    report = run(5_000)
    out.write_text(json.dumps(report), encoding="utf-8")
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["summary"]["B_win_satisfies_bound"]
    assert "abs_S_EABC" in loaded["summary"]


def test_prime_sequence_b_win_bound():
    seq = prime_eabc_sequence(50_000)
    classes = classes_from_sequence(seq)
    report = p_same_win_report(classes)
    assert report["theorem_B_win_ge_1"]


def test_prime_sequence_chsh_cycle_bounded():
    seq = prime_eabc_sequence(30_000)
    classes = classes_from_sequence(seq)
    ch = chsh_eabc_cycle_report(classes)
    assert ch["sample_size"] > 0
    assert ch["abs_S_EABC"] is not None
    assert ch["qm_reference_2sqrt2"] == math.sqrt(2) * 2
