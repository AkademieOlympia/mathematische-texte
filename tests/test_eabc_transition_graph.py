"""Tests für EABC-Übergangsgraph (collatz_eabc_transition_graph.py)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_transition_graph import (
    chi_from_bias,
    chi_transport,
    chi_transport_vs_quadruplet,
    count_t_cycle_windows,
    count_word_windows,
    edge_frequencies,
    markov_irreducible,
    prime_eabc_sequence,
    run,
    stationary_distribution,
    transition_counts,
    transition_matrix_report,
)
from eabc_from_lean import EClass, class_of, t


def test_t_rotation_definition():
    assert t(EClass.E) is EClass.A
    assert t(EClass.A) is EClass.B
    assert t(EClass.B) is EClass.C
    assert t(EClass.C) is EClass.E
    assert t(t(t(t(EClass.E)))) is EClass.E


def test_transition_counts_row_sum():
    seq = prime_eabc_sequence(500)
    classes = [r["class"] for r in seq]
    counts = transition_counts(classes)
    total_edges = sum(sum(row) for row in counts)
    assert total_edges == len(classes) - 1


def test_stationary_distribution_sums_to_one():
    seq = prime_eabc_sequence(2000)
    classes = [r["class"] for r in seq]
    counts = transition_counts(classes)
    stat = stationary_distribution(counts)
    pi = stat["pi_vector"]
    assert abs(sum(pi) - 1.0) < 1e-9
    assert all(x >= 0 for x in pi)


def test_word_window_abce_manual():
    classes = list("EABCEAB")
    assert count_word_windows(classes, "ABCE") == 1
    assert count_word_windows(classes, "CEAB") == 1


def test_t_cycle_window_manual():
    classes = ["E", "A", "B", "C", "E"]
    assert count_t_cycle_windows(classes, forward=True) == 2
    classes_rev = ["E", "C", "B", "A", "E"]
    assert count_t_cycle_windows(classes_rev, forward=False) == 2


def test_chi_from_bias_formula():
    assert chi_from_bias(3, 1) == 0.5
    assert chi_from_bias(0, 0) == 0.0


def test_edge_frequencies_t_bias_bounded():
    seq = prime_eabc_sequence(5000)
    classes = [r["class"] for r in seq]
    counts = transition_counts(classes)
    freqs = edge_frequencies(counts)
    assert -1.0 <= freqs["t_bias"] <= 1.0
    assert 0.0 <= freqs["t_forward_fraction"] <= 1.0
    assert 0.0 <= freqs["t_inverse_fraction"] <= 1.0
    assert freqs["total_transitions"] > 0


def test_chi_transport_bounded():
    seq = prime_eabc_sequence(5000)
    classes = [r["class"] for r in seq]
    trans = chi_transport(classes)
    assert -1.0 <= trans["chi_trans"] <= 1.0
    assert -1.0 <= trans["chi_t_cycle"] <= 1.0


def test_chi_transport_vs_quadruplet_honest():
    cmp = chi_transport_vs_quadruplet(5000, null_trials=100)
    assert "nicht identisch" in cmp["comparison"]["verdict"]
    assert "chi_E" in cmp
    assert "chi_trans" in cmp


def test_markov_irreducible_for_large_sample():
    report = transition_matrix_report(10_000)
    assert report["irreducible"]


def test_class_of_primes_gt3():
    for p in (5, 7, 11, 13, 17):
        assert class_of(p) is not None


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "trans.json"
    report = run(max_p=500, null_trials=50, output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["transition_matrix"]["prime_count"] > 0
    assert "cycle_holonomy" in loaded
    assert report["output_path"] == str(out)
