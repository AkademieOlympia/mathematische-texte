"""Tests für EABC-Übergangsgraph und Zyklus-Holonomie (collatz_eabc_transition_graph.py)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_transition_graph import (
    chi_E_sliding,
    chi_from_bias,
    chi_sliding_vs_quadruplet,
    chi_transport,
    chi_transport_vs_quadruplet,
    count_t_cycle_windows,
    count_word_windows,
    edge_frequencies,
    isotropy_null_chi_E,
    markov_irreducible,
    omega_window,
    prime_eabc_sequence,
    run,
    shuffle_null_chi_E,
    sliding_windows,
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


def test_omega_window_values():
    assert omega_window("ABCE") == 1
    assert omega_window("CEAB") == -1
    assert omega_window("EABC") == 0
    assert omega_window("ABCA") == 0


def test_chi_E_sliding_formula_manual():
    classes = list("EABCEAB")
    report = chi_E_sliding(classes)
    assert report["abce_windows"] == 1
    assert report["ceab_windows"] == 1
    assert report["omega_sum"] == 0
    assert report["chi_E"] == 0.0


def test_chi_E_sliding_bias_manual():
    classes = list("ABCEEABC")
    report = chi_E_sliding(classes)
    assert report["abce_windows"] == 1
    assert report["ceab_windows"] == 0
    assert report["chi_E"] == 1.0


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


def test_sliding_windows_count():
    classes = ["E", "A", "B", "C", "E"]
    assert len(sliding_windows(classes)) == 2


def test_edge_frequencies_t_bias_bounded():
    seq = prime_eabc_sequence(5000)
    classes = [r["class"] for r in seq]
    counts = transition_counts(classes)
    freqs = edge_frequencies(counts)
    assert -1.0 <= freqs["t_bias"] <= 1.0
    assert 0.0 <= freqs["t_forward_fraction"] <= 1.0
    assert 0.0 <= freqs["t_inverse_fraction"] <= 1.0
    assert freqs["total_transitions"] > 0


def test_chi_transport_equals_chi_E_sliding():
    seq = prime_eabc_sequence(5000)
    classes = [r["class"] for r in seq]
    trans = chi_transport(classes)
    sliding = chi_E_sliding(classes)
    assert trans["chi_E"] == sliding["chi_E"]
    assert trans["chi_trans"] == sliding["chi_E"]


def test_chi_sliding_vs_quadruplet_honest():
    cmp = chi_sliding_vs_quadruplet(5000, null_trials=100)
    assert "nicht identisch" in cmp["comparison"]["verdict"]
    assert "chi_E" in cmp
    assert "chi_E_quad" in cmp
    assert "null_models" in cmp


def test_shuffle_and_isotropy_null_bounded():
    seq = prime_eabc_sequence(3000)
    classes = [r["class"] for r in seq]
    shuffle = shuffle_null_chi_E(classes, trials=50)
    isotropy = isotropy_null_chi_E(classes, trials=50)
    assert shuffle["null_type"] == "marginal_shuffle"
    assert isotropy["null_type"] == "isotropy_relabel"
    assert -1.0 <= shuffle["observed_chi_E"] <= 1.0
    assert -1.0 <= isotropy["observed_chi_E"] <= 1.0


def test_chi_transport_vs_quadruplet_alias():
    a = chi_sliding_vs_quadruplet(2000, null_trials=20)
    b = chi_transport_vs_quadruplet(2000, null_trials=20)
    assert a["chi_E"] == b["chi_E"]


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
    assert loaded["meta"]["theory"] == "collatz_eabc_zyklus_holonomie.md"
    assert "hol_E_support" in loaded["hol_E_estimates"]
    assert report["output_path"] == str(out)
