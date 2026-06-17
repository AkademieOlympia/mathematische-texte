"""Tests für EABC-Übergangsgraph, Pfadorientierung (4-Block) und Holonomie (5-Block)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_transition_graph import (
    chi_E_sliding,
    chi_from_bias,
    chi_hol_sliding,
    chi_path_sliding,
    chi_path_vs_hol,
    chi_sliding_vs_quadruplet,
    chi_transport,
    chi_transport_vs_quadruplet,
    count_t_cycle_windows,
    count_word_windows,
    edge_frequencies,
    isotropy_null_chi_E,
    isotropy_null_chi_hol,
    isotropy_null_chi_path,
    markov_irreducible,
    omega_5,
    omega_path,
    omega_window,
    prime_eabc_sequence,
    run,
    shuffle_null_chi_E,
    shuffle_null_chi_hol,
    shuffle_null_chi_path,
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


def test_omega_path_values():
    assert omega_path("ABCE") == 1
    assert omega_path("CEAB") == -1
    assert omega_path("EABC") == 0
    assert omega_path("ABCA") == 0
    assert omega_window("ABCE") == 1  # legacy alias


def test_omega_5_values():
    assert omega_5("ABCEA") == 1
    assert omega_5("CEABC") == -1
    assert omega_5("ABCE") == 0
    assert omega_5("EABCE") == 0


def test_chi_path_sliding_formula_manual():
    classes = list("EABCEAB")
    report = chi_path_sliding(classes)
    assert report["abce_windows"] == 1
    assert report["ceab_windows"] == 1
    assert report["omega_sum"] == 0
    assert report["chi_path"] == 0.0
    assert report["chi_E"] == 0.0


def test_chi_path_sliding_bias_manual():
    classes = list("ABCEEABC")
    report = chi_path_sliding(classes)
    assert report["abce_windows"] == 1
    assert report["ceab_windows"] == 0
    assert report["chi_path"] == 1.0


def test_chi_hol_sliding_manual():
    classes = ["A", "B", "C", "E", "A"]
    report = chi_hol_sliding(classes)
    assert report["abcea_windows"] == 1
    assert report["ceabc_windows"] == 0
    assert report["chi_hol"] == 1.0


def test_chi_E_sliding_legacy_alias():
    classes = list("ABCEEABC")
    assert chi_E_sliding(classes)["chi_E"] == chi_path_sliding(classes)["chi_path"]


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
    assert len(sliding_windows(classes, width=4)) == 2
    assert len(sliding_windows(classes, width=5)) == 1


def test_edge_frequencies_t_bias_bounded():
    seq = prime_eabc_sequence(5000)
    classes = [r["class"] for r in seq]
    counts = transition_counts(classes)
    freqs = edge_frequencies(counts)
    assert -1.0 <= freqs["t_bias"] <= 1.0
    assert 0.0 <= freqs["t_forward_fraction"] <= 1.0
    assert 0.0 <= freqs["t_inverse_fraction"] <= 1.0
    assert freqs["total_transitions"] > 0


def test_chi_transport_equals_chi_path():
    seq = prime_eabc_sequence(5000)
    classes = [r["class"] for r in seq]
    trans = chi_transport(classes)
    path = chi_path_sliding(classes)
    assert trans["chi_path"] == path["chi_path"]
    assert trans["chi_trans"] == path["chi_path"]


def test_chi_path_vs_hol_structure():
    cmp = chi_path_vs_hol(5000, null_trials=100)
    assert "Pfad" in cmp["comparison"]["verdict"] or "Holonomie" in cmp["comparison"]["verdict"]
    assert "chi_path" in cmp
    assert "chi_hol" in cmp
    assert "null_models" in cmp
    assert "path" in cmp["null_models"]
    assert "hol" in cmp["null_models"]


def test_shuffle_and_isotropy_nulls_bounded():
    seq = prime_eabc_sequence(3000)
    classes = [r["class"] for r in seq]
    shuffle_path = shuffle_null_chi_path(classes, trials=50)
    isotropy_path = isotropy_null_chi_path(classes, trials=50)
    shuffle_hol = shuffle_null_chi_hol(classes, trials=50)
    isotropy_hol = isotropy_null_chi_hol(classes, trials=50)
    assert shuffle_path["null_type"] == "marginal_shuffle"
    assert isotropy_hol["null_type"] == "isotropy_relabel"
    assert -1.0 <= shuffle_path["observed_chi_path"] <= 1.0
    assert -1.0 <= shuffle_hol["observed_chi_hol"] <= 1.0
    # legacy aliases
    assert shuffle_null_chi_E(classes, trials=10)["observed_chi_path"] == shuffle_path["observed_chi_path"]
    assert isotropy_null_chi_E(classes, trials=10)["observed_chi_path"] == isotropy_path["observed_chi_path"]


def test_chi_transport_vs_quadruplet_alias():
    a = chi_path_vs_hol(2000, null_trials=20)
    b = chi_transport_vs_quadruplet(2000, null_trials=20)
    c = chi_sliding_vs_quadruplet(2000, null_trials=20)
    assert a["chi_path"] == b["chi_path"] == c["chi_path"]
    assert a["chi_hol"] == b["chi_hol"]


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
    assert "path_vs_holonomy" in loaded
    assert "chi_path" in loaded["path_vs_holonomy"]
    assert "chi_hol" in loaded["path_vs_holonomy"]
    assert loaded["meta"]["theory"] == "collatz_eabc_zyklus_holonomie.md"
    assert "hol_E_support" in loaded["hol_E_estimates"]
    assert report["output_path"] == str(out)
