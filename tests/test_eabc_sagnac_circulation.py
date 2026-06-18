"""Tests für EABC-Sagnac-Zirkulation C_E(X), Kantenorientierung ω und 1-Form."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_holonomie_fehlerterm import holonomy_counts
from collatz_eabc_sagnac_circulation import (
    CANONICAL_GAP_PATTERN,
    circulation_C_E,
    circulation_report,
    cycle_omega_graph,
    discrete_one_form,
    edge_omega,
    edge_omega_from_gap,
    gap_from_classes,
    line_integral_one_form,
    omega_cycle,
    run,
)
from collatz_eabc_transition_graph import ABCEA_WORD, CEABC_WORD


def test_edge_omega_canonical_cycle():
    assert edge_omega("A", "B") == 1
    assert edge_omega("B", "C") == 1
    assert edge_omega("C", "E") == 1
    assert edge_omega("E", "A") == 1
    assert edge_omega("B", "A") == -1
    assert edge_omega("A", "C") == 0


def test_gap_pattern_canonical_edges():
    for src, dst in (("A", "B"), ("B", "C"), ("C", "E"), ("E", "A")):
        assert gap_from_classes(src, dst) in CANONICAL_GAP_PATTERN
        assert edge_omega_from_gap(src, dst) == edge_omega(src, dst)


def test_discrete_one_form_line_integral():
    form = discrete_one_form()
    assert len(form) == 8
    assert abs(line_integral_one_form(ABCEA_WORD) - 1.0) < 1e-12
    assert abs(line_integral_one_form(CEABC_WORD) + 1.0) < 1e-12


def test_omega_cycle_words():
    assert omega_cycle(ABCEA_WORD) == 1
    assert omega_cycle(CEABC_WORD) == -1
    assert omega_cycle("ABCDE") == 0


def test_cycle_omega_graph_abcea():
    assert cycle_omega_graph(ABCEA_WORD) == 1


def test_circulation_equals_D_E():
    for max_p in (10_000, 50_000, 200_000):
        circ = circulation_C_E(max_p)
        row = holonomy_counts(max_p)
        assert circ["C_E"] == circ["D_E"]
        assert circ["C_E"] == row["Delta_E"]
        assert circ["N_plus"] == row["N_plus"]
        assert circ["N_minus"] == row["N_minus"]
        assert circ["S_E"] == row["S_E"]
        assert circ["C_E_equals_D_E"] is True


def test_circulation_report_structure():
    report = circulation_report(30_000)
    assert report["theory"] == "collatz_eabc_sagnac.md"
    assert "circulation" in report
    assert "graph_orientation" in report
    assert report["circulation"]["C_E_equals_D_E"] is True
    orient = report["graph_orientation"]
    assert orient["ABCEA_omega_cycle"] == 1
    assert orient["CEABC_omega_cycle"] == -1


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "circulation.json"
    report = run(max_p=10_000, output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["circulation"]["C_E"] == loaded["circulation"]["D_E"]
    assert report["output_path"] == str(out)
