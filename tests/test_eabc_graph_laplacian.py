"""Tests für EABC Graph-Laplace L_E und Spektralbericht."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_graph_laplacian import (
  adjacency_matrix,
  eigenvalues_sorted,
  laplacian_directed,
  laplacian_symmetrized,
  run,
  spectral_gap,
  spectral_report,
  spectral_series,
)
from collatz_eabc_sagnac_circulation import circulation_C_E
from collatz_eabc_transition_graph import LABELS, transition_counts


def test_adjacency_row_sums_positive():
  adj, meta = adjacency_matrix(50_000)
  assert adj.shape == (4, 4)
  assert meta["total_edges"] > 0
  assert adj.sum() == meta["total_edges"]


def test_laplacian_row_sums_zero():
  adj, _ = adjacency_matrix(20_000)
  l_dir = laplacian_directed(adj)
  l_sym = laplacian_symmetrized(adj)
  np.testing.assert_allclose(l_dir.sum(axis=1), 0.0, atol=1e-10)
  np.testing.assert_allclose(l_sym.sum(axis=1), 0.0, atol=1e-10)


def test_symmetrized_eigenvalues_real():
  adj, _ = adjacency_matrix(30_000)
  l_sym = laplacian_symmetrized(adj)
  eig = eigenvalues_sorted(l_sym)
  assert eig[0] >= -1e-10
  assert spectral_gap(eig) >= 0.0


def test_spectral_report_links_circulation():
  rep = spectral_report(100_000)
  circ = circulation_C_E(100_000)
  assert rep["circulation"]["C_E"] == circ["C_E"]
  assert rep["circulation"]["D_E"] == circ["D_E"]
  assert rep["relations"]["C_E_equals_D_E"] is True
  assert len(rep["eigenvalues_symmetrized"]) == 4


def test_spectral_series_monotone_limits():
  series = spectral_series([5_000, 20_000, 50_000])
  assert len(series["series"]) == 3
  for row in series["series"]:
    assert row["D_E"] == row["C_E"]


def test_adjacency_matches_transition_counts():
  from collatz_eabc_transition_graph import classes_from_sequence, prime_eabc_sequence

  max_p = 25_000
  seq = prime_eabc_sequence(max_p)
  classes = classes_from_sequence(seq)
  counts = transition_counts(classes)
  adj, _ = adjacency_matrix(max_p)
  np.testing.assert_array_equal(adj, np.array(counts, dtype=float))


def test_run_writes_json(tmp_path: Path):
  out = tmp_path / "laplacian.json"
  report = run(max_p=10_000, output=out)
  assert out.is_file()
  loaded = json.loads(out.read_text(encoding="utf-8"))
  assert loaded["theory"] == "collatz_eabc_zirkulation_spektral.md"
  assert len(loaded["eigenvalues_symmetrized"]) == len(LABELS)
  assert report["output_path"] == str(out)
