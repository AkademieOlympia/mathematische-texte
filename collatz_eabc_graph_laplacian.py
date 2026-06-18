#!/usr/bin/env python3
"""
EABC-Transportgraph G_E und Graph-Laplace L_E: Spektralgeometrie der Primfolge.

Theorie: collatz_eabc_zirkulation_spektral.md §6–§7

  G_E aus Primfolge-Übergängen κ(p_n) → κ(p_{n+1})
  A_E gewichtete Adjazenz, L_E = D_out - A_E (gerichtet)
  L_E^sym symmetrisiert für reelle Eigenwerte
  Spec(L_E), Spektrallücke λ_1 - λ_0
  Bezug zu C_E, D_E aus collatz_eabc_sagnac_circulation

Ausführung:
    python3 collatz_eabc_graph_laplacian.py
    python3 collatz_eabc_graph_laplacian.py --max-p 100000
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from collatz_eabc_sagnac_circulation import circulation_C_E
from collatz_eabc_transition_graph import (
    LABELS,
    classes_from_sequence,
    prime_eabc_sequence,
    transition_counts,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_graph_laplacian.json"
THEORY = "collatz_eabc_zirkulation_spektral.md"
N_VERTICES = len(LABELS)


def adjacency_matrix(max_p: int) -> tuple[np.ndarray, dict[str, Any]]:
  """Gewichtete Adjazenz A_E aus Primfolge-Übergängen bis max_p."""
  seq = prime_eabc_sequence(max_p)
  classes = classes_from_sequence(seq)
  counts = transition_counts(classes)
  mat = np.array(counts, dtype=float)
  meta = {
    "max_p": max_p,
    "prime_count": len(seq),
    "transition_counts": {
      LABELS[i]: {LABELS[j]: int(counts[i][j]) for j in range(N_VERTICES)}
      for i in range(N_VERTICES)
    },
    "total_edges": int(mat.sum()),
  }
  return mat, meta


def laplacian_directed(adj: np.ndarray) -> np.ndarray:
  """L_E = D_out - A_E (gerichteter Graph-Laplace)."""
  out_deg = adj.sum(axis=1)
  return np.diag(out_deg) - adj


def laplacian_symmetrized(adj: np.ndarray) -> np.ndarray:
  """L_E^sym = D - (A + A^T)/2 mit D = diag(A 1)."""
  sym_adj = 0.5 * (adj + adj.T)
  deg = sym_adj.sum(axis=1)
  return np.diag(deg) - sym_adj


def eigenvalues_sorted(laplacian: np.ndarray) -> np.ndarray:
  """Reelle Eigenwerte aufsteigend (symmetrische Matrix vorausgesetzt)."""
  vals = np.linalg.eigvalsh(laplacian)
  return np.sort(vals)


def spectral_gap(eigenvalues: np.ndarray) -> float:
  """λ_1 - λ_0 für sortierte Eigenwerte."""
  if len(eigenvalues) < 2:
    return 0.0
  return float(eigenvalues[1] - eigenvalues[0])


def spectral_report(max_p: int) -> dict[str, Any]:
  """Vollständiger Spektralbericht inkl. Zirkulations-Bezug."""
  adj, graph_meta = adjacency_matrix(max_p)
  l_dir = laplacian_directed(adj)
  l_sym = laplacian_symmetrized(adj)
  eig_dir = np.sort(np.linalg.eigvals(l_dir).real)
  eig_sym = eigenvalues_sorted(l_sym)
  circ = circulation_C_E(max_p)

  return {
    "X": max_p,
    "theory": THEORY,
    "graph": graph_meta,
    "adjacency": adj.tolist(),
    "laplacian_directed": l_dir.tolist(),
    "laplacian_symmetrized": l_sym.tolist(),
    "eigenvalues_directed": [float(x) for x in eig_dir],
    "eigenvalues_symmetrized": [float(x) for x in eig_sym],
    "spectral_gap_directed": spectral_gap(eig_dir),
    "spectral_gap_symmetrized": spectral_gap(eig_sym),
    "circulation": {
      "C_E": circ["C_E"],
      "D_E": circ["D_E"],
      "S_E": circ["S_E"],
      "N_plus": circ["N_plus"],
      "N_minus": circ["N_minus"],
      "detected_cycles": circ["detected_cycles"],
    },
    "relations": {
      "C_E_equals_D_E": circ["C_E_equals_D_E"],
      "laplacian_formula": "L_E = D_out - A_E",
      "sym_formula": "L_E^sym = D - (A + A^T)/2",
    },
    "epistemic": {
      "Spec_L_E": "Definition",
      "prim_spectral_anomaly": "Hypothese",
      "C_E_D_E_link": "Definition",
    },
  }


def spectral_series(limits: list[int]) -> dict[str, Any]:
  """Spektrallücke und D_E als Funktion von X."""
  rows: list[dict[str, Any]] = []
  for x in limits:
    rep = spectral_report(x)
    rows.append(
      {
        "X": x,
        "spectral_gap_sym": rep["spectral_gap_symmetrized"],
        "spectral_gap_dir": rep["spectral_gap_directed"],
        "lambda_min_sym": rep["eigenvalues_symmetrized"][0],
        "lambda_max_sym": rep["eigenvalues_symmetrized"][-1],
        "C_E": rep["circulation"]["C_E"],
        "D_E": rep["circulation"]["D_E"],
        "S_E": rep["circulation"]["S_E"],
      }
    )
  return {"limits": limits, "series": rows}


def run(max_p: int = 100_000, output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
  report = spectral_report(max_p)
  report["series_sample"] = spectral_series(
    [10**3, 10**4, min(max_p, 10**5)] if max_p >= 10**3 else [max_p]
  )
  output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
  report["output_path"] = str(output)
  return report


def main() -> None:
  parser = argparse.ArgumentParser(description="EABC Graph-Laplace Spec(L_E)")
  parser.add_argument("--max-p", type=int, default=100_000)
  parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
  args = parser.parse_args()
  report = run(max_p=args.max_p, output=args.output)
  eig = report["eigenvalues_symmetrized"]
  circ = report["circulation"]
  print("=== EABC Graph-Laplace Spec(L_E) ===")
  print(f"X={report['X']}: edges={report['graph']['total_edges']}")
  print(f"Spec(L_E^sym) = {[round(x, 4) for x in eig]}")
  print(f"gap_sym = {report['spectral_gap_symmetrized']:.4f}")
  print(f"C_E={circ['C_E']:+d}  D_E={circ['D_E']:+d}  S_E={circ['S_E']:+.4f}")
  print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
  main()
