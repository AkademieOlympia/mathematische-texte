#!/usr/bin/env python3
"""G_PC^(V4) audit: Dreifachstabilität über N × W × U.

Claim-Wall-Fokus:
  N-Stabilität, T-Fenster-Stabilität (W_k), U-Kernel-Stabilität
  mit Bonferroni über m = |N|·|W|·|U|.

Beispiele:
  # Smoke (schnell)
  python "Topologische Tests.py" --smoke

  # Lokaler Vorbereitungs-Lauf (m=48, reduzierte N / n_iter)
  python "Topologische Tests.py" --local-grid --seed 2026 \\
      --export gp_c_runs/local/triple_local.json \\
      --export-csv gp_c_runs/local/triple_local.csv

  # Voller Grid (CoCalc; m=48, α_Bonf≈1.04e-3)
  python "Topologische Tests.py" --full-grid --seed 2026 --n-iter 1000 \\
      --export gp_c_triple_stability.json --export-csv gp_c_triple_stability.csv

  # Produktive Ordinaten (lcalc-Export)
  python "Topologische Tests.py" --full-grid --zeros-json lcalc_zeros.json \\
      --n-iter 1000 --export gp_c_N5000_real.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import time
from pathlib import Path
from typing import Any

import numpy as np

CONDUCTORS = {"chi_0": 1, "chi_3": 3, "chi_4": 4, "chi_12": 12}

# Kanonisches Raster für Dreifachstabilität (m = 3·4·4 = 48)
DEFAULT_N_LIST = [500, 1000, 5000]
DEFAULT_U_LIST = [0.5, 1.0, 1.5, 2.0]
DEFAULT_HEIGHT_WINDOWS: dict[str, tuple[float, float]] = {
    "W1": (10.0, 200.0),
    "W2": (200.0, 600.0),
    "W3": (600.0, 1200.0),
    "W4": (1200.0, 2500.0),
}
# Lokaler Vorbereitungs-Grid: gleiche Topologie (m=48), kleinere N / n_iter
LOCAL_N_LIST = [200, 400, 800]
LOCAL_DEFAULT_N_ITER = 40


# ---------------------------------------------------------
# 1. Hilfsfunktionen & V4-Hadamard
# ---------------------------------------------------------
def get_hadamard_matrix() -> np.ndarray:
  return (
      np.array(
          [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]],
          dtype=float,
      )
      / 2.0
  )


def triangle_kernel_matrix(diff: np.ndarray, U: float) -> np.ndarray:
  """Vektorisierter Dreieckskern K_U(x) = max(0, 1 - |x|/U)."""
  abs_diff = np.abs(diff)
  out = 1.0 - abs_diff / U
  out[abs_diff > U] = 0.0
  return out


def frobenius_norm(M: np.ndarray) -> float:
  return float(np.linalg.norm(M, "fro"))


def hadamard_isometry_delta(C: np.ndarray, H: np.ndarray) -> dict[str, float]:
  """Prüft ‖C‖_F ≈ ‖H C Hᵀ‖_F (Hadamard orthonormal bis Skala 1)."""
  C_eabc = H @ C @ H.T
  n_c = frobenius_norm(C)
  n_eabc = frobenius_norm(C_eabc)
  return {
      "norm_C": n_c,
      "norm_C_eabc": n_eabc,
      "delta_frobenius": abs(n_c - n_eabc),
  }


# ---------------------------------------------------------
# 2. PC1: Unfolding
# ---------------------------------------------------------
def unfold_spectrum(gamma_list: list[float] | np.ndarray, q_chi: float) -> list[float]:
  """Rohe Ordinaten → entfaltete Einheitsdichte-Punkte (von-Mangoldt)."""
  u_list: list[float] = []
  for g in gamma_list:
    if g > 0:
      u = (g / (2 * math.pi)) * math.log(q_chi * g / (2 * math.pi * math.e))
      u_list.append(float(u))
  return sorted(u_list)


def filter_height_window(
    spectra_raw: dict[str, list[float]], t_lo: float, t_hi: float
) -> dict[str, list[float]]:
  return {
      ch: [g for g in gammas if t_lo <= g < t_hi] for ch, gammas in spectra_raw.items()
  }


def take_first_n(spectra: dict[str, list[float]], N: int) -> dict[str, list[float]]:
  return {ch: gammas[:N] for ch, gammas in spectra.items()}


# ---------------------------------------------------------
# 3. PC4–PC6: Connected Cross-Channel Matrix C
# ---------------------------------------------------------
def compute_connected_cross_matrix(
    spectra_dict: dict[str, list[float]], U: float, channels: list[str]
) -> np.ndarray:
  """4×4 Connected Cross-Correlation Matrix C[chi, psi] (vektorisiert)."""
  n_chan = len(channels)
  C_mat = np.zeros((n_chan, n_chan), dtype=float)
  arrays = {ch: np.asarray(spectra_dict[ch], dtype=float) for ch in channels}

  for i, chi in enumerate(channels):
    for j, psi in enumerate(channels):
      spec_i = arrays[chi]
      spec_j = arrays[psi]
      if len(spec_i) < 2 or len(spec_j) < 2:
        C_mat[i, j] = 0.0
        continue

      diff = spec_i[:, None] - spec_j[None, :]
      sum_val = float(np.sum(triangle_kernel_matrix(diff, U)))

      span = float(spec_i[-1] - spec_i[0])
      density_factor = (len(spec_i) * len(spec_j) / span) if span > 0 else 0.0
      background = U * density_factor

      if chi != psi:
        C_val = (sum_val - background) / max(1, len(spec_i))
      else:
        C_val = sum_val / max(1, len(spec_i))
      C_mat[i, j] = C_val

  return C_mat


# ---------------------------------------------------------
# 4. PC8: Metriken und Monte-Carlo Nullmodell
# ---------------------------------------------------------
def compute_eta_pc(C_mat: np.ndarray) -> float:
  """Relative Off-Diagonal-Energie eta_pc^2."""
  sq = C_mat**2
  den = float(np.sum(sq))
  if den == 0:
    return 0.0
  off = sq.copy()
  np.fill_diagonal(off, 0.0)
  return float(np.sum(off) / den)


def monte_carlo_null_test(
    spectra_dict: dict[str, list[float]],
    U: float,
    channels: list[str],
    n_iter: int = 100,
    rng: np.random.Generator | None = None,
    shift_scale: float = 10.0,
) -> dict[str, Any]:
  """Shift-Shuffling-Nulltest für Z_PC und zweiseitiges MC-p."""
  if rng is None:
    rng = np.random.default_rng()

  emp_C = compute_connected_cross_matrix(spectra_dict, U, channels)
  emp_eta_sq = compute_eta_pc(emp_C)

  null_etas: list[float] = []
  for _ in range(n_iter):
    shuffled: dict[str, list[float]] = {}
    for ch in channels:
      shift = float(rng.uniform(-shift_scale, shift_scale))
      shuffled[ch] = [x + shift for x in spectra_dict[ch]]
    null_C = compute_connected_cross_matrix(shuffled, U, channels)
    null_etas.append(compute_eta_pc(null_C))

  null_arr = np.asarray(null_etas, dtype=float)
  mean_null = float(np.mean(null_arr)) if len(null_arr) else 0.0
  std_null = float(np.std(null_arr, ddof=1)) if len(null_arr) > 1 else 1e-5
  if std_null <= 0:
    std_null = 1e-5
  z_score = (emp_eta_sq - mean_null) / std_null

  abs_dev_emp = abs(emp_eta_sq - mean_null)
  abs_dev_null = np.abs(null_arr - mean_null)
  p_mc = (1 + int(np.count_nonzero(abs_dev_null >= abs_dev_emp))) / (len(null_arr) + 1)

  return {
      "C": emp_C,
      "eta_pc2": emp_eta_sq,
      "null_mean": mean_null,
      "null_std": std_null,
      "z_pc": z_score,
      "p_mc": float(p_mc),
      "n_iter": n_iter,
  }


def bonferroni_correct(p_values: list[float], alpha: float = 0.05) -> dict[str, Any]:
  m = len(p_values)
  if m == 0:
    return {"m": 0, "alpha": alpha, "alpha_corr": alpha, "any_reject": False, "p_corr": []}
  alpha_corr = alpha / m
  p_corr = [min(1.0, m * p) for p in p_values]
  return {
      "m": m,
      "alpha": alpha,
      "alpha_corr": alpha_corr,
      "any_reject": any(p <= alpha_corr for p in p_values),
      "p_corr": p_corr,
  }


# ---------------------------------------------------------
# 5. Spektren laden / synthetisieren
# ---------------------------------------------------------
def synthesize_window_spectra(
    N: int,
    t_lo: float,
    t_hi: float,
    conductors: dict[str, int],
    rng: np.random.Generator,
) -> dict[str, list[float]]:
  """N synthetische Ordinaten pro Kanal im Fenster [t_lo, t_hi)."""
  spectra: dict[str, list[float]] = {}
  for ch in conductors:
    raw = rng.uniform(t_lo, t_hi, size=N)
    spectra[ch] = sorted(float(x) for x in raw)
  return spectra


def load_zeros_json(path: Path) -> dict[str, list[float]]:
  data = json.loads(path.read_text(encoding="utf-8"))
  out: dict[str, list[float]] = {}
  for ch in CONDUCTORS:
    if ch not in data:
      raise KeyError(f"Kanal {ch!r} fehlt in {path}")
    out[ch] = sorted(float(x) for x in data[ch])
  return out


def parse_float_list(s: str) -> list[float]:
  return [float(x.strip()) for x in s.split(",") if x.strip()]


def parse_int_list(s: str) -> list[int]:
  return [int(x.strip()) for x in s.split(",") if x.strip()]


def parse_named_windows(s: str) -> dict[str, tuple[float, float]]:
  """'W1:10:200,W2:200:600' oder '10:200,200:600' → benannte Fenster."""
  windows: dict[str, tuple[float, float]] = {}
  anon = 0
  for part in s.split(","):
    part = part.strip()
    if not part:
      continue
    bits = part.split(":")
    if len(bits) == 3:
      name, lo_s, hi_s = bits
    elif len(bits) == 2:
      anon += 1
      name, lo_s, hi_s = f"W{anon}", bits[0], bits[1]
    else:
      raise ValueError(f"Ungültiges Fenster: {part!r}")
    windows[name] = (float(lo_s), float(hi_s))
  return windows


# ---------------------------------------------------------
# 6. Dreifachstabilitäts-Raster: N × W × U
# ---------------------------------------------------------
def run_triple_stability_audit(
    N_list: list[int],
    height_windows: dict[str, tuple[float, float]],
    U_list: list[float],
    n_iter: int,
    seed: int,
    alpha: float,
    zeros_global: dict[str, list[float]] | None = None,
    min_counts: int = 20,
    store_matrices: bool = False,
) -> dict[str, Any]:
  """Vollständige Grid-Expansion mit Bonferroni über m = |N|·|W|·|U|."""
  channels = list(CONDUCTORS.keys())
  H = get_hadamard_matrix()
  rng = np.random.default_rng(seed)

  m_planned = len(N_list) * len(height_windows) * len(U_list)
  alpha_bonf_planned = alpha / m_planned if m_planned else alpha

  print("=== STARTE VOLLSTÄNDIGEN GRID-LAUF (N × W × U) ===")
  print(f"Testzellen geplant (m): {m_planned}")
  print(f"Bonferroni-Schwellenwert (α_Bonf): {alpha_bonf_planned:.6e}")
  print(f"N={N_list}  |  W={list(height_windows)}  |  U={U_list}")
  print(f"n_iter={n_iter}  |  seed={seed}\n")

  cells: list[dict[str, Any]] = []
  cell_counter = 0

  for N in N_list:
    for w_name, (t_lo, t_hi) in height_windows.items():
      if zeros_global is not None:
        raw_w = filter_height_window(zeros_global, t_lo, t_hi)
        # N-Stabilität: Prefix der Fenster-Ordinaten
        available = min(len(v) for v in raw_w.values())
        if available < min_counts:
          cell_counter += len(U_list)
          cells.append(
              {
                  "N": N,
                  "window": w_name,
                  "W": [t_lo, t_hi],
                  "U": None,
                  "status": "skip_too_few",
                  "counts": {ch: len(raw_w[ch]) for ch in channels},
              }
          )
          continue
        n_use = min(N, available)
        raw_w = take_first_n(raw_w, n_use)
      else:
        raw_w = synthesize_window_spectra(N, t_lo, t_hi, CONDUCTORS, rng)
        n_use = N

      unfolded = {
          ch: unfold_spectrum(raw_w[ch], CONDUCTORS[ch]) for ch in channels
      }
      counts = {ch: len(unfolded[ch]) for ch in channels}

      for U in U_list:
        cell_counter += 1
        cell_seed = int(rng.integers(0, 2**31 - 1))
        cell_rng = np.random.default_rng(cell_seed)
        t0 = time.perf_counter()
        result = monte_carlo_null_test(
            unfolded, U, channels, n_iter=n_iter, rng=cell_rng
        )
        elapsed = time.perf_counter() - t0

        C = result["C"]
        C_eabc = H @ C @ H.T
        iso = hadamard_isometry_delta(C, H)

        record: dict[str, Any] = {
            "N": N,
            "N_used": n_use,
            "window": w_name,
            "W": [t_lo, t_hi],
            "U": U,
            "status": "ok",
            "counts": counts,
            "eta_pc2": result["eta_pc2"],
            "null_mean": result["null_mean"],
            "null_std": result["null_std"],
            "z_pc": result["z_pc"],
            "p_mc": result["p_mc"],
            "delta_frobenius": iso["delta_frobenius"],
            "norm_C": iso["norm_C"],
            "norm_C_eabc": iso["norm_C_eabc"],
            "cell_seed": cell_seed,
            "elapsed_s": elapsed,
            "cell_index": cell_counter,
        }
        if store_matrices:
          record["C"] = C.tolist()
          record["C_eabc"] = C_eabc.tolist()
        cells.append(record)

        print(
            f"[{cell_counter}/{m_planned}] "
            f"N={N} | W={w_name} | U={U:g} | "
            f"η²={result['eta_pc2']:.4f} | "
            f"Z={result['z_pc']:+.2f} | "
            f"p={result['p_mc']:.4f} | "
            f"ΔF={iso['delta_frobenius']:.1e} | "
            f"{elapsed:.1f}s"
        )

  ok_cells = [c for c in cells if c.get("status") == "ok"]
  p_values = [c["p_mc"] for c in ok_cells]
  bonf = bonferroni_correct(p_values, alpha=alpha)
  for c, p_corr in zip(ok_cells, bonf["p_corr"]):
    c["p_bonferroni"] = p_corr
    c["reject_bonferroni"] = c["p_mc"] <= bonf["alpha_corr"]

  stability = summarize_triple_stability(ok_cells)

  print(f"\n=== GRID-LAUF VOLLSTÄNDIG: {len(ok_cells)} Zellen evaluiert ===")
  print(
      f"Bonferroni: m={bonf['m']}, α_Bonf={bonf['alpha_corr']:.6e}, "
      f"family-wise reject={bonf['any_reject']}"
  )
  print(
      f"Max |Z_PC|={stability['max_abs_z']:.2f}  |  "
      f"Max ΔF={stability['max_delta_frobenius']:.2e}  |  "
      f"Reject-Zellen={stability['n_reject']}"
  )

  return {
      "schema": "gp_c_v4_triple_stability.v1",
      "seed": seed,
      "n_iter": n_iter,
      "alpha": alpha,
      "N_list": N_list,
      "U_list": U_list,
      "windows": {k: list(v) for k, v in height_windows.items()},
      "m_planned": m_planned,
      "m_tests": bonf["m"],
      "bonferroni": {
          "alpha_corr": bonf["alpha_corr"],
          "any_reject": bonf["any_reject"],
      },
      "stability": stability,
      "cells": cells,
      "source": "zeros-json" if zeros_global is not None else "synthetic-per-window",
  }


def summarize_triple_stability(ok_cells: list[dict[str, Any]]) -> dict[str, Any]:
  """Kompakte Stabilitäts-Kennzahlen über N / W / U."""
  if not ok_cells:
    return {
        "max_abs_z": 0.0,
        "max_delta_frobenius": 0.0,
        "n_reject": 0,
        "by_N": {},
        "by_W": {},
        "by_U": {},
    }

  def group_stats(key: str) -> dict[str, Any]:
    groups: dict[Any, list[dict[str, Any]]] = {}
    for c in ok_cells:
      groups.setdefault(c[key], []).append(c)
    out: dict[str, Any] = {}
    for k, rows in groups.items():
      zs = [abs(r["z_pc"]) for r in rows]
      out[str(k)] = {
          "n": len(rows),
          "mean_abs_z": float(np.mean(zs)),
          "max_abs_z": float(np.max(zs)),
          "mean_eta": float(np.mean([r["eta_pc2"] for r in rows])),
          "n_reject": sum(1 for r in rows if r.get("reject_bonferroni")),
      }
    return out

  return {
      "max_abs_z": float(max(abs(c["z_pc"]) for c in ok_cells)),
      "max_delta_frobenius": float(max(c["delta_frobenius"] for c in ok_cells)),
      "n_reject": sum(1 for c in ok_cells if c.get("reject_bonferroni")),
      "by_N": group_stats("N"),
      "by_W": group_stats("window"),
      "by_U": group_stats("U"),
  }


def export_csv(path: Path, cells: list[dict[str, Any]]) -> None:
  fields = [
      "cell_index",
      "N",
      "N_used",
      "window",
      "W_lo",
      "W_hi",
      "U",
      "status",
      "eta_pc2",
      "null_mean",
      "null_std",
      "z_pc",
      "p_mc",
      "p_bonferroni",
      "reject_bonferroni",
      "delta_frobenius",
      "elapsed_s",
      "cell_seed",
  ]
  with path.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for c in cells:
      if c.get("status") != "ok":
        continue
      w = c.get("W") or [None, None]
      writer.writerow(
          {
              "cell_index": c.get("cell_index"),
              "N": c.get("N"),
              "N_used": c.get("N_used"),
              "window": c.get("window"),
              "W_lo": w[0],
              "W_hi": w[1],
              "U": c.get("U"),
              "status": c.get("status"),
              "eta_pc2": c.get("eta_pc2"),
              "null_mean": c.get("null_mean"),
              "null_std": c.get("null_std"),
              "z_pc": c.get("z_pc"),
              "p_mc": c.get("p_mc"),
              "p_bonferroni": c.get("p_bonferroni"),
              "reject_bonferroni": c.get("reject_bonferroni"),
              "delta_frobenius": c.get("delta_frobenius"),
              "elapsed_s": c.get("elapsed_s"),
              "cell_seed": c.get("cell_seed"),
          }
      )


def build_argparser() -> argparse.ArgumentParser:
  p = argparse.ArgumentParser(
      description="G_PC^(V4) Dreifachstabilität: N × W × U + Bonferroni"
  )
  p.add_argument(
      "--full-grid",
      action="store_true",
      help="Kanonisches Raster N∈{500,1000,5000} × W1..W4 × U∈{0.5,1,1.5,2}",
  )
  p.add_argument(
      "--local-grid",
      action="store_true",
      help=(
          "Lokaler Prep-Lauf: gleiche W×U-Topologie (m=48), "
          f"N∈{LOCAL_N_LIST}, Default n_iter={LOCAL_DEFAULT_N_ITER}"
      ),
  )
  p.add_argument(
      "--N-list",
      type=str,
      default=None,
      help="Kommaliste der Stichprobengrößen, z.B. 500,1000,5000",
  )
  p.add_argument("--N", type=int, default=300, help="Einzel-N (wenn kein --N-list/--full-grid)")
  p.add_argument(
      "--W",
      type=str,
      default=None,
      help="Fenster: W1:10:200,W2:200:600 oder 10:200,200:600",
  )
  p.add_argument("--U", type=str, default="1.5", help="Kernelbreiten als Kommaliste")
  p.add_argument("--seed", type=int, default=2026)
  p.add_argument(
      "--n-iter",
      type=int,
      default=None,
      help="MC-Shuffling-Iterationen (Defaults: smoke=8, local=40, sonst=50, full oft 200–1000)",
  )
  p.add_argument("--bonferroni", action="store_true", help="(immer aktiv; Flag dokumentarisch)")
  p.add_argument("--alpha", type=float, default=0.05)
  p.add_argument("--zeros-json", type=Path, default=None)
  p.add_argument("--export", type=Path, default=None, help="JSON-Report")
  p.add_argument("--export-csv", type=Path, default=None, help="CSV der Zellen")
  p.add_argument(
      "--store-matrices",
      action="store_true",
      help="C und C_EABC in JSON mitschreiben",
  )
  p.add_argument(
      "--smoke",
      action="store_true",
      help="Mini-Grid: N=60,80 × W1,W2 × U=1.0,1.5 × n_iter=8",
  )
  return p


def main(argv: list[str] | None = None) -> int:
  args = build_argparser().parse_args(argv)
  _ = args.bonferroni

  if args.smoke:
    N_list = [60, 80]
    height_windows = {"W1": (10.0, 120.0), "W2": (120.0, 250.0)}
    U_list = [1.0, 1.5]
    if args.n_iter is None:
      args.n_iter = 8
  elif args.local_grid:
    N_list = list(LOCAL_N_LIST)
    height_windows = dict(DEFAULT_HEIGHT_WINDOWS)
    U_list = list(DEFAULT_U_LIST)
    if args.n_iter is None:
      args.n_iter = LOCAL_DEFAULT_N_ITER
  elif args.full_grid:
    N_list = list(DEFAULT_N_LIST)
    height_windows = dict(DEFAULT_HEIGHT_WINDOWS)
    U_list = list(DEFAULT_U_LIST)
    if args.n_iter is None:
      args.n_iter = 200
  else:
    N_list = parse_int_list(args.N_list) if args.N_list else [args.N]
    if args.W:
      height_windows = parse_named_windows(args.W)
    else:
      height_windows = {"W1": (10.0, 500.0)}
    U_list = parse_float_list(args.U)
    if args.n_iter is None:
      args.n_iter = 50

  random.seed(args.seed)
  zeros_global = load_zeros_json(args.zeros_json) if args.zeros_json else None

  report = run_triple_stability_audit(
      N_list=N_list,
      height_windows=height_windows,
      U_list=U_list,
      n_iter=args.n_iter,
      seed=args.seed,
      alpha=args.alpha,
      zeros_global=zeros_global,
      store_matrices=args.store_matrices,
  )

  report["profile"] = (
      "smoke"
      if args.smoke
      else "local-grid"
      if args.local_grid
      else "full-grid"
      if args.full_grid
      else "custom"
  )

  if args.export is not None:
    args.export.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"\nExport JSON: {args.export}")

  if args.export_csv is not None:
    export_csv(args.export_csv, report["cells"])
    print(f"Export CSV:  {args.export_csv}")

  return 0


if __name__ == "__main__":
  raise SystemExit(main())
