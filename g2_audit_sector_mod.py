#!/usr/bin/env python3
"""V4-Sektor-Modulation der Assoziator-Sequenz α_n (64 Tripel-Kanäle).

Zerlegt den globalen Effekt in Kanäle (S_n, S_{n+1}, S_{n+2}) mit
S ∈ {E, A, B, C} ≅ (ℤ/12ℤ)^× ≅ V4  via p mod 12 ∈ {1,5,7,11}.

Voraussetzung: primes_107.npy

Beispiele:
  python g2_audit_sector_mod.py
  python g2_audit_sector_mod.py --min-count 500 --export scripts/gp_sector_v4_107.csv
"""

from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent

SECTOR_NAMES = ("E", "A", "B", "C")  # 1,5,7,11 mod 12
REM_TO_SECTOR = {1: 0, 5: 1, 7: 2, 11: 3}

# --- Fano ---
fano_table = np.zeros((8, 8), dtype=int)
_triples = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
for i, j, k in _triples:
  fano_table[i, j] = k
  fano_table[j, i] = -k
  fano_table[j, k] = i
  fano_table[k, j] = -i
  fano_table[k, i] = j
  fano_table[i, k] = -j
for i in range(1, 8):
  fano_table[i, i] = -1


def resolve_cache() -> Path:
  candidates = [
      SCRIPT_DIR / "primes_107.npy",
      Path.home() / "Desktop" / "Mathematische Texte" / "primes_107.npy",
      Path.home() / "Projects" / "Kepler-Hurrwitz" / "scripts" / "primes_107.npy",
      SCRIPT_DIR.parent / "scripts" / "primes_107.npy",
  ]
  for p in candidates:
    if p.exists():
      return p
  raise FileNotFoundError("Cache 'primes_107.npy' nicht gefunden.")


def get_primes_cached() -> np.ndarray:
  path = resolve_cache()
  print(f"[Cache] {path}")
  return np.load(path).astype(np.int64)


def oct_mul_batch(o1: np.ndarray, o2: np.ndarray) -> np.ndarray:
  res = np.zeros((len(o1), 8), dtype=float)
  res[:, 0] = o1[:, 0] * o2[:, 0] - np.sum(o1[:, 1:] * o2[:, 1:], axis=1)
  res[:, 1:] = o1[:, 0:1] * o2[:, 1:] + o2[:, 0:1] * o1[:, 1:]
  for i in range(1, 8):
    for j in range(1, 8):
      if i == j:
        continue
      target = int(fano_table[i, j])
      if target == 0:
        continue
      idx = abs(target)
      sgn = 1.0 if target > 0 else -1.0
      res[:, idx] += sgn * o1[:, i] * o2[:, j]
  return res


def fano_index(p_arr: np.ndarray) -> np.ndarray:
  ind = (p_arr % 7).astype(np.int64)
  ind[ind == 0] = 7
  return ind


def compute_alphas_batch(primes_arr: np.ndarray) -> np.ndarray:
  primes_arr = np.asarray(primes_arr, dtype=np.int64)
  n = len(primes_arr)
  if n < 4:
    return np.array([], dtype=float)

  g_all = np.diff(primes_arr).astype(float) / np.log(primes_arr[:-1].astype(float))
  n_triples = n - 3
  g_x = g_all[:n_triples]
  g_y = g_all[1 : n_triples + 1]
  g_z = g_all[2 : n_triples + 2]

  ind_x = fano_index(primes_arr[:n_triples])
  ind_y = fano_index(primes_arr[1 : n_triples + 1])
  ind_z = fano_index(primes_arr[2 : n_triples + 2])

  x_vec = np.zeros((n_triples, 8), dtype=float)
  y_vec = np.zeros((n_triples, 8), dtype=float)
  z_vec = np.zeros((n_triples, 8), dtype=float)
  for k in range(1, 8):
    mx = ind_x == k
    my = ind_y == k
    mz = ind_z == k
    x_vec[mx, k] = g_x[mx]
    y_vec[my, k] = g_y[my]
    z_vec[mz, k] = g_z[mz]

  xy = oct_mul_batch(x_vec, y_vec)
  yz = oct_mul_batch(y_vec, z_vec)
  ass = oct_mul_batch(xy, z_vec) - oct_mul_batch(x_vec, yz)
  return np.linalg.norm(ass, axis=1)


def v4_sectors(primes: np.ndarray) -> np.ndarray:
  """Sektor-Index 0..3 für p≡1,5,7,11 (mod 12); sonst −1 (z. B. 2,3)."""
  rem = primes % 12
  out = np.full(len(primes), -1, dtype=np.int8)
  for r, s in REM_TO_SECTOR.items():
    out[rem == r] = s
  return out


def channel_label(t: tuple[int, int, int]) -> str:
  return "→".join(SECTOR_NAMES[s] for s in t)


def run_sector_audit(
    primes: np.ndarray,
    alphas: np.ndarray,
    min_count: int,
    z_hit: float,
) -> tuple[list[dict[str, object]], dict[str, float | int]]:
  """Gruppiert α_i nach (S_i, S_{i+1}, S_{i+2})."""
  n = len(alphas)
  sectors = v4_sectors(primes)
  s0 = sectors[:n]
  s1 = sectors[1 : n + 1]
  s2 = sectors[2 : n + 2]

  valid = (s0 >= 0) & (s1 >= 0) & (s2 >= 0)
  # Packe Tripel als Base-4 Code 0..63
  code = (s0.astype(np.int32) * 16 + s1.astype(np.int32) * 4 + s2.astype(np.int32))
  code = code[valid]
  a_valid = alphas[valid]

  global_mean = float(np.mean(a_valid))
  global_std = float(np.std(a_valid, ddof=1)) if len(a_valid) > 1 else 0.0

  rows: list[dict[str, object]] = []
  for c in range(64):
    mask = code == c
    vals = a_valid[mask]
    count = int(len(vals))
    t = ((c // 16) % 4, (c // 4) % 4, c % 4)
    if count == 0:
      continue
    mean = float(np.mean(vals))
    std = float(np.std(vals, ddof=1)) if count > 1 else 0.0
    sem = std / np.sqrt(count) if count > 0 and std > 0 else 0.0
    # Z gegen globalen Mittelwert (SEM des Kanals)
    z = (mean - global_mean) / sem if sem > 0 else 0.0
    # Anteil assoziativer Tripel
    occ = float(np.mean(vals < 1e-5))
    rows.append(
        {
            "code": c,
            "channel": channel_label(t),
            "tripel": t,
            "N": count,
            "mean": mean,
            "std": std,
            "occ0": occ,
            "delta": mean - global_mean,
            "z": z,
            "pass_min": count >= min_count,
        }
    )

  # Bonferroni über Kanäle mit N ≥ min_count
  tested = [r for r in rows if r["pass_min"]]
  m = len(tested)
  alpha_bonf = 0.05 / m if m else 0.05
  from math import erfc, sqrt

  def normal_sf_two_sided(z: float) -> float:
    return float(erfc(abs(z) / sqrt(2.0)))

  for r in rows:
    if r["pass_min"]:
      p = normal_sf_two_sided(float(r["z"]))
      r["p"] = p
      r["p_bonf"] = min(1.0, m * p)
      r["reject_bonf"] = p <= alpha_bonf
      r["hit"] = abs(float(r["z"])) >= z_hit
    else:
      r["p"] = float("nan")
      r["p_bonf"] = float("nan")
      r["reject_bonf"] = False
      r["hit"] = False

  meta = {
      "global_mean": global_mean,
      "global_std": global_std,
      "n_valid": int(valid.sum()),
      "n_excluded": int((~valid).sum()),
      "m_tested": m,
      "alpha_bonf": alpha_bonf,
  }
  return rows, meta


def build_argparser() -> argparse.ArgumentParser:
  p = argparse.ArgumentParser(description="V4-Sektor-Modulation der α_n (64 Kanäle)")
  p.add_argument("--min-count", type=int, default=1000, help="Min. Tripel pro Kanal")
  p.add_argument("--z-hit", type=float, default=10.0, help="|Z|-Schwelle für 'Fund'-Markierung")
  p.add_argument("--top", type=int, default=20, help="Top-|Z| Kanäle in der Konsole")
  p.add_argument("--export", type=Path, default=None, help="CSV aller 64 Kanäle")
  return p


def main(argv: list[str] | None = None) -> int:
  args = build_argparser().parse_args(argv)
  primes = get_primes_cached()
  print(f"[Daten] {len(primes):,} Primzahlen")

  t0 = time.time()
  alphas = compute_alphas_batch(primes)
  print(f"[α] {len(alphas):,} Tripel in {time.time() - t0:.2f}s")

  rows, meta = run_sector_audit(primes, alphas, args.min_count, args.z_hit)
  tested = [r for r in rows if r["pass_min"]]
  tested.sort(key=lambda r: abs(float(r["z"])), reverse=True)

  print("\n" + "=" * 72)
  print(" V4-SEKTOR-MODULATION (64 Kanäle E,A,B,C)")
  print("=" * 72)
  print(
      f" Global ⟨α⟩ = {meta['global_mean']:.6f} ± {meta['global_std']:.6f}  |  "
      f"gültig={meta['n_valid']:,}  exkl.(2/3)={meta['n_excluded']}"
  )
  print(
      f" Getestet (N≥{args.min_count}): m={meta['m_tested']}  |  "
      f"α_Bonf={meta['alpha_bonf']:.3e}  |  |Z|-Hit≥{args.z_hit:g}"
  )
  print("-" * 72)
  print(f"{'Kanal (V4)':<14} {'N':>8} {'Mean':>10} {'Δ':>10} {'Z':>8} {'P(α≈0)':>8} {'Bonf':>6}")
  print("-" * 72)

  show = tested[: args.top]
  for r in show:
    bonf = "YES" if r["reject_bonf"] else ("HIT" if r["hit"] else "")
    print(
        f"{r['channel']:<14} {r['N']:>8} {r['mean']:>10.4f} {r['delta']:>+10.4f} "
        f"{r['z']:>+8.2f} {100*r['occ0']:>7.2f}% {bonf:>6}"
    )

  n_hit = sum(1 for r in tested if r["hit"])
  n_bonf = sum(1 for r in tested if r["reject_bonf"])
  print("-" * 72)
  print(f" Kanäle mit |Z|≥{args.z_hit:g}: {n_hit}  |  Bonferroni-Reject: {n_bonf}")
  if tested:
    best = tested[0]
    print(
        f" Top-|Z|: {best['channel']}  N={best['N']}  mean={best['mean']:.4f}  "
        f"Z={best['z']:+.2f}"
    )
  print("=" * 72)

  if args.export is not None:
    args.export.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "channel",
        "code",
        "N",
        "mean",
        "std",
        "delta",
        "z",
        "occ0",
        "p",
        "p_bonf",
        "reject_bonf",
        "hit",
        "pass_min",
    ]
    with args.export.open("w", encoding="utf-8", newline="") as f:
      w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
      w.writeheader()
      for r in sorted(rows, key=lambda x: -abs(float(x["z"]))):
        w.writerow(r)
    print(f"[Export] {args.export}")

  return 0


if __name__ == "__main__":
  raise SystemExit(main())
