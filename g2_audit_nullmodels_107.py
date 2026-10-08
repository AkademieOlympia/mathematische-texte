#!/usr/bin/env python3
"""G2-Assoziator-Audit bis 10^7 mit hierarchischen Nullmodellen N1/N2.

Voraussetzung: primes_107.npy (z. B. aus g2_audit_107.py).

Nullmodelle (Gap-Ebene):
  N1 / Markov-1: P(g_{n+1} | g_n) aus empirischen Gap-Übergängen
  N2 / Markov-2: P(g_{n+2} | g_n, g_{n+1}) aus empirischen Paar-Übergängen

Metriken:
  Δ_M1 = ⟨α⟩_emp − E_N1[⟨α⟩],  Z_M1 = Δ_M1 / σ_N1
  Δ_M2 = ⟨α⟩_emp − E_N2[⟨α⟩],  Z_M2 = Δ_M2 / σ_N2

Beispiele:
  python g2_audit_nullmodels_107.py
  python g2_audit_nullmodels_107.py --n-iter 50 --seed 2026
"""

from __future__ import annotations

import argparse
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
LIMIT = 10_000_000

# --- Fano-Struktur ---
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
  raise FileNotFoundError(
      "Cache 'primes_107.npy' nicht gefunden. Bitte zuerst g2_audit_107.py ausführen."
  )


def get_primes_cached() -> np.ndarray:
  path = resolve_cache()
  print(f"[Cache] Lade Primzahlen aus '{path}'...")
  t0 = time.time()
  primes = np.load(path)
  print(f"[Cache] {len(primes):,} Primzahlen in {time.time() - t0:.3f} s.")
  return primes.astype(np.int64)


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
  """α_n = ‖[x,y,z]‖ mit gap-gewichteter Fano-Einbettung (N = |P|−3)."""
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


def rebuild_primes(start: int, gaps: np.ndarray) -> np.ndarray:
  out = np.empty(len(gaps) + 1, dtype=np.int64)
  out[0] = start
  out[1:] = start + np.cumsum(gaps, dtype=np.int64)
  return out


def build_markov1_successors(gaps: np.ndarray) -> dict[int, np.ndarray]:
  succ: dict[int, list[int]] = defaultdict(list)
  for a, b in zip(gaps[:-1], gaps[1:]):
    succ[int(a)].append(int(b))
  return {k: np.asarray(v, dtype=np.int64) for k, v in succ.items()}


def build_markov2_successors(gaps: np.ndarray) -> dict[tuple[int, int], np.ndarray]:
  succ: dict[tuple[int, int], list[int]] = defaultdict(list)
  for a, b, c in zip(gaps[:-2], gaps[1:-1], gaps[2:]):
    succ[(int(a), int(b))].append(int(c))
  return {k: np.asarray(v, dtype=np.int64) for k, v in succ.items()}


def sample_markov1(
    gaps: np.ndarray,
    succ: dict[int, np.ndarray],
    n: int,
    rng: np.random.Generator,
) -> np.ndarray:
  out = np.empty(n, dtype=np.int64)
  out[0] = int(gaps[rng.integers(0, len(gaps))])
  for i in range(1, n):
    options = succ.get(int(out[i - 1]))
    if options is None or len(options) == 0:
      out[i] = int(gaps[rng.integers(0, len(gaps))])
    else:
      out[i] = int(options[rng.integers(0, len(options))])
  return out


def sample_markov2(
    gaps: np.ndarray,
    succ: dict[tuple[int, int], np.ndarray],
    n: int,
    rng: np.random.Generator,
) -> np.ndarray:
  out = np.empty(n, dtype=np.int64)
  j = int(rng.integers(0, len(gaps) - 1))
  out[0] = int(gaps[j])
  out[1] = int(gaps[j + 1])
  for i in range(2, n):
    options = succ.get((int(out[i - 2]), int(out[i - 1])))
    if options is None or len(options) == 0:
      out[i] = int(gaps[rng.integers(0, len(gaps))])
    else:
      out[i] = int(options[rng.integers(0, len(options))])
  return out


def monte_carlo_null(
    primes: np.ndarray,
    gaps: np.ndarray,
    model: str,
    n_iter: int,
    rng: np.random.Generator,
) -> dict[str, float | np.ndarray]:
  start_p = int(primes[0])
  n_gaps = len(gaps)
  means: list[float] = []

  if model == "N1":
    succ1 = build_markov1_successors(gaps)
    print(f"[Nullmodell N1] {n_iter} MC-Läufe (Markov-1, {len(succ1)} Zustände)...")
    sampler = lambda r: sample_markov1(gaps, succ1, n_gaps, r)
  elif model == "N2":
    succ2 = build_markov2_successors(gaps)
    print(f"[Nullmodell N2] {n_iter} MC-Läufe (Markov-2, {len(succ2)} Zustände)...")
    sampler = lambda r: sample_markov2(gaps, succ2, n_gaps, r)
  else:
    raise ValueError(model)

  t0 = time.time()
  for i in range(n_iter):
    child = rng.spawn(1)[0]
    synth_gaps = sampler(child)
    synth_primes = rebuild_primes(start_p, synth_gaps)
    alphas = compute_alphas_batch(synth_primes)
    means.append(float(np.mean(alphas)) if len(alphas) else 0.0)
    if (i + 1) % max(1, n_iter // 5) == 0 or i + 1 == n_iter:
      print(f"  [{model}] {i + 1}/{n_iter}  ⟨α⟩={means[-1]:.6f}")

  arr = np.asarray(means, dtype=float)
  mean_null = float(np.mean(arr))
  std_null = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
  print(f"  [{model}] fertig in {time.time() - t0:.2f}s")
  return {"mean": mean_null, "std": std_null, "means": arr}


def build_argparser() -> argparse.ArgumentParser:
  p = argparse.ArgumentParser(description="G2 α-Audit mit N1/N2-Nullmodellen (10^7)")
  p.add_argument("--n-iter", type=int, default=10, help="MC-Wiederholungen pro Nullmodell")
  p.add_argument("--seed", type=int, default=2026)
  p.add_argument("--skip-n2", action="store_true", help="Nur N1 berechnen")
  return p


def main(argv: list[str] | None = None) -> int:
  args = build_argparser().parse_args(argv)
  primes = get_primes_cached()
  print(f"[Daten] {len(primes):,} Primzahlen (Limit-Kontext: {LIMIT:,})")

  t0 = time.time()
  alphas_emp = compute_alphas_batch(primes)
  mean_emp = float(np.mean(alphas_emp))
  std_emp = float(np.std(alphas_emp))
  occ = float(np.mean(alphas_emp < 1e-5))
  print(
      f"[Empirie] ⟨α⟩_emp = {mean_emp:.6f} ± {std_emp:.6f}  |  "
      f"P(α≈0) = {occ * 100:.2f}%  |  {time.time() - t0:.2f}s"
  )

  gaps = np.diff(primes).astype(np.int64)
  rng = np.random.default_rng(args.seed)

  n1 = monte_carlo_null(primes, gaps, "N1", args.n_iter, rng)
  delta_m1 = mean_emp - float(n1["mean"])
  z_m1 = delta_m1 / float(n1["std"]) if float(n1["std"]) > 0 else 0.0

  if args.skip_n2:
    n2_mean = float("nan")
    n2_std = float("nan")
    delta_m2 = float("nan")
    z_m2 = float("nan")
  else:
    n2 = monte_carlo_null(primes, gaps, "N2", args.n_iter, rng)
    n2_mean = float(n2["mean"])
    n2_std = float(n2["std"])
    delta_m2 = mean_emp - n2_mean
    z_m2 = delta_m2 / n2_std if n2_std > 0 else 0.0

  print("\n" + "=" * 54)
  print(" ERGEBNISSE DER NULLMODELL-KONTROLLE (BIS 10^7)")
  print("=" * 54)
  print(f" Empirischer Mittelwert ⟨α⟩_emp : {mean_emp:.6f}")
  print(f" N1-Erwartungswert E_N1[⟨α⟩]    : {float(n1['mean']):.6f} ± {float(n1['std']):.6f}")
  print(f" Differenz-Metrik Δ_M1          : {delta_m1:.6f}")
  print(f" Standardisierter Z-Score Z_M1  : {z_m1:+.2f}")
  if not args.skip_n2:
    print(f" N2-Erwartungswert E_N2[⟨α⟩]    : {n2_mean:.6f} ± {n2_std:.6f}")
    print(f" Differenz-Metrik Δ_M2          : {delta_m2:.6f}")
    print(f" Standardisierter Z-Score Z_M2  : {z_m2:+.2f}")
  print(f" MC-Iterationen / Seed          : {args.n_iter} / {args.seed}")
  print("=" * 54)
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
