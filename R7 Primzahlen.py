#!/usr/bin/env python3
"""R7 / G2-Assoziator-Audit über Primzahl-Tripel (α_n vs. N1-Gap-Markov).

Einbettung: p ↦ e_{(p mod 7)+1} in Im(𝕆).
Observable: α_n = ‖[x_n, x_{n+1}, x_{n+2}]‖ (Oktonionen-Assoziator).
Nullmodell N1: Gap-Shuffle (lokale Gap-Struktur erhalten, Reihenfolge randomisiert).

Beispiele:
  python "R7 Primzahlen.py"
  python "R7 Primzahlen.py" --limit 50000 --n-iter 200 --seed 2026
  python "R7 Primzahlen.py" --primes-file primes.txt --n-iter 100
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

# --- 1. Fano-Ebene Multiplikationstabelle (e1..e7) ---
# Fano-Tripel (123, 145, 176, 246, 257, 347, 365)
fano_table = np.zeros((8, 8), dtype=int)
triples = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
for i, j, k in triples:
  fano_table[i, j] = k
  fano_table[j, i] = -k
  fano_table[j, k] = i
  fano_table[k, j] = -i
  fano_table[k, i] = j
  fano_table[i, k] = -j
# e_i^2 = -1 → Beitrag auf dem Skalar (Index 0); Marker -1 auf Diagonale
for i in range(1, 8):
  fano_table[i, i] = -1


def oct_mul(o1: np.ndarray, o2: np.ndarray) -> np.ndarray:
  """Oktonionen-Multiplikation in ℝ⊕Im(𝕆) mit Fano-Strukturkonstanten."""
  res = np.zeros(8, dtype=float)
  res[0] += o1[0] * o2[0]
  for j in range(1, 8):
    res[j] += o1[0] * o2[j]
    res[j] += o2[0] * o1[j]
  for i in range(1, 8):
    for j in range(1, 8):
      k = int(fano_table[i, j])
      if i == j:
        # e_i^2 = -1
        res[0] -= o1[i] * o2[j]
        continue
      if k == 0:
        continue
      res[abs(k)] += o1[i] * o2[j] * (1.0 if k > 0 else -1.0)
  return res


def get_fano_vec(p: int) -> np.ndarray:
  """Einbettung: Primzahl → Einheitsvektor e_{(p mod 7)+1}."""
  vec = np.zeros(8, dtype=float)
  vec[(p % 7) + 1] = 1.0
  return vec


def get_associator_norm(p1: int, p2: int, p3: int) -> float:
  """α = ‖(xy)z − x(yz)‖ für x,y,z ∈ Im(𝕆)."""
  x = get_fano_vec(p1)
  y = get_fano_vec(p2)
  z = get_fano_vec(p3)
  res = oct_mul(oct_mul(x, y), z) - oct_mul(x, oct_mul(y, z))
  return float(np.linalg.norm(res))


def run_audit(prime_list: list[int] | np.ndarray) -> tuple[float, float, np.ndarray]:
  alphas = [
      get_associator_norm(prime_list[i], prime_list[i + 1], prime_list[i + 2])
      for i in range(len(prime_list) - 2)
  ]
  arr = np.asarray(alphas, dtype=float)
  if len(arr) == 0:
    return 0.0, 0.0, arr
  return float(np.mean(arr)), float(np.std(arr, ddof=1) if len(arr) > 1 else 0.0), arr


def sieve_primes(limit: int) -> list[int]:
  """Primzahlen in [2, limit]."""
  if limit < 2:
    return []
  is_prime = np.ones(limit + 1, dtype=bool)
  is_prime[:2] = False
  for p in range(2, int(limit**0.5) + 1):
    if is_prime[p]:
      is_prime[p * p :: p] = False
  return [int(i) for i in np.nonzero(is_prime)[0]]


def load_primes_file(path: Path) -> list[int]:
  text = path.read_text(encoding="utf-8")
  out: list[int] = []
  for part in text.replace(",", " ").split():
    out.append(int(part))
  return out


def gap_shuffle_primes(
    primes: list[int], rng: np.random.Generator
) -> list[int]:
  """N1: gleiche Gaps, zufällige Reihenfolge; Start = erste Primzahl."""
  if len(primes) < 2:
    return list(primes)
  gaps = np.diff(primes)
  rng.shuffle(gaps)
  fake = np.empty(len(primes), dtype=int)
  fake[0] = primes[0]
  fake[1:] = primes[0] + np.cumsum(gaps)
  return fake.tolist()


def markov_n1_null(
    primes: list[int],
    emp_mean: float,
    n_iter: int,
    rng: np.random.Generator,
) -> dict[str, float]:
  null_means: list[float] = []
  for _ in range(n_iter):
    fake = gap_shuffle_primes(primes, rng)
    m, _, _ = run_audit(fake)
    null_means.append(m)
  arr = np.asarray(null_means, dtype=float)
  mean_null = float(np.mean(arr)) if len(arr) else 0.0
  std_null = float(np.std(arr, ddof=1)) if len(arr) > 1 else 1e-12
  if std_null <= 0:
    std_null = 1e-12
  z = (emp_mean - mean_null) / std_null
  # zweiseitiges MC-p
  abs_emp = abs(emp_mean - mean_null)
  p_mc = (1 + int(np.count_nonzero(np.abs(arr - mean_null) >= abs_emp))) / (
      len(arr) + 1
  )
  return {
      "null_mean": mean_null,
      "null_std": std_null,
      "z_n1": float(z),
      "p_mc": float(p_mc),
  }


def build_argparser() -> argparse.ArgumentParser:
  p = argparse.ArgumentParser(description="R7 Primzahl-Assoziator-Audit (α_n vs N1)")
  p.add_argument("--limit", type=int, default=10000, help="Sieb-Obergrenze (Default)")
  p.add_argument("--primes-file", type=Path, default=None, help="Datei mit Primzahlen")
  p.add_argument("--n-iter", type=int, default=100, help="N1-Shuffle-Iterationen")
  p.add_argument("--seed", type=int, default=2026)
  p.add_argument("--smoke", action="store_true", help="Schnell: limit=500, n_iter=20")
  return p


def main(argv: list[str] | None = None) -> int:
  args = build_argparser().parse_args(argv)
  if args.smoke:
    args.limit = 500
    args.n_iter = 20

  if args.primes_file is not None:
    primes = load_primes_file(args.primes_file)
    source = f"file:{args.primes_file}"
  else:
    primes = sieve_primes(args.limit)
    source = f"sieve:limit={args.limit}"

  if len(primes) < 3:
    raise SystemExit("Zu wenige Primzahlen (≥3 nötig).")

  rng = np.random.default_rng(args.seed)
  mean_val, std_val, alphas = run_audit(primes)
  null = markov_n1_null(primes, mean_val, args.n_iter, rng)
  occupancy = float(np.mean(alphas == 0.0)) if len(alphas) else 0.0

  print("=== R7 / G2 Assoziator-Audit ===")
  print(f"  Quelle     : {source}")
  print(f"  #Primzahlen: {len(primes)}  |  #Tripel: {len(alphas)}")
  print(f"  seed       : {args.seed}  |  n_iter: {args.n_iter}")
  print(f"  ᾱ_emp      : {mean_val:.6f} ± {std_val:.6f}")
  print(f"  I(α=0)     : {occupancy:.4f}  (associating-triplet occupancy)")
  print(f"  N1 null    : {null['null_mean']:.6f} ± {null['null_std']:.6f}")
  print(f"  Z_N1       : {null['z_n1']:+.2f}")
  print(f"  p_mc       : {null['p_mc']:.4f}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
