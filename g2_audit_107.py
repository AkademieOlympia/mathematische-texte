#!/usr/bin/env python3
"""G2-Assoziator-Audit bis 10^7 mit lokalem Primzahl-Cache.

Cache: primes_107.npy (einmaliges Sieb, danach Sofortladen).
Observable: α_n = ‖[x_n,x_{n+1},x_{n+2}]‖ mit gap-gewichteter Fano-Einbettung.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
CACHE_FILE = SCRIPT_DIR / "primes_107.npy"
LIMIT = 10_000_000  # 10^7 Obergrenze


def get_primes_cached(
    limit: int = LIMIT, cache_filename: Path | str = CACHE_FILE
) -> np.ndarray:
  """Lädt Primzahlen aus dem Cache oder generiert sie einmalig per NumPy-Sieb."""
  cache_path = Path(cache_filename)
  if cache_path.exists():
    print(f"[Cache] Lade Primzahlen aus '{cache_path}'...")
    start = time.time()
    primes = np.load(cache_path)
    print(f"[Cache] {len(primes):,} Primzahlen in {time.time() - start:.3f} s geladen.")
    return primes

  print(f"[Sieb] Generiere Primzahlen bis {limit:,} (Sieb von Eratosthenes)...")
  start = time.time()
  sieve = np.ones(limit + 1, dtype=bool)
  sieve[0] = sieve[1] = False
  for p in range(2, int(np.sqrt(limit)) + 1):
    if sieve[p]:
      sieve[p * p :: p] = False
  primes = np.nonzero(sieve)[0].astype(np.int64)
  print(f"[Sieb] {len(primes):,} Primzahlen in {time.time() - start:.2f} s ermittelt.")

  cache_path.parent.mkdir(parents=True, exist_ok=True)
  np.save(cache_path, primes)
  print(f"[Cache] Datei erfolgreich unter '{cache_path}' gespeichert.")
  return primes


# --- Fano-Ebene Multiplikationstabelle (e1..e7) ---
fano_table = np.zeros((8, 8), dtype=int)
triples = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
for i, j, k in triples:
  fano_table[i, j] = k
  fano_table[j, i] = -k
  fano_table[j, k] = i
  fano_table[k, j] = -i
  fano_table[k, i] = j
  fano_table[i, k] = -j
for i in range(1, 8):
  fano_table[i, i] = -1  # Marker; Diagonale wird in oct_mul_batch über -⟨·,·⟩ behandelt


def oct_mul_batch(o1: np.ndarray, o2: np.ndarray) -> np.ndarray:
  """Vektorisierte Oktonionen-Multiplikation für Arrays der Form (N, 8)."""
  res = np.zeros((len(o1), 8), dtype=float)
  # Skalar: a0 b0 − ⟨a⃗, b⃗⟩  (enthält e_i² = −1)
  res[:, 0] = o1[:, 0] * o2[:, 0] - np.sum(o1[:, 1:] * o2[:, 1:], axis=1)
  # Vektor: a0 b⃗ + b0 a⃗
  res[:, 1:] = o1[:, 0:1] * o2[:, 1:] + o2[:, 0:1] * o1[:, 1:]

  # Kreuzprodukt-Anteil über Fano (nur i ≠ j; Diagonale schon im Skalar)
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
  """p mod 7 → Index 1..7 (Rest 0 → 7)."""
  ind = (p_arr % 7).astype(np.int64)
  ind[ind == 0] = 7
  return ind


def build_gap_weighted_vecs(
    primes: np.ndarray, g_arr: np.ndarray, ind: np.ndarray
) -> np.ndarray:
  """Einbettung: x = g · e_{p mod 7} ∈ Im(𝕆)."""
  vecs = np.zeros((len(ind), 8), dtype=float)
  for k in range(1, 8):
    mask = ind == k
    vecs[mask, k] = g_arr[mask]
  return vecs


def run_full_audit() -> None:
  primes = get_primes_cached()

  print(f"Starte Assoziator-Audit für {len(primes):,} Primzahlen...")
  start_time = time.time()

  # Gaps g_n = (p_{n+1} − p_n) / log(p_n); für (x,y,z) brauchen wir g_n, g_{n+1}, g_{n+2}
  # → Tripel-Index n = 0 .. len(primes)-4  ⇒ N_triples = len(primes) - 3
  g_all = np.diff(primes).astype(float) / np.log(primes[:-1].astype(float))
  N_triples = len(primes) - 3
  if N_triples < 1:
    raise SystemExit("Zu wenige Primzahlen für Tripel mit drei Gaps.")

  g_x = g_all[:N_triples]
  g_y = g_all[1 : N_triples + 1]
  g_z = g_all[2 : N_triples + 2]

  ind_x = fano_index(primes[:N_triples])
  ind_y = fano_index(primes[1 : N_triples + 1])
  ind_z = fano_index(primes[2 : N_triples + 2])

  print(f"Analysiere {N_triples:,} Primzahl-Tripel...")

  X_vec = build_gap_weighted_vecs(primes, g_x, ind_x)
  Y_vec = build_gap_weighted_vecs(primes, g_y, ind_y)
  Z_vec = build_gap_weighted_vecs(primes, g_z, ind_z)

  # Assoziator: [X,Y,Z] = (XY)Z − X(YZ)
  XY = oct_mul_batch(X_vec, Y_vec)
  YZ = oct_mul_batch(Y_vec, Z_vec)
  Ass = oct_mul_batch(XY, Z_vec) - oct_mul_batch(X_vec, YZ)
  alphas = np.linalg.norm(Ass, axis=1)

  mean_alpha = float(np.mean(alphas))
  std_alpha = float(np.std(alphas))
  zero_occupancy = float(np.mean(alphas < 1e-5))

  print("\n" + "=" * 40)
  print(f" ERGEBNISSE (Limit: {LIMIT:,})")
  print("=" * 40)
  print(f" Ausgewertete Tripel : {N_triples:,}")
  print(f" Mittelwert alpha_n  : {mean_alpha:.6f} ± {std_alpha:.6f}")
  print(f" Assoziative Belegung: {zero_occupancy * 100:.2f}% (P(alpha≈0))")
  print(f" Rechenzeit          : {time.time() - start_time:.2f} Sekunden")
  print("=" * 40)


if __name__ == "__main__":
  run_full_audit()
