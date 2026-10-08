#!/usr/bin/env python3
"""Hypothesis 1: Langzeitgedächtnis der Assoziator-Sequenz α_n (ACF + PSD).

Voraussetzung: primes_107.npy (neben diesem Skript oder unter scripts/).

Beispiele:
  python g2_audit_hypothesis1.py
  python g2_audit_hypothesis1.py --max-lag 2000 --psd-fit-bins 2000
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent

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


def compute_acf(x: np.ndarray, max_lag: int) -> np.ndarray:
  """Normalisierte ACF für Lags 0..max_lag-1 (FFT-basiert, O(N log N))."""
  x = np.asarray(x, dtype=float)
  x = x - np.mean(x)
  n = len(x)
  if n < 2:
    return np.zeros(max_lag)
  variance = float(np.dot(x, x) / n)
  if variance <= 0:
    return np.zeros(max_lag)

  n_fft = 1 << (2 * n - 1).bit_length()
  fx = np.fft.rfft(x, n=n_fft)
  acov = np.fft.irfft(fx * np.conjugate(fx), n=n_fft)[:max_lag].real
  return acov / (n * variance)


def psd_loglog_slope(
    alphas: np.ndarray, n_fit: int = 1000
) -> tuple[float, float, np.ndarray, np.ndarray]:
  """Welch-freier Periodogramm-PSD; Log-Log-Fit auf den niedrigsten Frequenzen."""
  x = alphas - np.mean(alphas)
  fft_vals = np.fft.rfft(x)
  psd = np.abs(fft_vals) ** 2
  freqs = np.fft.rfftfreq(len(x))
  valid = freqs > 0
  log_f = np.log(freqs[valid])
  log_p = np.log(np.maximum(psd[valid], 1e-300))
  m = min(n_fit, len(log_f))
  slope, intercept = np.polyfit(log_f[:m], log_p[:m], 1)
  return float(slope), float(intercept), freqs, psd


def build_argparser() -> argparse.ArgumentParser:
  p = argparse.ArgumentParser(description="Hypothesis 1: ACF/PSD der α_n-Sequenz")
  p.add_argument(
      "--lags",
      type=str,
      default="1,2,5,10,20,50,100,500,1000",
      help="Kommaliste der ACF-Lags",
  )
  p.add_argument("--max-lag", type=int, default=None, help="Obere Lag-Grenze (Default: max(lags)+1)")
  p.add_argument("--psd-fit-bins", type=int, default=1000, help="Anzahl niedrigster Frequenzen für Slope-Fit")
  p.add_argument(
      "--export-acf",
      type=Path,
      default=None,
      help="Optional: ACF als .npz speichern",
  )
  return p


def main(argv: list[str] | None = None) -> int:
  args = build_argparser().parse_args(argv)
  lags = [int(x.strip()) for x in args.lags.split(",") if x.strip()]
  max_lag = args.max_lag if args.max_lag is not None else max(lags) + 1

  print("[Hypothese 1] Starte Autokorrelations-Analyse der Assoziator-Sequenz...")
  primes = get_primes_cached()

  t0 = time.time()
  alphas = compute_alphas_batch(primes)
  print(f"[Daten] {len(alphas):,} Assoziator-Werte in {time.time() - t0:.2f}s berechnet.")
  print(
      f"[Empirie] ⟨α⟩ = {float(np.mean(alphas)):.6f} ± {float(np.std(alphas)):.6f}  |  "
      f"P(α≈0) = {100.0 * float(np.mean(alphas < 1e-5)):.2f}%"
  )

  t1 = time.time()
  acf_vals = compute_acf(alphas, max_lag=max_lag)
  print(f"[ACF] berechnet bis Lag {max_lag - 1} in {time.time() - t1:.2f}s")

  # Grobe Signifikanz: Weißes-Rauschen-Band ≈ ±1.96/√N
  band = 1.96 / np.sqrt(len(alphas))

  print("\n" + "=" * 58)
  print(" ERGEBNISSE: AUTOKORRELATION (LONG-RANGE MEMORY)")
  print("=" * 58)
  print(f" 95%-Weißrausch-Band ≈ ±{band:.6f}")
  for lag in lags:
    if lag >= len(acf_vals):
      print(f" Lag k = {lag:4d}  |  ACF(k) = (außerhalb max_lag)")
      continue
    val = float(acf_vals[lag])
    flag = "*" if abs(val) > band else " "
    print(f" Lag k = {lag:4d}  |  ACF(k) = {val:+.6f}  {flag}")
  print("  (* = |ACF| > 95%-Weißrausch-Band)")
  print("=" * 58)

  slope, intercept, freqs, psd = psd_loglog_slope(alphas, n_fit=args.psd_fit_bins)
  print(f" Geschätzter PSD-Anstieg (Log-Log Slope): {slope:.3f}")
  print(" (nahe 0 ≈ weißes Rauschen; um −1 ≈ 1/f / Langzeitgedächtnis)")
  print(f" Fit über die {min(args.psd_fit_bins, max(0, (freqs > 0).sum()))} niedrigsten Frequenzen")
  print("=" * 58)

  if args.export_acf is not None:
    out = np.column_stack([np.arange(len(acf_vals)), acf_vals])
    args.export_acf.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        args.export_acf,
        out,
        delimiter=",",
        header="lag,acf",
        comments="",
    )
    print(f"[Export] ACF → {args.export_acf}")

  return 0


if __name__ == "__main__":
  raise SystemExit(main())
