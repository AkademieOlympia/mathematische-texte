#!/usr/bin/env python3
"""
Experiment: Lean-Kopplung des EABC-Bernoulli-Sensors Φ(n)=V_n.

EABC-Zerlegungsprinzip: N = (N_glatt, N_EABC). Dieses Experiment prüft die N_EABC-Schicht:
ob von-Staudt-Clausen-Primsignaturen P_n dieselbe mod-12-EABC-Klassifikation
tragen wie die formale Lean-Schicht (`EABC.lean` → `eabc_from_lean.py`) und wendet
Lean-observable Strukturen (T-Rotation, Chiralität, Vierlinge, Bernoulli-Uhr) auf P_n an.

Hinweis: „LEA-M“ im Gespräch = Fehlhörung von „Lean“; dieses Modul heißt bewusst *lean*.

Ausführung:
    python3 collatz_eabc_bernoulli_lean_test.py
    python3 collatz_eabc_bernoulli_lean_test.py --max-n 200 --output collatz_eabc_bernoulli_lean.json
"""

from __future__ import annotations

import argparse
import json
import math
import random
import statistics
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from collatz_eabc_bernoulli_sensor import (
    EabcVector,
    bernoulli_row,
    chi_eabc,
    i_chir,
    prime_sig,
    sigma_eabc,
    v_bernoulli,
    _sieve_primes,
)
from eabc_from_lean import (
    Chirality,
    EClass,
    chirality_order,
    class_of,
    is_prime_quadruplet,
    q,
    residue,
    t,
    t4,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_bernoulli_lean.json"

LEAN_SOURCES = [
    "EABC.lean (classOf, residue, T, T4, Chirality, Q, IsPrimeQuadruplet)",
    "eabc_from_lean.py (Python-Spiegel)",
    "CollatzEabc.Mod12Matrix.lean (EabcIndex, eabcResidue)",
    "CollatzEabc.BernoulliClock.lean (BernoulliCell, radialAt)",
    "CollatzEabc.PrefProjection.lean (theta, rho — Wortebene)",
]

CHIRALITY_WORDS = ("ABCE", "CEAB")
IOTA_E = {EClass.E: 0, EClass.A: 1, EClass.B: 2, EClass.C: 3}


@dataclass(frozen=True, slots=True)
class LeanClassCounts:
    e: int
    a: int
    b: int
    c: int

    def as_vector(self) -> EabcVector:
        return EabcVector(self.e, self.a, self.b, self.c)

    def as_tuple(self) -> tuple[int, int, int, int]:
        return (self.e, self.a, self.b, self.c)


def lean_class_counts(sig: list[int]) -> LeanClassCounts:
    """Unabhängige Zählung über eabc_from_lean.class_of (Lean classOf)."""
    counts = {EClass.E: 0, EClass.A: 0, EClass.B: 0, EClass.C: 0}
    for p in sig:
        cls = class_of(p)
        if cls is not None:
            counts[cls] += 1
    return LeanClassCounts(
        e=counts[EClass.E],
        a=counts[EClass.A],
        b=counts[EClass.B],
        c=counts[EClass.C],
    )


def prime_sig_word(sig: list[int], *, min_prime: int = 5) -> str:
  """EABC-Buchstabenwort aus P_n (aufsteigend, nur EABC-Primzahlen ≥ min_prime)."""
  letters: list[str] = []
  for p in sorted(sig):
    if p < min_prime:
      continue
    cls = class_of(p)
    if cls is not None:
      letters.append(cls.value)
  return "".join(letters)


def theta_word(w: str) -> int:
  return sum(IOTA_E[EClass(c)] for c in w)


def theta4_word(w: str) -> int:
  return theta_word(w) % 4


def rho_word(w: str) -> float:
  return 2.0 ** (-len(w)) if w else 0.0


def bernoulli_cell_indices(n: int) -> tuple[int, int, int]:
  """Lean BernoulliCell m=n: Indizes (2m-2, 2m, 2m+2)."""
  m = n
  return (2 * m - 2, 2 * m, 2 * m + 2)


def radial_at(m: int) -> float:
  """Lean radialAt m = 2^{-m}."""
  return 2.0 ** (-m)


def contains_chirality_subword(w: str) -> dict[str, bool]:
  return {cw: cw in w for cw in CHIRALITY_WORDS}


def quadruplet_witnesses(sig: list[int]) -> list[dict[str, Any]]:
  """Primvierlinge p∈P_n mit Lean IsPrimeQuadruplet."""
  out: list[dict[str, Any]] = []
  for p in sig:
    if not is_prime_quadruplet(p):
      continue
    classes = [class_of(x) for x in q(p)]
    out.append(
      {
        "p": p,
        "quadruplet": q(p),
        "center": p + 4,
        "classes": [c.value if c else None for c in classes],
        "word": "".join(c.value for c in classes if c is not None),
      }
    )
  return out


def t_rotate_word(w: str) -> str:
  return "".join(t(EClass(c)).value for c in w)


def classification_consistency(sig: list[int], sensor_v: EabcVector) -> dict[str, Any]:
  """Sensor V_n vs. unabhängige Lean-Zählung."""
  lean = lean_class_counts(sig)
  sensor_t = sensor_v.as_tuple()
  lean_t = lean.as_tuple()
  residue_ok = all(
    class_of(p) is not None and residue(class_of(p)) == p % 12
    for p in sig
    if class_of(p) is not None
  )
  t4_ok = all(t4(cls) is cls for cls in EClass)
  return {
    "sensor_V": sensor_t,
    "lean_V": lean_t,
    "match": sensor_t == lean_t,
    "residue_roundtrip_ok": residue_ok,
    "T4_identity_ok": t4_ok,
  }


def null_model_sigma_chi(
  sig: list[int],
  *,
  trials: int,
    seed: int,
) -> dict[str, float]:
  """Zufalls-Nullmodell: gleiche Anzahl EABC-Primzahlen, zufällige Klassenverteilung."""
  eabc_primes = [p for p in sig if class_of(p) is not None]
  k = len(eabc_primes)
  if k == 0:
    return {"sigma_mean": 0.0, "chi_mean": 0.0, "sigma_std": 0.0, "chi_std": 0.0}

  rng = random.Random(seed)
  sigmas: list[int] = []
  chis: list[int] = []
  for _ in range(trials):
    counts = {EClass.E: 0, EClass.A: 0, EClass.B: 0, EClass.C: 0}
    for _p in eabc_primes:
      counts[rng.choice(list(EClass))] += 1
    vec = EabcVector(counts[EClass.E], counts[EClass.A], counts[EClass.B], counts[EClass.C])
    sigmas.append(sigma_eabc(vec))
    chis.append(chi_eabc(vec))

  return {
    "sigma_mean": statistics.mean(sigmas),
    "chi_mean": statistics.mean(chis),
    "sigma_std": statistics.pstdev(sigmas) if len(sigmas) > 1 else 0.0,
    "chi_std": statistics.pstdev(chis) if len(chis) > 1 else 0.0,
  }


def analyze_row(
  n: int,
  sig: list[int],
  sensor_v: EabcVector,
  *,
  null_trials: int,
  null_seed: int,
) -> dict[str, Any]:
  word = prime_sig_word(sig)
  triplet_idx = bernoulli_cell_indices(n)
  cons = classification_consistency(sig, sensor_v)
  null = null_model_sigma_chi(sig, trials=null_trials, seed=null_seed + n)
  sigma = sigma_eabc(sensor_v)
  chi = chi_eabc(sensor_v)
  z_sigma = (sigma - null["sigma_mean"]) / null["sigma_std"] if null["sigma_std"] > 0 else 0.0
  z_chi = (chi - null["chi_mean"]) / null["chi_std"] if null["chi_std"] > 0 else 0.0

  return {
    "n": n,
    "two_n": 2 * n,
    "prime_sig_size": len(sig),
    "consistency": cons,
    "word": word,
    "word_len": len(word),
    "theta": theta_word(word),
    "theta4": theta4_word(word),
    "rho": rho_word(word),
    "T_word": t_rotate_word(word),
    "chirality_subwords": contains_chirality_subword(word),
    "quadruplet_witnesses": quadruplet_witnesses(sig),
    "bernoulli_cell": {
      "m": n,
      "triplet_indices": triplet_idx,
      "radial_at_m": radial_at(n),
    },
    "sigma": sigma,
    "chi": chi,
    "i_chir": i_chir(sensor_v),
    "null_model": null,
    "z_sigma_vs_null": z_sigma,
    "z_chi_vs_null": z_chi,
  }


def summarize(rows: list[dict[str, Any]], *, max_n: int) -> dict[str, Any]:
  n_rows = len(rows)
  match_all = all(r["consistency"]["match"] for r in rows)
  residue_all = all(r["consistency"]["residue_roundtrip_ok"] for r in rows)
  t4_all = all(r["consistency"]["T4_identity_ok"] for r in rows)

  words = [r["word"] for r in rows if r["word"]]
  abce_hits = sum(1 for r in rows if r["chirality_subwords"]["ABCE"])
  ceab_hits = sum(1 for r in rows if r["chirality_subwords"]["CEAB"])
  quad_total = sum(len(r["quadruplet_witnesses"]) for r in rows)

  z_sigmas = [r["z_sigma_vs_null"] for r in rows]
  z_chis = [r["z_chi_vs_null"] for r in rows]
  large_z_sigma = sum(1 for z in z_sigmas if abs(z) > 2.0)
  large_z_chi = sum(1 for z in z_chis if abs(z) > 2.0)

  return {
    "classification_match_all": match_all,
    "residue_roundtrip_all": residue_all,
    "T4_identity_all": t4_all,
    "rows_with_eabc_word": len(words),
    "mean_word_len": statistics.mean([len(w) for w in words]) if words else 0.0,
    "chirality_subword_hits": {"ABCE": abce_hits, "CEAB": ceab_hits},
    "quadruplet_witness_total": quad_total,
    "z_sigma_large_count": large_z_sigma,
    "z_chi_large_count": large_z_chi,
    "z_sigma_mean_abs": statistics.mean(abs(z) for z in z_sigmas) if z_sigmas else 0.0,
    "z_chi_mean_abs": statistics.mean(abs(z) for z in z_chis) if z_chis else 0.0,
    "max_n": max_n,
    "n_rows": n_rows,
  }


def run_lean_coupling(
  max_n: int = 100,
  *,
  null_trials: int = 500,
  null_seed: int = 20260617,
) -> dict[str, Any]:
  primes = _sieve_primes(2 * max_n + 1)
  rows: list[dict[str, Any]] = []

  for n in range(1, max_n + 1):
    row = bernoulli_row(n, primes)
    rows.append(
      analyze_row(
        n,
        row.prime_sig,
        row.v,
        null_trials=null_trials,
        null_seed=null_seed,
      )
    )

  summary = summarize(rows, max_n=max_n)

  verdict_parts: list[str] = []
  if summary["classification_match_all"] and summary["residue_roundtrip_all"]:
    verdict_parts.append(
      "Lean-Klassifikation und Sensor V_n stimmen punktweise überein (Konsistenzcheck bestanden)."
    )
  else:
    verdict_parts.append("Klassifikationsabweichung — Bridge oder Sensor prüfen.")

  if summary["z_sigma_large_count"] + summary["z_chi_large_count"] > max_n * 0.05:
    verdict_parts.append(
      "Mehr als 5 % der Stufen mit |z|>2 ggü. Zufalls-Nullmodell — explorative Strukturhinweise."
    )
  else:
    verdict_parts.append(
      "Chirale Bilanz σ/χ liegt überwiegend im Zufallsbereich des Nullmodells "
      "(keine starke Lean-Kopplungs-Resonanz)."
    )

  return {
    "framework": "EABC",
    "experiment": "Lean-Kopplung",
    "note": "LEA-M war Fehlhörung von Lean; Test gegen formale EABC-Schicht.",
    "epistemic_label": "Experiment",
    "lean_sources": LEAN_SOURCES,
    "python_bridge": "eabc_from_lean.py",
    "sensor_module": "collatz_eabc_bernoulli_sensor.py",
    "max_n": max_n,
    "null_model": {
      "trials_per_n": null_trials,
      "seed": null_seed,
      "description": "Zufällige EABC-Verteilung bei fester |P_n∩EABC|",
    },
    "summary": summary,
    "verdict": " ".join(verdict_parts),
    "rows": rows,
    "chirality_reference": {
      chir.value: "".join(c.value for c in chirality_order(chir)) for chir in Chirality
    },
  }


def main() -> None:
  parser = argparse.ArgumentParser(description="EABC-Bernoulli-Sensor vs. Lean EABC")
  parser.add_argument("--max-n", type=int, default=100)
  parser.add_argument("--null-trials", type=int, default=500)
  parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
  args = parser.parse_args()
  if args.max_n < 1:
    raise SystemExit("--max-n muss ≥ 1 sein")

  report = run_lean_coupling(args.max_n, null_trials=args.null_trials)
  args.output.write_text(
    json.dumps(report, indent=2, ensure_ascii=False) + "\n",
    encoding="utf-8",
  )
  s = report["summary"]
  print(
    f"Geschrieben: {args.output} | match_all={s['classification_match_all']} | "
    f"quad_witnesses={s['quadruplet_witness_total']} | "
    f"|z|>2: σ={s['z_sigma_large_count']}, χ={s['z_chi_large_count']}"
  )


if __name__ == "__main__":
  main()
