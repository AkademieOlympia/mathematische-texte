#!/usr/bin/env python3
"""Hadamard-Fourier-Analyse des V4-Sektor-Tensors: Δ → Δ̂ auf V̂₄^{⊗3}.

Liest gp_sector_v4_107.csv (Kanäle E→A→B, Spalte delta/mean/z),
baut den 4×4×4-Tensor und projiziert

  Δ̂_abc = Σ_{ijk} H_ai H_bj H_ck Δ_ijk

mit der kanonischen V4-Hadamard-Matrix aus dem Preprint.

Beispiele:
  python g2_audit_v4_tensor_fourier.py
  python g2_audit_v4_tensor_fourier.py --field delta --top 12 \\
      --export scripts/gp_v4_tensor_fourier_107.csv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent

SECTOR_MAP = {"E": 0, "A": 1, "B": 2, "C": 3}
CHAR_NAMES = ("χ₀", "χ₁", "χ₂", "χ₃")  # Zeilen von H

# V4-Hadamard (Preprint §4), orthonormal bis Skala 1
H = 0.5 * np.array(
    [
        [1, 1, 1, 1],
        [1, -1, 1, -1],
        [1, 1, -1, -1],
        [1, -1, -1, 1],
    ],
    dtype=float,
)


def resolve_csv(path: Path | None) -> Path:
  candidates = []
  if path is not None:
    candidates.append(path)
  candidates.extend(
      [
          SCRIPT_DIR / "gp_sector_v4_107.csv",
          Path.home()
          / "Desktop"
          / "Mathematische Texte"
          / "gp_sector_v4_107.csv",
          Path.home()
          / "Projects"
          / "Kepler-Hurrwitz"
          / "scripts"
          / "gp_sector_v4_107.csv",
      ]
  )
  for p in candidates:
    if p.exists():
      return p
  raise FileNotFoundError(
      "gp_sector_v4_107.csv nicht gefunden. Bitte zuerst g2_audit_sector_mod.py ausführen."
  )


def parse_channel(ch: str) -> tuple[int, int, int]:
  parts = str(ch).replace("->", "→").split("→")
  if len(parts) != 3:
    raise ValueError(f"Ungültiger Kanal: {ch!r}")
  return (
      SECTOR_MAP[parts[0].strip()],
      SECTOR_MAP[parts[1].strip()],
      SECTOR_MAP[parts[2].strip()],
  )


def load_tensor(csv_path: Path, field: str) -> tuple[np.ndarray, np.ndarray]:
  """Baut 4×4×4-Tensor aus CSV; gibt (T, N_counts) zurück."""
  T = np.zeros((4, 4, 4), dtype=float)
  N = np.zeros((4, 4, 4), dtype=float)
  with csv_path.open(encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f)
    if field not in (reader.fieldnames or []):
      raise KeyError(f"Spalte {field!r} fehlt in {csv_path}")
    for row in reader:
      i, j, k = parse_channel(row["channel"])
      T[i, j, k] = float(row[field])
      N[i, j, k] = float(row.get("N", 0) or 0)
  return T, N


def hadamard_tensor3(T: np.ndarray, Hmat: np.ndarray = H) -> np.ndarray:
  """Δ̂_abc = Σ_ijk H_ai H_bj H_ck Δ_ijk."""
  return np.einsum("ai,bj,ck,ijk->abc", Hmat, Hmat, Hmat, T, optimize=True)


def mode_energy(That: np.ndarray) -> np.ndarray:
  return That**2


def classify_mode(a: int, b: int, c: int) -> str:
  """Grobe Mode-Klasse für Claim-Wall-Lesbarkeit."""
  idxs = (a, b, c)
  if idxs == (0, 0, 0):
    return "trivial"
  if a == b == c and a > 0:
    return "diagonal-nontrivial"
  if len(set(idxs)) == 3 and 0 not in idxs:
    return "cyclic-nontrivial"
  if idxs.count(0) == 2:
    return "single-character"
  if idxs.count(0) == 1:
    return "two-character"
  if a == b or b == c or a == c:
    return "paired"
  return "mixed"


def build_argparser() -> argparse.ArgumentParser:
  p = argparse.ArgumentParser(description="V̂₄^{⊗3} Hadamard-Fourier des Sektor-Tensors")
  p.add_argument("--csv", type=Path, default=None)
  p.add_argument(
      "--field",
      choices=["delta", "mean", "z"],
      default="delta",
      help="Tensor-Einträge aus CSV-Spalte",
  )
  p.add_argument("--top", type=int, default=8)
  p.add_argument("--export", type=Path, default=None)
  return p


def main(argv: list[str] | None = None) -> int:
  args = build_argparser().parse_args(argv)
  csv_path = resolve_csv(args.csv)
  print(f"[Analyse] Lade Sektor-Daten aus '{csv_path}'...")
  print(f"[Feld] Tensor-Einträge = '{args.field}'")

  T, N = load_tensor(csv_path, args.field)
  filled = int(np.count_nonzero(N))
  print(f"[Tensor] 4×4×4 gefüllt: {filled}/64 Kanäle  |  ‖T‖_F = {np.linalg.norm(T):.6f}")

  That = hadamard_tensor3(T)
  delta_fro = abs(float(np.linalg.norm(That) - np.linalg.norm(T)))
  print(f"[Fourier] ‖T̂‖_F = {np.linalg.norm(That):.6f}  |  ΔF = {delta_fro:.2e}")

  E = mode_energy(That)
  E_tot = float(np.sum(E))
  flat = np.argsort(E.ravel())[::-1]

  print("\n" + "=" * 64)
  print(f" SPEKTRALE MODEN IN V̂₄^{{⊗3}}  (Top {args.top} nach Energie)")
  print("=" * 64)
  print(f"{'Rang':<5} {'Modus':<14} {'Chars':<16} {'Ampl.':>12} {'E-Anteil':>10} {'Klasse'}")
  print("-" * 64)

  cum = 0.0
  for rank, idx in enumerate(flat[: args.top], start=1):
    a, b, c = np.unravel_index(int(idx), (4, 4, 4))
    val = float(That[a, b, c])
    share = float(E[a, b, c] / E_tot) if E_tot > 0 else 0.0
    cum += share
    cls = classify_mode(a, b, c)
    chars = f"{CHAR_NAMES[a]},{CHAR_NAMES[b]},{CHAR_NAMES[c]}"
    print(
        f"{rank:<5} ({a},{b},{c}){'':<6} {chars:<16} {val:>+12.6f} {100 * share:>9.2f}% {cls}"
    )

  top1 = float(E.ravel()[flat[0]] / E_tot) if E_tot > 0 else 0.0
  top3 = float(np.sum(E.ravel()[flat[:3]]) / E_tot) if E_tot > 0 else 0.0
  top8 = float(np.sum(E.ravel()[flat[:8]]) / E_tot) if E_tot > 0 else 0.0
  triv = float(E[0, 0, 0] / E_tot) if E_tot > 0 else 0.0

  print("-" * 64)
  print(f" Energie Top-1 / Top-3 / Top-8 : {100*top1:.2f}% / {100*top3:.2f}% / {100*top8:.2f}%")
  print(f" Trivialmodus (0,0,0) Anteil   : {100*triv:.2f}%")
  print(f" Kumuliert (angezeigte Top)    : {100*cum:.2f}%")
  print("=" * 64)
  print(
      " Lesart: Konzentration auf wenige Moden + nichttriviale Klasse\n"
      " stützt die C₄^dyn–Kopplung an die V₄-Fourierbasis (Claim-Wall [B])."
  )

  if args.export is not None:
    args.export.parent.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for idx in flat:
      a, b, c = np.unravel_index(int(idx), (4, 4, 4))
      share = float(E[a, b, c] / E_tot) if E_tot > 0 else 0.0
      all_rows.append(
          {
              "a": a,
              "b": b,
              "c": c,
              "chars": f"{CHAR_NAMES[a]},{CHAR_NAMES[b]},{CHAR_NAMES[c]}",
              "amplitude": float(That[a, b, c]),
              "energy_share": share,
              "class": classify_mode(a, b, c),
          }
      )
    with args.export.open("w", encoding="utf-8", newline="") as f:
      w = csv.DictWriter(
          f,
          fieldnames=["a", "b", "c", "chars", "amplitude", "energy_share", "class"],
      )
      w.writeheader()
      w.writerows(all_rows)
    print(f"[Export] {args.export}")

  return 0


if __name__ == "__main__":
  raise SystemExit(main())
