#!/usr/bin/env python3
"""C₄^dyn-Symmetrie-Gate: Zeitumkehr-Chiralität C_F und Orbit-Struktur.

Analysen:
  1. C_F = ‖Ẑ − Ẑ^{rev}‖² / (2‖Ẑ‖²) mit Ẑ^{rev}_{abc} = Ẑ_{cba}
     (äquivalent zur Fourier-Transformation des zeitumgekehrten Tensors Z_kji).
  2. C₄^dyn-Orbit: T: E→A→B→C→E auf dem Positions-Tensor Z_ijk,
     Fourier beider Seiten, Orbit-Energie und Mode-Kopplung.

Voraussetzung:
  gp_sector_v4_107.csv  und/oder  gp_v4_tensor_fourier_z_107.csv

Beispiele:
  python g2_audit_symmetry_orbits.py
  python g2_audit_symmetry_orbits.py --field z --top 8
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent

SECTOR_MAP = {"E": 0, "A": 1, "B": 2, "C": 3}
SECTOR_NAMES = ("E", "A", "B", "C")
CHAR_NAMES = ("χ₀", "χ₁", "χ₂", "χ₃")

H = 0.5 * np.array(
    [
        [1, 1, 1, 1],
        [1, -1, 1, -1],
        [1, 1, -1, -1],
        [1, -1, -1, 1],
    ],
    dtype=float,
)


def resolve_path(name: str, explicit: Path | None = None) -> Path:
  cands: list[Path] = []
  if explicit is not None:
    cands.append(explicit)
  cands.extend(
      [
          SCRIPT_DIR / name,
          Path.home() / "Desktop" / "Mathematische Texte" / name,
          Path.home() / "Projects" / "Kepler-Hurrwitz" / "scripts" / name,
      ]
  )
  for p in cands:
    if p.exists():
      return p
  raise FileNotFoundError(f"{name} nicht gefunden.")


def parse_channel(ch: str) -> tuple[int, int, int]:
  parts = str(ch).replace("->", "→").split("→")
  return (
      SECTOR_MAP[parts[0].strip()],
      SECTOR_MAP[parts[1].strip()],
      SECTOR_MAP[parts[2].strip()],
  )


def load_sector_tensor(csv_path: Path, field: str) -> np.ndarray:
  T = np.zeros((4, 4, 4), dtype=float)
  with csv_path.open(encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
      i, j, k = parse_channel(row["channel"])
      T[i, j, k] = float(row[field])
  return T


def load_fourier_tensor(csv_path: Path) -> np.ndarray:
  That = np.zeros((4, 4, 4), dtype=float)
  with csv_path.open(encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f)
    amp_key = "amplitude" if "amplitude" in (reader.fieldnames or []) else "ampl"
    for row in reader:
      That[int(row["a"]), int(row["b"]), int(row["c"])] = float(row[amp_key])
  return That


def hadamard_tensor3(T: np.ndarray) -> np.ndarray:
  return np.einsum("ai,bj,ck,ijk->abc", H, H, H, T, optimize=True)


def time_reversal_chirality(That: np.ndarray) -> dict[str, float]:
  """Zeitumkehr-Chiralität auf dem Fourier-Tensor.

  Mit Ẑ^rev_abc = Ẑ_cba gilt:
    ‖Ẑ − Ẑ^rev‖² = 2‖Ẑ‖² (1 − ρ),  ρ = ⟨Ẑ,Ẑ^rev⟩/‖Ẑ‖² ∈ [−1,1].

  C_F^{raw} = ‖Ẑ − Ẑ^rev‖² / (2‖Ẑ‖²) = 1 − ρ   ∈ [0, 2]
  C_F       = (1 − ρ) / 2                     ∈ [0, 1]  (Manuskript-Norm)
  """
  That_rev = That.transpose(2, 1, 0)
  diff_sq = float(np.sum((That - That_rev) ** 2))
  total_sq = float(np.sum(That**2))
  if total_sq <= 0:
    return {"C_F": 0.0, "C_F_raw": 0.0, "overlap": 1.0, "diff_sq": 0.0, "total_sq": 0.0}
  overlap = float(np.sum(That * That_rev) / total_sq)
  C_F_raw = diff_sq / (2.0 * total_sq)  # = 1 - overlap
  C_F = 0.5 * (1.0 - overlap)  # ∈ [0,1]
  return {
      "C_F": C_F,
      "C_F_raw": C_F_raw,
      "overlap": overlap,
      "diff_sq": diff_sq,
      "total_sq": total_sq,
  }


def apply_C4_on_sectors(T: np.ndarray, power: int = 1) -> np.ndarray:
  """(T·Z)_{ijk} = Z_{T^{-p}i, T^{-p}j, T^{-p}k} mit T: s ↦ s+1 (mod 4)."""
  p = power % 4
  out = np.zeros_like(T)
  for i in range(4):
    for j in range(4):
      for k in range(4):
        out[i, j, k] = T[(i - p) % 4, (j - p) % 4, (k - p) % 4]
  return out


def mode_rank_map(That: np.ndarray) -> dict[tuple[int, int, int], int]:
  E = That**2
  order = np.argsort(E.ravel())[::-1]
  ranks: dict[tuple[int, int, int], int] = {}
  for r, idx in enumerate(order, start=1):
    a, b, c = np.unravel_index(int(idx), (4, 4, 4))
    ranks[(a, b, c)] = r
  return ranks


def top_modes(That: np.ndarray, n: int) -> list[tuple[int, int, int]]:
  E = That**2
  order = np.argsort(E.ravel())[::-1][:n]
  return [tuple(np.unravel_index(int(i), (4, 4, 4))) for i in order]  # type: ignore[misc]


def build_argparser() -> argparse.ArgumentParser:
  p = argparse.ArgumentParser(description="C_F und C₄^dyn-Orbit-Analyse des Z-Tensors")
  p.add_argument("--sector-csv", type=Path, default=None)
  p.add_argument("--fourier-csv", type=Path, default=None)
  p.add_argument("--field", choices=["z", "delta", "mean"], default="z")
  p.add_argument("--top", type=int, default=8)
  p.add_argument(
      "--export",
      type=Path,
      default=None,
      help="Orbit-Tabelle als CSV",
  )
  return p


def main(argv: list[str] | None = None) -> int:
  args = build_argparser().parse_args(argv)

  sector_csv = resolve_path("gp_sector_v4_107.csv", args.sector_csv)
  print(f"[Daten] Sektor-Tensor aus '{sector_csv}' (Feld={args.field})")
  Z = load_sector_tensor(sector_csv, args.field)
  That = hadamard_tensor3(Z)
  print(f"[Fourier] ‖Ẑ‖_F = {np.linalg.norm(That):.6f}")

  try:
    fourier_csv = resolve_path("gp_v4_tensor_fourier_z_107.csv", args.fourier_csv)
    That_file = load_fourier_tensor(fourier_csv)
    dF = float(np.linalg.norm(That - That_file))
    print(f"[Check] Abstand zu '{fourier_csv.name}': ΔF={dF:.2e}")
  except FileNotFoundError:
    print("[Check] Fourier-CSV optional nicht vorhanden — übersprungen.")

  cf = time_reversal_chirality(That)
  print("\n" + "=" * 64)
  print(" ZEITUMKEHR-CHIRALITÄT (C_F)")
  print("=" * 64)
  print(f" C_F (normiert ∈[0,1]) = {cf['C_F']:.6f}")
  print(f" C_F_raw (=1−ρ)         = {cf['C_F_raw']:.6f}")
  print(f" Overlap ρ=⟨Ẑ,Ẑ^rev⟩/‖Ẑ‖² = {cf['overlap']:+.6f}")
  print(" (C_F=0 ⇒ zeitumkehr-symmetrisch; C_F→1 ⇒ maximal antichiral)")
  if cf["C_F"] > 0.05:
    print(" → Signifikante Zeitrichtungs-Asymmetrie (Claim-Wall [B]).")
  else:
    print(" → Schwache/keine Zeitumkehr-Asymmetrie.")

  print("\n" + "=" * 64)
  print(" C₄^dyn-ORBIT (T: E→A→B→C→E auf Sektoren)")
  print("=" * 64)

  energies = []
  overlaps = []
  That_orbit = []
  for p in range(4):
    Zp = apply_C4_on_sectors(Z, power=p)
    Thp = hadamard_tensor3(Zp)
    That_orbit.append(Thp)
    energies.append(float(np.sum(Thp**2)))
    overlaps.append(float(np.sum(That * Thp) / np.sum(That**2)))
    print(
        f" T^{p}: ‖Ẑ‖²={energies[-1]:.4f}  |  "
        f"Overlap mit T^0: {overlaps[-1]:+.4f}"
    )

  e_std = float(np.std(energies))
  print(f" Energie-Isometrie unter ⟨T⟩: std(‖Ẑ‖²)={e_std:.2e}")

  ranks0 = mode_rank_map(That)
  tops = top_modes(That, args.top)

  print("\n[Orbit] Top-Moden und Amplituden unter T^p (gleicher Fourier-Index):")
  print(f"{'Modus':<10} {'Chars':<16} " + " ".join(f"{'T^'+str(p):>10}" for p in range(4)) + "  Orbit?")
  print("-" * 78)

  export_rows: list[dict[str, object]] = []
  closed = 0
  for m in tops:
    a, b, c = m
    chars = f"{CHAR_NAMES[a]},{CHAR_NAMES[b]},{CHAR_NAMES[c]}"
    amps = [float(That_orbit[p][a, b, c]) for p in range(4)]
    ranks = []
    for p in range(4):
      rp = mode_rank_map(That_orbit[p]).get(m, 64)
      ranks.append(rp)
    stable = max(ranks) <= max(2 * args.top, 16)
    if stable:
      closed += 1
    flag = "YES" if stable else "no"
    print(
        f"({a},{b},{c})  {chars:<16} "
        + " ".join(f"{v:>+10.2f}" for v in amps)
        + f"  {flag}"
    )
    export_rows.append(
        {
            "a": a,
            "b": b,
            "c": c,
            "chars": chars,
            "amp_T0": amps[0],
            "amp_T1": amps[1],
            "amp_T2": amps[2],
            "amp_T3": amps[3],
            "rank_T0": ranks[0],
            "rank_T1": ranks[1],
            "rank_T2": ranks[2],
            "rank_T3": ranks[3],
            "orbit_stable": stable,
        }
    )

  print("\n[Heuristik] Index-Shift T:(a,b,c)↦(a+1,b+1,c+1) auf Top-Moden:")
  for m in tops[:4]:
    mt = ((m[0] + 1) % 4, (m[1] + 1) % 4, (m[2] + 1) % 4)
    print(
        f"  {m} ({That[m]:+.2f}, Rang {ranks0[m]}) → "
        f"{mt} ({That[mt]:+.2f}, Rang {ranks0[mt]})"
    )

  print("\n" + "=" * 64)
  print(" FAZIT (Claim-Wall [B])")
  print("=" * 64)
  print(f" C_F = {cf['C_F']:.4f}  |  Orbit-stabile Top-Moden: {closed}/{len(tops)}")
  if cf["C_F"] > 0.05 and closed >= max(1, len(tops) // 3):
    print(" → C₄^dyn-Symmetrie-Gate: Chiralität + Orbit-Kopplung empirisch gestützt.")
  elif cf["C_F"] > 0.05:
    print(" → Chiralität ja; Orbit-Stabilität nur partiell — vorsichtig formulieren.")
  else:
    print(" → Kein klares Symmetrie-Gate aus C_F/Orbit.")
  print("=" * 64)

  if args.export is not None:
    args.export.parent.mkdir(parents=True, exist_ok=True)
    with args.export.open("w", encoding="utf-8", newline="") as f:
      w = csv.DictWriter(f, fieldnames=list(export_rows[0].keys()))
      w.writeheader()
      w.writerows(export_rows)
    print(f"[Export] {args.export}")

  return 0


if __name__ == "__main__":
  raise SystemExit(main())
