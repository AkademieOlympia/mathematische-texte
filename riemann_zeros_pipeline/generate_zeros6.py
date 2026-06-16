#!/usr/bin/env python3
"""Erzeugt zeros6.npy / zeros6.npz aus einer externen Nullstellenquelle.

Die im Projekt genutzte Datei ``zeros6.npy`` enthält N=2_001_052 aufsteigende
Imaginärteile γₙ der nicht-trivialen Nullstellen von ζ(s) auf Re(s)=½
(γ₁ ≈ 14.134725 … γ_N ≈ 1.13249×10⁶).

Laut ``zeros6_meta.json`` im Repository-Root wurde die lokale Kopie aus
``zeros6.gz`` (Text, eine Zahl pro Zeile) gebaut. Diese Gzip-Datei ist
nicht versioniert (>10 MB). Typische öffentliche Tabellen stammen von
Andrew Odlyzko (Montgomery-Odlyzko-Nullstellen):

  https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros1

Einzelne γₙ können zur Plausibilisierung auch über LMFDB abgefragt werden:

  https://www.lmfdb.org/riemann/

Dieses Skript beweist oder prüft **nicht** die Riemannsche Vermutung.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import time
from pathlib import Path

import numpy as np

_NUM_RE = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
PIPELINE_DIR = Path(__file__).resolve().parent
REPO_ROOT = PIPELINE_DIR.parent
DEFAULT_INPUT = PIPELINE_DIR / "zeros6.gz"
FIRST_GAMMA_REF = 14.134725142
EXPECTED_N = 2_001_052


def load_zeros_text(path: Path) -> np.ndarray:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    values: list[float] = []
    with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = _NUM_RE.search(line)
            if m:
                values.append(float(m.group(0)))
    if not values:
        raise ValueError(f"Keine γ-Werte in {path}")
    return np.asarray(values, dtype=np.float64)


def load_zeros_binary(path: Path) -> np.ndarray:
    """Rohe IEEE-754 float64-Folge (Odlyzko-Tabellenformat)."""
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    mode = "rb"
    with opener(path, mode) as f:
        raw = f.read()
    arr = np.frombuffer(raw, dtype=np.float64)
    if arr.size == 0:
        raise ValueError(f"Leere Binärdatei: {path}")
    return np.asarray(arr, dtype=np.float64)


def load_zeros(path: Path, *, binary: bool = False) -> np.ndarray:
    if binary:
        return load_zeros_binary(path)
    try:
        arr = load_zeros_text(path)
    except ValueError:
        arr = load_zeros_binary(path)
    return arr


def validate_gamma(arr: np.ndarray) -> None:
    if arr.ndim != 1:
        raise ValueError(f"erwartet 1D-Array, erhalten shape={arr.shape}")
    if not np.all(np.isfinite(arr)):
        raise ValueError("nicht-finite Werte in γ-Folge")
    if arr[0] < 14.0 or arr[0] > 14.2:
        raise ValueError(
            f"γ₁={arr[0]:.9f} passt nicht zur ersten bekannten Nullstelle (~14.134725)"
        )
    if not np.all(arr[1:] >= arr[:-1] - 1e-9):
        raise ValueError("γ-Folge ist nicht (nahezu) monoton nicht-absteigend")


def write_outputs(
    arr: np.ndarray,
    *,
    out_npy: Path,
    out_npz: Path,
    out_meta: Path,
    source: str,
    build_seconds: float,
) -> None:
    out_npy.parent.mkdir(parents=True, exist_ok=True)
    np.save(out_npy, arr)
    np.savez_compressed(out_npz, zeros=arr)
    meta = {
        "source": source,
        "n": int(arr.size),
        "dtype": "float64",
        "min": float(arr.min()),
        "max": float(arr.max()),
        "monotone_non_decreasing": bool(np.all(arr[1:] >= arr[:-1] - 1e-15)),
        "build_seconds": build_seconds,
        "npy_file": out_npy.name,
        "npz_file": out_npz.name,
    }
    out_meta.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description="zeros6.npy aus externer γ-Quelle bauen")
    ap.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Eingabe (Text-.gz oder Binär-float64), Standard: {DEFAULT_INPUT.name}",
    )
    ap.add_argument(
        "--binary",
        action="store_true",
        help="Eingabe als rohe float64-Bytes lesen (Odlyzko-Format)",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=REPO_ROOT,
        help="Zielverzeichnis für zeros6.npy / .npz / _meta.json (Standard: Repo-Root)",
    )
    args = ap.parse_args()

    if not args.input.is_file():
        print(
            f"Eingabedatei fehlt: {args.input}\n"
            "Bitte zeros6.gz (oder Odlyzko-Tabelle) nach riemann_zeros_pipeline/ legen\n"
            "oder --input /pfad/zur/datei setzen.\n"
            "Siehe README.md (Odlyzko-URL, LMFDB)."
        )
        return 2

    t0 = time.perf_counter()
    gamma = load_zeros(args.input, binary=args.binary)
    validate_gamma(gamma)
    elapsed = time.perf_counter() - t0

    out_npy = args.out_dir / "zeros6.npy"
    out_npz = args.out_dir / "zeros6.npz"
    out_meta = args.out_dir / "zeros6_meta.json"
    write_outputs(
        gamma,
        out_npy=out_npy,
        out_npz=out_npz,
        out_meta=out_meta,
        source=str(args.input.name),
        build_seconds=elapsed,
    )

    print(f"γ geladen: n={gamma.size:,}, γ₁={gamma[0]:.9f}, γ_N={gamma[-1]:.6f}")
    if gamma.size != EXPECTED_N:
        print(f"Hinweis: erwartete N={EXPECTED_N:,}, erhalten {gamma.size:,}")
    print(f"geschrieben: {out_npy}, {out_npz}, {out_meta} ({elapsed:.2f}s)")
    print("Verifikation: python verify_zeros6.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
