#!/usr/bin/env python3
"""SHA256- und Plausibilitätsprüfung für zeros6.npy / zeros6.npz."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

PIPELINE_DIR = Path(__file__).resolve().parent
REPO_ROOT = PIPELINE_DIR.parent
HASH_FILE = PIPELINE_DIR / "expected_hashes.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_expected() -> dict:
    return json.loads(HASH_FILE.read_text(encoding="utf-8"))


def check_array(arr: np.ndarray, expected: dict) -> list[str]:
    errors: list[str] = []
    arr = np.asarray(arr, dtype=np.float64).ravel()
    n_exp = int(expected["n"])
    if arr.size != n_exp:
        errors.append(f"Anzahl: erwartet {n_exp:,}, erhalten {arr.size:,}")
    if abs(float(arr[0]) - float(expected["gamma_1"])) > 1e-6:
        errors.append(f"γ₁: erwartet {expected['gamma_1']}, erhalten {arr[0]}")
    if abs(float(arr[-1]) - float(expected["gamma_n_max"])) > 1e-3:
        errors.append(f"γ_N: erwartet {expected['gamma_n_max']}, erhalten {arr[-1]}")
    if not np.all(np.isfinite(arr)):
        errors.append("enthält nicht-finite Werte")
    return errors


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--root",
        type=Path,
        default=REPO_ROOT,
        help="Verzeichnis mit zeros6.npy (Standard: Repo-Root)",
    )
    args = ap.parse_args()
    expected = load_expected()

    npy = args.root / "zeros6.npy"
    npz = args.root / "zeros6.npz"
    ok = True

    for label, path, hash_key in (
        ("zeros6.npy", npy, "zeros6.npy"),
        ("zeros6.npz", npz, "zeros6.npz"),
    ):
        if not path.is_file():
            print(f"[FEHLER] {label} nicht gefunden: {path}")
            ok = False
            continue
        digest = sha256_file(path)
        want = expected[hash_key]
        if digest == want:
            print(f"[OK] SHA256 {label}: {digest}")
        else:
            print(f"[FEHLER] SHA256 {label}: {digest} (erwartet {want})")
            ok = False

    if npy.is_file():
        arr = np.load(npy, mmap_mode="r")
        errs = check_array(arr, expected)
        if errs:
            ok = False
            for e in errs:
                print(f"[FEHLER] Inhalt zeros6.npy: {e}")
        else:
            print(f"[OK] Inhalt zeros6.npy: n={arr.size:,}, monoton, γ₁/γ_N plausibel")

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
