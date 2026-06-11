"""
Liest bm60_twin_mirror_w_scan.csv (oder per --csv) und druckt Dominanzen + Pol_11_47-Peak.

Erwartete Spalten: W, Pol_11_47, Pol_17_41, Pol_29_59, Spread (weitere Spalten ok).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path(__file__).resolve().parent / "bm60_twin_mirror_w_scan.csv",
        help="Pfad zur CSV (Default: neben diesem Skript)",
    )
    args = parser.parse_args()
    path = args.csv
    if not path.is_file():
        print(f"Datei nicht gefunden: {path}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(path)

    required = ["W", "Pol_11_47", "Pol_17_41", "Pol_29_59", "Spread"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"Fehlende Spalten: {missing}", file=sys.stderr)
        print("Vorhanden:", list(df.columns), file=sys.stderr)
        sys.exit(1)

    df = df.copy()
    df["Dominanz_11_47_vs_17_41"] = df["Pol_11_47"] - df["Pol_17_41"]
    df["Dominanz_11_47_vs_29_59"] = df["Pol_11_47"] - df["Pol_29_59"]

    print(
        df[
            [
                "W",
                "Pol_11_47",
                "Pol_17_41",
                "Pol_29_59",
                "Dominanz_11_47_vs_17_41",
                "Dominanz_11_47_vs_29_59",
                "Spread",
            ]
        ].to_string(index=False)
    )

    peak = df.loc[df["Pol_11_47"].idxmax()]
    print("\nPeak von Pol_11_47:")
    print(peak)


if __name__ == "__main__":
    main()
