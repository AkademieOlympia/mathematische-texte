"""
Paper-Figur: BM60 P_11/47(W) für zwei Siebgrenzen (z. B. 10^8 und 2·10^8) in einem Koordinatensystem.

Erwartete CSV: Ausgabe von bm60_p_11_47_wscan.py, Spalte Pol_11_47 (3- oder 4-Spalten-Format).

Für 420-Permutationstest siehe separates Skript; diese Figur dient der 60-Scan-Replikation in der Publikation.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _load(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "Pol_11_47" not in df.columns:
        raise ValueError(f"{path}: Spalte Pol_11_47 fehlt, Spalten: {list(df.columns)}")
    w = np.asarray(df["W"], dtype=float)
    p = np.asarray(df["Pol_11_47"], dtype=float)
    n = int(df["N"].iloc[0]) if "N" in df.columns and len(df) else 0
    return pd.DataFrame({"W": w, "Pol_11_47": p, "N": n})


def _label_n(n: int) -> str:
    if n == 100_000_000:
        return r"$N=10^8$"
    if n == 200_000_000:
        return r"$N=2\cdot 10^8$"
    return f"$N={n}$"


def main() -> None:
    p = argparse.ArgumentParser(description="P_11/47 vs W, zwei N-Kurven (Paper-Export).")
    base = Path(__file__).resolve().parent
    p.add_argument(
        "--n100m",
        type=Path,
        default=base / "bm60_p_11_47_wscan_curve_n100000000_p0.csv",
        help="CSV für N=100M (Kurven-Run: --perms 0)",
    )
    p.add_argument(
        "--n200m",
        type=Path,
        default=base / "bm60_p_11_47_wscan_curve_n200000000_p5000.csv",
        help="CSV für N=200M (laufende Kurve, Dateiname ggf. _p0)",
    )
    p.add_argument("--out-pdf", type=Path, default=base / "fig_bm60_p11_47_wscan_100m_200m.pdf")
    p.add_argument("--out-png", type=Path, default=base / "fig_bm60_p11_47_wscan_100m_200m.png")
    p.add_argument("--dpi", type=int, default=300)
    args = p.parse_args()

    d1 = _load(args.n100m)
    d2 = _load(args.n200m)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.labelsize": 11,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )
    fig, ax = plt.subplots(figsize=(5.6, 3.4), layout="constrained")
    ax.axhline(0.0, color="0.5", linewidth=0.6, zorder=0)
    ax.plot(
        d1["W"],
        d1["Pol_11_47"],
        "o-",
        color="#1f77b4",
        linewidth=1.2,
        markersize=4,
        label=_label_n(int(d1["N"].iloc[0])),
    )
    ax.plot(
        d2["W"],
        d2["Pol_11_47"],
        "s--",
        color="#d62728",
        linewidth=1.2,
        markersize=3.5,
        label=_label_n(int(d2["N"].iloc[0])),
    )

    ax.set_xlabel(r"Fensterbreite $W$")
    ax.set_ylabel(r"$P_{11/47}(W)=\Delta_{47}(W)-\Delta_{11}(W)$")
    ax.set_title("BM60: Zwillings-Scan (mod 60, Kanal 11 vs. 47)")

    ax.legend(frameon=True, loc="best")
    ax.set_xlim(350, 3150)
    for out_path, kind in ((args.out_pdf, "pdf"), (args.out_png, "png")):
        out = out_path.resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, format=kind, dpi=args.dpi if kind == "png" else None)
        print("geschrieben:", out)
    plt.close(fig)


if __name__ == "__main__":
    main()
