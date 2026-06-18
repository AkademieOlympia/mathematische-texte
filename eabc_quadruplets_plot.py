#!/usr/bin/env python3
"""Vierfeld-Diagnose-Plot: W_E, Z_E (= R_{1/2}), alpha_loc, R_beta."""

from pathlib import Path

import numpy as np
import pandas as pd

R_BETA_COLUMNS = {
    "R_1_2": 0.5,
    "R_2_3": 2 / 3,
    "R_3_4": 0.75,
    "R_9_10": 0.9,
    "R_1": 1.0,
}

LEGACY_BETA_COLUMNS = {
    "D_over_Q_half": 0.5,
    "D_over_Q_two_thirds": 2 / 3,
    "D_over_Q_three_fourths": 3 / 4,
    "D_over_Q_one": 1.0,
}


def ensure_beta_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    diff = out["diff"].to_numpy(dtype=np.float64)
    q = out["Q_total"].to_numpy(dtype=np.float64)
    for name, beta in R_BETA_COLUMNS.items():
        if name not in out.columns:
            legacy = next(
                (old for old, b in LEGACY_BETA_COLUMNS.items() if b == beta and old in out.columns),
                None,
            )
            if legacy is not None:
                out[name] = out[legacy]
            else:
                out[name] = diff / np.power(q, beta)
    if "Z_E" not in out.columns:
        out["Z_E"] = out["R_1_2"]
    if "W_E" not in out.columns:
        out["W_E"] = out["R_1"]
    return out


def alpha_loc_series(df: pd.DataFrame) -> pd.Series:
    if len(df) < 2:
        return pd.Series(np.nan, index=df.index, dtype=np.float64)

    x = np.log(df["Q_total"].to_numpy(dtype=np.float64))
    y = np.log(np.abs(df["diff"].to_numpy(dtype=np.float64)))
    dx = np.diff(x)
    dy = np.diff(y)
    valid = dx != 0
    loc = np.full(len(df), np.nan, dtype=np.float64)
    loc[1:][valid] = dy[valid] / dx[valid]
    return pd.Series(loc, index=df.index)


def make_diagnose_plot(df: pd.DataFrame, loc: pd.Series, out_path: Path) -> None:
    import matplotlib.pyplot as plt

    x = df["X"].to_numpy(dtype=np.float64)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    fig.suptitle("EABC-Quadruplets: Ebenen 0–3 / H0a–H3")

    ax = axes[0, 0]
    ax.plot(x, df["W_E"], "o-", color="C0", label=r"$W_E = R_1 = D/Q$")
    ax.axhline(0.0, color="gray", linewidth=0.8, linestyle="--")
    ax.set_ylabel(r"$W_E(X)$")
    ax.set_title("Panel 1: Ebene 0 — Orientierung (H0b / H3)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(x, df["R_1_2"], "s-", color="C1", label=r"$Z_E = R_{1/2} = D/\sqrt{Q}$")
    ax.axhline(0.0, color="gray", linewidth=0.8, linestyle="--")
    ax.set_ylabel(r"$Z_E(X) = R_{1/2}(X)$")
    ax.set_title("Panel 2: Ebene 1 — erste Testobservable (H0a / H1)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    loc_vals = loc.to_numpy(dtype=np.float64)
    ax.plot(x, loc_vals, "^-", color="C2", label=r"$\alpha_{\mathrm{loc}}$")
    ax.axhline(0.5, color="gray", linewidth=0.8, linestyle="--", label=r"$\alpha=1/2$")
    ax.set_ylabel(r"$\alpha_{\mathrm{loc}}$")
    ax.set_xlabel(r"$X$")
    ax.set_title("Panel 3: Ebene 2 — lokaler Exponent")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    labels = {
        "R_1_2": r"$R_{1/2}$",
        "R_2_3": r"$R_{2/3}$",
        "R_3_4": r"$R_{3/4}$",
        "R_9_10": r"$R_{0.9}$",
        "R_1": r"$R_1$",
    }
    for i, (col, label) in enumerate(labels.items()):
        ax.plot(x, df[col], "o-", markersize=4, label=label, color=f"C{i + 3}")
    ax.axhline(0.0, color="gray", linewidth=0.8, linestyle="--")
    ax.set_ylabel(r"$R_\beta(X)$")
    ax.set_xlabel(r"$X$")
    ax.set_title("Panel 4: Ebene 1 — Skalierungsdiagnostik")
    ax.legend(loc="best", fontsize=7)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="EABC-Quadruplets Diagnose-Plot")
    parser.add_argument("csv", nargs="?", default="eabc_quadruplets.csv")
    parser.add_argument(
        "--out",
        type=str,
        default="eabc_quadruplets_diagnose.png",
        help="Ausgabepfad (Standard: eabc_quadruplets_diagnose.png)",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    df = ensure_beta_columns(df)
    valid = df[(df["diff"] != 0) & (df["Q_total"] > 1)].copy()
    if valid.empty:
        raise SystemExit("Keine gültigen Zeilen (diff≠0, Q>1).")

    loc = alpha_loc_series(valid)
    out_path = Path(args.out)
    make_diagnose_plot(valid, loc, out_path)
    print(f"Diagnose-Plot gespeichert: {out_path.resolve()}")


if __name__ == "__main__":
    main()
