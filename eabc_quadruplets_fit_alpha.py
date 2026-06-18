#!/usr/bin/env python3
"""Log-log-Steigung alpha aus eabc_quadruplets.csv: |D(X)| ~ Q(X)^alpha."""

import argparse
import sys

import numpy as np
import pandas as pd


def load_valid(df: pd.DataFrame) -> pd.DataFrame:
    return df[(df["diff"] != 0) & (df["Q_total"] > 1)].copy()


def alpha_eff_series(df: pd.DataFrame) -> pd.Series:
    if "alpha_eff" in df.columns and df["alpha_eff"].notna().any():
        return df["alpha_eff"]
    log_q = np.log(df["Q_total"].to_numpy(dtype=np.float64))
    log_d = np.log(np.abs(df["diff"].to_numpy(dtype=np.float64)))
    return pd.Series(log_d / log_q, index=df.index)


def alpha_loc_series(df: pd.DataFrame) -> pd.Series:
    if len(df) < 2:
        return pd.Series(dtype=np.float64)

    x = np.log(df["Q_total"].to_numpy(dtype=np.float64))
    y = np.log(np.abs(df["diff"].to_numpy(dtype=np.float64)))
    dx = np.diff(x)
    dy = np.diff(y)
    valid = dx != 0
    loc = np.full(len(df), np.nan, dtype=np.float64)
    loc[1:][valid] = dy[valid] / dx[valid]
    return pd.Series(loc, index=df.index)


def fit_alpha(df: pd.DataFrame) -> tuple[float, float]:
    if len(df) < 2:
        raise SystemExit("Zu wenige Zeilen für Regression (diff≠0, Q>1).")

    x = np.log(df["Q_total"].to_numpy(dtype=np.float64))
    y = np.log(np.abs(df["diff"].to_numpy(dtype=np.float64)))
    alpha, beta = np.polyfit(x, y, 1)
    return float(alpha), float(beta)


def interpret(alpha: float, alpha_loc_max: float | None = None) -> str:
    hints = []
    if alpha_loc_max is not None and alpha_loc_max > 0.5:
        hints.append(f"alpha_loc max={alpha_loc_max:.4f} → H1 (empirischer Bias)")
    if abs(alpha - 0.5) < 0.15:
        hints.append(f"global alpha≈1/2 → H0 (Rausch)")
    elif 0.5 < alpha < 1.0:
        hints.append(f"global alpha={alpha:.4f} → H2-Kandidat (asymptotischer Bias, W_E→0)")
    elif alpha >= 0.85:
        hints.append(f"global alpha≈1 → H3-Kandidat (Holonomie, W_E→Φ_E≠0)")
    else:
        hints.append("Zwischenbereich — genauere Grenze bei größerem X")
    return "; ".join(hints)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", nargs="?", default="eabc_quadruplets.csv")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    valid = load_valid(df)
    if valid.empty:
        raise SystemExit("Keine gültigen Zeilen (diff≠0, Q>1).")

    eff = alpha_eff_series(valid)
    loc = alpha_loc_series(valid)

    print("alpha_eff pro Checkpoint:")
    for _, row in valid.iterrows():
        x = int(row["X"])
        ae = eff.loc[row.name]
        al = loc.loc[row.name]
        al_str = f"{al:.4f}" if not np.isnan(al) else "—"
        print(f"  X={x:>12}  alpha_eff={ae:.4f}  alpha_loc={al_str}")

    alpha_loc_vals = loc.dropna()
    alpha_loc_max = float(alpha_loc_vals.max()) if not alpha_loc_vals.empty else None
    if not alpha_loc_vals.empty:
        print("\nalpha_loc (lokale Steigung zwischen Checkpoints):")
        for _, row in valid.iloc[1:].iterrows():
            al = loc.loc[row.name]
            if not np.isnan(al):
                print(f"  bis X={int(row['X']):>12}: {al:.4f}")

    alpha, beta = fit_alpha(valid)
    print(f"\nglobal alpha (polyfit) = {alpha:.4f}")
    print(f"beta                   = {beta:.4f}")
    print(f"→ {interpret(alpha, alpha_loc_max)}")


if __name__ == "__main__":
    main()
