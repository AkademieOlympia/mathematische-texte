#!/usr/bin/env python3
"""Log-log-Steigung aus eabc_quadruplets.csv; H0a/H0b–H3-Einordnung und Diagnose-Plot.

Fünf Ebenen (vgl. collatz_eabc_zirkulationshypothese.md §4.2):
  Ebene 0 — Geometrie: G_E, gamma^+=ABCEA, gamma^-=CEABC
  Ebene 1 — Zirkulationsfehler: D_E = N_+ - N_- (primäre Observable)
  Ebene 2 — Skalierung: R_beta, alpha_loc (Numerik primär), alpha_E_hat
  Ebene 3 — Orientierung: W_E = D_E / Q
  Ebene 4 — Holonomie: Phi_E = lim W_E (Hypothese am Ende)

Reihenfolge Ausgabe: Ebene 0 → 1 → 2 → 3 → 4.
Vorwärtskette: D_E → alpha_E → W_E → Phi_E (keine Umkehrungen).
"""

import argparse
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


def load_valid(df: pd.DataFrame) -> pd.DataFrame:
    return df[(df["diff"] != 0) & (df["Q_total"] > 1)].copy()


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


def alpha_eff_series(df: pd.DataFrame) -> pd.Series:
    if "alpha_eff" in df.columns and df["alpha_eff"].notna().any():
        return df["alpha_eff"]
    log_q = np.log(df["Q_total"].to_numpy(dtype=np.float64))
    log_d = np.log(np.abs(df["diff"].to_numpy(dtype=np.float64)))
    return pd.Series(log_d / log_q, index=df.index)


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


def fit_alpha(df: pd.DataFrame) -> tuple[float, float]:
    if len(df) < 2:
        raise SystemExit("Zu wenige Zeilen für Regression (diff≠0, Q>1).")

    x = np.log(df["Q_total"].to_numpy(dtype=np.float64))
    y = np.log(np.abs(df["diff"].to_numpy(dtype=np.float64)))
    alpha, beta = np.polyfit(x, y, 1)
    return float(alpha), float(beta)


def coefficient_of_variation(values: np.ndarray) -> float:
    if values.size < 2:
        return float("inf")
    mean = float(np.mean(values))
    if abs(mean) < 1e-15:
        return float(np.std(values))
    return float(np.std(values) / abs(mean))


def r_beta_log_slope(df: pd.DataFrame, col: str) -> float:
    """Steigung von log|R_beta| vs log Q (≈ alpha_E - beta)."""
    if len(df) < 2:
        return float("nan")
    q = df["Q_total"].to_numpy(dtype=np.float64)
    r = np.abs(df[col].to_numpy(dtype=np.float64))
    r = np.maximum(r, 1e-15)
    slope, _ = np.polyfit(np.log(q), np.log(r), 1)
    return float(slope)


def estimate_alpha_E(df: pd.DataFrame) -> tuple[float, float, dict[float, float]]:
    """Heuristische alpha_E-Schaetzung aus R_beta-Plateaus (Experiment, kein Theorem).

    Fuer jedes beta: alpha_E-Kandidat = beta + Steigung(log|R_beta| vs log Q).
    Waehlt den Kandidaten bei minimalem |Steigung| (Plateau-Lesart).
    """
    slopes: dict[float, float] = {}
    candidates: dict[float, float] = {}
    for col, beta in R_BETA_COLUMNS.items():
        slope = r_beta_log_slope(df, col)
        slopes[beta] = slope
        candidates[beta] = beta + slope

    best_beta = min(slopes, key=lambda b: abs(slopes[b]))
    alpha_E_hat = candidates[best_beta]
    return alpha_E_hat, best_beta, slopes


def diagnose_hypothesis(alpha_E_hat: float, alpha_polyfit: float) -> str:
    if alpha_E_hat <= 0.55:
        band = "H0a (alpha_E <= 1/2, Wurzelrauschen)"
    elif alpha_E_hat < 0.85:
        band = f"H1/H2 (1/2 < alpha_E < 1, sublinearer Bias)"
    else:
        band = "H3-Kandidat (alpha_E ≈ 1, Holonomie-Grenzfall)"
    return (
        f"alpha_E_hat={alpha_E_hat:.3f} (heuristisch, Experiment) → {band}; "
        f"polyfit alpha={alpha_polyfit:.3f}"
    )


def diagnose_orientation(df: pd.DataFrame) -> str:
    w_e = df["W_E"].to_numpy(dtype=np.float64)
    if w_e.size < 2:
        return "zu wenige Checkpoints"
    w_last, w_first = abs(w_e[-1]), abs(w_e[0])
    if w_last < w_first * 0.5:
        trend = "W_E tendiert gegen 0 → H0b-Kandidat"
    elif w_last > w_first * 1.5:
        trend = "W_E waechst → H3-Kandidat (Orientierung)"
    else:
        trend = "W_E ohne klaren Trend"
    return f"W_E: {w_e[0]:+.4e} → {w_e[-1]:+.4e}; {trend}"


def diagnose_scaling(
    alpha_E_hat: float, alpha_polyfit: float, slopes: dict[float, float], df: pd.DataFrame
) -> str:
    parts: list[str] = [diagnose_hypothesis(alpha_E_hat, alpha_polyfit)]
    if slopes.get(0.5, 0.0) > 0.1:
        parts.append(f"R_{{1/2}} waechst (Steigung={slopes[0.5]:.3f}) → alpha_E > 1/2")
    z = df["R_1_2"].to_numpy(dtype=np.float64)
    r23 = df["R_2_3"].to_numpy(dtype=np.float64)
    if coefficient_of_variation(z) < coefficient_of_variation(r23):
        parts.append("R_{1/2} stabiler als R_{2/3} → eher H0a")
    elif coefficient_of_variation(r23) < coefficient_of_variation(df["R_1"]):
        parts.append("R_{2/3} relativ stabil → alpha_E≈2/3 (H2, nicht H3)")
    w_stable = coefficient_of_variation(df["R_1"].to_numpy(dtype=np.float64))
    if w_stable < 0.25 and alpha_E_hat >= 0.85:
        parts.append("R_1 stabil → Holonomie-Hinweis (H3-Kandidat)")
    return "; ".join(parts)


def diagnose_alpha_eff(alpha: float, alpha_loc_max: float | None) -> str:
    parts: list[str] = [f"polyfit alpha={alpha:.4f} (Experiment, kein Satz)"]
    if alpha_loc_max is not None:
        parts.append(f"alpha_loc max={alpha_loc_max:.4f}")
        if alpha_loc_max > 0.5:
            parts.append("alpha_loc > 1/2 → numerisches Signal, kein asymptotischer Satz")
    return "; ".join(parts)


def interpret(
    alpha: float,
    alpha_E_hat: float,
    alpha_loc_max: float | None,
    df: pd.DataFrame,
    slopes: dict[float, float],
) -> str:
    return (
        f"{diagnose_orientation(df)} | "
        f"{diagnose_scaling(alpha_E_hat, alpha, slopes, df)} | "
        f"{diagnose_alpha_eff(alpha, alpha_loc_max)} | "
        "H0b ⇏ alpha_E≤1/2 — Orientierung und Skalierung getrennt"
    )


def make_diagnose_plot(df: pd.DataFrame, loc: pd.Series, out_path: Path) -> None:
    import matplotlib.pyplot as plt

    x = df["X"].to_numpy(dtype=np.float64)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    fig.suptitle("EABC-Quadruplets: Diagnose Ebenen 0–4 (H0a/H0b–H3)")

    ax = axes[0, 0]
    ax.plot(x, df["W_E"], "o-", color="C0", label=r"$W_E = R_1 = D/Q$")
    ax.axhline(0.0, color="gray", linewidth=0.8, linestyle="--")
    ax.set_ylabel(r"$W_E(X)$")
    ax.set_title("Ebene 3/4: Orientierung / Holonomie (H0b / H3)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(x, df["R_1_2"], "s-", color="C1", label=r"$R_{1/2}$ (Heuristik/Diagnose)")
    ax.axhline(0.0, color="gray", linewidth=0.8, linestyle="--")
    ax.set_ylabel(r"$R_{1/2}(X)$")
    ax.set_title("Ebene 2: R_{1/2} Skalierung (H0a / H1)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    loc_vals = loc.to_numpy(dtype=np.float64)
    ax.plot(x, loc_vals, "^-", color="C2", label=r"$\alpha_{\mathrm{loc}}$")
    ax.axhline(0.5, color="gray", linewidth=0.8, linestyle="--", label=r"$\alpha=1/2$")
    ax.set_ylabel(r"$\alpha_{\mathrm{loc}}$")
    ax.set_xlabel(r"$X$")
    ax.set_title("Ebene 2: alpha_loc-Diagnose (primär; kein Satz)")
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
    ax.set_title("Ebene 2: Skalierungsdiagnostik R_beta")
    ax.legend(loc="best", fontsize=7)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", nargs="?", default="eabc_quadruplets.csv")
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Vierfeld-Diagnose-Plot als PNG erzeugen",
    )
    parser.add_argument(
        "--plot-out",
        type=str,
        default="eabc_quadruplets_diagnose.png",
        help="Ausgabepfad für --plot (Standard: eabc_quadruplets_diagnose.png)",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    df = ensure_beta_columns(df)
    valid = load_valid(df)
    if valid.empty:
        raise SystemExit("Keine gültigen Zeilen (diff≠0, Q>1).")

    eff = alpha_eff_series(valid)
    loc = alpha_loc_series(valid)

    print("=== Ebene 0: Geometrie (G_E, gamma^+=ABCEA, gamma^-=CEABC) ===")
    print("  (implizit in Vierlingszählung; keine numerische Ausgabe)")

    print("\n=== Ebene 1: Zirkulationsfehler D_E(X) = N_+ - N_- ===")
    for idx, row in valid.iterrows():
        x = int(row["X"])
        print(
            f"  X={x:>12}  D={row['diff']:+.0f}  Q={row['Q_total']:.0f}  "
            f"R_1_2={row['R_1_2']:+.4e}  R_1={row['R_1']:+.4e}"
        )

    print("\n=== Ebene 2: Skalierung — alpha_loc (primär) und alpha_eff (sekundär) ===")
    for idx, row in valid.iterrows():
        x = int(row["X"])
        ae = eff.loc[idx]
        al_str = "—"
        if idx in loc.index and not np.isnan(loc.loc[idx]):
            al_str = f"{loc.loc[idx]:.4f}"
        print(f"  X={x:>12}  alpha_loc={al_str}  alpha_eff={ae:.4f}")

    alpha_loc_vals = loc.dropna()
    alpha_loc_max = float(alpha_loc_vals.max()) if not alpha_loc_vals.empty else None

    if len(valid) < 2:
        print("\nZu wenige Checkpoints für globalen Log-log-Fit (mind. 2).")
        if args.plot and len(valid) >= 1:
            out_path = Path(args.plot_out)
            make_diagnose_plot(valid, loc.reindex(valid.index), out_path)
            print(f"\nDiagnose-Plot gespeichert: {out_path.resolve()}")
        return

    alpha, beta = fit_alpha(valid)
    alpha_E_hat, best_beta, slopes = estimate_alpha_E(valid)

    print("\n=== Ebene 2: Asymptotischer Exponent (alpha_E_hat; Vermutung/Diagnostik) ===")
    print("R_beta-Steigungen (log|R_beta| vs log Q; ≈ alpha_E - beta):")
    for b in sorted(slopes):
        col = next(k for k, v in R_BETA_COLUMNS.items() if v == b)
        print(f"  beta={b:.2f} ({col}): Steigung={slopes[b]:+.4f}")

    print(f"\nalpha_E_hat (heuristisch, Experiment) = {alpha_E_hat:.4f}")
    print(f"  (Plateau bei beta={best_beta:.2f}, kein Theorem)")
    print(f"global alpha (polyfit) = {alpha:.4f}")
    print(f"→ {diagnose_scaling(alpha_E_hat, alpha, slopes, valid)}")
    print(f"→ {diagnose_alpha_eff(alpha, alpha_loc_max)}")

    print("\n=== Ebene 3: Orientierung (W_E; H0b) ===")
    print(f"→ {diagnose_orientation(valid)}")

    print("\n=== Ebene 4: Holonomie (Phi_E; H3-Kandidat) ===")
    if alpha_E_hat >= 0.85:
        print("  alpha_E_hat ≈ 1 → Holonomie-Grenzfall möglich (H3-Kandidat)")
    else:
        print("  alpha_E_hat < 1 → Phi_E ≠ 0 nicht aus Skalierung folgbar")

    print(f"\n=== Gesamt ===")
    print(f"→ {interpret(alpha, alpha_E_hat, alpha_loc_max, valid, slopes)}")

    if args.plot:
        out_path = Path(args.plot_out)
        make_diagnose_plot(valid, loc, out_path)
        print(f"\nDiagnose-Plot gespeichert: {out_path.resolve()}")


if __name__ == "__main__":
    main()
