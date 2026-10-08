#!/usr/bin/env python3
"""
EABC Crossover Gate v2
======================

Zweikanaliger Audit:
1) kumulative Skalierung F_<= (X)
2) disjunkte Shell-Skalierung F_shell(X_lo, X_hi)

Ziel:
- finite-range Drift und asymptotische Kandidaten trennen
- Korrelation kumulativer Cutoffs nicht als unabhängige Information missbrauchen
- M1, M2, M3, M12 vergleichen
- Rolling-Window-, OOS-, Plateau- und Crossover-Diagnostik
- Gate-Matrix G0..G5 ausgeben

Erwartete CSVs:

cumulative.csv:
    X,F[,sigma,N_prime,N_transitions]

shell.csv:
    X_lo,X_hi,F_shell[,sigma_shell,N_shell]

Beispiel:
    python eabc_crossover_gate_v2.py \
        --cumulative cumulative.csv \
        --shell shell.csv \
        --out crossover_v2_out
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy.optimize import brentq
except Exception:
    brentq = None

MODELS = {
    "M1": {
        "k": 1,
        "basis_names": ["1/logX"],
        "basis": lambda X: np.column_stack([1.0 / np.log(X)]),
    },
    "M2": {
        "k": 2,
        "basis_names": ["loglogX/logX^2", "1/logX^2"],
        "basis": lambda X: np.column_stack([
            np.log(np.log(X)) / np.log(X)**2,
            1.0 / np.log(X)**2
        ]),
    },
    "M3": {
        "k": 3,
        "basis_names": ["loglogX^2/logX^3", "loglogX/logX^3", "1/logX^3"],
        "basis": lambda X: np.column_stack([
            np.log(np.log(X))**2 / np.log(X)**3,
            np.log(np.log(X)) / np.log(X)**3,
            1.0 / np.log(X)**3
        ]),
    },
    "M12": {
        "k": 3,
        "basis_names": ["1/logX", "loglogX/logX^2", "1/logX^2"],
        "basis": lambda X: np.column_stack([
            1.0 / np.log(X),
            np.log(np.log(X)) / np.log(X)**2,
            1.0 / np.log(X)**2
        ]),
    },
}

DEFAULT_GRID = [
    1e4, 2e4, 5e4,
    1e5, 2e5, 5e5,
    1e6, 2e6, 5e6,
    1e7, 2e7, 5e7,
    1e8
]

def fit_model(name, X, y, sigma=None):
    spec = MODELS[name]
    A = spec["basis"](X)
    n, k = len(y), spec["k"]

    if sigma is not None:
        w = 1.0 / np.asarray(sigma, float)
        Aw = A * w[:, None]
        yw = y * w
        coef, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    else:
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)

    pred = A @ coef
    resid = y - pred
    rss = float(np.sum(resid**2))

    if sigma is not None:
        chi2 = float(np.sum((resid / sigma)**2))
        score = chi2
    else:
        chi2 = None
        score = rss

    eps = np.finfo(float).tiny
    aic = float(n*np.log(max(score/max(n,1), eps)) + 2*k)
    bic = float(n*np.log(max(score/max(n,1), eps)) + k*np.log(max(n,1)))
    aicc = None
    if n > k + 1:
        aicc = float(aic + 2*k*(k+1)/(n-k-1))

    return {
        "model": name,
        "n": n,
        "k": k,
        "coef": coef,
        "pred": pred,
        "resid": resid,
        "rss": rss,
        "chi2": chi2,
        "aic": aic,
        "aicc": aicc,
        "bic": bic,
        "rank": int(np.linalg.matrix_rank(A)),
        "condition_number": float(np.linalg.cond(A)),
    }

def predict(name, coef, X):
    X = np.asarray(X, float)
    return MODELS[name]["basis"](X) @ np.asarray(coef, float)

def rescaled_observables(X, F):
    L = np.log(X)
    LL = np.log(L)
    return pd.DataFrame({
        "X": X,
        "F": F,
        "R1_F_logX": F * L,
        "R2_F_logX2_over_loglogX": F * L**2 / LL,
        "R3_F_logX3_over_loglogX2": F * L**3 / LL**2,
    })

def effective_exponents(X, F):
    rows = []
    for i in range(len(X)-1):
        if F[i] == 0 or F[i+1] == 0:
            continue
        den = np.log(np.log(X[i+1]) / np.log(X[i]))
        if den == 0:
            continue
        alpha = -np.log(abs(F[i+1]/F[i])) / den
        rows.append({
            "X_left": X[i],
            "X_right": X[i+1],
            "alpha_eff": float(alpha),
        })
    return pd.DataFrame(rows)

def rolling_suffix(name, X, y, sigma=None):
    k = MODELS[name]["k"]
    min_n = max(k + 2, 4)
    rows = []
    for start in range(0, max(0, len(X)-min_n+1)):
        Xi = X[start:]
        yi = y[start:]
        if len(Xi) < min_n:
            continue
        si = None if sigma is None else sigma[start:]
        fit = fit_model(name, Xi, yi, si)
        row = {
            "model": name,
            "start_index": start,
            "X_min": float(Xi[0]),
            "X_max": float(Xi[-1]),
            "n": len(Xi),
            "rss": fit["rss"],
            "aicc": fit["aicc"],
            "bic": fit["bic"],
            "condition_number": fit["condition_number"],
        }
        for j, v in enumerate(fit["coef"]):
            row[f"coef_{j}"] = float(v)
        rows.append(row)
    return pd.DataFrame(rows)

def expanding_oos(name, X, y, sigma=None):
    k = MODELS[name]["k"]
    min_train = max(k + 2, 4)
    rows = []
    if len(X) <= min_train:
        return pd.DataFrame()

    for stop in range(min_train, len(X)):
        fit = fit_model(
            name, X[:stop], y[:stop],
            None if sigma is None else sigma[:stop]
        )
        pred = float(predict(name, fit["coef"], np.array([X[stop]]))[0])
        err = float(y[stop] - pred)
        z = None
        if sigma is not None and sigma[stop] > 0:
            z = float(err / sigma[stop])

        rows.append({
            "model": name,
            "train_n": stop,
            "X_train_max": float(X[stop-1]),
            "X_test": float(X[stop]),
            "F_test": float(y[stop]),
            "F_pred": pred,
            "error": err,
            "squared_error": err**2,
            "z_error": z,
        })
    return pd.DataFrame(rows)

def find_crossovers(name, coef, xmin, xmax, grid_n=5000):
    if brentq is None or xmin <= math.e:
        return []

    grid = np.exp(np.linspace(np.log(xmin), np.log(xmax), grid_n))
    vals = predict(name, coef, grid)
    roots = []

    for i in range(len(grid)-1):
        a, b = grid[i], grid[i+1]
        fa, fb = vals[i], vals[i+1]
        if not (np.isfinite(fa) and np.isfinite(fb)):
            continue
        if fa == 0:
            roots.append(float(a))
        elif fa * fb < 0:
            try:
                r = brentq(
                    lambda z: float(predict(name, coef, np.array([z]))[0]),
                    a, b
                )
                roots.append(float(r))
            except Exception:
                pass

    dedup = []
    for r in roots:
        if not dedup or abs(np.log(r/dedup[-1])) > 1e-7:
            dedup.append(r)
    return dedup

def model_table(X, y, sigma=None, crossover_max=1e100):
    rows = []
    fits = {}
    for name in MODELS:
        fit = fit_model(name, X, y, sigma)
        fits[name] = fit
        roots = find_crossovers(
            name, fit["coef"],
            xmin=max(float(np.min(X)), math.e + 1e-9),
            xmax=crossover_max
        )
        row = {
            "model": name,
            "n": fit["n"],
            "k": fit["k"],
            "rank": fit["rank"],
            "rss": fit["rss"],
            "chi2": fit["chi2"],
            "aic": fit["aic"],
            "aicc": fit["aicc"],
            "bic": fit["bic"],
            "condition_number": fit["condition_number"],
            "crossovers": ";".join(f"{r:.12g}" for r in roots),
        }
        for j, val in enumerate(fit["coef"]):
            row[f"coef_{j}"] = float(val)
            row[f"basis_{j}"] = MODELS[name]["basis_names"][j]
        rows.append(row)
    return pd.DataFrame(rows), fits

def plateau_cv(values):
    values = np.asarray(values, float)
    if len(values) < 2:
        return None
    m = np.mean(values)
    if m == 0:
        return None
    return float(np.std(values, ddof=1) / abs(m))

def assess_scale_candidate(diag):
    cvs = {}
    for key in [
        "R1_F_logX",
        "R2_F_logX2_over_loglogX",
        "R3_F_logX3_over_loglogX2",
    ]:
        cvs[key] = plateau_cv(diag[key].to_numpy(float))

    finite = {k:v for k,v in cvs.items() if v is not None and np.isfinite(v)}
    if not finite:
        return "UNRESOLVED", None, cvs

    best = min(finite, key=finite.get)
    mapping = {
        "R1_F_logX": "M1-like",
        "R2_F_logX2_over_loglogX": "M2-like",
        "R3_F_logX3_over_loglogX2": "M3-leading-like",
    }
    return "PRELIMINARY", mapping[best], cvs

def gate_matrix(cum_df, shell_df, cumulative_diag, shell_diag, cumulative_fits, shell_fits):
    # G0: signal observed if all finite and not identically zero
    g0 = "OBSERVED" if np.any(np.abs(cum_df["F"].to_numpy(float)) > 0) else "UNRESOLVED"

    # G1 requires external N0/N1/N2 audit; this script must not invent it.
    g1 = "PENDING AUDIT"

    g2_status, g2_best, g2_cvs = assess_scale_candidate(cumulative_diag)

    # G3: needs enough shell + rolling/oos evidence.
    if shell_df is not None and len(shell_df) >= 8 and len(cum_df) >= 8:
        g3 = "AUDITABLE"
    else:
        g3 = "UNRESOLVED"

    # G4: conservative. A fitted root alone is not enough; need robustness.
    g4 = "UNRESOLVED"

    # G5 is analytic theory, never auto-passed from numerics.
    g5 = "OPEN"

    return {
        "G0_SIGNAL": g0,
        "G1_NULL_RESIDUAL": g1,
        "G2_SCALE_CANDIDATE": {
            "status": g2_status,
            "candidate": g2_best,
            "plateau_cv": g2_cvs,
        },
        "G3_SCALE_STABLE": g3,
        "G4_CROSSOVER_SAFE": g4,
        "G5_THEORY_MATCH": g5,
    }

def load_cumulative(path):
    df = pd.read_csv(path)
    if not {"X","F"}.issubset(df.columns):
        raise ValueError("cumulative CSV requires columns X,F")
    df = df.dropna(subset=["X","F"]).sort_values("X").reset_index(drop=True)
    return df

def load_shell(path):
    df = pd.read_csv(path)
    if not {"X_lo","X_hi","F_shell"}.issubset(df.columns):
        raise ValueError("shell CSV requires columns X_lo,X_hi,F_shell")
    df = df.dropna(subset=["X_lo","X_hi","F_shell"]).sort_values(["X_lo","X_hi"]).reset_index(drop=True)

    # Representative scale = geometric midpoint, natural for logarithmic shells.
    df["X_mid"] = np.sqrt(df["X_lo"].to_numpy(float) * df["X_hi"].to_numpy(float))
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cumulative", required=True)
    ap.add_argument("--shell", default=None)
    ap.add_argument("--out", default="crossover_v2_out")
    ap.add_argument("--crossover-max", type=float, default=1e100)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    cum = load_cumulative(args.cumulative)
    Xc = cum["X"].to_numpy(float)
    Fc = cum["F"].to_numpy(float)
    Sc = cum["sigma"].to_numpy(float) if "sigma" in cum.columns else None

    if np.any(Xc <= math.e):
        raise SystemExit("All cumulative X must be > e.")
    if Sc is not None and np.any(Sc <= 0):
        raise SystemExit("All cumulative sigma values must be positive.")

    cum_diag = rescaled_observables(Xc, Fc)
    for col in cum.columns:
        if col not in cum_diag.columns:
            cum_diag[col] = cum[col].values

    cum_fits, cum_fit_objects = model_table(Xc, Fc, Sc, args.crossover_max)

    roll_parts = []
    oos_parts = []
    for name in MODELS:
        r = rolling_suffix(name, Xc, Fc, Sc)
        if not r.empty:
            roll_parts.append(r)
        o = expanding_oos(name, Xc, Fc, Sc)
        if not o.empty:
            oos_parts.append(o)

    rolling = pd.concat(roll_parts, ignore_index=True) if roll_parts else pd.DataFrame()
    oos = pd.concat(oos_parts, ignore_index=True) if oos_parts else pd.DataFrame()
    eff = effective_exponents(Xc, Fc)

    shell = None
    shell_diag = pd.DataFrame()
    shell_fits = pd.DataFrame()

    if args.shell:
        shell = load_shell(args.shell)
        Xs = shell["X_mid"].to_numpy(float)
        Fs = shell["F_shell"].to_numpy(float)
        Ss = shell["sigma_shell"].to_numpy(float) if "sigma_shell" in shell.columns else None

        if np.any(Xs <= math.e):
            raise SystemExit("All shell X_mid must be > e.")
        if Ss is not None and np.any(Ss <= 0):
            raise SystemExit("All sigma_shell values must be positive.")

        shell_diag = rescaled_observables(Xs, Fs).rename(columns={"F":"F_shell"})
        shell_diag["X_lo"] = shell["X_lo"].values
        shell_diag["X_hi"] = shell["X_hi"].values
        shell_diag["X_mid"] = shell["X_mid"].values
        for col in shell.columns:
            if col not in shell_diag.columns:
                shell_diag[col] = shell[col].values

        shell_fits, _ = model_table(Xs, Fs, Ss, args.crossover_max)

    gates = gate_matrix(cum, shell, cum_diag, shell_diag, cum_fits, shell_fits)

    # Freeze the recommended 1-2-5 grid for reproducibility.
    grid_df = pd.DataFrame({"X": DEFAULT_GRID})

    cum.to_csv(out / "cumulative_input_sorted.csv", index=False)
    cum_diag.to_csv(out / "cumulative_rescaled.csv", index=False)
    cum_fits.to_csv(out / "cumulative_model_fits.csv", index=False)
    eff.to_csv(out / "cumulative_effective_exponents.csv", index=False)
    grid_df.to_csv(out / "frozen_1-2-5_grid.csv", index=False)

    if not rolling.empty:
        rolling.to_csv(out / "cumulative_rolling_windows.csv", index=False)
    if not oos.empty:
        oos.to_csv(out / "cumulative_oos.csv", index=False)

    if shell is not None:
        shell.to_csv(out / "shell_input_sorted.csv", index=False)
        shell_diag.to_csv(out / "shell_rescaled.csv", index=False)
        shell_fits.to_csv(out / "shell_model_fits.csv", index=False)

    report = {
        "gate": "G_Cross^chi_v2",
        "method": {
            "cumulative_role": "global drift / plateau visualization; points are correlated",
            "shell_role": "local logarithmic-scale dynamics; intended to reduce cumulative dependence",
            "warning": (
                "Shell observations are less cumulatively dependent but are not automatically "
                "statistically independent if the underlying observable couples across shell boundaries."
            ),
            "theory_warning": (
                "A_chi^perp = 0 can motivate a subleading candidate but does not by itself derive "
                "the exact loglog(X)/(log X)^2 order."
            ),
        },
        "gates": gates,
        "n_cumulative": int(len(cum)),
        "n_shell": int(len(shell)) if shell is not None else 0,
        "X_cumulative_min": float(Xc.min()),
        "X_cumulative_max": float(Xc.max()),
    }
    (out / "gate_matrix.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    print("="*88)
    print(" EABC Crossover Gate v2 — cumulative + disjoint shells")
    print("="*88)
    print(f"Cumulative points: {len(cum)}  X=[{Xc.min():.6g}, {Xc.max():.6g}]")
    print(f"Shell points:      {0 if shell is None else len(shell)}")
    print()
    print("Gate matrix:")
    for k, v in gates.items():
        print(f"  {k}: {v}")
    print()
    print("Cumulative model fits:")
    cols = ["model","n","k","rss","aicc","bic","condition_number","crossovers"]
    print(cum_fits[cols].to_string(index=False))
    if shell is not None:
        print()
        print("Shell model fits:")
        print(shell_fits[cols].to_string(index=False))
    print()
    print("Outputs:", out.resolve())

if __name__ == "__main__":
    main()
