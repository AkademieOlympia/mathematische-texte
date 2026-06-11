import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path

from math import lgamma
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import r2_score, mean_squared_error
import warnings

# ============================================================
# Erwartete Eingabedaten
# ------------------------------------------------------------
# DataFrame df mit Spalten:
#   sector_S     z.B. "(0,0,0)"
#   sector_support   z.B. 1
#   sector_type      z.B. "B"
#   log_mB           = log(m_B)
#   C_fam            = Familienentropie
#   kappa13          = Kleinträgerprofil
#   log_mult         = log(Hurwitz-Multiplizität)
# ============================================================

REQUIRED_COLUMNS = [
    "sector_S",
    "sector_support",
    "sector_type",
    "log_mB",
    "C_fam",
    "kappa13",
    "log_mult",
]

ARITH_STRUCTURE_COLUMNS = [
    "mB",
    "C_fam",
    "kappa_13",
    "fam_support",
    "rho_E",
    "rho_A",
    "rho_B",
    "rho_C",
]

DEFAULT_INPUT_CANDIDATES = [
    "napoleon_input.csv",
    "arith_struktur_v3_bis_100000.csv",
    "arith_struktur_bis_1000.csv",
]

DEFAULT_MAX_SECTORS = 3
DEFAULT_MAX_ROWS_PER_SECTOR = 15

# ----------------------------
# Hilfsfunktionen
# ----------------------------
def rmse(y_true, y_pred):
    return np.sqrt(mean_squared_error(y_true, y_pred))

def aic_bic(y_true, y_pred, n_params):
    n = len(y_true)
    rss = np.sum((y_true - y_pred) ** 2)
    rss = max(rss, 1e-12)
    aic = n * np.log(rss / n) + 2 * n_params
    bic = n * np.log(rss / n) + n_params * np.log(n)
    return aic, bic

def loocv_linear(X, y):
    loo = LeaveOneOut()
    preds = np.zeros_like(y, dtype=float)
    for train_idx, test_idx in loo.split(X):
        model = LinearRegression()
        model.fit(X[train_idx], y[train_idx])
        preds[test_idx] = model.predict(X[test_idx])
    return preds

def loocv_nonlinear(
    x, c, k, y, model_func, p0, bounds=(-np.inf, np.inf), maxfev=100000
):
    loo = LeaveOneOut()
    preds = np.full_like(y, np.nan, dtype=float)
    success = []
    for train_idx, test_idx in loo.split(y):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", OptimizeWarning)
                params, _ = curve_fit(
                    model_func,
                    (x[train_idx], c[train_idx], k[train_idx]),
                    y[train_idx],
                    p0=p0,
                    bounds=bounds,
                    maxfev=maxfev
                )
            preds[test_idx] = model_func(
                (x[test_idx], c[test_idx], k[test_idx]), *params
            )
            success.append(True)
        except Exception:
            success.append(False)
    mask = ~np.isnan(preds)
    return preds, np.mean(success), mask

# ----------------------------
# Modellfunktionen
# ----------------------------
def model_A_fit(X, y):
    model = LinearRegression()
    model.fit(X, y)
    yhat = model.predict(X)
    return model, yhat

def model_B_fit(X, y):
    model = LinearRegression()
    model.fit(X, y)
    yhat = model.predict(X)
    return model, yhat

def gamma_model(vars_, a0, a1, a2, a3, beta, l0, l1, l2, l3):
    x, c, k = vars_
    z = l0 + l1 * x + l2 * c + l3 * k
    z = np.maximum(z, 1e-8)
    lg = np.vectorize(lgamma)(z)
    return a0 + a1 * x + a2 * c + a3 * k + beta * lg

def gamma_model_exp(vars_, a0, a1, a2, a3, beta, m0, m1, m2, m3):
    x, c, k = vars_
    z = np.exp(m0 + m1 * x + m2 * c + m3 * k)
    lg = np.vectorize(lgamma)(z)
    return a0 + a1 * x + a2 * c + a3 * k + beta * lg


def model_A_matrix(x, c, k):
    return np.column_stack([x, c, k])


def model_B1_matrix(x, c, k, delta=1.0):
    return np.column_stack([x, c, k, 1.0 / (x + delta)])


def model_B2_matrix_fixed(x, c, k):
    X = 1.0 + x + c + k
    return np.column_stack([x, c, k, 1.0 / X])


def model_C1(vars_, a0, a1, a2, a3, beta, delta):
    x, c, k = vars_
    z = x + delta
    z = np.maximum(z, 1e-8)
    lg = np.vectorize(lgamma)(z)
    return a0 + a1 * x + a2 * c + a3 * k + beta * lg


def model_C2(vars_, a0, a1, a2, a3, beta, lam, delta):
    x, c, k = vars_
    z = lam * x + delta
    z = np.maximum(z, 1e-8)
    lg = np.vectorize(lgamma)(z)
    return a0 + a1 * x + a2 * c + a3 * k + beta * lg


def model_D1(vars_, a0, a1, a2, a3, betaA, betaB):
    x, c, k = vars_
    zA = np.maximum(x + 1.0, 1e-8)
    zB = np.maximum(k + 1.0, 1e-8)
    lgA = np.vectorize(lgamma)(zA)
    lgB = np.vectorize(lgamma)(zB)
    return a0 + a1 * x + a2 * c + a3 * k + betaA * lgA + betaB * lgB


def fit_linear_model(X, y):
    model = LinearRegression()
    model.fit(X, y)
    yhat = model.predict(X)
    return model, yhat


def fit_D1(x, c, k, y):
    p0 = [0.0, 1.0, 1.0, 1.0, 0.05, 0.05]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OptimizeWarning)
        params, _ = curve_fit(
            model_D1,
            (x, c, k),
            y,
            p0=p0,
            maxfev=100000,
        )
    yhat = model_D1((x, c, k), *params)
    return params, yhat


def collect_metrics(
    sector_name, model_name, y, yhat, loocv_pred, n_params, mask=None, extra=None
):
    if mask is None:
        mask = np.ones_like(y, dtype=bool)

    result = {
        "sector": sector_name,
        "model": model_name,
        "n": len(y),
        "R2": r2_score(y, yhat),
        "RMSE": rmse(y, yhat),
        "AIC": aic_bic(y, yhat, n_params)[0],
        "BIC": aic_bic(y, yhat, n_params)[1],
        "LOOCV_R2": r2_score(y[mask], loocv_pred[mask]) if np.sum(mask) > 2 else np.nan,
        "LOOCV_RMSE": rmse(y[mask], loocv_pred[mask]) if np.sum(mask) > 2 else np.nan,
    }
    if extra:
        result.update(extra)
    return result


def residual_structure_report(resid, x, c, k, label="model"):
    out = {}

    out["label"] = label
    out["var"] = float(np.var(resid, ddof=1))
    out["std"] = float(np.std(resid, ddof=1))

    def safe_corr(a, b):
        if np.std(a) < 1e-12 or np.std(b) < 1e-12:
            return np.nan
        return float(np.corrcoef(a, b)[0, 1])

    out["corr_x"] = safe_corr(resid, x)
    out["corr_c"] = safe_corr(resid, c)
    out["corr_k"] = safe_corr(resid, k)

    X_quad = np.column_stack([x, x**2])
    reg = LinearRegression().fit(X_quad, resid)
    resid_hat = reg.predict(X_quad)
    out["drift_R2_x_quad"] = float(r2_score(resid, resid_hat))

    Xc_quad = np.column_stack([c, c**2])
    regc = LinearRegression().fit(Xc_quad, resid)
    resid_hat_c = regc.predict(Xc_quad)
    out["drift_R2_c_quad"] = float(r2_score(resid, resid_hat_c))

    Xk_quad = np.column_stack([k, k**2])
    regk = LinearRegression().fit(Xk_quad, resid)
    resid_hat_k = regk.predict(Xk_quad)
    out["drift_R2_k_quad"] = float(r2_score(resid, resid_hat_k))

    signs = np.sign(resid)
    signs = signs[signs != 0]
    if len(signs) > 1:
        sign_changes = np.sum(signs[1:] != signs[:-1])
        out["sign_changes"] = int(sign_changes)
        out["sign_change_rate"] = float(sign_changes / (len(signs) - 1))
    else:
        out["sign_changes"] = np.nan
        out["sign_change_rate"] = np.nan

    return out


def compare_residuals_A_C1(
    df_sector, params_C1, delta_C1=1.0, sector_name="sector", plot=True
):
    d = df_sector.copy().sort_values("log_mB")

    x = d["log_mB"].to_numpy(dtype=float)
    c = d["C_fam"].to_numpy(dtype=float)
    k = d["kappa13"].to_numpy(dtype=float)
    y = d["log_mult"].to_numpy(dtype=float)

    XA = np.column_stack([x, c, k])
    regA = LinearRegression().fit(XA, y)
    yhatA = regA.predict(XA)
    RA = y - yhatA

    a0, a1, a2, a3, beta = params_C1
    z = np.maximum(x + delta_C1, 1e-8)
    yhatC1 = a0 + a1 * x + a2 * c + a3 * k + beta * np.vectorize(lgamma)(z)
    RC1 = y - yhatC1

    repA = residual_structure_report(RA, x, c, k, label="A")
    repC1 = residual_structure_report(RC1, x, c, k, label="C1")

    report_df = pd.DataFrame([repA, repC1])
    report_df["sector"] = sector_name

    var_ratio = repC1["var"] / repA["var"] if repA["var"] > 0 else np.nan
    drift_gain_x = repA["drift_R2_x_quad"] - repC1["drift_R2_x_quad"]
    drift_gain_c = repA["drift_R2_c_quad"] - repC1["drift_R2_c_quad"]
    drift_gain_k = repA["drift_R2_k_quad"] - repC1["drift_R2_k_quad"]

    summary = pd.DataFrame(
        [
            {
                "sector": sector_name,
                "var_ratio_C1_over_A": var_ratio,
                "drift_gain_x": drift_gain_x,
                "drift_gain_c": drift_gain_c,
                "drift_gain_k": drift_gain_k,
                "abs_corr_x_reduction": abs(repA["corr_x"]) - abs(repC1["corr_x"]),
                "abs_corr_c_reduction": abs(repA["corr_c"]) - abs(repC1["corr_c"]),
                "abs_corr_k_reduction": abs(repA["corr_k"]) - abs(repC1["corr_k"]),
            }
        ]
    )

    if plot:
        fig, axes = plt.subplots(2, 3, figsize=(13, 8))

        axes[0, 0].scatter(x, RA, s=20, label="A")
        axes[0, 0].scatter(x, RC1, s=20, label="C1")
        axes[0, 0].axhline(0.0, linestyle="--")
        axes[0, 0].set_title("Residuen vs. log(m_B)")
        axes[0, 0].set_xlabel("log(m_B)")
        axes[0, 0].set_ylabel("Residuum")
        axes[0, 0].legend()

        axes[0, 1].scatter(c, RA, s=20, label="A")
        axes[0, 1].scatter(c, RC1, s=20, label="C1")
        axes[0, 1].axhline(0.0, linestyle="--")
        axes[0, 1].set_title("Residuen vs. C_fam")
        axes[0, 1].set_xlabel("C_fam")
        axes[0, 1].set_ylabel("Residuum")
        axes[0, 1].legend()

        axes[0, 2].scatter(k, RA, s=20, label="A")
        axes[0, 2].scatter(k, RC1, s=20, label="C1")
        axes[0, 2].axhline(0.0, linestyle="--")
        axes[0, 2].set_title("Residuen vs. kappa13")
        axes[0, 2].set_xlabel("kappa13")
        axes[0, 2].set_ylabel("Residuum")
        axes[0, 2].legend()

        axes[1, 0].plot(np.arange(len(RA)), RA, marker="o", label="A")
        axes[1, 0].plot(np.arange(len(RC1)), RC1, marker="o", label="C1")
        axes[1, 0].axhline(0.0, linestyle="--")
        axes[1, 0].set_title("Sortierte Residuen")
        axes[1, 0].set_xlabel("Rang nach log(m_B)")
        axes[1, 0].set_ylabel("Residuum")
        axes[1, 0].legend()

        axes[1, 1].hist(RA, bins=10, alpha=0.6, label="A")
        axes[1, 1].hist(RC1, bins=10, alpha=0.6, label="C1")
        axes[1, 1].set_title("Residuenhistogramm")
        axes[1, 1].legend()

        axes[1, 2].boxplot([RA, RC1], tick_labels=["A", "C1"])
        axes[1, 2].set_title("Residuenvergleich")

        plt.suptitle(f"Residuenvergleich A vs. C1 - {sector_name}")
        plt.tight_layout()
        plt.show()

    return report_df, summary


def analyze_rest_C1(x, c, k, r, delta=1.0):
    rows = []

    for name, z in [("x", x), ("c", c), ("k", k)]:
        X_lin = z.reshape(-1, 1)
        X_quad = np.column_stack([z, z**2])

        m_lin = LinearRegression().fit(X_lin, r)
        rhat_lin = m_lin.predict(X_lin)

        m_quad = LinearRegression().fit(X_quad, r)
        rhat_quad = m_quad.predict(X_quad)

        rows.append(
            {
                "test": f"drift_{name}_linear",
                "R2": r2_score(r, rhat_lin),
                "RMSE": rmse(r, rhat_lin),
            }
        )
        rows.append(
            {
                "test": f"drift_{name}_quadratic",
                "R2": r2_score(r, rhat_quad),
                "RMSE": rmse(r, rhat_quad),
            }
        )

    X_inv1 = (1.0 / (x + delta)).reshape(-1, 1)
    m_inv1 = LinearRegression().fit(X_inv1, r)
    rhat_inv1 = m_inv1.predict(X_inv1)
    loocv_inv1 = loocv_linear(X_inv1, r)

    rows.append(
        {
            "test": f"inv1_x_shift_{delta}",
            "R2": r2_score(r, rhat_inv1),
            "RMSE": rmse(r, rhat_inv1),
            "LOOCV_R2": r2_score(r, loocv_inv1),
            "LOOCV_RMSE": rmse(r, loocv_inv1),
        }
    )

    X_inv13 = np.column_stack([1.0 / (x + delta), 1.0 / ((x + delta) ** 3)])
    m_inv13 = LinearRegression().fit(X_inv13, r)
    rhat_inv13 = m_inv13.predict(X_inv13)
    loocv_inv13 = loocv_linear(X_inv13, r)

    rows.append(
        {
            "test": f"inv13_x_shift_{delta}",
            "R2": r2_score(r, rhat_inv13),
            "RMSE": rmse(r, rhat_inv13),
            "LOOCV_R2": r2_score(r, loocv_inv13),
            "LOOCV_RMSE": rmse(r, loocv_inv13),
        }
    )

    Xcomb = 1.0 + x + c + k
    X_invcomb = (1.0 / Xcomb).reshape(-1, 1)
    m_invcomb = LinearRegression().fit(X_invcomb, r)
    rhat_invcomb = m_invcomb.predict(X_invcomb)
    loocv_invcomb = loocv_linear(X_invcomb, r)

    rows.append(
        {
            "test": "inv1_combined",
            "R2": r2_score(r, rhat_invcomb),
            "RMSE": rmse(r, rhat_invcomb),
            "LOOCV_R2": r2_score(r, loocv_invcomb),
            "LOOCV_RMSE": rmse(r, loocv_invcomb),
        }
    )

    signs = np.sign(r)
    signs = signs[signs != 0]
    if len(signs) > 1:
        sign_changes = np.sum(signs[1:] != signs[:-1])
        sign_change_rate = sign_changes / (len(signs) - 1)
    else:
        sign_changes = np.nan
        sign_change_rate = np.nan

    summary = {
        "resid_var": float(np.var(r, ddof=1)),
        "resid_std": float(np.std(r, ddof=1)),
        "sign_changes": sign_changes,
        "sign_change_rate": sign_change_rate,
    }

    return pd.DataFrame(rows), pd.DataFrame([summary])


def spectral_preanalysis_rest(
    x, r, sector_name="sector", top_k=5, detrend_linear=False, plot=True
):
    idx = np.argsort(x)
    x_sorted = np.asarray(x)[idx]
    r_sorted = np.asarray(r)[idx]

    r_work = r_sorted - np.mean(r_sorted)

    if detrend_linear:
        X = np.column_stack([np.ones(len(x_sorted)), x_sorted])
        beta = np.linalg.lstsq(X, r_work, rcond=None)[0]
        r_work = r_work - X @ beta

    n = len(r_work)

    fft_vals = np.fft.rfft(r_work)
    power = np.abs(fft_vals) ** 2
    freqs = np.fft.rfftfreq(n, d=1.0)

    if len(power) > 0:
        power_no0 = power.copy()
        power_no0[0] = 0.0
    else:
        power_no0 = power

    peak_idx = np.argsort(power_no0)[::-1][:top_k]
    peaks = pd.DataFrame(
        {
            "rank": np.arange(1, len(peak_idx) + 1),
            "freq_index": freqs[peak_idx],
            "power": power_no0[peak_idx],
            "detrended": detrend_linear,
        }
    )

    summary = pd.DataFrame(
        [
            {
                "sector": sector_name,
                "n": n,
                "resid_mean": float(np.mean(r_sorted)),
                "resid_std": float(np.std(r_sorted, ddof=1)),
                "max_power": float(np.max(power_no0)) if len(power_no0) else np.nan,
                "power_sum": float(np.sum(power_no0)) if len(power_no0) else np.nan,
                "peak_to_total_ratio": (
                    float(np.max(power_no0) / np.sum(power_no0))
                    if np.sum(power_no0) > 0
                    else np.nan
                ),
                "detrended": detrend_linear,
            }
        ]
    )

    sequence_df = pd.DataFrame(
        {
            "rank_index": np.arange(n),
            "x_sorted": x_sorted,
            "r_sorted": r_sorted,
            "r_used": r_work,
        }
    )

    if plot:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        axes[0].plot(np.arange(n), r_sorted, marker="o")
        axes[0].axhline(0.0, linestyle="--")
        axes[0].set_title(f"Restfolge {sector_name}")
        axes[0].set_xlabel("Rang nach log(m_B)")
        axes[0].set_ylabel("R_C1")

        axes[1].stem(freqs[1:], power_no0[1:], basefmt=" ")
        axes[1].set_title(f"Periodogramm {sector_name}")
        axes[1].set_xlabel("Index-Frequenz")
        axes[1].set_ylabel("Leistung")

        plt.tight_layout()
        plt.show()

    return summary, peaks, sequence_df


def compare_spectral_variants(summary_df, peaks_df):
    if summary_df.empty or peaks_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    summary_false = (
        summary_df[summary_df["detrended"] == False]
        .drop(columns=["detrended"])
        .rename(
            columns={
                "max_power": "max_power_raw",
                "power_sum": "power_sum_raw",
                "peak_to_total_ratio": "peak_to_total_ratio_raw",
            }
        )
    )
    summary_true = (
        summary_df[summary_df["detrended"] == True]
        .drop(columns=["detrended"])
        .rename(
            columns={
                "max_power": "max_power_detrended",
                "power_sum": "power_sum_detrended",
                "peak_to_total_ratio": "peak_to_total_ratio_detrended",
            }
        )
    )

    summary_compare = summary_false.merge(summary_true, on=["sector", "n", "resid_mean", "resid_std"], how="inner")
    if not summary_compare.empty:
        summary_compare["delta_max_power"] = (
            summary_compare["max_power_detrended"] - summary_compare["max_power_raw"]
        )
        summary_compare["delta_power_sum"] = (
            summary_compare["power_sum_detrended"] - summary_compare["power_sum_raw"]
        )
        summary_compare["delta_peak_to_total_ratio"] = (
            summary_compare["peak_to_total_ratio_detrended"]
            - summary_compare["peak_to_total_ratio_raw"]
        )

    peaks_false = (
        peaks_df[peaks_df["detrended"] == False]
        .drop(columns=["detrended"])
        .rename(
            columns={
                "freq_index": "freq_index_raw",
                "power": "power_raw",
            }
        )
    )
    peaks_true = (
        peaks_df[peaks_df["detrended"] == True]
        .drop(columns=["detrended"])
        .rename(
            columns={
                "freq_index": "freq_index_detrended",
                "power": "power_detrended",
            }
        )
    )

    peaks_compare = peaks_false.merge(peaks_true, on=["sector", "rank"], how="inner")
    if not peaks_compare.empty:
        peaks_compare["delta_freq_index"] = (
            peaks_compare["freq_index_detrended"] - peaks_compare["freq_index_raw"]
        )
        peaks_compare["delta_power"] = (
            peaks_compare["power_detrended"] - peaks_compare["power_raw"]
        )

    return summary_compare, peaks_compare


def validate_dataframe(df):
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(f"Fehlende Spalten in Eingabedaten: {missing}")


def is_arith_structure_dataframe(df):
    return all(col in df.columns for col in ARITH_STRUCTURE_COLUMNS)


def family_support_size(label):
    if label == "1" or not label:
        return 0
    return len(str(label).split("·"))


def build_sector_signature(row):
    active = tuple(
        int(row[col] > 1e-12)
        for col in ("rho_E", "rho_A", "rho_B", "rho_C")
    )
    return str(active)


def build_napoleon_input_from_arith(df):
    d = df.copy()
    d = d[d["mB"] > 1].copy()
    if d.empty:
        raise ValueError("Die Strukturdaten enthalten keine Restkerne mB > 1.")

    d["sector_type"] = d["fam_support"].astype(str)
    d["sector_support"] = d["sector_type"].map(family_support_size)
    d["sector_S"] = d.apply(build_sector_signature, axis=1)

    grouped = (
        d.groupby(["sector_S", "sector_support", "sector_type", "mB"], dropna=False)
        .agg(
            C_fam=("C_fam", "mean"),
            kappa13=("kappa_13", "mean"),
            multiplicity=("mB", "size"),
        )
        .reset_index()
    )

    grouped["log_mB"] = np.log(grouped["mB"].astype(float))
    grouped["log_mult"] = np.log(grouped["multiplicity"].astype(float))

    return grouped[
        [
            "sector_S",
            "sector_support",
            "sector_type",
            "log_mB",
            "C_fam",
            "kappa13",
            "log_mult",
        ]
    ]


def discover_sectors(df, min_rows=8, max_sectors=DEFAULT_MAX_SECTORS):
    grouped = (
        df.groupby(["sector_S", "sector_support", "sector_type"])
        .size()
        .reset_index(name="n_rows")
        .sort_values(["n_rows", "sector_S", "sector_support", "sector_type"], ascending=[False, True, True, True])
    )
    valid = grouped[grouped["n_rows"] >= min_rows]
    return [
        (row["sector_S"], row["sector_support"], row["sector_type"])
        for _, row in valid.head(max_sectors).iterrows()
    ]


def reduce_sector_dataframe(df_sector, max_rows=DEFAULT_MAX_ROWS_PER_SECTOR):
    if len(df_sector) <= max_rows:
        return df_sector.copy()

    d = df_sector.sort_values("log_mB").reset_index(drop=True)
    idx = np.linspace(0, len(d) - 1, num=max_rows, dtype=int)
    return d.iloc[idx].copy()

# ----------------------------
# Analyse eines Sektors
# ----------------------------
def analyze_sector(df_sector, sector_name="sector", use_exp_gamma=False, plot=True):
    d = reduce_sector_dataframe(df_sector).sort_values("log_mB")
    x = d["log_mB"].to_numpy(dtype=float)
    c = d["C_fam"].to_numpy(dtype=float)
    k = d["kappa13"].to_numpy(dtype=float)
    y = d["log_mult"].to_numpy(dtype=float)

    # Modell A
    X_A = np.column_stack([x, c, k])
    modelA, yhatA = model_A_fit(X_A, y)
    predA_loocv = loocv_linear(X_A, y)

    # Modell B
    X_B = np.column_stack([x, c, k, 1.0 / x, 1.0 / (x ** 3)])
    modelB, yhatB = model_B_fit(X_B, y)
    predB_loocv = loocv_linear(X_B, y)

    # Modell C
    if use_exp_gamma:
        modelC_func = gamma_model_exp
        p0 = [0.0, 1.0, 1.0, 1.0, 0.05, 0.0, 0.05, 0.05, 0.05]
        n_params_C = 9
    else:
        modelC_func = gamma_model
        # l0 bewusst positiv
        p0 = [0.0, 1.0, 1.0, 1.0, 0.05, 2.0, 0.05, 0.05, 0.05]
        n_params_C = 9

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OptimizeWarning)
        paramsC, _ = curve_fit(
            modelC_func, (x, c, k), y, p0=p0, maxfev=100000
        )
    yhatC = modelC_func((x, c, k), *paramsC)

    predC_loocv, success_rate_C, maskC = loocv_nonlinear(
        x, c, k, y, modelC_func, p0=p0, maxfev=100000
    )

    # Kennzahlen
    results = []

    def collect(name, yhat, loocv_pred, n_params, mask=None):
        if mask is None:
            mask = np.ones_like(y, dtype=bool)
        res = {
            "sector": sector_name,
            "model": name,
            "n": len(y),
            "R2": r2_score(y, yhat),
            "RMSE": rmse(y, yhat),
            "AIC": aic_bic(y, yhat, n_params)[0],
            "BIC": aic_bic(y, yhat, n_params)[1],
            "LOOCV_R2": r2_score(y[mask], loocv_pred[mask]) if np.sum(mask) > 2 else np.nan,
            "LOOCV_RMSE": rmse(y[mask], loocv_pred[mask]) if np.sum(mask) > 2 else np.nan,
        }
        return res

    results.append(collect("A_linear", yhatA, predA_loocv, n_params=4))
    results.append(collect("B_linear_invpow", yhatB, predB_loocv, n_params=6))
    results.append(collect("C_logGamma", yhatC, predC_loocv, n_params=n_params_C, mask=maskC))

    metrics_df = pd.DataFrame(results)

    # Residuen
    resid_df = pd.DataFrame({
        "log_mB": x,
        "C_fam": c,
        "kappa13": k,
        "y": y,
        "resid_A": y - yhatA,
        "resid_B": y - yhatB,
        "resid_C": y - yhatC,
    })

    if plot:
        fig, axes = plt.subplots(2, 2, figsize=(11, 8))

        # Fitvergleich
        axes[0, 0].scatter(x, y, s=20, label="Daten")
        axes[0, 0].plot(x, yhatA, label="A")
        axes[0, 0].plot(x, yhatB, label="B")
        axes[0, 0].plot(x, yhatC, label="C")
        axes[0, 0].set_title(f"{sector_name}: Fits vs. log(m_B)")
        axes[0, 0].set_xlabel("log(m_B)")
        axes[0, 0].set_ylabel("log Multiplizität")
        axes[0, 0].legend()

        # Residuen A
        axes[0, 1].scatter(x, y - yhatA, s=20)
        axes[0, 1].axhline(0.0, linestyle="--")
        axes[0, 1].set_title("Residuen Modell A")
        axes[0, 1].set_xlabel("log(m_B)")
        axes[0, 1].set_ylabel("Residuum")

        # Residuen B
        axes[1, 0].scatter(x, y - yhatB, s=20)
        axes[1, 0].axhline(0.0, linestyle="--")
        axes[1, 0].set_title("Residuen Modell B")
        axes[1, 0].set_xlabel("log(m_B)")
        axes[1, 0].set_ylabel("Residuum")

        # Residuen C
        axes[1, 1].scatter(x, y - yhatC, s=20)
        axes[1, 1].axhline(0.0, linestyle="--")
        axes[1, 1].set_title("Residuen Modell C")
        axes[1, 1].set_xlabel("log(m_B)")
        axes[1, 1].set_ylabel("Residuum")

        plt.tight_layout()
        plt.show()

    return metrics_df, resid_df, paramsC

# ----------------------------
# Mehrere Sektoren vergleichen
# ----------------------------
def compare_sectors(df, sectors, use_exp_gamma=False, plot=True):
    all_metrics = []
    all_params = []

    for sector_S, sector_support, sector_type in sectors:
        mask = (
            (df["sector_S"] == sector_S) &
            (df["sector_support"] == sector_support) &
            (df["sector_type"] == sector_type)
        )
        dsec = df.loc[mask].copy()
        if len(dsec) < 8:
            continue

        sector_name = f"({sector_S},{sector_support},{sector_type})"
        dsec_reduced = reduce_sector_dataframe(dsec)
        metrics_df, resid_df, paramsC = analyze_sector(
            dsec_reduced,
            sector_name=sector_name,
            use_exp_gamma=use_exp_gamma,
            plot=plot
        )
        all_metrics.append(metrics_df)

        param_row = {"sector": sector_name}
        for i, p in enumerate(paramsC):
            param_row[f"p{i}"] = p
        all_params.append(param_row)

    metrics_all = pd.concat(all_metrics, ignore_index=True) if all_metrics else pd.DataFrame()
    params_all = pd.DataFrame(all_params) if all_params else pd.DataFrame()

    return metrics_all, params_all


def analyze_sector_reduced(
    df_sector,
    sector_name="sector",
    delta_B1=1.0,
    delta_C1=1.0,
    delta_C2=1.0,
    plot=True,
):
    d = reduce_sector_dataframe(df_sector).sort_values("log_mB")

    x = d["log_mB"].to_numpy(dtype=float)
    c = d["C_fam"].to_numpy(dtype=float)
    k = d["kappa13"].to_numpy(dtype=float)
    y = d["log_mult"].to_numpy(dtype=float)

    results = []
    param_rows = []
    residual_report_df = pd.DataFrame()
    residual_summary_df = pd.DataFrame()
    rest_c1_tests_df = pd.DataFrame()
    rest_c1_summary_df = pd.DataFrame()
    spectral_c1_summary_df = pd.DataFrame()
    spectral_c1_peaks_df = pd.DataFrame()

    X_A = model_A_matrix(x, c, k)
    model_A, yhat_A = fit_linear_model(X_A, y)
    loocv_A = loocv_linear(X_A, y)
    results.append(
        collect_metrics(
            sector_name, "A_linear", y, yhat_A, loocv_A, n_params=4
        )
    )

    X_B1 = model_B1_matrix(x, c, k, delta=delta_B1)
    model_B1, yhat_B1 = fit_linear_model(X_B1, y)
    loocv_B1 = loocv_linear(X_B1, y)
    results.append(
        collect_metrics(
            sector_name,
            f"B1_inv_shifted_delta={delta_B1}",
            y,
            yhat_B1,
            loocv_B1,
            n_params=5,
        )
    )

    X_B2 = model_B2_matrix_fixed(x, c, k)
    model_B2, yhat_B2 = fit_linear_model(X_B2, y)
    loocv_B2 = loocv_linear(X_B2, y)
    results.append(
        collect_metrics(
            sector_name,
            "B2_inv_combined_fixed",
            y,
            yhat_B2,
            loocv_B2,
            n_params=5,
        )
    )

    def model_C1_fixed_delta(vars_, a0, a1, a2, a3, beta):
        return model_C1(vars_, a0, a1, a2, a3, beta, delta_C1)

    p0_C1 = [0.0, 1.0, 1.0, 1.0, 0.05]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OptimizeWarning)
        params_C1, _ = curve_fit(
            model_C1_fixed_delta,
            (x, c, k),
            y,
            p0=p0_C1,
            maxfev=100000,
        )
    yhat_C1 = model_C1_fixed_delta((x, c, k), *params_C1)
    loocv_C1, succ_C1, mask_C1 = loocv_nonlinear(
        x, c, k, y, model_C1_fixed_delta, p0=p0_C1, maxfev=100000
    )
    results.append(
        collect_metrics(
            sector_name,
            f"C1_logGamma_xplus_delta={delta_C1}",
            y,
            yhat_C1,
            loocv_C1,
            n_params=5,
            mask=mask_C1,
            extra={"LOOCV_success_rate": succ_C1},
        )
    )
    param_rows.append(
        {
            "sector": sector_name,
            "model": f"C1_logGamma_xplus_delta={delta_C1}",
            "a0": params_C1[0],
            "a1": params_C1[1],
            "a2": params_C1[2],
            "a3": params_C1[3],
            "beta": params_C1[4],
        }
    )

    residual_report_df, residual_summary_df = compare_residuals_A_C1(
        d,
        params_C1=params_C1,
        delta_C1=delta_C1,
        sector_name=sector_name,
        plot=plot,
    )
    R_C1 = y - yhat_C1
    rest_tests_A, rest_summary_A = analyze_rest_C1(
        x, c, k, R_C1, delta=delta_C1
    )
    rest_tests_A["sector"] = sector_name
    rest_summary_A["sector"] = sector_name
    spectral_c1_summary_df, spectral_c1_peaks_df, _ = spectral_preanalysis_rest(
        x,
        R_C1,
        sector_name=sector_name,
        top_k=5,
        detrend_linear=False,
        plot=plot,
    )
    spectral_c1_peaks_df["sector"] = sector_name
    spectral_c1_summary_dt_df, spectral_c1_peaks_dt_df, _ = spectral_preanalysis_rest(
        x,
        R_C1,
        sector_name=sector_name,
        top_k=5,
        detrend_linear=True,
        plot=plot,
    )
    spectral_c1_peaks_dt_df["sector"] = sector_name
    spectral_c1_summary_df = pd.concat(
        [spectral_c1_summary_df, spectral_c1_summary_dt_df], ignore_index=True
    )
    spectral_c1_peaks_df = pd.concat(
        [spectral_c1_peaks_df, spectral_c1_peaks_dt_df], ignore_index=True
    )

    def model_C2_fixed_delta(vars_, a0, a1, a2, a3, beta, lam):
        return model_C2(vars_, a0, a1, a2, a3, beta, lam, delta_C2)

    p0_C2 = [0.0, 1.0, 1.0, 1.0, 0.05, 1.0]
    bounds_C2 = (
        [-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, 1e-6],
        [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OptimizeWarning)
        params_C2, _ = curve_fit(
            model_C2_fixed_delta,
            (x, c, k),
            y,
            p0=p0_C2,
            bounds=bounds_C2,
            maxfev=100000,
        )
    yhat_C2 = model_C2_fixed_delta((x, c, k), *params_C2)
    loocv_C2, succ_C2, mask_C2 = loocv_nonlinear(
        x,
        c,
        k,
        y,
        model_C2_fixed_delta,
        p0=p0_C2,
        bounds=bounds_C2,
        maxfev=100000,
    )
    results.append(
        collect_metrics(
            sector_name,
            f"C2_logGamma_lamxplus_delta={delta_C2}",
            y,
            yhat_C2,
            loocv_C2,
            n_params=6,
            mask=mask_C2,
            extra={"LOOCV_success_rate": succ_C2},
        )
    )
    param_rows.append(
        {
            "sector": sector_name,
            "model": f"C2_logGamma_lamxplus_delta={delta_C2}",
            "a0": params_C2[0],
            "a1": params_C2[1],
            "a2": params_C2[2],
            "a3": params_C2[3],
            "beta": params_C2[4],
            "lam": params_C2[5],
        }
    )

    p0_D1 = [0.0, 1.0, 1.0, 1.0, 0.05, 0.05]
    params_D1, yhat_D1 = fit_D1(x, c, k, y)
    loocv_D1, succ_D1, mask_D1 = loocv_nonlinear(
        x, c, k, y, model_D1, p0=p0_D1, maxfev=100000
    )
    results.append(
        collect_metrics(
            sector_name,
            "D1_dual_logGamma",
            y,
            yhat_D1,
            loocv_D1,
            n_params=6,
            mask=mask_D1,
            extra={"LOOCV_success_rate": succ_D1},
        )
    )
    param_rows.append(
        {
            "sector": sector_name,
            "model": "D1_dual_logGamma",
            "a0": params_D1[0],
            "a1": params_D1[1],
            "a2": params_D1[2],
            "a3": params_D1[3],
            "betaA": params_D1[4],
            "betaB": params_D1[5],
        }
    )

    metrics_df = pd.DataFrame(results)
    params_df = pd.DataFrame(param_rows)

    resid_df = pd.DataFrame(
        {
            "log_mB": x,
            "C_fam": c,
            "kappa13": k,
            "y": y,
            "resid_A": y - yhat_A,
            "resid_B1": y - yhat_B1,
            "resid_B2": y - yhat_B2,
            "resid_C1": y - yhat_C1,
            "resid_C2": y - yhat_C2,
            "resid_D1": y - yhat_D1,
        }
    )

    if plot:
        fig, axes = plt.subplots(3, 2, figsize=(12, 11))

        axes[0, 0].scatter(x, y, s=18, label="Daten")
        axes[0, 0].plot(x, yhat_A, label="A")
        axes[0, 0].plot(x, yhat_B1, label="B1")
        axes[0, 0].plot(x, yhat_B2, label="B2")
        axes[0, 0].plot(x, yhat_C1, label="C1")
        axes[0, 0].plot(x, yhat_C2, label="C2")
        axes[0, 0].plot(x, yhat_D1, label="D1")
        axes[0, 0].set_title(f"{sector_name}: Fits")
        axes[0, 0].set_xlabel("log(m_B)")
        axes[0, 0].set_ylabel("log Multiplizität")
        axes[0, 0].legend()

        axes[0, 1].scatter(x, y - yhat_A, s=18, label="A")
        axes[0, 1].scatter(x, y - yhat_B1, s=18, label="B1")
        axes[0, 1].scatter(x, y - yhat_C1, s=18, label="C1")
        axes[0, 1].axhline(0.0, linestyle="--")
        axes[0, 1].set_title("Residuen vs. log(m_B)")
        axes[0, 1].set_xlabel("log(m_B)")
        axes[0, 1].set_ylabel("Residuum")
        axes[0, 1].legend()

        axes[1, 0].scatter(c, y - yhat_A, s=18)
        axes[1, 0].axhline(0.0, linestyle="--")
        axes[1, 0].set_title("Residuen A vs. C_fam")
        axes[1, 0].set_xlabel("C_fam")
        axes[1, 0].set_ylabel("Residuum")

        axes[1, 1].scatter(k, y - yhat_A, s=18)
        axes[1, 1].axhline(0.0, linestyle="--")
        axes[1, 1].set_title("Residuen A vs. kappa13")
        axes[1, 1].set_xlabel("kappa13")
        axes[1, 1].set_ylabel("Residuum")

        axes[2, 0].scatter(c, y - yhat_C2, s=18)
        axes[2, 0].axhline(0.0, linestyle="--")
        axes[2, 0].set_title("Residuen C2 vs. C_fam")
        axes[2, 0].set_xlabel("C_fam")
        axes[2, 0].set_ylabel("Residuum")

        axes[2, 1].scatter(k, y - yhat_C2, s=18)
        axes[2, 1].axhline(0.0, linestyle="--")
        axes[2, 1].set_title("Residuen C2 vs. kappa13")
        axes[2, 1].set_xlabel("kappa13")
        axes[2, 1].set_ylabel("Residuum")

        plt.tight_layout()
        plt.show()

    return (
        metrics_df,
        params_df,
        resid_df,
        residual_report_df,
        residual_summary_df,
        rest_tests_A,
        rest_summary_A,
        spectral_c1_summary_df,
        spectral_c1_peaks_df,
    )


def compare_sectors_reduced(
    df,
    sectors,
    delta_B1=1.0,
    delta_C1=1.0,
    delta_C2=1.0,
    plot=False,
):
    all_metrics = []
    all_params = []
    all_residual_reports = []
    all_residual_summaries = []
    all_rest_c1_tests = []
    all_rest_c1_summaries = []
    all_spectral_c1_summaries = []
    all_spectral_c1_peaks = []

    for sector_S, sector_support, sector_type in sectors:
        mask = (
            (df["sector_S"] == sector_S)
            & (df["sector_support"] == sector_support)
            & (df["sector_type"] == sector_type)
        )
        dsec = df.loc[mask].copy()

        if len(dsec) < 8:
            continue

        sector_name = f"({sector_S},{sector_support},{sector_type})"
        (
            metrics_df,
            params_df,
            _,
            residual_report_df,
            residual_summary_df,
            rest_c1_tests_df,
            rest_c1_summary_df,
            spectral_c1_summary_df,
            spectral_c1_peaks_df,
        ) = analyze_sector_reduced(
            dsec,
            sector_name=sector_name,
            delta_B1=delta_B1,
            delta_C1=delta_C1,
            delta_C2=delta_C2,
            plot=plot,
        )
        all_metrics.append(metrics_df)
        all_params.append(params_df)
        all_residual_reports.append(residual_report_df)
        all_residual_summaries.append(residual_summary_df)
        all_rest_c1_tests.append(rest_c1_tests_df)
        all_rest_c1_summaries.append(rest_c1_summary_df)
        all_spectral_c1_summaries.append(spectral_c1_summary_df)
        all_spectral_c1_peaks.append(spectral_c1_peaks_df)

    metrics_all = (
        pd.concat(all_metrics, ignore_index=True) if all_metrics else pd.DataFrame()
    )
    params_all = (
        pd.concat(all_params, ignore_index=True) if all_params else pd.DataFrame()
    )
    residual_reports_all = (
        pd.concat(all_residual_reports, ignore_index=True)
        if all_residual_reports
        else pd.DataFrame()
    )
    residual_summaries_all = (
        pd.concat(all_residual_summaries, ignore_index=True)
        if all_residual_summaries
        else pd.DataFrame()
    )
    rest_c1_tests_all = (
        pd.concat(all_rest_c1_tests, ignore_index=True)
        if all_rest_c1_tests
        else pd.DataFrame()
    )
    rest_c1_summaries_all = (
        pd.concat(all_rest_c1_summaries, ignore_index=True)
        if all_rest_c1_summaries
        else pd.DataFrame()
    )
    spectral_c1_summaries_all = (
        pd.concat(all_spectral_c1_summaries, ignore_index=True)
        if all_spectral_c1_summaries
        else pd.DataFrame()
    )
    spectral_c1_peaks_all = (
        pd.concat(all_spectral_c1_peaks, ignore_index=True)
        if all_spectral_c1_peaks
        else pd.DataFrame()
    )
    return (
        metrics_all,
        params_all,
        residual_reports_all,
        residual_summaries_all,
        rest_c1_tests_all,
        rest_c1_summaries_all,
        spectral_c1_summaries_all,
        spectral_c1_peaks_all,
    )


def resolve_input_path():
    if len(sys.argv) >= 2 and not sys.argv[1].startswith("--"):
        return Path(sys.argv[1])

    for candidate in DEFAULT_INPUT_CANDIDATES:
        path = Path(candidate)
        if path.exists():
            return path

    return None


def should_plot():
    return "--plot" in sys.argv[1:]

def main():
    input_path = resolve_input_path()
    if input_path is None:
        print("Keine Eingabedaten gefunden.")
        print("Aufruf: python3 \"Napoleon.py\" <daten.csv> [--plot]")
        print("Erwartete Spalten:")
        print(" - " + ", ".join(REQUIRED_COLUMNS))
        print("Alternativ wird automatisch eine arithmetische Struktur-CSV erkannt.")
        return

    df_raw = pd.read_csv(input_path)
    generated_input_path = None

    if all(col in df_raw.columns for col in REQUIRED_COLUMNS):
        df = df_raw
        print(f"Lade Napoleon-Eingabedaten aus: {input_path}")
    elif is_arith_structure_dataframe(df_raw):
        df = build_napoleon_input_from_arith(df_raw)
        generated_input_path = input_path.with_name(
            f"{input_path.stem}__napoleon_input.csv"
        )
        df.to_csv(generated_input_path, index=False)
        print(f"Lade arithmetische Strukturdaten aus: {input_path}")
        print("Erzeuge daraus automatisch Napoleon-Eingabedaten.")
        print(f"Abgeleitete Eingabedatei: {generated_input_path}")
    else:
        raise ValueError(
            "Die CSV passt weder zum Napoleon-Format noch zum arithmetischen "
            "Strukturformat."
        )

    validate_dataframe(df)

    sectors = discover_sectors(df)
    if not sectors:
        print("Keine auswertbaren Sektoren gefunden.")
        print("Benötigt werden mindestens 8 Zeilen pro Sektor.")
        return

    print(f"Gefundene Sektoren für die Analyse: {len(sectors)}")
    for sector in sectors:
        print(f" - {sector}")

    (
        metrics_all,
        params_all,
        residual_reports_all,
        residual_summaries_all,
        rest_c1_tests_all,
        rest_c1_summaries_all,
        spectral_c1_summaries_all,
        spectral_c1_peaks_all,
    ) = compare_sectors_reduced(
        df,
        sectors,
        delta_B1=1.0,
        delta_C1=1.0,
        delta_C2=1.0,
        plot=should_plot(),
    )

    if metrics_all.empty:
        print("Keine Modellmetriken erzeugt.")
    else:
        print("\nReduzierte Modellmetriken:\n")
        print(metrics_all.sort_values(["sector", "model"]).to_string(index=False))

    if params_all.empty:
        print("\nKeine Gamma-Parameter erzeugt.")
    else:
        print("\nReduzierte Gamma-Parameter:\n")
        print(params_all.to_string(index=False))

    if residual_reports_all.empty:
        print("\nKeine Residuen-Strukturberichte erzeugt.")
    else:
        print("\nResiduen-Strukturberichte A vs. C1:\n")
        print(
            residual_reports_all.sort_values(["sector", "label"]).to_string(index=False)
        )

    if residual_summaries_all.empty:
        print("\nKeine Residuen-Zusammenfassungen erzeugt.")
    else:
        print("\nResiduen-Zusammenfassung A vs. C1:\n")
        print(
            residual_summaries_all.sort_values("sector").to_string(index=False)
        )

    if rest_c1_tests_all.empty:
        print("\nKeine C1-Resttests erzeugt.")
    else:
        print("\nC1-Resttests:\n")
        print(
            rest_c1_tests_all.sort_values(["sector", "test"]).to_string(index=False)
        )

    if rest_c1_summaries_all.empty:
        print("\nKeine C1-Rest-Zusammenfassungen erzeugt.")
    else:
        print("\nC1-Rest-Zusammenfassung:\n")
        print(
            rest_c1_summaries_all.sort_values("sector").to_string(index=False)
        )

    if spectral_c1_summaries_all.empty:
        print("\nKeine spektralen C1-Voranalysen erzeugt.")
    else:
        print("\nSpektrale C1-Voranalyse:\n")
        print(
            spectral_c1_summaries_all.sort_values("sector").to_string(index=False)
        )

    if spectral_c1_peaks_all.empty:
        print("\nKeine spektralen C1-Peaks erzeugt.")
    else:
        print("\nSpektrale C1-Peaks:\n")
        print(
            spectral_c1_peaks_all.sort_values(["sector", "rank"]).to_string(index=False)
        )

    spectral_compare_summary_df, spectral_compare_peaks_df = compare_spectral_variants(
        spectral_c1_summaries_all, spectral_c1_peaks_all
    )

    if spectral_compare_summary_df.empty:
        print("\nKeine spektrale Vergleichszusammenfassung erzeugt.")
    else:
        print("\nVergleich Spektralvarianten:\n")
        print(
            spectral_compare_summary_df.sort_values("sector").to_string(index=False)
        )

    if spectral_compare_peaks_df.empty:
        print("\nKeine spektralen Peak-Differenzen erzeugt.")
    else:
        print("\nPeak-Differenzen detrended vs. roh:\n")
        print(
            spectral_compare_peaks_df.sort_values(["sector", "rank"]).to_string(index=False)
        )

    if generated_input_path is not None:
        print(f"\nGespeicherte Napoleon-Eingabe: {generated_input_path}")


if __name__ == "__main__":
    main()