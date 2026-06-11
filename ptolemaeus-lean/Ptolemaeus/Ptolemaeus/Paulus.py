import argparse
import math
import numpy as np
import sympy as sp
import pandas as pd
from mpmath import zetazero
from scipy.linalg import eigh


# ============================================================
# 1. Grunddaten
# ============================================================

def first_prime_quadruplets(n: int):
    quadruplets = []
    p = 5
    while len(quadruplets) < n:
        if sp.isprime(p) and sp.isprime(p + 2) and sp.isprime(p + 6) and sp.isprime(p + 8):
            quadruplets.append((p, p + 2, p + 6, p + 8))
        p += 2
    return np.array(quadruplets, dtype=float)


def zscore(x: np.ndarray, eps: float = 1e-12):
    mu = np.mean(x)
    sigma = np.std(x)
    return (x - mu) / (sigma + eps)


def first_riemann_gammas(n: int):
    return np.array([float(zetazero(k).imag) for k in range(1, n + 1)], dtype=float)


# ============================================================
# 2. Kontrollmodelle
# ============================================================

def permute_within_quadruplets(quadruplets: np.ndarray, seed: int = 1234):
    rng = np.random.default_rng(seed)
    out = quadruplets.copy()
    for k in range(len(out)):
        out[k] = rng.permutation(out[k])
    return out


def random_like_quadruplets(quadruplets: np.ndarray, seed: int = 1234):
    rng = np.random.default_rng(seed)
    centers = quadruplets.mean(axis=1)
    out = []
    for c in centers:
        offsets = np.sort(rng.integers(-6, 7, size=4))
        q = np.maximum(2, np.round(c + offsets)).astype(float)
        out.append(q)
    return np.array(out, dtype=float)


# ============================================================
# 3. Features
# ============================================================

def build_features(quadruplets: np.ndarray):
    logs = np.log(quadruplets)
    mu = logs.mean(axis=1, keepdims=True)
    delta = logs - mu

    centers = quadruplets.mean(axis=1)
    R = np.floor(np.log2(centers / 8.0)).astype(int)
    theta = centers / (2.0 ** (R + 3)) - 1.0

    if not np.all((theta >= -1e-12) & (theta < 1 + 1e-12)):
        raise ValueError("Einige theta_k liegen nicht in [0,1).")

    R_norm = zscore(R.astype(float))

    return {
        "delta": delta,
        "centers": centers,
        "R": R,
        "R_norm": R_norm,
        "theta": theta,
    }


# ============================================================
# 4. Matrixmodell mit Kubikterm und 3er-Phase
# ============================================================

def build_block_matrix(
    quadruplets: np.ndarray,
    alpha: float = 2.0,
    alpha3: float = 0.0,      # NEU: Kubikterm in delta
    beta: float = 0.08,
    eta1: float = 0.5,        # cos(2π theta)
    eta3: float = 0.0,        # NEU: cos(6π theta)
    rho: float = 0.8,
    gamma0: float = 0.6,
    gamma1: float = 0.08,
    gamma2: float = 0.15,
    use_242: bool = True,
):
    features = build_features(quadruplets)
    delta = features["delta"]
    R_norm = features["R_norm"]
    theta = features["theta"]

    N = len(quadruplets)
    dim = 4 * N
    M = np.zeros((dim, dim), dtype=float)

    if use_242:
        u12 = 1.0 + rho * (2.0 / 8.0)
        u23 = 1.0 + rho * (4.0 / 8.0)
        u34 = 1.0 + rho * (2.0 / 8.0)
    else:
        u12 = 1.0 + rho * 0.25
        u23 = 1.0 + rho * 0.25
        u34 = 1.0 + rho * 0.25

    for k in range(N):
        idx = 4 * k

        phase1 = math.cos(2.0 * math.pi * theta[k])
        phase3 = math.cos(6.0 * math.pi * theta[k])

        diag = (
            alpha * delta[k]
            + alpha3 * (delta[k] ** 3)
            + beta * R_norm[k]
            + eta1 * phase1
            + eta3 * phase3
        )

        B = np.array([
            [diag[0], u12,    0.0,    0.0],
            [u12,    diag[1], u23,    0.0],
            [0.0,    u23,    diag[2], u34],
            [0.0,    0.0,    u34,    diag[3]],
        ], dtype=float)

        M[idx:idx+4, idx:idx+4] = B

        if k < N - 1:
            gk = (
                gamma0
                + gamma1 * (R_norm[k + 1] - R_norm[k])
                + gamma2 * math.cos(2.0 * math.pi * (theta[k + 1] - theta[k]))
            )
            C = gk * np.eye(4)
            M[idx:idx+4, idx+4:idx+8] = C
            M[idx+4:idx+8, idx:idx+4] = C.T

    return M, features


# ============================================================
# 5. Eichungen
# ============================================================

FIT1_A = 149.654970
FIT1_B = 35.538153

FIT2_A = 273.234031
FIT2_B = 14.302539


def affine_fit(x: np.ndarray, y: np.ndarray):
    A = np.vstack([x, np.ones_like(x)]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(a), float(b)


def score_prediction(pred: np.ndarray, target: np.ndarray):
    rmse = float(np.sqrt(np.mean((pred - target) ** 2)))
    mae = float(np.mean(np.abs(pred - target)))

    d_pred = np.diff(pred)
    d_tar = np.diff(target)

    if np.std(d_pred) < 1e-12 or np.std(d_tar) < 1e-12:
        spacing_corr = float("nan")
    else:
        spacing_corr = float(np.corrcoef(d_pred, d_tar)[0, 1])

    return rmse, mae, spacing_corr


# ============================================================
# 6. Auswertung eines Modells
# ============================================================

def evaluate_model(
    quadruplets: np.ndarray,
    n_zeros: int = 100,
    alpha: float = 2.0,
    alpha3: float = 0.0,
    beta: float = 0.08,
    eta1: float = 0.5,
    eta3: float = 0.0,
    rho: float = 0.8,
    gamma0: float = 0.6,
    gamma1: float = 0.08,
    gamma2: float = 0.15,
    use_242: bool = True,
):
    M, features = build_block_matrix(
        quadruplets=quadruplets,
        alpha=alpha,
        alpha3=alpha3,
        beta=beta,
        eta1=eta1,
        eta3=eta3,
        rho=rho,
        gamma0=gamma0,
        gamma1=gamma1,
        gamma2=gamma2,
        use_242=use_242,
    )

    eigvals = eigh(M, eigvals_only=True)
    eigvals = np.sort(eigvals)
    pos = eigvals[eigvals > 1e-10]

    if len(pos) < n_zeros:
        raise ValueError(f"Zu wenige positive Eigenwerte: {len(pos)} < {n_zeros}")

    lam = pos[:n_zeros]
    gam = first_riemann_gammas(n_zeros)

    # freie affine Eichung
    a_free, b_free = affine_fit(lam, gam)
    pred_free = a_free * lam + b_free
    rmse_free, mae_free, corr_free = score_prediction(pred_free, gam)

    # externer Fit 1
    pred_fit1 = FIT1_A * lam + FIT1_B
    rmse_fit1, mae_fit1, corr_fit1 = score_prediction(pred_fit1, gam)

    # externer Fit 2
    pred_fit2 = FIT2_A * lam + FIT2_B
    rmse_fit2, mae_fit2, corr_fit2 = score_prediction(pred_fit2, gam)

    return {
        "features": features,
        "eigvals": eigvals,
        "lambda_used": lam,
        "gammas": gam,

        "free_a": a_free,
        "free_b": b_free,
        "pred_free": pred_free,
        "rmse_free": rmse_free,
        "mae_free": mae_free,
        "corr_free": corr_free,

        "pred_fit1": pred_fit1,
        "rmse_fit1": rmse_fit1,
        "mae_fit1": mae_fit1,
        "corr_fit1": corr_fit1,

        "pred_fit2": pred_fit2,
        "rmse_fit2": rmse_fit2,
        "mae_fit2": mae_fit2,
        "corr_fit2": corr_fit2,
    }


# ============================================================
# 7. Datensätze
# ============================================================

def build_dataset(kind: str, base_quadruplets: np.ndarray, seed: int):
    if kind == "original":
        return base_quadruplets
    elif kind == "permutiert":
        return permute_within_quadruplets(base_quadruplets, seed=seed)
    elif kind == "random_like":
        return random_like_quadruplets(base_quadruplets, seed=seed)
    else:
        raise ValueError(f"Unbekannter Datensatztyp: {kind}")


# ============================================================
# 8. Modellfamilien M0-M3
# ============================================================

def model_family_params(model_name: str):
    """
    M0: ohne Kubik, ohne 3er-Phase
    M1: nur Kubikterm
    M2: nur 3er-Phase
    M3: beides
    """
    base = dict(
        alpha=2.0,
        beta=0.08,
        eta1=0.5,
        rho=0.8,
        gamma0=0.6,
        gamma1=0.08,
        gamma2=0.15,
    )

    if model_name == "M0":
        return dict(base, alpha3=0.0, eta3=0.0)
    elif model_name == "M1":
        return dict(base, alpha3=1.2, eta3=0.0)
    elif model_name == "M2":
        return dict(base, alpha3=0.0, eta3=0.25)
    elif model_name == "M3":
        return dict(base, alpha3=1.2, eta3=0.25)
    else:
        raise ValueError(f"Unbekannte Modellfamilie: {model_name}")


# ============================================================
# 9. Robuster Vergleich
# ============================================================

def run_family_scan(
    N_list=(50, 100, 150),
    seeds=(1234, 2024, 7777),
):
    rows = []

    maxN = max(N_list)
    base_all = first_prime_quadruplets(maxN)

    n_steps = len(N_list) * len(seeds) * 3 * 2 * 4
    step_i = 0

    for N in N_list:
        base = base_all[:N]
        n_zeros = min(N, 100)

        for seed in seeds:
            for dataset_kind in ("original", "permutiert", "random_like"):
                quadruplets = build_dataset(dataset_kind, base, seed)

                for use_242 in (True, False):
                    for model_name in ("M0", "M1", "M2", "M3"):
                        step_i += 1
                        if step_i == 1 or step_i % 25 == 0 or step_i == n_steps:
                            print(
                                f"[Paulus] Scan {step_i}/{n_steps} (N={N}, {dataset_kind}, {model_name})",
                                flush=True,
                            )

                        pars = model_family_params(model_name)

                        row_base = {
                            "N": N,
                            "seed": seed,
                            "dataset": dataset_kind,
                            "use_242": use_242,
                            "model": model_name,
                            "alpha": pars["alpha"],
                            "alpha3": pars["alpha3"],
                            "beta": pars["beta"],
                            "eta1": pars["eta1"],
                            "eta3": pars["eta3"],
                            "rho": pars["rho"],
                            "gamma0": pars["gamma0"],
                            "gamma1": pars["gamma1"],
                            "gamma2": pars["gamma2"],
                        }

                        try:
                            res = evaluate_model(
                                quadruplets=quadruplets,
                                n_zeros=n_zeros,
                                use_242=use_242,
                                **pars,
                            )

                            rows.append({
                                **row_base,
                                "rmse_free": res["rmse_free"],
                                "mae_free": res["mae_free"],
                                "corr_free": res["corr_free"],
                                "rmse_fit1": res["rmse_fit1"],
                                "mae_fit1": res["mae_fit1"],
                                "corr_fit1": res["corr_fit1"],
                                "rmse_fit2": res["rmse_fit2"],
                                "mae_fit2": res["mae_fit2"],
                                "corr_fit2": res["corr_fit2"],
                                "free_a": res["free_a"],
                                "free_b": res["free_b"],
                                "error": None,
                            })

                        except Exception as e:
                            rows.append({
                                **row_base,
                                "rmse_free": np.nan,
                                "mae_free": np.nan,
                                "corr_free": np.nan,
                                "rmse_fit1": np.nan,
                                "mae_fit1": np.nan,
                                "corr_fit1": np.nan,
                                "rmse_fit2": np.nan,
                                "mae_fit2": np.nan,
                                "corr_fit2": np.nan,
                                "free_a": np.nan,
                                "free_b": np.nan,
                                "error": str(e),
                            })

    return pd.DataFrame(rows)


# ============================================================
# 10. Zusammenfassungen
# ============================================================

def summarize_family_scan(df: pd.DataFrame):
    if "error" not in df.columns:
        ok = df.copy()
    else:
        ok = df[df["error"].isna()].copy()

    if ok.empty:
        empty = pd.DataFrame()
        return empty, empty, empty, empty, empty, empty

    summary_free = (
        ok.groupby(["N", "dataset", "use_242", "model"])[["rmse_free", "mae_free", "corr_free"]]
        .mean()
        .reset_index()
        .sort_values(["N", "rmse_free", "mae_free"])
    )

    summary_fit1 = (
        ok.groupby(["N", "dataset", "use_242", "model"])[["rmse_fit1", "mae_fit1", "corr_fit1"]]
        .mean()
        .reset_index()
        .sort_values(["N", "rmse_fit1", "mae_fit1"])
    )

    summary_fit2 = (
        ok.groupby(["N", "dataset", "use_242", "model"])[["rmse_fit2", "mae_fit2", "corr_fit2"]]
        .mean()
        .reset_index()
        .sort_values(["N", "rmse_fit2", "mae_fit2"])
    )

    best_free = ok.sort_values("rmse_free").head(30)
    best_fit1 = ok.sort_values("rmse_fit1").head(30)
    best_fit2 = ok.sort_values("rmse_fit2").head(30)

    return summary_free, summary_fit1, summary_fit2, best_free, best_fit1, best_fit2


# ============================================================
# 11. Hauptlauf
# ============================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Paulus: robuster Vergleich der Modellfamilien M0–M3.",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help=(
            "Voller Scan: N ∈ {50,100,150}, drei Seeds — kann mehrere Minuten dauern. "
            "Ohne dieses Flag: Kurzscan (N ∈ {30,50}, ein Seed), nur zum Schnelltest."
        ),
    )
    args = parser.parse_args()

    if args.full:
        N_list = (50, 100, 150)
        seeds = (1234, 2024, 7777)
    else:
        N_list = (30, 50)
        seeds = (1234,)

    print(
        f"[Paulus] Modus: {'VOLL (--full)' if args.full else 'KURZ (Standard)'} — "
        f"N_list={N_list}, seeds={seeds}",
        flush=True,
    )
    print("[Paulus] Starte run_family_scan …", flush=True)
    df = run_family_scan(N_list=N_list, seeds=seeds)
    n_err = int(df["error"].notna().sum()) if "error" in df.columns else 0
    if n_err:
        print(f"\n[Paulus] Hinweis: {n_err} Läufe mit Fehler (siehe Spalte error).\n", flush=True)

    print("\nErste Zeilen der Gesamttabelle:\n")
    print(df.head(20).to_string(index=False))

    summary_free, summary_fit1, summary_fit2, best_free, best_fit1, best_fit2 = summarize_family_scan(df)

    print("\n\n===== FREIE AFFINE EICHUNG =====\n")
    print(summary_free.to_string(index=False))

    print("\n\n===== FESTER FIT 1 =====\n")
    print(summary_fit1.to_string(index=False))

    print("\n\n===== FESTER FIT 2 =====\n")
    print(summary_fit2.to_string(index=False))

    print("\n\n===== BESTE 30 (frei) =====\n")
    print(best_free.to_string(index=False))

    print("\n\n===== BESTE 30 (Fit 1) =====\n")
    print(best_fit1.to_string(index=False))

    print("\n\n===== BESTE 30 (Fit 2) =====\n")
    print(best_fit2.to_string(index=False))