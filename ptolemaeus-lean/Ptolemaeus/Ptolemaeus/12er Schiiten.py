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
# 2. Kontrollmodelle für Vierergruppen
# ============================================================

def permute_within_quadruplets(quadruplets: np.ndarray, seed: int = 1234):
    rng = np.random.default_rng(seed)
    out = quadruplets.copy()
    for k in range(len(out)):
        out[k] = rng.permutation(out[k])
    return out


def random_like_quadruplets(quadruplets: np.ndarray, seed: int = 1234):
    """
    Erzeugt zufällige 4er-Gruppen ähnlicher Größenordnung wie die echten Vierlinge.
    Dazu wird pro Vierling ein Zentrum ungefähr beibehalten und zufällige Offsets genommen.
    Keine Primzahlbedingung.
    """
    rng = np.random.default_rng(seed)
    centers = quadruplets.mean(axis=1)
    out = []
    for c in centers:
        # zufällige Offsets in ähnlicher lokaler Größenordnung
        offsets = np.sort(rng.integers(-6, 7, size=4))
        q = np.maximum(2, np.round(c + offsets)).astype(float)
        out.append(q)
    return np.array(out, dtype=float)


# ============================================================
# 3. Feature-Bau
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
# 4. Blockmatrix
# ============================================================

def build_block_matrix(
    quadruplets: np.ndarray,
    alpha: float = 2.0,
    beta: float = 0.08,
    eta: float = 0.5,
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
        # zerstörtes 2-4-2-Muster
        u12 = 1.0 + rho * 0.25
        u23 = 1.0 + rho * 0.25
        u34 = 1.0 + rho * 0.25

    for k in range(N):
        idx = 4 * k

        phase_term = math.cos(2.0 * math.pi * theta[k])
        diag = alpha * delta[k] + beta * R_norm[k] + eta * phase_term

        B = np.array([
            [diag[0], u12,    0.0,   0.0],
            [u12,    diag[1], u23,   0.0],
            [0.0,    u23,    diag[2], u34],
            [0.0,    0.0,    u34,   diag[3]],
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
# 5. Auswertung
# ============================================================

def affine_fit(x: np.ndarray, y: np.ndarray):
    A = np.vstack([x, np.ones_like(x)]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(a), float(b)


def evaluate_model(
    quadruplets: np.ndarray,
    n_zeros: int = 100,
    alpha: float = 2.0,
    beta: float = 0.08,
    eta: float = 0.5,
    rho: float = 0.8,
    gamma0: float = 0.6,
    gamma1: float = 0.08,
    gamma2: float = 0.15,
    use_242: bool = True,
):
    M, features = build_block_matrix(
        quadruplets=quadruplets,
        alpha=alpha,
        beta=beta,
        eta=eta,
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

    a, b = affine_fit(lam, gam)
    pred = a * lam + b

    rmse = float(np.sqrt(np.mean((pred - gam) ** 2)))
    mae = float(np.mean(np.abs(pred - gam)))

    d_pred = np.diff(pred)
    d_gam = np.diff(gam)

    if np.std(d_pred) < 1e-12 or np.std(d_gam) < 1e-12:
        spacing_corr = float("nan")
    else:
        spacing_corr = float(np.corrcoef(d_pred, d_gam)[0, 1])

    return {
        "features": features,
        "eigvals": eigvals,
        "positive_eigvals_used": lam,
        "riemann_gammas": gam,
        "fit_a": a,
        "fit_b": b,
        "predicted_gammas": pred,
        "rmse": rmse,
        "mae": mae,
        "spacing_corr": spacing_corr,
    }


# ============================================================
# 6. Benchmark über mehrere Modelle
# ============================================================

def run_benchmarks(n_quadruplets: int = 100, n_zeros: int = 100, seed: int = 1234):
    base = first_prime_quadruplets(n_quadruplets)
    perm = permute_within_quadruplets(base, seed=seed)
    rnd  = random_like_quadruplets(base, seed=seed)

    configs = [
        {
            "name": "original",
            "quadruplets": base,
            "pars": dict(alpha=2.0, beta=0.08, eta=0.5, rho=0.8, gamma0=0.6, gamma1=0.08, gamma2=0.15, use_242=True),
        },
        {
            "name": "permutiert",
            "quadruplets": perm,
            "pars": dict(alpha=2.0, beta=0.08, eta=0.5, rho=0.8, gamma0=0.6, gamma1=0.08, gamma2=0.15, use_242=True),
        },
        {
            "name": "random_like",
            "quadruplets": rnd,
            "pars": dict(alpha=2.0, beta=0.08, eta=0.5, rho=0.8, gamma0=0.6, gamma1=0.08, gamma2=0.15, use_242=True),
        },
        {
            "name": "destroy_242",
            "quadruplets": base,
            "pars": dict(alpha=2.0, beta=0.08, eta=0.5, rho=0.8, gamma0=0.6, gamma1=0.08, gamma2=0.15, use_242=False),
        },
        {
            "name": "ohne_achterrahmen",
            "quadruplets": base,
            "pars": dict(alpha=2.0, beta=0.0, eta=0.0, rho=0.8, gamma0=0.6, gamma1=0.0, gamma2=0.0, use_242=True),
        },
    ]

    rows = []
    detailed = {}

    for cfg in configs:
        res = evaluate_model(
            quadruplets=cfg["quadruplets"],
            n_zeros=n_zeros,
            **cfg["pars"],
        )
        rows.append({
            "modell": cfg["name"],
            "rmse": res["rmse"],
            "mae": res["mae"],
            "spacing_corr": res["spacing_corr"],
            "fit_a": res["fit_a"],
            "fit_b": res["fit_b"],
        })
        detailed[cfg["name"]] = res

    df = pd.DataFrame(rows).sort_values(["rmse", "mae"], ascending=True).reset_index(drop=True)
    return df, detailed


# ============================================================
# 7. Ausgabe
# ============================================================

def print_model_head(name: str, res: dict, n: int = 8):
    print("=" * 72)
    print(f"MODELL: {name}")
    print("=" * 72)
    print(f"RMSE         : {res['rmse']:.6f}")
    print(f"MAE          : {res['mae']:.6f}")
    print(f"Spacing Corr : {res['spacing_corr']:.6f}")
    print(f"Affine Fit   : gamma ≈ {res['fit_a']:.6f} * lambda + {res['fit_b']:.6f}")
    print()
    print(" n |   lambda_n   | affine(lambda_n) |   gamma_n")
    print("-" * 58)
    lam = res["positive_eigvals_used"]
    pred = res["predicted_gammas"]
    gam = res["riemann_gammas"]
    for i in range(min(n, len(lam))):
        print(f"{i+1:2d} | {lam[i]:11.6f} | {pred[i]:16.6f} | {gam[i]:10.6f}")


if __name__ == "__main__":
    df, detailed = run_benchmarks(n_quadruplets=100, n_zeros=100, seed=1234)

    print("\nBenchmark-Tabelle:\n")
    print(df.to_string(index=False))

    print("\n")
    for name in df["modell"]:
        print_model_head(name, detailed[name], n=8)
        print()