import math
import itertools
import numpy as np
import sympy as sp
import pandas as pd
from mpmath import zetazero
from scipy.linalg import eigh


# ============================================================
# 1. Daten
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


def first_riemann_gammas(n: int, path="zeros6.npy"):
    arr = np.load(path)
    return np.array(arr[:n], dtype=float)


# ============================================================
# 2. Features
# ============================================================

def build_features(quadruplets: np.ndarray):
    q = quadruplets
    logs = np.log(q)
    mu = logs.mean(axis=1, keepdims=True)
    delta = logs - mu

    centers = q.mean(axis=1)
    R = np.floor(np.log2(centers / 8.0)).astype(int)
    theta = centers / (2.0 ** (R + 3)) - 1.0

    # inverse S^3-Randdichte
    D = 4.0 * np.sum(q ** 0.75, axis=1) / np.sum(q, axis=1)
    invD = 1.0 / D

    # ternäre Tripelgröße
    t123 = (q[:, 0] * q[:, 1] * q[:, 2]) ** (1.0 / 3.0)
    t234 = (q[:, 1] * q[:, 2] * q[:, 3]) ** (1.0 / 3.0)
    chi = np.log((t123 + t234) / 2.0) - np.log(centers)

    return {
        "delta": delta,
        "centers": centers,
        "R": R,
        "R_norm": zscore(R.astype(float)),
        "theta": theta,
        "invD": invD,
        "invD_norm": zscore(invD),
        "chi": chi,
        "chi_norm": zscore(chi),
        "invDchi": invD * chi,
        "invDchi_norm": zscore(invD * chi),
    }


# ============================================================
# 3. Blockmatrix
# ============================================================

def build_block_matrix(
    quadruplets: np.ndarray,
    alpha: float = 2.0,
    alpha3: float = 1.2,
    beta: float = 0.08,
    eta1: float = 0.5,
    eta3: float = 0.25,
    rho: float = 0.8,
    gamma0: float = 0.6,
    gamma1: float = 0.08,
    gamma2: float = 0.15,
    zeta1: float = 0.0,   # invD
    zeta2: float = 0.0,   # chi
    zeta3: float = 0.0,   # invD * chi
    use_242: bool = True,
):
    feat = build_features(quadruplets)

    delta = feat["delta"]
    R_norm = feat["R_norm"]
    theta = feat["theta"]
    invD_norm = feat["invD_norm"]
    chi_norm = feat["chi_norm"]
    invDchi_norm = feat["invDchi_norm"]

    N = len(quadruplets)
    dim = 4 * N
    M = np.zeros((dim, dim), dtype=float)

    if use_242:
        u12 = 1.0 + rho * (2.0 / 8.0)
        u23 = 1.0 + rho * (4.0 / 8.0)
        u34 = 1.0 + rho * (2.0 / 8.0)
    else:
        u12 = u23 = u34 = 1.0 + rho * 0.25

    for k in range(N):
        idx = 4 * k

        phase1 = math.cos(2.0 * math.pi * theta[k])
        phase3 = math.cos(6.0 * math.pi * theta[k])

        extra = (
            zeta1 * invD_norm[k]
            + zeta2 * chi_norm[k]
            + zeta3 * invDchi_norm[k]
        )

        diag = (
            alpha * delta[k]
            + alpha3 * (delta[k] ** 3)
            + beta * R_norm[k]
            + eta1 * phase1
            + eta3 * phase3
            + extra
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

    return M


# ============================================================
# 4. Auswertung
# ============================================================

def affine_fit(x: np.ndarray, y: np.ndarray):
    A = np.vstack([x, np.ones_like(x)]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(a), float(b)


def evaluate_model(
    quadruplets: np.ndarray,
    n_zeros: int = 100,
    zeros_path: str = "zeros6.npy",
    **pars,
):
    M = build_block_matrix(quadruplets, **pars)

    eigvals = np.sort(eigh(M, eigvals_only=True))
    pos = eigvals[eigvals > 1e-10]
    lam = pos[:n_zeros]

    gam = first_riemann_gammas(n_zeros, path=zeros_path)

    a, b = affine_fit(lam, gam)
    pred = a * lam + b

    rmse = float(np.sqrt(np.mean((pred - gam) ** 2)))
    mae = float(np.mean(np.abs(pred - gam)))
    spacing_corr = float(np.corrcoef(np.diff(pred), np.diff(gam))[0, 1])

    return {
        "rmse": rmse,
        "mae": mae,
        "spacing_corr": spacing_corr,
        "a": a,
        "b": b,
    }


# ============================================================
# 5. Getrennter Scan
# ============================================================

def run_split_scan():
    quadruplets = first_prime_quadruplets(100)

    base = dict(
        alpha=2.0,
        alpha3=1.2,
        beta=0.08,
        eta1=0.5,
        eta3=0.25,
        rho=0.8,
        gamma0=0.6,
        gamma1=0.08,
        gamma2=0.15,
        use_242=True,
    )

    grid1 = [0.0, 0.1, 0.25, 0.5]   # invD
    grid2 = [0.0, 0.1, 0.25, 0.5]   # chi
    grid3 = [0.0, 0.1, 0.25, 0.5]   # invDchi

    rows = []

    # A) nur invD
    for z1 in grid1:
        pars = dict(base, zeta1=z1, zeta2=0.0, zeta3=0.0)
        res = evaluate_model(quadruplets, **pars)
        rows.append({
            "family": "invD_only",
            "zeta1": z1, "zeta2": 0.0, "zeta3": 0.0,
            **res
        })

    # B) nur chi
    for z2 in grid2:
        pars = dict(base, zeta1=0.0, zeta2=z2, zeta3=0.0)
        res = evaluate_model(quadruplets, **pars)
        rows.append({
            "family": "chi_only",
            "zeta1": 0.0, "zeta2": z2, "zeta3": 0.0,
            **res
        })

    # C) nur invDchi
    for z3 in grid3:
        pars = dict(base, zeta1=0.0, zeta2=0.0, zeta3=z3)
        res = evaluate_model(quadruplets, **pars)
        rows.append({
            "family": "invDchi_only",
            "zeta1": 0.0, "zeta2": 0.0, "zeta3": z3,
            **res
        })

    # D) kleiner gemeinsamer Scan
    combos = [
        (0.10, 0.10, 0.00),
        (0.25, 0.10, 0.00),
        (0.25, 0.25, 0.00),
        (0.25, 0.25, 0.10),
        (0.50, 0.25, 0.10),
        (0.25, 0.50, 0.10),
        (0.50, 0.50, 0.25),
    ]

    for z1, z2, z3 in combos:
        pars = dict(base, zeta1=z1, zeta2=z2, zeta3=z3)
        res = evaluate_model(quadruplets, **pars)
        rows.append({
            "family": "combined",
            "zeta1": z1, "zeta2": z2, "zeta3": z3,
            **res
        })

    df = pd.DataFrame(rows).sort_values(["rmse", "mae"])
    return df


if __name__ == "__main__":
    df = run_split_scan()
    print(df.to_string(index=False))