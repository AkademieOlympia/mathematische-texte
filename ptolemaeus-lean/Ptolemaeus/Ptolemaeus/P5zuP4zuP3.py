import math
import numpy as np
import pandas as pd
import sympy as sp
from scipy.linalg import eigh


# ============================================================
# 1. Primzahlkonstellationen
# ============================================================

def first_prime_quadruplets(n: int):
    """
    Erste n Primzahlvierlinge der Form
    (p, p+2, p+6, p+8)
    """
    quadruplets = []
    p = 5
    while len(quadruplets) < n:
        if sp.isprime(p) and sp.isprime(p + 2) and sp.isprime(p + 6) and sp.isprime(p + 8):
            quadruplets.append((p, p + 2, p + 6, p + 8))
        p += 2
    return np.array(quadruplets, dtype=float)


def first_prime_quintuplets_pattern_026812(n: int):
    """
    Erste n Primzahl-Fünflinge der Form
    (p, p+2, p+6, p+8, p+12)
    """
    quints = []
    p = 5
    while len(quints) < n:
        vals = [p, p + 2, p + 6, p + 8, p + 12]
        if all(sp.isprime(v) for v in vals):
            quints.append(tuple(vals))
        p += 2
    return np.array(quints, dtype=float)


def first_prime_quintuplets_pattern_0461012(n: int):
    """
    Erste n Primzahl-Fünflinge der Form
    (p, p+4, p+6, p+10, p+12)
    """
    quints = []
    p = 5
    while len(quints) < n:
        vals = [p, p + 4, p + 6, p + 10, p + 12]
        if all(sp.isprime(v) for v in vals):
            quints.append(tuple(vals))
        p += 2
    return np.array(quints, dtype=float)


# ============================================================
# 2. Hilfsfunktionen
# ============================================================

def zscore(x: np.ndarray, eps: float = 1e-12):
    mu = np.mean(x)
    sigma = np.std(x)
    return (x - mu) / (sigma + eps)


def load_riemann_gammas(n: int, path="zeros6.npy"):
    arr = np.load(path)
    return np.array(arr[:n], dtype=float)


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
# 3. Vierlings-Features
# ============================================================

def build_quadruplet_features(quadruplets: np.ndarray):
    q = quadruplets
    logs = np.log(q)
    mu = logs.mean(axis=1, keepdims=True)
    delta = logs - mu

    centers = q.mean(axis=1)

    # dyadischer Achterrahmen
    R = np.floor(np.log2(centers / 8.0)).astype(int)
    theta = centers / (2.0 ** (R + 3)) - 1.0

    # S^3-Rest aus dem 4D-Modell
    # D_sigma = 4 * sum(q^(3/4)) / sum(q)
    D_sigma = 4.0 * np.sum(q ** 0.75, axis=1) / np.sum(q, axis=1)
    invD4 = 1.0 / D_sigma

    return {
        "delta": delta,
        "centers": centers,
        "R": R,
        "R_norm": zscore(R.astype(float)),
        "theta": theta,
        "invD4": invD4,
        "invD4_norm": zscore(invD4),
    }


# ============================================================
# 4. Fünflings-Features
# ============================================================

def build_quintuplet_features(quintuplets: np.ndarray):
    q = quintuplets

    # 5D -> S^4
    # D_54 = 5 * sum(q^(4/5)) / sum(q)
    D_54 = 5.0 * np.sum(q ** (4.0 / 5.0), axis=1) / np.sum(q, axis=1)

    # 5D -> S^3
    # D_53 = (15/4) * sum(q^(3/5)) / sum(q)
    D_53 = (15.0 / 4.0) * np.sum(q ** (3.0 / 5.0), axis=1) / np.sum(q, axis=1)

    invD54 = 1.0 / D_54
    invD53 = 1.0 / D_53

    return {
        "D_54": D_54,
        "D_53": D_53,
        "invD54": invD54,
        "invD53": invD53,
        "invD54_norm": zscore(invD54),
        "invD53_norm": zscore(invD53),
    }


# ============================================================
# 5. Hybridmodell
# ============================================================

def build_hybrid_block_matrix(
    quadruplets: np.ndarray,
    quintuplets: np.ndarray,
    alpha: float = 2.0,
    alpha3: float = 1.2,
    beta: float = 0.05,
    rho: float = 0.8,
    gamma0: float = 0.6,
    gamma1: float = 0.05,
    gamma2: float = 0.10,
    zeta4: float = 0.0,   # Vierlings-S^3-Term: 1 / D_sigma
    zeta5: float = 0.0,   # Fünflings-Filter: 1 / D_53
    use_242: bool = True,
):
    """
    Diagonale:
        d_{k,j} = alpha * delta_{k,j}
                + alpha3 * delta_{k,j}^3
                + beta * R_norm_k
                + zeta4 * invD4_norm_k
                + zeta5 * invD53_norm_k

    Interblock-Kopplung:
        gamma_k = gamma0
                + gamma1 (R_{k+1} - R_k)
                + gamma2 (F_{k+1} - F_k)
    mit F_k = invD53_norm_k
    """
    if len(quadruplets) != len(quintuplets):
        raise ValueError("Vierlings- und Fünflingslisten müssen gleich lang sein.")

    qfeat = build_quadruplet_features(quadruplets)
    pfeat = build_quintuplet_features(quintuplets)

    delta = qfeat["delta"]
    R_norm = qfeat["R_norm"]
    invD4_norm = qfeat["invD4_norm"]

    invD53_norm = pfeat["invD53_norm"]

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

        extra = zeta4 * invD4_norm[k] + zeta5 * invD53_norm[k]

        diag = (
            alpha * delta[k]
            + alpha3 * (delta[k] ** 3)
            + beta * R_norm[k]
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
                + gamma2 * (invD53_norm[k + 1] - invD53_norm[k])
            )
            C = gk * np.eye(4)
            M[idx:idx+4, idx+4:idx+8] = C
            M[idx+4:idx+8, idx:idx+4] = C.T

    return M, qfeat, pfeat


# ============================================================
# 6. Auswertung
# ============================================================

def evaluate_hybrid_model(
    quadruplets: np.ndarray,
    quintuplets: np.ndarray,
    n_zeros: int = 100,
    zeros_path: str = "zeros6.npy",
    **pars,
):
    M, qfeat, pfeat = build_hybrid_block_matrix(quadruplets, quintuplets, **pars)

    eigvals = np.sort(eigh(M, eigvals_only=True))
    pos = eigvals[eigvals > 1e-10]

    if len(pos) < n_zeros:
        raise ValueError(f"Zu wenige positive Eigenwerte: {len(pos)} < {n_zeros}")

    lam = pos[:n_zeros]
    gam = load_riemann_gammas(n_zeros, path=zeros_path)

    a, b = affine_fit(lam, gam)
    pred = a * lam + b

    rmse, mae, spacing_corr = score_prediction(pred, gam)

    return {
        "eigvals": eigvals,
        "lambda_used": lam,
        "gammas": gam,
        "pred": pred,
        "a": a,
        "b": b,
        "rmse": rmse,
        "mae": mae,
        "spacing_corr": spacing_corr,
        "quadruplet_features": qfeat,
        "quintuplet_features": pfeat,
    }


# ============================================================
# 7. Modellfamilien H0-H3
# ============================================================

def hybrid_family_params(model_name: str):
    """
    H0: nur lokale Vierlingsstruktur
    H1: + Vierlings-S^3-Term
    H2: + Fünflings-Filter
    H3: + beides
    """
    base = dict(
        alpha=2.0,
        alpha3=1.2,
        beta=0.05,
        rho=0.8,
        gamma0=0.6,
        gamma1=0.05,
        gamma2=0.10,
        use_242=True,
    )

    if model_name == "H0":
        return dict(base, zeta4=0.0, zeta5=0.0)
    elif model_name == "H1":
        return dict(base, zeta4=0.20, zeta5=0.0)
    elif model_name == "H2":
        return dict(base, zeta4=0.0, zeta5=0.35)
    elif model_name == "H3":
        return dict(base, zeta4=0.20, zeta5=0.35)
    else:
        raise ValueError(f"Unbekanntes Hybridmodell: {model_name}")


# ============================================================
# 8. Hauptvergleich
# ============================================================

def run_hybrid_scan(
    n: int = 100,
    zeros_path: str = "zeros6.npy",
    quint_pattern: str = "026812",
):
    quadruplets = first_prime_quadruplets(n)

    if quint_pattern == "026812":
        quintuplets = first_prime_quintuplets_pattern_026812(n)
    elif quint_pattern == "0461012":
        quintuplets = first_prime_quintuplets_pattern_0461012(n)
    else:
        raise ValueError("quint_pattern muss '026812' oder '0461012' sein.")

    rows = []
    detailed = {}

    for model_name in ("H0", "H1", "H2", "H3"):
        pars = hybrid_family_params(model_name)
        res = evaluate_hybrid_model(
            quadruplets=quadruplets,
            quintuplets=quintuplets,
            n_zeros=min(n, 100),
            zeros_path=zeros_path,
            **pars,
        )

        rows.append({
            "model": model_name,
            "rmse": res["rmse"],
            "mae": res["mae"],
            "spacing_corr": res["spacing_corr"],
            "a": res["a"],
            "b": res["b"],
            "zeta4": pars["zeta4"],
            "zeta5": pars["zeta5"],
        })

        detailed[model_name] = res

    df = pd.DataFrame(rows).sort_values(["rmse", "mae"]).reset_index(drop=True)
    return df, detailed


# ============================================================
# 9. Ausgabe
# ============================================================

def print_hybrid_summary(df: pd.DataFrame, detailed: dict, show_first: int = 8):
    print("=" * 78)
    print("HYBRIDMODELL: VIERLINGE (lokal) + FÜNFLINGSFILTER (global)")
    print("=" * 78)
    print(df.to_string(index=False))
    print()

    for model_name in df["model"]:
        res = detailed[model_name]
        print("-" * 78)
        print(f"MODELL {model_name}")
        print("-" * 78)
        print(f"RMSE         : {res['rmse']:.6f}")
        print(f"MAE          : {res['mae']:.6f}")
        print(f"Spacing Corr : {res['spacing_corr']:.6f}")
        print(f"Affine Fit   : gamma ≈ {res['a']:.6f} * lambda + {res['b']:.6f}")
        print()
        print(" n |   lambda_n   | affine(lambda_n) |   gamma_n")
        print("-" * 62)
        lam = res["lambda_used"]
        pred = res["pred"]
        gam = res["gammas"]
        for i in range(min(show_first, len(lam))):
            print(f"{i+1:2d} | {lam[i]:11.6f} | {pred[i]:16.6f} | {gam[i]:10.6f}")
        print()


# ============================================================
# 10. Start
# ============================================================

if __name__ == "__main__":
    # Muster 1: (0,2,6,8,12)
    df1, detailed1 = run_hybrid_scan(
        n=100,
        zeros_path="zeros6.npy",
        quint_pattern="026812",
    )

    print("\n\n===== HYBRIDSCAN MIT FÜNFLINGSMUSTER (0,2,6,8,12) =====\n")
    print_hybrid_summary(df1, detailed1, show_first=8)

    # Muster 2: (0,4,6,10,12)
    df2, detailed2 = run_hybrid_scan(
        n=100,
        zeros_path="zeros6.npy",
        quint_pattern="0461012",
    )

    print("\n\n===== HYBRIDSCAN MIT FÜNFLINGSMUSTER (0,4,6,10,12) =====\n")
    print_hybrid_summary(df2, detailed2, show_first=8)