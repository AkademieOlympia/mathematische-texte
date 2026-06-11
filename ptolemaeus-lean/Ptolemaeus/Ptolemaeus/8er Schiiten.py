import math
import numpy as np
import sympy as sp
from mpmath import zetazero
from scipy.linalg import eigh


# ============================================================
# 1. Primzahlvierlinge erzeugen
# ============================================================

def first_prime_quadruplets(n: int):
    """
    Erzeugt die ersten n Primzahlvierlinge der Form
    (p, p+2, p+6, p+8).
    """
    quadruplets = []
    p = 5
    while len(quadruplets) < n:
        if sp.isprime(p) and sp.isprime(p + 2) and sp.isprime(p + 6) and sp.isprime(p + 8):
            quadruplets.append((p, p + 2, p + 6, p + 8))
        p += 2
    return np.array(quadruplets, dtype=float)


# ============================================================
# 2. Hilfsgrößen für das Modell
# ============================================================

def zscore(x: np.ndarray, eps: float = 1e-12):
    mu = np.mean(x)
    sigma = np.std(x)
    return (x - mu) / (sigma + eps)


def build_features(quadruplets: np.ndarray):
    """
    Baut die Modellgrößen:
    - delta_{k,j}: zentrierte Log-Abweichung innerhalb eines Vierlings
    - m_k: Vierlingszentrum
    - R_k: dyadischer Achterrahmenindex
    - theta_k: Lage im dyadischen Rahmen
    - normierte Versionen von R_k
    """
    logs = np.log(quadruplets)
    mu = logs.mean(axis=1, keepdims=True)
    delta = logs - mu

    centers = quadruplets.mean(axis=1)

    R = np.floor(np.log2(centers / 8.0)).astype(int)
    theta = centers / (2.0 ** (R + 3)) - 1.0

    # Sicherheitscheck
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
# 3. Dyadisch erweitertes 4-Block-Modell
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
):
    """
    Baut die symmetrische Blockmatrix M_N.

    Lokaler Block:
        d_{k,j} = alpha * delta_{k,j} + beta * R_norm_k + eta * cos(2π theta_k)

    Innere Kopplung:
        u12 = 1 + rho * 2/8
        u23 = 1 + rho * 4/8
        u34 = 1 + rho * 2/8

    Inter-Block-Kopplung:
        C_k = gamma_k I_4
        gamma_k = gamma0 + gamma1(R_norm_{k+1} - R_norm_k)
                        + gamma2 cos(2π(theta_{k+1} - theta_k))
    """
    features = build_features(quadruplets)
    delta = features["delta"]
    R_norm = features["R_norm"]
    theta = features["theta"]

    N = len(quadruplets)
    dim = 4 * N
    M = np.zeros((dim, dim), dtype=float)

    u12 = 1.0 + rho * (2.0 / 8.0)
    u23 = 1.0 + rho * (4.0 / 8.0)
    u34 = 1.0 + rho * (2.0 / 8.0)

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
# 4. Riemann-Nullstellen
# ============================================================

def first_riemann_gammas(n: int):
    """
    Erste n imaginäre Teile der nichttrivialen Riemann-Nullstellen.
    """
    return np.array([float(zetazero(k).imag) for k in range(1, n + 1)], dtype=float)


# ============================================================
# 5. Vergleich Modell <-> Riemann-Nullstellen
# ============================================================

def affine_fit(x: np.ndarray, y: np.ndarray):
    """
    Findet a, b mit y ≈ a x + b per kleinster Quadrate.
    """
    A = np.vstack([x, np.ones_like(x)]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(a), float(b)


def evaluate_model(
    n_quadruplets: int = 100,
    n_zeros: int = 100,
    alpha: float = 2.0,
    beta: float = 0.08,
    eta: float = 0.5,
    rho: float = 0.8,
    gamma0: float = 0.6,
    gamma1: float = 0.08,
    gamma2: float = 0.15,
):
    quadruplets = first_prime_quadruplets(n_quadruplets)
    M, features = build_block_matrix(
        quadruplets,
        alpha=alpha,
        beta=beta,
        eta=eta,
        rho=rho,
        gamma0=gamma0,
        gamma1=gamma1,
        gamma2=gamma2,
    )

    # alle Eigenwerte der symmetrischen Matrix
    eigvals = eigh(M, eigvals_only=True)
    eigvals = np.sort(eigvals)

    # positive Eigenwerte nehmen
    pos = eigvals[eigvals > 1e-10]
    if len(pos) < n_zeros:
        raise ValueError(
            f"Zu wenige positive Eigenwerte: {len(pos)} < {n_zeros}. "
            "Parameter oder Modell anpassen."
        )

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
        "quadruplets": quadruplets,
        "features": features,
        "matrix": M,
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
# 6. Ausgabe
# ============================================================

def print_summary(result, show_first=10):
    print("=" * 72)
    print("DYADISCH ERWEITERTES 4-BLOCK-MODELL")
    print("=" * 72)
    print(f"RMSE            : {result['rmse']:.6f}")
    print(f"MAE             : {result['mae']:.6f}")
    print(f"Spacing Corr    : {result['spacing_corr']:.6f}")
    print(f"Affine Fit      : gamma ≈ {result['fit_a']:.6f} * lambda + {result['fit_b']:.6f}")
    print()

    print("Erste 10 Primzahlvierlinge:")
    for i, q in enumerate(result["quadruplets"][:10], start=1):
        print(f"{i:2d}: {tuple(int(x) for x in q)}")
    print()

    print("Erste 10 dyadische Rahmenindizes R_k:")
    print(result["features"]["R"][:show_first])
    print()

    print("Erste 10 dyadische Phasen theta_k:")
    print(np.round(result["features"]["theta"][:show_first], 12))
    print()

    print("Vergleich der ersten Werte:")
    print(" n |   lambda_n   |  affine(lambda_n) |   gamma_n")
    print("-" * 56)
    lam = result["positive_eigvals_used"]
    pred = result["predicted_gammas"]
    gam = result["riemann_gammas"]
    for i in range(min(show_first, len(lam))):
        print(f"{i+1:2d} | {lam[i]:11.6f} | {pred[i]:16.6f} | {gam[i]:10.6f}")


if __name__ == "__main__":
    configs = [
        ("ohne_achterrahmen", dict(alpha=2.0, beta=0.0, eta=0.0, rho=0.8, gamma0=0.6, gamma1=0.0, gamma2=0.0)),
        ("nur_schale", dict(alpha=2.0, beta=0.08, eta=0.0, rho=0.8, gamma0=0.6, gamma1=0.08, gamma2=0.0)),
        ("voll", dict(alpha=2.0, beta=0.08, eta=0.5, rho=0.8, gamma0=0.6, gamma1=0.08, gamma2=0.15)),
    ]

    for name, pars in configs:
        res = evaluate_model(n_quadruplets=100, n_zeros=100, **pars)
        print(name, res["rmse"], res["mae"], res["spacing_corr"])