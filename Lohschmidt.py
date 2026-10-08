# ================================================================
# EABC S6 DUAL-AUDIT: ANTIPODAL (A) vs. LOSCHMIDT-REVERSAL (Theta)
# Compatible with SageMath 10.5 / Python 3.x
#
# Class A: Z2 central inversion A = -I  (commutes with T; A2 antipodal hit)
# Class T: true reversal involution Theta (anticommutes: Theta T Theta = T^{-1})
# ================================================================

from __future__ import annotations

import numpy as np


def norm(x) -> float:
    return float(np.linalg.norm(np.asarray(x, dtype=float)))


def rel_error(a, b, eps: float = 1e-15) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return norm(a - b) / max(norm(a), norm(b), eps)


def normalize_s6(x, eps: float = 1e-15) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    n = norm(x)
    return x / n if n > eps else x


# ----------------------------------------------------------------
# 1. Dihedral / D6 block-rotation kernel on R^7 ~ Im(O)
# ----------------------------------------------------------------

def build_dihedral_kernel() -> np.ndarray:
    """Three 2D planes with theta=pi/3 plus last-axis flip => R^3 = -I."""
    theta = np.pi / 3.0
    R = np.eye(7)
    c, s = np.cos(theta), np.sin(theta)
    for i in range(3):
        a, b = 2 * i, 2 * i + 1
        R[a, a] = c
        R[a, b] = -s
        R[b, a] = s
        R[b, b] = c
    R[6, 6] = -1.0
    return R


R_kernel = build_dihedral_kernel()
assert np.allclose(np.linalg.matrix_power(R_kernel, 3), -np.eye(7), atol=1e-12)


def T_sys(x) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = R_kernel @ x + 1e-4 * np.sin(x)
    return normalize_s6(y)


def build_theta_involution() -> np.ndarray:
    """
    Block reflection: flip the second axis of each rotational plane.
    Then Theta R Theta = R(-theta) = R^{-1} on each 2D block.
    Last axis: Theta_66 = 1 leaves R_66 = -1 invariant (and equals its inverse).
    """
    Theta = np.eye(7)
    for i in range(3):
        Theta[2 * i + 1, 2 * i + 1] = -1.0
    Theta[6, 6] = 1.0
    return Theta


Theta_mat = build_theta_involution()

# Algebraic check on the linear kernel: Theta R Theta = R^{-1}
_R_conj = Theta_mat @ R_kernel @ Theta_mat
_R_inv = np.linalg.inv(R_kernel)
assert np.allclose(_R_conj, _R_inv, atol=1e-12), "Theta must conjugate R to R^{-1}"


def Theta_op(x) -> np.ndarray:
    return Theta_mat @ np.asarray(x, dtype=float)


def A_op(x) -> np.ndarray:
    """Central / antipodal inversion A = -I."""
    return -np.asarray(x, dtype=float)


def Tinv_theta(x) -> np.ndarray:
    """Loschmidt backward map induced by Theta: T^{-1} := Theta T Theta."""
    return Theta_op(T_sys(Theta_op(x)))


def evolve_forward(X, n: int) -> np.ndarray:
    Y = np.asarray(X, dtype=float).copy()
    for _ in range(n):
        Y = T_sys(Y)
    return Y


# ----------------------------------------------------------------
# 2. Dual-Audit Execution
# ----------------------------------------------------------------

def run_dual_audit(n_ens: int = 300, seed: int = 2026) -> None:
    rng = np.random.default_rng(seed)
    ensemble = [normalize_s6(rng.standard_normal(7)) for _ in range(n_ens)]

    print("=" * 72)
    print(" EABC S6 DUAL AUDIT: ANTIPODAL (A) vs. LOSCHMIDT (Theta)")
    print("=" * 72)

    # --- Class A: Antipodal Audit ---
    err_A1 = max(rel_error(A_op(A_op(X)), X) for X in ensemble)

    # A2: specifically test designed resonance T^3 X ~= A X
    success_A2 = 0
    for X in ensemble:
        if rel_error(evolve_forward(X, 3), A_op(X)) < 1e-3:
            success_A2 += 1

    # A3: Z2 commutation T A = A T  (oddness of T)
    err_A3 = max(rel_error(T_sys(A_op(X)), A_op(T_sys(X))) for X in ensemble)

    print(f"[A1] Antipodal Involution A^2 = I err : {err_A1:.2e}")
    print(f"[A2] Dynamic Antipodal Hit (n*=3 rate): {success_A2}/{n_ens}")
    print(f"[A3] Z2 Commutation T A = A T err     : {err_A3:.2e}")

    # --- Class T: True Loschmidt Reversal ---
    err_T1 = max(rel_error(Theta_op(Theta_op(X)), X) for X in ensemble)

    # T2a: definitional covariance Theta T Theta == Tinv_theta (should be ~0)
    err_T2_def = max(
        rel_error(Theta_op(T_sys(Theta_op(X))), Tinv_theta(X)) for X in ensemble
    )
    # T2b: dynamical inverse property T o Tinv ~= id and Tinv o T ~= id
    err_T2_fwd = max(rel_error(T_sys(Tinv_theta(X)), X) for X in ensemble)
    err_T2_bwd = max(rel_error(Tinv_theta(T_sys(X)), X) for X in ensemble)
    err_T2 = max(err_T2_def, err_T2_fwd, err_T2_bwd)

    # T3: does the forward orbit ever hit the true reversed configuration Theta(X)?
    success_T3 = 0
    n_stars: list[int] = []
    for X in ensemble:
        curr = np.asarray(X, dtype=float)
        hit_n = None
        for n in range(1, 15):
            curr = T_sys(curr)
            if rel_error(curr, Theta_op(X)) < 1e-3:
                hit_n = n
                break
        if hit_n is not None:
            success_T3 += 1
            n_stars.append(hit_n)

    # T4: orbit retracing around a reference n_star (use modal T3 hit if any, else 3)
    if n_stars:
        n_ref = int(np.bincount(n_stars).argmax())
    else:
        n_ref = 3
    X0 = ensemble[0]
    err_T4 = 0.0
    for m in range(0, n_ref + 1):
        future = evolve_forward(X0, n_ref + m)
        # past via Theta-induced inverse
        past = np.asarray(X0, dtype=float)
        for _ in range(n_ref - m):
            past = Tinv_theta(past)
        # Loschmidt echo: future ~= Theta(past) when echo closes; here check
        # retracing identity under true inverse: T^{n+m} ~= Tinv^{n-m} wait —
        # Standard: T^{n*+m}(X) vs Theta( T^{n*-m}(X) ) for reversible flows.
        # With our Theta-covariance, compare future to Theta_op(evolve_forward(X0, n_ref-m))?
        # Clean T4: Tinv^{n_ref-m}(T^{n_ref+m}(X)) ~= X  (echo closure)
        echoed = future
        for _ in range(n_ref + m):
            echoed = Tinv_theta(echoed)
        err_T4 = max(err_T4, rel_error(echoed, X0))

    print("-" * 72)
    print(f"[T1] Reversal Involution Theta^2 = I err: {err_T1:.2e}")
    print(f"[T2] Loschmidt Covariance / Inverse err : {err_T2:.2e}")
    print(f"     (def={err_T2_def:.2e}, T Tinv={err_T2_fwd:.2e}, Tinv T={err_T2_bwd:.2e})")
    n_star_msg = (
        f"modal n*={int(np.bincount(n_stars).argmax())}" if n_stars else "n* = none"
    )
    print(f"[T3] Dynamic True Reversal Hit rate   : {success_T3}/{n_ens} ({n_star_msg})")
    print(f"[T4] Orbit Retracing / Echo Closure   : {err_T4:.2e} (ref n={n_ref})")
    print("=" * 72)

    print("Interpretation:")
    print("  A2 = antipodal transit (central Z2), NOT Loschmidt time reversal.")
    print("  T2/T4 = genuine reversal covariance under block reflection Theta.")
    if success_A2 == n_ens and success_T3 == 0:
        print("  STATUS: A2 PASS / T3 idle — dihedral antipode, not dynamical Theta-hit.")
    elif success_T3 > 0:
        print("  STATUS: Partial/full dynamical Theta generation detected.")
    else:
        print("  STATUS: Neither A2 nor T3 fully realized under current thresholds.")


if __name__ == "__main__":
    run_dual_audit()
