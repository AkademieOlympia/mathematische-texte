#!/usr/bin/env python3
# ================================================================
# EABC S6 — SYMMETRY BREAKING PARAMETER SCAN (epsilon)
# Compatible with SageMath 10.5 / Python 3.x
#
# Boxed claim:
#   Reversibility of the laws  =/=  spontaneous reachability of Theta(X)
#
# Observables:
#   D_A(eps)   = mean || T_eps^3 X - A X ||_rel     (D6 antipodal resonance)
#   E_Theta    = mean || Theta T Theta X - Tinv X || (algebraic covariance;
#                ~0 by construction if Tinv := Theta T Theta)
#   E_id(eps)  = mean || T(Tinv(X)) - X ||_rel      (dynamical inverse closure)
#   P_T3(eps)  = Prob[ min_{n<=N} ||T^n X - Theta X|| < delta ]
# ================================================================

from __future__ import annotations

import argparse
from pathlib import Path

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


def build_dihedral_kernel() -> np.ndarray:
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


def build_theta_involution() -> np.ndarray:
    Theta = np.eye(7)
    for i in range(3):
        Theta[2 * i + 1, 2 * i + 1] = -1.0
    Theta[6, 6] = 1.0
    return Theta


R_kernel = build_dihedral_kernel()
Theta_mat = build_theta_involution()


def Theta_op(x) -> np.ndarray:
    return Theta_mat @ np.asarray(x, dtype=float)


def A_op(x) -> np.ndarray:
    return -np.asarray(x, dtype=float)


def make_perturbed_system(eps: float):
    def T_eps(x):
        x = np.asarray(x, dtype=float)
        y = R_kernel @ x + eps * np.sin(x)  # odd perturbation
        return normalize_s6(y)

    def Tinv_eps(x):
        # Algebraic Loschmidt covariance: T^{-1} := Theta T Theta
        return Theta_op(T_eps(Theta_op(x)))

    return T_eps, Tinv_eps


def run_scan(
    n_ens: int = 150,
    seed: int = 2026,
    n_pt3: int = 10,
    delta: float = 1e-3,
    export: str | None = None,
) -> np.ndarray:
    print("=" * 88)
    print(" EABC S6 SYMMETRY BREAKING PHASE SCAN: eps in [1e-8, 1e-1]")
    print(" Claim: law-reversibility ≠ spontaneous reachability of Theta(X)")
    print("=" * 88)
    header = (
        f"{'Epsilon':<12} | {'DA (Antipode)':<15} | {'ETheta (alg)':<14} | "
        f"{'Eid (T Tinv)':<14} | {'PT3 (Spont. T)':<15}"
    )
    print(header)
    print("-" * 88)

    rng = np.random.default_rng(seed)
    ensemble = [normalize_s6(rng.standard_normal(7)) for _ in range(n_ens)]
    eps_values = np.logspace(-8, -1, 8)

    rows = []
    for eps in eps_values:
        T_eps, Tinv_eps = make_perturbed_system(float(eps))

        # D_A: antipodal resonance error at n=3
        da_list = []
        for X in ensemble:
            curr = np.asarray(X, dtype=float)
            for _ in range(3):
                curr = T_eps(curr)
            da_list.append(rel_error(curr, A_op(X)))
        DA_mean = float(np.mean(da_list))

        # E_Theta: algebraic covariance (definitional if Tinv := Theta T Theta)
        et_list = []
        eid_list = []
        for X in ensemble:
            lhs = Theta_op(T_eps(Theta_op(X)))
            rhs = Tinv_eps(X)
            et_list.append(rel_error(lhs, rhs))
            eid_list.append(rel_error(T_eps(Tinv_eps(X)), X))
        ET_mean = float(np.mean(et_list))
        EID_mean = float(np.mean(eid_list))

        # P_T3: spontaneous hit of Theta(X)
        success_pt3 = 0
        for X in ensemble:
            curr = np.asarray(X, dtype=float)
            target = Theta_op(X)
            reached = False
            for _ in range(n_pt3):
                curr = T_eps(curr)
                if rel_error(curr, target) < delta:
                    reached = True
                    break
            if reached:
                success_pt3 += 1
        PT3_rate = success_pt3 / len(ensemble)

        print(
            f"{eps:<12.1e} | {DA_mean:<15.2e} | {ET_mean:<14.2e} | "
            f"{EID_mean:<14.2e} | {PT3_rate:<15.4f}"
        )
        rows.append((eps, DA_mean, ET_mean, EID_mean, PT3_rate))

    print("=" * 88)
    print("Reading guide:")
    print("  DA   -> robustness of A2 / D6 antipodal transit under NL")
    print("  ETheta -> algebraic Theta-covariance (odd NL keeps this ~0)")
    print("  Eid  -> dynamical inverse closure T(Tinv(X))~X (true NL stress)")
    print("  PT3  -> spontaneous reachability of Theta(X) (expected ~0)")

    data = np.array(rows, dtype=float)
    if export:
        out = Path(export)
        out.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(
            out,
            data,
            delimiter=",",
            header="eps,DA,ETheta,Eid,PT3",
            comments="",
        )
        print(f"[Export] {out}")
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="EABC S6 epsilon symmetry-breaking scan")
    parser.add_argument("--n-ens", type=int, default=150)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument(
        "--export",
        default="eabc_epsilon_scan.csv",
        help="CSV output path (relative to CWD)",
    )
    parser.add_argument("--no-export", action="store_true")
    args = parser.parse_args()
    run_scan(
        n_ens=args.n_ens,
        seed=args.seed,
        export=None if args.no_export else args.export,
    )
