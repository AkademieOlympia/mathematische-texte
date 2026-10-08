#!/usr/bin/env python3
"""
e8_relative_entropy_curvature_audit.py

Reproducible numerical audit for an E8 shape-deformation relative-entropy
curvature.

Definition used:
    C_theta = sum_{r in E8 roots} (A_theta r)(A_theta r)^T
    rho_theta = C_theta / Tr(C_theta)
    rho_0     = C_0 / Tr(C_0)

and the finite-dimensional Umegaki relative entropy
    D(rho || sigma) = Tr[rho (log rho - log sigma)].

Default deformation:
    A_theta = I + theta E_{12}
which is an exactly unimodular elementary shear (det = 1).

IMPORTANT:
- The numerical curvature depends on the chosen deformation A_theta and on
  the definition/normalization of rho_theta.
- This script does NOT hard-code the target value 1.221436.
- If a manuscript claims 1.221436, the exact deformation family used there
  must reproduce it here independently.

Requirements:
    numpy
"""

from __future__ import annotations
import json
import math
from pathlib import Path
import argparse
import numpy as np


def e8_roots():
    """Return the standard 240 E8 roots in R^8, each of squared norm 2."""
    roots = []

    # Type I: permutations of (±1, ±1, 0,...,0), all sign choices.
    for i in range(8):
        for j in range(i + 1, 8):
            for si in (-1.0, 1.0):
                for sj in (-1.0, 1.0):
                    v = np.zeros(8)
                    v[i] = si
                    v[j] = sj
                    roots.append(v)

    # Type II: (±1/2,...,±1/2) with an even number of minus signs.
    for mask in range(1 << 8):
        minus_count = mask.bit_count()
        if minus_count % 2 == 0:
            v = np.array([
                -0.5 if ((mask >> k) & 1) else 0.5
                for k in range(8)
            ], dtype=float)
            roots.append(v)

    R = np.asarray(roots, dtype=float)
    assert R.shape == (240, 8)
    norms2 = np.sum(R * R, axis=1)
    if not np.allclose(norms2, 2.0, atol=1e-14):
        raise RuntimeError("E8 root normalization failed.")
    return R


def elementary_shear(theta: float, i: int = 0, j: int = 1):
    """
    A(theta) = I + theta E_ij.
    Since E_ij^2=0 for i != j, det A(theta)=1 exactly.
    """
    if i == j:
        raise ValueError("Need i != j for elementary shear.")
    A = np.eye(8)
    A[i, j] = theta
    return A


def density_from_roots(R, A):
    X = R @ A.T
    C = X.T @ X
    C = 0.5 * (C + C.T)
    tr = np.trace(C)
    if tr <= 0:
        raise RuntimeError("Non-positive trace.")
    rho = C / tr
    return rho


def log_psd(M, tol=1e-15):
    vals, vecs = np.linalg.eigh(0.5 * (M + M.T))
    if np.min(vals) <= tol:
        raise RuntimeError(
            f"Matrix is not strictly positive definite; min eigenvalue={np.min(vals)}"
        )
    return (vecs * np.log(vals)) @ vecs.T


def relative_entropy(rho, sigma):
    logrho = log_psd(rho)
    logsigma = log_psd(sigma)
    return float(np.trace(rho @ (logrho - logsigma)))


def curvature_central(f, h):
    f0 = f(0.0)
    fp = f(h)
    fm = f(-h)
    return (fp - 2.0 * f0 + fm) / (h * h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="artifacts/e8_relative_entropy_audit")
    ap.add_argument(
        "--hs",
        default="1e-1,5e-2,2e-2,1e-2,5e-3,2e-3,1e-3,5e-4",
        help="Comma-separated step sizes for central second derivative."
    )
    ap.add_argument("--shear-i", type=int, default=0)
    ap.add_argument("--shear-j", type=int, default=1)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    R = e8_roots()

    # Root-system isotropy check.
    C0_raw = R.T @ R
    isotropic_scalar = np.trace(C0_raw) / 8.0
    isotropy_residual = np.linalg.norm(
        C0_raw - isotropic_scalar * np.eye(8), ord="fro"
    )

    A0 = elementary_shear(0.0, args.shear_i, args.shear_j)
    rho0 = density_from_roots(R, A0)

    def D(theta):
        A = elementary_shear(theta, args.shear_i, args.shear_j)
        rho = density_from_roots(R, A)
        return relative_entropy(rho, rho0)

    hs = [float(x.strip()) for x in args.hs.split(",") if x.strip()]
    rows = []
    for h in hs:
        curv = curvature_central(D, h)
        rows.append({
            "h": h,
            "D_plus": D(h),
            "D_zero": D(0.0),
            "D_minus": D(-h),
            "curvature_central": curv,
        })

    # A simple Richardson estimate from the two smallest h values:
    # c(h)=c0 + O(h^2), so c0 ≈ (4 c(h/2)-c(h))/3 if h ratio is exactly 2.
    richardson = None
    if len(rows) >= 2:
        r1, r2 = rows[-2], rows[-1]
        h1, h2 = r1["h"], r2["h"]
        if abs(h1 / h2 - 2.0) < 1e-12:
            richardson = (
                4.0 * r2["curvature_central"] - r1["curvature_central"]
            ) / 3.0

    eig0 = np.linalg.eigvalsh(rho0)

    summary = {
        "status": "NUMERICAL_IMPLEMENTATION_AUDIT",
        "root_count": int(R.shape[0]),
        "root_norm_squared": 2.0,
        "deformation": {
            "type": "elementary_unimodular_shear",
            "A_theta": f"I + theta E_({args.shear_i+1},{args.shear_j+1})",
            "determinant": 1.0,
        },
        "density_definition": "rho_theta = C_theta / Tr(C_theta), C_theta = sum (A_theta r)(A_theta r)^T",
        "relative_entropy": "finite-dimensional Umegaki",
        "C0_isotropic_scalar": float(isotropic_scalar),
        "C0_isotropy_frobenius_residual": float(isotropy_residual),
        "rho0_eigenvalues": [float(x) for x in eig0],
        "D_at_zero": D(0.0),
        "curvature_rows": rows,
        "richardson_estimate_if_available": richardson,
        "warning": (
            "The curvature is deformation- and normalization-dependent. "
            "A claimed value such as 1.221436 is certified only if the exact "
            "manuscript deformation and density definition reproduce it."
        ),
    }

    out = outdir / "e8_relative_entropy_curvature_summary.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("=" * 82)
    print("E8 RELATIVE-ENTROPY CURVATURE AUDIT")
    print("=" * 82)
    print(f"E8 roots                    : {R.shape[0]}")
    print(f"C0 isotropy residual        : {isotropy_residual:.3e}")
    print(f"rho0 eigenvalues            : {eig0}")
    print(f"D(rho0 || rho0)             : {D(0.0):.3e}")
    print("-" * 82)
    print(" h              D(+h)              D(-h)              curvature")
    for row in rows:
        print(
            f"{row['h']:<12.3e} "
            f"{row['D_plus']:<18.10e} "
            f"{row['D_minus']:<18.10e} "
            f"{row['curvature_central']:.12f}"
        )
    if richardson is not None:
        print("-" * 82)
        print(f"Richardson curvature estimate: {richardson:.12f}")
    print(f"[Export] {out}")


if __name__ == "__main__":
    main()
