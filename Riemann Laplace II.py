#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt


# ----------------------------
# Riemann zeros (your file)
# ----------------------------
def load_riemann_zeros(path="zeros6.npy", n_take=200_000):
    zeros = np.load(path).astype(float)
    return zeros[:n_take]


def N_riemann_smooth(T):
    """
    Riemann-von Mangoldt counting function (smooth approximation):
      N(T) ~ (T/2π) log(T/2π) - (T/2π) + 7/8
    """
    T = np.asarray(T, dtype=float)
    return (T / (2.0 * np.pi)) * (np.log(T / (2.0 * np.pi)) - 1.0) + 7.0 / 8.0


def unfolded_spacings_from_times(t):
    """Given unfolded coordinates t_i, return normalized nearest-neighbor spacings."""
    s = np.diff(t)
    s = s / np.mean(s)
    return s


# ----------------------------
# Tetrahedron graph Laplacian
# ----------------------------
def tetra_points(n):
    pts = []
    idx = {}
    for i in range(n + 1):
        for j in range(n + 1 - i):
            for k in range(n + 1 - i - j):
                idx[(i, j, k)] = len(pts)
                pts.append((i, j, k))
    return pts, idx


def laplacian_tetra(n):
    pts, idx = tetra_points(n)
    N = len(pts)
    dirs = [(1, 0, 0), (-1, 0, 0),
            (0, 1, 0), (0, -1, 0),
            (0, 0, 1), (0, 0, -1)]
    rows, cols, data = [], [], []
    deg = np.zeros(N, dtype=float)

    for a, (i, j, k) in enumerate(pts):
        for di, dj, dk in dirs:
            p = (i + di, j + dj, k + dk)
            b = idx.get(p, None)
            if b is not None:
                rows.append(a)
                cols.append(b)
                data.append(-1.0)
                deg[a] += 1.0

    # diagonal degree
    rows.extend(range(N))
    cols.extend(range(N))
    data.extend(deg.tolist())

    L = sp.csr_matrix((data, (rows, cols)), shape=(N, N))
    return L


# ----------------------------
# Spectral preprocessing
# ----------------------------
def cluster_degeneracies(eigs, tol=1e-8):
    """
    Collapse numerically degenerate eigenvalues.
    Returns unique eigenvalues (sorted), and multiplicities.
    """
    eigs = np.sort(np.asarray(eigs, dtype=float))
    uniq = [eigs[0]]
    mult = [1]
    for x in eigs[1:]:
        if abs(x - uniq[-1]) <= tol:
            mult[-1] += 1
        else:
            uniq.append(x)
            mult.append(1)
    return np.array(uniq), np.array(mult, dtype=int)


def bulk_slice(arr, qlo=0.2, qhi=0.8):
    """Take middle quantiles of a sorted array."""
    n = len(arr)
    lo = int(np.floor(qlo * n))
    hi = int(np.floor(qhi * n))
    hi = max(hi, lo + 3)
    return arr[lo:hi]


def local_unfolding(eigs_sorted, window=101):
    """
    Robust local unfolding: map eigenvalues λ_i to a smoothed rank function t_i.
    We use a moving average of ranks over a sliding window.
    """
    lam = np.asarray(eigs_sorted, dtype=float)
    n = len(lam)
    if n < window + 2:
        raise ValueError(f"Need more eigenvalues for local unfolding (have {n}, need > {window+2}).")

    ranks = np.arange(1, n + 1, dtype=float)

    # Moving-average smoothing of ranks as a function of index
    # Equivalent to smoothing N(λ) locally without fitting a global polynomial.
    w = window
    kernel = np.ones(w, dtype=float) / w
    smooth = np.convolve(ranks, kernel, mode="same")

    # Ensure strictly increasing (numerical safety)
    # (If not strictly inc due to edge effects, enforce monotonicity.)
    t = np.maximum.accumulate(smooth)

    # Normalize so that mean spacing ~ 1 automatically after spacing normalization
    return t


# ----------------------------
# Main experiment
# ----------------------------
def main():
    # --- parameters you can tweak ---
    zeros_path = "zeros6.npy"
    n_zeros = 200_000

    n_tetra = 12       # increase if you want a larger graph (18 = sehr langsam)
    k_eigs = 400       # number of eigenvalues to request (1200 bei n=18 = extrem langsam)
    deg_tol = 1e-8

    bulk_lo, bulk_hi = 0.2, 0.8
    unfold_window = 101

    # --- Riemann side ---
    Z = load_riemann_zeros(zeros_path, n_take=n_zeros)
    tZ = N_riemann_smooth(Z)
    sZ = unfolded_spacings_from_times(tZ)
    sZ = sZ[(sZ >= 0) & (sZ <= 3)]

    # --- Tetra Laplacian side ---
    L = laplacian_tetra(n_tetra)

    # Compute eigenvalues near 0 (shift-invert). For spacing stats, we will take BULK later anyway.
    # Note: For bigger n/k this can be slow; adjust n_tetra/k_eigs as needed.
    vals = spla.eigsh(
        L,
        k=k_eigs,
        sigma=1e-8,
        which="LM",
        return_eigenvectors=False,
        tol=1e-7,
        maxiter=5000,
    )
    vals = np.sort(vals)

    # Drop the (near-)zero eigenvalue of the Laplacian
    vals = vals[vals > 1e-9]

    # Collapse degeneracies
    uniq, mult = cluster_degeneracies(vals, tol=deg_tol)

    # Bulk only
    uniq_bulk = bulk_slice(uniq, qlo=bulk_lo, qhi=bulk_hi)

    # Local unfolding (smoothed rank)
    tL = local_unfolding(uniq_bulk, window=unfold_window)
    sL = unfolded_spacings_from_times(tL)
    sL = sL[(sL >= 0) & (sL <= 3)]

    # --- diagnostics ---
    print("=== RIEMANN ===")
    print("n zeros used:", len(Z))
    print("mean spacing:", float(np.mean(np.diff(tZ) / np.mean(np.diff(tZ)))))
    print("var  spacing:", float(np.var(np.diff(tZ) / np.mean(np.diff(tZ)))))

    print("\n=== TETRA LAPLACIAN ===")
    print("n_tetra:", n_tetra)
    print("k_eigs requested:", k_eigs)
    print("eigs (nonzero) raw:", len(vals))
    print("unique after clustering:", len(uniq))
    print("bulk unique:", len(uniq_bulk))
    print("mean spacing:", float(np.mean(sL)))
    print("var  spacing:", float(np.var(sL)))
    print("min spacing:", float(np.min(sL)))
    print("max spacing:", float(np.max(sL)))
    print("top multiplicities (largest 10):", np.sort(mult)[-10:])

    # --- plot comparison ---
    def wigner_gue_pdf(s):
        return (32.0 / np.pi**2) * s**2 * np.exp(-4.0 * s**2 / np.pi)

    bins = np.linspace(0, 3, 80)
    centers = 0.5 * (bins[1:] + bins[:-1])

    hZ, _ = np.histogram(sZ, bins=bins, density=True)
    hL, _ = np.histogram(sL, bins=bins, density=True)
    wg = wigner_gue_pdf(centers)

    plt.figure()
    plt.plot(centers, hZ, label="Riemann zeros (unfolded)")
    plt.plot(centers, hL, label=f"Tetra Laplacian (clustered+bulk+local unfold), n={n_tetra}")
    plt.plot(centers, wg, label="GUE Wigner surmise")
    plt.xlabel("s")
    plt.ylabel("density")
    plt.legend()
    plt.title("Nearest-neighbor spacing comparison")
    plt.savefig("riemann_laplace_comparison.png", dpi=150)
    plt.close()
    print("\nPlot gespeichert: riemann_laplace_comparison.png")


if __name__ == "__main__":
    main()