#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RH + Planck: Darstellung über Frequenz UND Wellenlänge (beides)

- Lädt Riemann-Nullstellen γ aus zeros6.gz (oder .txt)
- Lädt ein Spektrum als CSV:
    Option A: Spalten: nu_Hz,S
    Option B: Spalten: lambda_m,S
    Option C: 2 Spalten ohne Header (wird heuristisch interpretiert)
- Erzeugt:
    1) Planckkurven (klassisch) über ν und über λ
    2) Optional: RH-modulierte Kurven (mult_u oder shift_y)
    3) Optional: Score(T, N)-Heatmap (Temperatur vs. Nullstellenanzahl) für ν und λ

Hinweis:
- Für reale Daten ist "S" nur proportional zu u_ν oder u_λ nötig; absolute Skalierung wird intern normalisiert.
- Für λ-Darstellung wird u_λ aus u_ν konsistent über Planck-Formel berechnet (nicht "einfach umbeschriftet").
"""

import argparse
import gzip
import math
import os
import re
from typing import Tuple, Optional, List

import numpy as np
import matplotlib.pyplot as plt

# -----------------------
# Konstanten (wie in deinen bisherigen Skripten)
# -----------------------
c = 299792458.0
kB = 1.380649e-23
e = 1.602176634e-19
R_K = 25812.80745
h = (e**2) * R_K
hc = h * c
LOG_H_OVER_KB = math.log(h / kB)
LOG_HC_OVER_KB = math.log(hc / kB)

_num_re = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")

def load_zeros(path: str, nmax: int = 10000) -> np.ndarray:
    opener = gzip.open if path.lower().endswith(".gz") else open
    z: List[float] = []
    with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = _num_re.search(line)
            if m:
                z.append(float(m.group(0)))
            if len(z) >= nmax:
                break
    if not z:
        raise ValueError(f"Keine Nullstellen in {path}")
    return np.array(z, dtype=np.float64)

# -----------------------
# Planck: u_nu und u_lambda (Energiedichte, nicht Radiance)
# -----------------------
def planck_u_nu(nu_hz: np.ndarray, T: float) -> np.ndarray:
    nu = np.asarray(nu_hz, dtype=np.float64)
    x = (h * nu) / (kB * T)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        factor = (8.0 * np.pi * h * nu**3) / (c**3)
        return factor / (np.exp(x) - 1.0)

def planck_u_lambda(lam_m: np.ndarray, T: float) -> np.ndarray:
    lam = np.asarray(lam_m, dtype=np.float64)
    x = (hc) / (lam * kB * T)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        return (8.0 * np.pi * hc) / (lam**5) * (1.0 / (np.exp(x) - 1.0))

# -----------------------
# RH template: m(z)=Σ w(γ) cos(γ z)
# -----------------------
def weights_inv_sqrt(gamma: np.ndarray) -> np.ndarray:
    g = np.asarray(gamma, dtype=np.float64)
    return 1.0 / np.sqrt(0.25 + g**2)

def build_template(z_grid: np.ndarray, gamma: np.ndarray, w: np.ndarray, chunk: int = 700) -> np.ndarray:
    z = np.asarray(z_grid, dtype=np.float64)
    M = np.zeros_like(z)
    gamma = np.asarray(gamma, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    for i in range(0, len(gamma), chunk):
        g = gamma[i:i+chunk][:, None]
        ww = w[i:i+chunk][:, None]
        M += (ww * np.cos(g * z[None, :])).sum(axis=0)
    return M

def standardize(v: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    v = np.asarray(v, dtype=np.float64)
    v0 = v - v.mean()
    return v0 / (np.std(v0) + eps)

def z_from_nu(nu_hz: np.ndarray, T: float) -> np.ndarray:
    return np.log(np.asarray(nu_hz, dtype=np.float64)) + LOG_H_OVER_KB - math.log(float(T))

def z_from_lambda(lam_m: np.ndarray, T: float) -> np.ndarray:
    # z = log(hc/(λ k T)) = log(hc/k) - log(λ) - log(T)
    return LOG_HC_OVER_KB - np.log(np.asarray(lam_m, dtype=np.float64)) - math.log(float(T))

# -----------------------
# CSV Loader (nu oder lambda)
# -----------------------
def load_spectrum(path: str) -> Tuple[str, np.ndarray, np.ndarray]:
    """
    Returns: axis_type ("nu" or "lambda"), x, S
    Accepts:
      - header with nu_Hz or lambda_m
      - 2 columns without header -> heuristic:
           if median(x) > 1e6 -> assume Hz (nu)
           else assume meters (lambda)
    """
    import csv
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        sample = f.read(4096)
    # Try detect header
    header = sample.splitlines()[0].lower() if sample else ""
    has_alpha = any(ch.isalpha() for ch in header)
    if has_alpha and (("nu" in header) or ("lambda" in header) or ("lam" in header)):
        # use genfromtxt with names
        data = np.genfromtxt(path, delimiter=",", names=True, dtype=np.float64)
        names = [n.lower() for n in data.dtype.names]
        if any("nu" in n for n in names):
            x = np.asarray(data[[n for n in data.dtype.names if "nu" in n.lower()][0]], dtype=np.float64)
            S = np.asarray(data[[n for n in data.dtype.names if n.lower() in ("s","spec","s_obs","intensity","power") or "s"==n.lower()][0]], dtype=np.float64) if any(n.lower() in ("s","spec","s_obs","intensity","power") or n.lower()=="s" for n in data.dtype.names) else np.asarray(data[data.dtype.names[1]], dtype=np.float64)
            axis="nu"
        else:
            x = np.asarray(data[[n for n in data.dtype.names if "lam" in n.lower() or "lambda" in n.lower()][0]], dtype=np.float64)
            S = np.asarray(data[data.dtype.names[1]], dtype=np.float64)
            axis="lambda"
    else:
        arr = np.genfromtxt(path, delimiter=",", dtype=np.float64)
        if arr.ndim != 2 or arr.shape[1] < 2:
            arr = np.genfromtxt(path, dtype=np.float64)
        x = np.asarray(arr[:,0], dtype=np.float64)
        S = np.asarray(arr[:,1], dtype=np.float64)
        axis = "nu" if np.nanmedian(x) > 1e6 else "lambda"

    m = np.isfinite(x) & np.isfinite(S) & (x > 0) & (S > 0)
    x, S = x[m], S[m]
    idx = np.argsort(x)
    return axis, x[idx], S[idx]

# -----------------------
# Detrend helper (poly)
# -----------------------
def detrend_poly(R: np.ndarray, deg: Optional[int]) -> np.ndarray:
    if deg is None:
        return R
    deg = int(deg)
    x = np.linspace(0.0, 1.0, len(R), dtype=np.float64)
    X = np.vstack([x**k for k in range(deg+1)]).T
    coef, *_ = np.linalg.lstsq(X, R, rcond=None)
    return R - (X @ coef)

# -----------------------
# Score(T,N) für mult_u (Template-Korrelation)
# -----------------------
def score_heatmap_mult_u(axis: str, x: np.ndarray, S: np.ndarray, zeros: np.ndarray,
                         T_grid: np.ndarray, N_list: List[int], detrend_deg: int = 5) -> np.ndarray:
    # Choose z mapping and Planck base
    if axis == "nu":
        wlog = np.log(x)
        z_min = (wlog.min() + LOG_H_OVER_KB - math.log(float(np.max(T_grid)))) - 0.25
        z_max = (wlog.max() + LOG_H_OVER_KB - math.log(float(np.min(T_grid)))) + 0.25
        z_grid = np.linspace(z_min, z_max, 2400)
        def u_base(T): return planck_u_nu(x, T)
        def z_of(T):  return wlog + LOG_H_OVER_KB - math.log(float(T))
    else:
        wlog = np.log(x)
        # z = log(hc/k) - log(lam) - log(T)
        z_min = (LOG_HC_OVER_KB - wlog.max() - math.log(float(np.max(T_grid)))) - 0.25
        z_max = (LOG_HC_OVER_KB - wlog.min() - math.log(float(np.min(T_grid)))) + 0.25
        z_grid = np.linspace(z_min, z_max, 2400)
        def u_base(T): return planck_u_lambda(x, T)
        def z_of(T):  return (LOG_HC_OVER_KB - wlog - math.log(float(T)))

    out = np.zeros((len(N_list), len(T_grid)), dtype=np.float64)

    for j, N in enumerate(N_list):
        g = zeros[:N]
        w = weights_inv_sqrt(g)
        M_grid = build_template(z_grid, g, w, chunk=800)

        for i, T in enumerate(T_grid):
            uT = u_base(float(T))
            R = (S / (uT + 1e-300)) - 1.0
            R = detrend_poly(R, detrend_deg)
            Rt = R - R.mean()

            z = z_of(float(T))
            Mt = np.interp(z, z_grid, M_grid)
            m = standardize(Mt)
            denom = (np.linalg.norm(Rt) * np.linalg.norm(m) + 1e-12)
            out[j, i] = float(np.dot(Rt, m) / denom)
    return out

# -----------------------
# Plots: Planckkurven (klassisch) & optional RH-moduliert
# -----------------------
def plot_planck_both(T_list: List[float], outdir: str) -> None:
    # ν-Plot
    nu = np.linspace(0.1e12, 1200e12, 3000)
    plt.figure(figsize=(11,6))
    for T in T_list:
        u = planck_u_nu(nu, T)
        plt.plot(nu/1e12, u, lw=1.6, label=f"T={T}K")
    plt.yscale("log")
    plt.xlabel("Frequenz ν [THz]")
    plt.ylabel("uν(ν,T) [J·s/m³] (log)")
    plt.title("Planckkurven über Frequenz (absolute Skala, log)")
    plt.grid(True, alpha=0.25, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "planck_over_nu_log.png"), dpi=170, bbox_inches="tight")
    plt.close()

    # λ-Plot (weiter Bereich, damit Abfall sichtbar)
    lam_um = np.linspace(0.5, 300.0, 6000)
    lam_m = lam_um * 1e-6
    plt.figure(figsize=(11,6))
    for T in T_list:
        u = planck_u_lambda(lam_m, T)
        plt.plot(lam_um, u, lw=1.6, label=f"T={T}K")
    plt.yscale("log")
    plt.xlabel("Wellenlänge λ [µm]")
    plt.ylabel("uλ(λ,T) [J·s/m⁴] (log)")
    plt.title("Planckkurven über Wellenlänge (absolute Skala, log)")
    plt.grid(True, alpha=0.25, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "planck_over_lambda_log.png"), dpi=170, bbox_inches="tight")
    plt.close()

    # λ-normalisiert (typische Lehrbuch-Form)
    lam_um2 = np.linspace(0.5, 30.0, 3000)
    lam_m2 = lam_um2 * 1e-6
    plt.figure(figsize=(11,6))
    for T in T_list:
        u = planck_u_lambda(lam_m2, T)
        u_norm = u / (np.max(u) + 1e-300)
        plt.plot(lam_um2, u_norm, lw=1.8, label=f"T={T}K")
    plt.xlabel("Wellenlänge λ [µm]")
    plt.ylabel("relative Intensität (max=1)")
    plt.title("Planckkurven über Wellenlänge (normalisiert, „klassische“ Darstellung)")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "planck_over_lambda_normalized.png"), dpi=170, bbox_inches="tight")
    plt.close()

def plot_heatmap(score_mat: np.ndarray, T_grid: np.ndarray, N_list: List[int], title: str, outpath: str) -> None:
    T_hat = T_grid[np.argmax(score_mat, axis=1)]
    plt.figure(figsize=(11,7))
    im = plt.imshow(score_mat, aspect="auto", origin="lower",
                    extent=[T_grid.min(), T_grid.max(), 0, len(N_list)-1])
    plt.yticks(np.arange(len(N_list)), [str(n) for n in N_list])
    plt.xlabel("T-Kandidaten [K]")
    plt.ylabel("Anzahl Nullstellen N")
    plt.title(title)
    cb = plt.colorbar(im)
    cb.set_label("Score (Korrelation)")
    plt.plot(T_hat, np.arange(len(N_list)), color="white", lw=2.0, alpha=0.9)
    plt.scatter(T_hat, np.arange(len(N_list)), color="white", s=18, alpha=0.9)
    plt.tight_layout()
    plt.savefig(outpath, dpi=170, bbox_inches="tight")
    plt.close()

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--zeros", default="zeros6.gz", help="Riemann-Nullstellen (.gz oder .txt)")
    ap.add_argument("--nzeros", type=int, default=10000, help="Wie viele Nullstellen verwenden (Default 10000)")
    ap.add_argument("--spectrum", default=None, help="CSV: entweder nu_Hz,S oder lambda_m,S")
    ap.add_argument("--out", default="out_rh_both", help="Ausgabeordner")
    ap.add_argument("--Tlist", default="273,400,1000", help="Komma-Liste für Planck-Overlays")
    ap.add_argument("--do_heatmap", action="store_true", help="Score(T,N)-Heatmaps erzeugen")
    ap.add_argument("--Tmin", type=float, default=100.0)
    ap.add_argument("--Tmax", type=float, default=2000.0)
    ap.add_argument("--Tpoints", type=int, default=260)
    ap.add_argument("--Nlist", default="200,500,1000,2000,5000,10000")
    ap.add_argument("--detrend", type=int, default=5)
    args = ap.parse_args()

    outdir = os.path.abspath(args.out)
    os.makedirs(outdir, exist_ok=True)

    zeros = load_zeros(args.zeros, nmax=args.nzeros)
    T_list = [float(x.strip()) for x in args.Tlist.split(",") if x.strip()]

    # Always produce classic Planck plots (both axes)
    plot_planck_both(T_list, outdir)

    # If spectrum provided: produce overlays and optional heatmaps
    if args.spectrum:
        axis, x, S = load_spectrum(args.spectrum)
        print(f"Spectrum loaded: axis={axis}, points={len(x)}")
        # Save a quick overlay in the same axis for sanity
        if axis == "nu":
            plt.figure(figsize=(11,6))
            plt.plot(x/1e12, S, lw=1.2, label="S_obs")
            for T in T_list:
                plt.plot(x/1e12, planck_u_nu(x, T), lw=1.0, ls="--", alpha=0.7, label=f"Planck(T={T}K)")
            plt.yscale("log")
            plt.xlabel("Frequenz ν [THz]")
            plt.ylabel("Spektrum (log)")
            plt.title("S_obs vs Planck (ν-Darstellung)")
            plt.grid(True, alpha=0.25, which="both")
            plt.legend(ncols=2, fontsize=9)
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, "spectrum_overlay_nu.png"), dpi=170, bbox_inches="tight")
            plt.close()
        else:
            lam_um = x * 1e6
            plt.figure(figsize=(11,6))
            plt.plot(lam_um, S, lw=1.2, label="S_obs")
            for T in T_list:
                plt.plot(lam_um, planck_u_lambda(x, T), lw=1.0, ls="--", alpha=0.7, label=f"Planck(T={T}K)")
            plt.yscale("log")
            plt.xlabel("Wellenlänge λ [µm]")
            plt.ylabel("Spektrum (log)")
            plt.title("S_obs vs Planck (λ-Darstellung)")
            plt.grid(True, alpha=0.25, which="both")
            plt.legend(ncols=2, fontsize=9)
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, "spectrum_overlay_lambda.png"), dpi=170, bbox_inches="tight")
            plt.close()

        if args.do_heatmap:
            T_grid = np.linspace(args.Tmin, args.Tmax, int(args.Tpoints), dtype=np.float64)
            N_list = [int(x.strip()) for x in args.Nlist.split(",") if x.strip()]
            S_mat = score_heatmap_mult_u(axis, x, S, zeros, T_grid, N_list, detrend_deg=int(args.detrend))
            outpath = os.path.join(outdir, f"score_heatmap_{axis}.png")
            plot_heatmap(S_mat, T_grid, N_list,
                         title=f"Score(T,N) – {axis}-Achse (mult_u, inv_sqrt, detrend={args.detrend})",
                         outpath=outpath)
            # Save matrix as csv
            np.savetxt(os.path.join(outdir, f"score_heatmap_{axis}.csv"),
                       np.column_stack([T_grid] + [S_mat[j,:] for j in range(S_mat.shape[0])]).T,
                       delimiter=",")
            print(f"Heatmap saved: {outpath}")

    print(f"Done. Output: {outdir}")

if __name__ == "__main__":
    main()
