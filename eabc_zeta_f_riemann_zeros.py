#!/usr/bin/env python3
"""
ζ_F an Riemann-Nullstellen-Ordinaten γ_n
=========================================

Zwei Einsetzungen:
  • reell:     s = γ_n           (trivial: F_1=F_2=1 ⇒ ζ_F(γ_n) ≈ 2)
  • kritisch:  s = 1/2 + i γ_n   (Beinahe-Auslöschungen, Resonanztest)

Reihe (stabil, ohne Resonatorformel):
  ζ_F(s) = Σ_{k≥1} F_k^{-s}

Optional: R_n = min_{t ∈ K_φ} |γ_n − t|,  K_φ = {2πk / log φ : k ∈ ℤ}.

Referenz: collatz_eabc_zeta_exponential_gedankenexperiment.md §9.12–9.13;
         collatz_eabc_zirkulationshypothese.md §4.8.6 (Gesamtbaustein)
Label: Experiment (Schicht C) — kein RH-Anspruch, kein Beweis.

CLI:
  python eabc_zeta_f_riemann_zeros.py --zeros-file zeros6.npz --limit 100000 -o out.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

PHI = (1.0 + math.sqrt(5.0)) / 2.0
LOG_PHI = math.log(PHI)

DEFAULT_ZEROS = Path("zeros6.npz")
DEFAULT_OUTPUT = Path("eabc_zeta_f_riemann_zeros.json")


# ---------------------------------------------------------------------------
# Nullstellen laden
# ---------------------------------------------------------------------------

def load_zeros(path: Path) -> np.ndarray:
    """Lädt γ_n (Imaginärteile auf Re(s)=1/2), aufsteigend sortiert."""
    if not path.exists():
        raise FileNotFoundError(f"Nullstellen-Datei nicht gefunden: {path}")
    raw = np.load(path, allow_pickle=False)
    if hasattr(raw, "files"):
        key = "zeros" if "zeros" in raw.files else raw.files[0]
        gamma = raw[key]
    else:
        gamma = raw
    gamma = np.asarray(gamma, dtype=np.float64).ravel()
    if gamma.size == 0:
        raise ValueError("leere Nullstellen-Datei")
    if not np.all(np.diff(gamma) >= 0):
        gamma = np.sort(gamma)
    return gamma


# ---------------------------------------------------------------------------
# ζ_F(s) = Σ F_k^{-s}
# ---------------------------------------------------------------------------

def _fibonacci_int(n: int) -> int:
    """Exaktes F_n (n ≥ 0), iterativ."""
    if n < 0:
        raise ValueError("n < 0")
    if n == 0:
        return 0
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def zeta_f_series(
    s: complex,
    *,
    rtol: float = 1e-14,
    atol: float = 1e-18,
    max_terms: int = 300,
) -> tuple[complex, int]:
    """
    Partialsumme Σ_{k≥1} F_k^{-s} mit adaptiver Trunkierung.

    Abbruch wenn |Term| ≤ max(rtol·|Summe|, atol).
    Für große |Im(s)| auf Re(s)=1/2 konvergiert |F_k^{-s}| = F_k^{-1/2} geometrisch.
    """
    total = 0.0 + 0.0j
    f_prev, f_cur = 1, 1  # F_1, F_2
    used = 0
    for k in range(1, max_terms + 1):
        if k == 1:
            fk = 1
        elif k == 2:
            fk = 1
        else:
            f_prev, f_cur = f_cur, f_prev + f_cur
            fk = f_cur
        term = fk ** (-s)
        total += term
        used = k
        bound = max(atol, rtol * abs(total))
        if abs(term) <= bound:
            break
    return total, used


def zeta_f_batch_real(gamma: np.ndarray, **kw: Any) -> tuple[np.ndarray, np.ndarray]:
    """|ζ_F(γ_n)| für reelle Einsetzung."""
    out = np.empty(gamma.size, dtype=np.float64)
    terms = np.empty(gamma.size, dtype=np.int32)
    for i, g in enumerate(gamma):
        val, n = zeta_f_series(float(g), **kw)
        out[i] = abs(val)
        terms[i] = n
    return out, terms


def zeta_f_batch_critical(gamma: np.ndarray, **kw: Any) -> tuple[np.ndarray, np.ndarray]:
    """|ζ_F(1/2 + i γ_n)|."""
    out = np.empty(gamma.size, dtype=np.float64)
    terms = np.empty(gamma.size, dtype=np.int32)
    for i, g in enumerate(gamma):
        val, n = zeta_f_series(0.5 + 1j * float(g), **kw)
        out[i] = abs(val)
        terms[i] = n
    return out, terms


# ---------------------------------------------------------------------------
# Goldener Kamm K_φ = {2πk / log φ}
# ---------------------------------------------------------------------------

def distance_to_phi_comb(gamma: float, k_span: int | None = None) -> tuple[float, int]:
    """
    R_n = min_{k∈ℤ} |γ − 2πk/log φ|.
    k_span: halbe Breite der k-Suche; default aus γ-Schätzung.
    """
    if k_span is None:
        k_center = gamma * LOG_PHI / (2.0 * math.pi)
        k_span = max(3, int(k_center) + 2)
    k_center = int(round(gamma * LOG_PHI / (2.0 * math.pi)))
    best_d = float("inf")
    best_k = 0
    for k in range(k_center - k_span, k_center + k_span + 1):
        t = 2.0 * math.pi * k / LOG_PHI
        d = abs(gamma - t)
        if d < best_d:
            best_d = d
            best_k = k
    return best_d, best_k


def phi_comb_distances(gamma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    d = np.empty(gamma.size, dtype=np.float64)
    k_idx = np.empty(gamma.size, dtype=np.int64)
    for i, g in enumerate(gamma):
        d[i], k_idx[i] = distance_to_phi_comb(float(g))
    return d, k_idx


# ---------------------------------------------------------------------------
# Statistik
# ---------------------------------------------------------------------------

def threshold_counts(vals: np.ndarray, thresholds: tuple[float, ...]) -> dict[str, int]:
    return {f"count_lt_{thr:g}": int(np.sum(vals < thr)) for thr in thresholds}


def threshold_comparison(
    riemann_vals: np.ndarray,
    random_vals: np.ndarray,
    thresholds: tuple[float, ...],
) -> dict[str, Any]:
    """
    Vergleich #{A_n < ε}: Riemann-γ vs. gleichverteilte Zufalls-γ (gleiche Dichte).
    """
    n_r = riemann_vals.size
    n_z = random_vals.size
    rows: list[dict[str, Any]] = []
    for thr in thresholds:
        c_r = int(np.sum(riemann_vals < thr))
        c_z = int(np.sum(random_vals < thr))
        rows.append(
            {
                "epsilon": thr,
                "count_riemann": c_r,
                "count_random": c_z,
                "excess_riemann_minus_random": c_r - c_z,
                "rate_riemann": c_r / n_r if n_r else 0.0,
                "rate_random": c_z / n_z if n_z else 0.0,
            }
        )
    return {
        "description": "#{|ζ_F(1/2+iγ)| < ε}: Riemann-γ_n vs. Zufalls-γ gleicher Dichte",
        "thresholds": rows,
    }


def quantile_summary(vals: np.ndarray) -> dict[str, float]:
    qs = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
    return {
        "min": float(np.min(vals)),
        "max": float(np.max(vals)),
        "mean": float(np.mean(vals)),
        "median": float(np.median(vals)),
        "std": float(np.std(vals)),
        **{f"q{int(100*q):02d}": float(np.quantile(vals, q)) for q in qs},
    }


def random_baseline(
    gamma_lo: float,
    gamma_hi: float,
    n: int,
    seed: int,
    **kw: Any,
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    t_rand = rng.uniform(gamma_lo, gamma_hi, size=n)
    abs_crit, _ = zeta_f_batch_critical(t_rand, **kw)
    return {
        "n": n,
        "gamma_range": [gamma_lo, gamma_hi],
        "seed": seed,
        "gamma_samples": t_rand,
        "abs_critical": abs_crit,
        "critical": quantile_summary(abs_crit),
        "thresholds": threshold_counts(abs_crit, (0.005, 0.01, 0.05, 0.1, 1.0)),
    }


# ---------------------------------------------------------------------------
# Hauptlauf
# ---------------------------------------------------------------------------

def run_analysis(
    zeros_path: Path,
    limit: int,
    *,
    seed: int = 42,
    rtol: float = 1e-14,
    max_terms: int = 300,
    compute_phi_comb: bool = True,
    random_baseline_n: int | None = None,
    detail_limit: int = 20,
) -> dict[str, Any]:
    gamma_full = load_zeros(zeros_path)
    n = min(limit, gamma_full.size)
    gamma = gamma_full[:n]
    kw = {"rtol": rtol, "max_terms": max_terms}

    abs_real, terms_real = zeta_f_batch_real(gamma, **kw)
    abs_crit, terms_crit = zeta_f_batch_critical(gamma, **kw)

    min_idx = int(np.argmin(abs_crit))
    thresholds = (0.005, 0.01, 0.05, 0.1, 1.0)

    detail_rows: list[dict[str, Any]] = []
    for i in range(min(detail_limit, n)):
        s_real = float(gamma[i])
        z_real, _ = zeta_f_series(s_real, **kw)
        z_crit, _ = zeta_f_series(0.5 + 1j * s_real, **kw)
        row: dict[str, Any] = {
            "n": i + 1,
            "gamma": s_real,
            "zeta_F_real": {"re": z_real.real, "im": z_real.imag, "abs": abs(z_real)},
            "zeta_F_critical": {"re": z_crit.real, "im": z_crit.imag, "abs": abs(z_crit)},
            "terms_real": int(terms_real[i]),
            "terms_critical": int(terms_crit[i]),
        }
        if compute_phi_comb:
            r_n, k_n = distance_to_phi_comb(s_real)
            row["R_n"] = r_n
            row["k_nearest"] = k_n
            row["t_k_nearest"] = 2.0 * math.pi * k_n / LOG_PHI
        detail_rows.append(row)

    result: dict[str, Any] = {
        "label": "Experiment (Schicht C) — Resonanztest, kein Beweis",
        "zeros_file": str(zeros_path),
        "n_zeros": n,
        "gamma_range": [float(gamma[0]), float(gamma[-1])],
        "series": "zeta_F(s) = sum_{k>=1} F_k^{-s}",
        "rtol": rtol,
        "max_terms": max_terms,
        "real_insertion": {
            "description": "s = gamma_n (reell); F_1=F_2=1 ⇒ Summe startet bei 2",
            "summary": quantile_summary(abs_real),
            "thresholds": threshold_counts(abs_real, thresholds),
            "max_terms_used": int(np.max(terms_real)),
        },
        "critical_insertion": {
            "description": "s = 1/2 + i gamma_n",
            "summary": quantile_summary(abs_crit),
            "thresholds": threshold_counts(abs_crit, thresholds),
            "minimum": {
                "index": min_idx + 1,
                "gamma": float(gamma[min_idx]),
                "abs_zeta_F": float(abs_crit[min_idx]),
            },
            "max_terms_used": int(np.max(terms_crit)),
        },
        "first_zeros_detail": detail_rows,
    }

    if compute_phi_comb:
        r_all, _ = phi_comb_distances(gamma)
        result["phi_comb"] = {
            "description": "K_phi = {2 pi k / log phi}; R_n = min |gamma_n - t|",
            "summary": quantile_summary(r_all),
        }

    n_rand = random_baseline_n if random_baseline_n is not None else n
    rand_result = random_baseline(
        float(gamma[0]),
        float(gamma[-1]),
        n_rand,
        seed,
        **kw,
    )
    # JSON-serialisierbar: große Arrays weglassen
    rand_json = {k: v for k, v in rand_result.items() if k not in ("gamma_samples", "abs_critical")}
    result["random_baseline"] = rand_json
    result["threshold_comparison"] = threshold_comparison(
        abs_crit,
        rand_result["abs_critical"],
        thresholds,
    )

    return result


def write_csv(result: dict[str, Any], path: Path) -> None:
    rows = result.get("first_zeros_detail", [])
    if not rows:
        path.write_text("n,gamma,abs_real,abs_critical,R_n\n")
        return
    fieldnames = [
        "n",
        "gamma",
        "abs_real",
        "re_real",
        "im_real",
        "abs_critical",
        "re_critical",
        "im_critical",
        "R_n",
        "k_nearest",
        "terms_real",
        "terms_critical",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(
                {
                    "n": r["n"],
                    "gamma": r["gamma"],
                    "abs_real": r["zeta_F_real"]["abs"],
                    "re_real": r["zeta_F_real"]["re"],
                    "im_real": r["zeta_F_real"]["im"],
                    "abs_critical": r["zeta_F_critical"]["abs"],
                    "re_critical": r["zeta_F_critical"]["re"],
                    "im_critical": r["zeta_F_critical"]["im"],
                    "R_n": r.get("R_n"),
                    "k_nearest": r.get("k_nearest"),
                    "terms_real": r["terms_real"],
                    "terms_critical": r["terms_critical"],
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="ζ_F an Riemann-Nullstellen-Ordinaten (reell vs. kritisch)"
    )
    parser.add_argument(
        "--zeros-file",
        type=Path,
        default=DEFAULT_ZEROS,
        help="zeros6.npz / zeros6.npy / Textdatei mit γ pro Zeile",
    )
    parser.add_argument("--limit", type=int, default=100_000, help="erste N Ordinaten")
    parser.add_argument("--detail-limit", type=int, default=20, help="Tabellenzeilen im Detail")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--rtol", type=float, default=1e-14)
    parser.add_argument("--max-terms", type=int, default=300)
    parser.add_argument(
        "--no-phi-comb",
        action="store_true",
        help="R_n zum goldenen Kamm weglassen",
    )
    parser.add_argument(
        "--random-n",
        type=int,
        default=None,
        help="Stichprobengröße Zufalls-Baseline (default: wie --limit)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Ausgabe .json oder .csv",
    )
    args = parser.parse_args()

    # Textdatei: eine γ pro Zeile
    zeros_path = args.zeros_file
    if zeros_path.suffix in {".txt", ".dat"}:
        gamma = np.loadtxt(zeros_path, dtype=np.float64)
        tmp = Path("_zeta_f_zeros_tmp.npy")
        np.save(tmp, gamma)
        zeros_path = tmp

    result = run_analysis(
        zeros_path,
        args.limit,
        seed=args.seed,
        rtol=args.rtol,
        max_terms=args.max_terms,
        compute_phi_comb=not args.no_phi_comb,
        random_baseline_n=args.random_n,
        detail_limit=args.detail_limit,
    )

    out = args.output
    if out.suffix.lower() == ".csv":
        write_csv(result, out)
    else:
        out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")

    cr = result["critical_insertion"]
    mn = cr["minimum"]
    sm = cr["summary"]
    rb = result["random_baseline"]["critical"]
    print("ζ_F an Riemann-Nullstellen (Schicht C)")
    print(f"  N={result['n_zeros']}  γ∈[{result['gamma_range'][0]:.2f}, {result['gamma_range'][1]:.2f}]")
    print(f"  reell:  median|ζ_F(γ)|={result['real_insertion']['summary']['median']:.6g}")
    print(
        f"  kritisch: min={mn['abs_zeta_F']:.8g} @ n={mn['index']} γ={mn['gamma']:.6g} "
        f"median={sm['median']:.4f}"
    )
    print(
        f"  Zufall: min={rb['min']:.8g} median={rb['median']:.4f} "
        f"(n={result['random_baseline']['n']})"
    )
    tc = result.get("threshold_comparison", {}).get("thresholds", [])
    if tc:
        print("  #{|ζ_F(1/2+iγ)|<ε}  Riemann vs. Zufall:")
        for row in tc:
            print(
                f"    ε={row['epsilon']:g}: "
                f"{row['count_riemann']} vs {row['count_random']} "
                f"(Δ={row['excess_riemann_minus_random']:+d})"
            )
    print(f"  → {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
