#!/usr/bin/env python3
"""
Fibonacci-Zeta: drei numerische Zeugen (kontrolliertes Gegenmodell)
====================================================================

Zeuge 1 — Goldene Resonanzpunkte t_k = 2πk/log φ:
  |ζ(1/2+it_k)| vs |ζ_F(1/2+it_k)| und zufällige Höhen auf der kritischen Linie.

Zeuge 2 — Prim-Frequenzen auf der goldenen Achse:
  θ_p = (log p)/(log φ) mod 1, Gleichverteilungstest.

Zeuge 3 — EABC-Fibonacci-Fenster F_k ≤ p < F_{k+1}:
  ABCE- vs. CEAB-Viererfolgen (mod 12) auf aufeinanderfolgenden Primzahlen, D_k.

Referenz: collatz_eabc_zeta_exponential_gedankenexperiment.md §9
Label: Experiment (Schicht C) — kein RH-Anspruch, kein EABC-Theorem.

Ausgabe: eabc_zeta_fibonacci_witnesses.json (+ stdout)
"""

from __future__ import annotations

import argparse
import json
import math
from math import isqrt
from pathlib import Path

import numpy as np

try:
    import mpmath as mp

    _HAS_MPMATH = True
except ImportError:
    mp = None  # type: ignore[assignment]
    _HAS_MPMATH = False

PHI = (1.0 + math.sqrt(5.0)) / 2.0
LOG_PHI = math.log(PHI)
PSI = -1.0 / PHI

RESIDUE_MAP: dict[int, str] = {1: "E", 5: "A", 7: "B", 11: "C"}
ABCE_PATTERN = ("A", "B", "C", "E")
CEAB_PATTERN = ("C", "E", "A", "B")

DEFAULT_OUTPUT = Path("eabc_zeta_fibonacci_witnesses.json")


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

def sieve_primes(limit: int) -> list[int]:
    if limit < 2:
        return []
    flags = bytearray([1]) * (limit + 1)
    flags[0] = flags[1] = 0
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p * p :: p] = bytearray(len(flags[p * p :: p]))
    return [i for i, ok in enumerate(flags) if ok]


def fibonacci_upto(n_max: int) -> list[int]:
    """F_0 … F_{n_max} (iterativ)."""
    if n_max < 0:
        return []
    fib = [0, 1]
    while len(fib) <= n_max:
        fib.append(fib[-1] + fib[-2])
    return fib[: n_max + 1]


def fibonacci_below(bound: int) -> list[int]:
    """Alle F_k mit F_k < bound (k ≥ 0)."""
    fib = [0, 1]
    while fib[-1] < bound:
        fib.append(fib[-1] + fib[-2])
    return fib


def mod12_label(p: int) -> str | None:
    return RESIDUE_MAP.get(p % 12)


def zeta_riemann(s: complex) -> complex:
    if not _HAS_MPMATH:
        raise RuntimeError("mpmath nicht verfügbar — ζ(s) nur mit mpmath")
    with mp.workdps(50):
        return complex(mp.zeta(s))


def zeta_fibonacci(s: complex, n_terms: int = 80) -> complex:
    """ζ_F(s) = Σ_{n≥1} F_n^{-s} (endliche Partialsumme)."""
    total = 0.0 + 0.0j
    for n in range(1, n_terms + 1):
        fn = round((PHI**n - PSI**n) / math.sqrt(5))
        if fn <= 0:
            continue
        total += fn ** (-s)
    return total


def fourgram_pattern(primes: list[int], i: int) -> str | None:
    labels = [mod12_label(primes[i + j]) for j in range(4)]
    if any(l is None for l in labels):
        return None
    tup = tuple(labels)  # type: ignore[arg-type]
    if tup == ABCE_PATTERN:
        return "ABCE"
    if tup == CEAB_PATTERN:
        return "CEAB"
    return "other"


# ---------------------------------------------------------------------------
# Zeuge 1
# ---------------------------------------------------------------------------

def witness_golden_resonance(
    k_max: int,
    n_random: int,
    t_max: float,
    seed: int,
    zeta_f_terms: int,
) -> dict:
    t_res = [2.0 * math.pi * k / LOG_PHI for k in range(1, k_max + 1)]
    rng = np.random.default_rng(seed)
    t_rand = rng.uniform(0.0, t_max, size=n_random)

    zeta_res: list[float] = []
    zeta_f_res: list[float] = []
    zeta_rand: list[float] = []
    zeta_f_rand: list[float] = []

    backend = "mpmath" if _HAS_MPMATH else "unavailable"
    if _HAS_MPMATH:
        for t in t_res:
            s = 0.5 + 1j * t
            zeta_res.append(abs(zeta_riemann(s)))
            zeta_f_res.append(abs(zeta_fibonacci(s, zeta_f_terms)))
        for t in t_rand:
            s = 0.5 + 1j * float(t)
            zeta_rand.append(abs(zeta_riemann(s)))
            zeta_f_rand.append(abs(zeta_fibonacci(s, zeta_f_terms)))

    def summary(vals: list[float], ref: list[float] | None = None) -> dict:
        if not vals:
            return {"count": 0}
        arr = np.array(vals)
        out: dict = {
            "count": len(vals),
            "mean": float(np.mean(arr)),
            "max": float(np.max(arr)),
            "min": float(np.min(arr)),
        }
        if ref:
            ref_arr = np.array(ref)
            out["mean_ratio_vs_random"] = float(np.mean(arr) / np.mean(ref_arr))
            out["max_percentile_among_random"] = float(
                np.mean(ref_arr <= np.max(arr))
            )
        return out

    rows = []
    for k, t in enumerate(t_res, start=1):
        row: dict = {"k": k, "t_k": t}
        if _HAS_MPMATH:
            row["abs_zeta"] = zeta_res[k - 1]
            row["abs_zeta_F"] = zeta_f_res[k - 1]
        rows.append(row)

    return {
        "label": "Experiment (Schicht C)",
        "description": "Goldene Resonanzpunkte t_k=2πk/log φ auf s=1/2+it",
        "backend_zeta": backend,
        "k_max": k_max,
        "n_random": n_random,
        "t_max": t_max,
        "resonance_points": rows,
        "zeta_summary": {
            "at_resonance": summary(zeta_res, zeta_rand),
            "at_random": summary(zeta_rand),
        },
        "zeta_F_summary": {
            "at_resonance": summary(zeta_f_res, zeta_f_rand),
            "at_random": summary(zeta_f_rand),
            "n_terms": zeta_f_terms,
        },
        "note": (
            "Resonanzen von ζ_F liegen bei s=2πik/log φ (imaginäre Achse); "
            "Vergleich auf Re(s)=1/2 ist diagnostisch, kein RH-Anspruch."
        ),
    }


# ---------------------------------------------------------------------------
# Zeuge 2
# ---------------------------------------------------------------------------

def witness_prime_phases(prime_bound: int, n_bins: int) -> dict:
    primes = [p for p in sieve_primes(prime_bound) if p > 3]
    theta = np.array([math.log(p) / LOG_PHI % 1.0 for p in primes])

    hist, edges = np.histogram(theta, bins=n_bins, range=(0.0, 1.0))
    expected = len(theta) / n_bins
    chi2 = float(np.sum((hist - expected) ** 2 / expected)) if expected > 0 else float("nan")
    chi2_p = None
    try:
        from scipy import stats

        chi2_p = float(stats.chisquare(hist).pvalue)
    except ImportError:
        pass

    return {
        "label": "Experiment (Schicht C)",
        "description": "θ_p = (log p)/(log φ) mod 1 — Fibonacci-Sampling der Prim-Frequenzen",
        "prime_bound": prime_bound,
        "n_primes": len(primes),
        "n_bins": n_bins,
        "histogram_counts": hist.tolist(),
        "bin_edges": edges.tolist(),
        "mean_theta": float(np.mean(theta)),
        "std_theta": float(np.std(theta)),
        "chi2_vs_uniform": chi2,
        "chi2_pvalue": chi2_p,
        "expected_per_bin": expected,
        "note": "Gleichverteilung wäre ein Nullbefund; Cluster wären arithmetische Struktur.",
    }


# ---------------------------------------------------------------------------
# Zeuge 3
# ---------------------------------------------------------------------------

def witness_fibonacci_windows(prime_bound: int) -> dict:
    primes = [p for p in sieve_primes(prime_bound) if p > 3]
    fib = fibonacci_below(prime_bound)
    # Fenster F_k ≤ p < F_{k+1} für k ≥ 4 (F_4=3)
    windows: list[dict] = []

    for k in range(4, len(fib) - 1):
        lo, hi = fib[k], fib[k + 1]
        if lo >= prime_bound:
            break
        in_win = [p for p in primes if lo <= p < hi]
        if len(in_win) < 4:
            windows.append(
                {
                    "k": k,
                    "F_k": lo,
                    "F_k1": hi,
                    "n_primes": len(in_win),
                    "ABCE": 0,
                    "CEAB": 0,
                    "other": 0,
                    "D_k": None,
                }
            )
            continue

        abce = ceab = other = 0
        for i in range(len(in_win) - 3):
            pat = fourgram_pattern(in_win, i)
            if pat == "ABCE":
                abce += 1
            elif pat == "CEAB":
                ceab += 1
            elif pat == "other":
                other += 1

        total = abce + ceab
        d_k = (abce - ceab) / total if total > 0 else None

        windows.append(
            {
                "k": k,
                "F_k": lo,
                "F_k1": hi,
                "n_primes": len(in_win),
                "ABCE": abce,
                "CEAB": ceab,
                "other": other,
                "D_k": d_k,
            }
        )

    d_vals = [w["D_k"] for w in windows if w["D_k"] is not None]
    return {
        "label": "Experiment (Schicht C)",
        "description": "ABCE vs CEAB in Fibonacci-Fenstern auf aufeinanderfolgenden Prim-4-Grammen (mod 12)",
        "prime_bound": prime_bound,
        "n_windows": len(windows),
        "windows": windows,
        "aggregate": {
            "total_ABCE": sum(w["ABCE"] for w in windows),
            "total_CEAB": sum(w["CEAB"] for w in windows),
            "mean_D_k": float(np.mean(d_vals)) if d_vals else None,
            "std_D_k": float(np.std(d_vals)) if d_vals else None,
        },
        "note": "D_k = (ABCE−CEAB)/(ABCE+CEAB); Wigner-Zellen optional — hier minimale 4-Gramm-Zählung.",
    }


# ---------------------------------------------------------------------------
# Hauptprogramm
# ---------------------------------------------------------------------------

def run_witnesses(
    k_max: int = 12,
    n_random: int = 200,
    t_max: float = 80.0,
    prime_bound: int = 50_000,
    n_bins: int = 20,
    seed: int = 42,
    zeta_f_terms: int = 80,
) -> dict:
    return {
        "epistemic_label": "Experiment (Schicht C) — kontrolliertes Gegenmodell, kein Theorem",
        "phi": PHI,
        "log_phi": LOG_PHI,
        "witness_1_golden_resonance": witness_golden_resonance(
            k_max, n_random, t_max, seed, zeta_f_terms
        ),
        "witness_2_prime_phases": witness_prime_phases(prime_bound, n_bins),
        "witness_3_fibonacci_windows": witness_fibonacci_windows(prime_bound),
        "references": [
            "collatz_eabc_zeta_exponential_gedankenexperiment.md §9",
            "collatz_eabc_zirkulationshypothese.md §4.8.2 (Stufe 0: regulärer Kamm)",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Fibonacci-Zeta numerische Zeugen")
    parser.add_argument("--k-max", type=int, default=12)
    parser.add_argument("--n-random", type=int, default=200)
    parser.add_argument("--t-max", type=float, default=80.0)
    parser.add_argument("--prime-bound", type=int, default=50_000)
    parser.add_argument("--n-bins", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--zeta-f-terms", type=int, default=80)
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    result = run_witnesses(
        k_max=args.k_max,
        n_random=args.n_random,
        t_max=args.t_max,
        prime_bound=args.prime_bound,
        n_bins=args.n_bins,
        seed=args.seed,
        zeta_f_terms=args.zeta_f_terms,
    )

    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")

    w1 = result["witness_1_golden_resonance"]
    w2 = result["witness_2_prime_phases"]
    w3 = result["witness_3_fibonacci_windows"]

    print("Fibonacci-Zeta Zeugen (Schicht C, Experiment)")
    print(f"  mpmath: {_HAS_MPMATH}")
    if _HAS_MPMATH:
        zr = w1["zeta_summary"]["at_resonance"]
        zz = w1["zeta_summary"]["at_random"]
        print(
            f"  Zeuge 1: mean|ζ| resonance={zr.get('mean', 'n/a'):.4f} "
            f"random={zz.get('mean', 'n/a'):.4f} "
            f"ratio={zr.get('mean_ratio_vs_random', 'n/a')}"
        )
    else:
        print("  Zeuge 1: ζ(s) übersprungen (mpmath fehlt)")

    print(
        f"  Zeuge 2: n_primes={w2['n_primes']} chi2={w2['chi2_vs_uniform']:.2f} "
        f"p={w2['chi2_pvalue']}"
    )
    agg = w3["aggregate"]
    print(
        f"  Zeuge 3: windows={w3['n_windows']} ABCE={agg['total_ABCE']} "
        f"CEAB={agg['total_CEAB']} mean D_k={agg['mean_D_k']}"
    )
    print(f"  → {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
