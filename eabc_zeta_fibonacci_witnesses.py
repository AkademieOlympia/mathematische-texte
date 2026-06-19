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

Zeuge 4 — Meromorphe Normalform von ζ_F(s) (klassisch):
  Binomialentwicklung von (1-u)^{-s}, u=(-1)^n φ^{-2n}; Resonanztürme.

Zeuge 5 — Goldene Fourier-Zeugen auf Primvierlingen:
  θ_φ(p)=(log p)/(log φ) mod 1, χ(p)=±1 (ABCE/CEAB), C_m(N), Z_m(N)=|C_m|/√Q.

Zeuge 6 — Fibonacci-Schalen mit mod-210-Triple (11,101,191) vs. lineares Fenster:
  Signatur (+,+,-) = (D_11>0, D_101>0, D_191<0) entlang F_k ≤ p < F_{k+1}.

Referenz: collatz_eabc_zeta_exponential_gedankenexperiment.md §9.8–9.9
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

# mod-210-Wigner-Zellen (p mod 420) — vgl. eabc_wigner_zellen.py, eabc_quadruplets.csv
MOD210_CELLS = (11, 101, 191)
ZELLEN_MOD420: dict[int, tuple[int, int]] = {
    11: (221, 11),
    101: (101, 311),
    191: (401, 191),
}
SIGNATURE_TARGET = ("+", "+", "-")  # (D_11, D_101, D_191)

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


def rising_factorial(s: complex, m: int) -> complex:
    """Pochhammer (s)_m = s(s+1)…(s+m-1)."""
    if m == 0:
        return 1.0 + 0.0j
    if _HAS_MPMATH:
        with mp.workdps(50):
            return complex(mp.rf(s, m))
    prod = 1.0 + 0.0j
    for j in range(m):
        prod *= s + j
    return prod


def zeta_f_meromorphic(s: complex, m_max: int = 60) -> complex:
    r"""
    Meromorphe Normalform (Partialsumme m=0..m_max):
    ζ_F(s) = 5^{s/2} Σ_m (s)_m/m! · (-1)^m φ^{-(s+2m)} / (1-(-1)^m φ^{-(s+2m)}).
    """
    total = 0.0 + 0.0j
    for m in range(m_max + 1):
        parity = (-1) ** m
        pow_term = PHI ** (-(s + 2 * m))
        denom = 1.0 - parity * pow_term
        if abs(denom) < 1e-14:
            continue
        coeff = rising_factorial(s, m) / math.factorial(m) * parity
        total += coeff * pow_term / denom
    return (5.0 ** (s / 2)) * total


def find_quadruplets(limit: int) -> np.ndarray:
    """Startprimzahlen p aller Primvierlinge (p,p+2,p+6,p+8) mit p ≤ limit."""
    if limit < 11:
        return np.array([], dtype=np.int64)
    flags = bytearray([1]) * (limit + 9)
    flags[0] = flags[1] = 0
    for p in range(2, isqrt(limit + 8) + 1):
        if flags[p]:
            flags[p * p :: p] = bytearray(len(flags[p * p :: p]))
    lim = limit - 8
    quads: list[int] = []
    for p in range(2, lim + 1):
        if flags[p] and flags[p + 2] and flags[p + 6] and flags[p + 8]:
            quads.append(p)
    return np.array(quads, dtype=np.int64)


def chi_quadruplet(p: int) -> int:
    """χ(p)=+1 (ABCE, p≡5 mod 12), -1 (CEAB, p≡11 mod 12), sonst 0."""
    r = p % 12
    if r == 5:
        return 1
    if r == 11:
        return -1
    return 0


def theta_phi(p: int) -> float:
    return math.log(p) / LOG_PHI % 1.0


def classify_mod210(p: int) -> tuple[int, str] | None:
    """Klassifiziert Vierling p nach mod-210-Zelle und ABCE/CEAB."""
    r = p % 210
    if r not in ZELLEN_MOD420:
        return None
    mod420 = p % 420
    abce_val, ceab_val = ZELLEN_MOD420[r]
    if mod420 == abce_val:
        return r, "ABCE"
    if mod420 == ceab_val:
        return r, "CEAB"
    return None


def mod210_triple(counts: dict[int, dict[str, int]]) -> tuple[int, int, int]:
    """(D_11, D_101, D_191) = ABCE−CEAB je Zelle."""
    return tuple(
        counts[r]["ABCE"] - counts[r]["CEAB"] for r in MOD210_CELLS
    )  # type: ignore[return-value]


def signature_label(d11: int, d101: int, d191: int) -> tuple[str, str, str]:
    def sgn(x: int) -> str:
        if x > 0:
            return "+"
        if x < 0:
            return "-"
        return "0"

    return sgn(d11), sgn(d101), sgn(d191)


def matches_signature(d11: int, d101: int, d191: int) -> bool:
    return signature_label(d11, d101, d191) == SIGNATURE_TARGET


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
# Zeuge 4 — meromorphe Normalform
# ---------------------------------------------------------------------------

def witness_meromorphic_normal_form(
    test_points: list[complex] | None = None,
    m_max: int = 50,
    n_terms: int = 80,
) -> dict:
    if test_points is None:
        test_points = [2.0 + 0.0j, 1.5 + 0.3j, 0.5 + 1.0j, 3.0 + 2.0j]

    rows: list[dict] = []
    for s in test_points:
        direct = zeta_fibonacci(s, n_terms)
        merom = zeta_f_meromorphic(s, m_max)
        rel_err = abs(direct - merom) / max(abs(direct), 1e-15)
        rows.append(
            {
                "s": {"re": s.real, "im": s.imag},
                "zeta_F_direct": {"re": direct.real, "im": direct.imag},
                "zeta_F_meromorphic": {"re": merom.real, "im": merom.imag},
                "relative_error": float(rel_err),
            }
        )

    log_phi = LOG_PHI
    resonance_towers = [
        {
            "tower": "main (m=0)",
            "poles": "s = 2πik/log φ, k ∈ ℤ\\{0}",
            "examples_im": [2 * math.pi * k / log_phi for k in (1, 2, 3)],
        },
        {
            "tower": "even m",
            "poles": "s = -2m + 2πik/log φ, m ≥ 1",
            "examples_im": [-2 * m + 2 * math.pi * k / log_phi for m in (1, 2) for k in (1,)],
        },
        {
            "tower": "odd m",
            "poles": "s = -2m + (2k+1)πi/log φ, m ≥ 1",
            "examples_im": [
                -2 * m + (2 * k + 1) * math.pi / log_phi for m in (1,) for k in (0, 1)
            ],
        },
    ]

    return {
        "label": "Theorem (klassisch) — numerische Verifikation",
        "description": "Meromorphe Normalform via Binomialentwicklung (1-u)^{-s}, u=(-1)^n φ^{-2n}",
        "formula": (
            "ζ_F(s)=5^{s/2} Σ_{m≥0} (s)_m/m! · (-1)^m φ^{-(s+2m)} / (1-(-1)^m φ^{-(s+2m)})"
        ),
        "m_max": m_max,
        "n_terms_direct": n_terms,
        "resonance_towers": resonance_towers,
        "verification_points": rows,
        "max_relative_error": max(r["relative_error"] for r in rows),
        "note": "Exakte Reihe = klassische Analysis; Abgleich direct vs. meromorph ist Sanity-Check.",
    }


# ---------------------------------------------------------------------------
# Zeuge 5 — goldene Fourier-Zeugen C_m, Z_m
# ---------------------------------------------------------------------------

def compute_Cm_Zm(
    quads: np.ndarray,
    m_max: int,
) -> tuple[list[complex], list[float], int]:
    """C_m und Z_m=|C_m|/√Q für m=0..m_max."""
    if len(quads) == 0:
        return [0.0 + 0.0j] * (m_max + 1), [0.0] * (m_max + 1), 0

    chi_vals = np.array([chi_quadruplet(int(p)) for p in quads], dtype=np.float64)
    theta = np.array([theta_phi(int(p)) for p in quads], dtype=np.float64)
    mask = chi_vals != 0
    chi_vals = chi_vals[mask]
    theta = theta[mask]
    q = int(len(chi_vals))
    sqrt_q = math.sqrt(q) if q > 0 else 1.0

    m_range = np.arange(m_max + 1, dtype=np.float64)
    phases = np.exp(2j * math.pi * np.outer(m_range, theta))  # (m_max+1, q)
    c_m = (phases * chi_vals).sum(axis=1)
    z_m = np.abs(c_m) / sqrt_q
    return c_m.tolist(), z_m.tolist(), q


def witness_golden_fourier(
    quadruplet_bound: int,
    m_max: int = 8,
    checkpoints: list[int] | None = None,
) -> dict:
    if checkpoints is None:
        checkpoints = sorted(
            {min(quadruplet_bound, x) for x in (10_000, 50_000, 200_000, 500_000, quadruplet_bound)}
        )

    quads_all = find_quadruplets(quadruplet_bound)
    checkpoint_rows: list[dict] = []

    for n_cp in checkpoints:
        if n_cp > quadruplet_bound:
            continue
        sub = quads_all[quads_all <= n_cp]
        c_m, z_m, q = compute_Cm_Zm(sub, m_max)
        checkpoint_rows.append(
            {
                "N": n_cp,
                "Q": q,
                "sqrt_Q": math.sqrt(q) if q > 0 else 0.0,
                "C_m": [{"re": z.real, "im": z.imag} for z in c_m],
                "Z_m": [float(z) for z in z_m],
            }
        )

    return {
        "label": "Experiment (Schicht C)",
        "description": (
            "θ_φ(p)=(log p)/(log φ) mod 1, χ(p)=+1 ABCE (p≡5 mod 12), "
            "-1 CEAB (p≡11 mod 12); C_m=Σ χ(p)e^{2πimθ_φ}, Z_m=|C_m|/√Q"
        ),
        "quadruplet_bound": quadruplet_bound,
        "m_max": m_max,
        "checkpoints": checkpoint_rows,
        "note": "EABC-Kopplung auf Primvierlingen = Experiment/Hypothese (Schicht C).",
    }


# ---------------------------------------------------------------------------
# Zeuge 6 — Fibonacci-Schalen mod-210 vs. lineares Fenster
# ---------------------------------------------------------------------------

def _shell_mod210_counts(quads: np.ndarray, lo: int, hi: int) -> dict[int, dict[str, int]]:
    counts = {r: {"ABCE": 0, "CEAB": 0} for r in MOD210_CELLS}
    sub = quads[(quads >= lo) & (quads < hi)]
    for p in sub:
        if int(p) == 5:
            continue
        cls = classify_mod210(int(p))
        if cls is None:
            continue
        r, pat = cls
        counts[r][pat] += 1
    return counts


def witness_fibonacci_mod210_shells(
    quadruplet_bound: int,
    min_quads_per_shell: int = 1,
) -> dict:
    quads = find_quadruplets(quadruplet_bound)
    if len(quads) == 0:
        return {"label": "Experiment (Schicht C)", "n_shells": 0, "shells": []}

    fib = fibonacci_below(int(quads[-1]) + 1)
    shells: list[dict] = []

    for k in range(4, len(fib) - 1):
        lo, hi = fib[k], fib[k + 1]
        if lo >= quadruplet_bound:
            break
        counts = _shell_mod210_counts(quads, lo, hi)
        d11, d101, d191 = mod210_triple(counts)
        q_shell = sum(counts[r]["ABCE"] + counts[r]["CEAB"] for r in MOD210_CELLS)
        if q_shell < min_quads_per_shell:
            continue
        sig = signature_label(d11, d101, d191)
        shells.append(
            {
                "k": k,
                "F_k": lo,
                "F_k1": hi,
                "Q_shell": q_shell,
                "D_11": d11,
                "D_101": d101,
                "D_191": d191,
                "signature": sig,
                "matches_target": matches_signature(d11, d101, d191),
            }
        )

    n_match = sum(1 for s in shells if s["matches_target"])
    frac_match = n_match / len(shells) if shells else None

    # Lineares Fenster-Baseline: gleiche Fensterbreiten wie Fibonacci-Schalen
    linear_windows: list[dict] = []
    p_min = int(quads[0])
    idx = 0
    for sh in shells:
        width = sh["F_k1"] - sh["F_k"]
        lo = p_min + idx * width
        hi = lo + width
        idx += 1
        if hi > quadruplet_bound:
            break
        counts = _shell_mod210_counts(quads, lo, hi)
        d11, d101, d191 = mod210_triple(counts)
        q_win = sum(counts[r]["ABCE"] + counts[r]["CEAB"] for r in MOD210_CELLS)
        linear_windows.append(
            {
                "window_index": idx - 1,
                "lo": lo,
                "hi": hi,
                "width": width,
                "Q_window": q_win,
                "D_11": d11,
                "D_101": d101,
                "D_191": d191,
                "signature": signature_label(d11, d101, d191),
                "matches_target": matches_signature(d11, d101, d191),
            }
        )

    n_lin_match = sum(1 for w in linear_windows if w["matches_target"])
    frac_lin = n_lin_match / len(linear_windows) if linear_windows else None

    return {
        "label": "Experiment (Schicht C)",
        "description": (
            "mod-210-Triple (D_11,D_101,D_191) in Fibonacci-Schalen F_k≤p<F_{k+1}; "
            f"Zielsignatur {SIGNATURE_TARGET}"
        ),
        "quadruplet_bound": quadruplet_bound,
        "signature_target": SIGNATURE_TARGET,
        "n_fibonacci_shells": len(shells),
        "fibonacci_shells": shells,
        "fibonacci_stability": {
            "n_matches": n_match,
            "fraction_matches": frac_match,
        },
        "linear_baseline": {
            "n_windows": len(linear_windows),
            "windows": linear_windows,
            "n_matches": n_lin_match,
            "fraction_matches": frac_lin,
        },
        "aggregate_mod210": _shell_mod210_counts(quads, 0, quadruplet_bound + 1),
        "note": (
            "Primvierlinge (p,p+2,p+6,p+8); Zellen r∈{11,101,191} mod 210 "
            "(vgl. eabc_wigner_zellen.py, eabc_quadruplets.csv)."
        ),
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
    m_max: int = 8,
    quadruplet_bound: int = 500_000,
    meromorphic_m_max: int = 50,
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
        "witness_4_meromorphic_normal_form": witness_meromorphic_normal_form(
            m_max=meromorphic_m_max, n_terms=zeta_f_terms
        ),
        "witness_5_golden_fourier": witness_golden_fourier(
            quadruplet_bound=quadruplet_bound, m_max=m_max
        ),
        "witness_6_fibonacci_mod210_shells": witness_fibonacci_mod210_shells(
            quadruplet_bound=quadruplet_bound
        ),
        "references": [
            "collatz_eabc_zeta_exponential_gedankenexperiment.md §9.8–9.9",
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
    parser.add_argument("--m-max", type=int, default=8, help="Fourier-Index m für C_m/Z_m")
    parser.add_argument(
        "--quadruplet-bound",
        type=int,
        default=500_000,
        help="Obergrenze p für Primvierlinge (Zeuge 5–6)",
    )
    parser.add_argument("--meromorphic-m-max", type=int, default=50)
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
        m_max=args.m_max,
        quadruplet_bound=args.quadruplet_bound,
        meromorphic_m_max=args.meromorphic_m_max,
    )

    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")

    w1 = result["witness_1_golden_resonance"]
    w2 = result["witness_2_prime_phases"]
    w3 = result["witness_3_fibonacci_windows"]
    w4 = result["witness_4_meromorphic_normal_form"]
    w5 = result["witness_5_golden_fourier"]
    w6 = result["witness_6_fibonacci_mod210_shells"]

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
    print(
        f"  Zeuge 4: meromorph max rel.err={w4['max_relative_error']:.2e} "
        f"(m_max={w4['m_max']})"
    )
    if w5["checkpoints"]:
        last = w5["checkpoints"][-1]
        z0 = last["Z_m"][0] if last["Z_m"] else float("nan")
        z1 = last["Z_m"][1] if len(last["Z_m"]) > 1 else float("nan")
        print(
            f"  Zeuge 5: N={last['N']} Q={last['Q']} Z_0={z0:.4f} Z_1={z1:.4f} "
            f"(m_max={w5['m_max']})"
        )
    stab = w6["fibonacci_stability"]
    lin = w6["linear_baseline"]
    print(
        f"  Zeuge 6: fib shells={w6['n_fibonacci_shells']} "
        f"sig (+,+,-) frac={stab.get('fraction_matches')} "
        f"linear frac={lin.get('fraction_matches')}"
    )
    print(f"  → {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
