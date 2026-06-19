#!/usr/bin/env python3
"""
EABC Kompressions-Robustheitstest (Test C, §4.8.8–§4.8.9)
=========================================================

Programm-Kette:
    Primdaten → Kompression C_k → Observablen M → Resonanztest

Observablenvektor:
    M = (D_E, Q_E, R_N, C_F)

Strukturstabilitätsindex:
    S(C) = 1 - ||V_C - V_raw|| / ||V_raw||
    S = 1: perfekte Erhaltung; S = 0: vollständiger Verlust (relativ zur Referenz).

Referenz V_raw: M auf dem kanonischen EABC-Wort (C_E auf dem Primstrom).

Kompressionsfamilie:
    C_E   — mod 12, EABC-Labels {E,A,B,C} ≡ {1,5,7,11}
    C_8   — mod 8,  Einheiten {1,3,5,7} → vier Klassen
    C_10  — mod 10, Einheiten {1,3,7,9} → vier Klassen
    C_30  — mod 30, acht Einheiten → vier Klassen (Rang-Quartile)
    C_210 — mod 210, 48 Einheiten → vier Klassen (Rang-Quartile)

Tests (Priorität):
    C — C_E vs C_8,C_10,C_30,C_210  (dieses Skript)
    B — Σ_A(N)-Stabilität           (eabc_level2_fluctuation.py)
    A — A_n Riemann-Detektor        (eabc_zeta_f_riemann_zeros.py)

Referenz: collatz_eabc_zirkulationshypothese.md §4.8.8–§4.8.9
"""

from __future__ import annotations

import argparse
import json
import math
from math import gcd, isqrt
from pathlib import Path
from typing import Any

import numpy as np

from eabc_level2_fluctuation import (
    LABEL_TO_CODE,
    RESIDUE_MAP,
    build_eabc_word_with_primes,
    class_of,
    collect_window_vectors,
    empirical_covariance,
)
from eabc_sigma_fibonacci_coupling import build_k_F_vector, coupling_C_F
from eabc_zeta_f_riemann_zeros import load_zeros

LOG_PHI = math.log((1.0 + math.sqrt(5.0)) / 2.0)
DEFAULT_OUTPUT = Path("eabc_compression_robustness.json")
DEFAULT_ZEROS = Path("zeros6.npz")
DEFAULT_N_PRIMES = 1_000_000
DEFAULT_WINDOW_M = 5_000
DEFAULT_R_N_CAP = 50_000
DEFAULT_SEED = 42

COMPRESSION_MODULI: dict[str, int] = {
    "C_E": 12,
    "C_8": 8,
    "C_10": 10,
    "C_30": 30,
    "C_210": 210,
}


def sieve_primes(limit: int) -> list[int]:
    if limit < 2:
        return []
    flags = bytearray([1]) * (limit + 1)
    flags[0] = flags[1] = 0
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p * p :: p] = bytearray(len(flags[p * p :: p]))
    return [i for i, ok in enumerate(flags) if ok]


def coprime_units(modulus: int) -> list[int]:
    return sorted(k for k in range(1, modulus) if gcd(k, modulus) == 1)


def build_residue_partition(modulus: int) -> dict[int, int]:
    """
    Deterministische 4-Klassen-Partition der Einheiten modulo `modulus`.

    C_E: feste EABC-Zuordnung {1→E, 5→A, 7→B, 11→C}.
    Sonst: sortierte Einheitenliste in vier Rang-Quartile.
    """
    if modulus == 12:
        return {r: LABEL_TO_CODE[RESIDUE_MAP[r]] for r in RESIDUE_MAP}

    units = coprime_units(modulus)
    n_units = len(units)
    mapping: dict[int, int] = {}
    for idx, residue in enumerate(units):
        cls = min(3, (4 * idx) // n_units)
        mapping[residue] = cls
    return mapping


def compress_primes(
    primes: np.ndarray,
    modulus: int,
    residue_map: dict[int, int],
) -> np.ndarray:
    """Wendet C_k auf Primarray an; nur Primzahlen mit zulässiger Restklasse."""
    codes: list[int] = []
    for p in primes:
        r = int(p % modulus)
        if r in residue_map:
            codes.append(residue_map[r])
    return np.array(codes, dtype=np.int8)


def counts_from_word(word: np.ndarray) -> dict[str, int]:
    inv = {0: "E", 1: "A", 2: "B", 3: "C"}
    counts = {label: 0 for label in "EABC"}
    for code in word:
        counts[inv[int(code)]] += 1
    return counts


def compute_D_E_Q_E(word: np.ndarray) -> tuple[float, float]:
    """D_E = N_A + N_C - N_B - N_E; Q_E = N_E + N_A + N_B + N_C."""
    c = counts_from_word(word)
    d_e = float(c["A"] + c["C"] - c["B"] - c["E"])
    q_e = float(len(word))
    return d_e, q_e


CODE_TO_RESIDUE = np.array([1, 5, 7, 11], dtype=np.float64)


def compute_R_N_abs(
    word: np.ndarray,
    gamma: np.ndarray,
    *,
    cap: int = DEFAULT_R_N_CAP,
) -> float:
    """
    |R_N| mit EABC-Phasen aus dem (komprimierten) Wort:
    Θ(n) = (π/2)·σ_EABC(n), σ_EABC ∈ {1,5,7,11} je Klasse; Φ(γ) = γ log φ.
    """
    n = min(len(word), len(gamma), cap)
    if n == 0:
        return 0.0
    sigma = CODE_TO_RESIDUE[word[:n].astype(np.intp)]
    g = gamma[:n].astype(np.float64)
    theta = 0.5 * math.pi * sigma
    phi = g * LOG_PHI
    delta = theta - phi
    r_n = np.mean(np.cos(delta) + 1j * np.sin(delta))
    return float(abs(r_n))


def compute_C_F(word: np.ndarray, m: int) -> float:
    """C_F aus Σ_A am Fenster-Checkpoint m."""
    if len(word) < 2 * m:
        return float("nan")
    vecs = collect_window_vectors(word, m)
    if vecs.shape[0] < 2:
        return float("nan")
    _, sigma = empirical_covariance(vecs)
    k_f = build_k_F_vector(dim=6)
    return float(coupling_C_F(sigma, k_f)["C_F"])


def feature_vector_for_stability(
    word: np.ndarray,
    gamma: np.ndarray,
    m: int,
    *,
    r_n_cap: int,
) -> tuple[np.ndarray, dict[str, float]]:
    """
  Skaleninvariante Features für S(C):
    (W_E, |R_N|, C_F) mit W_E = D_E/Q_E.
  Rohobservablen M=(D_E,Q_E,R_N,C_F) separat in Metadaten.
    """
    d_e, q_e = compute_D_E_Q_E(word)
    w_e = d_e / q_e if q_e else 0.0
    r_n = compute_R_N_abs(word, gamma, cap=r_n_cap)
    c_f = compute_C_F(word, m)
    if not np.isfinite(c_f):
        c_f = 0.0
    raw = {
        "D_E": d_e,
        "Q_E": q_e,
        "W_E": w_e,
        "R_N": r_n,
        "C_F": c_f,
    }
    return np.array([w_e, r_n, c_f], dtype=np.float64), raw


def observable_vector(
    word: np.ndarray,
    gamma: np.ndarray,
    m: int,
    *,
    r_n_cap: int,
) -> np.ndarray:
    """M = (D_E, Q_E, R_N, C_F) als reeller 4-Vektor (Rohwerte)."""
    _, raw = feature_vector_for_stability(word, gamma, m, r_n_cap=r_n_cap)
    return np.array([raw["D_E"], raw["Q_E"], raw["R_N"], raw["C_F"]], dtype=np.float64)


def stability_index(v_c: np.ndarray, v_raw: np.ndarray) -> float:
    """
    S(C) = 1 - ||V_C - V_raw|| / ||V_raw|| auf skaleninvarianten Features
    (W_E, |R_N|, C_F); Q_E als Normierung in W_E = D_E/Q_E enthalten.
    """
    norm_raw = float(np.linalg.norm(v_raw))
    if norm_raw < 1e-18 or not np.isfinite(norm_raw):
        return float("nan")
    diff = float(np.linalg.norm(v_c - v_raw))
    s = 1.0 - diff / norm_raw
    return max(0.0, min(1.0, s))


def classify_fall(s_c_e: float, s_k: float) -> str:
    """Grobe Fall-Lesart (1)–(3) aus §4.8.8."""
    if not np.isfinite(s_c_e) or not np.isfinite(s_k):
        return "unbestimmt"
    if s_k >= 0.95:
        return "1_signal_roh"
    if s_c_e > s_k + 0.05:
        return "2_nur_C_E"
    return "3_artifact"


def run_test_c(
    n_primes: int = DEFAULT_N_PRIMES,
    *,
    window_m: int = DEFAULT_WINDOW_M,
    zeros_path: Path = DEFAULT_ZEROS,
    r_n_cap: int = DEFAULT_R_N_CAP,
    compressions: list[str] | None = None,
) -> dict[str, Any]:
    """Test C: Strukturstabilität S(C) für alle Kompressionen."""
    if compressions is None:
        compressions = list(COMPRESSION_MODULI.keys())

    word_raw, primes = build_eabc_word_with_primes(n_primes)
    gamma = load_zeros(zeros_path)

    # Fenster m anpassen, falls n_primes nicht teilbar
    m = window_m
    while m > 100 and len(word_raw) % m != 0 and len(word_raw) // m < 2:
        m //= 2

    v_raw_feat, raw_meta = feature_vector_for_stability(word_raw, gamma, m, r_n_cap=r_n_cap)
    v_raw = observable_vector(word_raw, gamma, m, r_n_cap=r_n_cap)
    d_raw, q_raw = raw_meta["D_E"], raw_meta["Q_E"]

    results: dict[str, Any] = {}
    for name in compressions:
        modulus = COMPRESSION_MODULI[name]
        residue_map = build_residue_partition(modulus)
        word_c = compress_primes(primes, modulus, residue_map)
        v_c_feat, meta = feature_vector_for_stability(word_c, gamma, m, r_n_cap=r_n_cap)
        v_c = observable_vector(word_c, gamma, m, r_n_cap=r_n_cap)
        s_c = stability_index(v_c_feat, v_raw_feat)
        d_c, q_c = meta["D_E"], meta["Q_E"]
        results[name] = {
            "modulus": modulus,
            "word_length": int(len(word_c)),
            "D_E": d_c,
            "Q_E": q_c,
            "W_E": d_c / q_c if q_c else None,
            "R_N_abs": float(v_c[2]),
            "C_F": float(v_c[3]) if np.isfinite(v_c[3]) else None,
            "M": {
                "D_E": float(v_c[0]),
                "Q_E": float(v_c[1]),
                "R_N": float(v_c[2]),
                "C_F": float(v_c[3]) if np.isfinite(v_c[3]) else None,
            },
            "S": s_c,
            "delta_M_l2": float(np.linalg.norm(v_c - v_raw)),
            "relative_delta_M": float(np.linalg.norm(v_c - v_raw) / np.linalg.norm(v_raw))
            if np.linalg.norm(v_raw) > 0
            else None,
        }

    s_ce = results.get("C_E", {}).get("S")
    interpretation = "offen"
    if s_ce is not None and np.isfinite(s_ce):
        alt = [results[k]["S"] for k in results if k != "C_E" and np.isfinite(results[k]["S"])]
        if alt:
            mean_alt = float(np.mean(alt))
            if s_ce > mean_alt + 0.05:
                interpretation = "S(C_E) >> S(C_k) — EABC-spezifische Strukturerhaltung"
            elif abs(s_ce - mean_alt) < 0.05:
                interpretation = "S(C_E) ≈ S(C_k) — keine EABC-Spezifität"

    return {
        "label": "Test C — Kompressionsrobustheit (§4.8.8–§4.8.9)",
        "test": "C",
        "priority": 1,
        "n_primes_requested": n_primes,
        "n_primes_eabc": int(len(word_raw)),
        "window_m": m,
        "r_n_cap": r_n_cap,
        "zeros_file": str(zeros_path),
        "reference": "C_E",
        "V_raw": {
            "D_E": float(v_raw[0]),
            "Q_E": float(v_raw[1]),
            "R_N": float(v_raw[2]),
            "C_F": float(v_raw[3]) if np.isfinite(v_raw[3]) else None,
            "W_E": d_raw / q_raw if q_raw else None,
        },
        "S_definition": (
            "S(C) = 1 - ||V_C - V_raw|| / ||V_raw||; "
            "V = (W_E, |R_N|, C_F) mit W_E = D_E/Q_E (skaleninvariant zu M=(D_E,Q_E,R_N,C_F))"
        ),
        "compressions": results,
        "interpretation": interpretation,
        "fall_classification": {
            k: classify_fall(float(s_ce) if s_ce is not None else float("nan"), v["S"])
            for k, v in results.items()
            if k != "C_E"
        },
    }


def print_report(result: dict[str, Any]) -> None:
    print("=" * 72)
    print(result["label"])
    print("=" * 72)
    print(f"  n_primes (EABC) = {result['n_primes_eabc']:,}")
    print(f"  Fenster m         = {result['window_m']}")
    print(f"  R_N cap           = {result['r_n_cap']:,}")
    v = result["V_raw"]
    print(f"\n  V_raw: D_E={v['D_E']:.1f}  Q_E={v['Q_E']:.0f}  "
          f"R_N={v['R_N']:.6f}  C_F={v['C_F']}")
    print(f"\n  {'C_k':<8} {'S(C)':>8} {'D_E':>10} {'W_E':>12} {'|R_N|':>10} {'C_F':>12}")
    print("  " + "-" * 64)
    for name, row in result["compressions"].items():
        w_e = row["W_E"]
        w_str = f"{w_e:+.6f}" if w_e is not None else "n/a"
        c_f = row["C_F"]
        c_str = f"{c_f:.6f}" if c_f is not None else "n/a"
        print(
            f"  {name:<8} {row['S']:8.6f} {row['D_E']:10.1f} "
            f"{w_str:>12} {row['R_N_abs']:10.6f} {c_str:>12}"
        )
    print(f"\n  Interpretation: {result['interpretation']}")
    print("=" * 72)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="EABC Kompressions-Robustheitstest (Test C, §4.8.8)"
    )
    parser.add_argument(
        "--n-primes",
        type=int,
        default=DEFAULT_N_PRIMES,
        help=f"Anzahl Primzahlen im EABC-Strom (Default: {DEFAULT_N_PRIMES:,})",
    )
    parser.add_argument(
        "--window-m",
        type=int,
        default=DEFAULT_WINDOW_M,
        help=f"Fenstergröße m für Σ_A / C_F (Default: {DEFAULT_WINDOW_M})",
    )
    parser.add_argument(
        "--r-n-cap",
        type=int,
        default=DEFAULT_R_N_CAP,
        help=f"Max. Indizes für R_N (Default: {DEFAULT_R_N_CAP:,})",
    )
    parser.add_argument(
        "--zeros-file",
        type=Path,
        default=DEFAULT_ZEROS,
        help="Riemann-Nullstellen (zeros6.npz)",
    )
    parser.add_argument(
        "--compressions",
        nargs="+",
        default=list(COMPRESSION_MODULI.keys()),
        choices=list(COMPRESSION_MODULI.keys()),
        help="Kompressionsoperatoren",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="JSON-Ausgabedatei",
    )
    args = parser.parse_args()

    result = run_test_c(
        n_primes=args.n_primes,
        window_m=args.window_m,
        zeros_path=args.zeros_file,
        r_n_cap=args.r_n_cap,
        compressions=args.compressions,
    )
    print_report(result)

    out = args.output
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nJSON → {out}")


if __name__ == "__main__":
    main()
