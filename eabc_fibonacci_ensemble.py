#!/usr/bin/env python3
"""
EABC-Fibonacci-Ensemble: gewichtete Klassensumme unter φ-Kontraktion (§4.8.24)
==============================================================================

Komplementäre Perspektive zum δ(p)-Defektprogramm (§4.8.23): statt lokalem
Schalenabstand δ(p) aggregiert ein Ensemble über EABC-Klassen mit Fibonacci-
Kontraktionsgewicht κ(θ)=φ^{-θ}, θ=θ_φ(p)=(log p)/(log φ) mod 1.

Kerngröße (Charakter χ ∈ {χ_A, χ_B, χ_C}, default χ_A):
  Ψ_F = Σ_c w_c · Σ_{p: v(p)=c} κ(θ_φ(p)) · χ(v(p))

Vergleich:
  - Einzelklasse Ψ_c
  - Gesamtensemble Ψ_F
  - reine Fibonacci-Kontraktion Ψ_Fib = Σ_p κ(θ_φ(p))  (ohne Klassenasymmetrie)
  - Permutations-Null (Labels vertauscht, Marginalen erhalten)

Referenz: collatz_eabc_zirkulationshypothese.md §4.8.3–4.8.7, §4.8.24
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from eabc_level2_fluctuation import (
    DEFAULT_SEED,
    LABEL_TO_CODE,
    build_eabc_word_with_primes,
    theta_phi_from_primes,
)
from eabc_shell_prime_distribution import LABELS, V4_CHARACTERS

PHI = (1.0 + math.sqrt(5.0)) / 2.0
DEFAULT_OUTPUT = Path("eabc_fibonacci_ensemble.json")
DEFAULT_CLASS_WEIGHTS = {"E": 1.0, "A": 1.0, "B": 1.0, "C": 1.0}
CODE_TO_LABEL = {v: k for k, v in LABEL_TO_CODE.items()}


def fibonacci_contraction(theta: np.ndarray) -> np.ndarray:
    """κ(θ)=φ^{-θ} für θ∈[0,1) — Kontraktionsfaktor am goldenen Log-Gitter."""
    return PHI ** (-theta.astype(np.float64))


def character_value(label: str, character: str = "chi_A") -> int:
    return int(V4_CHARACTERS[label][character])


def class_counts_in_bins(
    word: np.ndarray,
    theta: np.ndarray,
    n_bins: int,
) -> np.ndarray:
    """N_c(k): 4×n_bins Zählmatrix (Klasse × φ-Bin)."""
    bins = np.minimum((theta * n_bins).astype(np.int64), n_bins - 1)
    counts = np.zeros((4, n_bins), dtype=np.int64)
    for code in range(4):
        mask = word == code
        if not np.any(mask):
            continue
        b = bins[mask]
        for k in range(n_bins):
            counts[code, k] = int(np.sum(b == k))
    return counts


def weighted_class_sums(
    word: np.ndarray,
    theta: np.ndarray,
    class_weights: dict[str, float],
    character: str,
) -> dict[str, float]:
    """Ψ_c und Ψ_F = Σ_c w_c Ψ_c mit κ(θ)-Gewicht und V₄-Charakter."""
    kappa = fibonacci_contraction(theta)
    per_class: dict[str, float] = {lab: 0.0 for lab in LABELS}
    for i, code in enumerate(word):
        lab = CODE_TO_LABEL[int(code)]
        per_class[lab] += float(class_weights[lab] * kappa[i] * character_value(lab, character))
    ensemble = float(sum(per_class.values()))
    return {"per_class": per_class, "ensemble": ensemble}


def pure_fibonacci_mass(theta: np.ndarray) -> float:
    """Ψ_Fib = Σ_p κ(θ_φ(p)) — reine Kontraktion ohne Klassenstruktur."""
    return float(np.sum(fibonacci_contraction(theta)))


def perm_null_ensemble(
    word: np.ndarray,
    theta: np.ndarray,
    class_weights: dict[str, float],
    character: str,
    n_draws: int,
    rng: np.random.Generator,
) -> dict:
    """Permutations-Null: EABC-Labels zufällig, κ(θ) und Marginalen fix."""
    values: list[float] = []
    shuffled = word.copy()
    for _ in range(n_draws):
        rng.shuffle(shuffled)
        values.append(
            weighted_class_sums(shuffled, theta, class_weights, character)["ensemble"]
        )
    arr = np.array(values, dtype=np.float64)
    return {
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
        "n_draws": n_draws,
        "percentile_2_5": float(np.percentile(arr, 2.5)),
        "percentile_97_5": float(np.percentile(arr, 97.5)),
        "values": [float(v) for v in values],
    }


def standardized_Z(observed: float, null_stats: dict) -> float | None:
    std = null_stats.get("std")
    if std is None or std == 0.0 or math.isnan(std):
        return None
    return float((observed - null_stats["mean"]) / std)


def bin_ensemble_asymmetry(
    counts: np.ndarray,
    class_weights: dict[str, float],
    character: str,
    n_bins: int,
) -> dict:
    """
    Diskrete Bin-Variante: Σ_k τ(k) Σ_c w_c N_c(k) σ_c
    mit τ(k)=φ^{-k/n_bins} (Fibonacci-Fenster).
    """
    tau = np.array([PHI ** (-k / n_bins) for k in range(n_bins)], dtype=np.float64)
    per_class_bins: dict[str, list[float]] = {lab: [] for lab in LABELS}
    ensemble_bins = np.zeros(n_bins, dtype=np.float64)
    for code, lab in enumerate(LABELS):
        sigma = character_value(lab, character)
        row = counts[code].astype(np.float64) * class_weights[lab] * sigma
        per_class_bins[lab] = [float(tau[k] * row[k]) for k in range(n_bins)]
        ensemble_bins += tau * row
    pure_fib = float(np.sum(tau * counts.sum(axis=0)))
    ensemble_total = float(np.sum(ensemble_bins))
    return {
        "tau": [float(x) for x in tau],
        "per_class_weighted": per_class_bins,
        "ensemble_by_bin": [float(x) for x in ensemble_bins],
        "ensemble_total": ensemble_total,
        "pure_fibonacci_bins": pure_fib,
        "asymmetry_ratio": ensemble_total / pure_fib if pure_fib != 0.0 else None,
    }


def run_ensemble_test(
    n_primes: int,
    n_bins: int = 8,
    class_weights: dict[str, float] | None = None,
    character: str = "chi_A",
    B_perm: int = 50,
    seed: int = DEFAULT_SEED,
) -> dict:
    weights = dict(DEFAULT_CLASS_WEIGHTS if class_weights is None else class_weights)
    rng = np.random.default_rng(seed)
    word, primes = build_eabc_word_with_primes(n_primes)
    theta = theta_phi_from_primes(primes)

    sums = weighted_class_sums(word, theta, weights, character)
    psi_fib = pure_fibonacci_mass(theta)
    counts = class_counts_in_bins(word, theta, n_bins)
    bin_stats = bin_ensemble_asymmetry(counts, weights, character, n_bins)
    perm_null = perm_null_ensemble(word, theta, weights, character, B_perm, rng)

    psi_ens = sums["ensemble"]
    asym_vs_fib = psi_ens / psi_fib if psi_fib != 0.0 else None
    z_perm = standardized_Z(psi_ens, perm_null)

    return {
        "epistemic_label": "Experiment (Schicht C) — Ensemble unter φ-Kontraktion, kein Theorem",
        "section": "4.8.24",
        "n_primes": n_primes,
        "n_bins": n_bins,
        "seed": seed,
        "B_perm": B_perm,
        "character": character,
        "class_weights": weights,
        "formula": {
            "kappa": "κ(θ)=φ^{-θ}, θ=(log p)/(log φ) mod 1",
            "Psi_F": "Σ_c w_c · Σ_{p: v(p)=c} κ(θ_φ(p)) · χ(v(p))",
            "Psi_Fib": "Σ_p κ(θ_φ(p))",
            "tau_bin": "τ(k)=φ^{-k/B} für Bin-Ensemble Σ_k τ(k) Σ_c w_c N_c(k) χ(v)",
        },
        "Psi_F": psi_ens,
        "Psi_Fib": psi_fib,
        "per_class": sums["per_class"],
        "asymmetry_vs_fibonacci": asym_vs_fib,
        "bin_ensemble": bin_stats,
        "perm_null": perm_null,
        "Z_perm": z_perm,
        "ratio_prime_over_perm_null": (
            psi_ens / perm_null["mean"] if perm_null["mean"] != 0.0 else None
        ),
        "class_counts_bins": {
            lab: [int(x) for x in counts[i]] for i, lab in enumerate(LABELS)
        },
    }


def print_summary(report: dict) -> None:
    print()
    print("=" * 72)
    print("EABC-FIBONACCI-ENSEMBLE (§4.8.24) — kein Theorem")
    print("=" * 72)
    print(
        f"  N_PRIMES = {report['n_primes']:,}  |  Bins = {report['n_bins']}  |  "
        f"χ = {report['character']}  |  seed = {report['seed']}"
    )
    print()
    print(f"  Ψ_F (Ensemble)     = {report['Psi_F']:.6f}")
    print(f"  Ψ_Fib (rein φ)     = {report['Psi_Fib']:.6f}")
    r_fib = report.get("asymmetry_vs_fibonacci")
    print(f"  Ψ_F / Ψ_Fib        = {r_fib if r_fib is not None else float('nan'):.6f}")
    print()
    print("  Einzelklassen Ψ_c:")
    for lab in LABELS:
        print(f"    {lab}: {report['per_class'][lab]:.6f}")
    print()
    perm = report["perm_null"]
    z = report.get("Z_perm")
    print(
        f"  Perm-Null: E[Ψ_F]={perm['mean']:.6f}, σ={perm['std']:.6f}, "
        f"Z={z if z is not None else float('nan'):.2f}"
    )
    r_pp = report.get("ratio_prime_over_perm_null")
    print(f"  Ψ_F / E[Perm]      = {r_pp if r_pp is not None else float('nan'):.6f}")
    be = report["bin_ensemble"]
    print(
        f"  Bin-Ensemble Σ_k τ N_c χ: total={be['ensemble_total']:.2f}, "
        f"rein φ={be['pure_fibonacci_bins']:.2f}, "
        f"ratio={be.get('asymmetry_ratio') or float('nan'):.6f}"
    )
    print()
    print("  Lesart: Komplement zu δ(p) — aggregierte Klassensumme unter κ(θ);")
    print("  „Atmen des Zahlenraums“ = benannte Analogie (Kontraktion/Expansion), kein Beweis.")
    print("=" * 72)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="EABC-Fibonacci-Ensemble: gewichtete Klassensumme unter φ-Kontraktion (§4.8.24)"
    )
    parser.add_argument(
        "--limit",
        "--n-primes",
        dest="n_primes",
        type=int,
        default=100_000,
        help="Anzahl Primzahlen > 3 (default: 100000)",
    )
    parser.add_argument("--bins", type=int, default=8, help="φ-Phasen-Bins (default: 8)")
    parser.add_argument(
        "--character",
        choices=("chi_A", "chi_B", "chi_C"),
        default="chi_A",
        help="V₄-Charakter für σ_c (default: chi_A)",
    )
    parser.add_argument("--B-perm", type=int, default=50, help="Permutations-Null-Züge")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--weight",
        nargs=4,
        type=float,
        metavar=("w_E", "w_A", "w_B", "w_C"),
        default=None,
        help="Klassengewichte w_c (default: alle 1)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="JSON-Ausgabe",
    )
    args = parser.parse_args()

    weights = None
    if args.weight is not None:
        weights = dict(zip(LABELS, args.weight))

    print("Berechne EABC-Fibonacci-Ensemble …")
    report = run_ensemble_test(
        n_primes=args.n_primes,
        n_bins=args.bins,
        class_weights=weights,
        character=args.character,
        B_perm=args.B_perm,
        seed=args.seed,
    )
    print_summary(report)
    args.output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"JSON gespeichert: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
