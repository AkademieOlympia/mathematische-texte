#!/usr/bin/env python3
"""
Level-2 ↔ Stufe-0: optionale Fibonacci-Kamm-Kopplung
=====================================================

Testet spektrale Kommensurabilität zwischen Σ_A^prime und dem φ-Resonatorkamm
(ζ_F-Stufe 0) — **kein Erklärungsanspruch**, nur diagnostische Projektion.

K_F(s) = Σ_{m≥0} w_m · (-1)^m φ^{-(s+2m)} / (1 - (-1)^m φ^{-(s+2m)})

k_F ∈ ℝ^6: Re(K_F) an goldenen Resonanzpunkten t_k = 2πk/log φ, k=1..6.

C_F(Σ, k_F) = Σ_i λ_i |⟨v_i, k_F⟩|²   (λ_i, v_i = Spektrum von Σ)

Nullvergleich: Σ_A^prime vs. golden / perm / Markov (wie Level-2).

Referenz: collatz_eabc_zirkulationshypothese.md §4.8.3
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from eabc_level2_fluctuation import (
    CHECKPOINTS,
    DEFAULT_B_RAND,
    DEFAULT_GOLDEN_BINS,
    DEFAULT_SEED,
    PAIR_NAMES,
    build_eabc_word_with_primes,
    collect_golden_null_vectors,
    collect_markov_null_vectors,
    collect_perm_null_vectors,
    collect_window_vectors,
    empirical_covariance,
)

PHI = (1.0 + math.sqrt(5.0)) / 2.0
LOG_PHI = math.log(PHI)
DEFAULT_M_MAX = 12
DEFAULT_OUTPUT = Path("eabc_sigma_fibonacci_coupling.json")


def resonator_term(s: complex, m: int) -> complex:
    """Einzelmodus (-1)^m φ^{-(s+2m)} / (1 - (-1)^m φ^{-(s+2m)})."""
    parity = (-1) ** m
    pow_term = PHI ** (-(s + 2 * m))
    denom = 1.0 - parity * pow_term
    if abs(denom) < 1e-14:
        return complex(float("inf"), float("inf"))
    return parity * pow_term / denom


def K_F(s: complex, m_max: int, weights: np.ndarray) -> complex:
    """Fibonacci-φ-Kamm (minimale Form ohne Pochhammer/5^{s/2})."""
    total = 0.0 + 0.0j
    n = min(m_max + 1, len(weights))
    for m in range(n):
        term = resonator_term(s, m)
        if math.isinf(term.real):
            continue
        total += weights[m] * term
    return total


def golden_resonance_t(k: int) -> float:
    """t_k = 2πk/log φ (imaginärer Resonanzpunkt auf Re(s)=1/2)."""
    return 2.0 * math.pi * k / LOG_PHI


def build_k_F_vector(
    dim: int = 6,
    m_max: int = DEFAULT_M_MAX,
    weights: np.ndarray | None = None,
    resonance_indices: tuple[int, ...] | None = None,
) -> np.ndarray:
    """
    k_F[j] = Re(K_F(1/2 + i t_{k_j})), normiert.

    Default: k_j an goldenen Resonanzpunkten k_j = j+1.
    """
    if weights is None:
        weights = np.ones(m_max + 1, dtype=np.float64)
    if resonance_indices is None:
        resonance_indices = tuple(range(1, dim + 1))

    k = np.empty(dim, dtype=np.float64)
    for j, rk in enumerate(resonance_indices[:dim]):
        s = 0.5 + 1j * golden_resonance_t(rk)
        k[j] = K_F(s, m_max, weights).real

    norm = float(np.linalg.norm(k))
    if norm > 0.0:
        k /= norm
    return k


def coupling_C_F(sigma: np.ndarray, k_F: np.ndarray) -> dict:
    """
    C_F = Σ_i λ_i |⟨v_i, k_F⟩|² mit λ_i ≥ 0 (numerisch geklemmt).
    """
    evals, evecs = np.linalg.eigh(sigma)
    k = k_F / np.linalg.norm(k_F) if np.linalg.norm(k_F) > 0 else k_F
    overlaps = np.array([abs(np.dot(evecs[:, i], k)) ** 2 for i in range(len(evals))])
    evals_pos = np.maximum(evals, 0.0)
    c_f = float(np.dot(evals_pos, overlaps))
    return {
        "C_F": c_f,
        "eigenvalues": [float(v) for v in evals],
        "overlaps": [float(v) for v in overlaps],
        "weighted_overlaps": [float(evals_pos[i] * overlaps[i]) for i in range(len(evals))],
    }


def random_k_F_baseline(
    sigma: np.ndarray,
    n_draws: int,
    rng: np.random.Generator,
    dim: int = 6,
) -> dict:
    """Mittelwert von C_F über zufällige Einheitsvektoren (Isotropie-Null)."""
    c_vals: list[float] = []
    for _ in range(n_draws):
        v = rng.normal(size=dim)
        v /= np.linalg.norm(v)
        c_vals.append(coupling_C_F(sigma, v)["C_F"])
    arr = np.array(c_vals)
    return {
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
        "n_draws": n_draws,
    }


def analyze_coupling_checkpoint(
    word: np.ndarray,
    primes: np.ndarray,
    m: int,
    B_rand: int,
    rng: np.random.Generator,
    m_max: int,
    weights: np.ndarray,
    golden_bins: int,
    random_baseline_draws: int = 20,
) -> dict:
    """C_F für Σ_A^prime und Nullmodelle am Checkpoint m."""
    eabc_vecs = collect_window_vectors(word, m)
    golden_vecs = collect_golden_null_vectors(
        word, primes, m, B_rand, rng, n_bins=golden_bins
    )
    perm_vecs = collect_perm_null_vectors(word, m, B_rand, rng)
    markov_vecs = collect_markov_null_vectors(word, m, B_rand, rng)

    _, sigma_prime = empirical_covariance(eabc_vecs)
    _, sigma_golden = empirical_covariance(golden_vecs)
    _, sigma_perm = empirical_covariance(perm_vecs)
    _, sigma_markov = empirical_covariance(markov_vecs)

    k_F = build_k_F_vector(dim=6, m_max=m_max, weights=weights)
    k_F_shuffled = k_F.copy()
    rng.shuffle(k_F_shuffled)
    k_F_shuffled /= np.linalg.norm(k_F_shuffled)

    c_prime = coupling_C_F(sigma_prime, k_F)
    c_golden = coupling_C_F(sigma_golden, k_F)
    c_perm = coupling_C_F(sigma_perm, k_F)
    c_markov = coupling_C_F(sigma_markov, k_F)
    c_shuffled = coupling_C_F(sigma_prime, k_F_shuffled)
    baseline = random_k_F_baseline(sigma_prime, random_baseline_draws, rng)

    def ratio(a: float, b: float) -> float | None:
        if b == 0.0 or math.isnan(b):
            return None
        return float(a / b)

    return {
        "m": m,
        "K_prime": int(eabc_vecs.shape[0]),
        "k_F": [float(x) for x in k_F],
        "k_F_labels": [f"t_{j+1}" for j in range(6)],
        "pair_names": PAIR_NAMES,
        "C_F_prime": c_prime["C_F"],
        "C_F_golden": c_golden["C_F"],
        "C_F_perm": c_perm["C_F"],
        "C_F_markov": c_markov["C_F"],
        "C_F_shuffled_k": c_shuffled["C_F"],
        "C_F_random_baseline": baseline,
        "ratio_prime_over_golden": ratio(c_prime["C_F"], c_golden["C_F"]),
        "ratio_prime_over_perm": ratio(c_prime["C_F"], c_perm["C_F"]),
        "ratio_prime_over_markov": ratio(c_prime["C_F"], c_markov["C_F"]),
        "ratio_prime_over_random": ratio(c_prime["C_F"], baseline["mean"]),
        "overlaps_prime": c_prime["overlaps"],
        "spec_prime": c_prime["eigenvalues"],
    }


def run_coupling_test(
    n_primes: int,
    checkpoints: list[int] | None = None,
    B_rand: int = DEFAULT_B_RAND,
    seed: int = DEFAULT_SEED,
    m_max: int = DEFAULT_M_MAX,
    weights: np.ndarray | None = None,
    golden_bins: int = DEFAULT_GOLDEN_BINS,
    random_baseline_draws: int = 20,
) -> dict:
    checkpoints = list(checkpoints or CHECKPOINTS)
    if weights is None:
        weights = np.ones(m_max + 1, dtype=np.float64)
    rng = np.random.default_rng(seed)
    word, primes = build_eabc_word_with_primes(n_primes)

    results: list[dict] = []
    for m in checkpoints:
        if m > len(word):
            continue
        results.append(
            analyze_coupling_checkpoint(
                word,
                primes,
                m,
                B_rand,
                rng,
                m_max,
                weights,
                golden_bins,
                random_baseline_draws,
            )
        )

    return {
        "epistemic_label": "Experiment (Schicht C) — optionale Stufe-0-Kopplung, kein Theorem",
        "n_primes": n_primes,
        "word_length": int(len(word)),
        "B_rand": B_rand,
        "seed": seed,
        "m_max": m_max,
        "weights": [float(w) for w in weights[: m_max + 1]],
        "golden_bins": golden_bins,
        "formula": {
            "K_F": "Σ_m w_m (-1)^m φ^{-(s+2m)} / (1-(-1)^m φ^{-(s+2m)})",
            "C_F": "Σ_i λ_i |⟨v_i, k_F⟩|²",
            "k_F": "Re(K_F(1/2 + i·2πk/log φ)), k=1..6, normiert",
        },
        "null_models": ["golden", "perm", "markov", "shuffled_k", "random_unit"],
        "checkpoints": checkpoints,
        "results": results,
    }


def print_summary(report: dict) -> None:
    print()
    print("=" * 88)
    print("LEVEL-2 ↔ STUFE-0: Fibonacci-Kamm-Kopplung C_F (optional, kein Theorem)")
    print("=" * 88)
    print(
        f"  N_PRIMES = {report['n_primes']:,}  |  m_max = {report['m_max']}  |  "
        f"B_RAND = {report['B_rand']}  |  seed = {report['seed']}"
    )
    print()
    print(
        f"  {'m':>8}  {'C_F^prime':>12}  {'C_F^gold':>12}  {'C_F^perm':>12}  "
        f"{'C_F^mark':>12}  {'prim/gold':>9}  {'prim/rand':>9}"
    )
    print("  " + "-" * 84)
    for row in report["results"]:
        r_pg = row.get("ratio_prime_over_golden")
        r_pr = row.get("ratio_prime_over_random")
        print(
            f"  {row['m']:>8,}  {row['C_F_prime']:>12.6e}  "
            f"{row['C_F_golden']:>12.6e}  {row['C_F_perm']:>12.6e}  "
            f"{row['C_F_markov']:>12.6e}  "
            f"{r_pg if r_pg is not None else float('nan'):>9.3f}  "
            f"{r_pr if r_pr is not None else float('nan'):>9.3f}"
        )
    print()
    print("  Lesart: C_F misst gewichtete Überlappung Spec(Σ_A) mit φ-Kamm — diagnostisch")
    print("=" * 88)


def export_json(report: dict, path: Path) -> None:
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"JSON gespeichert: {path}")


def export_csv(report: dict, path: Path) -> None:
    fields = [
        "m",
        "C_F_prime",
        "C_F_golden",
        "C_F_perm",
        "C_F_markov",
        "C_F_shuffled_k",
        "ratio_prime_over_golden",
        "ratio_prime_over_perm",
        "ratio_prime_over_markov",
        "ratio_prime_over_random",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in report["results"]:
            writer.writerow({k: row.get(k) for k in fields})
    print(f"CSV gespeichert: {path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Level-2 Σ_A ↔ Stufe-0 ζ_F: optionale Fibonacci-Kamm-Kopplung C_F"
    )
    parser.add_argument(
        "--limit",
        "--n-primes",
        dest="n_primes",
        type=int,
        default=100_000,
        help="Anzahl Primzahlen > 3 (default: 100000)",
    )
    parser.add_argument(
        "--checkpoints",
        type=int,
        nargs="+",
        default=CHECKPOINTS,
        help="Fenstergrößen m",
    )
    parser.add_argument("--B-rand", type=int, default=DEFAULT_B_RAND)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--m-max",
        type=int,
        default=DEFAULT_M_MAX,
        help="Maximaler Resonator-Index m in K_F",
    )
    parser.add_argument(
        "--weight",
        type=float,
        nargs="+",
        default=None,
        help="Gewichte w_m (default: alle 1); Länge ≤ m_max+1",
    )
    parser.add_argument("--golden-bins", type=int, default=DEFAULT_GOLDEN_BINS)
    parser.add_argument(
        "--random-draws",
        type=int,
        default=20,
        help="Zufallsvektor-Null für C_F-Baseline",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Ausgabedatei (.json oder .csv)",
    )
    parser.add_argument(
        "--format",
        choices=("json", "csv"),
        default=None,
        help="Ausgabeformat (default: aus Dateiendung)",
    )
    args = parser.parse_args()

    if args.weight is not None:
        weights = np.array(args.weight, dtype=np.float64)
        if len(weights) < args.m_max + 1:
            weights = np.pad(weights, (0, args.m_max + 1 - len(weights)), constant_values=1.0)
    else:
        weights = np.ones(args.m_max + 1, dtype=np.float64)

    print("Generiere EABC-Codewort und berechne C_F …")
    report = run_coupling_test(
        n_primes=args.n_primes,
        checkpoints=args.checkpoints,
        B_rand=args.B_rand,
        seed=args.seed,
        m_max=args.m_max,
        weights=weights,
        golden_bins=args.golden_bins,
        random_baseline_draws=args.random_draws,
    )
    print_summary(report)

    fmt = args.format or ("csv" if args.output.suffix.lower() == ".csv" else "json")
    if fmt == "csv":
        export_csv(report, args.output)
    else:
        export_json(report, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
