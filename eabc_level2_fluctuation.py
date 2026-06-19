#!/usr/bin/env python3
"""
EABC Level-2-Fluktuationsgeometrie auf Λ²(ℝ⁴)
=============================================

Observable: a = (a_EA, a_EB, a_EC, a_AB, a_AC, a_BC) ∈ ℝ⁶
  a_XY = (N_XY - N_YX) / (N_XY + N_YX) ∈ [-1, 1]

Kernmetrik (Pflicht-Checkpoints):
  Δ_F(m) = ||Σ_A^prime(m) - Σ_A^null(m)||_F / ||Σ_A^null(m)||_F

wobei Σ_A = E[(a - μ_A)(a - μ_A)^T] über disjunkte Fenster der Größe m
im EABC-Primstrom bzw. unter Nullmodell-Ensemble.

Nullmodell-Hierarchie (Stufe 1–3):
  Stufe 1 — perm:      volle Permutation (marginaltreu, Reihenfolge zerstört)
  Stufe 2 — markov:    Markov-erhaltend (lokale Übergangswahrscheinlichkeiten)
  Stufe 3 — hl:        Hardy-Littlewood-konsistent (Stub, noch nicht implementiert)

Ausgabe: JSON mit delta_F_perm, delta_F_markov (+ delta_F_hl=null) + stdout
"""

from __future__ import annotations

import argparse
import json
from math import isqrt
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Konstanten
# ---------------------------------------------------------------------------

CHECKPOINTS = [1_000, 2_000, 5_000, 10_000, 20_000]
DEFAULT_N_PRIMES = 500_000
DEFAULT_B_RAND = 50
DEFAULT_SEED = 42

PAIR_CODES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
PAIR_NAMES = ["EA", "EB", "EC", "AB", "AC", "BC"]
LABEL_TO_CODE: dict[str, int] = {"E": 0, "A": 1, "B": 2, "C": 3}
RESIDUE_MAP = {1: "E", 5: "A", 7: "B", 11: "C"}


def sieve_primes(limit: int) -> list[int]:
    """Sieb des Eratosthenes bis `limit` (inklusiv)."""
    if limit < 2:
        return []
    flags = bytearray([1]) * (limit + 1)
    flags[0] = flags[1] = 0
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p * p :: p] = bytearray(len(flags[p * p :: p]))
    return [i for i, ok in enumerate(flags) if ok]


def class_of(p: int) -> str | None:
    return RESIDUE_MAP.get(p % 12)


# ---------------------------------------------------------------------------
# EABC-Wort
# ---------------------------------------------------------------------------

def build_eabc_word(n_primes: int) -> np.ndarray:
    """Erstellt das EABC-Codewort (int8) der ersten n_primes Primzahlen > 3."""
    limit = max(10_000_000, n_primes * 25)
    while True:
        primes = sieve_primes(limit)
        codes: list[int] = []
        for p in primes:
            if p <= 3:
                continue
            c = class_of(p)
            if c is not None:
                codes.append(LABEL_TO_CODE[c])
            if len(codes) >= n_primes:
                break
        if len(codes) >= n_primes:
            return np.array(codes[:n_primes], dtype=np.int8)
        limit = int(limit * 1.5)


# ---------------------------------------------------------------------------
# Level-2-Vektor a ∈ Λ²(ℝ⁴)
# ---------------------------------------------------------------------------

def compute_nxy_nyx(word: np.ndarray, xc: int, yc: int) -> tuple[int, int]:
    is_x = (word == xc).view(np.uint8)
    is_y = (word == yc).view(np.uint8)
    cumx = np.cumsum(is_x, dtype=np.int64)
    cumy = np.cumsum(is_y, dtype=np.int64)
    n_xy = int(np.dot(cumx - is_x, is_y))
    n_yx = int(np.dot(cumy - is_y, is_x))
    return n_xy, n_yx


def compute_a_vector(window: np.ndarray) -> np.ndarray:
    """Normierter antisymmetrischer 6-Vektor für ein Fenster."""
    a = np.empty(6, dtype=np.float64)
    for i, (xc, yc) in enumerate(PAIR_CODES):
        n_xy, n_yx = compute_nxy_nyx(window, xc, yc)
        total = n_xy + n_yx
        a[i] = (n_xy - n_yx) / total if total > 0 else np.nan
    return a


def collect_window_vectors(word: np.ndarray, m: int) -> np.ndarray:
    """Sammelt a-Vektoren für alle disjunkten Fenster der Größe m."""
    n = len(word)
    rows: list[np.ndarray] = []
    for start in range(0, n - m + 1, m):
        a = compute_a_vector(word[start:start + m])
        if not np.any(np.isnan(a)):
            rows.append(a)
    return np.array(rows, dtype=np.float64)


def collect_perm_null_vectors(
    word: np.ndarray,
    m: int,
    B: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Stufe 1 — Permutations-Nullmodell: B Zufallspermutationen pro Fenster."""
    n = len(word)
    rows: list[np.ndarray] = []
    buf = np.empty(m, dtype=np.int8)
    for start in range(0, n - m + 1, m):
        win = word[start:start + m].copy()
        for _ in range(B):
            buf[:] = win
            rng.shuffle(buf)
            a = compute_a_vector(buf)
            if not np.any(np.isnan(a)):
                rows.append(a)
    return np.array(rows, dtype=np.float64)


def _markov_transition_matrix(window: np.ndarray) -> np.ndarray:
    """Zeilen-stochastische 4×4-Übergangsmatrix aus Adjazenzpaaren."""
    counts = np.ones((4, 4), dtype=np.float64)  # Laplace-Glättung
    for i in range(len(window) - 1):
        counts[int(window[i]), int(window[i + 1])] += 1.0
    row_sums = counts.sum(axis=1, keepdims=True)
    return counts / row_sums


def _markov_resample(
    window: np.ndarray,
    trans: np.ndarray,
    buf: np.ndarray,
    rng: np.random.Generator,
) -> None:
    """Erzeugt Markov-Nullfolge gleicher Länge mit erhaltener Startklasse."""
    buf[0] = window[0]
    for i in range(1, len(buf)):
        buf[i] = int(rng.choice(4, p=trans[int(buf[i - 1])]))


def collect_markov_null_vectors(
    word: np.ndarray,
    m: int,
    B: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Stufe 2 — Markov-erhaltendes Nullmodell (lokale Übergangswahrscheinlichkeiten)."""
    n = len(word)
    rows: list[np.ndarray] = []
    buf = np.empty(m, dtype=np.int8)
    for start in range(0, n - m + 1, m):
        win = word[start:start + m]
        trans = _markov_transition_matrix(win)
        for _ in range(B):
            _markov_resample(win, trans, buf, rng)
            a = compute_a_vector(buf)
            if not np.any(np.isnan(a)):
                rows.append(a)
    return np.array(rows, dtype=np.float64)


def collect_hl_null_vectors(
    word: np.ndarray,
    m: int,
    B: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Stufe 3 — Hardy-Littlewood-konsistentes Nullmodell (PLACEHOLDER).

    Geplant: Erzeugung von Nullfolgen mit HL-konsistenten Paar-/Mehrfachkorrelationen
    (Cramér-ähnlicher Primprozess + mod-12-Kanalrestriktion). Siehe §4.8.2.
    """
    raise NotImplementedError(
        "HL-Nullmodell (Stufe 3): Hardy-Littlewood-konsistente Ensemble-Generierung "
        "noch nicht implementiert — siehe collatz_eabc_zirkulationshypothese.md §4.8.2"
    )


# Rückwärtskompatibilität
collect_null_vectors = collect_perm_null_vectors


# ---------------------------------------------------------------------------
# Kovarianz und Δ_F
# ---------------------------------------------------------------------------

def empirical_covariance(vectors: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """μ_A und Σ_A aus Fenstervektoren."""
    mu = vectors.mean(axis=0)
    sigma = np.cov(vectors.T)
    return mu, sigma


def delta_F(Sigma_prime: np.ndarray, Sigma_rand: np.ndarray) -> float:
    """Relative Frobenius-Abweichung der Kovarianzmatrizen."""
    diff = Sigma_prime - Sigma_rand
    norm_rand = float(np.linalg.norm(Sigma_rand, "fro"))
    if norm_rand == 0.0:
        return float("nan")
    return float(np.linalg.norm(diff, "fro") / norm_rand)


def spectrum_summary(Sigma: np.ndarray) -> list[float]:
    evals = np.linalg.eigvalsh(Sigma)
    return [float(v) for v in sorted(evals, reverse=True)]


def analyze_checkpoint(
    word: np.ndarray,
    m: int,
    B_rand: int,
    rng: np.random.Generator,
) -> dict:
    """Berechnet Σ_A^prime und Δ_F(m) gegen perm- und Markov-Nullmodell."""
    eabc_vecs = collect_window_vectors(word, m)
    perm_vecs = collect_perm_null_vectors(word, m, B_rand, rng)
    markov_vecs = collect_markov_null_vectors(word, m, B_rand, rng)

    if eabc_vecs.shape[0] < 2:
        raise ValueError(f"Zu wenige Fenster für m={m}: K_prime={eabc_vecs.shape[0]}")
    for label, vecs in (("perm", perm_vecs), ("markov", markov_vecs)):
        if vecs.shape[0] < 2:
            raise ValueError(
                f"Zu wenige Null-Vektoren ({label}) für m={m}: K={vecs.shape[0]}"
            )

    mu_prime, sigma_prime = empirical_covariance(eabc_vecs)
    _, sigma_perm = empirical_covariance(perm_vecs)
    _, sigma_markov = empirical_covariance(markov_vecs)
    dF_perm = delta_F(sigma_prime, sigma_perm)
    dF_markov = delta_F(sigma_prime, sigma_markov)

    return {
        "m": m,
        "K_prime": int(eabc_vecs.shape[0]),
        "K_perm": int(perm_vecs.shape[0]),
        "K_markov": int(markov_vecs.shape[0]),
        "mu_A_prime": [float(x) for x in mu_prime],
        "mu_A_prime_norm": float(np.linalg.norm(mu_prime)),
        "Sigma_A_prime": sigma_prime.tolist(),
        "Sigma_A_perm": sigma_perm.tolist(),
        "Sigma_A_markov": sigma_markov.tolist(),
        "spec_prime": spectrum_summary(sigma_prime),
        "spec_perm": spectrum_summary(sigma_perm),
        "spec_markov": spectrum_summary(sigma_markov),
        "Delta_F": dF_perm,
        "delta_F_perm": dF_perm,
        "delta_F_markov": dF_markov,
        "delta_F_hl": None,
    }


def run_fluctuation_test(
    n_primes: int = DEFAULT_N_PRIMES,
    checkpoints: list[int] | None = None,
    B_rand: int = DEFAULT_B_RAND,
    seed: int = DEFAULT_SEED,
) -> dict:
    """Führt Δ_F(m)-Test an allen Checkpoints aus."""
    checkpoints = list(checkpoints or CHECKPOINTS)
    rng = np.random.default_rng(seed)
    word = build_eabc_word(n_primes)

    results: list[dict] = []
    for m in checkpoints:
        if m > len(word):
            continue
        results.append(analyze_checkpoint(word, m, B_rand, rng))

    return {
        "n_primes": n_primes,
        "word_length": int(len(word)),
        "B_rand": B_rand,
        "seed": seed,
        "checkpoints": checkpoints,
        "null_models": {
            "perm": "Stufe 1 — volle Permutation (marginaltreu)",
            "markov": "Stufe 2 — Markov-erhaltend (lokale Übergänge)",
            "hl": "Stufe 3 — HL-konsistent (Stub, delta_F_hl=null)",
        },
        "results": results,
    }


def print_summary(report: dict) -> None:
    print()
    print("=" * 78)
    print("LEVEL-2-FLUKTUATIONSGEOMETRIE: Δ_F(m) auf Λ²(ℝ⁴) — Multi-Nullmodell")
    print("=" * 78)
    print(f"  N_PRIMES = {report['n_primes']:,}  |  B_RAND = {report['B_rand']}")
    print("  Nullmodelle: Stufe 1 perm | Stufe 2 markov | Stufe 3 hl (Stub)")
    print()
    print(
        f"  {'m':>8}  {'K':>5}  {'|μ_A|':>10}  "
        f"{'Δ_F^perm':>10}  {'Δ_F^markov':>10}  Befund"
    )
    print("  " + "-" * 68)
    for row in report["results"]:
        mu_norm = row["mu_A_prime_norm"]
        dF_perm = row["delta_F_perm"]
        dF_markov = row["delta_F_markov"]
        ref = max(dF_perm, dF_markov)
        if ref > 0.10:
            verdict = "deutlich ≠ Null"
        elif ref > 0.03:
            verdict = "moderat ≠ Null"
        else:
            verdict = "≈ Null"
        print(
            f"  {row['m']:>8,}  {row['K_prime']:>5}  "
            f"{mu_norm:>10.6f}  {dF_perm:>10.6f}  {dF_markov:>10.6f}  {verdict}"
        )
    print()
    print("  Schlüssel: Gegner ist falsches Nullmodell — perm vs. markov vs. HL (§4.8.2)")
    print("  Interpretation: Δ_F(m) ↛ 0  ⇒  robuste Level-2-Fluktuationsstruktur")
    print("=" * 78)


def export_json(report: dict, out_path: Path) -> None:
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    print(f"JSON gespeichert: {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="EABC Level-2-Fluktuationsgeometrie: Δ_F(m)-Test"
    )
    parser.add_argument(
        "--n-primes", type=int, default=DEFAULT_N_PRIMES,
        help=f"Anzahl Primzahlen > 3 (default: {DEFAULT_N_PRIMES:,})",
    )
    parser.add_argument(
        "--checkpoints", type=int, nargs="+", default=CHECKPOINTS,
        help="Fenstergrößen m für Δ_F(m)",
    )
    parser.add_argument(
        "--B-rand", type=int, default=DEFAULT_B_RAND,
        help="Permutationen pro Fenster im Nullmodell",
    )
    parser.add_argument(
        "--seed", type=int, default=DEFAULT_SEED,
        help="Zufallsseed für Permutationsnull",
    )
    parser.add_argument(
        "--json", type=Path, default=None,
        help="JSON-Ausgabedatei (default: eabc_level2_fluctuation.json)",
    )
    args = parser.parse_args()

    print("Generiere EABC-Codewort …")
    report = run_fluctuation_test(
        n_primes=args.n_primes,
        checkpoints=args.checkpoints,
        B_rand=args.B_rand,
        seed=args.seed,
    )
    print_summary(report)

    out = args.json or Path(__file__).resolve().parent / "eabc_level2_fluctuation.json"
    export_json(report, out)


if __name__ == "__main__":
    main()
