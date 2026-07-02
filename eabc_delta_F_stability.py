#!/usr/bin/env python3
"""
EABC-Pfadsignatur: Δ_F(m)-Stabilitätstest
==========================================

Observable:
  Δ_F(m) = ||Σ_A^prime(m) - Σ_A^rand(m)||_F  /  ||Σ_A^rand(m)||_F

Wenn Δ_F(m) für wachsendes m nicht gegen 0 kollabiert, ist die
Kovarianzstruktur des Primstroms statistisch vom Nullmodell unterscheidbar —
ein belastbarer Level-2-Befund der EABC-Fluktuationsgeometrie.

Methoden:
  1. Σ_A^prime(m) aus disjunkten Fenstern des EABC-Primstroms (500 000 Prim.)
  2. Σ_A^rand(m)  aus B=50 unabhängigen Zufallswörtern gleicher Marginalen
  3. Unsicherheitsband: B=100 Zufallswörter, R=30× in Hälften splitten →
     Verteilung von Δ_F^null(m) unter H0
  4. Eigenwert-Spektrum: λ₁(Σ_A^prime) vs λ₁(Σ_A^rand) über alle m

Ausgabedateien:
  eabc_delta_F_stability.py   (dieses Skript)
  eabc_delta_F_results.csv
  eabc_delta_F_plot.png
"""

from __future__ import annotations

import csv
from pathlib import Path
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eabc_signature_analysis import sieve_primes, class_of

# ---------------------------------------------------------------------------
# Konstanten
# ---------------------------------------------------------------------------

WINDOW_SIZES = [500, 1000, 2000, 5000, 10_000, 20_000]
N_PRIMES     = 500_000
B_RAND       = 50       # Zufallswörter für Σ_A^rand
B_NULL_BAND  = 100      # Zufallswörter für Null-Unsicherheitsband
R_NULL_REPS  = 30       # Wiederholungen der Split-half-Prozedur
RANDOM_SEED  = 42

PAIR_CODES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
PAIR_NAMES = ["EA", "EB", "EC", "AB", "AC", "BC"]
LABEL_TO_CODE: dict[str, int] = {"E": 0, "A": 1, "B": 2, "C": 3}

EABC_COLOR  = "#1f77b4"
NULL_COLOR  = "#aaaaaa"
BAND_COLOR  = "#ff7f0e"


# ---------------------------------------------------------------------------
# Primzahlwort aufbauen
# ---------------------------------------------------------------------------

def build_eabc_word(n_primes: int) -> np.ndarray:
    """EABC-Codewort (int8) der ersten n_primes Primzahlen > 3."""
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
# Kern-Observable: normierter antisymmetrischer Vektor a ∈ ℝ⁶
# ---------------------------------------------------------------------------

def compute_nxy_nyx(word: np.ndarray, xc: int, yc: int) -> tuple[int, int]:
    """N_XY = #{(i,j): i<j, word[i]=X, word[j]=Y}  (O(m) vektorisiert)."""
    is_x = (word == xc).view(np.uint8)
    is_y = (word == yc).view(np.uint8)
    cumx = np.cumsum(is_x, dtype=np.int64)
    cumy = np.cumsum(is_y, dtype=np.int64)
    n_xy = int(np.dot(cumx - is_x, is_y))
    n_yx = int(np.dot(cumy - is_y, is_x))
    return n_xy, n_yx


def compute_a_vector(window: np.ndarray) -> np.ndarray:
    """Normierter 6-Vektor a ∈ Λ²(ℝ⁴). NaN wenn Paar gar nicht vorkommt."""
    a = np.empty(6, dtype=np.float64)
    for i, (xc, yc) in enumerate(PAIR_CODES):
        n_xy, n_yx = compute_nxy_nyx(window, xc, yc)
        total = n_xy + n_yx
        a[i] = (n_xy - n_yx) / total if total > 0 else np.nan
    return a


# ---------------------------------------------------------------------------
# Fensterdaten: EABC-Primstrom
# ---------------------------------------------------------------------------

def collect_prime_vectors(word: np.ndarray, m: int) -> np.ndarray:
    """
    a-Vektoren aus allen disjunkten Fenstern der Länge m.
    Rückgabe: (K, 6)-Array (NaN-Fenster werden verworfen).
    """
    n = len(word)
    starts = range(0, n - m + 1, m)
    K = len(starts)
    rows: list[np.ndarray] = []
    for k, start in enumerate(starts):
        a = compute_a_vector(word[start:start + m])
        if not np.any(np.isnan(a)):
            rows.append(a)
        if (k + 1) % 100 == 0 or k == K - 1:
            print(f"    EABC-Fenster: {k+1}/{K}", end="\r", flush=True)
    print()
    return np.array(rows, dtype=np.float64)


# ---------------------------------------------------------------------------
# Zufallswörter (Nullmodell)
# ---------------------------------------------------------------------------

def compute_random_a_vectors(
    m: int,
    B: int,
    marginals: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Erzeugt B unabhängige Zufallswörter der Länge m mit den gegebenen
    4-Klassen-Marginalverteilungen und berechnet jeweils den a-Vektor.
    Rückgabe: (B', 6)-Array (NaN-freie Zeilen).
    """
    rows: list[np.ndarray] = []
    for _ in range(B):
        word = rng.choice(4, size=m, p=marginals).astype(np.int8)
        a = compute_a_vector(word)
        if not np.any(np.isnan(a)):
            rows.append(a)
    return np.array(rows, dtype=np.float64)


# ---------------------------------------------------------------------------
# Kovarianzmatrix (robust: mindestens 7 Zeilen nötig)
# ---------------------------------------------------------------------------

def safe_cov(vecs: np.ndarray) -> np.ndarray | None:
    """np.cov(vecs.T), None wenn < 7 Zeilen."""
    if vecs.shape[0] < 7:
        return None
    return np.cov(vecs.T)


def frobenius_rel(A: np.ndarray, B: np.ndarray) -> float:
    """||A - B||_F / ||B||_F"""
    norm_B = float(np.linalg.norm(B, "fro"))
    if norm_B == 0.0:
        return np.nan
    return float(np.linalg.norm(A - B, "fro")) / norm_B


def top_eigenvalue(Sigma: np.ndarray) -> float:
    """Größter Eigenwert (reell, symmetrisch)."""
    evals = np.linalg.eigvalsh(Sigma)
    return float(evals.max())


# ---------------------------------------------------------------------------
# Null-Unsicherheitsband: Split-half über R Wiederholungen
# ---------------------------------------------------------------------------

def null_band(
    m: int,
    B_total: int,
    R_reps: int,
    marginals: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float, np.ndarray]:
    """
    Erzeugt R_reps × B_total/2 Zufalls-a-Vektoren-Paare und berechnet
    pro Paar Δ_F^null = ||Σ_rand_1 - Σ_rand_2||_F / ||Σ_rand_2||_F.

    Rückgabe: (mean_delta, std_delta, alle_delta_werte)
    """
    half = B_total // 2
    delta_vals: list[float] = []
    for _ in range(R_reps):
        vecs_all = compute_random_a_vectors(m, B_total, marginals, rng)
        if len(vecs_all) < 2 * 7:
            continue
        n_avail = len(vecs_all)
        split = n_avail // 2
        v1 = vecs_all[:split]
        v2 = vecs_all[split:2 * split]
        S1 = safe_cov(v1)
        S2 = safe_cov(v2)
        if S1 is None or S2 is None:
            continue
        delta_vals.append(frobenius_rel(S1, S2))
    arr = np.array(delta_vals, dtype=np.float64)
    return float(np.mean(arr)), float(np.std(arr, ddof=1)), arr


# ---------------------------------------------------------------------------
# Gesamtanalyse pro Fenstergröße
# ---------------------------------------------------------------------------

def analyse_window_size(
    word: np.ndarray,
    m: int,
    marginals: np.ndarray,
    rng: np.random.Generator,
) -> dict:
    """
    Führt die vollständige Δ_F-Analyse für eine Fenstergröße m durch.
    """
    t0 = time.perf_counter()

    # --- EABC-Primstrom ---
    print(f"  Sammle EABC-Fenstervektoren …")
    prime_vecs = collect_prime_vectors(word, m)
    K = len(prime_vecs)
    print(f"  → K = {K} Fenster")

    Sigma_prime = safe_cov(prime_vecs)
    if Sigma_prime is None:
        raise RuntimeError(f"Zu wenige EABC-Fenster für m={m}: K={K}")

    # --- Nullmodell: B=50 Zufallswörter ---
    print(f"  Erzeuge {B_RAND} Nullmodell-Wörter …", end=" ", flush=True)
    rand_vecs = compute_random_a_vectors(m, B_RAND, marginals, rng)
    print(f"→ {len(rand_vecs)} Vektoren")

    Sigma_rand = safe_cov(rand_vecs)
    if Sigma_rand is None:
        raise RuntimeError(f"Zu wenige Rand-Vektoren für m={m}: {len(rand_vecs)}")

    delta_F = frobenius_rel(Sigma_prime, Sigma_rand)

    # --- Eigenwerte ---
    lambda1_prime = top_eigenvalue(Sigma_prime)
    lambda1_rand  = top_eigenvalue(Sigma_rand)

    evals_prime = np.sort(np.linalg.eigvalsh(Sigma_prime))[::-1]
    evals_rand  = np.sort(np.linalg.eigvalsh(Sigma_rand))[::-1]

    # --- Null-Band (Split-half, R=30 Wiederholungen) ---
    print(
        f"  Berechne Null-Band "
        f"(B={B_NULL_BAND}, R={R_NULL_REPS} Splits) …",
        end=" ", flush=True,
    )
    null_mean, null_std, null_vals = null_band(
        m, B_NULL_BAND, R_NULL_REPS, marginals, rng
    )
    print(f"→ μ_null = {null_mean:.5f}, σ_null = {null_std:.5f}")

    signifikant = bool(delta_F > null_mean + 2.0 * null_std)

    elapsed = time.perf_counter() - t0
    print(
        f"  Δ_F = {delta_F:.5f}  |  Grenze = {null_mean + 2*null_std:.5f}"
        f"  |  {'SIGNIFIKANT ✓' if signifikant else 'n.s.'}  [{elapsed:.1f}s]"
    )

    return {
        "m"              : m,
        "K"              : K,
        "Sigma_prime"    : Sigma_prime,
        "Sigma_rand"     : Sigma_rand,
        "evals_prime"    : evals_prime,
        "evals_rand"     : evals_rand,
        "delta_F"        : delta_F,
        "null_mean"      : null_mean,
        "null_std"       : null_std,
        "null_vals"      : null_vals,
        "signifikant"    : signifikant,
        "lambda1_prime"  : lambda1_prime,
        "lambda1_rand"   : lambda1_rand,
    }


# ---------------------------------------------------------------------------
# Ausgabetabelle (Terminal)
# ---------------------------------------------------------------------------

def print_table(results: list[dict]) -> None:
    header = (
        f"{'m':>7} | {'K':>5} | {'Δ_F':>8} | "
        f"{'Δ_null_μ':>10} | {'Δ_null_σ':>10} | "
        f"{'sig?':^6} | {'λ₁_prime':>10} | {'λ₁_rand':>10}"
    )
    sep = "─" * len(header)
    print()
    print(sep)
    print("ΔELTA_F(m) STABILITÄTSTEST — Ergebnisübersicht")
    print(sep)
    print(header)
    print(sep)
    for r in results:
        sig_str = "YES" if r["signifikant"] else "no"
        print(
            f"{r['m']:>7} | {r['K']:>5} | {r['delta_F']:>8.5f} | "
            f"{r['null_mean']:>10.5f} | {r['null_std']:>10.5f} | "
            f"{sig_str:^6} | {r['lambda1_prime']:>10.6f} | {r['lambda1_rand']:>10.6f}"
        )
    print(sep)
    print()
    any_sig = any(r["signifikant"] for r in results)
    if any_sig:
        print(
            "BEFUND: Δ_F(m) überschreitet das 2σ-Null-Band für mindestens eine"
            " Fenstergröße.\n"
            "  → Kovarianzstruktur des Primstroms ist statistisch vom Nullmodell"
            " unterscheidbar.\n"
            "  → Erster belastbarer Level-2-Befund der EABC-Fluktuationsgeometrie."
        )
    else:
        print(
            "KEIN BEFUND: Δ_F(m) bleibt innerhalb des Null-Bandes für alle m.\n"
            "  → Keine statistisch belastbare Abweichung der Kovarianzstruktur."
        )
    print()


# ---------------------------------------------------------------------------
# CSV-Export
# ---------------------------------------------------------------------------

def export_csv(results: list[dict], out_path: Path) -> None:
    fields = [
        "m", "K", "delta_F",
        "delta_F_null_mean", "delta_F_null_std",
        "signifikant",
        "lambda1_prime", "lambda1_rand",
    ]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in results:
            writer.writerow({
                "m"                  : r["m"],
                "K"                  : r["K"],
                "delta_F"            : round(r["delta_F"], 8),
                "delta_F_null_mean"  : round(r["null_mean"], 8),
                "delta_F_null_std"   : round(r["null_std"], 8),
                "signifikant"        : int(r["signifikant"]),
                "lambda1_prime"      : round(r["lambda1_prime"], 10),
                "lambda1_rand"       : round(r["lambda1_rand"], 10),
            })
    print(f"CSV gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Visualisierung
# ---------------------------------------------------------------------------

def plot_results(results: list[dict], out_path: Path) -> None:
    """
    2 Subplots:
      1. Δ_F(m) als Funktion von m (log-Skala x-Achse), mit Null-Band
      2. λ₁(Σ_A^prime) vs λ₁(Σ_A^rand) über m
    """
    ms            = np.array([r["m"]             for r in results], dtype=float)
    delta_Fs      = np.array([r["delta_F"]        for r in results])
    null_means    = np.array([r["null_mean"]       for r in results])
    null_stds     = np.array([r["null_std"]        for r in results])
    lambda1_prime = np.array([r["lambda1_prime"]   for r in results])
    lambda1_rand  = np.array([r["lambda1_rand"]    for r in results])
    signifs       = [r["signifikant"] for r in results]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        r"$\Delta_F(m)$-Stabilitätstest: EABC-Primstrom vs. Nullmodell"
        "\n"
        r"$\Delta_F(m) = \|\Sigma_A^\mathrm{prime}(m) - "
        r"\Sigma_A^\mathrm{rand}(m)\|_F \;/\; \|\Sigma_A^\mathrm{rand}(m)\|_F$",
        fontsize=12,
        fontweight="bold",
    )

    # ── Subplot 1: Δ_F(m) ──────────────────────────────────────────────────
    ax = axes[0]
    band_lo = null_means - 2.0 * null_stds
    band_hi = null_means + 2.0 * null_stds
    band_lo = np.clip(band_lo, 0, None)

    ax.fill_between(
        ms, band_lo, band_hi,
        color=NULL_COLOR, alpha=0.35, label=r"Null-Band $\mu \pm 2\sigma$",
    )
    ax.plot(ms, null_means, color=NULL_COLOR, linewidth=1.4,
            linestyle="--", label=r"$\mu(\Delta_F^\mathrm{null})$")
    ax.plot(ms, delta_Fs, color=EABC_COLOR, linewidth=2.0,
            marker="o", markersize=6, label=r"$\Delta_F^\mathrm{prime}(m)$", zorder=3)

    # Signifikante Punkte markieren
    for i, (m_i, df_i, sig) in enumerate(zip(ms, delta_Fs, signifs)):
        if sig:
            ax.scatter(m_i, df_i, color="red", s=80, zorder=5,
                       marker="*", linewidths=0)

    ax.set_xscale("log")
    ax.set_xlabel("Fenstergröße $m$", fontsize=11)
    ax.set_ylabel(r"$\Delta_F(m)$", fontsize=11)
    ax.set_title(r"$\Delta_F(m)$ als Funktion der Fenstergröße", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.30, which="both")
    ax.set_xticks(ms)
    ax.set_xticklabels([str(int(v)) for v in ms], rotation=30, ha="right", fontsize=8)

    # Annotation: Trend
    if len(results) >= 3:
        log_ms = np.log10(ms)
        coeffs = np.polyfit(log_ms, delta_Fs, deg=1)
        trend_label = f"Trend (slope={coeffs[0]:+.4f})"
        trend_line  = np.polyval(coeffs, log_ms)
        ax.plot(ms, trend_line, color="darkblue", linewidth=1.0,
                linestyle=":", label=trend_label, zorder=2)
        ax.legend(fontsize=8.5)

    # ── Subplot 2: λ₁(m) ──────────────────────────────────────────────────
    ax2 = axes[1]
    ax2.plot(ms, lambda1_prime, color=EABC_COLOR, linewidth=2.0,
             marker="o", markersize=6, label=r"$\lambda_1(\Sigma_A^\mathrm{prime})$")
    ax2.plot(ms, lambda1_rand, color=NULL_COLOR, linewidth=2.0,
             marker="s", markersize=6, label=r"$\lambda_1(\Sigma_A^\mathrm{rand})$")

    # Verhältnis als zweite y-Achse
    ax2r = ax2.twinx()
    ratio = lambda1_prime / np.where(lambda1_rand > 0, lambda1_rand, np.nan)
    ax2r.plot(ms, ratio, color=BAND_COLOR, linewidth=1.4,
              linestyle="--", marker="^", markersize=5,
              label=r"$\lambda_1^\mathrm{prime}/\lambda_1^\mathrm{rand}$")
    ax2r.set_ylabel(r"Verhältnis $\lambda_1^\mathrm{prime}/\lambda_1^\mathrm{rand}$",
                    fontsize=9, color=BAND_COLOR)
    ax2r.tick_params(axis="y", labelcolor=BAND_COLOR, labelsize=8)
    ax2r.axhline(1.0, color=BAND_COLOR, linewidth=0.8, linestyle=":", alpha=0.6)

    ax2.set_xscale("log")
    ax2.set_xlabel("Fenstergröße $m$", fontsize=11)
    ax2.set_ylabel(r"Größter Eigenwert $\lambda_1$", fontsize=11)
    ax2.set_title(
        r"Größter Eigenwert $\lambda_1$ von $\Sigma_A$", fontsize=10
    )
    # Gemeinsame Legende
    lines_a, labs_a = ax2.get_legend_handles_labels()
    lines_b, labs_b = ax2r.get_legend_handles_labels()
    ax2.legend(lines_a + lines_b, labs_a + labs_b, fontsize=8.5)
    ax2.grid(True, alpha=0.30, which="both")
    ax2.set_xticks(ms)
    ax2.set_xticklabels([str(int(v)) for v in ms], rotation=30, ha="right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Hauptprogramm
# ---------------------------------------------------------------------------

def main() -> None:
    out_dir = Path(__file__).resolve().parent
    rng = np.random.default_rng(RANDOM_SEED)

    print("=" * 72)
    print("EABC-Pfadsignatur: Δ_F(m)-Stabilitätstest")
    print(f"  N_PRIMES    = {N_PRIMES:,}")
    print(f"  Fenstergrößen: {WINDOW_SIZES}")
    print(f"  B_RAND      = {B_RAND}  (Zufallswörter für Σ_rand)")
    print(f"  B_NULL_BAND = {B_NULL_BAND}, R_NULL_REPS = {R_NULL_REPS}")
    print("=" * 72)

    # Primzahlwort
    print(f"\nGeneriere EABC-Codewort der ersten {N_PRIMES:,} Primzahlen > 3 …")
    word = build_eabc_word(N_PRIMES)
    print(f"  Wortlänge: {len(word):,}")

    from collections import Counter
    counts = Counter(word.tolist())
    marginals = np.array(
        [counts[c] / len(word) for c in range(4)],
        dtype=np.float64,
    )
    for code, label in enumerate(("E", "A", "B", "C")):
        n = counts[code]
        print(f"  N_{label} = {n:>7}  ({100 * n / len(word):.3f}%)"
              f"  [Marginal: {marginals[code]:.5f}]")

    # Analyse pro Fenstergröße
    results: list[dict] = []

    for m in WINDOW_SIZES:
        print(f"\n{'─'*60}")
        print(f"  Fenstergröße m = {m:,}")
        print(f"{'─'*60}")
        res = analyse_window_size(word, m, marginals, rng)
        results.append(res)

    # Ausgabe
    print_table(results)

    # CSV
    csv_path = out_dir / "eabc_delta_F_results.csv"
    export_csv(results, csv_path)

    # Plot
    print("\nErstelle Δ_F-Plot …")
    plot_path = out_dir / "eabc_delta_F_plot.png"
    plot_results(results, plot_path)

    print("\nFertig.")


if __name__ == "__main__":
    main()
