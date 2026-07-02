#!/usr/bin/env python3
"""
EABC-Pfadsignatur: Kovarianzstruktur der antisymmetrischen Level-2-Signatur
===========================================================================

Observable: a = (a_EA, a_EB, a_EC, a_AB, a_AC, a_BC) ∈ ℝ⁶
  a_XY = (N_XY - N_YX) / (N_XY + N_YX) ∈ [-1, 1]
  lebt in Λ²(ℝ⁴) mit 6 Freiheitsgraden.

Kernfrage:
  Unterscheidet sich die Geometrie der Fluktuationen (Kovarianzmatrix Σ_A)
  im EABC-Primstrom vom Permutations-Nullmodell — auch wenn die Mittelwerte
  gegen Null gehen?

Methoden:
  1. Fensterdaten: normierter 6-Vektor a pro disjunktem Fenster
  2. Σ_A = empirische 6×6-Kovarianzmatrix (EABC)
  3. Σ_A^rand = Kovarianzmatrix unter Permutationsnullmodell
  4. Eigenwertzerlegung beider Matrizen
  5. Vergleich: Frobenius-Norm ||Σ_A - Σ_A^rand||_F / ||Σ_A^rand||_F
  6. Scatter der Fenstervektoren in der PC1-PC2-Ebene

Ausgabedateien:
  eabc_altsig2_covariance.py
  eabc_covariance_results.csv
  eabc_covariance_plot.png
"""

from __future__ import annotations

import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eabc_signature_analysis import sieve_primes, class_of

# ---------------------------------------------------------------------------
# Konstanten
# ---------------------------------------------------------------------------

WINDOW_SIZES = [1000, 5000]
N_PRIMES = 500_000
B_RAND = 50           # Permutationen pro Fenster für Nullmodell
RANDOM_SEED = 42

# Alle 6 Paare in Λ²(ℝ⁴): E=0, A=1, B=2, C=3
PAIR_CODES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
PAIR_NAMES = ["EA", "EB", "EC", "AB", "AC", "BC"]
PAIR_LABELS = [
    r"$a_{EA}$", r"$a_{EB}$", r"$a_{EC}$",
    r"$a_{AB}$", r"$a_{AC}$", r"$a_{BC}$",
]
LABEL_TO_CODE: dict[str, int] = {"E": 0, "A": 1, "B": 2, "C": 3}

EABC_COLOR = "#1f77b4"
NULL_COLOR = "#aaaaaa"
SCATTER_ALPHA = 0.45


# ---------------------------------------------------------------------------
# Primzahlwort erzeugen
# ---------------------------------------------------------------------------

def build_eabc_word(n_primes: int) -> np.ndarray:
    """
    Erstellt das EABC-Codewort (int8) der ersten n_primes Primzahlen > 3.
    Codes: E=0, A=1, B=2, C=3.
    """
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
    """
    N_XY = #{(i,j): i<j, word[i]=X, word[j]=Y}  (O(m) vectorized)
    """
    is_x = (word == xc).view(np.uint8)
    is_y = (word == yc).view(np.uint8)
    cumx = np.cumsum(is_x, dtype=np.int64)
    cumy = np.cumsum(is_y, dtype=np.int64)
    n_xy = int(np.dot(cumx - is_x, is_y))
    n_yx = int(np.dot(cumy - is_y, is_x))
    return n_xy, n_yx


def compute_a_vector(window: np.ndarray) -> np.ndarray:
    """
    Berechnet den normierten 6-Vektor a ∈ Λ²(ℝ⁴) für ein Fenster.
    Komponente i = (N_XY - N_YX) / (N_XY + N_YX), NaN wenn Nenner = 0.
    """
    a = np.empty(6, dtype=np.float64)
    for i, (xc, yc) in enumerate(PAIR_CODES):
        n_xy, n_yx = compute_nxy_nyx(window, xc, yc)
        total = n_xy + n_yx
        a[i] = (n_xy - n_yx) / total if total > 0 else np.nan
    return a


# ---------------------------------------------------------------------------
# Fensterdaten sammeln (EABC und Nullmodell)
# ---------------------------------------------------------------------------

def collect_window_vectors(word: np.ndarray, m: int) -> np.ndarray:
    """
    Sammelt den normierten a-Vektor für alle disjunkten Fenster.
    Rückgabe: (K, 6)-Array, nur vollständige (kein NaN) Zeilen.
    """
    n = len(word)
    starts = range(0, n - m + 1, m)
    rows: list[np.ndarray] = []
    K = len(starts)
    for k, start in enumerate(starts):
        a = compute_a_vector(word[start:start + m])
        if not np.any(np.isnan(a)):
            rows.append(a)
        if (k + 1) % 50 == 0 or k == K - 1:
            print(f"    EABC-Fenster: {k+1}/{K}", end="\r", flush=True)
    print()
    return np.array(rows, dtype=np.float64)


def collect_null_vectors(
    word: np.ndarray,
    m: int,
    B: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Nullmodell: Für jedes disjunkte Fenster werden B Zufallspermutationen
    des Fensters erzeugt. Rückgabe: (K*B, 6)-Array der Permutations-a-Vektoren.
    """
    n = len(word)
    starts = list(range(0, n - m + 1, m))
    K = len(starts)
    rows: list[np.ndarray] = []
    buf = np.empty(m, dtype=np.int8)
    for k, start in enumerate(starts):
        win = word[start:start + m].copy()
        for _ in range(B):
            buf[:] = win
            rng.shuffle(buf)
            a = compute_a_vector(buf)
            if not np.any(np.isnan(a)):
                rows.append(a)
        if (k + 1) % 20 == 0 or k == K - 1:
            print(f"    Null-Permutationen: Fenster {k+1}/{K}", end="\r", flush=True)
    print()
    return np.array(rows, dtype=np.float64)


# ---------------------------------------------------------------------------
# Kovarianzstruktur-Analyse
# ---------------------------------------------------------------------------

def covariance_analysis(
    eabc_vectors: np.ndarray,
    null_vectors: np.ndarray,
) -> dict:
    """
    Berechnet Kovarianzmatrix, Korrelationsmatrix, Eigenwertzerlegung und
    Frobenius-Abstand zum Nullmodell.

    Rückgabe: dict mit allen relevanten Größen.
    """
    mu_A = eabc_vectors.mean(axis=0)
    Sigma_A = np.cov(eabc_vectors.T)

    mu_rand = null_vectors.mean(axis=0)
    Sigma_rand = np.cov(null_vectors.T)

    # Korrelationsmatrizen
    diag_A = np.sqrt(np.diag(Sigma_A))
    diag_rand = np.sqrt(np.diag(Sigma_rand))

    with np.errstate(invalid="ignore"):
        Corr_A = Sigma_A / np.outer(diag_A, diag_A)
        Corr_rand = Sigma_rand / np.outer(diag_rand, diag_rand)

    # Eigenwertzerlegung (absteigende Reihenfolge)
    evals_A, evecs_A = np.linalg.eigh(Sigma_A)
    order_A = np.argsort(evals_A)[::-1]
    evals_A = evals_A[order_A]
    evecs_A = evecs_A[:, order_A]

    evals_rand, evecs_rand = np.linalg.eigh(Sigma_rand)
    order_rand = np.argsort(evals_rand)[::-1]
    evals_rand = evals_rand[order_rand]
    evecs_rand = evecs_rand[:, order_rand]

    # Frobenius-Abstand (normiert)
    diff = Sigma_A - Sigma_rand
    frob_abs = float(np.linalg.norm(diff, "fro"))
    frob_rel = frob_abs / float(np.linalg.norm(Sigma_rand, "fro"))

    # Projektion auf PC1/PC2 von Σ_A
    proj_eabc = (eabc_vectors - mu_A) @ evecs_A[:, :2]
    proj_null = (null_vectors - mu_A) @ evecs_A[:, :2]

    # Erklärte Varianz
    total_var_A = float(np.sum(evals_A[evals_A > 0]))
    var_explained = evals_A[:3] / total_var_A if total_var_A > 0 else np.zeros(3)

    return {
        "mu_A": mu_A,
        "Sigma_A": Sigma_A,
        "Corr_A": Corr_A,
        "evals_A": evals_A,
        "evecs_A": evecs_A,
        "var_explained": var_explained,
        "mu_rand": mu_rand,
        "Sigma_rand": Sigma_rand,
        "Corr_rand": Corr_rand,
        "evals_rand": evals_rand,
        "evecs_rand": evecs_rand,
        "frob_abs": frob_abs,
        "frob_rel": frob_rel,
        "proj_eabc": proj_eabc,
        "proj_null": proj_null,
    }


# ---------------------------------------------------------------------------
# CSV-Export
# ---------------------------------------------------------------------------

def export_csv(
    results: dict[int, dict],
    out_path: Path,
) -> None:
    """
    Exportiert Kovarianzmatrizen, Eigenwerte und Vergleichsmetriken nach CSV.
    """
    rows: list[dict] = []

    for m, res in results.items():
        # Kovarianzmatrix-Einträge
        for i in range(6):
            for j in range(i, 6):
                rows.append({
                    "type": "cov_eabc",
                    "m": m,
                    "row": PAIR_NAMES[i],
                    "col": PAIR_NAMES[j],
                    "value": round(float(res["Sigma_A"][i, j]), 10),
                })
                rows.append({
                    "type": "cov_rand",
                    "m": m,
                    "row": PAIR_NAMES[i],
                    "col": PAIR_NAMES[j],
                    "value": round(float(res["Sigma_rand"][i, j]), 10),
                })

        # Eigenwerte
        for k, (ev, evr) in enumerate(zip(res["evals_A"], res["evals_rand"])):
            rows.append({
                "type": "eigenvalue_eabc",
                "m": m,
                "row": f"PC{k+1}",
                "col": "",
                "value": round(float(ev), 10),
            })
            rows.append({
                "type": "eigenvalue_rand",
                "m": m,
                "row": f"PC{k+1}",
                "col": "",
                "value": round(float(evr), 10),
            })

        # Frobenius
        rows.append({
            "type": "frobenius_rel",
            "m": m,
            "row": "",
            "col": "",
            "value": round(res["frob_rel"], 8),
        })
        rows.append({
            "type": "frobenius_abs",
            "m": m,
            "row": "",
            "col": "",
            "value": round(res["frob_abs"], 10),
        })

        # Erklärte Varianz
        for k, v in enumerate(res["var_explained"]):
            rows.append({
                "type": "var_explained_eabc",
                "m": m,
                "row": f"PC{k+1}",
                "col": "",
                "value": round(float(v), 8),
            })

        # Mittelwerte
        for k, (mean_a, mean_r) in enumerate(zip(res["mu_A"], res["mu_rand"])):
            rows.append({
                "type": "mean_eabc",
                "m": m,
                "row": PAIR_NAMES[k],
                "col": "",
                "value": round(float(mean_a), 10),
            })
            rows.append({
                "type": "mean_rand",
                "m": m,
                "row": PAIR_NAMES[k],
                "col": "",
                "value": round(float(mean_r), 10),
            })

    fields = ["type", "m", "row", "col", "value"]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"CSV gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Visualisierung
# ---------------------------------------------------------------------------

def _heatmap(ax, matrix: np.ndarray, title: str, vmax: float | None = None) -> None:
    """Zeichnet eine symmetrische Kovarianz-Heatmap."""
    if vmax is None:
        vmax = np.abs(matrix).max()
    if vmax == 0:
        vmax = 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    im = ax.imshow(matrix, cmap="RdBu_r", norm=norm, aspect="equal")
    ax.set_xticks(range(6))
    ax.set_yticks(range(6))
    short = ["EA", "EB", "EC", "AB", "AC", "BC"]
    ax.set_xticklabels(short, fontsize=7)
    ax.set_yticklabels(short, fontsize=7)
    ax.set_title(title, fontsize=8, pad=3)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    # Beschriftung der Zellen
    for i in range(6):
        for j in range(6):
            val = matrix[i, j]
            color = "white" if abs(val) > 0.5 * vmax else "black"
            ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                    fontsize=5.5, color=color)


def plot_results(
    results: dict[int, dict],
    out_path: Path,
) -> None:
    """
    2×3-Layout:
      Spalte 0: Heatmap Σ_A (EABC)
      Spalte 1: Heatmap Σ_A^rand (Null)
      Spalte 2: Eigenwert-Spektrum
      Zeile 0: m = WINDOW_SIZES[0]
      Zeile 1: m = WINDOW_SIZES[1]
    Zusätzlich: Scatter-Plot PC1 vs PC2 als eigenständige Zeile unten.
    """
    n_m = len(WINDOW_SIZES)
    fig = plt.figure(figsize=(18, 6 * n_m + 5))
    fig.suptitle(
        r"EABC-Pfadsignatur: Kovarianzstruktur von $\Sigma_A$ in $\Lambda^2(\mathbb{R}^4)$"
        "\nVergleich EABC-Primstrom vs. Permutations-Nullmodell",
        fontsize=12,
        fontweight="bold",
        y=0.99,
    )

    # Outer grid: obere Reihen (Heatmaps + Eigenwerte) + untere Reihe (Scatter)
    outer = gridspec.GridSpec(
        n_m + 1, 1, hspace=0.55, figure=fig,
        height_ratios=[1.0] * n_m + [1.2],
    )

    for row_idx, m in enumerate(WINDOW_SIZES):
        res = results[m]
        K = res["proj_eabc"].shape[0]
        K_null = res["proj_null"].shape[0]

        inner = gridspec.GridSpecFromSubplotSpec(
            1, 3, subplot_spec=outer[row_idx], wspace=0.40,
        )
        ax_cov = fig.add_subplot(inner[0])
        ax_null = fig.add_subplot(inner[1])
        ax_ev = fig.add_subplot(inner[2])

        # Gemeinsamer Farbbereich für beide Heatmaps
        vmax = max(
            np.abs(res["Sigma_A"]).max(),
            np.abs(res["Sigma_rand"]).max(),
        )

        _heatmap(ax_cov, res["Sigma_A"],
                 rf"$\Sigma_A$ (EABC, $m={m}$, $K={K}$)", vmax=vmax)
        _heatmap(ax_null, res["Sigma_rand"],
                 rf"$\Sigma_A^\mathrm{{rand}}$ (Null, $m={m}$)", vmax=vmax)

        # Eigenwert-Spektrum
        ev_idx = np.arange(1, 7)
        ax_ev.bar(ev_idx - 0.2, res["evals_A"], width=0.38, color=EABC_COLOR,
                  alpha=0.75, label="EABC")
        ax_ev.bar(ev_idx + 0.2, res["evals_rand"], width=0.38, color=NULL_COLOR,
                  alpha=0.75, label="Null")
        ax_ev.set_xticks(ev_idx)
        ax_ev.set_xticklabels([f"PC{i}" for i in ev_idx], fontsize=8)
        ax_ev.set_ylabel("Eigenwert", fontsize=8)
        ax_ev.set_title(
            rf"Eigenwert-Spektrum von $\Sigma_A$ ($m={m}$)"
            f"\nFrobenius rel.: {res['frob_rel']:.4f}",
            fontsize=8,
        )
        ax_ev.legend(fontsize=8)
        ax_ev.grid(True, alpha=0.25, axis="y")
        ax_ev.tick_params(labelsize=8)

        # Erklärte Varianz annotieren
        ve = res["var_explained"]
        ax_ev.text(
            0.97, 0.97,
            f"PC1: {ve[0]:.1%}\nPC2: {ve[1]:.1%}\nPC3: {ve[2]:.1%}",
            transform=ax_ev.transAxes,
            fontsize=7.5,
            ha="right", va="top",
            fontfamily="monospace",
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none"),
        )

    # Untere Zeile: Scatter PC1 vs PC2 für beide Fenstergrößen
    scatter_inner = gridspec.GridSpecFromSubplotSpec(
        1, n_m, subplot_spec=outer[n_m], wspace=0.40,
    )
    for col_idx, m in enumerate(WINDOW_SIZES):
        res = results[m]
        ax_sc = fig.add_subplot(scatter_inner[col_idx])

        pn = res["proj_null"]
        pe = res["proj_eabc"]

        ax_sc.scatter(
            pn[:, 0], pn[:, 1],
            color=NULL_COLOR, alpha=SCATTER_ALPHA * 0.6,
            s=8, label=f"Null ($n={len(pn)}$)", rasterized=True,
        )
        ax_sc.scatter(
            pe[:, 0], pe[:, 1],
            color=EABC_COLOR, alpha=SCATTER_ALPHA,
            s=10, label=f"EABC ($n={len(pe)}$)", rasterized=True,
        )

        ax_sc.axhline(0, color="black", linewidth=0.7, alpha=0.4)
        ax_sc.axvline(0, color="black", linewidth=0.7, alpha=0.4)
        ax_sc.set_xlabel(
            rf"PC1 ({res['var_explained'][0]:.1%} Varianz)", fontsize=8
        )
        ax_sc.set_ylabel(
            rf"PC2 ({res['var_explained'][1]:.1%} Varianz)", fontsize=8
        )
        ax_sc.set_title(
            rf"Projektion auf Hauptachsen von $\Sigma_A$ ($m={m}$)"
            "\nBlau = EABC, Grau = Permutations-Null",
            fontsize=8,
        )
        ax_sc.legend(fontsize=7.5, markerscale=1.5)
        ax_sc.grid(True, alpha=0.20)
        ax_sc.tick_params(labelsize=8)

        # PC1-Ladungsvektor annotieren
        pc1 = res["evecs_A"][:, 0]
        loading_str = "PC1-Ladung:\n" + "\n".join(
            f"{PAIR_NAMES[i]}: {pc1[i]:+.3f}" for i in range(6)
        )
        ax_sc.text(
            0.02, 0.97, loading_str,
            transform=ax_sc.transAxes,
            fontsize=6.5, va="top",
            fontfamily="monospace",
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none"),
        )

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Ergebnistabelle (Terminal)
# ---------------------------------------------------------------------------

def print_summary(results: dict[int, dict]) -> None:
    print()
    print("=" * 80)
    print("KOVARIANZSTRUKTUR-ANALYSE: Λ²(ℝ⁴) — antisymmetrische Level-2-Signatur")
    print("=" * 80)

    for m, res in results.items():
        K = res["proj_eabc"].shape[0]
        K_null = res["proj_null"].shape[0]
        print(f"\n  Fenstergröße m = {m:>6}  (K_EABC = {K}, K_null = {K_null})")
        print()
        print("  Mittelwertvektor μ_A (EABC):")
        for i, (name, val) in enumerate(zip(PAIR_NAMES, res["mu_A"])):
            print(f"    a_{name} = {val:+.7f}")
        print()
        print("  Eigenwerte von Σ_A (EABC):      Eigenwerte Σ_A^rand:")
        for i, (ev, evr) in enumerate(zip(res["evals_A"], res["evals_rand"])):
            bar_eabc = "█" * int(ev / res["evals_A"][0] * 15)
            bar_rand = "█" * int(evr / max(res["evals_rand"][0], 1e-12) * 15)
            ve_str = f"({res['var_explained'][i]:.1%})" if i < 3 else ""
            print(
                f"    PC{i+1}: {ev:9.6f} {bar_eabc:<15}  {evr:9.6f} {bar_rand:<15}  {ve_str}"
            )
        print()
        print(f"  Frobenius-Abstand (absolut): {res['frob_abs']:.6f}")
        print(f"  Frobenius-Abstand (relativ): {res['frob_rel']:.6f}")
        if res["frob_rel"] > 0.10:
            print("  → Σ_A weicht DEUTLICH vom Nullmodell ab (>10% relativ)")
        elif res["frob_rel"] > 0.03:
            print("  → Σ_A weicht MODERAT vom Nullmodell ab (3–10% relativ)")
        else:
            print("  → Σ_A ist konsistent mit dem Nullmodell (<3% relativ)")

    print()
    print("=" * 80)
    print("GEOMETRISCHE INTERPRETATION:")
    print()
    for m, res in results.items():
        print(f"  m = {m}: Dominante Richtungen in Λ²(ℝ⁴):")
        for k in range(2):
            ev_k = res["evecs_A"][:, k]
            dominant = sorted(
                range(6), key=lambda i: abs(ev_k[i]), reverse=True
            )
            top2 = [f"{PAIR_NAMES[i]} ({ev_k[i]:+.3f})" for i in dominant[:3]]
            print(f"    PC{k+1} (λ={res['evals_A'][k]:.5f}): " + ", ".join(top2))
    print("=" * 80)


# ---------------------------------------------------------------------------
# Hauptprogramm
# ---------------------------------------------------------------------------

def main() -> None:
    out_dir = Path(__file__).resolve().parent
    rng = np.random.default_rng(RANDOM_SEED)

    print("=" * 80)
    print("EABC-Pfadsignatur: Kovarianzstruktur der antisymmetrischen Level-2-Signatur")
    print(f"N_PRIMES = {N_PRIMES:,}  |  Fenstergrößen: {WINDOW_SIZES}  |  B_RAND = {B_RAND}")
    print("=" * 80)

    # Primzahlwort
    print(f"\nGeneriere EABC-Codewort der ersten {N_PRIMES:,} Primzahlen > 3 …")
    word = build_eabc_word(N_PRIMES)
    print(f"  Wortlänge: {len(word):,}")
    from collections import Counter
    counts = Counter(word.tolist())
    for code, label in enumerate(("E", "A", "B", "C")):
        n = counts[code]
        print(f"  N_{label} = {n:>7}  ({100 * n / len(word):.3f}%)")

    # Analyse pro Fenstergröße
    results: dict[int, dict] = {}

    for m in WINDOW_SIZES:
        print(f"\n{'─'*60}")
        print(f"  Fenstergröße m = {m}")
        print(f"{'─'*60}")

        print(f"  Sammle EABC-Fenstervektoren (m={m}) …")
        eabc_vecs = collect_window_vectors(word, m)
        print(f"  → {eabc_vecs.shape[0]} Fenster, {eabc_vecs.shape[1]} Komponenten")

        print(f"  Sammle Nullmodell-Vektoren (B={B_RAND} Permutationen/Fenster) …")
        null_vecs = collect_null_vectors(word, m, B_RAND, rng)
        print(f"  → {null_vecs.shape[0]} Null-Vektoren")

        print(f"  Berechne Kovarianzstruktur …")
        res = covariance_analysis(eabc_vecs, null_vecs)
        results[m] = res

    # Ergebnisse ausgeben
    print_summary(results)

    # CSV
    csv_path = out_dir / "eabc_covariance_results.csv"
    export_csv(results, csv_path)

    # Plot
    print("\nErstelle Kovarianz-Plot …")
    plot_path = out_dir / "eabc_covariance_plot.png"
    plot_results(results, plot_path)

    print("\nFertig.")


if __name__ == "__main__":
    main()
