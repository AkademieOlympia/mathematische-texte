#!/usr/bin/env python3
"""
EABC-Pfadasymmetrie: Permutations-Nullmodell-Test (altsig2)
============================================================

Observable: a_XY = (N_XY - N_YX) / (N_XY + N_YX) ∈ [-1, 1]
  - skalierungsinvariant gegenüber Fenstergröße
  - NaN wenn N_XY + N_YX = 0

Nullhypothese H₀:
    a_XY verhält sich wie unter einer zufälligen Permutation des Fensters
    mit denselben Klassenhäufigkeiten (gleiche Marginalen, zufällige Reihenfolge).

Alternative H₁:
    a_EX besitzt eine systematische Verschiebung gegenüber H₀.

p-Wert:
    Permutationsbasiert, einseitig (Anteil der Permutationen mit a_XY(perm) ≤ a_XY(window)).
    Keine Normalverteilungsannahme erforderlich.

Stichprobenstrategie:
    - a_XY wird für alle disjunkten Fenster berechnet.
    - Permutationstest nur für jedes PERM_STRIDE-te Fenster (rechenintensiv).

Ausgabe:
    eabc_altsig2_permtest_results.csv
    eabc_altsig2_permtest_plot.png
"""

from __future__ import annotations

import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eabc_signature_analysis import (
    sieve_primes,
    class_of,
)

# ---------------------------------------------------------------------------
# Konstanten
# ---------------------------------------------------------------------------

WINDOW_SIZES = [1000, 5000, 10000]
E_PAIRS = [("E", "A"), ("E", "B"), ("E", "C")]
PAIR_NAMES = ["EA", "EB", "EC"]

N_PRIMES = 200_000
B_PERMS = 500
PERM_STRIDE = 5  # Permutationstest nur für jedes PERM_STRIDE-te Fenster

# Interne Buchstaben-Kodierung (int8 für schnelles numpy-Shuffling)
LABEL_TO_CODE: dict[str, int] = {"E": 0, "A": 1, "B": 2, "C": 3}
CODE_TO_LABEL: dict[int, str] = {v: k for k, v in LABEL_TO_CODE.items()}

COLORS_BY_M = {1000: "#1f77b4", 5000: "#2ca02c", 10000: "#d62728"}
NULL_COLOR = "#aaaaaa"


# ---------------------------------------------------------------------------
# Primzahlwort erzeugen
# ---------------------------------------------------------------------------

def build_eabc_word(n_primes: int) -> list[str]:
    """Erstellt das EABC-Wort der ersten n_primes Primzahlen > 3."""
    limit = 3_000_000
    while True:
        primes = sieve_primes(limit)
        word = [class_of(p) for p in primes if p > 3 and class_of(p) is not None]
        if len(word) >= n_primes:
            return word[:n_primes]
        limit = int(limit * 1.5)


def word_to_codes(word: list[str]) -> np.ndarray:
    """Konvertiert Buchstabenliste in int8-Array für schnelles Shuffling."""
    return np.array([LABEL_TO_CODE[c] for c in word], dtype=np.int8)


# ---------------------------------------------------------------------------
# Kern-Observable: normierte Antisymmetrie
# ---------------------------------------------------------------------------

def compute_nxy_nyx(word_codes: np.ndarray, x_code: int, y_code: int) -> tuple[int, int]:
    """
    Berechnet N_XY und N_YX in O(m) mit numpy-Vektorisierung.

    N_XY = #{(i,j) : i < j, word[i]=X, word[j]=Y}
    N_YX = #{(i,j) : i < j, word[i]=Y, word[j]=X}
    """
    is_x = (word_codes == x_code).view(np.uint8)
    is_y = (word_codes == y_code).view(np.uint8)
    # Für jede Y-Position j: Anzahl X's vor j = cumsum_x[j] - is_x[j]
    cumx = np.cumsum(is_x, dtype=np.int64)
    cumy = np.cumsum(is_y, dtype=np.int64)
    n_xy = int(np.dot(cumx - is_x, is_y))
    n_yx = int(np.dot(cumy - is_y, is_x))
    return n_xy, n_yx


def normalized_antisym(n_xy: int, n_yx: int) -> float:
    """a_XY = (N_XY - N_YX) / (N_XY + N_YX) ∈ [-1, 1], NaN wenn beide 0."""
    total = n_xy + n_yx
    if total == 0:
        return float("nan")
    return (n_xy - n_yx) / total


# ---------------------------------------------------------------------------
# Permutations-p-Wert
# ---------------------------------------------------------------------------

def permutation_pvalue(
    window_codes: np.ndarray,
    x_code: int,
    y_code: int,
    obs: float,
    B: int,
    rng: np.random.Generator,
) -> tuple[float, np.ndarray]:
    """
    Einseitiger Permutations-p-Wert.

    p = Anteil(a_XY(perm) ≤ obs) über B zufällige Permutationen.

    Rückgabe: (p_value, array der Permutations-a_XY-Werte)
    """
    if np.isnan(obs):
        return float("nan"), np.array([])

    perm = window_codes.copy()
    perm_vals = np.empty(B, dtype=np.float64)

    for b in range(B):
        rng.shuffle(perm)
        n_xy, n_yx = compute_nxy_nyx(perm, x_code, y_code)
        perm_vals[b] = normalized_antisym(n_xy, n_yx)

    valid = perm_vals[~np.isnan(perm_vals)]
    if len(valid) == 0:
        return float("nan"), np.array([])

    p_val = float(np.sum(valid <= obs) / len(valid))
    return p_val, valid


# ---------------------------------------------------------------------------
# Sliding-Window-Analyse
# ---------------------------------------------------------------------------

def run_analysis(
    word: list[str],
    rng: np.random.Generator,
) -> tuple[list[dict], dict]:
    """
    Führt die vollständige Sliding-Window-Analyse durch.

    Rückgabe:
        - summary_rows: Liste von Zusammenfassungsdicts pro (m, pair)
        - window_data: Dict für Visualisierung {(m, pair_name): {'a': [...], 'p': [...], 'null': [...]}}
    """
    word_codes = word_to_codes(word)
    n = len(word_codes)

    summary_rows: list[dict] = []
    window_data: dict[tuple[int, str], dict] = {}

    for m in WINDOW_SIZES:
        # Disjunkte Fensterpositionen
        starts = list(range(0, n - m + 1, m))
        K = len(starts)
        print(f"\n  m = {m:>6}, K = {K} disjunkte Fenster")

        for pair_name, (X, Y) in zip(PAIR_NAMES, E_PAIRS):
            x_code = LABEL_TO_CODE[X]
            y_code = LABEL_TO_CODE[Y]

            a_vals: list[float] = []
            p_vals: list[float] = []
            null_samples: list[float] = []

            for k, start in enumerate(starts):
                win = word_codes[start : start + m]
                n_xy, n_yx = compute_nxy_nyx(win, x_code, y_code)
                obs = normalized_antisym(n_xy, n_yx)
                a_vals.append(obs)

                # Permutationstest nur für jedes PERM_STRIDE-te Fenster
                if k % PERM_STRIDE == 0:
                    p, perms = permutation_pvalue(win, x_code, y_code, obs, B_PERMS, rng)
                    p_vals.append(p)
                    # Für Nullverteilungs-Plot: begrenzte Stichprobe
                    if len(null_samples) < 5000 and len(perms) > 0:
                        take = min(len(perms), 5000 - len(null_samples))
                        null_samples.extend(perms[:take].tolist())

                if (k + 1) % 20 == 0 or k == K - 1:
                    n_perm_done = sum(1 for kk in range(k + 1) if kk % PERM_STRIDE == 0)
                    print(
                        f"    A_{pair_name}: {k+1:>4}/{K} Fenster, "
                        f"{n_perm_done} Permutationstests",
                        end="\r",
                        flush=True,
                    )

            print()

            a_arr = np.array([v for v in a_vals if not np.isnan(v)], dtype=np.float64)
            p_arr = np.array([v for v in p_vals if not np.isnan(v)], dtype=np.float64)

            mean_a = float(np.mean(a_arr)) if len(a_arr) > 0 else float("nan")
            std_a = float(np.std(a_arr, ddof=1)) if len(a_arr) > 1 else float("nan")
            median_p = float(np.median(p_arr)) if len(p_arr) > 0 else float("nan")
            frac_sig = float(np.mean(p_arr < 0.05)) if len(p_arr) > 0 else float("nan")

            summary_rows.append({
                "m": m,
                "pair": pair_name,
                "K": K,
                "K_perm": len(p_arr),
                "mean_a": round(mean_a, 8),
                "std_a": round(std_a, 8),
                "median_p": round(median_p, 6),
                "frac_sig": round(frac_sig, 4),
            })

            window_data[(m, pair_name)] = {
                "a": a_vals,
                "p": p_vals,
                "null": null_samples,
            }

            print(
                f"    A_{pair_name}: mean_a={mean_a:+.6f}, std_a={std_a:.6f}, "
                f"median_p={median_p:.4f}, frac_sig={frac_sig:.3f}"
            )

    return summary_rows, window_data


# ---------------------------------------------------------------------------
# CSV-Export
# ---------------------------------------------------------------------------

CSV_FIELDS = [
    "m", "pair", "K", "K_perm",
    "mean_a", "std_a", "median_p", "frac_sig",
]


def export_csv(summary_rows: list[dict], out_path: Path) -> None:
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"CSV gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Ergebnistabelle
# ---------------------------------------------------------------------------

def print_results_table(summary_rows: list[dict]) -> None:
    print()
    print("=" * 100)
    print("ERGEBNISTABELLE: Normierte Antisymmetrie a_EX = (N_EX - N_XE) / (N_EX + N_XE)")
    print("Nullmodell: Permutation mit gleichen Klassenhäufigkeiten")
    print("=" * 100)

    header = (
        f"{'m':>7} | {'K':>5} | {'K_perm':>6} | "
        f"{'mean_a_EA':>11} | {'std_EA':>9} | {'med_p_EA':>9} | {'frac<0.05':>9} | "
        f"{'mean_a_EB':>11} | {'std_EB':>9} | {'med_p_EB':>9} | {'frac<0.05':>9} | "
        f"{'mean_a_EC':>11} | {'std_EC':>9} | {'med_p_EC':>9} | {'frac<0.05':>9}"
    )
    print(header)
    print("-" * 100)

    # Gruppiere nach m
    by_m: dict[int, dict[str, dict]] = {}
    for row in summary_rows:
        m = row["m"]
        if m not in by_m:
            by_m[m] = {}
        by_m[m][row["pair"]] = row

    for m in WINDOW_SIZES:
        pairs = by_m.get(m, {})
        ea = pairs.get("EA", {})
        eb = pairs.get("EB", {})
        ec = pairs.get("EC", {})

        def fmt_a(v):
            return f"{v:+11.6f}" if not (isinstance(v, float) and np.isnan(v)) else "        nan"

        def fmt_f(v):
            return f"{v:9.6f}" if not (isinstance(v, float) and np.isnan(v)) else "      nan"

        def fmt_frac(v):
            return f"{v:9.4f}" if not (isinstance(v, float) and np.isnan(v)) else "      nan"

        # Kernbefund-Flag
        flags = []
        for p_name, pr in [("EA", ea), ("EB", eb), ("EC", ec)]:
            mp = pr.get("median_p", float("nan"))
            fs = pr.get("frac_sig", float("nan"))
            if not np.isnan(mp) and mp < 0.05 and not np.isnan(fs) and fs > 0.5:
                flags.append(f"A_{p_name}✓")
        flag_str = "  ← " + ", ".join(flags) if flags else ""

        line = (
            f"{m:>7} | {ea.get('K','?'):>5} | {ea.get('K_perm','?'):>6} | "
            f"{fmt_a(ea.get('mean_a', float('nan')))} | "
            f"{fmt_f(ea.get('std_a', float('nan')))} | "
            f"{fmt_f(ea.get('median_p', float('nan')))} | "
            f"{fmt_frac(ea.get('frac_sig', float('nan')))} | "
            f"{fmt_a(eb.get('mean_a', float('nan')))} | "
            f"{fmt_f(eb.get('std_a', float('nan')))} | "
            f"{fmt_f(eb.get('median_p', float('nan')))} | "
            f"{fmt_frac(eb.get('frac_sig', float('nan')))} | "
            f"{fmt_a(ec.get('mean_a', float('nan')))} | "
            f"{fmt_f(ec.get('std_a', float('nan')))} | "
            f"{fmt_f(ec.get('median_p', float('nan')))} | "
            f"{fmt_frac(ec.get('frac_sig', float('nan')))}"
            f"{flag_str}"
        )
        print(line)

    print("=" * 100)


def print_core_finding(summary_rows: list[dict]) -> None:
    print()
    print("=" * 70)
    print("KERNBEFUND: Systematische Verschiebung gegenüber Permutations-H₀?")
    print("Kriterien: median_p < 0.05 UND frac_sig > 0.50")
    print("=" * 70)

    all_confirmed = True
    any_confirmed = False
    for pair_name in PAIR_NAMES:
        pair_rows = [r for r in summary_rows if r["pair"] == pair_name]
        confirmed_m = []
        for r in pair_rows:
            mp = r.get("median_p", float("nan"))
            fs = r.get("frac_sig", float("nan"))
            if not np.isnan(mp) and not np.isnan(fs) and mp < 0.05 and fs > 0.50:
                confirmed_m.append(r["m"])

        if len(confirmed_m) == len(pair_rows) and len(pair_rows) > 0:
            status = "BESTÄTIGT"
            any_confirmed = True
        elif len(confirmed_m) > 0:
            status = f"TEILWEISE (m={confirmed_m})"
            any_confirmed = True
            all_confirmed = False
        else:
            status = "NEIN"
            all_confirmed = False

        mp_vals = [r["median_p"] for r in pair_rows]
        fs_vals = [r["frac_sig"] for r in pair_rows]
        print(f"  A_{pair_name}: {status}")
        print(f"    median_p = {[round(v,4) for v in mp_vals]}")
        print(f"    frac_sig = {[round(v,4) for v in fs_vals]}")

    print()
    if all_confirmed and any_confirmed:
        print("  → POSITIVER BEFUND: Alle drei E-Paare weichen systematisch")
        print("    vom Permutations-Nullmodell ab.")
        print("    Erster belastbarer Befund einer arithmetischen")
        print("    Pfadasymmetrie zweiter Ordnung.")
    elif any_confirmed:
        print("  → TEILBEFUND: Einige Paare zeigen systematische Verschiebung.")
        print("    Ergebnis ist nicht konsistent über alle Paare/Fenstergrößen.")
    else:
        print("  → NEGATIVER BEFUND: Keine systematische Verschiebung")
        print("    gegenüber dem Permutations-Nullmodell nachweisbar.")

    print("=" * 70)


# ---------------------------------------------------------------------------
# Visualisierung
# ---------------------------------------------------------------------------

def plot_results(
    summary_rows: list[dict],
    window_data: dict[tuple[int, str], dict],
    out_path: Path,
) -> None:
    """
    2×2-Layout:
      [0,0..2]: Histogramm der a_EX-Verteilung + Nullmodell (je ein Subplot pro Paar)
      [1,0]:   a_EC
      [1,1]:   Zusammenfassung: frac_sig und median_p als Balkendiagramm
    """
    fig = plt.figure(figsize=(16, 12))
    fig.suptitle(
        r"EABC-Pfadasymmetrie: $a_{EX} = (N_{EX}-N_{XE})/(N_{EX}+N_{XE})$ — Permutations-Nullmodell"
        "\nBlau/Grün/Rot = EABC (m=1000/5000/10000), Grau = Permutations-Nullverteilung",
        fontsize=11,
        fontweight="bold",
    )

    gs = gridspec.GridSpec(2, 3, hspace=0.50, wspace=0.35)
    axes_hist = [
        fig.add_subplot(gs[0, 0]),  # EA
        fig.add_subplot(gs[0, 1]),  # EB
        fig.add_subplot(gs[0, 2]),  # EC
    ]
    ax_frac = fig.add_subplot(gs[1, 0:2])  # frac_sig Balkendiagramm
    ax_mp = fig.add_subplot(gs[1, 2])       # median_p

    # --- Histogramme der a_EX-Verteilung ---
    for ax, pair_name in zip(axes_hist, PAIR_NAMES):
        ax.axvline(0, color="black", linewidth=1.0, zorder=5)

        null_plotted = False
        for m, color in COLORS_BY_M.items():
            key = (m, pair_name)
            if key not in window_data:
                continue
            data = window_data[key]
            a_vals = [v for v in data["a"] if not np.isnan(v)]
            null_vals = data["null"]

            if not a_vals:
                continue

            bins = np.linspace(-1, 1, 51)

            # Nullverteilung (einmal, grau)
            if null_vals and not null_plotted:
                ax.hist(
                    null_vals,
                    bins=bins,
                    density=True,
                    color=NULL_COLOR,
                    alpha=0.5,
                    label="Null (Perm.)",
                    zorder=2,
                )
                null_plotted = True

            # EABC-Verteilung
            ax.hist(
                a_vals,
                bins=bins,
                density=True,
                color=color,
                alpha=0.55,
                label=f"m={m}",
                zorder=3,
            )

        ax.set_title(
            rf"$a_{{E{pair_name[1]}}}$ = $(N_{{E{pair_name[1]}}}-N_{{{pair_name[1]}E}})/(N_{{E{pair_name[1]}}}+N_{{{pair_name[1]}E}})$",
            fontsize=9,
        )
        ax.set_xlabel(rf"$a_{{E{pair_name[1]}}}$", fontsize=8)
        ax.set_ylabel("Dichte", fontsize=8)
        ax.legend(fontsize=7.5, loc="upper left")
        ax.grid(True, alpha=0.25)
        ax.tick_params(labelsize=8)

        # Annotiere mean und median_p
        annot_lines = []
        for m in WINDOW_SIZES:
            row = next(
                (r for r in summary_rows if r["m"] == m and r["pair"] == pair_name),
                None,
            )
            if row:
                ma = row["mean_a"]
                mp = row["median_p"]
                fs = row["frac_sig"]
                ma_str = f"{ma:+.5f}" if not np.isnan(ma) else "nan"
                mp_str = f"{mp:.4f}" if not np.isnan(mp) else "nan"
                fs_str = f"{fs:.3f}" if not np.isnan(fs) else "nan"
                annot_lines.append(f"m={m}: μ={ma_str}, med_p={mp_str}, sig={fs_str}")
        ax.text(
            0.02, 0.97,
            "\n".join(annot_lines),
            transform=ax.transAxes,
            fontsize=6.5,
            verticalalignment="top",
            fontfamily="monospace",
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none"),
        )

    # --- frac_sig Balkendiagramm ---
    x_base = np.arange(len(PAIR_NAMES))
    width = 0.25
    for i, (m, color) in enumerate(COLORS_BY_M.items()):
        frac_vals = []
        for pair_name in PAIR_NAMES:
            row = next(
                (r for r in summary_rows if r["m"] == m and r["pair"] == pair_name),
                None,
            )
            frac_vals.append(row["frac_sig"] if row else float("nan"))

        ax_frac.bar(
            x_base + (i - 1) * width,
            frac_vals,
            width=width * 0.9,
            color=color,
            alpha=0.75,
            label=f"m={m}",
        )

    ax_frac.axhline(0.05, color="black", linewidth=1.2, linestyle="--", label="α=0.05", zorder=5)
    ax_frac.set_xticks(x_base)
    ax_frac.set_xticklabels([f"A_E{p[1]}" for p in PAIR_NAMES], fontsize=9)
    ax_frac.set_ylabel("Anteil Fenster mit p < 0.05", fontsize=8)
    ax_frac.set_title("Anteil signifikanter Fenster (Permutationstest, einseitig)", fontsize=9)
    ax_frac.set_ylim(0, 1.05)
    ax_frac.legend(fontsize=8)
    ax_frac.grid(True, alpha=0.25, axis="y")
    ax_frac.tick_params(labelsize=8)

    # --- median_p Balkendiagramm ---
    for i, (m, color) in enumerate(COLORS_BY_M.items()):
        mp_vals = []
        for pair_name in PAIR_NAMES:
            row = next(
                (r for r in summary_rows if r["m"] == m and r["pair"] == pair_name),
                None,
            )
            mp_vals.append(row["median_p"] if row else float("nan"))

        ax_mp.bar(
            x_base + (i - 1) * width,
            mp_vals,
            width=width * 0.9,
            color=color,
            alpha=0.75,
            label=f"m={m}",
        )

    ax_mp.axhline(0.05, color="black", linewidth=1.2, linestyle="--", label="α=0.05", zorder=5)
    ax_mp.set_xticks(x_base)
    ax_mp.set_xticklabels([f"A_E{p[1]}" for p in PAIR_NAMES], fontsize=9)
    ax_mp.set_ylabel("Median p-Wert", fontsize=8)
    ax_mp.set_title("Median Permutations-p-Wert", fontsize=9)
    ax_mp.set_ylim(0, 1.05)
    ax_mp.legend(fontsize=8)
    ax_mp.grid(True, alpha=0.25, axis="y")
    ax_mp.tick_params(labelsize=8)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Hauptprogramm
# ---------------------------------------------------------------------------

def main() -> None:
    out_dir = Path(__file__).resolve().parent
    rng = np.random.default_rng(42)

    print("=" * 70)
    print("EABC-Pfadasymmetrie: Permutations-Nullmodell-Test (altsig2)")
    print(f"Observable: a_XY = (N_XY - N_YX) / (N_XY + N_YX) ∈ [-1, 1]")
    print(f"Fenstergrössen: {WINDOW_SIZES}")
    print(f"Permutationen pro Fenster: B={B_PERMS}, Schrittweite: {PERM_STRIDE}")
    print("=" * 70)

    # Primzahlwort erzeugen
    print(f"\nGeneriere EABC-Wort der ersten {N_PRIMES:,} Primzahlen > 3 …")
    word = build_eabc_word(N_PRIMES)
    print(f"  Wortlänge: {len(word):,}")

    from collections import Counter
    counts = Counter(word)
    for lab in ("E", "A", "B", "C"):
        print(f"  N_{lab} = {counts[lab]:>7}  ({100*counts[lab]/len(word):.3f}%)")

    # Analyse
    print("\nStarte Sliding-Window-Analyse …")
    summary_rows, window_data = run_analysis(word, rng)

    # Ergebnistabelle
    print_results_table(summary_rows)
    print_core_finding(summary_rows)

    # CSV
    csv_path = out_dir / "eabc_altsig2_permtest_results.csv"
    export_csv(summary_rows, csv_path)

    # Plot
    print("\nErstelle Visualisierung …")
    plot_path = out_dir / "eabc_altsig2_permtest_plot.png"
    plot_results(summary_rows, window_data, plot_path)

    print("\nFertig.")


if __name__ == "__main__":
    main()
