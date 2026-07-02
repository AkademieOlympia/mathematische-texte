#!/usr/bin/env python3
"""
EABC-Pfadasymmetrie-Test zweiter Ordnung
=========================================

Observable: A_XY = N_XY - N_YX  (antisymmetrische Level-2-Pfadsignatur)

Nullhypothese:  H₀: E[A_EX^(k,m)] = 0
Alternative:    H₁: E[A_EX^(k,m)] < 0   für X ∈ {A, B, C}

Zwei Modi werden explizit getrennt:
    Modus A (disjoint):  step = m   — z-Test statistisch valide
    Modus B (overlap):   step = m//10 — nur deskriptiv

Normierung: a_XY = A_XY / (Q*(Q-1)/2)  damit Fenstergrößen vergleichbar.

Kernfrage: Sind z_EA, z_EB, z_EC < −2 bei disjunkten Fenstern?

Ausgabe:
    eabc_asymmetry_results.csv
    eabc_asymmetry_plot.png
"""

from __future__ import annotations

import csv
import sys
from math import sqrt
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eabc_signature_analysis import (
    LABELS,
    eabc_signature,
    sieve_primes,
    class_of,
)

# ---------------------------------------------------------------------------
# Konstanten
# ---------------------------------------------------------------------------

WINDOW_SIZES = [1000, 5000, 10000]
E_PAIRS = [("E", "A"), ("E", "B"), ("E", "C")]
PAIR_KEYS = ["A_EA", "A_EB", "A_EC"]

COLORS_EABC = {1000: "#1f77b4", 5000: "#1a6fa8", 10000: "#155d8e"}
COLORS_EABC_OVERLAP = {1000: "#aec7e8", 5000: "#93b8da", 10000: "#7aaacb"}
COLORS_RAND = {1000: "#aaaaaa", 5000: "#888888", 10000: "#666666"}
COLORS_RAND_OVERLAP = {1000: "#dddddd", 5000: "#cccccc", 10000: "#bbbbbb"}


# ---------------------------------------------------------------------------
# Primzahlwort + Zufallsbaseline erzeugen
# ---------------------------------------------------------------------------

def build_eabc_word(n_primes: int) -> list[str]:
    """Erstellt EABC-Wort der ersten n_primes Primzahlen > 3."""
    limit = 3_000_000
    while True:
        primes = sieve_primes(limit)
        word = [class_of(p) for p in primes if p > 3 and class_of(p) is not None]
        if len(word) >= n_primes:
            return word[:n_primes]
        limit = int(limit * 1.5)


def marginal_probs(word: list[str]) -> dict[str, float]:
    """Berechnet empirische Häufigkeiten der EABC-Buchstaben."""
    n = len(word)
    return {lab: word.count(lab) / n for lab in LABELS}


def random_word_matched_marginals(
    word: list[str],
    seed: int = 42,
    rng: np.random.Generator | None = None,
) -> list[str]:
    """
    Zufallswort gleicher Länge mit gleichen Marginalwahrscheinlichkeiten
    wie das Originalwort (aus tatsächlichen Primzahlfrequenzen).
    Kein Shuffle — unabhängig gezogen, damit Fenster wirklich i.i.d. sind.
    """
    if rng is None:
        rng = np.random.default_rng(seed)
    probs = marginal_probs(word)
    labels_arr = list(LABELS)
    p_arr = np.array([probs[lab] for lab in labels_arr])
    indices = rng.choice(len(labels_arr), size=len(word), p=p_arr)
    return [labels_arr[i] for i in indices]


# ---------------------------------------------------------------------------
# Kern-Observable: antisymmetrische Level-2-Signatur
# ---------------------------------------------------------------------------

def antisymmetry_vals(window: list[str]) -> dict[str, float]:
    """
    Berechnet A_EX = N_EX - N_XE für X ∈ {A, B, C} sowie A_XY für alle Paare.
    Gibt rohe ganzzahlige Werte zurück.
    """
    sig = eabc_signature(window, max_level=2)
    l2 = sig.level2
    return {
        "A_EA": l2[("E", "A")] - l2[("A", "E")],
        "A_EB": l2[("E", "B")] - l2[("B", "E")],
        "A_EC": l2[("E", "C")] - l2[("C", "E")],
    }


def normalise(a_raw: float, m: int) -> float:
    """Normierung: a = A / (m*(m-1)/2) — macht Fenstergrößen vergleichbar."""
    denom = m * (m - 1) / 2.0
    return a_raw / denom if denom > 0 else float("nan")


# ---------------------------------------------------------------------------
# Sliding-Window-Berechnung
# ---------------------------------------------------------------------------

def sliding_window(
    word: list[str],
    window_size: int,
    step: int,
    label: str = "",
) -> list[dict]:
    """
    Berechnet die normierten A_EX-Werte für alle Fenster der gegebenen
    Schrittweite.  Rückgabe: Liste von Dicts mit start, end, rohen und
    normierten A_EX-Werten.
    """
    rows: list[dict] = []
    n = len(word)
    m = window_size
    for start in range(0, n - m + 1, step):
        w = word[start : start + m]
        raw = antisymmetry_vals(w)
        row: dict = {
            "start": start,
            "end": start + m,
            "label": label,
        }
        for pk in PAIR_KEYS:
            row[pk] = raw[pk]
            row[pk + "_norm"] = normalise(raw[pk], m)
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Statistische Auswertung
# ---------------------------------------------------------------------------

def compute_statistics(
    rows: list[dict],
    m: int,
    mode: str,
    source: str,
) -> list[dict]:
    """
    Für jeden Pair-Key:
        K, mean_A, sd_A, z (nur für disjoint valide), p (einseitig, H₁:<0),
        mean_a (normiert).

    z-Test: z = mean_A / (sd_A / sqrt(K))
    p-Wert: P(Z ≤ z) unter N(0,1)  (einseitig, linke Seite)
    """
    K = len(rows)
    results = []
    for pk in PAIR_KEYS:
        pair = pk[2:]  # 'EA', 'EB', 'EC'
        vals = np.array([r[pk] for r in rows], dtype=float)
        vals_norm = np.array([r[pk + "_norm"] for r in rows], dtype=float)

        mu = float(np.mean(vals))
        mu_norm = float(np.mean(vals_norm))
        sd = float(np.std(vals, ddof=1)) if K > 1 else float("nan")
        se = sd / sqrt(K) if (K > 1 and sd > 0) else float("nan")

        z = mu / se if not (np.isnan(se) or se == 0) else float("nan")
        # p-Wert (einseitig links: H₁: < 0)
        p = float(norm.cdf(z)) if not np.isnan(z) else float("nan")

        results.append({
            "m": m,
            "mode": mode,
            "source": source,
            "pair": pair,
            "K": K,
            "mean_A": round(mu, 6),
            "sd_A": round(sd, 4),
            "mean_a": round(mu_norm, 8),
            "z": round(z, 4) if not np.isnan(z) else float("nan"),
            "p": round(p, 6) if not np.isnan(p) else float("nan"),
            "z_lt_minus2": bool(z < -2) if not np.isnan(z) else False,
        })
    return results


# ---------------------------------------------------------------------------
# CSV-Export
# ---------------------------------------------------------------------------

RESULT_FIELDS = [
    "m", "mode", "source", "pair", "K",
    "mean_A", "sd_A", "mean_a", "z", "p", "z_lt_minus2",
]


def export_results_csv(all_stats: list[dict], out_path: Path) -> None:
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(all_stats)
    print(f"Ergebnis-CSV gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Ergebnistabelle (Print)
# ---------------------------------------------------------------------------

def print_results_table(all_stats: list[dict]) -> None:
    """Gibt die Haupttabelle in der angeforderten Form aus."""
    print()
    print("=" * 110)
    print("ERGEBNISTABELLE: A_EX-Pfadasymmetrie")
    print(f"  z-Test NUR für mode=disjoint statistisch valide.")
    print(f"  mode=overlap: nur deskriptiv (Fenster korreliert).")
    print("=" * 110)

    header = (
        f"{'m':>7} | {'mode':>10} | {'src':>6} | "
        f"{'K':>5} | "
        f"{'μ_EA':>11} | {'z_EA':>7} | {'p_EA':>7} | "
        f"{'μ_EB':>11} | {'z_EB':>7} | {'p_EB':>7} | "
        f"{'μ_EC':>11} | {'z_EC':>7} | {'p_EC':>7}"
    )
    print(header)
    print("-" * 110)

    # Gruppiere nach (m, mode, source)
    groups: dict[tuple, dict[str, dict]] = {}
    for row in all_stats:
        key = (row["m"], row["mode"], row["source"])
        if key not in groups:
            groups[key] = {}
        groups[key][row["pair"]] = row

    for (m, mode, src), pairs in sorted(groups.items()):
        ea = pairs.get("EA", {})
        eb = pairs.get("EB", {})
        ec = pairs.get("EC", {})

        def fmt_z(v):
            if isinstance(v, float) and np.isnan(v):
                return "     nan"
            return f"{v:+7.3f}"

        def fmt_p(v):
            if isinstance(v, float) and np.isnan(v):
                return "     nan"
            return f"{v:7.5f}"

        def fmt_mu(v):
            if isinstance(v, float) and np.isnan(v):
                return "        nan"
            return f"{v:+11.6f}"

        z_ea = ea.get("z", float("nan"))
        z_eb = eb.get("z", float("nan"))
        z_ec = ec.get("z", float("nan"))

        flag = ""
        if mode == "disjoint":
            sig_count = sum(1 for z in [z_ea, z_eb, z_ec]
                            if not np.isnan(z) and z < -2)
            if sig_count > 0:
                flag = f"  ← {sig_count}× z<−2"

        line = (
            f"{m:>7} | {mode:>10} | {src:>6} | "
            f"{ea.get('K', '?'):>5} | "
            f"{fmt_mu(ea.get('mean_a', float('nan')))} | "
            f"{fmt_z(z_ea)} | {fmt_p(ea.get('p', float('nan')))} | "
            f"{fmt_mu(eb.get('mean_a', float('nan')))} | "
            f"{fmt_z(z_eb)} | {fmt_p(eb.get('p', float('nan')))} | "
            f"{fmt_mu(ec.get('mean_a', float('nan')))} | "
            f"{fmt_z(z_ec)} | {fmt_p(ec.get('p', float('nan')))}"
            f"{flag}"
        )
        print(line)

    print("=" * 110)


# ---------------------------------------------------------------------------
# Kernbefund (Print)
# ---------------------------------------------------------------------------

def print_core_finding(all_stats: list[dict]) -> None:
    """Beantwortet die Kernfrage explizit."""
    print()
    print("=" * 70)
    print("KERNFRAGE: Sind z_EA, z_EB, z_EC < −2 (disjunkte Fenster)?")
    print("=" * 70)

    eabc_disjoint = [
        r for r in all_stats
        if r["mode"] == "disjoint" and r["source"] == "eabc"
    ]

    pair_summary: dict[str, list[float]] = {"EA": [], "EB": [], "EC": []}
    for row in eabc_disjoint:
        pair = row["pair"]
        if pair in pair_summary and not np.isnan(row["z"]):
            pair_summary[pair].append(row["z"])

    all_confirmed = True
    any_confirmed = False
    for pair in ["EA", "EB", "EC"]:
        z_vals = pair_summary[pair]
        below = [z for z in z_vals if z < -2]
        confirmed = len(below) == len(z_vals) and len(z_vals) > 0
        partial = 0 < len(below) < len(z_vals)
        status = "BESTÄTIGT" if confirmed else ("TEILWEISE" if partial else "NEIN")
        if confirmed:
            any_confirmed = True
        else:
            all_confirmed = False
        print(f"  A_{pair}: {status}   z-Werte = {[round(z,3) for z in z_vals]}")

    print()
    if all_confirmed:
        print("  → POSITIVER BEFUND: Alle drei E-Paare zeigen z < −2")
        print("    über alle disjunkten Fenstergrößen.")
        print("    Dies ist der erste belastbare Befund einer")
        print("    EABC-Pfadasymmetrie zweiter Ordnung.")
    elif any_confirmed:
        confirmed_pairs = [
            pair for pair in ["EA", "EB", "EC"]
            if pair_summary[pair] and all(z < -2 for z in pair_summary[pair])
        ]
        print(f"  → TEILBEFUND: {confirmed_pairs} bestätigt, aber nicht alle drei Paare.")
    else:
        print("  → NEGATIVER BEFUND: Kein Paar zeigt z < −2 konsistent.")
        print("    H₀: E[A_EX] = 0 kann nicht verworfen werden.")
        print("    Die beobachteten Tendenzen sind mit Zufallsrauschen konsistent.")

    print()
    print("  Hinweis: Für mode=overlap ist der z-Test NICHT valide")
    print("  (Fenster sind korreliert). Diese Werte sind nur deskriptiv.")
    print("=" * 70)


# ---------------------------------------------------------------------------
# Visualisierung
# ---------------------------------------------------------------------------

def plot_results(
    rows_eabc: dict[tuple[int, str], list[dict]],
    rows_rand: dict[tuple[int, str], list[dict]],
    all_stats: list[dict],
    out_path: Path,
) -> None:
    """
    4 Subplots:
        A_EA, A_EB, A_EC, Gesamtnorm (Mittel über alle drei E-Paare)

    Pro Subplot:
        EABC blau (disjoint durchgezogen, overlap gestrichelt)
        Zufall grau (disjoint durchgezogen, overlap gestrichelt)
    """
    fig = plt.figure(figsize=(16, 14))
    fig.suptitle(
        r"EABC-Pfadasymmetrie 2. Ordnung: $a_{EX} = A_{EX} / \binom{m}{2}$ pro Fenster"
        "\nBlau = EABC-Primstrom, Grau = Zufall (gleiche Marginalen)"
        "\nDurchgezogen = disjoint (statistisch valide), Gestrichelt = overlap (deskriptiv)",
        fontsize=11,
        fontweight="bold",
    )

    gs = gridspec.GridSpec(2, 2, hspace=0.45, wspace=0.35)
    ax_ea = fig.add_subplot(gs[0, 0])
    ax_eb = fig.add_subplot(gs[0, 1])
    ax_ec = fig.add_subplot(gs[1, 0])
    ax_tot = fig.add_subplot(gs[1, 1])

    pk_to_ax = {
        "A_EA": (ax_ea, r"$a_{EA} = A_{EA}/\binom{m}{2}$"),
        "A_EB": (ax_eb, r"$a_{EB} = A_{EB}/\binom{m}{2}$"),
        "A_EC": (ax_ec, r"$a_{EC} = A_{EC}/\binom{m}{2}$"),
    }

    # Farb-/Stil-Mapping
    style_map = {
        ("eabc", "disjoint"):  {"lw": 1.4, "alpha": 0.80, "ls": "-"},
        ("eabc", "overlap"):   {"lw": 0.7, "alpha": 0.40, "ls": "--"},
        ("rand", "disjoint"):  {"lw": 1.2, "alpha": 0.70, "ls": "-"},
        ("rand", "overlap"):   {"lw": 0.6, "alpha": 0.30, "ls": "--"},
    }
    color_map_eabc = {1000: "#1f77b4", 5000: "#1a6fa8", 10000: "#155d8e"}
    color_map_rand = {1000: "#aaaaaa", 5000: "#888888", 10000: "#555555"}

    # Legende-Handles sammeln
    legend_handles: dict[str, object] = {}

    def _plot_series(ax, pk, rows_dict, color_map, src_key):
        for (m, mode), rows in sorted(rows_dict.items()):
            if not rows:
                continue
            x = np.array([r["start"] + m / 2 for r in rows])
            y = np.array([r[pk + "_norm"] for r in rows])
            st = style_map[(src_key, mode)]
            col = color_map[m]
            lbl_key = f"{src_key}_m{m}_{mode}"
            line, = ax.plot(
                x, y,
                color=col,
                linewidth=st["lw"],
                alpha=st["alpha"],
                linestyle=st["ls"],
                label=lbl_key,
            )
            legend_handles[lbl_key] = line

    for pk, (ax, title) in pk_to_ax.items():
        ax.axhline(0, color="black", linewidth=1.0, zorder=5)
        _plot_series(ax, pk, rows_eabc, color_map_eabc, "eabc")
        _plot_series(ax, pk, rows_rand, color_map_rand, "rand")
        ax.set_title(title, fontsize=10, pad=4)
        ax.set_ylabel(r"$a_{EX}$", fontsize=9)
        ax.set_xlabel("Primzahl-Index", fontsize=8)
        ax.grid(True, alpha=0.25)
        ax.tick_params(labelsize=8)

        # z-Score-Annotation für disjoint EABC
        pair = pk[2:]  # EA, EB, EC
        stats_here = [r for r in all_stats
                      if r["mode"] == "disjoint" and r["source"] == "eabc"
                      and r["pair"] == pair]
        annot_lines = []
        for s in sorted(stats_here, key=lambda x: x["m"]):
            z_str = f"{s['z']:+.2f}" if not np.isnan(s["z"]) else "n/a"
            p_str = f"{s['p']:.4f}" if not np.isnan(s["p"]) else "n/a"
            annot_lines.append(f"m={s['m']:>5}: z={z_str}, p={p_str}")
        annot = "\n".join(annot_lines)
        ax.text(
            0.02, 0.02, annot,
            transform=ax.transAxes,
            fontsize=7, verticalalignment="bottom",
            fontfamily="monospace",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
        )

    # Subplot 4: Gesamtnorm (Mittelwert aller drei E-Paare)
    ax_tot.axhline(0, color="black", linewidth=1.0, zorder=5)

    def _plot_total(ax, rows_dict, color_map, src_key):
        for (m, mode), rows in sorted(rows_dict.items()):
            if not rows:
                continue
            x = np.array([r["start"] + m / 2 for r in rows])
            y = np.array([
                (r["A_EA_norm"] + r["A_EB_norm"] + r["A_EC_norm"]) / 3.0
                for r in rows
            ])
            st = style_map[(src_key, mode)]
            col = color_map[m]
            lbl_suffix = "disjoint" if mode == "disjoint" else "overlap"
            label = f"m={m} ({lbl_suffix})"
            ax.plot(
                x, y,
                color=col,
                linewidth=st["lw"],
                alpha=st["alpha"],
                linestyle=st["ls"],
                label=label if src_key == "eabc" else None,
            )

    _plot_total(ax_tot, rows_eabc, color_map_eabc, "eabc")
    _plot_total(ax_tot, rows_rand, color_map_rand, "rand")

    ax_tot.set_title(
        r"$\bar{a}_E = (a_{EA}+a_{EB}+a_{EC})/3$ (Gesamtasymmetrie)",
        fontsize=10, pad=4,
    )
    ax_tot.set_ylabel(r"$\bar{a}_E$", fontsize=9)
    ax_tot.set_xlabel("Primzahl-Index", fontsize=8)
    ax_tot.grid(True, alpha=0.25)
    ax_tot.tick_params(labelsize=8)
    ax_tot.legend(fontsize=7.5, loc="upper right", ncol=2)

    # Gemeinsame Legende für EABC/Zufall
    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0], [0], color="#1f77b4", lw=1.5, ls="-",  label="EABC (disjoint)"),
        Line2D([0], [0], color="#1f77b4", lw=1.0, ls="--", label="EABC (overlap)", alpha=0.5),
        Line2D([0], [0], color="#888888", lw=1.3, ls="-",  label="Zufall (disjoint)"),
        Line2D([0], [0], color="#888888", lw=0.8, ls="--", label="Zufall (overlap)", alpha=0.5),
        Line2D([0], [0], color="black",   lw=1.0, ls="-",  label="Null (A=0)"),
    ]
    fig.legend(
        handles=legend_elems,
        loc="lower center",
        ncol=5,
        fontsize=8.5,
        bbox_to_anchor=(0.5, 0.00),
        framealpha=0.9,
    )
    fig.subplots_adjust(bottom=0.07)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Hauptprogramm
# ---------------------------------------------------------------------------

def main() -> None:
    out_dir = Path(__file__).resolve().parent
    N_PRIMES = 200_000
    RAND_SEED = 42

    print("=" * 70)
    print("EABC-Pfadasymmetrie-Test zweiter Ordnung")
    print("Observable: A_EX = N_EX - N_XE  für X ∈ {A, B, C}")
    print("H₀: E[A_EX] = 0   H₁: E[A_EX] < 0")
    print("=" * 70)

    # -----------------------------------------------------------------------
    # Primzahlwort erzeugen
    # -----------------------------------------------------------------------
    print(f"\nGeneriere EABC-Wort der ersten {N_PRIMES:,} Primzahlen > 3 …")
    word = build_eabc_word(N_PRIMES)
    n_word = len(word)
    print(f"  Wortlänge: {n_word:,}")

    probs = marginal_probs(word)
    for lab in LABELS:
        cnt = word.count(lab)
        print(f"  N_{lab} = {cnt:>7}  ({100*probs[lab]:.3f}%)")

    # -----------------------------------------------------------------------
    # Zufallsbaseline erzeugen (gleiche Marginalwahrscheinlichkeiten)
    # -----------------------------------------------------------------------
    print(f"\nGeneriere Zufallsbaseline (gleiche Marginalen, seed={RAND_SEED}) …")
    rng = np.random.default_rng(RAND_SEED)
    rand_word = random_word_matched_marginals(word, rng=rng)
    print(f"  Zufallswortlänge: {len(rand_word):,}")

    # -----------------------------------------------------------------------
    # Sliding-Window für beide Wörter und beide Modi
    # -----------------------------------------------------------------------
    all_stats: list[dict] = []

    # Speichere Fenster-Rows für den Plot
    rows_eabc: dict[tuple[int, str], list[dict]] = {}
    rows_rand: dict[tuple[int, str], list[dict]] = {}

    for m in WINDOW_SIZES:
        step_disjoint = m
        step_overlap = max(1, m // 10)

        for mode, step in [("disjoint", step_disjoint), ("overlap", step_overlap)]:
            print(
                f"  m={m:>6}, mode={mode:>9}, step={step:>5} "
                f"→ EABC …", end="", flush=True
            )
            rows_e = sliding_window(word, m, step, label="eabc")
            print(f" {len(rows_e):>5} Fenster | Zufall …", end="", flush=True)
            rows_r = sliding_window(rand_word, m, step, label="rand")
            print(f" {len(rows_r):>5} Fenster")

            rows_eabc[(m, mode)] = rows_e
            rows_rand[(m, mode)] = rows_r

            stats_e = compute_statistics(rows_e, m, mode, "eabc")
            stats_r = compute_statistics(rows_r, m, mode, "rand")
            all_stats.extend(stats_e)
            all_stats.extend(stats_r)

    # -----------------------------------------------------------------------
    # Ausgabe
    # -----------------------------------------------------------------------
    print_results_table(all_stats)
    print_core_finding(all_stats)

    # Zusatz: Vergleich EABC vs Zufall (Vorzeichen)
    print()
    print("=" * 70)
    print("VERGLEICH EABC vs. ZUFALL (disjoint, normierter Mittelwert mean_a)")
    print("=" * 70)
    print(f"{'m':>7} | {'Paar':>5} | {'EABC mean_a':>13} | {'Rand mean_a':>13} | {'Diff':>13}")
    print("-" * 65)
    for m in WINDOW_SIZES:
        for pair in ["EA", "EB", "EC"]:
            e_row = next(
                (r for r in all_stats
                 if r["m"] == m and r["mode"] == "disjoint"
                 and r["source"] == "eabc" and r["pair"] == pair),
                None,
            )
            r_row = next(
                (r for r in all_stats
                 if r["m"] == m and r["mode"] == "disjoint"
                 and r["source"] == "rand" and r["pair"] == pair),
                None,
            )
            if e_row and r_row:
                diff = e_row["mean_a"] - r_row["mean_a"]
                print(
                    f"{m:>7} | A_{pair} | {e_row['mean_a']:>+13.8f} | "
                    f"{r_row['mean_a']:>+13.8f} | {diff:>+13.8f}"
                )
    print("=" * 70)

    # CSV
    csv_path = out_dir / "eabc_asymmetry_results.csv"
    export_results_csv(all_stats, csv_path)

    # Plot
    print("\nErstelle Visualisierung …")
    plot_path = out_dir / "eabc_asymmetry_plot.png"
    plot_results(rows_eabc, rows_rand, all_stats, plot_path)

    print("\nFertig.")


if __name__ == "__main__":
    main()
