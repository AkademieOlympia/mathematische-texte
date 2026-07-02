#!/usr/bin/env python3
"""
EABC Sliding-Window Bootstrap-Test (Referee-Formulierung)
==========================================================
Observable: AltSig²_E = (N_XY - N_YX)_{X,Y}  — erste ordnungssensitive
            Observable des EABC-Primstroms.

Offene Testfrage: E[AltSig²_E] = 0 ?

Wenn z_EA, z_EB, z_EC < -2 stabil für mehrere Fenstergrößen:
→ arithmetische Pfadasymmetrie zweiter Ordnung bestätigt.

Verwendete Implementierungsstruktur (Referee-Vorgabe):
    antisymmetry_matrix(M) = M - M.T
    sliding_window_A(word, window_size, step) → rows mit A_EX etc.
"""

from __future__ import annotations

import csv
import sys
from math import sqrt
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Importiere bestehende Infrastruktur aus eabc_signature_analysis.py
# ---------------------------------------------------------------------------
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

# Alle antisymmetrischen Off-Diagonal-Paare (X < Y lexikografisch):
# Vollständig: EA, EB, EC, AB, AC, BC
ALL_PAIRS = [("E", "A"), ("E", "B"), ("E", "C"), ("A", "B"), ("A", "C"), ("B", "C")]


# ---------------------------------------------------------------------------
# Primzahlwort erzeugen
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


# ---------------------------------------------------------------------------
# Kernfunktionen nach Referee-Vorgabe
# ---------------------------------------------------------------------------

def antisymmetry_matrix(level2_dict: dict) -> dict[str, dict[str, int]]:
    """
    A_XY = N_XY - N_YX  für alle X, Y ∈ {E, A, B, C}.
    Entspricht M - M.T (als verschachtelte dict-Darstellung).
    """
    A: dict[str, dict[str, int]] = {}
    for x in LABELS:
        A[x] = {}
        for y in LABELS:
            A[x][y] = level2_dict[(x, y)] - level2_dict[(y, x)]
    return A


def sliding_window_A(
    word: list[str],
    window_size: int,
    step: int,
) -> list[dict]:
    """
    Berechnet die antisymmetrische Level-2-Observable für alle
    nicht-überlappenden Fenster.

    Rückgabe: Liste von Dicts mit:
        start, end, A_EA, A_EB, A_EC, A_AB, A_AC, A_BC
    (rohe, nicht normierte Werte)
    """
    rows = []
    n = len(word)
    for start in range(0, n - window_size + 1, step):
        w = word[start : start + window_size]
        sig = eabc_signature(w, max_level=2)
        A = antisymmetry_matrix(sig.level2)
        rows.append({
            "start": start,
            "end": start + window_size,
            "A_EA": A["E"]["A"],
            "A_EB": A["E"]["B"],
            "A_EC": A["E"]["C"],
            "A_AB": A["A"]["B"],
            "A_AC": A["A"]["C"],
            "A_BC": A["B"]["C"],
        })
    return rows


# ---------------------------------------------------------------------------
# Statistische Auswertung pro Fenstergröße
# ---------------------------------------------------------------------------

def compute_statistics(rows: list[dict], m: int) -> dict[str, dict]:
    """
    Für jedes Paar XY berechne:
        μ_EX(m)  = Mittelwert von A_EX
        z_EX(m)  = μ_EX(m) / (sd(A_EX) / sqrt(K))
        (z ist der einseitige Z-Score unter H0: E[A_EX] = 0)
    """
    pair_keys = ["A_EA", "A_EB", "A_EC", "A_AB", "A_AC", "A_BC"]
    K = len(rows)
    stats: dict[str, dict] = {}

    for pk in pair_keys:
        vals = np.array([r[pk] for r in rows], dtype=float)
        mu = float(np.mean(vals))
        sd = float(np.std(vals, ddof=1))
        se = sd / sqrt(K) if K > 1 else float("nan")
        z = mu / se if se > 0 else float("nan")
        stats[pk] = {
            "m": m,
            "K": K,
            "mu": mu,
            "sd": sd,
            "se": se,
            "z": z,
        }
    return stats


# ---------------------------------------------------------------------------
# CSV-Export: eabc_signature_bootstrap.csv
# ---------------------------------------------------------------------------

def export_window_csv(
    all_rows: list[dict],   # rows haben zusätzlich Spalte 'm'
    out_path: Path,
) -> None:
    fieldnames = ["m", "start", "end", "A_EA", "A_EB", "A_EC", "A_AB", "A_AC", "A_BC"]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"Fenster-CSV gespeichert: {out_path}")


def export_stats_csv(
    all_stats: dict[int, dict[str, dict]],
    out_path: Path,
) -> None:
    fieldnames = ["pair", "m", "K", "mu", "sd", "se", "z",
                  "z_below_minus2", "sign_neg"]
    rows = []
    for m in WINDOW_SIZES:
        for pk, s in all_stats[m].items():
            rows.append({
                "pair": pk,
                "m": m,
                "K": s["K"],
                "mu": round(s["mu"], 6),
                "sd": round(s["sd"], 6),
                "se": round(s["se"], 6),
                "z": round(s["z"], 4),
                "z_below_minus2": s["z"] < -2,
                "sign_neg": s["mu"] < 0,
            })
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Statistik-CSV gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Visualisierung
# ---------------------------------------------------------------------------

COLORS = {1000: "#e74c3c", 5000: "#2980b9", 10000: "#27ae60"}
E_PAIRS = ["A_EA", "A_EB", "A_EC"]
E_PAIR_TITLES = {
    "A_EA": r"$A_{EA} = N_{EA} - N_{AE}$",
    "A_EB": r"$A_{EB} = N_{EB} - N_{BE}$",
    "A_EC": r"$A_{EC} = N_{EC} - N_{CE}$",
}


def plot_results(
    all_rows_by_m: dict[int, list[dict]],
    all_stats: dict[int, dict[str, dict]],
    out_path: Path,
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=False)
    fig.suptitle(
        r"EABC: $\mathrm{AltSig}^2_E = A_{EX}/m$ pro Fenster"
        "\n(Referee-Formulierung: Offene Testfrage E[AltSig²_E] = 0 ?)",
        fontsize=12,
        fontweight="bold",
    )

    for ax_idx, pk in enumerate(E_PAIRS):
        ax = axes[ax_idx]

        # Globaler Mittelwert (gewichtet über alle m)
        all_vals = []
        for m in WINDOW_SIZES:
            all_vals.extend([r[pk] / m for r in all_rows_by_m[m]])
        global_mean = float(np.mean(all_vals))

        ax.axhline(0, color="black", linewidth=1.0, linestyle="-", zorder=5)
        ax.axhline(
            global_mean,
            color="gray",
            linewidth=1.2,
            linestyle="--",
            label=f"Gesamtmittel = {global_mean:+.5f}",
            zorder=4,
        )

        for m in WINDOW_SIZES:
            rows = all_rows_by_m[m]
            s = all_stats[m][pk]
            x_pos = np.array([r["start"] + m / 2 for r in rows])
            y_vals = np.array([r[pk] / m for r in rows])

            label = (
                f"m={m}:  μ={s['mu']/m:+.4f},  "
                f"z={s['z']:+.3f}"
            )
            ax.plot(
                x_pos,
                y_vals,
                color=COLORS[m],
                linewidth=0.9 if m == 1000 else 1.3,
                alpha=0.65 if m == 1000 else 0.85,
                label=label,
                zorder=3,
            )

        ax.set_title(E_PAIR_TITLES[pk], fontsize=11, pad=5)
        ax.set_ylabel(r"$A_{EX}/m$", fontsize=9)
        ax.set_xlabel("Primzahl-Index", fontsize=9)
        ax.legend(fontsize=8.5, loc="upper right")
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=8)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Hauptprogramm
# ---------------------------------------------------------------------------

def main() -> None:
    out_dir = Path(__file__).resolve().parent
    N_PRIMES = 200_000

    print("=" * 65)
    print("EABC Sliding-Window Bootstrap-Test (Referee-Formulierung)")
    print("AltSig²_E = (N_XY - N_YX)_{X,Y}   Testfrage: E[AltSig²_E] = 0?")
    print("=" * 65)

    print(f"\nGeneriere EABC-Wort der ersten {N_PRIMES:,} Primzahlen > 3 …")
    word = build_eabc_word(N_PRIMES)
    print(f"  Wortlänge: {len(word):,}")

    for lab in LABELS:
        cnt = word.count(lab)
        print(f"  N_{lab} = {cnt:>7}  ({100*cnt/len(word):.3f}%)")

    # -----------------------------------------------------------------------
    # Sliding-Window-Analyse
    # -----------------------------------------------------------------------
    all_rows_by_m: dict[int, list[dict]] = {}
    all_stats: dict[int, dict[str, dict]] = {}
    all_csv_rows: list[dict] = []

    for m in WINDOW_SIZES:
        print(f"\nFenstergröße m={m} (step={m}, nicht-überlappend) …", flush=True)
        rows = sliding_window_A(word, window_size=m, step=m)
        print(f"  → {len(rows)} Fenster berechnet")

        all_rows_by_m[m] = rows
        all_stats[m] = compute_statistics(rows, m)

        for r in rows:
            all_csv_rows.append({"m": m, **r})

    # -----------------------------------------------------------------------
    # Ergebnistabelle
    # -----------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("ERGEBNISTABELLE  (rohe A_XY-Werte, normiert durch m für μ)")
    print(f"{'Paar':<8} {'m':>6} {'K':>5} {'μ/m':>11} {'sd':>10} "
          f"{'z':>8}  {'z < -2?':>8}  {'μ < 0?':>7}")
    print("-" * 75)

    for m in WINDOW_SIZES:
        for pk in ["A_EA", "A_EB", "A_EC", "A_AB", "A_AC", "A_BC"]:
            s = all_stats[m][pk]
            z_flag = "JA ***" if s["z"] < -2 else ("" if s["z"] > 2 else "")
            neg_flag = "JA" if s["mu"] < 0 else "NEIN"
            print(f"{pk:<8} {m:>6} {s['K']:>5} {s['mu']/m:>+11.5f} "
                  f"{s['sd']:>10.2f} {s['z']:>+8.3f}  {z_flag:<10} {neg_flag:>7}")
        print()

    # -----------------------------------------------------------------------
    # Referee-Kernbefund: z < -2 stabil?
    # -----------------------------------------------------------------------
    print("=" * 65)
    print("REFEREE-KERNFRAGE: z_EA, z_EB, z_EC < -2 für mehrere m?")
    print("=" * 65)

    for pk in ["A_EA", "A_EB", "A_EC"]:
        z_vals = {m: all_stats[m][pk]["z"] for m in WINDOW_SIZES}
        count_below = sum(1 for z in z_vals.values() if z < -2)
        print(f"\n  {pk}:")
        for m, z in z_vals.items():
            flag = " ← z < -2" if z < -2 else ""
            print(f"    m={m:>6}: z = {z:+.3f}{flag}")
        if count_below == len(WINDOW_SIZES):
            print(f"  → BESTÄTIGT: z < -2 für alle {len(WINDOW_SIZES)} Fenstergrößen")
        elif count_below > 0:
            print(f"  → TEILWEISE: z < -2 für {count_below}/{len(WINDOW_SIZES)} Fenstergrößen")
        else:
            print(f"  → NICHT BESTÄTIGT: kein z < -2")

    # -----------------------------------------------------------------------
    # Gesamtfazit
    # -----------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("GESAMTFAZIT")
    print("=" * 65)

    e_pairs_confirmed = []
    for pk in ["A_EA", "A_EB", "A_EC"]:
        all_below = all(all_stats[m][pk]["z"] < -2 for m in WINDOW_SIZES)
        if all_below:
            e_pairs_confirmed.append(pk)

    if len(e_pairs_confirmed) == 3:
        print("""
  POSITIV: Arithmetische Pfadasymmetrie 2. Ordnung bestätigt.
  z_EA, z_EB, z_EC < -2 stabil über alle Fenstergrößen.
  Dies wäre der erste statistisch gesicherte Level-2-Befund
  des EABC-Programms.
""")
    elif len(e_pairs_confirmed) > 0:
        print(f"""
  TEILBEFUND: {', '.join(e_pairs_confirmed)} zeigen z < -2,
  aber nicht alle drei E-Paare konsistent bestätigt.
""")
    else:
        print("""
  NEGATIV: Kein Paar zeigt z < -2 stabil für mehrere Fenstergrößen.
  H0: E[AltSig²_E] = 0 kann nicht verworfen werden.
  Die beobachteten Vorzeichen-Tendenzen sind mit Zufallsrauschen
  konsistent.
""")

    # -----------------------------------------------------------------------
    # Ausgabe
    # -----------------------------------------------------------------------
    csv_windows = out_dir / "eabc_signature_bootstrap.csv"
    csv_stats = out_dir / "eabc_bootstrap_results.csv"
    plot_path = out_dir / "eabc_bootstrap_plot.png"

    export_window_csv(all_csv_rows, csv_windows)
    export_stats_csv(all_stats, csv_stats)

    print("\nErstelle Visualisierung …")
    plot_results(all_rows_by_m, all_stats, plot_path)

    print("\nFertig.")


if __name__ == "__main__":
    main()
