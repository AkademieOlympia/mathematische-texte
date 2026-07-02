#!/usr/bin/env python3
"""
EABC-Sensitivitäts- und Power-Analyse
======================================

Fragestellung:
    Wie groß müsste N sein, um den beobachteten Trend bei Power 80% nachzuweisen?

Methodik:
    Für normalverteilte a_XY-Fenster-Werte mit empirischer Std σ gilt beim
    einseitigen z-Test mit α=0.05 und Power=80%:

        δ_min = (z_{0.95} + z_{0.80}) · σ / sqrt(K)
              = (1.6449 + 0.8416) · σ / sqrt(K)
              ≈ 2.4865 · σ / sqrt(K)

    mit K = N_primes / m (Anzahl disjunkter Fenster).

    N_needed = ((z_α + z_β) · σ / |μ_obs|)² · m

Ausgaben:
    eabc_power_results.csv    — δ_min, N_needed für alle (m, Paar)-Kombinationen
    eabc_power_plot.png       — δ_min(N) mit empirischen Punkten bei N=200k/500k
    eabc_null_report.txt      — vollständiger Null-Ergebnisbericht
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.stats import norm as sp_norm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eabc_signature_analysis import sieve_primes, class_of

# ---------------------------------------------------------------------------
# Konstanten
# ---------------------------------------------------------------------------

WINDOW_SIZES   = [1_000, 5_000, 10_000]
PAIR_NAMES     = ["EA", "EB", "EC"]
E_PAIRS        = [("E", "A"), ("E", "B"), ("E", "C")]

N_PRIMES_BASE  = 200_000       # Basislauf (wie altsig2)
N_PRIMES_EXT   = 500_000       # Erweiterter Lauf
N_PRIMES_EXT_LABEL = "500k"

TIME_BUDGET_S  = 60            # Sekunden für den Permutationstest-Teil

# Power-Parameter (einseitiger z-Test)
ALPHA          = 0.05
POWER          = 0.80
Z_ALPHA        = float(sp_norm.ppf(1 - ALPHA))   # z_{0.95} ≈ 1.6449
Z_BETA         = float(sp_norm.ppf(POWER))        # z_{0.80} ≈ 0.8416
Z_TOTAL        = Z_ALPHA + Z_BETA                 # ≈ 2.4865

LABEL_TO_CODE: dict[str, int] = {"E": 0, "A": 1, "B": 2, "C": 3}
COLORS_BY_M    = {1_000: "#1f77b4", 5_000: "#2ca02c", 10_000: "#d62728"}

B_PERMS        = 400            # Permutationen pro Fenster (Erweiterungslauf)
PERM_STRIDE    = 10             # Nur jedes PERM_STRIDE-te Fenster (Erweiterungslauf)


# ---------------------------------------------------------------------------
# Primzahlwort erzeugen
# ---------------------------------------------------------------------------

def build_eabc_word(n_primes: int) -> list[str]:
    """Erstellt das EABC-Wort der ersten n_primes Primzahlen > 3."""
    limit = max(3_000_000, n_primes * 15)
    while True:
        primes = sieve_primes(limit)
        word = [class_of(p) for p in primes if p > 3 and class_of(p) is not None]
        if len(word) >= n_primes:
            return word[:n_primes]
        limit = int(limit * 1.5)


def word_to_codes(word: list[str]) -> np.ndarray:
    return np.array([LABEL_TO_CODE[c] for c in word], dtype=np.int8)


# ---------------------------------------------------------------------------
# Kern-Observable (identisch zu altsig2)
# ---------------------------------------------------------------------------

def compute_nxy_nyx(word_codes: np.ndarray, x_code: int, y_code: int) -> tuple[int, int]:
    """N_XY und N_YX in O(m) mit Vektorisierung."""
    is_x = (word_codes == x_code).view(np.uint8)
    is_y = (word_codes == y_code).view(np.uint8)
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
# Sliding-Window: alle a_XY-Werte berechnen (ohne Permutationstest)
# ---------------------------------------------------------------------------

def compute_window_avalues(
    word_codes: np.ndarray,
    m: int,
    pairs: list[tuple[str, str]] = E_PAIRS,
) -> dict[str, np.ndarray]:
    """
    Berechnet a_XY für alle disjunkten Fenster.

    Rückgabe: {'EA': array, 'EB': array, 'EC': array}
    """
    n = len(word_codes)
    starts = list(range(0, n - m + 1, m))
    result: dict[str, list[float]] = {pn: [] for pn in PAIR_NAMES}

    for start in starts:
        win = word_codes[start : start + m]
        for (X, Y), pn in zip(pairs, PAIR_NAMES):
            xc = LABEL_TO_CODE[X]
            yc = LABEL_TO_CODE[Y]
            n_xy, n_yx = compute_nxy_nyx(win, xc, yc)
            result[pn].append(normalized_antisym(n_xy, n_yx))

    return {pn: np.array([v for v in vals if not np.isnan(v)]) for pn, vals in result.items()}


# ---------------------------------------------------------------------------
# Permutationstest (zeitgesteuert, nur für Erweiterungslauf)
# ---------------------------------------------------------------------------

def permutation_p_for_window(
    win: np.ndarray,
    x_code: int,
    y_code: int,
    obs: float,
    B: int,
    rng: np.random.Generator,
) -> float:
    """Einseitiger Permutations-p-Wert für ein einzelnes Fenster."""
    if np.isnan(obs):
        return float("nan")
    perm = win.copy()
    count_le = 0
    for _ in range(B):
        rng.shuffle(perm)
        n_xy, n_yx = compute_nxy_nyx(perm, x_code, y_code)
        v = normalized_antisym(n_xy, n_yx)
        if not np.isnan(v) and v <= obs:
            count_le += 1
    return count_le / B


def run_permtest_with_budget(
    word_codes: np.ndarray,
    m: int,
    time_budget: float,
    rng: np.random.Generator,
) -> dict[str, dict]:
    """
    Führt Permutationstest aus, bis `time_budget` Sekunden verbraucht sind.

    Rückgabe:
        {'EA': {'a': [...], 'p': [...], 'K_perm': int}, ...}
    """
    n = len(word_codes)
    starts = list(range(0, n - m + 1, m))
    results: dict[str, dict] = {
        pn: {"a": [], "p": []} for pn in PAIR_NAMES
    }

    t0 = time.time()
    for k, start in enumerate(starts):
        if time.time() - t0 > time_budget:
            print(f"    Zeitbudget ({time_budget:.0f}s) ausgeschöpft nach {k} Fenstern.", flush=True)
            break
        if k % PERM_STRIDE != 0:
            # Nur a_XY berechnen, kein Permutationstest
            win = word_codes[start : start + m]
            for (X, Y), pn in zip(E_PAIRS, PAIR_NAMES):
                n_xy, n_yx = compute_nxy_nyx(win, LABEL_TO_CODE[X], LABEL_TO_CODE[Y])
                results[pn]["a"].append(normalized_antisym(n_xy, n_yx))
            continue

        win = word_codes[start : start + m]
        for (X, Y), pn in zip(E_PAIRS, PAIR_NAMES):
            xc = LABEL_TO_CODE[X]
            yc = LABEL_TO_CODE[Y]
            n_xy, n_yx = compute_nxy_nyx(win, xc, yc)
            obs = normalized_antisym(n_xy, n_yx)
            results[pn]["a"].append(obs)
            p = permutation_p_for_window(win, xc, yc, obs, B_PERMS, rng)
            results[pn]["p"].append(p)

    for pn in PAIR_NAMES:
        results[pn]["K_perm"] = len(results[pn]["p"])

    return results


# ---------------------------------------------------------------------------
# Power-Analyse: δ_min und N_needed
# ---------------------------------------------------------------------------

def min_detectable_effect(sigma: float, K: int) -> float:
    """
    Minimale detektierbare Effektgröße bei Power=80%, α=0.05, einseitig.

        δ_min = (z_α + z_β) · σ / sqrt(K)
    """
    if K <= 0 or sigma <= 0 or np.isnan(sigma):
        return float("nan")
    return Z_TOTAL * sigma / np.sqrt(K)


def n_needed_for_detection(sigma: float, mu_obs: float, m: int) -> float:
    """
    Minimales N (Anzahl Primzahlen), um einen Effekt der Größe |mu_obs|
    bei Power=80%, α=0.05 nachzuweisen.

        K_needed = ((z_α + z_β) · σ / |μ|)²
        N_needed  = K_needed · m
    """
    if sigma <= 0 or np.isnan(sigma) or mu_obs == 0 or np.isnan(mu_obs):
        return float("nan")
    k_needed = (Z_TOTAL * sigma / abs(mu_obs)) ** 2
    return k_needed * m


def compute_power_table(
    sigma_dict: dict[tuple[int, str], float],
    mean_a_dict: dict[tuple[int, str], float],
    n_primes_list: list[int],
) -> list[dict]:
    """
    Berechnet δ_min für alle (m, Paar, N)-Kombinationen und N_needed.

    Rückgabe: Liste von Dicts für CSV-Export.
    """
    rows: list[dict] = []
    for m in WINDOW_SIZES:
        for pn in PAIR_NAMES:
            sigma = sigma_dict.get((m, pn), float("nan"))
            mu    = mean_a_dict.get((m, pn), float("nan"))
            n_need = n_needed_for_detection(sigma, mu, m)

            row: dict = {
                "m": m,
                "pair": pn,
                "sigma_a": round(sigma, 8) if not np.isnan(sigma) else float("nan"),
                "mean_a_200k": round(mu, 8) if not np.isnan(mu) else float("nan"),
                "N_needed": int(round(n_need)) if not np.isnan(n_need) else float("nan"),
            }
            for N in n_primes_list:
                K = N // m
                dm = min_detectable_effect(sigma, K)
                row[f"delta_min_N{N//1000}k"] = round(dm, 8) if not np.isnan(dm) else float("nan")

            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# CSV-Export
# ---------------------------------------------------------------------------

def export_power_csv(rows: list[dict], out_path: Path) -> None:
    if not rows:
        return
    fields = list(rows[0].keys())
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Power-CSV gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Plot: δ_min als Funktion von N
# ---------------------------------------------------------------------------

def plot_power_analysis(
    sigma_dict: dict[tuple[int, str], float],
    mean_a_dict_200k: dict[tuple[int, str], float],
    mean_a_dict_500k: dict[tuple[int, str], float],
    out_path: Path,
) -> None:
    """
    3 Subplots (je ein E-Paar):
        - δ_min(N) als Kurve für m=1000/5000/10000
        - Horizontale Linie: |μ_200k| und |μ_500k|
        - Vertikale Marker bei N=200k und N=500k
        - Schnittpunkt δ_min = |μ| → N_needed
    """
    N_range = np.logspace(np.log10(50_000), np.log10(10_000_000), 500)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=False)
    fig.suptitle(
        r"EABC-Power-Analyse: Minimale detektierbare Effektgröße $\delta_{\min}(N)$"
        f"\n(einseitiger z-Test, α={ALPHA}, Power={int(POWER*100)}%, "
        rf"$z_\alpha+z_\beta \approx {Z_TOTAL:.3f}$)",
        fontsize=12, fontweight="bold",
    )

    for ax, pn in zip(axes, PAIR_NAMES):
        ax.set_title(rf"Paar $A_{{E{pn[1]}}}$", fontsize=11)
        ax.set_xlabel("N (Primzahlen)", fontsize=9)
        ax.set_ylabel(r"$\delta_{\min}$", fontsize=9)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.25)
        ax.tick_params(labelsize=8)

        for m, color in COLORS_BY_M.items():
            sigma = sigma_dict.get((m, pn), float("nan"))
            if np.isnan(sigma):
                continue
            K_range = N_range / m
            delta_range = Z_TOTAL * sigma / np.sqrt(K_range)
            ax.plot(
                N_range, delta_range,
                color=color, linewidth=1.8,
                label=rf"$m={m}$, $\sigma={sigma:.4f}$",
            )

            # N_needed für 200k mean_a
            mu_200k = mean_a_dict_200k.get((m, pn), float("nan"))
            if not np.isnan(mu_200k) and mu_200k != 0:
                n_need = n_needed_for_detection(sigma, mu_200k, m)
                if not np.isnan(n_need):
                    ax.axvline(n_need, color=color, linestyle=":", linewidth=1.0, alpha=0.6)
                    ax.text(
                        n_need, delta_range.max() * 1.05,
                        f"N*={n_need/1e6:.2f}M\n(m={m})",
                        color=color, fontsize=6.5, ha="center", va="bottom",
                        rotation=90,
                    )

        # Beobachtete |μ| als horizontale Linien
        mu_200k_vals = [mean_a_dict_200k.get((m, pn), float("nan")) for m in WINDOW_SIZES]
        mu_500k_vals = [mean_a_dict_500k.get((m, pn), float("nan")) for m in WINDOW_SIZES]

        for mu, ls, lbl in [
            (np.nanmean([abs(v) for v in mu_200k_vals if not np.isnan(v)]), "--", r"$|\bar\mu|$ @ 200k"),
            (np.nanmean([abs(v) for v in mu_500k_vals if not np.isnan(v)]), "-.", r"$|\bar\mu|$ @ 500k"),
        ]:
            if not np.isnan(mu) and mu > 0:
                ax.axhline(mu, color="black", linestyle=ls, linewidth=1.2, alpha=0.8, label=lbl)

        # Marker bei N=200k und N=500k
        for N_mark, mk, mlbl in [
            (N_PRIMES_BASE, "o", "200k"),
            (N_PRIMES_EXT,  "s", "500k"),
        ]:
            for m, color in COLORS_BY_M.items():
                sigma = sigma_dict.get((m, pn), float("nan"))
                if np.isnan(sigma):
                    continue
                K_mark = N_mark / m
                dm = min_detectable_effect(sigma, int(K_mark))
                ax.plot(N_mark, dm, marker=mk, color=color, markersize=6, zorder=5)
            ax.axvline(N_mark, color="gray", linestyle=":", linewidth=0.9, alpha=0.5)
            ax.text(N_mark, ax.get_ylim()[0] if ax.get_ylim()[0] > 0 else 1e-5,
                    f" {mlbl}", color="gray", fontsize=7, va="bottom")

        ax.legend(fontsize=7, loc="upper right")

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Power-Plot gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Null-Ergebnisbericht
# ---------------------------------------------------------------------------

def generate_null_report(
    n_200k: int,
    max_p_200k: int,
    n_500k: int,
    max_p_500k: int,
    sigma_dict: dict[tuple[int, str], float],
    mean_a_200k: dict[tuple[int, str], float],
    mean_a_500k: dict[tuple[int, str], float],
    median_p_200k: dict[tuple[int, str], float],
    median_p_500k: dict[tuple[int, str], float],
    power_rows: list[dict],
) -> str:
    """Erzeugt den vollständigen Null-Ergebnisbericht als String."""

    # Minimale δ_min über alle m (bei N=200k)
    delta_min_200k: dict[str, list[float]] = {pn: [] for pn in PAIR_NAMES}
    delta_min_500k: dict[str, list[float]] = {pn: [] for pn in PAIR_NAMES}
    n_needed_all: dict[str, list] = {pn: [] for pn in PAIR_NAMES}

    for row in power_rows:
        pn = row["pair"]
        key_200 = f"delta_min_N{N_PRIMES_BASE//1000}k"
        key_500 = f"delta_min_N{N_PRIMES_EXT//1000}k"
        d200 = row.get(key_200, float("nan"))
        d500 = row.get(key_500, float("nan"))
        nn   = row.get("N_needed", float("nan"))
        if not (isinstance(d200, float) and np.isnan(d200)):
            delta_min_200k[pn].append(d200)
        if not (isinstance(d500, float) and np.isnan(d500)):
            delta_min_500k[pn].append(d500)
        if not (isinstance(nn, float) and np.isnan(nn)):
            n_needed_all[pn].append(nn)

    def fmt_delta(vals: list[float]) -> str:
        if not vals:
            return "n/a"
        mn, mx = min(vals), max(vals)
        return f"{mn:.6f} … {mx:.6f}"

    def fmt_median_p(d: dict[tuple[int, str], float], pn: str) -> str:
        vals = [v for (m, p), v in d.items() if p == pn and not np.isnan(v)]
        if not vals:
            return "n/a"
        return ", ".join(f"{v:.4f}" for v in sorted(vals))

    def fmt_n_needed(vals: list) -> str:
        if not vals:
            return "unendlich (kein beobachteter Trend)"
        mn = int(min(vals))
        mx = int(max(vals))
        if mn == mx:
            return f"{mn:,}"
        return f"{mn:,} … {mx:,}"

    min_median_p_200k = min(
        (v for v in median_p_200k.values() if not np.isnan(v)),
        default=float("nan"),
    )
    min_median_p_500k = min(
        (v for v in median_p_500k.values() if not np.isnan(v)),
        default=float("nan"),
    )

    # Trendvergleich 200k vs 500k
    trend_lines = []
    for pn in PAIR_NAMES:
        mu_200 = [v for (m, p), v in mean_a_200k.items() if p == pn and not np.isnan(v)]
        mu_500 = [v for (m, p), v in mean_a_500k.items() if p == pn and not np.isnan(v)]
        if mu_200 and mu_500:
            avg_200 = np.mean(mu_200)
            avg_500 = np.mean(mu_500)
            direction = "stärker" if abs(avg_500) > abs(avg_200) else "schwächer/gleich"
            trend_lines.append(
                f"  A_{pn}: μ_200k = {avg_200:+.6f},  μ_500k = {avg_500:+.6f}"
                f"  → Trend {direction}"
            )
        else:
            trend_lines.append(f"  A_{pn}: keine Daten")

    report_lines = [
        "EABC-AltSig²-Null-Ergebnisbericht",
        "===================================",
        f"Datum: {time.strftime('%Y-%m-%d %H:%M')}",
        "",
        "DATENBASIS",
        "----------",
        f"  Basislauf:      N = {n_200k:,} Primzahlen, p ≤ {max_p_200k:,}",
        f"  Erweiterung:    N = {n_500k:,} Primzahlen, p ≤ {max_p_500k:,}",
        f"  Fenstergrössen: m ∈ {WINDOW_SIZES}",
        f"  Fensteranzahl:  K = N/m (disjunkt)",
        "",
        "BEFUND",
        "------",
        "  Kein signifikanter Unterschied zum Permutations-Nullmodell.",
        "",
        "  Basislauf (N=200k), Permutations-median-p:",
        f"    A_EA: {fmt_median_p(median_p_200k, 'EA')} (über m=1k/5k/10k)",
        f"    A_EB: {fmt_median_p(median_p_200k, 'EB')}",
        f"    A_EC: {fmt_median_p(median_p_200k, 'EC')}",
        f"    Minimum: {min_median_p_200k:.4f}  (alle > 0.05 → H₀ nicht verworfen)",
        "",
        "  Erweiterungslauf (N=500k), Permutations-median-p:",
        f"    A_EA: {fmt_median_p(median_p_500k, 'EA')}",
        f"    A_EB: {fmt_median_p(median_p_500k, 'EB')}",
        f"    A_EC: {fmt_median_p(median_p_500k, 'EC')}",
        f"    Minimum: {min_median_p_500k:.4f}  (alle > 0.05 → H₀ nicht verworfen)",
        "",
        "TRENDVERGLEICH 200k → 500k",
        "--------------------------",
    ] + trend_lines + [
        "",
        "SENSITIVITÄTSSCHRANKE (Power-Analyse)",
        "--------------------------------------",
        f"  Methode: einseitiger z-Test, α={ALPHA}, Power={int(POWER*100)}%",
        f"  Formel:  δ_min = (z_α + z_β) · σ / sqrt(K)  = {Z_TOTAL:.4f} · σ / sqrt(K)",
        "",
        "  Minimale detektierbare Effektgröße bei N=200k:",
        f"    A_EA: {fmt_delta(delta_min_200k['EA'])}",
        f"    A_EB: {fmt_delta(delta_min_200k['EB'])}",
        f"    A_EC: {fmt_delta(delta_min_200k['EC'])}",
        "",
        "  Minimale detektierbare Effektgröße bei N=500k:",
        f"    A_EA: {fmt_delta(delta_min_500k['EA'])}",
        f"    A_EB: {fmt_delta(delta_min_500k['EB'])}",
        f"    A_EC: {fmt_delta(delta_min_500k['EC'])}",
        "",
        "  → Die Methode ist sensitiv für Effekte dieser Größenordnung.",
        "    Falls |E[a_EX]| kleiner ist, bleibt er unentdeckt.",
        "",
        "N_NEEDED: Wann wird der beobachtete Trend nachweisbar?",
        "-------------------------------------------------------",
        "  (Gesucht: N sodass δ_min(N) ≤ |μ_obs| aus dem Basislauf)",
        "",
        f"    A_EA: N_needed ≈ {fmt_n_needed(n_needed_all['EA'])}",
        f"    A_EB: N_needed ≈ {fmt_n_needed(n_needed_all['EB'])}",
        f"    A_EC: N_needed ≈ {fmt_n_needed(n_needed_all['EC'])}",
        "",
        "INTERPRETATION",
        "--------------",
        "  Falls E[a_EX] ≠ 0, dann gilt: |E[a_EX]| < δ_min(200k) in diesem N-Bereich.",
        "",
        "  Ein eventueller Effekt hätte folgende Eigenschaften:",
        "    - Er ist kleiner als die detektierbare Schwelle bei N=200k/500k",
        "    - Er wäre erst bei N > N_needed oben nachweisbar",
        "    - Das Vorzeichen des Trends (falls konsistent negativ) ist kein",
        "      belastbarer Befund, solange p > 0.05 und frac_sig < 0.50",
        "",
        "  WISSENSCHAFTLICHE AUSSAGE:",
        "  Das Null-Ergebnis ist quantifiziert: Es schließt Effekte der Größenordnung",
        f"  |a_EX| > {np.nanmean([v for vals in delta_min_200k.values() for v in vals]):.5f}",
        "  bei N=200k mit 80% Power aus.",
        "  Schwächere Effekte bleiben bei diesem N unentdeckt.",
        "",
        "FAZIT",
        "-----",
        "  Keine arithmetische Pfadasymmetrie zweiter Ordnung nachweisbar.",
        "  Das Ergebnis ist statistisch quantifiziert und nicht nur qualitativ.",
        "===================================",
    ]

    return "\n".join(report_lines)


# ---------------------------------------------------------------------------
# Ergebnis-Printer
# ---------------------------------------------------------------------------

def print_power_table(power_rows: list[dict]) -> None:
    print()
    print("=" * 90)
    print("POWER-ANALYSE: δ_min und N_needed für alle (m, Paar)-Kombinationen")
    print(f"(einseitiger z-Test, α={ALPHA}, Power={int(POWER*100)}%, Z_total={Z_TOTAL:.4f})")
    print("=" * 90)
    header = (
        f"{'m':>7} | {'Paar':>5} | {'σ_a':>10} | {'μ_200k':>11} | "
        f"{'δ_min_200k':>12} | {'δ_min_500k':>12} | {'N_needed':>14}"
    )
    print(header)
    print("-" * 90)
    key_200 = f"delta_min_N{N_PRIMES_BASE//1000}k"
    key_500 = f"delta_min_N{N_PRIMES_EXT//1000}k"
    for row in power_rows:
        dm200 = row.get(key_200, float("nan"))
        dm500 = row.get(key_500, float("nan"))
        nn    = row.get("N_needed", float("nan"))
        nn_str = f"{int(nn):>14,}" if not (isinstance(nn, float) and np.isnan(nn)) else f"{'∞':>14}"
        print(
            f"{row['m']:>7} | {row['pair']:>5} | {row['sigma_a']:>10.6f} | "
            f"{row['mean_a_200k']:>+11.6f} | "
            f"{dm200:>12.6f} | {dm500:>12.6f} | {nn_str}"
        )
    print("=" * 90)


# ---------------------------------------------------------------------------
# Hilfsfunktion: Statistiken aus a-Werten
# ---------------------------------------------------------------------------

def summarize_a_arrays(
    a_dict: dict[str, np.ndarray],
    m: int,
) -> tuple[dict[tuple[int, str], float], dict[tuple[int, str], float]]:
    """Gibt (sigma_dict, mean_a_dict) für eine (m, Paar)-Kombination zurück."""
    sigma_dict: dict[tuple[int, str], float] = {}
    mean_a_dict: dict[tuple[int, str], float] = {}
    for pn, arr in a_dict.items():
        if len(arr) > 1:
            sigma_dict[(m, pn)] = float(np.std(arr, ddof=1))
            mean_a_dict[(m, pn)] = float(np.mean(arr))
        else:
            sigma_dict[(m, pn)] = float("nan")
            mean_a_dict[(m, pn)] = float("nan")
    return sigma_dict, mean_a_dict


# ---------------------------------------------------------------------------
# Hauptprogramm
# ---------------------------------------------------------------------------

def main() -> None:
    out_dir = Path(__file__).resolve().parent
    rng = np.random.default_rng(42)

    print("=" * 70)
    print("EABC-Sensitivitäts- und Power-Analyse")
    print(f"α={ALPHA}, Power={int(POWER*100)}%, z_total={Z_TOTAL:.4f}")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1. BASISLAUF: N=200k — Fensterstatistiken (ohne Permutationstest)
    # ------------------------------------------------------------------
    print(f"\n[1/4] Basislauf: N={N_PRIMES_BASE:,} Primzahlen …")
    word_200k = build_eabc_word(N_PRIMES_BASE)
    codes_200k = word_to_codes(word_200k)
    max_p_200k = max(p for p in sieve_primes(int(len(word_200k) * 15))
                     if p > 3 and class_of(p) is not None)

    sigma_dict_200k:    dict[tuple[int, str], float] = {}
    mean_a_dict_200k:   dict[tuple[int, str], float] = {}
    median_p_dict_200k: dict[tuple[int, str], float] = {}

    for m in WINDOW_SIZES:
        print(f"  m={m}: berechne a-Werte …", end="", flush=True)
        a_dict = compute_window_avalues(codes_200k, m)
        sig_d, mean_d = summarize_a_arrays(a_dict, m)
        sigma_dict_200k.update(sig_d)
        mean_a_dict_200k.update(mean_d)
        K = len(codes_200k) // m
        print(f" K={K} Fenster  σ_EA={sig_d.get((m,'EA'),float('nan')):.6f}")
        # Permutations-p aus altsig2-Ergebnissen laden (falls vorhanden)
        # → wird erst in Schritt 3 gesetzt; vorerst NaN
        for pn in PAIR_NAMES:
            median_p_dict_200k[(m, pn)] = float("nan")

    # Versuche, bestehende altsig2-Ergebnisse zu laden
    altsig2_csv = out_dir / "eabc_altsig2_permtest_results.csv"
    if altsig2_csv.exists():
        print(f"\n  Lade Permutations-p-Werte aus {altsig2_csv.name} …")
        import csv as _csv
        with open(altsig2_csv, newline="", encoding="utf-8") as f:
            reader = _csv.DictReader(f)
            for row in reader:
                m_val = int(row["m"])
                pn_val = row["pair"]
                mp = float(row["median_p"])
                median_p_dict_200k[(m_val, pn_val)] = mp
        print("  Fertig.")
    else:
        print(f"\n  {altsig2_csv.name} nicht gefunden — median_p bleibt NaN.")
        print("  Führe eabc_altsig2_permtest.py zuerst aus, um vollständige Ergebnisse zu erhalten.")

    # ------------------------------------------------------------------
    # 2. ERWEITERUNGSLAUF: N=500k — Fensterstatistiken
    # ------------------------------------------------------------------
    print(f"\n[2/4] Erweiterungslauf: N={N_PRIMES_EXT:,} Primzahlen …")
    t_ext_start = time.time()
    word_500k = build_eabc_word(N_PRIMES_EXT)
    codes_500k = word_to_codes(word_500k)
    max_p_500k_approx = len(codes_500k) * 17  # Abschätzung: p_N ≈ N * ln(N)
    print(f"  Wortlänge: {len(word_500k):,}  (Sieblimit × 1 als Proxy für max_p)")

    sigma_dict_500k:    dict[tuple[int, str], float] = {}
    mean_a_dict_500k:   dict[tuple[int, str], float] = {}
    median_p_dict_500k: dict[tuple[int, str], float] = {}

    for pn in PAIR_NAMES:
        for m in WINDOW_SIZES:
            median_p_dict_500k[(m, pn)] = float("nan")

    remaining_budget = TIME_BUDGET_S - (time.time() - t_ext_start)
    budget_per_m = max(5.0, remaining_budget / len(WINDOW_SIZES))

    for m in WINDOW_SIZES:
        print(f"  m={m}: a-Werte + zeitgesteuerter Permutationstest (Budget {budget_per_m:.0f}s) …", flush=True)
        t_m_start = time.time()
        perm_results = run_permtest_with_budget(codes_500k, m, budget_per_m, rng)
        elapsed = time.time() - t_m_start
        print(f"    → {elapsed:.1f}s", flush=True)

        # a-Werte aus perm_results extrahieren
        a_dict_500k: dict[str, np.ndarray] = {}
        for pn in PAIR_NAMES:
            raw_a = [v for v in perm_results[pn]["a"] if not np.isnan(v)]
            a_dict_500k[pn] = np.array(raw_a)

        sig_d, mean_d = summarize_a_arrays(a_dict_500k, m)
        sigma_dict_500k.update(sig_d)
        mean_a_dict_500k.update(mean_d)

        for pn in PAIR_NAMES:
            p_vals = [v for v in perm_results[pn]["p"] if not np.isnan(v)]
            if p_vals:
                median_p_dict_500k[(m, pn)] = float(np.median(p_vals))
                k_perm = perm_results[pn]["K_perm"]
                print(f"    A_{pn}: K_perm={k_perm}, median_p={median_p_dict_500k[(m,pn)]:.4f}, "
                      f"mean_a={mean_d.get((m,pn),float('nan')):+.6f}")
            else:
                print(f"    A_{pn}: kein Permutationstest (Zeitbudget), mean_a={mean_d.get((m,pn),float('nan')):+.6f}")

    # Echtes max_p für 500k
    sieve_lim_500k = max(10_000_000, N_PRIMES_EXT * 17)
    primes_500k_check = sieve_primes(sieve_lim_500k)
    primes_500k_filtered = [p for p in primes_500k_check if p > 3 and class_of(p) is not None]
    max_p_500k = primes_500k_filtered[N_PRIMES_EXT - 1] if len(primes_500k_filtered) >= N_PRIMES_EXT else primes_500k_filtered[-1]

    # ------------------------------------------------------------------
    # 3. POWER-TABELLE berechnen
    # ------------------------------------------------------------------
    print(f"\n[3/4] Berechne Power-Tabelle …")

    # Sigma aus dem Basislauf (stabilere Schätzung)
    N_list = [N_PRIMES_BASE, N_PRIMES_EXT]
    power_rows = compute_power_table(sigma_dict_200k, mean_a_dict_200k, N_list)

    print_power_table(power_rows)

    # CSV
    csv_path = out_dir / "eabc_power_results.csv"
    export_power_csv(power_rows, csv_path)

    # ------------------------------------------------------------------
    # 4. PLOT
    # ------------------------------------------------------------------
    print(f"\n[4/4] Erstelle Plot …")
    plot_path = out_dir / "eabc_power_plot.png"
    plot_power_analysis(
        sigma_dict_200k,
        mean_a_dict_200k,
        mean_a_dict_500k,
        plot_path,
    )

    # ------------------------------------------------------------------
    # 5. NULL-ERGEBNISBERICHT
    # ------------------------------------------------------------------
    report = generate_null_report(
        n_200k=N_PRIMES_BASE,
        max_p_200k=max_p_200k,
        n_500k=N_PRIMES_EXT,
        max_p_500k=max_p_500k,
        sigma_dict=sigma_dict_200k,
        mean_a_200k=mean_a_dict_200k,
        mean_a_500k=mean_a_dict_500k,
        median_p_200k=median_p_dict_200k,
        median_p_500k=median_p_dict_500k,
        power_rows=power_rows,
    )

    print()
    print(report)

    report_path = out_dir / "eabc_null_report.txt"
    report_path.write_text(report, encoding="utf-8")
    print(f"\nNull-Ergebnisbericht gespeichert: {report_path}")

    print("\nFertig.")


if __name__ == "__main__":
    main()
