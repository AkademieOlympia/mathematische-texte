#!/usr/bin/env python3
"""
EABC-Pfadsignatur (Chen-Signaturen auf diskreten Worten)

Klassifikation: Primzahlen p > 3 werden nach p mod 12 eingeteilt:
  E: p ≡  1 (mod 12)
  A: p ≡  5 (mod 12)
  B: p ≡  7 (mod 12)
  C: p ≡ 11 (mod 12)

Pfadsignatur bis Stufe 3:
  Stufe 1: N_X         — Häufigkeiten (4 Werte)
  Stufe 2: N_XY        — geordnete Teilwörter: #{(i,j): i<j, w[i]=X, w[j]=Y}
  Stufe 3: N_XYZ       — geordnete Teilwörter: #{(i,j,k): i<j<k, ...}

Zentralfrage: Weicht die Level-2-Signatur signifikant von der
Unabhängigkeits-Nullhypothese ab (N_XY / (N_X * N_Y) ≠ 1)?
"""

from __future__ import annotations

import csv
import random
import sys
from math import isqrt
from pathlib import Path
from typing import NamedTuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# ---------------------------------------------------------------------------
# Konstanten
# ---------------------------------------------------------------------------

LABELS = ("E", "A", "B", "C")
LABEL_SET = set(LABELS)
IDX = {lab: i for i, lab in enumerate(LABELS)}

RESIDUE_MAP = {1: "E", 5: "A", 7: "B", 11: "C"}


# ---------------------------------------------------------------------------
# Primzahl-Sieb & EABC-Klassifikation
# ---------------------------------------------------------------------------

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


def eabc_word_first_n(n: int) -> list[str]:
    """Erstellt das EABC-Wort aus den ersten `n` Primzahlen > 3."""
    word: list[str] = []
    candidate = 5
    while len(word) < n:
        if is_prime_trial(candidate):
            c = class_of(candidate)
            if c is not None:
                word.append(c)
        candidate += 2 if candidate % 2 == 1 else 1
    return word


def is_prime_trial(n: int) -> bool:
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for d in range(3, isqrt(n) + 1, 2):
        if n % d == 0:
            return False
    return True


def eabc_word_up_to(limit: int) -> list[str]:
    """EABC-Wort aller Primzahlen p > 3 bis `limit`."""
    primes = sieve_primes(limit)
    return [class_of(p) for p in primes if p > 3 and class_of(p) is not None]


# ---------------------------------------------------------------------------
# Pfadsignatur — effiziente O(n·|Σ|²)-Implementierung
# ---------------------------------------------------------------------------

class Signature(NamedTuple):
    level1: dict[str, int]
    level2: dict[tuple[str, str], int]
    level3: dict[tuple[str, str, str], int]
    n: int


def eabc_signature(word: list[str], max_level: int = 3) -> Signature:
    """
    Berechnet die Pfadsignatur (Chen-Signatur) eines diskreten Wortes
    über dem Alphabet {E, A, B, C} bis zu `max_level`.

    Algorithmus: Iterativer Durchlauf, O(n·|Σ|²).
    Stufe 2: N_XY = #{(i,j): i<j, w[i]=X, w[j]=Y}
    Stufe 3: N_XYZ = #{(i,j,k): i<j<k, w[i]=X, w[j]=Y, w[k]=Z}
    """
    # Initialisierung
    l1: dict[str, int] = {x: 0 for x in LABELS}
    l2: dict[tuple[str, str], int] = {(x, y): 0 for x in LABELS for y in LABELS}
    l3: dict[tuple[str, str, str], int] = {}
    if max_level >= 3:
        l3 = {(x, y, z): 0 for x in LABELS for y in LABELS for z in LABELS}

    for c in word:
        if max_level >= 3:
            # Level-3-Update muss VOR Level-2-Update stattfinden
            for x in LABELS:
                for y in LABELS:
                    if l2[(x, y)] > 0:
                        l3[(x, y, c)] += l2[(x, y)]

        # Level-2-Update muss VOR Level-1-Update stattfinden
        for x in LABELS:
            if l1[x] > 0:
                l2[(x, c)] += l1[x]

        # Level-1-Update
        l1[c] += 1

    return Signature(level1=l1, level2=l2, level3=l3, n=len(word))


# ---------------------------------------------------------------------------
# Zufallswort-Nullmodell (gleiche Länge, gleichverteilte Buchstaben)
# ---------------------------------------------------------------------------

def random_word(n: int, seed: int | None = None) -> list[str]:
    rng = random.Random(seed)
    return [rng.choice(LABELS) for _ in range(n)]


def random_word_matched_marginals(word: list[str], seed: int | None = None) -> list[str]:
    """Zufallspermutation des tatsächlichen Wortes (marginaltreu)."""
    shuffled = list(word)
    rng = random.Random(seed)
    rng.shuffle(shuffled)
    return shuffled


# ---------------------------------------------------------------------------
# Auswertung & Statistiken
# ---------------------------------------------------------------------------

def level2_matrix(sig: Signature) -> np.ndarray:
    """4×4-Matrix: M[i,j] = N_{LABELS[i], LABELS[j]}."""
    M = np.zeros((4, 4), dtype=np.float64)
    for i, x in enumerate(LABELS):
        for j, y in enumerate(LABELS):
            M[i, j] = sig.level2[(x, y)]
    return M


def expected_level2_independence(sig: Signature) -> np.ndarray:
    """
    Erwartungswert unter Zufallspermutation (feste Marginalverteilung):
    Für X≠Y: E[N_XY] = N_X * N_Y / 2
      (Herleitung: für jedes Paar (i<j) ist P(w[i]=X, w[j]=Y) = N_X/n * N_Y/(n-1),
       summiert über C(n,2) Paare ergibt n*(n-1)/2 * N_X/n * N_Y/(n-1) = N_X*N_Y/2)
    Für X=Y: N_XX = C(N_X, 2) = N_X*(N_X-1)/2 (deterministisch aus Level-1)
    """
    M_exp = np.zeros((4, 4), dtype=np.float64)
    for i, x in enumerate(LABELS):
        for j, y in enumerate(LABELS):
            nx = sig.level1[x]
            ny = sig.level1[y]
            if x == y:
                M_exp[i, j] = nx * (nx - 1) / 2
            else:
                M_exp[i, j] = nx * ny / 2
    return M_exp


def relative_deviation(sig: Signature) -> np.ndarray:
    """
    Relative Abweichung von Unabhängigkeit:
    rho_XY = N_XY / E[N_XY] - 1
    Wert 0 bedeutet perfekte Unabhängigkeit.
    """
    M_obs = level2_matrix(sig)
    M_exp = expected_level2_independence(sig)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = np.where(M_exp > 0, M_obs / M_exp - 1, np.nan)
    return rho


def asymmetry_matrix(sig: Signature) -> np.ndarray:
    """
    A_XY = N_XY - N_YX (Nichtkommutativität)
    Misst, ob X→Y-Übergänge häufiger sind als Y→X.
    """
    M = level2_matrix(sig)
    return M - M.T


def compute_d_e(sig: Signature) -> int:
    """
    D_E = N_A + N_C - N_B - N_E
    Netto-Zirkulationsmaß (Level-1-Information).
    """
    l1 = sig.level1
    return l1["A"] + l1["C"] - l1["B"] - l1["E"]


def signature_summary(sig: Signature, label: str) -> dict:
    """Zusammenfassung aller Signaturdaten."""
    d_e = compute_d_e(sig)
    M = level2_matrix(sig)
    rho = relative_deviation(sig)
    asym = asymmetry_matrix(sig)

    # Frobenius-Norm der relativen Abweichung (ohne Diagonale, da diese
    # durch C(N_X,2) erklärt werden kann)
    mask_offdiag = ~np.eye(4, dtype=bool)
    rho_offdiag = rho[mask_offdiag]
    rho_rms = float(np.sqrt(np.nanmean(rho_offdiag**2)))

    return {
        "label": label,
        "n": sig.n,
        "level1": {k: sig.level1[k] for k in LABELS},
        "D_E": d_e,
        "level2_matrix": M.tolist(),
        "rho_matrix": rho.tolist(),
        "asymmetry_matrix": asym.tolist(),
        "rho_rms_offdiag": rho_rms,
        "max_abs_rho": float(np.nanmax(np.abs(rho))),
        "noncommutativity_norm": float(np.linalg.norm(asym, "fro")),
    }


# ---------------------------------------------------------------------------
# Visualisierung
# ---------------------------------------------------------------------------

def plot_heatmap_grid(signatures: list[tuple[str, Signature]],
                      random_sig: Signature,
                      out_path: Path) -> None:
    """
    Erstellt ein Heatmap-Raster für Level-2-Matrizen verschiedener Primzahlbereiche.
    Zeile 1: N_XY (normalisiert)
    Zeile 2: Relative Abweichung rho_XY
    """
    n_sigs = len(signatures)
    fig, axes = plt.subplots(
        3, n_sigs + 1,
        figsize=(4.0 * (n_sigs + 1), 11),
        gridspec_kw={"hspace": 0.45, "wspace": 0.35},
    )
    fig.suptitle(
        "EABC-Pfadsignatur Stufe 2: Level-2-Heatmaps\n"
        r"$N_{XY}$ = geordnete Teilwörter (i<j, w[i]=X, w[j]=Y)",
        fontsize=13,
    )

    def _draw_heatmap(ax, data, title, fmt=".2f", cmap="RdYlGn", center_val=None,
                      vmin=None, vmax=None, annot_size=8):
        if vmin is None:
            vmin = np.nanmin(data)
        if vmax is None:
            vmax = np.nanmax(data)
        if center_val is not None:
            abs_max = max(abs(vmin - center_val), abs(vmax - center_val))
            vmin = center_val - abs_max
            vmax = center_val + abs_max
        im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        plt.colorbar(im, ax=ax, shrink=0.8)
        ax.set_xticks(range(4))
        ax.set_yticks(range(4))
        ax.set_xticklabels(LABELS)
        ax.set_yticklabels(LABELS)
        ax.set_xlabel("Y", fontsize=8)
        ax.set_ylabel("X", fontsize=8)
        ax.set_title(title, fontsize=8, pad=4)
        for i in range(4):
            for j in range(4):
                val = data[i, j]
                if np.isnan(val):
                    txt = "—"
                else:
                    txt = f"{val:{fmt}}"
                ax.text(j, i, txt, ha="center", va="center",
                        fontsize=annot_size, color="black")

    # Zufalls-Nullmodell immer in der letzten Spalte
    all_sigs_with_rand = list(signatures) + [("Zufall\n(Shuffle)", random_sig)]

    for col_idx, (lbl, sig) in enumerate(all_sigs_with_rand):
        M = level2_matrix(sig)
        M_norm = M / sig.n**2  # normiert auf [0,1]-Skala
        rho = relative_deviation(sig)
        asym = asymmetry_matrix(sig)
        asym_norm = asym / max(sig.n, 1)

        _draw_heatmap(axes[0, col_idx], M_norm,
                      f"{lbl}\n$N_{{XY}}/n^2$",
                      fmt=".4f", cmap="Blues", vmin=0)
        _draw_heatmap(axes[1, col_idx], rho,
                      f"{lbl}\n$\\rho_{{XY}}$ (Abw. v. Unabh.)",
                      fmt=".3f", cmap="RdBu_r", center_val=0)
        _draw_heatmap(axes[2, col_idx], asym_norm,
                      f"{lbl}\n$(N_{{XY}}-N_{{YX}})/n$",
                      fmt=".4f", cmap="PRGn", center_val=0)

    axes[0, 0].set_ylabel("$N_{XY}/n^2$\nX", fontsize=8)
    axes[1, 0].set_ylabel("$\\rho_{XY}$\nX", fontsize=8)
    axes[2, 0].set_ylabel("Asymmetrie\nX", fontsize=8)

    fig.subplots_adjust(top=0.93, hspace=0.55, wspace=0.40)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Heatmap gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Hauptanalyse
# ---------------------------------------------------------------------------

def run_analysis(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    sizes = [1_000, 10_000, 100_000]
    size_labels = ["1 000", "10 000", "100 000"]

    print("=" * 65)
    print("EABC-Pfadsignatur-Analyse")
    print("Klassifikation: p mod 12 → {E=1, A=5, B=7, C=11}")
    print("=" * 65)

    # -----------------------------------------------------------------------
    # Erzeuge Primzahlwörter
    # -----------------------------------------------------------------------
    # Wir brauchen die ersten 100 000 Primzahlen > 3.
    # Die 100 000. Primzahl > 3 ≈ 1 299 709 → Sieb bis 1 500 000 sicher genug.
    print("\nErzeuge EABC-Wörter aus Primzahlsieb …")
    primes = sieve_primes(1_600_000)
    full_word = [class_of(p) for p in primes if p > 3 and class_of(p) is not None]
    print(f"  Primzahlen > 3 bis 1 600 000: {len(full_word)}")

    if len(full_word) < 100_000:
        print("  Warnung: Nicht genug Primzahlen, erhöhe Siebgrenze …")
        primes = sieve_primes(2_000_000)
        full_word = [class_of(p) for p in primes if p > 3 and class_of(p) is not None]

    # -----------------------------------------------------------------------
    # Signaturen berechnen
    # -----------------------------------------------------------------------
    signatures: list[tuple[str, Signature]] = []
    results_rows: list[dict] = []

    for n, lbl in zip(sizes, size_labels):
        word = full_word[:n]
        print(f"\nBerechne Signatur für n={lbl} …")
        sig = eabc_signature(word, max_level=3)
        signatures.append((f"Prim\nn={lbl}", sig))

        # Nullmodell: gleiche Marginalverteilung, gemischt
        rand_word = random_word_matched_marginals(word, seed=42)
        rand_sig = eabc_signature(rand_word, max_level=3)

        summary_prim = signature_summary(sig, f"EABC Prim n={n}")
        summary_rand = signature_summary(rand_sig, f"Zufall n={n}")

        print(f"  Level-1: {sig.level1}")
        print(f"  D_E = {compute_d_e(sig)}")
        print(f"  rho_RMS (Off-Diag): EABC={summary_prim['rho_rms_offdiag']:.5f}  "
              f"Zufall={summary_rand['rho_rms_offdiag']:.5f}")
        print(f"  Max |rho|:          EABC={summary_prim['max_abs_rho']:.5f}  "
              f"Zufall={summary_rand['max_abs_rho']:.5f}")
        print(f"  Asymmetrie-Norm:    EABC={summary_prim['noncommutativity_norm']:.1f}  "
              f"Zufall={summary_rand['noncommutativity_norm']:.1f}")

        # Zeige Level-2-Relative-Abweichung
        rho = np.array(summary_prim["rho_matrix"])
        print(f"  rho-Matrix (EABC n={n}):")
        header = "      " + "  ".join(f"{y:>8}" for y in LABELS)
        print(header)
        for i, x in enumerate(LABELS):
            row_vals = "  ".join(
                f"{rho[i, j]:>+8.4f}" if not np.isnan(rho[i, j]) else "     nan" 
                for j in range(4)
            )
            print(f"  {x}:  {row_vals}")

        for row in [summary_prim, summary_rand]:
            flat = {
                "label": row["label"],
                "n": row["n"],
                "D_E": row["D_E"],
                "rho_rms_offdiag": row["rho_rms_offdiag"],
                "max_abs_rho": row["max_abs_rho"],
                "noncommutativity_norm": row["noncommutativity_norm"],
            }
            for labi in LABELS:
                flat[f"N_{labi}"] = row["level1"][labi]
            for i, xi in enumerate(LABELS):
                for j, yj in enumerate(LABELS):
                    flat[f"N_{xi}{yj}"] = round(np.array(row["level2_matrix"])[i, j])
                    flat[f"rho_{xi}{yj}"] = round(np.array(row["rho_matrix"])[i, j], 6)
            results_rows.append(flat)

    # -----------------------------------------------------------------------
    # Zufallsmodell für Heatmap (n=100 000, gemischt)
    # -----------------------------------------------------------------------
    rand_word_100k = random_word_matched_marginals(full_word[:100_000], seed=99)
    rand_sig_100k = eabc_signature(rand_word_100k, max_level=3)

    # -----------------------------------------------------------------------
    # Zusatz: Level-1-Analyse (D_E)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("Level-1-Analyse: D_E = N_A + N_C - N_B - N_E")
    print("(Dies ist Level-1-Information — nur von Häufigkeiten abhängig)")
    for (sig_label, sig), n, lbl_n in zip(signatures, sizes, size_labels):
        d = compute_d_e(sig)
        q = sig.n
        print(f"  n={lbl_n:>8}:  D_E={d:+5d}  W_E=D_E/n={d/q:+.5f}")

    # -----------------------------------------------------------------------
    # Level-2-Nichtkommutativitäts-Test: M_AB vs M_BA
    # -----------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("Nichtkommutativitäts-Test: N_XY vs N_YX für n=100 000")
    sig_100k = signatures[2][1]
    print(f"  {'Paar':>8}  {'N_XY':>10}  {'N_YX':>10}  {'Diff':>10}  {'Asymm%':>9}")
    for i, x in enumerate(LABELS):
        for j, y in enumerate(LABELS):
            if j <= i:
                continue
            n_xy = sig_100k.level2[(x, y)]
            n_yx = sig_100k.level2[(y, x)]
            diff = n_xy - n_yx
            sym = (n_xy + n_yx)
            asym_pct = 100 * diff / sym if sym > 0 else 0
            print(f"  {x}-{y}:  {n_xy:>10}  {n_yx:>10}  {diff:>+10}  {asym_pct:>+8.3f}%")

    # -----------------------------------------------------------------------
    # Level-3 Beispiel (ausgewählte Tripel)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("Level-3-Pfadsignatur (ausgewählte Tripel, n=100 000):")
    sig_100k = signatures[2][1]
    triplets_to_show = [
        ("A", "B", "C"), ("A", "B", "E"), ("E", "A", "B"), ("C", "E", "A"),
        ("A", "B", "C"), ("C", "E", "A"), ("B", "C", "E"), ("E", "A", "B"),
    ]
    # Deduplizieren, Reihenfolge behalten
    seen: set = set()
    for trip in triplets_to_show:
        if trip not in seen:
            seen.add(trip)
            val = sig_100k.level3.get(trip, 0)
            print(f"  N_{trip[0]}{trip[1]}{trip[2]} = {val:>12}")

    # -----------------------------------------------------------------------
    # Heatmap
    # -----------------------------------------------------------------------
    print("\nErstelle Heatmap …")
    heatmap_path = out_dir / "eabc_signature_heatmap.png"
    plot_heatmap_grid(signatures, rand_sig_100k, heatmap_path)

    # -----------------------------------------------------------------------
    # CSV-Export
    # -----------------------------------------------------------------------
    csv_path = out_dir / "eabc_signature_results.csv"
    if results_rows:
        fieldnames = list(results_rows[0].keys())
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results_rows)
        print(f"CSV gespeichert: {csv_path}")

    # -----------------------------------------------------------------------
    # Abschlusskfazit
    # -----------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("FAZIT: Ist die Level-2-Pfadsignatur arithmetisch nichttrivial?")
    print("=" * 65)

    sig_100k = signatures[2][1]
    summary = signature_summary(sig_100k, "EABC n=100000")
    rho = np.array(summary["rho_matrix"])

    rand_sum = signature_summary(rand_sig_100k, "Zufall n=100000")
    rho_rand = np.array(rand_sum["rho_matrix"])

    mask = ~np.eye(4, dtype=bool)
    rho_max = float(np.nanmax(np.abs(rho[mask])))
    rho_rand_max = float(np.nanmax(np.abs(rho_rand[mask])))

    print(f"\n  Max |rho_XY| (Off-Diag) für EABC-Strom:     {rho_max:.5f}")
    print(f"  Max |rho_XY| (Off-Diag) für Zufallspermut.: {rho_rand_max:.5f}")
    print()

    if rho_max > 10 * rho_rand_max:
        verdict = "JA — Level-2-Signatur ist stark nichttrivial"
    elif rho_max > 3 * rho_rand_max:
        verdict = "WAHRSCHEINLICH JA — Level-2-Signatur weicht signifikant ab"
    elif rho_max > 1.5 * rho_rand_max:
        verdict = "SCHWACH — moderate Abweichung, aber nicht eindeutig"
    else:
        verdict = "NEIN — Level-2-Signatur ist konsistent mit Zufallserwartung"

    print(f"  Befund: {verdict}")
    print()
    print("  Interpretation:")
    print("  • rho_XY = N_XY / E[N_XY] - 1 misst relative Abweichung von Unabhängigkeit")
    print("  • rho_XY > 0: Paar (X,Y) erscheint häufiger als erwartet")
    print("  • rho_XY < 0: Paar (X,Y) erscheint seltener als erwartet")
    print()
    print("  Stärkste Abweichungen (n=100 000):")
    flat_rho = [(LABELS[i], LABELS[j], rho[i, j])
                for i in range(4) for j in range(4) if not np.isnan(rho[i, j])]
    flat_rho.sort(key=lambda x: -abs(x[2]))
    for x, y, r in flat_rho[:8]:
        n_xy = sig_100k.level2[(x, y)]
        print(f"    N_{x}{y}={n_xy:>10}  rho_{x}{y}={r:+.5f}")

    print()
    print("  D_E ist Level-1-Information:")
    print("  D_E hängt nur von N_A, N_B, N_C, N_E ab — wird durch die")
    print("  Level-1-Marginalverteilung vollständig bestimmt. Die Level-2-")
    print("  Signatur enthält zusätzliche Information über geordnete Paare.")
    print("  Konkret: Selbst wenn D_E=0, kann N_AB ≠ N_BA sein.")


# ---------------------------------------------------------------------------
# Einstiegspunkt
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    out_dir = Path(__file__).resolve().parent
    run_analysis(out_dir)
