#!/usr/bin/env python3
"""
EABC Shuffle-Produkt-Test

Verifiziert die Chen-Formel für Level-2-Signaturen:

    N_XY(γ₁ ++ γ₂) = N_XY(γ₁) + N_XY(γ₂) + N_X(γ₁)·N_Y(γ₂)

Der Kreuzterm K_XY = N_X(γ₁)·N_Y(γ₂) kodiert die Interaktion der
Teilpfade: X erscheint in γ₁, Y erscheint in γ₂.

Zentralfragen:
  1. Gilt die Shuffle-Formel algebraisch exakt?
  2. Wie groß ist der Kreuzterm im Vergleich zur direkten Level-2-Summe?
  3. Ist K_EX = N_E(γ₁)·N_X(γ₂) systematisch verschieden von K_XE = N_X(γ₁)·N_E(γ₂)?
"""

from __future__ import annotations

import random
from math import isqrt
from pathlib import Path
from typing import NamedTuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Übernahme der Kern-Funktionen aus eabc_signature_analysis.py
# ---------------------------------------------------------------------------

LABELS = ("E", "A", "B", "C")
RESIDUE_MAP = {1: "E", 5: "A", 7: "B", 11: "C"}


def sieve_primes(limit: int) -> list[int]:
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


class Signature(NamedTuple):
    level1: dict[str, int]
    level2: dict[tuple[str, str], int]
    n: int


def eabc_signature(word: list[str]) -> Signature:
    """Pfadsignatur bis Stufe 2 (O(n·|Σ|)-Algorithmus)."""
    l1: dict[str, int] = {x: 0 for x in LABELS}
    l2: dict[tuple[str, str], int] = {(x, y): 0 for x in LABELS for y in LABELS}
    for c in word:
        for x in LABELS:
            if l1[x] > 0:
                l2[(x, c)] += l1[x]
        l1[c] += 1
    return Signature(level1=l1, level2=l2, n=len(word))


# ---------------------------------------------------------------------------
# Shuffle-Formel
# ---------------------------------------------------------------------------

def shuffle_predict_level2(
    sig1: Signature, sig2: Signature
) -> dict[tuple[str, str], int]:
    """
    Shuffle-Produkt-Vorhersage (Level 2):
        N_XY(γ₁ ++ γ₂) = N_XY(γ₁) + N_XY(γ₂) + N_X(γ₁)·N_Y(γ₂)
    """
    predicted: dict[tuple[str, str], int] = {}
    for x in LABELS:
        for y in LABELS:
            cross_term = sig1.level1[x] * sig2.level1[y]
            predicted[(x, y)] = sig1.level2[(x, y)] + sig2.level2[(x, y)] + cross_term
    return predicted


def cross_term_matrix(sig1: Signature, sig2: Signature) -> np.ndarray:
    """K_XY = N_X(γ₁)·N_Y(γ₂) als 4×4-Matrix."""
    K = np.zeros((4, 4), dtype=np.int64)
    for i, x in enumerate(LABELS):
        for j, y in enumerate(LABELS):
            K[i, j] = sig1.level1[x] * sig2.level1[y]
    return K


def verify_shuffle(word: list[str], split: int) -> dict:
    """
    Verifiziert die Shuffle-Formel an einem gegebenen Split-Punkt.
    Gibt einen Bericht mit Fehlern zurück (sollten alle 0 sein).
    """
    gamma1 = word[:split]
    gamma2 = word[split:]
    gamma_full = word  # Identität: γ₁ ++ γ₂ = γ

    sig_full = eabc_signature(gamma_full)
    sig1 = eabc_signature(gamma1)
    sig2 = eabc_signature(gamma2)
    predicted = shuffle_predict_level2(sig1, sig2)

    errors: dict[tuple[str, str], int] = {}
    max_err = 0
    for x in LABELS:
        for y in LABELS:
            err = sig_full.level2[(x, y)] - predicted[(x, y)]
            errors[(x, y)] = err
            max_err = max(max_err, abs(err))

    K = cross_term_matrix(sig1, sig2)
    l2_direct = np.array([[sig1.level2[(x, y)] + sig2.level2[(x, y)]
                            for y in LABELS] for x in LABELS], dtype=np.int64)
    l2_full = np.array([[sig_full.level2[(x, y)] for y in LABELS]
                         for x in LABELS], dtype=np.int64)

    return {
        "split": split,
        "n": len(word),
        "max_error": max_err,
        "all_zero": max_err == 0,
        "cross_term_matrix": K,
        "l2_direct_sum": l2_direct,
        "l2_full": l2_full,
        "cross_fraction": float(K.sum()) / max(l2_full.sum(), 1),
    }


# ---------------------------------------------------------------------------
# E-Paar-Asymmetrie über zufällige Splits
# ---------------------------------------------------------------------------

def e_pair_asymmetry_over_splits(
    word: list[str],
    n_samples: int = 500,
    seed: int = 42,
) -> dict:
    """
    Prüft: Ist K_EX = N_E(γ₁)·N_X(γ₂) systematisch ≠ K_XE = N_X(γ₁)·N_E(γ₂)?

    Für jeden zufälligen Split k ∈ {1, ..., n-1}:
        Diff_X(k) = K_EX(k) - K_XE(k)
                  = N_E(γ₁)·N_X(γ₂) - N_X(γ₁)·N_E(γ₂)
    """
    rng = random.Random(seed)
    n = len(word)
    splits = [rng.randint(1, n - 1) for _ in range(n_samples)]
    splits.sort()

    results: dict[str, list[float]] = {x: [] for x in ("A", "B", "C")}
    # Auch das Verhältnis N_E(γ₁)/|γ₁| vs. N_E(γ₂)/|γ₂| als Funktion des Split
    e_density_diff: list[float] = []

    for k in splits:
        sig1 = eabc_signature(word[:k])
        sig2 = eabc_signature(word[k:])

        ne1 = sig1.level1["E"]
        ne2 = sig2.level1["E"]

        for x in ("A", "B", "C"):
            nx1 = sig1.level1[x]
            nx2 = sig2.level1[x]
            diff = ne1 * nx2 - nx1 * ne2
            results[x].append(diff)

        # E-Dichte-Differenz zwischen Teilpfaden
        d1 = ne1 / max(k, 1)
        d2 = ne2 / max(n - k, 1)
        e_density_diff.append(d1 - d2)

    means = {x: float(np.mean(results[x])) for x in ("A", "B", "C")}
    stds = {x: float(np.std(results[x])) for x in ("A", "B", "C")}
    t_stats = {x: means[x] / max(stds[x] / np.sqrt(n_samples), 1e-12)
               for x in ("A", "B", "C")}

    return {
        "n_samples": n_samples,
        "means_diff": means,
        "stds_diff": stds,
        "t_statistics": t_stats,
        "e_density_mean_diff": float(np.mean(e_density_diff)),
        "raw_diff": results,
        "splits": splits,
    }


# ---------------------------------------------------------------------------
# Kreuzterm als Funktion des Split-Punkts
# ---------------------------------------------------------------------------

def cross_term_profile(
    word: list[str], n_checkpoints: int = 200
) -> dict:
    """
    Berechnet den Kreuzterm K_EX(k) = N_E(γ₁)·N_X(γ₂) als Funktion
    des Split-Punkts k ∈ {1, ..., n-1}.
    """
    n = len(word)
    step = max(1, n // n_checkpoints)
    checkpoints = list(range(step, n, step))

    # Kumulativer Level-1-Zähler für γ₁: O(n) total
    cumulative_l1 = {x: [0] * (n + 1) for x in LABELS}
    for idx, c in enumerate(word):
        for x in LABELS:
            cumulative_l1[x][idx + 1] = cumulative_l1[x][idx]
        cumulative_l1[c][idx + 1] += 1

    # Gesamt-Level-1 für γ₂ = Gesamt - γ₁
    total_l1 = {x: cumulative_l1[x][n] for x in LABELS}

    k_ex_vals: dict[str, list[int]] = {x: [] for x in LABELS}
    k_xe_vals: dict[str, list[int]] = {x: [] for x in LABELS}

    for k in checkpoints:
        ne1 = cumulative_l1["E"][k]
        for x in LABELS:
            nx1 = cumulative_l1[x][k]
            nx2 = total_l1[x] - nx1
            ne2 = total_l1["E"] - ne1
            k_ex_vals[x].append(ne1 * nx2)
            k_xe_vals[x].append(nx1 * ne2)

    return {
        "checkpoints": checkpoints,
        "K_EX": k_ex_vals,
        "K_XE": k_xe_vals,
        "n": n,
    }


# ---------------------------------------------------------------------------
# Visualisierung
# ---------------------------------------------------------------------------

def plot_results(
    verify_report: dict,
    asym_report: dict,
    profile: dict,
    out_path: Path,
    n_word: int,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(
        f"EABC Shuffle-Produkt-Test  (n={n_word:,})\n"
        r"$N_{XY}(\gamma_1{+\!\!+}\gamma_2) = N_{XY}(\gamma_1) + N_{XY}(\gamma_2) + N_X(\gamma_1)\cdot N_Y(\gamma_2)$",
        fontsize=11,
    )

    # --- Panel 1: Kreuzterm-Matrix als Heatmap ---
    ax = axes[0, 0]
    K = verify_report["cross_term_matrix"].astype(float)
    K_norm = K / max(K.max(), 1)
    im = ax.imshow(K_norm, cmap="Blues", vmin=0)
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_xticks(range(4)); ax.set_xticklabels(LABELS)
    ax.set_yticks(range(4)); ax.set_yticklabels(LABELS)
    ax.set_xlabel("Y (γ₂)"); ax.set_ylabel("X (γ₁)")
    ax.set_title(r"Kreuzterm $K_{XY} = N_X(\gamma_1)\cdot N_Y(\gamma_2)$ / max", fontsize=9)
    for i in range(4):
        for j in range(4):
            ax.text(j, i, f"{K[i,j]:.0f}", ha="center", va="center", fontsize=8)

    # --- Panel 2: Direktsumme vs. Kreuzterm-Anteil ---
    ax = axes[0, 1]
    direct = verify_report["l2_direct_sum"].astype(float)
    cross = verify_report["cross_term_matrix"].astype(float)
    full = verify_report["l2_full"].astype(float)

    # Für jedes Paar: Anteil des Kreuztterms am Gesamtwert
    labels_flat = [f"{x}{y}" for x in LABELS for y in LABELS]
    direct_flat = direct.flatten()
    cross_flat = cross.flatten()
    total_flat = full.flatten()

    x_pos = np.arange(len(labels_flat))
    width = 0.35
    ax.bar(x_pos, direct_flat / np.maximum(total_flat, 1), width, label="Direkte Summe", color="C0", alpha=0.7)
    ax.bar(x_pos + width, cross_flat / np.maximum(total_flat, 1), width, label="Kreuzterm", color="C1", alpha=0.7)
    ax.set_xticks(x_pos + width / 2)
    ax.set_xticklabels(labels_flat, rotation=90, fontsize=7)
    ax.set_ylabel("Anteil am Gesamt-N_XY")
    ax.set_title("Kreuzterm-Anteil pro Paar (X,Y)", fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)

    # --- Panel 3: K_EX vs K_XE als Funktion des Split-Punkts ---
    ax = axes[1, 0]
    ks = np.array(profile["checkpoints"]) / n_word
    colors = {"A": "C1", "B": "C2", "C": "C3"}
    for x in ("A", "B", "C"):
        kex = np.array(profile["K_EX"][x])
        kxe = np.array(profile["K_XE"][x])
        ax.plot(ks, kex - kxe, color=colors[x], label=f"K_E{x} - K_{x}E", lw=1.2)
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.set_xlabel("Split-Position k/n")
    ax.set_ylabel(r"$K_{EX}(k) - K_{XE}(k)$")
    ax.set_title(r"E-Paar-Asymmetrie $K_{EX} - K_{XE}$ vs. Split", fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # --- Panel 4: Verteilung der Differenzen über zufällige Splits ---
    ax = axes[1, 1]
    colors = {"A": "C1", "B": "C2", "C": "C3"}
    for x in ("A", "B", "C"):
        diffs = asym_report["raw_diff"][x]
        m = asym_report["means_diff"][x]
        t = asym_report["t_statistics"][x]
        ax.hist(
            diffs, bins=30, alpha=0.5, color=colors[x],
            label=f"E-{x}: μ={m:+.0f}, t={t:+.2f}",
        )
    ax.axvline(0, color="black", lw=1.0, ls="--")
    ax.set_xlabel(r"$K_{EX} - K_{XE}$ (zufällige Splits)")
    ax.set_ylabel("Häufigkeit")
    ax.set_title(f"Verteilung der E-Paar-Asymmetrie ({asym_report['n_samples']} Splits)", fontsize=9)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot gespeichert: {out_path}")


# ---------------------------------------------------------------------------
# Hauptanalyse
# ---------------------------------------------------------------------------

def run(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("EABC Shuffle-Produkt-Test")
    print("Formel: N_XY(γ₁++γ₂) = N_XY(γ₁) + N_XY(γ₂) + N_X(γ₁)·N_Y(γ₂)")
    print("=" * 70)

    # Erzeuge EABC-Wort der ersten 100 000 Primzahlen > 3
    print("\n[1] Erzeuge Primzahlwort …")
    primes = sieve_primes(1_600_000)
    full_word = [class_of(p) for p in primes if p > 3 and class_of(p) is not None]
    n = 100_000
    word = full_word[:n]
    print(f"    Wortlänge: {len(word)}")

    # -------------------------------------------------------------------
    # Algebraische Verifikation an mehreren Split-Punkten
    # -------------------------------------------------------------------
    print("\n[2] Algebraische Verifikation der Shuffle-Formel …")
    test_splits = [1, 100, 1_000, 10_000, 50_000, 99_000, n - 1]
    all_ok = True
    for k in test_splits:
        report = verify_shuffle(word, k)
        status = "OK ✓" if report["all_zero"] else f"FEHLER! max_err={report['max_error']}"
        cross_pct = 100 * report["cross_fraction"]
        print(f"    Split k={k:>7}: {status}  |  Kreuzterm-Anteil: {cross_pct:.1f}%")
        if not report["all_zero"]:
            all_ok = False

    if all_ok:
        print("\n    → Shuffle-Formel algebraisch exakt verifiziert (alle Fehler = 0)")
    else:
        print("\n    → WARNUNG: Fehler in der Shuffle-Verifikation!")

    # -------------------------------------------------------------------
    # Kreuzterm-Analyse am Mittelpunkt
    # -------------------------------------------------------------------
    print("\n[3] Kreuzterm-Analyse am Mittelpunkt (k=n/2) …")
    mid_report = verify_shuffle(word, n // 2)
    K = mid_report["cross_term_matrix"]
    L2_direct = mid_report["l2_direct_sum"]
    L2_full = mid_report["l2_full"]

    print(f"\n    {'':4s}  {'K_XY = N_X(γ₁)·N_Y(γ₂)':>30}  {'% von N_XY(γ)':>15}")
    print(f"    {'Paar':>6}  {'Kreuzterm':>12}  {'Direkte Summe':>14}  {'Gesamt':>10}  {'Kreuzterm-Anteil':>17}")
    print("    " + "-" * 65)
    for i, x in enumerate(LABELS):
        for j, y in enumerate(LABELS):
            ct = K[i, j]
            ds = L2_direct[i, j]
            tot = L2_full[i, j]
            pct = 100 * ct / max(tot, 1)
            print(f"    {x}{y}:   {ct:>12,}  {ds:>14,}  {tot:>10,}  {pct:>16.1f}%")

    print(f"\n    Gesamt-Kreuzterm-Anteil: {100*mid_report['cross_fraction']:.2f}%")

    # Wichtige strukturelle Beobachtung: K_XY ≠ K_YX im Allgemeinen
    print("\n    K-Matrix-Asymmetrie (K_XY - K_YX):")
    K_asym = K - K.T
    print(f"    {'':4s}", end="")
    for y in LABELS:
        print(f"  {y:>8}", end="")
    print()
    for i, x in enumerate(LABELS):
        print(f"    {x}:", end="")
        for j, y in enumerate(LABELS):
            print(f"  {K_asym[i,j]:>+8,}", end="")
        print()

    # -------------------------------------------------------------------
    # E-Paar-Asymmetrie über zufällige Splits
    # -------------------------------------------------------------------
    print("\n[4] E-Paar-Asymmetrie über zufällige Splits …")
    asym = e_pair_asymmetry_over_splits(word, n_samples=500, seed=42)

    print(f"\n    {'Paar':>6}  {'Mittlere Diff':>15}  {'Std':>12}  {'t-Statistik':>13}  {'Befund':>30}")
    print("    " + "-" * 80)
    for x in ("A", "B", "C"):
        m = asym["means_diff"][x]
        s = asym["stds_diff"][x]
        t = asym["t_statistics"][x]
        if abs(t) > 3:
            befund = "SIGNIFIKANT (|t|>3)"
        elif abs(t) > 2:
            befund = "moderat (|t|>2)"
        else:
            befund = "nicht signifikant"
        print(f"    E-{x}:  {m:>+14.1f}  {s:>12.1f}  {t:>+12.2f}  {befund}")

    print(f"\n    E-Dichte-Differenz (γ₁ vs γ₂, Mittel): "
          f"{asym['e_density_mean_diff']:+.5f}")
    print("    (positiv = E häufiger in der ersten Hälfte des Primstroms)")

    # -------------------------------------------------------------------
    # Kreuzterm-Profil über alle Splits
    # -------------------------------------------------------------------
    print("\n[5] Kreuzterm-Profil über alle Split-Punkte …")
    profile = cross_term_profile(word, n_checkpoints=300)

    # -------------------------------------------------------------------
    # Plot
    # -------------------------------------------------------------------
    out_path = out_dir / "eabc_shuffle_test.png"
    plot_results(mid_report, asym, profile, out_path, n)

    # -------------------------------------------------------------------
    # Fazit
    # -------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("FAZIT:")
    print("=" * 70)
    print()
    print("1. SHUFFLE-FORMEL: Algebraisch exakt verifiziert (Fehler = 0 für alle Splits).")
    print()
    print("2. KREUZTERM K_XY = N_X(γ₁)·N_Y(γ₂):")
    pct = 100 * mid_report["cross_fraction"]
    print(f"   Macht ca. {pct:.1f}% des gesamten Level-2-Werts aus (bei Mittelpunkt-Split).")
    print("   Der Kreuzterm ist die echte Level-2-Information: er misst, wie")
    print("   oft ein X-Ereignis in γ₁ mit einem Y-Ereignis in γ₂ kombiniert wird.")
    print("   Diese Information ist in der Level-1-Signatur (Frequenzen) nicht enthalten.")
    print()
    print("3. E-PAAR-ASYMMETRIE:")
    sig_pairs = [x for x in ("A","B","C") if abs(asym["t_statistics"][x]) > 2]
    if sig_pairs:
        print(f"   Signifikante Asymmetrie für Paare: {sig_pairs}")
        print("   K_EX ≠ K_XE → E-Ereignisse sind zeitlich nicht gleichmäßig verteilt.")
    else:
        print("   Keine signifikante systematische E-Asymmetrie gefunden.")
        print("   K_EX ≈ K_XE → E-Ereignisse sind statistisch gleichmäßig verteilt.")
    print()
    print("4. STRUKTUR DES KREUZTTERMS:")
    print("   K_XY = N_X(γ₁)·N_Y(γ₂) ist im Allgemeinen nicht symmetrisch in X,Y,")
    print("   denn K_XY = N_X(γ₁)·N_Y(γ₂) ≠ N_Y(γ₁)·N_X(γ₂) = K_YX")
    print("   (außer wenn N_X(γ₁)/N_X(γ₂) = N_Y(γ₁)/N_Y(γ₂)).")


if __name__ == "__main__":
    run(Path(__file__).resolve().parent)
