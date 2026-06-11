#!/usr/bin/env python3
"""Erzeugt die Abbildungen für Random Vierlinge.tex → Random_Vierlinge_figures/*.pdf"""
from __future__ import annotations

import random
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parent
OUT = SCRIPT / "Random_Vierlinge_figures"
CSV_NAME = "primzahlvierlinge_1000.csv"

classes36 = [1, 5, 7, 11, 13, 17, 19, 23, 25, 29, 31, 35]
index36 = {r: i for i, r in enumerate(classes36)}


def class36(p: int) -> int:
    return int(p) % 36


def transition36(quadruplets: np.ndarray, max_gap: int) -> np.ndarray:
    M = np.zeros((12, 12), dtype=int)
    for i in range(len(quadruplets) - 1):
        p = int(quadruplets[i])
        q = int(quadruplets[i + 1])
        if q - p <= max_gap:
            r1, r2 = class36(p), class36(q)
            if r1 in index36 and r2 in index36:
                M[index36[r1], index36[r2]] += 1
    return M


def gap_class(d: int, *, s_max: int = 100, m_max: int = 500) -> str:
    if d <= s_max:
        return "S"
    if d <= m_max:
        return "M"
    return "L"


def seq_label(t: tuple[str, str, str]) -> str:
    return "".join(t)


def v_p(n: int, p: int) -> int:
    if n == 0:
        return 0
    c = 0
    while n % p == 0:
        n //= p
        c += 1
    return c


def shuffled_sequences(
    gaps: np.ndarray,
    *,
    trials: int = 1000,
    s_max: int = 100,
    m_max: int = 500,
    seed: int | None = 42,
) -> Counter[tuple[str, str, str]]:
    gaps = np.asarray(gaps, dtype=np.int64).ravel()
    if len(gaps) < 3:
        return Counter()
    rng = random.Random(seed) if seed is not None else random.Random()
    counts: Counter[tuple[str, str, str]] = Counter()
    g_list = [int(x) for x in gaps]
    for _ in range(trials):
        g = g_list.copy()
        rng.shuffle(g)
        for i in range(len(g) - 2):
            seq = (
                gap_class(g[i], s_max=s_max, m_max=m_max),
                gap_class(g[i + 1], s_max=s_max, m_max=m_max),
                gap_class(g[i + 2], s_max=s_max, m_max=m_max),
            )
            counts[seq] += 1
    return counts


def load_quadruplets(csv_path: Path) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(csv_path)
    keep = [c for c in df.columns if not str(c).startswith("Unnamed")]
    df = df[keep].copy()
    df["Mittel"] = pd.to_numeric(df["Mittel"], errors="raise")
    df = df.sort_values("Mittel").reset_index(drop=True)
    p_start = (df["Mittel"] - 4).astype(np.int64).to_numpy()
    centers = df["Mittel"].to_numpy(dtype=np.int64)
    return p_start, centers


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    csv = SCRIPT / CSV_NAME
    if not csv.is_file():
        raise SystemExit(f"CSV fehlt: {csv}")

    quadruplets, centers = load_quadruplets(csv)
    mod = 210
    counts210 = Counter(int(p) % mod for p in quadruplets)

    # --- (A) mod 210: alle Reste mit Treffer > 0 ---
    items = sorted(counts210.items(), key=lambda x: (-x[1], x[0]))
    labels_a = [str(r) for r, _ in items]
    vals_a = [c for _, c in items]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(labels_a, vals_a, color="steelblue", edgecolor="black", linewidth=0.5)
    ax.set_title("Verteilung der Vierling-Starts $p$ mod $210$")
    ax.set_xlabel("Restklasse")
    ax.set_ylabel("Anzahl")
    plt.xticks(rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(OUT / "figA_mod210.pdf")
    plt.close(fig)

    # --- (B) nur 11, 101, 191 normiert (Erwartung 1/3) ---
    main3 = [11, 101, 191]
    vals_b = [counts210[r] for r in main3]
    labels_b = [str(r) for r in main3]
    probs = np.array(vals_b, dtype=float) / sum(vals_b)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(labels_b, probs, color="darkseagreen", edgecolor="black", linewidth=0.5)
    ax.axhline(1 / 3, color="crimson", linestyle="--", linewidth=1.5, label="$1/3$ (gleichverteilt)")
    ax.set_ylim(0, max(probs.max() * 1.15, 1 / 3 * 1.15))
    ax.set_title("Normiert auf die drei Einheitsklassen mod $210$")
    ax.set_ylabel("Anteil")
    ax.set_xlabel("Restklasse")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "figB_normiert.pdf")
    plt.close(fig)

    # --- (C) Top-4 Gap-Tripel vs. Permutations-Baseline ---
    gaps = np.diff(quadruplets)
    seq_counter: Counter[tuple[str, str, str]] = Counter()
    for i in range(len(gaps) - 2):
        seq = (
            gap_class(int(gaps[i])),
            gap_class(int(gaps[i + 1])),
            gap_class(int(gaps[i + 2])),
        )
        seq_counter[seq] += 1
    top4 = seq_counter.most_common(4)
    trials = 1000
    baseline = shuffled_sequences(gaps, trials=trials, seed=42)
    labels_c = [seq_label(seq) for seq, _ in top4]
    obs = [cnt for _, cnt in top4]
    exp = [baseline[seq] / trials for seq, _ in top4]

    x = np.arange(len(obs))
    w = 0.35
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(x - w / 2, obs, width=w, label="beobachtet", color="steelblue", edgecolor="black", lw=0.5)
    ax.bar(
        x + w / 2,
        exp,
        width=w,
        label=f"Baseline (Mittel / Trial, $n={trials}$)",
        color="orange",
        edgecolor="black",
        lw=0.5,
        alpha=0.85,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels_c)
    ax.set_title("Häufigste Gap-Klassen-Tripel (S/M/L) vs. Shuffle")
    ax.set_ylabel("Anzahl bzw. Mittel")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "figC_gap_sequences.pdf")
    plt.close(fig)

    # --- (D) v3,v5,v7 Mittel small vs large gaps ---
    small_thr, large_thr = 100, 1000
    q = quadruplets
    c = centers
    data = []
    for i in range(len(q) - 1):
        m = int(c[i])
        gap = int(q[i + 1] - q[i])
        data.append({"gap": gap, "v3": v_p(m, 3), "v5": v_p(m, 5), "v7": v_p(m, 7)})
    small = [d for d in data if d["gap"] <= small_thr]
    large = [d for d in data if d["gap"] > large_thr]

    def avg(lst: list, key: str) -> float:
        return float(np.mean([x[key] for x in lst])) if lst else 0.0

    lab = ["$v_3$", "$v_5$", "$v_7$"]
    small_y = [avg(small, "v3"), avg(small, "v5"), avg(small, "v7")]
    large_y = [avg(large, "v3"), avg(large, "v5"), avg(large, "v7")]

    x = np.arange(3)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(x - 0.2, small_y, width=0.4, label=rf"kleiner Gap ($\Delta \leq {small_thr}$)", color="#4C72B0")
    ax.bar(x + 0.2, large_y, width=0.4, label=rf"großer Gap ($\Delta > {large_thr}$)", color="#DD8452")
    ax.set_xticks(x)
    ax.set_xticklabels(lab)
    ax.set_title("Mittlere $p$-adische Bewertungen des Zentrums $m=p+4$")
    ax.set_ylabel("Mittelwert")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "figD_glattheit.pdf")
    plt.close(fig)

    # --- (E) M36 heatmap ---
    max_gap = 1000
    M36 = transition36(quadruplets, max_gap=max_gap)
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    im = ax.imshow(M36, cmap="viridis", aspect="auto")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title(rf"Übergangsmatrix mod $36$ (Gap $\leq {max_gap}$)")
    ax.set_xticks(range(12))
    ax.set_xticklabels([str(c) for c in classes36], rotation=45, ha="right")
    ax.set_yticks(range(12))
    ax.set_yticklabels([str(c) for c in classes36])
    ax.set_xlabel("Ziel: $p$ mod 36")
    ax.set_ylabel("Start: $p$ mod 36")
    fig.tight_layout()
    fig.savefig(OUT / "figE_M36_heatmap.pdf")
    plt.close(fig)

    print("Gespeichert unter:", OUT)


if __name__ == "__main__":
    main()
