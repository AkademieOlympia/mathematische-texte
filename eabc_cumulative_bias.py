#!/usr/bin/env python3
"""
Kumulative Kanalzählungen N_i(x) und B_{j,i}(x) für mod-420-Vierlingskanäle.

B_{101,11}(x) = N_{101}(x) - N_{11}(x) entlang wachsender Vierlingssequenz.
Zeitliche Tests: Vorzeichenwechsel, Drift in B/√x, LIL-Vergleich, Random-Walk-Null.
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np

# Stiefel.py im gleichen Verzeichnis
sys.path.insert(0, str(Path(__file__).resolve().parent))
from Stiefel import HALF_CHANNELS, INDEX, channel_mod420, quadruplets  # noqa: E402

CHANNELS = HALF_CHANNELS
N_CHANNELS = len(CHANNELS)
IDX_11 = CHANNELS.index(11)
IDX_101 = CHANNELS.index(101)


def collect_channel_events(max_prime: int) -> tuple[np.ndarray, np.ndarray]:
    """Vierlingsstarts p und zugehörige Kanalrestklassen (nur gültige Kanäle)."""
    t0 = time.perf_counter()
    starts = quadruplets(max_prime)
    t_sieve = time.perf_counter() - t0

    primes: list[int] = []
    channels: list[int] = []
    for p in starts:
        r = channel_mod420(int(p))
        if r is not None:
            primes.append(int(p))
            channels.append(r)

    t_total = time.perf_counter() - t0
    return (
        np.array(primes, dtype=np.int64),
        np.array(channels, dtype=np.int64),
        t_sieve,
        t_total,
    )


def cumulative_counts(channels: np.ndarray) -> np.ndarray:
    """Shape (K, 6): kumulative N_i nach k Vierlingen."""
    k = len(channels)
    idx = np.array([INDEX[c] for c in channels], dtype=np.int64)
    one_hot = np.zeros((k, N_CHANNELS), dtype=np.int64)
    one_hot[np.arange(k), idx] = 1
    return np.cumsum(one_hot, axis=0)


def all_pair_biases(cum: np.ndarray) -> dict[tuple[int, int], np.ndarray]:
    """B_{j,i}(k) = N_j(k) - N_i(k) für alle i < j."""
    out: dict[tuple[int, int], np.ndarray] = {}
    for i, ci in enumerate(CHANNELS):
        for j in range(i + 1, N_CHANNELS):
            cj = CHANNELS[j]
            out[(cj, ci)] = cum[:, j] - cum[:, i]
    return out


def sign_changes(series: np.ndarray) -> tuple[int, list[int]]:
    """Anzahl und Indizes (1-basiert) der Vorzeichenwechsel; Nullen behalten Vorzeichen."""
    changes: list[int] = []
    prev_sign = 0
    for i, val in enumerate(series):
        if val == 0:
            continue
        s = 1 if val > 0 else -1
        if prev_sign != 0 and s != prev_sign:
            changes.append(i + 1)  # 1-basierter k-Index
        prev_sign = s
    return len(changes), changes


def linear_drift(y: np.ndarray, x: np.ndarray) -> tuple[float, float, float]:
    """OLS y ~ a*x + b; Rückgabe (a, b, r²)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2:
        return 0.0, float(y[0]) if len(y) else 0.0, 0.0
    xm, ym = x.mean(), y.mean()
    ss_xx = np.sum((x - xm) ** 2)
    if ss_xx == 0:
        return 0.0, ym, 0.0
    a = np.sum((x - xm) * (y - ym)) / ss_xx
    b = ym - a * xm
    y_hat = a * x + b
    ss_res = np.sum((y - y_hat) ** 2)
    ss_tot = np.sum((y - ym) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return float(a), float(b), float(r2)


def lil_bound(k: int) -> float:
    """Heuristische LIL-Grenze √(2k log log k), k≥3."""
    if k < 3:
        return math.sqrt(max(k, 1))
    return math.sqrt(2.0 * k * math.log(max(math.log(k), math.e)))


def random_walk_null(
    n_events: int,
    n_trials: int,
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    """
    Nullhypothese: gleichverteilte Kanäle (1/6 je Kanal).
    Schritt für B_{101,11}: +1 (101), -1 (11), 0 sonst.
    """
    steps = np.zeros(N_CHANNELS, dtype=float)
    steps[IDX_101] = 1.0
    steps[IDX_11] = -1.0

    max_abs_ratio = np.empty(n_trials)
    final_ratio = np.empty(n_trials)
    sign_changes_rw = np.empty(n_trials, dtype=int)

    for t in range(n_trials):
        draws = rng.integers(0, N_CHANNELS, size=n_events)
        walk = steps[draws].cumsum()
        k = np.arange(1, n_events + 1, dtype=float)
        ratio = walk / np.sqrt(k)
        max_abs_ratio[t] = np.max(np.abs(ratio))
        final_ratio[t] = ratio[-1]
        sc, _ = sign_changes(walk)
        sign_changes_rw[t] = sc

    return {
        "max_abs_ratio": max_abs_ratio,
        "final_ratio": final_ratio,
        "sign_changes": sign_changes_rw,
    }


def ascii_sparkline(values: np.ndarray, width: int = 60) -> str:
    """Einfache ASCII-Sparkline für normalisierte Werte in [-1,1]."""
    if len(values) == 0:
        return ""
    chars = " .:-=+*#%@"
    v = np.asarray(values, dtype=float)
    vmin, vmax = v.min(), v.max()
    if abs(vmax - vmin) < 1e-15:
        mid = chars[len(chars) // 2]
        return mid * min(width, len(v))
    idx = np.linspace(0, len(v) - 1, width).astype(int)
    sampled = v[idx]
    norm = (sampled - vmin) / (vmax - vmin)
    return "".join(chars[min(int(x * (len(chars) - 1)), len(chars) - 1)] for x in norm)


def save_csv(
    path: Path,
    primes: np.ndarray,
    cum: np.ndarray,
    b_101_11: np.ndarray,
) -> None:
    header = "k,p," + ",".join(f"N_{c}" for c in CHANNELS)
    header += ",B_101_11,B_101_11_over_sqrt_k"
    lines = [header]
    k_arr = np.arange(1, len(primes) + 1)
    sqrt_k = np.sqrt(k_arr.astype(float))
    for i in range(len(primes)):
        row = [str(i + 1), str(primes[i])]
        row.extend(str(int(cum[i, j])) for j in range(N_CHANNELS))
        row.append(str(int(b_101_11[i])))
        row.append(f"{b_101_11[i] / sqrt_k[i]:.8f}")
        lines.append(",".join(row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def maybe_plot(out_dir: Path, k: np.ndarray, b: np.ndarray, ratio: np.ndarray) -> str:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return "matplotlib nicht verfügbar — nur ASCII/CSV"

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    axes[0].plot(k, b, lw=0.8, color="C0")
    axes[0].axhline(0, color="gray", ls="--", lw=0.5)
    axes[0].set_ylabel("B_{101,11}(k)")
    axes[0].set_title("Kumulativer Kanal-Bias B_{101,11}")

    axes[1].plot(k, ratio, lw=0.8, color="C1", label="Beobachtet")
    axes[1].axhline(0, color="gray", ls="--", lw=0.5)
    lil = np.array([lil_bound(int(kk)) for kk in k])
    axes[1].plot(k, lil / np.sqrt(k), "g--", lw=0.6, alpha=0.7, label="LIL/√k (heuristisch)")
    axes[1].plot(k, -lil / np.sqrt(k), "g--", lw=0.6, alpha=0.7)
    axes[1].set_xlabel("k (Vierlingsindex)")
    axes[1].set_ylabel("B / √k")
    axes[1].legend(fontsize=8)
    axes[1].set_title("√k-Skalierung")

    fig.tight_layout()
    png = out_dir / "eabc_cumulative_bias_B101_11.png"
    fig.savefig(png, dpi=120)
    plt.close(fig)
    return str(png)


def run(max_prime: int, n_null_trials: int, seed: int) -> str:
    out_dir = Path(__file__).resolve().parent
    report_lines: list[str] = []

    def log(s: str = "") -> None:
        print(s)
        report_lines.append(s)

    log("=" * 72)
    log("EABC kumulative Kanal-Bias-Analyse mod 420")
    log("=" * 72)
    log(f"Prim-Obergrenze: {max_prime:,}")
    log(f"Kanäle: {CHANNELS}")

    primes, channels, t_sieve, t_collect = collect_channel_events(max_prime)
    log(f"\nVierlingssuche: {t_sieve:.2f}s")
    log(f"Kanal-Extraktion gesamt: {t_collect:.2f}s")
    log(f"Gültige Kanal-Ereignisse K = {len(primes):,}")

    if len(primes) == 0:
        log("Keine Ereignisse — Abbruch.")
        return "\n".join(report_lines)

    cum = cumulative_counts(channels)
    final_counts = cum[-1].tolist()
    log(f"Endzählungen N_i(K): {final_counts}")
    log(f"Summe: {sum(final_counts)}")

    k = np.arange(1, len(primes) + 1, dtype=np.int64)
    sqrt_k = np.sqrt(k.astype(float))
    b_101_11 = cum[:, IDX_101] - cum[:, IDX_11]
    ratio = b_101_11 / sqrt_k

    # --- B_{101,11} Verlauf ---
    log("\n--- B_{101,11}(k) ---")
    log(f"  Start: B(1) = {b_101_11[0]:+d}")
    log(f"  Ende:  B(K) = {b_101_11[-1]:+d}")
    log(f"  min/max: {b_101_11.min():+d} / {b_101_11.max():+d}")

    milestones = [10, 50, 100, 500, 1000, 2000, 3000, 4000, len(primes)]
    log("\n  Meilensteine (k, p, B, B/√k):")
    for m in milestones:
        if m > len(primes):
            continue
        i = m - 1
        log(
            f"    k={m:5d}  p={primes[i]:12,d}  "
            f"B={b_101_11[i]:+5d}  B/√k={ratio[i]:+.4f}"
        )

    log(f"\n  ASCII B/√k: {ascii_sparkline(ratio)}")

    # --- Vorzeichenwechsel ---
    n_sc, sc_idx = sign_changes(b_101_11)
    log("\n--- Vorzeichenwechsel B_{101,11} ---")
    log(f"  Anzahl: {n_sc}")
    if sc_idx:
        show = sc_idx[:15]
        details = ", ".join(
            f"k={idx} (B={b_101_11[idx-2]:+d}→{b_101_11[idx-1]:+d})"
            for idx in show
            if idx >= 2
        )
        log(f"  Erste Wechsel: {details}")
        if len(sc_idx) > 15:
            log(f"  ... und {len(sc_idx) - 15} weitere")

    # --- Drift ---
    log("\n--- Drift / linearer Trend ---")
    a_r, b_r, r2_r = linear_drift(ratio, k.astype(float))
    log(f"  B/√k ~ {a_r:.6e}·k + {b_r:.6f}  (R²={r2_r:.6f})")

    a_b, b_b, r2_b = linear_drift(b_101_11.astype(float), sqrt_k)
    log(f"  B ~ {a_b:.6f}·√k + {b_b:.4f}  (R²={r2_b:.6f})  [erwartet ~ konstant×√k bei RW]")

    # Drift in letzter Hälfte
    half = len(k) // 2
    a2, b2, r2_2 = linear_drift(ratio[half:], k[half:].astype(float))
    log(f"  B/√k (zweite Hälfte): Steigung {a2:.6e}, R²={r2_2:.6f}")

    # --- √x-Skalierung & LIL ---
    max_abs_ratio = float(np.max(np.abs(ratio)))
    k_at_max = int(k[np.argmax(np.abs(ratio))])
    lil_at_k = lil_bound(len(primes))
    lil_ratio = abs(b_101_11[-1]) / lil_at_k

    log("\n--- √k-Skalierung & LIL ---")
    log(f"  max |B|/√k = {max_abs_ratio:.4f}  bei k={k_at_max}")
    log(f"  final B/√K = {ratio[-1]:+.4f}")
    log(f"  LIL-Grenze √(2K log log K) ≈ {lil_at_k:.2f}")
    log(f"  |B(K)| / LIL ≈ {lil_ratio:.4f}")
    log(f"  max |B| / LIL ≈ {abs(b_101_11).max() / lil_at_k:.4f}")

    bounded = max_abs_ratio < 5.0  # grobe Heuristik
    log(f"  B/√k bounded (max < 5): {'ja' if bounded else 'nein'}")

    # --- Alle Paare (Endwerte) ---
    log("\n--- B_{{j,i}}(K) für alle Kanalpaare ---")
    pairs = all_pair_biases(cum)
    for (cj, ci), series in sorted(pairs.items(), key=lambda x: -abs(x[1][-1])):
        log(f"  B_{{{cj},{ci}}} = {series[-1]:+d}  (B/√K = {series[-1]/sqrt_k[-1]:+.4f})")

    # --- Random-Walk-Null ---
    log(f"\n--- Nullhypothese: symmetrischer 6-Kanal-Random-Walk ({n_null_trials} Trials) ---")
    rng = np.random.default_rng(seed)
    null = random_walk_null(len(primes), n_null_trials, rng)

    obs_max = max_abs_ratio
    obs_final = ratio[-1]
    obs_sc = n_sc

    pct_max = 100.0 * np.mean(null["max_abs_ratio"] >= obs_max)
    pct_final = 100.0 * np.mean(np.abs(null["final_ratio"]) >= abs(obs_final))
    pct_sc = 100.0 * np.mean(null["sign_changes"] >= obs_sc)

    log(f"  Beobachtet max|B|/√k = {obs_max:.4f}")
    log(
        f"  Null: max|B|/√k  Mittel={null['max_abs_ratio'].mean():.4f}, "
        f"95%-Quantil={np.percentile(null['max_abs_ratio'], 95):.4f}, "
        f"99%-Quantil={np.percentile(null['max_abs_ratio'], 99):.4f}"
    )
    log(f"  p-Wert (max ≥ beobachtet): {pct_max:.1f}%")

    log(f"  Beobachtet |B/√K| = {abs(obs_final):.4f}")
    log(
        f"  Null |B/√K|: Mittel={np.abs(null['final_ratio']).mean():.4f}, "
        f"95%-Quantil={np.percentile(np.abs(null['final_ratio']), 95):.4f}"
    )
    log(f"  p-Wert (|final| ≥ beobachtet): {pct_final:.1f}%")

    log(f"  Beobachtet Vorzeichenwechsel: {obs_sc}")
    log(
        f"  Null Vorzeichenwechsel: Mittel={null['sign_changes'].mean():.1f}, "
        f"Median={np.median(null['sign_changes']):.0f}"
    )
    log(f"  p-Wert (Wechsel ≥ beobachtet): {pct_sc:.1f}%")

    # --- Speichern ---
    csv_path = out_dir / "eabc_cumulative_bias.csv"
    save_csv(csv_path, primes, cum, b_101_11)
    log(f"\nCSV: {csv_path}")

    plot_msg = maybe_plot(out_dir, k, b_101_11, ratio)
    log(f"Plot: {plot_msg}")

    report_path = out_dir / "eabc_cumulative_bias_report.txt"
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    log(f"Report: {report_path}")

    log("\n" + "=" * 72)
    log("Fertig.")
    log("=" * 72)

    return "\n".join(report_lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Kumulative mod-420 Kanal-Bias-Analyse")
    parser.add_argument(
        "--max-prime",
        type=int,
        default=100_000_000,
        help="Prim-Obergrenze für Vierlingssuche (Default: 10^8)",
    )
    parser.add_argument(
        "--null-trials",
        type=int,
        default=5000,
        help="Monte-Carlo-Trials für Random-Walk-Null",
    )
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    run(args.max_prime, args.null_trials, args.seed)


if __name__ == "__main__":
    main()
