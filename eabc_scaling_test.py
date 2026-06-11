#!/usr/bin/env python3
"""
Skalierungstest für Chebyshev-artigen Kanal-Bias B_{101,11} mod 420.

Zählt Vierlinge (p, p+2, p+6, p+8) bis zu mehreren Prim-Obergrenzen und
vergleicht, ob B/√K stabil bleibt (√x-Drift) oder gegen 0 kollabiert.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eabc_cumulative_bias import (  # noqa: E402
    CHANNELS,
    IDX_101,
    IDX_11,
    collect_channel_events,
    cumulative_counts,
    random_walk_null,
    sign_changes,
)


DEFAULT_LIMITS = [100_000_000, 300_000_000, 1_000_000_000]
TIME_BUDGET_SEC = 30 * 60  # 30 min pro Grenze


@dataclass
class ScalingResult:
    n_max: int
    runtime_sieve_s: float
    runtime_total_s: float
    k: int
    n_11: int
    n_101: int
    n_191: int
    n_221: int
    n_311: int
    n_401: int
    b_101_11_final: int
    b_over_sqrt_k: float
    max_abs_b_over_sqrt_k: float
    k_at_max_ratio: int
    frac_b_positive: float
    sign_changes: int
    max_pair_abs_ratio: float
    max_pair_label: str
    null_p_final: float
    null_p_max: float
    extrapolated: bool = False


def analyze_limit(n_max: int, n_null_trials: int, seed: int) -> ScalingResult:
    t0 = time.perf_counter()
    primes, channels, t_sieve, t_collect = collect_channel_events(n_max)
    runtime_total = time.perf_counter() - t0

    k = len(primes)
    if k == 0:
        return ScalingResult(
            n_max=n_max,
            runtime_sieve_s=t_sieve,
            runtime_total_s=runtime_total,
            k=0,
            n_11=0,
            n_101=0,
            n_191=0,
            n_221=0,
            n_311=0,
            n_401=0,
            b_101_11_final=0,
            b_over_sqrt_k=0.0,
            max_abs_b_over_sqrt_k=0.0,
            k_at_max_ratio=0,
            frac_b_positive=0.0,
            sign_changes=0,
            max_pair_abs_ratio=0.0,
            max_pair_label="",
            null_p_final=1.0,
            null_p_max=1.0,
        )

    cum = cumulative_counts(channels)
    final = cum[-1]
    sqrt_k = np.sqrt(np.arange(1, k + 1, dtype=float))
    b = cum[:, IDX_101] - cum[:, IDX_11]
    ratio = b / sqrt_k

    n_sc, _ = sign_changes(b)
    frac_pos = float(np.mean(b > 0))

    # max |B_{j,i}|/√k über alle Paare
    max_pair_ratio = 0.0
    max_pair_label = ""
    n_ch = len(CHANNELS)
    for i in range(n_ch):
        for j in range(i + 1, n_ch):
            pair_b = cum[:, j] - cum[:, i]
            pair_ratio = pair_b / sqrt_k
            m = float(np.max(np.abs(pair_ratio)))
            if m > max_pair_ratio:
                max_pair_ratio = m
                max_pair_label = f"B_{{{CHANNELS[j]},{CHANNELS[i]}}}"

    rng = np.random.default_rng(seed + n_max)
    null = random_walk_null(k, n_null_trials, rng)
    obs_max = float(np.max(np.abs(ratio)))
    obs_final = float(ratio[-1])
    null_p_max = 100.0 * float(np.mean(null["max_abs_ratio"] >= obs_max))
    null_p_final = 100.0 * float(np.mean(np.abs(null["final_ratio"]) >= abs(obs_final)))

    idx_max = int(np.argmax(np.abs(ratio)))

    return ScalingResult(
        n_max=n_max,
        runtime_sieve_s=t_sieve,
        runtime_total_s=t_collect,
        k=k,
        n_11=int(final[0]),
        n_101=int(final[1]),
        n_191=int(final[2]),
        n_221=int(final[3]),
        n_311=int(final[4]),
        n_401=int(final[5]),
        b_101_11_final=int(b[-1]),
        b_over_sqrt_k=float(ratio[-1]),
        max_abs_b_over_sqrt_k=obs_max,
        k_at_max_ratio=idx_max + 1,
        frac_b_positive=frac_pos,
        sign_changes=n_sc,
        max_pair_abs_ratio=max_pair_ratio,
        max_pair_label=max_pair_label,
        null_p_final=null_p_final,
        null_p_max=null_p_max,
    )


def extrapolate_from(results: list[ScalingResult], target: int) -> ScalingResult:
    """Lineare Extrapolation K ~ c·N_max und B ~ d·√K."""
    if len(results) < 2:
        raise ValueError("Mindestens zwei Messpunkte für Extrapolation nötig")

    xs = np.array([r.n_max for r in results], dtype=float)
    ks = np.array([r.k for r in results], dtype=float)
    bs = np.array([r.b_101_11_final for r in results], dtype=float)

    # K linear in N_max (Primzahltheorem-artig für Vierlinge)
    k_slope = np.polyfit(xs, ks, 1)
    k_est = float(np.polyval(k_slope, target))

    # B ~ c·√K wenn √x-Skalierung; sonst B/K → 0
    sqrt_ks = np.sqrt(ks)
    b_over_sqrt = bs / sqrt_ks
    mean_ratio = float(np.mean(b_over_sqrt))
    b_est = mean_ratio * math.sqrt(k_est)

    last = results[-1]
    runtime_per = last.runtime_total_s / last.n_max

    return ScalingResult(
        n_max=target,
        runtime_sieve_s=runtime_per * target,
        runtime_total_s=runtime_per * target,
        k=int(round(k_est)),
        n_11=0,
        n_101=0,
        n_191=0,
        n_221=0,
        n_311=0,
        n_401=0,
        b_101_11_final=int(round(b_est)),
        b_over_sqrt_k=mean_ratio,
        max_abs_b_over_sqrt_k=mean_ratio * 1.1,
        k_at_max_ratio=0,
        frac_b_positive=last.frac_b_positive,
        sign_changes=0,
        max_pair_abs_ratio=0.0,
        max_pair_label="extrapoliert",
        null_p_final=float("nan"),
        null_p_max=float("nan"),
        extrapolated=True,
    )


def write_csv(path: Path, results: list[ScalingResult]) -> None:
    fields = list(asdict(results[0]).keys()) if results else []
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in results:
            w.writerow(asdict(r))


def format_table(results: list[ScalingResult]) -> str:
    hdr = (
        f"{'N_max':>14} {'K':>8} {'B_101,11':>9} {'B/√K':>8} "
        f"{'max|B|/√k':>10} {'δ(B>0)':>8} {'VZW':>5} {'Zeit[s]':>9} {'extr.':>6}"
    )
    lines = [hdr, "-" * len(hdr)]
    for r in results:
        ext = "ja" if r.extrapolated else "nein"
        lines.append(
            f"{r.n_max:>14,} {r.k:>8,} {r.b_101_11_final:>+9d} {r.b_over_sqrt_k:>+8.4f} "
            f"{r.max_abs_b_over_sqrt_k:>10.4f} {r.frac_b_positive:>8.4f} {r.sign_changes:>5d} "
            f"{r.runtime_total_s:>9.2f} {ext:>6}"
        )
    return "\n".join(lines)


def scaling_verdict(results: list[ScalingResult]) -> str:
    measured = [r for r in results if not r.extrapolated]
    if len(measured) < 2:
        return "Zu wenige Messpunkte für Skalierungsurteil."

    ratios = [r.b_over_sqrt_k for r in measured]
    ks = [r.k for r in measured]
    bs = [r.b_101_11_final for r in measured]

    # Trend B/√K
    slope_ratio = (ratios[-1] - ratios[0]) / max(len(ratios) - 1, 1)
    collapse = all(abs(r) < 0.3 for r in ratios[1:]) and abs(ratios[-1]) < abs(ratios[0]) * 0.5

    # B vs √K linear fit
    sqrt_ks = np.sqrt(np.array(ks, dtype=float))
    if len(sqrt_ks) >= 2:
        coef = np.polyfit(sqrt_ks, bs, 1)
        b_pred = np.polyval(coef, sqrt_ks)
        ss_res = np.sum((bs - b_pred) ** 2)
        ss_tot = np.sum((np.array(bs) - np.mean(bs)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    else:
        r2 = 0.0
        coef = [0.0, 0.0]

    lines = []
    if collapse:
        lines.append(
            "KOLLAPSE: B/√K fällt mit wachsendem K — kein persistenter Chebyshev-Drift."
        )
    elif abs(ratios[-1]) > 0.5 and abs(slope_ratio) < 0.15:
        lines.append(
            f"STABIL/√x: B/√K bleibt O(1) (Endwerte {ratios[0]:+.3f} → {ratios[-1]:+.3f}). "
            f"Linearer Fit B ~ {coef[0]:.3f}·√K + {coef[1]:.1f} (R²={r2:.3f}) "
            "spricht für Random-Walk-/Dirichlet-artige Fluktuation, nicht für Kollaps."
        )
    else:
        lines.append(
            f"GEMISCHT: B/√K von {ratios[0]:+.3f} nach {ratios[-1]:+.3f}; "
            f"Trend pro Stufe {slope_ratio:+.4f}."
        )

    # Null-p-Werte
    for r in measured:
        lines.append(
            f"  N_max={r.n_max:,}: p(|B/√K|≥beob.)={r.null_p_final:.1f}%, "
            f"p(max|B|/√k≥beob.)={r.null_p_max:.1f}%"
        )

    return "\n".join(lines)


def run(limits: list[int], n_null_trials: int, seed: int, time_budget: float) -> str:
    out_dir = Path(__file__).resolve().parent
    report: list[str] = []
    results: list[ScalingResult] = []

    def log(s: str = "") -> None:
        print(s, flush=True)
        report.append(s)

    log("=" * 72)
    log("EABC Skalierungstest — Chebyshev-artiger Bias B_{101,11}")
    log("=" * 72)
    log(f"Grenzen: {[f'{x:,}' for x in limits]}")
    log(f"Kanäle mod 420: {CHANNELS}")
    log(f"Random-Walk-Null: {n_null_trials} Trials pro Grenze")
    log(f"Zeitbudget pro Grenze: {time_budget/60:.0f} min")
    log("")

    for n_max in limits:
        log(f"--- N_max = {n_max:,} ---")
        t_start = time.perf_counter()
        res = analyze_limit(n_max, n_null_trials, seed)
        elapsed = time.perf_counter() - t_start
        results.append(res)

        log(f"  Laufzeit: {res.runtime_total_s:.2f}s (Sieb: {res.runtime_sieve_s:.2f}s)")
        log(f"  K = {res.k:,}")
        log(
            f"  N_i = [{res.n_11}, {res.n_101}, {res.n_191}, "
            f"{res.n_221}, {res.n_311}, {res.n_401}]"
        )
        log(f"  B_101,11 = {res.b_101_11_final:+d}, B/√K = {res.b_over_sqrt_k:+.4f}")
        log(f"  max|B|/√k = {res.max_abs_b_over_sqrt_k:.4f} bei k={res.k_at_max_ratio}")
        log(f"  δ(B>0) = {res.frac_b_positive:.4f}, Vorzeichenwechsel = {res.sign_changes}")
        log(
            f"  max Paar: {res.max_pair_label} mit max|B|/√k = {res.max_pair_abs_ratio:.4f}"
        )
        log(
            f"  Null p-Werte: |final|={res.null_p_final:.1f}%, max={res.null_p_max:.1f}%"
        )
        log("")

        # Nächste Grenze überspringen/extrapolieren wenn Budget überschritten
        remaining = [x for x in limits if x > n_max and x not in [r.n_max for r in results]]
        if remaining and elapsed > time_budget:
            log(f"  WARNUNG: Laufzeit {elapsed:.1f}s > Budget — extrapoliere restliche Grenzen.")
            for target in remaining:
                ext = extrapolate_from(results, target)
                results.append(ext)
                log(f"  [EXTRAPOLIERT] N_max={target:,}: K≈{ext.k:,}, B/√K≈{ext.b_over_sqrt_k:+.4f}")
            break

        # Extrapolation 10^9 wenn geschätzte Zeit zu lang
        if n_max == 300_000_000 and 1_000_000_000 in limits:
            rate = res.runtime_total_s / n_max
            est_1e9 = rate * 1_000_000_000
            log(f"  Geschätzte Zeit für 10^9: {est_1e9:.1f}s ({est_1e9/60:.1f} min)")
            if est_1e9 > time_budget and 1_000_000_000 not in [r.n_max for r in results]:
                ext = extrapolate_from(results, 1_000_000_000)
                results.append(ext)
                limits_done = {r.n_max for r in results}
                if 1_000_000_000 in [x for x in limits if x not in limits_done]:
                    pass
                log(
                    f"  [EXTRAPOLIERT] N_max=10^9: K≈{ext.k:,}, "
                    f"B/√K≈{ext.b_over_sqrt_k:+.4f} (geschätzt {est_1e9/60:.1f} min)"
                )
                # Entferne 10^9 aus weiterer Ausführung
                limits = [x for x in limits if x != 1_000_000_000 or ext.n_max == 1_000_000_000]

    log("\n" + format_table(results))
    log("\n--- Skalierungsurteil ---")
    log(scaling_verdict(results))

    # Random-Walk-Vergleich bei größtem K
    measured = [r for r in results if not r.extrapolated]
    if measured:
        best = max(measured, key=lambda r: r.k)
        log(f"\n--- Random-Walk bei K={best.k:,} ({best.n_max:,}) ---")
        log(f"  Beobachtet: B/√K={best.b_over_sqrt_k:+.4f}, max|B|/√k={best.max_abs_b_over_sqrt_k:.4f}")
        log(f"  Null p(|final|): {best.null_p_final:.1f}%, p(max): {best.null_p_max:.1f}%")

    csv_path = out_dir / "eabc_scaling_report.csv"
    write_csv(csv_path, results)
    log(f"\nCSV: {csv_path}")

    report_path = out_dir / "eabc_scaling_report.txt"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    log(f"Report: {report_path}")

    return "\n".join(report)


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC Skalierungstest B_{101,11}")
    parser.add_argument(
        "--limits",
        type=int,
        nargs="+",
        default=DEFAULT_LIMITS,
        help="Prim-Obergrenzen (Default: 10^8 3·10^8 10^9)",
    )
    parser.add_argument("--null-trials", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--time-budget-min",
        type=float,
        default=30.0,
        help="Max. Minuten pro Grenze",
    )
    args = parser.parse_args()
    run(args.limits, args.null_trials, args.seed, args.time_budget_min * 60)


if __name__ == "__main__":
    main()
