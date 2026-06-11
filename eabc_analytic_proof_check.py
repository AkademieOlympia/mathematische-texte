#!/usr/bin/env python3
"""
Numerische Verifikation der analytischen Skalierung
  ⟨δ_K, χ₃⟩ ∼ N(0, 1/(36K))   bzw.   |⟨δ_K, χ₃⟩| = O_p(1/√K)

Vergleicht beobachtete χ₃-Projektionen bei K = 4767, 10972, 28387
mit der Multinomial-Null und dem Hardy-Littlewood-artigen Referenzmodell.
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
from eabc_cumulative_bias import CHANNELS, collect_channel_events  # noqa: E402
from Stiefel import INDEX  # noqa: E402

R_VALS = np.array(CHANNELS)
OMEGA = np.exp(2j * np.pi / 3)
CHI3_TABLE = {1: 1.0, 2: OMEGA**2, 3: OMEGA, 4: 1.0, 5: OMEGA, 6: OMEGA**2}
CHI3_R = np.array([CHI3_TABLE[int(r % 7)] for r in R_VALS], dtype=complex)
M = 6  # Kanäle


def inner_product_delta(counts: np.ndarray, chi: np.ndarray) -> complex:
    """⟨δ, χ⟩ mit δ_r = n_r/K - 1/6."""
    k = counts.sum()
    delta = counts / k - 1.0 / M
    return np.sum(delta * np.conj(chi)) / M


def theoretical_std(k: int) -> float:
    """Std(⟨δ_K, χ₃⟩) unter H (Multinomial 1/6)."""
    return 1.0 / (6.0 * math.sqrt(k))


def theoretical_var(k: int) -> float:
    return 1.0 / (36.0 * k)


def z_score(obs: complex, k: int) -> float:
    return abs(obs) / theoretical_std(k)


@dataclass
class Chi3ScalingRow:
    n_max: int
    k: int
    counts: list[int]
    delta_prob: list[float]
    chi3_proj: complex
    chi3_abs: float
    chi3_arg_deg: float
    var_chi3_pair_pct: float
    std_theory: float
    z_score: float
    ratio_obs_theory: float
    b_101_11: int
    b_over_sqrt_k: float
    runtime_s: float


def analyze_limit(n_max: int) -> Chi3ScalingRow:
    t0 = time.perf_counter()
    _, channels, _, _ = collect_channel_events(n_max)
    k = len(channels)
    counts = np.zeros(M, dtype=np.int64)
    for c in channels:
        counts[INDEX[int(c)]] += 1

    chi3 = inner_product_delta(counts.astype(float), CHI3_R)
    delta = counts.astype(float) / k - 1.0 / M
    var_total = float(np.sum(delta**2) / M)
    var_chi3 = 2 * abs(chi3) ** 2  # χ₃ + χ̄₃ Paar
    var_pct = 100.0 * var_chi3 / var_total if var_total > 0 else 0.0

    idx_11 = CHANNELS.index(11)
    idx_101 = CHANNELS.index(101)
    b = int(counts[idx_101] - counts[idx_11])
    b_sqrt = b / math.sqrt(k)
    std_th = theoretical_std(k)

    return Chi3ScalingRow(
        n_max=n_max,
        k=k,
        counts=[int(x) for x in counts],
        delta_prob=[float(x) for x in counts / k - 1.0 / M],
        chi3_proj=chi3,
        chi3_abs=abs(chi3),
        chi3_arg_deg=float(np.degrees(np.angle(chi3))),
        var_chi3_pair_pct=var_pct,
        std_theory=std_th,
        z_score=z_score(chi3, k),
        ratio_obs_theory=abs(chi3) / std_th if std_th > 0 else 0.0,
        b_101_11=b,
        b_over_sqrt_k=b_sqrt,
        runtime_s=time.perf_counter() - t0,
    )


def multinomial_null_sim(
    k: int, n_trials: int, seed: int
) -> dict[str, np.ndarray]:
    """Null: Multinomial(K, (1/6,...,1/6)), ⟨δ,χ₃⟩."""
    rng = np.random.default_rng(seed)
    probs = np.full(M, 1.0 / M)
    abs_proj = np.empty(n_trials)
    for t in range(n_trials):
        counts = rng.multinomial(k, probs)
        abs_proj[t] = abs(inner_product_delta(counts.astype(float), CHI3_R))
    return {"abs_proj": abs_proj}


def format_report(rows: list[Chi3ScalingRow], null_stats: dict[int, dict]) -> str:
    lines = [
        "=" * 72,
        "EABC Analytischer Beweisversuch — χ₃-Skalierung",
        "=" * 72,
        "",
        "Theorie (Hypothese H, Multinomial 1/6):",
        "  E[⟨δ_K,χ₃⟩] = 0",
        "  Var(⟨δ_K,χ₃⟩) = 1/(36K)",
        "  Std(|⟨δ_K,χ₃⟩|) ≈ √(π/(72K))  (komplexe Normalverteilung)",
        "  Typische Größenordnung: |⟨δ_K,χ₃⟩| = O_p(1/√K)",
        "",
        f"{'N_max':>14} {'K':>8} {'|⟨δ,χ₃⟩|':>12} {'1/(6√K)':>10} "
        f"{'z':>6} {'Var%χ₃':>8} {'B/√K':>8}",
        "-" * 72,
    ]
    for r in rows:
        lines.append(
            f"{r.n_max:>14,} {r.k:>8,} {r.chi3_abs:>12.6f} {r.std_theory:>10.6f} "
            f"{r.z_score:>6.2f} {r.var_chi3_pair_pct:>8.2f} {r.b_over_sqrt_k:>+8.4f}"
        )

    lines.append("")
    lines.append("--- Skalierungsfit: |⟨δ,χ₃⟩| · √K ---")
    for r in rows:
        lines.append(f"  K={r.k:,}: |⟨δ,χ₃⟩|·√K = {r.chi3_abs * math.sqrt(r.k):.4f}")

    lines.append("")
    lines.append("--- Multinomial-Null (|⟨δ,χ₃⟩|) ---")
    for k, stats in null_stats.items():
        obs = next(x.chi3_abs for x in rows if x.k == k)
        p = 100.0 * float(np.mean(stats["abs_proj"] >= obs))
        lines.append(
            f"  K={k:,}: beob.={obs:.6f}, Null-Median={stats['median']:.6f}, "
            f"p(≥beob.)={p:.1f}%"
        )

    lines.append("")
    lines.append("--- Urteil ---")
    scaled = [r.chi3_abs * math.sqrt(r.k) for r in rows]
    mean_scaled = float(np.mean(scaled))
    std_scaled = float(np.std(scaled))
    if std_scaled / max(mean_scaled, 1e-9) < 0.35:
        lines.append(
            f"KONSISTENT: |⟨δ,χ₃⟩|·√K stabil ({mean_scaled:.3f} ± {std_scaled:.3f}) "
            f"— passt zu f(K) = Θ(1/√K)."
        )
    else:
        lines.append(
            f"GEMISCHT: |⟨δ,χ₃⟩|·√K variiert ({mean_scaled:.3f} ± {std_scaled:.3f})."
        )

    b_ratios = [r.b_over_sqrt_k for r in rows]
    lines.append(
        f"B/√K fällt ({b_ratios[0]:+.3f} → {b_ratios[-1]:+.3f}): "
        "erwartbar unter Random-Walk-Null (ein Pfad, kein Limes)."
    )

    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--limits",
        type=int,
        nargs="+",
        default=[100_000_000, 300_000_000, 1_000_000_000],
    )
    parser.add_argument("--null-trials", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    rows = [analyze_limit(n) for n in args.limits]

    null_stats: dict[int, dict] = {}
    for r in rows:
        sim = multinomial_null_sim(r.k, args.null_trials, args.seed + r.k)
        null_stats[r.k] = {
            "abs_proj": sim["abs_proj"],
            "median": float(np.median(sim["abs_proj"])),
        }

    report = format_report(rows, null_stats)
    print(report)

    csv_path = out_dir / "eabc_analytic_proof_check.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        fields = list(asdict(rows[0]).keys())
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            d = asdict(r)
            d["chi3_proj"] = str(d["chi3_proj"])
            w.writerow(d)

    report_path = out_dir / "eabc_analytic_proof_check_report.txt"
    report_path.write_text(report, encoding="utf-8")
    print(f"CSV: {csv_path}")
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
