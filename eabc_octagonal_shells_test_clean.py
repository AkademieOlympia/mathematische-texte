#!/usr/bin/env python3
"""
Minimal reproduzierbares Numerik-Experiment:

"Division" von Prime-Vierlingen (p, p+2, p+6, p+8) in einer oktagonalen Darstellung,
mit ABCE vs CEAB als zwei "Hurwitz half-shells" (operational),
anschließend Test auf "layers of nine" und E-Achsen-Analyse.

Da im Repo (bei den relevanten Suchen) keine expliziten Implementierungen für
octagonal representation / Hurwitz half-shells / R8-balls / E-axis / layers of nine
vorliegen, definiert dieses Skript eine *klare* und *deterministische* Abbildung.

Operationalisierung "layers of nine":
Für jeden Vierling erzeugt das Skript 9 Punkte (8 Oktagonpunkte + Mittelpunkt),
zentrert die 8 Oktagonpunkte auf Ursprung (Mittelpunkt -> (0,0,0)) und zählt
diskrete Norm-Schichten als Anzahl der unterschiedlichen Normabstände unter den 9 Punkten
(gerundet auf feste Dezimalstellen). "layers of nine" heißt hier: genau 9 Schichten.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


CHANNELS_MOD420 = [11, 101, 191, 221, 311, 401]
CHANNEL_TO_INDEX = {r: i for i, r in enumerate(CHANNELS_MOD420)}


def sieve_bool(n: int) -> np.ndarray:
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for k in range(2, int(math.isqrt(n)) + 1):
        if sieve[k]:
            sieve[k * k : n + 1 : k] = False
    return sieve


def find_prime_quadruplets(max_p: int) -> list[int]:
    sieve = sieve_bool(max_p + 8)
    out: list[int] = []
    for p in range(5, max_p + 1, 2):
        if sieve[p] and sieve[p + 2] and sieve[p + 6] and sieve[p + 8]:
            out.append(p)
    return out


def chirality_from_p(p: int) -> str | None:
    r12 = p % 12
    if r12 == 5:
        return "ABCE"
    if r12 == 11:
        return "CEAB"
    return None


def channel_index_mod420(p: int) -> int | None:
    return CHANNEL_TO_INDEX.get(p % 420, None)


def quad_octagonal_points_centered(p: int, channel_idx: int, chirality: str) -> np.ndarray:
    """
    9 Punkte in 3D:
    - 8 Oktagonpunkte auf Winkel k*pi/4
    - 1 Mittelpunktpunkt = centroid der 8 Punkte

    ABCE => z = +radius, CEAB => z = -radius ("up/down").
    Zwei "half-shells" (oben/unten) werden über half=k//4 eingeführt:
      radius = log(p+offset) + half*log(channel_idx+1)
    """
    radii4 = np.log(np.array([p, p + 2, p + 6, p + 8], dtype=float))
    boost = math.log(channel_idx + 1.0)
    sign = 1.0 if chirality == "ABCE" else -1.0
    theta_offset = math.pi / 4 if chirality == "ABCE" else -math.pi / 4

    pts8 = []
    base_angles = np.arange(8, dtype=float) * (math.pi / 4) + theta_offset
    for k in range(8):
        idx = k % 4
        half = k // 4
        r = radii4[idx] + half * boost
        x = r * math.cos(base_angles[k])
        y = r * math.sin(base_angles[k])
        z = sign * r
        pts8.append((x, y, z))
    pts8 = np.asarray(pts8, dtype=float)  # (8,3)
    center = pts8.mean(axis=0, keepdims=True)  # (1,3)
    pts9_centered = np.vstack([pts8 - center, np.zeros((1, 3), dtype=float)])
    return pts9_centered


def unique_shell_count(pts9: np.ndarray, *, decimals: int) -> tuple[int, np.ndarray]:
    norms = np.linalg.norm(pts9, axis=1)
    norms_r = np.round(norms, decimals=decimals)
    uniq = np.unique(norms_r)
    return int(len(uniq)), uniq


def dft_magnitudes(counts_by_channel: np.ndarray) -> list[float]:
    v = np.asarray(counts_by_channel, dtype=float)
    if v.sum() > 0:
        v = v / v.sum()
    F = np.fft.fft(v)
    return [float(abs(z)) for z in F]


def z_peaks(z: np.ndarray, *, decimals: int) -> dict[str, int]:
    z_r = np.round(z, decimals=decimals)
    u, c = np.unique(z_r, return_counts=True)
    # sort descending by count
    order = np.argsort(-c)
    return {str(float(u[i])): int(c[i]) for i in order[:15]}


@dataclass
class Trial:
    trial_index: int
    sample_size: int
    always_all9: bool
    frac_quads_shell9: float
    min_shells: int
    max_shells: int
    mean_shells: float
    shell_counts: list[int]


def wilson_interval(k: int, n: int, *, alpha: float = 0.05) -> dict[str, float]:
    if n == 0:
        return {"p_hat": float("nan"), "lo": float("nan"), "hi": float("nan")}
    # z_{1-alpha/2} for alpha=0.05
    z = 1.959963984540054
    p_hat = k / n
    denom = 1 + z * z / n
    center = p_hat + z * z / (2 * n)
    half = z * math.sqrt((p_hat * (1 - p_hat) + z * z / (4 * n)) / n)
    lo = (center - half) / denom
    hi = (center + half) / denom
    return {"p_hat": p_hat, "lo": lo, "hi": hi}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max_p", type=int, default=1_000_000, help="Suche bis p<=max_p")
    ap.add_argument("--sample_size", type=int, default=120, help="Vierlinge pro Trial")
    ap.add_argument("--trials", type=int, default=60, help="Anzahl random Trials")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--shell_decimals", type=int, default=7, help="Rundung für Norm-Schichten")
    ap.add_argument("--eaxis_decimals", type=int, default=6, help="Rundung für E-Achsen-Peaks")
    ap.add_argument("--report_name", type=str, default="eabc_octagonal_shells_test_clean")
    args = ap.parse_args()

    out_dir = Path(__file__).resolve().parent
    report_txt = out_dir / f"{args.report_name}.txt"
    report_json = out_dir / f"{args.report_name}.json"

    t0 = time.time()
    print(f"max_p={args.max_p:,} ...")
    quadruplets = find_prime_quadruplets(args.max_p)

    # Filter: mod-420 Kanäle + Chiralität
    items: list[tuple[int, str, int]] = []
    for p in quadruplets:
        ch = chirality_from_p(p)
        if ch is None:
            continue
        cidx = channel_index_mod420(p)
        if cidx is None:
            continue
        items.append((p, ch, int(cidx)))

    n_items = len(items)
    if n_items == 0:
        raise SystemExit("Keine Vierlinge in mod-420 Kanälen gefunden. max_p erhöhen.")

    rng = np.random.default_rng(args.seed)
    sample_size = int(min(args.sample_size, n_items))
    trials = int(max(1, args.trials))

    print(f"mod-420 gefiltert: {n_items} Vierlinge; sample_size={sample_size}, trials={trials}")

    # Precompute shells==count for each quadruplet
    shells_counts = np.empty(n_items, dtype=int)
    shells_uniq_example: dict[int, list[float]] = {}
    for i, (p, ch, cidx) in enumerate(items):
        pts9 = quad_octagonal_points_centered(p, cidx, ch)
        cnt, uniq = unique_shell_count(pts9, decimals=args.shell_decimals)
        shells_counts[i] = cnt
        if i < 3:
            shells_uniq_example[i] = [float(x) for x in uniq.tolist()]

    n_shell9 = int(np.sum(shells_counts == 9))
    frac_shell9 = n_shell9 / n_items
    print(f"shelldebug: shells==9 for {n_shell9}/{n_items} = {frac_shell9:.4f}")

    # Trials
    trial_list: list[Trial] = []
    trials_all9 = 0
    quad_scored = 0
    quad_shell9 = 0

    for tr in range(trials):
        idxs = rng.choice(n_items, size=sample_size, replace=False) if sample_size < n_items else np.arange(n_items)
        sc = shells_counts[idxs]
        quad_scored += int(sc.size)
        quad_shell9 += int(np.sum(sc == 9))
        always_all9 = bool(np.all(sc == 9))
        if always_all9:
            trials_all9 += 1
        trial_list.append(
            Trial(
                trial_index=tr,
                sample_size=sample_size,
                always_all9=always_all9,
                frac_quads_shell9=float(np.mean(sc == 9)),
                min_shells=int(sc.min()),
                max_shells=int(sc.max()),
                mean_shells=float(sc.mean()),
                shell_counts=[int(x) for x in sc.tolist()],
            )
        )

    frac_trials_all9 = trials_all9 / trials
    wilson_quad = wilson_interval(quad_shell9, quad_scored)
    wilson_trials = wilson_interval(trials_all9, trials)

    # Transformation / E-axis analysis on a larger sample:
    take_n = int(min(max(200, sample_size), n_items))
    idxs_big = rng.choice(n_items, size=take_n, replace=False) if take_n < n_items else np.arange(n_items)

    pts_blocks_by_ch = {"ABCE": [], "CEAB": []}
    channel_counts_by_ch = {"ABCE": np.zeros(6, dtype=int), "CEAB": np.zeros(6, dtype=int)}

    for i in idxs_big:
        p, ch, cidx = items[i]
        pts9 = quad_octagonal_points_centered(p, cidx, ch)
        pts_blocks_by_ch[ch].append(pts9)
        channel_counts_by_ch[ch][cidx] += int(pts9.shape[0])  # 9 points per quad

    # "move midpoints of respective R8-balls to zero point":
    # operational: shift z so that mean z of the 8 outer points per chirality is 0.
    z_shift = {}
    for ch in ["ABCE", "CEAB"]:
        if not pts_blocks_by_ch[ch]:
            z_shift[ch] = 0.0
            continue
        blocks = np.array(pts_blocks_by_ch[ch], dtype=float)  # (n_quads,9,3)
        outer = blocks[:, :8, :].reshape(-1, 3)
        z_shift[ch] = float(np.mean(outer[:, 2]))

    eaxis_peaks = {}
    eaxis_summary = {}
    for ch in ["ABCE", "CEAB"]:
        if not pts_blocks_by_ch[ch]:
            continue
        blocks = np.array(pts_blocks_by_ch[ch], dtype=float)
        outer_z = blocks[:, :8, 2].reshape(-1) - z_shift[ch]
        eaxis_peaks[ch] = z_peaks(outer_z, decimals=args.eaxis_decimals)
        eaxis_summary[ch] = {
            "z_shift_used": z_shift[ch],
            "n_outer_points": int(outer_z.size),
            "z_mean_after_shift": float(np.mean(outer_z)),
            "z_std_after_shift": float(np.std(outer_z)),
        }

    # Fourier analog: ring C6 using counts_by_channel
    fourier_by_ch = {ch: dft_magnitudes(channel_counts_by_ch[ch]) for ch in ["ABCE", "CEAB"]}

    # Global shells visible (outer points) after per-quad centering (not using z_shift for norms)
    global_shells = {}
    for ch in ["ABCE", "CEAB"]:
        if not pts_blocks_by_ch[ch]:
            continue
        blocks = np.array(pts_blocks_by_ch[ch], dtype=float)
        outer_norms = np.linalg.norm(blocks[:, :8, :], axis=2).reshape(-1)
        shells = np.unique(np.round(outer_norms, decimals=args.shell_decimals))
        global_shells[ch] = shells.tolist()

    elapsed = time.time() - t0
    print(f"done in {elapsed:.2f}s")

    # Write report
    lines: list[str] = []
    lines.append("EABC Octagonal Shells Test (clean minimal operational definition)")
    lines.append("")
    lines.append(f"out_dir: {out_dir}")
    lines.append(f"max_p: {args.max_p:,}")
    lines.append(f"quadruplets found total: {len(quadruplets):,}")
    lines.append(f"mod-420 filtered + chirality: {n_items:,}")
    lines.append(f"sample_size per trial: {sample_size:,} ; trials: {trials}")
    lines.append("")
    lines.append("layers-of-nine (operational):")
    lines.append(f"  count distinct rounded norms among 9 points after centering -> equals 9")
    lines.append(f"shell_decimals: {args.shell_decimals}")
    lines.append("")
    lines.append(f"shells==9 per quadruplet: {n_shell9}/{n_items} = {frac_shell9:.4f}")
    lines.append(f"Wilson (quad-level): {wilson_quad}")
    lines.append(f"Trials with 'all quads shell==9': {trials_all9}/{trials} = {frac_trials_all9:.4f}")
    lines.append(f"Wilson (trial-level): {wilson_trials}")
    lines.append("")
    lines.append("example unique shell sets (first 3 quads):")
    for k, v in shells_uniq_example.items():
        lines.append(f"  idx={k}: shells uniq (rounded) = {v}")
    lines.append("")
    lines.append("E-axis analysis (z-axis), after chirality-specific z-shift:")
    lines.append(f"eaxis_decimals: {args.eaxis_decimals}")
    for ch in ["ABCE", "CEAB"]:
        if ch not in eaxis_summary:
            continue
        s = eaxis_summary[ch]
        lines.append(
            f"  {ch}: z_shift={s['z_shift_used']:.6f}, mean(z)={s['z_mean_after_shift']:.6f}, std(z)={s['z_std_after_shift']:.6f}, n_outer={s['n_outer_points']}"
        )
        lines.append(f"    top z-peaks: {eaxis_peaks[ch]}")
    lines.append("")
    lines.append("mod-420 ring-Fourier analog using point-counts per channel (C6, k=0..5):")
    for ch in ["ABCE", "CEAB"]:
        lines.append(f"  {ch}: magnitudes k=0..5 = {[f'{m:.6f}' for m in fourier_by_ch[ch]]}")
    lines.append("")
    lines.append("Global visible norm-shells among outer points (after per-quad centering):")
    for ch in ["ABCE", "CEAB"]:
        if ch not in global_shells:
            continue
        lines.append(f"  {ch}: #unique_norm_shells={len(global_shells[ch])}, first={global_shells[ch][:10]}")
    lines.append("")
    lines.append(f"Runtime: {elapsed:.2f}s")

    report_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")

    payload = {
        "settings": asdict(
            {
                "max_p": args.max_p,
                "sample_size": sample_size,
                "trials": trials,
                "seed": args.seed,
                "shell_decimals": args.shell_decimals,
                "eaxis_decimals": args.eaxis_decimals,
                "channels_mod420": CHANNELS_MOD420,
            }
        )
        if False
        else {
            "max_p": args.max_p,
            "sample_size": sample_size,
            "trials": trials,
            "seed": args.seed,
            "shell_decimals": args.shell_decimals,
            "eaxis_decimals": args.eaxis_decimals,
            "channels_mod420": CHANNELS_MOD420,
        },
        "counts": {
            "quadruplets_total_found": len(quadruplets),
            "quadruplets_mod420_filtered": n_items,
            "quadruplets_shell9": n_shell9,
        },
        "layers_of_nine": {
            "fraction_quadruplet_shell9": frac_shell9,
            "wilson_quad_level": wilson_quad,
            "trials_all9": trials_all9,
            "fraction_trials_all9": frac_trials_all9,
            "wilson_trial_level": wilson_trials,
        },
        "eaxis": {
            "z_shift_used": z_shift,
            "summary": eaxis_summary,
            "peaks": eaxis_peaks,
        },
        "fourier_mod420_ring": fourier_by_ch,
        "global_norm_shells_outer": global_shells,
        "trials": [asdict(t) for t in trial_list],
    }
    report_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(f"Report: {report_txt}")
    print(f"JSON  : {report_json}")


if __name__ == "__main__":
    main()

