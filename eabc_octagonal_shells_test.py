#!/usr/bin/env python3
"""
Numerisches Minimal-Experiment zu
"Division" von Prime-Vierlingen (Primzahlvierlinge) in einer oktagonalen Darstellung,
inkl. ABCE vs CEAB (zwei "Hurwitz half-shells") sowie einer E-Achsen-Analyse.

Wichtige Hinweise / Annahmen:
- In diesem Projekt existieren (in den gesuchten Quellen) keine expliziten Implementierungen
  für "octagonal representation", "Hurwitz half-shells", "R8-balls", "E-axis" oder "layers of nine".
- Deshalb definiert dieses Skript eine *reproduzierbare* Abbildung von
  (p, p+2, p+6, p+8) -> 9 Punkte (8 Oktagonpunkte + Mittelpunkt) in 3D,
  wobei ABCE "nach oben" und CEAB "nach unten" gespiegelt wird.
- "layers of nine" wird hier operationalisiert als: nach Zentrierung (Mittelpunkt -> Ursprung)
  gibt es genau 9 verschiedene Norm-Abstände (Norm-Schichten) unter den 9 Punkten.

Ausgabe:
- Text-Report (.txt)
- Rohresultate (.json)

Laufzeit:
- Für max_p=1e7 und einige Dutzend Trials ist es typischerweise noch gut handhabbar.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np


CHANNELS_MOD420 = [11, 101, 191, 221, 311, 401]
CHANNEL_TO_INDEX = {r: i for i, r in enumerate(CHANNELS_MOD420)}


def isprime_sieve(limit: int) -> np.ndarray:
    """Sieve of Eratosthenes in numpy-boolean Form."""
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for k in range(2, int(math.isqrt(limit)) + 1):
        if sieve[k]:
            sieve[k * k : limit + 1 : k] = False
    return sieve


def find_prime_quadruplets(max_p: int) -> list[int]:
    """
    Findet Starts p so dass p, p+2, p+6, p+8 alle prim sind.
    """
    sieve = isprime_sieve(max_p + 8)
    out: list[int] = []
    # Starts p bei Vierlingen sind i.d.R. ungerade; p=2..4 liefern ohnehin keine Vierlinge.
    for p in range(5, max_p + 1, 2):
        if not (sieve[p] and sieve[p + 2] and sieve[p + 6] and sieve[p + 8]):
            continue
        out.append(p)
    return out


def chirality_from_p(p: int) -> str | None:
    """
    Chiralität analog zu den Projekt-Skripten:
    - ABCE wenn p % 12 == 5
    - CEAB wenn p % 12 == 11
    """
    r12 = p % 12
    if r12 == 5:
        return "ABCE"
    if r12 == 11:
        return "CEAB"
    return None


def channel_index_mod420(p: int) -> int | None:
    r = p % 420
    return CHANNEL_TO_INDEX.get(r, None)


def phi_for_family_code(code: str) -> float:
    """
    Phasenwinkel für {E,A,B,C}.
    (Wir übernehmen bewusst die gleiche Winkelwahl wie in den EABC-Riemann-Skripten üblich.)
    """
    mapping = {"E": 0.0, "A": math.pi / 2, "B": math.pi, "C": 3 * math.pi / 2}
    if code not in mapping:
        raise KeyError(code)
    return mapping[code]


def component_radii_logs(p: int) -> np.ndarray:
    comps = np.array([p, p + 2, p + 6, p + 8], dtype=float)
    return np.log(comps)


def quad_points_octagonal_3d(
    p: int,
    channel_idx: int,
    chirality: str,
    *,
    eps: float = 1e-12,
) -> np.ndarray:
    """
    Erzeuge genau 9 Punkte in 3D:
      - 8 Oktagonpunkte (k=0..7)
      - 1 Mittelpunktpunkt (Mittelpunkt der 8 Oktagonpunkte; wird später explizit als 9. Punkt genutzt)

    Operationale Konstruktion:
    - Jede Komponente log(p+offset) wird mit einem Oktagon-Radius verknüpft.
    - Zwei "Halbschalen" werden über half=k//4 unterschieden (Radius-Offset via mod-420 Kanalindex).
    - ABCE spiegelt nach oben, CEAB nach unten (E-Achse ~ z).
    - Zusätzlich wird die Oktagon-Phasenlage (theta_offset) chiraler gewählt.
    """
    radii4 = component_radii_logs(p)  # 4 radii basierend auf (p,p+2,p+6,p+8)
    boost = math.log(channel_idx + 1.0 + eps)
    chirality_sign = +1.0 if chirality == "ABCE" else -1.0
    theta_offset = +math.pi / 4 if chirality == "ABCE" else -math.pi / 4

    # Acht Oktagonpunkte
    pts = []
    base_angles = np.arange(8, dtype=float) * (math.pi / 4)
    for k in range(8):
        idx = k % 4
        half = k // 4  # 0 oder 1
        r = radii4[idx] + half * boost
        ang = base_angles[k] + theta_offset
        x = r * math.cos(ang)
        y = r * math.sin(ang)
        z = chirality_sign * r
        pts.append((x, y, z))

    pts8 = np.array(pts, dtype=float)  # shape (8,3)
    center = pts8.mean(axis=0, keepdims=True)  # 9. Punkt
    pts9 = np.vstack([pts8, center])  # shape (9,3)
    return pts9


def unique_shell_count_from_points(
    pts9: np.ndarray,
    *,
    decimals: int = 7,
) -> tuple[int, list[float]]:
    """
    zählt diskrete Norm-Schichten:
    - Normabstände der 9 Punkte bzgl. Ursprung
    - nach Zentrierung sollte der Mittelpunktpunkt nahe 0 liegen
    """
    norms = np.linalg.norm(pts9, axis=1)
    rounded = np.round(norms, decimals=decimals)
    uniq = np.unique(rounded)
    uniq_list = [float(x) for x in uniq.tolist()]
    return len(uniq_list), uniq_list


def z_histogram_peaks(z_values: np.ndarray, *, decimals: int = 6) -> dict[str, int]:
    """Rundet z-Werte und zählt Häufigkeit der Peaks."""
    z_round = np.round(z_values, decimals=decimals)
    u, counts = np.unique(z_round, return_counts=True)
    return {str(float(x)): int(c) for x, c in zip(u, counts)}


@dataclass
class TrialResult:
    trial_index: int
    sample_size: int
    shells_all_9: bool
    quad_all_9_fraction: float
    mean_shells: float
    min_shells: int
    max_shells: int
    shells_counts: list[int]


def bootstrap_wilson_interval(k: int, n: int, *, alpha: float = 0.05) -> dict[str, float]:
    """
    Wilson score interval für eine Bernoulli-Quote.
    """
    if n == 0:
        return {"p_hat": float("nan"), "lo": float("nan"), "hi": float("nan")}
    z = 1.959963984540054  # ~ 97.5% quantile for alpha=0.05
    p_hat = k / n
    denom = 1 + z * z / n
    center = p_hat + z * z / (2 * n)
    half = z * math.sqrt((p_hat * (1 - p_hat) + z * z / (4 * n)) / n)
    lo = (center - half) / denom
    hi = (center + half) / denom
    return {"p_hat": p_hat, "lo": lo, "hi": hi}


def dft_magnitudes_on_channels(counts_by_channel: np.ndarray) -> dict[str, list[float]]:
    """
    Analog zu mod-420 Ring-Fourier auf C6:
    counts_by_channel shape (6,)
    """
    v = np.asarray(counts_by_channel, dtype=float)
    if v.sum() > 0:
        v = v / v.sum()
    F = np.fft.fft(v)
    return {
        "magnitudes": [float(abs(Fk)) for Fk in F],
        "phases": [float(np.angle(Fk)) for Fk in F],
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max_p", type=int, default=1_000_000, help="Suche bis p<=max_p")
    ap.add_argument("--sample_size", type=int, default=120, help="Vierlinge pro Trial")
    ap.add_argument("--trials", type=int, default=60, help="Anzahl random Trials")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--shell_decimals", type=int, default=7, help="Rundung für Shell-Zählung")
    ap.add_argument("--eaxis_decimals", type=int, default=6, help="Rundung für E-Achsen-Peaks")
    ap.add_argument("--report_name", type=str, default="eabc_octagonal_shells_test")
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()

    t0 = time.time()

    max_p: int = args.max_p
    sample_size: int = args.sample_size
    trials: int = args.trials

    out_dir = Path(__file__).resolve().parent
    report_txt = out_dir / f"{args.report_name}.txt"
    report_json = out_dir / f"{args.report_name}.json"

    # 1) Vierlinge finden
    print(f"Finde Prime-Vierlinge bis max_p={max_p:,} ...")
    quadruplets = find_prime_quadruplets(max_p)
    # 2) mod-420 Kanäle & Chiralität filtern
    items = []
    for p in quadruplets:
        ch = chirality_from_p(p)
        if ch is None:
            continue
        cidx = channel_index_mod420(p)
        if cidx is None:
            continue
        items.append((p, ch, cidx))

    rng = np.random.default_rng(args.seed)
    n_items = len(items)
    print(f"Gefunden (nach mod-420 Channels): {n_items}")

    if n_items == 0:
        raise SystemExit("Keine mod-420 gefilterten Vierlinge gefunden. max_p erhöhen.")

    sample_size_eff = int(min(sample_size, n_items))
    trials_eff = int(max(1, trials))

    # Precompute representations + shell counts (schnell für Trials)
    # Wir speichern pro p eine Shell-Zahl (9? ggf. kleiner).
    precomputed = []
    for p, ch, cidx in items:
        pts9 = quad_points_octagonal_3d(p, cidx, ch)
        # Zentriere pro Quad (Mittelpunkt -> Ursprung), damit "layers of nine" konsistent ist.
        center = pts9[-1, :]  # centroid
        pts9_centered = pts9.copy()
        pts9_centered[:-1, :] -= center
        pts9_centered[-1, :] = 0.0
        shells_count, shells_uniq = unique_shell_count_from_points(
            pts9_centered, decimals=args.shell_decimals
        )
        precomputed.append(
            {
                "p": p,
                "chirality": ch,
                "channel_idx": int(cidx),
                "shells_count": int(shells_count),
                "shells_uniq": shells_uniq,
                "pts9_centered": pts9_centered.tolist(),  # klein genug: 9 Punkte
            }
        )

    shells_counts_all = np.array([x["shells_count"] for x in precomputed], dtype=int)
    eq9_total = int(np.sum(shells_counts_all == 9))
    frac_eq9_total = eq9_total / len(shells_counts_all)
    print(f"gesamt: shells==9 für {eq9_total}/{len(shells_counts_all)} Vierlinge ({frac_eq9_total:.4f})")

    # Trials
    trial_results: list[TrialResult] = []
    all_trials_all9 = 0
    quad_total_scored = 0
    quad_total_eq9 = 0

    for tr in range(trials_eff):
        idxs = rng.choice(len(precomputed), size=sample_size_eff, replace=False)
        shells_counts = shells_counts_all[idxs]
        quad_total_scored += int(shells_counts.size)
        quad_total_eq9 += int(np.sum(shells_counts == 9))

        shells_all_9 = bool(np.all(shells_counts == 9))
        if shells_all_9:
            all_trials_all9 += 1

        res = TrialResult(
            trial_index=tr,
            sample_size=sample_size_eff,
            shells_all_9=shells_all_9,
            quad_all_9_fraction=float(np.mean(shells_counts == 9)),
            mean_shells=float(np.mean(shells_counts)),
            min_shells=int(np.min(shells_counts)),
            max_shells=int(np.max(shells_counts)),
            shells_counts=[int(x) for x in shells_counts.tolist()],
        )
        trial_results.append(res)

    frac_trials_all9 = all_trials_all9 / trials_eff
    # Anteil der einzelnen Vierlinge (Bernoulli)
    wilson_quad = bootstrap_wilson_interval(quad_total_eq9, quad_total_scored)
    wilson_trials = bootstrap_wilson_interval(all_trials_all9, trials_eff)

    # 4) Transformation + E-axis Analyse (auf einem großen Sample)
    # Nimm den größten Sample: sample_size_eff oder alle (wenn klein)
    take_n = int(min(max(200, sample_size_eff), n_items))
    idxs_big = rng.choice(len(precomputed), size=take_n, replace=False) if take_n < n_items else np.arange(n_items)

    # Sammle alle z-Werte und Norm-Schalen (auf Quad-zentrierten Punkten)
    # Außerdem: Punkte nach Chirality gruppieren für "R8-ball midpoint -> 0".
    pts_by_ch = {"ABCE": [], "CEAB": []}
    shell_r_by_ch = {"ABCE": [], "CEAB": []}
    channel_counts_by_ch = {"ABCE": np.zeros(6, dtype=int), "CEAB": np.zeros(6, dtype=int)}

    all_point_z = {"ABCE": [], "CEAB": []}
    for i in idxs_big:
        item = precomputed[i]
        ch = item["chirality"]
        cidx = item["channel_idx"]
        pts9_centered = np.array(item["pts9_centered"], dtype=float)  # shape (9,3)
        # Distanz der 9 Punkte zum Ursprung: diese sind die "Norm-Schichten" innerhalb dieses Quad.
        norms = np.linalg.norm(pts9_centered, axis=1)
        pts_by_ch[ch].append(pts9_centered)
        shell_r_by_ch[ch].append(norms)
        # counts "pro Punkt" (8 äußere + 1 Mitte) im Kanal ring
        channel_counts_by_ch[ch][cidx] += int(pts9_centered.shape[0])
        # E-axis (z)
        all_point_z[ch].append(pts9_centered[:, 2])

    # R8-ball Mitte pro Chirality: als centroid der äußeren Punkte (ignoriert Mitte, die bei 0 liegt).
    # Wir holen die 8 äußeren Punkte aller Quads zusammen.
    z_shift = {}
    for ch in ["ABCE", "CEAB"]:
        if len(pts_by_ch[ch]) == 0:
            z_shift[ch] = 0.0
            continue
        pts_concat = np.vstack(pts_by_ch[ch])  # (Npts,3)
        # Mitte-Punkte liegen i.d.R. nahe 0, deshalb beschränken wir Shift auf äußere Punkte:
        # äußere Punkte entsprechen in unserem Konstruktions-pts9_centered den ersten 8 Einträgen pro Quad.
        # Wir approximieren: alle Punkte außer die letzten in jedem Block von 9.
        # -> einfach: n=9 und reshape.
        n_quads = len(pts_by_ch[ch])
        pts_blocks = np.array(pts_by_ch[ch], dtype=float)  # (n_quads,9,3)
        outer = pts_blocks[:, :8, :].reshape(-1, 3)
        z_shift[ch] = float(np.mean(outer[:, 2]))

    # Rekonstruieren z-Werte nach z_shift und Peaks/Histo.
    eaxis_peaks = {}
    eaxis_hist_summary = {}
    for ch in ["ABCE", "CEAB"]:
        if len(pts_by_ch[ch]) == 0:
            continue
        pts_blocks = np.array(pts_by_ch[ch], dtype=float)  # (n_quads,9,3)
        outer = pts_blocks[:, :8, :]
        outer_z = outer[:, :, 2].reshape(-1) - z_shift[ch]
        peaks = z_histogram_peaks(outer_z, decimals=args.eaxis_decimals)
        eaxis_peaks[ch] = peaks
        eaxis_hist_summary[ch] = {
            "z_shift_used": z_shift[ch],
            "n_outer_points": int(outer_z.size),
            "z_mean_after_shift": float(np.mean(outer_z)),
            "z_std_after_shift": float(np.std(outer_z)),
            "top_peaks": dict(list(peaks.items())[:10]),
        }

    # Fourier/Charakter-Projektion analog auf dem mod-420 Kanalring:
    # Nutzen counts über alle Punkte des großen Samples pro Kanal.
    fourier_by_ch = {}
    for ch in ["ABCE", "CEAB"]:
        fourier_by_ch[ch] = dft_magnitudes_on_channels(channel_counts_by_ch[ch])

    # Norm-Schalen Sicht: Clusterung über Normen aller Punkte (nach chirality Shift).
    # (Hier nehmen wir die Normen der 9 Punkte in jedem quad-zentrierten System => Norm(shell)=Distanz zu Quad-Zentrum.)
    # Für ein globales Bild, sammeln wir die äußeren Punkte Normen und zählen unique Rundungen.
    global_shells_by_ch = {}
    for ch in ["ABCE", "CEAB"]:
        norms_all_outer = []
        for block in shell_r_by_ch[ch]:
            norms_all_outer.append(block[:8])  # nur äußere Punkte
        if len(norms_all_outer) == 0:
            continue
        norms_all_outer = np.concatenate(norms_all_outer, axis=0)
        norms_r = np.unique(np.round(norms_all_outer, decimals=args.shell_decimals))
        global_shells_by_ch[ch] = norms_r.tolist()

    # Report schreiben
    elapsed = time.time() - t0
    report_lines = []
    report_lines.append("EABC Octagonal Shells Test (minimal, operational definitions)")
    report_lines.append("")
    report_lines.append(f"Projektpfad (root): {out_dir}")
    report_lines.append(f"max_p: {max_p:,}")
    report_lines.append(f"Gefundene Vierlinge (mod-420 filtered): {n_items}")
    report_lines.append("")
    report_lines.append("layers of nine (operational):")
    report_lines.append(
        "  Pro Vierling werden 8 Oktagonpunkte + Mittelpunkt gebildet, danach wird zentriert (Ball-Mitte -> Ursprung)."
    )
    report_lines.append("  layers-of-nine = Anzahl der *diskreten* Normabstände unter den 9 Punkten == 9.")
    report_lines.append("")
    report_lines.append(f"Shells==9 pro Vierling: {eq9_total}/{len(shells_counts_all)} ({frac_eq9_total:.4f})")
    report_lines.append(
        f"Wison-Intervall (Quad-Quote): p={wilson_quad['p_hat']:.4f}  CI=[{wilson_quad['lo']:.4f},{wilson_quad['hi']:.4f}]"
    )
    report_lines.append("")
    report_lines.append("Trials:")
    report_lines.append(f"  trials_eff={trials_eff}, sample_size_eff={sample_size_eff}")
    report_lines.append(
        f"  Trials mit 'alle Quads shells==9': {all_trials_all9}/{trials_eff} ({frac_trials_all9:.4f})"
    )
    report_lines.append(
        f"  Wilson-Intervall (Trial-Quote): p={wilson_trials['p_hat']:.4f}  CI=[{wilson_trials['lo']:.4f},{wilson_trials['hi']:.4f}]"
    )
    report_lines.append("")
    report_lines.append("E-axis Analyse nach 'move R8-ball midpoint -> 0' (per Chirality via z-Mittel):")
    report_lines.append(f"  z_shift ABCE: {eaxis_hist_summary.get('ABCE', {}).get('z_shift_used', None)}")
    report_lines.append(f"  z_shift CEAB: {eaxis_hist_summary.get('CEAB', {}).get('z_shift_used', None)}")
    report_lines.append("")
    for ch in ["ABCE", "CEAB"]:
        if ch not in eaxis_hist_summary:
            continue
        s = eaxis_hist_summary[ch]
        report_lines.append(
            f"  {ch}: n_outer_points={s['n_outer_points']}, z_mean_after_shift={s['z_mean_after_shift']:.6f}, z_std_after_shift={s['z_std_after_shift']:.6f}"
        )
        report_lines.append(f"    top_peaks (z_round->{args.eaxis_decimals}d): {s['top_peaks']}")
    report_lines.append("")
    report_lines.append("Fourier (DFT) der mod-420 Kanal-Punktecounts (Ring-C6):")
    report_lines.append("  Magnitudes over k=0..5:")
    for ch in ["ABCE", "CEAB"]:
        mags = fourier_by_ch[ch]["magnitudes"]
        report_lines.append(f"  {ch}: {['%.6f'%m for m in mags]}")
    report_lines.append("")
    report_lines.append("Globale Norm-Schichten (äußere Punkte) innerhalb der Quad-Zentrierung:")
    for ch in ["ABCE", "CEAB"]:
        if ch not in global_shells_by_ch:
            continue
        shells = global_shells_by_ch[ch]
        report_lines.append(f"  {ch}: #unique_norm_shells={len(shells)}")
        report_lines.append(f"    first few shells: {shells[:12]}")

    report_lines.append("")
    report_lines.append(f"Runtime: {elapsed:.2f}s")

    report_txt.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    payload = {
        "settings": {
            "max_p": max_p,
            "sample_size": sample_size_eff,
            "trials": trials_eff,
            "seed": args.seed,
            "shell_decimals": args.shell_decimals,
            "eaxis_decimals": args.eaxis_decimals,
            "channels_mod420": CHANNELS_MOD420,
        },
        "counts": {
            "quadruplets_total_found": len(quadruplets),
            "quadruplets_mod420_filtered": n_items,
        },
        "layers_of_nine": {
            "eq9_total": eq9_total,
            "eq9_frac_total": frac_eq9_total,
            "wilson_quad": wilson_quad,
            "trial_all9_total": all_trials_all9,
            "trial_all9_frac": frac_trials_all9,
            "wilson_trials": wilson_trials,
        },
        "eaxis": {
            "z_shift": z_shift,
            "peaks": eaxis_peaks,
            "hist_summary": eaxis_hist_summary,
            "global_norm_shells_outer_by_chirality": global_shells_by_ch,
            "fourier_by_chirality": fourier_by_ch,
        },
        "trial_results": [asdict(x) for x in trial_results],
        "precomputed_shell_counts": {
            "min": int(shells_counts_all.min()),
            "max": int(shells_counts_all.max()),
            "unique": {int(k): int(v) for k, v in zip(*np.unique(shells_counts_all, return_counts=True))},
        },
    }
    report_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(f"\nReport geschrieben: {report_txt}")
    print(f"JSON geschrieben: {report_json}")


if __name__ == "__main__":
    main()
