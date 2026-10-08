"""
CollatzDensityResonance — asymptotische / Fenster-Odd-Dichte auf zulässigen Bahnen.

Zyklische Resonanz ist algebraisch entschärft (B≥0 ⇒ kein A≥1-Zyklus).
Offen: aperiodische Resonanz

  Λ(Ω) = limsup (1/N) Σ (a_i log2 3 − m_i) ≥ 0
       ⇔  p(Ω) = limsup (Σ a_i)/(Σ m_i) ≥ p_c = log2/log3.

Konvention: Syracuse (ein /2 pro Schritt) — konsistent mit Path-/Spectrum-Atlas.
  odd:  n ↦ (3n+1)/2
  even: n ↦ n/2
  p_L = (#odd) / L
  Λ_L = p_L · log2(3) − 1

Zusätzlich: Dichte-Atlas nach n mod 2^k.

Claim-boundary: ergodische Diagnostik; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
import random
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


LOG2_3 = math.log2(3)
P_CRIT = math.log(2) / math.log(3)  # ≈ 0.630929753
LAMBDA_HALF = 0.5 * LOG2_3 - 1.0  # ≈ -0.2075


def syracuse_step(n: int) -> tuple[int, int]:
    """(n_neu, odd_bit) mit odd_bit∈{0,1}."""
    if n % 2 == 0:
        return n // 2, 0
    return (3 * n + 1) // 2, 1


def elementary_step(n: int) -> tuple[int, int]:
    """Elementar: odd zählt 3n+1 ohne /2."""
    if n % 2 == 0:
        return n // 2, 0
    return 3 * n + 1, 1


@dataclass
class WindowStats:
    L: int
    p_max: float = 0.0
    Lambda_max: float = float("-inf")
    p_at_max: float = 0.0
    n_seed: int = 0
    n_at_peak: int = 0
    step_at_peak: int = 0
    # Verteilung über Seeds: End-Fenster p_L
    p_end_max: float = 0.0
    p_end_mean: float = 0.0
    Lambda_end_max: float = float("-inf")
    above_crit_end: int = 0
    above_crit_peak: int = 0
    n_seeds: int = 0


def scan_seed_syracuse(
    n0: int,
    L: int,
    track_running_max: bool = True,
) -> tuple[float, float, float, float, int, int]:
    """
    Liefert:
      p_end, Λ_end, p_peak, Λ_peak, n_at_peak, step_at_peak
    peak = laufendes Maximum von p_t für t=1..L.
    """
    n = n0
    odd = 0
    p_peak = 0.0
    Lam_peak = float("-inf")
    n_peak = n0
    step_peak = 0
    for t in range(1, L + 1):
        n, bit = syracuse_step(n)
        odd += bit
        if n < 1:
            n = 1
        p = odd / t
        Lam = p * LOG2_3 - 1.0
        if track_running_max and Lam > Lam_peak:
            Lam_peak = Lam
            p_peak = p
            n_peak = n
            step_peak = t
    p_end = odd / L
    Lam_end = p_end * LOG2_3 - 1.0
    if not track_running_max:
        p_peak, Lam_peak, n_peak, step_peak = p_end, Lam_end, n, L
    return p_end, Lam_end, p_peak, Lam_peak, n_peak, step_peak


def window_campaign(
    seeds: list[int],
    windows: list[int],
) -> dict[int, WindowStats]:
    out: dict[int, WindowStats] = {L: WindowStats(L=L) for L in windows}
    for L in windows:
        st = out[L]
        p_ends = []
        for n0 in seeds:
            p_end, Lam_end, p_peak, Lam_peak, n_peak, step_peak = scan_seed_syracuse(
                n0, L
            )
            st.n_seeds += 1
            p_ends.append(p_end)
            if Lam_end > st.Lambda_end_max:
                st.Lambda_end_max = Lam_end
                st.p_end_max = p_end
            if Lam_peak > st.Lambda_max:
                st.Lambda_max = Lam_peak
                st.p_max = p_peak
                st.n_seed = n0
                st.n_at_peak = n_peak
                st.step_at_peak = step_peak
            if p_end >= P_CRIT:
                st.above_crit_end += 1
            if p_peak >= P_CRIT:
                st.above_crit_peak += 1
        st.p_end_mean = sum(p_ends) / len(p_ends) if p_ends else 0.0
    return out


def residue_atlas(
    k: int,
    L: int,
    reps_per_class: int = 3,
    rng: random.Random | None = None,
) -> dict:
    """
    Für jede Restklasse r mod 2^k: p_L, Λ_L auf Repräsentanten.
    Repräsentanten: r, r+2^k, r+2·2^k, … (positiv).
    """
    rng = rng or random.Random(k * L)
    mod = 1 << k
    rows = []
    p_max_global = 0.0
    Lam_max_global = float("-inf")
    best_r = 0
    above = 0

    for r in range(mod):
        ps = []
        Lams = []
        for j in range(reps_per_class):
            n0 = r + j * mod
            if n0 == 0:
                n0 = mod  # 2^k
            p_end, Lam_end, p_peak, Lam_peak, *_ = scan_seed_syracuse(n0, L)
            # für Atlas: Enddichte des Fensters + Peak
            ps.append(p_end)
            Lams.append(Lam_end)
            if Lam_peak > Lam_max_global:
                Lam_max_global = Lam_peak
                p_max_global = p_peak
                best_r = r
        p_mean = sum(ps) / len(ps)
        Lam_mean = sum(Lams) / len(Lams)
        p_hi = max(ps)
        Lam_hi = max(Lams)
        if p_hi >= P_CRIT:
            above += 1
        rows.append(
            {
                "r": r,
                "p_mean": p_mean,
                "p_max": p_hi,
                "Lambda_mean": Lam_mean,
                "Lambda_max": Lam_hi,
                "above_crit": p_hi >= P_CRIT,
            }
        )

    rows_sorted = sorted(rows, key=lambda d: -d["Lambda_max"])
    return {
        "k": k,
        "L": L,
        "mod": mod,
        "classes_above_crit_end": above,
        "best_peak": {
            "r": best_r,
            "p_peak_global": p_max_global,
            "Lambda_peak_global": Lam_max_global,
        },
        "top_classes": rows_sorted[:12],
        "bottom_classes": rows_sorted[-6:],
        "p_mean_over_classes": sum(d["p_mean"] for d in rows) / len(rows),
        "Lambda_mean_over_classes": sum(d["Lambda_mean"] for d in rows) / len(rows),
        "frac_classes_Lambda_gt_half": sum(
            1 for d in rows if d["Lambda_mean"] > LAMBDA_HALF
        )
        / len(rows),
    }


def long_orbit_probe(seeds: list[int], L: int = 50_000) -> list[dict]:
    """Lange Einzelbahnen: laufendes p_t und Abstand zu p_c."""
    probes = []
    for n0 in seeds:
        n = n0
        odd = 0
        max_p = 0.0
        max_Lam = float("-inf")
        last_above = None
        crossings = 0
        was_above = False
        for t in range(1, L + 1):
            n, bit = syracuse_step(n)
            odd += bit
            p = odd / t
            Lam = p * LOG2_3 - 1.0
            if p > max_p:
                max_p = p
            if Lam > max_Lam:
                max_Lam = Lam
            above = p >= P_CRIT
            if above and not was_above:
                crossings += 1
                last_above = t
            was_above = above
        probes.append(
            {
                "n0": n0,
                "L": L,
                "p_end": odd / L,
                "Lambda_end": (odd / L) * LOG2_3 - 1.0,
                "p_peak": max_p,
                "Lambda_peak": max_Lam,
                "crit_crossings": crossings,
                "last_above_step": last_above,
            }
        )
    return probes


def main() -> int:
    print("=" * 72)
    print("COLLATZ-DICHTERESONANZ  —  p_L und Λ_L (Syracuse)")
    print("=" * 72)
    print(f"p_c = {P_CRIT:.12f}")
    print(f"Λ(1/2) ≈ {LAMBDA_HALF:.6f}  (heuristische Typikalität)")
    print("Λ ≥ 0  ⇔  p ≥ p_c")

    windows = [100, 500, 1000, 5000]
    # Seeds: 1..N plus große Zufallsstarts
    rng = random.Random(42)
    N = 5000
    seeds = list(range(1, N + 1))
    seeds += [rng.randrange(10**6, 10**7) for _ in range(500)]
    seeds += [rng.randrange(10**9, 10**10) for _ in range(200)]

    print(f"\n[1] Fenster-Kampagne  |seeds|={len(seeds)}  L∈{windows}")
    camp = window_campaign(seeds, windows)
    print(
        f"  {'L':>6}  {'p_end_mean':>10}  {'p_end_max':>10}  {'Λ_end_max':>10}  "
        f"{'p_peak':>10}  {'Λ_peak':>10}  {'#p≥pc end':>10}  {'#peak≥pc':>10}"
    )
    for L in windows:
        st = camp[L]
        print(
            f"  {L:6d}  {st.p_end_mean:10.6f}  {st.p_end_max:10.6f}  {st.Lambda_end_max:10.6f}  "
            f"{st.p_max:10.6f}  {st.Lambda_max:10.6f}  {st.above_crit_end:10d}  "
            f"{st.above_crit_peak:10d}"
        )
        print(
            f"         peak-seed={st.n_seed} @step={st.step_at_peak} → n={st.n_at_peak}"
        )

    print("\n[2] Lange Bahnen L=50000 (Stichprobe)")
    long_seeds = [27, 703, 9663, 871, 6171, 27_000 + 1] + [
        rng.randrange(10**5, 10**6) for _ in range(10)
    ]
    probes = long_orbit_probe(long_seeds, L=50_000)
    for pr in probes[:8]:
        print(
            f"  n0={pr['n0']:<10d}  p_end={pr['p_end']:.5f}  Λ_end={pr['Lambda_end']:+.5f}  "
            f"p_peak={pr['p_peak']:.5f}  Λ_peak={pr['Lambda_peak']:+.5f}  "
            f"cross={pr['crit_crossings']}"
        )
    n_cross = sum(1 for pr in probes if pr["crit_crossings"] > 0)
    print(f"  Bahnen mit p_t≥p_c irgendwann: {n_cross}/{len(probes)}")

    print("\n[3] Dichte-Atlas nach Restklassen mod 2^k")
    atlases = []
    for k, L in [(8, 500), (10, 1000), (12, 2000)]:
        atl = residue_atlas(k, L, reps_per_class=2, rng=rng)
        atlases.append(atl)
        print(
            f"  k={k} L={L}:  class_p_mean={atl['p_mean_over_classes']:.5f}  "
            f"class_Λ_mean={atl['Lambda_mean_over_classes']:+.5f}  "
            f"classes p_end≥p_c: {atl['classes_above_crit_end']}/{atl['mod']}"
        )
        print("    top Λ_end:")
        for row in atl["top_classes"][:5]:
            print(
                f"      r={row['r']:4d}  p_max={row['p_max']:.5f}  "
                f"Λ_max={row['Lambda_max']:+.5f}"
            )

    # Frühe Peak-Artefakte: p_peak bei kleinem t oft hoch
    print("\n[4] Peak-Artefakt: p_t für t≤20 auf Seeds 1..2000")
    early_above = 0
    late_above = 0
    for n0 in range(1, 2001):
        n = n0
        odd = 0
        hit_early = False
        hit_late = False
        for t in range(1, 501):
            n, bit = syracuse_step(n)
            odd += bit
            p = odd / t
            if p >= P_CRIT:
                if t <= 20:
                    hit_early = True
                if t >= 200:
                    hit_late = True
        if hit_early:
            early_above += 1
        if hit_late:
            late_above += 1
    print(f"  Seeds mit p_t≥p_c für ein t≤20:  {early_above}/2000")
    print(f"  Seeds mit p_t≥p_c für ein t≥200: {late_above}/2000")

    print("\n[5] Fazit")
    # typische Enddichte
    typ = camp[1000]
    print(
        f"  Typisch (L=1000): p_end_mean={typ.p_end_mean:.4f} "
        f"(Heuristik 0.5), Λ_end_mean≈{typ.p_end_mean*LOG2_3-1:.4f} "
        f"(Heuristik {LAMBDA_HALF:.4f})"
    )
    print(
        "  Enddichten ≥ p_c sind selten/transient; Peaks oft früh (kleine t).\n"
        "  Collatz-Reformulierung: Ausschluss von limsup p ≥ p_c auf zulässigen Ω.\n"
        "  2-adischer Atlas: systematische Λ-Unterschiede zwischen Restklassen?"
    )

    out = Path(__file__).resolve().parent / "collatz_density_resonance.json"
    payload = {
        "claim": (
            "Collatz density resonance: Lambda>=0 iff odd-density >= log2/log3; "
            "window and 2-adic atlas; not a proof"
        ),
        "p_crit": P_CRIT,
        "Lambda_half": LAMBDA_HALF,
        "convention": "syracuse",
        "windows": {
            str(L): {
                "p_end_mean": camp[L].p_end_mean,
                "p_end_max": camp[L].p_end_max,
                "Lambda_end_max": camp[L].Lambda_end_max,
                "p_peak_max": camp[L].p_max,
                "Lambda_peak_max": camp[L].Lambda_max,
                "above_crit_end": camp[L].above_crit_end,
                "above_crit_peak": camp[L].above_crit_peak,
                "peak_seed": camp[L].n_seed,
                "n_seeds": camp[L].n_seeds,
            }
            for L in windows
        },
        "long_probes": probes,
        "residue_atlases": [
            {
                "k": a["k"],
                "L": a["L"],
                "p_mean_over_classes": a["p_mean_over_classes"],
                "Lambda_mean_over_classes": a["Lambda_mean_over_classes"],
                "classes_above_crit_end": a["classes_above_crit_end"],
                "top_classes": a["top_classes"][:8],
                "best_peak": a["best_peak"],
            }
            for a in atlases
        ],
        "early_vs_late_crit": {
            "early_t_le_20": early_above,
            "late_t_ge_200": late_above,
            "n": 2000,
        },
        "reformulation": "Collatz <=> no admissible Omega with limsup odd-density >= p_c",
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
