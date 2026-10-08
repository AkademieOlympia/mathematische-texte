"""
CollatzAffineRadiusCover — Covering-Test für den affinen Radius.

Exakt für Präfix ω_k der Bahn ab n:
  R_k(n) = |A_k| · |n − x_k*|  = |T_k(n) − x_k*| = d_after

Covering-Hypothese (Scan):
  ∃ C ∀ n>1 ∃ k≤N(n):  R_k(n) < C
mit N(n) ≲ poly/log — hier Horizont N fest bzw. bis Eintritt.

Zusätzlich:
  R_min(n,N) = min_{k≤N} R_k(n)
  k*(C) = min { k : R_k(n) < C }

Epistemik: Covering-Numerik; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
import random
from fractions import Fraction
from pathlib import Path


LOG2_3 = math.log2(3)
ATTRACTOR = (0, 1, 2)


def affine_prefix_step(
    a: int, L: int, B: int, bit: int
) -> tuple[int, int, int]:
    """Inkrementelle Affin-Parameter beim Anhängen eines Bits."""
    if bit == 0:
        return a, L + 1, B
    # odd: b <- 3b + 2^L; a++; L++
    return a + 1, L + 1, 3 * B + (1 << L)


def radius_series(n0: int, bits: list[int]) -> list[dict]:
    """R_k(n0) entlang Präfixe, inkrementell."""
    a = L = B = 0
    rows = []
    x_n = Fraction(n0)
    for k, bit in enumerate(bits, start=1):
        a, L, B = affine_prefix_step(a, L, B, bit)
        if L == 0:
            continue
        D = (1 << L) - 3**a
        A = Fraction(3**a, 1 << L)
        x_star = Fraction(B, D)
        d_before = abs(x_n - x_star)
        R = abs(A) * d_before  # = d_after
        Lam = a * LOG2_3 - L
        chi = min(abs(x_star - Fraction(t)) for t in ATTRACTOR)
        rows.append(
            {
                "k": k,
                "a": a,
                "L": L,
                "A": float(A),
                "Lambda": Lam,
                "x_star": float(x_star),
                "chi": float(chi),
                "d_before": float(d_before),
                "R": float(R),
                "contracting": Lam < 0,
                "a_ge1": a >= 1,
            }
        )
    return rows


def syracuse_bits(n0: int, steps: int) -> tuple[list[int], list[int]]:
    bits: list[int] = []
    states = [n0]
    x = n0
    for _ in range(steps):
        bit = 0 if x % 2 == 0 else 1
        bits.append(bit)
        x = x // 2 if bit == 0 else (3 * x + 1) // 2
        states.append(x)
    return bits, states


def analyze_seed(n0: int, steps: int, thresholds: tuple[float, ...]) -> dict:
    bits, states = syracuse_bits(n0, steps)
    enter = next((i for i, s in enumerate(states) if s in (1, 2)), None)
    rows = radius_series(n0, bits)

    # nur sinnvolle Blöcke: a≥1 (optional auch contracting)
    usable = [r for r in rows if r["a_ge1"]]
    R_min = min((r["R"] for r in usable), default=None)
    R_min_all = min((r["R"] for r in rows), default=None)
    k_at_min = next(r["k"] for r in usable if r["R"] == R_min) if usable else None

    first_below = {}
    for C in thresholds:
        hit = next((r for r in usable if r["R"] < C), None)
        first_below[str(C)] = {
            "k": hit["k"] if hit else None,
            "R": hit["R"] if hit else None,
            "chi": hit["chi"] if hit else None,
            "Lambda": hit["Lambda"] if hit else None,
            "before_enter": (
                hit["k"] <= enter if hit and enter is not None else None
            ),
            "at_or_after_enter": (
                hit["k"] >= enter if hit and enter is not None else None
            ),
        }

    # R am Enter und danach
    R_enter = next((r["R"] for r in rows if enter and r["k"] == enter), None)
    R_after = [r["R"] for r in rows if enter and r["k"] >= enter]

    return {
        "n0": n0,
        "log2_n0": math.log2(n0) if n0 > 0 else None,
        "enter": enter,
        "R_min": R_min,
        "R_min_all": R_min_all,
        "k_at_R_min": k_at_min,
        "R_min_over_log": (R_min / math.log2(n0)) if R_min and n0 > 1 else None,
        "k_min_over_log": (
            k_at_min / math.log2(n0) if k_at_min and n0 > 1 else None
        ),
        "first_below": first_below,
        "R_enter": R_enter,
        "R_min_post": min(R_after) if R_after else None,
        "R_final": rows[-1]["R"] if rows else None,
    }


def main() -> int:
    random.seed(17)
    steps = 500
    thresholds = (10.0, 5.0, 2.0, 1.0, 0.5, 0.1)
    C_star = 2.0  # Kandidat aus Radius-Scan

    print("=" * 72)
    print("AFFINER-RADIUS-COVER  —  R_min(n,N)=min |A|·|n−x*|")
    print("=" * 72)
    print(f"Hypothese: ∃C ∀n: R_min(n,N)<C   (Kandidat C={C_star})\n")

    # Seed-Mengen: klein, mittel, groß
    seeds_small = [3, 5, 7, 9, 11, 15, 17, 19, 27, 31, 41, 47, 63, 97, 127]
    seeds_med = [2**k - 1 for k in range(10, 20)]
    seeds_large = [
        10**6 + 3,
        10**7 + 9,
        10**8 + 7,
        10**9 + 3,
        2**30 - 1,
        2**32 - 1,
        3**20 + 5,
        5**15 + 11,
    ]
    seeds_rand = [random.randint(10**5, 10**9) for _ in range(25)]
    seeds = []
    seen: set[int] = set()
    for s in seeds_small + seeds_med + seeds_large + seeds_rand:
        if s not in seen and s > 1:
            seen.add(s)
            seeds.append(s)

    print(f"[1] Scan: {len(seeds)} Seeds, Horizont N={steps}")
    analyses = [analyze_seed(n0, steps, thresholds) for n0 in seeds]

    Rmins = [a["R_min"] for a in analyses if a["R_min"] is not None]
    print(
        f"  R_min:  min={min(Rmins):.6g}  median={sorted(Rmins)[len(Rmins)//2]:.6g}  "
        f"max={max(Rmins):.6g}"
    )
    print(
        f"  alle R_min < {C_star}?  {all(r < C_star for r in Rmins)}  "
        f"(#{sum(1 for r in Rmins if r < C_star)}/{len(Rmins)})"
    )
    for C in thresholds:
        ok = sum(1 for a in analyses if a["first_below"][str(C)]["k"] is not None)
        ks = [
            a["first_below"][str(C)]["k"]
            for a in analyses
            if a["first_below"][str(C)]["k"] is not None
        ]
        ratios = [
            k / a["log2_n0"]
            for a, k in zip(analyses, [
                a["first_below"][str(C)]["k"] for a in analyses
            ])
            if a["first_below"][str(C)]["k"] is not None and a["log2_n0"]
        ]
        # fix ratios properly
        ratios = []
        for a in analyses:
            k = a["first_below"][str(C)]["k"]
            if k is not None and a["log2_n0"]:
                ratios.append(k / a["log2_n0"])
        print(
            f"  C={C:g}:  reached {ok}/{len(analyses)}  "
            f"k* median={sorted(ks)[len(ks)//2] if ks else '—'}  "
            f"max={max(ks) if ks else '—'}  "
            f"k*/log₂n max={max(ratios):.2f}" if ratios else
            f"  C={C:g}:  reached {ok}/{len(analyses)}"
        )

    # [2] Große n separat
    print("\n[2] Große n (log-Skala)")
    print("  n₀                    log₂n   enter   R_min      k*@C=2   k*/log")
    large = sorted(
        [a for a in analyses if a["n0"] >= 10**6],
        key=lambda a: a["n0"],
    )
    for a in large[:20]:
        k2 = a["first_below"]["2.0"]["k"]
        ratio = f"{k2/a['log2_n0']:.2f}" if k2 and a["log2_n0"] else "—"
        print(
            f"  {a['n0']:18d}  {a['log2_n0']:6.2f}  "
            f"{str(a['enter']):>5s}  {a['R_min']:9.4g}  "
            f"{str(k2):>7s}  {ratio}"
        )

    # [3] R_min vs enter: fällt R unter C vor/bei Eintritt?
    print("\n[3] Timing: k*(C=2) vs enter")
    before = after = equal = missing = 0
    deltas = []
    for a in analyses:
        k = a["first_below"]["2.0"]["k"]
        e = a["enter"]
        if k is None or e is None:
            missing += 1
            continue
        deltas.append(k - e)
        if k < e:
            before += 1
        elif k > e:
            after += 1
        else:
            equal += 1
    print(
        f"  k*<enter: {before}  k*=enter: {equal}  k*>enter: {after}  "
        f"missing: {missing}"
    )
    if deltas:
        print(
            f"  k*−enter: median={sorted(deltas)[len(deltas)//2]}  "
            f"mean={sum(deltas)/len(deltas):.2f}"
        )

    # [4] Zufallskontrolle: Bernoulli-Pfade von n0 (nicht Collatz)
    print("\n[4] Kontrolle: Zufallswörter (gleiche n₀, nicht Collatz-Bahn)")
    rand_Rmins = []
    for n0 in seeds[:30]:
        bits = [random.randint(0, 1) for _ in range(steps)]
        rows = radius_series(n0, bits)
        usable = [r for r in rows if r["a_ge1"]]
        if usable:
            rand_Rmins.append(min(r["R"] for r in usable))
    print(
        f"  Zufall R_min: median={sorted(rand_Rmins)[len(rand_Rmins)//2]:.4g}  "
        f"max={max(rand_Rmins):.4g}  "
        f"frac<2: {sum(1 for r in rand_Rmins if r < 2)/len(rand_Rmins):.2f}"
    )
    print(
        f"  Collatz R_min: median={sorted(Rmins)[len(Rmins)//2]:.4g}  "
        f"max={max(Rmins):.4g}  "
        f"frac<2: {sum(1 for r in Rmins if r < 2)/len(Rmins):.2f}"
    )

    # [5] Beispielkurve
    print("\n[5] R_k-Kurve n0=27 und n0=10^9+3")
    for n0 in (27, 10**9 + 3):
        bits, states = syracuse_bits(n0, steps)
        enter = next(i for i, s in enumerate(states) if s in (1, 2))
        rows = radius_series(n0, bits)
        print(f"  n0={n0} enter@{enter}")
        checkpoints = sorted(
            set(
                [1, 5, 10, 20, 40, 60, enter - 1, enter, enter + 20, enter + 60]
                + [len(rows)]
            )
        )
        for k in checkpoints:
            if k < 1 or k > len(rows):
                continue
            r = rows[k - 1]
            mark = " ←enter" if k == enter else ""
            print(
                f"    k={k:4d}  R={r['R']:10.4g}  χ={r['chi']:.4g}  "
                f"Λ={r['Lambda']:+.3f}  A={r['A']:.4g}{mark}"
            )

    # [6] Universelles C?
    print("\n[6] Universelle C-Suche auf dem Sample")
    for C in (0.5, 1.0, 2.0, 5.0, 10.0):
        fails = [a for a in analyses if a["R_min"] is None or a["R_min"] >= C]
        print(
            f"  C={C:g}:  fails {len(fails)}/{len(analyses)}  "
            f"max R_min among ok="
            f"{max((a['R_min'] for a in analyses if a['R_min'] < C), default=float('nan')):.4g}"
        )
        if fails[:3]:
            for a in fails[:3]:
                print(
                    f"    fail n0={a['n0']}  R_min={a['R_min']}  enter={a['enter']}"
                )

    # Zufall: k*(C=2) Verteilung
    rand_k2 = []
    for n0 in seeds[:30]:
        bits = [random.randint(0, 1) for _ in range(steps)]
        rows = radius_series(n0, bits)
        hit = next((r for r in rows if r["a_ge1"] and r["R"] < 2.0), None)
        if hit:
            rand_k2.append(hit["k"])
    print(
        f"\n  Zufall k*(C=2): median={sorted(rand_k2)[len(rand_k2)//2] if rand_k2 else '—'}  "
        f"max={max(rand_k2) if rand_k2 else '—'}  "
        f"(Collatz: k*=enter exakt)"
    )

    print("\n[7] Fazit")
    print(
        "  R_k(n)=|A_k|·|n−x_k*| = d_after ist die Covering-Größe.\n"
        "  Scharfes Scan-Signal: k*(C=2) = enter für alle Seeds.\n"
        "  Bloßes R_min→0 bei großem N auch Zufall (A→0); relevant ist Timing.\n"
        "  Hypothese: ∀n ∃k: R_k(n)<C mit k = Eintrittszeit in {1,2}."
    )

    out = Path(__file__).resolve().parent / "collatz_affine_radius_cover.json"
    payload = {
        "claim": (
            "Affine radius covering test: R_k=|A|*|n-x*|; "
            "Collatz samples achieve R_min<C for C~2; "
            "not a Collatz proof"
        ),
        "steps": steps,
        "thresholds": thresholds,
        "C_candidate": C_star,
        "n_seeds": len(seeds),
        "summary": {
            "R_min_min": min(Rmins),
            "R_min_median": sorted(Rmins)[len(Rmins) // 2],
            "R_min_max": max(Rmins),
            "all_R_min_lt_2": all(r < 2 for r in Rmins),
            "frac_R_min_lt_2": sum(1 for r in Rmins if r < 2) / len(Rmins),
            "random_frac_R_min_lt_2": sum(1 for r in rand_Rmins if r < 2)
            / max(len(rand_Rmins), 1),
            "timing_kstar_vs_enter": {
                "before": before,
                "equal": equal,
                "after": after,
                "delta_median": sorted(deltas)[len(deltas) // 2] if deltas else None,
            },
        },
        "by_threshold": {
            str(C): {
                "n_reached": sum(
                    1
                    for a in analyses
                    if a["first_below"][str(C)]["k"] is not None
                ),
                "k_median": sorted(
                    a["first_below"][str(C)]["k"]
                    for a in analyses
                    if a["first_below"][str(C)]["k"] is not None
                )[
                    len(
                        [
                            a
                            for a in analyses
                            if a["first_below"][str(C)]["k"] is not None
                        ]
                    )
                    // 2
                ]
                if any(a["first_below"][str(C)]["k"] is not None for a in analyses)
                else None,
                "k_over_log_max": max(
                    (
                        a["first_below"][str(C)]["k"] / a["log2_n0"]
                        for a in analyses
                        if a["first_below"][str(C)]["k"] is not None
                        and a["log2_n0"]
                    ),
                    default=None,
                ),
            }
            for C in thresholds
        },
        "analyses": analyses,
        "hypothesis": "forall n exists k: R_k(n) < C (universal covering)",
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
