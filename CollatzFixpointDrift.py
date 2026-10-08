"""
CollatzFixpointDrift — Driftgesetz und Lyapunov-Kandidat.

Entlang eines Pfades ω₁…ωₙ:
  x_n* = b_n/(1−A_n)
  χ_n  = dist(x_n*, {0,1,2})
  Λ_n  = log₂ A_n
  D_n  = |x_n* − x_{n+1}^*|
  V_n  = χ_n + α max(Λ_n, 0)

Fragen (Scan):
  1) Σ D_n < ∞ auf Collatz-Bahnen?
  2) x_n* → {0,1,2} (χ_n→0) — wann?
  3) Ist V_{n+1}<V_n außerhalb des Attraktors? (globale Monotonie)

Epistemik: trennt lokalen Kanal von Erreichbarkeit; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
import random
from fractions import Fraction
from pathlib import Path


LOG2_3 = math.log2(3)
ATTRACTOR = (0, 1, 2)


def affine_of_word(omega: tuple[int, ...]) -> tuple[int, int, int]:
    a = L = b = 0
    for e in omega:
        if e == 0:
            L += 1
        else:
            b = 3 * b + (1 << L)
            a += 1
            L += 1
    return a, L, b


def invariants(omega: tuple[int, ...]) -> dict:
    a, L, B = affine_of_word(omega)
    if L == 0:
        return {
            "Lambda": 0.0,
            "x": Fraction(0),
            "chi": Fraction(0),
            "A": Fraction(0),
        }
    A = Fraction(3**a, 1 << L)
    x = Fraction(B, (1 << L) - 3**a)
    chi = min(abs(x - Fraction(t)) for t in ATTRACTOR)
    return {
        "a": a,
        "L": L,
        "B": B,
        "A": A,
        "Lambda": a * LOG2_3 - L,
        "x": x,
        "chi": chi,
    }


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


def path_series(bits: list[int]) -> list[dict]:
    """Serie (x_n*, χ_n, Λ_n, D_n) für Präfixe n=1..len."""
    rows = []
    prev_x = None
    for n in range(1, len(bits) + 1):
        inv = invariants(tuple(bits[:n]))
        D = None if prev_x is None else abs(float(inv["x"] - prev_x))
        rows.append(
            {
                "n": n,
                "x": inv["x"],
                "x_f": float(inv["x"]),
                "chi": float(inv["chi"]),
                "Lambda": inv["Lambda"],
                "Lambda_pos": max(inv["Lambda"], 0.0),
                "D": D,
            }
        )
        prev_x = inv["x"]
    return rows


def V(row: dict, alpha: float) -> float:
    return row["chi"] + alpha * row["Lambda_pos"]


def analyze_series(rows: list[dict], enter: int | None, alphas: tuple[float, ...]) -> dict:
    Ds = [r["D"] for r in rows if r["D"] is not None]
    chis = [r["chi"] for r in rows]
    pre = [r for r in rows if enter is None or r["n"] < enter]
    post = [r for r in rows if enter is not None and r["n"] >= enter]

    def sum_D(part: list[dict]) -> float:
        return sum(r["D"] for r in part if r["D"] is not None)

    out: dict = {
        "sum_D_all": sum(Ds),
        "sum_D_pre": sum_D(pre),
        "sum_D_post": sum_D(post),
        "chi_final": chis[-1] if chis else None,
        "chi_min": min(chis) if chis else None,
        "chi_mean_pre": (
            sum(r["chi"] for r in pre) / len(pre) if pre else None
        ),
        "chi_mean_post": (
            sum(r["chi"] for r in post) / len(post) if post else None
        ),
        "chi_enter": next((r["chi"] for r in rows if enter and r["n"] == enter), None),
        "chi_final_over_enter": None,
        "enter": enter,
        "lyapunov": {},
    }
    if out["chi_enter"] and out["chi_enter"] > 0 and out["chi_final"] is not None:
        out["chi_final_over_enter"] = out["chi_final"] / out["chi_enter"]

    for alpha in alphas:
        Vs = [V(r, alpha) for r in rows]
        dec = sum(1 for i in range(len(Vs) - 1) if Vs[i + 1] < Vs[i] - 1e-15)
        strict = sum(1 for i in range(len(Vs) - 1) if Vs[i + 1] < Vs[i] - 1e-15)
        noninc = sum(1 for i in range(len(Vs) - 1) if Vs[i + 1] <= Vs[i] + 1e-12)
        # nur Nicht-Attraktor-Phase: n < enter
        if enter and enter > 2:
            Vp = [V(r, alpha) for r in rows if r["n"] < enter]
            noninc_pre = sum(
                1 for i in range(len(Vp) - 1) if Vp[i + 1] <= Vp[i] + 1e-12
            )
            frac_pre = noninc_pre / max(len(Vp) - 1, 1)
            # Anzahl Verletzungen V_{n+1}>=V_n
            viol_pre = (len(Vp) - 1) - noninc_pre
        else:
            frac_pre = None
            viol_pre = None
        out["lyapunov"][str(alpha)] = {
            "frac_decrease": dec / max(len(Vs) - 1, 1),
            "frac_nonincrease": noninc / max(len(Vs) - 1, 1),
            "frac_nonincrease_pre": frac_pre,
            "violations_pre": viol_pre,
            "V_start": Vs[0],
            "V_end": Vs[-1],
            "V_min": min(Vs),
            "V_max": max(Vs),
            "mono_global": noninc == len(Vs) - 1,
            "mono_pre": frac_pre == 1.0 if frac_pre is not None else None,
        }
    return out


def main() -> int:
    random.seed(7)
    steps = 240
    alphas = (0.0, 0.25, 0.5, 1.0, 2.0, 5.0)
    n_seeds = 30

    print("=" * 72)
    print("FIXPUNKT-DRIFT  —  Dₙ=|xₙ*−xₙ₊₁*|  und  Vₙ=χₙ+α Λₙ⁺")
    print("=" * 72)
    print("Lokal: Kanal via (Λ,x*).  Global: existiert monotone V?\n")

    seeds = [3, 5, 7, 9, 11, 15, 17, 19, 27, 31, 41, 47, 63, 97, 127, 255, 703]
    seeds += [random.randint(10**3, 10**6) for _ in range(n_seeds - len(seeds))]
    seeds = seeds[:n_seeds]

    # [1] Kontroll: reine Attraktorpfade
    print("[1] Kontrolle (10)^∞ / Zufall")
    bits10 = ([1, 0] * (steps // 2))
    rows10 = path_series(bits10)
    a10 = analyze_series(rows10, enter=1, alphas=alphas)
    print(
        f"  (10)^∞:  ΣD={a10['sum_D_all']:.6g}  χ_final={a10['chi_final']:.3g}  "
        f"mono V(α=1)? {a10['lyapunov']['1.0']['mono_global']}"
    )

    bits_rand = [random.randint(0, 1) for _ in range(steps)]
    rows_r = path_series(bits_rand)
    a_r = analyze_series(rows_r, enter=None, alphas=alphas)
    print(
        f"  Zufall:  ΣD={a_r['sum_D_all']:.6g}  χ_final={a_r['chi_final']:.3g}  "
        f"meanχ={sum(r['chi'] for r in rows_r)/len(rows_r):.3f}  "
        f"frac ΔV≤0 (α=1)={a_r['lyapunov']['1.0']['frac_nonincrease']:.3f}"
    )

    # [2] Collatz-Bahnen Drift
    print(f"\n[2] Collatz-Bahnen: ΣD und χ-Konvergenz ({n_seeds} Seeds)")
    analyses = []
    for n0 in seeds:
        bits, states = syracuse_bits(n0, steps)
        enter = next((i for i, s in enumerate(states) if s in (1, 2)), None)
        # enter index = state index = prefix length to arrive
        rows = path_series(bits)
        an = analyze_series(rows, enter=enter, alphas=alphas)
        an["n0"] = n0
        analyses.append(an)

    sumD_pre = [a["sum_D_pre"] for a in analyses]
    sumD_post = [a["sum_D_post"] for a in analyses]
    ratios = [
        a["chi_final_over_enter"]
        for a in analyses
        if a["chi_final_over_enter"] is not None
    ]
    print(
        f"  ΣD_pre:  median={sorted(sumD_pre)[len(sumD_pre)//2]:.4g}  "
        f"mean={sum(sumD_pre)/len(sumD_pre):.4g}  max={max(sumD_pre):.4g}"
    )
    print(
        f"  ΣD_post: median={sorted(sumD_post)[len(sumD_post)//2]:.4g}  "
        f"mean={sum(sumD_post)/len(sumD_post):.4g}  max={max(sumD_post):.4g}"
    )
    print(
        f"  χ_final/χ_enter: median={sorted(ratios)[len(ratios)//2]:.4g}  "
        f"(≪1 ⇒ Kontraktion nach Eintritt)"
    )

    # Beispiel 27
    print("\n  Drift-Beispiel n0=27 (Auswahl):")
    bits27, st27 = syracuse_bits(27, steps)
    enter27 = next(i for i, s in enumerate(st27) if s in (1, 2))
    rows27 = path_series(bits27)
    for n in [10, 20, 40, 60, enter27 - 1, enter27, enter27 + 20, enter27 + 60]:
        if n < 1 or n >= len(rows27):
            continue
        r = rows27[n - 1]
        mark = " enter" if n == enter27 else ""
        print(
            f"    n={n:3d}  state={st27[n]:8d}  χ={r['chi']:.4g}  "
            f"Λ={r['Lambda']:+.3f}  D={r['D'] if r['D'] is not None else float('nan'):.4g}{mark}"
        )

    # [3] Lyapunov-Monotonie
    print("\n[3] Lyapunov-Kandidat V=χ+α Λ⁺ — Monotonie?")
    print("  (global und nur Transient-Phase n < enter)")
    for alpha in alphas:
        key = str(alpha)
        mono_g = sum(1 for a in analyses if a["lyapunov"][key]["mono_global"])
        mono_p = sum(
            1
            for a in analyses
            if a["lyapunov"][key]["mono_pre"] is True
        )
        frac_pre = [
            a["lyapunov"][key]["frac_nonincrease_pre"]
            for a in analyses
            if a["lyapunov"][key]["frac_nonincrease_pre"] is not None
        ]
        viol = [
            a["lyapunov"][key]["violations_pre"]
            for a in analyses
            if a["lyapunov"][key]["violations_pre"] is not None
        ]
        print(
            f"  α={alpha:4g}:  mono_global {mono_g}/{n_seeds}  "
            f"mono_pre {mono_p}/{n_seeds}  "
            f"mean frac(ΔV≤0)_pre={sum(frac_pre)/len(frac_pre):.3f}  "
            f"mean viol_pre={sum(viol)/len(viol):.1f}"
        )

    # Wo bricht V?
    print("\n[4] Struktur der Verletzungen (α=1, n0=27)")
    a_focus = 1.0
    V27 = [V(r, a_focus) for r in rows27]
    viol_idx = [
        i + 1
        for i in range(len(V27) - 1)
        if V27[i + 1] > V27[i] + 1e-12 and (i + 1) < enter27
    ]
    print(f"  #Verletzungen vor enter: {len(viol_idx)}")
    if viol_idx:
        # korreliere mit Odd-Schritten / Λ-Anstieg
        sample = viol_idx[:8] + viol_idx[len(viol_idx) // 2 : len(viol_idx) // 2 + 3]
        for i in sample[:12]:
            r0, r1 = rows27[i - 1], rows27[i]
            bit = bits27[i]  # step that produced prefix i+1? 
            # prefix length i -> i+1 uses bits[i]
            step_bit = bits27[i] if i < len(bits27) else -1
            print(
                f"    n={i}→{i+1}  bit={bits27[i]}  "
                f"Δχ={r1['chi']-r0['chi']:+.4g}  "
                f"ΔΛ={r1['Lambda']-r0['Lambda']:+.4g}  "
                f"ΔV={V27[i]-V27[i-1]:+.4g}  "
                f"state={st27[i+1]}"
            )

    # Odd steps often increase Lambda by log2(3)-1 > 0
    print("\n[5] Mechanismus: Odd-Schritte zerstören Monotonie")
    odd_up = even_up = odd_n = even_n = 0
    for n0 in seeds[:15]:
        bits, states = syracuse_bits(n0, steps)
        enter = next(i for i, s in enumerate(states) if s in (1, 2))
        rows = path_series(bits)
        Vs = [V(r, 1.0) for r in rows]
        for i in range(min(enter, len(Vs)) - 1):
            bit = bits[i]  # extending prefix i to i+1? 
            # rows[i] is prefix length i+1, step bits[i] was last
            # transition rows[i] -> rows[i+1] uses bits[i+1]
            if i + 1 >= enter:
                break
            b = bits[i + 1]
            up = Vs[i + 1] > Vs[i] + 1e-12
            if b == 1:
                odd_n += 1
                odd_up += int(up)
            else:
                even_n += 1
                even_up += int(up)
    print(
        f"  vor enter:  Odd-Schritte mit ΔV>0:  {odd_up}/{odd_n} "
        f"({odd_up/max(odd_n,1):.3f})"
    )
    print(
        f"              Even-Schritte mit ΔV>0: {even_up}/{even_n} "
        f"({even_up/max(even_n,1):.3f})"
    )
    print(
        "  Odd erhöht Λ um log₂3−1≈0.585 und kann χ stark ändern —\n"
        "  daher ist V=χ+αΛ⁺ auf der Transient-Phase nicht Lyapunov."
    )

    # [6] Nach Eintritt: Drift klein, Summe konvergent-ähnlich
    print("\n[6] Post-Enter: kontraktives Driftgesetz")
    print(
        f"  mean ΣD_post={sum(sumD_post)/len(sumD_post):.4g}  "
        f"≪ mean ΣD_pre={sum(sumD_pre)/len(sumD_pre):.4g}"
    )
    # geometrische Abnahme von χ nach enter
    geo_ok = 0
    for n0 in seeds:
        bits, states = syracuse_bits(n0, steps)
        enter = next(i for i, s in enumerate(states) if s in (1, 2))
        if enter + 60 > len(bits):
            continue
        rows = path_series(bits)
        c0 = rows[enter - 1]["chi"]
        c1 = rows[enter + 59]["chi"]
        # erwarte Faktor ~(3/4)^30 für 30 Perioden (10)
        if c0 > 0 and c1 < c0 * 0.01:
            geo_ok += 1
    print(f"  χ(enter+60)<0.01·χ(enter): {geo_ok}/{n_seeds}")

    print("\n[7] Fazit")
    print(
        "  ΣD_post klein / χ→0 nach Eintritt: lokaler Kanal kontrahiert.\n"
        "  ΣD_pre groß und V nicht monoton: keine globale Lyapunov-Form\n"
        "  V=χ+αΛ⁺ auf der Transient-Phase (Odd-Schritte).\n"
        "  Driftgesetz erklärt Stabilität IM Kanal, nicht Erreichbarkeit."
    )

    out = Path(__file__).resolve().parent / "collatz_fixpoint_drift.json"
    payload = {
        "claim": (
            "Fixpoint drift: post-enter sum D small and chi contracts; "
            "V=chi+alpha*Lambda_+ is not monotone on transients (odd steps); "
            "no global Lyapunov found; does not prove Collatz"
        ),
        "steps": steps,
        "alphas": alphas,
        "control_10": a10,
        "control_random": {
            "sum_D": a_r["sum_D_all"],
            "chi_final": a_r["chi_final"],
            "lyapunov_1": a_r["lyapunov"]["1.0"],
        },
        "collatz": {
            "n_seeds": n_seeds,
            "sum_D_pre": {
                "mean": sum(sumD_pre) / len(sumD_pre),
                "median": sorted(sumD_pre)[len(sumD_pre) // 2],
                "max": max(sumD_pre),
            },
            "sum_D_post": {
                "mean": sum(sumD_post) / len(sumD_post),
                "median": sorted(sumD_post)[len(sumD_post) // 2],
                "max": max(sumD_post),
            },
            "chi_final_over_enter_median": sorted(ratios)[len(ratios) // 2]
            if ratios
            else None,
            "lyapunov": {
                str(alpha): {
                    "mono_global_count": sum(
                        1
                        for a in analyses
                        if a["lyapunov"][str(alpha)]["mono_global"]
                    ),
                    "mono_pre_count": sum(
                        1
                        for a in analyses
                        if a["lyapunov"][str(alpha)]["mono_pre"] is True
                    ),
                    "mean_frac_nonincrease_pre": sum(
                        a["lyapunov"][str(alpha)]["frac_nonincrease_pre"]
                        for a in analyses
                        if a["lyapunov"][str(alpha)]["frac_nonincrease_pre"]
                        is not None
                    )
                    / n_seeds,
                }
                for alpha in alphas
            },
            "odd_even_violations_alpha1": {
                "odd_up_frac": odd_up / max(odd_n, 1),
                "even_up_frac": even_up / max(even_n, 1),
            },
        },
        "separation": {
            "local": "channel stability via (Lambda, x*) after entry",
            "global": "reachability of {1,2} — open; no monotone V found",
        },
        "examples": [
            {"n0": a["n0"], "enter": a["enter"], "sum_D_pre": a["sum_D_pre"],
             "sum_D_post": a["sum_D_post"], "chi_final": a["chi_final"]}
            for a in analyses[:12]
        ],
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
