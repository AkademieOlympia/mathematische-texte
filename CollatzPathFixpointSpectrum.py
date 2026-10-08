"""
CollatzPathFixpointSpectrum — Pfadprodukte und Fixpunktfolge x_n*.

Für ein Wortprodukt / Präfix ω^{(n)} der Länge n:
  T_n(x) = A_n x + b_n
  Λ_n = log₂ A_n
  x_n* = b_n / (1 − A_n)
  χ_n = dist(x_n*, {0,1,2})

Fragen (Scan, kein Beweis):
  1) Periodische Attraktorpfade: χ_n = 0 für alle n?
  2) Echte Collatz-Bahnen: wird χ_n → 0 erst beim Eintritt in {1,2}?
  3) Zufällige / nichtperiodische Wörter: bleibt χ_n von {0,1,2} weg?
  4) Kann χ_n klein sein ohne dass das Präfix ein Attraktorwort ist?

Epistemik: dynamische Erweiterung von R_C=(Λ,x*); kein Collatz-Beweis.
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


def path_invariants(omega: tuple[int, ...]) -> dict:
    a, L, B = affine_of_word(omega)
    if L == 0:
        return {
            "n": 0,
            "a": 0,
            "L": 0,
            "Lambda": 0.0,
            "x_star": None,
            "chi": None,
            "chi_float": None,
            "A": None,
        }
    A = Fraction(3**a, 1 << L)
    D = (1 << L) - 3**a
    xQ = Fraction(B, D)
    chi = min(abs(xQ - Fraction(t)) for t in ATTRACTOR)
    return {
        "n": L,
        "a": a,
        "L": L,
        "B": B,
        "A": str(A),
        "Lambda": a * LOG2_3 - L,
        "x_star": str(xQ),
        "x_star_float": float(xQ),
        "chi": str(chi),
        "chi_float": float(chi),
        "chi_zero": chi == 0,
        "is_attractor_word": chi == 0,
    }


def syracuse_parity_path(n0: int, steps: int) -> tuple[list[int], list[int]]:
    """Paritätswort und Zustandsspur der echten Syracuse-Bahn."""
    bits: list[int] = []
    states = [n0]
    x = n0
    for _ in range(steps):
        bit = 0 if x % 2 == 0 else 1
        bits.append(bit)
        if bit == 0:
            x //= 2
        else:
            x = (3 * x + 1) // 2
        states.append(x)
    return bits, states


def prefix_spectrum(bits: list[int], stride: int = 1) -> list[dict]:
    """χ entlang wachsender Präfixe."""
    rows = []
    for n in range(stride, len(bits) + 1, stride):
        inv = path_invariants(tuple(bits[:n]))
        rows.append(inv)
    return rows


def block_composition_spectrum(bits: list[int], block: int) -> list[dict]:
    """
    Produkte aufeinanderfolgender Blöcke fester Länge:
    ω1, ω1ω2, ω1ω2ω3, … (Präfixe in Blockschritten).
    """
    rows = []
    for k in range(1, len(bits) // block + 1):
        inv = path_invariants(tuple(bits[: k * block]))
        inv["n_blocks"] = k
        rows.append(inv)
    return rows


def summarize_curve(rows: list[dict]) -> dict:
    if not rows:
        return {}
    chis = [r["chi_float"] for r in rows if r["chi_float"] is not None]
    lams = [r["Lambda"] for r in rows]
    first_zero = next((r["n"] for r in rows if r["chi_zero"]), None)
    last = rows[-1]
    # wie oft χ unter Schwellen
    return {
        "n_max": last["n"],
        "chi_final": last["chi"],
        "chi_final_f": last["chi_float"],
        "Lambda_final": last["Lambda"],
        "chi_min": min(chis) if chis else None,
        "chi_max": max(chis) if chis else None,
        "chi_mean": sum(chis) / len(chis) if chis else None,
        "frac_chi_lt_0_1": sum(1 for c in chis if c < 0.1) / len(chis),
        "frac_chi_zero": sum(1 for r in rows if r["chi_zero"]) / len(rows),
        "first_chi_zero_n": first_zero,
        "ended_chi_zero": bool(last["chi_zero"]),
    }


def main() -> int:
    random.seed(42)
    steps = 200
    n_seeds = 40

    print("=" * 72)
    print("PFAD-FIXPUNKTSPEKTRUM  —  ω₁…ωₙ ↦ (Λₙ, xₙ*, χₙ)")
    print("=" * 72)
    print("Frage: wird xₙ* bei langen nichtperiodischen Pfaden von {0,1,2}")
    print("weggedrückt, oder erst beim Attraktor-Eintritt χ=0?\n")

    # [1] Periodische Kontrollpfade
    print("[1] Periodische Attraktorpfade")
    controls = {}
    for name, period in (("10", (1, 0)), ("01", (0, 1)), ("0", (0,))):
        bits = list(period * (steps // len(period)))
        rows = prefix_spectrum(bits, stride=len(period))
        s = summarize_curve(rows)
        controls[name] = {"summary": s, "sample": rows[:: max(1, len(rows) // 8)]}
        print(
            f"  ({name})^∞  n≤{s['n_max']}:  "
            f"frac(χ=0)={s['frac_chi_zero']:.3f}  "
            f"χ_final={s['chi_final']}  Λ_final={s['Lambda_final']:+.3f}"
        )

    # [2] Echte Collatz-Bahnen
    print(f"\n[2] Collatz-Bahnen: {n_seeds} Seeds, je {steps} Schritte")
    seeds = [3, 5, 7, 9, 11, 15, 17, 19, 27, 31, 41, 47, 63, 97, 127]
    seeds += [random.randint(10**3, 10**6) for _ in range(n_seeds - len(seeds))]
    seeds = seeds[:n_seeds]

    orbit_rows = []
    entered = 0
    chi0_before_enter = 0
    chi0_only_after = 0
    for n0 in seeds:
        bits, states = syracuse_parity_path(n0, steps)
        # Eintritt in {1,2}: erster Index i mit states[i]∈{1,2}
        enter = next((i for i, s in enumerate(states) if s in (1, 2)), None)
        rows = prefix_spectrum(bits, stride=1)
        # optional subsample for storage
        s = summarize_curve(rows)
        first0 = s["first_chi_zero_n"]
        if enter is not None:
            entered += 1
            # χ=0 vor Eintritt? (Präfixlänge = Schritte; Eintritt bei state index)
            # Präfix der Länge k entspricht states[k]
            if first0 is not None and enter is not None:
                if first0 < enter:
                    chi0_before_enter += 1
                if first0 >= enter:
                    chi0_only_after += 1
        orbit_rows.append(
            {
                "n0": n0,
                "enter_idx": enter,
                "summary": s,
                "chi_at_enter": (
                    path_invariants(tuple(bits[:enter]))["chi"]
                    if enter and enter <= len(bits)
                    else None
                ),
                "chi_before_enter": (
                    path_invariants(tuple(bits[: max(enter - 1, 1)]))["chi_float"]
                    if enter and enter > 1
                    else None
                ),
                "curve_stride": [
                    rows[i]
                    for i in range(0, len(rows), max(1, len(rows) // 20))
                ],
            }
        )

    # aggregierte Statistik vor Attraktor
    pre_enter_chi = []
    post_enter_chi = []
    for o in orbit_rows:
        n0 = o["n0"]
        bits, states = syracuse_parity_path(n0, steps)
        enter = o["enter_idx"]
        if enter is None:
            for k in range(1, len(bits) + 1):
                pre_enter_chi.append(path_invariants(tuple(bits[:k]))["chi_float"])
            continue
        for k in range(1, enter):
            pre_enter_chi.append(path_invariants(tuple(bits[:k]))["chi_float"])
        for k in range(enter, len(bits) + 1):
            post_enter_chi.append(path_invariants(tuple(bits[:k]))["chi_float"])

    def stats(xs: list[float]) -> dict:
        if not xs:
            return {"n": 0}
        return {
            "n": len(xs),
            "min": min(xs),
            "max": max(xs),
            "mean": sum(xs) / len(xs),
            "frac_zero": sum(1 for c in xs if c == 0.0) / len(xs),
            "frac_lt_0_01": sum(1 for c in xs if c < 0.01) / len(xs),
            "frac_lt_0_1": sum(1 for c in xs if c < 0.1) / len(xs),
        }

    print(f"  Seeds mit Eintritt in {{1,2}} innerhalb {steps}: {entered}/{n_seeds}")
    print(
        f"  Seeds mit early χ=0 (kurzes Attraktor-Präfix): {chi0_before_enter}  "
        f"— lokal, nicht Endkanal"
    )
    sp, spo = stats(pre_enter_chi), stats(post_enter_chi)
    print(
        f"  Präfixe vor Eintritt:  n={sp.get('n',0)}  "
        f"mean χ={sp.get('mean',0):.4f}  frac(χ=0)={sp.get('frac_zero',0):.4f}  "
        f"frac(χ<0.1)={sp.get('frac_lt_0_1',0):.4f}"
    )
    print(
        f"  Präfixe nach Eintritt: n={spo.get('n',0)}  "
        f"mean χ={spo.get('mean',0):.4f}  frac(χ=0)={spo.get('frac_zero',0):.4f}  "
        f"frac(χ<0.1)={spo.get('frac_lt_0_1',0):.4f}  "
        f"(χ→0 asymptotisch, exakt 0 selten/nie mit Transient)"
    )

    # Beispielkurve
    print("\n  Beispiel n0=27 (χ entlang Präfix, Stride 10):")
    bits27, st27 = syracuse_parity_path(27, 120)
    enter27 = next(i for i, s in enumerate(st27) if s in (1, 2))
    for n in list(range(10, min(121, enter27 + 40), 10)) + [enter27]:
        if n < 1 or n > len(bits27):
            continue
        inv = path_invariants(tuple(bits27[:n]))
        mark = " ← enter" if n == enter27 else ""
        print(
            f"    n={n:3d}  state={st27[n]:8d}  Λ={inv['Lambda']:+7.3f}  "
            f"x*={inv['x_star']:>16s}  χ={inv['chi']}{mark}"
        )

    # [3] Zufällige nichtperiodische Wörter
    print(f"\n[3] Zufällige Paritätswörter (Bernoulli 1/2), {n_seeds} Pfade × {steps}")
    rand_summaries = []
    rand_chi_all = []
    n_ever_zero = 0
    for _ in range(n_seeds):
        bits = [random.randint(0, 1) for _ in range(steps)]
        # vermeide reine Attraktorwörter zufällig — ok, selten
        rows = prefix_spectrum(bits, stride=1)
        s = summarize_curve(rows)
        rand_summaries.append(s)
        if s["frac_chi_zero"] > 0:
            n_ever_zero += 1
        rand_chi_all.extend(r["chi_float"] for r in rows)

    sr = stats(rand_chi_all)
    print(
        f"  alle Präfixe: mean χ={sr['mean']:.4f}  "
        f"frac(χ=0)={sr['frac_zero']:.6f}  frac(χ<0.1)={sr['frac_lt_0_1']:.4f}"
    )
    print(f"  Pfade mit jemals χ=0: {n_ever_zero}/{n_seeds}")
    print(
        f"  finale χ: mean={sum(s['chi_final_f'] for s in rand_summaries)/len(rand_summaries):.4f}  "
        f"min={min(s['chi_final_f'] for s in rand_summaries):.4g}"
    )

    # [4] Blockprodukte entlang Collatz
    print("\n[4] Blockprodukte (Länge 4) entlang Collatz vs Zufall")
    blk = 4
    # Collatz block
    c_chi = []
    for n0 in seeds[:20]:
        bits, _ = syracuse_parity_path(n0, steps)
        for row in block_composition_spectrum(bits, blk):
            c_chi.append(row["chi_float"])
    r_chi = []
    for _ in range(20):
        bits = [random.randint(0, 1) for _ in range(steps)]
        for row in block_composition_spectrum(bits, blk):
            r_chi.append(row["chi_float"])
    sc, srr = stats(c_chi), stats(r_chi)
    print(
        f"  Collatz-Blöcke: mean χ={sc['mean']:.4f}  frac(χ=0)={sc['frac_zero']:.4f}"
    )
    print(
        f"  Zufall-Blöcke:  mean χ={srr['mean']:.4f}  frac(χ=0)={srr['frac_zero']:.6f}"
    )

    # [5] Asymptotik nach Eintritt: χ_n → 0?
    print("\n[5] Asymptotik: nach Eintritt in {1,2} gilt χₙ→0 (nicht exakt 0)")
    print(
        "  Mit transientem Kopf σ·(10)^k ist χ=0 unmöglich (nur reine Attraktorwörter).\n"
        "  Aber A_n→0 und x_n* wird an den Kanal gedrückt:"
    )
    post_curves = []
    for n0 in seeds[:15]:
        bits, states = syracuse_parity_path(n0, steps)
        enter = next(i for i, s in enumerate(states) if s in (1, 2))
        curve = []
        for dn in (0, 10, 20, 40, 80, 120):
            n = enter + dn
            if n > len(bits):
                break
            inv = path_invariants(tuple(bits[:n]))
            curve.append((dn, inv["chi_float"], inv["Lambda"]))
        post_curves.append({"n0": n0, "enter": enter, "curve": curve})
        if n0 in (27, 7, 31, seeds[0]):
            parts = " ".join(f"Δ{dn}:χ={chi:.3g}" for dn, chi, _ in curve)
            print(f"  n0={n0:6d} enter@{enter:3d}  {parts}")

    # Monotonie / Abnahme: χ(enter+80) < χ(enter) ?
    shrink = 0
    shrink_n = 0
    ratios = []
    for n0 in seeds:
        bits, states = syracuse_parity_path(n0, steps)
        enter = next(i for i, s in enumerate(states) if s in (1, 2))
        if enter + 80 > len(bits):
            continue
        c0 = path_invariants(tuple(bits[:enter]))["chi_float"]
        c1 = path_invariants(tuple(bits[: enter + 80]))["chi_float"]
        shrink_n += 1
        if c1 < c0:
            shrink += 1
        if c0 > 0:
            ratios.append(c1 / c0)
    print(
        f"  χ(enter+80)<χ(enter): {shrink}/{shrink_n}  "
        f"median ratio={sorted(ratios)[len(ratios)//2]:.4g}"
        if ratios
        else f"  χ(enter+80)<χ(enter): {shrink}/{shrink_n}"
    )

    # [6] Vor Eintritt vs Zufall: bleibt χ weg von 0?
    print("\n[6] Vor Eintritt / Zufall: Abstand zum Attraktorspiegel")
    epsilons = (1e-3, 1e-2, 0.1, 0.5)
    for eps in epsilons:
        frac = sum(1 for c in pre_enter_chi if c < eps) / max(len(pre_enter_chi), 1)
        fr = sum(1 for c in rand_chi_all if c < eps) / max(len(rand_chi_all), 1)
        print(f"  frac(χ<{eps:g}):  Collatz-pre={frac:.4f}  random={fr:.4f}")
    print(
        f"  mean χ:  Collatz-pre={sp['mean']:.3f}  "
        f"Collatz-post={spo['mean']:.3f}  random={sr['mean']:.3f}"
    )

    # kurze χ=0-Präfixe = lokaler Zufallstreffer reiner Attraktorwörter
    early_zeros = sum(
        1
        for o in orbit_rows
        if o["summary"]["first_chi_zero_n"] is not None
        and o["enter_idx"] is not None
        and o["summary"]["first_chi_zero_n"] < o["enter_idx"]
    )
    print(
        f"  Seeds mit early χ=0 (kurzes Attraktor-Präfix vor Eintritt): "
        f"{early_zeros}/{n_seeds}  — nicht der Endkanal"
    )

    print("\n[7] Fazit")
    print(
        "  Reine Attraktorpfade: χₙ≡0.\n"
        "  Collatz mit Transient σ·(10)^k: χₙ→0 asymptotisch nach Eintritt,\n"
        "  aber χₙ=0 nur für reine Attraktorwörter (nie mit Transient-Kopf).\n"
        "  Vor Eintritt / Zufall: χ typisch O(1), kein stabiler Kanal.\n"
        "  Offene Collatz-Hürde: Erreichen von {1,2}, nicht die χ-Asymptotik danach."
    )

    out = Path(__file__).resolve().parent / "collatz_path_fixpoint_spectrum.json"
    payload = {
        "claim": (
            "Path fixpoint spectrum: pure attractor paths have chi=0; "
            "Collatz paths after entry to {1,2} have chi->0 asymptotically "
            "under (10)-suffix but chi!=0 with transient head; "
            "pre-entry/random chi typically O(1); does not prove Collatz"
        ),
        "steps": steps,
        "n_seeds": n_seeds,
        "controls": {
            k: {"summary": v["summary"]} for k, v in controls.items()
        },
        "collatz_orbits": {
            "n_entered": entered,
            "pre_enter_chi": sp,
            "post_enter_chi": spo,
            "chi0_before_enter_seeds": chi0_before_enter,
            "note": (
                "After enter, chi->0 asymptotically under (10)-suffix, "
                "but chi=0 exactly never for prefixes with nontrivial transient"
            ),
            "post_enter_shrink": {
                "n": shrink_n,
                "chi_enter80_lt_chi_enter": shrink,
                "median_ratio": sorted(ratios)[len(ratios) // 2] if ratios else None,
            },
            "post_curves_sample": post_curves[:8],
            "examples": [
                {
                    "n0": o["n0"],
                    "enter_idx": o["enter_idx"],
                    "summary": o["summary"],
                }
                for o in orbit_rows[:15]
            ],
        },
        "random_paths": {
            "chi_all": sr,
            "n_paths_ever_chi_zero": n_ever_zero,
            "final_chi_mean": sum(s["chi_final_f"] for s in rand_summaries)
            / len(rand_summaries),
        },
        "block_len_4": {"collatz": sc, "random": srr},
        "open_question": (
            "Can an infinite admissible path keep chi_n>0 without divergence? "
            "Scan shows Collatz finite prefixes before entry have chi>0; "
            "entry forces approach to attractor class."
        ),
        "example_27_curve": [
            path_invariants(tuple(bits27[:n]))
            | {"state": st27[n], "enter": n == enter27}
            for n in list(range(10, min(121, enter27 + 40), 10)) + [enter27]
            if 1 <= n <= len(bits27)
        ],
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
