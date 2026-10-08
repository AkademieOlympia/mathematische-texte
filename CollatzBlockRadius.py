"""
CollatzBlockRadius — Block-Radius-Spektrum.

Hierarchie:
  Wort → (Λ, x*) → Block → Radius → Attraktor

Für einen Block ω mit Startzustand n (ganzzahlig):
  T(x) = A x + b,  x* = b/(1−A)
  χ = dist(x*, {0,1,2})
  d_before = |n − x*|
  d_after  = |T(n) − x*| = |A| · d_before   (affine Identität)
  T(n) = S^L(n) falls word(n)=ω

Vermutung (zu prüfen / ggf. falsifizieren):
  χ<ε  ⟹  d_before < r_ω
ist für Präfixe großer n₀ falsch.

Präzisere Kandidatenbedingung:
  χ<ε  und  d_after < r_safe  ⟹  Eintritt in {1,2}
bzw.
  χ<ε  und  |A|·n₀ ≲ O(1)  ⟹  Eintritt.

Epistemik: Radius-/Zylinderbrücke; kein Collatz-Beweis.
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


def block_data(omega: tuple[int, ...]) -> dict:
    a, L, B = affine_of_word(omega)
    A = Fraction(3**a, 1 << L)
    D = (1 << L) - 3**a
    x_star = Fraction(B, D)
    chi = min(abs(x_star - Fraction(t)) for t in ATTRACTOR)
    return {
        "a": a,
        "L": L,
        "B": B,
        "A": A,
        "A_float": float(A),
        "Lambda": a * LOG2_3 - L,
        "x_star": x_star,
        "chi": float(chi),
    }


def syracuse_apply(n: int, omega: tuple[int, ...]) -> tuple[int, bool]:
    """Wende ω an; ok=False falls Paritätswort nicht passt."""
    x = n
    for bit in omega:
        actual = 0 if x % 2 == 0 else 1
        if actual != bit:
            return x, False
        x = x // 2 if bit == 0 else (3 * x + 1) // 2
    return x, True


def syracuse_path(n0: int, steps: int) -> tuple[list[int], list[int]]:
    bits: list[int] = []
    states = [n0]
    x = n0
    for _ in range(steps):
        bit = 0 if x % 2 == 0 else 1
        bits.append(bit)
        x = x // 2 if bit == 0 else (3 * x + 1) // 2
        states.append(x)
    return bits, states


def enter_index(states: list[int]) -> int | None:
    return next((i for i, s in enumerate(states) if s in (1, 2)), None)


def prefix_radius_row(n0: int, bits: list[int], L: int) -> dict:
    omega = tuple(bits[:L])
    bd = block_data(omega)
    x_star = bd["x_star"]
    A = bd["A"]
    d_before = abs(float(Fraction(n0) - x_star))
    # echte Syracuse-Anwendung
    t_n, ok = syracuse_apply(n0, omega)
    d_after = abs(float(Fraction(t_n) - x_star))
    d_after_pred = abs(float(A)) * d_before
    return {
        "n0": n0,
        "L": L,
        "a": bd["a"],
        "Lambda": bd["Lambda"],
        "chi": bd["chi"],
        "A": bd["A_float"],
        "x_star": str(x_star),
        "d_before": d_before,
        "d_after": d_after,
        "d_after_pred": d_after_pred,
        "T_n": t_n,
        "word_ok": ok,
        "T_in_attractor": t_n in (1, 2),
        "A_times_n0": abs(bd["A_float"]) * n0,
    }


def find_prefix_block(bits: list[int], eps: float) -> int | None:
    for L in range(1, len(bits) + 1):
        bd = block_data(tuple(bits[:L]))
        if bd["a"] >= 1 and bd["Lambda"] < 0 and bd["chi"] < eps:
            return L
    return None


def estimate_r_for_word(
    omega: tuple[int, ...],
    r_candidates: list[float],
    sample_mod_bits: int = 12,
    max_iter: int = 200,
) -> dict:
    """
    Empirischer Radius: unter x im Zylinder C_ω (word match),
    für welche r gilt |x−x*|<r ⇒ Eintritt in {1,2} in max_iter Schritten.
    """
    bd = block_data(omega)
    x_star = float(bd["x_star"])
    # sampmle Residuen mod 2^m die ω als Anfangswort erlauben
    m = max(sample_mod_bits, len(omega) + 2)
    mod = 1 << m
    cylinder = []
    for r in range(mod):
        _, ok = syracuse_apply(r, omega)
        if ok:
            # hebe zu Repräsentanten nahe x* und darüber
            cylinder.append(r)
            cylinder.append(r + mod)
            if r + 2 * mod < 10**6:
                cylinder.append(r + 2 * mod)

    def enters(x: int) -> bool:
        y = x
        for _ in range(max_iter):
            if y in (1, 2):
                return True
            y = y // 2 if y % 2 == 0 else (3 * y + 1) // 2
            if y <= 0:
                return False
        return y in (1, 2)

    safe = {}
    for r in r_candidates:
        xs = [x for x in cylinder if abs(x - x_star) < r and x > 0]
        if not xs:
            safe[r] = {"n": 0, "frac_enter": None}
            continue
        # cap sample
        if len(xs) > 80:
            xs = xs[:: max(1, len(xs) // 80)][:80]
        ok = sum(1 for x in xs if enters(x))
        safe[r] = {"n": len(xs), "frac_enter": ok / len(xs), "all": ok == len(xs)}

    # größtes r mit frac_enter==1
    r_max = 0.0
    for r in sorted(r_candidates):
        info = safe[r]
        if info["n"] and info.get("all"):
            r_max = r
    return {
        "omega": "".join(map(str, omega)),
        "x_star": str(bd["x_star"]),
        "A": bd["A_float"],
        "chi": bd["chi"],
        "r_max_safe_scan": r_max,
        "radius_table": {str(r): safe[r] for r in r_candidates},
    }


def main() -> int:
    random.seed(13)
    steps = 400
    epsilons = (0.5, 0.1, 0.01, 1e-3)
    n_seeds = 35

    print("=" * 72)
    print("BLOCK-RADIUS-SPEKTRUM  —  χ, d=|n−x*|, |A|d, Eintritt")
    print("=" * 72)
    print("Frage: folgt aus χ<ε schon d<r — oder braucht man d_after?\n")

    seeds = [3, 5, 7, 9, 11, 15, 17, 19, 27, 31, 41, 47, 63, 97, 127, 255, 703]
    seeds += [2**k - 1 for k in (10, 12, 14, 16)]
    seeds += [random.randint(10**3, 10**7) for _ in range(n_seeds)]
    seen: set[int] = set()
    seeds = [s for s in seeds if s not in seen and not seen.add(s)][:n_seeds]  # type: ignore

    # [1] Empirischer Radius für reine Attraktorwörter
    print("[1] Empirischer Radius r_ω für Attraktorwörter")
    r_cands = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    pure = {
        "10": (1, 0),
        "01": (0, 1),
        "1010": (1, 0, 1, 0),
        "0101": (0, 1, 0, 1),
    }
    radius_atlas = {}
    for name, omega in pure.items():
        rad = estimate_r_for_word(omega, r_cands, sample_mod_bits=10)
        radius_atlas[name] = rad
        print(
            f"  ω={name:4s}  x*={rad['x_star']}  A={rad['A']:.4f}  "
            f"r_max_safe≈{rad['r_max_safe_scan']}"
        )

    # [2] Präfixblöcke: χ vs d_before vs d_after
    print("\n[2] Präfixblöcke auf Collatz-Bahnen: χ ⇏ d_before < r")
    rows_by_eps: dict[float, list[dict]] = {e: [] for e in epsilons}
    for n0 in seeds:
        bits, states = syracuse_path(n0, steps)
        enter = enter_index(states)
        for eps in epsilons:
            L = find_prefix_block(bits, eps)
            if L is None:
                continue
            row = prefix_radius_row(n0, bits, L)
            row["eps"] = eps
            row["enter"] = enter
            row["entered_by_T"] = row["T_in_attractor"]
            row["entered_eventually"] = enter is not None and enter <= steps
            # naives Lemma χ<ε ⇒ d < r_ref
            r_ref = 2.0  # großzügig um {0,1,2}
            row["naive_lemma"] = row["d_before"] < r_ref
            row["after_lemma"] = row["d_after"] < r_ref
            row["A_n_lemma"] = row["A_times_n0"] < r_ref
            rows_by_eps[eps].append(row)

    for eps in epsilons:
        rows = rows_by_eps[eps]
        if not rows:
            print(f"  ε={eps}: keine Blöcke")
            continue
        naive = sum(1 for r in rows if r["naive_lemma"])
        after = sum(1 for r in rows if r["after_lemma"])
        an = sum(1 for r in rows if r["A_n_lemma"])
        tin = sum(1 for r in rows if r["T_in_attractor"])
        print(
            f"  ε={eps:g}: n={len(rows)}  "
            f"d_before<2: {naive}/{len(rows)}  "
            f"d_after<2: {after}/{len(rows)}  "
            f"|A|n₀<2: {an}/{len(rows)}  "
            f"T(n)∈{{1,2}}: {tin}/{len(rows)}"
        )
        db = sorted(r["d_before"] for r in rows)
        da = sorted(r["d_after"] for r in rows)
        print(
            f"         d_before median={db[len(db)//2]:.4g}  "
            f"d_after median={da[len(da)//2]:.4g}  "
            f"median |A|={sorted(r['A'] for r in rows)[len(rows)//2]:.4g}"
        )

    # [3] Beispielkurve: wachsendes Präfix für n0=27
    print("\n[3] Beispiel n0=27: wachsendes Präfix bis χ klein")
    bits27, st27 = syracuse_path(27, steps)
    enter27 = enter_index(st27)
    print(f"  enter@{enter27}")
    for eps in (0.5, 0.1, 0.01, 1e-3):
        L = find_prefix_block(bits27, eps)
        if L is None:
            continue
        row = prefix_radius_row(27, bits27, L)
        print(
            f"  ε={eps:g}: L={L}  χ={row['chi']:.4g}  A={row['A']:.4g}  "
            f"d_bef={row['d_before']:.4g}  d_aft={row['d_after']:.4g}  "
            f"T={row['T_n']}  ∈Attr={row['T_in_attractor']}"
        )

    # [4] Korrelation: Eintritt ⇔ d_after klein
    print("\n[4] Brückenbedingung: d_after < r  vs  T(n)∈{1,2}")
    rows = rows_by_eps[0.1]
    for r_safe in (0.5, 1.0, 2.0, 5.0, 10.0):
        tp = sum(1 for r in rows if r["d_after"] < r_safe and r["T_in_attractor"])
        fp = sum(1 for r in rows if r["d_after"] < r_safe and not r["T_in_attractor"])
        fn = sum(1 for r in rows if r["d_after"] >= r_safe and r["T_in_attractor"])
        tn = sum(1 for r in rows if r["d_after"] >= r_safe and not r["T_in_attractor"])
        print(
            f"  r={r_safe:g}:  TP={tp} FP={fp} FN={fn} TN={tn}  "
            f"(d_after<r ⇒ Attraktor? prec={tp/max(tp+fp,1):.2f})"
        )

    # [5] χ<ε ∧ |A|n₀ < R
    print("\n[5] Alternative: χ<ε und |A|·n₀ < R")
    for R in (1, 2, 5, 10, 50, 100):
        subset = [r for r in rows if r["A_times_n0"] < R]
        if not subset:
            print(f"  R={R}: keine")
            continue
        tin = sum(1 for r in subset if r["T_in_attractor"])
        print(
            f"  R={R}: n={len(subset)}  T∈{{1,2}}: {tin}/{len(subset)}  "
            f"median d_after={sorted(r['d_after'] for r in subset)[len(subset)//2]:.4g}"
        )

    # [6] Fenster am aktuellen Zustand (nicht Präfix ab n0)
    print("\n[6] Lokale Fenster-Blöcke am Zustand s (ω=10/01)")
    local_rows = []
    for n0 in seeds:
        bits, states = syracuse_path(n0, min(steps, 200))
        for i in range(len(bits) - 1):
            omega = tuple(bits[i : i + 2])
            if omega not in ((1, 0), (0, 1)):
                continue
            s = states[i]
            bd = block_data(omega)
            d_b = abs(s - float(bd["x_star"]))
            t, ok = syracuse_apply(s, omega)
            d_a = abs(t - float(bd["x_star"]))
            local_rows.append(
                {
                    "n0": n0,
                    "i": i,
                    "s": s,
                    "omega": "".join(map(str, omega)),
                    "chi": bd["chi"],
                    "d_before": d_b,
                    "d_after": d_a,
                    "T": t,
                    "T_in_attr": t in (1, 2),
                    "s_in_attr": s in (1, 2),
                }
            )
            break  # erstes Digramm pro Seed
    if local_rows:
        near = [r for r in local_rows if r["d_before"] < 3]
        far = [r for r in local_rows if r["d_before"] >= 3]
        print(
            f"  erstes 10/01-Fenster: n={len(local_rows)}  "
            f"d_before<3: {len(near)} (davon T∈Attr "
            f"{sum(1 for r in near if r['T_in_attr'])})  "
            f"d_before≥3: {len(far)} (T∈Attr "
            f"{sum(1 for r in far if r['T_in_attr'])})"
        )

    print("\n[7] Fazit")
    print(
        "  χ<ε allein ⇒ nicht d_before < r  (Präfixe: d_before∼n₀).\n"
        "  Affine Identität: d_after = |A| d_before.\n"
        "  Eintritt korreliert mit d_after ≲ O(1), d.h. genug Kontraktion\n"
        "  |A|·n₀ ≲ O(1), nicht mit bloßer Fixpunktnähe χ.\n"
        "  Offene Hürde: Überdeckung durch Blockzylinder mit d_after < r."
    )

    out = Path(__file__).resolve().parent / "collatz_block_radius.json"
    payload = {
        "claim": (
            "Block radius spectrum: chi<eps does not imply d_before<r; "
            "entry correlates with d_after=|A|*d_before small; "
            "not a Collatz proof"
        ),
        "radius_atlas_pure": radius_atlas,
        "prefix_by_eps": {
            str(eps): {
                "n": len(rows_by_eps[eps]),
                "frac_naive_d_before_lt_2": (
                    sum(1 for r in rows_by_eps[eps] if r["naive_lemma"])
                    / max(len(rows_by_eps[eps]), 1)
                ),
                "frac_d_after_lt_2": (
                    sum(1 for r in rows_by_eps[eps] if r["after_lemma"])
                    / max(len(rows_by_eps[eps]), 1)
                ),
                "frac_T_in_attractor": (
                    sum(1 for r in rows_by_eps[eps] if r["T_in_attractor"])
                    / max(len(rows_by_eps[eps]), 1)
                ),
                "median_d_before": sorted(r["d_before"] for r in rows_by_eps[eps])[
                    len(rows_by_eps[eps]) // 2
                ]
                if rows_by_eps[eps]
                else None,
                "median_d_after": sorted(r["d_after"] for r in rows_by_eps[eps])[
                    len(rows_by_eps[eps]) // 2
                ]
                if rows_by_eps[eps]
                else None,
                "rows": rows_by_eps[eps],
            }
            for eps in epsilons
        },
        "bridge": {
            "note": "chi measures fixpoint class; d_after measures cylinder landing",
            "refined_conjecture": (
                "chi<eps and |A|*n0 = O(1) (equivalently d_after small) "
                "=> attractor entry for the prefix block"
            ),
        },
        "local_digram_sample": local_rows[:20],
        "hierarchy": "word -> (Lambda,x*) -> block -> radius/d_after -> attractor",
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
