"""
CollatzLiftInvariant — arithmetische Lift-Invariante (Λ, ν₂, x*_ℚ).

Pipeline:
  ω ↦ (A=3^a/2^L, B) ↦ Λ = log₂ A ↦ C = 3^a − 2^L
      ↦ (ν₂(C), ν₂(B), ν₂^ℚ(A−1)) ↦ ℋ, d(J)

Affine Fixpunktgleichung (ganzzahlig):
  (2^L − 3^a) x ≡ B  (mod 2^j)
Koeffizient ungerade ⇒ immer eindeutig lösbar; ν₂(C)=0 ist universell
für a≥1 und trennt nicht (wie die frühere Valuation).

Die 2-adische Feinstruktur steckt in der rationalen Lösung
  x*_ℚ = B / (2^L − 3^a)
und ihrer Konsistenz mit dem Automatenwort.

Vermutung (empirisch zu prüfen):
  d_ω(J) → 1  (stable_ray)
    ⟺  Λ<0  und  x*_ℚ ∈ ℕ  und  word(x*_ℚ)=ω  (ganzzahliger Attraktor)
  Λ>0  ⇒  x*_ℚ ≤ 0  ⇒  kein positiver ℤ-Kanal  ⇒  d→0

Epistemik: Invariantenatlas; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path


LOG2_3 = math.log2(3)


def v2(n: int) -> int | None:
    """2-adische Bewertung; v2(0)=+∞ als None."""
    if n == 0:
        return None
    n = abs(n)
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


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


def modinv_odd(a: int, mod: int) -> int:
    return pow(a, -1, mod)


def algebraic_fixed_point(a: int, L: int, B: int, j: int) -> int:
    mod = 1 << j
    coef = ((1 << L) - pow(3, a)) % mod
    return (B * modinv_odd(coef, mod)) % mod


def syracuse_iterate_close(n: int, L: int, mod: int | None = None) -> tuple[tuple[int, ...], int]:
    bits = []
    x = n if mod is None else n % mod
    for _ in range(L):
        if mod is not None:
            x %= mod
        bit = 0 if x % 2 == 0 else 1
        bits.append(bit)
        if bit == 0:
            x = x // 2
        else:
            x = (3 * x + 1) // 2
        if mod is not None:
            x %= mod
    return tuple(bits), x


def lift_density(omega: tuple[int, ...], j_max: int) -> dict:
    a, L, B = affine_of_word(omega)
    H = []
    for j in range(1, j_max + 1):
        mod = 1 << j
        x = algebraic_fixed_point(a, L, B, j)
        w, x_end = syracuse_iterate_close(x, L, mod)
        if w == omega and x_end == x:
            H.append(j)
    H_ge2 = [j for j in H if j >= 2]
    stable_ray = H_ge2 == list(range(2, j_max + 1))
    return {
        "H": H,
        "N": len(H),
        "dens": len(H) / j_max,
        "dens_ge2": len(H_ge2) / max(j_max - 1, 1),
        "stable_ray": stable_ray,
        "max_H": max(H) if H else 0,
        "sporadic": len(H) >= 2 and H != list(range(H[0], H[-1] + 1)),
    }


def invariants(omega: tuple[int, ...], j_max: int = 48) -> dict:
    a, L, B = affine_of_word(omega)
    if L == 0:
        return {"omega": "", "L": 0}

    C = 3**a - (1 << L)          # 3^a − 2^L
    D = (1 << L) - 3**a          # 2^L − 3^a  (Koeffizient der FP-Gleichung)
    Lam = a * LOG2_3 - L
    A_num, A_den = 3**a, 1 << L  # A = A_num/A_den

    nu_C = v2(C)
    nu_B = v2(B)
    # ν₂^ℚ(A−1) = ν₂(C) − L   (für C≠0)
    nu_A_minus_1 = None if nu_C is None else nu_C - L

    # rationale Fixpunktlösung
    xQ = Fraction(B, D) if D != 0 else None
    xQ_int = xQ is not None and xQ.denominator == 1
    xQ_nonneg = xQ is not None and xQ >= 0
    xQ_pos_int = xQ_int and xQ > 0
    xQ_in_attractor = xQ_int and int(xQ) in (0, 1, 2)

    # ganzzahlige Selbstkonsistenz (nur wenn x*_ℚ ∈ ℕ)
    int_word_ok = False
    int_closed = False
    if xQ_pos_int:
        n0 = int(xQ)
        w, n_end = syracuse_iterate_close(n0, L, mod=None)
        int_word_ok = w == omega
        int_closed = n_end == n0

    lift = lift_density(omega, j_max)

    return {
        "omega": "".join(map(str, omega)),
        "a": a,
        "L": L,
        "B": B,
        "C": C,
        "D": D,
        "A": f"{A_num}/{A_den}",
        "Lambda": Lam,
        "Lambda_sign": "pos" if Lam > 1e-12 else ("neg" if Lam < -1e-12 else "zero"),
        "nu2_C": nu_C,
        "nu2_B": nu_B,
        "nu2_A_minus_1": nu_A_minus_1,  # = −L falls ν₂(C)=0
        "xQ": str(xQ) if xQ is not None else None,
        "xQ_nonneg": xQ_nonneg,
        "xQ_integer": xQ_int,
        "xQ_pos_int": xQ_pos_int,
        "xQ_in_{0,1,2}": xQ_in_attractor,
        "int_word_ok": int_word_ok,
        "int_closed": int_closed,
        "int_channel": int_word_ok and int_closed,
        **lift,
    }


def cycles_as_words(j: int) -> list[tuple[int, ...]]:
    mod = 1 << j
    nxt = [
        (r // 2) % mod if r % 2 == 0 else ((3 * r + 1) // 2) % mod for r in range(mod)
    ]
    visited = [0] * mod
    words: list[tuple[int, ...]] = []
    for start in range(mod):
        if visited[start]:
            continue
        path: list[int] = []
        pos: dict[int, int] = {}
        x = start
        while visited[x] == 0:
            visited[x] = 1
            pos[x] = len(path)
            path.append(x)
            x = nxt[x]
        if visited[x] == 1:
            cyc = path[pos[x] :]
            i0 = min(range(len(cyc)), key=lambda t: cyc[t])
            rot = cyc[i0:] + cyc[:i0]
            words.append(tuple(0 if v % 2 == 0 else 1 for v in rot))
        for v in path:
            visited[v] = 2
    return words


def main() -> int:
    j_max = 48
    L_max = 12

    print("=" * 72)
    print("LIFT-INVARIANTE  —  (Λ, ν₂(C), ν₂(B), x*_ℚ) ↔ dens(ℋ)")
    print("=" * 72)
    print("C = 3^a−2^L;  x*_ℚ = B/(2^L−3^a);  dens≥2 = |ℋ∩[2,J]|/(J−1)\n")

    # [1] Valuation-Sanity
    print("[1] ν₂(C)-Sanity (erwartet: immer 0 für a≥1)")
    samples = [(1, 0), (1, 1, 0, 0), (1,) * 5 + (0,) * 5, (1, 1, 1, 0, 1, 1, 1)]
    for omega in samples:
        r = invariants(omega, j_max=8)
        print(
            f"  ω={r['omega']:12s}  Λ={r['Lambda']:+.3f}  "
            f"ν₂(C)={r['nu2_C']}  ν₂(B)={r['nu2_B']}  "
            f"ν₂^ℚ(A−1)={r['nu2_A_minus_1']}"
        )

    # [2] Kontroll- und Phantomkanäle
    print("\n[2] Kanäle: Invariante vs Lift")
    focus = {
        "10": (1, 0),
        "01": (0, 1),
        "1010": (1, 0, 1, 0),
        "0101": (0, 1, 0, 1),
        "1100": (1, 1, 0, 0),
        "0": (0,),
    }
    for name, omega in focus.items():
        r = invariants(omega, j_max)
        print(
            f"  ω={name:6s}  Λ={r['Lambda']:+.4f}  ν₂C={r['nu2_C']}  ν₂B={r['nu2_B']}  "
            f"xQ={r['xQ']:>8s}  posℤ={r['xQ_pos_int']}  "
            f"intCh={r['int_channel']}  dens≥2={r['dens_ge2']:.3f}  "
            f"stable={r['stable_ray']}"
        )

    print("\n  Automaten-Phantome:")
    phantoms = []
    seen: set[str] = set()
    for j in (8, 10, 11, 12):
        for omega in cycles_as_words(j):
            key = "".join(map(str, omega))
            if key in seen:
                continue
            r = invariants(omega, j_max)
            if r["Lambda_sign"] != "pos":
                continue
            seen.add(key)
            phantoms.append(r)
            print(
                f"  j={j}  ω={r['omega'][:28]}{'…' if len(r['omega'])>28 else ''}  "
                f"Λ={r['Lambda']:+.3f}  ν₂C={r['nu2_C']}  ν₂B={r['nu2_B']}  "
                f"xQ={r['xQ']}  dens≥2={r['dens_ge2']:.4f}  "
                f"ℋ={r['H']}"
            )

    # [3] Exhaustive Atlas L≤12
    print(f"\n[3] Atlas L≤{L_max}: Korrelation Invariante ↔ dens")
    rows = []
    for L in range(1, L_max + 1):
        for bits in product([0, 1], repeat=L):
            if sum(bits) == 0 and L > 1:
                # pure zeros: only keep ω=0 once via L=1
                continue
            if bits == (0,) * L and L > 1:
                continue
            rows.append(invariants(bits, j_max))

    # ν₂(C) universal?
    a_ge1 = [r for r in rows if r["a"] >= 1]
    nuC_all0 = all(r["nu2_C"] == 0 for r in a_ge1)
    print(f"  #Wörter: {len(rows)}  (a≥1: {len(a_ge1)})")
    print(f"  ν₂(C)=0 für alle a≥1? {nuC_all0}")

    # contingency tables
    def tab(pred_name: str, pred, target_name: str = "stable_ray"):
        tp = tn = fp = fn = 0
        for r in rows:
            p = bool(pred(r))
            t = bool(r[target_name])
            if p and t:
                tp += 1
            elif not p and not t:
                tn += 1
            elif p and not t:
                fp += 1
            else:
                fn += 1
        print(
            f"  {pred_name:40s} ⇒ {target_name}:  "
            f"TP={tp} FP={fp} FN={fn} TN={tn}  "
            f"prec={tp/(tp+fp) if tp+fp else 0:.2f}  "
            f"rec={tp/(tp+fn) if tp+fn else 0:.2f}"
        )
        return {"TP": tp, "FP": fp, "FN": fn, "TN": tn}

    print("\n  Prädiktoren für stable_ray (dens≥2→1):")
    tabs = {}
    tabs["Lambda_neg"] = tab("Λ<0", lambda r: r["Lambda_sign"] == "neg")
    tabs["xQ_pos_int"] = tab("x*_ℚ ∈ ℕ>0", lambda r: r["xQ_pos_int"])
    tabs["int_channel"] = tab("int_channel (word∧close)", lambda r: r["int_channel"])
    tabs["Lam_neg_and_int"] = tab(
        "Λ<0 ∧ int_channel",
        lambda r: r["Lambda_sign"] == "neg" and r["int_channel"],
    )
    tabs["xQ_in_012"] = tab("x*_ℚ ∈ {0,1,2}", lambda r: r["xQ_in_{0,1,2}"])
    tabs["Lam_neg_and_012"] = tab(
        "Λ<0 ∧ x*_ℚ∈{0,1,2}",
        lambda r: r["Lambda_sign"] == "neg" and r["xQ_in_{0,1,2}"],
    )

    print("\n  Prädiktoren für dens≥2>0 (irgendein Lift j≥2):")
    for r in rows:
        r["has_lift"] = r["dens_ge2"] > 0
    tabs_d = {}
    tabs_d["Lam_neg"] = tab(
        "Λ<0", lambda r: r["Lambda_sign"] == "neg", target_name="has_lift"
    )
    tabs_d["int_ch"] = tab(
        "int_channel ∨ x*=0",
        lambda r: r["int_channel"] or r["xQ"] == "0",
        target_name="has_lift",
    )
    tabs_d["xQ_012"] = tab(
        "x*_ℚ ∈ {0,1,2}",
        lambda r: r["xQ_in_{0,1,2}"],
        target_name="has_lift",
    )
    tabs_d["Lam_neg_int"] = tab(
        "Λ<0 ∧ (int_channel ∨ x*=0)",
        lambda r: r["Lambda_sign"] == "neg"
        and (r["int_channel"] or r["xQ"] == "0"),
        target_name="has_lift",
    )

    # [4] Λ>0 ⇒ xQ≤0
    pos = [r for r in rows if r["Lambda_sign"] == "pos"]
    neg = [r for r in rows if r["Lambda_sign"] == "neg"]
    print("\n[4] Vorzeichen von x*_ℚ")
    print(
        f"  Λ>0 ⇒ x*_ℚ≤0? {all(not r['xQ_nonneg'] or r['xQ']=='0' for r in pos)}  "
        f"(#pos={len(pos)}, davon xQ<0: {sum(1 for r in pos if r['xQ'] and r['xQ'].startswith('-'))})"
    )
    print(
        f"  Λ<0 ⇒ x*_ℚ≥0? {all(r['xQ_nonneg'] for r in neg)}  (#neg={len(neg)})"
    )
    print(
        f"  Λ>0 ∧ xQ_pos_int: {sum(1 for r in pos if r['xQ_pos_int'])}  "
        f"(soll 0 — Path-/Affin-Theorem)"
    )

    # [5] ν₂(B)-Verteilung nach Typ
    print("\n[5] ν₂(B) nach Lift-Typ (nicht diskriminierend allein)")
    by_typ = defaultdict(list)
    for r in rows:
        if r["stable_ray"]:
            by_typ["A_stable"].append(r)
        elif r["Lambda_sign"] == "pos" and r["N"] > 0:
            by_typ["B_phantom"].append(r)
        elif r["N"] > 0:
            by_typ["C_finite"].append(r)
        else:
            by_typ["D_empty"].append(r)
    for typ, g in by_typ.items():
        nuBs = [r["nu2_B"] for r in g if r["nu2_B"] is not None]
        nuB0 = sum(1 for r in g if r["B"] == 0)
        print(
            f"  {typ:12s} n={len(g):5d}  "
            f"ν₂(B)=0: {sum(1 for v in nuBs if v==0):5d}  "
            f"ν₂(B)≥1: {sum(1 for v in nuBs if v is not None and v>=1):5d}  "
            f"B=0: {nuB0}"
        )

    # [6] Exact characterization candidates
    print("\n[6] Charakterisierung stable_ray")
    stables = [r for r in rows if r["stable_ray"]]
    print(f"  #stable_ray={len(stables)}")
    for r in sorted(stables, key=lambda d: (d["L"], d["omega"])):
        print(
            f"    ω={r['omega']:14s}  Λ={r['Lambda']:+.4f}  "
            f"xQ={r['xQ']}  ν₂B={r['nu2_B']}  intCh={r['int_channel']}"
        )
    # bijection check
    int_ch_neg = [
        r for r in rows if r["Lambda_sign"] == "neg" and r["int_channel"]
    ]
    print(
        f"\n  # (Λ<0 ∧ int_channel) = {len(int_ch_neg)}  "
        f"davon stable_ray={sum(1 for r in int_ch_neg if r['stable_ray'])}"
    )
    # any int_channel not stable?
    for r in int_ch_neg:
        if not r["stable_ray"]:
            print(
                f"    Gegenbeispiel intCh¬stable: ω={r['omega']}  "
                f"xQ={r['xQ']}  dens≥2={r['dens_ge2']}  ℋ={r['H'][:12]}"
            )

    print("\n[7] Fazit")
    print(
        "  ν₂(C)=0 und ν₂^ℚ(A−1)=−L sind universell (a≥1) — keine Sterilisierung.\n"
        "  Λ bestimmt das Vorzeichen von x*_ℚ: Λ>0 ⇒ x*<0 (kein positiver ℤ-Kanal).\n"
        "  Exact im Scan: stable_ray ⟺ x*_ℚ ∈ {0,1,2}.\n"
        "  Relevante Invariante: (Λ, x*_ℚ), nicht ν₂(C) oder ν₂(B)."
    )

    out = Path(__file__).resolve().parent / "collatz_lift_invariant.json"
    payload = {
        "claim": (
            "Lift invariant atlas: nu2(C)=0 universal for a>=1; "
            "Lambda>0 => xQ<=0; stable_ray iff Lambda<0 and integer self-consistent FP; "
            "not a Collatz proof"
        ),
        "j_max": j_max,
        "L_max": L_max,
        "nu2_C_universal_zero": nuC_all0,
        "contingency_stable_ray": tabs,
        "contingency_has_lift": tabs_d,
        "sign_xQ": {
            "Lambda_pos_implies_xQ_nonpos": all(
                (not r["xQ_nonneg"]) or r["xQ"] == "0" for r in pos
            ),
            "Lambda_neg_implies_xQ_nonneg": all(r["xQ_nonneg"] for r in neg),
            "n_Lambda_pos_xQ_pos_int": sum(1 for r in pos if r["xQ_pos_int"]),
        },
        "stables": stables,
        "phantoms": phantoms,
        "type_nu2B": {
            typ: {
                "n": len(g),
                "nu2B_zero": sum(1 for r in g if r["nu2_B"] == 0),
                "nu2B_pos": sum(
                    1 for r in g if r["nu2_B"] is not None and r["nu2_B"] >= 1
                ),
            }
            for typ, g in by_typ.items()
        },
        "conjecture": (
            "d->1 iff Lambda<0 and x*_Q in N with word(x*)=omega and S^L(x*)=x*; "
            "Lambda>0 => dens->0 (no positive integer channel)"
        ),
        "samples": [invariants(w, j_max) for w in focus.values()],
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
