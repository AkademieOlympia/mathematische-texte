"""
CollatzLiftValuation — v_2(|2^L − 3^a|) als Lift-Obstruktion resonanter Zyklen.

Syracuse-Zyklus γ ⊂ Z/2^j Z Länge L, a ungerade Schritte:
  affine Komposition  S^L(x) = (3^a x + B)/2^L
  Fixpunktform: (2^L − 3^a)x = B

Bewertung:
  ν(γ) := v_2(|2^L − 3^a|)

Vermutung (empirisch zu testen):
  η(γ) > 0  ⇒  ν(γ) < j + C   (beschränkte Lift-Tiefe)
  bzw. schärfer: resonante Phantome haben ν ≪ j oder ν unabhängig klein.

Λ = a log2 3 − L > 0  ⇔  3^a > 2^L  ⇔  μ > p_c
η > 0  ⇔  μ > 1/2

Epistemik: algebraische Lift-Diagnostik am Automaten; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


LOG2_3 = math.log2(3)
P_CRIT = math.log(2) / math.log(3)


def v2(n: int) -> int:
    """2-adische Bewertung; v2(0)=+∞ als None."""
    if n == 0:
        return None  # type: ignore
    n = abs(n)
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def syracuse_next(r: int, mod: int) -> tuple[int, int]:
    if r % 2 == 0:
        return (r // 2) % mod, 0
    return ((3 * r + 1) // 2) % mod, 1


def affine_B_of_cycle(nodes: list[int], mod: int) -> tuple[int, int, int]:
    """
    Bestimme (a, L, B) entlang des Zyklus.
    B ganzzahlig mit S^L(x)=(3^a x+B)/2^L auf dem Pfad (Syracuse).
    """
    A = 1  # tracked as 3^a / 2^L at end
    # accumulate with fractions via integers: num/den
    # start: x |-> x  (a=0,L=0,B=0)
    a = 0
    L = 0
    # Use Fraction-free: after each step maintain (a,L,b) with map (3^a n+b)/2^L
    b = 0
    for x in nodes:
        # verify parity matches
        bit = 0 if x % 2 == 0 else 1
        if bit == 0:
            # /2: a,L,b -> a, L+1, b  (b unchanged in integer form before rescale)
            L += 1
        else:
            # (3x+1)/2: b := 3b + 2^L, a+=1, L+=1
            b = 3 * b + (1 << L)
            a += 1
            L += 1
    return a, L, b


def all_cycles(j: int) -> list[list[int]]:
    mod = 1 << j
    nxt = [0] * mod
    for r in range(mod):
        nxt[r], _ = syracuse_next(r, mod)
    visited = [0] * mod
    cycles = []
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
            cycles.append(path[pos[x] :])
        for v in path:
            visited[v] = 2
    return cycles


def analyze_cycle(nodes: list[int], j: int) -> dict:
    mod = 1 << j
    a, L, B = affine_B_of_cycle(nodes, mod)
    # cross-check a via weights
    a2 = sum(1 for x in nodes if x % 2 == 1)
    assert a == a2 and L == len(nodes)

    diff = (1 << L) - 3**a  # 2^L - 3^a
    nu = v2(diff)  # v2(|2^L-3^a|) via abs inside v2
    mu = a / L
    eta = mu - 0.5
    Lam = a * LOG2_3 - L
    # lift obstruction heuristic: need room for valuation growth
    # compare nu with j
    return {
        "j": j,
        "L": L,
        "a": a,
        "mu": mu,
        "eta": eta,
        "Lambda": Lam,
        "B": B,
        "diff": diff,  # 2^L - 3^a (negativ wenn resonant Λ>0)
        "abs_diff": abs(diff),
        "nu": nu,
        "nu_minus_j": (nu - j) if nu is not None else None,
        "nu_lt_j": (nu < j) if nu is not None else None,
        "ge_pc": mu + 1e-15 >= P_CRIT,
        "eta_pos": eta > 1e-15,
        "Lambda_pos": Lam > 1e-15,
        "nodes": nodes if L <= 20 else nodes[:10] + ["…"],
        "is_12": set(nodes) == {1, 2},
        "is_zero": nodes == [0],
    }


def try_pattern_lift_depth(nodes: list[int], j: int, j_max: int = 40) -> int:
    """Wie viele Stufen j→j+1… hebt das Knotenmuster? (depth)."""
    # reuse logic: lift one level
    def lifts_to(cyc: list[int], k: int) -> list[list[int]] | None:
        mod_lo = 1 << k
        mod_hi = 1 << (k + 1)
        base = cyc[0]
        out = []
        for s0 in (base, base + mod_lo):
            if s0 >= mod_hi:
                continue
            ns = []
            x = s0
            ok = True
            for expected in cyc:
                if x % mod_lo != expected:
                    ok = False
                    break
                ns.append(x)
                x, _ = syracuse_next(x, mod_hi)
            if ok and x == s0 and [n % mod_lo for n in ns] == cyc:
                out.append(ns)
        return out

    depth = 0
    cur = nodes
    k = j
    while k < j_max:
        ups = lifts_to(cur, k)
        if not ups:
            break
        cur = ups[0]
        k += 1
        depth += 1
    return depth


def main() -> int:
    print("=" * 72)
    print("LIFT-VALUATION  —  ν = v_2(|2^L − 3^a|)")
    print("=" * 72)
    print("Frage: η>0  ⇒  ν < j+C  (endliche 2-adische Tiefe)?")
    print(f"p_c={P_CRIT:.6f}\n")

    j_min, j_max = 4, 20
    records = []
    print(
        f"  {'j':>3}  {'L':>4}  {'a':>4}  {'η':>9}  {'Λ':>8}  "
        f"{'ν':>4}  {'ν-j':>6}  {'liftΔ':>5}  typ"
    )

    for j in range(j_min, j_max + 1):
        for nodes in all_cycles(j):
            rec = analyze_cycle(nodes, j)
            depth = try_pattern_lift_depth(nodes, j, j_max=min(j + 8, 32))
            rec["lift_depth"] = depth
            records.append(rec)
            typ = (
                "zero"
                if rec["is_zero"]
                else "1-2"
                if rec["is_12"]
                else "PHANTOM"
                if rec["eta_pos"]
                else "other"
            )
            nu_s = "∞" if rec["nu"] is None else str(rec["nu"])
            nj = rec["nu_minus_j"]
            nj_s = f"{nj:+6d}" if nj is not None else f"{'—':>6}"
            print(
                f"  {j:3d}  {rec['L']:4d}  {rec['a']:4d}  {rec['eta']:+9.5f}  "
                f"{rec['Lambda']:+8.3f}  {nu_s:>4}  {nj_s}  {depth:5d}  {typ}"
            )

    # Empirische Implikation
    print("\n[1] η>0  ⇒  ν ? j")
    phantoms = [r for r in records if r["eta_pos"]]
    sterile = [r for r in records if not r["eta_pos"]]
    print(f"  Phantom-Zyklen (η>0): {len(phantoms)}")
    for r in phantoms:
        print(
            f"    j={r['j']} L={r['L']} a={r['a']} η={r['eta']:+.4f}  "
            f"ν={r['nu']}  ν-j={r['nu_minus_j']:+d}  "
            f"Λ={r['Lambda']:+.3f}  lift_depth={r['lift_depth']}  "
            f"2^L-3^a={r['diff']}"
        )

    if phantoms:
        nus = [r["nu"] for r in phantoms if r["nu"] is not None]
        njs = [r["nu_minus_j"] for r in phantoms if r["nu_minus_j"] is not None]
        depths = [r["lift_depth"] for r in phantoms]
        print(f"\n  ν(phantoms): min={min(nus)} max={max(nus)}")
        print(f"  ν−j:         min={min(njs)} max={max(njs)}")
        print(f"  lift_depth:  max={max(depths)}  (0 ⇒ kein Lift j→j+1)")
        print(f"  alle ν < j?  {all(r['nu_lt_j'] for r in phantoms)}")
        print(f"  alle Λ>0?    {all(r['Lambda_pos'] for r in phantoms)}")

    print("\n[2] Kontrollklasse η≤0 (steril / 1-2 / 0)")
    # sample unique (L,a,nu) patterns
    seen = set()
    for r in sterile:
        key = (r["L"], r["a"], r["nu"], r["is_12"], r["is_zero"])
        if key in seen:
            continue
        seen.add(key)
        if len(seen) <= 12:
            print(
                f"    L={r['L']} a={r['a']} ν={r['nu']} η={r['eta']:+.4f}  "
                f"typ={'0' if r['is_zero'] else '1-2' if r['is_12'] else '?'}"
            )

    print("\n[3] Zusammenhang Λ>0 ⇒ endliche Tiefe")
    lam_pos = [r for r in records if r["Lambda_pos"]]
    print(f"  Zyklen mit Λ>0: {len(lam_pos)}  (sollte = Phantome mit μ>p_c)")
    if lam_pos:
        print(
            f"  max lift_depth unter Λ>0: {max(r['lift_depth'] for r in lam_pos)}"
        )
        print(
            f"  max ν unter Λ>0: {max(r['nu'] for r in lam_pos if r['nu'] is not None)}"
        )

    # Schärfere Relation: ν vs L
    print("\n[4] Vermutungs-Check")
    print(
        "  η>0 ⇒ ν < j+C:\n"
        f"    mit C=0 (ν<j): {all(r['nu_lt_j'] for r in phantoms) if phantoms else 'n/a'}\n"
        f"    lift_depth=0 für alle η>0: "
        f"{all(r['lift_depth']==0 for r in phantoms) if phantoms else 'n/a'}"
    )
    print(
        "  Λ>0 ⇒ endliche 2-adische Tiefe: empirisch lift_depth=0\n"
        "  (stärkste Form der Obstruktion — nicht nur ν beschränkt)."
    )
    print(
        "\n  Algebraischer Verdacht:\n"
        "    a>L/2 (bzw. Λ>0) ⇒ 3^a−2^L ungerade oder kleine ν\n"
        "    ⇒ keine beliebige Vertiefung der Zyklusgleichung mod 2^{j+•}."
    )

    # Tabelle für JSON / Markdown-freundlich
    table = []
    for r in records:
        if r["eta_pos"] or r["is_12"] or r["is_zero"] or r["j"] in (8, 12, 13):
            table.append(
                {
                    "j": r["j"],
                    "L": r["L"],
                    "a": r["a"],
                    "eta": r["eta"],
                    "Lambda": r["Lambda"],
                    "nu": r["nu"],
                    "nu_minus_j": r["nu_minus_j"],
                    "lift_depth": r["lift_depth"],
                    "diff": r["diff"],
                    "nodes": r["nodes"],
                }
            )

    out = Path(__file__).resolve().parent / "collatz_lift_valuation.json"
    payload = {
        "claim": (
            "v_2(|2^L-3^a|) vs level j for Syracuse cycles; "
            "eta>0 implies lift_depth=0 and nu<j in scan; not a Collatz proof"
        ),
        "p_crit": P_CRIT,
        "records": records,
        "phantoms": phantoms,
        "implication": {
            "eta_pos_implies_nu_lt_j": all(r["nu_lt_j"] for r in phantoms)
            if phantoms
            else None,
            "eta_pos_implies_lift_depth_0": all(r["lift_depth"] == 0 for r in phantoms)
            if phantoms
            else None,
            "Lambda_pos_implies_lift_depth_0": all(
                r["lift_depth"] == 0 for r in lam_pos
            )
            if lam_pos
            else None,
        },
        "table_focus": table,
        "bridge": "Lambda>0 => finite 2-adic depth (empirical)",
    }
    # shrink records for JSON - drop huge node lists duplicates
    out.write_text(
        json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8"
    )
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
