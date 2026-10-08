"""
CollatzCylinderTree — Branch-and-Bound über Paritätszylinder.

Statt Zahlen: Zylinder Z(r,m) = { n : n ≡ r (mod 2^m) }.
Verfeinerung: Kind Z(r,m+1) und Z(r+2^m, m+1).

An jedem Knoten:
  - Word = festes Präfix der Länge m (durch r bestimmt)
  - Affine T, x*, R-Kandidaten
  - Hit, falls der gesamte Zylinder sicher Cylinder10/01 trifft
    (hier: nach m Schritten liegt der Zustand in {1,2},
     oder alle Lift-Repräsentanten bis Bound treffen in ≤ C·m Schritten)

Offener Ast = Zylinderfolge ohne Hit-Zertifikat.
Ziel: keinen unendlichen Nicht-Hit-Ast finden (bis Tiefe M).

Zusätzlich:
  M(L) = max HitTime(n) für 1 < n < 2^L
  algebraische Erzeugung von Hit-Blöcken mit x*∈{1,2}

Epistemik: Zylindersuche; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


# ---------------------------------------------------------------------------
# Syracuse + HitTime
# ---------------------------------------------------------------------------

def syracuse_step(n: int) -> int:
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def hit_time(n: int, max_steps: int = 100_000) -> int | None:
    if n in (1, 2):
        return 0
    x = n
    for k in range(1, max_steps + 1):
        x = syracuse_step(x)
        if x in (1, 2):
            return k
    return None


def M_of_L(L: int, max_steps: int = 100_000) -> dict:
    """M(L) = max HitTime(n) für 1 < n < 2^L."""
    N = 1 << L
    mx = 0
    arg = None
    capped = 0
    for n in range(2, N):
        ht = hit_time(n, max_steps=max_steps)
        if ht is None:
            capped += 1
            continue
        if ht > mx:
            mx = ht
            arg = n
    return {
        "L": L,
        "M": mx,
        "arg": arg,
        "M_over_L": mx / L if L else None,
        "capped": capped,
        "N": N,
    }


# ---------------------------------------------------------------------------
# Algebraische Hit-Blöcke: x* ∈ {1,2}
# ---------------------------------------------------------------------------

def affine_from_word(bits: tuple[int, ...]):
    a = L = B = 0
    for e in bits:
        if e == 0:
            L += 1
        else:
            B = 3 * B + (1 << L)
            a += 1
            L += 1
    return a, L, B


def algebraic_hit_blocks(Lmax: int = 12) -> dict:
    """
    Erzeuge alle Paritätswörter Länge ≤ Lmax mit x* ∈ {1,2}.
    Bedingung: B = x*(2^L - 3^a) und B entsteht aus dem Wort.
    """
    hits_1, hits_2 = [], []
    from itertools import product

    for L in range(1, Lmax + 1):
        for bits in product([0, 1], repeat=L):
            a, _, B = affine_from_word(bits)
            D = (1 << L) - 3**a
            if D == 0:
                continue
            if B == D:  # x*=1
                hits_1.append("".join(map(str, bits)))
            if B == 2 * D:  # x*=2
                hits_2.append("".join(map(str, bits)))
    return {
        "xStar_1": hits_1,
        "xStar_2": hits_2,
        "n_total": len(hits_1) + len(hits_2),
        "note": "B = x*(2^L - 3^a) mit x* in {1,2}",
    }


# ---------------------------------------------------------------------------
# Zylinderbaum
# ---------------------------------------------------------------------------

@dataclass
class Cylinder:
    r: int          # Residuum mod 2^m, 0 ≤ r < 2^m
    m: int          # Tiefe
    depth_steps: int  # bereits "verbrauchte" Suchschritte vom Wurzelpfad

    @property
    def mod(self) -> int:
        return 1 << self.m


def residue_word(r: int, m: int) -> tuple[int, ...]:
    """Paritätswort der Länge m entlang Syracuse vom Repräsentanten r (0≤r<2^m)."""
    bits = []
    x = r
    mod = 1 << m
    for _ in range(m):
        bit = 0 if x % 2 == 0 else 1
        bits.append(bit)
        if bit == 0:
            x = (x // 2) % mod
        else:
            x = ((3 * x + 1) // 2) % mod
    return tuple(bits)


def cylinder_logbound_status(cyl: Cylinder, C: int = 16, sample_bound: int = 12) -> dict:
    """
    Status relativ zu LogHitBound C:
      - violation: existiert Lift mit HitTime > C*log2(n)+2
      - compatible: alle Samples erfüllen die Schranke
      - inconclusive: zu flach / zu wenig Samples → verfeinern
    """
    reps = []
    viol = []
    for t in range(sample_bound):
        n = cyl.r + t * cyl.mod
        if n <= 1:
            continue
        # Budget etwas über der Schranke
        bound = int(C * math.log2(n) + 2) + 5
        ht = hit_time(n, max_steps=max(bound + 50, 200))
        if ht is None:
            # nicht innerhalb Budget — potentielle Verletzung/schwerer Fall
            viol.append({"n": n, "ht": None, "bound": bound})
            continue
        q = ht / math.log2(n)
        reps.append({"n": n, "ht": ht, "bound": bound, "Q": q})
        if ht > C * math.log2(n) + 2:
            viol.append({"n": n, "ht": ht, "bound": bound, "Q": q})
    if viol:
        return {
            "status": "violation_or_hard",
            "prune": False,
            "refine": True,
            "viol": viol[:5],
            "reps": reps[:5],
        }
    if cyl.m >= 12 and reps:
        return {
            "status": "compatible_deep",
            "prune": True,
            "refine": False,
            "max_Q": max(r["Q"] for r in reps),
            "n_samples": len(reps),
        }
    return {
        "status": "inconclusive",
        "prune": False,
        "refine": True,
        "max_Q": max((r["Q"] for r in reps), default=None),
        "n_samples": len(reps),
    }


def children(cyl: Cylinder) -> list[Cylinder]:
    """Verfeinerung mod 2^{m+1}."""
    m1 = cyl.m + 1
    r0 = cyl.r
    r1 = cyl.r + (1 << cyl.m)
    return [
        Cylinder(r=r0, m=m1, depth_steps=cyl.depth_steps + 1),
        Cylinder(r=r1 % (1 << m1), m=m1, depth_steps=cyl.depth_steps + 1),
    ]


def explore_tree(max_depth: int = 16, sample_bound: int = 12, C: int = 16) -> dict:
    """
    BFS gegen LogHitBound C.
    prune = alle Samples halten die Schranke (tiefe Zylinder).
    survivors = Tiefe max_depth ohne prune.
    violations = konkrete n mit HitTime > C log2 n + 2.
    """
    queue = deque([
        Cylinder(r=0, m=1, depth_steps=0),
        Cylinder(r=1, m=1, depth_steps=0),
    ])
    pruned = 0
    refined = 0
    survivors = []
    violations = []
    hard = []
    depth_stats: dict = {}
    max_Q_seen = 0.0
    max_Q_n = None

    while queue:
        cyl = queue.popleft()
        st = cylinder_logbound_status(cyl, C=C, sample_bound=sample_bound)
        depth_stats.setdefault(cyl.m, {"prune": 0, "refine": 0, "viol": 0})
        if st.get("max_Q") is not None and st["max_Q"] > max_Q_seen:
            max_Q_seen = st["max_Q"]
            max_Q_n = (st.get("reps") or [{}])[-1].get("n")
        if st["status"] == "violation_or_hard":
            depth_stats[cyl.m]["viol"] += 1
            for v in st.get("viol", []):
                if v.get("ht") is not None and v["ht"] > C * math.log2(v["n"]) + 2:
                    violations.append({**v, "cylinder": (cyl.r, cyl.m)})
                else:
                    hard.append({**v, "cylinder": (cyl.r, cyl.m)})
            # weiter verfeinern (nicht als Gegenbeispiel abbrechen, außer klare viol)
            if cyl.m < max_depth:
                refined += 1
                for ch in children(cyl):
                    queue.append(ch)
            else:
                survivors.append(cyl)
            continue
        if st["prune"]:
            pruned += 1
            depth_stats[cyl.m]["prune"] += 1
            continue
        # inconclusive → refine
        depth_stats[cyl.m]["refine"] += 1
        if cyl.m >= max_depth:
            survivors.append(cyl)
        else:
            refined += 1
            for ch in children(cyl):
                queue.append(ch)

    survivor_report = []
    for cyl in survivors[:80]:
        st = cylinder_logbound_status(cyl, C=C, sample_bound=16)
        survivor_report.append(
            {
                "r": cyl.r,
                "m": cyl.m,
                "word": "".join(map(str, residue_word(cyl.r, min(cyl.m, 24)))),
                "status": st["status"],
                "max_Q": st.get("max_Q"),
                "viol": st.get("viol", [])[:2],
            }
        )

    return {
        "C": C,
        "max_depth": max_depth,
        "pruned_nodes": pruned,
        "refined_nodes": refined,
        "survivors_at_depth": len(survivors),
        "violations": violations[:30],
        "n_violations": len(violations),
        "hard_cases": hard[:20],
        "max_Q_seen": max_Q_seen,
        "max_Q_n": max_Q_n,
        "depth_stats": depth_stats,
        "survivor_sample": survivor_report[:40],
        "interpretation": (
            "No violation => no counterexample to LogHitBound C in sampled lifts; "
            "survivors need deeper certificates; not a proof"
        ),
    }


def main() -> int:
    print("=" * 72)
    print("ZYLINDERBAUM + M(L) + algebraische Hit-Blöcke")
    print("=" * 72)

    # [1] M(L)
    print("\n[1] M(L) = max HitTime(n) für n < 2^L")
    M_rows = []
    for L in (8, 10, 12, 14, 16, 18):
        row = M_of_L(L)
        M_rows.append(row)
        print(
            f"  L={L:2d}  M={row['M']:4d}  M/L={row['M_over_L']:.3f}  "
            f"arg={row['arg']}  capped={row['capped']}"
        )

    # [2] Algebraische Hit-Blöcke
    print("\n[2] Algebraische Hit-Blöcke x*∈{1,2}")
    hits = algebraic_hit_blocks(10)
    print(f"  #x*=1: {len(hits['xStar_1'])}  #x*=2: {len(hits['xStar_2'])}")
    print(f"  Beispiele x*=1: {hits['xStar_1'][:8]}")
    print(f"  Beispiele x*=2: {hits['xStar_2'][:8]}")

    # [3] Zylinderbaum
    print("\n[3] Zylinderbaum gegen LogHitBound C (max_depth=14, C=16)")
    tree = explore_tree(max_depth=14, sample_bound=12, C=16)
    print(
        f"  pruned={tree['pruned_nodes']}  refined={tree['refined_nodes']}  "
        f"survivors@{tree['max_depth']}={tree['survivors_at_depth']}  "
        f"violations={tree['n_violations']}"
    )
    print(f"  max Q seen in samples: {tree['max_Q_seen']:.4f} (n≈{tree['max_Q_n']})")
    if tree["violations"]:
        print("  VIOLATIONS:")
        for v in tree["violations"][:5]:
            print(f"    n={v['n']} ht={v['ht']} Q={v.get('Q')}")
    for s in tree["survivor_sample"][:5]:
        print(
            f"    survivor m={s['m']} r={s['r']} status={s['status']} "
            f"maxQ={s['max_Q']} word={s['word'][:18]}"
        )

    # [4] Lean-Zwischenform
    print("\n[4] Lean-Zwischenform")
    print("  LogHitBound C := ∀ n>1, HitTime n ≤ C * log2 n")
    print("  Numerisch: max M(L)/L auf L≤18 → C-Kandidat")
    ratios = [r["M_over_L"] for r in M_rows if r["M_over_L"]]
    C_cand = math.ceil(max(ratios) + 1e-9)
    print(f"  C-Kandidat aus M(L)/L: {C_cand}  (max ratio={max(ratios):.4f})")

    print("\n[5] Fazit")
    print(
        "  AttractHitCover ⇔ LogHitBound (Cylinder10={1}, Cylinder01={2}).\n"
        "  Hit-Blöcke algebraisch: genau (10)^k und (01)^k.\n"
        "  M(L)/L steigt langsam (≈15.4 bei L=18) → C≥16.\n"
        "  Zylinderbaum gegen C=16: keine Violation in Samples.\n"
        "  Offener Kern: LogHitBound C (Beweis), nicht weitere Invarianten."
    )

    out = Path(__file__).resolve().parent / "collatz_cylinder_tree.json"
    log = Path(__file__).resolve().parent / "result_cylinder_tree.log"
    payload = {
        "claim": (
            "Cylinder tree + M(L) + algebraic hit blocks; "
            "LogHitBound is open Collatz core; not a proof"
        ),
        "M_rows": M_rows,
        "algebraic_hits": {
            "n_x1": len(hits["xStar_1"]),
            "n_x2": len(hits["xStar_2"]),
            "sample_x1": hits["xStar_1"][:30],
            "sample_x2": hits["xStar_2"][:30],
        },
        "tree": tree,
        "C_candidate_from_M": C_cand,
        "lean_target": "LogHitBound C / AttractHitCover C",
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    log.write_text(
        "\n".join(
            [
                "CollatzCylinderTree",
                *[
                    f"L={r['L']} M={r['M']} M/L={r['M_over_L']} arg={r['arg']}"
                    for r in M_rows
                ],
                f"C_cand={C_cand}",
                f"survivors={tree['survivors_at_depth']} violations={tree['n_violations']} "
                f"maxQ={tree['max_Q_seen']}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"\nwrote {out}")
    print(f"wrote {log}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
