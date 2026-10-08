"""
CollatzAffineAtlasClose — Schlussstein des affinen Collatz-Atlas.

Parameter:
  Λ_k = a_k log₂3 − k
  D(n,k) = −Λ_k − log₂ n     (Kompressionsreserve)
  R_k(n) = |A_k| |n − x_k*|
  k* = min{k: R_k ≤ 2} (= enter)
  Q(n) = k*/log₂ n
  S(n) = max_{k≤k*} |x_k*|

Vermutung (affine Covering-Form):
  ∃C ∀n ∃k≤C log n:
      −Λ_k ≥ log₂ n − C  und  |x_k*| ≤ C
  ⇒ R_k(n) < 2.

Epistemik: Atlas-Freeze; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
import random
from collections import defaultdict
from fractions import Fraction
from pathlib import Path


LOG2_3 = math.log2(3)


def path_to_enter(n0: int, max_steps: int = 200_000) -> dict | None:
    bits: list[int] = []
    x = n0
    for k in range(1, max_steps + 1):
        bit = 0 if x % 2 == 0 else 1
        bits.append(bit)
        x = x // 2 if bit == 0 else (3 * x + 1) // 2
        if x in (1, 2):
            return {"k_star": k, "bits": bits, "enter_state": x}
    return None


def prefix_series(n0: int, bits: list[int]) -> list[dict]:
    a = L = B = 0
    rows = []
    logn = math.log2(n0)
    for k, bit in enumerate(bits, start=1):
        if bit == 0:
            L += 1
        else:
            B = 3 * B + (1 << L)
            a += 1
            L += 1
        A = Fraction(3**a, 1 << L)
        x_star = Fraction(B, (1 << L) - 3**a)
        Lam = a * LOG2_3 - k  # L==k for Syracuse bit length
        R = float(abs(A) * abs(Fraction(n0) - x_star))
        D = -Lam - logn
        rows.append(
            {
                "k": k,
                "a": a,
                "Lambda": Lam,
                "A": float(A),
                "x_star": float(x_star),
                "abs_x": abs(float(x_star)),
                "R": R,
                "D": D,
                "neg_Lam": -Lam,
                "logn": logn,
            }
        )
    return rows


def analyze(n0: int) -> dict | None:
    hit = path_to_enter(n0)
    if hit is None:
        return None
    k_star = hit["k_star"]
    rows = prefix_series(n0, hit["bits"])
    assert len(rows) == k_star
    star = rows[-1]
    abs_x = [r["abs_x"] for r in rows]
    S = max(abs_x)
    S_at = next(r["k"] for r in rows if r["abs_x"] == S)
    # Kompression am Eintritt
    D_star = star["D"]
    # −Λ ≥ log n − C  ⇔  D ≥ −C; prüfe D_star
    return {
        "n0": n0,
        "log2_n": math.log2(n0),
        "k_star": k_star,
        "Q": k_star / math.log2(n0),
        "R_star": star["R"],
        "R_prev": rows[-2]["R"] if k_star > 1 else None,
        "Lambda_star": star["Lambda"],
        "neg_Lambda_star": -star["Lambda"],
        "D_star": D_star,
        "A_star": star["A"],
        "x_star": star["x_star"],
        "abs_x_star": star["abs_x"],
        "S": S,
        "S_at_k": S_at,
        "S_over_log": S / math.log2(n0),
        "max_abs_x_pre": max(abs_x[:-1]) if k_star > 1 else S,
        "enter_state": hit["enter_state"],
        # Checks für Covering-Form am k*
        "negLam_ge_logn": -star["Lambda"] >= math.log2(n0) - 5,  # soft C=5 probe
    }


def main() -> int:
    random.seed(23)

    print("=" * 72)
    print("AFFINER COLLATZ-ATLAS — SCHLUSSSTEIN  (S, D, Q)")
    print("=" * 72)
    print("S(n)=max|x_k*|;  D(n,k*)=−Λ_{k*}−log₂n;  Q=k*/log₂n\n")

    seeds: list[int] = list(range(3, 150))
    seeds += [2**k - 1 for k in range(8, 36)]
    seeds += [10**k + 3 for k in range(3, 12)]
    for e in range(3, 12):
        for _ in range(25):
            seeds.append(random.randint(10**e, 10 ** (e + 1) - 1))
    for _ in range(20):
        seeds.append(random.randint(10**12, 10**16))

    seen: set[int] = set()
    seeds = [s for s in seeds if s > 1 and s not in seen and not seen.add(s)]  # type: ignore

    rows = []
    for n0 in seeds:
        rec = analyze(n0)
        if rec:
            rows.append(rec)

    print(f"[1] Sample: {len(rows)} Bahnen bis Eintritt\n")

    def stats(key: str) -> dict:
        xs = [r[key] for r in rows]
        return {
            "min": min(xs),
            "median": sorted(xs)[len(xs) // 2],
            "max": max(xs),
            "mean": sum(xs) / len(xs),
        }

    for key, label in (
        ("Q", "Q(n)"),
        ("S", "S(n)=max|x*|"),
        ("S_over_log", "S/log₂n"),
        ("D_star", "D(n,k*)"),
        ("abs_x_star", "|x*|_{k*}"),
        ("R_star", "R_{k*}"),
    ):
        st = stats(key)
        print(
            f"  {label:16s}  min={st['min']:.4g}  med={st['median']:.4g}  "
            f"max={st['max']:.4g}  mean={st['mean']:.4g}"
        )

    # [2] Covering-Form am k*
    print("\n[2] Affine Covering-Form am k*")
    print("  Ziel: −Λ_{k*} ≥ log₂n − C  und  |x*| ≤ C'")
    for C in (1, 2, 3, 5, 8, 10):
        ok_lam = sum(
            1 for r in rows if r["neg_Lambda_star"] >= r["log2_n"] - C
        )
        print(f"  −Λ ≥ log n − {C}:  {ok_lam}/{len(rows)}")
    for Cp in (1, 2, 5, 10, 20, 50, 100, 1e3, 1e6):
        ok_x = sum(1 for r in rows if r["abs_x_star"] <= Cp)
        print(f"  |x*|_{{k*}} ≤ {Cp:g}:  {ok_x}/{len(rows)}")

    # S-Schranken
    print("\n[3] S(n)-Schranken")
    for B in (10, 50, 100, 1e3, 1e4, 1e6, 1e9, 1e12):
        fails = sum(1 for r in rows if r["S"] > B)
        print(f"  S≤{B:g}: fails {fails}/{len(rows)}")

    # [4] Buckets
    print("\n[4] S und D nach log10-Buckets")
    buckets: dict[int, list[dict]] = defaultdict(list)
    for r in rows:
        buckets[int(math.floor(math.log10(r["n0"])))].append(r)
    print("  bucket    #   S_med      S_max     D_med    Q_med")
    for b in sorted(buckets):
        g = buckets[b]
        Sm = sorted(r["S"] for r in g)
        Dm = sorted(r["D_star"] for r in g)
        Qm = sorted(r["Q"] for r in g)
        print(
            f"  1e{b:<2d}  {len(g):4d}  "
            f"{Sm[len(Sm)//2]:9.4g}  {Sm[-1]:9.4g}  "
            f"{Dm[len(Dm)//2]:7.3f}  {Qm[len(Qm)//2]:6.3f}"
        )

    # [5] Extremale S
    print("\n[5] Größte S(n) (Fixpunktverschiebung entlang Transient)")
    for r in sorted(rows, key=lambda d: -d["S"])[:12]:
        print(
            f"  S={r['S']:12.4g}  @k={r['S_at_k']:4d}/{r['k_star']:4d}  "
            f"|x*|_*={r['abs_x_star']:.4g}  D_*={r['D_star']:+.3f}  "
            f"Q={r['Q']:.2f}  n={r['n0']}"
        )

    # [6] Relation R ≈ |A| n wenn |x*|≪n
    print("\n[6] Konsistenz: R_* ≤ 2 und D_*")
    r_ok = sum(1 for r in rows if r["R_star"] <= 2 + 1e-9)
    r_prev = sum(
        1
        for r in rows
        if r["R_prev"] is None or r["R_prev"] >= 2 - 1e-9
    )
    print(f"  R_{{k*}}≤2: {r_ok}/{len(rows)}")
    print(f"  R_{{k*-1}}≥2: {r_prev}/{len(rows)}")
    # D_* typisch nahe 0?
    print(
        f"  D_* ∈ [-2,5]: "
        f"{sum(1 for r in rows if -2 <= r['D_star'] <= 5)}/{len(rows)}"
    )

    # [7] Freeze-Thesen
    print("\n[7] ATLAS-FREEZE — Thesen")
    print(
        """
  E1  Wortklassifikation:  x*∈{0,1,2} ⟺ stabiler Kanal (lokal).
  E2  Radius:  R_k=|A_k||n−x_k*| = d_after;  R_k≤2 ⟺ Eintritt.
  E3  Skalierung:  Q(n)=k*/log₂n typisch O(1)–O(15) im Scan.
  E4  Kompression:  D(n,k*)=−Λ_{k*}−log₂n = O(1) am Eintritt.
  E5  Fixpunktlage:  |x*|_{k*} = O(1);  S(n)=max_{k≤k*}|x_k*| kann
      während Transient groß sein — entscheidend ist Kontrolle bei k*.
  E6  Offene Hürde = universelle Covering-Aussage:
        ∀n ∃k≤C log n:  R_k(n)≤2
      äquivalent Collatz; äquivalente Affine-Form:
        −Λ_k ≥ log₂n−C  und  |x_k*|≤C'  (bei geeignetem C).
  E7  Kein Schritt-Lyapunov V(Λ,x*); Block/Radius/Covering statt Abstieg.
""".rstrip()
    )

    out = Path(__file__).resolve().parent / "collatz_affine_atlas_close.json"
    manifesto = Path(__file__).resolve().parent / "collatz_affine_atlas_freeze.json"

    payload = {
        "claim": (
            "Affine Collatz atlas freeze: R_k<=2 iff enter; "
            "Q=O(1)-O(15); D(n,k*)=O(1); |x*|_{k*}=O(1); "
            "S(n) may be large on transient; covering open; not a proof"
        ),
        "n_seeds": len(rows),
        "summary": {key: stats(key) for key in (
            "Q", "S", "S_over_log", "D_star", "abs_x_star", "R_star", "neg_Lambda_star"
        )},
        "covering_checks": {
            "R_star_le_2": r_ok,
            "R_prev_ge_2": r_prev,
            "D_star_in_[-2,5]": sum(1 for r in rows if -2 <= r["D_star"] <= 5),
            "abs_x_star_le_2": sum(1 for r in rows if r["abs_x_star"] <= 2),
            "abs_x_star_le_5": sum(1 for r in rows if r["abs_x_star"] <= 5),
            "S_le_1e6": sum(1 for r in rows if r["S"] <= 1e6),
            "S_le_1e12": sum(1 for r in rows if r["S"] <= 1e12),
        },
        "buckets": {
            str(b): {
                "n": len(buckets[b]),
                "S_median": sorted(r["S"] for r in buckets[b])[len(buckets[b]) // 2],
                "S_max": max(r["S"] for r in buckets[b]),
                "D_median": sorted(r["D_star"] for r in buckets[b])[
                    len(buckets[b]) // 2
                ],
                "Q_median": sorted(r["Q"] for r in buckets[b])[len(buckets[b]) // 2],
            }
            for b in sorted(buckets)
        },
        "worst_S": sorted(rows, key=lambda d: -d["S"])[:20],
        "worst_Q": sorted(rows, key=lambda d: -d["Q"])[:15],
        "theses": {
            "E1_local_channel": "x* in {0,1,2} iff stable channel",
            "E2_radius": "R_k<=2 iff enter",
            "E3_scaling": "Q(n)=k*/log2 n typically O(1)-O(15)",
            "E4_compression": "D(n,k*)=O(1) at entry",
            "E5_fixpoint": "|x*|_{k*}=O(1); S may spike on transient",
            "E6_open": "forall n exists k<=C log n: R_k(n)<=2",
            "E7_no_step_lyapunov": "block/radius/covering, not step V",
        },
        "rows": rows,
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")

    freeze = {
        "title": "Affine Collatz Operator Atlas — Freeze",
        "status": "local mechanism closed; global covering open",
        "pipeline": "word -> (A,x*) -> R_k(n) -> k* -> Q,D,S",
        "theses": payload["theses"],
        "summary": payload["summary"],
        "covering_checks": payload["covering_checks"],
        "scripts": [
            "CollatzFixpointSpectrum.py",
            "CollatzPathFixpointSpectrum.py",
            "CollatzFixpointDrift.py",
            "CollatzBlockSelector.py",
            "CollatzBlockRadius.py",
            "CollatzAffineRadiusCover.py",
            "CollatzRadiusScaling.py",
            "CollatzAffineAtlasClose.py",
        ],
        "open_equivalent_to_collatz": (
            "exists C forall n exists k<=C log2(n): R_k(n)<=2"
        ),
        "not_claimed": "proof of Collatz conjecture",
    }
    manifesto.write_text(json.dumps(freeze, indent=2) + "\n", encoding="utf-8")

    print(f"\nwrote {out}")
    print(f"wrote {manifesto}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
