"""
CollatzCoverGrowth — Wachstumsrate N_B(L) des skalenabhängigen Covers.

Cover(L): für alle ungeraden n < 2^L existiert ω mit
  |ω| ≤ C·L  und  R_ω ≤ 2  (Operatoratlas).

Drei konstruktive Atlanten (obere Schranken an N_B):

  1) Prefix-Atlas:
       ω = Bahnpräfix Länge k*(n)  (First Passage R≤2)
       N_B^pref(L) = #{ solche ω }

  2) Hit-Block-Atlas (kurz):
       ω = erster kurzer Operator-Treffer aus Bibliothek Länge ≤ L_blk
       N_B^hit(L) = #{ getroffene Bibliotheksblöcke }

  3) Attraktor-Atlas χ=0:
       ω ∈ {01,10,(01)^m,(10)^m,…} Länge ≤ L_blk
       N_B^attr(L) = #{ genutzte Attraktorwörter }

Szenarien:
  Fall 1: N_B = O(1)
  Fall 2: N_B = O(L^c)
  Fall 3: N_B ∼ 2^L

Zusätzlich: max k*(n) / L für n < 2^L.

Epistemik: Cover-Wachstum; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
import time
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path


def affine(omega: tuple[int, ...]):
    a = L = B = 0
    for e in omega:
        if e == 0:
            L += 1
        else:
            B = 3 * B + (1 << L)
            a += 1
            L += 1
    A = Fraction(3**a, 1 << L) if L else Fraction(1)
    D = (1 << L) - 3**a
    x = Fraction(B, D) if D else Fraction(0)
    return a, L, B, A, x


def radius_prefix(n0: int, bits: list[int]) -> float:
    a = L = B = 0
    for bit in bits:
        if bit == 0:
            L += 1
        else:
            B = 3 * B + (1 << L)
            a += 1
            L += 1
    A = Fraction(3**a, 1 << L)
    x = Fraction(B, (1 << L) - 3**a)
    return float(abs(A * (Fraction(n0) - x)))


def path_with_post(n0: int, post_enter: int = 24, max_steps: int = 200_000):
    bits: list[int] = []
    states = [n0]
    x = n0
    enter = None
    for _ in range(max_steps):
        b = 0 if x % 2 == 0 else 1
        bits.append(b)
        x = x // 2 if b == 0 else (3 * x + 1) // 2
        states.append(x)
        if enter is None and x in (1, 2):
            enter = len(bits)
            for _ in range(post_enter):
                b2 = 0 if x % 2 == 0 else 1
                bits.append(b2)
                x = x // 2 if b2 == 0 else (3 * x + 1) // 2
                states.append(x)
            break
    return bits, states, enter


def k_star_prefix(n0: int, bits: list[int]) -> int | None:
    """Min k mit R(prefix_k) ≤ 2."""
    for k in range(1, len(bits) + 1):
        if radius_prefix(n0, bits[:k]) <= 2 + 1e-12:
            return k
    return None


def build_attr_library(Lmax: int) -> list[tuple[int, ...]]:
    lib = []
    for m in range(1, Lmax // 2 + 1):
        lib.append(tuple([1, 0] * m))
        lib.append(tuple([0, 1] * m))
    # unique
    seen = set()
    out = []
    for w in lib:
        if w not in seen and len(w) <= Lmax:
            seen.add(w)
            out.append(w)
    return out


def build_contract_library(Lmax: int, x_bound: float = 50.0) -> list[tuple[int, ...]]:
    lib = []
    for L in range(2, Lmax + 1):
        for bits in product([0, 1], repeat=L):
            if sum(bits) == 0:
                continue
            _, _, _, A, x = affine(bits)
            if float(A) >= 1:
                continue
            if abs(float(x)) > x_bound:
                continue
            lib.append(bits)
    return lib


def first_hit_block(
    bits: list[int],
    states: list[int],
    library: list[tuple[int, ...]],
) -> tuple[int, tuple[int, ...]] | None:
    by_L: dict[int, set[tuple[int, ...]]] = {}
    for w in library:
        by_L.setdefault(len(w), set()).add(w)
    lengths = sorted(by_L)
    for i in range(len(bits)):
        for L in lengths:
            if i + L > len(bits):
                continue
            seg = tuple(bits[i : i + L])
            if seg not in by_L[L]:
                continue
            s = states[i]
            _, _, _, A, x = affine(seg)
            R = abs(float(A * (Fraction(s) - x)))
            if R <= 2 + 1e-12:
                return i, seg
    return None


def scan_level(L: int, C_len: int = 20) -> dict:
    """Alle ungeraden n < 2^L."""
    N = 1 << L
    L_blk = min(max(C_len, 2), max(2 * L, 8))  # Blocklänge-Budget ~ O(L)
    attr_lib = build_attr_library(min(L_blk, 16))
    # kontraktive Lib nur bis begrenzter Tiefe (sonst Explosion)
    contr_L = min(L_blk, 8)
    contr_lib = build_contract_library(contr_L)

    pref_words: set[str] = set()
    attr_used: set[str] = set()
    contr_used: set[str] = set()
    k_stars = []
    Qs = []
    misses_pref = 0
    misses_attr = 0
    misses_contr = 0
    t0 = time.time()
    n_tested = 0

    for n in range(3, N, 2):
        n_tested += 1
        bits, states, enter = path_with_post(n)
        # 1) Prefix-Atlas
        ks = k_star_prefix(n, bits)
        if ks is None:
            misses_pref += 1
        else:
            pref_words.add("".join(map(str, bits[:ks])))
            k_stars.append(ks)
            Qs.append(ks / L)  # k*/L da log2 n < L für n<2^L

        # 2) Attraktor-Hit
        hit_a = first_hit_block(bits, states, attr_lib)
        if hit_a is None:
            misses_attr += 1
        else:
            attr_used.add("".join(map(str, hit_a[1])))

        # 3) Kontraktions-Hit
        hit_c = first_hit_block(bits, states, contr_lib)
        if hit_c is None:
            misses_contr += 1
        else:
            contr_used.add("".join(map(str, hit_c[1])))

    elapsed = time.time() - t0
    n_odds = n_tested

    def growth_class(nb: int) -> str:
        if nb <= 20:
            return "O(1)-ish"
        # vergleich mit L^c und 2^L
        if nb <= 5 * (L**3):
            return "O(L^c)-plausible"
        if nb >= 0.1 * n_odds:
            return "~2^L-ish"
        return "intermediate"

    return {
        "L": L,
        "N": N,
        "n_tested": n_odds,
        "L_blk": L_blk,
        "contr_L": contr_L,
        "N_B_prefix": len(pref_words),
        "N_B_attr": len(attr_used),
        "N_B_contract": len(contr_used),
        "N_B_prefix_over_odds": len(pref_words) / max(n_odds, 1),
        "miss_prefix": misses_pref,
        "miss_attr": misses_attr,
        "miss_contract": misses_contr,
        "max_k_star": max(k_stars) if k_stars else None,
        "max_k_star_over_L": (max(k_stars) / L) if k_stars else None,
        "med_k_star_over_L": (
            sorted(Qs)[len(Qs) // 2] if Qs else None
        ),
        "class_prefix": growth_class(len(pref_words)),
        "class_attr": growth_class(len(attr_used)),
        "class_contract": growth_class(len(contr_used)),
        "attr_lib_size": len(attr_lib),
        "contr_lib_size": len(contr_lib),
        "top_attr": Counter(attr_used).most_common(8),
        "elapsed_s": elapsed,
    }


def main() -> int:
    # L=8..16: 2^16=65536, odd ~32k — machbar
    levels = [8, 10, 12, 14, 16]
    print("=" * 72)
    print("COVER-WACHSTUM  —  N_B(L) für n < 2^L")
    print("=" * 72)
    print("Fall 1: O(1)  |  Fall 2: O(L^c)  |  Fall 3: ∼2^L\n")

    rows = []
    for L in levels:
        rec = scan_level(L)
        rows.append(rec)
        print(
            f"L={L:2d}  n_odd={rec['n_tested']:6d}  "
            f"N_pref={rec['N_B_prefix']:6d} ({rec['class_prefix']})  "
            f"N_attr={rec['N_B_attr']:3d} ({rec['class_attr']})  "
            f"N_contr={rec['N_B_contract']:4d} ({rec['class_contract']})  "
            f"max k*/L={rec['max_k_star_over_L']:.3f}  "
            f"miss_attr={rec['miss_attr']}  t={rec['elapsed_s']:.1f}s"
        )

    print("\n[Wachstum N_B(L)]")
    print("  L   N_pref   ratio_pref   N_attr  N_contr  max_k*/L")
    prev_p = None
    for r in rows:
        ratio = (
            r["N_B_prefix"] / prev_p if prev_p else float("nan")
        )
        print(
            f"  {r['L']:2d}  {r['N_B_prefix']:7d}  "
            f"{ratio if prev_p else float('nan'):10.3f}  "
            f"{r['N_B_attr']:6d}  {r['N_B_contract']:7d}  "
            f"{r['max_k_star_over_L']:8.3f}"
        )
        prev_p = r["N_B_prefix"]

    # Fit-Hinweise
    print("\n[Interpretation]")
    attr_sizes = [r["N_B_attr"] for r in rows]
    pref_sizes = [r["N_B_prefix"] for r in rows]
    odds = [r["n_tested"] for r in rows]
    print(
        f"  Attraktor-Atlas: N_B ∈ [{min(attr_sizes)}, {max(attr_sizes)}]  "
        f"→ Fall 1 (fixes/kurzes Cover), miss_attr alle "
        f"{[r['miss_attr'] for r in rows]}"
    )
    print(
        f"  Prefix-Atlas: N_B/odds → "
        f"{[round(r['N_B_prefix_over_odds'], 3) for r in rows]}  "
        f"(∼1 ⇒ Fall 3 für volle Präfixe — erwartet, keine Kompression)"
    )
    print(
        f"  max k*/L: {[round(r['max_k_star_over_L'], 3) for r in rows]}  "
        f"(soll ≤ C universell)"
    )

    print("\n[Hypothese]")
    print(
        "  Präzise Form: N_B^hit(L)=O(L^c) oder O(1) für Hit-Blöcke,\n"
        "  und max k*(n) ≤ C·L für n<2^L.\n"
        "  Volle Präfixe komprimieren nicht (jedes n eigenes Wort).\n"
        "  Operator-Kompression steckt in kurzen Hit-Blöcken + log-Zeit."
    )

    out = Path(__file__).resolve().parent / "collatz_cover_growth.json"
    log = Path(__file__).resolve().parent / "result_cover_growth.log"
    payload = {
        "claim": (
            "Scale-dependent cover growth N_B(L); "
            "prefix atlas ~2^L, attractor hit atlas O(1); "
            "not a Collatz proof"
        ),
        "levels": rows,
        "hypothesis": {
            "fall1": "N_B_hit = O(1)",
            "fall2": "N_B_hit = O(L^c)",
            "fall3": "N_B_prefix ~ 2^L (no compression)",
            "open": "exists c,C: N_B_hit(L)<=C L^c and max k* <= C L",
        },
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    lines = [
        "CollatzCoverGrowth",
        *[
            f"L={r['L']} N_pref={r['N_B_prefix']} N_attr={r['N_B_attr']} "
            f"N_contr={r['N_B_contract']} max_k*/L={r['max_k_star_over_L']} "
            f"miss_attr={r['miss_attr']}"
            for r in rows
        ],
    ]
    log.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    print(f"wrote {log}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
