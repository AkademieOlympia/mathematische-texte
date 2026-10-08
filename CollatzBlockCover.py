"""
CollatzBlockCover — naechster Test: endliche Blockmenge B.

Ziel: endliche Menge kontraktiver Bloecke, die Bahnen mit R<=2 ueberdeckt.
"""

from __future__ import annotations

import json
import math
import random
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path


def affine(omega):
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


def path_bits(n0, steps):
    bits, states = [], [n0]
    x = n0
    for _ in range(steps):
        b = 0 if x % 2 == 0 else 1
        bits.append(b)
        x = x // 2 if b == 0 else (3 * x + 1) // 2
        states.append(x)
    return bits, states


def build_library(Lmax=8):
    lib = []
    for L in range(2, Lmax + 1):
        for bits in product([0, 1], repeat=L):
            if sum(bits) == 0:
                continue
            a, _, B, A, x = affine(bits)
            if A >= 1:
                continue
            if abs(x) > 20:
                continue
            lib.append({
                "omega": "".join(map(str, bits)),
                "bits": bits,
                "L": L,
                "A": A,
                "x": x,
                "chi": min(abs(x - t) for t in (0, 1, 2)),
            })
    lib.sort(key=lambda r: (r["chi"], r["L"], float(r["A"])))
    return lib


def first_hit(bits, states, library, R_thresh=2.0):
    by_L = {}
    for b in library:
        by_L.setdefault(b["L"], []).append(b)
    n = len(bits)
    for i in range(n):
        for L, blocks in by_L.items():
            if i + L > n:
                continue
            seg = tuple(bits[i:i+L])
            for blk in blocks:
                if blk["bits"] != seg:
                    continue
                s = states[i]
                _, _, _, A, x = affine(seg)
                R = abs(float(A * (s - x)))
                if R <= R_thresh:
                    return {"i": i, "L": L, "omega": blk["omega"], "s": s, "R": R, "x": float(x)}
    return None


def main():
    random.seed(29)
    lib = build_library(8)
    core = [b for b in lib if b["chi"] == 0][:30]
    near = [b for b in lib if 0 < b["chi"] < 0.5][:50]
    library = core + near
    print("=" * 72)
    print("BLOCK-COVER  —  endliche Bibliothek")
    print("=" * 72)
    print(f"Bibliothek: {len(library)} (chi=0: {len(core)}, nahe: {len(near)})\n")

    seeds = list(range(3, 80))
    seeds += [2**k - 1 for k in range(8, 24)]
    seeds += [random.randint(10**4, 10**8) for _ in range(40)]
    seen = set()
    seeds = [s for s in seeds if s not in seen and not seen.add(s)]

    hits, miss = [], []
    for n0 in seeds:
        bits, states = path_bits(n0, 400)
        enter = next((i for i, s in enumerate(states) if s in (1, 2)), None)
        h = first_hit(bits, states, library, 2.0)
        if h is None:
            miss.append(n0)
            continue
        h["n0"] = n0
        h["enter"] = enter
        h["Q"] = h["i"] / math.log2(n0)
        hits.append(h)

    print(f"[1] Cover-Treffer: {len(hits)}/{len(seeds)}  miss={len(miss)}")
    if hits:
        Is = [h["i"] for h in hits]
        Qs = [h["Q"] for h in hits]
        print(f"  i_hit med={sorted(Is)[len(Is)//2]} max={max(Is)} Q_med={sorted(Qs)[len(Qs)//2]:.3f} Q_max={max(Qs):.3f}")
        eq = sum(1 for h in hits if h["enter"] is not None and h["i"] <= h["enter"])
        print(f"  i_hit <= enter: {eq}/{len(hits)}")
        print("  Top-Bloecke:")
        for w, c in Counter(h["omega"] for h in hits).most_common(8):
            print(f"    {w}: {c}")

    hits2 = 0
    for n0 in seeds:
        bits, states = path_bits(n0, 400)
        if first_hit(bits, states, core, 2.0):
            hits2 += 1
    print(f"\n[2] Nur chi=0: {hits2}/{len(seeds)}")
    print("\n[3] Naechster Beweisschritt: endliche B -> AffineLogCovering")

    Path(__file__).resolve().parent.joinpath("collatz_block_cover.json").write_text(
        json.dumps({
            "library_size": len(library),
            "n_hit": len(hits),
            "n_miss": len(miss),
            "core_only_hit": hits2,
            "top_blocks": Counter(h["omega"] for h in hits).most_common(20),
        }, indent=2) + "\n"
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
