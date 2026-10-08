"""
CollatzL — Top-Down-Messung der Ballot-Tiefe L(n).

L(n) = max { K : die ersten K Odd-to-Odd-Blöcke von n erfüllen
                 2^{Q_j} ≤ 3^{R_j} für alle j ≤ K }.

Ein hypothetischer unendlicher No-Descent-Integerpfad im Regime A≤1
ist genau L(n)=∞. Empirische Schranken an L(n) sind Diagnose, kein Beweis.

Epistemik: Integer-Census. Kein Collatz-Beweis.
Oszillation A_K ≶ 1 mit n_K ≥ n_0 bleibt separat.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass

from CollatzBallotTree import actual_block, prefix_safe


@dataclass
class LRecord:
    n: int
    L: int
    reason: str
    path: list[tuple[int, int]]
    n_exit: int
    n_next: int
    break_block: tuple[int, int] | None
    descended: bool
    hit_one: bool
    R: int
    S: int
    Q: int


def measure_L(n_start: int, max_blocks: int = 200) -> LRecord:
    """Maximale Präfixtiefe im gefährlichen Regime. n_start ungerade > 1."""
    if n_start % 2 == 0 or n_start <= 1:
        raise ValueError("n_start must be odd and > 1")
    curr = n_start
    R = S = 0
    path: list[tuple[int, int]] = []
    for k in range(1, max_blocks + 1):
        if curr == 1:
            return LRecord(
                n=n_start,
                L=k - 1,
                reason="Reached 1",
                path=path,
                n_exit=curr,
                n_next=curr,
                break_block=None,
                descended=curr < n_start,
                hit_one=True,
                R=R,
                S=S,
                Q=R + S,
            )
        r, s, nxt = actual_block(curr)
        if not prefix_safe(R + r, S + s):
            return LRecord(
                n=n_start,
                L=k - 1,
                reason="Exited dangerous regime",
                path=path,
                n_exit=curr,
                n_next=nxt,
                break_block=(r, s),
                descended=nxt < n_start,
                hit_one=False,
                R=R,
                S=S,
                Q=R + S,
            )
        path.append((r, s))
        R += r
        S += s
        curr = nxt
    return LRecord(
        n=n_start,
        L=max_blocks,
        reason="Max depth reached",
        path=path,
        n_exit=curr,
        n_next=curr,
        break_block=None,
        descended=curr < n_start,
        hit_one=curr == 1,
        R=R,
        S=S,
        Q=R + S,
    )


def census_L(X: int, max_blocks: int = 200) -> dict:
    """Alle ungeraden n in (1, X]."""
    recs: dict[int, LRecord] = {}
    hist: Counter[int] = Counter()
    records: list[LRecord] = []
    best = -1
    max_ratio = 0.0
    arg_ratio = None
    n_capped = 0
    n_hit_one_inside = 0
    n_exit_descend = 0
    n_exit_rise = 0

    for n in range(3, X + 1, 2):
        rec = measure_L(n, max_blocks=max_blocks)
        recs[n] = rec
        hist[rec.L] += 1
        if rec.L > best:
            best = rec.L
            records.append(rec)
        if rec.reason == "Max depth reached":
            n_capped += 1
        if rec.hit_one:
            n_hit_one_inside += 1
        if rec.reason == "Exited dangerous regime":
            if rec.descended:
                n_exit_descend += 1
            else:
                n_exit_rise += 1
        ratio = rec.L / math.log2(n)
        if ratio > max_ratio:
            max_ratio = ratio
            arg_ratio = n

    # m(K)
    mK = []
    maxL = max(hist) if hist else 0
    for K in range(0, maxL + 1):
        cands = [n for n, rec in recs.items() if rec.L >= K]
        if not cands:
            break
        m = min(cands)
        mK.append((K, m, recs[m].L))

    top = sorted(recs.values(), key=lambda r: (-r.L, r.n))[:20]
    return {
        "X": X,
        "n_odd": (X - 1) // 2,
        "hist": dict(sorted(hist.items())),
        "max_L": best,
        "record_setters": records,
        "top": top,
        "mK": mK,
        "max_L_over_log2": max_ratio,
        "arg_max_ratio": arg_ratio,
        "n_capped": n_capped,
        "n_hit_one_inside": n_hit_one_inside,
        "n_exit_descend": n_exit_descend,
        "n_exit_rise": n_exit_rise,
        "mean_L": sum(rec.L for rec in recs.values()) / len(recs),
    }


def fmt_path(path: list[tuple[int, int]], limit: int = 12) -> str:
    if len(path) <= limit:
        return str(path)
    return str(path[:limit]) + f" … (+{len(path) - limit})"


def main() -> None:
    print("=== Spot-Checks ===", flush=True)
    for n in (7, 27, 31, 71, 127, 255, 511, 703, 1023, 7527, 10087, 35655):
        rec = measure_L(n)
        print(
            f"  n={n}: L={rec.L}  {rec.reason}  break={rec.break_block}  "
            f"descended={rec.descended}  path={fmt_path(rec.path)}",
            flush=True,
        )

    for X in (10**5, 10**6):
        print(f"\n=== Census L(n), ungerade n ≤ {X} ===", flush=True)
        out = census_L(X)
        print(f"  odd n = {out['n_odd']}, max L = {out['max_L']}, mean L = {out['mean_L']:.3f}", flush=True)
        print(
            f"  max L/log2(n) = {out['max_L_over_log2']:.4f} at n={out['arg_max_ratio']}",
            flush=True,
        )
        print(
            f"  exit+descend = {out['n_exit_descend']}, exit+rise = {out['n_exit_rise']}, "
            f"hit 1 inside = {out['n_hit_one_inside']}, capped = {out['n_capped']}",
            flush=True,
        )
        print("  Histogram L:", flush=True)
        for L, c in out["hist"].items():
            if c:
                print(f"    L={L:2d}: {c}", flush=True)
        print("  Record-Setter (neues Max-L bei wachsendem n):", flush=True)
        for rec in out["record_setters"]:
            print(
                f"    n={rec.n} L={rec.L}  {fmt_path(rec.path, 8)}",
                flush=True,
            )
        print("  Top-20 nach L:", flush=True)
        for rec in out["top"]:
            print(
                f"    n={rec.n} L={rec.L}  L/log2={rec.L/math.log2(rec.n):.3f}  "
                f"break={rec.break_block} desc={rec.descended}",
                flush=True,
            )
        print("  m(K):", flush=True)
        last_m = None
        for K, m, d in out["mK"]:
            if m != last_m:
                print(f"    m({K}) = {m} (L={d})", flush=True)
                last_m = m
            elif K == out["mK"][-1][0]:
                print(f"    … bis K={K}, m={m}", flush=True)


if __name__ == "__main__":
    main()
