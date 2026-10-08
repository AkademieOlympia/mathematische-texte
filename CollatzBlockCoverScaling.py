"""
CollatzBlockCoverScaling — Q_B(n) = k_hit(n) / log₂ n.

Operatoratlas (nicht Wortatlas):
  Treffer erst, wenn ein Block ω ∈ 𝓑 am Zustand s die Parität matcht
  UND R_ω(s) = |A|·|s − x*| ≤ 2.

Skalierungs-Audit über N ∈ {10³,10⁴,10⁵,…}:
  Szenario A: max Q_B bleibt beschränkt, Misses = 0 bei fester 𝓑
  Szenario B: max Q_B wächst oder 𝓑 muss wachsen

Epistemik: kombinatorisches Cover; kein Collatz-Beweis.
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


def build_B(mode: str = "chi0", Lmax: int = 8, x_bound: float = 20.0) -> list[dict]:
    """Endliche Blockmenge 𝓑."""
    lib = []
    for L in range(2, Lmax + 1):
        for bits in product([0, 1], repeat=L):
            if sum(bits) == 0:
                continue
            a, _, _, A, x = affine(bits)
            if float(A) >= 1:
                continue
            if abs(float(x)) > x_bound:
                continue
            chi = min(abs(float(x) - t) for t in (0.0, 1.0, 2.0))
            if mode == "chi0" and chi != 0:
                continue
            if mode == "near" and not (chi == 0 or chi < 0.5):
                continue
            lib.append(
                {
                    "omega": "".join(map(str, bits)),
                    "bits": bits,
                    "L": L,
                    "A": A,
                    "x": x,
                    "chi": chi,
                }
            )
    lib.sort(key=lambda r: (r["chi"], r["L"], float(r["A"])))
    return lib


def path_until_enter(n0: int, max_steps: int = 100_000, post_enter: int = 16):
    """
    Bahn bis Eintritt in {1,2}, dann noch `post_enter` Schritte im Attraktor.
    Ohne Post-Zyklus sind Blöcke 10/01 am Eintrittszustand nie anwendbar.
    """
    bits: list[int] = []
    states = [n0]
    x = n0
    entered_at = None
    for _ in range(max_steps):
        b = 0 if x % 2 == 0 else 1
        bits.append(b)
        x = x // 2 if b == 0 else (3 * x + 1) // 2
        states.append(x)
        if entered_at is None and x in (1, 2):
            entered_at = len(bits)  # Index des Zustands x in states
            # noch post_enter weitere Schritte aufnehmen
            for _ in range(post_enter):
                b2 = 0 if x % 2 == 0 else 1
                bits.append(b2)
                x = x // 2 if b2 == 0 else (3 * x + 1) // 2
                states.append(x)
            break
    return bits, states, entered_at


def first_operator_hit(
    bits: list[int],
    states: list[int],
    library: list[dict],
    R_thresh: float = 2.0,
) -> dict | None:
    by_L: dict[int, list[dict]] = {}
    for blk in library:
        by_L.setdefault(blk["L"], []).append(blk)
    # Längen absteigend: längere Blöcke zuerst optional; hier kürzeste zuerst
    for i in range(len(bits)):
        for L in sorted(by_L):
            if i + L > len(bits):
                continue
            seg = tuple(bits[i : i + L])
            for blk in by_L[L]:
                if blk["bits"] != seg:
                    continue
                s = states[i]
                _, _, _, A, x = affine(seg)
                R = abs(float(A * (Fraction(s) - x)))
                if R <= R_thresh:
                    return {
                        "k_hit": i,  # Schritte vom Start bis Blockbeginn
                        "L": L,
                        "omega": blk["omega"],
                        "s": s,
                        "R": R,
                        "x": float(x),
                        "A": float(A),
                        "k_end": i + L,  # inkl. Blocklänge
                    }
    return None


def scan_scale(N: int, library: list[dict], odd_only: bool = True) -> dict:
    t0 = time.time()
    records = []
    misses = []
    n_tested = 0
    start = 3
    step = 2 if odd_only else 1
    for n in range(start, N + 1, step):
        n_tested += 1
        bits, states, enter = path_until_enter(n)
        hit = first_operator_hit(bits, states, library)
        if hit is None:
            misses.append(n)
            continue
        log2n = math.log2(n)
        # Q_B über Trefferbeginn (First passage in Cover)
        q = hit["k_hit"] / log2n if log2n > 0 else 0.0
        q_end = hit["k_end"] / log2n
        records.append(
            {
                "n": n,
                "k_hit": hit["k_hit"],
                "k_end": hit["k_end"],
                "Q_B": q,
                "Q_B_end": q_end,
                "omega": hit["omega"],
                "R": hit["R"],
                "enter": enter,
                "before_enter": enter is None or hit["k_hit"] <= enter,
                "at_or_after_enter": enter is not None and hit["k_hit"] >= enter,
            }
        )
    elapsed = time.time() - t0
    if not records:
        return {
            "N": N,
            "n_tested": n_tested,
            "n_hit": 0,
            "n_miss": len(misses),
            "misses_sample": misses[:20],
            "elapsed_s": elapsed,
        }
    Qs = [r["Q_B"] for r in records]
    Qe = [r["Q_B_end"] for r in records]
    return {
        "N": N,
        "n_tested": n_tested,
        "n_hit": len(records),
        "n_miss": len(misses),
        "misses_sample": misses[:20],
        "max_Q_B": max(Qs),
        "avg_Q_B": sum(Qs) / len(Qs),
        "med_Q_B": sorted(Qs)[len(Qs) // 2],
        "max_Q_B_end": max(Qe),
        "avg_Q_B_end": sum(Qe) / len(Qe),
        "max_k_hit": max(r["k_hit"] for r in records),
        "frac_before_enter": sum(1 for r in records if r["before_enter"]) / len(records),
        "frac_at_or_after_enter": sum(1 for r in records if r["at_or_after_enter"])
        / len(records),
        "top_blocks": Counter(r["omega"] for r in records).most_common(10),
        "worst_Q": sorted(records, key=lambda r: -r["Q_B"])[:8],
        "elapsed_s": elapsed,
        "library_size": len(library),
    }


def main() -> int:
    scales = [10**3, 10**4, 10**5]
    # optional: 10**6 is ~50x slower than 10**5 for odd-only; keep optional
    modes = {
        "B_chi0": build_B("chi0", Lmax=8),
        "B_near": build_B("near", Lmax=8),
    }

    print("=" * 72)
    print("BLOCK-COVER SKALIERUNG  —  Q_B(n)=k_hit/log₂n  (Operatoratlas)")
    print("=" * 72)
    print("Treffer ⇔ Wortmatch ∧ R≤2 am Zustand (nicht nur Suffix-Liste)\n")

    all_results = {}
    log_lines = []
    log_lines.append("CollatzBlockCoverScaling — Operatoratlas Q_B audit")
    log_lines.append("Treffer: parity match + R=|A|*|s-x*| <= 2")
    log_lines.append("")

    for mode, library in modes.items():
        print(f"--- {mode}: |𝓑|={len(library)} ---")
        log_lines.append(f"=== {mode} |B|={len(library)} ===")
        series = []
        for N in scales:
            res = scan_scale(N, library, odd_only=True)
            series.append(res)
            if res["n_hit"] == 0:
                msg = f"N={N:7d}  HITLOS  miss={res['n_miss']}"
            else:
                msg = (
                    f"N={N:7d} | max Q_B={res['max_Q_B']:.4f} | "
                    f"avg={res['avg_Q_B']:.4f} | med={res['med_Q_B']:.4f} | "
                    f"miss={res['n_miss']} | max_k={res['max_k_hit']} | "
                    f"t={res['elapsed_s']:.1f}s"
                )
            print(f"  {msg}")
            log_lines.append(msg)
            if res.get("worst_Q"):
                w = res["worst_Q"][0]
                detail = (
                    f"    worst: n={w['n']} Q_B={w['Q_B']:.4f} "
                    f"k_hit={w['k_hit']} ω={w['omega']} R={w['R']:.4g}"
                )
                print(detail)
                log_lines.append(detail)

        # Stabilität von max Q_B?
        maxes = [r["max_Q_B"] for r in series if r.get("max_Q_B") is not None]
        misses = [r["n_miss"] for r in series]
        if maxes:
            growth = maxes[-1] / maxes[0] if maxes[0] > 0 else float("inf")
            scenario = (
                "A_stable"
                if all(m == 0 for m in misses) and growth <= 1.5
                else (
                    "A_soft"
                    if all(m == 0 for m in misses)
                    else "B_miss_or_expand"
                )
            )
        else:
            scenario = "B_miss_or_expand"
            growth = None
        print(f"  Szenario: {scenario}  maxQ-growth={growth}")
        log_lines.append(f"Szenario: {scenario}  maxQ_growth={growth}")
        log_lines.append("")
        all_results[mode] = {
            "library_size": len(library),
            "library_omegas": [b["omega"] for b in library[:40]],
            "scales": series,
            "scenario": scenario,
            "maxQ_growth": growth,
        }

    # Lean-Schranken-Vorschlag
    chi0 = all_results["B_chi0"]["scales"]
    if chi0 and chi0[-1].get("max_Q_B") is not None:
        C_cand = math.ceil(max(r["max_Q_B"] for r in chi0) + 1e-9)
        lean_note = (
            f"Lean-Kandidat: C = {C_cand} fuer AffineLogCovering / "
            f"blockCover_implies_affineLogCovering (aus max Q_B auf N<={scales[-1]})"
        )
    else:
        C_cand = None
        lean_note = "Kein C-Kandidat (Misses oder leere Hits)"
    print(f"\n{lean_note}")
    log_lines.append(lean_note)

    print("\nFazit:")
    print(
        "  Wortatlas allein ist zu schwach; Operatoratlas verlangt R≤2.\n"
        "  Stabiles max Q_B + 0 Misses ⇒ endliches Cover numerisch plausibel.\n"
        "  Wachsendes max Q_B / Misses ⇒ Cover muss wachsen (Szenario B)."
    )

    out_dir = Path(__file__).resolve().parent
    json_path = out_dir / "collatz_block_cover_scaling.json"
    log_path = out_dir / "result.log"
    payload = {
        "claim": (
            "Operator-atlas block cover scaling Q_B=k_hit/log2 n; "
            "scenario A vs B; not a Collatz proof"
        ),
        "scales": scales,
        "results": all_results,
        "lean_C_candidate": C_cand,
        "lean_target": "blockCover_implies_affineLogCovering",
    }
    json_path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print(f"\nwrote {json_path}")
    print(f"wrote {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
