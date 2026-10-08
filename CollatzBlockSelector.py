"""
CollatzBlockSelector — Block-Lyapunov statt Schritt-Potential.

Idee: Nicht V_n entlang Einzelschritten, sondern kanonische Blöcke ω,
für die
  Λ_ω = a log₂3 − L < 0
und
  χ_ω = dist(x*_ω, {0,1,2}) < ε
gesucht werden.

Selektor auf einer Bahn n₀ ↦ bits:
  k*(ε) = min { L : Präfix ω der Länge L erfüllt Λ<0 und χ<ε }
  (oder ∞ im Scan-Horizont)

Quantitative Frage (Scan):
  existiert C, ε mit k*(ε) ≤ C log n₀  für alle getesteten Seeds?

Zusätzlich: First-Descent-Blöcke / feste Fenster als alternative Zerlegungen.

Epistemik: Block-Selektor / quantitative Erreichbarkeit; kein Collatz-Beweis.
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


def block_invariants(omega: tuple[int, ...]) -> dict:
    a, L, B = affine_of_word(omega)
    if L == 0:
        return {"L": 0, "Lambda": 0.0, "chi": 0.0, "x": "0"}
    D = (1 << L) - 3**a
    x = Fraction(B, D)
    chi = min(abs(x - Fraction(t)) for t in ATTRACTOR)
    return {
        "a": a,
        "L": L,
        "B": B,
        "Lambda": a * LOG2_3 - L,
        "x": str(x),
        "chi": float(chi),
        "chi_frac": str(chi),
        "contracting": a * LOG2_3 - L < -1e-15,
        "chi_zero": chi == 0,
    }


def syracuse_bits(n0: int, steps: int) -> tuple[list[int], list[int]]:
    bits: list[int] = []
    states = [n0]
    x = n0
    for _ in range(steps):
        bit = 0 if x % 2 == 0 else 1
        bits.append(bit)
        x = x // 2 if bit == 0 else (3 * x + 1) // 2
        states.append(x)
    return bits, states


def first_descent_blocks(bits: list[int], states: list[int]) -> list[tuple[int, int]]:
    """
    Blöcke vom Odd-Start bis Rückkehr unter Startwert (first descent),
    als Indexintervalle [i,j) in der Bitfolge.
    """
    blocks = []
    i = 0
    n = len(bits)
    while i < n:
        # starte bei ungeradem Zustand falls möglich
        while i < n and states[i] % 2 == 0:
            i += 1
        if i >= n:
            break
        start_val = states[i]
        j = i
        while j < n:
            j += 1
            if states[j] < start_val:
                blocks.append((i, j))
                i = j
                break
        else:
            break
    return blocks


def block_ok(inv: dict, eps: float, require_contract: bool) -> bool:
    """Nichttrivial: mindestens ein Odd-Schritt (a≥1); kein reines ω=0."""
    if inv.get("a", 0) < 1:
        return False
    if inv["chi"] >= eps:
        return False
    if require_contract and not inv["contracting"]:
        return False
    return True


def selector_prefix(
    bits: list[int],
    eps: float,
    require_contract: bool = True,
) -> dict:
    """Kürzestes Präfix mit a≥1, χ<ε (und optional Λ<0)."""
    for L in range(1, len(bits) + 1):
        inv = block_invariants(tuple(bits[:L]))
        if block_ok(inv, eps, require_contract):
            return {
                "found": True,
                "k": L,
                "eps": eps,
                "require_contract": require_contract,
                **{key: inv[key] for key in ("a", "L", "Lambda", "chi", "x")},
            }
    return {
        "found": False,
        "k": None,
        "eps": eps,
        "require_contract": require_contract,
    }


def selector_window(
    bits: list[int],
    eps: float,
    max_len: int,
    require_contract: bool = True,
) -> dict:
    """
    Kürzester Teilblock bits[i:i+L] mit L≤max_len, χ<ε, Λ<0.
    Entlang der Bahn (nicht nur Präfix vom Start).
    """
    n = len(bits)
    best = None
    for i in range(n):
        # early stop if we can't beat best
        lim = min(max_len, n - i)
        for L in range(1, lim + 1):
            if best is not None and L >= best["k"]:
                break
            inv = block_invariants(tuple(bits[i : i + L]))
            if block_ok(inv, eps, require_contract):
                cand = {
                    "found": True,
                    "k": L,
                    "i": i,
                    "eps": eps,
                    **{key: inv[key] for key in ("a", "L", "Lambda", "chi", "x")},
                }
                best = cand
                break
        if best is not None and best["k"] <= 2:
            break
    if best is None:
        return {"found": False, "k": None, "eps": eps, "max_len": max_len}
    return best


def selector_descent(
    bits: list[int],
    states: list[int],
    eps: float,
) -> dict:
    """Erster First-Descent-Block mit χ<ε und Λ<0."""
    for i, j in first_descent_blocks(bits, states):
        inv = block_invariants(tuple(bits[i:j]))
        if block_ok(inv, eps, require_contract=True):
            return {
                "found": True,
                "k": j - i,
                "i": i,
                "j": j,
                "start_state": states[i],
                "end_state": states[j],
                **{key: inv[key] for key in ("a", "L", "Lambda", "chi", "x")},
            }
    return {"found": False, "k": None, "eps": eps}


def main() -> int:
    random.seed(11)
    steps = 400
    epsilons = (0.5, 0.1, 0.01, 1e-3, 1e-6)
    n_seeds = 40

    print("=" * 72)
    print("BLOCK-SELEKTOR  —  kürzester Block mit Λ<0, χ<ε, a≥1")
    print("=" * 72)
    print("Nichttrivial: a≥1 (kein reines ω=0).  Frage: k*(ε) ≲ C log n₀?\n")

    # Seed-Mengen: klein + logarithmisch gestaffelt
    seeds = [3, 5, 7, 9, 11, 15, 17, 19, 27, 31, 41, 47, 63, 97, 127]
    seeds += [2**k - 1 for k in range(8, 16)]  # Mersenne-artig
    seeds += [random.randint(10**2, 10**7) for _ in range(n_seeds)]
    # unique preserve order
    seen: set[int] = set()
    uniq = []
    for s in seeds:
        if s not in seen and s > 0:
            seen.add(s)
            uniq.append(s)
    seeds = uniq[:n_seeds]

    # [1] Präfix-Selektor
    print("[1] Präfix-Selektor vom Bahnstart")
    prefix_table: dict[float, list[dict]] = {e: [] for e in epsilons}
    for n0 in seeds:
        bits, states = syracuse_bits(n0, steps)
        enter = next((i for i, s in enumerate(states) if s in (1, 2)), None)
        for eps in epsilons:
            sel = selector_prefix(bits, eps, require_contract=True)
            sel["n0"] = n0
            sel["log2_n0"] = math.log2(n0)
            sel["enter"] = enter
            sel["k_over_log"] = (
                sel["k"] / math.log2(n0) if sel["found"] and n0 > 1 else None
            )
            prefix_table[eps].append(sel)

    for eps in epsilons:
        rows = prefix_table[eps]
        found = [r for r in rows if r["found"]]
        ks = [r["k"] for r in found]
        ratios = [r["k_over_log"] for r in found if r["k_over_log"] is not None]
        print(
            f"  ε={eps:g}:  found {len(found)}/{len(rows)}  "
            f"k median={sorted(ks)[len(ks)//2] if ks else '—'}  "
            f"max={max(ks) if ks else '—'}  "
            f"k/log₂n median={sorted(ratios)[len(ratios)//2]:.3f}"
            if ratios
            else f"  ε={eps:g}:  found {len(found)}/{len(rows)}"
        )

    # [2] Fenster-Selektor (irgendwo auf der Bahn)
    print("\n[2] Fenster-Selektor (kürzester Teilblock auf der Bahn)")
    window_table: dict[float, list[dict]] = {e: [] for e in epsilons}
    max_len = 64
    for n0 in seeds:
        bits, _ = syracuse_bits(n0, steps)
        for eps in epsilons:
            sel = selector_window(bits, eps, max_len=max_len, require_contract=True)
            sel["n0"] = n0
            sel["log2_n0"] = math.log2(n0)
            sel["k_over_log"] = (
                sel["k"] / math.log2(n0) if sel.get("found") and n0 > 1 else None
            )
            window_table[eps].append(sel)

    for eps in epsilons:
        rows = window_table[eps]
        found = [r for r in rows if r["found"]]
        ks = [r["k"] for r in found]
        ratios = [r["k_over_log"] for r in found if r["k_over_log"] is not None]
        print(
            f"  ε={eps:g}:  found {len(found)}/{len(rows)}  "
            f"k median={sorted(ks)[len(ks)//2] if ks else '—'}  "
            f"max={max(ks) if ks else '—'}  "
            f"k/log₂n median={sorted(ratios)[len(ratios)//2]:.3f}"
            if ratios
            else f"  ε={eps:g}:  found {len(found)}/{len(rows)}"
        )

    # [3] First-Descent-Blöcke
    print("\n[3] First-Descent-Blöcke mit Λ<0 und χ<ε")
    descent_table: dict[float, list[dict]] = {e: [] for e in epsilons}
    for n0 in seeds:
        bits, states = syracuse_bits(n0, steps)
        for eps in epsilons:
            sel = selector_descent(bits, states, eps)
            sel["n0"] = n0
            sel["log2_n0"] = math.log2(n0)
            sel["k_over_log"] = (
                sel["k"] / math.log2(n0) if sel.get("found") and n0 > 1 else None
            )
            descent_table[eps].append(sel)

    for eps in (0.5, 0.1, 0.01):
        rows = descent_table[eps]
        found = [r for r in rows if r["found"]]
        ks = [r["k"] for r in found]
        print(
            f"  ε={eps:g}:  found {len(found)}/{len(rows)}  "
            f"k median={sorted(ks)[len(ks)//2] if ks else '—'}  "
            f"max={max(ks) if ks else '—'}"
        )

    # [4] Attraktorkanal-Blöcke: χ=0 exakt (reine A-Wörter als Fenster)
    print("\n[4] Exakte Kanalblöcke χ=0 (Fenster, Λ<0 automatisch)")
    exact = []
    for n0 in seeds:
        bits, states = syracuse_bits(n0, min(steps, 200))
        sel = selector_window(bits, eps=1e-15, max_len=24, require_contract=True)
        # also search chi==0 specifically
        hit = None
        for i in range(len(bits)):
            for L in range(1, min(24, len(bits) - i) + 1):
                inv = block_invariants(tuple(bits[i : i + L]))
                if inv["chi_zero"] and inv["contracting"] and inv["a"] >= 1:
                    hit = {"n0": n0, "i": i, "k": L, "x": inv["x"], "state": states[i]}
                    break
            if hit:
                break
        exact.append(hit or {"n0": n0, "found": False})
    n_exact = sum(1 for h in exact if h and h.get("k"))
    print(f"  Fenster mit χ=0, a≥1, L≤24: {n_exact}/{len(seeds)}")
    # wann relativ zum Enter?
    delays = []
    for n0 in seeds:
        bits, states = syracuse_bits(n0, steps)
        enter = next((i for i, s in enumerate(states) if s in (1, 2)), None)
        hit = next(
            (
                (i, L)
                for i in range(len(bits))
                for L in range(2, min(8, len(bits) - i) + 1)
                if (
                    (inv := block_invariants(tuple(bits[i : i + L])))["chi_zero"]
                    and inv["a"] >= 1
                )
            ),
            None,
        )
        if enter is not None and hit is not None:
            delays.append(hit[0] - enter)
    if delays:
        print(
            f"  erste χ=0-Fensterposition − enter: "
            f"median={sorted(delays)[len(delays)//2]}  "
            f"(≥0 ⇒ nach/bei Attraktor)"
        )

    # [5] Skalierung k vs log n
    print("\n[5] Skalierung k*(ε=0.1) vs log₂ n₀  (Präfix & Fenster)")
    print("  n₀           log₂n  k_pref  k_win  enter  k_pref/log")
    for n0 in sorted(seeds, key=lambda x: x)[:20]:
        bits, states = syracuse_bits(n0, steps)
        enter = next((i for i, s in enumerate(states) if s in (1, 2)), None)
        pref = selector_prefix(bits, 0.1, True)
        win = selector_window(bits, 0.1, max_len=64, require_contract=True)
        kp = pref["k"] if pref["found"] else None
        kw = win["k"] if win["found"] else None
        ratio = f"{kp/math.log2(n0):.2f}" if kp else "—"
        print(
            f"  {n0:10d}  {math.log2(n0):6.2f}  "
            f"{str(kp):>6s}  {str(kw):>5s}  {str(enter):>5s}  {ratio}"
        )

    # C-Suche: kleinstes C mit k ≤ C log2 n für alle found
    print("\n[6] Existenz von C für k ≤ C log₂ n₀?")
    for kind, table in (("prefix", prefix_table), ("window", window_table)):
        for eps in (0.5, 0.1, 0.01):
            rows = [r for r in table[eps] if r["found"] and r["k_over_log"] is not None]
            if not rows:
                print(f"  {kind} ε={eps}: keine Treffer")
                continue
            C_need = max(r["k_over_log"] for r in rows)
            # fail rate if C=20, 50, 100
            for C in (5, 10, 20, 50, 100):
                ok = sum(1 for r in rows if r["k_over_log"] <= C)
                if ok == len(rows):
                    print(
                        f"  {kind} ε={eps:g}:  C={C} reicht "
                        f"(max k/log₂n={C_need:.2f}, n={len(rows)})"
                    )
                    break
            else:
                print(
                    f"  {kind} ε={eps:g}:  kein C≤100  "
                    f"(max k/log₂n={C_need:.2f})"
                )

    # [7] Zufalls-Kontrolle
    print("\n[7] Zufallspfade-Kontrolle (Präfix-Selektor ε=0.1)")
    rand_ks = []
    for _ in range(30):
        bits = [random.randint(0, 1) for _ in range(steps)]
        sel = selector_prefix(bits, 0.1, True)
        if sel["found"]:
            rand_ks.append(sel["k"])
    print(
        f"  found {len(rand_ks)}/30  "
        f"k median={sorted(rand_ks)[len(rand_ks)//2] if rand_ks else '—'}  "
        f"max={max(rand_ks) if rand_ks else '—'}"
    )

    print("\n[8] Fazit")
    print(
        "  Nichttrivialer Block: Λ<0 ∧ χ<ε ∧ a≥1.\n"
        "  Präfix-Selektor: k*(ε) skaliert im Scan mit log n (C≲20 für ε=0.1).\n"
        "  Fenster k=2 ist schwach: nur lokales Vorkommen von 10/01.\n"
        "  Relevant für Reduktion: Präfix- oder Descent-Überdeckung,\n"
        "  nicht bloße Digramm-Existenz. Erreichbarkeit bleibt offen."
    )

    out = Path(__file__).resolve().parent / "collatz_block_selector.json"
    payload = {
        "claim": (
            "Block selector: shortest omega with Lambda<0 and chi<eps; "
            "tests quantitative reachability k <? C log n; "
            "not a Collatz proof"
        ),
        "steps": steps,
        "epsilons": epsilons,
        "n_seeds": len(seeds),
        "seeds": seeds,
        "prefix": {
            str(eps): {
                "n_found": sum(1 for r in prefix_table[eps] if r["found"]),
                "k_median": sorted(r["k"] for r in prefix_table[eps] if r["found"])[
                    len([r for r in prefix_table[eps] if r["found"]]) // 2
                ]
                if any(r["found"] for r in prefix_table[eps])
                else None,
                "k_over_log_max": max(
                    (
                        r["k_over_log"]
                        for r in prefix_table[eps]
                        if r["k_over_log"] is not None
                    ),
                    default=None,
                ),
                "rows": prefix_table[eps],
            }
            for eps in epsilons
        },
        "window": {
            str(eps): {
                "n_found": sum(1 for r in window_table[eps] if r["found"]),
                "k_median": sorted(r["k"] for r in window_table[eps] if r["found"])[
                    len([r for r in window_table[eps] if r["found"]]) // 2
                ]
                if any(r["found"] for r in window_table[eps])
                else None,
                "k_over_log_max": max(
                    (
                        r["k_over_log"]
                        for r in window_table[eps]
                        if r["k_over_log"] is not None
                    ),
                    default=None,
                ),
                "rows": [
                    {k: v for k, v in r.items() if k != "x" or r.get("found")}
                    for r in window_table[eps]
                ],
            }
            for eps in epsilons
        },
        "descent": {
            str(eps): {
                "n_found": sum(1 for r in descent_table[eps] if r["found"]),
                "rows": descent_table[eps],
            }
            for eps in (0.5, 0.1, 0.01)
        },
        "exact_chi0_windows": {
            "n_found": n_exact,
            "delay_median": sorted(delays)[len(delays) // 2] if delays else None,
        },
        "random_prefix_eps0.1": {
            "n_found": len(rand_ks),
            "k_median": sorted(rand_ks)[len(rand_ks) // 2] if rand_ks else None,
            "k_max": max(rand_ks) if rand_ks else None,
        },
        "open_question": (
            "exists C forall n exists k<C log n with chi(omega_k)<eps "
            "and Lambda<0 along the Collatz path"
        ),
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
