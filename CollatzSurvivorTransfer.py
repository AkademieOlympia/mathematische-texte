"""
CollatzSurvivorTransfer — Survivor-Zylindergraph, Φ, Blockdrift, Spektralradius.

Zielstruktur (asymptotische Drift, nicht punktweise Λ ≤ −εL):

  ∀ ω ∈ Survivor:  liminf Λ_L < 0
  oder stärker:     liminf Λ_L/L ≤ −ε

Operator auf Survivor-Knoten ω (Wortlänge L):
  Φ(ω) = max{ Λ_{ω0}, Λ_{ω1} }   (über existierende Survivor-Kinder)
  ΔΛ(ω→ωb) = Λ_{ωb} − Λ_ω

Blockmittel:
  (1/m) Σ ΔΛ_i  ≤ −δ   auf Survivor-Pfaden / Kanten

Übergangsgraph G_L:
  Knoten = Survivor-Wörter Länge L (mode Clog, Schranke C)
  Kante ω → ω' wenn ω' = (ω·b)[1:]  (Shift) oder Prefix-Extension ω→ωb

Epistemik: Diagnostik der Übergangsmatrix; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def syracuse_step(n: int) -> int:
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def hit_time(n: int, max_steps: int) -> int | None:
    if n in (1, 2):
        return 0
    x = n
    for k in range(1, max_steps + 1):
        x = syracuse_step(x)
        if x in (1, 2):
            return k
    return None


def affine_from_word(word: str) -> tuple[int, float, float, float]:
    """a, Λ, A, und B für T(n)=(3^a n + B)/2^L; x* = B/(2^L - 3^a) falls A≠1."""
    L = len(word)
    a = word.count("1")
    # B aus der affinen Akkumulation: start T=id, apply bits left-to-right
    # Even: n ↦ n/2; Odd: n ↦ (3n+1)/2
    # In Integer-Lift: A_num/A_den, B_num/B_den with den=2^L
    A_num, B_num = 1, 0  # T = (A_num n + B_num) / 2^k accumulating
    k = 0
    for b in word:
        if b == "0":
            # n ↦ n/2: A' = A/2, B' = B/2
            k += 1
        else:
            # n ↦ (3n+1)/2: A' = 3A/2, B' = (3B+1)/2
            A_num *= 3
            B_num = 3 * B_num + (1 << k)
            k += 1
    assert k == L
    Lam = a * math.log2(3) - L
    A = (3**a) / (2**L)
    B = B_num / (2**L)
    return a, Lam, A, B


def x_star(A: float, B: float) -> float | None:
    if abs(1 - A) < 1e-15:
        return None
    return B / (1 - A)


def residue_for_word(word: str) -> int:
    """Kleinstes r ≥ 0 mit Paritätswort = word (mod 2^L)."""
    L = len(word)
    # Löse das System: Residue aus Bits vorwärts
    # Standard: r bestimmt Bits via Syracuse mod 2^L
    # Brute für L≤20 ok: scan
    mod = 1 << L
    target = word
    for r in range(mod):
        x = r
        bits = []
        for _ in range(L):
            if x % 2 == 0:
                bits.append("0")
                x = (x // 2) % mod
            else:
                bits.append("1")
                x = ((3 * x + 1) // 2) % mod
        if "".join(bits) == target:
            return r
    raise ValueError(f"no residue for {word}")


def is_survivor_residue(r: int, L: int, C: float) -> tuple[bool, int | None, float]:
    n = r if r > 1 else r + (1 << L)
    budget = int(max(C * math.log2(n) * 3 + 50, 500))
    ht = hit_time(n, budget)
    thresh = C * math.log2(n) + 2
    return (ht is None or ht > thresh), ht, thresh


def enumerate_survivors(L: int, C: float) -> list[dict]:
    mod = 1 << L
    out = []
    for r in range(mod):
        ok, ht, thresh = is_survivor_residue(r, L, C)
        if not ok:
            continue
        # word
        x, bits = r, []
        a = 0
        for _ in range(L):
            if x % 2 == 0:
                bits.append("0")
                x = (x // 2) % mod
            else:
                bits.append("1")
                a += 1
                x = ((3 * x + 1) // 2) % mod
        word = "".join(bits)
        Lam = a * math.log2(3) - L
        A = (3**a) / (2**L)
        _, _, _, B = affine_from_word(word)
        xs = x_star(A, B)
        out.append(
            {
                "r": r,
                "n": r if r > 1 else r + mod,
                "word": word,
                "a": a,
                "Lambda": Lam,
                "A": A,
                "B": B,
                "x_star": xs,
                "ht": ht,
                "thresh": thresh,
            }
        )
    return out


def extension_graph(survivors_L: list[dict], survivors_Lp1: list[dict]) -> dict:
    """Prefix-Extension: ω → ωb wenn beide Survivor."""
    by_word = {s["word"]: s for s in survivors_Lp1}
    nodes = {s["word"]: s for s in survivors_L}
    edges = []
    Phi = {}
    for w, s in nodes.items():
        kids = []
        dLams = []
        for b in "01":
            wb = w + b
            if wb in by_word:
                child = by_word[wb]
                dLam = child["Lambda"] - s["Lambda"]
                kids.append(wb)
                dLams.append(dLam)
                edges.append(
                    {
                        "from": w,
                        "to": wb,
                        "bit": b,
                        "Lambda_from": s["Lambda"],
                        "Lambda_to": child["Lambda"],
                        "dLambda": dLam,
                    }
                )
        Phi[w] = max(dLams) + s["Lambda"] if dLams else None  # max Λ of children
        Phi[w] = {
            "Phi": max((by_word[k]["Lambda"] for k in kids), default=None),
            "n_children": len(kids),
            "mean_dLambda": float(np.mean(dLams)) if dLams else None,
            "max_dLambda": max(dLams) if dLams else None,
            "min_dLambda": min(dLams) if dLams else None,
        }
    return {"edges": edges, "Phi": Phi, "n_nodes": len(nodes), "n_edges": len(edges)}


def shift_adjacency(survivors: list[dict]) -> tuple[np.ndarray, list[str]]:
    """
    De-Bruijn/Shift-Graph auf Survivor-Wörtern fester Länge L:
    Kante u→v wenn u[1:]==v[:-1] und uv-kompatibel (u + v[-1] hat Präfix u).
    """
    words = [s["word"] for s in survivors]
    idx = {w: i for i, w in enumerate(words)}
    n = len(words)
    M = np.zeros((n, n), dtype=float)
    for u in words:
        for b in "01":
            v = u[1:] + b
            if v in idx:
                M[idx[u], idx[v]] = 1.0
    return M, words


def spectral_radius(M: np.ndarray) -> float:
    if M.size == 0:
        return 0.0
    # nichtnegativ → Perron
    vals = np.linalg.eigvals(M)
    return float(np.max(np.abs(vals)))


def block_drift_stats(edges: list[dict]) -> dict:
    if not edges:
        return {"n": 0}
    d = np.array([e["dLambda"] for e in edges], dtype=float)
    return {
        "n_edges": len(d),
        "mean_dLambda": float(d.mean()),
        "median_dLambda": float(np.median(d)),
        "frac_neg": float(np.mean(d < 0)),
        "p10": float(np.percentile(d, 10)),
        "p90": float(np.percentile(d, 90)),
        "min": float(d.min()),
        "max": float(d.max()),
    }


def path_block_means(survivors_by_L: dict[int, list[dict]], L0: int, m: int) -> dict:
    """
    Sample extension paths of length m from survivors at L0;
    report mean ΔΛ along each path and fraction with mean ≤ −δ candidates.
    """
    # Build nested maps word -> children words
    maps = {}
    for L in range(L0, L0 + m):
        S = survivors_by_L.get(L, [])
        Sp = survivors_by_L.get(L + 1, [])
        by = {s["word"]: s for s in Sp}
        ch = {}
        for s in S:
            kids = []
            for b in "01":
                wb = s["word"] + b
                if wb in by:
                    kids.append(by[wb])
            ch[s["word"]] = kids
        maps[L] = ch

    roots = survivors_by_L.get(L0, [])
    means = []
    worst = []
    for s0 in roots:
        # DFS limited: all survivor extensions of length m
        stack = [(s0["word"], s0["Lambda"], 0, 0.0)]  # word, Lam, depth, sum_d
        while stack:
            w, Lam, depth, sum_d = stack.pop()
            if depth == m:
                means.append(sum_d / m)
                continue
            for child in maps.get(L0 + depth, {}).get(w, []):
                d = child["Lambda"] - Lam
                stack.append((child["word"], child["Lambda"], depth + 1, sum_d + d))
    arr = np.array(means) if means else np.array([])
    if arr.size == 0:
        return {"n_paths": 0}
    return {
        "n_paths": int(arr.size),
        "m": m,
        "mean_of_block_means": float(arr.mean()),
        "frac_block_mean_neg": float(np.mean(arr < 0)),
        "frac_block_mean_le_-0.1": float(np.mean(arr <= -0.1)),
        "frac_block_mean_le_-0.2": float(np.mean(arr <= -0.2)),
        "min_block_mean": float(arr.min()),
        "max_block_mean": float(arr.max()),
        "p10_block_mean": float(np.percentile(arr, 10)),
        "p90_block_mean": float(np.percentile(arr, 90)),
    }


def xstar_control(survivors: list[dict]) -> dict:
    xs = [abs(s["x_star"]) for s in survivors if s["x_star"] is not None]
    if not xs:
        return {"n": 0}
    arr = np.array(xs, dtype=float)
    L = len(survivors[0]["word"])
    return {
        "n": len(arr),
        "L": L,
        "max_|x*|": float(arr.max()),
        "mean_|x*|": float(arr.mean()),
        "p99_|x*|": float(np.percentile(arr, 99)),
        "max_|x*|/2^L": float(arr.max() / (1 << L)),
        "max_|x*|/L**3": float(arr.max() / (L**3)),
        "frac_|x*|_le_2^L": float(np.mean(arr <= (1 << L))),
    }


def main() -> int:
    # Soft C damit Survivor-Graph nicht leer; kritisch nahe C=12..14
    Cs = [8.0, 10.0, 12.0]
    levels = [8, 9, 10, 11, 12, 13, 14]
    print("=" * 72)
    print("SURVIVOR-TRANSFER — Φ, Blockdrift ΔΛ, Spektralradius")
    print("=" * 72)
    print(
        "Ziel: liminf Λ<0 entlang Survivor-Pfaden; "
        "kein punktweises Λ≤−εL.\n"
    )

    all_payload = {}
    for C in Cs:
        print(f"\n### C = {C}")
        by_L = {L: enumerate_survivors(L, C) for L in levels}
        for L in levels:
            print(f"  |S_{L}| = {len(by_L[L])}")

        # Extension edges L → L+1
        ext_stats = []
        for L in levels[:-1]:
            g = extension_graph(by_L[L], by_L[L + 1])
            bd = block_drift_stats(g["edges"])
            # Φ-Rekurrenz: mean (Φ(child)-Λ) already in edges
            Phi_vals = [v["Phi"] for v in g["Phi"].values() if v["Phi"] is not None]
            mean_d = [
                v["mean_dLambda"]
                for v in g["Phi"].values()
                if v["mean_dLambda"] is not None
            ]
            rec = {
                "L": L,
                **bd,
                "mean_Phi": float(np.mean(Phi_vals)) if Phi_vals else None,
                "mean_node_mean_dLambda": float(np.mean(mean_d)) if mean_d else None,
                "frac_nodes_all_children_neg_d": (
                    float(
                        np.mean(
                            [
                                v["max_dLambda"] is not None and v["max_dLambda"] < 0
                                for v in g["Phi"].values()
                            ]
                        )
                    )
                    if g["Phi"]
                    else None
                ),
            }
            ext_stats.append(rec)
            print(
                f"  ext L={L}→{L+1}: edges={bd.get('n_edges',0)}  "
                f"meanΔΛ={bd.get('mean_dLambda')}  "
                f"fracΔΛ<0={bd.get('frac_neg')}  "
                f"frac nodes maxΔΛ<0={rec['frac_nodes_all_children_neg_d']}"
            )

        # Shift spectral radius at each L
        spec = []
        for L in levels:
            S = by_L[L]
            if len(S) == 0:
                spec.append({"L": L, "n": 0, "rho": 0.0, "rho_per_2": 0.0})
                continue
            M, words = shift_adjacency(S)
            rho = spectral_radius(M)
            spec.append(
                {
                    "L": L,
                    "n": len(S),
                    "rho": rho,
                    "rho_per_2": rho / 2.0,  # vs volle Verzweigung 2
                    "log2_rho_over_nbits": math.log2(rho) / L if rho > 0 else None,
                }
            )
            print(
                f"  shift L={L}: |V|={len(S)}  ρ(M)={rho:.4f}  "
                f"ρ/2={rho/2:.4f}"
            )

        # Block means m=3,4 from L0=8,10
        blocks = {}
        for L0 in (8, 10):
            for m in (3, 4):
                if L0 + m not in by_L:
                    continue
                key = f"L0={L0}_m={m}"
                blocks[key] = path_block_means(by_L, L0, m)
                b = blocks[key]
                print(
                    f"  block {key}: paths={b.get('n_paths')}  "
                    f"mean={b.get('mean_of_block_means')}  "
                    f"frac<0={b.get('frac_block_mean_neg')}  "
                    f"frac≤-0.1={b.get('frac_block_mean_le_-0.1')}  "
                    f"max={b.get('max_block_mean')}"
                )

        # x* control on survivors
        xctrl = {L: xstar_control(by_L[L]) for L in levels}
        print("  |x*| control (Survivor):")
        for L, xc in xctrl.items():
            if xc.get("n", 0) == 0:
                continue
            print(
                f"    L={L}: max|x*|={xc['max_|x*|']:.3g}  "
                f"max/2^L={xc['max_|x*|/2^L']:.3g}  "
                f"max/L^3={xc['max_|x*|/L**3']:.3g}"
            )

        all_payload[str(C)] = {
            "n_S": {str(L): len(by_L[L]) for L in levels},
            "extension": ext_stats,
            "shift_spectral": spec,
            "block_means": blocks,
            "x_star": xctrl,
        }

    print("\n[Lesart]")
    print(
        "  • mean ΔΛ < 0 und Blockmittel < 0 ⇒ Drifteliminierung (dynamisch)\n"
        "  • ρ(M)/2 < 1 ⇒ Survivor-Shift verdünnt vs volle Entropie\n"
        "  • |x*| ≲ 2^L typisch; für R→0 braucht man Λ→−∞ stärker als log|x*|\n"
        "  • punktweises Φ(ω)≤Φ−δ pro Schritt erwartet NICHT; Blockmittel ja"
    )

    out = Path(
        "/Users/thomashoffbauer/Projects/Kepler-Hurrwitz/docs/exports/"
        "collatz_survivor_transfer.json"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "claim": (
            "Survivor transfer graph: block-averaged dLambda and spectral radius; "
            "target liminf Lambda<0 on infinite survivor paths; not Collatz proof"
        ),
        "target_theorems": [
            "forall survivor omega: liminf Lambda_L < 0",
            "stronger: liminf Lambda_L/L <= -eps",
            "aux: |x_star_L| subexponential / poly(L) on survivor prefixes",
        ],
        "gap_analogy": (
            "no single step need contract; every long survivor block has "
            "negative mean resonance — spectral gap of survivor transfer"
        ),
        "data": all_payload,
    }
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")

    # also desktop copy
    desk = Path(__file__).resolve().parent / "collatz_survivor_transfer.json"
    try:
        desk.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {desk}")
    except OSError as e:
        print(f"(desktop write skipped: {e})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
