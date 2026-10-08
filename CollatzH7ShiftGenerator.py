"""
CollatzH7ShiftGenerator — Shift-Graph M_w auf Clog-Survivors.

Kante:  u → v  ⇔  u[1:] = v[:-1]
Gewicht: w(b) = 1/2 (b=0), 3/2 (b=1)

Vorbedingung: |S_L(C)| > 0, sonst Abbruch ohne Matrix.

Export: JSON (+ CSV Knotenliste) für Lean-SurvivorNode-Mapping.
Epistemik: numerische Basis; ρ<1 ist nicht garantiert.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


W = {"0": 0.5, "1": 1.5}


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


def word_of_residue(r: int, L: int) -> str:
    mod = 1 << L
    x, bits = r, []
    for _ in range(L):
        if x % 2 == 0:
            bits.append("0")
            x = (x // 2) % mod
        else:
            bits.append("1")
            x = ((3 * x + 1) // 2) % mod
    return "".join(bits)


def is_survivor_clog(r: int, L: int, C: float) -> bool:
    n = r if r > 1 else r + (1 << L)
    budget = int(max(C * math.log2(n) * 3 + 50, 500))
    ht = hit_time(n, budget)
    return ht is None or ht > C * math.log2(n) + 2


def enumerate_survivors(L: int, C: float) -> list[dict]:
    """Ein Knoten pro Survivor-Wort (kanonisch: kleinstes Residuum)."""
    best: dict[str, dict] = {}
    for r in range(1 << L):
        if not is_survivor_clog(r, L, C):
            continue
        w = word_of_residue(r, L)
        n = r if r > 1 else r + (1 << L)
        a = w.count("1")
        Lam = a * math.log2(3) - L
        rec = {
            "word": w,
            "residue": r,
            "n_rep": n,
            "a": a,
            "Lambda": Lam,
            "Lambda_over_L": Lam / L,
        }
        if w not in best or r < best[w]["residue"]:
            best[w] = rec
    # stabile Indizierung: lexikographisch nach Wort
    nodes = [best[w] for w in sorted(best)]
    for i, node in enumerate(nodes):
        node["index"] = i
    return nodes


def build_Mw(nodes: list[dict]) -> tuple[np.ndarray, list[dict]]:
    words = [n["word"] for n in nodes]
    idx = {w: i for i, w in enumerate(words)}
    n = len(words)
    M = np.zeros((n, n), dtype=float)
    edges = []
    for u in words:
        for b, wgt in W.items():
            v = u[1:] + b
            if v in idx:
                i, j = idx[u], idx[v]
                M[i, j] = wgt
                edges.append(
                    {
                        "i": i,
                        "j": j,
                        "u": u,
                        "v": v,
                        "bit": b,
                        "weight": wgt,
                        "weight_rat": "1/2" if b == "0" else "3/2",
                    }
                )
    return M, edges


def spectral_audit(M: np.ndarray) -> dict:
    if M.size == 0:
        return {"n": 0, "rho": 0.0}
    vals, vecs = np.linalg.eig(M)
    k = int(np.argmax(np.abs(vals)))
    rho = float(np.abs(vals[k]))
    v = np.real(vecs[:, k])
    if np.min(v) < -1e-12:
        v = -v
    if np.max(np.abs(v)) > 0:
        v = v / np.max(np.abs(v))
    # prüfen Mv ≈ ρ v
    Mv = M @ v
    resid = float(np.max(np.abs(Mv - rho * v))) if rho > 0 else float(np.max(np.abs(Mv)))
    row_sums = M.sum(axis=1)
    return {
        "n": int(M.shape[0]),
        "rho": rho,
        "rho_lt_1": rho < 1.0 - 1e-12,
        "row_sum_max": float(row_sums.max()) if len(row_sums) else None,
        "row_dissipative": bool(np.all(row_sums < 1.0 - 1e-12)) if len(row_sums) else False,
        "dom_eigenvalue_real": float(np.real(vals[k])),
        "dom_eigenvalue_imag": float(np.imag(vals[k])),
        "right_evec_max1": v.tolist(),
        "Mv_minus_rho_v_maxabs": resid,
    }


def export_bundle(
    L: int,
    C: float,
    nodes: list[dict],
    M: np.ndarray,
    edges: list[dict],
    audit: dict,
    out_dir: Path,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"h7_shift_L{L}_C{str(C).replace('.', 'p')}"
    json_path = out_dir / f"{stem}.json"
    csv_path = out_dir / f"{stem}_nodes.csv"
    mat_csv = out_dir / f"{stem}_Mw.csv"

    payload = {
        "claim": (
            "H7/Clog survivor shift graph M_w; edge u->v iff u[1:]==v[:-1]; "
            "weights 1/2 and 3/2; not a Collatz proof; rho<1 not guaranteed"
        ),
        "L": L,
        "C": C,
        "rule": "clog_HitTime > C*log2(n)+2",
        "edge_rule": "u[1:] == v[:-1]",
        "weights": {"0": "1/2", "1": "3/2"},
        "n_nodes": len(nodes),
        "nodes": nodes,
        "edges": edges,
        "M_w": M.tolist(),
        "spectral": audit,
        "lean_indexing_recommendation": (
            "Abstract: SurvivorNode L := {w : BitVec L // isSurvivor C w}; "
            "Concrete audit Fin n with explicit survivors list (this JSON); "
            "matrix as Matrix (Fin n) (Fin n) Rat"
        ),
    }
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["index", "word", "residue", "n_rep", "a", "Lambda", "Lambda_over_L"],
        )
        w.writeheader()
        for node in nodes:
            w.writerow(node)

    np.savetxt(mat_csv, M, delimiter=",", fmt="%.10g")
    return {
        "json": str(json_path),
        "nodes_csv": str(csv_path),
        "Mw_csv": str(mat_csv),
        **audit,
        "n_nodes": len(nodes),
    }




def emit_lean_fragment(nodes, M, L, C) -> str:
    """Optional Lean Matrix dump for small |V| (<= 8)."""
    n = len(nodes)
    if n == 0 or n > 8:
        return f"-- skip Lean emit: n={n}\n"
    def rat(x: float) -> str:
        if abs(x) < 1e-15:
            return "0"
        if abs(x - 0.5) < 1e-12:
            return "(1 : Rat) / 2"
        if abs(x - 1.5) < 1e-12:
            return "(3 : Rat) / 2"
        return repr(x)
    lines = [
        f"-- Auto from CollatzH7ShiftGenerator L={L} C={C} n={n}",
        f"def Mw_L{L}_C{str(C).replace('.', 'p')} : Matrix (Fin {n}) (Fin {n}) Rat :=",
        "  !![",
    ]
    for i, row in enumerate(M):
        ents = ", ".join(rat(float(x)) for x in row)
        semi = ";" if i + 1 < n else ""
        lines.append(f"    {ents}{semi}")
    lines.append("  ]")
    for node in nodes:
        lines.append(f"-- {node['index']}: {node['word']} residue={node['residue']}")
    return "\n".join(lines) + "\n"

def main() -> int:
    ap = argparse.ArgumentParser(description="H7 shift M_w generator")
    ap.add_argument("--L", type=int, default=8)
    ap.add_argument("--C", type=float, default=8.0)
    ap.add_argument(
        "--out",
        type=str,
        default=str(
            Path("/Users/thomashoffbauer/Projects/Kepler-Hurrwitz/docs/exports")
        ),
    )
    args = ap.parse_args()
    L, C = args.L, args.C
    out_dir = Path(args.out)

    print("=" * 72)
    print(f"H7 SHIFT-GENERATOR  L={L}  C={C}")
    print("  Kante: u[1:]==v[:-1]   w∈{1/2, 3/2}")
    print("=" * 72)

    nodes = enumerate_survivors(L, C)
    print(f"|S_L| (unique words) = {len(nodes)}")
    if len(nodes) == 0:
        print("VORBEDINGUNG VERLETZT: V_L = ∅ — keine Matrix, kein Lean-Dump.")
        # minimal manifest
        man = {
            "L": L,
            "C": C,
            "n_nodes": 0,
            "status": "empty_survivor_set",
            "claim": "no matrix exported",
        }
        out_dir.mkdir(parents=True, exist_ok=True)
        p = out_dir / f"h7_shift_L{L}_C{str(C).replace('.', 'p')}_EMPTY.json"
        p.write_text(json.dumps(man, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {p}")
        return 2

    M, edges = build_Mw(nodes)
    audit = spectral_audit(M)
    meta = export_bundle(L, C, nodes, M, edges, audit, out_dir)

    print(f"edges = {len(edges)}")
    print(f"ρ(M_w) = {audit['rho']:.6f}  (<1? {audit['rho_lt_1']})")
    print(f"row_sum_max = {audit['row_sum_max']}  row_dissipative? {audit['row_dissipative']}")
    print(f"Mv−ρv max|·| = {audit['Mv_minus_rho_v_maxabs']:.3e}")
    print(f"wrote {meta['json']}")
    print(f"wrote {meta['nodes_csv']}")
    print(f"wrote {meta['Mw_csv']}")
    if len(nodes) <= 8:
        lean = emit_lean_fragment(nodes, M, L, C)
        lean_path = out_dir / f"h7_shift_L{L}_C{str(C).replace('.', 'p')}.leanfrag"
        lean_path.write_text(lean, encoding="utf-8")
        print(f"wrote {lean_path}")
        print(lean)
    if not audit["rho_lt_1"]:
        print(
            "\nHINWEIS: ρ≥1 — numerische Basis für Lean-Struktur ja, "
            "Perron-Zertifikat θ<1 nein. Nicht als h7_dissipativ behaupten."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
