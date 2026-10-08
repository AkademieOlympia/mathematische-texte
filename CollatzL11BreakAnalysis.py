"""
CollatzL11BreakAnalysis — Bruch der Komposit-Kontraktion bei L=11, C=8.

Vergleicht L=10 (Gap hält) vs L=11 (ρ(M3 M2), ρ(M3 M5) > 1):
  * Graph-Statistik V_L, nnz, SCC
  * Zyklen in M_comp = M3 @ M2 (und optional M3 @ M5)
  * Instabilitäts-Knoten: größte Verletzung von Mv ≤ θ v
    für θ∈{1/2, 19/20, 1} und Perron-/Eins-Vektoren

Epistemik: Strukturdiagnose; kein Collatz-Beweis.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from CollatzBlockWeightedTransfer import (
    adjacency_list,
    enumerate_survivor_words,
    get_block_matrix,
    hit_time,
    nilpotency_index,
    spectral,
    survivor_word_set,
    word_of_residue,
)


def digraph_edges(M: np.ndarray) -> list[tuple[int, int, float]]:
    ii, jj = np.where(M > 1e-15)
    return [(int(i), int(j), float(M[i, j])) for i, j in zip(ii, jj)]


def scc_count(M: np.ndarray) -> dict:
    """Tarjan SCC on support of M>0."""
    n = M.shape[0]
    adj = [[] for _ in range(n)]
    for i, j, _ in digraph_edges(M):
        adj[i].append(j)
    index = 0
    stack: list[int] = []
    onstack = [False] * n
    indices = [-1] * n
    lowlink = [0] * n
    sccs: list[list[int]] = []

    def strongconnect(v: int) -> None:
        nonlocal index
        indices[v] = index
        lowlink[v] = index
        index += 1
        stack.append(v)
        onstack[v] = True
        for w in adj[v]:
            if indices[w] == -1:
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif onstack[w]:
                lowlink[v] = min(lowlink[v], indices[w])
        if lowlink[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                onstack[w] = False
                comp.append(w)
                if w == v:
                    break
            sccs.append(comp)

    for v in range(n):
        if indices[v] == -1 and adj[v]:
            strongconnect(v)
        elif indices[v] == -1:
            # isolated
            pass
    nontrivial = [c for c in sccs if len(c) > 1 or (len(c) == 1 and M[c[0], c[0]] > 0)]
    return {
        "n_scc_total": len(sccs),
        "n_scc_nontrivial": len(nontrivial),
        "largest_scc": max((len(c) for c in sccs), default=0),
        "nontrivial_sizes": sorted((len(c) for c in nontrivial), reverse=True)[:10],
    }


def find_cycles(
    M: np.ndarray, nodes: list[str], max_len: int = 8, max_report: int = 15
) -> list[dict]:
    """DFS cycles; report top by geometric mean edge weight."""
    n = M.shape[0]
    adj = [[] for _ in range(n)]
    for i, j, w in digraph_edges(M):
        adj[i].append((j, w))
    best: list[tuple[float, float, tuple[int, ...]]] = []

    # length-1: genuine self-loops (avoid DFS double-count s→s→s)
    for i in range(n):
        if M[i, i] > 1e-15:
            w = float(M[i, i])
            best.append((w, w, (i, i)))

    def dfs(start: int, u: int, prod: float, path: tuple[int, ...], seen: set[int]) -> None:
        if len(path) > max_len:
            return
        for v, w in adj[u]:
            nw = prod * w
            if v == start and len(path) >= 2:
                L = len(path)
                geom = nw ** (1 / L)
                best.append((geom, nw, path + (v,)))
            elif v != start and v not in seen:
                dfs(start, v, nw, path + (v,), seen | {v})

    starts = [i for i in range(n) if adj[i]]
    for s in starts:
        dfs(s, s, 1.0, (s,), {s})

    best_sorted = sorted(best, key=lambda t: -t[0])
    # unique by frozenset of nodes + length
    seen_sig = set()
    out = []
    for geom, prod, path in best_sorted:
        sig = (len(path) - 1, frozenset(path[:-1]))
        if sig in seen_sig:
            continue
        seen_sig.add(sig)
        out.append(
            {
                "length": len(path) - 1,
                "geom_mean_weight": float(geom),
                "product_weight": float(prod),
                "indices": list(path[:-1]),
                "words": [nodes[i] for i in path[:-1]],
            }
        )
        if len(out) >= max_report:
            break
    return out


def instability_nodes(
    M: np.ndarray,
    nodes: list[str],
    theta: float,
    v: np.ndarray,
    top_k: int = 20,
) -> list[dict]:
    """Nodes maximizing (Mv)_i - theta v_i (positive = certificate violation)."""
    Mv = M @ v
    slack = Mv - theta * v
    order = np.argsort(-slack)
    rows = []
    for i in order[:top_k]:
        if slack[i] <= 1e-12 and len(rows) >= 5:
            break
        rows.append(
            {
                "index": int(i),
                "word": nodes[i],
                "v_i": float(v[i]),
                "Mv_i": float(Mv[i]),
                "theta_v_i": float(theta * v[i]),
                "slack": float(slack[i]),
                "ratio_Mv_over_v": float(Mv[i] / v[i]) if v[i] > 1e-15 else None,
            }
        )
    return rows


def build_composite(nodes, C, L, mB, mA, cache):
    MB = get_block_matrix(nodes, C, L, mB, cache)
    MA = get_block_matrix(nodes, C, L, mA, cache)
    return MB @ MA, MB, MA


def analyze_level(C: float, L: int, pairs: list[tuple[int, int]]) -> dict:
    nodes = enumerate_survivor_words(L, C)
    cache: dict[int, np.ndarray] = {}
    Mw = get_block_matrix(nodes, C, L, 1, cache)
    spw = spectral(Mw)
    level = {
        "C": C,
        "L": L,
        "n_V": len(nodes),
        "rho_Mw": spw["rho"],
        "nnz_Mw": spw["nnz"],
        "Mw_scc": scc_count(Mw),
        "composites": {},
    }
    for mB, mA in pairs:
        key = f"M{mB}_M{mA}"
        Mcomp, _, _ = build_composite(nodes, C, L, mB, mA, cache)
        sp = spectral(Mcomp)
        nil_k = nilpotency_index(Mcomp)
        v_per = np.array(sp["v"], float)
        if v_per.max() > 0:
            v_per = v_per / v_per.max()
        v_ones = np.ones(len(nodes))
        cycles = find_cycles(Mcomp, nodes, max_len=6, max_report=12)
        diag = np.diag(Mcomp)
        self_loops = [
            {"index": int(i), "word": nodes[i], "weight": float(diag[i])}
            for i in np.where(diag > 1e-12)[0]
        ]
        # instability vs several thetas
        inst = {}
        for th, v, tag in (
            (0.5, v_ones, "theta_half_ones"),
            (0.95, v_per, "theta_095_perron"),
            (1.0, v_per, "theta_1_perron"),
        ):
            inst[tag] = {
                "theta": th,
                "max_slack": float((Mcomp @ v - th * v).max()),
                "n_violations": int(np.sum(Mcomp @ v > th * v + 1e-12)),
                "top_nodes": instability_nodes(Mcomp, nodes, th, v, top_k=15),
            }
        level["composites"][key] = {
            "mB": mB,
            "mA": mA,
            "rho": sp["rho"],
            "rho_lt_1": sp["rho"] < 1,
            "nilpotency_index": nil_k,
            "nnz": int(np.count_nonzero(Mcomp)),
            "scc": scc_count(Mcomp),
            "self_loops": self_loops,
            "top_cycles": cycles,
            "instability": inst,
            "perron_support": int((v_per > 1e-10).sum()),
        }
    return level


def compare_break(C: float, L_ok: int, L_break: int, out_dir: Path) -> dict:
    pairs = [(3, 2), (3, 5)]
    print("=" * 72)
    print(f"L-BREAK ANALYSIS  C={C}  L_ok={L_ok} vs L_break={L_break}")
    print("=" * 72)
    a = analyze_level(C, L_ok, pairs)
    b = analyze_level(C, L_break, pairs)

    def summarize(level: dict, label: str) -> None:
        print(f"\n--- {label} L={level['L']} |V|={level['n_V']} ρ(Mw)={level['rho_Mw']:.4f} ---")
        print(f"  Mw SCC nontrivial: {level['Mw_scc']['n_scc_nontrivial']}")
        for key, comp in level["composites"].items():
            print(
                f"  {key}: ρ={comp['rho']:.4f} nnz={comp['nnz']} nil={comp['nilpotency_index']} "
                f"SCC_nt={comp['scc']['n_scc_nontrivial']} cyc0_geom="
                f"{comp['top_cycles'][0]['geom_mean_weight'] if comp['top_cycles'] else '—'}"
            )
            top = comp["instability"]["theta_1_perron"]["top_nodes"][:5]
            print("    top slack vs θ=1 (Perron v):")
            for row in top:
                print(
                    f"      {row['word']}  slack={row['slack']:.4f}  "
                    f"Mv/v={row['ratio_Mv_over_v']}"
                )

    summarize(a, "OK")
    summarize(b, "BREAK")

    # word-level: new nodes at L_break that appear in top cycles / top slack
    payload = {
        "claim": (
            "L=11 break diagnosis vs L=10: cycles and instability nodes of M3M2/M3M5; "
            "not a Collatz proof"
        ),
        "C": C,
        "L_ok": L_ok,
        "L_break": L_break,
        "level_ok": a,
        "level_break": b,
        "diff_highlights": {},
    }
    for key in ("M3_M2", "M3_M5"):
        okc = a["composites"][key]
        brc = b["composites"][key]
        payload["diff_highlights"][key] = {
            "rho_ok": okc["rho"],
            "rho_break": brc["rho"],
            "nnz_ok": okc["nnz"],
            "nnz_break": brc["nnz"],
            "scc_nt_ok": okc["scc"]["n_scc_nontrivial"],
            "scc_nt_break": brc["scc"]["n_scc_nontrivial"],
            "break_top_cycle": brc["top_cycles"][:3],
            "break_top_instability_theta1": brc["instability"]["theta_1_perron"]["top_nodes"][:10],
            "break_top_instability_theta_half_ones": brc["instability"]["theta_half_ones"][
                "top_nodes"
            ][:10],
        }
        print(f"\n[DIFF {key}] ρ {okc['rho']:.4f} → {brc['rho']:.4f}  "
              f"nnz {okc['nnz']}→{brc['nnz']}  SCC_nt {okc['scc']['n_scc_nontrivial']}→"
              f"{brc['scc']['n_scc_nontrivial']}")

    out = out_dir / f"h7_L{L_break}_break_vs_L{L_ok}_C{str(C).replace('.','p')}.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    # also short node list CSV-like
    rows = []
    for key in ("M3_M2", "M3_M5"):
        for row in b["composites"][key]["instability"]["theta_1_perron"]["top_nodes"]:
            rows.append({"composite": key, **row})
    nodes_out = out_dir / f"h7_L{L_break}_instability_nodes_C{str(C).replace('.','p')}.json"
    nodes_out.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    print(f"wrote {nodes_out}")
    return payload


def clog_probe(word: str, C: float) -> dict:
    """Clog-Status für ein Bitwort: Residuen, die dieses Wort erzeugen, + HitTime."""
    L = len(word)
    hits = []
    for r in range(1 << L):
        if word_of_residue(r, L) != word:
            continue
        n = r if r > 1 else r + (1 << L)
        budget = int(max(C * __import__("math").log2(n) * 3 + 50, 500))
        ht = hit_time(n, budget)
        thr = C * __import__("math").log2(n) + 2
        surv = ht is None or ht > thr
        hits.append(
            {
                "r": r,
                "n": n,
                "hit_time": ht,
                "threshold": float(thr),
                "is_survivor": surv,
                "margin": None if ht is None else float(ht - thr),
            }
        )
    any_surv = any(h["is_survivor"] for h in hits)
    return {
        "word": word,
        "L": L,
        "n_residues": len(hits),
        "in_survivor_set": any_surv,
        "residues": hits,
        "tightest_survivor": min(
            (h for h in hits if h["is_survivor"] and h["margin"] is not None),
            key=lambda h: h["margin"],
            default=None,
        ),
    }


def enumerate_block_walks(
    nodes: list[str], start: str, m: int, C: float, L: int
) -> list[dict]:
    """Alle m-Walks von start, die den L+m-Survivor-Filter bestehen."""
    idx = {w: i for i, w in enumerate(nodes)}
    if start not in idx:
        return []
    adj = adjacency_list(nodes)
    surv_Lm = survivor_word_set(L + m, C)
    out: list[dict] = []

    def rec(i: int, depth: int, path: str, wprod: float, chain: list[str]) -> None:
        if depth == m:
            full = start + path
            kept = full in surv_Lm
            out.append(
                {
                    "path_bits": path,
                    "weight": float(wprod),
                    "end": chain[-1],
                    "chain": list(chain),
                    "filter_word": full,
                    "filter_kept": kept,
                    "filter_clog": clog_probe(full, C) if kept or depth <= 8 else None,
                }
            )
            return
        for j, b, wgt in adj[i]:
            rec(j, depth + 1, path + b, wprod * wgt, chain + [nodes[j]])

    rec(idx[start], 0, "", 1.0, [start])
    return out


def extract_path_details(
    start_node: str,
    nodes: list[str],
    C: float,
    L: int,
    mA: int,
    mB: int,
) -> dict:
    """
    Dekonstruktion der Selbstschleife (M_mB @ M_mA)_{ωω}.

    Matrix: (MB MA)_{ii} = Σ_ν MB[i,ν] MA[ν,i]
      ⇒ Pfad ω --mB--> ν --mA--> ω
      (Operator auf Spalten: zuerst MA=M2/M5, dann MB=M3).

    Für Diagnose der Filter-Lücke zählen wir trotzdem beide Lesarten
    der Walk-Kette; Beiträge matchen exakt die Diagonale von MB@MA.
    """
    idx = {w: i for i, w in enumerate(nodes)}
    if start_node not in idx:
        return {"error": "start not in V_L", "start": start_node}

    cache: dict[int, np.ndarray] = {}
    MA = get_block_matrix(nodes, C, L, mA, cache)
    MB = get_block_matrix(nodes, C, L, mB, cache)
    i0 = idx[start_node]
    Mcomp = MB @ MA
    diag = float(Mcomp[i0, i0])

    # Intermediate ν mit MB[i0,ν] MA[ν,i0] > 0
    mids = [
        nodes[k]
        for k in range(len(nodes))
        if MB[i0, k] > 1e-15 and MA[k, i0] > 1e-15
    ]

    contributions = []
    for mid in mids:
        im = idx[mid]
        walks_B = [
            w
            for w in enumerate_block_walks(nodes, start_node, mB, C, L)
            if w["filter_kept"] and w["end"] == mid
        ]
        walks_A = [
            w
            for w in enumerate_block_walks(nodes, mid, mA, C, L)
            if w["filter_kept"] and w["end"] == start_node
        ]
        for wb in walks_B:
            for wa in walks_A:
                contrib = float(wb["weight"] * wa["weight"])
                # Operator-Reihenfolge (zuerst MA, dann MB): Bits A dann B
                # Matrix-Faktorisierung MB@MA: Walk B dann A am Zustand
                total_bits_state = wb["path_bits"] + wa["path_bits"]
                total_bits_op = wa["path_bits"] + wb["path_bits"]
                prefix_clog = []
                for t in range(1, len(total_bits_state) + 1):
                    window = (start_node + total_bits_state[:t])[-L:]
                    prefix_clog.append(
                        {
                            "t": t,
                            "bit": total_bits_state[t - 1],
                            "window_L": window,
                            "window_in_VL": window in idx,
                            "clog": clog_probe(window, C),
                        }
                    )
                contributions.append(
                    {
                        "intermediate": mid,
                        "mB_path": wb["path_bits"],
                        "mA_path": wa["path_bits"],
                        "total_path_bits_state_order_BmA": total_bits_state,
                        "total_path_bits_operator_AmB": total_bits_op,
                        "weight_A": wa["weight"],
                        "weight_B": wb["weight"],
                        "weight_product": contrib,
                        "block_weights": {
                            "MA_mid_start": float(MA[im, i0]),
                            "MB_start_mid": float(MB[i0, im]),
                        },
                        "chain_B": wb["chain"],
                        "chain_A": wa["chain"],
                        "filter_word_B": wb["filter_word"],
                        "filter_word_A": wa["filter_word"],
                        "filter_clog_A": clog_probe(wa["filter_word"], C),
                        "filter_clog_B": clog_probe(wb["filter_word"], C),
                        "prefix_windows_clog": prefix_clog,
                    }
                )

    # L=10-Vergleich: Fenster / Trunkierungen
    L_ok = L - 1
    V_ok = set(enumerate_survivor_words(L_ok, C)) if L_ok >= 2 else set()
    trunc = {
        "prefix_L_ok": start_node[:L_ok],
        "suffix_L_ok": start_node[-L_ok:],
        "shift_windows_L_ok": [start_node[i : i + L_ok] for i in range(L - L_ok + 1)],
    }
    trunc["in_V_ok"] = {w: w in V_ok for w in trunc["shift_windows_L_ok"]}
    trunc["clog_windows"] = {w: clog_probe(w, C) for w in trunc["shift_windows_L_ok"]}

    # fehlender Schnitt: Filterwörter bei L+mA / L+mB, deren L_ok-Fenster nicht Survivor sind
    missing_cuts = []
    for c in contributions:
        for tag, fword in (("A", c["filter_word_A"]), ("B", c["filter_word_B"])):
            for i in range(len(fword) - L_ok + 1):
                win = fword[i : i + L_ok]
                if win not in V_ok:
                    missing_cuts.append(
                        {
                            "block": tag,
                            "filter_word": fword,
                            "window_L_ok": win,
                            "offset": i,
                            "clog_L_ok": clog_probe(win, C),
                        }
                    )

    return {
        "start": start_node,
        "L": L,
        "C": C,
        "mA": mA,
        "mB": mB,
        "composite_diag": diag,
        "n_contributions": len(contributions),
        "weight_sum": float(sum(c["weight_product"] for c in contributions)),
        "contributions": contributions,
        "truncation_vs_L_ok": trunc,
        "missing_cuts_vs_L_ok": missing_cuts,
        "start_clog": clog_probe(start_node, C),
    }


# Die vier Resonanz-Selbstschleifen aus der L=11-Diagnose
RESONANCE_TARGETS = {
    "M3_M2": {
        "mB": 3,
        "mA": 2,
        "words": ["10111101111", "11011110111"],
    },
    "M3_M5": {
        "mB": 3,
        "mA": 5,
        "words": ["10111110101", "11110101111"],
    },
}


def deconstruct_resonance(C: float, L: int, out_dir: Path) -> dict:
    nodes = enumerate_survivor_words(L, C)
    print("=" * 72)
    print(f"PATH DECONSTRUCT  C={C} L={L}  |V|={len(nodes)}")
    print("=" * 72)
    payload = {
        "claim": (
            "Pfad-Dekonstruktion der vier L=11-Selbstschleifen in M3M2/M3M5; "
            "Clog-Schnitt vs L=10-Fenster; kein Collatz-Beweis"
        ),
        "C": C,
        "L": L,
        "targets": {},
    }
    for key, spec in RESONANCE_TARGETS.items():
        mB, mA = spec["mB"], spec["mA"]
        print(f"\n### {key}  (M{mB} @ M{mA})")
        block = {}
        for w in spec["words"]:
            det = extract_path_details(w, nodes, C, L, mA, mB)
            block[w] = det
            print(
                f"  {w}: diag={det.get('composite_diag')}  "
                f"contribs={det.get('n_contributions')}  "
                f"wsum={det.get('weight_sum')}  "
                f"missing_cuts={len(det.get('missing_cuts_vs_L_ok', []))}"
            )
            for c in det.get("contributions", [])[:5]:
                print(
                    f"    ν={c['intermediate']}  "
                    f"bits(B|A)={c['total_path_bits_state_order_BmA']}  "
                    f"w={c['weight_product']:.6f}  "
                    f"B={c['mB_path']} A={c['mA_path']}"
                )
                print(f"      chainB={c['chain_B']}")
                print(f"      chainA={c['chain_A']}")
            # highlight missing cuts unique
            cuts = det.get("missing_cuts_vs_L_ok", [])
            uniq = {}
            for mc in cuts:
                uniq.setdefault(mc["window_L_ok"], mc)
            if uniq:
                print(f"    fehlende L={L-1}-Fenster ({len(uniq)}):")
                for win, mc in list(uniq.items())[:8]:
                    cl = mc["clog_L_ok"]
                    ht = None
                    if cl["residues"]:
                        # prefer non-survivor residue info
                        bad = [r for r in cl["residues"] if not r["is_survivor"]]
                        sample = (bad or cl["residues"])[0]
                        ht = sample["hit_time"]
                    print(
                        f"      {win}  in_surv_set={cl['in_survivor_set']}  "
                        f"sample_ht={ht}"
                    )
        payload["targets"][key] = block

    out = out_dir / f"h7_L{L}_path_deconstruct_C{str(C).replace('.','p')}.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return payload


def contains_forbidden(word: str, forbidden: set[str]) -> bool:
    for f in forbidden:
        if len(f) > len(word):
            continue
        for i in range(len(word) - len(f) + 1):
            if word[i : i + len(f)] == f:
                return True
    return False


def filter_nodes(nodes: list[str], forbidden: set[str]) -> tuple[list[str], list[str]]:
    kept, killed = [], []
    for w in nodes:
        if contains_forbidden(w, forbidden):
            killed.append(w)
        else:
            kept.append(w)
    return kept, killed


def block_matrix_with_forbidden(
    nodes: list[str], C: float, L: int, m: int, forbidden: set[str]
) -> tuple[np.ndarray, dict]:
    """
    Wie block_matrix_walk_filter, plus Virtual-Cut:
      full = start+path muss Survivor bei L+m sein UND kein Forbidden-Fenster enthalten.
    Für m=1: Basis-Mw auf dem bereits beschnittenen V_L.
    """
    from CollatzBlockWeightedTransfer import adjacency_list, build_base_Mw

    if m == 1:
        M, _ = build_base_Mw(nodes)
        return M, {"mode": "base", "forbidden": sorted(forbidden)}

    n = len(nodes)
    adj = adjacency_list(nodes)
    M = np.zeros((n, n), dtype=float)
    surv_Lm = survivor_word_set(L + m, C)
    stats = {
        "walks_total": 0,
        "walks_kept": 0,
        "walks_rejected_clog": 0,
        "walks_rejected_forbidden": 0,
        "n_survivors_Lpm": len(surv_Lm),
    }

    def rec(i0: int, i: int, depth: int, path: str, wprod: float) -> None:
        if depth == m:
            stats["walks_total"] += 1
            full = nodes[i0] + path
            if full not in surv_Lm:
                stats["walks_rejected_clog"] += 1
                return
            if contains_forbidden(full, forbidden):
                stats["walks_rejected_forbidden"] += 1
                return
            stats["walks_kept"] += 1
            M[i0, i] += wprod
            return
        for j, b, wgt in adj[i]:
            rec(i0, j, depth + 1, path + b, wprod * wgt)

    for i0 in range(n):
        rec(i0, i0, 0, "", 1.0)
    return M, stats


def composite_spectrum_cut(
    C: float, L: int, forbidden: set[str], pairs: list[tuple[int, int]]
) -> dict:
    nodes0 = enumerate_survivor_words(L, C)
    nodes, killed = filter_nodes(nodes0, forbidden)
    killed_resonance = {
        key: [w for w in spec["words"] if w in killed or contains_forbidden(w, forbidden)]
        for key, spec in RESONANCE_TARGETS.items()
    }
    out = {
        "L": L,
        "C": C,
        "forbidden": sorted(forbidden),
        "n_V_raw": len(nodes0),
        "n_V_cut": len(nodes),
        "n_killed": len(killed),
        "killed_words": killed,
        "killed_resonance_targets": killed_resonance,
        "composites": {},
        "block_stats": {},
    }
    if not nodes:
        out["error"] = "V_cut empty"
        return out

    cache: dict[int, np.ndarray] = {}
    needed = sorted({m for pair in pairs for m in pair})
    for m in needed:
        M, st = block_matrix_with_forbidden(nodes, C, L, m, forbidden)
        cache[m] = M
        out["block_stats"][m] = {
            **st,
            "rho": spectral(M)["rho"],
            "nnz": int(np.count_nonzero(M)),
            "nilpotency_index": nilpotency_index(M),
        }

    for mB, mA in pairs:
        key = f"M{mB}_M{mA}"
        Mcomp = cache[mB] @ cache[mA]
        sp = spectral(Mcomp)
        diag = np.diag(Mcomp)
        loops = [
            {"word": nodes[i], "weight": float(diag[i])}
            for i in np.where(diag > 1e-12)[0]
        ]
        out["composites"][key] = {
            "rho": sp["rho"],
            "rho_lt_1": sp["rho"] < 1.0 - 1e-12,
            "nilpotency_index": nilpotency_index(Mcomp),
            "nnz": int(np.count_nonzero(Mcomp)),
            "self_loops": loops,
            "scc": scc_count(Mcomp),
        }
    return out


# Cuts aus der L=11-Pfad-Dekonstruktion
CUT_PRIMARY = "1110111101"  # knappe Margin ≈ 0.5
CUT_M3M2 = {"1110111101", "0111101111"}
CUT_M3M5 = {"1111010111"}
CUT_ALL = CUT_M3M2 | CUT_M3M5


def virtual_cut_surgery(C: float, out_dir: Path) -> dict:
    """
    Chirurgie: Virtual-Cut auf Forbidden-Fenstern.
    Prüft Eliminierung in V_11 und ρ(M3M2)/ρ(M3M5); Kontrolle L=8..10.
    """
    pairs = [(3, 2), (3, 5)]
    scenarios = {
        "baseline_no_cut": set(),
        "cut_primary_1110111101": {CUT_PRIMARY},
        "cut_M3M2_pair": set(CUT_M3M2),
        "cut_M3M5": set(CUT_M3M5),
        "cut_all_deconstruct": set(CUT_ALL),
    }
    print("=" * 72)
    print(f"VIRTUAL-CUT SURGERY  C={C}")
    print("=" * 72)

    payload = {
        "claim": (
            "Virtual-Cut / Forbidden-Fenster gegen L=11-Selbstschleifen; "
            "minimale Filter-Erweiterung; kein Collatz-Beweis"
        ),
        "C": C,
        "scenarios": {},
        "scale_control_L8_to_L10": {},
    }

    # --- L=11 Szenarien ---
    for name, forb in scenarios.items():
        print(f"\n### scenario={name}  forbidden={sorted(forb) or '∅'}")
        row = composite_spectrum_cut(C, 11, forb, pairs)
        payload["scenarios"][name] = row
        print(
            f"  |V| {row['n_V_raw']} → {row['n_V_cut']}  "
            f"(killed {row['n_killed']})"
        )
        if row["killed_words"] and row["n_killed"] <= 20:
            print(f"  killed: {row['killed_words']}")
        elif row["n_killed"] > 20:
            print(f"  killed (sample): {row['killed_words'][:12]} ...")
        print(f"  resonance targets hit: {row['killed_resonance_targets']}")
        for key, comp in row["composites"].items():
            print(
                f"  {key}: ρ={comp['rho']:.6f} nil={comp['nilpotency_index']} "
                f"loops={len(comp['self_loops'])} SCC_nt={comp['scc']['n_scc_nontrivial']}"
            )
            for lp in comp["self_loops"][:4]:
                print(f"    loop {lp['word']} w={lp['weight']}")

    # --- Kontrolle: bester Cut auf L=8..10 (soll Kontraktion nicht zerstören) ---
    # Wähle den kleinsten Cut, der bei L=11 beide Komposite nilpotent macht; sonst cut_all
    best_name = None
    for name in (
        "cut_primary_1110111101",
        "cut_M3M2_pair",
        "cut_M3M5",
        "cut_all_deconstruct",
    ):
        comps = payload["scenarios"][name]["composites"]
        if all(
            comps[k]["nilpotency_index"] is not None or comps[k]["rho"] < 1e-12
            for k in ("M3_M2", "M3_M5")
        ):
            best_name = name
            break
    if best_name is None:
        # fallback: kleinster Cut mit ρ<1 für beide
        for name in (
            "cut_primary_1110111101",
            "cut_M3M2_pair",
            "cut_M3M5",
            "cut_all_deconstruct",
        ):
            comps = payload["scenarios"][name]["composites"]
            if all(comps[k]["rho"] < 1 for k in ("M3_M2", "M3_M5")):
                best_name = name
                break
    control_forb = scenarios.get(best_name or "cut_all_deconstruct", CUT_ALL)
    payload["selected_control_cut"] = {
        "name": best_name or "cut_all_deconstruct",
        "forbidden": sorted(control_forb),
    }
    print(f"\n### scale control with cut={payload['selected_control_cut']['name']}")
    for L in (8, 9, 10, 11):
        base = composite_spectrum_cut(C, L, set(), pairs)
        cut = composite_spectrum_cut(C, L, set(control_forb), pairs)
        payload["scale_control_L8_to_L10"][str(L)] = {
            "baseline": {
                k: {"rho": v["rho"], "nil": v["nilpotency_index"]}
                for k, v in base["composites"].items()
            },
            "with_cut": {
                k: {"rho": v["rho"], "nil": v["nilpotency_index"]}
                for k, v in cut["composites"].items()
            },
            "n_V_raw": base["n_V_raw"],
            "n_V_cut": cut["n_V_cut"],
            "n_killed": cut["n_killed"],
        }
        print(
            f"  L={L}: |V| {base['n_V_raw']}→{cut['n_V_cut']}  "
            f"M3M2 ρ {base['composites']['M3_M2']['rho']:.4f}→"
            f"{cut['composites']['M3_M2']['rho']:.4f}  "
            f"M3M5 ρ {base['composites']['M3_M5']['rho']:.4f}→"
            f"{cut['composites']['M3_M5']['rho']:.4f}"
        )

    out = out_dir / f"h7_virtual_cut_surgery_C{str(C).replace('.','p')}.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return payload


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--mode",
        choices=("compare", "path_deconstruct", "virtual_cut"),
        default="virtual_cut",
    )
    ap.add_argument("--C", type=float, default=8.0)
    ap.add_argument("--L-ok", type=int, default=10)
    ap.add_argument("--L-break", type=int, default=11)
    ap.add_argument(
        "--out",
        type=str,
        default=str(Path("/Users/thomashoffbauer/Projects/Kepler-Hurrwitz/docs/exports")),
    )
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.mode == "compare":
        compare_break(args.C, args.L_ok, args.L_break, out_dir)
    elif args.mode == "path_deconstruct":
        deconstruct_resonance(args.C, args.L_break, out_dir)
    else:
        virtual_cut_surgery(args.C, out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


