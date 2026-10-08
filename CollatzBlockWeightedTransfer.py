"""
CollatzBlockWeightedTransfer — Blockoperator M_{w,m} vs M_w^m.

Kernpunkt:
  M_{w,m} ≠ M_w^m, weil Pfade der Länge m zusätzlich gegen die
  Survivor-Regel auf Wörtern der Länge L+m gefiltert werden.
  Diese Asymmetrie ist die mögliche Quelle spektraler Kontraktion.

Setup (Default): C=8, L=8, m∈{2,3,4}
  V_L  = Clog-Survivor-Wörter Länge L (fix, lex)
  Kante u→v auf V_L: Shift u[1:]==v[:-1], Gewicht w(v[-1])∈{1/2,3/2}

Blockvariante B (Graph-Walk + Tiefenfilter):
  Summe über Walks u=v0→…→vm=v in V_L der Länge m,
  deren Bitfolge path die Erweiterung word(u)+path als
  Clog-Survivor der Länge L+m besteht;
  Beitrag = Produkt der Kanten-Gewichte.

Export: JSON + Spektren + Vergleich zu M_w^m.
Epistemik: Diagnostik; θ<1 nicht garantiert; kein Collatz-Beweis.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

WBIT = {"0": 0.5, "1": 1.5}


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


def enumerate_survivor_words(L: int, C: float) -> list[str]:
    best: dict[str, int] = {}
    for r in range(1 << L):
        if not is_survivor_clog(r, L, C):
            continue
        w = word_of_residue(r, L)
        if w not in best or r < best[w]:
            best[w] = r
    return sorted(best)


def survivor_word_set(L: int, C: float) -> set[str]:
    return set(enumerate_survivor_words(L, C))


def build_base_Mw(nodes: list[str]) -> tuple[np.ndarray, list[tuple[int, int, str, float]]]:
    idx = {w: i for i, w in enumerate(nodes)}
    n = len(nodes)
    M = np.zeros((n, n), dtype=float)
    edges = []
    for u in nodes:
        for b, wgt in WBIT.items():
            v = u[1:] + b
            if v in idx:
                i, j = idx[u], idx[v]
                M[i, j] = wgt
                edges.append((i, j, b, wgt))
    return M, edges


def adjacency_list(nodes: list[str]) -> list[list[tuple[int, str, float]]]:
    idx = {w: i for i, w in enumerate(nodes)}
    adj: list[list[tuple[int, str, float]]] = [[] for _ in nodes]
    for i, u in enumerate(nodes):
        for b, wgt in WBIT.items():
            v = u[1:] + b
            if v in idx:
                adj[i].append((idx[v], b, wgt))
    return adj


def block_matrix_walk_filter(
    nodes: list[str], C: float, L: int, m: int
) -> tuple[np.ndarray, dict]:
    """
    Walks on V_L of length m; keep only if u+path is survivor at L+m;
    endpoint automatically in V_L.
    """
    n = len(nodes)
    adj = adjacency_list(nodes)
    M = np.zeros((n, n), dtype=float)
    surv_Lm = survivor_word_set(L + m, C)
    stats = {
        "walks_total": 0,
        "walks_kept": 0,
        "walks_rejected_filter": 0,
        "n_survivors_Lpm": len(surv_Lm),
    }

    def rec(i0: int, i: int, depth: int, path: str, wprod: float) -> None:
        if depth == m:
            stats["walks_total"] += 1
            full = nodes[i0] + path
            assert len(full) == L + m
            if full in surv_Lm:
                stats["walks_kept"] += 1
                M[i0, i] += wprod
            else:
                stats["walks_rejected_filter"] += 1
            return
        for j, b, wgt in adj[i]:
            rec(i0, j, depth + 1, path + b, wprod * wgt)

    for i0 in range(n):
        rec(i0, i0, 0, "", 1.0)

    return M, stats


def spectral(M: np.ndarray) -> dict:
    if M.size == 0:
        return {"rho": 0.0, "n": 0}
    vals, vecs = np.linalg.eig(M)
    k = int(np.argmax(np.abs(vals)))
    rho = float(np.abs(vals[k]))
    v = np.real(vecs[:, k])
    if v.min() < 0:
        v = -v
    v = np.maximum(v, 0.0)
    if v.max() > 0:
        v = v / v.max()
    return {
        "rho": rho,
        "n": int(M.shape[0]),
        "row_sum_max": float(M.sum(axis=1).max()) if M.shape[0] else 0.0,
        "nnz": int(np.count_nonzero(M)),
        "v_support": int((v > 1e-12).sum()),
        "v": v.tolist(),
    }


def find_subinvariant(M: np.ndarray, v_hint: np.ndarray) -> dict | None:
    supp = v_hint > 1e-10
    if supp.sum() == 0:
        return None
    # try on full with padded positives
    best = None
    x = v_hint.copy()
    x = np.maximum(x, 0)
    if x.max() > 0:
        x = x / x.max()
    for scale in (50, 100, 200, 500, 1000, 2000, 5000):
        vi = np.maximum(np.rint(x * scale), 0).astype(int)
        # restrict search to support; set zeros to 0 and only certify on support indices
        if (vi > 0).sum() == 0:
            continue
        Mv = M @ vi.astype(float)
        # only require inequality where vi>0; where vi=0 require Mv≈0 or skip
        pos = vi > 0
        if not np.any(pos):
            continue
        th = float((Mv[pos] / vi[pos]).max())
        for den in (10, 20, 50, 100, 200, 500, 1000):
            num = int(np.ceil(th * den - 1e-15))
            theta = num / den
            ok_pos = np.all(Mv[pos] <= theta * vi[pos] + 1e-8)
            ok_zero = np.all(Mv[~pos] <= 1e-8)
            if ok_pos and ok_zero:
                if best is None or theta < best["theta"]:
                    best = {
                        "theta_num": num,
                        "theta_den": den,
                        "theta": theta,
                        "theta_lt_1": theta < 1,
                        "v": vi.tolist(),
                        "scale": scale,
                    }
    return best


def nilpotency_index(M: np.ndarray, max_k: int = 40) -> int | None:
    P = np.eye(M.shape[0])
    for k in range(1, max_k + 1):
        P = P @ M
        if np.allclose(P, 0):
            return k
    return None


def power_iteration(M: np.ndarray, iters: int = 200) -> tuple[float, np.ndarray]:
    n = M.shape[0]
    v = np.ones(n)
    for _ in range(iters):
        w = M @ v
        nrm = np.linalg.norm(w, 1)
        if nrm < 1e-30:
            return 0.0, v
        v = w / nrm
    w = M @ v
    # Rayleigh-like for nonnegative: max component ratio on support
    pos = v > 1e-15
    if not np.any(pos):
        return 0.0, v
    rho_est = float(np.max(w[pos] / v[pos]))
    if v.max() > 0:
        v = v / v.max()
    return rho_est, v


def rational_subinvariant_from_perron(
    M: np.ndarray, v_float: np.ndarray, rho: float, margin: float = 1.01
) -> dict | None:
    """Integer v≈Perron, θ ≥ margin·ρ with Mv ≤ θv (upper certificate)."""
    if rho < 1e-14:
        return {
            "regime": "nilpotent",
            "rho": 0.0,
            "theta_num": 1,
            "theta_den": 2,
            "theta": 0.5,
            "theta_lt_1": True,
            "note": "use constructive nilpotent certificate separately (m even / some odd)",
        }
    x = np.maximum(v_float, 0.0)
    if x.max() <= 0:
        return None
    x = x / x.max()
    best = None
    target = rho * margin
    for scale in (50, 100, 200, 500, 1000, 2000, 5000, 10000):
        vi = np.maximum(np.rint(x * scale), 0).astype(int)
        pos = vi > 0
        if pos.sum() == 0:
            continue
        # lift exact zeros of evec slightly for strict positivity attempts
        vi2 = vi.copy()
        vi2[vi2 == 0] = 1
        for cand_v, tag in ((vi, "support"), (vi2, "strict_pos")):
            Mv = M @ cand_v.astype(float)
            p = cand_v > 0
            th_need = float((Mv[p] / cand_v[p]).max())
            th_need = max(th_need, target)
            for den in (10, 20, 50, 100, 200, 500, 1000, 2000):
                num = int(np.ceil(th_need * den - 1e-15))
                theta = num / den
                if np.all(Mv[p] <= theta * cand_v[p] + 1e-8) and np.all(
                    Mv[~p] <= 1e-8
                ):
                    rec = {
                        "regime": "recurrent",
                        "v_mode": tag,
                        "rho": rho,
                        "theta_num": num,
                        "theta_den": den,
                        "theta": theta,
                        "theta_lt_1": theta < 1,
                        "gap_rho_minus_1": rho - 1.0,
                        "gap_theta_minus_1": theta - 1.0,
                        "v_int": cand_v.tolist(),
                        "v_support": int((cand_v > 0).sum()),
                        "scale": scale,
                    }
                    if best is None or theta < best["theta"]:
                        best = rec
    return best


def analyze_odd_m(C: float, L: int, odd_ms: list[int], out_dir: Path) -> dict:
    """
    Perron-Konvergenz-Analyse für ungerade Blocklängen.
    Speichert h7_perron_chain_L{L}_C….json (ohne volle Matrizen).
    """
    nodes = enumerate_survivor_words(L, C)
    Mw, _ = build_base_Mw(nodes)
    print("=" * 72)
    print(f"PERRON-CHAIN (odd m)  C={C} L={L} m={odd_ms}")
    print(f"|V_L|={len(nodes)}")
    print("=" * 72)

    chain = []
    print(
        f"{'m':>3} {'ρ':>10} {'Δ=ρ-1':>10} {'ρ_pow':>10} {'supp':>5} "
        f"{'nil?':>5} {'θ_up':>10} {'θ-1':>10} {'kept':>6}"
    )
    for m in odd_ms:
        assert m % 2 == 1
        Mb, stats = block_matrix_walk_filter(nodes, C, L, m)
        sp = spectral(Mb)
        rho_eig = sp["rho"]
        rho_pow_it, v_pow = power_iteration(Mb)
        rho = max(rho_eig, rho_pow_it) if rho_eig > 1e-14 else rho_eig
        # prefer eig if recurrent
        if rho_eig > 1e-14:
            rho = rho_eig
            v = np.array(sp["v"], float)
        else:
            rho = 0.0
            v = v_pow
        nil_k = nilpotency_index(Mb)
        Mw_pow = np.linalg.matrix_power(Mw, m)
        rho_power = float(np.max(np.abs(np.linalg.eigvals(Mw_pow))))
        cert = rational_subinvariant_from_perron(Mb, v, rho, margin=1.01)
        delta = abs(rho - 1.0)
        theta_s = f"{cert['theta']:.6f}" if cert and "theta" in cert else "—"
        gap_t = cert["gap_theta_minus_1"] if cert and cert.get("regime") == "recurrent" else float("nan")
        print(
            f"{m:3d} {rho:10.6f} {rho-1:10.6f} {rho_power:10.6f} "
            f"{int((v>1e-12).sum()):5d} {str(nil_k):>5} {theta_s:>10} "
            f"{gap_t:10.6f} {stats['walks_kept']:6d}"
        )
        chain.append(
            {
                "m": m,
                "rho": rho,
                "delta_abs_rho_minus_1": delta,
                "gap_rho_minus_1": rho - 1.0,
                "rho_power_Mw_m": rho_power,
                "rho_block_over_power": (rho / rho_power) if rho_power > 1e-15 else None,
                "v_support": int((v > 1e-12).sum()),
                "nnz": int(np.count_nonzero(Mb)),
                "nilpotency_index": nil_k,
                "is_nilpotent": nil_k is not None,
                "rho_lt_1": rho < 1.0 - 1e-12,
                "stats": stats,
                "v_perron_float": v.tolist(),
                "subinvariant_upper": cert,
            }
        )

    # monotonicity among recurrent (rho>1) entries
    recurrent = [c for c in chain if c["rho"] > 1 + 1e-12]
    mono = None
    if len(recurrent) >= 2:
        rhos = [c["rho"] for c in recurrent]
        mono = all(rhos[i] >= rhos[i + 1] - 1e-9 for i in range(len(rhos) - 1))

    payload = {
        "claim": (
            "Odd-m Perron chain for block operator M_w,m; "
            "tests convergence of rho toward 1 vs resonance barrier; "
            "odd m are NOT always non-nilpotent; not a Collatz proof"
        ),
        "C": C,
        "L": L,
        "n": len(nodes),
        "nodes": nodes,
        "odd_m": odd_ms,
        "chain": chain,
        "analysis": {
            "recurrent_count": len(recurrent),
            "nilpotent_odd_m": [c["m"] for c in chain if c["is_nilpotent"]],
            "recurrent_odd_m": [c["m"] for c in chain if c["rho"] > 1 + 1e-12],
            "rho_monotone_nonincreasing_on_recurrent": mono,
            "suggests_limit_to_1": (
                mono is True
                and len(recurrent) >= 2
                and recurrent[-1]["rho"] < recurrent[0]["rho"]
            ),
            "suggests_resonance_barrier_above_1": (
                len(recurrent) >= 1 and all(c["rho"] > 1 + 1e-6 for c in recurrent)
            ),
        },
    }
    out = out_dir / f"h7_perron_chain_L{L}_C{str(C).replace('.', 'p')}.json"
    # also alias requested name
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    alias = out_dir / "h7_perron_chain_L8_C8.json"
    if out != alias:
        alias.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    if alias != out:
        print(f"wrote {alias}")
    print("analysis:", json.dumps(payload["analysis"], indent=2))
    return payload


def get_block_matrix(
    nodes: list[str], C: float, L: int, m: int, cache: dict[int, np.ndarray]
) -> np.ndarray:
    if m not in cache:
        if m == 1:
            cache[m], _ = build_base_Mw(nodes)
        else:
            cache[m], _ = block_matrix_walk_filter(nodes, C, L, m)
    return cache[m]


def analyze_composite(
    C: float,
    L: int,
    mA_list: list[int],
    mB_list: list[int],
    out_dir: Path,
) -> dict:
    """
    Komposition Klasse-B (rekurrent-Kandidat) × Klasse-A (Dämpfer/nilpotent).

    Misst ρ(M_B M_A) und ρ(M_A M_B). Caveat: ρ(A)=0 ⇏ ρ(BA)<1.
    """
    nodes = enumerate_survivor_words(L, C)
    cache: dict[int, np.ndarray] = {}
    print("=" * 72)
    print(f"COMPOSITE FILTER  C={C} L={L}")
    print(f"  m_A (damper)={mA_list}  m_B (resonant)={mB_list}")
    print(f"|V_L|={len(nodes)}")
    print("=" * 72)

    # classify singles
    singles = {}
    for m in sorted(set(mA_list + mB_list)):
        M = get_block_matrix(nodes, C, L, m, cache)
        sp = spectral(M)
        nil_k = nilpotency_index(M)
        singles[m] = {
            "rho": sp["rho"],
            "nilpotency_index": nil_k,
            "class": "A_damper" if (nil_k is not None or sp["rho"] < 1e-12) else "B_recurrent",
            "nnz": sp["nnz"],
        }
        print(
            f"  single m={m}: ρ={sp['rho']:.6f}  nil={nil_k}  → {singles[m]['class']}"
        )

    rows = []
    print(
        f"\n{'mB':>3} {'mA':>3} {'order':>8} {'ρ':>10} {'ρ<1':>5} "
        f"{'nnz':>5} {'nil?':>5} {'θ_up':>10}"
    )
    for mB in mB_list:
        for mA in mA_list:
            MB = get_block_matrix(nodes, C, L, mB, cache)
            MA = get_block_matrix(nodes, C, L, mA, cache)
            for order, Mcomp in (
                ("B*A", MB @ MA),
                ("A*B", MA @ MB),
            ):
                sp = spectral(Mcomp)
                nil_k = nilpotency_index(Mcomp)
                cert = rational_subinvariant_from_perron(
                    Mcomp, np.array(sp["v"], float), sp["rho"], margin=1.01
                )
                # constructive nilpotent cert marker
                contractive = sp["rho"] < 1.0 - 1e-12
                if contractive and nil_k is not None:
                    # α=2 ⇒ θ=1/2 constructive exists for nilpotent nonnegative
                    theta_up = 0.5
                else:
                    theta_up = cert["theta"] if cert and "theta" in cert else None
                print(
                    f"{mB:3d} {mA:3d} {order:>8} {sp['rho']:10.6f} "
                    f"{str(contractive):>5} {sp['nnz']:5d} {str(nil_k):>5} "
                    f"{(f'{theta_up:.4f}' if theta_up is not None else '—'):>10}"
                )
                rows.append(
                    {
                        "mB": mB,
                        "mA": mA,
                        "order": order,
                        "formula": f"M_{mB} @ M_{mA}" if order == "B*A" else f"M_{mA} @ M_{mB}",
                        "rho": sp["rho"],
                        "rho_lt_1": contractive,
                        "nnz": sp["nnz"],
                        "nilpotency_index": nil_k,
                        "single_A": singles.get(mA),
                        "single_B": singles.get(mB),
                        "subinvariant_upper": cert,
                        "theta_contractive_candidate": theta_up if contractive else None,
                    }
                )

    winners = [r for r in rows if r["rho_lt_1"]]
    payload = {
        "claim": (
            "Composite M_B M_A / M_A M_B block products; "
            "tests damper×resonance hypothesis; "
            "rho(A)=0 does not imply rho(BA)<1 a priori; not Collatz proof"
        ),
        "C": C,
        "L": L,
        "n": len(nodes),
        "mA_list": mA_list,
        "mB_list": mB_list,
        "singles": {str(k): v for k, v in singles.items()},
        "composites": rows,
        "summary": {
            "n_tested": len(rows),
            "n_rho_lt_1": len(winners),
            "first_rho_lt_1": winners[0] if winners else None,
            "all_winners": [
                {"mB": w["mB"], "mA": w["mA"], "order": w["order"], "rho": w["rho"]}
                for w in winners
            ],
        },
    }
    out = out_dir / f"h7_composite_filter_L{L}_C{str(C).replace('.', 'p')}.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    alias = out_dir / "h7_composite_filter_L8_C8.json"
    alias.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    print("winners (ρ<1):", payload["summary"]["all_winners"])
    return payload


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--C", type=float, default=8.0)
    ap.add_argument("--L", type=int, default=8)
    ap.add_argument("--m", type=str, default="1,2,3,4")
    ap.add_argument(
        "--mode",
        type=str,
        default="block",
        choices=["block", "analyze_odd_m", "composite_filter"],
        help="block | analyze_odd_m | composite_filter",
    )
    ap.add_argument(
        "--odd-m",
        type=str,
        default="3,5,7,9",
        help="odd m list for analyze_odd_m (default 3,5,7,9)",
    )
    ap.add_argument(
        "--mA",
        type=str,
        default="2,4,5,7",
        help="damper class-A block lengths for composite_filter",
    )
    ap.add_argument(
        "--mB",
        type=str,
        default="3,9",
        help="resonant class-B block lengths for composite_filter",
    )
    ap.add_argument(
        "--out",
        type=str,
        default=str(Path("/Users/thomashoffbauer/Projects/Kepler-Hurrwitz/docs/exports")),
    )
    args = ap.parse_args()
    C, L = args.C, args.L
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.mode == "analyze_odd_m":
        odd_ms = [int(x) for x in args.odd_m.split(",") if x.strip()]
        for m in odd_ms:
            if m % 2 == 0:
                raise SystemExit(f"analyze_odd_m expects odd m, got {m}")
        analyze_odd_m(C, L, odd_ms, out_dir)
        return 0

    if args.mode == "composite_filter":
        mA_list = [int(x) for x in args.mA.split(",") if x.strip()]
        mB_list = [int(x) for x in args.mB.split(",") if x.strip()]
        analyze_composite(C, L, mA_list, mB_list, out_dir)
        return 0

    ms = [int(x) for x in args.m.split(",") if x.strip()]

    print("=" * 72)
    print(f"BLOCK WEIGHTED TRANSFER  C={C} L={L} m={ms}")
    print("  M_{w,m} = Walks on V_L + Survivor-Filter auf Wörtern Länge L+m")
    print("  Asymmetrie zu M_w^m ist die erwartete Kontraktionsquelle")
    print("=" * 72)

    nodes = enumerate_survivor_words(L, C)
    print(f"|V_L| = {len(nodes)}")
    if not nodes:
        print("V_L leer — Abbruch")
        return 2

    Mw, edges = build_base_Mw(nodes)
    spec1 = spectral(Mw)
    print(f"m=1 (base): ρ(M_w)={spec1['rho']:.6f}  nnz={spec1['nnz']}  edges={len(edges)}")

    payload = {
        "claim": (
            "Block operator M_w,m with walk-on-V_L plus L+m survivor filter; "
            "asymmetry vs M_w^m is contraction candidate; not Collatz proof"
        ),
        "C": C,
        "L": L,
        "nodes": nodes,
        "n": len(nodes),
        "base": {"rho": spec1["rho"], "nnz": spec1["nnz"], "M": Mw.tolist()},
        "blocks": {},
    }

    print(f"\n{'m':>3} {'ρ(M_w,m)':>12} {'ρ(M_w^m)':>12} {'ratio':>10} {'nnz':>6} {'kept':>8} {'rej':>8} {'θ?':>8}")
    for m in ms:
        if m == 1:
            Mb = Mw.copy()
            stats = {"walks_total": len(edges), "walks_kept": len(edges), "walks_rejected_filter": 0}
        else:
            Mb, stats = block_matrix_walk_filter(nodes, C, L, m)
        sp = spectral(Mb)
        Mw_pow = np.linalg.matrix_power(Mw, m)
        rho_pow = float(np.max(np.abs(np.linalg.eigvals(Mw_pow))))
        # asymmetry measure
        asym = float(np.max(np.abs(Mb - Mw_pow)))
        ratio = sp["rho"] / rho_pow if rho_pow > 1e-15 else float("inf")
        cert = find_subinvariant(Mb, np.array(sp["v"], float))
        theta_s = f"{cert['theta']:.4f}" if cert else "—"
        print(
            f"{m:3d} {sp['rho']:12.6f} {rho_pow:12.6f} {ratio:10.4f} "
            f"{sp['nnz']:6d} {stats['walks_kept']:8d} {stats['walks_rejected_filter']:8d} {theta_s:>8}"
        )
        payload["blocks"][str(m)] = {
            "m": m,
            "rho_block": sp["rho"],
            "rho_power": rho_pow,
            "rho_block_over_rho_power": ratio,
            "asym_maxabs_Mb_minus_Mwpow": asym,
            "nnz": sp["nnz"],
            "stats": stats,
            "rho_lt_1": sp["rho"] < 1.0 - 1e-12,
            "subinvariant": cert,
            "M": Mb.tolist(),
            "v_perron_float": sp["v"],
        }

    # summary
    rhos = {int(k): v["rho_block"] for k, v in payload["blocks"].items()}
    payload["summary"] = {
        "rho_by_m": rhos,
        "first_m_with_rho_lt_1": next((m for m in ms if rhos.get(m, 2) < 1), None),
        "note": (
            "If rho(M_w,m)<1 for some m, block spectral gap holds for this C,L. "
            "If all rho>1, need larger m or different filter/weight."
        ),
    }

    out = out_dir / f"h7_block_transfer_L{L}_C{str(C).replace('.', 'p')}.json"
    # omit huge matrices from pretty print? keep them for Lean export
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    print("first m with ρ<1:", payload["summary"]["first_m_with_rho_lt_1"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
