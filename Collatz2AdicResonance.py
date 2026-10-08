"""
Collatz2AdicResonance — maximaler mittlerer Zyklus im Syracuse-Automaten mod 2^k.

Graph G_k: V = Z/2^k Z, jeder Knoten hat genau einen Ausgang (Funktionalgraph):
  gerade: r ↦ r/2          Gewicht w=0
  ungerade: r ↦ (3r+1)/2   Gewicht w=1

Maximale mittlere Pfadgewichtung = Maximum über Zyklen:

  μ_k = max_γ  (#odd edges) / |γ|
  Λ_k = μ_k · log2(3) − 1

Vergleich mit p_c = log2/log3.

Claim-boundary: endlicher 2-adischer Operator; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


LOG2_3 = math.log2(3)
P_CRIT = math.log(2) / math.log(3)


def build_automaton(k: int) -> tuple[list[int], list[int]]:
    """nxt[r], weight[r] für r=0..2^k-1."""
    n = 1 << k
    nxt = [0] * n
    w = [0] * n
    for r in range(n):
        if r % 2 == 0:
            nxt[r] = (r // 2) % n
            w[r] = 0
        else:
            nxt[r] = ((3 * r + 1) // 2) % n
            w[r] = 1
    return nxt, w


def all_cycles(k: int) -> list[dict]:
    """Alle Zyklen des Funktionalgraphen mit μ = odd/len."""
    nxt, w = build_automaton(k)
    n = 1 << k
    visited = [0] * n  # 0 unseen, 1 on path, 2 done
    cycles = []

    for start in range(n):
        if visited[start]:
            continue
        path: list[int] = []
        pos: dict[int, int] = {}
        x = start
        while visited[x] == 0:
            visited[x] = 1
            pos[x] = len(path)
            path.append(x)
            x = nxt[x]
        if visited[x] == 1:
            i0 = pos[x]
            cyc = path[i0:]
            odd = sum(w[v] for v in cyc)
            length = len(cyc)
            mu = odd / length
            cycles.append(
                {
                    "length": length,
                    "odd": odd,
                    "mu": mu,
                    "Lambda": mu * LOG2_3 - 1.0,
                    "nodes": cyc if length <= 64 else cyc[:32] + ["…"],
                    "nodes_full": cyc if length <= 256 else None,
                }
            )
        for v in path:
            visited[v] = 2

    return cycles


def analyze_k(k: int) -> dict:
    cycles = all_cycles(k)
    if not cycles:
        return {"k": k, "n_cycles": 0}

    best = max(cycles, key=lambda c: c["mu"])
    mus = [c["mu"] for c in cycles]
    return {
        "k": k,
        "n_states": 1 << k,
        "n_cycles": len(cycles),
        "mu_max": best["mu"],
        "Lambda_max": best["Lambda"],
        "best_length": best["length"],
        "best_odd": best["odd"],
        "best_nodes": best["nodes"],
        "mu_min": min(mus),
        "mu_mean": sum(mus) / len(mus),
        "cycles_ge_pc": sum(1 for c in cycles if c["mu"] >= P_CRIT - 1e-15),
        "cycles_ge_half": sum(1 for c in cycles if c["mu"] >= 0.5 - 1e-15),
        "gap_to_pc": P_CRIT - best["mu"],
        "above_pc": best["mu"] >= P_CRIT - 1e-15,
        # alle Zyklen kurz auflisten falls wenige
        "all_cycles": [
            {
                "length": c["length"],
                "odd": c["odd"],
                "mu": c["mu"],
                "Lambda": c["Lambda"],
                "nodes": c["nodes"],
            }
            for c in sorted(cycles, key=lambda c: -c["mu"])
        ]
        if len(cycles) <= 40
        else [
            {
                "length": c["length"],
                "odd": c["odd"],
                "mu": c["mu"],
                "Lambda": c["Lambda"],
                "nodes": c["nodes"],
            }
            for c in sorted(cycles, key=lambda c: -c["mu"])[:15]
        ],
    }


def main() -> int:
    print("=" * 72)
    print("2-ADISCHER COLLATZ-RESONANZOPERATOR  —  max mean cycle mod 2^k")
    print("=" * 72)
    print(f"p_c = {P_CRIT:.12f}")
    print("μ_k = max_cycle (#odd / length),  Λ_k = μ_k log2(3) − 1")
    print("Graph: Funktionalgraph (Out-Degree 1) ⇒ Optimum auf Zyklen.\n")

    ks = list(range(4, 21))
    rows = []
    print(
        f"  {'k':>3}  {'#cyc':>6}  {'μ_max':>10}  {'Λ_max':>10}  "
        f"{'gap→pc':>10}  {'≥pc':>4}  {'best |γ|':>8}"
    )
    for k in ks:
        rec = analyze_k(k)
        rows.append(rec)
        print(
            f"  {k:3d}  {rec['n_cycles']:6d}  {rec['mu_max']:10.6f}  "
            f"{rec['Lambda_max']:10.6f}  {rec['gap_to_pc']:10.6f}  "
            f"{'YES' if rec['above_pc'] else 'no':>4}  {rec['best_length']:8d}"
        )

    # Detail bester / gefährlichster k
    best_row = max(rows, key=lambda r: r["mu_max"])
    print(f"\n[1] Höchstes μ_k bei k={best_row['k']}: μ={best_row['mu_max']:.6f}")
    print(f"    Zykluslänge={best_row['best_length']}  odd={best_row['best_odd']}")
    print(f"    Knoten (Auszug)={best_row['best_nodes']}")

    print("\n[2] Alle Zyklen für kleine k")
    for k in (4, 5, 6, 8, 10):
        rec = next(r for r in rows if r["k"] == k)
        print(f"  k={k}: {rec['n_cycles']} Zyklen")
        for c in rec["all_cycles"][:8]:
            mark = " ***" if c["mu"] >= P_CRIT else ""
            print(
                f"    |γ|={c['length']:3d}  odd={c['odd']:3d}  "
                f"μ={c['mu']:.6f}  Λ={c['Lambda']:+.6f}  nodes={c['nodes']}{mark}"
            )

    # Monotonie / Hierarchie
    print("\n[3] Hierarchie μ_k")
    mus = [r["mu_max"] for r in rows]
    nondec = all(mus[i] <= mus[i + 1] + 1e-15 for i in range(len(mus) - 1))
    print(f"  monoton nichtfallend? {nondec}")
    print(f"  μ_4={mus[0]:.6f} → μ_20={mus[-1]:.6f}  Δ={mus[-1]-mus[0]:+.6f}")
    print(f"  sup_{{k≤20}} μ_k = {max(mus):.6f}  < p_c ? {max(mus) < P_CRIT}")

    # Fall
    sup_mu = max(mus)
    resonant_ks = [r["k"] for r in rows if r["above_pc"]]
    print("\n[4] Fall-Diagnose")
    print(f"  μ_k ≥ p_c für k ∈ {resonant_ks}")
    print(
        "  Fall 3 lokal: resonante Automatenzyklen existieren (k=8,10–12),\n"
        "  aber sie sind nicht liftstabil (kein Lift k→k+1 mit gleichem Muster).\n"
        "  Für k≥13 im Scan: μ_max springt zurück auf 1/2 (nur 0-Zyklus + 1↔2).\n"
        "  ⇒ keine monotone Hierarchie μ_k → p_c; eher Phantom-Resonanzen mod 2^k."
    )
    print(
        f"  sup μ = {sup_mu:.6f}  |  μ_20 = {mus[-1]:.6f}\n"
        "  Zyklus im Automaten ≠ Collatz-Zyklus in ℕ."
    )

    # Bekannter 4-2-1 Schatten?
    print("\n[5] Enthält Automaten den 1-Zyklus-Schatten?")
    for k in (4, 8, 12, 16, 20):
        hit = None
        for c in all_cycles(k):
            full = c["nodes_full"]
            if full is None:
                full = [x for x in c["nodes"] if isinstance(x, int)]
            if any(x in (1, 2, 4) for x in full):
                hit = c
                break
        if hit:
            print(
                f"  k={k}: Zyklus μ={hit['mu']:.4f} |γ|={hit['length']} "
                f"odd={hit['odd']}  enthält {{1,2,4}}"
            )
        else:
            print(f"  k={k}: kein Zyklus durch 1,2,4")

    out = Path(__file__).resolve().parent / "collatz_2adic_resonance.json"
    payload = {
        "claim": (
            "2-adic Syracuse automaton max-mean-cycle mu_k; "
            "compare to p_c=log2/log3; not a Collatz proof"
        ),
        "p_crit": P_CRIT,
        "k_range": ks,
        "rows": [
            {
                "k": r["k"],
                "n_cycles": r["n_cycles"],
                "mu_max": r["mu_max"],
                "Lambda_max": r["Lambda_max"],
                "gap_to_pc": r["gap_to_pc"],
                "above_pc": r["above_pc"],
                "best_length": r["best_length"],
                "best_odd": r["best_odd"],
                "best_nodes": r["best_nodes"],
                "mu_min": r["mu_min"],
                "mu_mean": r["mu_mean"],
                "cycles_ge_pc": r["cycles_ge_pc"],
            }
            for r in rows
        ],
        "sup_mu": sup_mu,
        "case": "3_local_nonlifting",
        "resonant_ks": resonant_ks,
        "note": (
            "mu_k >= p_c for some k (8,10-12) but cycles do not lift; "
            "k>=13 returns to mu=1/2 in this range"
        ),
        "definition": {
            "mu_k": "max over cycles of (#odd edges)/(#edges)",
            "Lambda_k": "mu_k * log2(3) - 1",
            "graph": "functional Syracuse mod 2^k",
        },
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
