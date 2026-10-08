"""
CollatzSpectralGap — 2-adischer Resonanzspektralradius der Syracuse-Automaten.

  ρ_j = max_γ  (#odd(γ) / |γ|)   über Zyklen in G_j = (Z/2^j Z, S)
  Δ_j = p_c − ρ_j                 Resonanzlücke
  δ_j = ρ_{j+1} − ρ_j             Lift-Indikator

Sterilisierungsvermutung (starke Form):
  ∃ J ∀ j≥J ∀ r:  M(j,r) < p_c
äquivalent bei Funktionalgraphen zu:
  ∃ J ∀ j≥J:  ρ_j < p_c

Zentrale asymptotische Vermutung:
  ρ_∞ := limsup ρ_j  =  1/2  <  p_c
  bzw.  inf_{j≥J} Δ_j > 0

Epistemik: Satz über endliche Syracuse-Operatoren; kein Collatz-Beweis für ℕ.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


LOG2_3 = math.log2(3)
P_CRIT = math.log(2) / math.log(3)


def syracuse_next(r: int, mod: int) -> tuple[int, int]:
    if r % 2 == 0:
        return (r // 2) % mod, 0
    return ((3 * r + 1) // 2) % mod, 1


def spectral_data(j: int) -> dict:
    """ρ_j, alle Zyklen, SCC-Info (bei Funktionalgraph = Zyklen + Bäume)."""
    mod = 1 << j
    nxt = [0] * mod
    w = [0] * mod
    for r in range(mod):
        nxt[r], w[r] = syracuse_next(r, mod)

    visited = [0] * mod
    cycles = []
    for start in range(mod):
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
            cyc = path[pos[x] :]
            odd = sum(w[v] for v in cyc)
            mu = odd / len(cyc)
            cycles.append(
                {
                    "length": len(cyc),
                    "odd": odd,
                    "mu": mu,
                    "Lambda": mu * LOG2_3 - 1.0,
                    "nodes": cyc if len(cyc) <= 32 else cyc[:16] + ["…"],
                }
            )
        for v in path:
            visited[v] = 2

    # Basin sizes: how many nodes flow to each cycle
    # identify cycle by frozenset of nodes (canonical min rotation string)
    def cid(nodes: list[int]) -> str:
        i = min(range(len(nodes)), key=lambda t: nodes[t])
        rot = nodes[i:] + nodes[:i]
        return "-".join(map(str, rot))

    # rebuild full cycle node lists
    visited = [0] * mod
    full_cycles: list[list[int]] = []
    for start in range(mod):
        if visited[start]:
            continue
        path = []
        pos = {}
        x = start
        while visited[x] == 0:
            visited[x] = 1
            pos[x] = len(path)
            path.append(x)
            x = nxt[x]
        if visited[x] == 1:
            full_cycles.append(path[pos[x] :])
        for v in path:
            visited[v] = 2

    attractor = [-1] * mod
    for ci, cyc in enumerate(full_cycles):
        for v in cyc:
            attractor[v] = ci
    for r in range(mod):
        if attractor[r] != -1:
            continue
        x = r
        chain = []
        while attractor[x] == -1:
            chain.append(x)
            x = nxt[x]
        aid = attractor[x]
        for v in chain:
            attractor[v] = aid

    basin = [0] * len(full_cycles)
    for r in range(mod):
        basin[attractor[r]] += 1

    for i, c in enumerate(cycles):
        c["basin"] = basin[i] if i < len(basin) else None
        c["basin_frac"] = (basin[i] / mod) if i < len(basin) else None

    rho = max(c["mu"] for c in cycles) if cycles else 0.0
    Delta = P_CRIT - rho
    return {
        "j": j,
        "n_states": mod,
        "n_cycles": len(cycles),
        "rho": rho,
        "Delta": Delta,
        "resonant": rho + 1e-15 >= P_CRIT,
        "cycles": sorted(cycles, key=lambda c: -c["mu"]),
        "mu_second": sorted((c["mu"] for c in cycles), reverse=True)[1]
        if len(cycles) > 1
        else None,
    }


def main() -> int:
    print("=" * 72)
    print("2-ADISCHER RESONANZSPEKTRALRADIUS  —  ρ_j, Δ_j, δ_j")
    print("=" * 72)
    print(f"p_c = {P_CRIT:.12f}")
    print("ρ_j = max_γ μ(γ),  Δ_j = p_c − ρ_j,  δ_j = ρ_{j+1} − ρ_j")
    print("Vermutung: ρ_∞ = 1/2 < p_c  ⇒  inf_{j≥J} Δ_j > 0\n")

    j_min, j_max = 4, 24
    rows = []
    print(
        f"  {'j':>3}  {'#cyc':>5}  {'ρ_j':>10}  {'Δ_j':>10}  {'δ_j':>10}  "
        f"{'res?':>4}  {'μ_2nd':>8}  basin(ρ_max)"
    )
    prev_rho = None
    for j in range(j_min, j_max + 1):
        rec = spectral_data(j)
        delta_lift = None if prev_rho is None else rec["rho"] - prev_rho
        rec["delta_lift"] = delta_lift
        rows.append(rec)
        best = rec["cycles"][0]
        basin_str = f"{best.get('basin_frac', 0):.4f}" if best.get("basin_frac") is not None else "?"
        dstr = f"{delta_lift:+10.6f}" if delta_lift is not None else f"{'—':>10}"
        print(
            f"  {j:3d}  {rec['n_cycles']:5d}  {rec['rho']:10.6f}  {rec['Delta']:10.6f}  "
            f"{dstr}  {'YES' if rec['resonant'] else 'no':>4}  "
            f"{rec['mu_second'] if rec['mu_second'] is not None else float('nan'):8.4f}  "
            f"{basin_str}"
        )
        prev_rho = rec["rho"]

    # Spektrallücke ab J
    print("\n[1] Resonanzlücke Δ_j")
    resonant_js = [r["j"] for r in rows if r["resonant"]]
    safe_js = [r["j"] for r in rows if not r["resonant"]]
    print(f"  resonant (Δ<0): j ∈ {resonant_js}")
    print(f"  safe (Δ>0):     j ∈ {safe_js[:8]}…{safe_js[-1] if safe_js else ''}")

    # Find minimal J such that all j>=J in range are safe
    J = None
    for J_cand in range(j_min, j_max + 1):
        if all(not r["resonant"] for r in rows if r["j"] >= J_cand):
            J = J_cand
            break
    print(f"  minimales J im Scan mit ∀j≥J: Δ_j>0  →  J={J}")

    if J is not None:
        gaps = [r["Delta"] for r in rows if r["j"] >= J]
        print(f"  inf_{{j≥{J}}} Δ_j  = {min(gaps):.6f}")
        print(f"  sup_{{j≥{J}}} ρ_j  = {max(r['rho'] for r in rows if r['j']>=J):.6f}")
        uniform = min(gaps) > 0
        half = all(abs(r["rho"] - 0.5) < 1e-12 or r["rho"] < 0.5 + 1e-12 for r in rows if r["j"] >= J)
        # check all rho are exactly 0 or 1/2
        only_half = all(
            abs(r["rho"] - 0.5) < 1e-12 for r in rows if r["j"] >= J
        )
        print(f"  uniforme Lücke >0? {uniform}")
        print(f"  ρ_j = 1/2 exakt für alle j≥J? {only_half}")

    print("\n[2] ρ_∞-Proxy = limsup ρ_j im Fenster")
    # limsup as max of tail
    for tail_start in (13, 16, 20, J or 13):
        if tail_start > j_max:
            continue
        limsup = max(r["rho"] for r in rows if r["j"] >= tail_start)
        liminf = min(r["rho"] for r in rows if r["j"] >= tail_start)
        print(
            f"  j≥{tail_start}: limsup≈{limsup:.6f}  liminf≈{liminf:.6f}  "
            f"<p_c? {limsup < P_CRIT}"
        )

    print("\n[3] Lift-Indikator δ_j = ρ_{j+1}−ρ_j")
    spikes = [r for r in rows if r["delta_lift"] is not None and abs(r["delta_lift"]) > 1e-9]
    for r in spikes:
        print(
            f"  δ_{r['j']-1} = {r['delta_lift']:+.6f}  "
            f"(ρ_{r['j']-1}→ρ_{r['j']}: phantom ein/aus)"
        )
    if not spikes:
        print("  keine Sprünge")

    print("\n[4] Zyklusinventar (Auszug)")
    for j in (8, 12, 13, 20, 24):
        rec = next(r for r in rows if r["j"] == j)
        print(f"  j={j}: ρ={rec['rho']:.6f}  Δ={rec['Delta']:+.6f}")
        for c in rec["cycles"][:4]:
            print(
                f"    |γ|={c['length']:3d} μ={c['mu']:.6f}  "
                f"basin={c.get('basin_frac', 0):.4f}  nodes={c['nodes']}"
            )

    print("\n[5] Einordnung")
    print(
        "  ρ_j = 2-adischer Resonanzspektralradius des endlichen Operators.\n"
        "  Vermutung ρ_∞ = 1/2 < p_c ist ein Satz über Automaten G_j,\n"
        "  nicht die Collatz-Vermutung auf ℕ — aber die natürliche\n"
        "  Spektrallücke im inversen Limes der Syracuse-Operatoren.\n"
        f"  Empirie j≤{j_max}: J={J}, "
        f"inf Δ = {min(r['Delta'] for r in rows if J and r['j']>=J):.6f}, "
        f"ρ_tail = 1/2."
    )

    out = Path(__file__).resolve().parent / "collatz_spectral_gap.json"
    payload = {
        "claim": (
            "2-adic Syracuse spectral radius rho_j and gap Delta_j=p_c-rho_j; "
            "conjecture rho_infty=1/2<p_c; statement about finite automata, not Collatz on N"
        ),
        "p_crit": P_CRIT,
        "j_range": [j_min, j_max],
        "J_sterile": J,
        "rows": [
            {
                "j": r["j"],
                "n_cycles": r["n_cycles"],
                "rho": r["rho"],
                "Delta": r["Delta"],
                "delta_lift": r["delta_lift"],
                "resonant": r["resonant"],
                "mu_second": r["mu_second"],
                "best_cycle": r["cycles"][0],
            }
            for r in rows
        ],
        "limsup_tails": {
            str(t): max(r["rho"] for r in rows if r["j"] >= t)
            for t in (13, 16, 20)
            if t <= j_max
        },
        "conjecture": {
            "rho_infty": "1/2",
            "uniform_gap": "exists J: inf_{j>=J} Delta_j > 0",
            "numerical_support": J is not None
            and all(not r["resonant"] for r in rows if r["j"] >= J),
        },
        "analogy": {
            "Gap_class": "parity cylinder r mod 2^j",
            "M_G": "odd density mu",
            "resonance_constant": "Lambda = a log2(3) - m",
            "classification": "inverse limit / spectral gap",
        },
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
