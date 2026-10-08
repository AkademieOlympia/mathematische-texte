"""
CollatzLiftResonance — Liftstabilität resonanter 2-adischer Zyklen.

Lokale Resonanz μ(γ_k) ≥ p_c ist nötig, aber nicht hinreichend.
Entscheidend: Liftstabilität entlang π_{j,k}: Z/2^j → Z/2^k.

  L(γ_k) = max { j ≥ k : γ_k besitzt Lift bis mod 2^j }
  R_k    = max { L(γ) : γ resonant bei Level k }

Vermutung: Der Resonanzbaum hat keine unendlichen Äste
  (jedes resonante γ_k hat endliche Lift-Tiefe).

Claim-boundary: Inverse-Limit-Diagnostik; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path


LOG2_3 = math.log2(3)
P_CRIT = math.log(2) / math.log(3)


def syracuse_next(r: int, mod: int) -> tuple[int, int]:
    """(next, odd_weight) mod `mod`."""
    if r % 2 == 0:
        return (r // 2) % mod, 0
    return ((3 * r + 1) // 2) % mod, 1


def build_automaton(k: int) -> tuple[list[int], list[int]]:
    n = 1 << k
    nxt = [0] * n
    w = [0] * n
    for r in range(n):
        nxt[r], w[r] = syracuse_next(r, n)
    return nxt, w


@dataclass(frozen=True)
class Cycle:
    k: int
    nodes: tuple[int, ...]
    mu: float
    odd: int

    @property
    def length(self) -> int:
        return len(self.nodes)

    @property
    def Lambda(self) -> float:
        return self.mu * LOG2_3 - 1.0

    @property
    def resonant(self) -> bool:
        return self.mu + 1e-15 >= P_CRIT

    @property
    def cid(self) -> str:
        # kanonische Rotation
        ns = self.nodes
        i = min(range(len(ns)), key=lambda t: ns[t])
        rot = ns[i:] + ns[:i]
        return f"k{self.k}:" + "-".join(str(x) for x in rot)


def all_cycles(k: int) -> list[Cycle]:
    nxt, w = build_automaton(k)
    n = 1 << k
    visited = [0] * n
    out: list[Cycle] = []
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
            cyc = path[pos[x] :]
            odd = sum(w[v] for v in cyc)
            mu = odd / len(cyc)
            out.append(Cycle(k=k, nodes=tuple(cyc), mu=mu, odd=odd))
        for v in path:
            visited[v] = 2
    return out


def try_lift_one_level(cyc: Cycle) -> list[Cycle]:
    """
    Alle Lifts von γ_k nach mod 2^{k+1}.
    Startlift s ≡ nodes[0] (mod 2^k), zwei Kandidaten; Orbit muss
    projizieren und schließen.
    """
    k = cyc.k
    mod_lo = 1 << k
    mod_hi = 1 << (k + 1)
    base = cyc.nodes[0]
    L = cyc.length
    lifts: list[Cycle] = []

    for s0 in (base, base + mod_lo):
        if s0 >= mod_hi:
            continue
        nodes: list[int] = []
        odd = 0
        x = s0
        ok = True
        for expected in cyc.nodes:
            if x % mod_lo != expected:
                ok = False
                break
            nodes.append(x)
            x2, bit = syracuse_next(x, mod_hi)
            odd += bit
            x = x2
        if not ok:
            continue
        if x != s0:
            continue
        # Zusatz: gesamte Projektion muss exakt γ sein (Reihenfolge)
        if tuple(n % mod_lo for n in nodes) != cyc.nodes:
            continue
        mu = odd / L
        lifts.append(Cycle(k=k + 1, nodes=tuple(nodes), mu=mu, odd=odd))
    return lifts


def lift_tree(
    cyc: Cycle,
    j_max: int,
    require_resonant_lifts: bool = False,
) -> dict:
    """
    BFS/DFS der Lift-Nachfolger bis j_max.
    L(γ) = max erreichte Ebene (oder k-1 falls kein Lift).
    Wenn require_resonant_lifts: nur Lifts mit μ≥p_c zählen für R_k-Pfad.
    """
    # depth = j - k (Anzahl erfolgreicher Lift-Schritte)
    frontier = [cyc]
    levels: dict[int, list[Cycle]] = {cyc.k: [cyc]}
    max_j = cyc.k
    total_nodes = 1
    dead_ends = 0

    for j in range(cyc.k, j_max):
        nxt_level: list[Cycle] = []
        for g in frontier:
            lifts = try_lift_one_level(g)
            if not lifts:
                dead_ends += 1
                continue
            for h in lifts:
                if require_resonant_lifts and not h.resonant:
                    continue
                nxt_level.append(h)
        if not nxt_level:
            break
        # Dedup by cid
        uniq: dict[str, Cycle] = {}
        for h in nxt_level:
            uniq[h.cid] = h
        frontier = list(uniq.values())
        levels[j + 1] = frontier
        max_j = j + 1
        total_nodes += len(frontier)

    return {
        "root": _cyc_pub(cyc),
        "L": max_j,  # max modulus exponent reached
        "depth": max_j - cyc.k,  # number of successful lifts
        "j_max_attempted": j_max,
        "total_lift_nodes": total_nodes,
        "dead_ends_encountered": dead_ends,
        "levels_summary": {
            str(j): [
                {"cid": g.cid, "mu": g.mu, "Lambda": g.Lambda, "length": g.length}
                for g in gs
            ]
            for j, gs in sorted(levels.items())
        },
        "survives_to_jmax": max_j >= j_max,
    }


def _cyc_pub(c: Cycle) -> dict:
    return {
        "cid": c.cid,
        "k": c.k,
        "length": c.length,
        "odd": c.odd,
        "mu": c.mu,
        "Lambda": c.Lambda,
        "resonant": c.resonant,
        "nodes": list(c.nodes) if c.length <= 48 else list(c.nodes[:24]) + ["…"],
    }


def scan_resonant_levels(
    k_min: int,
    k_max: int,
    j_max: int,
) -> dict:
    """Alle resonanten Zyklen k∈[k_min,k_max], ihre Lift-Tiefen."""
    records = []
    R = {}
    any_infinite_candidate = False

    for k in range(k_min, k_max + 1):
        cycs = all_cycles(k)
        resonant = [c for c in cycs if c.resonant]
        depths = []
        for c in resonant:
            # absolute Lift-Tiefe (Muster-Lift, μ darf fallen)
            tree_all = lift_tree(c, j_max=j_max, require_resonant_lifts=False)
            # resonanter Lift-Pfad
            tree_res = lift_tree(c, j_max=j_max, require_resonant_lifts=True)
            depths.append(tree_res["depth"])
            if tree_res["survives_to_jmax"] or tree_all["survives_to_jmax"]:
                any_infinite_candidate = True
            records.append(
                {
                    "k": k,
                    "cycle": _cyc_pub(c),
                    "lift_depth_pattern": tree_all["depth"],
                    "lift_max_j_pattern": tree_all["L"],
                    "lift_depth_resonant": tree_res["depth"],
                    "lift_max_j_resonant": tree_res["L"],
                    "survives_pattern_to_jmax": tree_all["survives_to_jmax"],
                    "survives_resonant_to_jmax": tree_res["survives_to_jmax"],
                    "tree_resonant_levels": {
                        j: len(v) for j, v in tree_res["levels_summary"].items()
                    },
                }
            )
        R[k] = max(depths) if depths else None
        print(
            f"  k={k:2d}: cycles={len(cycs):3d}  resonant={len(resonant):3d}  "
            f"R_k(resonant-lifts)={R[k]}"
        )
        for rec in records:
            if rec["k"] != k:
                continue
            c = rec["cycle"]
            print(
                f"    μ={c['mu']:.6f} |γ|={c['length']:3d}  "
                f"depth_pattern={rec['lift_depth_pattern']}  "
                f"depth_resonant={rec['lift_depth_resonant']}  "
                f"cid…{c['cid'][:60]}"
            )

    return {
        "records": records,
        "R_k": R,
        "any_survives_to_jmax": any_infinite_candidate,
    }


def main() -> int:
    print("=" * 72)
    print("COLLATZ-LIFT-RESONANZ  —  Resonanzbaum & Lift-Tiefe")
    print("=" * 72)
    print(f"p_c = {P_CRIT:.12f}")
    print("resonant ⇔ μ ≥ p_c")
    print("L(γ)=max j mit Lift bis 2^j;  R_k=max L über resonante γ_k\n")

    j_max = 40
    k_min, k_max = 4, 20

    print(f"[1] Resonante Zyklen k={k_min}..{k_max}, Lift-Versuch bis j={j_max}")
    scan = scan_resonant_levels(k_min, k_max, j_max)

    print("\n[2] R_k = max resonante Lift-Tiefe (nur Lifts mit μ≥p_c)")
    for k, rk in scan["R_k"].items():
        if rk is not None:
            print(f"  R_{k} = {rk}")
        else:
            print(f"  R_{k} = — (keine resonanten Zyklen)")

    print("\n[3] Pattern-Lifts (μ darf fallen) vs. resonante Lifts")
    # aggregiere
    max_pat = max((r["lift_depth_pattern"] for r in scan["records"]), default=0)
    max_res = max((r["lift_depth_resonant"] for r in scan["records"]), default=0)
    print(f"  max depth_pattern  = {max_pat}")
    print(f"  max depth_resonant = {max_res}")
    print(f"  irgendein Ast bis j_max={j_max}? {scan['any_survives_to_jmax']}")

    # Detail: bekannte gefährliche k
    print("\n[4] Detail gefährliche Levels")
    for k in (8, 10, 11, 12):
        recs = [r for r in scan["records"] if r["k"] == k]
        for r in recs:
            print(
                f"  k={k} μ={r['cycle']['mu']:.6f}:  "
                f"pattern→j={r['lift_max_j_pattern']} (depth {r['lift_depth_pattern']}), "
                f"resonant→j={r['lift_max_j_resonant']} (depth {r['lift_depth_resonant']})"
            )

    print("\n[5] Lift-Resonanz-Vermutung")
    finite = (not scan["any_survives_to_jmax"]) and max_res == 0 and max_pat == 0
    # actually max_pat might be 0 if no lift at all - depth 0 means L=k, no lift to k+1
    all_finite = not scan["any_survives_to_jmax"]
    print(
        "  Vermutung: Resonanzbaum hat keine unendlichen Äste.\n"
        f"  Empirie bis j={j_max}: alle resonanten Äste endlich? {all_finite}\n"
        f"  Schärfer: kein resonanter Lift um auch nur 1 Stufe "
        f"(max depth_resonant={max_res})? {max_res == 0}"
    )
    print(
        "\n  Trennung:\n"
        "    lokale Resonanz  ≠  2-adische (liftstabile) Resonanz  ≠  Collatz-Bahn"
    )

    out = Path(__file__).resolve().parent / "collatz_lift_resonance.json"
    payload = {
        "claim": (
            "Lift depth of resonant 2-adic Syracuse cycles; "
            "resonance tree has no infinite branches (empirical); not a proof"
        ),
        "p_crit": P_CRIT,
        "j_max": j_max,
        "k_range": [k_min, k_max],
        "R_k": scan["R_k"],
        "max_depth_pattern": max_pat,
        "max_depth_resonant": max_res,
        "any_survives_to_jmax": scan["any_survives_to_jmax"],
        "conjecture_supported": all_finite,
        "records": scan["records"],
        "definitions": {
            "resonant": "mu >= p_c",
            "liftstable": "exists lifts gamma_j for all j>k with pi(gamma_j)=gamma_k and mu>=p_c",
            "L": "max j with a (pattern) lift to mod 2^j",
            "R_k": "max resonant-lift depth over resonant cycles at level k",
        },
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
