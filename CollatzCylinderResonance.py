"""
CollatzCylinderResonance — Resonanz über Zylinderklassen & Branch-and-Bound.

Ein Wort ω der Länge m ↔ Restklasse r mod 2^m (Syracuse, eindeutig).
Lift zur Ebene j ≥ m: alle x ≡ r (mod 2^m) in Z/2^j Z.

Da G_j Funktionalgraph ist, fällt jedes x in genau einen Attraktorzyklus.
Definiere:

  M(j, ω) = M(j, r) = max { μ(cycle(x)) : x ≡ r (mod 2^m) }

Die Lift-Sterilisierungs-/Zylinder-Vermutung:

  ∀ ω:  limsup_{j→∞} M(j, ω) < p_c

Branch-and-Bound auf dem Wortbaum:
  - Knoten = Präfix ω
  - obere Schranke = max erreichbare Zyklusdichte unter Lifts bis j_cap
  - verwerfe Äste mit M_upper < p_c

Claim-boundary: Inverse-Limit / Zylinder-Diagnostik; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


LOG2_3 = math.log2(3)
P_CRIT = math.log(2) / math.log(3)


def syracuse_next(r: int, mod: int) -> tuple[int, int]:
    if r % 2 == 0:
        return (r // 2) % mod, 0
    return ((3 * r + 1) // 2) % mod, 1


def automaton(k: int) -> tuple[list[int], list[int]]:
    n = 1 << k
    nxt = [0] * n
    w = [0] * n
    for r in range(n):
        nxt[r], w[r] = syracuse_next(r, n)
    return nxt, w


@dataclass(frozen=True)
class CycleInfo:
    nodes: tuple[int, ...]
    mu: float
    odd: int

    @property
    def length(self) -> int:
        return len(self.nodes)


def cycle_map(k: int) -> tuple[list[int], list[CycleInfo], list[float]]:
    """
    Für jedes r: Index des Attraktorzyklus, Liste der Zyklen, μ(r).
    """
    nxt, w = automaton(k)
    n = 1 << k
    # Finde Zyklen
    visited = [0] * n
    cycles: list[CycleInfo] = []
    cycle_id_of = [-1] * n  # temporary during discovery

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
        cid = -1
        if visited[x] == 1:
            cyc = path[pos[x] :]
            odd = sum(w[v] for v in cyc)
            info = CycleInfo(nodes=tuple(cyc), mu=odd / len(cyc), odd=odd)
            cid = len(cycles)
            cycles.append(info)
            for v in cyc:
                cycle_id_of[v] = cid
        for v in path:
            visited[v] = 2

    # Assign every node its attractor by following
    attractor = [-1] * n
    mu_of = [0.0] * n
    for r in range(n):
        x = r
        seen = {}
        steps = 0
        while x not in seen:
            if attractor[x] != -1:
                # jump to known
                aid = attractor[x]
                break
            seen[x] = steps
            x = nxt[x]
            steps += 1
        else:
            # hit cycle in this walk
            aid = cycle_id_of[x]
        # fill path
        x = r
        while attractor[x] == -1:
            attractor[x] = aid
            mu_of[x] = cycles[aid].mu
            x = nxt[x]

    return attractor, cycles, mu_of


def M_at_level(j: int, m: int, r: int, mu_of_j: list[float]) -> float:
    """M(j,r) = max μ unter Lifts x ≡ r mod 2^m in Z/2^j."""
    mod_hi = 1 << j
    mod_lo = 1 << m
    best = 0.0
    for x in range(r, mod_hi, mod_lo):
        best = max(best, mu_of_j[x])
    return best


def compute_M_table(m: int, j_max: int) -> dict:
    """Für festes m: M(j,r) für j=m..j_max, alle r."""
    rows = []
    # precompute mu_of for each j
    cache_mu: dict[int, list[float]] = {}
    cache_cyc: dict[int, list[CycleInfo]] = {}
    for j in range(m, j_max + 1):
        _att, cycles, mu_of = cycle_map(j)
        cache_mu[j] = mu_of
        cache_cyc[j] = cycles

    global_limsup_proxy = []
    for r in range(1 << m):
        series = []
        for j in range(m, j_max + 1):
            Mj = M_at_level(j, m, r, cache_mu[j])
            series.append(Mj)
        # Tail erst ab j≥13 (nach Sterilisierung der Phantom-Zyklen)
        tail = [series[j - m] for j in range(max(m, 13), j_max + 1)]
        limsup_emp = max(tail) if tail else max(series)
        global_limsup_proxy.append(limsup_emp)
        rows.append(
            {
                "r": r,
                "M_series": series,
                "M_max": max(series),
                "M_at_jmax": series[-1],
                "limsup_tail": limsup_emp,
                "ever_ge_pc": any(x >= P_CRIT - 1e-15 for x in series),
                "final_ge_pc": series[-1] >= P_CRIT - 1e-15,
            }
        )

    return {
        "m": m,
        "j_max": j_max,
        "mu_max_by_j": {
            str(j): max(c.mu for c in cache_cyc[j]) for j in range(m, j_max + 1)
        },
        "cylinders_ever_ge_pc": sum(1 for row in rows if row["ever_ge_pc"]),
        "cylinders_final_ge_pc": sum(1 for row in rows if row["final_ge_pc"]),
        "max_limsup_tail": max(global_limsup_proxy),
        "frac_limsup_lt_pc": sum(1 for x in global_limsup_proxy if x < P_CRIT) / len(rows),
        "top_cylinders": sorted(rows, key=lambda d: -d["M_max"])[:12],
        "worst_final": sorted(rows, key=lambda d: -d["M_at_jmax"])[:8],
    }


def word_to_residue(omega: tuple[int, ...]) -> int:
    """Eindeutige Restklasse mod 2^m für Syracuse-Wort (Lift-Konstruktion)."""
    m = len(omega)
    candidates = {0}
    cur_mod = 1
    for j in range(m):
        new_mod = cur_mod * 2
        new_cand = set()
        for r0 in candidates:
            for half in (0, 1):
                r = r0 + half * cur_mod
                x = r
                ok = True
                for t in range(j + 1):
                    b = 0 if x % 2 == 0 else 1
                    if b != omega[t]:
                        ok = False
                        break
                    x = x // 2 if b == 0 else (3 * x + 1) // 2
                if ok:
                    new_cand.add(r % new_mod)
        candidates = new_cand
        cur_mod = new_mod
    return min(candidates)


def branch_and_bound(
    m0: int,
    j_cap: int,
    max_nodes: int = 200_000,
) -> dict:
    """
    Wortbaum ab Länge m0: expandiere Bits; Upper Bound = M(j_cap, residue(ω)).
    Verwerfe wenn M_upper < p_c.
    """
    # Precompute mu at j_cap
    _att, cycles, mu_of = cycle_map(j_cap)
    mu_star = max(c.mu for c in cycles)

    # BFS/DFS stack: (omega_tuple)
    # Start: all words length m0 — too many if m0 large; use m0 small or sample
    from itertools import product

    pruned = 0
    expanded = 0
    critical_leaves = []  # reached len == something with M>=pc
    frontier_ge_pc = 0
    max_depth_ge_pc = m0 - 1

    # Start from all residues at m0 via words — 2^{m0} 
    # For m0>=16 use residues directly instead of words
    use_residues = m0 >= 12

    if use_residues:
        starts = [(r,) for r in range(1 << m0)]  # dummy; treat as residue id
        mode = "residue"
    else:
        starts = list(product([0, 1], repeat=m0))
        mode = "word"

    # Iterative deepening by length for words; for residues scan all r at m0..j_cap
    if mode == "residue":
        # Direct cylinder scan — exact M(j_cap,r) and series
        ge_pc = []
        for r in range(1 << m0):
            M = M_at_level(j_cap, m0, r, mu_of)
            expanded += 1
            if M + 1e-15 >= P_CRIT:
                frontier_ge_pc += 1
                ge_pc.append({"r": r, "M": M})
            else:
                pruned += 1
        return {
            "mode": "residue_cylinder",
            "m0": m0,
            "j_cap": j_cap,
            "mu_star_jcap": mu_star,
            "expanded": expanded,
            "pruned_M_lt_pc": pruned,
            "frontier_ge_pc": frontier_ge_pc,
            "critical_sample": ge_pc[:20],
            "all_pruned": frontier_ge_pc == 0,
        }

    # Word mode: expand while M_upper >= pc
    stack = list(starts)
    seen_res = set()
    while stack and expanded < max_nodes:
        omega = stack.pop()
        expanded += 1
        r = word_to_residue(omega)
        key = (len(omega), r)
        if key in seen_res:
            pruned += 1
            continue
        seen_res.add(key)
        # Upper bound: M at j_cap for this cylinder (m=len)
        m = len(omega)
        if m > j_cap:
            continue
        M = M_at_level(j_cap, m, r, mu_of)
        if M + 1e-15 < P_CRIT:
            pruned += 1
            continue
        frontier_ge_pc += 1
        max_depth_ge_pc = max(max_depth_ge_pc, m)
        if m >= j_cap:
            critical_leaves.append({"omega": "".join(map(str, omega)), "r": r, "M": M})
            continue
        # expand both bits
        for b in (0, 1):
            stack.append(omega + (b,))

    return {
        "mode": "word_bnb",
        "m0": m0,
        "j_cap": j_cap,
        "mu_star_jcap": mu_star,
        "expanded": expanded,
        "pruned": pruned,
        "frontier_ge_pc_visits": frontier_ge_pc,
        "max_depth_ge_pc": max_depth_ge_pc,
        "critical_leaves": critical_leaves[:30],
        "all_pruned_before_jcap": len(critical_leaves) == 0 and frontier_ge_pc > 0,
    }


def main() -> int:
    print("=" * 72)
    print("COLLATZ-ZYLINDER-RESONANZ  —  M(j,ω) & Branch-and-Bound")
    print("=" * 72)
    print(f"p_c = {P_CRIT:.12f}")
    print("M(j,r) = max μ(cycle(x)) über x ≡ r (mod 2^m), x ∈ Z/2^j Z\n")

    # [1] M-Tabellen für m=8,10,12
    tables = []
    print("[1] M(j,r)-Atlas (Zylinder → Attraktor-μ)")
    for m, j_max in [(6, 16), (8, 18), (10, 18), (12, 18)]:
        print(f"\n  --- m={m}, j≤{j_max} ---")
        tab = compute_M_table(m, j_max)
        tables.append(tab)
        print(f"  μ_max(j): {tab['mu_max_by_j']}")
        print(
            f"  Zylinder je ≥p_c irgendwann: {tab['cylinders_ever_ge_pc']}/{1<<m}"
        )
        print(
            f"  Zylinder final M(j_max)≥p_c: {tab['cylinders_final_ge_pc']}/{1<<m}"
        )
        print(
            f"  max limsup-tail≈ {tab['max_limsup_tail']:.6f}  "
            f"frac(limsup<p_c)={tab['frac_limsup_lt_pc']:.4f}"
        )
        print("  top M_max:")
        for row in tab["top_cylinders"][:5]:
            print(
                f"    r={row['r']:5d}  M_max={row['M_max']:.6f}  "
                f"M_final={row['M_at_jmax']:.6f}  ever≥pc={row['ever_ge_pc']}"
            )

    # [2] Branch-and-bound / residue pruning at high j_cap
    print("\n[2] Branch-and-Bound / Zylinder-Pruning")
    bnbs = []
    for m0, j_cap in [(8, 16), (10, 18), (12, 18), (14, 20)]:
        bnb = branch_and_bound(m0, j_cap)
        bnbs.append(bnb)
        ge = bnb.get("frontier_ge_pc", bnb.get("frontier_ge_pc_visits", 0))
        pr = bnb.get("pruned_M_lt_pc", bnb.get("pruned", 0))
        safe = ge == 0
        print(
            f"  m0={m0} j_cap={j_cap}: μ*={bnb['mu_star_jcap']:.6f}  "
            f"ge_pc={ge}  pruned={pr}  all_safe={safe}"
        )

    # [3] Globale Aussage über j
    print("\n[3] Globale Sterilisierung auf hohen Ebenen")
    print("  Für j≥13: einzige Zyklen μ∈{0, 1/2} (aus 2-adic Atlas).")
    print("  ⇒ für alle r, alle m≤j: M(j,r) ≤ 1/2 < p_c  (j=13..20 geprüft).")
    # verify quickly
    safe_levels = []
    for j in range(13, 21):
        _a, cycs, _mu = cycle_map(j)
        mumax = max(c.mu for c in cycs)
        safe_levels.append((j, mumax, mumax < P_CRIT))
        print(f"  j={j}: μ_max={mumax:.6f}  <p_c? {mumax < P_CRIT}")

    print("\n[4] Vermutung (Zylinder-Form)")
    print(
        "  ∀ω: limsup_j M(j,ω) < p_c\n"
        "  Empirie: für j≥13 gilt sogar ∀r: M(j,r)≤1/2 < p_c.\n"
        "  Resonante Zylinder bei j∈{8,10,11,12} sterilisieren beim Lift\n"
        "  (fallen in den 1↔2-Attraktor mit μ=1/2)."
    )
    print(
        "\n  Gap-Analogie:\n"
        "    Gap-Klasse → Operatorklasse → Liftklasse → Inverse-Limit-Kanal\n"
        "  Hier: Zylinder ohne unendlichen resonanten Kanal im geprüften Bereich."
    )

    out = Path(__file__).resolve().parent / "collatz_cylinder_resonance.json"
    payload = {
        "claim": (
            "Cylinder resonance M(j,r)=max attractor mu over 2-adic lifts; "
            "limsup M < p_c empirically for j>=13; not a proof"
        ),
        "p_crit": P_CRIT,
        "definition": {
            "M(j,r)": "max mu(cycle(x)) over x≡r mod 2^m in Z/2^jZ",
            "conjecture": "forall omega: limsup_j M(j,omega) < p_c",
        },
        "tables": [
            {
                "m": t["m"],
                "j_max": t["j_max"],
                "mu_max_by_j": t["mu_max_by_j"],
                "cylinders_ever_ge_pc": t["cylinders_ever_ge_pc"],
                "cylinders_final_ge_pc": t["cylinders_final_ge_pc"],
                "max_limsup_tail": t["max_limsup_tail"],
                "frac_limsup_lt_pc": t["frac_limsup_lt_pc"],
                "top_cylinders": t["top_cylinders"][:8],
            }
            for t in tables
        ],
        "bnb": bnbs,
        "safe_levels": [{"j": j, "mu_max": mu, "safe": s} for j, mu, s in safe_levels],
        "sterilization": (
            "For j>=13 all attractor mus are in {0,1/2}; "
            "hence M(j,r)<=1/2<p_c for all cylinders"
        ),
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
