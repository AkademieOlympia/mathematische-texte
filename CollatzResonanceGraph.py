"""
CollatzResonanceGraph — Klassen-Dynamik Γ auf Erstabstiegsklassen.

Äquivalenz:
  n ∼ m  ⇔  (a_τ, τ, c) gleich
  G = ℕ / ∼

Resonanz (Analog zu R=(ρ,Σ)):
  R_C = (λ, ρ, σ)
    λ = log2(A) = log2(3^a / 2^τ)
    ρ = |c| / (2^τ + 3^a)
    σ = |B| / |A − 1|   (affiner Fixpunkt-Maßstab)

Klassenabbildung (Relation):
  Γ([n]) ∋ [T(n)]   (elementarer Collatz-Schritt)
  Δ([n]) ∋ [T^τ(n)] (Bild nach Erstabstieg)

Frage: Ist Γ/Δ funktionswertig? Gibt es gerichtete Zyklen außer 4-2-1?
Zielbild: jede Klasse hat λ<0 und erreicht (0,1,0).

Claim-boundary: Graph-Diagnostik; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


# ---------------------------------------------------------------------------
# Syracuse first-descent (wie CollatzResonanceAtlas)
# ---------------------------------------------------------------------------

def collatz_step(x: int) -> int:
    return x // 2 if x % 2 == 0 else 3 * x + 1


def first_descent(n: int, max_steps: int = 100_000) -> tuple[int, int, int, int]:
    """(a, τ, c, T^τ(n)) in Syracuse-Konvention."""
    if n <= 1:
        return 0, 0, 0, n
    a = 0
    tau = 0
    c = 0
    x = n
    for _ in range(max_steps):
        if x % 2 == 0:
            tau += 1
            x //= 2
        else:
            c = 3 * c + (1 << tau)
            a += 1
            tau += 1
            x = (3 * x + 1) // 2
        if x < n:
            return a, tau, c, x
    raise RuntimeError(f"no descent for n={n}")


ClassKey = tuple[int, int, int]  # (a, tau, c)


def class_of(n: int) -> ClassKey:
    a, tau, c, _ = first_descent(n)
    return a, tau, c


@dataclass(frozen=True)
class CollatzClass:
    a: int
    tau: int
    c: int

    @property
    def key(self) -> ClassKey:
        return self.a, self.tau, self.c

    @property
    def A(self) -> Fraction:
        if self.tau == 0:
            return Fraction(3**self.a, 1)
        return Fraction(3**self.a, 1 << self.tau)

    @property
    def B(self) -> Fraction:
        if self.tau == 0:
            return Fraction(self.c, 1)
        return Fraction(self.c, 1 << self.tau)

    def resonance(self) -> tuple[float, Fraction, float]:
        """R_C = (λ, ρ, σ)."""
        A = self.A
        B = self.B
        if A == 0:
            lam = float("-inf")
        else:
            lam = math.log2(float(A))
        rho = Fraction(abs(self.c), (1 << self.tau) + 3**self.a) if self.tau or self.a else Fraction(0)
        if A == 1:
            sigma = float("inf")
        else:
            sigma = float(abs(B) / abs(A - 1))
        return lam, rho, sigma


# ---------------------------------------------------------------------------
# Atlas + Graph
# ---------------------------------------------------------------------------

TRIVIAL = (0, 1, 0)  # alle Geraden: ein /2


def build_members(N: int) -> dict[ClassKey, list[int]]:
    members: dict[ClassKey, list[int]] = defaultdict(list)
    for n in range(2, N):
        members[class_of(n)].append(n)
    return members


def build_relation(
    members: dict[ClassKey, list[int]],
    mode: str,
) -> dict[ClassKey, set[ClassKey]]:
    """
    mode='step':   Γ via T (elementarer Schritt)
    mode='descent': Δ via T^τ (Erstabstiegsbild)
    """
    rel: dict[ClassKey, set[ClassKey]] = defaultdict(set)
    for G, ns in members.items():
        for n in ns:
            if mode == "step":
                m = collatz_step(n)
            else:
                _a, _t, _c, m = first_descent(n)
            if m <= 1:
                # 1 als Senke: Klasse (0,0,0) Marker
                rel[G].add((0, 0, 0))
            else:
                rel[G].add(class_of(m))
    return rel


def well_defined_stats(rel: dict[ClassKey, set[ClassKey]]) -> dict:
    out_deg = [len(v) for v in rel.values()]
    functional = sum(1 for d in out_deg if d == 1)
    return {
        "nodes": len(rel),
        "functional_nodes": functional,
        "functional_rate": functional / len(rel) if rel else None,
        "max_out_degree": max(out_deg) if out_deg else 0,
        "mean_out_degree": sum(out_deg) / len(out_deg) if out_deg else 0,
        "multi_valued": sum(1 for d in out_deg if d > 1),
    }


def find_directed_cycles(
    rel: dict[ClassKey, set[ClassKey]],
    max_cycles: int = 50,
    max_len: int = 20,
) -> list[list[ClassKey]]:
    """Einfache gerichtete Zyklen via DFS (begrenzt)."""
    cycles: list[list[ClassKey]] = []
    nodes = list(rel.keys())

    def dfs(start: ClassKey, node: ClassKey, path: list[ClassKey], stack_set: set[ClassKey]):
        if len(cycles) >= max_cycles or len(path) > max_len:
            return
        for nxt in rel.get(node, ()):
            if nxt == start and len(path) >= 1:
                cycles.append(path + [start])
                if len(cycles) >= max_cycles:
                    return
            elif nxt not in stack_set and nxt in rel:
                stack_set.add(nxt)
                path.append(nxt)
                dfs(start, nxt, path, stack_set)
                path.pop()
                stack_set.remove(nxt)

    for s in nodes:
        dfs(s, s, [s], {s})
        if len(cycles) >= max_cycles:
            break

    # normalisiere Rotation / Richtung
    uniq = []
    seen = set()
    for cyc in cycles:
        body = cyc[:-1]
        if not body:
            continue
        # kanonische Rotation
        i = body.index(min(body))
        norm = tuple(body[i:] + body[:i])
        if norm not in seen:
            seen.add(norm)
            uniq.append(list(norm) + [norm[0]])
    return uniq


def scc_tarjan(rel: dict[ClassKey, set[ClassKey]]) -> list[list[ClassKey]]:
    index = 0
    stack: list[ClassKey] = []
    onstack: set[ClassKey] = set()
    indices: dict[ClassKey, int] = {}
    lowlink: dict[ClassKey, int] = {}
    sccs: list[list[ClassKey]] = []

    def strongconnect(v: ClassKey):
        nonlocal index
        indices[v] = index
        lowlink[v] = index
        index += 1
        stack.append(v)
        onstack.add(v)
        for w in rel.get(v, ()):
            if w not in rel and w not in indices:
                # Senke ohne eigene Out-Edges: als Blatt behandeln
                indices.setdefault(w, index)
                lowlink.setdefault(w, index)
                index += 1
                stack.append(w)
                onstack.add(w)
                # pop singleton later via fallthrough
            if w not in indices:
                # ensure w is in graph view
                rel.setdefault(w, set())
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif w in onstack:
                lowlink[v] = min(lowlink[v], indices[w])
        if lowlink[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                onstack.remove(w)
                comp.append(w)
                if w == v:
                    break
            sccs.append(comp)

    # include sinks
    nodes = set(rel.keys())
    for vs in rel.values():
        nodes |= vs
    for v in nodes:
        rel.setdefault(v, set())
    for v in list(rel.keys()):
        if v not in indices:
            strongconnect(v)
    return sccs


def reaches_trivial(
    rel: dict[ClassKey, set[ClassKey]],
    members: dict[ClassKey, list[int]],
    max_depth: int = 10_000,
) -> dict:
    """Erreichbarkeit von TRIVIAL=(0,1,0) oder Senke (0,0,0) im Klassengraphen."""
    # reverse graph for ancestors of trivial
    rev: dict[ClassKey, set[ClassKey]] = defaultdict(set)
    nodes = set(rel.keys())
    for u, vs in rel.items():
        for v in vs:
            rev[v].add(u)
            nodes.add(v)

    targets = {TRIVIAL, (0, 0, 0)}
    reach: set[ClassKey] = set()
    dq = deque(targets)
    while dq:
        v = dq.popleft()
        if v in reach:
            continue
        reach.add(v)
        for u in rev.get(v, ()):
            if u not in reach:
                dq.append(u)

    all_nodes = set(members.keys())
    missing = sorted(all_nodes - reach, key=lambda g: (-len(members[g]), g))
    return {
        "nodes": len(all_nodes),
        "reach_trivial": len(all_nodes & reach),
        "miss_trivial": len(missing),
        "miss_samples": [
            {"G": g, "size": len(members[g]), "sample": members[g][:6]} for g in missing[:12]
        ],
    }


def analyze_lambda(members: dict[ClassKey, list[int]]) -> dict:
    pos = []
    zero = []
    neg = 0
    for g, ns in members.items():
        cls = CollatzClass(*g)
        lam, rho, sigma = cls.resonance()
        if lam > 0:
            pos.append({"G": g, "lam": lam, "size": len(ns), "sample": ns[:5]})
        elif lam == 0:
            zero.append(g)
        else:
            neg += 1
    return {
        "lambda_neg": neg,
        "lambda_zero": len(zero),
        "lambda_pos": len(pos),
        "lambda_pos_samples": pos[:15],
    }


def check_gamma_welldefined_examples(
    members: dict[ClassKey, list[int]], mode: str, limit: int = 8
) -> list[dict]:
    """Beispiele mehrwertiger Übergänge."""
    rel = build_relation(members, mode)
    ex = []
    for G, outs in rel.items():
        if len(outs) > 1:
            # zeige Partner
            pairs = []
            for n in members[G][:20]:
                m = collatz_step(n) if mode == "step" else first_descent(n)[3]
                tgt = (0, 0, 0) if m <= 1 else class_of(m)
                pairs.append((n, m, tgt))
            ex.append(
                {
                    "G": G,
                    "out_degree": len(outs),
                    "targets": sorted(outs),
                    "pairs": pairs[:10],
                    "size": len(members[G]),
                }
            )
            if len(ex) >= limit:
                break
    return ex


def main() -> int:
    N = 10_000
    print("=" * 72)
    print("COLLATZ-RESONANZGRAPH  —  Γ : G → P(G)")
    print("=" * 72)
    print(f"N={N},  R_C=(λ,ρ,σ),  λ=log2(A),  trivial=(0,1,0)")

    members = build_members(N)
    print(f"\n[1] Klassen |G| = {len(members)}")

    lam = analyze_lambda(members)
    print("\n[2] Expansion λ auf Klassen")
    print(f"  λ<0 : {lam['lambda_neg']}")
    print(f"  λ=0 : {lam['lambda_zero']}")
    print(f"  λ>0 : {lam['lambda_pos']}")
    if lam["lambda_pos_samples"]:
        for s in lam["lambda_pos_samples"][:5]:
            print(f"    G={s['G']} lam={s['lam']:.4f} size={s['size']} sample={s['sample']}")

    # Resonanz-Stichprobe
    print("\n[3] R_C Stichprobe (größte Klassen)")
    largest = sorted(members.items(), key=lambda kv: -len(kv[1]))[:8]
    for G, ns in largest:
        cls = CollatzClass(*G)
        L, rho, sigma = cls.resonance()
        print(
            f"  size={len(ns):4d}  G={G}  λ={L:.4f}  ρ={rho}  σ={sigma:.4f}  "
            f"sample={ns[:5]}"
        )

    rel_step = build_relation(members, "step")
    rel_desc = build_relation(members, "descent")
    st_step = well_defined_stats(rel_step)
    st_desc = well_defined_stats(rel_desc)

    print("\n[4] Wohldefiniertheit")
    print(
        f"  Γ (ein Schritt):   funktional={st_step['functional_nodes']}/{st_step['nodes']}  "
        f"({st_step['functional_rate']:.3f})  "
        f"max_out={st_step['max_out_degree']}  mean_out={st_step['mean_out_degree']:.2f}"
    )
    print(
        f"  Δ (Erstabstieg):   funktional={st_desc['functional_nodes']}/{st_desc['nodes']}  "
        f"({st_desc['functional_rate']:.3f})  "
        f"max_out={st_desc['max_out_degree']}  mean_out={st_desc['mean_out_degree']:.2f}"
    )

    print("\n[5] Mehrwertige Γ-Beispiele")
    for ex in check_gamma_welldefined_examples(members, "step", limit=5):
        print(
            f"  G={ex['G']} out={ex['out_degree']} size={ex['size']}  "
            f"targets={ex['targets'][:6]}"
        )
        print(f"    pairs={ex['pairs'][:4]}")

    print("\n[6] Mehrwertige Δ-Beispiele")
    for ex in check_gamma_welldefined_examples(members, "descent", limit=5):
        print(
            f"  G={ex['G']} out={ex['out_degree']} size={ex['size']}  "
            f"targets={ex['targets'][:6]}"
        )
        print(f"    pairs={ex['pairs'][:4]}")

    reach_g = reaches_trivial(rel_step, members)
    reach_d = reaches_trivial(rel_desc, members)
    print("\n[7] Erreichbarkeit von (0,1,0) / Senke")
    print(
        f"  via Γ: reach={reach_g['reach_trivial']}/{reach_g['nodes']}  "
        f"miss={reach_g['miss_trivial']}"
    )
    print(
        f"  via Δ: reach={reach_d['reach_trivial']}/{reach_d['nodes']}  "
        f"miss={reach_d['miss_trivial']}"
    )
    if reach_d["miss_samples"]:
        print("  Δ-Miss-Samples:")
        for s in reach_d["miss_samples"][:5]:
            print(f"    G={s['G']} size={s['size']} sample={s['sample']}")

    print("\n[8] Gerichtete Zyklen (Δ-Graph)")
    cycles = find_directed_cycles(rel_desc, max_cycles=80, max_len=12)
    # filter self-loops and trivial-ish
    nontrivial = []
    for cyc in cycles:
        body = cyc[:-1]
        if TRIVIAL in body and len(set(body)) <= 2:
            continue
        if set(body) <= {TRIVIAL, (0, 0, 0)}:
            continue
        nontrivial.append(cyc)
    print(f"  Zyklen gefunden (roh): {len(cycles)}")
    print(f"  nicht-triviale Auszüge: {len(nontrivial)}")
    for cyc in nontrivial[:10]:
        print(f"    {[list(g) for g in cyc]}")

    print("\n[9] SCCs (Δ-Graph)")
    sccs = scc_tarjan({k: set(v) for k, v in rel_desc.items()})
    multi = [c for c in sccs if len(c) > 1]
    print(f"  SCCs total={len(sccs)}  multi-node={len(multi)}")
    for comp in sorted(multi, key=len, reverse=True)[:8]:
        print(f"    size={len(comp)}  sample={comp[:5]}")

    # Funktionaler Kern: nur Knoten mit out-degree 1
    print("\n[10] Funktionaler Kern (out-degree = 1)")
    func_nodes = {g for g, outs in rel_desc.items() if len(outs) == 1}
    func_rel = {g: rel_desc[g] for g in func_nodes}
    # Ziel der funktionalen Kante muss nicht funktional sein
    func_cycles = find_directed_cycles(func_rel, max_cycles=40, max_len=15)
    print(f"  funktionale Knoten: {len(func_nodes)}/{len(rel_desc)}")
    print(f"  Zyklen nur im funktionalen Kern: {len(func_cycles)}")
    for cyc in func_cycles[:6]:
        print(f"    {[list(g) for g in cyc]}")

    # Echte Zahlenzyklen in Klassenmitglieder?
    print("\n[11] Echte Collatz-Zyklen in Mitgliedern (N-Fenster)")
    real_cycles = []
    for n in range(2, min(N, 2000)):
        x = n
        seen = {}
        for i in range(50):
            if x in seen:
                if x == n and i - seen[x] >= 1:
                    cyc_len = i - seen[x]
                    if n not in (1, 2, 4) and x in (1, 2, 4):
                        break
                    real_cycles.append((n, cyc_len))
                break
            seen[x] = i
            x = collatz_step(x)
            if x < 1:
                break
    # nur Zyklen die nicht in 4-2-1 fallen
    nontrivial_real = [rc for rc in real_cycles if rc[0] not in (1, 2, 4)]
    print(f"  detektierte Rückkehren: {len(real_cycles)}  (meist Eintritt in 4-2-1)")
    print(f"  nichttriviale geschlossene n→…→n: {len(nontrivial_real)}")

    # Zielsatz-Diagnostik: λ<0 und Pfad nach trivial
    print("\n[12] Zielsatz-Diagnostik (endlich, N-beschränkt)")
    all_neg = lam["lambda_pos"] == 0
    all_reach = reach_d["miss_trivial"] == 0
    print(f"  alle Klassen λ<0?           {all_neg}")
    print(f"  alle Klassen → trivial via Δ? {all_reach}")
    print(
        "  Bruch der Gap-Analogie: Γ/Δ sind Relationen, keine Funktionen.\n"
        "  Mehrwertigkeit vor allem bei (0,1,0): T(2n)=n streut in viele Klassen.\n"
        "  SCC/Zyklen im Mengengraphen ≠ Collatz-Zyklen (Artefakt der Relation)."
    )

    out = Path(__file__).resolve().parent / "collatz_resonance_graph.json"
    payload = {
        "claim": (
            "Collatz class graph on first-descent classes; R_C=(lambda,rho,sigma); "
            "Gamma is set-valued; search cycles / absorption; not a proof"
        ),
        "N": N,
        "num_classes": len(members),
        "lambda": lam,
        "well_defined": {"Gamma_step": st_step, "Delta_descent": st_desc},
        "reachability": {"Gamma": reach_g, "Delta": reach_d},
        "cycles_delta_nontrivial": [
            [list(g) for g in cyc] for cyc in nontrivial[:20]
        ],
        "scc_multi_count": len(multi),
        "scc_multi_sample": [[list(g) for g in c[:8]] for c in multi[:5]],
        "functional_core": {
            "nodes": len(func_nodes),
            "cycles": [[list(g) for g in cyc] for cyc in func_cycles[:15]],
        },
        "nontrivial_real_cycles": nontrivial_real[:20],
        "largest_resonances": [
            {
                "G": list(G),
                "size": len(ns),
                "R_C": {
                    "lambda": CollatzClass(*G).resonance()[0],
                    "rho": str(CollatzClass(*G).resonance()[1]),
                    "sigma": CollatzClass(*G).resonance()[2],
                },
                "sample": ns[:8],
            }
            for G, ns in largest
        ],
        "target_theorem_finite": {
            "all_lambda_neg": all_neg,
            "all_reach_trivial_via_Delta": all_reach,
        },
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
