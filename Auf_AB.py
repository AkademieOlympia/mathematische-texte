from __future__ import annotations

from collections import defaultdict
import math
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


# ---------------------------------------------------------
# Daten
# ---------------------------------------------------------

beta_ladder = {
    6:   [6],
    12:  [6, 8, 10, 12],
    18:  [6, 8, 10, 12, 18],
    102: [12, 18, 102],
    108: [18, 102],
    192: [102, 108],
    198: [12, 18, 102],
}

alpha_nodes = {
    "A024": (0, 2, 4),
    "A026": (0, 2, 6),
    "A246": (2, 4, 6),
}

gamma_nodes = {
    "G115": (1, 1, 5),
    "G133": (1, 3, 3),
    "G135": (1, 3, 5),
    "G157": (1, 5, 7),
}


# ---------------------------------------------------------
# Kombinierter ungerichteter Graph
# ---------------------------------------------------------

def make_beta_nodes(ladder):
    nodes = []
    for rho in sorted(ladder):
        for ell in sorted(set(ladder[rho])):
            nodes.append((rho, ell))
    return nodes


def add_vertical_edges(G, ladder):
    for rho in sorted(ladder):
        ells = sorted(set(ladder[rho]))
        for i in range(len(ells) - 1):
            G.add_edge((rho, ells[i]), (rho, ells[i + 1]), kind="vertical")


def add_diagonal_edges_same_ell(G, ladder):
    ell_to_rhos = defaultdict(list)
    for rho in sorted(ladder):
        for ell in sorted(set(ladder[rho])):
            ell_to_rhos[ell].append(rho)

    for ell, rhos in ell_to_rhos.items():
        rhos = sorted(rhos)
        for i in range(len(rhos) - 1):
            G.add_edge((rhos[i], ell), (rhos[i + 1], ell), kind="diagonal")


def add_cascade_edges(G, ladder):
    nodes = set(G.nodes())
    rhos = sorted(ladder)
    for rho1 in rhos:
        source = (rho1, rho1)
        if source not in nodes:
            continue
        for rho2 in rhos:
            if rho2 <= rho1:
                continue
            target = (rho2, rho1)
            if target in nodes:
                G.add_edge(source, target, kind="cascade")


def build_combined_graph():
    G = nx.Graph()

    for node in make_beta_nodes(beta_ladder):
        rho, ell = node
        G.add_node(node, family="beta", rho=rho, ell=ell)

    add_vertical_edges(G, beta_ladder)
    add_diagonal_edges_same_ell(G, beta_ladder)
    add_cascade_edges(G, beta_ladder)

    for name, pat in alpha_nodes.items():
        G.add_node(name, family="alpha", pattern=pat)

    for name, pat in gamma_nodes.items():
        G.add_node(name, family="gamma", pattern=pat)

    G.add_edge((6, 6), "A024", kind="beta_to_alpha")
    G.add_edge((6, 6), "A246", kind="beta_to_alpha")
    G.add_edge((12, 6), "A026", kind="beta_to_alpha")
    G.add_edge((12, 12), "A246", kind="beta_to_alpha")

    G.add_edge((6, 6), "G115", kind="beta_to_gamma")
    G.add_edge((6, 6), "G135", kind="beta_to_gamma")
    G.add_edge((12, 6), "G157", kind="beta_to_gamma")
    G.add_edge((18, 6), "G157", kind="beta_to_gamma")
    G.add_edge((12, 12), "G135", kind="beta_to_gamma")

    return G


# ---------------------------------------------------------
# Energie
# ---------------------------------------------------------

def energy(node, G, model="additive", lam=1.0):
    fam = G.nodes[node]["family"]

    if fam == "beta":
        rho = G.nodes[node]["rho"]
        ell = G.nodes[node]["ell"]

        if model == "additive":
            return rho + lam * ell
        elif model == "multiplicative":
            return rho * ell
        elif model == "sqrt":
            return rho + lam * math.sqrt(ell)
        else:
            raise ValueError("Unbekanntes Modell")

    elif fam == "alpha":
        pat = G.nodes[node]["pattern"]
        return sum(pat) / 2.0

    elif fam == "gamma":
        pat = G.nodes[node]["pattern"]
        return sum(pat)

    return 0.0


# ---------------------------------------------------------
# Gerichteter Graph
# ---------------------------------------------------------

def build_directed_graph(G, model="additive", lam=1.0, eps=1e-12):
    DG = nx.DiGraph()

    for node, data in G.nodes(data=True):
        E = energy(node, G, model=model, lam=lam)
        DG.add_node(node, **data, E=E)

    for u, v, data in G.edges(data=True):
        Eu = DG.nodes[u]["E"]
        Ev = DG.nodes[v]["E"]

        if Eu > Ev + eps:
            DG.add_edge(u, v, kind=data["kind"])
        elif Ev > Eu + eps:
            DG.add_edge(v, u, kind=data["kind"])
        else:
            DG.add_edge(u, v, kind=data["kind"])
            DG.add_edge(v, u, kind=data["kind"])

    return DG


# ---------------------------------------------------------
# Gewichtung
# ---------------------------------------------------------

TYPE_FACTOR = {
    "vertical": 1.00,
    "diagonal": 1.10,
    "cascade": 1.40,
    "beta_to_alpha": 1.60,
    "beta_to_gamma": 1.20,
}

def assign_transition_weights(DG, beta=0.035):
    """
    beta steuert, wie stark Energiedifferenzen bevorzugt werden.
    """
    for u, v, data in DG.edges(data=True):
        dE = DG.nodes[u]["E"] - DG.nodes[v]["E"]
        kind = data["kind"]
        factor = TYPE_FACTOR.get(kind, 1.0)

        raw = factor * math.exp(beta * dE)
        data["dE"] = dE
        data["raw_weight"] = raw

    # normiere pro Startknoten zu Wahrscheinlichkeiten
    for u in DG.nodes():
        outs = list(DG.out_edges(u, data=True))
        if not outs:
            continue
        total = sum(d["raw_weight"] for _, _, d in outs)
        for _, _, d in outs:
            d["prob"] = d["raw_weight"] / total


# ---------------------------------------------------------
# Analyse
# ---------------------------------------------------------

def sinks_by_family(DG):
    sinks = [n for n in DG.nodes() if DG.out_degree(n) == 0]
    by_family = defaultdict(list)
    for s in sinks:
        by_family[DG.nodes[s]["family"]].append(s)
    return sinks, by_family


def print_local_probabilities(DG, node):
    print(f"Lokale Übergänge von {node}:")
    outs = list(DG.out_edges(node, data=True))
    if not outs:
        print("  keine")
        return
    outs.sort(key=lambda x: x[2]["prob"], reverse=True)
    for _, v, d in outs:
        print(
            f"  -> {v:>10} | kind={d['kind']:>13} | "
            f"dE={d['dE']:>7.3f} | p={d['prob']:.4f}"
        )
    print()


def absorption_probabilities(DG):
    """
    Berechnet Absorptionswahrscheinlichkeiten für alle Knoten
    auf alle Senken mittels Rückwärtseinsetzen, da DAG-artig.
    """
    sinks = [n for n in DG.nodes() if DG.out_degree(n) == 0]
    sink_index = {s: i for i, s in enumerate(sinks)}

    order = list(nx.topological_sort(DG))
    order.reverse()

    probs = {n: np.zeros(len(sinks)) for n in DG.nodes()}

    for s, i in sink_index.items():
        probs[s][i] = 1.0

    for u in order:
        if DG.out_degree(u) == 0:
            continue
        vec = np.zeros(len(sinks))
        for _, v, d in DG.out_edges(u, data=True):
            vec += d["prob"] * probs[v]
        probs[u] = vec

    return sinks, probs


def summarize_absorption(DG, start_nodes=None):
    sinks, probs = absorption_probabilities(DG)
    sink_families = {s: DG.nodes[s]["family"] for s in sinks}

    if start_nodes is None:
        start_nodes = [n for n in DG.nodes() if DG.nodes[n]["family"] == "beta"]

    print("Absorptionswahrscheinlichkeiten")
    print("-" * 70)
    for u in start_nodes:
        vec = probs[u]
        p_alpha = sum(vec[i] for i, s in enumerate(sinks) if sink_families[s] == "alpha")
        p_gamma = sum(vec[i] for i, s in enumerate(sinks) if sink_families[s] == "gamma")

        print(f"{str(u):>10} | P(alpha)={p_alpha:.4f} | P(gamma)={p_gamma:.4f}")

    print()
    print("Senken:")
    for i, s in enumerate(sinks):
        print(f"  {i}: {s} ({DG.nodes[s]['family']})")
    print()


def most_probable_path(DG, start):
    """
    Wahrscheinlichster Pfad zu einer Senke über -log(prob).
    """
    H = nx.DiGraph()
    for u, v, d in DG.edges(data=True):
        p = d.get("prob", 0.0)
        if p <= 0:
            continue
        H.add_edge(u, v, cost=-math.log(p))

    sinks = [n for n in DG.nodes() if DG.out_degree(n) == 0]
    best_path = None
    best_cost = float("inf")
    best_sink = None

    for s in sinks:
        try:
            path = nx.shortest_path(H, source=start, target=s, weight="cost")
            cost = nx.shortest_path_length(H, source=start, target=s, weight="cost")
            if cost < best_cost:
                best_cost = cost
                best_path = path
                best_sink = s
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            pass

    if best_path is None:
        return None, None, None

    return best_path, best_sink, math.exp(-best_cost)


# ---------------------------------------------------------
# Hauptprogramm
# ---------------------------------------------------------

def main():
    G = build_combined_graph()
    DG = build_directed_graph(G, model="additive", lam=1.0)
    assign_transition_weights(DG, beta=0.035)

    print_local_probabilities(DG, (198, 102))
    print_local_probabilities(DG, (102, 102))
    print_local_probabilities(DG, (12, 12))
    print_local_probabilities(DG, (6, 6))

    summarize_absorption(DG, start_nodes=[
        (198, 102),
        (192, 108),
        (102, 102),
        (18, 18),
        (12, 12),
        (6, 6),
    ])

    for start in [(198, 102), (192, 108), (102, 102), (18, 18), (12, 12), (6, 6)]:
        path, sink, p = most_probable_path(DG, start)
        print(f"Start: {start}")
        print(f"  wahrscheinlichster Pfad: {path}")
        print(f"  endet in: {sink}")
        print(f"  Pfadstärke ~ {p:.6f}")
        print()


if __name__ == "__main__":
    main()