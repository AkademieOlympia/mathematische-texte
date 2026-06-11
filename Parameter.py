from __future__ import annotations

from collections import defaultdict
import math
import itertools
import numpy as np
import networkx as nx


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
# Graphaufbau
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
        return sum(G.nodes[node]["pattern"]) / 2.0

    elif fam == "gamma":
        return sum(G.nodes[node]["pattern"])

    return 0.0


def build_directed_graph(G, model="additive", lam=1.0, eps=1e-12):
    DG = nx.DiGraph()

    for node, data in G.nodes(data=True):
        DG.add_node(node, **data, E=energy(node, G, model=model, lam=lam))

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
# Gewichtung und Absorption
# ---------------------------------------------------------

def assign_transition_weights(DG, beta_strength=0.035, type_factor=None):
    if type_factor is None:
        type_factor = {
            "vertical": 1.00,
            "diagonal": 1.10,
            "cascade": 1.40,
            "beta_to_alpha": 1.60,
            "beta_to_gamma": 1.20,
        }

    for u, v, data in DG.edges(data=True):
        dE = DG.nodes[u]["E"] - DG.nodes[v]["E"]
        factor = type_factor.get(data["kind"], 1.0)
        raw = factor * math.exp(beta_strength * dE)
        data["dE"] = dE
        data["raw_weight"] = raw

    for u in DG.nodes():
        outs = list(DG.out_edges(u, data=True))
        if not outs:
            continue
        total = sum(d["raw_weight"] for _, _, d in outs)
        for _, _, d in outs:
            d["prob"] = d["raw_weight"] / total


def absorption_probabilities(DG):
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


def summarize_channel_probabilities(DG, start_nodes):
    sinks, probs = absorption_probabilities(DG)
    sink_families = {s: DG.nodes[s]["family"] for s in sinks}

    rows = []
    for u in start_nodes:
        vec = probs[u]
        p_alpha = sum(vec[i] for i, s in enumerate(sinks) if sink_families[s] == "alpha")
        p_gamma = sum(vec[i] for i, s in enumerate(sinks) if sink_families[s] == "gamma")
        rows.append((u, p_alpha, p_gamma))

    mean_alpha = sum(r[1] for r in rows) / len(rows)
    mean_gamma = sum(r[2] for r in rows) / len(rows)
    return rows, mean_alpha, mean_gamma


# ---------------------------------------------------------
# Parameterstudie
# ---------------------------------------------------------

def parameter_study():
    G = build_combined_graph()
    DG = build_directed_graph(G, model="additive", lam=1.0)

    start_nodes = [
        (198, 102),
        (192, 108),
        (102, 102),
        (18, 18),
        (12, 12),
        (6, 6),
    ]

    beta_grid = [0.02, 0.035, 0.05]
    c_alpha_grid = [1.2, 1.6, 2.0]
    c_gamma_grid = [0.8, 1.2, 1.6]
    c_cascade_grid = [1.1, 1.4, 1.8]

    results = []

    for beta_strength, c_alpha, c_gamma, c_cascade in itertools.product(
        beta_grid, c_alpha_grid, c_gamma_grid, c_cascade_grid
    ):
        type_factor = {
            "vertical": 1.00,
            "diagonal": 1.10,
            "cascade": c_cascade,
            "beta_to_alpha": c_alpha,
            "beta_to_gamma": c_gamma,
        }

        DG_local = DG.copy()
        assign_transition_weights(DG_local, beta_strength=beta_strength, type_factor=type_factor)
        _, mean_alpha, mean_gamma = summarize_channel_probabilities(DG_local, start_nodes)

        results.append({
            "beta_strength": beta_strength,
            "c_alpha": c_alpha,
            "c_gamma": c_gamma,
            "c_cascade": c_cascade,
            "mean_alpha": mean_alpha,
            "mean_gamma": mean_gamma,
            "delta": mean_alpha - mean_gamma,
        })

    results.sort(key=lambda r: (r["delta"], r["mean_alpha"]))

    print("Parameterstudie: schwächste alpha-Dominanz zuerst")
    print("-" * 110)
    for r in results[:20]:
        print(
            f"beta={r['beta_strength']:.3f} | "
            f"cA={r['c_alpha']:.2f} | cG={r['c_gamma']:.2f} | cC={r['c_cascade']:.2f} || "
            f"P(alpha)={r['mean_alpha']:.4f} | P(gamma)={r['mean_gamma']:.4f} | "
            f"Delta={r['delta']:.4f}"
        )

    print()
    print("Parameterstudie: stärkste alpha-Dominanz zuerst")
    print("-" * 110)
    for r in results[-20:]:
        print(
            f"beta={r['beta_strength']:.3f} | "
            f"cA={r['c_alpha']:.2f} | cG={r['c_gamma']:.2f} | cC={r['c_cascade']:.2f} || "
            f"P(alpha)={r['mean_alpha']:.4f} | P(gamma)={r['mean_gamma']:.4f} | "
            f"Delta={r['delta']:.4f}"
        )

    print()
    near_balanced = [r for r in results if abs(r["delta"]) < 0.05]
    print(f"Nahezu balancierte Parameterpunkte (|Delta| < 0.05): {len(near_balanced)}")
    for r in near_balanced[:20]:
        print(
            f"beta={r['beta_strength']:.3f} | "
            f"cA={r['c_alpha']:.2f} | cG={r['c_gamma']:.2f} | cC={r['c_cascade']:.2f} || "
            f"P(alpha)={r['mean_alpha']:.4f} | P(gamma)={r['mean_gamma']:.4f}"
        )


if __name__ == "__main__":
    parameter_study()