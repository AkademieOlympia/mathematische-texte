from __future__ import annotations

from collections import defaultdict
import math
import numpy as np
import matplotlib.pyplot as plt
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
# Energie und gerichteter Graph
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
# Gewichtung
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
        data["raw_weight"] = raw

    for u in DG.nodes():
        outs = list(DG.out_edges(u, data=True))
        if not outs:
            continue
        total = sum(d["raw_weight"] for _, _, d in outs)
        for _, _, d in outs:
            d["prob"] = d["raw_weight"] / total


# ---------------------------------------------------------
# Absorption
# ---------------------------------------------------------

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


def mean_channel_bias(DG, start_nodes):
    sinks, probs = absorption_probabilities(DG)
    sink_families = {s: DG.nodes[s]["family"] for s in sinks}

    p_alpha_list = []
    p_gamma_list = []

    for u in start_nodes:
        vec = probs[u]
        p_alpha = sum(vec[i] for i, s in enumerate(sinks) if sink_families[s] == "alpha")
        p_gamma = sum(vec[i] for i, s in enumerate(sinks) if sink_families[s] == "gamma")
        p_alpha_list.append(p_alpha)
        p_gamma_list.append(p_gamma)

    mean_alpha = float(np.mean(p_alpha_list))
    mean_gamma = float(np.mean(p_gamma_list))
    return mean_alpha, mean_gamma, mean_alpha - mean_gamma


# ---------------------------------------------------------
# Phasenkarte
# ---------------------------------------------------------

def phase_map(beta_strength=0.035, c_cascade=1.40, resolution=41):
    G = build_combined_graph()
    DG_base = build_directed_graph(G, model="additive", lam=1.0)

    start_nodes = [
        (198, 102),
        (192, 108),
        (102, 102),
        (18, 18),
        (12, 12),
        (6, 6),
    ]

    cA_vals = np.linspace(0.8, 2.2, resolution)
    cG_vals = np.linspace(0.8, 2.2, resolution)

    Delta = np.zeros((resolution, resolution))
    Alpha = np.zeros((resolution, resolution))
    Gamma = np.zeros((resolution, resolution))

    for i, cG in enumerate(cG_vals):
        for j, cA in enumerate(cA_vals):
            DG = DG_base.copy()

            type_factor = {
                "vertical": 1.00,
                "diagonal": 1.10,
                "cascade": c_cascade,
                "beta_to_alpha": float(cA),
                "beta_to_gamma": float(cG),
            }

            assign_transition_weights(DG, beta_strength=beta_strength, type_factor=type_factor)
            mean_alpha, mean_gamma, delta = mean_channel_bias(DG, start_nodes)

            Delta[i, j] = delta
            Alpha[i, j] = mean_alpha
            Gamma[i, j] = mean_gamma

    return cA_vals, cG_vals, Delta, Alpha, Gamma


def plot_phase_map(cA_vals, cG_vals, Delta, savepath=None):
    plt.figure(figsize=(8, 6))

    extent = [cA_vals.min(), cA_vals.max(), cG_vals.min(), cG_vals.max()]
    im = plt.imshow(
        Delta,
        origin="lower",
        extent=extent,
        aspect="auto",
        interpolation="nearest",
    )

    plt.colorbar(im, label=r"$\Delta = P(\alpha)-P(\gamma)$")
    plt.contour(
        cA_vals,
        cG_vals,
        Delta,
        levels=[0.0],
        linewidths=2,
    )

    plt.xlabel(r"$c_{\beta\to\alpha}$")
    plt.ylabel(r"$c_{\beta\to\gamma}$")
    plt.title(r"Phasenkarte der bifurkativen Relaxation")

    if savepath:
        plt.savefig(savepath, bbox_inches="tight", dpi=160)
        print(f"Phasenkarte gespeichert nach: {savepath}")

    plt.show()


def main():
    beta_strength = 0.035
    c_cascade = 1.40
    resolution = 41

    cA_vals, cG_vals, Delta, Alpha, Gamma = phase_map(
        beta_strength=beta_strength,
        c_cascade=c_cascade,
        resolution=resolution,
    )

    # einige Referenzpunkte ausgeben
    test_points = [
        (1.2, 1.6),
        (1.6, 1.2),
        (1.6, 1.6),
        (1.2, 1.2),
    ]

    print("Referenzpunkte:")
    for cA, cG in test_points:
        j = np.argmin(np.abs(cA_vals - cA))
        i = np.argmin(np.abs(cG_vals - cG))
        print(
            f"cA={cA:.2f}, cG={cG:.2f} | "
            f"P(alpha)={Alpha[i,j]:.4f} | P(gamma)={Gamma[i,j]:.4f} | Delta={Delta[i,j]:.4f}"
        )

    plot_phase_map(cA_vals, cG_vals, Delta, savepath="phasenkarte_abg.png")


if __name__ == "__main__":
    main()