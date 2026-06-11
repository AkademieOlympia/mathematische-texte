from __future__ import annotations

from collections import defaultdict
import matplotlib.pyplot as plt
import networkx as nx


# ---------------------------------------------------------
# Beobachtete beta-Schalenleiter
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


# ---------------------------------------------------------
# Knotenmenge
# ---------------------------------------------------------

def make_nodes(ladder: dict[int, list[int]]) -> list[tuple[int, int]]:
    nodes = []
    for rho in sorted(ladder):
        for ell in sorted(set(ladder[rho])):
            nodes.append((rho, ell))
    return nodes


# ---------------------------------------------------------
# Übergangsregeln
# ---------------------------------------------------------

def add_vertical_edges(G: nx.Graph, ladder: dict[int, list[int]]) -> None:
    """
    Gleiche rho-Schale, benachbarte ell-Werte verbinden.
    """
    for rho in sorted(ladder):
        ells = sorted(set(ladder[rho]))
        for i in range(len(ells) - 1):
            a = (rho, ells[i])
            b = (rho, ells[i + 1])
            G.add_edge(a, b, kind="vertical")


def add_diagonal_edges_same_ell(G: nx.Graph, ladder: dict[int, list[int]]) -> None:
    """
    Gleicher ell-Wert über verschiedene rho-Schalen.
    """
    ell_to_rhos = defaultdict(list)
    for rho in sorted(ladder):
        for ell in sorted(set(ladder[rho])):
            ell_to_rhos[ell].append(rho)

    for ell, rhos in ell_to_rhos.items():
        rhos = sorted(rhos)
        for i in range(len(rhos) - 1):
            a = (rhos[i], ell)
            b = (rhos[i + 1], ell)
            G.add_edge(a, b, kind="diagonal")


def add_cascade_edges(G: nx.Graph, ladder: dict[int, list[int]]) -> None:
    """
    Kaskadenregel:
    Wenn rho_1 als ell in einer höheren Schale rho_2 auftaucht,
    verbinde (rho_1, rho_1) mit (rho_2, rho_1), falls beide existieren.
    """
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


# ---------------------------------------------------------
# Graph aufbauen
# ---------------------------------------------------------

def build_transition_graph(ladder: dict[int, list[int]]) -> nx.Graph:
    G = nx.Graph()

    for node in make_nodes(ladder):
        rho, ell = node
        G.add_node(node, rho=rho, ell=ell)

    add_vertical_edges(G, ladder)
    add_diagonal_edges_same_ell(G, ladder)
    add_cascade_edges(G, ladder)

    return G


# ---------------------------------------------------------
# Textausgabe
# ---------------------------------------------------------

def print_graph_summary(G: nx.Graph) -> None:
    print("Übergangsgraph der beta-Schalen")
    print("-" * 50)
    print(f"Anzahl Knoten : {G.number_of_nodes()}")
    print(f"Anzahl Kanten : {G.number_of_edges()}")
    print()

    kind_counter = defaultdict(int)
    for _, _, data in G.edges(data=True):
        kind_counter[data["kind"]] += 1

    print("Kantentypen:")
    for kind, cnt in sorted(kind_counter.items()):
        print(f"  {kind:>10}: {cnt}")
    print()

    print("Knoten:")
    for node in sorted(G.nodes()):
        print(f"  {node}")
    print()

    print("Kanten:")
    for a, b, data in sorted(G.edges(data=True)):
        print(f"  {a} -- {b} [{data['kind']}]")
    print()


# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------

def plot_transition_graph(G: nx.Graph, savepath: str | None = None) -> None:
    pos = {}
    for node, data in G.nodes(data=True):
        rho = data["rho"]
        ell = data["ell"]
        pos[node] = (rho, ell)

    plt.figure(figsize=(10, 6))

    # Knoten
    nx.draw_networkx_nodes(G, pos, node_size=500)

    # Kanten nach Typ
    edge_styles = {
        "vertical": {"style": "solid", "width": 2},
        "diagonal": {"style": "dashed", "width": 2},
        "cascade":  {"style": "dotted", "width": 2},
    }

    for kind, style in edge_styles.items():
        edgelist = [(u, v) for u, v, d in G.edges(data=True) if d["kind"] == kind]
        nx.draw_networkx_edges(
            G,
            pos,
            edgelist=edgelist,
            style=style["style"],
            width=style["width"],
        )

    labels = {node: f"({node[0]},{node[1]})" for node in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)

    plt.xlabel(r"Resonanzzahl $\rho$")
    plt.ylabel(r"Kopplungszahl $\ell$")
    plt.title(r"Übergangsgraph der $\beta$-Schalen")
    plt.grid(True)

    # Legende improvisiert
    plt.plot([], [], linestyle='solid', label='vertikal')
    plt.plot([], [], linestyle='dashed', label='diagonal')
    plt.plot([], [], linestyle='dotted', label='Kaskade')
    plt.legend()

    if savepath:
        plt.savefig(savepath, bbox_inches="tight", dpi=150)
        print(f"Plot gespeichert nach: {savepath}")

    plt.show()


# ---------------------------------------------------------
# Hauptprogramm
# ---------------------------------------------------------

def main() -> None:
    G = build_transition_graph(beta_ladder)
    print_graph_summary(G)
    plot_transition_graph(G, savepath="beta_uebergangsgraph.png")


if __name__ == "__main__":
    main()