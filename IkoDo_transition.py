"""
Vier Matrizen (Ikosaeder/Dodekaeder) und Matrixprodukte zum 60-Ecken-Raum (Halbkanten).

Erzeugt:
  - ikodo_four_matrices.png   (Spy-Plots A_I, A_D, B_I, B_D)
  - ikodo_transition_products.png (H, H B_I^T, Gram H^T H)
"""

from __future__ import annotations

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


def planar_faces(emb: nx.PlanarEmbedding) -> list[list]:
    """Alle begrenzten Facetten als geschlossene Kantenfolgen (Knotenlisten)."""
    faces: list[list] = []
    visited: set[tuple[int, int]] = set()
    for u in emb.nodes():
        for v in emb.neighbors_cw_order(u):
            if (u, v) in visited:
                continue
            face = list(emb.traverse_face(u, v))
            for i in range(len(face)):
                visited.add((face[i], face[(i + 1) % len(face)]))
            faces.append(face)
    return faces


def icosahedron_face_vertex_matrix() -> tuple[np.ndarray, nx.PlanarEmbedding]:
    """Phi in R^{20 x 12}: Phi_{f,v}=1, wenn Ecke v zur Dreiecksflaeche f gehoert."""
    G = nx.icosahedral_graph()
    _, emb = nx.check_planarity(G, True)
    faces = planar_faces(emb)
    if len(faces) != 20 or any(len(f) != 3 for f in faces):
        raise RuntimeError("Ikosaeder: 20 Dreiecksflaechen erwartet.")
    Phi = np.zeros((20, 12), dtype=float)
    for fi, cycle in enumerate(faces):
        for v in set(cycle):
            Phi[fi, v] = 1.0
    return Phi, emb


def edge_vertex_incidence(G: nx.Graph, n: int) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """B in R^{30 x n} mit sortierten Kanten (u,v), u<v."""
    edges = sorted(tuple(sorted(e)) for e in G.edges())
    B = np.zeros((len(edges), n), dtype=float)
    for ei, (u, v) in enumerate(edges):
        B[ei, u] = 1.0
        B[ei, v] = 1.0
    return B, edges


def half_edge_incidence_from_icosa_edges(edges: list[tuple[int, int]]) -> np.ndarray:
    """
    H in R^{60 x 12}: je zwei Zeilen pro Kante (u,v), u<v —
    Zeile 2e traegt 1 an u, Zeile 2e+1 traegt 1 an v (Halbkanten / C60-Ecken).
    """
    H = np.zeros((2 * len(edges), 12), dtype=float)
    for e, (u, v) in enumerate(edges):
        H[2 * e, u] = 1.0
        H[2 * e + 1, v] = 1.0
    return H


def build_four_matrices():
    G_i = nx.icosahedral_graph()
    G_d = nx.dodecahedral_graph()
    A_I = nx.to_numpy_array(G_i, nodelist=range(12), dtype=float)
    A_D = nx.to_numpy_array(G_d, nodelist=range(20), dtype=float)
    B_I, edges_i = edge_vertex_incidence(G_i, 12)
    B_D, edges_d = edge_vertex_incidence(G_d, 20)
    if B_I.shape != (30, 12) or B_D.shape != (30, 20):
        raise RuntimeError("Platonische Graphen: je 30 Kanten erwartet.")
    return A_I, A_D, B_I, B_D, edges_i, edges_d


def verify_gram_identities(A_I: np.ndarray, A_D: np.ndarray, B_I: np.ndarray, B_D: np.ndarray):
    """B^T B = diag(deg) + A fuer einfache Graphen."""
    np.testing.assert_allclose(B_I.T @ B_I, 5 * np.eye(12) + A_I, atol=1e-10)
    np.testing.assert_allclose(B_D.T @ B_D, 3 * np.eye(20) + A_D, atol=1e-10)


def plot_four_and_transition(
    A_I: np.ndarray,
    A_D: np.ndarray,
    B_I: np.ndarray,
    B_D: np.ndarray,
    H: np.ndarray,
    out_four: str = "ikodo_four_matrices.png",
    out_prod: str = "ikodo_transition_products.png",
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    for ax, M, title in zip(
        axes.flat,
        [A_I, A_D, B_I, B_D],
        [
            r"$A_I \in \{0,1\}^{12 \times 12}$ (Ikosaeder)",
            r"$A_D \in \{0,1\}^{20 \times 20}$ (Dodekaeder)",
            r"$B_I \in \{0,1\}^{30 \times 12}$ (Kanten–Ecken, Ikosaeder)",
            r"$B_D \in \{0,1\}^{30 \times 20}$ (Kanten–Ecken, Dodekaeder)",
        ],
    ):
        ax.spy(M, markersize=2, color="black")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Spalte")
        ax.set_ylabel("Zeile")
    plt.tight_layout()
    plt.savefig(out_four, dpi=150)
    plt.close()

    Mhb = H @ B_I.T
    Ggram = H @ H.T
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    axes[0].spy(H, markersize=1, color="C0")
    axes[0].set_title(r"$H \in \{0,1\}^{60 \times 12}$ (Halbkanten $\to$ Ikosaeder-Ecken)")
    axes[1].spy(Mhb, markersize=0.8, color="C1")
    axes[1].set_title(
        r"$H B_I^{\top} \in \{0,1\}^{60 \times 30}$ (Halbkante $\to$ Ikosaeder-Kanten)"
    )
    axes[2].spy(Ggram, markersize=0.5, color="C2")
    axes[2].set_title(r"$H H^{\top}$ (Gram; $12$ Blöcke je $5\times 5$ Einsen)")
    for ax in axes:
        ax.set_xlabel("Spalte")
        ax.set_ylabel("Zeile")
    plt.tight_layout()
    plt.savefig(out_prod, dpi=150)
    plt.close()


def main() -> None:
    A_I, A_D, B_I, B_D, edges_i, _ = build_four_matrices()
    verify_gram_identities(A_I, A_D, B_I, B_D)
    Phi, _ = icosahedron_face_vertex_matrix()
    H = half_edge_incidence_from_icosa_edges(edges_i)

    np.testing.assert_allclose(H.T @ H, 5 * np.eye(12), atol=1e-10)
    # Zwei Dreiecksflaechen pro Kante des Ikosaeders
    np.testing.assert_allclose(Phi.T @ Phi, 5 * np.eye(12) + 2 * A_I, atol=1e-10)

    plot_four_and_transition(A_I, A_D, B_I, B_D, H)
    print("Gespeichert: ikodo_four_matrices.png, ikodo_transition_products.png")
    print("Identitaeten B_I^T B_I = 5I + A_I, B_D^T B_D = 3I + A_D, H^T H = 5I verifiziert.")
    print("Phi^T Phi = 5I + 2 A_I (Flaechen-Inzidenz Ikosaeder) verifiziert.")


if __name__ == "__main__":
    main()
