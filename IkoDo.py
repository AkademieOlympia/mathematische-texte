import networkx as nx
import numpy as np


def _even_cyclic_perms(a, b, c):
    return np.array([[a, b, c], [b, c, a], [c, a, b]], dtype=float)


def truncated_icosahedron_graph():
    """Platonischer Fußball / C60-Skelett: 60 Ecken, 3-regulär, 90 Kanten."""
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    raw = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            raw.extend(_even_cyclic_perms(0.0, float(s1), float(s2) * 3.0 * phi))
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            for s3 in (-1, 1):
                raw.extend(
                    _even_cyclic_perms(
                        float(s1), float(s2) * (2.0 + phi), float(s3) * 2.0 * phi
                    )
                )
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            for s3 in (-1, 1):
                raw.extend(
                    _even_cyclic_perms(
                        float(s1) * phi,
                        float(s2) * 2.0,
                        float(s3) * (2.0 * phi + 1.0),
                    )
                )
    V = np.unique(np.round(np.array(raw), decimals=10), axis=0)
    if V.shape[0] != 60:
        raise RuntimeError(f"erwartet 60 Ecken, erhalten {V.shape[0]}")

    diff = V[:, np.newaxis, :] - V[np.newaxis, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=-1))
    d_min = np.min(dist[np.triu_indices(60, k=1)])
    tol = 1e-5 * max(d_min, 1.0)
    G = nx.Graph()
    G.add_nodes_from(range(60))
    for i in range(60):
        for j in range(i + 1, 60):
            if dist[i, j] <= d_min + tol:
                G.add_edge(i, j)
    if G.number_of_edges() != 90:
        raise RuntimeError(f"erwartet 90 Kanten, erhalten {G.number_of_edges()}")
    return G


def tight_binding_spectrum(G, epsilon=0, t=1, decimals=5):
    """H = epsilon*I - t*A, sortierte reelle Eigenwerte."""
    A = nx.to_numpy_array(G, dtype=float)
    n = G.number_of_nodes()
    H = epsilon * np.eye(n) - t * A
    energies = np.linalg.eigvalsh(H)
    return sorted(np.round(energies, decimals).tolist())


# Parameter (Beispielwerte)
epsilon = 0
t = 1

# Dodekaeder (20 Knoten, 3-regulär)
print("Dodekaeder:", tight_binding_spectrum(nx.dodecahedral_graph(), epsilon, t))

# Ikosaeder (12 Knoten, 5-regulär) — dual zum Dodekaeder
print("Ikosaeder:", tight_binding_spectrum(nx.icosahedral_graph(), epsilon, t))

# Abgestumpftes Ikosaeder / C60-Gitter (60 Knoten, 3-regulär)
print(
    "C60 (abgestumpftes Ikosaeder, 60 Ecken):",
    tight_binding_spectrum(truncated_icosahedron_graph(), epsilon, t),
)
