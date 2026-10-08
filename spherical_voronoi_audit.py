#!/usr/bin/env python3
"""Sphaerische Voronoi-Zellen auf S^2 ohne SciPy.

Die Zellenecken der Zelle um den Samen s_i entstehen als Schnittpunkte
zweier Bisektorebenen

    x . (s_i - s_j) = 0,   x . (s_i - s_k) = 0,   ||x|| = 1.

Ihre beiden Kandidaten sind die normierten Vektoren

    +/- ((s_i-s_j) x (s_i-s_k)).

Ein Kandidat gehoert genau dann zur Zelle i, wenn

    x . s_i >= x . s_l

fuer alle Samen s_l gilt. Das Skript verwendet die 12 normierten Ecken
eines Ikosaeders. Die duale Voronoi-Zerlegung ist ein Dodekaeder mit
20 Ecken, 30 Kanten und 12 fuenfeckigen Zellen.

Abhaengigkeiten: numpy, matplotlib
"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


TOL = 1.0e-9


def normalize(v: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(v)
    if norm <= TOL:
        raise ValueError("Nullvektor kann nicht normiert werden.")
    return v / norm


def icosahedron_seeds() -> np.ndarray:
    """Gibt die 12 Einheitsrichtungen eines regulaeren Ikosaeders zurueck."""
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    seeds = []
    for a in (-1.0, 1.0):
        for b in (-phi, phi):
            seeds.extend(((0.0, a, b), (a, b, 0.0), (b, 0.0, a)))
    return np.array([normalize(np.asarray(v, dtype=float)) for v in seeds])


def deduplicate(points: list[np.ndarray], tol: float = TOL) -> list[np.ndarray]:
    unique: list[np.ndarray] = []
    for point in points:
        if not any(np.linalg.norm(point - old) <= tol for old in unique):
            unique.append(point)
    return unique


def cell_vertices(seeds: np.ndarray, i: int, tol: float = TOL) -> list[np.ndarray]:
    """Berechnet alle Ecken der sphaerischen Voronoi-Zelle von seeds[i]."""
    candidates: list[np.ndarray] = []
    other = [j for j in range(len(seeds)) if j != i]

    for j, k in combinations(other, 2):
        normal_ij = seeds[i] - seeds[j]
        normal_ik = seeds[i] - seeds[k]
        cross = np.cross(normal_ij, normal_ik)
        cross_norm = np.linalg.norm(cross)
        if cross_norm <= tol:
            continue

        direction = cross / cross_norm
        for vertex in (direction, -direction):
            dominance = seeds @ vertex
            if np.all(dominance[i] >= dominance - tol):
                candidates.append(vertex)

    return order_around_site(deduplicate(candidates, 50.0 * tol), seeds[i])


def order_around_site(vertices: list[np.ndarray], site: np.ndarray) -> list[np.ndarray]:
    """Ordnet Zellenecken zyklisch in der Tangentialebene des Samens."""
    if len(vertices) < 3:
        return vertices

    reference = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(reference, site)) > 0.9:
        reference = np.array([1.0, 0.0, 0.0])
    tangent_1 = normalize(reference - np.dot(reference, site) * site)
    tangent_2 = np.cross(site, tangent_1)

    def angle(v: np.ndarray) -> float:
        projected = v - np.dot(v, site) * site
        return float(np.arctan2(np.dot(projected, tangent_2),
                                np.dot(projected, tangent_1)))

    return sorted(vertices, key=angle)


def global_vertex_index(cells: list[list[np.ndarray]], tol: float = TOL):
    vertices: list[np.ndarray] = []
    indexed_cells: list[list[int]] = []
    for cell in cells:
        indices = []
        for point in cell:
            match = next((idx for idx, old in enumerate(vertices)
                          if np.linalg.norm(point - old) <= 50.0 * tol), None)
            if match is None:
                vertices.append(point)
                match = len(vertices) - 1
            indices.append(match)
        indexed_cells.append(indices)
    return np.asarray(vertices), indexed_cells


def cell_edges(indexed_cells: list[list[int]]) -> set[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    for cell in indexed_cells:
        for left, right in zip(cell, cell[1:] + cell[:1]):
            edges.add(tuple(sorted((left, right))))
    return edges


def slerp(a: np.ndarray, b: np.ndarray, samples: int = 40) -> np.ndarray:
    """Kurzer Grosskreisbogen zwischen zwei nicht-antipodalen Punkten."""
    omega = np.arccos(np.clip(np.dot(a, b), -1.0, 1.0))
    if omega <= TOL:
        return np.repeat(a[None, :], samples, axis=0)
    t = np.linspace(0.0, 1.0, samples)
    return (np.sin((1.0 - t) * omega)[:, None] * a
            + np.sin(t * omega)[:, None] * b) / np.sin(omega)


def audit(seeds: np.ndarray, cells: list[list[np.ndarray]],
          vertices: np.ndarray, indexed_cells: list[list[int]]) -> None:
    edges = cell_edges(indexed_cells)
    unit_error = float(np.max(np.abs(np.linalg.norm(vertices, axis=1) - 1.0)))
    dominance_error = 0.0
    bisector_error = 0.0

    for i, cell in enumerate(cells):
        for vertex in cell:
            scores = seeds @ vertex
            dominance_error = max(dominance_error,
                                  float(np.max(scores - scores[i])))
            tied = np.sort(np.abs(scores - scores[i]))
            # Neben dem Samen i muessen mindestens zwei weitere Samen
            # am Voronoi-Eckpunkt denselben Abstand besitzen.
            bisector_error = max(bisector_error, float(tied[2]))

    face_count = len(cells)
    edge_count = len(edges)
    vertex_count = len(vertices)
    euler = vertex_count - edge_count + face_count

    assert len(seeds) == 12
    assert all(len(cell) == 5 for cell in cells)
    assert vertex_count == 20
    assert edge_count == 30
    assert euler == 2
    assert unit_error < 1.0e-12
    assert dominance_error < 1.0e-9
    assert bisector_error < 1.0e-9

    print("SPHAERISCHES VORONOI-AUDIT")
    print("Status                       : BESTANDEN")
    print(f"Samen / Zellen               : {face_count}")
    print(f"Globale Zellenecken          : {vertex_count}")
    print(f"Globale Kanten               : {edge_count}")
    print(f"Euler-Charakteristik V-E+F   : {euler}")
    print(f"Zellenecken je Zelle         : {[len(c) for c in cells]}")
    print(f"Maximaler Einheitskugelfehler: {unit_error:.3e}")
    print(f"Maximale Halbraumverletzung  : {dominance_error:.3e}")
    print(f"Maximaler Bisektorfehler     : {bisector_error:.3e}")
    print("Interpretation               : endliches numerisches Zertifikat")


def plot_diagram(seeds: np.ndarray, vertices: np.ndarray,
                 indexed_cells: list[list[int]], output: Path) -> None:
    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_subplot(111, projection="3d")

    u = np.linspace(0.0, 2.0 * np.pi, 80)
    v = np.linspace(0.0, np.pi, 40)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(x, y, z, rstride=5, cstride=5,
                      color="#cbd5e1", linewidth=0.35, alpha=0.35)

    for left, right in sorted(cell_edges(indexed_cells)):
        arc = slerp(vertices[left], vertices[right])
        ax.plot(arc[:, 0], arc[:, 1], arc[:, 2],
                color="#155e75", linewidth=1.8)

    ax.scatter(seeds[:, 0], seeds[:, 1], seeds[:, 2],
               color="#dc2626", s=45, label="Ikosaeder-Samen")
    ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2],
               color="#f59e0b", s=28, label="Voronoi-Ecken")
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlim(-1.08, 1.08)
    ax.set_ylim(-1.08, 1.08)
    ax.set_zlim(-1.08, 1.08)
    ax.set_title("Sphaerische Voronoi-Zerlegung auf $S^2$")
    ax.legend(loc="upper left")
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Berechnet und auditiert eine sphaerische Voronoi-Zerlegung."
    )
    parser.add_argument("--output", type=Path,
                        default=Path("spherical_voronoi_audit.png"),
                        help="Ausgabedatei der Visualisierung")
    parser.add_argument("--no-plot", action="store_true",
                        help="Nur Audit ausfuehren, keine Grafik schreiben")
    args = parser.parse_args()

    seeds = icosahedron_seeds()
    cells = [cell_vertices(seeds, i) for i in range(len(seeds))]
    vertices, indexed_cells = global_vertex_index(cells)
    audit(seeds, cells, vertices, indexed_cells)

    if not args.no_plot:
        plot_diagram(seeds, vertices, indexed_cells, args.output)
        print(f"Grafik                       : {args.output.resolve()}")


if __name__ == "__main__":
    main()
