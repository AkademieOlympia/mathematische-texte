# -*- coding: utf-8 -*-
#
# Primzahlvierlinge -> Simplizialer Komplex -> Homologie
#
# Homologie über die Bibliothek **Dionysus** (persistente Homologie, alle
# Simplizes bei Filtration 0 = gewöhnliche simpliziale Homologie, Koeff. Z2).
#
# Installation:  pip install dionysus

from __future__ import annotations

import math
import sys
from collections import Counter, defaultdict
from itertools import combinations

try:
    import dionysus as dion
except ImportError as e:
    raise ImportError(
        "Homology.py benoetigt 'dionysus'. Installieren mit:\n"
        "  pip install dionysus\n"
    ) from e

# ------------------------------------------------------------
# 0. Sieb + Dionysus-gestützter Komplex
# ------------------------------------------------------------


def _sieve_primes_upto(limit: int) -> set[int]:
    if limit < 2:
        return set()
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i * i : limit + 1 : i] = b"\x00" * (((limit - i * i) // i) + 1)
    return {i for i in range(2, limit + 1) if sieve[i]}


class DionysusSimplicialComplex:
    """
    Simplizes als Tupel von Ecken (beliebige hashbare Labels, z. B. Primzahlen).
    Homologie: Dionysus-Filtration, alle Geburtszeiten 0.
    """

    def __init__(self, simplices: list[tuple]) -> None:
        self._simplices = [tuple(sorted(s)) for s in simplices]
        self._faces: dict[int, set[frozenset]] = defaultdict(set)
        for s in self._simplices:
            self._faces[len(s) - 1].add(frozenset(s))

    def dimension(self) -> int:
        return max(self._faces.keys()) if self._faces else -1

    def faces(self) -> dict[int, set[frozenset]]:
        return {d: set(s) for d, s in self._faces.items()}

    def f_vector(self) -> list[int]:
        if not self._faces:
            return [1]
        d = self.dimension()
        fv = [1]
        for i in range(d + 1):
            fv.append(len(self._faces.get(i, set())))
        return fv

    def betti_numbers(self, max_dim: int) -> dict[int, int]:
        verts = sorted({v for s in self._simplices for v in s})
        vm = {v: i for i, v in enumerate(verts)}
        filt = dion.Filtration()
        for s in self._simplices:
            idx = sorted(vm[v] for v in s)
            filt.append(dion.Simplex(idx, 0.0))
        filt.sort()
        m = dion.homology_persistence(filt)
        dgms = dion.init_diagrams(m, filt)
        beta: dict[int, int] = {}
        for k in range(max_dim + 1):
            if k < len(dgms):
                beta[k] = sum(1 for p in dgms[k] if math.isinf(p.death))
            else:
                beta[k] = 0
        return beta


# ------------------------------------------------------------
# 1. Hilfsfunktionen
# ------------------------------------------------------------


def family_mod12(p: int) -> str:
    if p == 2:
        return "P2"
    if p == 3:
        return "P3"
    r = p % 12
    if r == 1:
        return "E"
    if r == 5:
        return "A"
    if r == 7:
        return "B"
    if r == 11:
        return "C"
    return "?"


def prime_quadruplets(limit: int) -> list[tuple[int, int, int, int]]:
    primes = _sieve_primes_upto(limit)
    quads: list[tuple[int, int, int, int]] = []
    for p in range(2, max(2, limit - 7)):
        q = (p, p + 2, p + 6, p + 8)
        if q[3] <= limit and all(x in primes for x in q):
            quads.append(q)
    return quads


def _iter_k_subsets(quad: tuple, k: int):
    for e in combinations(sorted(quad), k):
        yield tuple(e)


# ------------------------------------------------------------
# 2. Komplex-Konstruktionen
# ------------------------------------------------------------


def build_quadruplet_complex(quadruplets, mode: str = "tetra"):
    simplices: set[tuple] = set()
    vertices: set[int] = set()

    for quad in quadruplets:
        quad = tuple(sorted(quad))
        for v in quad:
            vertices.add(v)
            simplices.add((v,))

        for e in _iter_k_subsets(quad, 2):
            simplices.add(tuple(sorted(e)))

        if mode in ("triangles", "tetra"):
            for f in _iter_k_subsets(quad, 3):
                simplices.add(tuple(sorted(f)))

        if mode == "tetra":
            simplices.add(tuple(sorted(quad)))

    return sorted(simplices, key=lambda s: (len(s), s)), sorted(vertices)


def build_overlap_complex(quadruplets, overlap: int = 1):
    quads = [tuple(q) for q in quadruplets]
    simplices: set[tuple] = {(i,) for i in range(len(quads))}

    for i in range(len(quads)):
        for j in range(i + 1, len(quads)):
            if len(set(quads[i]).intersection(quads[j])) >= overlap:
                simplices.add((i, j))

    return sorted(simplices, key=lambda s: (len(s), s)), quads


# ------------------------------------------------------------
# 3. Homologie (Dionysus)
# ------------------------------------------------------------


def make_simplicial_complex(simplices):
    return DionysusSimplicialComplex(simplices)


def betti_data(sc: DionysusSimplicialComplex, maxdim: int = 3, base_ring=None) -> dict:
    del base_ring  # früher Sage-QQ; Dionysus nutzt Z2
    b = sc.betti_numbers(maxdim)
    data = {k: b.get(k, 0) for k in range(maxdim + 1)}
    data["_engine"] = "dionysus (Z2, Filtration 0)"
    return data


def euler_characteristic_from_betti(betti: dict) -> int:
    chi = 0
    for k, b in betti.items():
        if isinstance(k, str):
            continue
        if isinstance(b, int):
            chi += ((-1) ** k) * b
    return chi


# ------------------------------------------------------------
# 4. Statistische Zusatzanalyse
# ------------------------------------------------------------


def simplex_dimension_counts(sc: DionysusSimplicialComplex) -> dict:
    return {d: len(s) for d, s in sc.faces().items()}


def family_signature_of_quad(quad):
    return tuple(family_mod12(p) for p in quad)


def analyze_quadruplet_families(quadruplets):
    cnt = Counter()
    for q in quadruplets:
        cnt[family_signature_of_quad(q)] += 1
    return cnt


# ------------------------------------------------------------
# 5. Hauptanalyse
# ------------------------------------------------------------


def run_analysis(limit: int = 1000) -> None:
    print("=" * 72)
    print(f"Primzahlvierlinge bis {limit}")
    print("=" * 72)
    print("[Homology] Homologie-Engine: Dionysus (Z2).", file=sys.stderr)

    quads = prime_quadruplets(limit)
    print(f"Anzahl Primzahlvierlinge: {len(quads)}")
    if len(quads) == 0:
        print("Keine Vierlinge gefunden.")
        return

    print("\nErste 10 Vierlinge:")
    for q in quads[:10]:
        print(q, "Familien:", family_signature_of_quad(q))

    fam_counts = analyze_quadruplet_families(quads)
    print("\nHäufigste Familiensignaturen:")
    for sig, n in fam_counts.most_common(10):
        print(sig, "->", n)

    # --------------------------------------------------------
    # A) Graph-Modus
    # --------------------------------------------------------
    print("\n" + "-" * 72)
    print("A) 1-Skelett / Graph-Komplex")
    print("-" * 72)
    simplices_graph, vertices_graph = build_quadruplet_complex(quads, mode="graph")
    SC_graph = make_simplicial_complex(simplices_graph)

    print("Anzahl Vertices:", len(vertices_graph))
    print("Dimension des Komplexes:", SC_graph.dimension())
    print("f-Vektor:", SC_graph.f_vector())
    print("Simplex-Anzahlen:", simplex_dimension_counts(SC_graph))

    betti_graph = betti_data(SC_graph, maxdim=2)
    print("Bettizahlen:", betti_graph)
    print("Euler-Charakteristik:", euler_characteristic_from_betti(betti_graph))

    # --------------------------------------------------------
    # B) 2-Skelett
    # --------------------------------------------------------
    print("\n" + "-" * 72)
    print("B) 2-Skelett / Dreiecks-Komplex")
    print("-" * 72)
    simplices_tri, vertices_tri = build_quadruplet_complex(quads, mode="triangles")
    SC_tri = make_simplicial_complex(simplices_tri)

    print("Anzahl Vertices:", len(vertices_tri))
    print("Dimension des Komplexes:", SC_tri.dimension())
    print("f-Vektor:", SC_tri.f_vector())
    print("Simplex-Anzahlen:", simplex_dimension_counts(SC_tri))

    betti_tri = betti_data(SC_tri, maxdim=3)
    print("Bettizahlen:", betti_tri)
    print("Euler-Charakteristik:", euler_characteristic_from_betti(betti_tri))

    # --------------------------------------------------------
    # C) Voller Tetraeder-Komplex
    # --------------------------------------------------------
    print("\n" + "-" * 72)
    print("C) Voller Tetraeder-Komplex")
    print("-" * 72)
    simplices_tetra, vertices_tetra = build_quadruplet_complex(quads, mode="tetra")
    SC_tetra = make_simplicial_complex(simplices_tetra)

    print("Anzahl Vertices:", len(vertices_tetra))
    print("Dimension des Komplexes:", SC_tetra.dimension())
    print("f-Vektor:", SC_tetra.f_vector())
    print("Simplex-Anzahlen:", simplex_dimension_counts(SC_tetra))

    betti_tetra = betti_data(SC_tetra, maxdim=3)
    print("Bettizahlen:", betti_tetra)
    print("Euler-Charakteristik:", euler_characteristic_from_betti(betti_tetra))

    # --------------------------------------------------------
    # D) Meta-Komplex der Vierlinge
    # --------------------------------------------------------
    print("\n" + "-" * 72)
    print("D) Meta-Komplex der Vierlinge (Überlappung als Relation)")
    print("-" * 72)
    simplices_meta, meta_vertices = build_overlap_complex(quads, overlap=1)
    SC_meta = make_simplicial_complex(simplices_meta)

    print("Anzahl Vierling-Vertices:", len(meta_vertices))
    print("Dimension des Meta-Komplexes:", SC_meta.dimension())
    print("f-Vektor:", SC_meta.f_vector())
    print("Simplex-Anzahlen:", simplex_dimension_counts(SC_meta))

    betti_meta = betti_data(SC_meta, maxdim=2)
    print("Bettizahlen:", betti_meta)
    print("Euler-Charakteristik:", euler_characteristic_from_betti(betti_meta))

    print("\n" + "=" * 72)
    print("Interpretation")
    print("=" * 72)
    print("beta_0 = Anzahl zusammenhängender Komponenten")
    print("beta_1 = Anzahl unabhängiger 1-Zyklen")
    print("beta_2 = Anzahl unabhängiger 2-Hohlräume")
    print()
    print("Graph-Modus zeigt eher rohe Zyklik.")
    print("2-Skelett zeigt, welche Schleifen durch Dreiecke schon gefüllt werden.")
    print("Tetra-Modus zeigt die echte 3D-Topologie der Vierlingsstruktur.")
    print()
    print("Für das Bamberg-Modell ist besonders interessant,")
    print("ob stabile beta_1- oder beta_2-Muster unter wachsendem cutoff auftreten.")
    print()
    print("Hinweis: Bettizahlen stammen aus Dionysus (Koeffizientenkörper Z2).")


# ------------------------------------------------------------
# 6. Start
# ------------------------------------------------------------

if __name__ == "__main__":
    run_analysis(limit=1000)
