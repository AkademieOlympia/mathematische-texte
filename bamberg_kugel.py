#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
bamberg_kugel.py

Parametrisierbare SageMath-10.5/Python-Klasse zur Erzeugung von
Bamberg-Kugel-Normschalen auf Basis des Icosianen-/2I-Modells.

WICHTIG:
- Die Klasse ist für SageMath 10.5 ausgelegt.
- Die vollständig zertifizierte BK-I1/BK-I2-Struktur bezieht sich auf die
  120 Einheitsicosianen (Q_I = 1).
- Für beliebige Normschalen wird eine allgemeine Shell-Engine bereitgestellt.
  Die konkrete mathematische Definition der "BK-Kanten" höherer Schalen muss
  separat parametrisiert bzw. erforscht werden.

Aufruf:
    sage -python bamberg_kugel.py

oder importierbar:
    from bamberg_kugel import BambergKugel
"""

from dataclasses import dataclass, asdict
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

try:
    from sage.all import (
        QQ, ZZ, PolynomialRing, NumberField, QuaternionAlgebra,
        Graph, matrix, vector
    )
    SAGE_AVAILABLE = True
except Exception:
    SAGE_AVAILABLE = False


@dataclass
class BambergKugelConfig:
    norm_shell: Any = 1
    edge_norm: Any = 2
    use_antipodal_quotient: bool = False
    require_connected: bool = True
    require_regular: bool = False
    expected_degree: Optional[int] = None
    max_search_coefficient: int = 2
    include_zero: bool = False
    deduplicate_signs: bool = False
    exact_arithmetic: bool = True


class BambergKugel:
    """
    Parametrisierbare Engine für Bamberg-Kugel-Normschalen.

    Kernbegriffe:
        N_H(q) = reduzierte Quaternionennorm in Q(sqrt(5))
        Q_I(q) = x + y für N_H(q) = x + sqrt(5)*y

    Hauptmodi:
    1. unit_shell():
       Liefert die 120 binär-ikosaedrischen Einheiten (Q_I = 1),
       sofern der klassische 2I-Generatoraufbau verwendet wird.

    2. generate_shell():
       Enumeriert Icosianen aus einem endlich gewählten Z-Modulfenster
       und filtert nach Q_I(q) = norm_shell.

    3. build_distance_graph():
       Baut den vollständigen Distanzgraphen innerhalb einer Schale
       über Q_I(r-q) = edge_norm.

    4. build_bk_graph():
       Verwendet entweder
       - eine benutzerdefinierte edge_selector-Funktion, oder
       - für norm_shell == 1 optional eine spezielle, extern implementierbare
         BK-I1-Selektionsregel.

    Die Klasse ist absichtlich so gebaut, dass höhere Normschalen nicht
    fälschlich als bereits mathematisch zertifizierte BK-Strukturen ausgegeben
    werden.
    """

    def __init__(self, config: Optional[BambergKugelConfig] = None):
        self._config = config or BambergKugelConfig()

        self._K = None
        self._sqrt5 = None
        self._H = None
        self._one = None
        self._minus_one = None

        self._group_2I: Optional[Set[Any]] = None
        self._icosian_basis: Optional[List[Any]] = None
        self._shell_points: Optional[List[Any]] = None
        self._distance_graph = None
        self._bk_graph = None

        self._edge_selector: Optional[Callable[[Any, Any, "BambergKugel"], bool]] = None
        self._metadata: Dict[str, Any] = {}

        if SAGE_AVAILABLE:
            self._initialize_algebra()

    # ------------------------------------------------------------------
    # Getter / Setter: Konfiguration
    # ------------------------------------------------------------------

    def get_config(self) -> BambergKugelConfig:
        return self._config

    def set_config(self, config: BambergKugelConfig) -> None:
        self._config = config
        self._invalidate_computed_state()

    def get_norm_shell(self):
        return self._config.norm_shell

    def set_norm_shell(self, value) -> None:
        self._config.norm_shell = value
        self._invalidate_shell_state()

    def get_edge_norm(self):
        return self._config.edge_norm

    def set_edge_norm(self, value) -> None:
        self._config.edge_norm = value
        self._distance_graph = None
        self._bk_graph = None

    def get_max_search_coefficient(self) -> int:
        return self._config.max_search_coefficient

    def set_max_search_coefficient(self, value: int) -> None:
        if int(value) < 0:
            raise ValueError("max_search_coefficient muss >= 0 sein.")
        self._config.max_search_coefficient = int(value)
        self._invalidate_shell_state()

    def get_expected_degree(self) -> Optional[int]:
        return self._config.expected_degree

    def set_expected_degree(self, value: Optional[int]) -> None:
        if value is not None and int(value) < 0:
            raise ValueError("expected_degree muss >= 0 oder None sein.")
        self._config.expected_degree = None if value is None else int(value)

    def get_require_connected(self) -> bool:
        return self._config.require_connected

    def set_require_connected(self, value: bool) -> None:
        self._config.require_connected = bool(value)

    def get_require_regular(self) -> bool:
        return self._config.require_regular

    def set_require_regular(self, value: bool) -> None:
        self._config.require_regular = bool(value)

    def get_use_antipodal_quotient(self) -> bool:
        return self._config.use_antipodal_quotient

    def set_use_antipodal_quotient(self, value: bool) -> None:
        self._config.use_antipodal_quotient = bool(value)
        self._invalidate_shell_state()

    def get_deduplicate_signs(self) -> bool:
        return self._config.deduplicate_signs

    def set_deduplicate_signs(self, value: bool) -> None:
        self._config.deduplicate_signs = bool(value)
        self._invalidate_shell_state()

    def get_edge_selector(self):
        return self._edge_selector

    def set_edge_selector(
        self,
        selector: Optional[Callable[[Any, Any, "BambergKugel"], bool]]
    ) -> None:
        """
        selector(q, r, self) -> bool

        Wird für build_bk_graph() benutzt.
        """
        self._edge_selector = selector
        self._bk_graph = None

    # ------------------------------------------------------------------
    # Getter: mathematische Objekte
    # ------------------------------------------------------------------

    def get_number_field(self):
        self._require_sage()
        return self._K

    def get_quaternion_algebra(self):
        self._require_sage()
        return self._H

    def get_group_2I(self) -> Set[Any]:
        self._require_sage()
        if self._group_2I is None:
            self._group_2I = self._build_binary_icosahedral_group()
        return set(self._group_2I)

    def get_icosian_basis(self) -> List[Any]:
        self._require_sage()
        if self._icosian_basis is None:
            self._icosian_basis = self._default_icosian_basis()
        return list(self._icosian_basis)

    def set_icosian_basis(self, basis: Sequence[Any]) -> None:
        self._require_sage()
        if len(basis) != 8:
            raise ValueError("Die Icosianenbasis muss als Z-Modulbasis Rang 8 haben.")
        self._icosian_basis = list(basis)
        self._invalidate_shell_state()

    def get_shell_points(self) -> List[Any]:
        if self._shell_points is None:
            return self.generate_shell()
        return list(self._shell_points)

    def set_shell_points(self, points: Sequence[Any]) -> None:
        self._shell_points = list(points)
        self._distance_graph = None
        self._bk_graph = None

    def get_distance_graph(self):
        return self._distance_graph

    def get_bk_graph(self):
        return self._bk_graph

    def get_metadata(self) -> Dict[str, Any]:
        return dict(self._metadata)

    # ------------------------------------------------------------------
    # Algebra
    # ------------------------------------------------------------------

    def _require_sage(self) -> None:
        if not SAGE_AVAILABLE:
            raise RuntimeError(
                "SageMath ist nicht verfügbar. Bitte mit "
                "`sage -python bamberg_kugel.py` oder innerhalb SageMath 10.5 starten."
            )

    def _initialize_algebra(self) -> None:
        R = PolynomialRing(QQ, 't')
        t = R.gen()
        self._K = NumberField(t**2 - 5, 'sqrt5')
        self._sqrt5 = self._K.gen()
        self._H = QuaternionAlgebra(self._K, -1, -1)
        self._one = self._H(1)
        self._minus_one = self._H(-1)

    def _coeffs_quaternion(self, q) -> Tuple[Any, Any, Any, Any]:
        """
        Robust quaternion coefficient extraction.
        """
        try:
            return tuple(q.coefficient_tuple())
        except Exception:
            try:
                return tuple(q)
            except Exception as exc:
                raise TypeError(f"Kann Quaternionenkoeffizienten nicht lesen: {q}") from exc

    def reduced_norm(self, q):
        self._require_sage()
        try:
            return q.reduced_norm()
        except Exception:
            a, b, c, d = self._coeffs_quaternion(q)
            return a*a + b*b + c*c + d*d

    def split_Qsqrt5(self, z) -> Tuple[Any, Any]:
        """
        Liefert x,y mit z = x + y*sqrt(5), x,y in QQ.
        """
        self._require_sage()
        z = self._K(z)
        coeffs = list(z.list())
        if len(coeffs) == 0:
            return QQ(0), QQ(0)
        if len(coeffs) == 1:
            return QQ(coeffs[0]), QQ(0)
        return QQ(coeffs[0]), QQ(coeffs[1])

    def icosian_norm(self, q):
        """
        Q_I(q) = x + y, wenn N_H(q) = x + sqrt(5)*y.
        """
        n = self.reduced_norm(q)
        x, y = self.split_Qsqrt5(n)
        return x + y

    # ------------------------------------------------------------------
    # 2I / Einheitsicosianen
    # ------------------------------------------------------------------

    def _hurwitz_units(self) -> Set[Any]:
        """
        24 Hurwitz-Einheiten im Quaternionenraum.
        """
        self._require_sage()
        H = self._H
        K = self._K

        elems = set()

        # ±1, ±i, ±j, ±k
        basis = [
            H(1),
            H.gen(0),
            H.gen(1),
            H.gen(0) * H.gen(1),
        ]
        for u in basis:
            elems.add(u)
            elems.add(-u)

        # (±1 ± i ± j ± k)/2
        one = H(1)
        i = H.gen(0)
        j = H.gen(1)
        k = i*j
        half = K(1)/2

        for s0 in (-1, 1):
            for s1 in (-1, 1):
                for s2 in (-1, 1):
                    for s3 in (-1, 1):
                        elems.add(
                            half * (s0*one + s1*i + s2*j + s3*k)
                        )

        assert len(elems) == 24
        return elems

    def _icosahedral_seed(self):
        """
        Ein Standard-Seed mit Goldener-Zahl-Koeffizienten.
        Der Gruppenabschluss mit den Hurwitz-Einheiten erzeugt 2I.
        """
        self._require_sage()
        K = self._K
        H = self._H
        sqrt5 = self._sqrt5
        phi = (K(1) + sqrt5) / 2
        invphi = 1 / phi

        one = H(1)
        i = H.gen(0)
        j = H.gen(1)
        k = i*j

        # Standard H4/600-cell-artiges Seed.
        return (phi*one + i + invphi*j) / 2

    def _group_closure(self, generators: Iterable[Any]) -> Set[Any]:
        gens = list(generators)
        gens = gens + [g.inverse() for g in gens]

        G = {self._one}
        changed = True

        while changed:
            changed = False
            old = list(G)
            for x in old:
                for g in gens:
                    z = x * g
                    if z not in G:
                        G.add(z)
                        changed = True

        return G

    def _build_binary_icosahedral_group(self) -> Set[Any]:
        hurwitz = self._hurwitz_units()
        seed = self._icosahedral_seed()

        generators = list(hurwitz) + [seed, -seed]
        G = self._group_closure(generators)

        if len(G) != 120:
            raise AssertionError(
                f"2I-Konstruktion ergab {len(G)} statt 120 Elemente. "
                "Seed/Basis-Konvention prüfen."
            )

        # Exakte Gruppenchecks
        assert self._one in G
        assert self._minus_one in G
        assert all(g.inverse() in G for g in G)
        assert all((a*b) in G for a in G for b in G)

        return G

    def unit_shell(self) -> List[Any]:
        """
        Die zertifizierbare Q_I=1-Schale: 120 Einheiten von 2I.
        """
        G = sorted(self.get_group_2I(), key=self.exact_key)
        assert all(self.icosian_norm(g) == 1 for g in G)
        return G

    # ------------------------------------------------------------------
    # Icosianenbasis und allgemeine Normschalen
    # ------------------------------------------------------------------

    def _default_icosian_basis(self) -> List[Any]:
        """
        Liefert eine Z-Modulbasis aus 8 linear unabhängigen Icosianen.

        Die Basis wird aus 2I-Elementen automatisch gewählt, indem deren
        Q-Koordinaten in der 8-dimensionalen Basis
            1, sqrt5
        pro Quaternionenkoordinate auf Rang 8 geprüft werden.

        Diese Wahl ist deterministisch nach exact_key.
        """
        self._require_sage()
        candidates = sorted(self.get_group_2I(), key=self.exact_key)

        chosen = []
        rows = []

        for q in candidates:
            row = self.q_vector_8(q)
            trial = rows + [row]
            M = matrix(QQ, trial)
            if M.rank() > len(rows):
                chosen.append(q)
                rows.append(row)
            if len(chosen) == 8:
                break

        if len(chosen) != 8:
            raise AssertionError("Konnte keine Rang-8-Icosianenbasis bestimmen.")

        assert matrix(QQ, rows).rank() == 8
        return chosen

    def q_vector_8(self, q):
        """
        Quaternion q -> 8-dim QQ-Vektor:
        pro Koeffizient a = x + y*sqrt5 werden x,y abgelegt.
        """
        self._require_sage()
        out = []
        for c in self._coeffs_quaternion(q):
            x, y = self.split_Qsqrt5(c)
            out.extend([x, y])
        return vector(QQ, out)

    def exact_key(self, q):
        """
        Deterministischer Sortierschlüssel für Quaternionen.
        """
        v = self.q_vector_8(q)
        return tuple(v)

    def _canonical_sign(self, q):
        """
        Wählt deterministisch einen Repräsentanten aus {q,-q}.
        """
        return min((q, -q), key=self.exact_key)

    def generate_shell(self, norm_shell=None, coefficient_bound=None) -> List[Any]:
        """
        Enumeriert Punkte der gewählten Normschale in einem endlichen
        Koeffizientenfenster der 8-dimensionalen Icosianenbasis.

        ACHTUNG:
        Für beliebige Normschalen ist dies eine ENDLICHE SUCHENGINE.
        Vollständigkeit ist nur garantiert, wenn coefficient_bound groß
        genug gewählt und mathematisch separat begründet wird.

        Für norm_shell == 1 wird direkt die exakte 120er Einheitsschale benutzt.
        """
        self._require_sage()

        shell = self._config.norm_shell if norm_shell is None else norm_shell
        bound = (
            self._config.max_search_coefficient
            if coefficient_bound is None
            else int(coefficient_bound)
        )

        if shell == 1:
            pts = self.unit_shell()
            self._shell_points = pts
            self._metadata["shell_complete"] = True
            self._metadata["shell_generation_mode"] = "exact_2I_unit_shell"
            return list(pts)

        basis = self.get_icosian_basis()

        # Für Rang 8 wächst (2B+1)^8 sehr schnell.
        count_candidates = (2*bound + 1)**8
        if count_candidates > 5_000_000:
            raise ValueError(
                f"Suchraum zu groß: {(2*bound + 1)}^8 = {count_candidates}. "
                "Bitte coefficient_bound verkleinern oder eine spezialisierte "
                "Shell-Enumeration implementieren."
            )

        points = set()

        # Iterative kartesische Produkterzeugung ohne itertools.product-Materialisierung
        def rec(idx, current):
            if idx == 8:
                if not self._config.include_zero and current == 0:
                    return
                if self.icosian_norm(current) == shell:
                    q = current
                    if self._config.deduplicate_signs or self._config.use_antipodal_quotient:
                        q = self._canonical_sign(q)
                    points.add(q)
                return

            for c in range(-bound, bound + 1):
                rec(idx + 1, current + c*basis[idx])

        rec(0, self._H(0))

        pts = sorted(points, key=self.exact_key)
        self._shell_points = pts
        self._distance_graph = None
        self._bk_graph = None

        self._metadata["shell_complete"] = False
        self._metadata["shell_generation_mode"] = "bounded_icosian_basis_enumeration"
        self._metadata["coefficient_bound"] = bound
        self._metadata["candidate_count"] = count_candidates

        return list(pts)

    # ------------------------------------------------------------------
    # Graphen
    # ------------------------------------------------------------------

    def build_distance_graph(
        self,
        points: Optional[Sequence[Any]] = None,
        edge_norm=None
    ):
        """
        Vollständiger Distanzgraph innerhalb der Schale:
            q -- r  <=>  Q_I(r-q) == edge_norm
        """
        self._require_sage()
        pts = list(points) if points is not None else self.get_shell_points()
        en = self._config.edge_norm if edge_norm is None else edge_norm

        G = Graph()
        G.add_vertices(range(len(pts)))

        for i in range(len(pts)):
            qi = pts[i]
            for j in range(i + 1, len(pts)):
                qj = pts[j]
                if self.icosian_norm(qj - qi) == en:
                    G.add_edge(i, j)

        self._distance_graph = G
        self._metadata["distance_graph_vertices"] = G.order()
        self._metadata["distance_graph_edges"] = G.size()
        self._metadata["distance_graph_edge_norm"] = en

        self._validate_graph(G)
        return G

    def build_bk_graph(
        self,
        points: Optional[Sequence[Any]] = None,
        selector: Optional[Callable[[Any, Any, "BambergKugel"], bool]] = None
    ):
        """
        Baut einen parametrisierten BK-Teilgraphen.

        Für höhere Schalen wird KEINE unbegründete Standardselektion angenommen.
        Ein selector muss spezifiziert werden:

            selector(q, r, bk) -> bool

        Beispiel:
            lambda q, r, bk: bk.icosian_norm(r-q) == bk.get_edge_norm()

        Das wäre dann allerdings nur der vollständige Distanzgraph,
        nicht die spezielle kubische BK-I1-Selektion.
        """
        self._require_sage()
        pts = list(points) if points is not None else self.get_shell_points()
        sel = selector or self._edge_selector

        if sel is None:
            raise ValueError(
                "Kein BK-edge_selector gesetzt. Für beliebige Normschalen "
                "muss die Kantenregel explizit angegeben werden."
            )

        G = Graph()
        G.add_vertices(range(len(pts)))

        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                if sel(pts[i], pts[j], self):
                    G.add_edge(i, j)

        self._bk_graph = G
        self._metadata["bk_graph_vertices"] = G.order()
        self._metadata["bk_graph_edges"] = G.size()

        self._validate_graph(G)
        return G

    def _validate_graph(self, G) -> None:
        if self._config.require_connected and G.order() > 0:
            if not G.is_connected():
                raise AssertionError("Graph ist nicht zusammenhängend.")

        if self._config.require_regular and G.order() > 0:
            degrees = [G.degree(v) for v in G.vertices()]
            if len(set(degrees)) != 1:
                raise AssertionError(
                    f"Graph ist nicht regulär; Grade={sorted(set(degrees))}"
                )

        if self._config.expected_degree is not None and G.order() > 0:
            expected = self._config.expected_degree
            if any(G.degree(v) != expected for v in G.vertices()):
                raise AssertionError(
                    f"Nicht alle Vertices haben erwarteten Grad {expected}."
                )

    # ------------------------------------------------------------------
    # Schalenanalyse
    # ------------------------------------------------------------------

    def pair_norm_distribution(
        self,
        points: Optional[Sequence[Any]] = None
    ) -> Dict[Any, int]:
        pts = list(points) if points is not None else self.get_shell_points()
        dist: Dict[Any, int] = {}

        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                n = self.icosian_norm(pts[j] - pts[i])
                dist[n] = dist.get(n, 0) + 1

        return dict(sorted(dist.items(), key=lambda kv: kv[0]))

    def local_shell_distribution(
        self,
        vertex_index: int = 0,
        points: Optional[Sequence[Any]] = None
    ) -> Dict[Any, int]:
        pts = list(points) if points is not None else self.get_shell_points()

        if not (0 <= vertex_index < len(pts)):
            raise IndexError("vertex_index außerhalb des gültigen Bereichs.")

        q = pts[vertex_index]
        dist: Dict[Any, int] = {}

        for j, r in enumerate(pts):
            if j == vertex_index:
                continue
            n = self.icosian_norm(r - q)
            dist[n] = dist.get(n, 0) + 1

        return dict(sorted(dist.items(), key=lambda kv: kv[0]))

    def shell_summary(self) -> Dict[str, Any]:
        pts = self.get_shell_points()
        summary = {
            "norm_shell": self._config.norm_shell,
            "number_of_points": len(pts),
            "edge_norm": self._config.edge_norm,
            "pair_norm_distribution": self.pair_norm_distribution(pts),
            "config": asdict(self._config),
            "metadata": self.get_metadata(),
        }

        if self._distance_graph is not None:
            summary["distance_graph"] = {
                "vertices": self._distance_graph.order(),
                "edges": self._distance_graph.size(),
                "degree_set": sorted(
                    set(self._distance_graph.degree(v)
                        for v in self._distance_graph.vertices())
                ),
            }

        if self._bk_graph is not None:
            summary["bk_graph"] = {
                "vertices": self._bk_graph.order(),
                "edges": self._bk_graph.size(),
                "degree_set": sorted(
                    set(self._bk_graph.degree(v)
                        for v in self._bk_graph.vertices())
                ),
            }

        return summary

    # ------------------------------------------------------------------
    # Export / Externalisierung
    # ------------------------------------------------------------------

    def export_json(self, filename: str) -> None:
        """
        Exportiert Parameter und Schalenstatistik.
        Quaternionen selbst werden als Strings serialisiert.
        """
        import json

        pts = self.get_shell_points()
        payload = self.shell_summary()
        payload["points"] = [str(q) for q in pts]

        with open(filename, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False, default=str)

    def export_graphml(self, filename: str, which: str = "distance") -> None:
        self._require_sage()

        if which == "distance":
            G = self._distance_graph or self.build_distance_graph()
        elif which == "bk":
            if self._bk_graph is None:
                raise ValueError("BK-Graph wurde noch nicht erzeugt.")
            G = self._bk_graph
        else:
            raise ValueError("which muss 'distance' oder 'bk' sein.")

        G.export_to_file(filename, format="graphml")

    def clone(self) -> "BambergKugel":
        """
        Erzeugt eine neue parametrisierte Instanz mit derselben Konfiguration
        und derselben benutzerdefinierten Kantenregel.
        """
        cfg = BambergKugelConfig(**asdict(self._config))
        other = BambergKugel(cfg)
        other.set_edge_selector(self._edge_selector)
        if self._icosian_basis is not None:
            other.set_icosian_basis(self._icosian_basis)
        return other

    # ------------------------------------------------------------------
    # State Management
    # ------------------------------------------------------------------

    def _invalidate_shell_state(self) -> None:
        self._shell_points = None
        self._distance_graph = None
        self._bk_graph = None
        self._metadata = {}

    def _invalidate_computed_state(self) -> None:
        self._group_2I = None
        self._icosian_basis = None
        self._invalidate_shell_state()


def interactive_dialog():
    """
    Einfacher CLI-Dialog für SageMath 10.5.
    """
    if not SAGE_AVAILABLE:
        print(
            "FEHLER: SageMath ist nicht verfügbar.\n"
            "Starte mit: sage -python bamberg_kugel.py"
        )
        return

    print("=" * 72)
    print("BAMBERG-KUGEL — parametrisierbare Normschalen-Engine")
    print("=" * 72)

    bk = BambergKugel()

    raw_shell = input("Normschale Q_I(q) eingeben [Default 1]: ").strip()
    if raw_shell:
        # Sage verarbeitet rationale Strings sauber über QQ
        try:
            shell = QQ(raw_shell)
        except Exception:
            raise ValueError(
                "Normschale bitte als ganze/rationale Zahl eingeben, z.B. 1, 2, 3/2."
            )
        bk.set_norm_shell(shell)

    if bk.get_norm_shell() != 1:
        raw_bound = input(
            "Suchschranke für 8 Icosianenkoeffizienten "
            f"[Default {bk.get_max_search_coefficient()}]: "
        ).strip()
        if raw_bound:
            bk.set_max_search_coefficient(int(raw_bound))

    raw_edge = input("Kantennorm Q_I(r-q) [Default 2]: ").strip()
    if raw_edge:
        bk.set_edge_norm(QQ(raw_edge))

    print("\nErzeuge Normschale ...")
    pts = bk.generate_shell()

    print(f"Anzahl gefundener Punkte: {len(pts)}")
    print(f"Erzeugungsmodus: {bk.get_metadata().get('shell_generation_mode')}")
    print(
        "Vollständigkeit zertifiziert:",
        bk.get_metadata().get("shell_complete", False)
    )

    print("\nBerechne Paar-Normverteilung ...")
    dist = bk.pair_norm_distribution()
    for norm, count in dist.items():
        print(f"  Q_I(delta) = {norm}: {count} ungeordnete Paare")

    create_graph = input(
        "\nVollständigen Distanzgraphen für die gewählte Kantennorm erzeugen? [j/N]: "
    ).strip().lower()

    if create_graph in {"j", "ja", "y", "yes"}:
        # Ein allgemeiner Distanzgraph ist nicht automatisch zusammenhängend.
        bk.set_require_connected(False)
        G = bk.build_distance_graph()
        print(f"Graph: |V|={G.order()}, |E|={G.size()}")
        print(
            "Gradmenge:",
            sorted(set(G.degree(v) for v in G.vertices()))
        )

    export = input("\nJSON-Zusammenfassung exportieren? Dateiname oder Enter: ").strip()
    if export:
        bk.export_json(export)
        print(f"Exportiert: {export}")

    print("\nFertig.")
    return bk


if __name__ == "__main__":
    interactive_dialog()
