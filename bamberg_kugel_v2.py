#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
bamberg_kugel_v2.py
====================

Parametrisierbare SageMath-10.5-Klasse zur exakten Erzeugung und Analyse
von Icosianen-Normschalen und daraus abgeleiteten Bamberg-Kugel-Graphen.

Zentrale Trennung
-----------------
(1) Icosianen-Normschale:
        S_N = {q in I : Q_I(q) = N}

(2) vollständiger Distanzgraph:
        Gamma_{N,d}: q ~ r  <=>  Q_I(r-q) = d

(3) spezieller Bamberg-Kugel-Teilgraph:
        BK_N <= Gamma_{N,d}
    mit einer expliziten, austauschbaren Selektionsregel.

Für N=1 steht zusätzlich die klassische 120-elementige Einheitsschale 2I
zur Verfügung. Eine spezielle BK-I1-C2-Voltage-Lift-Konstruktion ist NICHT
stillschweigend in höhere Normschalen fortgesetzt; höhere BK-Selektionen
müssen explizit als Regel angegeben werden.

Wesentliche Neuerung gegenüber v1
---------------------------------
Die exakte Shell-Engine enumeriert NICHT mehr einen willkürlichen Würfel
[-B,B]^8. Stattdessen wird

    * die Z-Hülle aller 120 Einheitsicosianen bestimmt,
    * daraus eine echte Rang-8-Z-Basis des Icosianengitters konstruiert,
    * die exakte positive definite Quadratikform Q_I in dieser Basis gebildet,
    * und SageMaths exakte Kurzvektor-Enumeration für Q_I(v)=N verwendet.

Damit ist generate_shell_exact(N) eine vollständige endliche Enumeration
für nichtnegative ganzzahlige Normschalen N, vorausgesetzt die in dieser
Datei verwendete 2I-Realisierung ist erfolgreich zertifiziert.

Ausführung:
    sage -python bamberg_kugel_v2.py

Import:
    from bamberg_kugel_v2 import (
        BambergKugel, BambergKugelConfig,
        BKSelectionRule, MetricSelectionRule
    )
"""

from dataclasses import dataclass, asdict
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple
import json

try:
    from sage.all import (
        QQ, ZZ, PolynomialRing, NumberField, QuaternionAlgebra,
        Graph, matrix, vector, lcm, QuadraticForm
    )
    SAGE_AVAILABLE = True
except Exception:
    SAGE_AVAILABLE = False


# ======================================================================
# Konfiguration
# ======================================================================

@dataclass
class BambergKugelConfig:
    norm_shell: int = 1
    edge_norm: int = 2

    # Graphprüfungen
    require_connected: bool = False
    require_regular: bool = False
    expected_degree: Optional[int] = None

    # Darstellung / Quotienten
    use_antipodal_quotient: bool = False
    deduplicate_signs: bool = False
    include_zero: bool = False

    # Exakte Rechnung ist im mathematischen Kern verpflichtend.
    exact_arithmetic: bool = True

    # Nur Legacy/Fallback; generate_shell_exact() benötigt dies NICHT.
    max_search_coefficient: int = 2


# ======================================================================
# Externalisierbare Selektionsregeln
# ======================================================================

class BKSelectionRule:
    """
    Abstrakte Schnittstelle für BK-Teilgraphregeln.

    Subklassen implementieren:
        admissible(q, r, context) -> bool
    """

    name = "abstract"

    def admissible(self, q, r, context: "BambergKugel") -> bool:
        raise NotImplementedError

    def __call__(self, q, r, context: "BambergKugel") -> bool:
        return self.admissible(q, r, context)


class MetricSelectionRule(BKSelectionRule):
    """
    Wählt alle Paare mit Q_I(r-q) = edge_norm.
    Das ergibt den vollständigen Distanzgraphen, nicht automatisch
    die spezielle kubische BK-I1-Struktur.
    """

    name = "metric"

    def __init__(self, edge_norm=None):
        self._edge_norm = edge_norm

    def get_edge_norm(self):
        return self._edge_norm

    def set_edge_norm(self, value):
        self._edge_norm = value

    def admissible(self, q, r, context: "BambergKugel") -> bool:
        d = context.get_edge_norm() if self._edge_norm is None else self._edge_norm
        return context.icosian_norm(r - q) == d


class CallableSelectionRule(BKSelectionRule):
    """Adapter für beliebige Python-Funktionen selector(q,r,context)."""

    name = "callable"

    def __init__(self, function: Callable[[Any, Any, "BambergKugel"], bool]):
        if not callable(function):
            raise TypeError("function muss aufrufbar sein.")
        self._function = function

    def get_function(self):
        return self._function

    def set_function(self, function):
        if not callable(function):
            raise TypeError("function muss aufrufbar sein.")
        self._function = function

    def admissible(self, q, r, context: "BambergKugel") -> bool:
        return bool(self._function(q, r, context))


# ======================================================================
# Hauptklasse
# ======================================================================

class BambergKugel:
    """
    Exakte Normschalen-/Graph-Engine für das Icosianenmodell.
    """

    def __init__(self, config: Optional[BambergKugelConfig] = None):
        self._config = config or BambergKugelConfig()

        self._K = None
        self._sqrt5 = None
        self._H = None
        self._one = None
        self._minus_one = None

        self._group_2I: Optional[Set[Any]] = None

        # Icosianengitter
        self._icosian_basis: Optional[List[Any]] = None
        self._icosian_basis_vectors = None
        self._gram_matrix = None
        self._quadratic_form = None

        # Schalen und Graphen
        self._shell_cache: Dict[int, Tuple[Any, ...]] = {}
        self._shell_vector_cache: Dict[int, Tuple[Any, ...]] = {}
        self._shell_points: Optional[List[Any]] = None
        self._distance_graph = None
        self._bk_graph = None

        self._selection_rule: Optional[BKSelectionRule] = None
        self._metadata: Dict[str, Any] = {}

        if SAGE_AVAILABLE:
            self._initialize_algebra()

    # ------------------------------------------------------------------
    # Konfigurations-Getter / Setter
    # ------------------------------------------------------------------

    def get_config(self) -> BambergKugelConfig:
        return self._config

    def set_config(self, config: BambergKugelConfig) -> None:
        if not isinstance(config, BambergKugelConfig):
            raise TypeError("config muss BambergKugelConfig sein.")
        self._config = config
        self._invalidate_shell_state()

    def get_norm_shell(self) -> int:
        return self._config.norm_shell

    def set_norm_shell(self, value) -> None:
        value = self._validate_shell_norm(value)
        self._config.norm_shell = int(value)
        self._shell_points = None
        self._distance_graph = None
        self._bk_graph = None

    def get_edge_norm(self):
        return self._config.edge_norm

    def set_edge_norm(self, value) -> None:
        self._config.edge_norm = value
        self._distance_graph = None
        self._bk_graph = None

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
        self._shell_points = None
        self._distance_graph = None
        self._bk_graph = None

    def get_deduplicate_signs(self) -> bool:
        return self._config.deduplicate_signs

    def set_deduplicate_signs(self, value: bool) -> None:
        self._config.deduplicate_signs = bool(value)
        self._shell_points = None
        self._distance_graph = None
        self._bk_graph = None

    def get_include_zero(self) -> bool:
        return self._config.include_zero

    def set_include_zero(self, value: bool) -> None:
        self._config.include_zero = bool(value)
        self._shell_points = None

    def get_max_search_coefficient(self) -> int:
        """Legacy/Fallback-Parameter."""
        return self._config.max_search_coefficient

    def set_max_search_coefficient(self, value: int) -> None:
        value = int(value)
        if value < 0:
            raise ValueError("max_search_coefficient muss >= 0 sein.")
        self._config.max_search_coefficient = value

    def get_selection_rule(self) -> Optional[BKSelectionRule]:
        return self._selection_rule

    def set_selection_rule(self, rule) -> None:
        if rule is None:
            self._selection_rule = None
        elif isinstance(rule, BKSelectionRule):
            self._selection_rule = rule
        elif callable(rule):
            self._selection_rule = CallableSelectionRule(rule)
        else:
            raise TypeError("rule muss BKSelectionRule, callable oder None sein.")
        self._bk_graph = None

    # Kompatibilität mit v1
    def get_edge_selector(self):
        return self._selection_rule

    def set_edge_selector(self, selector) -> None:
        self.set_selection_rule(selector)

    # ------------------------------------------------------------------
    # Algebra-Getter
    # ------------------------------------------------------------------

    def get_number_field(self):
        self._require_sage()
        return self._K

    def get_quaternion_algebra(self):
        self._require_sage()
        return self._H

    def get_metadata(self) -> Dict[str, Any]:
        return dict(self._metadata)

    def get_group_2I(self) -> Set[Any]:
        self._require_sage()
        if self._group_2I is None:
            self._group_2I = self._build_binary_icosahedral_group()
        return set(self._group_2I)

    # ------------------------------------------------------------------
    # Initialisierung / Grundalgebra
    # ------------------------------------------------------------------

    def _require_sage(self) -> None:
        if not SAGE_AVAILABLE:
            raise RuntimeError(
                "SageMath nicht verfügbar. Start z.B. mit:\n"
                "    sage -python bamberg_kugel_v2.py"
            )

    def _initialize_algebra(self) -> None:
        R = PolynomialRing(QQ, "t")
        t = R.gen()
        self._K = NumberField(t**2 - 5, "sqrt5")
        self._sqrt5 = self._K.gen()
        self._H = QuaternionAlgebra(self._K, -1, -1)
        self._one = self._H(1)
        self._minus_one = self._H(-1)

    def _coeffs_quaternion(self, q) -> Tuple[Any, Any, Any, Any]:
        try:
            c = tuple(q.coefficient_tuple())
        except Exception:
            try:
                c = tuple(q)
            except Exception as exc:
                raise TypeError(
                    f"Quaternionenkoeffizienten nicht lesbar: {q}"
                ) from exc

        if len(c) != 4:
            raise ValueError(f"Erwartet 4 Quaternionenkoeffizienten, erhalten: {c}")
        return c

    def reduced_norm(self, q):
        self._require_sage()
        try:
            return q.reduced_norm()
        except Exception:
            a, b, c, d = self._coeffs_quaternion(q)
            return a*a + b*b + c*c + d*d

    def split_Qsqrt5(self, z) -> Tuple[Any, Any]:
        """z = x + y*sqrt(5) -> (x,y) in QQ."""
        self._require_sage()
        z = self._K(z)
        coeffs = list(z.list())
        x = QQ(coeffs[0]) if len(coeffs) >= 1 else QQ(0)
        y = QQ(coeffs[1]) if len(coeffs) >= 2 else QQ(0)
        return x, y

    def icosian_norm(self, q):
        """
        Q_I(q) = x+y, falls N_H(q)=x+sqrt(5)y.
        """
        x, y = self.split_Qsqrt5(self.reduced_norm(q))
        return x + y

    def bilinear_form(self, q, r):
        """
        Polarisierung:
            B_I(q,r) = (Q_I(q+r)-Q_I(q)-Q_I(r))/2.
        """
        return (
            self.icosian_norm(q + r)
            - self.icosian_norm(q)
            - self.icosian_norm(r)
        ) / 2

    def get_bilinear_form(self):
        """Externalisierbarer Funktionszugriff."""
        return self.bilinear_form

    # ------------------------------------------------------------------
    # 8D-Q-Koordinaten
    # ------------------------------------------------------------------

    def q_vector_8(self, q):
        """
        q -> QQ^8 in der Feldbasis (1,sqrt5) je Quaternionenkoordinate.
        """
        self._require_sage()
        out = []
        for c in self._coeffs_quaternion(q):
            x, y = self.split_Qsqrt5(c)
            out.extend([x, y])
        return vector(QQ, out)

    def vector8_to_quaternion(self, v):
        """
        Inverse Abbildung QQ^8 -> H.
        """
        self._require_sage()
        if len(v) != 8:
            raise ValueError("8 Koordinaten erwartet.")

        s5 = self._sqrt5
        coeff = [
            self._K(v[0]) + self._K(v[1])*s5,
            self._K(v[2]) + self._K(v[3])*s5,
            self._K(v[4]) + self._K(v[5])*s5,
            self._K(v[6]) + self._K(v[7])*s5,
        ]

        i = self._H.gen(0)
        j = self._H.gen(1)
        k = i*j

        return coeff[0] + coeff[1]*i + coeff[2]*j + coeff[3]*k

    def exact_key(self, q):
        return tuple(self.q_vector_8(q))

    def _canonical_sign(self, q):
        return min((q, -q), key=self.exact_key)

    # ------------------------------------------------------------------
    # Binäre Ikosaedergruppe 2I
    # ------------------------------------------------------------------

    def _hurwitz_units(self) -> Set[Any]:
        self._require_sage()

        H = self._H
        K = self._K
        one = H(1)
        i = H.gen(0)
        j = H.gen(1)
        k = i*j

        elems = set()

        for u in (one, i, j, k):
            elems.add(u)
            elems.add(-u)

        half = K(1)/2
        for s0 in (-1, 1):
            for s1 in (-1, 1):
                for s2 in (-1, 1):
                    for s3 in (-1, 1):
                        elems.add(
                            half*(s0*one + s1*i + s2*j + s3*k)
                        )

        assert len(elems) == 24
        return elems

    def _icosahedral_seed(self):
        """
        Goldener Standard-Seed für die H4/2I-Realisierung.
        """
        self._require_sage()

        phi = (self._K(1) + self._sqrt5)/2
        invphi = 1/phi

        one = self._H(1)
        i = self._H.gen(0)
        j = self._H.gen(1)

        return (phi*one + i + invphi*j)/2

    def _group_closure(self, generators: Iterable[Any]) -> Set[Any]:
        gens = list(generators)
        gens += [g.inverse() for g in gens]

        G = {self._one}
        changed = True

        while changed:
            changed = False
            old = list(G)
            for x in old:
                for g in gens:
                    z = x*g
                    if z not in G:
                        G.add(z)
                        changed = True

        return G

    def _build_binary_icosahedral_group(self) -> Set[Any]:
        seed = self._icosahedral_seed()
        generators = list(self._hurwitz_units()) + [seed, -seed]

        G = self._group_closure(generators)

        assert len(G) == 120, (
            f"2I-Konstruktion ergab {len(G)} statt 120 Elemente."
        )
        assert self._one in G
        assert self._minus_one in G
        assert all(g.inverse() in G for g in G)
        assert all(a*b in G for a in G for b in G)
        assert all(self.icosian_norm(g) == 1 for g in G)

        self._metadata["2I_order"] = 120
        self._metadata["2I_group_checks"] = "PASS"
        return G

    def unit_shell(self) -> List[Any]:
        return sorted(self.get_group_2I(), key=self.exact_key)

    # ------------------------------------------------------------------
    # Echte Z-Basis des Icosianengitters
    # ------------------------------------------------------------------

    def _common_denominator(self, M) -> Any:
        dens = []
        for x in M.list():
            try:
                dens.append(QQ(x).denominator())
            except Exception:
                pass
        return lcm(dens) if dens else ZZ(1)

    def _build_icosian_lattice_basis(self) -> List[Any]:
        """
        Konstruktion der Z-Hülle aller 120 Einheitsicosianen.

        Sei R die 120x8-Matrix ihrer QQ^8-Koordinaten.
        Nach Multiplikation mit einem gemeinsamen Nenner D ist D*R integral.
        Die Zeilenbasis des von D*R erzeugten ZZ-Moduls liefert nach Division
        durch D eine Z-Basis der ursprünglichen Icosianen-Z-Hülle.
        """
        self._require_sage()

        units = sorted(self.get_group_2I(), key=self.exact_key)
        R = matrix(QQ, [self.q_vector_8(q) for q in units])

        assert R.rank() == 8

        D = self._common_denominator(R)
        RZ = matrix(ZZ, [[ZZ(D*x) for x in row] for row in R.rows()])

        # Exakte Zeilen-Z-Hülle, unabhängig von einer willkürlichen
        # Auswahl rational unabhängiger 2I-Elemente.
        L = RZ.row_module()
        BZ = L.basis_matrix()

        assert BZ.nrows() == 8
        assert BZ.ncols() == 8
        assert BZ.rank() == 8

        BQ = matrix(QQ, [[QQ(x)/D for x in row] for row in BZ.rows()])
        basis = [self.vector8_to_quaternion(row) for row in BQ.rows()]

        # Jede Einheit muss ganzzahlig in dieser Basis liegen.
        # Mit Zeilenkonvention: c * BQ = unit_vector.
        Binv = BQ.inverse()
        for q in units:
            v = self.q_vector_8(q)
            c = v * Binv
            assert all(x in ZZ for x in c), (
                "Basis deckt Z-Hülle der 2I-Einheiten nicht exakt ab."
            )

        self._icosian_basis_vectors = BQ
        self._metadata["icosian_basis_rank"] = 8
        self._metadata["icosian_basis_common_denominator"] = str(D)
        self._metadata["icosian_basis_mode"] = "Z-span(2I)-row-lattice"

        return basis

    def get_icosian_basis(self) -> List[Any]:
        self._require_sage()
        if self._icosian_basis is None:
            self._icosian_basis = self._build_icosian_lattice_basis()
        return list(self._icosian_basis)

    def get_icosian_basis_vectors(self):
        self.get_icosian_basis()
        return self._icosian_basis_vectors.copy()

    def set_icosian_basis(self, basis: Sequence[Any]) -> None:
        """
        Setzt eine externe Rang-8-Z-Basis.

        ACHTUNG:
        Die Klasse kann nur prüfen, dass die Basis rational Rang 8 besitzt
        und dass die zugehörige Q_I-Form ganzzahlig ist. Ob die Basis exakt
        das vollständige Icosianengitter und nicht einen Unter-/Übergitter
        beschreibt, liegt bei extern gesetzten Basen in der Verantwortung
        des Aufrufers.
        """
        self._require_sage()

        if len(basis) != 8:
            raise ValueError("Genau 8 Basiselemente erwartet.")

        BQ = matrix(QQ, [self.q_vector_8(q) for q in basis])
        if BQ.rank() != 8:
            raise ValueError("Externe Basis hat nicht Rang 8.")

        self._icosian_basis = list(basis)
        self._icosian_basis_vectors = BQ
        self._gram_matrix = None
        self._quadratic_form = None
        self.clear_shell_cache()

        self._metadata["icosian_basis_mode"] = "external"
        self._metadata["external_basis_full_lattice_certified"] = False

    # ------------------------------------------------------------------
    # Gram-Matrix / ganzzahlige Quadratikform
    # ------------------------------------------------------------------

    def get_gram_matrix(self):
        """
        Gram-Matrix G mit
            Q_I(sum n_i b_i) = n^T G n.
        """
        self._require_sage()

        if self._gram_matrix is None:
            basis = self.get_icosian_basis()
            n = len(basis)

            G = matrix(
                QQ, n, n,
                lambda i, j: self.bilinear_form(basis[i], basis[j])
            )

            assert G == G.transpose()
            assert G.rank() == 8

            # Exakte Positivitätskontrolle via Hauptminoren (Sylvester).
            # Keine Floating-Point-Eigenwerte als Beweis.
            for k in range(1, 9):
                assert G[:k, :k].det() > 0, (
                    f"Gram-Matrix nicht positiv definit: Hauptminor {k}."
                )

            self._gram_matrix = G
            self._metadata["gram_determinant"] = str(G.det())
            self._metadata["gram_positive_definite"] = True

        return self._gram_matrix.copy()

    def get_quadratic_form(self):
        """
        Ganzzahlige Sage-QuadraticForm Q mit Q(n)=n^T G n.

        Für eine Integralform müssen gelten:
            G_ii in ZZ,
            2 G_ij in ZZ (i<j).
        """
        self._require_sage()

        if self._quadratic_form is None:
            G = self.get_gram_matrix()

            coeffs = []
            for i in range(8):
                diag = G[i, i]
                assert diag in ZZ, (
                    f"Nichtintegraler Diagonaleintrag G[{i},{i}]={diag}."
                )
                coeffs.append(ZZ(diag))

                for j in range(i + 1, 8):
                    cross = 2*G[i, j]
                    assert cross in ZZ, (
                        f"Nichtintegraler Kreuzkoeffizient 2G[{i},{j}]={cross}."
                    )
                    coeffs.append(ZZ(cross))

            Q = QuadraticForm(ZZ, 8, coeffs)

            # Matrixkonvention verifizieren, statt sie vorauszusetzen.
            assert Q.Gram_matrix_rational() == G

            # Basisvektoren müssen direkt übereinstimmen.
            for i, b in enumerate(self.get_icosian_basis()):
                e = [ZZ(0)]*8
                e[i] = ZZ(1)
                assert Q(e) == self.icosian_norm(b)

            self._quadratic_form = Q
            self._metadata["integral_quadratic_form"] = True

        return self._quadratic_form

    # ------------------------------------------------------------------
    # Exakte Normschalen-Enumeration
    # ------------------------------------------------------------------

    def _validate_shell_norm(self, N) -> int:
        if SAGE_AVAILABLE:
            q = QQ(N)
            if q.denominator() != 1:
                raise ValueError(
                    "Die aktuelle exakte Icosianen-Shell-Engine erwartet "
                    "nichtnegative ganzzahlige Q_I-Normen."
                )
            N = ZZ(q)
        else:
            N = int(N)

        if N < 0:
            raise ValueError("Q_I ist positiv definit; Normschale muss >= 0 sein.")
        return int(N)

    def clear_shell_cache(self) -> None:
        self._shell_cache.clear()
        self._shell_vector_cache.clear()
        self._shell_points = None
        self._distance_graph = None
        self._bk_graph = None

    def get_shell_cache_keys(self) -> List[int]:
        return sorted(self._shell_cache.keys())

    def lattice_vector_to_icosian(self, coeff_vector):
        basis = self.get_icosian_basis()
        if len(coeff_vector) != 8:
            raise ValueError("8 Gitterkoeffizienten erwartet.")

        q = self._H(0)
        for c, b in zip(coeff_vector, basis):
            q += ZZ(c)*b
        return q

    def generate_shell_vectors_exact(self, N=None) -> List[Any]:
        """
        Liefert ALLE Z^8-Koeffizientenvektoren mit Q_I=N.

        SageMath:
            QuadraticForm.short_vector_list_up_to_length(N+1)

        Die Methode liefert Listen für exakte ganzzahlige Normwerte
        0,...,N; wir entnehmen Index N.
        """
        self._require_sage()

        if N is None:
            N = self.get_norm_shell()
        N = self._validate_shell_norm(N)

        if N in self._shell_vector_cache:
            return [vector(ZZ, v) for v in self._shell_vector_cache[N]]

        Q = self.get_quadratic_form()

        # Default ist nicht "up to sign"; beide Vorzeichen werden benötigt.
        lists = Q.short_vector_list_up_to_length(ZZ(N + 1))
        vecs = [vector(ZZ, v) for v in lists[N]]

        # Harte Vollständigkeits-/Normkontrollen auf dem zurückgegebenen Satz.
        assert all(Q(v) == N for v in vecs)
        assert len(set(tuple(v) for v in vecs)) == len(vecs)

        self._shell_vector_cache[N] = tuple(tuple(ZZ(x) for x in v) for v in vecs)
        return vecs

    def generate_shell_exact(self, N=None) -> List[Any]:
        """
        Vollständige Icosianen-Normschale S_N.

        Für N=1 wird zusätzlich gegen die 120-elementige 2I-Schale geprüft.
        """
        self._require_sage()

        if N is None:
            N = self.get_norm_shell()
        N = self._validate_shell_norm(N)

        # Cache enthält die volle, nicht quotientierte Schale.
        if N not in self._shell_cache:
            vecs = self.generate_shell_vectors_exact(N)
            pts = [self.lattice_vector_to_icosian(v) for v in vecs]

            assert all(self.icosian_norm(q) == N for q in pts)
            assert len(set(pts)) == len(pts)

            if N == 0:
                expected = {self._H(0)}
                assert set(pts) == expected

            if N == 1:
                units = self.get_group_2I()
                assert set(pts) == units, (
                    "Norm-1-Schale stimmt nicht exakt mit den 120 Einheiten 2I überein."
                )

            self._shell_cache[N] = tuple(sorted(pts, key=self.exact_key))

        pts = list(self._shell_cache[N])

        if N == 0 and not self._config.include_zero:
            pts = []

        if self._config.deduplicate_signs or self._config.use_antipodal_quotient:
            pts = sorted(
                {self._canonical_sign(q) for q in pts},
                key=self.exact_key
            )

        self._shell_points = list(pts)
        self._metadata["shell_generation_mode"] = "exact_positive_definite_QF"
        self._metadata["shell_complete"] = True
        self._metadata["current_shell"] = N
        self._metadata["current_shell_cardinality"] = len(pts)

        return list(pts)

    # v1-kompatibel: standardmäßig jetzt EXAKT
    def generate_shell(self, norm_shell=None, coefficient_bound=None) -> List[Any]:
        if coefficient_bound is not None:
            return self.generate_shell_bounded(
                norm_shell=norm_shell,
                coefficient_bound=coefficient_bound
            )
        return self.generate_shell_exact(norm_shell)

    def get_shell_vectors(self, N=None) -> List[Any]:
        return self.generate_shell_vectors_exact(N)

    def get_shell_quaternions(self, N=None) -> List[Any]:
        return self.generate_shell_exact(N)

    def get_shell_cardinality(self, N=None, antipodal=False) -> int:
        old = self._config.use_antipodal_quotient
        try:
            self._config.use_antipodal_quotient = bool(antipodal)
            return len(self.generate_shell_exact(N))
        finally:
            self._config.use_antipodal_quotient = old

    # ------------------------------------------------------------------
    # Legacy: beschränkte Enumeration nur als Kontroll-/Fallbackmethode
    # ------------------------------------------------------------------

    def generate_shell_bounded(
        self,
        norm_shell=None,
        coefficient_bound=None
    ) -> List[Any]:
        """
        Alte endliche Würfelsuche. NICHT als Vollständigkeitsbeweis verwenden.
        Nützlich nur als unabhängiger Cross-Check kleiner Schalen.
        """
        self._require_sage()

        N = self.get_norm_shell() if norm_shell is None else norm_shell
        N = self._validate_shell_norm(N)

        B = (
            self._config.max_search_coefficient
            if coefficient_bound is None
            else int(coefficient_bound)
        )
        if B < 0:
            raise ValueError("coefficient_bound muss >= 0 sein.")

        candidate_count = (2*B + 1)**8
        if candidate_count > 5_000_000:
            raise ValueError(
                f"Legacy-Suchraum zu groß: (2*{B}+1)^8={candidate_count}."
            )

        basis = self.get_icosian_basis()
        points = set()

        def rec(i, q):
            if i == 8:
                if self.icosian_norm(q) == N:
                    points.add(q)
                return
            for c in range(-B, B + 1):
                rec(i + 1, q + c*basis[i])

        rec(0, self._H(0))

        self._metadata["bounded_crosscheck_complete"] = False
        self._metadata["bounded_crosscheck_bound"] = B
        self._metadata["bounded_crosscheck_candidates"] = candidate_count

        return sorted(points, key=self.exact_key)

    # ------------------------------------------------------------------
    # Shell-Zugriff
    # ------------------------------------------------------------------

    def get_shell_points(self) -> List[Any]:
        if self._shell_points is None:
            return self.generate_shell_exact(self.get_norm_shell())
        return list(self._shell_points)

    def set_shell_points(self, points: Sequence[Any]) -> None:
        """
        Setzt externe Punkte für Graphanalysen. Dies zertifiziert NICHT,
        dass die Punkte eine vollständige Icosianen-Normschale bilden.
        """
        self._shell_points = list(points)
        self._distance_graph = None
        self._bk_graph = None
        self._metadata["external_shell_points"] = True

    # ------------------------------------------------------------------
    # Distanz- und BK-Graphen
    # ------------------------------------------------------------------

    def build_distance_graph(
        self,
        points: Optional[Sequence[Any]] = None,
        edge_norm=None
    ):
        self._require_sage()

        pts = list(points) if points is not None else self.get_shell_points()
        d = self.get_edge_norm() if edge_norm is None else edge_norm

        G = Graph()
        G.add_vertices(range(len(pts)))

        for i in range(len(pts)):
            qi = pts[i]
            for j in range(i + 1, len(pts)):
                if self.icosian_norm(pts[j] - qi) == d:
                    G.add_edge(i, j)

        self._distance_graph = G
        self._metadata["distance_graph_vertices"] = G.order()
        self._metadata["distance_graph_edges"] = G.size()
        self._metadata["distance_graph_edge_norm"] = str(d)

        self._validate_graph(G)
        return G

    def build_bk_graph(
        self,
        points: Optional[Sequence[Any]] = None,
        rule: Optional[BKSelectionRule] = None
    ):
        """
        Baut einen explizit selektierten BK-Teilgraphen.

        Keine Default-Regel wird erfunden. Das ist bei höheren Schalen
        methodisch wichtig.
        """
        self._require_sage()

        pts = list(points) if points is not None else self.get_shell_points()
        active = rule or self._selection_rule

        if active is None:
            raise ValueError(
                "Keine BK-Selektionsregel gesetzt. "
                "Nutze set_selection_rule(...) oder rule=..."
            )

        if callable(active) and not isinstance(active, BKSelectionRule):
            active = CallableSelectionRule(active)

        G = Graph()
        G.add_vertices(range(len(pts)))

        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                if active(pts[i], pts[j], self):
                    G.add_edge(i, j)

        self._bk_graph = G
        self._metadata["bk_graph_vertices"] = G.order()
        self._metadata["bk_graph_edges"] = G.size()
        self._metadata["bk_selection_rule"] = getattr(active, "name", type(active).__name__)

        self._validate_graph(G)
        return G

    def get_distance_graph(self):
        return self._distance_graph

    def get_bk_graph(self):
        return self._bk_graph

    def _validate_graph(self, G) -> None:
        if self._config.require_connected and G.order() > 0:
            assert G.is_connected(), "Graph ist nicht zusammenhängend."

        if self._config.require_regular and G.order() > 0:
            degree_set = {G.degree(v) for v in G.vertices()}
            assert len(degree_set) == 1, (
                f"Graph ist nicht regulär; Grade={sorted(degree_set)}"
            )

        if self._config.expected_degree is not None and G.order() > 0:
            d = self._config.expected_degree
            assert all(G.degree(v) == d for v in G.vertices()), (
                f"Nicht alle Vertices haben den erwarteten Grad {d}."
            )

    # ------------------------------------------------------------------
    # Analysen
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

        if not 0 <= vertex_index < len(pts):
            raise IndexError("vertex_index außerhalb des Bereichs.")

        q = pts[vertex_index]
        dist: Dict[Any, int] = {}

        for j, r in enumerate(pts):
            if j == vertex_index:
                continue
            n = self.icosian_norm(r - q)
            dist[n] = dist.get(n, 0) + 1

        return dict(sorted(dist.items(), key=lambda kv: kv[0]))

    def shell_summary(self, N=None) -> Dict[str, Any]:
        if N is not None:
            points = self.generate_shell_exact(N)
        else:
            points = self.get_shell_points()

        summary = {
            "norm_shell": self._metadata.get("current_shell", self.get_norm_shell()),
            "number_of_points": len(points),
            "edge_norm": str(self.get_edge_norm()),
            "pair_norm_distribution": {
                str(k): v for k, v in self.pair_norm_distribution(points).items()
            },
            "config": asdict(self._config),
            "metadata": self.get_metadata(),
        }

        if self._distance_graph is not None:
            summary["distance_graph"] = {
                "vertices": self._distance_graph.order(),
                "edges": self._distance_graph.size(),
                "degree_set": sorted({
                    self._distance_graph.degree(v)
                    for v in self._distance_graph.vertices()
                }),
            }

        if self._bk_graph is not None:
            summary["bk_graph"] = {
                "vertices": self._bk_graph.order(),
                "edges": self._bk_graph.size(),
                "degree_set": sorted({
                    self._bk_graph.degree(v)
                    for v in self._bk_graph.vertices()
                }),
            }

        return summary

    def enumerate_shell_range(self, start: int, stop: int) -> Dict[int, int]:
        """
        Vollständige Kardinalitäten S_start,...,S_stop.
        """
        start = self._validate_shell_norm(start)
        stop = self._validate_shell_norm(stop)
        if stop < start:
            raise ValueError("stop muss >= start sein.")

        return {
            N: len(self.generate_shell_exact(N))
            for N in range(start, stop + 1)
        }

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def export_json(self, filename: str, N=None) -> None:
        points = self.generate_shell_exact(N) if N is not None else self.get_shell_points()
        payload = self.shell_summary(N)
        payload["points"] = [str(q) for q in points]

        with open(filename, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False, default=str)

    def export_graphml(self, filename: str, which: str = "distance") -> None:
        self._require_sage()

        if which == "distance":
            G = self._distance_graph
            if G is None:
                G = self.build_distance_graph()
        elif which == "bk":
            G = self._bk_graph
            if G is None:
                raise ValueError("BK-Graph wurde noch nicht erzeugt.")
        else:
            raise ValueError("which muss 'distance' oder 'bk' sein.")

        G.export_to_file(filename, format="graphml")

    def clone(self) -> "BambergKugel":
        cfg = BambergKugelConfig(**asdict(self._config))
        other = BambergKugel(cfg)
        other.set_selection_rule(self._selection_rule)

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


# ======================================================================
# Interaktiver Dialog
# ======================================================================

def interactive_dialog():
    if not SAGE_AVAILABLE:
        print(
            "FEHLER: SageMath nicht verfügbar.\n"
            "Bitte starten mit:\n"
            "    sage -python bamberg_kugel_v2.py"
        )
        return None

    print("=" * 78)
    print("BAMBERG-KUGEL v2 — EXAKTE ICOSIANEN-NORMSCHALEN")
    print("=" * 78)
    print("Die Schale wird über die positive definite Rang-8-Quadratikform")
    print("vollständig enumeriert; keine künstliche Koeffizientenschranke nötig.")
    print()

    bk = BambergKugel()

    raw = input("Normschale N = Q_I(q) [1]: ").strip()
    N = 1 if not raw else int(raw)
    bk.set_norm_shell(N)

    raw_edge = input("Kantennorm d = Q_I(r-q) [2]: ").strip()
    if raw_edge:
        bk.set_edge_norm(QQ(raw_edge))

    print("\n[1/4] Konstruiere/zertifiziere 2I und Icosianen-Z-Basis ...")
    basis = bk.get_icosian_basis()
    print(f"      Rang der Icosianenbasis: {len(basis)}")

    print("[2/4] Konstruiere exakte Gram-Matrix und Quadratikform ...")
    G = bk.get_gram_matrix()
    print(f"      det(G) = {G.det()}")

    print(f"[3/4] Enumeriere vollständige Normschale S_{N} ...")
    shell = bk.generate_shell_exact(N)
    print(f"      |S_{N}| = {len(shell)}")
    print("      Vollständigkeit: JA (exakte Kurzvektor-Enumeration)")

    do_pairs = input(
        "[4/4] Vollständige Paar-Normverteilung berechnen? [j/N]: "
    ).strip().lower()

    if do_pairs in {"j", "ja", "y", "yes"}:
        dist = bk.pair_norm_distribution(shell)
        print("\nPaar-Normverteilung:")
        for norm, count in dist.items():
            print(f"      Q_I(delta)={norm}: {count}")

    do_graph = input(
        "\nVollständigen Distanzgraphen Gamma_{N,d} erzeugen? [j/N]: "
    ).strip().lower()

    if do_graph in {"j", "ja", "y", "yes"}:
        graph = bk.build_distance_graph(shell)
        degree_set = sorted({
            graph.degree(v) for v in graph.vertices()
        })
        print(
            f"      |V|={graph.order()}, |E|={graph.size()}, "
            f"Gradmenge={degree_set}"
        )

    raw_range = input(
        "\nWeitere Schalen als Kardinalitätstabelle berechnen "
        "(z.B. 1:10) oder Enter: "
    ).strip()

    if raw_range:
        a, b = raw_range.split(":")
        counts = bk.enumerate_shell_range(int(a), int(b))
        print("\nSchalenkardinalitäten:")
        for n, c in counts.items():
            print(f"      |S_{n}| = {c}")

    export = input(
        "\nJSON-Exportdatei oder Enter zum Beenden: "
    ).strip()

    if export:
        bk.export_json(export, N=N)
        print(f"Exportiert: {export}")

    print("\nFertig.")
    return bk


if __name__ == "__main__":
    interactive_dialog()
