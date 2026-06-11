from __future__ import annotations

import math
import os
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Callable, List, Optional, Sequence, Tuple, Union

import numpy as np


Number = Union[int, float, Fraction]


@dataclass(frozen=True)
class Quaternion:
    """Hamilton-Quaternion (i^2 = j^2 = k^2 = -1, i*j = k)."""

    a: Fraction
    b: Fraction
    c: Fraction
    d: Fraction

    def reduced_norm(self) -> Fraction:
        return self.a**2 + self.b**2 + self.c**2 + self.d**2

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0 and self.c == 0 and self.d == 0

    @staticmethod
    def _as_frac(x: Number) -> Fraction:
        if isinstance(x, Fraction):
            return x
        if isinstance(x, int):
            return Fraction(x, 1)
        return Fraction(x).limit_denominator() if x != int(x) else Fraction(int(x), 1)

    def __add__(self, o: object) -> Quaternion:
        if not isinstance(o, Quaternion):
            return NotImplemented
        return Quaternion(self.a + o.a, self.b + o.b, self.c + o.c, self.d + o.d)

    def __neg__(self) -> Quaternion:
        return Quaternion(-self.a, -self.b, -self.c, -self.d)

    def conjugate(self) -> Quaternion:
        """Hamilton-Konjugation: a - bi - cj - dk."""
        return Quaternion(self.a, -self.b, -self.c, -self.d)

    def __sub__(self, o: object) -> Quaternion:
        if not isinstance(o, Quaternion):
            return NotImplemented
        return self + (-o)

    def __mul__(self, o: object) -> Quaternion:
        if isinstance(o, (int, float, Fraction)):
            t = self._as_frac(o)
            return Quaternion(
                t * self.a, t * self.b, t * self.c, t * self.d
            )
        if not isinstance(o, Quaternion):
            return NotImplemented
        a1, b1, c1, d1 = self.a, self.b, self.c, self.d
        a2, b2, c2, d2 = o.a, o.b, o.c, o.d
        return Quaternion(
            a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
            a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
            a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2,
        )

    def __rmul__(self, o: object) -> Quaternion:
        if isinstance(o, (int, float, Fraction)):
            t = self._as_frac(o)
            return Quaternion(
                t * self.a, t * self.b, t * self.c, t * self.d
            )
        return NotImplemented

    def __str__(self) -> str:
        return f"{self.a} + {self.b}*i + {self.c}*j + {self.d}*k"


def quat(a: Number = 0, b: Number = 0, c: Number = 0, d: Number = 0) -> Quaternion:
    return Quaternion(
        Quaternion._as_frac(a),
        Quaternion._as_frac(b),
        Quaternion._as_frac(c),
        Quaternion._as_frac(d),
    )


def _left_regular_4x4(q: Quaternion) -> List[List[float]]:
    """Linksmultiplikation q * v im R^4-Basis(1, i, j, k)."""
    a, b, c, d = (float(x) for x in (q.a, q.b, q.c, q.d))
    return [
        [a, -b, -c, -d],
        [b, a, -d, c],
        [c, d, a, -b],
        [d, -c, b, a],
    ]


def _det_real_square(a: List[List[float]]) -> float:
    n = len(a)
    if n == 0:
        return 1.0
    if any(len(r) != n for r in a):
        raise ValueError("determinant: quadratische Matrix nötig")
    try:
        import numpy as _np  # type: ignore[import-not-found]

        return float(_np.linalg.det(_np.array(a, dtype=float)))
    except Exception:
        pass
    m = [row[:] for row in a]
    sign = 1.0
    eps = 1e-15
    for c in range(n):
        p = c
        for r in range(c + 1, n):
            if abs(m[r][c]) > abs(m[p][c]):
                p = r
        if abs(m[p][c]) < eps:
            return 0.0
        if p != c:
            m[c], m[p] = m[p], m[c]
            sign = -sign
        piv = m[c][c]
        for r in range(c + 1, n):
            f = m[r][c] / piv
            for j in range(c, n):
                m[r][j] -= f * m[c][j]
    d = sign
    for i in range(n):
        d *= m[i][i]
    return d


def _qmatrix_to_real_left_regular(M: "QMatrix") -> np.ndarray:
    """Linkregular-Darstellung m×n über H → (4m)×(4n) reell."""
    m, n = M.nrows(), M.ncols()
    R = np.zeros((4 * m, 4 * n), dtype=np.float64)
    for i in range(m):
        for j in range(n):
            blk = _left_regular_4x4(M._rows[i][j])
            for bi in range(4):
                for bj in range(4):
                    R[4 * i + bi, 4 * j + bj] = blk[bi][bj]
    return R


def _det_qmatrix_left_regular(m: "QMatrix") -> float:
    n = m.nrows()
    if m.ncols() != n or n == 0:
        raise ValueError("determinant: quadratische Matrix n>0 nötig")
    return _det_real_square(_qmatrix_to_real_left_regular(m).tolist())


class QMatrix:
    def __init__(self, rows: List[List[Quaternion]]):
        self._rows: List[List[Quaternion]] = [list(r) for r in rows]
        if not self._rows:
            return
        w = len(self._rows[0])
        if any(len(r) != w for r in self._rows):
            raise ValueError("Matrix ist nicht rechteckig")

    def nrows(self) -> int:
        return len(self._rows)

    def ncols(self) -> int:
        return len(self._rows[0]) if self._rows else 0

    def __mul__(self, other: object) -> QMatrix:
        if isinstance(other, QMatrix):
            r, n, m = self.nrows(), self.ncols(), other.ncols()
            if n != other.nrows():
                raise ValueError("Dimensionen für Coend-Zusammensetzung inkompatibel.")
            z = quat(0, 0, 0, 0)
            out: List[List[Quaternion]] = []
            for i in range(r):
                row: List[Quaternion] = []
                for j in range(m):
                    s = z
                    for k in range(n):
                        s = s + self._rows[i][k] * other._rows[k][j]
                    row.append(s)
                out.append(row)
            return QMatrix(out)
        if isinstance(other, (int, float, Fraction)):
            s = Quaternion._as_frac(other)
            return QMatrix(
                [[q * s for q in row] for row in self._rows]
            )
        if isinstance(other, Quaternion):
            return QMatrix(
                [[q * other for q in row] for row in self._rows]
            )
        return NotImplemented

    def __rmul__(self, other: object) -> QMatrix:
        if isinstance(other, Quaternion):
            r, c = self.nrows(), self.ncols()
            return QMatrix(
                [
                    [other * self._rows[i][j] for j in range(c)]
                    for i in range(r)
                ]
            )
        return NotImplemented

    def __add__(self, other: object) -> QMatrix:
        if not isinstance(other, QMatrix):
            return NotImplemented
        if self.nrows() != other.nrows() or self.ncols() != other.ncols():
            raise ValueError("Matrixaddition: gleiche Zeilen- und Spaltenzahl nötig")
        r, c = self.nrows(), self.ncols()
        return QMatrix(
            [
                [self._rows[i][j] + other._rows[i][j] for j in range(c)]
                for i in range(r)
            ]
        )

    def __sub__(self, other: object) -> QMatrix:
        if not isinstance(other, QMatrix):
            return NotImplemented
        if self.nrows() != other.nrows() or self.ncols() != other.ncols():
            raise ValueError("Matrixsubtraktion: gleiche Form nötig")
        r, c = self.nrows(), self.ncols()
        return QMatrix(
            [
                [self._rows[i][j] - other._rows[i][j] for j in range(c)]
                for i in range(r)
            ]
        )

    def frobenius_norm(self) -> float:
        """Euklidische Frobenius-Norm: sqrt(Σ reduzierte Norm je Eintrag)."""
        s = 0.0
        for row in self._rows:
            for q in row:
                s += float(q.reduced_norm())
        return math.sqrt(s)

    def apply_map(
        self, f: Callable[[Quaternion], Quaternion]
    ) -> QMatrix:
        """Eintragweise Anwendung wie Sage ``matrix.apply_map``."""
        c = self.ncols()
        return QMatrix(
            [[f(self._rows[i][j]) for j in range(c)] for i in range(self.nrows())]
        )

    def determinant(self) -> float:
        """
        Determinante der 4n×4n-Linkregular-Darstellung (jede Quaternion-Position
        als 4×4-Block) — Skalarproxi für die Zustandsdichte-Visualisierung;
        kein vollwertiger Quaterninonen-Dieudonné-Determinant.
        """
        return _det_qmatrix_left_regular(self)

    def __str__(self) -> str:
        return "\n".join("[" + "  ".join(str(q) for q in row) + "]" for row in self._rows)

    def is_block_diagonal(
        self, block_sizes: Optional[Sequence[int]] = None
    ) -> bool:
        """
        True, wenn die Matrix quadratisch ist und
        (a) ohne block_sizes: alle Einträge außerhalb der Diagonalen 0, oder
        (b) mit block_sizes [n1, ..., nk], Summe n: alle Blöcke außerhalb
            der diagonalen Blöcke sind 0-Quaternione.
        """
        n = self.nrows()
        if self.ncols() != n:
            return False
        if not block_sizes:
            for i in range(n):
                for j in range(n):
                    if i != j and not self._rows[i][j].is_zero():
                        return False
            return True
        if sum(block_sizes) != n:
            return False
        starts: List[int] = []
        s = 0
        for b in block_sizes:
            starts.append(s)
            s += b
        for bi, bi_size in enumerate(block_sizes):
            for bj, bj_size in enumerate(block_sizes):
                if bi == bj:
                    continue
                r0, r1 = starts[bi], starts[bi] + bi_size
                c0, c1 = starts[bj], starts[bj] + bj_size
                for i in range(r0, r1):
                    for j in range(c0, c1):
                        if not self._rows[i][j].is_zero():
                            return False
        return True

    def submatrix(self, row0: int, col0: int, nrows: int, ncols: int) -> QMatrix:
        """Rechteckige Untermatrix (Sage-API: ab ``(row0, col0)`` mit Größe nrows×ncols)."""
        rmax, cmax = self.nrows(), self.ncols()
        if (
            row0 < 0
            or col0 < 0
            or nrows < 0
            or ncols < 0
            or row0 + nrows > rmax
            or col0 + ncols > cmax
        ):
            raise ValueError("submatrix: Indizes oder Größe ungültig")
        return QMatrix(
            [
                [self._rows[row0 + i][col0 + j] for j in range(ncols)]
                for i in range(nrows)
            ]
        )

    def conjugate(self) -> QMatrix:
        """Eintragsweise Quaternion-Konjugation (ohne Transposition)."""
        c = self.ncols()
        return QMatrix(
            [
                [self._rows[i][j].conjugate() for j in range(c)]
                for i in range(self.nrows())
            ]
        )

    def conjugate_transpose(self) -> QMatrix:
        """Hermitesche Adjungierte M* mit (M*)_{ij} = conj(M_{ji})."""
        m, n = self.nrows(), self.ncols()
        return QMatrix(
            [
                [self._rows[j][i].conjugate() for j in range(m)]
                for i in range(n)
            ]
        )

    def tensor_product(self, other: QMatrix) -> QMatrix:
        """Kronecker-/Tensorprodukt über H (Einträge A_ij · B_kl)."""
        return kronecker_quaternion(self, other)


def kronecker_quaternion(A: QMatrix, B: QMatrix) -> QMatrix:
    """Kroneckerprodukt A ⊗ B mit (A⊗B)_{i p+u, j q+v} = A_ij · B_uv."""
    m, n = A.nrows(), A.ncols()
    p, q = B.nrows(), B.ncols()
    rows: List[List[Quaternion]] = []
    for i in range(m):
        for u in range(p):
            line: List[Quaternion] = []
            for j in range(n):
                for v in range(q):
                    line.append(A._rows[i][j] * B._rows[u][v])
            rows.append(line)
    return QMatrix(rows)


def quaternion_dual(M: QMatrix) -> QMatrix:
    """
    Dualitäts-Funktor (-)* für quaternionische Prozesse: adjungierte Matrix
    (Konjugation + Transposition), analog dualer Kegel c* (Sektion 2.2).
    """
    return M.conjugate_transpose()


def apply_par_product(op_A: QMatrix, op_B: QMatrix) -> QMatrix:
    """
    Par-Operator (A ℘ B) als (A* ⊗ B*)* (Definition 5 / BV-Logik im Manuskript).
    """
    A_dual = quaternion_dual(op_A)
    B_dual = quaternion_dual(op_B)
    tensor_dual = A_dual.tensor_product(B_dual)
    return quaternion_dual(tensor_dual)


def calculate_pair_coherence(
    p1: int, p2: int, product_type: str = "tensor"
) -> float:
    """
    Kohärenz zweier Primzahlen als 1×1-quaternionische Prozess-Objekte: Norm
    (Frobenius) des Ergebnisses.

    ``product_type``: ``'tensor'`` = unabhängiges Kroneckerprodukt; ``'par'`` =
    = ``(A* ⊗ B*)*`` (Par-Produkt, Sekt. 2.2).  Für rein *reelle* 1×1-Einträge
    fällen Tensor und Par mathematisch zusammen, ``Diff`` ist dann 0; mit echten
    Quaternion-Größen in größeren Matrizen entstehen Unterschiede.
    """
    a = quat(p1, 0, 0, 0)
    b = quat(p2, 0, 0, 0)
    A = QMatrix([[a]])
    B = QMatrix([[b]])
    if product_type == "tensor":
        res = A.tensor_product(B)
    elif product_type == "par":
        res = apply_par_product(A, B)
    else:
        raise ValueError("product_type muss 'tensor' oder 'par' sein")
    return res.frobenius_norm()


def get_prime_operator(p: int, t: Optional[float] = None) -> QMatrix:
    """
    Erzeugt einen nicht-trivialen quaternionischen Operator zur Primzahl ``p``:
    Phasen-Signatur in den Nebendiagonaleinträgen, Diagonale ``p + 0i+0j+0k`` / ``1``.

    Wenn ``t`` gesetzt ist (z. B. Riemann-Nullstellen-Höhe), die Oszillation
    ``sin, cos`` mit ``t`` statt ``p`` (kritische Linie / #Energiedoku).
    """
    fp = float(p)
    if t is None:
        ph = (math.sin(fp), math.cos(fp))
    else:
        ft = float(t)
        ph = (math.sin(ft), math.cos(ft))
    phase = quat(p % 4, ph[0], ph[1], 1.0 / fp)
    return QMatrix(
        [
            [quat(p, 0, 0, 0), phase],
            [phase.conjugate(), quat(1, 0, 0, 0)],
        ]
    )


def analyze_bv_coherence(p1: int, p2: int) -> Tuple[float, float]:
    """
    Vergleicht Frobenius-Normen von ``Op1 ⊗ Op2`` (lokal) und Par-Produkt
    ``(Op1* ⊗ Op2*)*`` (verschränkt, BV-Logik).
    """
    op1 = get_prime_operator(p1)
    op2 = get_prime_operator(p2)
    tensor_prod = op1.tensor_product(op2)
    par_prod = apply_par_product(op1, op2)
    return tensor_prod.frobenius_norm(), par_prod.frobenius_norm()


def get_causal_capacity(M: QMatrix) -> float:
    """
    Kausales Kapazitäts-Maß: Shannon-Entropie der auf |λ| normierten Gewichte,
    wobei λ die Eigenwerte der reellen Linkregular-4n×4n-Darstellung von ``M``
    sind (Spektral-Stub analog #Energiedoku / Eigenwert-Log-Normen).
    H = -Σ p log p (natürlicher Logarithmus).
    """
    if M.nrows() != M.ncols() or M.nrows() == 0:
        return 0.0
    R = _qmatrix_to_real_left_regular(M)
    evs = np.linalg.eigvals(R)
    norms = [abs(float(z)) for z in evs if abs(z) > 1e-10]
    total = float(sum(norms))
    if total == 0.0:
        return 0.0
    probs = [n / total for n in norms]
    h = 0.0
    for p in probs:
        if p > 1e-30:
            h -= p * math.log(p)
    return float(h)


def get_signalling_flux(m_4x4: QMatrix) -> float:
    """
    Kausale Korrelation: Frobenius-Norm von (A - D), wobei M = [[A, B], [C, D]]
    in 2×2-Blöcken (Partial-Trace-Stub, semi-causality / Sektion 1 im Paper).
    """
    if m_4x4.nrows() != 4 or m_4x4.ncols() != 4:
        raise ValueError("get_signalling_flux: 4×4-Supermap nötig")
    a = m_4x4.submatrix(0, 0, 2, 2)
    d = m_4x4.submatrix(2, 2, 2, 2)
    return (a - d).frobenius_norm()


def analyze_flux_at_zeros(
    zeros_subset: Union[np.ndarray, List[float], Sequence[float]],
) -> List[Tuple[float, float]]:
    """
    Pro Höhe ``t`` aus den Nullstellen: Par-Produkt der Operatoren
    ``get_prime_operator(17, t)`` und ``get_prime_operator(19, t)``,
    dann ``get_signalling_flux`` auf der 4×4-Supermap.
    """
    results: List[Tuple[float, float]] = []
    for t in zeros_subset:
        tf = float(t)
        op1 = get_prime_operator(17, tf)
        op2 = get_prime_operator(19, tf)
        par_res = apply_par_product(op1, op2)
        results.append((tf, get_signalling_flux(par_res)))
    return results


def analyze_causal_svd(t_instable: float) -> Tuple[List[float], float]:
    """
    SVD der Supermap (Par-Produkt) an der ''instabilen'' Höhe ``t``.

    Baut die gleiche 2×2-Choi-Struktur wie in ``Starke Profunktoren``:
    ``[[p, phase], [phase*, 1]]`` mit
    ``phase = (0.5, t, t^2, sin t)``; Par-Produkt = ``apply_par_product``;

    Führt SVD an der *reellen* 16×16-Linkregular-Matrix aus (4×4 über H
    statt 8×8 komplex, wie in Sage) — singuläre Werte = kausale
    Kanalstärken in dieser Doppelung.

    ``purity`` = ``(Σ σ_k²) / (Σ σ_k)²`` wie im Skizzen-Code (0/1-Interpretation
    nur heuristisch).
    """
    u = float(t_instable)
    phase = quat(0.5, u, u**2, math.sin(u))
    op1 = QMatrix(
        [
            [quat(17, 0, 0, 0), phase],
            [phase.conjugate(), quat(1, 0, 0, 0)],
        ]
    )
    op2 = QMatrix(
        [
            [quat(19, 0, 0, 0), phase],
            [phase.conjugate(), quat(1, 0, 0, 0)],
        ]
    )
    m = apply_par_product(op1, op2)
    r = _qmatrix_to_real_left_regular(m)
    s = np.linalg.svd(r, compute_uv=False)
    s_list = [float(x) for x in s]
    sm = float(sum(s_list))
    if sm < 1e-30:
        return s_list, 0.0
    purity = sum(x * x for x in s_list) / (sm * sm)
    return s_list, float(purity)


def check_no_signalling(
    total_matrix: QMatrix, dim_A: int, dim_B: int
) -> bool:
    """
    Prüft (vereinfacht), ob eine Supermap No-Signalling im Sinne
    teilsystem-entkoppelter Blöcke erlaubt (Tensor-/Block-Diagonal-Analogon;
    vgl. Abschnitt 2.2 / Semi-Lokalisierbarkeit in Ihrem Manuskript).

    total_matrix: Darstellung der Supermap (eindeutige Koordinatenbasis vorausgesetzt).
    dim_A, dim_B: Blockgrößen entlang der Diagonal-Partition, mit dim_A + dim_B = n
        (2×2-Beispiel: dim_A=dim_B=1). Ist n ≠ dim_A + dim_B, fällt die Prüfung auf
        strikte Diagonalität (nur Nichtdiagonal-Einträge 0) zurück.
    """
    n = total_matrix.nrows()
    if total_matrix.ncols() != n or n == 0:
        is_separable = False
    elif dim_A + dim_B == n and dim_A >= 0 and dim_B >= 0:
        is_separable = total_matrix.is_block_diagonal([dim_A, dim_B])
    else:
        is_separable = total_matrix.is_block_diagonal()

    if is_separable:
        print("Bedingung erfüllt: No-Signalling (Raumartige Trennung / Tensor)")
    else:
        print("Bedingung nicht erfüllt: Potenzielles Signalling (Zeitartig / Sequencer)")

    return is_separable


def block_matrix_2x2(
    top_left: QMatrix,
    top_right: QMatrix,
    bottom_left: QMatrix,
    bottom_right: QMatrix,
) -> QMatrix:
    """Fügt vier gleich große quadratische n×n-Blöcke zu einer 2n×2n-Matrix."""
    n = top_left.nrows()
    for name, b in (
        ("top_left", top_left),
        ("top_right", top_right),
        ("bottom_left", bottom_left),
        ("bottom_right", bottom_right),
    ):
        if b.nrows() != b.ncols() or b.nrows() != n:
            raise ValueError(
                f"{name} muss n×n sein und dieselbe Größe wie top_left (n={n})"
            )
    rows: List[List[Quaternion]] = []
    for i in range(n):
        rows.append(
            [top_left._rows[i][j] for j in range(n)]
            + [top_right._rows[i][j] for j in range(n)]
        )
    for i in range(n):
        rows.append(
            [bottom_left._rows[i][j] for j in range(n)]
            + [bottom_right._rows[i][j] for j in range(n)]
        )
    return QMatrix(rows)


def create_quantum_switch_profunctor(
    op_A: QMatrix, op_B: QMatrix, control_qubit_norm: Number
) -> QMatrix:
    """
    Simuliert einen Quanten-Switch als Profunktor.
    Modelliert indefinite Kausalität nach Abschnitt 2.2: zwei Kausalzweige
    (op_A ∘ op_B bzw. op_B ∘ op_A) liegen schematisch blockdiagonal; die
    Nebendiagonal-Blöcke bleiben 0 (Vereinfachung, keine Hadamart-Kreuzkoppelung
    in diesem Stub).

    control_qubit_norm: positive skalare Norm des Kontrollsystems (z. B.
        p_q.reduced_norm() statt p_q.norm()).
    """
    c = (
        control_qubit_norm
        if isinstance(control_qubit_norm, Fraction)
        else Quaternion._as_frac(control_qubit_norm)
    )
    if c == 0:
        raise ValueError("control_qubit_norm ist 0")
    ab = op_A * op_B
    ba = op_B * op_A
    n = ab.nrows()
    if ab.ncols() != n or ba.nrows() != n or ba.ncols() != n:
        raise ValueError("Operatormatrizen inkompatibel (erwartet n×n, passende Multiplikation).")
    z0 = quat(0, 0, 0, 0)
    zero = QMatrix([[z0] * n for _ in range(n)])
    switch_matrix = block_matrix_2x2(ab, zero, zero, ba)
    return switch_matrix * (Fraction(1) / c)


def apply_quaternion_interference(
    switch_matrix: QMatrix, phase_q: Quaternion
) -> QMatrix:
    """
    Mischt die Pfade des Switches über eine eintragweise nicht-kommutative
    „Sliding“-Korrektur: x ↦ x + (x*phase - phase*x) = x + [x, phase]
    (vereinfachter Interferenzterm gemäß Ihrem Modell).
    """
    def comm(x: Quaternion) -> Quaternion:
        return x * phase_q - phase_q * x

    interference = switch_matrix.apply_map(comm)
    return switch_matrix + interference


def analyze_critical_line(
    op_P: QMatrix,
    op_Q: QMatrix,
    phases: Optional[object] = None,
    step_size: float = 0.1,
) -> List[Tuple[float, float]]:
    """
    Untersucht (heuristisch) die reelle Determinante des gestörten Switches
    entlang eines Parameters t (Bild: Bewegung entlang des Imaginärteils);
    reeller Block aus der 4n-Linksdarstellung, Betrag als Dichte-Analogon.

    ``phases`` ist optionaler Platzhalter (API wie Sage) und wird im Stub nicht
    benutzt; die „kritische“ Phase wird pro t als ``(1/2, t, t^2, sin(t))`` in
    die vier Quaternion-Koordinaten gelegt (Modell, nicht 1:1 C mit H).
    """
    _ = phases
    results: List[Tuple[float, float]] = []
    k = 0
    while True:
        t = k * step_size
        if t >= 5.0:
            break
        s_phase = quat(0.5, t, t**2, math.sin(t))
        nrm = s_phase.reduced_norm()
        if nrm == 0:
            k += 1
            continue
        sw = create_quantum_switch_profunctor(op_P, op_Q, nrm)
        sw_int = apply_quaternion_interference(sw, s_phase)
        det = sw_int.determinant()
        results.append((t, abs(det)))
        k += 1
    return results


def quaternion_coend_composition(
    P_matrix: QMatrix, Q_matrix: QMatrix, p_prime: Quaternion
) -> QMatrix:
    """
    Berechnet die Komposition zweier Prozesse (Profunktoren) P und Q
    unter Nutzung des Coend-Kalküls (Definition 7)[cite: 127].

    P_matrix: Repräsentiert P(A, X)
    Q_matrix: Repräsentiert Q(X, B)
    p_prime: Quaternion zur Normierung (reduzierte Norm).
    """
    scale = p_prime.reduced_norm()
    if scale == 0:
        raise ValueError("reduzierte Norm von p_prime ist 0")
    if P_matrix.ncols() != Q_matrix.nrows():
        raise ValueError("Dimensionen für Coend-Zusammensetzung inkompatibel.")

    composed_matrix = P_matrix * Q_matrix
    return composed_matrix * (Fraction(1) / scale)


def _riemann_zeros_npy_path() -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "zeros6.npy")


_riemann_zeros_memmap: Any = None


def get_riemann_zeros() -> np.ndarray:
    """
    Riemann-Nullstellen t > 0 auf der kritischen Linie, aus ``zeros6.npy``
    (read-only memmap, große Datei).
    """
    global _riemann_zeros_memmap
    if _riemann_zeros_memmap is None:
        _riemann_zeros_memmap = np.load(
            _riemann_zeros_npy_path(), mmap_mode="r", allow_pickle=False
        )
    return _riemann_zeros_memmap


def calculate_causal_signature(zero_index: int) -> float:
    """
    Signatur der Supermap-Phase an der Nullstelle ``zero_index`` in ``zeros6.npy``.

    t = γ_{zero_index+1} (Index in der NPY-Datei) speist die Konstruktion
    ``(1/2, t, t^2, sin(t))`` in die Quaternion-Koordinaten; die
    (reduzierte) Norm a^2+b^2+c^2+d^2 dient als skalare Kausal-Signatur
    (Lax–Lax-Duoidalität, vgl. Ihr Manuskript [cite: 2386]).
    """
    t = float(get_riemann_zeros()[zero_index])
    s_phase = quat(0.5, t, t**2, math.sin(t))
    return float(s_phase.reduced_norm())


# Nicht-kommutative Störung (Sage: matrix(H,2,2,[H([0,1,0,0]), H(0), H(0), H([0,0,1,0])]))
P_perturbed = QMatrix(
    [
        [quat(0, 1, 0, 0), quat(0, 0, 0, 0)],
        [quat(0, 0, 0, 0), quat(0, 0, 1, 0)],
    ]
)


def get_signalling_strength(
    t: float, p_op: Optional[QMatrix] = None
) -> float:
    """
    Stärke der kausalen Interferenz: Norm([P, phase]) / ||phase||.

    ``phase`` auf der kritischen Linie als H([0.5, t, t^2, sin(t)]);
    ``P`` = ``P_perturbed`` (gestörter Operator), sofern ``p_op`` nicht
    gesetzt ist. Maß für die Abweichung von reiner No-Signalling-/Lax-Situation
    (Sektion 2.2, E#Energiedoku).
    """
    phase = quat(0.5, t, t**2, math.sin(t))
    P_p = P_perturbed if p_op is None else p_op
    comm = P_p * phase - phase * P_p
    ph_n = math.sqrt(float(phase.reduced_norm()))
    if ph_n < 1e-30:
        return 0.0
    return comm.frobenius_norm() / ph_n


def get_causal_entropy(
    t: float, p_op: Optional[QMatrix] = None
) -> float:
    """
    Heuristisches 'kausales Entropie'-Analog: ln(||[P, phase]||_F + ε).

    ``phase`` wird auf Einheits-Euklidnorm in R^4 gebracht (Sättigungs-Stub);
    ``P`` = ``P_perturbed`` wie im Paper (Kap. 2.2, Abweichung No-Signalling).
    """
    raw_phase = quat(0.5, t, t**2, math.sin(t))
    n = math.sqrt(float(raw_phase.reduced_norm()))
    if n < 1e-30:
        return math.log(1e-20)
    phase = raw_phase * (1.0 / n)
    P_p = P_perturbed if p_op is None else p_op
    comm = P_p * phase - phase * P_p
    c_norm = comm.frobenius_norm()
    return math.log(c_norm + 1e-20)


def get_pure_interference_signal(
    t: float, p_op: Optional[QMatrix] = None
) -> float:
    """
    Isoliertes Interferenz-Analog: ln(||[P, phase]||_F / ||phase||) mit
    phase = (0.5, t, t^2, sin t) in H, P = P_perturbed
    (Choi/Jamiołkowski-Kontext in Ihrem Manuskript, Sekt. 6).
    """
    phase = quat(0.5, t, t**2, math.sin(t))
    phase_norm = math.sqrt(float(phase.reduced_norm()))
    if phase_norm < 1e-30:
        return math.log(1e-20)
    P_p = P_perturbed if p_op is None else p_op
    comm = P_p * phase - phase * P_p
    c_n = comm.frobenius_norm()
    return math.log(c_n / phase_norm + 1e-20)


def run_supermap_interference_spectrum(
    n: int = 2048,
    out_name: str = "causal_supermap_spectrum_riemann.png",
) -> None:
    """
    FFT des entrendeten Log-Interferenzsignals längs der ersten ``n`` Nullstellen;
    markiert ``log(p)/(2π)`` für kleine Primzahlen. Speichert Plot als PNG.
    """
    try:
        from scipy.fft import fft, fftfreq
    except ImportError:
        from numpy.fft import fft, fftfreq

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "matplotlib nicht installiert — Spektrum übersprungen "
            "(pip install matplotlib)."
        )
        return

    zeros = get_riemann_zeros()[:n]
    t_zeros = np.asarray(zeros, dtype=np.float64)
    signal = np.array(
        [get_pure_interference_signal(float(t)) for t in t_zeros],
        dtype=np.float64,
    )
    N = len(signal)
    if N < 2:
        print("Zu wenig Nullstellen für Spektrum.")
        return
    x = np.arange(N, dtype=np.float64)
    trend = np.polyval(np.polyfit(x, signal, 1), x)
    signal_detrended = signal - trend

    yf = fft(signal_detrended)
    xf = fftfreq(N, 1.0)[: N // 2]
    amplitudes = 2.0 / N * np.abs(yf[0 : N // 2])

    primes_to_check = [2, 3, 5, 7, 11]
    log_primes = [float(math.log(p)) for p in primes_to_check]
    two_pi = 2.0 * math.pi

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(xf, amplitudes, label="Kausales Spektrum (Supermap)")
    for i, lp in enumerate(log_primes):
        ax.axvline(
            x=lp / two_pi,
            color="r",
            linestyle="--",
            alpha=0.5,
            label=(f"log({primes_to_check[i]})" if i == 0 else "_nolegend_"),
        )
    ax.set_title("Spektralanalyse der Supermap-Interferenz (#Energiedoku)")
    ax.set_xlabel("Frequenz (entspricht log p / 2pi)")
    ax.set_ylabel("Amplitude (Spektrale Dichte)")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), out_name)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Supermap-Spektrum gespeichert: {out_path}")


def run_signalling_midpoint_histogram(
    sample_size: int = 1000,
    bins: int = 50,
    out_name: str = "causal_interference_riemann.png",
) -> None:
    """
    Vergleicht ``get_signalling_strength`` an den ersten ``sample_size``
    Nullstellen vs. Mittelpunkte aufeinanderfolgender Nullstellen (Kontrolle).
    Speichert Histogramm als PNG neben dieser Datei.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "matplotlib nicht installiert — Histogramm übersprungen "
            "(pip install matplotlib)."
        )
        return

    zeros = get_riemann_zeros()
    test_zeros = zeros[:sample_size]
    mid_points = [
        float((zeros[i] + zeros[i + 1]) / 2.0) for i in range(sample_size - 1)
    ]

    zero_strengths = [get_signalling_strength(float(z)) for z in test_zeros]
    mid_strengths = [get_signalling_strength(m) for m in mid_points]

    avg_zero = float(np.mean(zero_strengths))
    avg_mid = float(np.mean(mid_strengths))
    print(f"Durchschnittliche Signalling-Stärke (Nullstellen): {avg_zero:.4e}")
    print(
        f"Durchschnittliche Signalling-Stärke (Zwischenräume): {avg_mid:.4e}"
    )

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(zero_strengths, bins=bins, alpha=0.5, label="An Nullstellen")
    ax.hist(mid_strengths, bins=bins, alpha=0.5, label="Zwischen Nullstellen")
    ax.set_title("Verteilung der kausalen Interferenz (Riemann #Energiedoku)")
    ax.set_xlabel("Norm-Verhältnis [cite: 10, 118]")
    ax.legend()
    fig.tight_layout()
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), out_name)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Histogramm gespeichert: {out_path}")


def run_spacing_causal_entropy_scatter(
    sample_size: int = 5000,
    out_name: str = "causal_spacing_entropy_riemann.png",
) -> None:
    """
    Normierte Spacings ``delta_n / mean(delta)`` vs. ``get_causal_entropy`` an
    den ersten ``sample_size`` Nullstellen (letzte Nullstelle wegen Koppelung
    an ``diff`` weggelassen). Speichert Streudiagramm als PNG.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "matplotlib nicht installiert — Scatter übersprungen "
            "(pip install matplotlib)."
        )
        return

    zeros = get_riemann_zeros()[:sample_size]
    test_zeros = np.asarray(zeros, dtype=np.float64)
    spacings = np.diff(test_zeros)
    avg_spacing = float(np.mean(spacings))
    if avg_spacing < 1e-30:
        print("Mittlerer Spacing zu klein, Abbruch.")
        return
    norm_spacings = spacings / avg_spacing
    entropies = [get_causal_entropy(float(z)) for z in test_zeros[:-1]]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.scatter(
        norm_spacings,
        entropies,
        alpha=0.3,
        c="purple",
        s=10,
    )
    ax.set_title(
        "Zusammenhang: Riemann-Spacing vs. Kausale Entropie (#Energiedoku)"
    )
    ax.set_xlabel("Normalisierter Abstand (GUE-Analogon)")
    ax.set_ylabel("Kausale Entropie (Log-Norm des Kommutators)")
    ax.grid(True, linestyle="--")
    fig.tight_layout()
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), out_name)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Spacing–Entropie-Scatter gespeichert: {out_path}")


# Beispielanwendung
if __name__ == "__main__":
    p_q = quat(1, 1, 1, 1)
    P = QMatrix(
        [
            [quat(1, 0, 0, 0), quat(0, 0, 0, 0)],
            [quat(0, 0, 0, 0), quat(1, 0, 0, 0)],
        ]
    )
    Q = QMatrix(
        [
            [p_q, quat(0, 0, 0, 0)],
            [quat(0, 0, 0, 0), p_q],
        ]
    )
    
    result = quaternion_coend_composition(P, Q, p_q)
    print(f"Resultierende Supermap (Coend):\n{result}")
    
    # Wie in Sage: H([1/4,1/4,1/4,1/4]) auf der Diagonale, sonst 0
    entry = quat(1, 1, 1, 1) * Fraction(1, 4)
    z = quat(0, 0, 0, 0)
    res_matrix = QMatrix([[entry, z], [z, entry]])
    check_no_signalling(res_matrix, 1, 1)
    
    switch_res = create_quantum_switch_profunctor(P, Q, p_q.reduced_norm())
    print(f"Indefinite Kausalstruktur (Switch):\n{switch_res}")
    
    # Phase i+j+k (Sage: H([0,1,1,1])); führt zu nicht-skalaren Eintragskommutatoren
    phase = quat(0, 1, 1, 1)
    critical_data = analyze_critical_line(P_perturbed, Q, phase)
    print(
        f"Determinanten-Normen entlang der kritischen Linie (Auszug):\n"
        f"{critical_data[:3]}"
    )
    switch_interfered = apply_quaternion_interference(switch_res, phase)
    print(f"Switch mit quaternionischer Interferenz:\n{switch_interfered}")
    # Mit P_perturbed (i, j) auf der Diagonalen ist [Eintrag, phase] i.d.R. ≠ 0
    switch_pert = create_quantum_switch_profunctor(
        P_perturbed, Q, p_q.reduced_norm()
    )
    switch_interfered_pert = apply_quaternion_interference(switch_pert, phase)
    print(
        "Switch (P_perturbed, Q) mit Interferenz (Eintragskommutatoren, hier Diagonale "
        "gegenüber reiner Skalierung geändert):\n"
        f"{switch_interfered_pert}"
    )
    
    signatures = [calculate_causal_signature(i) for i in range(5)]
    print(f"Kausale Signaturen für die ersten Nullstellen: {signatures}")
    
    run_signalling_midpoint_histogram(sample_size=1000, bins=50)
    
    run_spacing_causal_entropy_scatter(sample_size=5000)
    
    run_supermap_interference_spectrum(n=2048)
    
    P_par = QMatrix(
        [
            [quat(1, 0, 0, 0), quat(0, 0, 0, 0)],
            [quat(0, 0, 0, 0), quat(1, 0, 0, 0)],
        ]
    )
    Q_par = QMatrix(
        [
            [quat(0, 1, 1, 1), quat(0, 0, 0, 0)],
            [quat(0, 0, 0, 0), quat(0, 1, 1, 1)],
        ]
    )
    par_result = apply_par_product(P_par, Q_par)
    print(f"Par-Produkt (Maximale Verschränkung):\n{par_result}")
    
    twins = [(3, 5), (5, 7), (11, 13), (17, 19)]
    random_pairs = [(3, 11), (5, 13), (7, 17), (11, 23)]
    print("--- Kohärenz-Check (BV-Logik des Papers) ---")
    for p1, p2 in twins:
        t_coh = calculate_pair_coherence(p1, p2, "tensor")
        p_coh = calculate_pair_coherence(p1, p2, "par")
        print(
            f"Zwilling ({p1}, {p2}): Tensor={t_coh:.2f}, Par={p_coh:.2f} -> "
            f"Diff={p_coh - t_coh:.2f}"
        )
    print("\n--- Vergleich mit Nicht-Zwillingen ---")
    for p1, p2 in random_pairs:
        t_coh = calculate_pair_coherence(p1, p2, "tensor")
        p_coh = calculate_pair_coherence(p1, p2, "par")
        print(
            f"Paar ({p1}, {p2}): Tensor={t_coh:.2f}, Par={p_coh:.2f} -> "
            f"Diff={p_coh - t_coh:.2f}"
        )
    
    bv_p1, bv_p2 = 17, 19
    op1_bv = get_prime_operator(bv_p1)
    op2_bv = get_prime_operator(bv_p2)
    tensor_bv = op1_bv.tensor_product(op2_bv)
    par_bv = apply_par_product(op1_bv, op2_bv)
    t_norm = tensor_bv.frobenius_norm()
    p_norm = par_bv.frobenius_norm()
    print(f"\nBV-Analyse für Zwilling ({bv_p1}, {bv_p2}):")
    print(f"Tensor-Norm (Lokal): {t_norm:.4f}")
    print(f"Par-Norm (Verschränkt): {p_norm:.4f}")
    print(f"Kausale Differenz: {abs(p_norm - t_norm):.4e}")
    t_cap = get_causal_capacity(tensor_bv)
    p_cap = get_causal_capacity(par_bv)
    print("\n--- #Energiedoku Bamberg: Erweitertes Funktional ---")
    print(f"Zwilling ({bv_p1}, {bv_p2})")
    print(f"Kausale Kapazität (Tensor): {t_cap:.6f}")
    print(f"Kausale Kapazität (Par):    {p_cap:.6f}")
    print(f"BV-Differenz (Delta_H):     {abs(p_cap - t_cap):.4e}")
    
    flux_data = analyze_flux_at_zeros(get_riemann_zeros()[:5])
    print("\n--- #Energiedoku: Signalling-Flux an Nullstellen ---")
    for t, f in flux_data:
        print(f"t = {t:.4f} -> Flux: {f:.6e}")
    
    t_svd = 40.9187
    s_vals, p_val = analyze_causal_svd(t_svd)
    print(f"\n--- SVD Analyse (t = {t_svd:.4f}) ---")
    print(
        f"Singulärwerte (erste 4 von {len(s_vals)}; reelle 16×16-Darstellung): "
        f"{s_vals[:4]}"
    )
    print(
        f"Kausale Reinheit (0=indefinit, 1=lokal, heuristisch): {p_val:.6f}"
    )
