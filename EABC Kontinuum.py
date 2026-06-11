"""
EABC-Hamiltonian im Kontinuum: -t Δ ⊗ I_4 + H_eabc (4×4 Kopplungsmatrix).

kwant.continuum.discretize kann keine symbolischen exp(I*phi)-Einträge lambdifizieren.
Für die Gitterdiskretisierung werden daher konkrete Phasen eingesetzt (siehe DEFAULT_PHASES).

Benötigt: numpy, scipy, matplotlib, sympy, kwant.

Zusätzlich (optional): **Minimal-TB** nur EABC-Onsite + Hopping ``-t I`` (ohne Laplace-Onsite),
Ensemble + ``scipy.sparse.linalg.eigsh``: ``finite_square_minimal_random_eabc``,
``compute_eigenvalues_eigsh``, ``ensemble_minimal_random_eabc_beta``.

Riemann-Nullstellen: ``zeros6.npy`` (γ_n = Im(ρ_n)) im selben Ordner wie dieses Skript —
API: ``get_riemann_zeros()`` bzw. ``load_riemann_imag_parts()``.

Räumliche Phasen φ_A(x,y), φ_B(x,y), φ_C(x,y): Python-Funktionen ``(x,y)->float``
(Radiant) mit ``finite_square_eabc_spatial``; nur φ_A zufällig: ``random_spatial_phi_A``;
**φ_A/φ_B/φ_C aus Rastercodierung:** ``spatial_phi_from_delta_field`` /
``spatial_phases_from_three_delta_fields``; gleiche Bedeutung als **phi_K = Delta_K(Feld,q=…)**
(alias ``delta_map_A`` / ``_B`` / ``_C``);
alle drei zufällig: ``random_spatial_eabc_phases``; symbolisch in SymPy über
``SPATIAL_X``/``SPATIAL_Y`` und ``build_hamiltonian_spatial`` (ohne zuverlässiges
``continuum.discretize`` bei ``exp(I φ(x,y))`` in dieser Kwant-Version).
"""

from __future__ import annotations

import os
from collections.abc import Callable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import kwant
import scipy.sparse.linalg as sla
from scipy.interpolate import UnivariateSpline
import kwant.continuum as continuum
import kwant.wraparound
import sympy as sp

# SymPy-Ortskoordinaten für φ_A(x,y), φ_B(x,y), φ_C(x,y) (siehe build_hamiltonian_spatial)
SPATIAL_X, SPATIAL_Y = sp.symbols("x y")


def build_hamiltonian():
    """Sympy-Matrix H(k_x, k_y; t, phi_A, phi_B, phi_C)."""
    k_x, k_y = sp.symbols("k_x k_y")
    t = sp.Symbol("t", positive=True)
    phi_A, phi_B, phi_C = sp.symbols("phi_A phi_B phi_C")

    alpha = sp.exp(sp.I * phi_A)
    beta = sp.exp(sp.I * phi_B)
    gamma = sp.exp(sp.I * phi_C)

    h_eabc = sp.Matrix(
        [
            [0, alpha, beta, gamma],
            [alpha, 0, gamma, beta],
            [beta, gamma, 0, alpha],
            [gamma, beta, alpha, 0],
        ]
    )
    h_kin = -t * (k_x**2 + k_y**2) * sp.eye(4)
    return h_kin + h_eabc, (k_x, k_y, t, phi_A, phi_B, phi_C)


def build_dirac_hamiltonian():
    """
    Massenloses Dirac-ähnliches 8×8-Modell im Impulsraum:

    ``H = k_x (σ_x ⊗ I_4) + k_y (σ_y ⊗ I_4) + (I_2 ⊗ H_eabc)``.

    Die ersten beiden Terme koppeln zwei vierkomponentige „Fermion“-Blöcke linear
    in ``k``; ``H_eabc`` ist dieselbe 4×4-EABC-Matrix wie in ``build_hamiltonian``
    (ohne Schrödinger-Laplace ``-t k^2``).
    """
    k_x, k_y = sp.symbols("k_x k_y")
    phi_A, phi_B, phi_C = sp.symbols("phi_A phi_B phi_C")

    alpha = sp.exp(sp.I * phi_A)
    beta = sp.exp(sp.I * phi_B)
    gamma = sp.exp(sp.I * phi_C)
    h_eabc = sp.Matrix(
        [
            [0, alpha, beta, gamma],
            [alpha, 0, gamma, beta],
            [beta, gamma, 0, alpha],
            [gamma, beta, alpha, 0],
        ]
    )

    sigma_x = sp.Matrix([[0, 1], [1, 0]])
    sigma_y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    i2 = sp.eye(2)
    i4 = sp.eye(4)

    h_dirac = (
        sp.kronecker_product(sigma_x, i4) * k_x
        + sp.kronecker_product(sigma_y, i4) * k_y
        + sp.kronecker_product(i2, h_eabc)
    )
    return h_dirac, (k_x, k_y, phi_A, phi_B, phi_C)


def build_hamiltonian_spatial(
    phi_A: sp.Expr,
    phi_B: sp.Expr,
    phi_C: sp.Expr,
) -> tuple[sp.Matrix, tuple]:
    """
    Schrödinger + EABC mit **örtlich variierenden** Phasen φ_A(x,y), φ_B(x,y), φ_C(x,y).

    ``phi_*`` sind Sympy-Ausdrücke in ``SPATIAL_X``, ``SPATIAL_Y``. Nutzbar z.B. für
    algebraische Reduktion oder punktweise Auswertung mit ``.subs({x: x0, y: y0})``.

    Für **Kwant-Gitter** mit variablen Phasen siehe ``finite_square_eabc_spatial``
    (Callables ``(x,y) -> φ``); ``continuum.discretize`` liefert hier oft fehlerhafte
    Onsite-Funktionen, sobald ``exp(I·φ(x,y))`` vorkommt.
    """
    k_x, k_y = sp.symbols("k_x k_y")
    t = sp.Symbol("t", positive=True)
    alpha = sp.exp(sp.I * phi_A)
    beta = sp.exp(sp.I * phi_B)
    gamma = sp.exp(sp.I * phi_C)
    h_eabc = sp.Matrix(
        [
            [0, alpha, beta, gamma],
            [alpha, 0, gamma, beta],
            [beta, gamma, 0, alpha],
            [gamma, beta, alpha, 0],
        ]
    )
    h_kin = -t * (k_x**2 + k_y**2) * sp.eye(4)
    return h_kin + h_eabc, (k_x, k_y, t, SPATIAL_X, SPATIAL_Y)


def build_dirac_hamiltonian_spatial(
    phi_A: sp.Expr,
    phi_B: sp.Expr,
    phi_C: sp.Expr,
) -> tuple[sp.Matrix, tuple]:
    """Dirac-8×8 mit φ(x,y); rein symbolisch, siehe ``build_hamiltonian_spatial``."""
    k_x, k_y = sp.symbols("k_x k_y")
    alpha = sp.exp(sp.I * phi_A)
    beta = sp.exp(sp.I * phi_B)
    gamma = sp.exp(sp.I * phi_C)
    h_eabc = sp.Matrix(
        [
            [0, alpha, beta, gamma],
            [alpha, 0, gamma, beta],
            [beta, gamma, 0, alpha],
            [gamma, beta, alpha, 0],
        ]
    )
    sigma_x = sp.Matrix([[0, 1], [1, 0]])
    sigma_y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    i2 = sp.eye(2)
    i4 = sp.eye(4)
    h_dirac = (
        sp.kronecker_product(sigma_x, i4) * k_x
        + sp.kronecker_product(sigma_y, i4) * k_y
        + sp.kronecker_product(i2, h_eabc)
    )
    return h_dirac, (k_x, k_y, SPATIAL_X, SPATIAL_Y)


def dirac_discretized_builder(*, grid: float = 1.0, **phases: float):
    """
    Entspricht ``builder = kwant.continuum.discretize(H_dirac)`` nach Einsetzen
    numerischer ``phi_A, phi_B, phi_C`` (symbolische ``exp(I phi)`` lambdifiziert Kwant nicht).

    Rückgabe: ``DiscretizedBuilder`` mit **8 NORBs** pro Gitterplatz. Für den
    2D-Bulk: ``kwant.wraparound.wraparound(builder).finalized()``.
    """
    h_dirac_sym, _ = build_dirac_hamiltonian()
    h_num = sp.N(_subs_named(h_dirac_sym, **phases))
    return continuum.discretize(h_num, grid=grid)


DEFAULT_T = 1.0
DEFAULT_PHASES = dict(phi_A=0.2, phi_B=0.5, phi_C=-0.3)

# Die vier Eck-/Kanal-Labels E,A,B,C in Viertelaquivalenz φ ∈ {0, π/2, π, 3π/2} (Radiant).
# Für codierte Raster können z.B. Indices 0..3 auf diese Zeichen dann auf ``PHASE_MAP`` gelegt werden.
PHASE_MAP: dict[str, float] = {
    "E": 0.0,
    "A": 0.5 * np.pi,
    "B": 1.0 * np.pi,
    "C": 1.5 * np.pi,
}

_DELTA_SYMBOL_LOOKUP: tuple[str, str, str, str] = ("E", "A", "B", "C")


def delta_symbol_at(i: int, j: int, channel: int) -> str:
    """
    Platzhalter-Δ für ein Symbol aus ``{'E','A','B','C'}`` am Gitterpunkt ``(i,j)``
    zur Kopplungslinie „channel“. **Deterministisch**, vorläufig nichttrivial zerhackt –
    gegen echte Δ-Physik austauschen.

    ``channel``: ``0``, ``1`` oder ``2`` — Zuordnung zu den drei Kantenphasen phi_A,
    phi_B, phi_C.
    """
    ii = int(i)
    jj = int(j)
    ch = int(channel)
    # Hash-Mischkonstanten (Platzhalter; nicht kryptografisch gedacht).
    key = (ii * 73856093) ^ (jj * 19349663) ^ (ch * 83492791)
    return _DELTA_SYMBOL_LOOKUP[key % 4]


def _delta_phase_fields_from_sym(
    Lx: int,
    Ly: int,
    sym: Callable[[int, int, int], str],
) -> tuple[
    Callable[[float, float], float],
    Callable[[float, float], float],
    Callable[[float, float], float],
]:
    """
    Intern: dreifaches Phasenraster aus Funktion ``sym(i, j, c) -> 'E'|'A'|'B'|'C'``.
    """
    lx = max(0, int(Lx))
    ly = max(0, int(Ly))
    field_a = np.zeros((lx, ly))
    field_b = np.zeros((lx, ly))
    field_c = np.zeros((lx, ly))

    for ii in range(lx):
        for jj in range(ly):
            sa = sym(ii, jj, 0)
            sb = sym(ii, jj, 1)
            sc = sym(ii, jj, 2)
            field_a[ii, jj] = PHASE_MAP[sa]
            field_b[ii, jj] = PHASE_MAP[sb]
            field_c[ii, jj] = PHASE_MAP[sc]

    def _reader(field: np.ndarray) -> Callable[[float, float], float]:
        def fn(x: float, y: float) -> float:
            i = int(round(float(x)))
            j = int(round(float(y)))
            if 0 <= i < lx and 0 <= j < ly:
                return float(field[i, j])
            return 0.0

        return fn

    return _reader(field_a), _reader(field_b), _reader(field_c)


def delta_phase_fields(
    Lx: int,
    Ly: int,
) -> tuple[
    Callable[[float, float], float],
    Callable[[float, float], float],
    Callable[[float, float], float],
]:
    """
    Raster der drei Kantenphasen (Radiant) aus ``delta_symbol_at(i, j, ch)`` und
    ``PHASE_MAP`` über alle ``i``, ``j``. Gibt drei Callables für
    ``finite_square_eabc_spatial`` zurück (wie andere räumliche Phasen ausgewertet).

    Pro Platz: Kanäle ``0``, ``1``, ``2`` füllen die drei Felder (phi_A-, phi_B- und phi_C-Eintrag).

    Außerhalb ``[0, Lx) x [0, Ly)``: die Callables liefern ``0.0``.
    """
    return _delta_phase_fields_from_sym(
        Lx,
        Ly,
        lambda i, j, c: delta_symbol_at(int(i), int(j), int(c)),
    )


def delta_phase_fields_with_offset(
    Lx: int,
    Ly: int,
    offx: int,
    offy: int,
) -> tuple[
    Callable[[float, float], float],
    Callable[[float, float], float],
    Callable[[float, float], float],
]:
    """
    Wie ``delta_phase_fields``, Zuordnung mit ``delta_symbol_at(i + offx, j + offy, c)``.

    Negative ``offx`` / ``offy`` sind erlaubt.
    """
    ox = int(offx)
    oy = int(offy)

    def sym(ii: int, jj: int, ch: int) -> str:
        return delta_symbol_at(ii + ox, jj + oy, int(ch))

    return _delta_phase_fields_from_sym(Lx, Ly, sym)


# Unabhängige Spektren / Realisierungen für Levelstatistik: typisch **50–100**;
# Standard = Mitte dieses Bereichs (Ensemble-Funktionen).
DEFAULT_ENSEMBLE_SPECTRA: int = 75
ENSEMBLE_SPECTRA_RECOMMENDED: tuple[int, int] = (50, 100)

# Riemann-Nullstellen (Projekt „zeros6“): γ_n entlang Re(s)=½
_PACKAGE_DIR = Path(__file__).resolve().parent
ZEROS6_NPY: Path = _PACKAGE_DIR / "zeros6.npy"
ZEROS6_NPZ: Path = _PACKAGE_DIR / "zeros6.npz"

_RIEMANN_ZEROS_MMAP: np.ndarray | None = None


def _subs_named(expr: sp.Expr, **name_to_value: float) -> sp.Expr:
    """Ersetzt Symbole nach Namen (vermeidet Kollisionen zweier Symbole 't')."""
    mapping = {}
    for s in expr.free_symbols:
        if s.name in name_to_value:
            mapping[s] = name_to_value[s.name]
    return expr.subs(mapping)


def make_square_shape(Lx: int, Ly: int):
    """
    shape-Funktion für builder.fill(...).

    WICHTIG: fill übergibt hier ein Site-Objekt s, nicht das Tupel pos.
    (lat.shape(shape, origin) dagegen ruft shape(pos) mit Koordinaten auf.)
    """

    def shape(site):
        x, y = site.pos
        return 0 <= x < Lx and 0 <= y < Ly

    return shape


def eabc_matrix(*, phi_A: float, phi_B: float, phi_C: float) -> np.ndarray:
    """Numerische 4×4 EABC-Kopplung (gleiche Konvention wie build_hamiltonian)."""
    alpha = np.exp(1j * phi_A)
    beta = np.exp(1j * phi_B)
    gamma = np.exp(1j * phi_C)
    return np.array(
        [
            [0, alpha, beta, gamma],
            [alpha, 0, gamma, beta],
            [beta, gamma, 0, alpha],
            [gamma, beta, alpha, 0],
        ],
        dtype=complex,
    )


def random_eabc_onsite_consistent(
    Lx: int,
    Ly: int,
    *,
    seed: int | None = None,
) -> Callable[..., np.ndarray]:
    """
    Liefert eine Funktion ``onsite(site, *args, **kwargs) -> H_eabc`` (nur der 4×4-Block),
    wie in typischen Kwant-Snippets — aber mit **festen** Zufallsphasen pro Gitterplatz.

    **Anti-Pattern** (nicht so bauen)::

        def onsite(site, params):
            phi_A = np.random.uniform(0, 2 * np.pi)  # bei jedem Aufruf neu!

    ``onsite`` wird beim Finalisieren / bei ``hamiltonian_submatrix`` ggf. mehrfach
    aufgerufen; dann sind die Phasen nicht mehr pro Site konstant und der Operator
    ist nicht wohldefiniert. Hier werden drei ``(Lx, Ly)``-Felder einmal mit
    ``numpy.random.default_rng(seed)`` gefüllt (gleiche Konvention wie
    ``random_spatial_eabc_phases``).

    Für das **vollständige** TB (Laplace + EABC): ``finite_square_eabc_random`` /
    ``finite_square_eabc_spatial``.
    """
    rng = np.random.default_rng(seed)
    field_A = rng.uniform(0.0, 2.0 * np.pi, size=(Lx, Ly))
    field_B = rng.uniform(0.0, 2.0 * np.pi, size=(Lx, Ly))
    field_C = rng.uniform(0.0, 2.0 * np.pi, size=(Lx, Ly))

    def onsite(site, *args, **kwargs) -> np.ndarray:
        del args, kwargs  # Kwant übergibt ggf. zusätzliche Argumente
        x, y = site.pos
        i = int(round(float(x)))
        j = int(round(float(y)))
        if not (0 <= i < Lx and 0 <= j < Ly):
            raise ValueError(f"Site außerhalb des {Lx}×{Ly}-Rechtecks: pos={site.pos!r}")
        return eabc_matrix(
            phi_A=float(field_A[i, j]),
            phi_B=float(field_B[i, j]),
            phi_C=float(field_C[i, j]),
        )

    return onsite


def finite_square_manual_lat(
    *,
    Lx: int,
    Ly: int,
    t: float = 1.0,
    a: float = 1.0,
    **phases: float,
) -> kwant.system.FiniteSystem:
    """
    Gleiches TB wie ``continuum.discretize(-t(k_x²+k_y²)I + H_eabc, grid=a)``,
    aber mit explizitem ``kwant.lattice.square(a, norbs=4)``.

    - ``lat.shape(shape, origin)`` übergibt an ``shape`` ein Koordinaten-**pos**-Tupel.
    - **Kein** ``sys[site] = builder``: On-Site ist eine **4×4-Matrix**.
    - **Nachbarn** sind nicht „automatisch“: ``sys[lat.neighbors()] = …`` setzt
      die kinetischen NN-Hoppings (hier +t/a² auf I₄).

    Häufiger Fehler: Onsite +4t/a² und Hop -t/a² — das ist **vorzeichenverkehrt**
    zu Kwant's Standarddiskretisierung von ``-(k_x²+k_y²)``.
    """
    lat = kwant.lattice.square(a, norbs=4)
    kin_onsite = (-4.0 * t / a**2) * np.eye(4, dtype=complex)
    kin_hop = (t / a**2) * np.eye(4, dtype=complex)
    onsite = kin_onsite + eabc_matrix(
        phi_A=phases["phi_A"],
        phi_B=phases["phi_B"],
        phi_C=phases["phi_C"],
    )

    def shape(pos):
        x, y = pos
        return 0 <= x < Lx and 0 <= y < Ly

    sys = kwant.Builder()
    sys[lat.shape(shape, (0, 0))] = onsite
    sys[lat.neighbors()] = kin_hop
    return sys.finalized()


def random_spatial_eabc_phases(
    Lx: int,
    Ly: int,
    *,
    seed: int | None = None,
) -> tuple[
    Callable[[float, float], float],
    Callable[[float, float], float],
    Callable[[float, float], float],
]:
    """
    Pro Platz unabhängig **φ_A, φ_B, φ_C ~ Uniform(0, 2π)** (drei Felder ``(Lx, Ly)``).

    Auswertung an Gitterkoordinaten über
    ``int(round(x))``, ``int(round(y))``. Außerhalb des Rechtecks wird ``0.0`` geliefert.

    ``seed``: Zahl → reproduzierbar mit ``numpy.random.default_rng``; ``None`` → jedes Mal neu.
    """
    rng = np.random.default_rng(seed)
    field_A = rng.uniform(0.0, 2.0 * np.pi, size=(Lx, Ly))
    field_B = rng.uniform(0.0, 2.0 * np.pi, size=(Lx, Ly))
    field_C = rng.uniform(0.0, 2.0 * np.pi, size=(Lx, Ly))

    def _from_field(field: np.ndarray) -> Callable[[float, float], float]:
        def phi(x: float, y: float) -> float:
            i = int(round(float(x)))
            j = int(round(float(y)))
            if 0 <= i < Lx and 0 <= j < Ly:
                return float(field[i, j])
            return 0.0

        return phi

    return _from_field(field_A), _from_field(field_B), _from_field(field_C)


def random_spatial_phi_A(
    Lx: int,
    Ly: int,
    *,
    seed: int | None = None,
) -> Callable[[float, float], float]:
    """
    **φ_A(x,y) ~ Uniform(0, 2π)** i.i.d. pro Gitterplatz (``numpy.random.default_rng.uniform``);
    φ_B und φ_C bleiben frei wählbar (z.B. ``spatial_phase_constant`` oder eigene Callables).

    Beispiel::

        phA = random_spatial_phi_A(30, 30, seed=1)
        fin = finite_square_eabc_spatial(
            Lx=30, Ly=30,
            phi_A=phA,
            phi_B=spatial_phase_constant(0.5),
            phi_C=spatial_phase_constant(DEFAULT_PHASES["phi_C"]),
        )
    """
    rng = np.random.default_rng(seed)
    field = rng.uniform(0.0, 2.0 * np.pi, size=(Lx, Ly))

    def phi_A(x: float, y: float) -> float:
        i = int(round(float(x)))
        j = int(round(float(y)))
        if 0 <= i < Lx and 0 <= j < Ly:
            return float(field[i, j])
        return 0.0

    return phi_A


def spatial_phase_constant(phi: float) -> Callable[[float, float], float]:
    """Callable ``(x,y) -> φ`` (ortsunabhängige Phase im Bogenmaß)."""

    def f(_x: float, _y: float) -> float:
        return float(phi)

    return f


def spatial_phi_from_delta(
    delta: Callable[[float, float], float],
    phase_map: Callable[[float], float],
) -> Callable[[float, float], float]:
    """
    **φ(x,y) = phase_map(delta(x,y))** für eine beliebige Kante (**φ_A**, **φ_B** oder **φ_C**).

    ``delta`` liefert eine codierte Größe (Index, …), ausgewertet an den Gitterkoordinaten
    wie bei ``random_spatial_phi_A`` (Kwant-``site.pos``). ``phase_map`` bildet nach
    **Radiant** ab — typisch ``phase_map_uniform_q(q)``.
    """

    def phi_xy(x: float, y: float) -> float:
        return float(phase_map(delta(x, y)))

    return phi_xy


def phase_map_uniform_q(q: int) -> Callable[[float], float]:
    """
    Diskrete Klassen ``k ∈ {0,…,q-1}`` (als ``float`` übergeben) → gleichmäßige Phasen
    ``2π·k/q``. Nichtganzzahlige ``delta`` werden auf die nächste Klasse gerundet.
    """
    q = int(q)
    if q < 1:
        raise ValueError("phase_map_uniform_q: q muss mindestens 1 sein")

    def f(delta: float) -> float:
        k = int(round(float(delta))) % q
        return (2.0 * np.pi / q) * float(k)

    return f


def spatial_phi_from_delta_field(
    delta_field: np.ndarray,
    *,
    q: int = 4,
) -> Callable[[float, float], float]:
    """
    **φ(x,y) = phase_map(Δ(x,y))** mit Raster ``delta_field[i,j]``, gleiche Auswertung
    wie bei ``random_spatial_phi_A`` (``round`` auf Gitterplatzzentren). Gilt für
    **φ_A**, **φ_B** und **φ_C** — nur andere Eingabetensoren Δ.

    ``delta_field.shape`` ist ``(Lx, Ly)``.

    Für Lesarten ``phi_K = delta_map_K(...)`` oder **``phi_K = Delta_K(...)``** siehe
    ``delta_map_A`` / ``_B`` / ``_C`` bzw. Aliasse ``Delta_A``, ``Delta_B``, ``Delta_C``.
    """
    arr = np.asarray(delta_field, dtype=float)
    if arr.ndim != 2:
        raise ValueError("spatial_phi_from_delta_field: delta_field muss 2-dimensional sein")
    Lx_f, Ly_f = int(arr.shape[0]), int(arr.shape[1])
    phase_map = phase_map_uniform_q(q)

    def delta_xy(x: float, y: float) -> float:
        i = int(round(float(x)))
        j = int(round(float(y)))
        if 0 <= i < Lx_f and 0 <= j < Ly_f:
            return float(arr[i, j])
        return 0.0

    return spatial_phi_from_delta(delta_xy, phase_map)


def spatial_phi_A_from_delta_field(
    delta_field: np.ndarray,
    *,
    q: int = 4,
) -> Callable[[float, float], float]:
    """
    Alias für ``spatial_phi_from_delta_field`` (Kompatibilität: nur Delta_A gekoppelt).

    Für **alle drei Kanten** aus Rasterfeldern: ``spatial_phi_from_delta_field``
    drei Mal aufrufen oder ``spatial_phases_from_three_delta_fields``.
    """
    return spatial_phi_from_delta_field(delta_field, q=q)


def delta_map_A(
    delta_field: np.ndarray,
    *,
    q: int = 4,
) -> Callable[[float, float], float]:
    """phi_A(x,y) aus dem Rasterfeld Delta_A; Kurzschreibweise für ``spatial_phi_from_delta_field``."""
    return spatial_phi_from_delta_field(delta_field, q=q)


def delta_map_B(
    delta_field: np.ndarray,
    *,
    q: int = 4,
) -> Callable[[float, float], float]:
    """phi_B(x,y) aus Delta_B (gleiche Logik wie ``delta_map_A``)."""
    return spatial_phi_from_delta_field(delta_field, q=q)


def delta_map_C(
    delta_field: np.ndarray,
    *,
    q: int = 4,
) -> Callable[[float, float], float]:
    """phi_C(x,y) aus Delta_C (gleiche Logik wie ``delta_map_A``)."""
    return spatial_phi_from_delta_field(delta_field, q=q)


# Notation **phi_K = Delta_K(Raster, q=…)** (mathematisch oft ``\Delta_K``); identisch mit ``delta_map_*``.
Delta_A = delta_map_A
Delta_B = delta_map_B
Delta_C = delta_map_C


def spatial_phases_from_three_delta_fields(
    delta_A: np.ndarray,
    delta_B: np.ndarray,
    delta_C: np.ndarray,
    *,
    q_A: int = 4,
    q_B: int = 4,
    q_C: int = 4,
) -> tuple[
    Callable[[float, float], float],
    Callable[[float, float], float],
    Callable[[float, float], float],
]:
    """
    (phi_A, phi_B, phi_C) aus drei Feldern Delta_A, Delta_B, Delta_C, jeweils mit
    ``spatial_phi_from_delta_field``.

    Alle Felder haben dieselbe Shape ``(Lx, Ly)``.
    """
    a = np.asarray(delta_A)
    b = np.asarray(delta_B)
    c = np.asarray(delta_C)
    if a.shape != b.shape or a.shape != c.shape:
        raise ValueError(
            "spatial_phases_from_three_delta_fields: delta_A/B/C gleiche Shape nötig, "
            f"erhalten {a.shape}, {b.shape}, {c.shape}"
        )
    return (
        spatial_phi_from_delta_field(delta_A, q=q_A),
        spatial_phi_from_delta_field(delta_B, q=q_B),
        spatial_phi_from_delta_field(delta_C, q=q_C),
    )


def finite_square_eabc_random(
    *,
    Lx: int,
    Ly: int,
    seed: int | None = None,
    t: float = 1.0,
    a: float = 1.0,
) -> kwant.system.FiniteSystem:
    """Kurzform: ``finite_square_eabc_spatial`` mit ``random_spatial_eabc_phases(Lx, Ly, seed=...)``."""
    phA, phB, phC = random_spatial_eabc_phases(Lx, Ly, seed=seed)
    return finite_square_eabc_spatial(
        Lx=Lx, Ly=Ly, phi_A=phA, phi_B=phB, phi_C=phC, t=t, a=a
    )


def eabc_matrix_xy(
    x: float,
    y: float,
    *,
    phi_A: Callable[[float, float], float],
    phi_B: Callable[[float, float], float],
    phi_C: Callable[[float, float], float],
) -> np.ndarray:
    """4×4 EABC-Kopplung am Punkt ``(x,y)``; ``phi_*`` liefern Phasen im **Bogenmaß**."""
    return eabc_matrix(
        phi_A=phi_A(x, y),
        phi_B=phi_B(x, y),
        phi_C=phi_C(x, y),
    )


def finite_square_eabc_spatial(
    *,
    Lx: int,
    Ly: int,
    phi_A: Callable[[float, float], float],
    phi_B: Callable[[float, float], float],
    phi_C: Callable[[float, float], float],
    t: float = 1.0,
    a: float = 1.0,
) -> kwant.system.FiniteSystem:
    """
    Quadratgitter Lx×Ly mit **φ_A(x,y), φ_B(x,y), φ_C(x,y)** (reelle Callables).

    Entspricht dem gleichen kinetischen TB wie ``finite_square_manual_lat``, aber
    ortsabhängigem EABC-Onsite. ``(x,y)`` sind die Kwant-Gitterkoordinaten von
    ``site.pos`` / ``shape(pos)``.

    Beispiel::

        import numpy as np
        fin = finite_square_eabc_spatial(
            Lx=20, Ly=20,
            phi_A=lambda x, y: 0.1 * x,
            phi_B=lambda x, y: 0.1 * y,
            phi_C=lambda x, y: 0.05 * (x + y),
            t=1.0,
        )

        # nur φ_A zufällig; φ_B, φ_C konstant:
        # phA = random_spatial_phi_A(20, 20, seed=42)
        # fin = finite_square_eabc_spatial(Lx=20, Ly=20, phi_A=phA,
        #     phi_B=spatial_phase_constant(0.5), phi_C=spatial_phase_constant(-0.3))
        # alle drei zufällig (reproduzierbar mit seed=42):
        # phA, phB, phC = random_spatial_eabc_phases(20, 20, seed=42)
        # fin = finite_square_eabc_spatial(Lx=20, Ly=20, phi_A=phA, phi_B=phB, phi_C=phC)
        # oder: fin = finite_square_eabc_random(Lx=20, Ly=20, seed=42)
    """
    lat = kwant.lattice.square(a, norbs=4)
    kin_onsite = (-4.0 * t / a**2) * np.eye(4, dtype=complex)
    kin_hop = (t / a**2) * np.eye(4, dtype=complex)

    def onsite(site):
        x, y = site.pos
        return kin_onsite + eabc_matrix_xy(
            x, y, phi_A=phi_A, phi_B=phi_B, phi_C=phi_C
        )

    def shape(pos):
        x, y = pos
        return 0 <= x < Lx and 0 <= y < Ly

    sys = kwant.Builder()
    sys[lat.shape(shape, (0, 0))] = onsite
    sys[lat.neighbors()] = kin_hop
    return sys.finalized()


def finite_square_minimal_random_eabc(
    Lx: int,
    Ly: int,
    *,
    t: float = 1.0,
    norbs: int = 4,
    a: float = 1.0,
    seed: int | None = None,
) -> kwant.system.FiniteSystem:
    """
    **Minimal-Modell** (wie typisches Nutzer-Skript): pro Platz **nur** zufällige
    4×4-EABC-Matrix (φ_A,φ_B,φ_C je einmal pro Site gezogen), NN-Hopping ``-t I_4``.

    Kein kinetischer Onsite ``-4t/a²`` — das unterscheidet sich bewusst von
    ``finite_square_eabc_spatial`` / ``finite_square_eabc_random`` (dort Laplace+EABC).

    ``seed``: Realisierung reproduzierbar; ``None``: neue Störung.
    """
    rng = np.random.default_rng(seed)
    lat = kwant.lattice.square(a=a, norbs=norbs)
    sys = kwant.Builder()
    for x in range(int(Lx)):
        for y in range(int(Ly)):
            sys[lat(x, y)] = eabc_matrix(
                phi_A=float(rng.uniform(0.0, 2.0 * np.pi)),
                phi_B=float(rng.uniform(0.0, 2.0 * np.pi)),
                phi_C=float(rng.uniform(0.0, 2.0 * np.pi)),
            )

    hop = -float(t) * np.eye(norbs, dtype=complex)
    sys[lat.neighbors()] = hop
    return sys.finalized()


def compute_eigenvalues_eigsh(
    sys: kwant.system.FiniteSystem,
    *,
    k: int = 200,
    sigma: float = 0.0,
) -> np.ndarray:
    """
    ``k`` Eigenwerte nahe ``sigma`` (shift-invert), sortiert reell — für große Gitter,
    wenn nicht das volle Spektrum nötig ist.

    ``k`` wird bei kleiner Matrix automatisch auf ``n-1`` begrenzt.
    """
    H = sys.hamiltonian_submatrix(sparse=True)
    n = int(H.shape[0])
    k_eff = min(max(1, int(k)), max(1, n - 1))
    ev = sla.eigsh(
        H,
        k=k_eff,
        sigma=float(sigma),
        return_eigenvectors=False,
    )
    return np.sort(np.asarray(ev.real, dtype=float))


def ensemble_minimal_random_eabc_beta(
    n_runs: int = DEFAULT_ENSEMBLE_SPECTRA,
    *,
    L: int = 20,
    t: float = 1.0,
    k_eigs: int = 200,
    sigma: float = 0.0,
    spline_s: float | None = None,
    run_seed0: int = 0,
    s_max: float = 0.5,
    bins: int = 30,
    verbose: bool = False,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    """
    Ensemble wie im Nutzer-Skript: pro Lauf ``finite_square_minimal_random_eabc``,
    ``compute_eigenvalues_eigsh``, ``level_spacings``, ``estimate_beta``.

    - Erstes Tuple-Element: alle Abstände aller Läufe (wie ``all_spacings`` mit
      ``all_spacings.extend(s_run)`` pro Lauf).
    - ``betas``: ein ``estimate_beta``-Wert pro Lauf.
    - ``stats``: u.a. ``mean_beta``, ``std_beta``, ``beta_pooled`` (Fit auf gepoolte Abstände).

    **Hinweis:** Nur ``k_eigs`` Modi nahe ``sigma`` — Entfaltung bezieht sich dann
    auf diesen Ausschnitt, nicht auf das vollständige Spektrum.

    Für die Statistik oft **50–100** Spektren (``ENSEMBLE_SPECTRA_RECOMMENDED``);
    Voreinstellung ``n_runs=DEFAULT_ENSEMBLE_SPECTRA`` (75).
    """
    n_runs = int(n_runs)
    if n_runs < 1:
        raise ValueError("n_runs muss mindestens 1 sein")
    all_spacings: list[float] = []
    betas: list[float] = []
    for run in range(n_runs):
        syst = finite_square_minimal_random_eabc(
            Lx=L, Ly=L, t=t, seed=run_seed0 + run
        )
        E = compute_eigenvalues_eigsh(syst, k=k_eigs, sigma=sigma)
        s_run = level_spacings(E, spline_s=spline_s)
        b = estimate_beta(s_run, s_max=s_max, bins=bins)
        betas.append(b)
        all_spacings.extend(s_run)
        if verbose:
            print(f"Run {run:02d}  beta ≈ {b:.3f}")
    arr_s = np.asarray(all_spacings, dtype=float)
    arr_b = np.asarray(betas, dtype=float)
    stats: dict[str, float] = {
        "mean_beta": float(np.nanmean(arr_b)),
        "std_beta": float(np.nanstd(arr_b)),
        "beta_pooled": float(estimate_beta(arr_s, s_max=s_max, bins=bins)),
    }
    return arr_s, arr_b, stats


def finite_square_from_continuum(
    h_sym: sp.Matrix,
    *,
    Lx: int,
    Ly: int,
    t: float,
    grid: float = 1.0,
    **phases: float,
) -> kwant.system.FiniteSystem:
    """
    Endliches Lx×Ly-Gitter aus dem diskretisierten Kontinuum (4 NORBs pro Platz).

    Entspricht faktisch lat = quadratisches Gitter mit norbs=4, aber das Gitter
    und die nächsten Nachbar-Sprünge kommen aus continuum.discretize(H), nicht
    aus einem freien „builder“ pro Site.
    """
    h_num = sp.N(_subs_named(h_sym, t=t, **phases))
    template = continuum.discretize(h_num, grid=grid)
    syst = kwant.Builder()
    syst.fill(template, make_square_shape(Lx, Ly), (0, 0))
    return syst.finalized()


def compute_eigenvalues(
    params: dict[str, float],
    *,
    Lx: int = 30,
    Ly: int = 30,
    grid: float = 1.0,
    h_sym: sp.Matrix | None = None,
) -> np.ndarray:
    """
    Sortierte Eigenwerte des endlichen Lx×Ly-EABC-Gitters (``continuum.discretize`` + ``fill``).

    ``params`` enthält mindestens ``t``, ``phi_A``, ``phi_B``, ``phi_C`` (wie bei
    ``finite_square_from_continuum``). Optional kann ``h_sym`` wiederverwendet werden,
    um SymPy-``build_hamiltonian()`` bei vielen Läufen zu sparen.
    """
    if h_sym is None:
        h_sym, _ = build_hamiltonian()
    t = float(params["t"])
    phases = {k: float(params[k]) for k in ("phi_A", "phi_B", "phi_C")}
    fin = finite_square_from_continuum(
        h_sym, Lx=Lx, Ly=Ly, t=t, grid=grid, **phases
    )
    return np.linalg.eigvalsh(fin.hamiltonian_submatrix())


def compute_eigenvalues_randomized(
    *,
    Lx: int = 22,
    Ly: int = 22,
    t: float = 1.0,
    a: float = 1.0,
    seed: int | None = None,
) -> np.ndarray:
    """
    Sortierte Eigenwerte für **eine** Realisierung **örtlich zufälliger** EABC-Phasen
    (``finite_square_eabc_random``): pro Platz unabhängig φ_A, φ_B, φ_C ∈ [0, 2π).

    - ``seed=None``: jeder Aufruf neue Störung (nicht reproduzierbar).
    - ``seed`` ganzzahlig: dieselbe Realisierung bei wiederholtem Aufruf mit gleichem ``seed``.

    Beispiel (entfaltete Abstände + Histogramm)::

        E = compute_eigenvalues_randomized(Lx=30, Ly=30, seed=42)
        s = level_spacings(E)  # optional: spline_s=... bei Spline-Warnungen
        plot_spacing_histogram(E, unfolded=True)  # dieselben Abstände bei gleichem spline_s
    """
    fin = finite_square_eabc_random(Lx=Lx, Ly=Ly, t=t, a=a, seed=seed)
    return np.linalg.eigvalsh(fin.hamiltonian_submatrix())


def compute_eigenvalues_delta_phi_a_field(
    *,
    Lx: int,
    Ly: int,
    delta_field: np.ndarray,
    q: int = 4,
    t: float = 1.0,
    a: float = 1.0,
    phi_B: float | None = None,
    phi_C: float | None = None,
) -> np.ndarray:
    """
    Wie ``compute_eigenvalues_randomized``, aber **nur φ_A ortsabhängig** aus
    ``spatial_phi_A_from_delta_field(delta_field, q=q)``; φ_B und φ_C **konstant**
    (Vorgaben ``DEFAULT_PHASES``, überschreibbar).

    Minimalbeispiel (unterscheidet sich vom rein zufälligen Ensemble)::

        L = 30
        delta = np.fromfunction(lambda i, j: (i + j) % 4, (L, L), dtype=float)
        E = compute_eigenvalues_delta_phi_a_field(Lx=L, Ly=L, delta_field=delta, q=4)
        plot_spacing_histogram(E, unfolded=True)
    """
    arr = np.asarray(delta_field, dtype=float)
    if arr.shape != (int(Lx), int(Ly)):
        raise ValueError(
            "compute_eigenvalues_delta_phi_a_field: delta_field.shape muss (Lx, Ly) sein, "
            f"erhalten {arr.shape}"
        )
    phi_b_val = float(DEFAULT_PHASES["phi_B"] if phi_B is None else phi_B)
    phi_c_val = float(DEFAULT_PHASES["phi_C"] if phi_C is None else phi_C)
    ph_a = spatial_phi_A_from_delta_field(arr, q=int(q))
    fin = finite_square_eabc_spatial(
        Lx=int(Lx),
        Ly=int(Ly),
        phi_A=ph_a,
        phi_B=spatial_phase_constant(phi_b_val),
        phi_C=spatial_phase_constant(phi_c_val),
        t=t,
        a=a,
    )
    return np.linalg.eigvalsh(fin.hamiltonian_submatrix())


def compute_eigenvalues_delta_three_fields(
    *,
    Lx: int,
    Ly: int,
    delta_A: np.ndarray,
    delta_B: np.ndarray,
    delta_C: np.ndarray,
    q_A: int = 4,
    q_B: int = 4,
    q_C: int = 4,
    t: float = 1.0,
    a: float = 1.0,
) -> np.ndarray:
    """
    Eigenwerte mit örtlichen Kantenphasen aus den drei Rasterfeldern Delta_A, Delta_B, Delta_C
    (wie ``spatial_phases_from_three_delta_fields`` dann ``finite_square_eabc_spatial``).

    Alle Felder haben Shape ``(Lx, Ly)``.
    """
    sh = (int(Lx), int(Ly))
    for name, fld in ("delta_A", delta_A), ("delta_B", delta_B), ("delta_C", delta_C):
        ff = np.asarray(fld)
        if ff.shape != sh:
            raise ValueError(
                f"compute_eigenvalues_delta_three_fields: {name}.shape muss {sh} sein, "
                f"erhalten {ff.shape}"
            )

    phi_A, phi_B, phi_C = spatial_phases_from_three_delta_fields(
        delta_A, delta_B, delta_C, q_A=q_A, q_B=q_B, q_C=q_C
    )
    fin = finite_square_eabc_spatial(
        Lx=sh[0],
        Ly=sh[1],
        phi_A=phi_A,
        phi_B=phi_B,
        phi_C=phi_C,
        t=t,
        a=a,
    )
    return np.linalg.eigvalsh(fin.hamiltonian_submatrix())


def hamiltonian_k_numeric(
    h_sym: sp.Matrix,
    *,
    t: float,
    kx: float,
    ky: float,
    phi_A: float,
    phi_B: float,
    phi_C: float,
) -> np.ndarray:
    """Wertet H an (k_x, k_y) mit gegebenen Parametern numerisch aus."""
    m = _subs_named(
        h_sym,
        t=t,
        phi_A=phi_A,
        phi_B=phi_B,
        phi_C=phi_C,
        k_x=kx,
        k_y=ky,
    )
    return np.array(m.evalf(), dtype=complex)


def unfold_spectrum(
    E: np.ndarray,
    *,
    spline_s: float | None = None,
) -> np.ndarray:
    """
    Glatte empirische Zählfunktion Ñ(E): Spline durch (E_i, i) bei sortiertem Spektrum.

    ``spline_s`` ist der Glättungsparameter ``s`` von ``UnivariateSpline`` (größer =
    glatter). Standard: ``0.1 * N`` wie in gängigen Skizzen; bei zu wenig Punkten
    kleiner wählen.
    """
    E = np.sort(np.asarray(E, dtype=float))
    N = E.size
    if N < 3:
        raise ValueError("unfold_spectrum: mindestens drei Eigenwerte nötig")
    # Streng monoton in E, sonst scheitert der Spline als Funktion von E
    for i in range(1, N):
        if E[i] <= E[i - 1]:
            E[i] = np.nextafter(E[i - 1], np.inf)

    x = np.arange(N, dtype=float)
    if spline_s is None:
        spline_s = N * 0.1
    spline = UnivariateSpline(E, x, s=float(spline_s))
    return np.asarray(spline(E), dtype=float)


def level_spacings(E: np.ndarray, *, spline_s: float | None = None) -> np.ndarray:
    """
    Nächst-Nachbar-Abstände in **entfalteten** Einheiten, normiert auf ⟨s⟩ = 1.

    Entspricht: ``xi = unfold_spectrum(E)``, ``s = diff(xi)``, ``s / mean(s)``.
    """
    xi = unfold_spectrum(E, spline_s=spline_s)
    s = np.diff(xi)
    m = np.mean(s)
    if not np.isfinite(m) or m <= 0:
        raise ValueError(f"level_spacings: ungültiger mittlerer Abstand {m!r}")
    return s / m


def level_spacings_hermitian(
    H: np.ndarray,
    *,
    spline_s: float | None = None,
) -> np.ndarray:
    """
    Kurzform: ``eigs = np.linalg.eigvalsh(H)``, dann ``level_spacings(eigs)``.

    ``H`` quadratisch, reell-symmetrisch oder hermitesch (``eigvalsh`` liefert reelle
    Eigenwerte).
    """
    H = np.asarray(H)
    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError("level_spacings_hermitian: H muss eine quadratische Matrix sein")
    eigs = np.linalg.eigvalsh(H)
    return level_spacings(np.asarray(eigs, dtype=float), spline_s=spline_s)


def pooled_level_spacings_spatial_random(
    n_runs: int = DEFAULT_ENSEMBLE_SPECTRA,
    *,
    Lx: int = 22,
    Ly: int = 22,
    t: float = 1.0,
    a: float = 1.0,
    spline_s: float | None = None,
    seed: int | None = None,
) -> np.ndarray:
    """
    Wiederholt ``n_runs``-mal: neue räumliche Zufallsphasen → ``level_spacings``;
    alle Abstände werden **aneinandergehängt** (gemeinsames Histogramm, ``estimate_beta``, …).

    Entspricht::

        all_spacings = []
        for _ in range(n_runs):
            E = compute_eigenvalues_randomized(...)
            s_run = level_spacings(E, spline_s=spline_s)
            all_spacings.extend(s_run)

    ``seed``: gesetzt → Folge der Einzel-Realisierungen ist reproduzierbar; ``None`` → nicht.

    Empfohlen: **50–100** Spektren (``ENSEMBLE_SPECTRA_RECOMMENDED``); Standard ``n_runs=75``.
    """
    n = int(n_runs)
    if n < 1:
        raise ValueError("n_runs muss mindestens 1 sein")
    all_spacings: list[float] = []
    if seed is None:
        for _ in range(n):
            E = compute_eigenvalues_randomized(
                Lx=Lx, Ly=Ly, t=t, a=a, seed=None
            )
            s_run = level_spacings(E, spline_s=spline_s)
            all_spacings.extend(s_run)
    else:
        rng = np.random.default_rng(seed)
        for _ in range(n):
            run_seed = int(rng.integers(0, 2**31))
            E = compute_eigenvalues_randomized(
                Lx=Lx, Ly=Ly, t=t, a=a, seed=run_seed
            )
            s_run = level_spacings(E, spline_s=spline_s)
            all_spacings.extend(s_run)
    return np.asarray(all_spacings, dtype=float)


def pooled_level_spacings_delta_phase_offset_scan(
    n_runs: int,
    *,
    Lx: int,
    Ly: int,
    off_step_x: int = 7,
    off_step_y: int = 11,
    base_offx: int = 0,
    base_offy: int = 0,
    t: float = 1.0,
    a: float = 1.0,
    spline_s: float | None = None,
) -> np.ndarray:
    """
    Ensemble aus ``delta_phase_fields_with_offset`` mit Rasterverschiebung pro Lauf ``k``:

    ``offx = base_offx + k * off_step_x``, ``offy = base_offy + k * off_step_y``

    Für jedes ``k`` ein Spektrum, ``level_spacings`` (entfaltet), dann **alle Abstände
    konkateniert** (wie ``pooled_level_spacings_spatial_random`` für Histogramme und
    ``estimate_beta``).

    Beispiel::

        all_s = pooled_level_spacings_delta_phase_offset_scan(
            20, Lx=30, Ly=30, off_step_x=7, off_step_y=11,
        )
        beta_pooled = estimate_beta(all_s, s_max=0.5, bins=30)
        print("beta_full_delta (pooled) =", beta_pooled)
    """
    n = int(n_runs)
    if n < 1:
        raise ValueError("n_runs muss mindestens 1 sein")
    dx = int(off_step_x)
    dy = int(off_step_y)

    all_spacings: list[float] = []
    for k in range(n):
        ox = int(base_offx) + k * dx
        oy = int(base_offy) + k * dy
        phi_a, phi_b, phi_c = delta_phase_fields_with_offset(Lx, Ly, ox, oy)
        fin = finite_square_eabc_spatial(
            Lx=int(Lx),
            Ly=int(Ly),
            phi_A=phi_a,
            phi_B=phi_b,
            phi_C=phi_c,
            t=float(t),
            a=float(a),
        )
        evec = np.linalg.eigvalsh(fin.hamiltonian_submatrix())
        s_run = level_spacings(evec, spline_s=spline_s)
        all_spacings.extend(s_run.ravel())

    return np.asarray(all_spacings, dtype=float)


def normalized_nn_spacings(eigenvalues: np.ndarray) -> np.ndarray:
    """
    Aufeinanderfolgende Eigenwertabstände, normiert auf mittleren Abstand ⟨s⟩ = 1.

    Häufige Nutzung: Histogramm der s/⟨s⟩ (ungefaltet); bei chaotischen Systemen
    oft Vergleich mit Wigner-Dyson (GOE/GUE) oder Poisson exp(-s).
    """
    E = np.sort(np.asarray(eigenvalues, dtype=float))
    if E.size < 2:
        raise ValueError("mindestens zwei Eigenwerte nötig")
    s = np.diff(E)
    m = np.mean(s)
    if not np.isfinite(m) or m <= 0:
        raise ValueError(f"ungültiger mittlerer Abstand: {m!r}")
    return s / m


def estimate_beta(
    s: np.ndarray,
    *,
    s_max: float = 0.5,
    bins: int = 30,
) -> float:
    """
    Heuristischer Exponent im **kleinen-s**-Bereich: Histogramm-Dichte ``p(s)`` (density=True),
    dann **log-log-Regression** ``log p ≈ β log s + const`` nur für ``s < s_max``.

    Interpretation (grob): Wenn ``p(s) ∝ s^α`` nahe 0, ist ``β ≈ α`` — z.B. GOE
    ``p ∝ s`` → ``β ≈ 1``, GUE ``p ∝ s²`` → ``β ≈ 2``. Ergebnis hängt von ``bins``,
    ``s_max`` und der Stichprobengröße ab.

    Gibt ``nan`` zurück, wenn zu wenige positive Histogrammwerte im Fitbereich liegen.
    """
    s = np.asarray(s, dtype=float).ravel()
    s_small = s[s < float(s_max)]
    if s_small.size < 2:
        return float("nan")

    hist, edges = np.histogram(s_small, bins=int(bins), density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])

    mask = (hist > 0) & (centers > 0)
    if np.count_nonzero(mask) < 2:
        return float("nan")

    x = np.log(centers[mask])
    y = np.log(hist[mask])
    beta_fit, _ = np.polyfit(x, y, 1)
    return float(beta_fit)


def estimate_betas(
    *,
    s_max: float = 0.5,
    bins: int = 30,
    **named_spacings: np.ndarray,
) -> dict[str, float]:
    """
    Wie ``estimate_beta``, aber für mehrere Abstands-Samples auf einmal.

    Beispiel::

        out = estimate_betas(E=s_E, A=s_A, B=s_B, C=s_C)
        beta_E, beta_A = out["E"], out["A"]

    Alle Arrays werden mit denselben ``s_max`` / ``bins`` ausgewertet.
    """
    kw = dict(s_max=s_max, bins=bins)
    return {
        name: estimate_beta(np.asarray(s, dtype=float).ravel(), **kw)
        for name, s in named_spacings.items()
    }


def beta_monte_carlo_uniform_phases(
    n_draws: int = DEFAULT_ENSEMBLE_SPECTRA,
    *,
    Lx: int = 22,
    Ly: int = 22,
    t: float = 1.0,
    seed: int | None = None,
    spline_s: float | None = None,
    s_max: float = 0.5,
    bins: int = 30,
) -> tuple[np.ndarray, dict[str, float | int]]:
    """
    Zieht ``n_draws``-mal unabhängig uniforme Phasen φ_A,φ_B,φ_C ∈ [0,2π),
    jeweils ``level_spacings`` + ``estimate_beta`` auf dem endlichen Gitter.

    Rückgabe: ``(betas, stats)`` mit ``stats`` = ``mean``, ``std`` (``nanmean``/``nanstd``),
    ``n_nan``, ``n_draws``.

    Üblich: **50–100** unabhängige Spektren (``ENSEMBLE_SPECTRA_RECOMMENDED``); Standard ``n_draws=75``.
    """
    rng = np.random.default_rng(seed)
    h_sym, _ = build_hamiltonian()
    betas: list[float] = []
    for _ in range(int(n_draws)):
        params = dict(
            t=float(t),
            phi_A=float(rng.uniform(0.0, 2.0 * np.pi)),
            phi_B=float(rng.uniform(0.0, 2.0 * np.pi)),
            phi_C=float(rng.uniform(0.0, 2.0 * np.pi)),
        )
        E = compute_eigenvalues(params, Lx=Lx, Ly=Ly, h_sym=h_sym)
        s = level_spacings(E, spline_s=spline_s)
        betas.append(estimate_beta(s, s_max=s_max, bins=bins))
    arr = np.asarray(betas, dtype=float)
    stats: dict[str, float | int] = {
        "mean": float(np.nanmean(arr)),
        "std": float(np.nanstd(arr)),
        "n_nan": int(np.isnan(arr).sum()),
        "n_draws": int(n_draws),
    }
    return arr, stats


def measure_beta(
    lam: float,
    *,
    phi0: tuple[float, float, float] | None = None,
    t: float = 1.0,
    Lx: int = 22,
    Ly: int = 22,
    grid: float = 1.0,
    h_sym: sp.Matrix | None = None,
    spline_s: float | None = None,
    s_max: float = 0.5,
    bins: int = 30,
) -> float:
    """
    Eigenwerte → ``level_spacings`` → ``estimate_beta`` für ein Gitter mit

    ``φ_A = lam·φ_A0``, ``φ_B = lam·φ_B0``, ``φ_C = lam·φ_C0``.

    Damit lässt sich λ ∈ [0, 1] als **Einschalten** der EABC-Phasen von 0 bis zu
    festen Zielwerten ``phi0`` interpretieren. Standard-``phi0``: ``DEFAULT_PHASES``.

    Für einen Sweep: ``phi0`` einmal setzen (oder zufällig ziehen) und ``lam`` variieren;
    ``h_sym`` zwischen Aufrufen wiederverwenden.
    """
    if phi0 is None:
        phi0 = (
            float(DEFAULT_PHASES["phi_A"]),
            float(DEFAULT_PHASES["phi_B"]),
            float(DEFAULT_PHASES["phi_C"]),
        )
    lam = float(lam)
    params = dict(
        t=float(t),
        phi_A=lam * float(phi0[0]),
        phi_B=lam * float(phi0[1]),
        phi_C=lam * float(phi0[2]),
    )
    if h_sym is None:
        h_sym, _ = build_hamiltonian()
    E = compute_eigenvalues(params, Lx=Lx, Ly=Ly, grid=grid, h_sym=h_sym)
    s = level_spacings(E, spline_s=spline_s)
    return estimate_beta(s, s_max=s_max, bins=bins)


def beta_vs_lambda(
    lambdas: np.ndarray,
    *,
    phi0: tuple[float, float, float] | None = None,
    t: float = 1.0,
    Lx: int = 22,
    Ly: int = 22,
    grid: float = 1.0,
    spline_s: float | None = None,
    s_max: float = 0.5,
    bins: int = 30,
) -> np.ndarray:
    """Wie ``measure_beta`` für jedes λ in ``lambdas``; baut ``h_sym`` nur einmal."""
    h_sym, _ = build_hamiltonian()
    lam_arr = np.asarray(lambdas, dtype=float).ravel()
    return np.asarray(
        [
            measure_beta(
                lam,
                phi0=phi0,
                t=t,
                Lx=Lx,
                Ly=Ly,
                grid=grid,
                h_sym=h_sym,
                spline_s=spline_s,
                s_max=s_max,
                bins=bins,
            )
            for lam in lam_arr
        ],
        dtype=float,
    )


def p_goe(s):
    """Wigner-Surmise (GOE)."""
    return (np.pi / 2) * s * np.exp(-np.pi * s**2 / 4)


def p_gue(s):
    """Wigner-Surmise (GUE)."""
    return (32 / np.pi**2) * s**2 * np.exp(-4 * s**2 / np.pi)


def p_poisson(s):
    """NN-Abstands-Dichte bei unkorrelierten Niveaus."""
    return np.exp(-s)


def load_riemann_imag_parts(path: str | Path | None = None) -> np.ndarray:
    """
    Lädt die Folge γₙ (Imaginärteile der Nullstellen auf Re(s)=½), wie in ``zeros6.npy``.

    Nutzt Speicher-Mapping für große Dateien. Die übliche ``zeros6.npy`` ist bereits
    aufsteigend sortiert.
    """
    if path is not None:
        p = Path(path)
        if not p.is_file():
            raise FileNotFoundError(p)
    else:
        if ZEROS6_NPY.is_file():
            p = ZEROS6_NPY
        elif ZEROS6_NPZ.is_file():
            p = ZEROS6_NPZ
        else:
            raise FileNotFoundError(
                f"Weder {ZEROS6_NPY} noch {ZEROS6_NPZ} gefunden."
            )
    raw = np.load(p, mmap_mode="r", allow_pickle=False)
    if hasattr(raw, "files"):
        key = "zeros" if "zeros" in raw.files else raw.files[0]
        raw = raw[key]
    return raw


def get_riemann_zeros(
    path: str | Path | None = None,
    *,
    copy: bool = False,
) -> np.ndarray:
    """
    Zugriff auf die in ``zeros6.npy`` gespeicherten γ-Werte (Standard: ``ZEROS6_NPY``).

    Ohne ``path``: einmaliges Laden als **Memmap** (speicherschonend), danach gecacht.
    ``copy=True``: volle ``ndarray``-Kopie im RAM.
    """
    global _RIEMANN_ZEROS_MMAP
    if path is None:
        if _RIEMANN_ZEROS_MMAP is None:
            _RIEMANN_ZEROS_MMAP = load_riemann_imag_parts(None)
        src = _RIEMANN_ZEROS_MMAP
    else:
        src = load_riemann_imag_parts(path)
    if copy:
        return np.asarray(src, dtype=float).copy()
    return src


def normalized_zeta_spacings(gamma: np.ndarray) -> np.ndarray:
    """Aufeinanderfolgende γ-Abstände, normiert mit ``/ mean`` (wie in deinem Snippet)."""
    g = np.asarray(gamma, dtype=float).ravel()
    if g.size < 2:
        raise ValueError("mindestens zwei γ-Werte")
    d = np.diff(g)
    m = np.mean(d)
    if not np.isfinite(m) or m <= 0:
        raise ValueError("ungültiger mittlerer Zeta-Abstand")
    return d / m


def plot_eabc_vs_riemann_spacings(
    eabc_eigenvalues: np.ndarray,
    *,
    bins: int = 50,
    density: bool = True,
    eabc_unfolded: bool = False,
    spline_s: float | None = None,
    zeta_gamma: np.ndarray | None = None,
    zeta_path: str | Path | None = None,
    match_spacing_count: bool = True,
    outfile: str | Path | None = None,
    show: bool | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Überlagerung: normierte EABC-Abstände vs. normierte Riemann-γ-Abstände.

    ``zeta_s = zeta_spacings / mean(zeta_spacings)`` wie in deinem Code; EABC-Seite
    analog ``normalized_nn_spacings`` oder bei ``eabc_unfolded=True`` ``level_spacings``.

    Wenn ``match_spacing_count=True``, werden so viele aufeinanderfolgende γₙ
    verwendet, wie EABC-Abstände erzeugen (gleiche Histogramm-Stichprobengröße).
    """
    if eabc_unfolded:
        s = level_spacings(eabc_eigenvalues, spline_s=spline_s)
    else:
        s = normalized_nn_spacings(eabc_eigenvalues)

    if zeta_gamma is None:
        gamma_full = get_riemann_zeros(zeta_path)
    else:
        gamma_full = np.asarray(zeta_gamma, dtype=float).ravel()

    if match_spacing_count:
        n = int(s.size) + 1
        if gamma_full.size < n:
            raise ValueError(
                f"zu wenige Zeta-Nullstellen: brauche {n}, habe {gamma_full.size}"
            )
        gamma_use = np.asarray(gamma_full[:n], dtype=float)
    else:
        gamma_use = np.asarray(gamma_full, dtype=float)

    zeta_spacings = np.diff(gamma_use)
    # zeta spacings (normiert)
    zeta_s = zeta_spacings / np.mean(zeta_spacings)

    if outfile is None:
        outfile = Path(__file__).resolve().parent / "EABC_Kontinuum_spacings_vs_riemann.png"
    outfile = Path(outfile)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(s, bins=bins, density=density, alpha=0.5, color="steelblue", label="EABC")
    ax.hist(zeta_s, bins=bins, density=density, alpha=0.5, color="orangered", label="Riemann")
    ax.set_xlabel(r"Abstand / $\langle s\rangle$")
    ax.set_ylabel("Dichte" if density else "Anzahl")
    ax.set_title("Levelabstände: EABC vs. Riemann-γ")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    print(f"Vergleichs-Histogramm EABC/Riemann gespeichert: {outfile}")

    if show is None:
        show = os.environ.get("MPLBACKEND", "").lower() in (
            "tkagg",
            "qt5agg",
            "qtagg",
            "macosx",
        )
    if show:
        plt.show()
    plt.close(fig)
    return s, zeta_s


def plot_unfolded_spacings_classical_curves(
    eigenvalues: np.ndarray,
    *,
    bins: int = 50,
    xmax: float = 3.0,
    n_curve: int = 200,
    hist_alpha: float = 0.5,
    spline_s: float | None = None,
    outfile: str | Path | None = None,
    show: bool | None = None,
) -> np.ndarray:
    """
    Entfaltete Abstände (``level_spacings``) vs. GOE/GUE/Poisson — wie übliches RMT-Plot.

    Entspricht dem Muster: ``hist`` + ``plot`` auf ``np.linspace(0, 3, 200)``.
    """
    s = level_spacings(eigenvalues, spline_s=spline_s)
    x = np.linspace(0.0, xmax, n_curve)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(s, bins=bins, density=True, alpha=hist_alpha, label="data")
    ax.plot(x, p_goe(x), label="GOE")
    ax.plot(x, p_gue(x), label="GUE")
    ax.plot(x, p_poisson(x), label="Poisson")
    ax.set_xlabel(r"$(\tilde N_{i+1}-\tilde N_i)/\langle s\rangle$")
    ax.set_ylabel("Dichte")
    ax.set_title("Entfaltet: Daten vs. GOE / GUE / Poisson")
    ax.legend(loc="upper right")
    fig.tight_layout()

    if outfile is None:
        outfile = Path(__file__).resolve().parent / "EABC_Kontinuum_unfolded_classical.png"
    outfile = Path(outfile)
    fig.savefig(outfile, dpi=150)
    print(f"Entfaltungsplot (klassische Kurven) gespeichert: {outfile}")

    if show is None:
        show = os.environ.get("MPLBACKEND", "").lower() in (
            "tkagg",
            "qt5agg",
            "qtagg",
            "macosx",
        )
    if show:
        plt.show()
    plt.close(fig)
    return s


def plot_spacing_histogram(
    eigenvalues: np.ndarray,
    *,
    bins: int = 50,
    density: bool = True,
    unfolded: bool = False,
    spline_s: float | None = None,
    outfile: str | Path | None = None,
    show: bool | None = None,
    reference_curves: bool | tuple[str, ...] = True,
) -> np.ndarray:
    """
    Histogramm der normierten Nächst-Nachbar-Abstände.

    ``unfolded=True``: Abstände der **entfalteten** Zählfunktion (Spline), siehe
    ``level_spacings``; sonst rohe Energieabstände ``normalized_nn_spacings``.
    Wenn du ``s`` separat mit ``level_spacings(E, spline_s=σ)`` bildest, dieselbe
    ``spline_s=σ`` hier übergeben, damit Histogramm und ``s`` übereinstimmen.

    ``reference_curves``: Wenn ``True``, zeichnet GOE/GUE/Poisson-Vergleichskurven
    (nur sinnvoll bei ``density=True``). Oder Tupel aus ``"goe"``, ``"gue"``, ``"poisson"``.

    ``outfile``: PNG-Pfad (Standard: ``EABC_Kontinuum_level_spacings.png`` bzw.
    ``..._unfolded.png``).
    ``show``: bei ``True`` ``plt.show()``; bei ``None`` nur anzeigen, wenn
    ``MPLBACKEND`` interaktiv ist (sonst nur speichern).
    """
    if unfolded:
        s = level_spacings(eigenvalues, spline_s=spline_s)
    else:
        s = normalized_nn_spacings(eigenvalues)

    if outfile is None:
        base = Path(__file__).resolve().parent / (
            "EABC_Kontinuum_level_spacings_unfolded.png"
            if unfolded
            else "EABC_Kontinuum_level_spacings.png"
        )
        outfile = base
    outfile = Path(outfile)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(
        s,
        bins=bins,
        density=density,
        color="steelblue",
        edgecolor="white",
        linewidth=0.5,
        label="Daten",
    )

    if reference_curves and density:
        if reference_curves is True:
            keys = ("goe", "gue", "poisson")
        else:
            keys = tuple(reference_curves)
        smax = float(max(3.0, np.max(s) * 1.05))
        xs = np.linspace(0.0, smax, 300)
        styles = {
            "goe": ("GOE (Wigner)", "C1", "-"),
            "gue": ("GUE (Wigner)", "C2", "--"),
            "poisson": ("Poisson", "C3", ":"),
        }
        funcs = {"goe": p_goe, "gue": p_gue, "poisson": p_poisson}
        for k in keys:
            k = k.lower()
            if k not in funcs:
                continue
            lab, col, ls = styles[k]
            ax.plot(xs, funcs[k](xs), color=col, ls=ls, lw=2.0, label=lab)
        ax.legend(loc="upper right", fontsize=8)

    if unfolded:
        ax.set_xlabel(r"$(\tilde N_{i+1}-\tilde N_i)/\langle s\rangle$")
        ax.set_title("Normierte Abstände (entfaltet, Spline)")
    else:
        ax.set_xlabel(r"$(E_{i+1}-E_i)/\langle s\rangle$")
        ax.set_title("Normierte Levelabstände (NN)")
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    print(f"Levelabstands-Histogramm gespeichert: {outfile}")

    if show is None:
        show = os.environ.get("MPLBACKEND", "").lower() in (
            "tkagg",
            "qt5agg",
            "qtagg",
            "macosx",
        )
    if show:
        plt.show()
    plt.close(fig)
    return s


def plot_beta_vs_lambda(
    lambdas: np.ndarray,
    betas: np.ndarray,
    *,
    outfile: str | Path | None = None,
    show: bool | None = None,
    title: str | None = None,
) -> None:
    """
    ``β(λ)`` aus ``measure_beta`` / ``beta_vs_lambda`` mit Referenzlinien GOE (1) und GUE (2).
    """
    if outfile is None:
        outfile = Path(__file__).resolve().parent / "EABC_Kontinuum_beta_vs_lambda.png"
    outfile = Path(outfile)

    lam = np.asarray(lambdas, dtype=float).ravel()
    b = np.asarray(betas, dtype=float).ravel()
    if lam.size != b.size:
        raise ValueError("lambdas und betas müssen gleiche Länge haben")

    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(lam, b, "o-", color="steelblue", lw=1.5, markersize=5, label=r"$\hat\beta$")
    ax.axhline(1.0, color="C1", ls="--", lw=1.5, label="GOE (klein-$s$-Exponent ~1)")
    ax.axhline(2.0, color="C2", ls="--", lw=1.5, label="GUE (klein-$s$-Exponent ~2)")
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$\hat\beta$ (log-log Histogramm)")
    ax.set_title(title or r"Kleiner-$s$-Exponent $\hat\beta$ über $\lambda$")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    print(f"β–λ-Plot gespeichert: {outfile}")

    if show is None:
        show = os.environ.get("MPLBACKEND", "").lower() in (
            "tkagg",
            "qt5agg",
            "qtagg",
            "macosx",
        )
    if show:
        plt.show()
    plt.close(fig)


def main() -> None:
    h_sym, _ = build_hamiltonian()

    t = DEFAULT_T
    ph = DEFAULT_PHASES

    # Spektrum am Γ-Punkt (kontinuierliches Modell, ohne Gitter)
    h_gamma = hamiltonian_k_numeric(
        h_sym, t=t, kx=0.0, ky=0.0, phi_A=ph["phi_A"], phi_B=ph["phi_B"], phi_C=ph["phi_C"]
    )
    evals_gamma = np.linalg.eigvalsh(h_gamma)
    print("Kontinuum: Eigenwerte H(k=0) [t={}, Phasen φ_A,φ_B,φ_C = {:.2f},{:.2f},{:.2f}]".format(
        t, ph["phi_A"], ph["phi_B"], ph["phi_C"]
    ))
    print(np.array2string(evals_gamma, precision=4, suppress_small=True))

    h_dirac_sym, _ = build_dirac_hamiltonian()
    h_d_k0 = np.array(
        _subs_named(h_dirac_sym, k_x=0.0, k_y=0.0, **ph).evalf(),
        dtype=complex,
    )
    evals_dirac_k0 = np.linalg.eigvalsh(h_d_k0)
    print(
        f"\nDirac 8×8 bei k=0 (I₂⊗H_eabc): erste 4 Eigenwerte "
        f"{np.array2string(evals_dirac_k0[:4], precision=4)}"
    )

    builder_dirac = dirac_discretized_builder(grid=1.0, **ph)
    syst_dirac = kwant.wraparound.wraparound(builder_dirac).finalized()
    print("\nDirac diskretisiert (wraparound, 8 NORBs): Eigenwertspanne k_y=0:")
    for kx in np.linspace(-0.5, 0.5, 5):
        mat_d = syst_dirac.hamiltonian_submatrix(params=dict(k_x=kx, k_y=0.0))
        ev_d = np.linalg.eigvalsh(mat_d)
        print(f"  k_x = {kx:+.3f}  E in [{ev_d.min():.4f}, {ev_d.max():.4f}]")

    # Diskretisierung + 2D-Periodizität (Wraparound); nur k_x, k_y bleiben symbolisch
    h_num = sp.N(_subs_named(h_sym, t=t, **ph))
    tb = continuum.discretize(h_num, grid=1.0)
    syst = kwant.wraparound.wraparound(tb).finalized()

    print("\nDiskretisiert (wraparound): Eigenwertspanne entlang k_x, k_y=0:")
    for kx in np.linspace(-0.5, 0.5, 5):
        mat = syst.hamiltonian_submatrix(params=dict(k_x=kx, k_y=0.0))
        ev = np.linalg.eigvalsh(mat)
        print(f"  k_x = {kx:+.3f}  E in [{ev.min():.4f}, {ev.max():.4f}]")

    # Endliches 30×30-Gitter (4 NORBs ≙ EABC+kin aus discretize); kein wraparound
    L = 30
    fin = finite_square_from_continuum(h_sym, Lx=L, Ly=L, t=t, **ph)
    h_fin = fin.hamiltonian_submatrix()
    evals_fin = np.linalg.eigvalsh(h_fin)
    print(f"\nEndlich {L}×{L} Sites, 4 NORBs: Matrix {h_fin.shape[0]}×{h_fin.shape[1]}")
    print(f"  E_min = {evals_fin.min():.4f},  E_max = {evals_fin.max():.4f}")
    print(f"  erste 6 Eigenwerte: {np.array2string(evals_fin[:6], precision=4)}")
    print(f"  letzte 6 Eigenwerte: {np.array2string(evals_fin[-6:], precision=4)}")

    plot_spacing_histogram(evals_fin, bins=50, density=True, show=False)
    plot_unfolded_spacings_classical_curves(evals_fin, show=False)

    # Zweites Histogramm: φ_A(x,y) = phase_map(Δ(i,j)); φ_B/φ_C konstant wie oben –
    # sonst wären die Hilfsfunktionen „spatial_phi_from_delta“ ohne diesen Pfad nicht sichtbar.
    delta_demo = np.fromfunction(
        lambda i, j: (i + j) % 4, (L, L), dtype=float
    )
    evals_dp = compute_eigenvalues_delta_phi_a_field(
        Lx=L, Ly=L, delta_field=delta_demo, q=4, t=t
    )
    print(
        "\nΔ-gemapptes φ_A: E_min = {:.4f}, E_max = {:.4f} "
        "(unterscheidet sich strukturell vom Block mit konstanten Phasen oben)".format(
            evals_dp.min(), evals_dp.max()
        )
    )
    plot_spacing_histogram(
        evals_dp,
        unfolded=True,
        bins=50,
        density=True,
        show=False,
        outfile=_PACKAGE_DIR / "EABC_Kontinuum_level_spacings_unfolded_delta_phi.png",
    )

    # Alle drei Kanten aus delta_symbol_at + PHASE_MAP (delta_phase_fields)
    phi_da, phi_db, phi_dc = delta_phase_fields(L, L)
    fin_delta_full = finite_square_eabc_spatial(
        Lx=L,
        Ly=L,
        phi_A=phi_da,
        phi_B=phi_db,
        phi_C=phi_dc,
        t=t,
    )
    evals_delta_full = np.linalg.eigvalsh(fin_delta_full.hamiltonian_submatrix())
    print(
        "\nΔ voll Kanäle (delta_phase_fields): E_min = {:.4f}, E_max = {:.4f}".format(
            evals_delta_full.min(), evals_delta_full.max()
        )
    )
    plot_spacing_histogram(
        evals_delta_full,
        unfolded=True,
        bins=50,
        density=True,
        show=False,
        outfile=_PACKAGE_DIR / "EABC_Kontinuum_level_spacings_unfolded_delta_phi_full.png",
    )

    # Entfaltete Abstände → heuristischer klein-s-Exponent (wie Dokumentation zu estimate_beta)
    _spl = None  # konsistent bei Bedarf spline_s=N*0.06 o.Ä. gemeinsam setzen
    s_const = level_spacings(evals_fin, spline_s=_spl)
    s_delta = level_spacings(evals_dp, spline_s=_spl)
    s_delta_full = level_spacings(evals_delta_full, spline_s=_spl)
    beta_const = estimate_beta(s_const, s_max=0.5, bins=50)
    beta_delta = estimate_beta(s_delta, s_max=0.5, bins=50)
    beta_delta_full = estimate_beta(s_delta_full, s_max=0.5, bins=50)

    def _fmt_beta(bx: float) -> str:
        return "nan (zu wenige Bin-/Fitpunkte — Heuristik greift nicht)" if not np.isfinite(
            bx
        ) else f"{float(bx):.4f}"

    print(
        "\nestimate_beta auf entfalteten Levelabständen (ein 30×30-Spektrum; nur Orientierung):\n"
        f"  β̂ (konstante φ)                ≈ {_fmt_beta(beta_const)}\n"
        f"  β̂ (Δ nur φ_A, q·phase_map)     ≈ {_fmt_beta(beta_delta)}\n"
        f"  β̂ (Δ phi_A,B,C delta_phase…)  ≈ {_fmt_beta(beta_delta_full)}"
    )

    try:
        plot_eabc_vs_riemann_spacings(
            evals_fin, bins=50, density=True, eabc_unfolded=False, show=False
        )
    except FileNotFoundError as exc:
        print(f"\nRiemann-Vergleich übersprungen ({exc})")

    # Exakt gleiche Matrix wie fill(), wenn man lat.shape + neighbors wie oben setzt
    Lchk = 3
    h_man = finite_square_manual_lat(Lx=Lchk, Ly=Lchk, t=t, a=1.0, **ph).hamiltonian_submatrix()
    h_fil = finite_square_from_continuum(h_sym, Lx=Lchk, Ly=Lchk, t=t, **ph).hamiltonian_submatrix()
    dmax = float(np.max(np.abs(h_man - h_fil)))
    print(f"\nAbgleich manuelles Gitter vs. fill() ({Lchk}×{Lchk}): max|H_man−H_fill| = {dmax:.3e}")

    n_beta_mc = int(os.environ.get("EABC_BETA_MC", "0"))
    if n_beta_mc > 0:
        Lmc = int(os.environ.get("EABC_BETA_L", "22"))
        seed_mc = os.environ.get("EABC_BETA_SEED")
        seed_i = int(seed_mc) if seed_mc is not None else None
        print(
            f"\nMonte Carlo β: {n_beta_mc} Ziehungen, Gitter {Lmc}×{Lmc} "
            f"(EABC_BETA_MC / EABC_BETA_L / EABC_BETA_SEED)"
        )
        _, st = beta_monte_carlo_uniform_phases(
            n_beta_mc, Lx=Lmc, Ly=Lmc, t=t, seed=seed_i
        )
        print(f"  mean beta = {st['mean']:.6f}")
        print(f"  std beta  = {st['std']:.6f}")
        if st["n_nan"]:
            print(f"  nan-Zähler: {st['n_nan']}")


if __name__ == "__main__":
    main()
