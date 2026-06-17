#!/usr/bin/env python3
"""
Geodätischer Morley-Operator T_M^(g) — Stufen M1 → M2 → M3.

Reihenfolge (Reviewer-Vorgabe, Tao-Stil)
----------------------------------------
M1  Operator-Konsistenz: vier Realisierungen von T_M^(g) auf schrumpfenden
    geodätischen Dreiecken Δ_ε; paarweise Abstände d_ij(ε); Skalierung O(ε) vs O(ε²).
M2  Modellräume: R² (K_G=0), S² (K_G>0), H² (K_G<0) — nach bestandener M1.
M3  Morley-Sensor: F_M(Δ), Tests F_M ∝ K_G A bzw. K_G² A² — erst nach M1+M2.

Epistemischer Rahmen
--------------------
Morley ist ein geometrisches Analogie-/Testmodul — kein Collatz-Beweisbaustein.
Möglicherweise eigenständiges Projekt „Morley-Operatoren auf riemannschen Flächen“.

Grenzen (ehrlich)
-----------------
- Vier Varianten sind numerische Approximationen, keine Levi-Civita-Implementierung
  in voller Allgemeinheit.
- Hyperbolische Fläche (H²): in M1/M2 nur als Platzhalter dokumentiert.
- Paralleltransport: diskret entlang Großkreisen auf S²; lokale Karten / exp: Log-Exp
  am Dreieckschwerpunkt (nicht isotherm).
- Kein Beweis; nur numerische Plausibilität / Definitionsvergleich.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Callable, Iterable, Sequence

import numpy as np

Array = np.ndarray
Vec2 = Array
Vec3 = Array

VARIANT_NAMES = (
    "geodesic_angles",
    "local_chart",
    "exp_euclidean",
    "parallel_transport",
)


class Variant(str, Enum):
    GEODESIC_ANGLES = "geodesic_angles"
    LOCAL_CHART = "local_chart"
    EXP_EUCLIDEAN = "exp_euclidean"
    PARALLEL_TRANSPORT = "parallel_transport"


# ---------------------------------------------------------------------------
# Grundwerkzeuge
# ---------------------------------------------------------------------------


def _unit(v: Array, eps: float = 1e-12) -> Array:
    n = float(np.linalg.norm(v))
    if n < eps:
        raise ValueError("degeneriertes Objekt")
    return v / n


def _signed_area2(a: Vec2, b: Vec2, c: Vec2) -> float:
    return float((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def _interior_angle(p: Array, q: Array, r: Array) -> float:
    v1 = _unit(p - q)
    v2 = _unit(r - q)
    return float(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)))


def _rotate_2d(v: Vec2, angle: float) -> Vec2:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1]], dtype=float)


def _rotate_toward_2d(v1: Vec2, v2: Vec2, angle: float) -> Vec2:
    v1 = _unit(v1)
    v2 = _unit(v2)
    cross = v1[0] * v2[1] - v1[1] * v2[0]
    sign = 1.0 if cross >= 0.0 else -1.0
    return _unit(_rotate_2d(v1, sign * angle))


def _rotate_toward_3d(v1: Array, v2: Array, angle: float) -> Array:
    """Dreht v1 im von v1×v2 bestimmten Sinne um angle Richtung v2."""
    v1 = _unit(v1)
    v2 = _unit(v2)
    axis = np.cross(v1, v2)
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm < 1e-12:
        return v1
    axis = axis / axis_norm
    return _unit(
        v1 * math.cos(angle)
        + np.cross(axis, v1) * math.sin(angle)
        + axis * np.dot(axis, v1) * (1.0 - math.cos(angle))
    )


def _rodrigues_rotate(v: Array, axis: Array, angle: float) -> Array:
    axis = _unit(axis)
    return (
        v * math.cos(angle)
        + np.cross(axis, v) * math.sin(angle)
        + axis * np.dot(axis, v) * (1.0 - math.cos(angle))
    )


def _line_intersection_2d(p1: Vec2, d1: Vec2, p2: Vec2, d2: Vec2) -> Vec2:
    cross = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(cross) < 1e-12:
        raise ValueError("parallele Trisektoren")
    diff = p2 - p1
    t = (diff[0] * d2[1] - diff[1] * d2[0]) / cross
    return p1 + t * d1


def _triangle_angles(vertices: Sequence[Array]) -> tuple[float, float, float]:
    a, b, c = vertices
    return (
        _interior_angle(c, a, b),
        _interior_angle(a, b, c),
        _interior_angle(b, c, a),
    )


def _triangle_area(vertices: Sequence[Array]) -> float:
    a, b, c = vertices
    return 0.5 * abs(_signed_area2(a, b, c))


# ---------------------------------------------------------------------------
# Euklidischer Morley (Referenz)
# ---------------------------------------------------------------------------


def morley_vertices_euclidean(vertices: Sequence[Vec2]) -> tuple[Vec2, Vec2, Vec2]:
    """Klassische Morley-Konstruktion in R^2 (Morley-Satz)."""
    a, b, c = [np.asarray(v, dtype=float) for v in vertices]
    if _signed_area2(a, b, c) < 0:
        a, c = c, a

    def trisector_from(q: Vec2, toward: Vec2, away: Vec2) -> Vec2:
        ang = _interior_angle(toward, q, away)
        base = _unit(toward - q)
        return _rotate_toward_2d(base, _unit(away - q), ang / 3.0)

    d_a = trisector_from(a, b, c)
    d_b = trisector_from(b, c, a)
    d_c = trisector_from(c, a, b)
    d_a2 = trisector_from(a, c, b)
    d_b2 = trisector_from(b, a, c)
    d_c2 = trisector_from(c, b, a)

    p_ab = _line_intersection_2d(a, d_a, b, d_b2)
    p_bc = _line_intersection_2d(b, d_b, c, d_c2)
    p_ca = _line_intersection_2d(c, d_c, a, d_a2)
    return p_ab, p_bc, p_ca


def morley_form_fm(vertices: Sequence[Array]) -> float:
    """F_M(Δ) = Σ (θ_i^M - π/3)^2 für Morley-Dreieck (euklidisch)."""
    mor = morley_vertices_euclidean(vertices)
    angles = _triangle_angles(mor)
    target = math.pi / 3.0
    return sum((ang - target) ** 2 for ang in angles)


# ---------------------------------------------------------------------------
# Sphärische Hilfsfunktionen
# ---------------------------------------------------------------------------


def _tangent_basis_at_sphere(q: Vec3) -> tuple[Vec3, Vec3]:
    q = _unit(q)
    ref = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(ref, q)) > 0.9:
        ref = np.array([1.0, 0.0, 0.0])
    e1 = _unit(np.cross(ref, q))
    e2 = _unit(np.cross(q, e1))
    return e1, e2


def _tangent_at_sphere(q: Vec3, p: Vec3) -> Vec3:
    q = _unit(q)
    v = p - q
    return _unit(v - np.dot(v, q) * q)


def _project_tangent(q: Vec3, v: Vec3) -> Vec3:
    q = _unit(q)
    return v - np.dot(v, q) * q


def _great_circle_plane_normal(q: Vec3, tangent_dir: Vec3) -> Vec3:
    return _unit(np.cross(q, tangent_dir))


def _great_circle_intersection(n1: Vec3, n2: Vec3) -> tuple[Vec3, Vec3]:
    v = np.cross(n1, n2)
    n = float(np.linalg.norm(v))
    if n < 1e-12:
        raise ValueError("parallele Großkreise")
    p = v / n
    return _unit(p), _unit(-p)


def _spherical_angle(p: Vec3, q: Vec3, r: Vec3) -> float:
    t1 = _tangent_at_sphere(q, p)
    t2 = _tangent_at_sphere(q, r)
    return float(np.arccos(np.clip(np.dot(t1, t2), -1.0, 1.0)))


def _spherical_triangle_area(angles: Iterable[float]) -> float:
    return float(sum(angles) - math.pi)


def _geodesic_distance_sphere(p: Vec3, q: Vec3) -> float:
    return float(np.arccos(np.clip(np.dot(_unit(p), _unit(q)), -1.0, 1.0)))


def _triangle_centroid_sphere(a: Vec3, b: Vec3, c: Vec3) -> Vec3:
    return _unit(a + b + c)


def _parallel_transport_sphere(v: Vec3, p: Vec3, q: Vec3) -> Vec3:
    """Paralleltransport von v ∈ T_p S² entlang des Großkreises p → q."""
    p, q = _unit(p), _unit(q)
    v = _project_tangent(p, v)
    axis = np.cross(p, q)
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm < 1e-12:
        return _unit(v)
    angle = math.acos(np.clip(np.dot(p, q), -1.0, 1.0))
    w = _rodrigues_rotate(v, axis / axis_norm, angle)
    return _unit(_project_tangent(q, w))


def _log_map_sphere(p: Vec3, x: Vec3, e1: Vec3, e2: Vec3) -> Vec2:
    """Log_p(x) in der Basis (e1,e2) am Punkt p."""
    p, x = _unit(p), _unit(x)
    t = _project_tangent(p, x)
    t_norm = float(np.linalg.norm(t))
    if t_norm < 1e-14:
        return np.zeros(2)
    dist = _geodesic_distance_sphere(p, x)
    t_hat = t / t_norm
    return np.array([np.dot(t_hat, e1), np.dot(t_hat, e2)]) * dist


def _exp_map_sphere(p: Vec3, v2: Vec2, e1: Vec3, e2: Vec3) -> Vec3:
    """exp_p(v) mit v in der Basis (e1,e2)."""
    p = _unit(p)
    v = v2[0] * e1 + v2[1] * e2
    v_norm = float(np.linalg.norm(v))
    if v_norm < 1e-14:
        return p
    return _unit(math.cos(v_norm) * p + math.sin(v_norm) * (v / v_norm))


def _choose_branch_sphere(p_pos: Vec3, p_neg: Vec3, reference: Vec3) -> Vec3:
    return p_pos if np.dot(p_pos, reference) >= np.dot(p_neg, reference) else p_neg


# ---------------------------------------------------------------------------
# Vier T_M^(g)-Varianten auf S²
# ---------------------------------------------------------------------------


def _trisector_direction_sphere(q: Vec3, p: Vec3, r: Vec3) -> Vec3:
    t_p = _tangent_at_sphere(q, p)
    t_r = _tangent_at_sphere(q, r)
    ang = _spherical_angle(p, q, r)
    return _rotate_toward_3d(t_p, t_r, ang / 3.0)


def morley_vertices_sphere_geodesic(vertices: Sequence[Vec3]) -> tuple[Vec3, Vec3, Vec3]:
    """
    Variante 1: geodätische Winkel — Trisektion im Tangentialraum, Kanten als Großkreise.
    """
    a, b, c = [_unit(np.asarray(v, dtype=float)) for v in vertices]
    centroid = _triangle_centroid_sphere(a, b, c)

    pairs = (
        (a, b, c, _trisector_direction_sphere(a, b, c), _trisector_direction_sphere(b, a, c)),
        (b, c, a, _trisector_direction_sphere(b, c, a), _trisector_direction_sphere(c, b, a)),
        (c, a, b, _trisector_direction_sphere(c, a, b), _trisector_direction_sphere(a, c, b)),
    )
    morley_pts: list[Vec3] = []
    for q1, q2, _, d1, d2 in pairs:
        n1 = _great_circle_plane_normal(q1, d1)
        n2 = _great_circle_plane_normal(q2, d2)
        p_pos, p_neg = _great_circle_intersection(n1, n2)
        morley_pts.append(_choose_branch_sphere(p_pos, p_neg, centroid))
    return morley_pts[0], morley_pts[1], morley_pts[2]


def morley_vertices_sphere_local_chart(vertices: Sequence[Vec3]) -> tuple[Vec3, Vec3, Vec3]:
    """
    Variante 2: lokale Karte am Schwerpunkt — Log-Exp (Approximation konformer Karte).
    """
    a, b, c = [_unit(np.asarray(v, dtype=float)) for v in vertices]
    p = _triangle_centroid_sphere(a, b, c)
    e1, e2 = _tangent_basis_at_sphere(p)
    tri2d = (
        _log_map_sphere(p, a, e1, e2),
        _log_map_sphere(p, b, e1, e2),
        _log_map_sphere(p, c, e1, e2),
    )
    m_ab, m_bc, m_ca = morley_vertices_euclidean(tri2d)
    return (
        _exp_map_sphere(p, m_ab, e1, e2),
        _exp_map_sphere(p, m_bc, e1, e2),
        _exp_map_sphere(p, m_ca, e1, e2),
    )


def morley_vertices_sphere_exp_euclidean(vertices: Sequence[Vec3]) -> tuple[Vec3, Vec3, Vec3]:
    """
    Variante 3: tangentialeuklidisches Modell — jedes Eck via exp_p aus R².

    Unterschied zu Variante 2: Morley-Schnittpunkte werden in der Tangentialebene
    konstruiert und einzeln zurücktransportiert (nicht identisch mit chart-Morley
    bei großem ε).
    """
    a, b, c = [_unit(np.asarray(v, dtype=float)) for v in vertices]
    p = _triangle_centroid_sphere(a, b, c)
    e1, e2 = _tangent_basis_at_sphere(p)
    tri2d = (
        _log_map_sphere(p, a, e1, e2),
        _log_map_sphere(p, b, e1, e2),
        _log_map_sphere(p, c, e1, e2),
    )
    return tuple(
        _exp_map_sphere(p, np.asarray(v, dtype=float), e1, e2) for v in morley_vertices_euclidean(tri2d)
    )


def morley_vertices_sphere_parallel_transport(vertices: Sequence[Vec3]) -> tuple[Vec3, Vec3, Vec3]:
    """
    Variante 4: Paralleltransport — Trisektionsrichtung der einen Ecke zum
    Großkreismittelpunkt der Kante transportiert, dann Schnitt (Approximation).
    """
    a, b, c = [_unit(np.asarray(v, dtype=float)) for v in vertices]
    centroid = _triangle_centroid_sphere(a, b, c)

    def morley_pair(q1: Vec3, q2: Vec3, r_other: Vec3, d1_raw: Vec3, d2_raw: Vec3) -> Vec3:
        mid = _unit(q1 + q2)
        d1 = _unit(_parallel_transport_sphere(d1_raw, q1, mid))
        d2 = _unit(_parallel_transport_sphere(d2_raw, q2, mid))
        n1 = _great_circle_plane_normal(q1, d1)
        n2 = _great_circle_plane_normal(q2, d2)
        p_pos, p_neg = _great_circle_intersection(n1, n2)
        ref = _unit(centroid + 0.1 * _unit(r_other))
        return _choose_branch_sphere(p_pos, p_neg, ref)

    return (
        morley_pair(a, b, c, _trisector_direction_sphere(a, b, c), _trisector_direction_sphere(b, a, c)),
        morley_pair(b, c, a, _trisector_direction_sphere(b, c, a), _trisector_direction_sphere(c, b, a)),
        morley_pair(c, a, b, _trisector_direction_sphere(c, a, b), _trisector_direction_sphere(a, c, b)),
    )


MORLEY_VARIANTS: dict[str, Callable[[Sequence[Vec3]], tuple[Vec3, Vec3, Vec3]]] = {
    Variant.GEODESIC_ANGLES.value: morley_vertices_sphere_geodesic,
    Variant.LOCAL_CHART.value: morley_vertices_sphere_local_chart,
    Variant.EXP_EUCLIDEAN.value: morley_vertices_sphere_exp_euclidean,
    Variant.PARALLEL_TRANSPORT.value: morley_vertices_sphere_parallel_transport,
}


def morley_vertices_sphere(vertices: Sequence[Vec3], variant: str) -> tuple[Vec3, Vec3, Vec3]:
    if variant not in MORLEY_VARIANTS:
        raise ValueError(f"unbekannte Variante: {variant}")
    return MORLEY_VARIANTS[variant](vertices)


def morley_form_fm_sphere(vertices: Sequence[Vec3], variant: str = "geodesic_angles") -> float:
    mor = morley_vertices_sphere(vertices, variant)
    angles = tuple(_spherical_angle(mor[(i + 2) % 3], mor[i], mor[(i + 1) % 3]) for i in range(3))
    target = math.pi / 3.0
    return sum((ang - target) ** 2 for ang in angles)


def _point_distance(a: Array, b: Array) -> float:
    """Euklidisch (R²/R³) oder geodätisch (S²-Einbettung), je nach Dimension."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.shape == (2,) and b.shape == (2,):
        return float(np.linalg.norm(a - b))
    return _geodesic_distance_sphere(a, b)


def triangle_morley_distance(
    mor_a: Sequence[Array],
    mor_b: Sequence[Array],
) -> float:
    """Maximale Eckpunkt-Abweichung zwischen zwei Morley-Dreiecken."""
    dists = [_point_distance(p, q) for p, q in zip(mor_a, mor_b, strict=True)]
    return float(max(dists))


def pairwise_variant_distances(
    vertices: Sequence[Vec3],
) -> dict[str, float]:
    """Paarweise d_ij für alle Variantenpaare."""
    outputs = {name: MORLEY_VARIANTS[name](vertices) for name in VARIANT_NAMES}
    distances: dict[str, float] = {}
    names = list(VARIANT_NAMES)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            key = f"{names[i]}__{names[j]}"
            distances[key] = triangle_morley_distance(outputs[names[i]], outputs[names[j]])
    return distances


# ---------------------------------------------------------------------------
# Dreiecksfamilien Δ_ε
# ---------------------------------------------------------------------------


def sphere_patch_triangle(
    center: Vec3,
    side_angle: float,
    orientation: float = 0.0,
) -> tuple[Vec3, Vec3, Vec3]:
    """Kleines sphärisches Dreieck um center (Einheitskugel, K_G=1)."""
    c = _unit(np.asarray(center, dtype=float))
    e1, e2 = _tangent_basis_at_sphere(c)
    rot = math.cos(orientation) * e1 + math.sin(orientation) * e2
    pts = []
    for phi in (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0):
        dir2d = math.cos(phi) * rot + math.sin(phi) * (
            math.cos(orientation) * e2 - math.sin(orientation) * e1
        )
        tangent = _unit(dir2d)
        p = math.cos(side_angle) * c + math.sin(side_angle) * tangent
        pts.append(_unit(p))
    return pts[0], pts[1], pts[2]


def euclidean_patch_triangle(
    side_length: float,
    orientation: float = 0.37,
) -> tuple[Vec2, Vec2, Vec2]:
    """Kleines euklidisches Dreieck (Kontrolle K_G=0), skaliert mit side_length."""
    base = np.array([1.0, 0.0])
    rot = _rotate_2d(base, orientation)
    pts = []
    for phi in (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0):
        direction = _rotate_2d(rot, phi)
        pts.append(side_length * direction)
    return pts[0], pts[1], pts[2]


def euclidean_variant_distances_plane(
    vertices: Sequence[Vec2],
) -> dict[str, float]:
    """Auf R²: alle Varianten degenerieren zum euklidischen Morley (Kontrolle)."""
    mor = morley_vertices_euclidean(vertices)
    zero = {f"{a}__{b}": 0.0 for i, a in enumerate(VARIANT_NAMES) for b in VARIANT_NAMES[i + 1 :]}
    # numerische Abweichung der euklidischen Realisierung von sich selbst
    for key in zero:
        zero[key] = triangle_morley_distance(mor, mor)
    return zero


# ---------------------------------------------------------------------------
# Regression / Skalungsfit
# ---------------------------------------------------------------------------


def linear_fit(x: Sequence[float], y: Sequence[float]) -> tuple[float, float]:
    xs = np.asarray(x, dtype=float)
    ys = np.asarray(y, dtype=float)
    if len(xs) < 2:
        return float("nan"), float("nan")
    A = np.vstack([xs, np.ones_like(xs)]).T
    slope, _intercept = np.linalg.lstsq(A, ys, rcond=None)[0]
    y_hat = slope * xs + _intercept
    ss_res = float(np.sum((ys - y_hat) ** 2))
    ss_tot = float(np.sum((ys - np.mean(ys)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(r2)


def loglog_fit(epsilons: Sequence[float], distances: Sequence[float]) -> tuple[float, float]:
    """log(d) ~ slope * log(ε) + c; slope ≈ 1 → O(ε), ≈ 2 → O(ε²)."""
    xs, ys = [], []
    for e, d in zip(epsilons, distances, strict=True):
        if e > 0 and d > 0:
            xs.append(math.log(e))
            ys.append(math.log(d))
    return linear_fit(xs, ys)


@dataclass
class ScalingFit:
    pair: str
    slope_loglog: float
    r2_loglog: float
    slope_linear_eps: float
    r2_linear_eps: float
    slope_linear_eps2: float
    r2_linear_eps2: float
    order_hint: str


@dataclass
class M1EpsilonSample:
    epsilon: float
    area: float
    distances: dict[str, float]


@dataclass
class M1Report:
    stage: str = "M1"
    surface: str = "S2"
    epsilons: list[float] = field(default_factory=list)
    samples: list[M1EpsilonSample] = field(default_factory=list)
    fits: list[ScalingFit] = field(default_factory=list)
    euclidean_control_max_distance: float = 0.0
    notes: list[str] = field(default_factory=list)


@dataclass
class ScalingSample:
    side_angle: float
    area: float
    kg_area: float
    f_m: float


@dataclass
class M3NumerikReport:
    stage: str = "M3"
    realization: str = "geodesic_angles_great_circles"
    euclidean_fm_max: float = 0.0
    sphere_scaling: list[ScalingSample] = field(default_factory=list)
    fit_slope: float | None = None
    fit_r2: float | None = None
    notes: list[str] = field(default_factory=list)


@dataclass
class M2SurfaceSample:
    surface: str
    kg: float
    epsilon: float
    area: float
    f_m: float
    f_m_over_a: float


@dataclass
class M2ExponentFit:
    log_c: float
    c: float
    alpha: float
    beta: float
    r2: float
    n_samples: int


@dataclass
class M2Report:
    stage: str = "M2"
    variant: str = "local_chart"
    epsilons: list[float] = field(default_factory=list)
    samples: list[M2SurfaceSample] = field(default_factory=list)
    euclidean_fm_max: float = 0.0
    naive_c: float | None = None
    naive_r2: float | None = None
    exponent_fit: M2ExponentFit | None = None
    m1_gate_passed: bool = False
    notes: list[str] = field(default_factory=list)


def _triangle_angles_sphere(vertices: Sequence[Vec3]) -> tuple[float, float, float]:
    a, b, c = vertices
    return (
        _spherical_angle(c, a, b),
        _spherical_angle(a, b, c),
        _spherical_angle(b, c, a),
    )


def _classify_order(slope: float) -> str:
    if not math.isfinite(slope):
        return "identisch (d≡0)"
    if slope >= 2.75:
        return "O(epsilon^3) — Konsistenz bis 2. Ordnung"
    if slope >= 1.75:
        return "O(epsilon^2) oder besser"
    if abs(slope - 1.0) < 0.35:
        return "O(epsilon) — Definitionskonflikt"
    return f"Zwischenordnung (slope≈{slope:.2f})"


def run_m1_konsistenz(
    center: Vec3 | None = None,
    epsilons: Sequence[float] | None = None,
    orientation: float = 0.37,
) -> M1Report:
    """
  M1: Operator-Konsistenz auf S² — paarweise Variantenabstände vs. ε.
    """
    center = np.array([0.0, 0.0, 1.0]) if center is None else _unit(np.asarray(center, float))
    if epsilons is None:
        epsilons = [0.04, 0.06, 0.08, 0.10, 0.14, 0.18, 0.24, 0.32]

    notes = [
        "M1 vor M2/M3: Definitionsvergleich, nicht F_M-Conjecture.",
        "Vier Varianten sind numerische Approximationen (s. Modul-Docstring).",
        "log-log-Steigung ≈ 3 spricht für O(ε³)-Abweichung → Konsistenz bis 2. Ordnung.",
        "H² (K_G<0) in M2 implementiert — M2a Modellraum-Vergleich.",
    ]

    samples: list[M1EpsilonSample] = []
    pair_series: dict[str, list[float]] = {f"{a}__{b}": [] for i, a in enumerate(VARIANT_NAMES) for b in VARIANT_NAMES[i + 1 :]}

    for eps in epsilons:
        tri = sphere_patch_triangle(center, float(eps), orientation=orientation)
        angles = _triangle_angles_sphere(tri)
        area = _spherical_triangle_area(angles)
        dists = pairwise_variant_distances(tri)
        samples.append(M1EpsilonSample(epsilon=float(eps), area=area, distances=dists))
        for key, val in dists.items():
            pair_series[key].append(val)

    eps_list = [s.epsilon for s in samples]
    fits: list[ScalingFit] = []
    for pair, dists in pair_series.items():
        slope_ll, r2_ll = loglog_fit(eps_list, dists)
        slope_e, r2_e = linear_fit(eps_list, dists)
        slope_e2, r2_e2 = linear_fit([e * e for e in eps_list], dists)
        fits.append(
            ScalingFit(
                pair=pair,
                slope_loglog=slope_ll,
                r2_loglog=r2_ll,
                slope_linear_eps=slope_e,
                r2_linear_eps=r2_e,
                slope_linear_eps2=slope_e2,
                r2_linear_eps2=r2_e2,
                order_hint=_classify_order(slope_ll),
            )
        )

    # Euklidische Kontrolle: Variantenabstand auf R² sollte ~0 sein
    euc_max = 0.0
    for eps in epsilons[:3]:
        tri_e = euclidean_patch_triangle(float(eps), orientation=orientation)
        mor = morley_vertices_euclidean(tri_e)
        euc_max = max(euc_max, triangle_morley_distance(mor, mor))

    return M1Report(
        epsilons=list(eps_list),
        samples=samples,
        fits=fits,
        euclidean_control_max_distance=euc_max,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Hyperbolische Hilfsfunktionen (Hyperboloid-Modell, K_G = -1)
# ---------------------------------------------------------------------------


def _minkowski_dot(u: Array, v: Array) -> float:
    return float(-u[0] * v[0] + u[1] * v[1] + u[2] * v[2])


def _hyperboloid_normalize(p: Array) -> Vec3:
    p = np.asarray(p, dtype=float)
    norm_sq = -_minkowski_dot(p, p)
    if norm_sq <= 0:
        raise ValueError("kein Hyperboloid-Punkt")
    return p / math.sqrt(norm_sq)


def _tangent_basis_at_hyperboloid(q: Vec3) -> tuple[Vec3, Vec3]:
    q = _hyperboloid_normalize(q)
    ref = np.array([0.0, 0.0, 1.0], dtype=float)
    if abs(-ref[0] * q[0] + ref[2] * q[2]) > 0.9:
        ref = np.array([0.0, 1.0, 0.0], dtype=float)
    e1 = _project_tangent_hyperboloid(q, ref)
    n1 = float(np.linalg.norm(e1))
    if n1 < 1e-12:
        ref = np.array([0.0, 1.0, 0.0], dtype=float)
        e1 = _project_tangent_hyperboloid(q, ref)
        n1 = float(np.linalg.norm(e1))
    e1 = e1 / n1
    raw_e2 = np.cross(np.array([q[0], -q[1], -q[2]], dtype=float), e1)
    e2 = _project_tangent_hyperboloid(q, raw_e2)
    e2 = e2 / float(np.linalg.norm(e2))
    return e1, e2


def _project_tangent_hyperboloid(q: Vec3, v: Vec3) -> Vec3:
    q = _hyperboloid_normalize(q)
    return v + q * _minkowski_dot(v, q)


def _tangent_at_hyperboloid(q: Vec3, p: Vec3) -> Vec3:
    q, p = _hyperboloid_normalize(q), _hyperboloid_normalize(p)
    e1, e2 = _tangent_basis_at_hyperboloid(q)
    v2 = _log_map_hyperboloid(q, p, e1, e2)
    v = v2[0] * e1 + v2[1] * e2
    n = float(np.linalg.norm(v))
    if n < 1e-14:
        raise ValueError("degenerierter Tangentialvektor")
    return v / n


def _geodesic_distance_hyperboloid(p: Vec3, q: Vec3) -> float:
    p, q = _hyperboloid_normalize(p), _hyperboloid_normalize(q)
    return float(math.acosh(max(1.0, -_minkowski_dot(p, q))))


def _hyperbolic_angle(p: Vec3, q: Vec3, r: Vec3) -> float:
    t1 = _tangent_at_hyperboloid(q, p)
    t2 = _tangent_at_hyperboloid(q, r)
    m1 = math.sqrt(max(1e-16, _minkowski_dot(t1, t1)))
    m2 = math.sqrt(max(1e-16, _minkowski_dot(t2, t2)))
    cos_th = _minkowski_dot(t1, t2) / (m1 * m2)
    return float(math.acos(np.clip(cos_th, -1.0, 1.0)))


def _hyperbolic_triangle_area(angles: Iterable[float]) -> float:
    return float(math.pi - sum(angles))


def _triangle_centroid_hyperboloid(a: Vec3, b: Vec3, c: Vec3) -> Vec3:
    raw = a + b + c
    return _hyperboloid_normalize(raw)


def _log_map_hyperboloid(p: Vec3, x: Vec3, e1: Vec3, e2: Vec3) -> Vec2:
    p, x = _hyperboloid_normalize(p), _hyperboloid_normalize(x)
    t = _project_tangent_hyperboloid(p, x)
    t_norm = float(np.linalg.norm(t))
    if t_norm < 1e-14:
        return np.zeros(2)
    dist = _geodesic_distance_hyperboloid(p, x)
    t_hat = t / t_norm
    return np.array([np.dot(t_hat, e1), np.dot(t_hat, e2)]) * dist


def _exp_map_hyperboloid(p: Vec3, v2: Vec2, e1: Vec3, e2: Vec3) -> Vec3:
    p = _hyperboloid_normalize(p)
    v = v2[0] * e1 + v2[1] * e2
    v_norm = float(np.linalg.norm(v))
    if v_norm < 1e-14:
        return p
    return _hyperboloid_normalize(math.cosh(v_norm) * p + math.sinh(v_norm) * (v / v_norm))


def hyperboloid_patch_triangle(
    center: Vec3 | None = None,
    side_angle: float = 0.1,
    orientation: float = 0.0,
) -> tuple[Vec3, Vec3, Vec3]:
    """Kleines hyperbolisches Dreieck (Hyperboloid, K_G = -1)."""
    c = _hyperboloid_normalize(
        np.array([1.0, 0.0, 0.0]) if center is None else np.asarray(center, dtype=float)
    )
    e1, e2 = _tangent_basis_at_hyperboloid(c)
    rot = math.cos(orientation) * e1 + math.sin(orientation) * e2
    pts = []
    for phi in (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0):
        dir3d = math.cos(phi) * rot + math.sin(phi) * (
            math.cos(orientation) * e2 - math.sin(orientation) * e1
        )
        tangent = dir3d / float(np.linalg.norm(dir3d))
        p = math.cosh(side_angle) * c + math.sinh(side_angle) * tangent
        pts.append(_hyperboloid_normalize(p))
    return pts[0], pts[1], pts[2]


def morley_vertices_hyperbolic_local_chart(vertices: Sequence[Vec3]) -> tuple[Vec3, Vec3, Vec3]:
    """Lokale Karte am Schwerpunkt — analog zu S²-Variante 2."""
    a, b, c = [_hyperboloid_normalize(np.asarray(v, dtype=float)) for v in vertices]
    p = _triangle_centroid_hyperboloid(a, b, c)
    e1, e2 = _tangent_basis_at_hyperboloid(p)
    tri2d = (
        _log_map_hyperboloid(p, a, e1, e2),
        _log_map_hyperboloid(p, b, e1, e2),
        _log_map_hyperboloid(p, c, e1, e2),
    )
    m_ab, m_bc, m_ca = morley_vertices_euclidean(tri2d)
    return (
        _exp_map_hyperboloid(p, m_ab, e1, e2),
        _exp_map_hyperboloid(p, m_bc, e1, e2),
        _exp_map_hyperboloid(p, m_ca, e1, e2),
    )


def morley_form_fm_hyperbolic(vertices: Sequence[Vec3]) -> float:
    mor = morley_vertices_hyperbolic_local_chart(vertices)
    angles = tuple(
        _hyperbolic_angle(mor[(i + 2) % 3], mor[i], mor[(i + 1) % 3]) for i in range(3)
    )
    target = math.pi / 3.0
    return sum((ang - target) ** 2 for ang in angles)


def _triangle_angles_hyperboloid(vertices: Sequence[Vec3]) -> tuple[float, float, float]:
    a, b, c = vertices
    return (
        _hyperbolic_angle(c, a, b),
        _hyperbolic_angle(a, b, c),
        _hyperbolic_angle(b, c, a),
    )


# ---------------------------------------------------------------------------
# M2 / M3 (hinter Flags — erst nach M1)
# ---------------------------------------------------------------------------


def euclidean_smoke_triangles() -> list[tuple[Vec2, Vec2, Vec2]]:
    return [
        (np.array([0.0, 0.0]), np.array([5.0, 0.0]), np.array([1.5, 4.0])),
        (np.array([0.0, 0.0]), np.array([3.0, 0.0]), np.array([0.8, 2.2])),
        (np.array([-1.0, 0.0]), np.array([2.0, 0.0]), np.array([0.5, 3.0])),
    ]


def probe_sphere_scaling(
    center: Vec3 | None = None,
    side_angles: Sequence[float] | None = None,
    variant: str = "geodesic_angles",
) -> list[ScalingSample]:
    """M3-Hilfe: F_M vs. K_G·Area auf S²."""
    center = np.array([0.0, 0.0, 1.0]) if center is None else _unit(np.asarray(center, float))
    if side_angles is None:
        side_angles = [0.05, 0.08, 0.12, 0.18, 0.25, 0.35]

    samples: list[ScalingSample] = []
    for s in side_angles:
        tri = sphere_patch_triangle(center, s, orientation=0.37)
        angles = _triangle_angles_sphere(tri)
        area = _spherical_triangle_area(angles)
        fm = morley_form_fm_sphere(tri, variant=variant)
        samples.append(ScalingSample(side_angle=float(s), area=area, kg_area=area, f_m=fm))
    return samples


def _m2_fm_for_surface(surface: str, tri: Sequence[Array], variant: str) -> float:
    if surface == "R2":
        return morley_form_fm(tri)
    if surface == "S2":
        return morley_form_fm_sphere(tri, variant=variant)
    if surface == "H2":
        return morley_form_fm_hyperbolic(tri)
    raise ValueError(f"unbekannte Fläche: {surface}")


def _m2_triangle_area(surface: str, tri: Sequence[Array]) -> float:
    if surface == "R2":
        return _triangle_area(tri)
    if surface == "S2":
        return _spherical_triangle_area(_triangle_angles_sphere(tri))
    if surface == "H2":
        return _hyperbolic_triangle_area(_triangle_angles_hyperboloid(tri))
    raise ValueError(f"unbekannte Fläche: {surface}")


def _multilinear_log_fit(
    log_abs_kg: Sequence[float],
    log_a: Sequence[float],
    log_fm: Sequence[float],
) -> tuple[float, float, float, float]:
    """log F_M ≈ log_c + α log|K_G| + β log A."""
    xs = np.asarray(log_abs_kg, dtype=float)
    ys_a = np.asarray(log_a, dtype=float)
    zs = np.asarray(log_fm, dtype=float)
    if len(xs) < 3:
        return float("nan"), float("nan"), float("nan"), float("nan")
    design = np.column_stack([xs, ys_a, np.ones_like(xs)])
    alpha, beta, log_c = np.linalg.lstsq(design, zs, rcond=None)[0]
    pred = design @ np.array([alpha, beta, log_c])
    ss_res = float(np.sum((zs - pred) ** 2))
    ss_tot = float(np.sum((zs - np.mean(zs)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(log_c), float(alpha), float(beta), float(r2)


def _m1_gate_passed(m1: M1Report, min_slope: float = 2.75) -> bool:
    for fit in m1.fits:
        if fit.pair == "local_chart__exp_euclidean":
            continue
        if not math.isfinite(fit.slope_loglog) or fit.slope_loglog < min_slope:
            return False
    return True


def run_m2_sensor(
    epsilons: Sequence[float] | None = None,
    orientation: float = 0.37,
    variant: str = "local_chart",
) -> M2Report:
    """
    M2a: Konstant-Krümmungsflächen (R², S², H²) — gleiche Δ_ε-Familie.
    M2b: Exponentenfit F_M = c |K_G|^α A^β (α, β aus Daten).
    """
    if epsilons is None:
        epsilons = [0.05, 0.08, 0.12, 0.18, 0.25, 0.35]

    m1 = run_m1_konsistenz(epsilons=list(epsilons)[:4])
    gate = _m1_gate_passed(m1)

    surfaces = (
        ("R2", 0.0),
        ("S2", 1.0),
        ("H2", -1.0),
    )
    samples: list[M2SurfaceSample] = []
    euc_max = 0.0

    for surf, kg in surfaces:
        for eps in epsilons:
            if surf == "R2":
                tri = euclidean_patch_triangle(float(eps), orientation=orientation)
            elif surf == "S2":
                tri = sphere_patch_triangle(np.array([0.0, 0.0, 1.0]), float(eps), orientation=orientation)
            else:
                tri = hyperboloid_patch_triangle(side_angle=float(eps), orientation=orientation)
            area = _m2_triangle_area(surf, tri)
            fm = _m2_fm_for_surface(surf, tri, variant=variant)
            ratio = fm / area if area > 1e-16 else float("nan")
            samples.append(
                M2SurfaceSample(
                    surface=surf,
                    kg=kg,
                    epsilon=float(eps),
                    area=area,
                    f_m=fm,
                    f_m_over_a=ratio,
                )
            )
            if surf == "R2":
                euc_max = max(euc_max, fm)

    curved = [s for s in samples if s.kg != 0.0 and s.f_m > 1e-18 and s.area > 1e-18]
    naive_c, naive_r2 = float("nan"), float("nan")
    if len(curved) >= 2:
        x_naive = [s.kg * s.area for s in curved]
        y_naive = [s.f_m for s in curved]
        naive_c, naive_r2 = linear_fit(x_naive, y_naive)

    log_abs_kg = [math.log(abs(s.kg)) for s in curved]
    log_a = [math.log(s.area) for s in curved]
    log_fm = [math.log(s.f_m) for s in curved]
    log_c, alpha, beta, r2_exp = _multilinear_log_fit(log_abs_kg, log_a, log_fm)
    exp_fit = M2ExponentFit(
        log_c=log_c,
        c=math.exp(log_c) if math.isfinite(log_c) else float("nan"),
        alpha=alpha,
        beta=beta,
        r2=r2_exp,
        n_samples=len(curved),
    )

    notes = [
        "M2a: gleiche ε-Familie auf R² (K_G=0), S² (+1), H² (−1).",
        "M2b: F_M = c |K_G|^α A^β — α, β nicht vorausgesetzt.",
        "Kanoniche Variante: local_chart (M1: O(ε³)-Konsistenz aller Varianten).",
        "K_G=0-Punkte nur Kontrolle (F_M≈0); Exponentenfit nur K_G≠0.",
        "M3 erst nach bestandenem M2-Gate (m1_gate_passed).",
    ]

    return M2Report(
        variant=variant,
        epsilons=list(epsilons),
        samples=samples,
        euclidean_fm_max=euc_max,
        naive_c=naive_c,
        naive_r2=naive_r2,
        exponent_fit=exp_fit,
        m1_gate_passed=gate,
        notes=notes,
    )


def run_m3_sensor(
    variant: str = "geodesic_angles",
) -> M3NumerikReport:
    """M3: Morley-Sensor F_M — nur nach M1+M2 sinnvoll."""
    notes = [
        "M3: F_M ~ c * K_G * A ist Conjecture, kein Theorem.",
        "K_G=1 auf der Einheitskugel; kg_area = K_G * Area.",
        "Vor M3 sollten M1-Varianten konsistent (O(ε²)) sein.",
    ]
    euc_max = max(morley_form_fm(t) for t in euclidean_smoke_triangles())
    scaling = probe_sphere_scaling(variant=variant)
    slope, r2 = linear_fit([s.kg_area for s in scaling], [s.f_m for s in scaling])
    return M3NumerikReport(
        realization=variant,
        euclidean_fm_max=euc_max,
        sphere_scaling=scaling,
        fit_slope=slope,
        fit_r2=r2,
        notes=notes,
    )


def run_numerik() -> M3NumerikReport:
    """Rückwärtskompatibler Alias für M3 (ältere Tests)."""
    return run_m3_sensor()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _print_m1_report(report: M1Report) -> None:
    print("=== Morley M1: Operator-Konsistenz (S²) ===")
    print(f"ε-Werte: {', '.join(f'{e:.3f}' for e in report.epsilons)}")
    print(f"Euklidische Kontrolle max d (soll ~0): {report.euclidean_control_max_distance:.3e}")
    print("\nStichproben (max paarweise d):")
    for s in report.samples:
        d_max = max(s.distances.values()) if s.distances else 0.0
        print(f"  ε={s.epsilon:.3f}  Area={s.area:.5f}  max d_ij={d_max:.5e}")
    print("\nSkalungsfits (log-log):")
    for fit in report.fits:
        print(
            f"  {fit.pair}: slope≈{fit.slope_loglog:.3f}  R²≈{fit.r2_loglog:.4f}  → {fit.order_hint}"
        )
    for note in report.notes:
        print(f"  • {note}")


def _print_m3_report(report: M3NumerikReport) -> None:
    print("=== Morley M3: Sensor F_M (S²) ===")
    print(f"Realisierung: {report.realization}")
    print(f"Euklid max F_M (soll ~0): {report.euclidean_fm_max:.3e}")
    for s in report.sphere_scaling:
        ratio = s.f_m / s.kg_area if s.kg_area > 0 else float("nan")
        print(
            f"  side={s.side_angle:.3f} rad  Area={s.area:.5f}  "
            f"F_M={s.f_m:.5e}  F_M/Area={ratio:.5e}"
        )
    if report.fit_slope is not None:
        print(f"Lineare Anpassung F_M ~ c * (K_G A): c ≈ {report.fit_slope:.5f}, R² ≈ {report.fit_r2:.4f}")


def _print_m2_report(report: M2Report) -> None:
    print("=== Morley M2: Modellräume + Exponentenfit ===")
    print(f"Variante: {report.variant}  |  M1-Gate: {'bestanden' if report.m1_gate_passed else 'offen'}")
    print(f"Euklid max F_M (K_G=0): {report.euclidean_fm_max:.3e}")
    print("\nM2a — F_M/A nach Fläche:")
    for surf in ("R2", "S2", "H2"):
        pts = [s for s in report.samples if s.surface == surf]
        if not pts:
            continue
        ratios = [s.f_m_over_a for s in pts if math.isfinite(s.f_m_over_a)]
        kg = pts[0].kg
        med = float(np.median(ratios)) if ratios else float("nan")
        print(f"  {surf} (K_G={kg:+.0f}): median F_M/A ≈ {med:.5e}  (n={len(pts)})")
    if report.naive_c is not None and math.isfinite(report.naive_c):
        print(f"\nNaiv F_M ~ c·K_G·A: c ≈ {report.naive_c:.5f}, R² ≈ {report.naive_r2:.4f}")
    if report.exponent_fit is not None:
        ef = report.exponent_fit
        print(
            f"\nM2b — F_M ≈ c |K_G|^α A^β:  c≈{ef.c:.5f}, α≈{ef.alpha:.3f}, β≈{ef.beta:.3f}, "
            f"R²≈{ef.r2:.4f}  (n={ef.n_samples})"
        )
    for note in report.notes:
        print(f"  • {note}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Morley T_M^(g): M1 Konsistenz (Default) → M2 Modelle → M3 Sensor"
    )
    parser.add_argument(
        "mode",
        nargs="?",
        default="m1",
        choices=("m1", "m2", "m3"),
        help="m1=Varianten-Konsistenz (Default), m2=Modellräume, m3=F_M-Sensor",
    )
    parser.add_argument("--json", type=str, default="", help="JSON-Ausgabepfad")
    parser.add_argument(
        "--variant",
        type=str,
        default="geodesic_angles",
        choices=VARIANT_NAMES,
        help="T_M^(g)-Variante für M3",
    )
    args = parser.parse_args()

    if args.mode == "m1":
        report = run_m1_konsistenz()
        _print_m1_report(report)
        payload = asdict(report)
    elif args.mode == "m2":
        report = run_m2_sensor()
        _print_m2_report(report)
        payload = asdict(report)
    else:
        report = run_m3_sensor(variant=args.variant)
        _print_m3_report(report)
        payload = asdict(report)

    out = args.json
    if not out:
        out = {
            "m1": "collatz_morley_m1_konsistenz.json",
            "m2": "collatz_morley_m2_sensor.json",
            "m3": "",
        }[args.mode]
    if out:
        with open(out, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)
        print(f"\nJSON geschrieben: {out}")


if __name__ == "__main__":
    main()
