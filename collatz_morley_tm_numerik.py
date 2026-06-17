#!/usr/bin/env python3
"""
Erster numerischer Test des geodätischen Morley-Operators T_M^(g).

Ziel: die offene Vermutung F_M(Δ) ~ c * K_G(p) * Area(Δ) für kleine geodätische
Dreiecke empirisch sondieren (collatz_morley_metrik_erweiterung.md).

Epistemischer Rahmen
--------------------
Morley ist ein geometrisches Analogie-/Testmodul — kein Collatz-Beweisbaustein.

Grenzen (ehrlich)
-----------------
- Realisierung: Variante „geodätische Winkel“ (Winkeltrisektion im Tangentialraum,
  Fortsetzung entlang Großkreisen auf der Einheitskugel).
- Keine Levi-Civita-Paralleltransport- oder Exponentialkarten-Variante.
- Kein Beweis; nur numerische Plausibilität / Gegenbeispiel-Suche.
- Hyperbolische Flächen: nicht implementiert (nur S^2 und R^2).
- Konstante c hängt von der Normalisierung von F_M und der gewählten T_M^(g)-Variante ab.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from typing import Iterable, Sequence

import numpy as np

Array = np.ndarray
Vec2 = Array
Vec3 = Array


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
    """Dreht v1 im von v1×v2 bestimmten Sinne um angle Richtung v2 (2D oder 3D-Tangente)."""
    v1 = _unit(v1)
    v2 = _unit(v2)
    axis = np.cross(v1, v2)
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm < 1e-12:
        return v1
    axis = axis / axis_norm
    # Rodrigues
    return _unit(
        v1 * math.cos(angle)
        + np.cross(axis, v1) * math.sin(angle)
        + axis * np.dot(axis, v1) * (1.0 - math.cos(angle))
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


def morley_vertices_euclidean(vertices: Sequence[Vec2]) -> tuple[Vec2, Vec2, Vec2]:
    """Klassische Morley-Konstruktion in R^2 (Morley-Satz)."""
    a, b, c = [np.asarray(v, dtype=float) for v in vertices]
    if _signed_area2(a, b, c) < 0:
        a, c = c, a

    def trisector_from(q: Vec2, toward: Vec2, away: Vec2) -> Vec2:
        ang = _interior_angle(toward, q, away)
        base = _unit(toward - q)
        return _rotate_toward_2d(base, _unit(away - q), ang / 3.0)

    # Trisektoren „nah“ an der gemeinsamen Kante (Standard-Morley-Wahl).
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
    """F_M(Δ) = Σ (θ_i^M - π/3)^2 für Morley-Dreieck."""
    mor = morley_vertices_euclidean(vertices)
    angles = _triangle_angles(mor)
    target = math.pi / 3.0
    return sum((ang - target) ** 2 for ang in angles)


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


def morley_vertices_sphere_geodesic(vertices: Sequence[Vec3]) -> tuple[Vec3, Vec3, Vec3]:
    """
    T_M^(g) via geodätische Winkel: Trisektion im Tangentialraum, Kanten als Großkreise.
    """
    a, b, c = [_unit(np.asarray(v, dtype=float)) for v in vertices]

    def trisector_direction(q: Vec3, p: Vec3, r: Vec3) -> Vec3:
        t_p = _tangent_at_sphere(q, p)
        t_r = _tangent_at_sphere(q, r)
        ang = _spherical_angle(p, q, r)
        return _rotate_toward_3d(t_p, t_r, ang / 3.0)

    pairs = (
        (a, b, c, trisector_direction(a, b, c), trisector_direction(b, a, c)),
        (b, c, a, trisector_direction(b, c, a), trisector_direction(c, b, a)),
        (c, a, b, trisector_direction(c, a, b), trisector_direction(a, c, b)),
    )
    morley_pts: list[Vec3] = []
    for q1, q2, _, d1, d2 in pairs:
        n1 = _great_circle_plane_normal(q1, d1)
        n2 = _great_circle_plane_normal(q2, d2)
        p_pos, p_neg = _great_circle_intersection(n1, n2)
        # Wähle den Schnittpunkt näher am sphärischen Schwerpunkt.
        centroid = _unit(a + b + c)
        morley_pts.append(p_pos if np.dot(p_pos, centroid) >= np.dot(p_neg, centroid) else p_neg)
    return morley_pts[0], morley_pts[1], morley_pts[2]


def morley_form_fm_sphere(vertices: Sequence[Vec3]) -> float:
    mor = morley_vertices_sphere_geodesic(vertices)
    angles = tuple(_spherical_angle(mor[(i + 2) % 3], mor[i], mor[(i + 1) % 3]) for i in range(3))
    target = math.pi / 3.0
    return sum((ang - target) ** 2 for ang in angles)


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
    for k, phi in enumerate((0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0)):
        dir2d = math.cos(phi) * rot + math.sin(phi) * (
            math.cos(orientation) * e2 - math.sin(orientation) * e1
        )
        tangent = _unit(dir2d)
        # Exponentialabbildung: exp_c(side_angle * tangent)
        p = math.cos(side_angle) * c + math.sin(side_angle) * tangent
        pts.append(_unit(p))
    return pts[0], pts[1], pts[2]


@dataclass
class ScalingSample:
    side_angle: float
    area: float
    kg_area: float
    f_m: float


@dataclass
class NumerikReport:
    realization: str
    euclidean_fm_max: float
    sphere_scaling: list[ScalingSample]
    fit_slope: float | None
    fit_r2: float | None
    notes: list[str]


def euclidean_smoke_triangles() -> list[tuple[Vec2, Vec2, Vec2]]:
    return [
        (np.array([0.0, 0.0]), np.array([5.0, 0.0]), np.array([1.5, 4.0])),
        (np.array([0.0, 0.0]), np.array([3.0, 0.0]), np.array([0.8, 2.2])),
        (np.array([-1.0, 0.0]), np.array([2.0, 0.0]), np.array([0.5, 3.0])),
    ]


def probe_sphere_scaling(
    center: Vec3 | None = None,
    side_angles: Sequence[float] | None = None,
) -> list[ScalingSample]:
    center = np.array([0.0, 0.0, 1.0]) if center is None else _unit(np.asarray(center, float))
    if side_angles is None:
        side_angles = [0.05, 0.08, 0.12, 0.18, 0.25, 0.35]

    samples: list[ScalingSample] = []
    for s in side_angles:
        tri = sphere_patch_triangle(center, s, orientation=0.37)
        angles = _triangle_angles_sphere(tri)
        area = _spherical_triangle_area(angles)
        fm = morley_form_fm_sphere(tri)
        samples.append(ScalingSample(side_angle=float(s), area=area, kg_area=area, f_m=fm))
    return samples


def _triangle_angles_sphere(vertices: Sequence[Vec3]) -> tuple[float, float, float]:
    a, b, c = vertices
    return (
        _spherical_angle(c, a, b),
        _spherical_angle(a, b, c),
        _spherical_angle(b, c, a),
    )


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


def run_numerik() -> NumerikReport:
    notes = [
        "K_G=1 auf der Einheitskugel; kg_area = K_G * Area.",
        "Lineare Regression F_M ~ slope * (K_G * Area) nur heuristisch.",
    ]
    euc_max = max(morley_form_fm(t) for t in euclidean_smoke_triangles())
    scaling = probe_sphere_scaling()
    slope, r2 = linear_fit([s.kg_area for s in scaling], [s.f_m for s in scaling])
    return NumerikReport(
        realization="geodetic_angles_great_circles",
        euclidean_fm_max=euc_max,
        sphere_scaling=scaling,
        fit_slope=slope,
        fit_r2=r2,
        notes=notes,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Numerischer T_M^(g)-Test (Morley-Form F_M)")
    parser.add_argument("--json", type=str, default="", help="Optionaler JSON-Ausgabepfad")
    args = parser.parse_args()

    report = run_numerik()
    print("=== Morley T_M^(g) Numerik ===")
    print(f"Realisierung: {report.realization}")
    print(f"Euklid max F_M (soll ~0): {report.euclidean_fm_max:.3e}")
    print("Sphäre (kleine Dreiecke):")
    for s in report.sphere_scaling:
        ratio = s.f_m / s.kg_area if s.kg_area > 0 else float("nan")
        print(
            f"  side={s.side_angle:.3f} rad  Area={s.area:.5f}  "
            f"F_M={s.f_m:.5e}  F_M/Area={ratio:.5e}"
        )
    if report.fit_slope is not None:
        print(f"Lineare Anpassung F_M ~ c * (K_G A): c ≈ {report.fit_slope:.5f}, R² ≈ {report.fit_r2:.4f}")

    if args.json:
        payload = asdict(report)
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)
        print(f"JSON geschrieben: {args.json}")


if __name__ == "__main__":
    main()
