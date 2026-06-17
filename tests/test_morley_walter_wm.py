from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_morley_tm_numerik import (
    WALTER_EUCLIDEAN_AREA_RATIO,
    euclidean_patch_triangle,
    euclidean_smoke_triangles,
    hyperboloid_patch_triangle,
    run_m2_sensor,
    sphere_patch_triangle,
    walter_form_wm,
    walter_hexagon_area,
    walter_hexagon_area_hyperbolic,
    walter_hexagon_vertices_euclidean,
)
import numpy as np


def test_walter_euclidean_area_ratio_smoke():
    """Marion-Walter-Satz: Area(H_W)/Area(Δ) = 1/10 auf R²."""
    for tri in euclidean_smoke_triangles():
        area = 0.5 * abs(
            (tri[1][0] - tri[0][0]) * (tri[2][1] - tri[0][1])
            - (tri[1][1] - tri[0][1]) * (tri[2][0] - tri[0][0])
        )
        hex_area = walter_hexagon_area(tri)
        assert abs(hex_area / area - WALTER_EUCLIDEAN_AREA_RATIO) < 1e-12
        assert abs(walter_form_wm(tri, triangle_area=area)) < 1e-12


def test_walter_euclidean_patch_near_zero():
    for eps in (0.05, 0.12, 0.25, 0.35):
        tri = euclidean_patch_triangle(eps)
        assert abs(walter_form_wm(tri)) < 1e-10


def test_walter_hexagon_six_vertices():
    tri = euclidean_patch_triangle(0.12)
    hex_v = walter_hexagon_vertices_euclidean(tri)
    assert len(hex_v) == 6
    hex_area = walter_hexagon_area(tri)
    assert hex_area > 0.0


def test_walter_sphere_hyperbolic_finite():
    tri_s = sphere_patch_triangle(np.array([0.0, 0.0, 1.0]), 0.12)
    tri_h = hyperboloid_patch_triangle(side_angle=0.12)
    from collatz_morley_tm_numerik import _m2_wm_for_surface

    wm_s = _m2_wm_for_surface("S2", tri_s)
    wm_h = _m2_wm_for_surface("H2", tri_h)
    assert math.isfinite(wm_s)
    assert math.isfinite(wm_h)
    assert abs(wm_s) > 1e-6
    assert abs(wm_h) > 1e-6
    assert wm_s * wm_h < 0.0
    assert walter_hexagon_area(tri_s) > 0.0
    assert walter_hexagon_area_hyperbolic(tri_h) > 0.0


def test_m2_wm_plane_control():
    report = run_m2_sensor(epsilons=[0.08, 0.16])
    st = report.sign_test
    assert st is not None
    assert st.plane_wm_near_zero is True
    plane = [s for s in report.samples if s.surface == "R2"]
    assert all(abs(s.w_m) < 1e-10 for s in plane)


def test_m2_wm_opposite_sign_curved():
    """W_M trennt S²/H², aber anti-parallel zu G_M (chart-Näherung)."""
    report = run_m2_sensor(epsilons=[0.06, 0.10, 0.14, 0.18])
    st = report.sign_test
    assert st is not None
    assert st.wm_sign_detected is True
    assert st.sphere_wm_median < 0.0
    assert st.hyperbolic_wm_median > 0.0
    assert st.gm_wm_sign_agreement is False
