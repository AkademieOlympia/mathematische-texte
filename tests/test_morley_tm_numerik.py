from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_morley_tm_numerik import (
    morley_form_fm,
    morley_form_fm_sphere,
    probe_sphere_scaling,
    run_numerik,
    sphere_patch_triangle,
)


def test_euclidean_morley_fm_near_zero():
    tri = (
        np.array([0.0, 0.0]),
        np.array([4.0, 0.0]),
        np.array([1.0, 3.5]),
    )
    assert morley_form_fm(tri) < 1e-10


def test_sphere_small_triangle_positive_fm():
    tri = sphere_patch_triangle(np.array([0.0, 0.0, 1.0]), side_angle=0.15)
    fm = morley_form_fm_sphere(tri)
    assert fm > 0.0
    assert math.isfinite(fm)


def test_sphere_scaling_correlates_with_area():
    samples = probe_sphere_scaling(side_angles=[0.06, 0.10, 0.14, 0.20])
    areas = np.array([s.kg_area for s in samples])
    fms = np.array([s.f_m for s in samples])
    corr = float(np.corrcoef(areas, fms)[0, 1])
    assert corr > 0.95
    assert fms[-1] > fms[0]


def test_run_numerik_smoke():
    report = run_numerik()
    assert report.euclidean_fm_max < 1e-8
    assert len(report.sphere_scaling) >= 4
    assert report.fit_slope is not None
    assert report.fit_r2 is not None
