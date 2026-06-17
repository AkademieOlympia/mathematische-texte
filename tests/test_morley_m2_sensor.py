from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_morley_tm_numerik import (
    hyperboloid_patch_triangle,
    morley_form_fm,
    morley_form_fm_hyperbolic,
    morley_form_gm,
    morley_form_gm_hyperbolic,
    run_m2_sensor,
)


def test_hyperbolic_patch_finite_area():
    tri = hyperboloid_patch_triangle(side_angle=0.12)
    fm = morley_form_fm_hyperbolic(tri)
    assert math.isfinite(fm)
    assert fm >= 0.0


def test_m2_three_surfaces_present():
    report = run_m2_sensor(epsilons=[0.06, 0.12, 0.18])
    surfaces = {s.surface for s in report.samples}
    assert surfaces == {"R2", "S2", "H2"}
    assert len(report.samples) == 9


def test_m2_euclidean_control_near_zero():
    report = run_m2_sensor(epsilons=[0.08, 0.16])
    assert report.euclidean_fm_max < 1e-20
    assert report.sign_test is not None
    assert report.sign_test.plane_near_zero is True


def test_m2_sign_structure_positive_on_curved():
    """M2a: F_M>0 auf S² und H², F_M≈0 auf R²."""
    report = run_m2_sensor(epsilons=[0.08, 0.14, 0.20])
    st = report.sign_test
    assert st is not None
    assert st.sphere_positive is True
    assert st.hyperbolic_positive is True
    assert st.plane_near_zero is True
    assert st.sphere_fm_median > 0.0
    assert st.hyperbolic_fm_median > 0.0


def test_m2_sphere_hyperbolic_near_equal_fm():
    """M2a: bei gleichem ε ist F_M(S²)≈F_M(H²) — kein Vorzeichen-Sensor."""
    report = run_m2_sensor(epsilons=[0.06, 0.10, 0.14, 0.18])
    st = report.sign_test
    assert st is not None
    assert len(st.ratio_sphere_over_hyperbolic) >= 3
    for ratio in st.ratio_sphere_over_hyperbolic:
        assert 0.99 < ratio < 1.01
    assert st.curvature_sign_detected is False
    assert 0.99 < st.median_ratio_sphere_over_hyperbolic < 1.01


def test_m2_exponent_fit_structure():
    report = run_m2_sensor(epsilons=[0.05, 0.08, 0.12, 0.18])
    assert report.exponent_fit is not None
    ef = report.exponent_fit
    assert math.isfinite(ef.alpha)
    assert math.isfinite(ef.beta)
    assert ef.n_samples >= 8
    # erste Daten: F_M ∝ A^2, schwache |K_G|-Kopplung
    assert ef.beta > 1.5


def test_m2_m1_gate_passed():
    report = run_m2_sensor(epsilons=[0.06, 0.10, 0.14])
    assert report.m1_gate_passed is True


def test_m2_sphere_hyperbolic_same_order_fm_over_a():
    report = run_m2_sensor(epsilons=[0.06, 0.10, 0.14, 0.18])
    s_ratios = [s.f_m_over_a for s in report.samples if s.surface == "S2"]
    h_ratios = [s.f_m_over_a for s in report.samples if s.surface == "H2"]
    assert len(s_ratios) == len(h_ratios) >= 3
    med_s = float(np.median(s_ratios))
    med_h = float(np.median(h_ratios))
    assert 0.1 < med_h / med_s < 10.0


def test_m2_geometry_table_present():
    report = run_m2_sensor(epsilons=[0.06, 0.12, 0.18])
    assert len(report.geometry_table) == 3
    for row in report.geometry_table:
        assert math.isfinite(row.f_m_median)
        assert math.isfinite(row.g_m_median)
        assert math.isfinite(row.f_m_over_a_median)
        assert math.isfinite(row.f_m_over_a2_median)
    plane = next(r for r in report.geometry_table if r.kg == 0.0)
    sphere = next(r for r in report.geometry_table if r.kg == 1.0)
    assert plane.f_m_median < 1e-20
    assert sphere.f_m_median > 0.0
    assert math.isfinite(plane.w_m_median)
    assert abs(plane.w_m_median) < 1e-10


def test_m2_euclidean_gm_near_zero():
    """M2b: G_M≈0 auf R² (Morley-Satz — gleichseitiger Kern)."""
    report = run_m2_sensor(epsilons=[0.08, 0.16])
    st = report.sign_test
    assert st is not None
    assert st.plane_gm_near_zero is True
    assert abs(st.plane_gm_median) < 1e-12


def test_m2_gm_definition_euclidean():
    """G_M = Σ(θ-π/3), F_M = Σ(θ-π/3)² auf R²."""
    from collatz_morley_tm_numerik import euclidean_patch_triangle

    tri = euclidean_patch_triangle(0.12)
    assert abs(morley_form_gm(tri)) < 1e-12
    assert morley_form_fm(tri) < 1e-20


def test_m2_gm_hyperbolic_finite():
    tri = hyperboloid_patch_triangle(side_angle=0.12)
    gm = morley_form_gm_hyperbolic(tri)
    assert math.isfinite(gm)


def test_m2_gm_opposite_signs_sphere_hyperbolic():
    """M2b: G_M(S²)>0, G_M(H²)<0 — Vorzeichen trennt Krümmung."""
    report = run_m2_sensor(epsilons=[0.06, 0.10, 0.14, 0.18])
    st = report.sign_test
    assert st is not None
    assert st.gm_sign_detected is True
    assert st.sphere_gm_median > 0.0
    assert st.hyperbolic_gm_median < 0.0
    for gs, gh in zip(st.sphere_gm_at_eps, st.hyperbolic_gm_at_eps, strict=True):
        assert gs > 0.0
        assert gh < 0.0


def test_m2_json_cli(tmp_path):
    out = tmp_path / "m2_test.json"
    subprocess.run(
        [sys.executable, "collatz_morley_tm_numerik.py", "m2", "--json", str(out)],
        cwd=Path(__file__).resolve().parents[1],
        check=True,
    )
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["stage"] == "M2"
    assert len(data["samples"]) >= 9
    assert data["sign_test"] is not None
    assert data["exponent_fit"] is not None
    st = data["sign_test"]
    assert st["plane_near_zero"] is True
    assert st["sphere_positive"] is True
    assert st["hyperbolic_positive"] is True
    assert "plane_gm_median" in st
    assert "sphere_gm_median" in st
    assert "hyperbolic_gm_median" in st
    assert "gm_sign_detected" in st
    assert "plane_wm_median" in st
    assert "wm_sign_detected" in st
    assert st["plane_wm_near_zero"] is True
