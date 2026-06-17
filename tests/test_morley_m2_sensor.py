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
    morley_form_fm_hyperbolic,
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


def test_m2_curved_surfaces_positive_fm():
    report = run_m2_sensor(epsilons=[0.08, 0.14, 0.20])
    for surf in ("S2", "H2"):
        vals = [s.f_m for s in report.samples if s.surface == surf]
        assert all(v > 0.0 for v in vals)


def test_m2_exponent_fit_structure():
    report = run_m2_sensor(epsilons=[0.05, 0.08, 0.12, 0.18])
    assert report.exponent_fit is not None
    ef = report.exponent_fit
    assert math.isfinite(ef.alpha)
    assert math.isfinite(ef.beta)
    assert ef.n_samples >= 8


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
    # gleiche ε-Familie: F_M/A auf S² und H² sollte dieselbe Größenordnung haben
    assert 0.1 < med_h / med_s < 10.0


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
    assert data["exponent_fit"] is not None
