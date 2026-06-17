from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_morley_tm_numerik import run_m3_beweisversuch


def test_m3_gm_sign_sphere_hyperbolic():
    """G_M trennt Vorzeichen: S²>0, H²<0."""
    report = run_m3_beweisversuch(
        epsilons=[0.06, 0.10, 0.14, 0.18],
        check_variants=False,
    )
    sphere = [s for s in report.samples if s.surface == "S2"]
    hyper = [s for s in report.samples if s.surface == "H2"]
    assert all(s.g_m > 0.0 for s in sphere)
    assert all(s.g_m < 0.0 for s in hyper)


def test_m3_fm_near_zero_plane():
    """F_M≈0 auf der Ebene (Morley-Satz)."""
    report = run_m3_beweisversuch(epsilons=[0.08, 0.16], check_variants=False)
    pc = report.plane_control
    assert pc is not None
    assert pc.f_m_near_zero is True
    assert pc.f_m_max < 1e-20
    plane = [s for s in report.samples if s.surface == "R2"]
    assert all(s.f_m < 1e-20 for s in plane)


def test_m3_dual_fit_finite_exponents():
    """Dualer Fit liefert endliche α, β mit hohem R²."""
    report = run_m3_beweisversuch(
        epsilons=[0.05, 0.08, 0.12, 0.18],
        check_variants=False,
    )
    assert report.f_m_fit is not None
    assert report.g_m_fit is not None
    ff = report.f_m_fit
    gf = report.g_m_fit
    assert math.isfinite(ff.alpha)
    assert math.isfinite(ff.beta)
    assert math.isfinite(gf.alpha)
    assert math.isfinite(gf.beta)
    assert ff.n_samples >= 8
    assert gf.n_samples >= 8
    assert ff.r2 > 0.99
    assert gf.r2 > 0.99
    assert ff.beta > 1.5
    assert gf.beta > 0.5


def test_m3_variant_sign_cross_check():
    """sign(G_M)=sign(K_G) über alle M1-Varianten."""
    report = run_m3_beweisversuch(
        epsilons=[0.06, 0.12, 0.18],
        check_variants=True,
    )
    sc = report.sign_cross_check
    assert sc is not None
    assert sc.all_variants_sign_ok is True
    assert len(sc.variant_rows) == 4


def test_m3_m2_gate_passed():
    report = run_m3_beweisversuch(epsilons=[0.06, 0.10, 0.14], check_variants=False)
    assert report.m2_gate_passed is True


def test_m3_json_cli(tmp_path):
    out = tmp_path / "m3_test.json"
    subprocess.run(
        [
            sys.executable,
            "collatz_morley_tm_numerik.py",
            "m3",
            "--no-variant-check",
            "--json",
            str(out),
        ],
        cwd=Path(__file__).resolve().parents[1],
        check=True,
    )
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["stage"] == "M3"
    assert data["f_m_fit"] is not None
    assert data["g_m_fit"] is not None
    assert data["plane_control"]["f_m_near_zero"] is True
    assert len(data["samples"]) >= 9
