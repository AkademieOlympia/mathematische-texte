from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_morley_tm_numerik import (
    babylon_orthogonalize,
    morley_form_gm,
    run_m2_babylon,
    walter_form_wm,
)
from collatz_morley_tm_numerik import euclidean_patch_triangle, _triangle_area


def test_babylon_orthogonal_euclidean_plane_near_zero():
    """R²: G_M≈0, W_M≈0 → Orthogonalitätsdefekt ≈0."""
    tri = euclidean_patch_triangle(0.12)
    gm = morley_form_gm(tri)
    wm = walter_form_wm(tri, triangle_area=_triangle_area(tri))
    res = babylon_orthogonalize(gm, wm)
    assert res.orthogonality_defect < 1e-20
    assert res.hypotenuse_345 < 1e-10
    assert abs(gm) < 1e-12
    assert abs(wm) < 1e-10


def test_babylon_norm_ratio_3_4_on_synthetic():
    """Synthetisch: Projektionen α=β=1 → ‖leg3‖/‖leg4‖ = 3/4, Hypotenuse = 5."""
    res = babylon_orthogonalize(1.0, 1.0)
    assert abs(res.u3_norm / res.u4_norm - 3.0 / 4.0) < 1e-12
    assert abs(res.hypotenuse_345 - 5.0) < 1e-12
    assert res.orthogonality_defect < 1e-12


def test_babylon_pythagoras_defect_zero():
    """3²+4²=5²: Orthogonalitätsdefekt verschwindet für Babylon-Tripel."""
    res = babylon_orthogonalize(1.0, 1.0)
    assert res.orthogonality_defect < 1e-12
    assert abs(res.hypotenuse_345 - 5.0) < 1e-12


def test_babylon_batch_matches_scalar():
    gms = [0.001, -0.002, 0.003]
    wms = [-0.0008, 0.0016, -0.0024]
    batch = babylon_orthogonalize(gms, wms)
    assert isinstance(batch, list)
    assert len(batch) == 3
    for g, w, r in zip(gms, wms, batch, strict=True):
        single = babylon_orthogonalize(g, w)
        assert abs(single.hypotenuse_345 - r.hypotenuse_345) < 1e-15


def test_m2_babylon_report_structure():
    report = run_m2_babylon(epsilons=[0.06, 0.12, 0.18])
    assert report.stage == "M2-babylon"
    assert report.n_samples == 9
    assert report.plane_orthogonality_defect_median < 1e-20
    assert math.isfinite(report.curved_hypotenuse_median)
    assert report.pca_comparison is not None
    # Anti-korrelierte Daten: PCA fast 1D, Babylon erzwingt 2D-Achsen
    assert report.pca_comparison.explained_variance_ratio[0] > 0.99
    assert report.pca_comparison.angle_pc1_u3_deg > 25.0


def test_m2_babylon_json_cli(tmp_path):
    out = tmp_path / "babylon.json"
    subprocess.run(
        [
            sys.executable,
            "collatz_morley_tm_numerik.py",
            "m2-babylon",
            "--json",
            str(out),
        ],
        cwd=Path(__file__).resolve().parents[1],
        check=True,
    )
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["stage"] == "M2-babylon"
    assert len(data["per_sample"]) >= 6
    assert "gram_schmidt_u3" in data
    assert data["pca_comparison"]["babylon_u3"] == [1.0, 0.0]


def test_babylon_fm_ratio_finite_on_curved():
    report = run_m2_babylon(epsilons=[0.08, 0.14])
    curved = [
        r for r in report.per_sample if abs(r.gm) > 1e-12 or abs(r.wm) > 1e-12
    ]
    ratios = [r.fm_babylon_ratio for r in curved if r.fm_babylon_ratio is not None]
    assert ratios
    assert all(math.isfinite(x) and x > 0 for x in ratios)
