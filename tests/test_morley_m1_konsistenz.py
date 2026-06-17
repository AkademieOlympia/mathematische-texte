from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_morley_tm_numerik import (
    MORLEY_VARIANTS,
    VARIANT_NAMES,
    euclidean_patch_triangle,
    morley_vertices_euclidean,
    morley_vertices_sphere,
    pairwise_variant_distances,
    run_m1_konsistenz,
    sphere_patch_triangle,
    triangle_morley_distance,
)


def test_euclidean_plane_variants_agree():
    """Auf R²: eine Morley-Realisierung — paarweise d_ij = 0."""
    tri = euclidean_patch_triangle(0.5, orientation=0.2)
    mor = morley_vertices_euclidean(tri)
    for i, a in enumerate(VARIANT_NAMES):
        for b in VARIANT_NAMES[i + 1 :]:
            d = triangle_morley_distance(mor, mor)
            assert d < 1e-12, f"Euklid {a} vs {b}: d={d}"


def test_sphere_variants_produce_finite_outputs():
    tri = sphere_patch_triangle(np.array([0.0, 0.0, 1.0]), side_angle=0.12)
    for name in VARIANT_NAMES:
        verts = morley_vertices_sphere(tri, name)
        for v in verts:
            assert abs(float(np.linalg.norm(v)) - 1.0) < 1e-10
            assert np.all(np.isfinite(v))


def test_pairwise_distances_positive_for_dissimilar_variants():
    tri = sphere_patch_triangle(np.array([0.0, 0.0, 1.0]), side_angle=0.20)
    dists = pairwise_variant_distances(tri)
    assert len(dists) == 6  # C(4,2)
    # chart und exp_euclidean sind identisch implementiert → d=0
    assert dists["local_chart__exp_euclidean"] < 1e-12
    # geodätisch vs. chart sollte bei nicht-trivialem Dreieck sichtbar sein
    assert dists["geodesic_angles__local_chart"] > 0.0


def test_m1_scaling_data_structure():
    report = run_m1_konsistenz(epsilons=[0.06, 0.10, 0.16, 0.22])
    assert report.stage == "M1"
    assert len(report.samples) == 4
    assert len(report.fits) == 6
    for fit in report.fits:
        if fit.pair == "local_chart__exp_euclidean":
            assert fit.order_hint == "identisch (d≡0)"
        else:
            assert math.isfinite(fit.slope_loglog)
            assert fit.slope_loglog >= 1.5


def test_m1_distances_shrink_with_epsilon():
    report = run_m1_konsistenz(epsilons=[0.08, 0.24])
    small = report.samples[0].distances["geodesic_angles__parallel_transport"]
    large = report.samples[1].distances["geodesic_angles__parallel_transport"]
    assert large > small


def test_m1_json_cli(tmp_path):
    out = tmp_path / "m1_test.json"
    subprocess.run(
        [sys.executable, "collatz_morley_tm_numerik.py", "m1", "--json", str(out)],
        cwd=Path(__file__).resolve().parents[1],
        check=True,
    )
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["stage"] == "M1"
    assert len(data["samples"]) >= 4
    assert len(data["fits"]) == 6
