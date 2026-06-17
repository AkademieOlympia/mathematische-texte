from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_morley_tm_numerik import (
    _pearson_r,
    run_m2_correlations,
    walter_form_wm_oriented,
)


def _find_corr(report, x: str, y: str, domain: str):
    for c in report.correlations:
        if c.x == x and c.y == y and c.domain == domain:
            return c
    raise AssertionError(f"missing {x},{y} in {domain}")


def test_pearson_r_basic():
    r, n = _pearson_r([1.0, 2.0, 3.0, 4.0], [2.0, 4.0, 6.0, 8.0])
    assert n == 4
    assert abs(r - 1.0) < 1e-12


def test_pearson_r_degenerate():
    r, n = _pearson_r([1.0, 1.0, 1.0], [2.0, 3.0, 4.0])
    assert n == 3
    assert math.isnan(r)


def test_m2_correlations_structure():
    report = run_m2_correlations(epsilons=[0.06, 0.12, 0.18])
    assert report.stage == "M2-correlations"
    assert len(report.epsilons) == 3
    domains = {c.domain for c in report.correlations}
    assert domains == {"all_pooled", "curved_only"}
    for domain in domains:
        rows = [c for c in report.correlations if c.domain == domain]
        assert len(rows) == 6
        for row in rows:
            assert math.isfinite(row.pearson_r)
            assert row.n >= 6


def test_m2_correlations_sign_pattern_curved():
    """G_M ∥ K_G, W_M anti-parallel zu K_G und G_M (chart-nah)."""
    report = run_m2_correlations(epsilons=[0.06, 0.10, 0.14, 0.18])
    curved = "curved_only"
    r_gk = _find_corr(report, "G_M", "K_G", curved).pearson_r
    r_wk = _find_corr(report, "W_M", "K_G", curved).pearson_r
    r_gw = _find_corr(report, "G_M", "W_M", curved).pearson_r
    assert r_gk > 0.75
    assert r_wk < -0.75
    assert r_gw < -0.95


def test_m2_correlations_fm_kg_weak_on_curved():
    """F_M hängt von |K_G| ab — schwache ρ(F_M, K_G) bei konstantem |K|=1."""
    report = run_m2_correlations(epsilons=[0.06, 0.10, 0.14, 0.18])
    r_fk = _find_corr(report, "F_M", "K_G", "curved_only").pearson_r
    assert abs(r_fk) < 0.5


def test_m2_correlations_json_cli(tmp_path):
    out = tmp_path / "m2_corr.json"
    subprocess.run(
        [
            sys.executable,
            "collatz_morley_tm_numerik.py",
            "m2-correlations",
            "--json",
            str(out),
        ],
        cwd=Path(__file__).resolve().parents[1],
        check=True,
    )
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["stage"] == "M2-correlations"
    assert len(data["correlations"]) == 12
    keys = {(c["x"], c["y"], c["domain"]) for c in data["correlations"]}
    assert ("G_M", "K_G", "curved_only") in keys
    assert ("W_M", "K_G", "all_pooled") in keys


def test_walter_form_wm_oriented_stub():
    from collatz_morley_tm_numerik import euclidean_patch_triangle

    try:
        walter_form_wm_oriented(euclidean_patch_triangle(0.12))
    except NotImplementedError as exc:
        assert "W_M^or" in str(exc) or "orientierte" in str(exc).lower()
    else:
        raise AssertionError("expected NotImplementedError")
