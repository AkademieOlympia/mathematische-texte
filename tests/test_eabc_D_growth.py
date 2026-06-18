"""Tests für EABC D_E-Wachstumsdiagnostik und Dirichlet-Stub."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_D_growth import (
    THEORY_EVOLUTION,
    THEORY_ZIRKULATION,
    c4_laplacian_spectrum,
    classify_growth,
    d_e_at_x,
    d_e_series,
    default_x_grid,
    dirichlet_decomposition_stub,
    growth_report,
    run,
)
from collatz_eabc_holonomie_fehlerterm import holonomy_counts


def test_d_e_matches_holonomy_counts():
    for x in (10_000, 50_000, 100_000):
        row = d_e_at_x(x)
        ref = holonomy_counts(x)
        assert row["D_E"] == ref["D_E"]
        assert row["N_plus"] == ref["N_plus"]
        assert row["N_minus"] == ref["N_minus"]
        assert row["C_E"] == ref["C_E"]


def test_default_grid_bounds():
    grid = default_x_grid(1_000_000)
    assert grid[0] >= 1_000
    assert grid[-1] == 1_000_000
    assert all(grid[i] < grid[i + 1] for i in range(len(grid) - 1))


def test_growth_classification_structure():
    series = d_e_series((10_000, 100_000, 1_000_000))
    gc = classify_growth(series)
    assert gc["best_fit_model"] in ("O(1)", "O(log X)", "O(sqrt X)", "power_law")
    assert gc["scenario_letter"] in ("A", "B", "C", "D", "?")
    assert "models" in gc
    assert "diagnostics_at_max_X" in gc


def test_growth_at_1e6_not_O1_or_sqrt():
    """Bei X=10^6: D_E wächst, liegt aber weit unter sqrt(X)."""
    row = d_e_at_x(1_000_000)
    assert row["D_E"] > 10
    assert row["D_E"] < 0.2 * (1_000_000**0.5)
    gc = classify_growth(d_e_series(default_x_grid(1_000_000)))
    assert gc["preferred_scenario"] != "O(1)"
    assert gc["preferred_scenario"] != "O(sqrt X)"


def test_c4_spectrum():
    spec = c4_laplacian_spectrum()
    assert spec["eigenvalues_symmetrized"] == [0.0, 2.0, 2.0, 4.0]
    assert spec["spectral_gap"] == 2.0


def test_dirichlet_stub():
    stub = dirichlet_decomposition_stub(100_000)
    assert stub["D_E"] == holonomy_counts(100_000)["D_E"]
    assert set(stub["dirichlet_coefficients"]) == {"chi_4", "chi_3", "chi_12"}
    assert "character_sums" in stub
    assert stub["status"].startswith("stub")


def test_growth_report_meta():
    report = growth_report(max_x=50_000, grid=(1_000, 10_000, 50_000))
    assert report["meta"]["theory_evolution"] == THEORY_EVOLUTION
    assert report["meta"]["theory_zirkulation"] == THEORY_ZIRKULATION
    assert len(report["D_E_series"]) == 3
    assert report["boxed_conclusion"]["evolution"].startswith("Bell")


def test_run_writes_json(tmp_path):
    out = tmp_path / "growth.json"
    report = run(max_x=10_000, output=out)
    assert out.exists()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["meta"]["max_x"] == 10_000
    assert report["output_path"] == str(out)
