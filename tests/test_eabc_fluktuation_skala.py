from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_fluktuation_skala_test import (
    DEFAULT_GRID,
    fit_scaling_hypotheses,
    run_scaling_experiment,
    snapshots_at_grid,
)
from collatz_eabc_invarianzprogramm import snapshot_at_x


def test_snapshots_match_snapshot_at_x():
    grid = [100, 500, 1000]
    rows = snapshots_at_grid(grid)
    for r in rows:
        ref = snapshot_at_x(r.x)
        assert abs(ref.h - r.h) < 1e-9
        assert abs(ref.chi - r.chi) < 1e-12
        assert ref.pi == r.pi


def test_scaling_rows_monotone_pi():
    rows = snapshots_at_grid([100, 500, 1000, 5000])
    pis = [r.pi for r in rows]
    assert pis == sorted(pis)
    assert all(r.pi > 0 for r in rows)


def test_run_scaling_experiment_structure():
    result = run_scaling_experiment([100, 500, 1000])
    assert result["label"] == "Experiment"
    assert len(result["rows"]) == 3
    assert "scaling_analysis" in result
    assert "best_fit" in result["scaling_analysis"]
    row = result["rows"][-1]
    assert "H_over_pi" in row
    assert "chi_over_sqrt_pi" in row
    assert "mode_c" in row


def test_fit_scaling_hypotheses_returns_scores():
    rows = snapshots_at_grid([100, 500, 1000, 5000, 10000])
    analysis = fit_scaling_hypotheses(rows)
    assert len(analysis["hypothesis_scores"]) >= 3
    assert analysis["H_scaling_interpretation"]


def test_default_grid_endpoints():
    assert DEFAULT_GRID[0] == 100
    assert DEFAULT_GRID[-1] == 1_000_000
