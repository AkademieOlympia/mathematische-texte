from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_kappa_test import check_dynamics, kappa_prefix, kappa_prefix_word, run_test


def test_dynamics_small():
    for n in range(1, 51, 2):
        assert check_dynamics(n, 4)


def test_kappa_prefix_length():
    assert len(kappa_prefix(27, 6)) == 6


def test_kappa_undefined_mod12():
    # n ≡ 3 (mod 12) hat keine EABC-Klasse am Start
    assert kappa_prefix(3, 1)[0] is None


def test_run_test_smoke():
    detail, summary = run_test(500, 4)
    assert summary.dynamics_holds
    assert summary.odd_starts == 250
    assert "summary" in detail
