"""Tests für Hurwitz–EABC-Orbit-Experiment (collatz_eabc_hurwitz_orbit_test.py)."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_hurwitz_orbit_test import (
    aggregate_orbit_report,
    chirality_score,
    enum_integer_solutions,
    entropy_from_counts,
    hurwitz_norm,
    hurwitz_orbit_elements,
    kappa_leg,
    prime_orbit_report,
    run,
    shannon_entropy,
)


def test_hurwitz_norm_unit():
    assert hurwitz_norm(1, 0, 0, 0) == 1
    assert hurwitz_norm(1, 1, 1, 1) == 4


def test_enum_integer_solutions_p5():
    sols = enum_integer_solutions(5)
    assert (1, 2, 0, 0) in sols or (-1, -2, 0, 0) in sols
    assert all(hurwitz_norm(*s) == 5 for s in sols)


def test_orbit_size_p2_matches_hurwitz_csv():
    elems = hurwitz_orbit_elements(2)
    assert len(elems) == 24


def test_orbit_size_p5():
    elems = hurwitz_orbit_elements(5)
    assert len(elems) == 144


def test_kappa_leg_zero():
    assert kappa_leg(0) == "0"


def test_kappa_leg_odd():
    assert kappa_leg(1) == "E"
    assert kappa_leg(5) == "A"


def test_chirality_balanced_four_e():
    assert chirality_score(("E", "E", "E", "E")) == 4
    assert chirality_score(("A", "B", "A", "B")) == -4
    assert chirality_score(("0", "E", "0", "C")) == 2


def test_prime_orbit_report_structure():
    report = prime_orbit_report(7)
    assert report["orbit_size"] == 192
    assert "mean_chi" in report
    assert "H_p" in report
    assert "mu_p_top" in report
    assert "channel_correlation" in report
    assert report["distinct_gamma"] >= 1
    assert report["H_p"] > 0


def test_shannon_entropy_uniform():
    mu = {"a": 0.5, "b": 0.5}
    assert abs(shannon_entropy(mu) - math.log(2)) < 1e-9


def test_entropy_from_counts():
    from collections import Counter

    counts = Counter({("E", "E", "E", "E"): 2, ("A", "B", "C", "0"): 2})
    H = entropy_from_counts(counts, 4)
    assert abs(H - math.log(2)) < 1e-9


def test_aggregate_has_H_stats():
    scale = aggregate_orbit_report(500)
    assert "H_stats" in scale
    assert scale["H_stats"]["mean"] > 0
    assert "H_by_p_mod_12" in scale
    assert "by_residue" in scale["H_by_p_mod_12"]


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "hurwitz_orbit.json"
    report = run(max_ps=[200], output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert "scales" in loaded
    assert (
        loaded["scales"]["200"]["prime_count"] == report["scales"]["200"]["prime_count"]
    )


def test_aggregate_has_orbit_stats(tmp_path: Path):
    out = tmp_path / "hurwitz_orbit.json"
    report = run(max_ps=[500], output=out)
    scale = report["scales"]["500"]
    assert scale["orbit_size_stats"]["min"] >= 24
    assert scale["total_orbit_points"] > 0


def test_half_integer_contributes_for_some_primes():
    report = prime_orbit_report(3)
    assert report["half_integer_coords"] > 0
