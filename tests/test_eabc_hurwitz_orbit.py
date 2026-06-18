"""Tests für Hurwitz–EABC-Schalen- und Orbit-Experiment (collatz_eabc_hurwitz_orbit_test.py)."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_hurwitz_orbit_test import (
    aggregate_shell_report,
    chirality_score,
    enum_integer_solutions,
    entropy_from_counts,
    full_gamma4,
    hurwitz_norm,
    hurwitz_orbit_elements,
    hurwitz_shell_elements,
    kappa_leg,
    leg_signature,
    prime_orbit_report,
    prime_shell_report,
    r4_theorem,
    run,
    shannon_entropy,
    uh_orbit_partition,
    verify_r4_theorem,
)


def test_hurwitz_norm_unit():
    assert hurwitz_norm(1, 0, 0, 0) == 1
    assert hurwitz_norm(1, 1, 1, 1) == 4


def test_enum_integer_solutions_p5():
    sols = enum_integer_solutions(5)
    assert (1, 2, 0, 0) in sols or (-1, -2, 0, 0) in sols
    assert all(hurwitz_norm(*s) == 5 for s in sols)


def test_r4_theorem_sample_primes():
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 97]:
        check = verify_r4_theorem(p)
        assert check["r4_matches"], f"r_4({p}) != 8({p}+1)"
        assert check["integer_reps"] == r4_theorem(p)


def test_shell_size_vs_r4():
    assert len(hurwitz_shell_elements(2)) == r4_theorem(2)
    for p in [3, 5, 7, 11, 13]:
        shell = len(hurwitz_shell_elements(p))
        r4 = r4_theorem(p)
        assert shell == 3 * r4, f"|Σ_{p}|={shell}, erwartet 3*r_4={3*r4}"


def test_shell_size_p2_matches_hurwitz_csv():
    assert len(hurwitz_shell_elements(2)) == 24
    assert len(hurwitz_orbit_elements(2)) == 24


def test_shell_size_p5():
    assert len(hurwitz_shell_elements(5)) == 144


def test_uh_orbit_partition_p7():
    orbits = uh_orbit_partition(7)
    assert len(orbits) == 2
    sizes = sorted(len(o) for o in orbits)
    assert sizes == [96, 96]


def test_uh_orbit_partition_p5():
    orbits = uh_orbit_partition(5)
    assert len(orbits) == 1
    assert len(orbits[0]) == 144


def test_kappa_leg_zero():
    assert kappa_leg(0) == "0"


def test_kappa_leg_odd():
    assert kappa_leg(1) == "E"
    assert kappa_leg(5) == "A"


def test_leg_signature_full():
    assert leg_signature(0) == (0, 0, "0")
    alpha, beta, ka = leg_signature(12)
    assert alpha >= 2 and beta >= 1 and ka in ("E", "A", "B", "C")


def test_full_gamma4_structure():
    fg = full_gamma4(1, 2, 0, 0)
    assert len(fg) == 4
    assert fg[2] == (0, 0, "0")
    assert all(len(leg) == 3 for leg in fg)


def test_prime_shell_report_has_full_gamma():
    report = prime_shell_report(7)
    assert "distinct_full_gamma" in report
    assert report["distinct_full_gamma"] >= report["distinct_gamma"]
    assert "M_p_full_gamma_top" in report
    sp = report["sample_points"][0]
    assert "full_gamma" in sp
    assert "full_gamma_compact" in sp


def test_chirality_balanced_four_e():
    assert chirality_score(("E", "E", "E", "E")) == 4
    assert chirality_score(("A", "B", "A", "B")) == -4
    assert chirality_score(("0", "E", "0", "C")) == 2


def test_prime_shell_report_structure():
    report = prime_shell_report(7)
    assert report["shell_size"] == 192
    assert report["orbit_size"] == 192
    assert report["r4_check"]["r4_matches"]
    assert "chi_p" in report
    assert "H_p" in report
    assert "K_p" in report
    assert "M_p" in report
    assert "uh_orbits" in report
    assert report["uh_orbits"]["orbit_count"] == 2
    assert "shuffle_null" in report
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
    scale = aggregate_shell_report(500)
    assert "H_stats" in scale
    assert scale["H_stats"]["mean"] > 0
    assert "H_trend" in scale
    assert "H_by_p_mod_12" in scale
    assert "r4_verification" in scale
    assert scale["r4_verification"]["all_primes_match"]
    assert "asymptotic_buckets" in scale
    assert "mu_infinity_hint" in scale["asymptotic_buckets"]


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "hurwitz_orbit.json"
    report = run(max_ps=[200], output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert "scales" in loaded
    assert loaded["scales"]["200"]["prime_count"] == report["scales"]["200"]["prime_count"]


def test_aggregate_has_shell_stats(tmp_path: Path):
    out = tmp_path / "hurwitz_orbit.json"
    report = run(max_ps=[500], output=out)
    scale = report["scales"]["500"]
    assert scale["shell_size_stats"]["min"] >= 24
    assert scale["total_shell_points"] > 0
    assert scale["orbit_count_stats"]["mean"] >= 1


def test_half_integer_contributes_for_some_primes():
    report = prime_orbit_report(3)
    assert report["half_integer_coords"] > 0


def test_prime_orbit_report_alias():
    assert prime_orbit_report(5)["shell_size"] == prime_shell_report(5)["shell_size"]
