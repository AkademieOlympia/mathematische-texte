from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from eabc_from_lean import EClass, class_of
from collatz_eabc_invarianzprogramm import (
    chi_at_x,
    kappa,
    pi_eabc,
    quadruplet_signature_frequencies,
    s_at_x,
    sigma_quadruplet,
    snapshot_at_x,
    v_at_x,
)


def test_kappa_matches_class_of_for_primes_above_3():
    for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        assert kappa(p) is class_of(p)
        assert kappa(p) is not None


def test_kappa_none_for_two_and_three():
    assert kappa(2) is None
    assert kappa(3) is None


def test_s_at_x_sums_to_one():
    for x in [20, 100, 500, 1000, 5000]:
        s = s_at_x(x)
        total = s.e + s.a + s.b + s.c
        assert abs(total - 1.0) < 1e-12


def test_v_at_x_known_small_values():
    # p ≤ 10 in P_{>3}: 5→A, 7→B
    v10 = v_at_x(10)
    assert v10.as_tuple() == (0, 1, 1, 0)
    assert pi_eabc(10, v10) == 2

    # p ≤ 20: +11→C, 13→E, 17→A, 19→B
    v20 = v_at_x(20)
    assert v20.as_tuple() == (1, 2, 2, 1)
    assert pi_eabc(20, v20) == 6

    # p ≤ 100: π_{>3}(100)=23, Verteilung (5,6,6,6)
    v100 = v_at_x(100)
    assert v100.as_tuple() == (5, 6, 6, 6)
    assert pi_eabc(100, v100) == 23


def test_chi_at_x_known_values():
    v20 = v_at_x(20)
    # χ = ((E+C)-(A+B))/π = ((1+1)-(2+2))/6 = -2/6
    assert abs(chi_at_x(20, v20) - (-1 / 3)) < 1e-12

    v100 = v_at_x(100)
    # χ = ((5+6)-(6+6))/23 = -1/23
    assert abs(chi_at_x(100, v100) - (-1 / 23)) < 1e-12


def test_sigma_quadruplet_first_example():
    assert sigma_quadruplet(5) == "ABCE"
    assert sigma_quadruplet(7) is None


def test_quadruplet_frequencies_include_abce():
    stats = quadruplet_signature_frequencies(100)
    assert stats["quadruplet_count"] >= 1
    assert "ABCE" in stats["signature_counts"]


def test_snapshot_structure():
    snap = snapshot_at_x(100)
    assert snap.x == 100
    assert snap.v.total == snap.pi
    assert abs(snap.s.e + snap.s.a + snap.s.b + snap.s.c - 1.0) < 1e-12
    assert snap.chi == chi_at_x(100, snap.v)


def test_kappa_residue_classes():
    assert kappa(13) is EClass.E
    assert kappa(17) is EClass.A
    assert kappa(19) is EClass.B
    assert kappa(23) is EClass.C
