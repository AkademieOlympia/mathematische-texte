from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_bernoulli_sensor import (
    EabcVector,
    bernoulli_row,
    chi_eabc,
    delta_q4,
    i_chir,
    non_eabc_primes,
    prime_sig,
    q4_vector,
    run_sensor,
    sigma_eabc,
    staudt_denominator,
    v_bernoulli,
    verify_staudt_sympy,
    _sieve_primes,
)


def test_staudt_denominator_small_n():
    """von Staudt--Clausen: den(B_{2n}) = ∏_{p-1|2n} p."""
    expected = {
        1: 6,  # B_2 = 1/6, 2n=2
        2: 30,  # B_4, 2n=4
        3: 42,  # B_6, 2n=6
        4: 30,  # B_8, 2n=8
        5: 66,  # B_10, 2n=10
    }
    primes = _sieve_primes(20)
    for n, den in expected.items():
        assert staudt_denominator(2 * n, primes) == den


def test_staudt_matches_sympy_up_to_20():
    checks = verify_staudt_sympy(20)
    assert checks
    assert all(c["match"] for c in checks)


def test_prime_sig_counts():
    primes = _sieve_primes(20)
    assert prime_sig(2, primes) == [2, 3]
    assert prime_sig(4, primes) == [2, 3, 5]
    assert prime_sig(6, primes) == [2, 3, 7]
    assert prime_sig(12, primes) == [2, 3, 5, 7, 13]


def test_v_sums_match_prime_sig():
    """e+a+b+c + |non-EABC| = |PrimeSig| für n=1..50."""
    primes = _sieve_primes(120)
    for n in range(1, 51):
        two_n = 2 * n
        sig = prime_sig(two_n, primes)
        vec = v_bernoulli(sig)
        extra = len(non_eabc_primes(sig))
        assert vec.total + extra == len(sig)


def test_v_known_values_n1_to_5():
    primes = _sieve_primes(20)
    expected = {
        1: (0, 0, 0, 0),
        2: (0, 1, 0, 0),
        3: (0, 0, 1, 0),
        4: (0, 1, 0, 0),
        5: (0, 0, 0, 1),
    }
    for n, tup in expected.items():
        row = bernoulli_row(n, primes)
        assert row.v.as_tuple() == tup


def test_sigma_chi_i_chir_n2():
    v = EabcVector(0, 1, 0, 0)
    assert sigma_eabc(v) == 1
    assert chi_eabc(v) == -1
    assert i_chir(v) == -1


def test_run_sensor_structure(tmp_path: Path):
    out = tmp_path / "sensor.json"
    report = run_sensor(10)
    assert report["framework"] == "EABC"
    assert report["sensor"] == "Phi"
    assert report["hypothesis"] == "EABC-Resonanzhypothese der Zetafunktion"
    assert len(report["samples"]) == 10
    sample0 = report["samples"][0]
    assert sample0["E_n"] == sample0["V"]["E"]
    assert "stats" in report
    out.write_text(json.dumps(report), encoding="utf-8")
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["max_n"] == 10


def test_q4_vector_counts_primes_up_to_n():
    """Q_4(N) zählt EABC-Klassen über alle p ≤ N."""
    primes = _sieve_primes(50)
    v = q4_vector(30, primes)
    # p≤30 in EABC: 5(A),7(B),11(C),13(A),17(E),19(A),23(E),29(A)
    assert v.as_tuple() == (1, 3, 2, 2)
    d = delta_q4(v)
    assert d["sigma"] == sigma_eabc(v)
    assert d["chi"] == chi_eabc(v)
    assert d["a_minus_c"] == v.a - v.c
