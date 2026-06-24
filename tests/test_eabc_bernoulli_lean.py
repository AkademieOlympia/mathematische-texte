from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_bernoulli_lean_test import (
    CHIRALITY_WORDS,
    IOTA_E,
    LEAN_SOURCES,
    bernoulli_cell_indices,
    classification_consistency,
    lean_class_counts,
    prime_sig_word,
    quadruplet_witnesses,
    radial_at,
    run_lean_coupling,
    t_rotate_word,
    theta4_word,
    theta_word,
)
from collatz_eabc_bernoulli_sensor import _sieve_primes, bernoulli_row, prime_sig, v_bernoulli
from eabc_from_lean import Chirality, EClass, class_of, is_prime_quadruplet, q, residue, t4


def test_lean_sources_reference_eabc_lean():
    assert any("EABC.lean" in s for s in LEAN_SOURCES)
    assert any("eabc_from_lean" in s for s in LEAN_SOURCES)


def test_lean_class_counts_matches_sensor_up_to_500():
    primes = _sieve_primes(1001)
    for n in range(1, 501):
        sig = prime_sig(2 * n, primes)
        assert lean_class_counts(sig).as_tuple() == v_bernoulli(sig).as_tuple()


def test_classification_consistency_small_n():
    primes = _sieve_primes(120)
    for n in range(1, 61):
        row = bernoulli_row(n, primes)
        cons = classification_consistency(row.prime_sig, row.v)
        assert cons["match"]
        assert cons["residue_roundtrip_ok"]
        assert cons["T4_identity_ok"]


def test_residue_and_t4_on_eabc_primes():
    for p in [5, 7, 11, 13, 17, 19, 23]:
        cls = class_of(p)
        assert cls is not None
        assert residue(cls) == p % 12
        assert t4(cls) is cls


def test_quadruplet_witness_at_5():
    sig = [2, 3, 5, 7, 11, 13]
    quads = quadruplet_witnesses(sig)
    assert len(quads) == 2  # p=5 und p=11 sind Vierlingsanfänge
    q5 = next(qd for qd in quads if qd["p"] == 5)
    assert q5["quadruplet"] == q(5)
    assert q5["classes"] == ["A", "B", "C", "E"]
    assert q5["word"] == "ABCE"


def test_prime_sig_word_skips_2_and_3():
    sig = [2, 3, 5, 7]
    assert prime_sig_word(sig) == "AB"


def test_chirality_words_constants():
    assert CHIRALITY_WORDS == ("ABCE", "CEAB")


def test_iota_e_aligns_with_residue():
    for cls in EClass:
        assert class_of(residue(cls)) is cls
        assert IOTA_E[cls] == {"E": 0, "A": 1, "B": 2, "C": 3}[cls.value]


def test_pref_projection_word_observables():
    w = "ABCE"
    assert theta_word(w) == sum(IOTA_E[c] for c in EClass if c.value in w)
    assert theta4_word(w) == theta_word(w) % 4
    assert t_rotate_word(w) == "BCEA"


def test_bernoulli_cell_and_radial_at():
    assert bernoulli_cell_indices(3) == (4, 6, 8)
    assert radial_at(3) == 0.125


def test_run_lean_coupling_summary():
    report = run_lean_coupling(100, null_trials=50)
    s = report["summary"]
    assert s["classification_match_all"]
    assert s["residue_roundtrip_all"]
    assert s["T4_identity_all"]
    assert report["experiment"] == "Lean-Kopplung"
    assert report["python_bridge"] == "eabc_from_lean.py"
    assert len(report["lean_sources"]) >= 3
    assert "note" in report and "Lean" in report["note"]


def test_run_lean_coupling_json_serializable(tmp_path: Path):
    report = run_lean_coupling(20, null_trials=20)
    out = tmp_path / "lean.json"
    out.write_text(json.dumps(report), encoding="utf-8")
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["max_n"] == 20
    assert loaded["summary"]["classification_match_all"]


def test_is_prime_quadruplet_5_in_sig_of_n6():
    primes = _sieve_primes(20)
    sig = prime_sig(12, primes)
    assert 5 in sig
    assert is_prime_quadruplet(5)
    quads = quadruplet_witnesses(sig)
    assert any(qd["p"] == 5 for qd in quads)


def test_chirality_reference_in_report():
    report = run_lean_coupling(5, null_trials=10)
    ref = report["chirality_reference"]
    assert ref["ABCE"] == "ABCE"
    assert ref["CEAB"] == "CEAB"
    for chir in Chirality:
        assert chir.value in ref
