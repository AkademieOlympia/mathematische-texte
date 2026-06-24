"""Tests für EABC-Oktanion-Associator (collatz_eabc_oktonion_associator.py)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_oktonion_associator import (
    associator,
    associator_norm,
    canonical_triples_test,
    eabc_associator_norm,
    o_mul,
    o_norm_sq,
    run,
    shell_size_z8,
    triple_factorizations,
)


def test_norm_multiplicativity():
    e1 = (0, 1, 0, 0, 0, 0, 0, 0)
    e2 = (0, 0, 1, 0, 0, 0, 0, 0)
    assert o_norm_sq(o_mul(e1, e2)) == o_norm_sq(e1) * o_norm_sq(e2)


def test_associator_zero_quaternion_subalgebra():
    e1 = (0, 1, 0, 0, 0, 0, 0, 0)
    e2 = (0, 0, 1, 0, 0, 0, 0, 0)
    e3 = (0, 0, 0, 1, 0, 0, 0, 0)
    assert associator_norm(e1, e2, e3) < 1e-9
    assert associator(e1, e2, e3) == (0,) * 8


def test_associator_nonzero_generic_octonion():
    e1 = (0, 1, 0, 0, 0, 0, 0, 0)
    e2 = (0, 0, 1, 0, 0, 0, 0, 0)
    e4 = (0, 0, 0, 0, 1, 0, 0, 0)
    assert associator_norm(e1, e2, e4) > 0.0
    canon = canonical_triples_test()
    assert canon["generic_associator_nonzero"]
    assert canon["quaternion_associator_zero"]


def test_eabc_associator_nonzero_on_shell_sample():
    """Glatt-Γ-Differenz: auf Einheitsbasis oft 0; auf Σ_2×Σ_2×Σ_3 nichttrivial."""
    x = (-1, -1, 0, 0, 0, 0, 0, 0)
    y = (-1, 0, 0, 0, -1, 0, 0, 0)
    z = (-1, 0, -1, -1, 0, 0, 0, 0)
    assert o_norm_sq(x) == 2 and o_norm_sq(y) == 2 and o_norm_sq(z) == 3
    assert eabc_associator_norm(x, y, z) > 0.0


def test_shell_size_matches_stub():
    assert shell_size_z8(1) == 16
    assert shell_size_z8(2) == 112


def test_triple_factorizations_12():
    facs = triple_factorizations(12)
    assert (2, 2, 3) in facs


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "assoc.json"
    report = run(max_n=20, samples=5, seed=1, output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["meta"]["max_n"] == 20
    assert report["canonical_tests"]["generic_associator_nonzero"]
    assert report["output_path"] == str(out)
