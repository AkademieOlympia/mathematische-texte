"""Tests für EABC-Oktanion-Assoziator-Spektrum (collatz_eabc_oktonion_spectrum.py)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_oktonion_spectrum import (
    alpha_algebraic,
    build_histogram,
    kl_divergence,
    normalize_hist,
    run,
    spectrum_for_n,
    triple_norm_factorizations,
)


def test_alpha_algebraic_unit_triple():
    e1 = (0, 1, 0, 0, 0, 0, 0, 0)
    e2 = (0, 0, 1, 0, 0, 0, 0, 0)
    e4 = (0, 0, 0, 0, 1, 0, 0, 0)
    assert alpha_algebraic(e1, e2, e4) == 4  # ||[e1,e2,e4]||^2 = 2^2


def test_triple_factorizations_6():
    facs = triple_norm_factorizations(6)
    assert (1, 1, 6) in facs
    assert (1, 2, 3) in facs


def test_kl_self_zero():
    h = {"0": 5, "4": 3}
    p = normalize_hist(h)
    assert kl_divergence(p, p) is not None
    assert abs(kl_divergence(p, p) or 1.0) < 1e-9


def test_spectrum_for_n_structure():
    import random

    row = spectrum_for_n(6, samples=10, rng=random.Random(1))
    assert row["n"] == 6
    assert not row["is_prime"]
    assert "M_n" in row
    assert "M_n_E" in row
    assert row["sample_triples"] > 0


def test_run_writes_json(tmp_path: Path):
    out = tmp_path / "spectrum.json"
    report = run(max_n=12, samples=8, seed=2, output=out)
    assert out.is_file()
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["meta"]["max_n"] == 12
    assert "falsification" in loaded
    assert "prime_vs_composite_kl" in loaded
    assert report["output_path"] == str(out)


def test_build_histogram():
    assert build_histogram([0, 0, 4]) == {"0": 2, "4": 1}
