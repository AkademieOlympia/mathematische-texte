"""Tests für diskreten EABC-Assoziator (collatz_eabc_discrete_associator.py)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz_eabc_discrete_associator import (
    ALL_CLASSES,
    ABC_CLASSES,
    associator_sign,
    bracket_preference,
    check_associativity,
    phi,
    phi_left,
    phi_right,
    prime_associator_mean,
    prove_v4_klein_associativity,
    quadruplet_chirality,
    run,
)
from eabc_from_lean import Chirality, EClass, class_of, is_prime_quadruplet, q, residue


def test_phi_matches_global_lokal_v4_table():
  """V₄-Tabelle wie Global Lokal.py."""
  expected = {
      ("E", "E"): "E", ("E", "A"): "A", ("E", "B"): "B", ("E", "C"): "C",
      ("A", "E"): "A", ("A", "A"): "E", ("A", "B"): "C", ("A", "C"): "B",
      ("B", "E"): "B", ("B", "A"): "C", ("B", "B"): "E", ("B", "C"): "A",
      ("C", "E"): "C", ("C", "A"): "B", ("C", "B"): "A", ("C", "C"): "E",
  }
  for (xs, ys), want in expected.items():
      x = EClass(xs)
      y = EClass(ys)
      assert phi(x, y).value == want


def test_phi_via_residue_mod12():
  assert phi(EClass.A, EClass.B) is EClass.C
  assert (residue(EClass.A) * residue(EClass.B)) % 12 == residue(EClass.C)


def test_associativity_full_v4():
  report = check_associativity(ALL_CLASSES)
  assert report["associative"]
  assert report["triples_tested"] == 64
  assert report["counterexample_count"] == 0


def test_associator_zero_on_abc_subtriple():
  report = check_associativity(ABC_CLASSES)
  assert report["associative"]
  for x in ABC_CLASSES:
      for y in ABC_CLASSES:
          for z in ABC_CLASSES:
              assert associator_sign(x, y, z) == 0
              assert phi_left(x, y, z) is phi_right(x, y, z)


def test_bracket_preference_chirality():
  assert bracket_preference(Chirality.ABCE) == "left"
  assert bracket_preference(Chirality.CEAB) == "right"


def test_lean_quadruplet_q5_abce():
  assert is_prime_quadruplet(5)
  assert quadruplet_chirality(5) is Chirality.ABCE
  word = "".join(class_of(n).value for n in q(5))
  assert word == "ABCE"


def test_v4_klein_proof():
  proof = prove_v4_klein_associativity()
  assert proof["associative"]
  assert proof["isomorphism_to_Z2_squared"]


def test_v4_klein_proof_in_run(tmp_path: Path):
  out = tmp_path / "disc.json"
  report = run(prime_limit=200, n_primes=20, output=out)
  assert report["v4_klein_proof"]["associative"]


def test_prime_associator_mean_zero():
  pa = prime_associator_mean(50)
  assert pa["A_N"] == 0.0
  assert pa["nonzero_count"] == 0


def test_run_writes_json(tmp_path: Path):
  out = tmp_path / "disc.json"
  report = run(prime_limit=500, n_primes=30, output=out)
  assert out.is_file()
  loaded = json.loads(out.read_text(encoding="utf-8"))
  assert loaded["associativity"]["full_V4"]["associative"]
  assert loaded["prime_associator"]["A_N"] == 0.0
  assert report["output_path"] == str(out)
