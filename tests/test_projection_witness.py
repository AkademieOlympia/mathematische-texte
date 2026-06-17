from __future__ import annotations

import math
import unittest

from susy_fourlinge.witness import (
    Vec3,
    chirality_word,
    classify_by_projection,
    delta_weight,
    edge_projection,
    edge_weight,
    monotone_transition_weights,
    monotone_transition_witnesses,
    projection_sum_witness,
    quadruplet_witness,
    sieve_statistics,
)


class ProjectionWitnessTests(unittest.TestCase):
  def test_weights_at_t0(self) -> None:
    self.assertTrue(math.isclose(edge_weight("A", "B", 0.0, sign=1, start_edge=True), 2.5))
    self.assertTrue(math.isclose(edge_weight("C", "E", 0.0, sign=-1, start_edge=True), 1.5))
    self.assertTrue(math.isclose(delta_weight(0.0), 1.0))

  def test_projections_at_t0(self) -> None:
    p_ab = edge_projection("A", "B", 0.0, sign=1, start_edge=True)
    p_ce = edge_projection("C", "E", 0.0, sign=-1, start_edge=True)

    self.assertTrue(math.isclose(p_ab.x, -0.346410, rel_tol=0, abs_tol=1e-6))
    self.assertTrue(math.isclose(p_ab.y, 0.200000, rel_tol=0, abs_tol=1e-6))
    self.assertTrue(math.isclose(p_ab.z, -0.346410, rel_tol=0, abs_tol=1e-6))

    self.assertTrue(math.isclose(p_ce.x, 0.577350, rel_tol=0, abs_tol=1e-6))
    self.assertTrue(math.isclose(p_ce.y, -0.333333, rel_tol=0, abs_tol=1e-6))
    self.assertTrue(math.isclose(p_ce.z, 0.577350, rel_tol=0, abs_tol=1e-6))

    self.assertTrue(math.isclose(p_ce.norm() / p_ab.norm(), 5.0 / 3.0, rel_tol=0, abs_tol=1e-9))

  def test_sum_xz_witness_at_t0(self) -> None:
    p_ab = edge_projection("A", "B", 0.0, sign=1, start_edge=True)
    p_ce = edge_projection("C", "E", 0.0, sign=-1, start_edge=True)

    self.assertLess(projection_sum_witness(p_ab), 0.0)
    self.assertGreater(projection_sum_witness(p_ce), 0.0)
    self.assertEqual(classify_by_projection(p_ab), "ABCE")
    self.assertEqual(classify_by_projection(p_ce), "CEAB")

  def test_monotone_transition_at_pi_over_4(self) -> None:
    weights = monotone_transition_weights(math.pi / 4)
    expected = {
        "W_EA": 1.292893,
        "W_AB": 2.000000,
        "W_BC": 2.707107,
        "W_CE": 2.000000,
    }
    for key, value in expected.items():
      self.assertTrue(math.isclose(weights[key], value, rel_tol=0, abs_tol=1e-6), key)

  def test_monotone_projection_gate_at_pi_over_4(self) -> None:
    witnesses = monotone_transition_witnesses(math.pi / 4)
    norms = {label: w.projection.norm() for label, w in witnesses.items()}
    self.assertGreater(norms["EA"], norms["AB"])
    self.assertGreater(norms["EA"], norms["CE"])
    self.assertTrue(math.isclose(norms["AB"], norms["CE"], rel_tol=0, abs_tol=1e-9))
    self.assertLess(norms["BC"], norms["AB"])
    self.assertLess(norms["BC"], norms["CE"])

  def test_first_quadruplet_witnesses(self) -> None:
    cases = [
        (101, "ABCE"),
        (191, "CEAB"),
        (821, "ABCE"),
        (1871, "CEAB"),
    ]
    for p, word in cases:
      witness = quadruplet_witness(p)
      self.assertEqual(chirality_word(p), word)
      self.assertEqual(witness.word, word)
      self.assertEqual(classify_by_projection(witness.start_witness.projection), word)

  def test_sieve_statistics_1e6(self) -> None:
    stats = sieve_statistics(10**6)
    self.assertEqual(stats.total, 166)
    self.assertEqual(stats.abce, 84)
    self.assertEqual(stats.ceab, 82)
    self.assertTrue(math.isclose(stats.bias, 0.01205, rel_tol=0, abs_tol=1e-5))
    self.assertTrue(math.isclose(stats.signed_bias, 0.01205, rel_tol=0, abs_tol=1e-5))

  def test_sieve_statistics_1e7(self) -> None:
    stats = sieve_statistics(10**7)
    self.assertEqual(stats.total, 899)
    self.assertEqual(stats.abce, 450)
    self.assertEqual(stats.ceab, 449)
    self.assertTrue(math.isclose(stats.bias, 0.00111, rel_tol=0, abs_tol=1e-5))
    self.assertTrue(math.isclose(stats.signed_bias, 0.00111, rel_tol=0, abs_tol=1e-5))

  def test_sieve_statistics_1e8(self) -> None:
    stats = sieve_statistics(10**8)
    self.assertEqual(stats.total, 4768)
    self.assertEqual(stats.abce, 2408)
    self.assertEqual(stats.ceab, 2360)
    self.assertTrue(math.isclose(stats.bias, 0.01007, rel_tol=0, abs_tol=1e-5))
    self.assertTrue(math.isclose(stats.signed_bias, 0.01007, rel_tol=0, abs_tol=1e-5))


if __name__ == "__main__":
  unittest.main()
