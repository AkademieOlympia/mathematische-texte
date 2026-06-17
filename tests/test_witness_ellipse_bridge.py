from __future__ import annotations

import math
import unittest

from susy_fourlinge.witness import (
    PrimeSieve,
    RHO_PV,
    STEP_PATTERN,
    QUADRUPLET_OFFSETS,
    canonical_ellipse_params,
    chirality_word,
    classify_by_projection,
    quadruplet_center,
    quadruplet_members,
    quadruplet_normal_form,
    witness_ellipse_bridge,
)


class WitnessEllipseBridgeTests(unittest.TestCase):
  def test_canonical_ellipse_params(self) -> None:
    params = canonical_ellipse_params()
    self.assertEqual(params.a, 4.0)
    self.assertEqual(params.b, 2.0)
    self.assertTrue(math.isclose(params.f, 2.0 * math.sqrt(3.0)))
    self.assertTrue(math.isclose(params.e, math.sqrt(3.0) / 2.0))
    self.assertLess(params.e, 1.0)
    self.assertEqual(params.step_pattern, STEP_PATTERN)
    self.assertEqual(params.offsets, QUADRUPLET_OFFSETS)

  def test_normal_form_and_members(self) -> None:
    center, offsets = quadruplet_normal_form(101)
    self.assertEqual(center, 105)
    self.assertEqual(offsets, (-4, -2, 2, 4))
    members = quadruplet_members(101)
    self.assertEqual(members, (101, 103, 107, 109))
    self.assertEqual(tuple(center + o for o in offsets), members)

  def test_bridge_at_t0(self) -> None:
    bridge = witness_ellipse_bridge(101, t=0.0)
    self.assertEqual(bridge.start, 101)
    self.assertEqual(bridge.center, quadruplet_center(101))
    self.assertEqual(bridge.word, "ABCE")
    self.assertTrue(math.isclose(bridge.rho_pv, RHO_PV))
    self.assertTrue(math.isclose(bridge.rho_pv, 1.5))
    self.assertLess(bridge.sum_xz, 0.0)
    self.assertTrue(math.isclose(bridge.ellipse.e, math.sqrt(3.0) / 2.0))
    self.assertNotEqual(bridge.rho_pv, bridge.ellipse.e)
    self.assertEqual(
        classify_by_projection(bridge.witness.start_witness.projection),
        bridge.word,
    )

  def test_bridge_orientation_ceab(self) -> None:
    bridge = witness_ellipse_bridge(191, t=0.0)
    self.assertEqual(bridge.word, "CEAB")
    self.assertGreater(bridge.sum_xz, 0.0)
    self.assertEqual(canonical_ellipse_params(), bridge.ellipse)

  def test_canonical_params_invariant_up_to_1e6(self) -> None:
    canonical = canonical_ellipse_params()
    sieve = PrimeSieve(10**6)
    for p in sieve.iter_quadruplet_starts():
      bridge = witness_ellipse_bridge(p)
      self.assertEqual(bridge.ellipse, canonical)
      self.assertEqual(bridge.center, p + 4)
      self.assertEqual(classify_by_projection(bridge.witness.start_witness.projection), chirality_word(p))

  def test_first_examples(self) -> None:
    cases = [
        (101, "ABCE"),
        (191, "CEAB"),
        (821, "ABCE"),
        (1871, "CEAB"),
    ]
    for p, word in cases:
      bridge = witness_ellipse_bridge(p)
      self.assertEqual(bridge.word, word)
      self.assertTrue(math.isclose(bridge.rho_pv, 1.5))


if __name__ == "__main__":
  unittest.main()
