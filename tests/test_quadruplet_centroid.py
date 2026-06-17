from __future__ import annotations

import unittest

from susy_fourlinge.witness import (
    centered_prime_quadruplet,
    centroid_error,
    centroid_statistics,
    prime_quadruplet_ellipse_parameters,
    quadruplet_center,
)


class QuadrupletCentroidTests(unittest.TestCase):
  def test_center_is_p_plus_four(self) -> None:
    self.assertEqual(quadruplet_center(5), 9)
    self.assertEqual(quadruplet_center(191), 195)

  def test_centroid_error_zero_for_first_quadruplets(self) -> None:
    for p in (5, 11, 101, 191, 821, 1871):
      self.assertEqual(centroid_error(p), 0.0)

  def test_centroid_statistics_1e6(self) -> None:
    stats = centroid_statistics(10**6)
    self.assertEqual(stats.total, 166)
    self.assertEqual(stats.abce, 84)
    self.assertEqual(stats.ceab, 82)
    self.assertEqual(stats.max_centroid_error, 0.0)

  def test_centered_prime_quadruplet_normal_form(self) -> None:
    self.assertEqual(centered_prime_quadruplet(5), (-4, -2, 2, 4))
    self.assertEqual(centered_prime_quadruplet(101), (-4, -2, 2, 4))

  def test_prime_quadruplet_ellipse_parameters(self) -> None:
    params = prime_quadruplet_ellipse_parameters()
    self.assertEqual(params["a_pv"], 4.0)
    self.assertEqual(params["b_pv"], 2.0)
    self.assertLess(abs(params["e_pv"] - (3 ** 0.5) / 2), 1e-12)
    self.assertEqual(params["rho_pv"], 1.5)


if __name__ == "__main__":
  unittest.main()
