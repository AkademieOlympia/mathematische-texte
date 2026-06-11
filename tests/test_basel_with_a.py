from __future__ import annotations

import math
import unittest

from basel_with_a.core import analyze_first_primes, first_n_primes, nth_prime


class BaselWithATests(unittest.TestCase):
    def test_first_n_primes(self) -> None:
        self.assertEqual(first_n_primes(5), [2, 3, 5, 7, 11])

    def test_nth_prime(self) -> None:
        self.assertEqual(nth_prime(100), 541)

    def test_analysis(self) -> None:
        summary = analyze_first_primes(5)

        self.assertEqual(summary.prime_count, 5)
        self.assertEqual(summary.nth_prime, 11)
        self.assertEqual(summary.prime_square_sum, 208)
        self.assertTrue(math.isclose(summary.difference, 208 - (math.pi**2) / 6.0))
        self.assertTrue(math.isclose(summary.index_value_a, math.sqrt(208)))


if __name__ == "__main__":
    unittest.main()
