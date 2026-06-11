from __future__ import annotations

import math
import unittest

from susy_fourlinge.core import check_quadruplet_geometry


class SusyFourlingeTests(unittest.TestCase):
    def test_default_quadruplet_summary(self) -> None:
        summary = check_quadruplet_geometry()

        self.assertEqual(summary.quadruplet, (5, 7, 11, 13))
        self.assertEqual(summary.gaps, (2, 4, 2))
        self.assertTrue(math.isclose(summary.ratio_mid_to_outer_gap, 2.0))
        self.assertEqual(len(summary.points), 4)

    def test_requires_sorted_values(self) -> None:
        with self.assertRaises(ValueError):
            check_quadruplet_geometry((7, 5, 11, 13))


if __name__ == "__main__":
    unittest.main()
