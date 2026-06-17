from __future__ import annotations

import math
import unittest

from susy_fourlinge.witness import (
    ideal_phase,
    perihel_shift,
    perihel_shift_uniform_epsilon,
    test_perihel_real_vs_random,
)


class PerihelDriftTests(unittest.TestCase):
    def test_ideal_perihel_shift_zero(self) -> None:
        for m in range(128):
            pi_val = perihel_shift(ideal_phase(m), ideal_phase(m + 4))
            self.assertTrue(math.isclose(pi_val, 0.0, abs_tol=1e-12), m)

    def test_uniform_perihel_shift_equals_four_epsilon(self) -> None:
        epsilon = 0.017
        self.assertTrue(
            math.isclose(perihel_shift_uniform_epsilon(epsilon), 4 * epsilon, abs_tol=1e-15)
        )

    def test_real_vs_random_witness_runs(self) -> None:
        stats = test_perihel_real_vs_random(10**4)
        self.assertGreater(stats.sample_count, 0)
        self.assertLess(stats.ideal_zero_max_abs, 1e-12)
        self.assertGreater(stats.real_mean_pi, 0.0)
        self.assertTrue(math.isfinite(stats.random_mean_pi))


if __name__ == "__main__":
    unittest.main()
