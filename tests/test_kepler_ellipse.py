from __future__ import annotations

import math
import unittest

from susy_fourlinge.witness import (
    kepler_ellipse_point,
    kepler_phase_tick,
    quadruplet_center,
)


class KeplerEllipseTests(unittest.TestCase):
    def test_abce_vs_ceab_at_zero_differ_by_pi_shift(self) -> None:
        p_abce, p_ceab = 101, 191
        m_abce = quadruplet_center(p_abce)
        m_ceab = quadruplet_center(p_ceab)

        e_plus_0 = kepler_ellipse_point(p_abce, 0.0)
        e_minus_0 = kepler_ellipse_point(p_ceab, 0.0)

        self.assertTrue(math.isclose(e_plus_0[0], m_abce + 4.0))
        self.assertTrue(math.isclose(e_plus_0[1], 0.0, abs_tol=1e-9))
        self.assertTrue(math.isclose(e_minus_0[0], m_ceab - 4.0))
        self.assertTrue(math.isclose(e_minus_0[1], 0.0, abs_tol=1e-9))

        # CEAB bei θ=0 entspricht der ABCE-Parametrisierung bei θ=π (gleicher Anker M).
        m = m_ceab
        expected_pi = (m + 4.0 * math.cos(math.pi), 2.0 * math.sin(math.pi))
        self.assertTrue(math.isclose(e_minus_0[0], expected_pi[0]))
        self.assertTrue(math.isclose(e_minus_0[1], expected_pi[1], abs_tol=1e-9))

    def test_discrete_phase_clock(self) -> None:
        expected = [0.0, math.pi / 2.0, math.pi, 3.0 * math.pi / 2.0]
        for t, theta in enumerate(expected):
            self.assertTrue(math.isclose(kepler_phase_tick(t), theta))

    def test_phase_clock_on_ellipse(self) -> None:
        p = 101
        flavors = ("E", "A", "B", "C")
        for t, flavor in enumerate(flavors):
            theta = kepler_phase_tick(t)
            x, y = kepler_ellipse_point(p, theta)
            if flavor == "E":
                self.assertTrue(math.isclose(x, quadruplet_center(p) + 4.0))
                self.assertTrue(math.isclose(y, 0.0, abs_tol=1e-9))
            elif flavor == "A":
                self.assertTrue(math.isclose(x, quadruplet_center(p)))
                self.assertTrue(math.isclose(y, 2.0))
            elif flavor == "B":
                self.assertTrue(math.isclose(x, quadruplet_center(p) - 4.0))
                self.assertTrue(math.isclose(y, 0.0, abs_tol=1e-9))
            elif flavor == "C":
                self.assertTrue(math.isclose(x, quadruplet_center(p)))
                self.assertTrue(math.isclose(y, -2.0))


if __name__ == "__main__":
    unittest.main()
