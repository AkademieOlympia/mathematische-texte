from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from collatz_praezession_test import (
    g_from_span,
    information_content,
    pearson,
    perihel_proxy_uniform,
    spearman,
)


def test_g_from_span_minimal():
    assert g_from_span(10) == 0
    assert g_from_span(22) == 1
    assert g_from_span(34) == 2


def test_information_content_monotone():
    assert information_content(0) == 0.0
    assert information_content(3) > information_content(1)


def test_perihel_proxy_decreases_with_p():
    pi_small = perihel_proxy_uniform(13)
    pi_large = perihel_proxy_uniform(1000)
    assert pi_small > pi_large
    assert math.isclose(pi_small, 4.0 / math.log(13))


def test_pearson_perfect_correlation():
    xs = [1.0, 2.0, 3.0, 4.0]
    assert math.isclose(pearson(xs, xs) or 0.0, 1.0)
    assert math.isclose(pearson(xs, [-x for x in xs]) or 0.0, -1.0)


def test_spearman_monotone():
    xs = [1.0, 2.0, 3.0, 5.0]
    ys = [10.0, 20.0, 30.0, 50.0]
    r = spearman(xs, ys)
    assert r is not None
    assert abs(r - 1.0) < 1e-12
