"""Pure Python replacement for the old Sage-based phase experiment."""

from __future__ import annotations

import cmath
from dataclasses import dataclass
from math import sqrt
from statistics import mean
from typing import Iterable


DEFAULT_QUADRUPLET = (5, 7, 11, 13)


@dataclass(frozen=True, slots=True)
class PhasePoint:
    prime: int
    angle: float
    phase: complex


@dataclass(frozen=True, slots=True)
class GeometrySummary:
    quadruplet: tuple[int, ...]
    phi: float
    gaps: tuple[int, ...]
    ratio_mid_to_outer_gap: float
    phi_squared: float
    points: tuple[PhasePoint, ...]
    mean_real: float
    mean_imag: float


def _validate_quadruplet(values: Iterable[int]) -> tuple[int, ...]:
    quadruplet = tuple(int(v) for v in values)
    if len(quadruplet) < 2:
        raise ValueError("Es werden mindestens zwei Werte benoetigt")
    if tuple(sorted(quadruplet)) != quadruplet:
        raise ValueError("Die Werte muessen aufsteigend sortiert sein")
    return quadruplet


def check_quadruplet_geometry(quadruplet: Iterable[int] = DEFAULT_QUADRUPLET) -> GeometrySummary:
    quad = _validate_quadruplet(quadruplet)
    phi = (1.0 + sqrt(5.0)) / 2.0
    gaps = tuple(right - left for left, right in zip(quad, quad[1:]))
    outer_gap = gaps[0]
    if outer_gap == 0:
        raise ValueError("Der erste Gap darf nicht 0 sein")

    points = tuple(
        PhasePoint(prime=prime, angle=prime / phi, phase=cmath.exp(1j * (prime / phi)))
        for prime in quad
    )

    return GeometrySummary(
        quadruplet=quad,
        phi=phi,
        gaps=gaps,
        ratio_mid_to_outer_gap=(gaps[1] / outer_gap) if len(gaps) > 1 else 0.0,
        phi_squared=phi**2,
        points=points,
        mean_real=mean(point.phase.real for point in points),
        mean_imag=mean(point.phase.imag for point in points),
    )
