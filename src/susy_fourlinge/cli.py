"""CLI for the SUSY quadruplet phase experiment."""

from __future__ import annotations

import argparse

from .core import check_quadruplet_geometry


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Analysiert ein Primzahl-Vierling-Phasenmodell ohne SageMath."
    )
    parser.add_argument(
        "quadruplet",
        nargs="*",
        type=int,
        help="Optional vier aufsteigende Zahlen, z. B. 5 7 11 13.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    quadruplet = args.quadruplet or None
    summary = check_quadruplet_geometry(quadruplet or (5, 7, 11, 13))

    print(f"Vierling: {summary.quadruplet}")
    print(f"phi: {summary.phi:.12f}")
    print(f"Gaps: {summary.gaps}")
    print(f"mittlerer/aeusserer Gap: {summary.ratio_mid_to_outer_gap:.12f}")
    print(f"phi^2: {summary.phi_squared:.12f}")
    print(f"Mittelwert Phase (Re): {summary.mean_real:.12f}")
    print(f"Mittelwert Phase (Im): {summary.mean_imag:.12f}")
    print()
    print("Einzelphasen")
    for point in summary.points:
        print(
            f"p={point.prime:>4} angle={point.angle:>12.6f} "
            f"phase=({point.phase.real:>10.6f}, {point.phase.imag:>10.6f})"
        )
