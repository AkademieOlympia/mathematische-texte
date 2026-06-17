"""CLI for the SUSY quadruplet phase experiment."""

from __future__ import annotations

import argparse

from .core import check_quadruplet_geometry
from .witness import (
    centered_prime_quadruplet,
    prime_quadruplet_ellipse_parameters,
    sieve_statistics,
    test_perihel_real_vs_random,
    witness_ellipse_bridge,
)


def _cmd_bridge(args: argparse.Namespace) -> None:
    bridge = witness_ellipse_bridge(args.p, t=args.t)
    e = bridge.ellipse
    centered = centered_prime_quadruplet(bridge.start)
    params = prime_quadruplet_ellipse_parameters()
    print(f"Start p={bridge.start}  M={bridge.center}  word={bridge.word}")
    print(f"Zentriert: {centered}")
    print(f"Ellipse: a={e.a} b={e.b} f={e.f:.6f} e={e.e:.6f} steps={e.step_pattern}")
    print(
        f"Kanonisch: a_pv={params['a_pv']} b_pv={params['b_pv']} "
        f"e_pv={params['e_pv']:.6f} rho_pv={params['rho_pv']}"
    )
    print(f"Witness: S=x+z={bridge.sum_xz:.6f}  rho_PV={bridge.rho_pv:.6f}  t={args.t}")
    print(f"Startkante: {bridge.witness.start_edge}  W={bridge.witness.start_witness.weight:.6f}")


def _cmd_sieve(args: argparse.Namespace) -> None:
    stats = sieve_statistics(args.limit)
    print(f"N={stats.limit}  total={stats.total}  ABCE={stats.abce}  CEAB={stats.ceab}")
    print(f"signed_bias={stats.signed_bias:.5f}  bias={stats.bias:.5f}")


def _cmd_perihel(args: argparse.Namespace) -> None:
    stats = test_perihel_real_vs_random(args.limit)
    print(
        f"N={stats.limit}  samples={stats.sample_count}  "
        f"ideal_|Pi|_max={stats.ideal_zero_max_abs:.2e}"
    )
    print(
        f"real_mean_Pi={stats.real_mean_pi:.6f}  real_std={stats.real_std_pi:.6f}"
    )
    print(
        f"random_mean_Pi={stats.random_mean_pi:.6f}  random_std={stats.random_std_pi:.6f}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Analysiert ein Primzahl-Vierling-Phasenmodell ohne SageMath."
    )
    sub = parser.add_subparsers(dest="command")

    bridge = sub.add_parser("bridge", help="Witness-Ellipse-Brücke für Startpunkt p")
    bridge.add_argument("p", type=int)
    bridge.add_argument("--t", type=float, default=0.0)
    bridge.set_defaults(handler=_cmd_bridge)

    sieve = sub.add_parser("sieve", help="ABCE/CEAB-Siebstatistik bis N")
    sieve.add_argument("limit", type=int)
    sieve.set_defaults(handler=_cmd_sieve)

    perihel = sub.add_parser(
        "perihel",
        help="Numerischer Vergleich Π_real vs Π_random (Evidenz, kein Beweis)",
    )
    perihel.add_argument("limit", type=int, nargs="?", default=10**5)
    perihel.set_defaults(handler=_cmd_perihel)

    parser.add_argument(
        "quadruplet",
        nargs="*",
        type=int,
        help="Optional vier aufsteigende Zahlen, z. B. 5 7 11 13 (Legacy-Modus).",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command is not None:
        args.handler(args)
        return

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


if __name__ == "__main__":
    main()
