#!/usr/bin/env python3
"""
Reproduzierbare Analyse des chiralen Flusses S(N) = #ABCE(N) - #CEAB(N)
für markierte Primzahlvierlinge (#Energiedoku).

Ausführung:
    python3 energiedoku_chiraler_fluss.py
    python3 energiedoku_chiraler_fluss.py --limits 1e5 1e6 1e7 --plot
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from src.susy_fourlinge.witness import PrimeSieve, chirality_word, sieve_statistics  # noqa: E402


@dataclass(frozen=True, slots=True)
class ChiralRow:
    limit: int
    abce: int
    ceab: int
    total: int
    flux: int
    signed_bias: float
    abs_flux_ratio: float
    z_score: float


def eabc_operator() -> np.ndarray:
    """Quaternionischer EABC-Operator (4×4) aus Durchbruch.py."""
    i = 1j
    return np.array(
        [
            [0, 1, 0, -i],
            [1, 0, 1, 0],
            [0, 1, 0, 1],
            [i, 0, 1, 0],
        ],
        dtype=complex,
    )


def verify_eabc_spectrum() -> dict[str, float]:
    """Prüft ±2cos(π/8) und ±2cos(3π/8) im 4×4-EABC-Spektrum."""
    vals = np.linalg.eigvals(eabc_operator())
    targets = {
        "+2cos(pi/8)": 2 * math.cos(math.pi / 8),
        "-2cos(pi/8)": -2 * math.cos(math.pi / 8),
        "+2cos(3pi/8)": 2 * math.cos(3 * math.pi / 8),
        "-2cos(3pi/8)": -2 * math.cos(3 * math.pi / 8),
    }
    return {name: min(abs(v - t) for v in vals) for name, t in targets.items()}


def _liouville_table(limit: int) -> np.ndarray:
    spf = np.arange(limit + 10, dtype=np.int32)
    for i in range(2, int(math.isqrt(limit)) + 1):
        if spf[i] == i:
            for j in range(i * i, limit + 10, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def liouville(n: int, spf: np.ndarray) -> int:
    omega = 0
    while n > 1:
        p = int(spf[n])
        while n % p == 0:
            n //= p
            omega += 1
    return -1 if omega % 2 else 1


def liouville_orthogonality(limit: int) -> dict[str, complex | float | int]:
    """Σχ_H und Σλ(p)χ_H(p) mit χ_H(ABCE)=i, χ_H(CEAB)=-i."""
    phase = {"ABCE": 1j, "CEAB": -1j}
    spf = _liouville_table(limit)
    sieve = PrimeSieve(limit)
    starts = sieve.iter_quadruplet_starts()
    sum_chi = 0j
    sum_lam_chi = 0j
    for p in starts:
        chi = phase[chirality_word(p)]
        sum_chi += chi
        sum_lam_chi += liouville(p, spf) * chi
    total = len(starts)
    return {
        "Q": total,
        "sum_chi": sum_chi,
        "sum_lam_chi": sum_lam_chi,
        "sum_chi_over_Q": sum_chi / total if total else 0j,
        "sum_lam_chi_over_Q": sum_lam_chi / total if total else 0j,
    }


def chiral_row(limit: int) -> ChiralRow:
    stats = sieve_statistics(limit)
    flux = stats.abce - stats.ceab
    total = stats.total
    z = flux / math.sqrt(total) if total else 0.0
    return ChiralRow(
        limit=limit,
        abce=stats.abce,
        ceab=stats.ceab,
        total=total,
        flux=flux,
        signed_bias=stats.signed_bias,
        abs_flux_ratio=abs(flux) / total if total else float("nan"),
        z_score=z,
    )


def default_limits() -> list[int]:
    return [10**5, 10**6, 5 * 10**6, 10**7, 10**8]


def analyze_limits(limits: list[int]) -> list[ChiralRow]:
    return [chiral_row(n) for n in limits]


def write_csv(rows: list[ChiralRow], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "N",
                "Q",
                "ABCE",
                "CEAB",
                "S",
                "B_sgn",
                "|S|/Q",
                "z",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.limit,
                    row.total,
                    row.abce,
                    row.ceab,
                    row.flux,
                    f"{row.signed_bias:.8f}",
                    f"{row.abs_flux_ratio:.8f}",
                    f"{row.z_score:.4f}",
                ]
            )


def maybe_plot(rows: list[ChiralRow], path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib nicht verfügbar — Plot übersprungen.")
        return

    ns = [r.limit for r in rows]
    ratios = [r.abs_flux_ratio for r in rows]
    signed = [r.signed_bias for r in rows]

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.loglog(ns, ratios, "o-", color="#204080", label=r"$|S(N)|/Q(N)$")
    ax.loglog(ns, [abs(s) for s in signed], "s--", color="#c06010", label=r"$|B_{\mathrm{sgn}}(N)|$")
    ax.set_xlabel(r"Siebgrenze $N$")
    ax.set_ylabel("normierter Fluss")
    ax.set_title(r"Chiraler Fluss markierter Primvierlinge (#Energiedoku)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Plot geschrieben: {path}")


def print_report(rows: list[ChiralRow], liouville_at: int | None) -> None:
    print("=== EABC-Spektrum (4×4) ===")
    for name, err in verify_eabc_spectrum().items():
        print(f"  {name}: max Abweichung {err:.2e}")

    print("\n=== Chiraler Fluss S(N) = ABCE - CEAB ===")
    print(f"{'N':>12} {'Q':>6} {'ABCE':>6} {'CEAB':>6} {'S':>5} {'|S|/Q':>10} {'z':>6}")
    for row in rows:
        print(
            f"{row.limit:12d} {row.total:6d} {row.abce:6d} {row.ceab:6d} "
            f"{row.flux:5d} {row.abs_flux_ratio:10.6f} {row.z_score:6.2f}"
        )

    if liouville_at is not None:
        lio = liouville_orthogonality(liouville_at)
        print(f"\n=== Liouville-Orthogonalität bei N={liouville_at}, Q={lio['Q']} ===")
        print(f"  Σ χ_H       = {lio['sum_chi']}")
        print(f"  Σ λ(p)χ_H   = {lio['sum_lam_chi']}")
        print(f"  Σ χ_H / Q   = {lio['sum_chi_over_Q']}")
        print(f"  Σ λχ / Q    = {lio['sum_lam_chi_over_Q']}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Chiraler Fluss S(N) für Primvierlinge")
    parser.add_argument(
        "--limits",
        nargs="+",
        type=float,
        default=[float(x) for x in default_limits()],
        help="Siebgrenzen (z. B. 1e5 1e6 1e7)",
    )
    parser.add_argument(
        "--liouville-at",
        type=int,
        default=10**7,
        help="N für Liouville-Orthogonalität (0 = aus)",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=ROOT / "energiedoku_chiraler_fluss.csv",
        help="CSV-Ausgabedatei",
    )
    parser.add_argument("--plot", action="store_true", help="PNG-Plot erzeugen")
    parser.add_argument(
        "--plot-path",
        type=Path,
        default=ROOT / "energiedoku_chiraler_fluss.png",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    limits = sorted({int(x) for x in args.limits})
    rows = analyze_limits(limits)
    write_csv(rows, args.csv)
    liouville_limit = args.liouville_at if args.liouville_at > 0 else None
    print_report(rows, liouville_limit)
    print(f"\nCSV geschrieben: {args.csv}")
    if args.plot:
        maybe_plot(rows, args.plot_path)


if __name__ == "__main__":
    main()
