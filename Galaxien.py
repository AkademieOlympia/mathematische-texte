"""
Galaxienspin-Simulation mit leichtem Orientierungs-Bias.
Läuft mit Standard-Python (ohne SageMath). Eine voll oktave Algebraische
Dynamik würde SageMath (Octonions) benötigen; die hier zählt nur gezupfte ±1.
"""

import argparse
import random


def check_octonionic_bias(n_galaxies: int, bias_factor: float = 0.01) -> float:
    """
    Simuliert Galaxienspins als ±1 mit leichter Anisotropie (bias_factor).

    bias_factor: Verschiebt die Schwelle 0.5, sodass die mittlere
    Richtung leicht von 0 wegdriftet (grobe Proxy für „abweichende
    Orientierung“; keine echte S^7-Verteilung).
    """
    if n_galaxies <= 0:
        return float("nan")
    results = []
    for _ in range(n_galaxies):
        # Leichte Bevorzugung +1, wenn Zufall + Bias oberhalb 0.5
        spin = 1 if (random.random() + bias_factor) > 0.5 else -1
        results.append(spin)
    return sum(results) / n_galaxies


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mittlerer Spin-Bias (oktonisch motiviert, numerisch simuliert)."
    )
    parser.add_argument(
        "-n",
        type=int,
        default=200_000,
        help="Anzahl simulierter Galaxien (Standard: 200000)",
    )
    parser.add_argument(
        "--bias",
        type=float,
        default=0.01,
        help="Bias-Faktor (Standard: 0.01)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Zufalls-Seed (Reproduzierbarkeit)",
    )
    args = parser.parse_args()
    if args.seed is not None:
        random.seed(args.seed)

    m = check_octonionic_bias(args.n, args.bias)
    print(f"n = {args.n}, bias_factor = {args.bias}", flush=True)
    print(f"Mittlerer Bias (Durchschnitt ±1): {m}", flush=True)


# __main__ = direkter Start; "<run_path>" = z. B. runpy/Cursor/IDE-„Datei ausführen“
if __name__ in ("__main__", "<run_path>"):
    main()
