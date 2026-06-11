import math

import mpmath
import numpy as np
from scipy.optimize import root_scalar


def riemann_z(t):
    """Vollständige Riemann-Siegel-Z-Funktion inkl. R(t)-Korrektur."""
    return float(mpmath.siegelz(t))


def find_riemann_zeros_adaptive(start_t, count):
    """Sucht lückenlos nach Nullstellen mit adaptiver Schrittweite."""
    zeros = []
    t = start_t

    while len(zeros) < count:
        erwartete_dichte = (1.0 / (2 * math.pi)) * math.log(t / (2 * math.pi))
        if erwartete_dichte <= 0:
            t += 0.1
            continue
        mittlerer_abstand = 1.0 / erwartete_dichte
        step = mittlerer_abstand / 4.0

        t_next = t + step
        try:
            val_current = riemann_z(t)
            val_next = riemann_z(t_next)

            if val_current * val_next < 0:
                sol = root_scalar(riemann_z, bracket=[t, t_next], method="brentq")
                if not zeros or abs(sol.root - zeros[-1]) > 1e-5:
                    zeros.append(sol.root)
        except (ValueError, ZeroDivisionError):
            pass
        t = t_next
    return zeros


def unfold_zeros(zeros):
    unfolded = []
    for gamma in zeros:
        N_t = (gamma / (2 * math.pi)) * math.log(gamma / (2 * math.pi * math.e)) + 7 / 8
        unfolded.append(N_t)
    return np.array(unfolded)


GUE_SIGMA2 = 0.10419


def compute_statistics(zeros):
    """Berechnet normierte Abstandsstatistiken und Monopol-Druck."""
    unfolded = unfold_zeros(zeros)
    spacings = np.diff(unfolded)
    monopol_factor = 1796 / math.sqrt(3)
    spectral_pressure = np.sum(np.sin(np.array(zeros) * math.log(monopol_factor) / 150))
    return {
        "n_zeros": len(zeros),
        "n_spacings": len(spacings),
        "t_min": zeros[0],
        "t_max": zeros[-1],
        "mean_spacing": float(np.mean(spacings)),
        "variance": float(np.var(spacings)),
        "std_spacing": float(np.std(spacings)),
        "delta_sigma2": float(np.var(spacings) - GUE_SIGMA2),
        "spectral_pressure": float(spectral_pressure),
    }


if __name__ == "__main__":
    import argparse
    import json
    import time
    from pathlib import Path

    parser = argparse.ArgumentParser(description="Monopol-Test (mpmath.siegelz)")
    parser.add_argument(
        "--preprint",
        action="store_true",
        help="Preprint-Validierung: START_T=5000, COUNT=500, Konvergenztabelle",
    )
    parser.add_argument("--start-t", type=float, default=5000.0)
    parser.add_argument("--count", type=int, default=500)
    args = parser.parse_args()

    if args.preprint or (args.start_t == 5000.0 and args.count == 500):
        START_T = 5000.0
        COUNT = 500
    else:
        START_T = args.start_t
        COUNT = args.count

    print("--- Korrigierter arithmetischer Monopol-Test ---")

    print(
        f"Suche {COUNT} lückenlose Riemann-Nullstellen ab t = {START_T} "
        "mit adaptiver Schrittweite..."
    )
    t0 = time.perf_counter()
    zeros = find_riemann_zeros_adaptive(START_T, COUNT)
    runtime_s = time.perf_counter() - t0

    stats = compute_statistics(zeros)

    print("\n--- Analytische Ergebnisse (Preprint-Validierung, N=500) ---")
    print(f"Gefundene Nullstellen:       {stats['n_zeros']}")
    print(f"Anzahl Abstände:             {stats['n_spacings']}")
    print(f"t-Bereich:                   [{stats['t_min']:.4f}, {stats['t_max']:.4f}]")
    print(f"Mittlerer normierter Abstand: {stats['mean_spacing']:.6f} (Soll: ~1.0)")
    print(f"Varianz σ²:                  {stats['variance']:.6f} (GUE-Soll: {GUE_SIGMA2})")
    print(f"Δσ² = σ² − σ²_GUE:           {stats['delta_sigma2']:.6f}")
    print(f"Std.-Abw. der Abstände:      {stats['std_spacing']:.6f}")
    print(f"Arithmetischer Monopol-Druck: {stats['spectral_pressure']:.6f}")
    print(f"Laufzeit:                    {runtime_s:.1f} s")

    print("\n--- Konvergenz der Varianz (Teilstichproben) ---")
    convergence = {}
    for n in (50, 100, 200, 500):
        if n <= len(zeros):
            sub = compute_statistics(zeros[:n])
            convergence[n] = sub
            print(
                f"N={n:3d}: σ²={sub['variance']:.6f}, "
                f"Δσ²={sub['delta_sigma2']:.6f}, "
                f"mean={sub['mean_spacing']:.6f}"
            )

    out_dir = Path(__file__).parent
    results_path = out_dir / "monopol_preprint_results.txt"
    json_path = out_dir / "monopol_preprint_results.json"

    payload = {
        "method": "mpmath.siegelz + adaptive step + brentq",
        "start_t": START_T,
        "count": COUNT,
        "runtime_s": round(runtime_s, 1),
        "gue_sigma2": GUE_SIGMA2,
        "n500": stats,
        "convergence": {
            str(n): {
                "variance": sub["variance"],
                "delta_sigma2": sub["delta_sigma2"],
                "mean_spacing": sub["mean_spacing"],
            }
            for n, sub in convergence.items()
        },
    }
    with json_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
        f.write("\n")

    with results_path.open("w", encoding="utf-8") as f:
        f.write("% Monopol-Test Ergebnisse (mpmath.siegelz, START_T=5000)\n")
        f.write(f"% Laufzeit: {runtime_s:.1f} s\n")
        f.write(f"% Unfolding: N(t) = (t/2π) log(t/2πe) + 7/8\n\n")
        f.write(f"\\newcommand{{\\monopolN}}{{{stats['n_zeros']}}}\n")
        f.write(f"\\newcommand{{\\monopolMeanSpacing}}{{{stats['mean_spacing']:.6f}}}\n")
        f.write(f"\\newcommand{{\\monopolVariance}}{{{stats['variance']:.6f}}}\n")
        f.write(f"\\newcommand{{\\monopolDeltaSigma}}{{{stats['delta_sigma2']:.6f}}}\n")
        f.write(f"\\newcommand{{\\monopolStdSpacing}}{{{stats['std_spacing']:.6f}}}\n")
        f.write(f"\\newcommand{{\\monopolPressure}}{{{stats['spectral_pressure']:.6f}}}\n")
        f.write(f"\\newcommand{{\\monopolTmin}}{{{stats['t_min']:.4f}}}\n")
        f.write(f"\\newcommand{{\\monopolTmax}}{{{stats['t_max']:.4f}}}\n")
        f.write(f"\\newcommand{{\\monopolGUE}}{{{GUE_SIGMA2}}}\n")
        f.write("\n% Konvergenz:\n")
        for n, sub in convergence.items():
            f.write(
                f"% N={n}: sigma2={sub['variance']:.6f}, "
                f"delta={sub['delta_sigma2']:.6f}, mean={sub['mean_spacing']:.6f}\n"
            )
    print(f"\nErgebnisse gespeichert: {results_path}")
    print(f"JSON gespeichert:         {json_path}")
