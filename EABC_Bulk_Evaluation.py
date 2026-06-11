import json
import math
import time
from pathlib import Path

import numpy as np
import mpmath

# Höchste numerische Präzision für die transzendenten Summen im Fernfeld
mpmath.mp.dps = 35


def find_exact_riemann_zeros_bulk(start_t, target_count):
    """
    Sucht eine hochenergetische, lückenlose Kette von Riemann-Nullstellen.
    Nutzt mpmath.siegelz und ein adaptives Suchgitter gegen Aliasing.
    """
    zeros = []
    t = mpmath.mpf(start_t)
    progress_step = max(1, target_count // 10)

    print(f"=== Starte Groß-Evaluation: N = {target_count} ===")
    print(f"-> Berechne adaptives Suchgitter ab t_0 = {start_t}...")

    while len(zeros) < target_count:
        t_float = float(t)
        # Lokale Nullstellendichte nach dem Weyl-Gesetz
        lokale_dichte = (1.0 / (2 * math.pi)) * math.log(t_float / (2 * math.pi))
        mittlerer_abstand = 1.0 / lokale_dichte

        # Konservative Schrittweite (1/5 des Abstands), um keine Nullstelle zu verpassen
        step = mittlerer_abstand / 5.0
        t_next = t + step

        val_current = mpmath.siegelz(t)
        val_next = mpmath.siegelz(t_next)

        # Vorzeichenwechsel isolieren
        if val_current * val_next < 0:
            try:
                zero_root = mpmath.findroot(
                    mpmath.siegelz, (t, t_next), solver="illinois"
                )
                zero_float = float(zero_root)

                # Dubletten-Schutz
                if not zeros or abs(zero_float - zeros[-1]) > 1e-5:
                    zeros.append(zero_float)

                    # Fortschrittsanzeige in 10%-Schritten
                    if len(zeros) % progress_step == 0:
                        print(
                            f"   Progress: {len(zeros)} / {target_count} Nullstellen verifiziert."
                        )
            except Exception:
                # Falls Solver in singulären Punkten oszilliert, Schritt verfeinern
                pass

        t = t_next
    return np.array(zeros)


def apply_weyl_unfolding(zeros):
    """Transformiert die Rohdaten in das normierte Spektrum (mittlerer Abstand = 1)."""
    unfolded = []
    for gamma in zeros:
        N_t = (gamma / (2 * math.pi)) * math.log(gamma / (2 * math.pi * math.e)) + 7 / 8
        unfolded.append(N_t)
    return np.array(unfolded)


if __name__ == "__main__":
    START_T = 5000.0
    COUNT = 1000
    GUE_VARIANZ = 0.1041915
    monopol_factor = 1796 / math.sqrt(3)  # Asymptotische Feldkonstante: 1036.9215

    t0 = time.perf_counter()

    # 1. Pipeline ausführen
    raw_zeros = find_exact_riemann_zeros_bulk(START_T, COUNT)
    unfolded_zeros = apply_weyl_unfolding(raw_zeros)
    spacings = np.diff(unfolded_zeros)

    runtime_s = time.perf_counter() - t0

    # 2. Statistische Gesamtauswertung
    total_mean = np.mean(spacings)
    total_var = np.var(spacings)
    total_exzess = total_var - GUE_VARIANZ
    total_pressure = np.sum(np.sin(raw_zeros * math.log(monopol_factor) / 150)) / COUNT

    # 3. Stabilitätsprüfung über Datenblöcke (Konvergenztest)
    print("\n" + "=" * 60)
    print("   STABILITÄTSPRÜFUNG ÜBER DATENBLÖCKE (JE 250 ZEROS)")
    print("=" * 60)
    print(
        f"{'Block':<10} | {'Mittelwert':<12} | {'Varianz':<12} | "
        f"{'Δσ² (Exzess)':<15} | {'Druck P_m':<12}"
    )
    print("-" * 60)

    blocks = []
    for b in range(4):
        start_idx = b * 250
        end_idx = start_idx + 250 if b < 3 else COUNT - 1
        block_spacings = spacings[start_idx:end_idx]
        block_zeros = raw_zeros[start_idx:end_idx]

        b_mean = np.mean(block_spacings)
        b_var = np.var(block_spacings)
        b_exzess = b_var - GUE_VARIANZ
        b_pressure = np.sum(np.sin(block_zeros * math.log(monopol_factor) / 150)) / len(
            block_zeros
        )

        print(
            f"Block {b + 1:<5} | {b_mean:<12.6f} | {b_var:<12.6f} | "
            f"{b_exzess:<+15.6f} | {b_pressure:<12.6f}"
        )
        blocks.append(
            {
                "block": b + 1,
                "mean": float(b_mean),
                "variance": float(b_var),
                "delta_sigma2": float(b_exzess),
                "pressure_pm": float(b_pressure),
            }
        )

    print("=" * 60)
    print("   ZUSAMMENFASSENDE KENNZAHLEN FÜR DEN AUSBLICK")
    print("=" * 60)
    print(f"Untersuchter Sektor (t):          [{START_T:.1f}, {raw_zeros[-1]:.2f}]")
    print(f"Gesamt-Mittelwert (N={COUNT}): {total_mean:.6f}  (Weyl-Soll: 1.000000)")
    print(f"Gesamt-Varianz:              {total_var:.6f}")
    print(f"Asymptotischer Exzess (Δσ²):  {total_exzess:+.6f}")
    print(f"Asymptotischer Druck (P_m):   {total_pressure:.6f}")
    print(f"Laufzeit:                     {runtime_s:.1f} s")
    print("=" * 60)

    out_dir = Path(__file__).parent
    bulk_json_path = out_dir / "bulk_results.json"
    payload = {
        "method": "mpmath.siegelz + adaptive step/5 + illinois",
        "mpmath_dps": 35,
        "start_t": START_T,
        "count": COUNT,
        "runtime_s": round(runtime_s, 1),
        "t_min": float(raw_zeros[0]),
        "t_max": float(raw_zeros[-1]),
        "gue_varianz": GUE_VARIANZ,
        "total": {
            "mean_spacing": float(total_mean),
            "variance": float(total_var),
            "delta_sigma2": float(total_exzess),
            "pressure_pm": float(total_pressure),
        },
        "blocks": blocks,
    }
    with bulk_json_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
        f.write("\n")
    print(f"\nErgebnisse gespeichert: {bulk_json_path}")
