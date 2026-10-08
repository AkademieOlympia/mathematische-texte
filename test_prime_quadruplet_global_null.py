#!/usr/bin/env python3
"""
Globaler Nullhypothesentest für Abstände konsekutiver Primzahlvierlinge.

Primzahlvierlinge:
    Q(p) = (p, p+2, p+6, p+8)

Globale Nullhypothese H0:
    Nach Entfernung der lokalen Primzahldichte sind die Abstände konsekutiver
    Vierlinge unabhängig von
      1. der Zahlenskala,
      2. der Größe der Vorgängerlücke,
      3. einem Überschuss oder Defizit kleiner Lücken,
      4. der relativen EABC-Orientierung benachbarter Vierlinge.

Teilstatistiken:
    T1 trend_slope
       negativ = relative Annäherung
       positiv = relative Entfernung

    T2 small_gap_excess
       positiv = mehr kleine Lücken / Clusterung
       negativ = Kleinlückendefizit / Abstoßung

    T3 lag1_spearman
       positiv = Persistenz kleiner bzw. großer Lücken
       negativ = alternierende Kompensation

    T4 orientation_gap_effect
       positiv = gleichphasige Vierlinge liegen weiter auseinander
       negativ = gleichphasige Vierlinge liegen näher zusammen

Globaler Test:
    Max-T-Bootstrap über die vier standardisierten Teilstatistiken.

Siebtreues Nullmodell:
    - reguläres Kandidatengitter p ≡ 11 (mod 30),
    - Entfernung aller Kandidaten, die durch kleine Primzahlen lokal
      ausgeschlossen sind,
    - Stichprobengewicht proportional zu 1/log(p)^4,
    - Konditionierung auf die beobachtete Zahl der Vierlinge,
    - standardmäßig Erhaltung der beobachteten ABCE/CEAB-Anzahlen.

Der Sondervierling (5,7,11,13) wird ausgeschlossen, da p=5 nicht zum
regulären Kandidatengitter p ≡ 11 (mod 30) gehört.

Benötigte Pakete:
    numpy pandas scipy matplotlib

Installation:
    python3 -m pip install numpy pandas scipy matplotlib

Beispiel:
    python3 test_prime_quadruplet_global_null.py \
        --limit 20000000 \
        --bootstrap 5000 \
        --sieve-prime-max 97 \
        --seed 20260719 \
        --output-dir results_quadruplets

Schneller Probelauf:
    python3 test_prime_quadruplet_global_null.py \
        --limit 1000000 \
        --bootstrap 200 \
        --no-plots
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

OFFSETS = np.array([0, 2, 6, 8], dtype=np.int64)
STAT_NAMES = (
    "trend_slope",
    "small_gap_excess",
    "lag1_spearman",
    "orientation_gap_effect",
)


@dataclass(frozen=True)
class TestStatistics:
    trend_slope: float
    small_gap_excess: float
    lag1_spearman: float
    orientation_gap_effect: float

    def as_array(self) -> np.ndarray:
        return np.array(
            [
                self.trend_slope,
                self.small_gap_excess,
                self.lag1_spearman,
                self.orientation_gap_effect,
            ],
            dtype=float,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Globaler Nulltest für konsekutive Primzahlvierlinge."
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=20_000_000,
        help="Maximale reguläre Startprimzahl p (Standard: 20.000.000).",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=1_000,
        help="Anzahl Bootstrap-Surrogate (Standard: 1000).",
    )
    parser.add_argument(
        "--sieve-prime-max",
        type=int,
        default=97,
        help="Kleine Primzahlen bis zu dieser Grenze im Nullmodell (Standard: 97).",
    )
    parser.add_argument(
        "--small-gap-threshold",
        type=float,
        default=0.5,
        help="Schwelle für kleine entfaltete Lücken (Standard: 0.5).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260719,
        help="Zufallsseed.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("quadruplet_null_test"),
        help="Ausgabeordner.",
    )
    parser.add_argument(
        "--free-orientation-counts",
        action="store_true",
        help=(
            "ABCE/CEAB-Anzahlen im Bootstrap nicht fixieren. "
            "Standardmäßig werden die beobachteten Anzahlen erhalten."
        ),
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Keine PNG-Diagramme erzeugen.",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.limit < 100_000:
        raise ValueError("--limit muss mindestens 100000 sein.")
    if args.bootstrap < 100:
        raise ValueError("--bootstrap muss mindestens 100 sein.")
    if args.sieve_prime_max < 7:
        raise ValueError("--sieve-prime-max muss mindestens 7 sein.")
    if not (0.0 < args.small_gap_threshold < 5.0):
        raise ValueError("--small-gap-threshold muss zwischen 0 und 5 liegen.")


def prime_sieve(limit: int) -> np.ndarray:
    """Eratosthenes-Sieb als Bool-Array."""
    is_prime = np.ones(limit + 1, dtype=np.bool_)
    is_prime[:2] = False
    if limit >= 4:
        is_prime[4::2] = False
    for p in range(3, math.isqrt(limit) + 1, 2):
        if is_prime[p]:
            is_prime[p * p :: 2 * p] = False
    return is_prime


def small_primes_up_to(limit: int) -> np.ndarray:
    return np.flatnonzero(prime_sieve(limit))


def find_quadruplet_starts(limit: int, is_prime: np.ndarray) -> np.ndarray:
    """Exakte reguläre Starts p ≡ 11 mod 30; der Sonderfall p=5 bleibt draußen."""
    candidates = np.arange(11, limit + 1, 30, dtype=np.int64)
    mask = np.ones(candidates.size, dtype=bool)
    for h in OFFSETS:
        mask &= is_prime[candidates + h]
    return candidates[mask]


def orientation(starts: np.ndarray) -> np.ndarray:
    """
    +1 = ABCE für p ≡ 41 mod 60
    -1 = CEAB für p ≡ 11 mod 60
    """
    residues = starts % 60
    valid = (residues == 11) | (residues == 41)
    if not np.all(valid):
        bad = np.unique(residues[~valid])
        raise ValueError(f"Unerwartete Startrestklassen modulo 60: {bad.tolist()}")
    return np.where(residues == 41, 1, -1).astype(np.int8)


def build_sieve_faithful_candidates(
    limit: int,
    sieve_prime_max: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Reguläre Kandidaten p ≡ 11 mod 30, die für keine kleine Primzahl q > 5
    bereits lokal ausgeschlossen sind.

    Nach dem kleinen Sieb wird als glatte Restintensität 1/log(p)^4 verwendet.
    """
    candidates = np.arange(11, limit + 1, 30, dtype=np.int64)
    survive = np.ones(candidates.size, dtype=bool)

    for q in small_primes_up_to(sieve_prime_max):
        if q in (2, 3, 5):
            continue
        forbidden = np.zeros(candidates.size, dtype=bool)
        for h in OFFSETS:
            values = candidates + h
            # Ein Wert, der exakt der Siebprimzahl q entspricht, ist prim und
            # darf nicht als lokales Hindernis entfernt werden. Ausgeschlossen
            # werden nur echte Vielfache 2q, 3q, ...
            forbidden |= (((values % q) == 0) & (values != q))
        survive &= ~forbidden

    valid = candidates[survive]
    if valid.size == 0:
        raise RuntimeError("Keine Kandidaten überleben das lokale Sieb.")

    weights = 1.0 / np.log(valid.astype(float)) ** 4
    if np.any(~np.isfinite(weights)) or np.any(weights <= 0):
        raise RuntimeError("Ungültige Kandidatengewichte.")
    return valid, weights


def event_indices(valid_candidates: np.ndarray, starts: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(valid_candidates, starts)
    if np.any(idx >= valid_candidates.size):
        raise RuntimeError("Mindestens ein Vierling liegt außerhalb des Kandidatenraums.")

    recovered = valid_candidates[idx]
    if not np.array_equal(recovered, starts):
        mismatch = recovered != starts
        missing = starts[mismatch]
        raise RuntimeError(
            "Das kleine Primzahlsieb hat echte Vierlinge entfernt. "
            f"Beispiele: {missing[:10].tolist()}"
        )
    return idx


def normalized_gaps(
    valid_candidates: np.ndarray,
    weights: np.ndarray,
    starts: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Diskretes Unfolding:
        hazard(p_j) ∝ 1/log(p_j)^4 auf dem siebtreuen Kandidatenraum.

    Rückgabe:
        gaps     entfaltete Nachbarlücken mit empirischem Mittelwert 1
        midpoints geometrische Mittelpunkte der konsekutiven Startzahlen
    """
    if starts.size < 20:
        raise ValueError("Zu wenige Ereignisse für einen stabilen Lückentest.")

    scaled_hazard = weights * (starts.size / weights.sum())
    cumulative_hazard = np.cumsum(scaled_hazard)
    idx = event_indices(valid_candidates, starts)
    tau = cumulative_hazard[idx]

    gaps = np.diff(tau)
    if np.any(gaps <= 0):
        raise RuntimeError("Nichtpositive entfaltete Lücke gefunden.")

    gaps = gaps / gaps.mean()
    midpoints = np.sqrt(starts[:-1].astype(float) * starts[1:].astype(float))
    return gaps, midpoints


def linear_slope(x: np.ndarray, y: np.ndarray) -> float:
    """OLS-Steigung pro einer Standardabweichung von x."""
    x_sd = x.std(ddof=1)
    if not np.isfinite(x_sd) or x_sd <= 0:
        return float("nan")
    z = (x - x.mean()) / x_sd
    denominator = np.dot(z, z)
    if denominator <= 0:
        return float("nan")
    return float(np.dot(z, y - y.mean()) / denominator)


def compute_statistics(
    starts: np.ndarray,
    valid_candidates: np.ndarray,
    weights: np.ndarray,
    small_gap_threshold: float,
) -> TestStatistics:
    gaps, midpoints = normalized_gaps(valid_candidates, weights, starts)
    log_gaps = np.log(gaps)

    # T1: Trend der normalisierten Lücken entlang der Zahlenskala.
    trend = linear_slope(np.log(midpoints), log_gaps)

    # T2: Kleinlückenüberschuss relativ zur Exp(1)-Referenz.
    # Die Signifikanz wird nicht aus Exp(1), sondern aus dem Bootstrap bestimmt.
    exp_reference = 1.0 - math.exp(-small_gap_threshold)
    small_excess = float(np.mean(gaps <= small_gap_threshold) - exp_reference)

    # T3: serielle Abhängigkeit konsekutiver Lücken.
    rho_result = spearmanr(log_gaps[:-1], log_gaps[1:])
    lag1 = float(rho_result.statistic)

    # T4: Abstand bei gleicher vs. verschiedener EABC-Orientierung.
    ori = orientation(starts)
    same_orientation = ori[:-1] == ori[1:]
    if not np.any(same_orientation) or np.all(same_orientation):
        orient_effect = float("nan")
    else:
        orient_effect = float(
            log_gaps[same_orientation].mean()
            - log_gaps[~same_orientation].mean()
        )

    return TestStatistics(
        trend_slope=trend,
        small_gap_excess=small_excess,
        lag1_spearman=lag1,
        orientation_gap_effect=orient_effect,
    )


def weighted_sample_without_replacement(
    rng: np.random.Generator,
    values: np.ndarray,
    weights: np.ndarray,
    size: int,
) -> np.ndarray:
    if size < 0 or size > values.size:
        raise ValueError(
            f"Stichprobengröße {size} ist für {values.size} Kandidaten unzulässig."
        )
    probabilities = weights / weights.sum()
    chosen_indices = rng.choice(
        values.size,
        size=size,
        replace=False,
        p=probabilities,
    )
    return values[chosen_indices]


def draw_surrogate(
    rng: np.random.Generator,
    valid_candidates: np.ndarray,
    weights: np.ndarray,
    n_total: int,
    n_abce: int,
    preserve_orientation_counts: bool,
) -> np.ndarray:
    if preserve_orientation_counts:
        mask_abce = (valid_candidates % 60) == 41
        mask_ceab = (valid_candidates % 60) == 11
        n_ceab = n_total - n_abce

        sample_abce = weighted_sample_without_replacement(
            rng,
            valid_candidates[mask_abce],
            weights[mask_abce],
            n_abce,
        )
        sample_ceab = weighted_sample_without_replacement(
            rng,
            valid_candidates[mask_ceab],
            weights[mask_ceab],
            n_ceab,
        )
        sample = np.concatenate((sample_abce, sample_ceab))
    else:
        sample = weighted_sample_without_replacement(
            rng,
            valid_candidates,
            weights,
            n_total,
        )

    sample.sort()
    return sample


def bootstrap_null(
    rng: np.random.Generator,
    repetitions: int,
    observed_starts: np.ndarray,
    valid_candidates: np.ndarray,
    weights: np.ndarray,
    small_gap_threshold: float,
    preserve_orientation_counts: bool,
) -> np.ndarray:
    n_total = observed_starts.size
    n_abce = int(np.sum(orientation(observed_starts) == 1))
    output = np.empty((repetitions, len(STAT_NAMES)), dtype=float)

    report_every = max(1, repetitions // 20)

    for b in range(repetitions):
        surrogate = draw_surrogate(
            rng=rng,
            valid_candidates=valid_candidates,
            weights=weights,
            n_total=n_total,
            n_abce=n_abce,
            preserve_orientation_counts=preserve_orientation_counts,
        )
        stats = compute_statistics(
            starts=surrogate,
            valid_candidates=valid_candidates,
            weights=weights,
            small_gap_threshold=small_gap_threshold,
        )
        output[b] = stats.as_array()

        if (b + 1) % report_every == 0 or b + 1 == repetitions:
            print(f"\rBootstrap: {b + 1}/{repetitions}", end="", flush=True)

    print()
    if np.any(~np.isfinite(output)):
        raise RuntimeError(
            "Nichtendliche Bootstrap-Statistik gefunden. "
            "Erhöhe --limit oder überprüfe die Orientierungsklassen."
        )
    return output


def empirical_two_sided_p(observed: float, null_values: np.ndarray) -> float:
    """
    Zweiseitiger empirischer p-Wert um das simulierte Nullzentrum.
    """
    center = float(null_values.mean())
    obs_distance = abs(observed - center)
    null_distances = np.abs(null_values - center)
    return float(
        (1 + np.sum(null_distances >= obs_distance))
        / (null_values.size + 1)
    )


def summarize_results(
    observed: np.ndarray,
    null_matrix: np.ndarray,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    means = null_matrix.mean(axis=0)
    sds = null_matrix.std(axis=0, ddof=1)

    if np.any(~np.isfinite(sds)) or np.any(sds <= 0):
        raise RuntimeError("Mindestens eine Nullverteilung besitzt Varianz null.")

    z_observed = (observed - means) / sds
    z_null = (null_matrix - means) / sds

    max_t_observed = float(np.max(np.abs(z_observed)))
    max_t_null = np.max(np.abs(z_null), axis=1)
    global_p = float(
        (1 + np.sum(max_t_null >= max_t_observed))
        / (null_matrix.shape[0] + 1)
    )

    interpretations = {
        "trend_slope":
            "negativ: relative Annäherung; positiv: relative Entfernung",
        "small_gap_excess":
            "positiv: Clusterung; negativ: Kleinlückendefizit/Abstoßung",
        "lag1_spearman":
            "positiv: Persistenz; negativ: alternierende Kompensation",
        "orientation_gap_effect":
            "positiv: gleichphasige Abstoßung; negativ: gleichphasige Nähe",
    }

    rows = []
    for j, name in enumerate(STAT_NAMES):
        rows.append(
            {
                "statistik": name,
                "beobachtet": observed[j],
                "null_mittel": means[j],
                "null_sd": sds[j],
                "z_wert": z_observed[j],
                "p_empirisch_zweiseitig": empirical_two_sided_p(
                    observed[j], null_matrix[:, j]
                ),
                "interpretation_vorzeichen": interpretations[name],
            }
        )

    summary = pd.DataFrame(rows)
    global_result = {
        "max_abs_z_observed": max_t_observed,
        "global_maxT_p_value": global_p,
    }
    return summary, global_result


def make_plots(
    output_dir: Path,
    starts: np.ndarray,
    valid_candidates: np.ndarray,
    weights: np.ndarray,
    observed: np.ndarray,
    null_matrix: np.ndarray,
    global_result: Dict[str, float],
) -> None:
    import matplotlib.pyplot as plt

    gaps, midpoints = normalized_gaps(valid_candidates, weights, starts)

    plt.figure(figsize=(9, 5))
    plt.scatter(midpoints, gaps, s=10, alpha=0.55)
    plt.axhline(1.0, linewidth=1)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Geometrischer Mittelpunkt konsekutiver Startprimzahlen")
    plt.ylabel("Entfaltete Lücke s_n")
    plt.title("Entfaltete Abstände konsekutiver Primzahlvierlinge")
    plt.tight_layout()
    plt.savefig(output_dir / "unfolded_gaps.png", dpi=180)
    plt.close()

    window = max(25, min(200, gaps.size // 10))
    rolling = pd.Series(gaps).rolling(window, center=True).median()
    plt.figure(figsize=(9, 5))
    plt.plot(midpoints, rolling)
    plt.axhline(1.0, linewidth=1)
    plt.xscale("log")
    plt.xlabel("Geometrischer Mittelpunkt")
    plt.ylabel(f"Rollender Median der s_n, Fenster={window}")
    plt.title("Lokale Entwicklung der entfalteten Vierlingslücken")
    plt.tight_layout()
    plt.savefig(output_dir / "rolling_median_gaps.png", dpi=180)
    plt.close()

    means = null_matrix.mean(axis=0)
    sds = null_matrix.std(axis=0, ddof=1)
    max_t_null = np.max(np.abs((null_matrix - means) / sds), axis=1)

    plt.figure(figsize=(8, 5))
    plt.hist(max_t_null, bins=40)
    plt.axvline(global_result["max_abs_z_observed"], linewidth=2)
    plt.xlabel("max_j |Z_j|")
    plt.ylabel("Bootstrap-Häufigkeit")
    plt.title("Globale Max-T-Nullverteilung")
    plt.tight_layout()
    plt.savefig(output_dir / "global_maxT_null.png", dpi=180)
    plt.close()

    for j, name in enumerate(STAT_NAMES):
        plt.figure(figsize=(8, 5))
        plt.hist(null_matrix[:, j], bins=40)
        plt.axvline(observed[j], linewidth=2)
        plt.xlabel(name)
        plt.ylabel("Bootstrap-Häufigkeit")
        plt.title(f"Nullverteilung: {name}")
        plt.tight_layout()
        plt.savefig(output_dir / f"null_{name}.png", dpi=180)
        plt.close()


def main() -> int:
    args = parse_args()

    try:
        validate_args(args)
    except ValueError as exc:
        print(f"Parameterfehler: {exc}", file=sys.stderr)
        return 2

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Erzeuge Primzahlsieb bis {args.limit + 8:,} ...")
    is_prime = prime_sieve(args.limit + 8)

    starts = find_quadruplet_starts(args.limit, is_prime)
    if starts.size < 30:
        print(
            "Zu wenige Primzahlvierlinge für einen stabilen Test. "
            "Erhöhe --limit.",
            file=sys.stderr,
        )
        return 3

    ori = orientation(starts)
    n_abce = int(np.sum(ori == 1))
    n_ceab = int(np.sum(ori == -1))

    print(
        f"Gefunden: {starts.size} reguläre Primzahlvierlinge "
        f"(ABCE={n_abce}, CEAB={n_ceab})."
    )

    print(
        f"Baue siebtreuen Kandidatenraum bis Primzahl "
        f"{args.sieve_prime_max} ..."
    )
    valid_candidates, weights = build_sieve_faithful_candidates(
        limit=args.limit,
        sieve_prime_max=args.sieve_prime_max,
    )
    event_indices(valid_candidates, starts)

    total_regular_candidates = ((args.limit - 11) // 30) + 1
    print(
        f"Siebtreue Kandidaten: {valid_candidates.size:,} "
        f"von {total_regular_candidates:,} regulären Kandidaten."
    )

    observed_stats = compute_statistics(
        starts=starts,
        valid_candidates=valid_candidates,
        weights=weights,
        small_gap_threshold=args.small_gap_threshold,
    )
    observed = observed_stats.as_array()

    print("Beobachtete Teilstatistiken:")
    for name, value in zip(STAT_NAMES, observed):
        print(f"  {name:28s} {value:+.8f}")

    rng = np.random.default_rng(args.seed)
    null_matrix = bootstrap_null(
        rng=rng,
        repetitions=args.bootstrap,
        observed_starts=starts,
        valid_candidates=valid_candidates,
        weights=weights,
        small_gap_threshold=args.small_gap_threshold,
        preserve_orientation_counts=not args.free_orientation_counts,
    )

    summary, global_result = summarize_results(
        observed=observed,
        null_matrix=null_matrix,
    )

    starts_df = pd.DataFrame(
        {
            "p": starts,
            "orientation": np.where(ori == 1, "ABCE", "CEAB"),
        }
    )
    starts_df.to_csv(args.output_dir / "quadruplet_starts.csv", index=False)
    summary.to_csv(args.output_dir / "component_tests.csv", index=False)
    pd.DataFrame(
        null_matrix,
        columns=STAT_NAMES,
    ).to_csv(
        args.output_dir / "bootstrap_null_statistics.csv",
        index=False,
    )

    metadata = {
        "limit": args.limit,
        "bootstrap_repetitions": args.bootstrap,
        "sieve_prime_max": args.sieve_prime_max,
        "small_gap_threshold": args.small_gap_threshold,
        "seed": args.seed,
        "preserve_orientation_counts": not args.free_orientation_counts,
        "regular_quadruplets": int(starts.size),
        "ABCE": n_abce,
        "CEAB": n_ceab,
        "valid_null_candidates": int(valid_candidates.size),
        **global_result,
    }
    with open(
        args.output_dir / "global_result.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)

    if not args.no_plots:
        make_plots(
            output_dir=args.output_dir,
            starts=starts,
            valid_candidates=valid_candidates,
            weights=weights,
            observed=observed,
            null_matrix=null_matrix,
            global_result=global_result,
        )

    print("\nKomponententests:")
    print(
        summary.to_string(
            index=False,
            float_format=lambda x: f"{x:.8g}",
        )
    )

    print(
        "\nGlobaler Max-T-Test:\n"
        f"  max |Z| = {global_result['max_abs_z_observed']:.6f}\n"
        f"  globaler p-Wert = "
        f"{global_result['global_maxT_p_value']:.6f}"
    )

    if global_result["global_maxT_p_value"] < 0.05:
        print(
            "\nEntscheidung bei alpha=0.05: "
            "Globale Nullhypothese verwerfen."
        )
    else:
        print(
            "\nEntscheidung bei alpha=0.05: "
            "Globale Nullhypothese nicht verwerfen."
        )

    print(f"\nErgebnisse gespeichert in: {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
