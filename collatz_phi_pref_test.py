#!/usr/bin/env python3
"""Φ_pref-Datentest: Primvierlings-EABC-Wörter vs. Zufallsnullmodell.

Referenz: collatz_kepler_gedankenexperiment.tex (Abschnitt Φ_pref),
CollatzEabc.PrefProjection.lean.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import statistics
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal

from eabc_from_lean import EClass, class_of, is_prime_quadruplet, q

Z0 = complex(-1.0, -1.0 / 12.0)
CHIRALITY_WORDS = frozenset({"ABCE", "CEAB"})
IOTA_E: dict[EClass, int] = {
    EClass.E: 0,
    EClass.A: 1,
    EClass.B: 2,
    EClass.C: 3,
}
LETTERS = "EABC"
Source = Literal["quadruplet", "random"]


@dataclass(frozen=True, slots=True)
class PrefPoint:
    word: str
    source: Source
    chirality: str | None
    z_re: float
    z_im: float
    t: float
    theta: int
    theta4: int
    rho: float
    phase_re: float
    phase_im: float


@dataclass(frozen=True, slots=True)
class TubeStats:
    label: str
    count: int
    z_re_mean: float
    z_im_mean: float
    z_re_std: float
    z_im_std: float
    t_mean: float


@dataclass(frozen=True, slots=True)
class SeparationSummary:
    limit: int
    quadruplet_count: int
    abce_count: int
    ceab_count: int
    random_count: int
    tube_abce: TubeStats
    tube_ceab: TubeStats
    tube_random: TubeStats
    inter_tube_distance: float
    random_mean_dist_to_nearest_tube: float
    random_phase_shell_hits: dict[str, int]
    quadruplet_phase_shell_hits: dict[str, int]
    real_clusters_in_two_tubes: bool
    random_more_diffuse_than_real: bool
    verdict: str


def theta(w: str) -> int:
    return sum(IOTA_E[EClass(c)] for c in w)


def theta4(w: str) -> int:
    return theta(w) % 4


def phase_unit(k: int) -> complex:
    return (1.0, 1.0j, -1.0, -1.0j)[k % 4]


def unit_phase(w: str) -> complex:
    """Chirale Normalisierung für ABCE/CEAB; sonst phaseUnit ∘ Θ₄."""
    if w == "ABCE":
        return 1.0j
    if w == "CEAB":
        return -1.0j
    return phase_unit(theta4(w))


def rho(w: str) -> float:
    return 2.0 ** (-len(w))


def z_complex(w: str) -> complex:
    return Z0 + rho(w) * unit_phase(w)


def phi_pref(w: str) -> tuple[complex, float]:
    return z_complex(w), float(len(w))


def quadruplet_word(p: int) -> str:
    return "".join(class_of(n).value for n in q(p))


def chirality_word(p: int) -> str:
    residue = p % 12
    if residue == 5:
        return "ABCE"
    if residue == 11:
        return "CEAB"
    raise ValueError(f"Ungültiger Vierlingsstart {p} (mod 12 = {residue})")


def iter_quadruplet_starts(limit: int) -> list[int]:
    starts: list[int] = []
    for p in range(5, limit - 7, 2):
        if p % 12 not in (5, 11):
            continue
        if is_prime_quadruplet(p):
            starts.append(p)
    return starts


def random_eabc_word(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(LETTERS) for _ in range(length))


def _phase_shell_label(u: complex) -> str:
    for name, value in (
        ("+1", 1.0 + 0.0j),
        ("+i", 1.0j),
        ("-1", -1.0 + 0.0j),
        ("-i", -1.0j),
    ):
        if abs(u - value) < 1e-12:
            return name
    return "other"


def _point_record(
    word: str,
    source: Source,
    chirality: str | None,
) -> PrefPoint:
    z, t = phi_pref(word)
    u = unit_phase(word)
    return PrefPoint(
        word=word,
        source=source,
        chirality=chirality,
        z_re=z.real,
        z_im=z.imag,
        t=t,
        theta=theta(word),
        theta4=theta4(word),
        rho=rho(word),
        phase_re=u.real,
        phase_im=u.imag,
    )


def _tube_stats(label: str, points: list[PrefPoint]) -> TubeStats:
    if not points:
        return TubeStats(label, 0, 0.0, 0.0, 0.0, 0.0, 0.0)
    z_re = [p.z_re for p in points]
    z_im = [p.z_im for p in points]
    t_vals = [p.t for p in points]
    return TubeStats(
        label=label,
        count=len(points),
        z_re_mean=statistics.fmean(z_re),
        z_im_mean=statistics.fmean(z_im),
        z_re_std=statistics.pstdev(z_re) if len(z_re) > 1 else 0.0,
        z_im_std=statistics.pstdev(z_im) if len(z_im) > 1 else 0.0,
        t_mean=statistics.fmean(t_vals),
    )


def _dist(a: complex, b: complex) -> float:
    return abs(a - b)


def run_test(
    limit: int,
    random_count: int,
    *,
    seed: int = 20260617,
) -> tuple[list[PrefPoint], SeparationSummary]:
    starts = iter_quadruplet_starts(limit)
    rng = random.Random(seed)

    quadruplet_points: list[PrefPoint] = []
    abce_points: list[PrefPoint] = []
    ceab_points: list[PrefPoint] = []
    quad_shell_hits: dict[str, int] = {"+1": 0, "+i": 0, "-1": 0, "-i": 0, "other": 0}

    for p in starts:
        word = quadruplet_word(p)
        chir = chirality_word(p)
        pt = _point_record(word, "quadruplet", chir)
        quadruplet_points.append(pt)
        shell = _phase_shell_label(unit_phase(word))
        quad_shell_hits[shell] = quad_shell_hits.get(shell, 0) + 1
        if chir == "ABCE":
            abce_points.append(pt)
        else:
            ceab_points.append(pt)

    lengths = [len(p.word) for p in quadruplet_points] or [4]
    random_points: list[PrefPoint] = []
    rand_shell_hits: dict[str, int] = {"+1": 0, "+i": 0, "-1": 0, "-i": 0, "other": 0}
    for _ in range(random_count):
        length = rng.choice(lengths)
        word = random_eabc_word(length, rng)
        pt = _point_record(word, "random", None)
        random_points.append(pt)
        shell = _phase_shell_label(unit_phase(word))
        rand_shell_hits[shell] = rand_shell_hits.get(shell, 0) + 1

    tube_abce = _tube_stats("T_ABCE", abce_points)
    tube_ceab = _tube_stats("T_CEAB", ceab_points)
    tube_random = _tube_stats("random", random_points)

    z_abce = complex(tube_abce.z_re_mean, tube_abce.z_im_mean)
    z_ceab = complex(tube_ceab.z_re_mean, tube_ceab.z_im_mean)
    inter_tube = _dist(z_abce, z_ceab)

    tube_centers = {"ABCE": z_abce, "CEAB": z_ceab}
    random_dists = [
        min(_dist(complex(p.z_re, p.z_im), center) for center in tube_centers.values())
        for p in random_points
    ]
    random_mean_dist = statistics.fmean(random_dists) if random_dists else 0.0

    real_tight = (
        tube_abce.z_re_std == 0.0
        and tube_abce.z_im_std == 0.0
        and tube_ceab.z_re_std == 0.0
        and tube_ceab.z_im_std == 0.0
        and inter_tube > 1e-9
    )
    real_two_shells = (
        quad_shell_hits.get("+i", 0) == tube_abce.count
        and quad_shell_hits.get("-i", 0) == tube_ceab.count
        and quad_shell_hits.get("+1", 0) == 0
        and quad_shell_hits.get("-1", 0) == 0
    )
    random_spread = (
        tube_random.z_re_std > 0.0 or tube_random.z_im_std > 0.0 or len(set(
            (_phase_shell_label(complex(p.phase_re, p.phase_im)) for p in random_points)
        )) > 2
    )

    if real_tight and real_two_shells and random_spread:
        verdict = (
            "Ja: reale Primvierlinge liegen in zwei degenerierten Röhren "
            "T_ABCE/T_CEAB (+i/-i); Zufall nutzt alle vier Phasenschalen."
        )
    elif real_tight and real_two_shells:
        verdict = (
            "Teilweise: zwei chirale Röhren für reale Daten, aber Zufallsstreuung "
            "in diesem Lauf nicht deutlich breiter (Länge 4 → nur 4 Phasenschalen)."
        )
    else:
        verdict = "Nein: keine klare Zwei-Röhren-Trennung unter den gewählten Kriterien."

    summary = SeparationSummary(
        limit=limit,
        quadruplet_count=len(quadruplet_points),
        abce_count=len(abce_points),
        ceab_count=len(ceab_points),
        random_count=len(random_points),
        tube_abce=tube_abce,
        tube_ceab=tube_ceab,
        tube_random=tube_random,
        inter_tube_distance=inter_tube,
        random_mean_dist_to_nearest_tube=random_mean_dist,
        random_phase_shell_hits=rand_shell_hits,
        quadruplet_phase_shell_hits=quad_shell_hits,
        real_clusters_in_two_tubes=real_tight and real_two_shells,
        random_more_diffuse_than_real=random_spread,
        verdict=verdict,
    )
    return quadruplet_points + random_points, summary


def write_json(path: Path, points: list[PrefPoint], summary: SeparationSummary) -> None:
    payload = {
        "z0": {"re": Z0.real, "im": Z0.imag},
        "summary": asdict(summary),
        "points": [asdict(p) for p in points],
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def print_summary(summary: SeparationSummary) -> None:
    print("=== Φ_pref-Datentest (Primvierlinge vs. Zufall) ===")
    print(f"Limit: {summary.limit:,}  |  Vierlinge: {summary.quadruplet_count}")
    print(f"  ABCE: {summary.abce_count}  CEAB: {summary.ceab_count}")
    print(f"  Zufall: {summary.random_count} Wörter (gleiche Längenverteilung)")
    print()
    for tube in (summary.tube_abce, summary.tube_ceab, summary.tube_random):
        print(
            f"{tube.label:8s}  n={tube.count:5d}  "
            f"z̄=({tube.z_re_mean:+.6f}, {tube.z_im_mean:+.6f})  "
            f"σ_z=({tube.z_re_std:.2e}, {tube.z_im_std:.2e})  t̄={tube.t_mean:.2f}"
        )
    print()
    print(f"Abstand T_ABCE ↔ T_CEAB: {summary.inter_tube_distance:.6f}")
    print(f"Zufall: mittlere Distanz zur nächsten Röhre: {summary.random_mean_dist_to_nearest_tube:.6f}")
    print(f"Phasenschalen reale Daten: {summary.quadruplet_phase_shell_hits}")
    print(f"Phasenschalen Zufall:      {summary.random_phase_shell_hits}")
    print()
    print(f"Zwei chirale Röhren (real): {summary.real_clusters_in_two_tubes}")
    print(f"Zufall diffuser:             {summary.random_more_diffuse_than_real}")
    print(f"Urteil: {summary.verdict}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Φ_pref-Datentest für EABC-Wörter")
    parser.add_argument("--limit", type=int, default=1_000_000, help="Obergrenze für Primzahlsuche")
    parser.add_argument("--random-count", type=int, default=5000, help="Anzahl Zufallswörter")
    parser.add_argument("--seed", type=int, default=20260617)
    parser.add_argument(
        "--json",
        type=Path,
        default=Path("collatz_phi_pref_test.json"),
        help="JSON-Ausgabedatei",
    )
    args = parser.parse_args()

    points, summary = run_test(args.limit, args.random_count, seed=args.seed)
    write_json(args.json, points, summary)
    print_summary(summary)
    print(f"\nJSON: {args.json}")


if __name__ == "__main__":
    main()
