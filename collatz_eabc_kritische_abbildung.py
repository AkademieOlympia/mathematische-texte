#!/usr/bin/env python3
"""
EABC kritische Abbildung: Geschwindigkeitsmodell s_v(x) entlang der kritischen Linie.

Theorie: collatz_eabc_kritische_abbildung.md (Modellabbildung, kein Physikanspruch)
Epistemik: collatz_eabc_epistemik_physik.md (Wegfunktion T, Nicht-SRT, Physik-vs.-EABC)
Verknüpfung: collatz_eabc_zirkulationshypothese.md (Holonomie ±1, ABCEA/CEABC)

  s_v(x) = 1/2 + i v (x - 1/2)
  x_{n,v} = 1/2 + γ_n / v
  EABC C4-Schaltkreis: Kantenlängen (2,4,2,4) mod 12

Ausführung:
    python3 collatz_eabc_kritische_abbildung.py
    python3 collatz_eabc_kritische_abbildung.py --v 2 --n-zeros 10
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_kritische_abbildung.json"
THEORY = "collatz_eabc_kritische_abbildung.md"
THEORY_ZIRKULATION = "collatz_eabc_zirkulationshypothese.md"
THEORY_ZYKLUS = "collatz_eabc_zyklus_holonomie.md"
THEORY_GENERALANGRIFF = "collatz_generalangriff_2026.md"

SOURCE_X = 0.5
CRITICAL_REAL = 0.5

# Kanonisches Lückenmuster mod 12 auf A→B→C→E→A
CANONICAL_GAP_PATTERN = (2, 4, 2, 4)
UNIT_GAP_PATTERN = (1, 1, 1, 1)

# Vorwärtsorientierung A→B→C→E→A
EDGE_AB = ("A", "B", 2)
EDGE_BC = ("B", "C", 4)
EDGE_CE = ("C", "E", 2)
EDGE_EA = ("E", "A", 4)

ABCEA_EDGES: tuple[tuple[str, str, int], ...] = (EDGE_AB, EDGE_BC, EDGE_CE, EDGE_EA)
CEABC_EDGES: tuple[tuple[str, str, int], ...] = (EDGE_CE, EDGE_EA, EDGE_AB, EDGE_BC)

# Gerichtete Kanten in zyklischer EAABC-Reihenfolge (E→A, A→B, B→C, C→E)
EAABC_EDGE_KEYS: tuple[str, ...] = ("EA", "AB", "BC", "CE")
ABCEA_EDGE_KEYS: tuple[str, ...] = ("AB", "BC", "CE", "EA")

# Erste 20 Imaginärteile nicht-trivialer ζ-Nullstellen (mpmath-kompatibel, Fallback)
ZETA_GAMMA_HARDCODED: tuple[float, ...] = (
    14.134725141734693,
    21.022039638771554,
    25.010857580145688,
    30.424876125859513,
    32.93506158773919,
    37.58617815882567,
    40.918719012147495,
    43.32707328091499,
    48.00515088116716,
    49.7738324776723,
    52.97032147771446,
    56.446247697063394,
    59.34704400260235,
    60.83177852460981,
    65.1125440480816,
    67.07981052949417,
    69.54640171117398,
    72.0671576744819,
    75.70469069908398,
    77.1448400688748,
)

GAMMA_1_APPROX = ZETA_GAMMA_HARDCODED[0]


def s_v(x: float, v: float) -> complex:
    """s_v(x) = 1/2 + i v (x - 1/2)."""
    return complex(CRITICAL_REAL, v * (x - SOURCE_X))


def gamma_v(x: float, v: float) -> float:
    """γ_v(x) = v (x - 1/2) — Imaginärteil auf der kritischen Linie."""
    return v * (x - SOURCE_X)


def x_from_gamma(gamma: float, v: float) -> float:
    """Inverse: x = 1/2 + γ / v."""
    if v == 0:
        raise ValueError("v must be positive")
    return SOURCE_X + gamma / v


def x_n_v(n: int, v: float, gammas: tuple[float, ...] | list[float] | None = None) -> float:
    """x_{n,v} = 1/2 + γ_n / v für 1-basiertes n."""
    if n < 1:
        raise ValueError("n must be >= 1")
    table = gammas if gammas is not None else zeta_imaginary_parts(len(ZETA_GAMMA_HARDCODED))
    if n > len(table):
        raise IndexError(f"n={n} exceeds available zeros ({len(table)})")
    return x_from_gamma(table[n - 1], v)


def zeta_imaginary_parts(n: int) -> list[float]:
    """Erste n Imaginärteile γ_n der ζ-Nullstellen (mpmath falls verfügbar)."""
    if n <= 0:
        return []
    try:
        from mpmath import zetazero  # type: ignore[import-untyped]

        return [float(zetazero(k).imag) for k in range(1, n + 1)]
    except ImportError:
        if n <= len(ZETA_GAMMA_HARDCODED):
            return list(ZETA_GAMMA_HARDCODED[:n])
        raise RuntimeError(
            f"mpmath not available and n={n} exceeds hardcoded table ({len(ZETA_GAMMA_HARDCODED)})"
        ) from None


def zero_mapping_table(
    n: int,
    velocities: tuple[float, ...] | list[float] = (1.0, 2.0, 10.0),
    gammas: list[float] | None = None,
) -> list[dict[str, Any]]:
    """Tabelle x_{n,v} für die ersten n Nullstellen und gegebene v."""
    if gammas is None:
        gammas = zeta_imaginary_parts(n)
    rows: list[dict[str, Any]] = []
    for idx, gamma in enumerate(gammas[:n], start=1):
        row: dict[str, Any] = {
            "n": idx,
            "gamma_n": gamma,
            "rho_n": f"1/2 + i*{gamma}",
        }
        for v in velocities:
            xv = x_from_gamma(gamma, v)
            row[f"x_n_v_{v:g}"] = xv
            row[f"s_v_{v:g}"] = {"re": CRITICAL_REAL, "im": gamma}
        rows.append(row)
    return rows


def example_gamma_1(velocities: tuple[float, ...] = (1.0, 2.0, 10.0)) -> dict[str, float]:
    """Beispiel γ_1 ≈ 14.134725 für v ∈ {1,2,10}."""
    return {f"v={v:g}": x_from_gamma(GAMMA_1_APPROX, v) for v in velocities}


def edge_velocities_uniform(v: float) -> dict[str, float]:
    """Ein-Parameter-Modell: v_E = v_A = v_B = v_C = v."""
    return {"v_E": v, "v_A": v, "v_B": v, "v_C": v}


def gaps_eaabc_to_abcea(gaps: tuple[int, ...] | list[int]) -> tuple[int, ...]:
    """(ℓ_EA, ℓ_AB, ℓ_BC, ℓ_CE) → Traversierungsreihenfolge ABCEA (ℓ_AB, ℓ_BC, ℓ_CE, ℓ_EA)."""
    g = tuple(gaps)
    if len(g) != 4:
        raise ValueError(f"expected 4 gaps, got {len(g)}")
    ell_ea, ell_ab, ell_bc, ell_ce = g
    return (ell_ab, ell_bc, ell_ce, ell_ea)


def gaps_abcea_to_eaabc(gaps: tuple[int, ...] | list[int]) -> tuple[int, ...]:
    """(ℓ_AB, ℓ_BC, ℓ_CE, ℓ_EA) → zyklische EAABC-Liste (ℓ_EA, ℓ_AB, ℓ_BC, ℓ_CE)."""
    g = tuple(gaps)
    if len(g) != 4:
        raise ValueError(f"expected 4 gaps, got {len(g)}")
    ell_ab, ell_bc, ell_ce, ell_ea = g
    return (ell_ea, ell_ab, ell_bc, ell_ce)


def normalize_gaps(
    gaps: tuple[int, ...] | list[int],
    gap_order: str = "ABCEA",
) -> tuple[int, ...]:
    """Normalisiere Lücken auf ABCEA-Traversierungsreihenfolge."""
    g = tuple(gaps)
    if gap_order == "ABCEA":
        return g
    if gap_order == "EAABC":
        return gaps_eaabc_to_abcea(g)
    raise ValueError(f"unknown gap_order: {gap_order}")


def edges_for_orientation(
    orientation: str,
    gaps_abcea: tuple[int, ...],
) -> tuple[tuple[str, str, int], ...]:
    """Gerichtete Kanten mit Längen ℓ_j für ABCEA oder CEABC."""
    length_by_edge = {
        f"{src}{dst}": ell for (src, dst, _), ell in zip(ABCEA_EDGES, gaps_abcea)
    }
    template = ABCEA_EDGES if orientation == "ABCEA" else CEABC_EDGES
    return tuple((src, dst, length_by_edge[f"{src}{dst}"]) for src, dst, _ in template)


def edge_velocities_from_gaps(
    gaps: tuple[int, ...] | list[int],
    gamma_ref: float,
    *,
    gap_order: str = "ABCEA",
) -> dict[str, Any]:
    """
    Kantengeschwindigkeiten aus Primlücken mod 12: v_j = γ_ref / ℓ_j.

    Liefert Schlüssel v_EA, v_AB, v_BC, v_CE sowie v_E, v_A, v_B, v_C
    (Geschwindigkeit am Startknoten jeder Kante für traverse_circuit).

    gap_order:
      'ABCEA' — (ℓ_AB, ℓ_BC, ℓ_CE, ℓ_EA), kanonisch (2,4,2,4)
      'EAABC' — (ℓ_EA, ℓ_AB, ℓ_BC, ℓ_CE), zyklisch (4,2,4,2)
    """
    if gamma_ref <= 0:
        raise ValueError("gamma_ref must be positive")
    gaps_ab = normalize_gaps(gaps, gap_order)
    result: dict[str, Any] = {"gamma_ref": gamma_ref}
    for (src, dst, _), ell in zip(ABCEA_EDGES, gaps_ab):
        if ell <= 0:
            raise ValueError(f"edge length must be positive, got {ell}")
        v_edge = gamma_ref / ell
        result[f"v_{src}{dst}"] = v_edge
        result[f"v_{src}"] = v_edge
    result["gaps_abcea"] = list(gaps_ab)
    return result


def path_time_T(
    gaps: tuple[int, ...] | list[int],
    gamma_ref: float,
    *,
    gap_order: str = "ABCEA",
) -> float:
    """
    Wegfunktion T = Σ_j ℓ_j / v_j mit v_j = γ_ref / ℓ_j  ⇒  T = Σ_j ℓ_j² / γ_ref.

    Rein euklidisch — keine Zeitdilatation (collatz_eabc_epistemik_physik.md §1).
    Verschiedene Pfadgeometrien auf denselben Knoten (gerade vs. Halbkreis) ändern
    die komplexe Weglänge via compare_path_times, nicht diese Kanten-Wegfunktion.

    Siehe collatz_eabc_kritische_abbildung.md §7.
    """
    if gamma_ref <= 0:
        raise ValueError("gamma_ref must be positive")
    gaps_norm = normalize_gaps(gaps, gap_order)
    for ell in gaps_norm:
        if ell <= 0:
            raise ValueError(f"edge length must be positive, got {ell}")
    return sum(ell * ell for ell in gaps_norm) / gamma_ref


def holonomy_sensor_trajectory(
    orientation: str = "ABCEA",
    gaps: tuple[int, ...] | list[int] | None = None,
    *,
    gamma_ref: float | None = None,
    v_base: float | None = None,
    gamma_n_index: int = 1,
    gap_order: str = "ABCEA",
) -> dict[str, Any]:
    """
    EABC-Holonomie-Sensor: Ray-Mapping mit aus Lücken abgeleiteten v_j.

    Pro Kante: Δγ_j = v_j · ℓ_j = γ_ref (konstante Höheninkremente).
    gamma_ref: Referenzskala (typisch γ_n); falls None: v_base oder γ_1.
    """
    gaps_in = tuple(gaps) if gaps is not None else CANONICAL_GAP_PATTERN
    if gamma_ref is None:
        if v_base is not None:
            gamma_ref = v_base
        else:
            gamma_ref = zeta_imaginary_parts(gamma_n_index)[gamma_n_index - 1]
    gaps_ab = normalize_gaps(gaps_in, gap_order)
    edge_v = edge_velocities_from_gaps(gaps_ab, gamma_ref, gap_order="ABCEA")
    edges = edges_for_orientation(orientation, gaps_ab)
    segments = traverse_circuit(edges, edge_v, x0=SOURCE_X, length_model="explicit")
    total_length = sum(s.get("length", 0) for s in segments if "length" in s)
    total_gamma = segments[-1]["gamma_cumulative"] if segments else 0.0
    return {
        "sensor": "EABC_holonomy",
        "epistemic": "Modellabbildung — kein Physikanspruch",
        "orientation": orientation,
        "holonomy_sign": holonomy_sign(orientation),
        "gamma_ref": gamma_ref,
        "gap_order_input": gap_order,
        "gaps_abcea": list(gaps_ab),
        "gaps_eaabc": list(gaps_abcea_to_eaabc(gaps_ab)),
        "edge_velocities": {k: edge_v[k] for k in ("v_EA", "v_AB", "v_BC", "v_CE")},
        "edge_velocities_vertex": {k: edge_v[k] for k in ("v_E", "v_A", "v_B", "v_C")},
        "segments": segments,
        "total_length": total_length,
        "total_gamma": total_gamma,
        "delta_gamma_per_edge": gamma_ref,
        "endpoint": segments[-1] if segments else None,
    }


def compare_holonomy_sensor_trajectories(
    gaps: tuple[int, ...] | list[int] | None = None,
    *,
    gamma_ref: float | None = None,
    v_base: float | None = None,
    gamma_n_index: int = 1,
    gap_order: str = "ABCEA",
) -> dict[str, Any]:
    """Vergleich ABCEA (+1) vs. CEABC (-1) mit gleichen abgeleiteten Kantengeschwindigkeiten."""
    kwargs: dict[str, Any] = {
        "gaps": gaps,
        "gamma_ref": gamma_ref,
        "v_base": v_base,
        "gamma_n_index": gamma_n_index,
        "gap_order": gap_order,
    }
    abcea = holonomy_sensor_trajectory("ABCEA", **kwargs)
    ceabc = holonomy_sensor_trajectory("CEABC", **kwargs)
    gref = abcea["gamma_ref"]
    return {
        "sensor": "EABC_holonomy_compare",
        "gamma_ref": gref,
        "gaps_abcea": abcea["gaps_abcea"],
        "gaps_eaabc": abcea["gaps_eaabc"],
        "edge_velocities": abcea["edge_velocities"],
        "ABCEA": abcea,
        "CEABC": ceabc,
        "holonomy_contrast": {
            "sign_ABCEA": abcea["holonomy_sign"],
            "sign_CEABC": ceabc["holonomy_sign"],
            "same_total_length": abcea["total_length"] == ceabc["total_length"],
            "same_total_gamma": math.isclose(
                abcea["total_gamma"], ceabc["total_gamma"], rel_tol=0, abs_tol=1e-12
            ),
            "same_edge_velocities": abcea["edge_velocities"] == ceabc["edge_velocities"],
        },
    }


def prime_window_gap_samples(
    max_p: int = 10_000,
    limit: int = 5,
) -> list[dict[str, Any]]:
    """
    Beispiel-Fenster aus der Primfolge: ABCEA/CEABC mit mod-12-Lücken und abgeleiteten v_j.
    """
    from collatz_eabc_holonomie_fehlerterm import EABC_RESIDUES, gap_pattern_mod12
    from collatz_eabc_transition_graph import ABCEA_WORD, CEABC_WORD, prime_eabc_sequence
    from eabc_from_lean import EClass

    seq = prime_eabc_sequence(max_p)
    classes = [row["class"] for row in seq]
    samples: list[dict[str, Any]] = []
    for word in (ABCEA_WORD, CEABC_WORD):
        count = 0
        orientation = "ABCEA" if word == ABCEA_WORD else "CEABC"
        for i in range(len(classes) - 4):
            window = "".join(classes[i : i + 5])
            if window != word:
                continue
            residues = tuple(EABC_RESIDUES[EClass(c)] for c in window)
            gaps = gap_pattern_mod12(residues)
            traj = holonomy_sensor_trajectory(
                orientation,
                gaps=gaps,
                gamma_ref=GAMMA_1_APPROX,
            )
            samples.append(
                {
                    "word": word,
                    "orientation": orientation,
                    "prime_index_start": i,
                    "primes": [seq[j]["p"] for j in range(i, i + 5)],
                    "residues_mod12": list(residues),
                    "gaps_abcea": list(gaps),
                    "edge_velocities": traj["edge_velocities"],
                    "total_gamma": traj["total_gamma"],
                }
            )
            count += 1
            if count >= limit:
                break
    return samples


def traverse_circuit(
    edges: tuple[tuple[str, str, int], ...],
    edge_velocities: dict[str, float],
    x0: float = SOURCE_X,
    length_model: str = "canonical",
) -> list[dict[str, Any]]:
    """
    Traversiere gerichtete Kanten; kumulative x und s_v an jedem Knoten.

    edge_velocities keys: v_A, v_B, v_C, v_E (Geschwindigkeit am Startknoten der Kante).
    length_model: 'canonical' → mod-12-Lücken aus Kantendefinition;
                  'explicit' → ℓ aus Kantentupel;
                  'unit' → alle ℓ=1.
    """
    segments: list[dict[str, Any]] = []
    x = x0
    gamma_cum = 0.0
    segments.append(
        {
            "step": 0,
            "vertex": edges[0][0] if edges else "P",
            "x": x,
            "gamma_cumulative": gamma_cum,
            "s_v": {"re": CRITICAL_REAL, "im": gamma_cum},
        }
    )
    for step, (src, dst, ell_canonical) in enumerate(edges, start=1):
        if length_model == "unit":
            ell = 1
        elif length_model == "explicit":
            ell = ell_canonical
        else:
            ell = ell_canonical
        v_key = f"v_{src}"
        v_seg = edge_velocities.get(v_key, edge_velocities.get("v", 1.0))
        x_prev = x
        x = x + ell
        delta_gamma = v_seg * ell
        gamma_cum += delta_gamma
        s_point = s_v(x, v_seg)
        segments.append(
            {
                "step": step,
                "edge": f"{src}->{dst}",
                "src": src,
                "dst": dst,
                "length": ell,
                "velocity": v_seg,
                "x_prev": x_prev,
                "x": x,
                "delta_gamma": delta_gamma,
                "gamma_cumulative": gamma_cum,
                "s_v": {"re": s_point.real, "im": s_point.imag},
            }
        )
    return segments


def holonomy_sign(orientation: str) -> int:
    """ABCEA → +1, CEABC → -1 (collatz_eabc_zirkulationshypothese.md)."""
    if orientation == "ABCEA":
        return 1
    if orientation == "CEABC":
        return -1
    raise ValueError(f"unknown orientation: {orientation}")


def chord_length(z1: complex, z2: complex) -> float:
    """Euklidische Sehnenlänge |z2 - z1|."""
    return abs(z2 - z1)


def semicircle_arc_length(z1: complex, z2: complex) -> float:
    """Halbkreis in der oberen Halbebene (Durchmesser = Sehne z1–z2)."""
    return math.pi * abs(z2 - z1) / 2


def linear_round_trip_time(trajectory_points: list[complex] | tuple[complex, ...]) -> float:
    """
    „Soldaten“-Photon: geradlinige Segmente P → s_v(x_1) → … → zurück zu P.

    Einheitsgeschwindigkeit |dz/dt| = 1 ⇒ Rückkehrzeit = Gesamtlänge.
  """
    points = list(trajectory_points)
    if len(points) < 2:
        return 0.0
    total = sum(chord_length(points[i], points[i + 1]) for i in range(len(points) - 1))
    total += chord_length(points[-1], points[0])
    return total


def semicircle_chain_time(trajectory_points: list[complex] | tuple[complex, ...]) -> float:
    """
    Verkettete Halbkreisbögen in der oberen Halbebene zwischen aufeinanderfolgenden Besuchspunkten.
    """
    points = list(trajectory_points)
    if len(points) < 2:
        return 0.0
    total = sum(
        semicircle_arc_length(points[i], points[i + 1]) for i in range(len(points) - 1)
    )
    total += semicircle_arc_length(points[-1], points[0])
    return total


def sensor_trajectory_points(
    orientation: str = "ABCEA",
    gaps: tuple[int, ...] | list[int] | None = None,
    *,
    gamma_ref: float | None = None,
    v_base: float | None = None,
    gamma_n_index: int = 1,
    gap_order: str = "ABCEA",
) -> list[complex]:
    """Besuchspunkte des Holonomie-Sensors als komplexe Zahlen (inkl. Start P)."""
    traj = holonomy_sensor_trajectory(
        orientation,
        gaps=gaps,
        gamma_ref=gamma_ref,
        v_base=v_base,
        gamma_n_index=gamma_n_index,
        gap_order=gap_order,
    )
    return [
        complex(seg["s_v"]["re"], seg["s_v"]["im"]) for seg in traj["segments"]
    ]


def compare_path_times(
    orientation: str = "ABCEA",
    gaps: tuple[int, ...] | list[int] | None = None,
    *,
    gamma_ref: float | None = None,
    v_base: float | None = None,
    gamma_n_index: int = 1,
    gap_order: str = "ABCEA",
) -> dict[str, Any]:
    """
    Vergleich Halbkreis-Kette vs. gerader Polygonzug auf dem Holonomie-Sensor.

    Kein Einstein-Zwillingsparadoxon — nur euklidische Weglängenvergleich (Modellabbildung).
    """
    if gamma_ref is None:
        if v_base is not None:
            gamma_ref = v_base
        else:
            gamma_ref = zeta_imaginary_parts(gamma_n_index)[gamma_n_index - 1]
    gaps_in = tuple(gaps) if gaps is not None else CANONICAL_GAP_PATTERN
    points = sensor_trajectory_points(
        orientation,
        gaps=gaps_in,
        gamma_ref=gamma_ref,
        gap_order=gap_order,
    )
    t_linear = linear_round_trip_time(points)
    t_semi = semicircle_chain_time(points)
    ratio = t_semi / t_linear if t_linear > 0 else float("nan")

    other = "CEABC" if orientation == "ABCEA" else "ABCEA"
    points_other = sensor_trajectory_points(
        other,
        gaps=gaps_in,
        gamma_ref=gamma_ref,
        gap_order=gap_order,
    )
    t_linear_other = linear_round_trip_time(points_other)
    t_semi_other = semicircle_chain_time(points_other)

    return {
        "epistemic": "Geometrischer Weglängenvergleich — kein SRT-Zwillingsparadoxon",
        "orientation": orientation,
        "holonomy_sign": holonomy_sign(orientation),
        "gamma_ref": gamma_ref,
        "gaps_abcea": list(normalize_gaps(gaps_in, gap_order)),
        "trajectory_points": [{"re": z.real, "im": z.imag} for z in points],
        "T_linear": t_linear,
        "T_semicircle": t_semi,
        "ratio_semi_over_linear": ratio,
        "other_orientation": other,
        "other_T_linear": t_linear_other,
        "other_T_semicircle": t_semi_other,
        "same_linear_time_both_orientations": math.isclose(
            t_linear, t_linear_other, rel_tol=0, abs_tol=1e-12
        ),
        "same_semicircle_time_both_orientations": math.isclose(
            t_semi, t_semi_other, rel_tol=0, abs_tol=1e-12
        ),
    }


def eabc_circuit_report(
    v: float = 1.0,
    orientation: str = "ABCEA",
    length_model: str = "canonical",
    edge_v: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Schaltkreisbericht für ABCEA oder CEABC."""
    if edge_v is None:
        edge_v = edge_velocities_uniform(v)
    edges = ABCEA_EDGES if orientation == "ABCEA" else CEABC_EDGES
    segments = traverse_circuit(edges, edge_v, x0=SOURCE_X, length_model=length_model)
    total_length = sum(s.get("length", 0) for s in segments if "length" in s)
    total_gamma = segments[-1]["gamma_cumulative"] if segments else 0.0
    return {
        "orientation": orientation,
        "holonomy_sign": holonomy_sign(orientation),
        "length_model": length_model,
        "gap_pattern": list(UNIT_GAP_PATTERN if length_model == "unit" else CANONICAL_GAP_PATTERN),
        "edge_velocities": edge_v,
        "source_P": {"x": SOURCE_X, "s_v": {"re": CRITICAL_REAL, "im": 0.0}},
        "segments": segments,
        "total_x_span": segments[-1]["x"] - SOURCE_X if segments else 0.0,
        "total_length": total_length,
        "total_gamma": total_gamma,
        "endpoint": segments[-1] if segments else None,
    }


def dual_circuit_report(
    v: float = 1.0,
    length_model: str = "canonical",
) -> dict[str, Any]:
    """ABCEA und CEABC bei gleichem v."""
    edge_v = edge_velocities_uniform(v)
    abcea = eabc_circuit_report(v=v, orientation="ABCEA", length_model=length_model, edge_v=edge_v)
    ceabc = eabc_circuit_report(v=v, orientation="CEABC", length_model=length_model, edge_v=edge_v)
    return {
        "v": v,
        "length_model": length_model,
        "ABCEA": abcea,
        "CEABC": ceabc,
        "holonomy_contrast": {
            "sign_ABCEA": abcea["holonomy_sign"],
            "sign_CEABC": ceabc["holonomy_sign"],
            "same_total_length": abcea["total_length"] == ceabc["total_length"],
            "same_total_gamma_magnitude": math.isclose(
                abs(abcea["total_gamma"]), abs(ceabc["total_gamma"]), rel_tol=0, abs_tol=1e-12
            ),
        },
    }


def run(
    n_zeros: int = 5,
    velocities: tuple[float, ...] = (1.0, 2.0, 10.0),
    v_circuit: float = 1.0,
    output: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    gammas = zeta_imaginary_parts(max(n_zeros, 1))
    report: dict[str, Any] = {
        "meta": {
            "module": "collatz_eabc_kritische_abbildung.py",
            "theory": THEORY,
            "theory_zirkulation": THEORY_ZIRKULATION,
            "theory_zyklus": THEORY_ZYKLUS,
            "theory_generalangriff": THEORY_GENERALANGRIFF,
            "epistemic": "Modellabbildung — kein Physikanspruch",
        },
        "formulas": {
            "source_P": f"({SOURCE_X}, 0)",
            "gamma_v": "v * (x - 1/2)",
            "s_v": "1/2 + i * v * (x - 1/2)",
            "inverse": "x = 1/2 + gamma / v",
            "x_n_v": "1/2 + gamma_n / v",
        },
        "example_gamma_1": example_gamma_1(velocities),
        "zeta_zero_mappings": zero_mapping_table(n_zeros, velocities, gammas),
        "eabc_circuits": dual_circuit_report(v=v_circuit, length_model="canonical"),
        "eabc_circuits_unit": dual_circuit_report(v=v_circuit, length_model="unit"),
        "holonomy_sensor": compare_holonomy_sensor_trajectories(
            gaps=CANONICAL_GAP_PATTERN, gamma_ref=GAMMA_1_APPROX
        ),
        "path_time_compare": compare_path_times(
            orientation="ABCEA", gaps=CANONICAL_GAP_PATTERN, gamma_ref=GAMMA_1_APPROX
        ),
        "boxed": {
            "mapping": "x ↦ 1/2 + i v (x - 1/2)",
            "holonomy_link": "ABCEA (+1) vs CEABC (-1) auf C4 mit Lücken (2,4,2,4)",
            "holonomy_sensor": "v_j = γ_ref / ℓ_j — Ray-Mapping als EABC-Holonomie-Sensor",
            "path_compare": "T_semi/T_linear = π/2 bei vertikalen Segmenten — kein Zwillingsparadoxon",
            "zirkulation": "D_E = N_plus - N_minus (collatz_eabc_zirkulationshypothese.md)",
        },
    }
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC kritische Abbildung s_v(x)")
    parser.add_argument("--n-zeros", type=int, default=5)
    parser.add_argument("--v", type=float, default=1.0, help="Schaltkreis-Geschwindigkeit")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(n_zeros=args.n_zeros, v_circuit=args.v, output=args.output)
    circ = report["eabc_circuits"]["ABCEA"]
    print("=== EABC kritische Abbildung ===")
    print(f"s_v(x) = 1/2 + i v (x - 1/2), Quelle P = ({SOURCE_X}, 0)")
    print()
    print("Beispiel γ_1:")
    for k, val in report["example_gamma_1"].items():
        print(f"  {k}: x ≈ {val:.7f}")
    print()
    print(f"ABCEA @ v={args.v} (kanonische Lücken {list(CANONICAL_GAP_PATTERN)}):")
    for seg in circ["segments"]:
        if seg["step"] == 0:
            print(f"  P: x={seg['x']}, s={seg['s_v']}")
        else:
            print(
                f"  {seg['edge']} ℓ={seg['length']} v={seg['velocity']}: "
                f"x={seg['x']}, γ_cum={seg['gamma_cumulative']}, s={seg['s_v']}"
            )
    print()
    print(f"Holonomie: ABCEA={report['eabc_circuits']['ABCEA']['holonomy_sign']}, "
          f"CEABC={report['eabc_circuits']['CEABC']['holonomy_sign']}")
    sensor = report["holonomy_sensor"]
    print()
    print(f"Holonomie-Sensor (γ_ref=γ_1, Lücken {sensor['gaps_abcea']}):")
    for ek, val in sensor["edge_velocities"].items():
        print(f"  {ek} = {val:.7f}")
    print(f"  Σγ = {sensor['ABCEA']['total_gamma']:.7f} (4·γ_ref)")
    paths = report["path_time_compare"]
    print()
    print(f"Weglängenvergleich ABCEA (γ_ref=γ_1, Lücken {paths['gaps_abcea']}):")
    print(f"  T_linear     = {paths['T_linear']:.7f}")
    print(f"  T_semicircle = {paths['T_semicircle']:.7f}")
    print(f"  T_semi/T_lin = {paths['ratio_semi_over_linear']:.7f}")
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
