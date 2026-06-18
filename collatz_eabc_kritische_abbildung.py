#!/usr/bin/env python3
"""
EABC kritische Abbildung: Geschwindigkeitsmodell s_v(x) entlang der kritischen Linie.

Theorie: collatz_eabc_kritische_abbildung.md (Modellabbildung, kein Physikanspruch)
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


def traverse_circuit(
    edges: tuple[tuple[str, str, int], ...],
    edge_velocities: dict[str, float],
    x0: float = SOURCE_X,
    length_model: str = "canonical",
) -> list[dict[str, Any]]:
    """
    Traversiere gerichtete Kanten; kumulative x und s_v an jedem Knoten.

    edge_velocities keys: v_A, v_B, v_C, v_E (Geschwindigkeit am Startknoten der Kante).
    length_model: 'canonical' → mod-12-Lücken; 'unit' → alle ℓ=1.
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
        ell = 1 if length_model == "unit" else ell_canonical
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
        "boxed": {
            "mapping": "x ↦ 1/2 + i v (x - 1/2)",
            "holonomy_link": "ABCEA (+1) vs CEABC (-1) auf C4 mit Lücken (2,4,2,4)",
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
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
