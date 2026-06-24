#!/usr/bin/env python3
"""
EABC-Invarianzprogramm: strikte Zählvektoren, Simplex-Anteile und Quadrupel-Signaturen.

Kanonsiche Definitionen: collatz_eabc_invarianzprogramm.md
Philosophischer Kontext: collatz_eabc_bernoulli_uebersetzung.md §22 (Querverweis).

Ausführung:
    python3 collatz_eabc_invarianzprogramm.py
    python3 collatz_eabc_invarianzprogramm.py --max-x 50000 --output collatz_eabc_invarianzprogramm.json
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from math import isqrt, sqrt
from pathlib import Path
from typing import Any

from eabc_from_lean import EClass, class_of, is_prime_quadruplet, q

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_invarianzprogramm.json"

EABC_ORDER = (EClass.E, EClass.A, EClass.B, EClass.C)
SIGMA4_ALPHABET = "EABC"

# Klein-V_4-Fouriermoden auf (E,A,B,C); Φ_0 aus Summen-Nebenbedingung ausgeschlossen.
PHI_0 = (1.0, 1.0, 1.0, 1.0)
PHI_1 = (1.0, -1.0, 1.0, -1.0)
PHI_2 = (1.0, 1.0, -1.0, -1.0)
PHI_3 = (1.0, -1.0, -1.0, 1.0)
V4_MODES = (PHI_1, PHI_2, PHI_3)
V4_MODE_NORM_SQ = 4.0


def _sieve_primes(limit: int) -> list[int]:
    if limit < 2:
        return []
    flags = [True] * (limit + 1)
    flags[0] = flags[1] = False
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            step = p
            start = p * p
            flags[start : limit + 1 : step] = [False] * (((limit - start) // step) + 1)
    return [i for i, ok in enumerate(flags) if ok]


def kappa(p: int) -> EClass | None:
    """κ: P_{>3} → {E,A,B,C}; für p ≤ 3 oder nicht-prim: None."""
    if p <= 3:
        return None
    return class_of(p)


@dataclass(frozen=True, slots=True)
class CountVector:
    """V(x) = (E(x), A(x), B(x), C(x)) — Zählvektor über p ≤ x in P_{>3}."""

    e: int
    a: int
    b: int
    c: int

    def as_tuple(self) -> tuple[int, int, int, int]:
        return (self.e, self.a, self.b, self.c)

    @property
    def total(self) -> int:
        return self.e + self.a + self.b + self.c

    def to_dict(self) -> dict[str, int]:
        return {"E": self.e, "A": self.a, "B": self.b, "C": self.c}


@dataclass(frozen=True, slots=True)
class SimplexPoint:
    """S(x) = (1/π(x))·V(x) ∈ Δ_3; π(x) hier = |{p ≤ x : p ∈ P_{>3}}|."""

    e: float
    a: float
    b: float
    c: float

    def as_tuple(self) -> tuple[float, float, float, float]:
        return (self.e, self.a, self.b, self.c)

    def to_dict(self) -> dict[str, float]:
        return {"E": self.e, "A": self.a, "B": self.b, "C": self.c}


def v_at_x(x: int, primes: list[int] | None = None) -> CountVector:
    """V(x): zählt Primzahlen p ≤ x, p > 3, nach EABC-Klasse."""
    if x < 5:
        return CountVector(0, 0, 0, 0)
    if primes is None:
        primes = _sieve_primes(x)
    counts = {cls: 0 for cls in EABC_ORDER}
    for p in primes:
        if p > x:
            break
        cls = kappa(p)
        if cls is not None:
            counts[cls] += 1
    return CountVector(
        e=counts[EClass.E],
        a=counts[EClass.A],
        b=counts[EClass.B],
        c=counts[EClass.C],
    )


def pi_eabc(x: int, v: CountVector | None = None) -> int:
    """π_{>3}(x) = E(x)+A(x)+B(x)+C(x); Nenner für S(x) und χ(x)."""
    if v is None:
        v = v_at_x(x)
    return v.total


def s_at_x(x: int, v: CountVector | None = None) -> SimplexPoint:
    """S(x) = (1/π(x))·V(x) auf dem Standard-3-Simplex (Summe 1)."""
    if v is None:
        v = v_at_x(x)
    denom = v.total
    if denom == 0:
        return SimplexPoint(0.0, 0.0, 0.0, 0.0)
    return SimplexPoint(
        e=v.e / denom,
        a=v.a / denom,
        b=v.b / denom,
        c=v.c / denom,
    )


def chi_at_x(x: int, v: CountVector | None = None) -> float:
    """χ(x) = ((E+C)−(A+B)) / π(x) — Beispiel-EABC-Invariante (Def. 4)."""
    if v is None:
        v = v_at_x(x)
    denom = v.total
    if denom == 0:
        return 0.0
    return ((v.e + v.c) - (v.a + v.b)) / denom


@dataclass(frozen=True, slots=True)
class FluctuationVector:
    """δ(x) = V(x) − (π(x)/4)·1 ∈ ℰ (Summe 0)."""

    e: float
    a: float
    b: float
    c: float

    def as_tuple(self) -> tuple[float, float, float, float]:
        return (self.e, self.a, self.b, self.c)

    def to_dict(self) -> dict[str, float]:
        return {"E": self.e, "A": self.a, "B": self.b, "C": self.c}

    @property
    def sum(self) -> float:
        return self.e + self.a + self.b + self.c

    @property
    def norm_sq(self) -> float:
        return self.e**2 + self.a**2 + self.b**2 + self.c**2


def delta_at_x(x: int, v: CountVector | None = None) -> FluctuationVector:
    """δ_i(x) = i(x) − π_{>3}(x)/4 für i ∈ {E,A,B,C}."""
    if v is None:
        v = v_at_x(x)
    pi = v.total
    if pi == 0:
        return FluctuationVector(0.0, 0.0, 0.0, 0.0)
    quarter = pi / 4.0
    return FluctuationVector(
        e=v.e - quarter,
        a=v.a - quarter,
        b=v.b - quarter,
        c=v.c - quarter,
    )


def chi_fluct_at_x(
    x: int,
    v: CountVector | None = None,
    delta: FluctuationVector | None = None,
) -> float:
    """χ_fluct(x) = (δ_E+δ_C)−(δ_A+δ_B) = (E+C)−(A+B); χ(x) = χ_fluct/π."""
    if delta is None:
        if v is None:
            v = v_at_x(x)
        delta = delta_at_x(x, v)
    return (delta.e + delta.c) - (delta.a + delta.b)


def h_at_x(
    x: int,
    v: CountVector | None = None,
    delta: FluctuationVector | None = None,
) -> float:
    """H(x) = ||δ(x)||² — EABC-Energie (perfekte Gleichverteilung ⟺ H=0)."""
    if delta is None:
        if v is None:
            v = v_at_x(x)
        delta = delta_at_x(x, v)
    return delta.norm_sq


def mode_coefficients(delta: FluctuationVector) -> tuple[float, float, float]:
    """Orthogonale Projektion δ = c₁Φ₁ + c₂Φ₂ + c₃Φ₃ mit c_i = (δ·Φ_i)/||Φ_i||²."""
    d = delta.as_tuple()
    return tuple(
        sum(d[j] * phi[j] for j in range(4)) / V4_MODE_NORM_SQ for phi in V4_MODES
    )


def _covariance_matrix(vectors: list[tuple[float, float, float, float]]) -> list[list[float]]:
    """Stichproben-Kovarianzmatrix (ddof=1) für 4-Tupel."""
    n = len(vectors)
    if n < 2:
        return [[0.0] * 4 for _ in range(4)]
    mean = [0.0] * 4
    for vec in vectors:
        for i in range(4):
            mean[i] += vec[i]
    mean = [m / n for m in mean]
    k = [[0.0] * 4 for _ in range(4)]
    for vec in vectors:
        for i in range(4):
            for j in range(4):
                k[i][j] += (vec[i] - mean[i]) * (vec[j] - mean[j])
    denom = n - 1
    return [[k[i][j] / denom for j in range(4)] for i in range(4)]


def _symmetric_eigh(matrix: list[list[float]]) -> tuple[list[float], list[list[float]]]:
    """Eigenwerte/-vektoren einer symmetrischen 4×4-Matrix (Spalten = Eigenvektoren)."""
    try:
        import numpy as np

        arr = np.array(matrix, dtype=float)
        eigenvalues, eigenvectors = np.linalg.eigh(arr)
        order = np.argsort(eigenvalues)[::-1]
        vals = [float(eigenvalues[i]) for i in order]
        vecs = [
            [float(eigenvectors[j, i]) for j in range(4)] for i in order
        ]
        return vals, vecs
    except ImportError:
        pass

    # Fallback: Potenziteration für dominante Richtung + Deflation (reicht für Tests).
    n = 4
    mat = [row[:] for row in matrix]
    eigenvalues: list[float] = []
    eigenvectors: list[list[float]] = []
    work = [row[:] for row in mat]
    for _ in range(n):
        v = [1.0 / sqrt(n)] * n
        for _ in range(80):
            w = [sum(work[i][j] * v[j] for j in range(n)) for i in range(n)]
            norm = sqrt(sum(x * x for x in w))
            if norm < 1e-15:
                break
            v = [x / norm for x in w]
        rayleigh = sum(v[i] * sum(work[i][j] * v[j] for j in range(n)) for i in range(n))
        eigenvalues.append(rayleigh)
        eigenvectors.append(v[:])
        for i in range(n):
            for j in range(n):
                work[i][j] -= rayleigh * v[i] * v[j]
    pairs = sorted(zip(eigenvalues, eigenvectors), key=lambda p: p[0], reverse=True)
    return [p[0] for p in pairs], [p[1] for p in pairs]


def fluctuation_covariance_at_grid(
    max_x: int, primes: list[int] | None = None
) -> dict[str, Any]:
    """Kovarianz K von δ entlang des Gitters x = 5..max_x (Primzähl-Punkte)."""
    if max_x < 5:
        return {
            "sample_count": 0,
            "matrix": [[0.0] * 4 for _ in range(4)],
            "eigenvalues": [],
            "eigenvectors": [],
            "top_eigenvector": None,
            "labels": ["E", "A", "B", "C"],
        }
    if primes is None:
        primes = _sieve_primes(max_x)
    deltas: list[tuple[float, float, float, float]] = []
    for x in range(5, max_x + 1):
        deltas.append(delta_at_x(x, v_at_x(x, primes)).as_tuple())
    k = _covariance_matrix(deltas)
    eigenvalues, eigenvectors = _symmetric_eigh(k)
    top = eigenvectors[0] if eigenvectors else None
    return {
        "sample_count": len(deltas),
        "matrix": k,
        "eigenvalues": eigenvalues,
        "eigenvectors": [
            {"lambda": eigenvalues[i], "vector": dict(zip("EABC", eigenvectors[i]))}
            for i in range(len(eigenvectors))
        ],
        "top_eigenvector": dict(zip("EABC", top)) if top else None,
        "labels": ["E", "A", "B", "C"],
    }


def sigma_quadruplet(p: int) -> str | None:
    """σ(Q) = (κ(p), κ(p+2), κ(p+6), κ(p+8)) als String über Σ_4."""
    if not is_prime_quadruplet(p):
        return None
    parts: list[str] = []
    for val in q(p):
        cls = kappa(val)
        if cls is None:
            return None
        parts.append(cls.value)
    return "".join(parts)


def quadruplet_signature_frequencies(
    max_p: int, primes: list[int] | None = None
) -> dict[str, Any]:
    """Schätzt μ auf Σ_4 durch Häufigkeiten der σ(Q) für Prim-Vierlinge mit p ≤ max_p."""
    if primes is None:
        primes = _sieve_primes(max_p + 8)
    counts: dict[str, int] = {}
    total = 0
    for p in primes:
        if p > max_p:
            break
        sig = sigma_quadruplet(p)
        if sig is None:
            continue
        counts[sig] = counts.get(sig, 0) + 1
        total += 1
    frequencies = {sig: c / total for sig, c in sorted(counts.items())} if total else {}
    return {
        "max_p": max_p,
        "quadruplet_count": total,
        "signature_counts": dict(sorted(counts.items())),
        "signature_frequencies": frequencies,
        "sigma4_alphabet": SIGMA4_ALPHABET,
        "sigma4_cardinality": len(SIGMA4_ALPHABET) ** 4,
    }


@dataclass(frozen=True, slots=True)
class Snapshot:
    x: int
    pi: int
    v: CountVector
    s: SimplexPoint
    chi: float
    delta: FluctuationVector
    chi_fluct: float
    h: float
    mode_c: tuple[float, float, float]


def snapshot_at_x(x: int, primes: list[int] | None = None) -> Snapshot:
    v = v_at_x(x, primes)
    delta = delta_at_x(x, v)
    return Snapshot(
        x=x,
        pi=v.total,
        v=v,
        s=s_at_x(x, v),
        chi=chi_at_x(x, v),
        delta=delta,
        chi_fluct=chi_fluct_at_x(x, v, delta),
        h=h_at_x(x, v, delta),
        mode_c=mode_coefficients(delta),
    )


def run_program(
    max_x: int,
    sample_points: list[int] | None = None,
) -> dict[str, Any]:
    primes = _sieve_primes(max_x + 8)
    if sample_points is None:
        sample_points = sorted(
            {
                20,
                100,
                500,
                1000,
                5000,
                max_x,
            }
            & {x for x in range(5, max_x + 1)}
        )
        if not sample_points:
            sample_points = [max_x] if max_x >= 5 else []

    samples = []
    chi_series: list[dict[str, float | int]] = []
    fluct_series: list[dict[str, float | int]] = []
    for x in range(5, max_x + 1):
        v = v_at_x(x, primes)
        delta = delta_at_x(x, v)
        chi_val = chi_at_x(x, v)
        chi_fluct = chi_fluct_at_x(x, v, delta)
        h_val = h_at_x(x, v, delta)
        pi = v.total
        c1, c2, c3 = mode_coefficients(delta)
        chi_series.append({"x": x, "chi": chi_val, "pi": pi})
        fluct_series.append(
            {
                "x": x,
                "pi": pi,
                "H": h_val,
                "H_over_pi": h_val / pi if pi else 0.0,
                "chi_fluct": chi_fluct,
                "chi_over_sqrt_pi": chi_val / sqrt(pi) if pi else 0.0,
                "c1": c1,
                "c2": c2,
                "c3": c3,
            }
        )
        if x in sample_points:
            snap = snapshot_at_x(x, primes)
            samples.append(
                {
                    "x": snap.x,
                    "pi_eabc": snap.pi,
                    "V": snap.v.to_dict(),
                    "S": snap.s.to_dict(),
                    "chi": snap.chi,
                    "delta": snap.delta.to_dict(),
                    "chi_fluct": snap.chi_fluct,
                    "H": snap.h,
                    "H_over_pi": snap.h / snap.pi if snap.pi else 0.0,
                    "chi_over_sqrt_pi": snap.chi / sqrt(snap.pi) if snap.pi else 0.0,
                    "mode_c": {
                        "c1": snap.mode_c[0],
                        "c2": snap.mode_c[1],
                        "c3": snap.mode_c[2],
                    },
                }
            )

    quad_stats = quadruplet_signature_frequencies(max_x, primes)
    cov_k = fluctuation_covariance_at_grid(max_x, primes)

    chi_values = [row["chi"] for row in chi_series]
    chi_trend = {
        "min": min(chi_values) if chi_values else 0.0,
        "max": max(chi_values) if chi_values else 0.0,
        "mean": sum(chi_values) / len(chi_values) if chi_values else 0.0,
        "last": chi_values[-1] if chi_values else 0.0,
        "first": chi_values[0] if chi_values else 0.0,
    }

    h_over_pi = [row["H_over_pi"] for row in fluct_series]
    chi_sqrt_pi = [row["chi_over_sqrt_pi"] for row in fluct_series]
    c1_vals = [row["c1"] for row in fluct_series]
    c2_vals = [row["c2"] for row in fluct_series]
    c3_vals = [row["c3"] for row in fluct_series]
    fluct_trend = {
        "H_over_pi": {
            "min": min(h_over_pi) if h_over_pi else 0.0,
            "max": max(h_over_pi) if h_over_pi else 0.0,
            "mean": sum(h_over_pi) / len(h_over_pi) if h_over_pi else 0.0,
            "last": h_over_pi[-1] if h_over_pi else 0.0,
        },
        "chi_over_sqrt_pi": {
            "min": min(chi_sqrt_pi) if chi_sqrt_pi else 0.0,
            "max": max(chi_sqrt_pi) if chi_sqrt_pi else 0.0,
            "mean": sum(chi_sqrt_pi) / len(chi_sqrt_pi) if chi_sqrt_pi else 0.0,
            "last": chi_sqrt_pi[-1] if chi_sqrt_pi else 0.0,
        },
        "mode_c1": {
            "min": min(c1_vals) if c1_vals else 0.0,
            "max": max(c1_vals) if c1_vals else 0.0,
            "last": c1_vals[-1] if c1_vals else 0.0,
        },
        "mode_c2": {
            "min": min(c2_vals) if c2_vals else 0.0,
            "max": max(c2_vals) if c2_vals else 0.0,
            "last": c2_vals[-1] if c2_vals else 0.0,
        },
        "mode_c3": {
            "min": min(c3_vals) if c3_vals else 0.0,
            "max": max(c3_vals) if c3_vals else 0.0,
            "last": c3_vals[-1] if c3_vals else 0.0,
        },
    }

    last_fluct = fluct_series[-1] if fluct_series else {}

    return {
        "program": "EABC-Invarianzprogramm",
        "canonical_doc": "collatz_eabc_invarianzprogramm.md",
        "max_x": max_x,
        "definitions": {
            "kappa": "P_{>3} → {E,A,B,C}, E≡1, A≡5, B≡7, C≡11 (mod 12)",
            "V": "Zählvektor (E,A,B,C) für p ≤ x, p > 3",
            "S": "V / π_{>3}(x) ∈ Δ_3",
            "chi": "((E+C)-(A+B)) / π_{>3}(x)",
            "delta": "δ(x)=V(x)−(π_{>3}(x)/4)·1 ∈ ℰ, Summe 0",
            "chi_fluct": "(δ_E+δ_C)−(δ_A+δ_B); χ=chi_fluct/π",
            "H": "||δ(x)||²",
            "K": "Stichproben-Kovarianz von δ auf Gitter x=5..max_x",
            "v4_modes": "Φ_1=(1,-1,1,-1), Φ_2=(1,1,-1,-1), Φ_3=(1,-1,-1,1)",
            "sigma_quadruplet": "σ(Q)=(κ(p),κ(p+2),κ(p+6),κ(p+8)) für Q=(p,p+2,p+6,p+8)",
        },
        "samples": samples,
        "chi_trend": chi_trend,
        "chi_series_tail": chi_series[-20:] if len(chi_series) > 20 else chi_series,
        "fluctuation_field": {
            "at_max_x": last_fluct,
            "trend": fluct_trend,
            "covariance_K": cov_k,
            "research_A": {
                "H_over_pi_limit_question": "lim sup / lim inf von H(x)/π_{>3}(x)",
                "chi_over_sqrt_pi_limit_question": "lim sup / lim inf von χ(x)/√π_{>3}(x)",
            },
            "research_B": "Spektralzerlegung von K; Eigenvektoren vs. E,A,B,C",
            "research_C": "Modenkoeffizienten c_i via Φ_1,Φ_2,Φ_3 (Klein V_4)",
            "spectral_conjecture": (
                "Statistik der c_i(x) kodiert Information zu L-Funktionen mod 12 / ζ-Nullstellen"
            ),
        },
        "quadruplet_signatures": quad_stats,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC-Invarianzprogramm numerisch")
    parser.add_argument("--max-x", type=int, default=10_000, help="Obergrenze für x")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="JSON-Ausgabedatei",
    )
    args = parser.parse_args()
    result = run_program(args.max_x)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(f"Wrote {args.output} (max_x={args.max_x})")
    for row in result["samples"]:
        print(
            f"x={row['x']:>6}  S={tuple(round(row['S'][k], 4) for k in 'EABC')}  "
            f"χ={row['chi']:+.4f}  H/π={row['H_over_pi']:.4f}  "
            f"χ/√π={row['chi_over_sqrt_pi']:+.4f}  "
            f"c=({row['mode_c']['c1']:+.3f},{row['mode_c']['c2']:+.3f},{row['mode_c']['c3']:+.3f})"
        )
    cov = result["fluctuation_field"]["covariance_K"]
    if cov.get("top_eigenvector"):
        ev = cov["top_eigenvector"]
        lam = cov["eigenvalues"][0] if cov["eigenvalues"] else 0.0
        print(
            f"K top-λ={lam:.6g}  v*={tuple(round(ev[k], 4) for k in 'EABC')}"
        )


if __name__ == "__main__":
    main()
