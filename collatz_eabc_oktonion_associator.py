#!/usr/bin/env python3
"""
EABC-Associator-Observable für Oktanionen (Γ-Differenz bei Klammerung).

Kanonsiche Theorie: collatz_eabc_oktonion_singularitaet.md §3.7–§3.8
                   collatz_eabc_plattenuebergang.md §2.6–§2.7

Fundamentale Größen:
  [x,y,z] = (xy)z - x(yz)                    (algebraischer Assoziator)
  A(T1,T2) = ||P_T1(...) - P_T2(...)||       (Baum-Abhängigkeit)
  Γ_EABC-Assoziator = smooth_Γ((xy)z) - smooth_Γ(x(yz))  (16 glatt-EABC-Koordinaten)
  𝔞_E(n) = mittlere EABC-Assoziator-Norm über Stichproben zu n = abc (a,b,c≥2)

Spektrum M_n(t): siehe collatz_eabc_oktonion_spectrum.py — Constraint N(x)N(y)N(z)=n,
Histogramm der Assoziator-Norm α=N([x,y,z]) über 𝔞_n (nicht nur Schalen-Mittel).

Sampling (n ≤ 50): repräsentative Faktorisierungen n = abc, a,b,c ≥ 2;
keine volle Σ_n-Enumeration (r_8(n) ~ n^{3/2}); Z^8-Stub (Hurwitz-Limitation).

Ausführung:
    python3 collatz_eabc_oktonion_associator.py
    python3 collatz_eabc_oktonion_associator.py --max-n 50 --samples 30
    python3 collatz_eabc_oktonion_spectrum.py   # M_n(t)-Histogramme
"""

from __future__ import annotations

import argparse
import json
import math
import random
import statistics
from functools import lru_cache
from math import isqrt
from pathlib import Path
from typing import Any

from collatz_eabc_gauss_spaltung_test import kappa_glatt

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_oktonion_associator.json"

Oct = tuple[int, ...]  # 8-tuple, Z^8 stub

# Fano-Zyklen (1,2,3), (1,4,5), (1,7,6), (2,4,6), (2,5,7), (3,4,7), (3,5,6)
_FANO_CYCLES: tuple[tuple[int, int, int], ...] = (
    (1, 2, 3),
    (1, 4, 5),
    (1, 7, 6),
    (2, 4, 6),
    (2, 5, 7),
    (3, 4, 7),
    (3, 5, 6),
)


def _build_basis_mul_table() -> list[list[tuple[int, int]]]:
    """e_i * e_j = sign * e_k; Index 0 = Einheit."""
    table: list[list[tuple[int, int] | None]] = [[None] * 8 for _ in range(8)]
    for i in range(8):
        table[0][i] = (i, 1)
        table[i][0] = (i, 1)
        if i == 0:
            table[0][0] = (0, 1)
        else:
            table[i][i] = (0, -1)
    for i, j, k in _FANO_CYCLES:
        table[i][j] = (k, 1)
        table[j][i] = (k, -1)
        table[j][k] = (i, 1)
        table[k][j] = (i, -1)
        table[k][i] = (j, 1)
        table[i][k] = (j, -1)
    out: list[list[tuple[int, int]]] = []
    for row in table:
        out.append([(0, 1) if x is None else x for x in row])
    return out


_BASIS_MUL = _build_basis_mul_table()


def o_add(x: Oct, y: Oct) -> Oct:
    return tuple(a + b for a, b in zip(x, y))


def o_scale(s: int, x: Oct) -> Oct:
    return tuple(s * a for a in x)


def o_neg(x: Oct) -> Oct:
    return o_scale(-1, x)


def o_mul(x: Oct, y: Oct) -> Oct:
    """Oktanionen-Multiplikation auf Z^8 (Standard-Cayley-Tabelle)."""
    acc = [0] * 8
    for i in range(8):
        if x[i] == 0:
            continue
        for j in range(8):
            if y[j] == 0:
                continue
            k, sign = _BASIS_MUL[i][j]
            acc[k] += sign * x[i] * y[j]
    return tuple(acc)


def o_norm_sq(x: Oct) -> int:
    return sum(c * c for c in x)


def o_norm(x: Oct) -> float:
    return math.sqrt(float(o_norm_sq(x)))


def associator(x: Oct, y: Oct, z: Oct) -> Oct:
    """[x,y,z] = (xy)z - x(yz)."""
    left = o_mul(o_mul(x, y), z)
    right = o_mul(x, o_mul(y, z))
    return o_add(left, o_neg(right))


def associator_norm(x: Oct, y: Oct, z: Oct) -> float:
    return o_norm(associator(x, y, z))


def leg_smooth(coord: int) -> tuple[int, int]:
    """Glatt-EABC (α, β) pro Koordinate; (0,0) bei coord=0."""
    if coord == 0:
        return (0, 0)
    alpha, beta, _core, _ka = kappa_glatt(abs(coord))
    return (alpha, beta)


def smooth_gamma8(x: Oct) -> tuple[int, ...]:
    """16-Tupel (α_1,β_1,…,α_8,β_8) — glatt-gestrippte EABC-Koordinaten."""
    out: list[int] = []
    for c in x:
        a, b = leg_smooth(c)
        out.extend((a, b))
    return tuple(out)


def eabc_associator_vector(x: Oct, y: Oct, z: Oct) -> tuple[int, ...]:
    """Γ((xy)z) - Γ(x(yz)) auf glatt-EABC-Koordinaten."""
    gl = smooth_gamma8(o_mul(o_mul(x, y), z))
    gr = smooth_gamma8(o_mul(x, o_mul(y, z)))
    return tuple(a - b for a, b in zip(gl, gr))


def eabc_associator_norm(x: Oct, y: Oct, z: Oct) -> float:
    v = eabc_associator_vector(x, y, z)
    return math.sqrt(float(sum(c * c for c in v)))


def enumerate_shell_z8(n: int) -> list[Oct]:
    """Alle x ∈ Z^8 mit Σ x_i² = n (Stub-Gitter)."""
    if n < 1:
        return []
    out: list[Oct] = []
    bound = isqrt(n)
    for x1 in range(-bound, bound + 1):
        s1 = x1 * x1
        if s1 > n:
            continue
        for x2 in range(-bound, bound + 1):
            s2 = s1 + x2 * x2
            if s2 > n:
                continue
            for x3 in range(-bound, bound + 1):
                s3 = s2 + x3 * x3
                if s3 > n:
                    continue
                for x4 in range(-bound, bound + 1):
                    s4 = s3 + x4 * x4
                    if s4 > n:
                        continue
                    for x5 in range(-bound, bound + 1):
                        s5 = s4 + x5 * x5
                        if s5 > n:
                            continue
                        for x6 in range(-bound, bound + 1):
                            s6 = s5 + x6 * x6
                            if s6 > n:
                                continue
                            for x7 in range(-bound, bound + 1):
                                s7 = s6 + x7 * x7
                                if s7 > n:
                                    continue
                                rem = n - s7
                                d = isqrt(rem)
                                if d * d != rem:
                                    continue
                                if d == 0:
                                    out.append((x1, x2, x3, x4, x5, x6, x7, 0))
                                else:
                                    for x8 in (d, -d):
                                        out.append((x1, x2, x3, x4, x5, x6, x7, x8))
    return out


@lru_cache(maxsize=None)
def shell_z8(n: int) -> tuple[Oct, ...]:
    return tuple(enumerate_shell_z8(n))


def shell_size_z8(n: int) -> int:
    return len(shell_z8(n))


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def triple_factorizations(n: int, min_factor: int = 2) -> list[tuple[int, int, int]]:
    """Ungeordnete Tripel (a,b,c) mit a≤b≤c, a,b,c≥min_factor, abc=n."""
    if n < min_factor**3:
        return []
    seen: set[tuple[int, int, int]] = set()
    for a in range(min_factor, int(n ** (1 / 3)) + 2):
        if n % a != 0:
            continue
        rest = n // a
        for b in range(max(a, min_factor), isqrt(rest) + 1):
            if rest % b != 0:
                continue
            c = rest // b
            if c < min_factor:
                continue
            key = tuple(sorted((a, b, c)))
            seen.add(key)  # type: ignore[arg-type]
    return sorted(seen)


def _sample_shell(n: int, k: int, rng: random.Random) -> list[Oct]:
    pts = list(shell_z8(n))
    if not pts:
        return []
    if len(pts) <= k:
        return pts
    return rng.sample(pts, k)


def bracketing_magnitudes(
    x: Oct, y: Oct, z: Oct
) -> dict[str, float]:
    """Algebraische und EABC-Assoziator-Normen für (xy)z vs x(yz)."""
    return {
        "algebraic_associator_norm": associator_norm(x, y, z),
        "eabc_associator_norm": eabc_associator_norm(x, y, z),
        "tree_product_norm_diff": o_norm(
            o_add(o_mul(o_mul(x, y), z), o_neg(o_mul(x, o_mul(y, z))))
        ),
    }


def a_E_for_factorization(
    a: int,
    b: int,
    c: int,
    samples: int,
    rng: random.Random,
    max_shell_points: int = 400,
) -> dict[str, Any]:
    """Mittelwerte über Stichproben x∈Σ_a, y∈Σ_b, z∈Σ_c."""
    xs = _sample_shell(a, min(samples, max_shell_points), rng)
    ys = _sample_shell(b, min(samples, max_shell_points), rng)
    zs = _sample_shell(c, min(samples, max_shell_points), rng)
    if not xs or not ys or not zs:
        return {
            "factors": [a, b, c],
            "sample_count": 0,
            "error": "empty_shell",
        }

    alg: list[float] = []
    eabc: list[float] = []
    trials = 0
    cap = samples
    for _ in range(cap):
        x = rng.choice(xs)
        y = rng.choice(ys)
        z = rng.choice(zs)
        m = bracketing_magnitudes(x, y, z)
        alg.append(m["algebraic_associator_norm"])
        eabc.append(m["eabc_associator_norm"])
        trials += 1

    return {
        "factors": [a, b, c],
        "shell_sizes": [shell_size_z8(a), shell_size_z8(b), shell_size_z8(c)],
        "sample_count": trials,
        "mean_algebraic_associator_norm": statistics.mean(alg) if alg else 0.0,
        "mean_eabc_associator_norm": statistics.mean(eabc) if eabc else 0.0,
        "max_eabc_associator_norm": max(eabc) if eabc else 0.0,
    }


def a_E_for_n(
    n: int,
    samples: int,
    rng: random.Random,
) -> dict[str, Any]:
    """𝔞_E(n): Mittel über repräsentative Tripel-Faktorisierungen."""
    if _is_prime(n):
        return {
            "n": n,
            "is_prime": True,
            "triple_factorizations": 0,
            "a_E_algebraic": None,
            "a_E_eabc": None,
            "note": "Prim: keine nichttriviale abc-Zerlegung; 𝔞_E nicht definiert (0 by convention)",
        }

    facs = triple_factorizations(n)
    if not facs:
        return {
            "n": n,
            "is_prime": False,
            "triple_factorizations": 0,
            "a_E_algebraic": None,
            "a_E_eabc": None,
            "note": "zusammengesetzt, aber kein abc mit a,b,c≥2",
        }

    per_fac = [a_E_for_factorization(a, b, c, samples, rng) for a, b, c in facs]
    valid = [p for p in per_fac if p.get("sample_count", 0) > 0]
    eabc_vals = [p["mean_eabc_associator_norm"] for p in valid]
    alg_vals = [p["mean_algebraic_associator_norm"] for p in valid]

    return {
        "n": n,
        "is_prime": False,
        "triple_factorizations": len(facs),
        "factorization_details": per_fac,
        "a_E_algebraic": statistics.mean(alg_vals) if alg_vals else None,
        "a_E_eabc": statistics.mean(eabc_vals) if eabc_vals else None,
        "max_factor_eabc": max(eabc_vals) if eabc_vals else None,
    }


def canonical_triples_test() -> dict[str, Any]:
    """Referenz: Quaternion-Teilalgebra vs. generisches O-Tripel."""
    e1 = (0, 1, 0, 0, 0, 0, 0, 0)
    e2 = (0, 0, 1, 0, 0, 0, 0, 0)
    e3 = (0, 0, 0, 1, 0, 0, 0, 0)
    e4 = (0, 0, 0, 0, 1, 0, 0, 0)
    quat = bracketing_magnitudes(e1, e2, e3)
    generic = bracketing_magnitudes(e1, e2, e4)
    # EABC auf Einheitsbasis: glatt-Γ oft identisch (nur Vorzeichen) → 0
    x = (-1, -1, 0, 0, 0, 0, 0, 0)
    y = (-1, 0, 0, 0, -1, 0, 0, 0)
    z = (-1, 0, -1, -1, 0, 0, 0, 0)
    shell_eabc = bracketing_magnitudes(x, y, z)
    return {
        "quaternion_subalgebra_e1_e2_e3": quat,
        "generic_o_e1_e2_e4": generic,
        "shell_sample_sigma2_sigma2_sigma3": shell_eabc,
        "quaternion_associator_zero": quat["algebraic_associator_norm"] < 1e-9,
        "generic_associator_nonzero": generic["algebraic_associator_norm"] > 0.0,
        "shell_eabc_associator_nonzero": shell_eabc["eabc_associator_norm"] > 0.0,
    }


def prime_vs_composite_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    prime_eabc: list[float] = []
    comp_eabc: list[float] = []
    comp_alg: list[float] = []
    for r in rows:
        if r.get("a_E_eabc") is None:
            continue
        if r.get("is_prime"):
            prime_eabc.append(0.0)  # convention: undefined → 0
        else:
            comp_eabc.append(float(r["a_E_eabc"]))
            if r.get("a_E_algebraic") is not None:
                comp_alg.append(float(r["a_E_algebraic"]))

    def _safe_mean(xs: list[float]) -> float | None:
        return statistics.mean(xs) if xs else None

    ratio = None
    if comp_eabc and _safe_mean(comp_eabc):
        m_comp = _safe_mean(comp_eabc) or 0.0
        m_prime = 0.0
        ratio = m_prime / m_comp if m_comp > 0 else None

    return {
        "composite_count_with_a_E": len(comp_eabc),
        "prime_count": sum(1 for r in rows if r.get("is_prime")),
        "mean_a_E_eabc_composite": _safe_mean(comp_eabc),
        "mean_a_E_eabc_prime_convention_zero": 0.0,
        "mean_a_E_algebraic_composite": _safe_mean(comp_alg),
        "ratio_prime_over_composite_eabc": ratio,
        "bias_visible": (
            "Primschalen minimieren 𝔞_E (by definition: keine abc-Zerlegung); "
            "zusammengesetzte n zeigen positive mittlere EABC-Assoziator-Norm "
            f"({len(comp_eabc)} n mit Daten, Mittel={_safe_mean(comp_eabc):.4f})"
            if comp_eabc
            else "keine zusammengesetzten abc-Daten"
        ),
        "sample_a_E_composite_first_10": [
            {"n": r["n"], "a_E_eabc": r["a_E_eabc"]}
            for r in rows
            if r.get("a_E_eabc") is not None
        ][:10],
    }


def associator_report(
    max_n: int = 50,
    samples: int = 25,
    seed: int = 54,
) -> dict[str, Any]:
    rng = random.Random(seed)
    rows = [a_E_for_n(n, samples, rng) for n in range(2, max_n + 1)]
    summary = prime_vs_composite_summary(rows)
    canon = canonical_triples_test()

    return {
        "meta": {
            "hypothesis_doc": [
                "collatz_eabc_oktonion_assoziator_spektralhypothese.md",
                "collatz_eabc_oktonion_singularitaet.md §3.7–§3.8",
                "collatz_eabc_plattenuebergang.md §2.6–§2.7",
            ],
            "script": "collatz_eabc_oktonion_associator.py",
            "lattice": "Z^8 stub (nicht Hurwitz O_H)",
            "max_n": max_n,
            "samples_per_factorization": samples,
            "seed": seed,
            "computed": [
                "algebraic associator ||(xy)z - x(yz)||",
                "EABC associator ||Γ((xy)z) - Γ(x(yz))||_smooth",
                "𝔞_E(n) over representative n=abc factorizations",
            ],
            "not_computed": [
                "full Σ_n enumeration for large n",
                "μ_n shell measure averaging",
                "Hurwitz O_H lattice",
                "all Catalan trees beyond single (xy)z vs x(yz) split",
            ],
            "epistemic_note": (
                "Sampling auf repräsentativen abc-Faktorisierungen; "
                "𝔞_E(p)=0 konventionell für Prim (keine nichttriviale Klammerung). "
                "Sichtbarer Bias: Primschalen minimieren 𝔞_E per Definition; "
                "Frage ist ob zusammengesetzte Schalen 𝔞_E maximieren oder Profil tragen."
            ),
        },
        "canonical_tests": canon,
        "rows": rows,
        "prime_vs_composite": summary,
    }


def format_summary(report: dict[str, Any]) -> str:
    s = report["prime_vs_composite"]
    c = report["canonical_tests"]
    lines = [
        "EABC-Oktanion-Associator (𝔞_E)",
        "=" * 40,
        f"n ≤ {report['meta']['max_n']}, samples={report['meta']['samples_per_factorization']}",
        f"Quaternion e1,e2,e3: assoc≈0? {c['quaternion_associator_zero']}",
        f"Generisch e1,e2,e4: assoc>0? {c['generic_associator_nonzero']}",
        f"  alg={c['generic_o_e1_e2_e4']['algebraic_associator_norm']:.4f}",
        f"  EABC(unit)={c['generic_o_e1_e2_e4']['eabc_associator_norm']:.4f}",
        f"Shell Σ2×Σ2×Σ3 EABC>0? {c['shell_eabc_associator_nonzero']}",
        "",
        f"Zusammengesetzt mit 𝔞_E: {s['composite_count_with_a_E']}",
        f"Mittel 𝔞_E(composite) = {s['mean_a_E_eabc_composite']}",
        f"Mittel 𝔞_E(prime conv.) = {s['mean_a_E_eabc_prime_convention_zero']}",
        "",
        s["bias_visible"],
    ]
    return "\n".join(lines)


def run(
    max_n: int = 50,
    samples: int = 25,
    seed: int = 54,
    output: Path | None = None,
) -> dict[str, Any]:
    report = associator_report(max_n=max_n, samples=samples, seed=seed)
    out = output or DEFAULT_OUTPUT
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report["output_path"] = str(out)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC-Oktanion-Associator 𝔞_E(n)")
    parser.add_argument("--max-n", type=int, default=50)
    parser.add_argument("--samples", type=int, default=25)
    parser.add_argument("--seed", type=int, default=54)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(
        max_n=args.max_n,
        samples=args.samples,
        seed=args.seed,
        output=args.output,
    )
    print(format_summary(report))
    print(f"\nJSON: {report['output_path']}")


if __name__ == "__main__":
    main()
