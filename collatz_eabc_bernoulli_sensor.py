#!/usr/bin/env python3
"""
EABC-Resonanzhypothese-Sensor: von-Staudt-Clausen-Primsignaturen als EABC-Zustände.

EABC-Zerlegungsprinzip: N = (N_glatt, N_EABC). Dieser Sensor liest N_EABC entlang der
Bernoulli-Brücke s=-2n → ζ(1-2n) → PrimeSig(B_{2n}) → V_n.

Φ(n) = V_n = (E_n, A_n, B_n, C_n) zählt Primzahlen p ∈ P_n mit p-1 | 2n nach nativer
EABC-Restklasse (E≡1, A≡5, B≡7, C≡11 mod 12). Bernoulli = Übersetzungsobjekt, nicht Endobjekt.

Notation: Python-Felder e,a,b,c entsprechen Dokumentation E_n, A_n, B_n, C_n.

Ausführung:
    python3 collatz_eabc_bernoulli_sensor.py
    python3 collatz_eabc_bernoulli_sensor.py --max-n 200 --output collatz_eabc_bernoulli_sensor.json
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from math import isqrt
from pathlib import Path
from typing import Any

from eabc_from_lean import EClass, class_of

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_bernoulli_sensor.json"

EABC_RESIDUES = {1, 5, 7, 11}


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


def prime_sig(two_n: int, primes: list[int] | None = None) -> list[int]:
    """PrimeSig(B_{2n}) = {p Primzahl : p-1 teilt 2n}."""
    if two_n < 2:
        return []
    if primes is None:
        primes = _sieve_primes(two_n + 1)
    return [p for p in primes if p <= two_n + 1 and (two_n % (p - 1) == 0)]


def staudt_denominator(two_n: int, primes: list[int] | None = None) -> int:
    """den(B_{2n}) = ∏_{p-1|2n} p (von Staudt--Clausen)."""
    sig = prime_sig(two_n, primes)
    out = 1
    for p in sig:
        out *= p
    return out


@dataclass(frozen=True, slots=True)
class EabcVector:
    """EABC-Zustandsvektor V_n = (E_n, A_n, B_n, C_n); Felder e,a,b,c = E_n,…,C_n."""

    e: int  # E_n: Anzahl p ≡ 1 (mod 12) in P_n
    a: int  # A_n: Anzahl p ≡ 5 (mod 12)
    b: int  # B_n: Anzahl p ≡ 7 (mod 12)
    c: int  # C_n: Anzahl p ≡ 11 (mod 12)

    def as_tuple(self) -> tuple[int, int, int, int]:
        return (self.e, self.a, self.b, self.c)

    def as_doc_dict(self) -> dict[str, int]:
        """Alias-Schlüssel E_n, A_n, B_n, C_n (Dokumentationsnotation)."""
        return {"E_n": self.e, "A_n": self.a, "B_n": self.b, "C_n": self.c}

    @property
    def total(self) -> int:
        return self.e + self.a + self.b + self.c


def v_bernoulli(sig: list[int]) -> EabcVector:
    """V_n = (E_n, A_n, B_n, C_n): EABC-Zählvektor über P_n (p ≡ 1,5,7,11 mod 12)."""
    counts = {EClass.E: 0, EClass.A: 0, EClass.B: 0, EClass.C: 0}
    for p in sig:
        cls = class_of(p)
        if cls is not None:
            counts[cls] += 1
    return EabcVector(
        e=counts[EClass.E],
        a=counts[EClass.A],
        b=counts[EClass.B],
        c=counts[EClass.C],
    )


def q4_vector(n: int, primes: list[int] | None = None) -> EabcVector:
    """Q_4(N): EABC-Zählvektor über alle Primzahlen p ≤ N in den vier Klassen (§2.3)."""
    if n < 2:
        return EabcVector(0, 0, 0, 0)
    if primes is None:
        primes = _sieve_primes(n)
    return v_bernoulli([p for p in primes if p <= n])


def delta_q4(v: EabcVector) -> dict[str, int]:
    """ΔQ_4(N): chirale Asymmetrie-Observablen aus Q_4(N) — σ, χ, A−C (§2.3)."""
    return {
        "sigma": sigma_eabc(v),
        "chi": chi_eabc(v),
        "a_minus_c": v.a - v.c,
        "i_chir": i_chir(v),
    }


def non_eabc_primes(sig: list[int]) -> list[int]:
    """Primzahlen in PrimeSig außerhalb der vier EABC-Klassen (typisch 2, 3)."""
    return [p for p in sig if class_of(p) is None]


def sigma_eabc(v: EabcVector) -> int:
    """σ: EA- vs. BC-Bilanz (grobe EABC-Symmetrie)."""
    return (v.e + v.a) - (v.b + v.c)


def chi_eabc(v: EabcVector) -> int:
    """χ: Diagonalen-Bilanz (E+B) − (A+C)."""
    return (v.e + v.b) - (v.a + v.c)


def i_chir(v: EabcVector) -> int:
    """ι_chir: chiral Vorzeichen-Index sgn(σ·χ), 0 falls σ=0 oder χ=0."""
    s, c = sigma_eabc(v), chi_eabc(v)
    prod = s * c
    if prod > 0:
        return 1
    if prod < 0:
        return -1
    return 0


@dataclass(frozen=True, slots=True)
class BernoulliRow:
    n: int
    two_n: int
    prime_sig: list[int]
    non_eabc: list[int]
    staudt_denominator: int
    v: EabcVector
    sigma: int
    chi: int
    i_chir: int


def bernoulli_row(n: int, primes: list[int]) -> BernoulliRow:
    two_n = 2 * n
    sig = prime_sig(two_n, primes)
    vec = v_bernoulli(sig)
    return BernoulliRow(
        n=n,
        two_n=two_n,
        prime_sig=sig,
        non_eabc=non_eabc_primes(sig),
        staudt_denominator=staudt_denominator(two_n, primes),
        v=vec,
        sigma=sigma_eabc(vec),
        chi=chi_eabc(vec),
        i_chir=i_chir(vec),
    )


def _asymmetry_stats(rows: list[BernoulliRow]) -> dict[str, Any]:
    sigmas = [r.sigma for r in rows]
    chis = [r.chi for r in rows]
    chirals = [r.i_chir for r in rows]
    nonzero_sigma = sum(1 for s in sigmas if s != 0)
    nonzero_chi = sum(1 for c in chis if c != 0)
    return {
        "sigma_mean": sum(sigmas) / len(sigmas) if rows else 0.0,
        "chi_mean": sum(chis) / len(chis) if rows else 0.0,
        "sigma_nonzero_fraction": nonzero_sigma / len(rows) if rows else 0.0,
        "chi_nonzero_fraction": nonzero_chi / len(rows) if rows else 0.0,
        "i_chir_counts": {
            "+1": sum(1 for x in chirals if x == 1),
            "0": sum(1 for x in chirals if x == 0),
            "-1": sum(1 for x in chirals if x == -1),
        },
        "max_abs_sigma": max((abs(s) for s in sigmas), default=0),
        "max_abs_chi": max((abs(c) for c in chis), default=0),
    }


def run_sensor(max_n: int = 100) -> dict[str, Any]:
    two_n_max = 2 * max_n
    primes = _sieve_primes(two_n_max + 1)
    rows = [bernoulli_row(n, primes) for n in range(1, max_n + 1)]

    samples = []
    for r in rows:
        samples.append(
            {
                "n": r.n,
                "two_n": r.two_n,
                "prime_sig": r.prime_sig,
                "non_eabc_primes": r.non_eabc,
                "staudt_denominator": r.staudt_denominator,
                "V": {"E": r.v.e, "A": r.v.a, "B": r.v.b, "C": r.v.c},
                **r.v.as_doc_dict(),
                "sigma": r.sigma,
                "chi": r.chi,
                "i_chir": r.i_chir,
            }
        )

    return {
        "framework": "EABC",
        "hypothesis": "EABC-Resonanzhypothese der Zetafunktion",
        "sensor": "Phi",
        "description": (
            "EABC-Zustände V_n=(E_n,A_n,B_n,C_n) aus von-Staudt-Clausen-Signaturen P_n; "
            "mod-12: E≡1, A≡5, B≡7, C≡11. Falsifikation: V_n vs. Δt_k."
        ),
        "epistemic_label": "Experiment",
        "max_n": max_n,
        "samples": samples,
        "stats": _asymmetry_stats(rows),
        "future_tests": {
            "zeta_coupling": "V_n vs. Δt_k (Resonanzhypothese, benötigt Nullstellendaten)",
            "curvature": "K_B(n) vs. π(x)-Li(x) (zukünftig)",
        },
    }


def verify_staudt_sympy(max_n: int = 20) -> list[dict[str, Any]]:
    """Optionaler Kreuzcheck mit sympy.bernoulli (nur für Tests/Diagnose)."""
    try:
        from sympy import Rational, bernoulli
    except ImportError:
        return []
    primes = _sieve_primes(2 * max_n + 1)
    checks = []
    for n in range(1, max_n + 1):
        two_n = 2 * n
        den_sc = staudt_denominator(two_n, primes)
        den_sym = int(Rational(bernoulli(two_n)).q)
        checks.append(
            {
                "n": n,
                "staudt": den_sc,
                "sympy": den_sym,
                "match": den_sc == den_sym,
            }
        )
    return checks


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC-Resonanz-Sensor Φ(n)=V_n")
    parser.add_argument("--max-n", type=int, default=100, help="Obergrenze für n (B_{2n})")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="JSON-Ausgabedatei",
    )
    args = parser.parse_args()
    if args.max_n < 1:
        raise SystemExit("--max-n muss ≥ 1 sein")

    report = run_sensor(args.max_n)
    args.output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Geschrieben: {args.output} ({len(report['samples'])} Zeilen, n=1..{args.max_n})")


if __name__ == "__main__":
    main()
