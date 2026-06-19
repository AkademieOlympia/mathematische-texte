#!/usr/bin/env python3
"""
A-Absorption und A-Transfer (§4.8.26)

Mehrschichtschale: R² = 2^n + 3^m + Σ_j A_j  mit A_j ≡ 5 (mod 12).
Operator A_comp: S_mehr → S_einfach,  2^n + 3^m + Σ A_j  ↦  2^{n'} + 3^{m'}.

Epistemische Trennung:
  - Multiplikative FTA: q^α eindeutig (Schicht A)
  - Additive Brücke 5 = 2 + 3 (≠ 2·3): Analogie / log-Umverteilung (Schicht C)
  - A_comp als definierter Operator (Schicht B)

Referenz: collatz_eabc_zirkulationshypothese.md §4.8.20–4.8.26
"""

from __future__ import annotations

import argparse
import json
import math
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np

from eabc_shell_prime_distribution import (
    LABEL_TO_IDX,
    LABELS,
    charge_observables,
    class_of,
    defect_delta,
    sieve_primes,
)

DEFAULT_OUTPUT = Path("eabc_a_absorption.json")

# Primachsen q ≡ 5 (mod 12) in A-Klasse (EABC RESIDUE_MAP: 5 → A)
A_CLASS_RESIDUE = 5


def is_a_class_prime(q: int) -> bool:
    """q ist Primachse in EABC-Klasse A (q ≡ 5 mod 12)."""
    return q >= 5 and q % 12 == A_CLASS_RESIDUE


def a_term_class(q: int) -> str | None:
    """EABC-Klasse eines A-Achsen-Terms q^α (Prim q → mod 12)."""
    if q == 2:
        return "E"
    if q == 3:
        return "E"
    return class_of(q) if q > 3 else None


def multilayer_r_squared(
    n: int,
    m: int,
    a_terms: list[tuple[int, int]] | None = None,
) -> int:
    """R² = 2^n + 3^m + Σ q^α für A-Terme (q, α)."""
    total = 2**n + 3**m
    for q, alpha in a_terms or []:
        if alpha < 1:
            raise ValueError("Exponent α ≥ 1 erforderlich")
        if not is_a_class_prime(q):
            raise ValueError(f"Achse q={q} nicht in A-Klasse (q ≡ 5 mod 12)")
        total += q**alpha
    return int(total)


def simple_shell_r_squared(n: int, m: int) -> int:
    return 2**n + 3**m


def absorption_defect(r_sq: int, r_tilde_sq: int) -> int:
    return abs(int(r_sq) - int(r_tilde_sq))


def compression_grade(k_multi: int) -> float:
    """Kompressionsgrad: k Mehrschicht-Terme → 2 Einfachschalen-Terme."""
    return float(k_multi) / 2.0


def count_multilayer_terms(n: int, m: int, a_terms: list[tuple[int, int]] | None) -> int:
    k = 0
    if n > 0:
        k += 1
    if m > 0:
        k += 1
    k += len(a_terms or [])
    return max(k, 1)


def absorb_exact_search(
    r_sq: int,
    n_max: int,
    m_max: int,
) -> tuple[int, int, int, int]:
    """
    A_abs via exakte Suche: minimales (n', m') mit minimalem |R² − R̃²|.
    Liefert (n', m', R̃², Defekt).
    """
    best_defect = r_sq + 1
    best = (0, 0, simple_shell_r_squared(0, 0))
    for n in range(n_max + 1):
        for m in range(m_max + 1):
            rt = simple_shell_r_squared(n, m)
            d = abs(r_sq - rt)
            if d < best_defect or (d == best_defect and (n + m) < (best[0] + best[1])):
                best_defect = d
                best = (n, m, rt)
    return best[0], best[1], best[2], best_defect


def absorb_greedy_log(
    n: int,
    m: int,
    a_terms: list[tuple[int, int]] | None,
) -> tuple[int, int, int, int, str]:
    """
    Greedy/log-Umverteilung: additive Brücke 5 = 2 + 3 (NICHT multiplikativ).

    Für jeden A-Term q^α: log-Beitrag ℓ·log(q) wird nach Gewichten
    w_2 = log(2)/(log 2 + log 3), w_3 = log(3)/(log 2 + log 3) auf Exponenten
    umverteilt (gerundet). Basis n,m bleibt; A-Beiträge addieren sich in n', m'.

    Liefert (n', m', R̃², Defekt, method='greedy_log').
    """
    log2 = math.log(2)
    log3 = math.log(3)
    w2 = log2 / (log2 + log3)
    w3 = log3 / (log2 + log3)
    n_prime = int(n)
    m_prime = int(m)
    for q, alpha in a_terms or []:
        contrib = alpha * math.log(q)
        n_prime += max(0, int(round(contrib * w2 / log2)))
        m_prime += max(0, int(round(contrib * w3 / log3)))
    r_tilde = simple_shell_r_squared(n_prime, m_prime)
    r_sq = multilayer_r_squared(n, m, a_terms)
    return n_prime, m_prime, r_tilde, absorption_defect(r_sq, r_tilde), "greedy_log"


def absorb_a_term_additive_bridge(q: int, alpha: int, exp_max: int) -> tuple[int, int, int, int]:
    """
    Absorption eines isolierten A-Terms q^α → (α', β') mit 2^{α'}+3^{β'} ≈ q^α.
    Additive Brücke: Referenz 5 = 2 + 3 (Analogie), Suche minimaler Defekt.
    """
    target = q**alpha
    return absorb_exact_search(target, exp_max, exp_max)[:3] + (
        absorption_defect(target, absorb_exact_search(target, exp_max, exp_max)[2]),
    )


def apply_a_comp(
    n: int,
    m: int,
    a_terms: list[tuple[int, int]] | None = None,
    *,
    method: str = "exact",
    exp_max: int = 12,
) -> dict[str, Any]:
    """
    Operator A_comp: Mehrschichtschale → komprimierte (2,3)-Einfachschale.
    """
    a_terms = list(a_terms or [])
    r_sq = multilayer_r_squared(n, m, a_terms)
    r = math.sqrt(r_sq)
    k = count_multilayer_terms(n, m, a_terms)
    a_classes = sorted({a_term_class(q) for q, _ in a_terms if a_term_class(q)})

    if method == "greedy_log":
        n_p, m_p, r_tilde_sq, defect, _ = absorb_greedy_log(n, m, a_terms)
        search_method = "greedy_log"
    else:
        n_p, m_p, r_tilde_sq, defect = absorb_exact_search(r_sq, exp_max, exp_max)
        search_method = "exact_search"

    r_tilde = math.sqrt(r_tilde_sq)
    grade = compression_grade(k)

    per_term: list[dict[str, Any]] = []
    for q, alpha in a_terms:
        tp = q**alpha
        tn, tm, _, td = absorb_exact_search(tp, min(exp_max, 10), min(exp_max, 10))
        per_term.append(
            {
                "term": f"{q}^{alpha}",
                "value": int(tp),
                "eabc_class": a_term_class(q),
                "absorbed_to": f"2^{tn}+3^{tm}",
                "term_defect": int(td),
            }
        )

    return {
        "input": {
            "n": n,
            "m": m,
            "a_terms": [(int(q), int(a)) for q, a in a_terms],
            "formula": _formula_str(n, m, a_terms),
        },
        "r_squared": int(r_sq),
        "R": round(r, 6),
        "k_multi": k,
        "a_eabc_classes": a_classes,
        "a_comp": {
            "method": search_method,
            "n_prime": int(n_p),
            "m_prime": int(m_p),
            "r_tilde_squared": int(r_tilde_sq),
            "R_tilde": round(r_tilde, 6),
            "absorption_defect": int(defect),
            "compression_grade": round(grade, 4),
            "formula": f"2^{n_p}+3^{m_p}",
        },
        "per_a_term_absorption": per_term,
    }


def _formula_str(n: int, m: int, a_terms: list[tuple[int, int]]) -> str:
    parts = []
    if n > 0:
        parts.append(f"2^{n}")
    if m > 0:
        parts.append(f"3^{m}")
    for q, alpha in a_terms:
        parts.append(f"{q}^{alpha}")
    return "+".join(parts) if parts else "0"


def primes_in_band(r: float, primes: list[int], epsilon: float) -> list[int]:
    lo, hi = r - epsilon, r + epsilon
    return [p for p in primes if lo <= p <= hi]


def eabc_band_observables(
    r: float,
    primes: list[int],
    epsilon: float,
    n_max: int,
    m_max: int,
) -> dict[str, Any]:
    """δ(p) und D_AB im ε-Band um Radius R."""
    band = primes_in_band(r, primes, epsilon)
    counts = np.zeros(4, dtype=np.int64)
    deltas: list[float] = []
    d_ab_vals: list[float] = []
    for p in band:
        lab = class_of(p)
        if lab is None:
            continue
        counts[LABEL_TO_IDX[lab]] += 1
        delta, _, _, _ = defect_delta(p, n_max, m_max)
        deltas.append(delta)
        obs = charge_observables(np.array([0, 1 if lab == "A" else 0, 1 if lab == "B" else 0, 1 if lab == "C" else 0]))
        d_ab_vals.append(obs["D_AB"])
    charge = charge_observables(counts) if counts.sum() > 0 else {
        "N_E": 0.0, "N_A": 0.0, "N_B": 0.0, "N_C": 0.0,
        "D_AB": 0.0, "D_AC": 0.0, "D_E": 0.0,
    }
    return {
        "n_primes": int(counts.sum()),
        "mean_delta": round(float(np.mean(deltas)), 4) if deltas else None,
        "mean_D_AB": round(float(np.mean(d_ab_vals)), 4) if d_ab_vals else None,
        "charge": {k: round(v, 4) if isinstance(v, float) else v for k, v in charge.items()},
        "primes_sample": band[:6],
    }


def compare_absorption_observables(
    shell: dict[str, Any],
    primes: list[int],
    epsilon: float,
    n_max: int,
    m_max: int,
) -> dict[str, Any]:
    """Vergleich δ(p), D_AB vor/nach Absorption im Prim-Umfeld."""
    r = float(shell["R"])
    r_t = float(shell["a_comp"]["R_tilde"])
    before = eabc_band_observables(r, primes, epsilon, n_max, m_max)
    after = eabc_band_observables(r_t, primes, epsilon, n_max, m_max)
    delta_shift = None
    if before["mean_delta"] is not None and after["mean_delta"] is not None:
        delta_shift = round(after["mean_delta"] - before["mean_delta"], 4)
    d_ab_shift = None
    if before["charge"]["D_AB"] != 0 or after["charge"]["D_AB"] != 0:
        d_ab_shift = round(after["charge"]["D_AB"] - before["charge"]["D_AB"], 4)
    return {
        "epsilon": epsilon,
        "before_absorption": before,
        "after_absorption": after,
        "delta_mean_shift": delta_shift,
        "D_AB_shift": d_ab_shift,
    }


def generate_multilayer_shells(
    n_max: int,
    m_max: int,
    l_max: int,
    *,
    a_axes: tuple[int, ...] = (5, 17),
) -> list[dict[str, Any]]:
    """Beispielfamilie R² = 2^n + 3^m + Σ q^ℓ mit q ≡ 5 mod 12."""
    shells: list[dict[str, Any]] = []
    for n in range(1, n_max + 1):
        for m in range(1, m_max + 1):
            for q in a_axes:
                for ell in range(1, l_max + 1):
                    a_terms = [(q, ell)]
                    try:
                        r_sq = multilayer_r_squared(n, m, a_terms)
                    except ValueError:
                        continue
                    shells.append(
                        {
                            "n": n,
                            "m": m,
                            "a_terms": a_terms,
                            "r_squared": r_sq,
                            "R": float(math.sqrt(r_sq)),
                        }
                    )
    return shells


def run_a_absorption_analysis(
    *,
    n_max: int = 3,
    m_max: int = 3,
    l_max: int = 2,
    exp_max: int = 10,
    prime_limit: int = 5_000,
    epsilon: float = 2.0,
    method: str = "exact",
    seed: int = 42,
) -> dict[str, Any]:
    """Smoke-Analyse: A-Absorption + δ/D_AB-Vergleich."""
    _ = seed  # reserviert für Perm-Null-Erweiterung
    primes = [p for p in sieve_primes(prime_limit) if p > 3]
    shells_raw = generate_multilayer_shells(n_max, m_max, l_max)
    results: list[dict[str, Any]] = []
    defects: list[int] = []
    grades: list[float] = []
    delta_shifts: list[float] = []

    for sh in shells_raw:
        comp = apply_a_comp(
            sh["n"], sh["m"], sh["a_terms"],
            method=method, exp_max=exp_max,
        )
        comp["observables"] = compare_absorption_observables(
            comp, primes, epsilon, exp_max, exp_max,
        )
        results.append(comp)
        defects.append(comp["a_comp"]["absorption_defect"])
        grades.append(comp["a_comp"]["compression_grade"])
        ds = comp["observables"]["delta_mean_shift"]
        if ds is not None:
            delta_shifts.append(ds)

    # Beispiel 2^n + 3^m + 5^ℓ explizit
    examples = [
        apply_a_comp(2, 2, [(5, 1)], method=method, exp_max=exp_max),
        apply_a_comp(3, 1, [(5, 2)], method=method, exp_max=exp_max),
        apply_a_comp(1, 2, [(5, 1), (17, 1)], method=method, exp_max=exp_max),
    ]

    return {
        "section": "4.8.26",
        "method": method,
        "n_shells": len(results),
        "prime_limit": prime_limit,
        "epsilon": epsilon,
        "summary": {
            "absorption_defect_mean": round(float(np.mean(defects)), 2) if defects else None,
            "absorption_defect_max": int(max(defects)) if defects else None,
            "compression_grade_mean": round(float(np.mean(grades)), 4) if grades else None,
            "delta_mean_shift_mean": round(float(np.mean(delta_shifts)), 4) if delta_shifts else None,
            "all_a_terms_in_A_class": all(
                a_term_class(q) == "A" for r in results for q, _ in r["input"]["a_terms"]
            ),
        },
        "examples": examples,
        "shells": results[:24],
    }


def print_smoke_summary(report: dict[str, Any]) -> None:
    s = report["summary"]
    print(f"§{report['section']} A-Absorption — {report['n_shells']} Schalen ({report['method']})")
    print(f"  Defekt: mean={s['absorption_defect_mean']}, max={s['absorption_defect_max']}")
    print(f"  Kompressionsgrad mean={s['compression_grade_mean']}")
    print(f"  δ-Verschiebung mean={s['delta_mean_shift_mean']}")
    print(f"  A-Terme EABC-Klasse A: {s['all_a_terms_in_A_class']}")
    for ex in report["examples"][:2]:
        inp = ex["input"]["formula"]
        out = ex["a_comp"]["formula"]
        d = ex["a_comp"]["absorption_defect"]
        print(f"  Beispiel: {inp} → {out}  (Defekt={d})")


def main() -> None:
    parser = argparse.ArgumentParser(description="A-Absorption / A-Transfer (§4.8.26)")
    parser.add_argument("--smoke", action="store_true", help="begrenzter Smoke-Lauf")
    parser.add_argument("--n-max", type=int, default=3)
    parser.add_argument("--m-max", type=int, default=3)
    parser.add_argument("--l-max", type=int, default=2)
    parser.add_argument("--exp-max", type=int, default=10)
    parser.add_argument("--prime-limit", type=int, default=5_000)
    parser.add_argument("--epsilon", type=float, default=2.0)
    parser.add_argument("--method", choices=("exact", "greedy_log"), default="exact")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if args.smoke:
        report = run_a_absorption_analysis(
            n_max=min(args.n_max, 3),
            m_max=min(args.m_max, 3),
            l_max=min(args.l_max, 2),
            exp_max=min(args.exp_max, 10),
            prime_limit=min(args.prime_limit, 5_000),
            epsilon=args.epsilon,
            method=args.method,
            seed=args.seed,
        )
    else:
        report = run_a_absorption_analysis(
            n_max=args.n_max,
            m_max=args.m_max,
            l_max=args.l_max,
            exp_max=args.exp_max,
            prime_limit=args.prime_limit,
            epsilon=args.epsilon,
            method=args.method,
            seed=args.seed,
        )

    print_smoke_summary(report)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"→ {args.output}")


if __name__ == "__main__":
    main()
