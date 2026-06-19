#!/usr/bin/env python3
"""
EABC-Primverteilung auf Schalen R_{n,m} (§4.8.21)

Dreistufige Schalenhierarchie:
  Ebene 1 — glatte Geometrie: R_{n,m} = sqrt(2^n + 3^m), Log-Gitter (n,m)
  Ebene 2 — Primprojektion: P(p) = (R_{n,m}, v(p)), v(p) ∈ {E,A,B,C} via p mod 12
  Ebene 3 — Impuls/Masse: M(N) = Σ e_i v(p_i) log(p_i) (allgemeine N; hier Prim-Stub)

Schalen-Zuordnung (dokumentierte Wahl):
  factor (Default, kanonisch §4.8.20):
      N = p als Prim → α,β aus 2^α 3^β | p.
      Für p > 3 gilt α = β = 0 (keine Faktoren 2,3) → Schale S_{0,0}, R_{0,0} = √5.
  log-grid (empirische Stratifikation):
      (n,m) minimiert |log p − log R_{n,m}| auf Gitter [0,n_max]×[0,m_max].

Nullmodell:
  perm — Permutation der EABC-Labels bei festen Marginalen (global).
  uniform — erwartete Gleichverteilung 1/4 je Klasse pro Schale.

CLI: --limit, --output json, --shell-mode {factor,log-grid}
"""

from __future__ import annotations

import argparse
import json
import math
from math import isqrt
from pathlib import Path
from typing import Any

import numpy as np

LABELS = ("E", "A", "B", "C")
LABEL_TO_IDX = {lab: i for i, lab in enumerate(LABELS)}
RESIDUE_MAP = {1: "E", 5: "A", 7: "B", 11: "C"}


def sieve_primes(limit: int) -> list[int]:
    if limit < 2:
        return []
    flags = bytearray([1]) * (limit + 1)
    flags[0] = flags[1] = 0
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p * p :: p] = bytearray(len(flags[p * p :: p]))
    return [i for i, ok in enumerate(flags) if ok]


def class_of(p: int) -> str | None:
    return RESIDUE_MAP.get(p % 12)


def shell_radius_squared(n: int, m: int) -> int:
    if n < 0 or m < 0:
        raise ValueError("n,m ≥ 0 erforderlich")
    return 2**n + 3**m


def shell_radius(n: int, m: int) -> float:
    return math.sqrt(shell_radius_squared(n, m))


def smooth_shell_exponents(n: int) -> tuple[int, int]:
    """α, β aus N = 2^α 3^β × m (glatte (2,3)-Schale)."""
    alpha = 0
    while n % 2 == 0:
        alpha += 1
        n //= 2
    beta = 0
    while n % 3 == 0:
        beta += 1
        n //= 3
    return alpha, beta


def assign_shell_factor(p: int) -> tuple[int, int]:
    """Kanonisch: Schale aus Faktorisierung von N = p."""
    return smooth_shell_exponents(p)


def assign_shell_log_grid(p: int, n_max: int, m_max: int) -> tuple[int, int]:
    """Empirisch: (n,m) mit minimalem |log p − log R_{n,m}|."""
    log_p = math.log(p)
    best = (0, 0)
    best_dist = float("inf")
    for n in range(n_max + 1):
        for m in range(m_max + 1):
            r = shell_radius(n, m)
            if r <= 0:
                continue
            d = abs(log_p - math.log(r))
            if d < best_dist:
                best_dist = d
                best = (n, m)
    return best


def primes_above_three(limit: int) -> list[int]:
    return [p for p in sieve_primes(limit) if p > 3]


def bucket_key(n: int, m: int) -> str:
    return f"{n},{m}"


def counts_to_dict(counts: np.ndarray) -> dict[str, int]:
    return {lab: int(counts[i]) for i, lab in enumerate(LABELS)}


def chi2_vs_uniform(counts: np.ndarray) -> float:
    n = int(counts.sum())
    if n == 0:
        return float("nan")
    expected = n / 4.0
    return float(np.sum((counts - expected) ** 2 / expected))


def perm_p_value_chi2(observed: np.ndarray, n_perm: int, rng: np.random.Generator) -> float:
    """Einseitig: Anteil Permutationen mit χ² ≥ χ²_obs (Marginalen erhalten)."""
    n = int(observed.sum())
    if n == 0:
        return float("nan")
    labels = np.repeat(np.arange(4, dtype=np.int64), observed.astype(np.int64))
    chi2_obs = chi2_vs_uniform(observed)
    ge = 0
    for _ in range(n_perm):
        rng.shuffle(labels)
        perm_counts = np.bincount(labels, minlength=4)
        if chi2_vs_uniform(perm_counts) >= chi2_obs:
            ge += 1
    return (ge + 1) / (n_perm + 1)


def analyze_distribution(
    limit: int = 10_000,
    shell_mode: str = "factor",
    n_max: int = 8,
    m_max: int = 8,
    n_perm: int = 499,
    seed: int = 42,
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    primes = primes_above_three(limit)

    buckets: dict[str, dict[str, Any]] = {}
    global_counts = np.zeros(4, dtype=np.int64)

    for p in primes:
        v = class_of(p)
        if v is None:
            continue
        if shell_mode == "factor":
            n, m = assign_shell_factor(p)
        elif shell_mode == "log-grid":
            n, m = assign_shell_log_grid(p, n_max, m_max)
        else:
            raise ValueError(f"unbekannter shell_mode: {shell_mode}")

        key = bucket_key(n, m)
        if key not in buckets:
            buckets[key] = {
                "n": n,
                "m": m,
                "R_nm": shell_radius(n, m),
                "counts": np.zeros(4, dtype=np.int64),
                "primes_sample": [],
            }
        buckets[key]["counts"][LABEL_TO_IDX[v]] += 1
        global_counts[LABEL_TO_IDX[v]] += 1
        if len(buckets[key]["primes_sample"]) < 5:
            buckets[key]["primes_sample"].append(int(p))

    bucket_rows: list[dict[str, Any]] = []
    for key in sorted(buckets, key=lambda k: (buckets[k]["n"], buckets[k]["m"])):
        b = buckets[key]
        counts = b["counts"]
        n_total = int(counts.sum())
        chi2 = chi2_vs_uniform(counts)
        p_perm = perm_p_value_chi2(counts, n_perm, rng) if n_total >= 8 else None
        bucket_rows.append(
            {
                "shell": key,
                "n": b["n"],
                "m": b["m"],
                "R_nm": round(b["R_nm"], 6),
                "n_primes": n_total,
                "counts": counts_to_dict(counts),
                "fractions": {
                    lab: round(int(counts[i]) / n_total, 4) if n_total else 0.0
                    for i, lab in enumerate(LABELS)
                },
                "chi2_vs_uniform": round(chi2, 4) if n_total else None,
                "perm_p_chi2": round(p_perm, 4) if p_perm is not None else None,
                "primes_sample": b["primes_sample"],
            }
        )

    n_all = int(global_counts.sum())
    global_chi2 = chi2_vs_uniform(global_counts)
    global_p = perm_p_value_chi2(global_counts, n_perm, rng)

    shell_assignment_note = (
        "Für p>3 und mode=factor: α=β=0, alle Primzahlen auf S_{0,0} "
        "(R_{0,0}=√5); Schalenstratifizierung trivial, globaler EABC-Test."
        if shell_mode == "factor"
        else f"log-grid auf [0,{n_max}]×[0,{m_max}]: Primzahlen nach log-Abstand zu R_{{n,m}}."
    )

    return {
        "section": "4.8.21",
        "script": "eabc_shell_prime_distribution.py",
        "limit": limit,
        "shell_mode": shell_mode,
        "shell_assignment_note": shell_assignment_note,
        "n_max": n_max if shell_mode == "log-grid" else None,
        "m_max": m_max if shell_mode == "log-grid" else None,
        "n_primes": n_all,
        "n_buckets": len(bucket_rows),
        "global": {
            "counts": counts_to_dict(global_counts),
            "fractions": {
                lab: round(int(global_counts[i]) / n_all, 4) for i, lab in enumerate(LABELS)
            },
            "chi2_vs_uniform": round(global_chi2, 4),
            "perm_p_chi2": round(global_p, 4),
            "expected_fraction": 0.25,
        },
        "buckets": bucket_rows,
        "null_model": "perm (Marginalen erhalten, χ²-Statistik)",
        "n_perm": n_perm,
        "seed": seed,
    }


def format_summary(result: dict[str, Any]) -> str:
    g = result["global"]
    fr = g["fractions"]
    return (
        f"n={result['n_primes']} Primzahlen>p≤{result['limit']}, "
        f"mode={result['shell_mode']}, buckets={result['n_buckets']}; "
        f"EABC={fr['E']:.3f}/{fr['A']:.3f}/{fr['B']:.3f}/{fr['C']:.3f}, "
        f"χ²={g['chi2_vs_uniform']}, perm_p={g['perm_p_chi2']}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="EABC-Primverteilung auf Schalen R_{n,m} (§4.8.21)"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=50_000,
        help="Obergrenze für Primzahl-Sieb (default: 50000)",
    )
    parser.add_argument(
        "--shell-mode",
        choices=("factor", "log-grid"),
        default="factor",
        help="Schalen-Zuordnung: factor (kanonisch) oder log-grid (Stratifizierung)",
    )
    parser.add_argument("--n-max", type=int, default=8, help="log-grid: max n (default: 8)")
    parser.add_argument("--m-max", type=int, default=8, help="log-grid: max m (default: 8)")
    parser.add_argument("--n-perm", type=int, default=499, help="Permutationen (default: 499)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=str, default="", help="JSON-Ausgabedatei")
    args = parser.parse_args()

    result = analyze_distribution(
        limit=args.limit,
        shell_mode=args.shell_mode,
        n_max=args.n_max,
        m_max=args.m_max,
        n_perm=args.n_perm,
        seed=args.seed,
    )

    print(format_summary(result))

    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", encoding="utf-8") as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        print(f"JSON → {out}")


if __name__ == "__main__":
    main()
