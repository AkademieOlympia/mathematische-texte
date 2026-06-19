#!/usr/bin/env python3
"""
EABC-Primverteilung auf Schalen R_{n,m} (§4.8.21–4.8.23)

Dreistufige Schalenhierarchie (§4.8.21):
  Ebene 1 — glatte Geometrie: R_{n,m} = sqrt(2^n + 3^m), Log-Gitter (n,m)
  Ebene 2 — Primprojektion: P(p) = (R_{n,m}, v(p)), v(p) ∈ {E,A,B,C} via p mod 12
  Ebene 3 — Impuls/Masse: M(N) = Σ e_i v(p_i) log(p_i) (allgemeine N; hier Prim-Stub)

Schalen-Defekte (§4.8.22):
  δ(p) = min_{n,m} |p − R_{n,m}|
  Band-Zuordnung: nearest R_{n,m}, ε-Band |p−R_{n,m}|<ε, oder Intervall R_k ≤ p < R_{k+1}
  Observablen: N_A, N_B, N_C, N_E; D_AB, D_AC, D_E = N_E − (N_A+N_B+N_C)/3

Defektprogramm (§4.8.23):
  log N = S(N) + P(N),  S(N)=α log 2 + β log 3,  P(N)=Σ e_i log p_i
  V₄-Charaktere χ_A, χ_B, χ_C auf {E,A,B,C} ≅ C₂×C₂
  δ-Buckets: N_A(δ), …; asymptotisches D_AB(X) vs. X

Nullmodell:
  perm — Permutation der EABC-Labels bei festen Marginalen (global).
  uniform — erwartete Gleichverteilung 1/4 je Klasse pro Schale.

CLI: --mode {shell,delta-program}, --limit, --output json
      --shell-mode {factor,log-grid,defect}; --epsilon, --band-mode (defect)
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

# V₄ ≅ C₂×C₂: E Identität; χ_A, χ_B, χ_C nichttriviale Charaktere (§4.8.23)
V4_CHARACTERS: dict[str, dict[str, int]] = {
    "E": {"chi_A": 1, "chi_B": 1, "chi_C": 1},
    "A": {"chi_A": -1, "chi_B": 1, "chi_C": -1},
    "B": {"chi_A": 1, "chi_B": -1, "chi_C": -1},
    "C": {"chi_A": -1, "chi_B": -1, "chi_C": 1},
}


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


def shell_grid(n_max: int, m_max: int) -> list[tuple[int, int, float]]:
    """Alle Schalen (n,m,R) auf endlichem Gitter, sortiert nach R."""
    cells = [(n, m, shell_radius(n, m)) for n in range(n_max + 1) for m in range(m_max + 1)]
    return sorted(cells, key=lambda t: t[2])


def defect_delta(p: int, n_max: int, m_max: int) -> tuple[float, int, int, float]:
    """δ(p) = min_{n,m} |p − R_{n,m}|; liefert (δ, n*, m*, R*)."""
    best_delta = float("inf")
    best = (0, 0, shell_radius(0, 0))
    pf = float(p)
    for n in range(n_max + 1):
        for m in range(m_max + 1):
            r = shell_radius(n, m)
            d = abs(pf - r)
            if d < best_delta:
                best_delta = d
                best = (n, m, r)
    return best_delta, best[0], best[1], best[2]


def assign_shell_defect_nearest(p: int, n_max: int, m_max: int) -> tuple[int, int, float, float]:
    """Nächste Schale nach |p − R_{n,m}|."""
    delta, n, m, r = defect_delta(p, n_max, m_max)
    return n, m, r, delta


def assign_shell_defect_interval(
    p: int, grid: list[tuple[int, int, float]]
) -> tuple[int, int, float, float] | None:
    """Intervallband R_k ≤ p < R_{k+1}; unteres Schalenpaar (n_k,m_k)."""
    pf = float(p)
    for i in range(len(grid) - 1):
        n0, m0, r0 = grid[i]
        _n1, _m1, r1 = grid[i + 1]
        if r0 <= pf < r1:
            return n0, m0, r0, abs(pf - r0)
    return None


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


def factorize_radical(n: int) -> list[tuple[int, int]]:
    """Primfaktorisierung des Radikalanteils (nach Abzug von 2 und 3)."""
    if n < 1:
        raise ValueError("n ≥ 1 erforderlich")
    m = n
    while m % 2 == 0:
        m //= 2
    while m % 3 == 0:
        m //= 3
    if m <= 1:
        return []
    factors: list[tuple[int, int]] = []
    d = 5
    while d * d <= m:
        if m % d == 0:
            e = 0
            while m % d == 0:
                e += 1
                m //= d
            factors.append((d, e))
        d += 2
    if m > 1:
        factors.append((m, 1))
    return factors


def log_decomposition(n: int) -> dict[str, Any]:
    """log N = S(N) + P(N) mit S(N)=α log 2 + β log 3, P(N)=Σ e_i log p_i."""
    alpha, beta = smooth_shell_exponents(n)
    prime_factors = factorize_radical(n)
    s_part = alpha * math.log(2) + beta * math.log(3)
    p_part = sum(e * math.log(p) for p, e in prime_factors)
    log_n = math.log(n)
    return {
        "N": n,
        "alpha": alpha,
        "beta": beta,
        "S": s_part,
        "P": p_part,
        "log_N": log_n,
        "residual": log_n - s_part - p_part,
        "prime_factors": [(int(p), int(e)) for p, e in prime_factors],
    }


def v4_characters(label: str) -> dict[str, int]:
    """Charakterwerte χ_A, χ_B, χ_C für EABC-Label."""
    if label not in V4_CHARACTERS:
        raise ValueError(f"unbekanntes Label: {label}")
    return dict(V4_CHARACTERS[label])


def v4_character_vector_prime(p: int) -> dict[str, Any]:
    """Charaktervektor (χ_A, χ_B, χ_C) für Primzahl p>3."""
    lab = class_of(p)
    if lab is None:
        raise ValueError(f"p={p} hat keine EABC-Klasse mod 12")
    chi = v4_characters(lab)
    return {"p": int(p), "label": lab, **chi}


def M_representation(n: int) -> dict[str, Any]:
    """M(N) als V₄-Gewicht und Charakter-Summen Σ e_i χ_k(v(p_i)) log p_i."""
    dec = log_decomposition(n)
    m_weight = {lab: 0.0 for lab in LABELS}
    m_chi = {k: 0.0 for k in ("chi_A", "chi_B", "chi_C")}
    for p, e in dec["prime_factors"]:
        lab = class_of(p)
        if lab is None:
            continue
        w = e * math.log(p)
        m_weight[lab] += w
        for ck, cv in v4_characters(lab).items():
            m_chi[ck] += e * cv * math.log(p)
    return {
        "N": n,
        "S": dec["S"],
        "P": dec["P"],
        "M_weight": {k: round(v, 6) for k, v in m_weight.items()},
        "M_chi": {k: round(v, 6) for k, v in m_chi.items()},
    }


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


def charge_observables(counts: np.ndarray) -> dict[str, float]:
    """D_AB, D_AC, D_E aus EABC-Zählvektor."""
    n_e, n_a, n_b, n_c = (int(counts[i]) for i in range(4))
    abc_third = (n_a + n_b + n_c) / 3.0
    return {
        "N_E": n_e,
        "N_A": n_a,
        "N_B": n_b,
        "N_C": n_c,
        "D_AB": float(n_a - n_b),
        "D_AC": float(n_a - n_c),
        "D_E": float(n_e - abc_third),
    }


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


def perm_p_value_D_abs(
    observed: np.ndarray, stat_key: str, n_perm: int, rng: np.random.Generator
) -> float:
    """Einseitig: |D| unter Label-Permutation vs. beobachtet."""
    n = int(observed.sum())
    if n == 0:
        return float("nan")
    obs = charge_observables(observed)
    stat_obs = abs(obs[stat_key])
    labels = np.repeat(np.arange(4, dtype=np.int64), observed.astype(np.int64))
    ge = 0
    for _ in range(n_perm):
        rng.shuffle(labels)
        perm_counts = np.bincount(labels, minlength=4)
        stat_p = abs(charge_observables(perm_counts)[stat_key])
        if stat_p >= stat_obs:
            ge += 1
    return (ge + 1) / (n_perm + 1)


def analyze_distribution(
    limit: int = 10_000,
    shell_mode: str = "factor",
    n_max: int = 8,
    m_max: int = 8,
    n_perm: int = 499,
    seed: int = 42,
    epsilon: float | None = None,
    band_mode: str = "nearest",
    null_model: str = "perm",
) -> dict[str, Any]:
    if shell_mode == "defect":
        return analyze_defect_distribution(
            limit=limit,
            n_max=n_max,
            m_max=m_max,
            n_perm=n_perm,
            seed=seed,
            epsilon=epsilon,
            band_mode=band_mode,
            null_model=null_model,
        )

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


def analyze_defect_distribution(
    limit: int = 10_000,
    n_max: int = 8,
    m_max: int = 8,
    n_perm: int = 499,
    seed: int = 42,
    epsilon: float | None = None,
    band_mode: str = "nearest",
    null_model: str = "perm",
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    primes = primes_above_three(limit)
    grid = shell_grid(n_max, m_max)

    if band_mode == "epsilon" and epsilon is None:
        raise ValueError("band_mode=epsilon erfordert --epsilon > 0")

    buckets: dict[str, dict[str, Any]] = {}
    global_counts = np.zeros(4, dtype=np.int64)
    delta_samples: list[float] = []
    skipped_interval = 0

    for p in primes:
        v = class_of(p)
        if v is None:
            continue

        if band_mode == "interval":
            assigned = assign_shell_defect_interval(p, grid)
            if assigned is None:
                skipped_interval += 1
                continue
            n, m, r_nm, delta = assigned
        else:
            n, m, r_nm, delta = assign_shell_defect_nearest(p, n_max, m_max)
            if band_mode == "epsilon":
                assert epsilon is not None
                if delta >= epsilon:
                    continue

        delta_samples.append(delta)
        key = bucket_key(n, m)
        if key not in buckets:
            buckets[key] = {
                "n": n,
                "m": m,
                "R_nm": r_nm,
                "counts": np.zeros(4, dtype=np.int64),
                "delta_sum": 0.0,
                "primes_sample": [],
            }
        buckets[key]["counts"][LABEL_TO_IDX[v]] += 1
        buckets[key]["delta_sum"] += delta
        global_counts[LABEL_TO_IDX[v]] += 1
        if len(buckets[key]["primes_sample"]) < 5:
            buckets[key]["primes_sample"].append(int(p))

    bucket_rows: list[dict[str, Any]] = []
    for key in sorted(buckets, key=lambda k: (buckets[k]["n"], buckets[k]["m"])):
        b = buckets[key]
        counts = b["counts"]
        n_total = int(counts.sum())
        obs = charge_observables(counts)
        chi2 = chi2_vs_uniform(counts)
        p_perm = perm_p_value_chi2(counts, n_perm, rng) if n_total >= 8 else None
        row: dict[str, Any] = {
            "shell": key,
            "n": b["n"],
            "m": b["m"],
            "R_nm": round(b["R_nm"], 6),
            "n_primes": n_total,
            "counts": counts_to_dict(counts),
            "charge": {k: round(v, 4) if isinstance(v, float) else v for k, v in obs.items()},
            "mean_delta": round(b["delta_sum"] / n_total, 6) if n_total else None,
            "fractions": {
                lab: round(int(counts[i]) / n_total, 4) if n_total else 0.0
                for i, lab in enumerate(LABELS)
            },
            "chi2_vs_uniform": round(chi2, 4) if n_total else None,
            "perm_p_chi2": round(p_perm, 4) if p_perm is not None else None,
            "primes_sample": b["primes_sample"],
        }
        if null_model == "perm" and n_total >= 4:
            row["perm_p_D_AB"] = round(perm_p_value_D_abs(counts, "D_AB", n_perm, rng), 4)
            row["perm_p_D_AC"] = round(perm_p_value_D_abs(counts, "D_AC", n_perm, rng), 4)
            row["perm_p_D_E"] = round(perm_p_value_D_abs(counts, "D_E", n_perm, rng), 4)
        bucket_rows.append(row)

    n_all = int(global_counts.sum())
    global_obs = charge_observables(global_counts)
    global_chi2 = chi2_vs_uniform(global_counts)
    global_p = perm_p_value_chi2(global_counts, n_perm, rng)

    band_note = {
        "nearest": f"δ(p)-Minimierer auf Gitter [0,{n_max}]×[0,{m_max}]",
        "epsilon": f"ε-Band: |p−R_{{n,m}}|<{epsilon} um nächste Schale",
        "interval": f"Intervallband R_k ≤ p < R_{{k+1}} auf sortiertem Gitter",
    }[band_mode]

    result: dict[str, Any] = {
        "section": "4.8.22",
        "script": "eabc_shell_prime_distribution.py",
        "limit": limit,
        "shell_mode": "defect",
        "band_mode": band_mode,
        "epsilon": epsilon,
        "shell_assignment_note": band_note,
        "n_max": n_max,
        "m_max": m_max,
        "n_primes": n_all,
        "n_primes_input": len(primes),
        "skipped_interval": skipped_interval if band_mode == "interval" else 0,
        "n_buckets": len(bucket_rows),
        "delta_stats": {
            "mean": round(float(np.mean(delta_samples)), 6) if delta_samples else None,
            "median": round(float(np.median(delta_samples)), 6) if delta_samples else None,
            "max": round(float(np.max(delta_samples)), 6) if delta_samples else None,
        },
        "global": {
            "counts": counts_to_dict(global_counts),
            "fractions": {
                lab: round(int(global_counts[i]) / n_all, 4) for i, lab in enumerate(LABELS)
            },
            "charge": {
                k: round(v, 4) if isinstance(v, float) else v for k, v in global_obs.items()
            },
            "chi2_vs_uniform": round(global_chi2, 4),
            "perm_p_chi2": round(global_p, 4),
            "expected_fraction": 0.25,
        },
        "buckets": bucket_rows,
        "null_model": null_model,
        "n_perm": n_perm,
        "seed": seed,
    }

    if null_model == "perm":
        result["global"]["perm_p_D_AB"] = round(
            perm_p_value_D_abs(global_counts, "D_AB", n_perm, rng), 4
        )
        result["global"]["perm_p_D_AC"] = round(
            perm_p_value_D_abs(global_counts, "D_AC", n_perm, rng), 4
        )
        result["global"]["perm_p_D_E"] = round(
            perm_p_value_D_abs(global_counts, "D_E", n_perm, rng), 4
        )
    elif null_model == "uniform":
        result["global"]["uniform_note"] = (
            "Erwartung: N_i=n/4, D_AB=D_AC=D_E=0; χ²-Test vs. Gleichverteilung"
        )

    return result


def delta_bin_index(delta: float, edges: list[float]) -> int:
    for i in range(len(edges) - 1):
        if edges[i] <= delta < edges[i + 1]:
            return i
    return len(edges) - 2


def asymptotic_D_AB_series(
    limit: int, checkpoints: int = 12
) -> list[dict[str, Any]]:
    """D_AB(X) und D_AB(X)/X für wachsende Schwellen X ≤ limit."""
    primes = primes_above_three(limit)
    if not primes:
        return []
    xs = np.unique(
        np.round(np.geomspace(max(100, primes[0]), limit, num=checkpoints)).astype(int)
    )
    rows: list[dict[str, Any]] = []
    idx = 0
    counts = np.zeros(4, dtype=np.int64)
    for x in xs:
        while idx < len(primes) and primes[idx] <= x:
            lab = class_of(primes[idx])
            if lab is not None:
                counts[LABEL_TO_IDX[lab]] += 1
            idx += 1
        obs = charge_observables(counts)
        n_tot = int(counts.sum())
        d_ab = obs["D_AB"]
        ratio = d_ab / float(x) if x > 0 else float("nan")
        rows.append(
            {
                "X": int(x),
                "n_primes": n_tot,
                "D_AB": round(d_ab, 4),
                "D_AC": round(obs["D_AC"], 4),
                "D_E": round(obs["D_E"], 4),
                "D_AB_over_X": round(ratio, 8),
            }
        )
    return rows


def assess_o_X(asymptotic_rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Heuristik: |D_AB(X)|/X fallend über letzte Hälfte der Checkpoints."""
    if len(asymptotic_rows) < 3:
        return {"o_X_likely": None, "note": "zu wenige Checkpoints"}
    ratios = [abs(r["D_AB_over_X"]) for r in asymptotic_rows]
    half = len(ratios) // 2
    early = float(np.mean(ratios[:half]))
    late = float(np.mean(ratios[half:]))
    decreasing = late < early
    return {
        "o_X_likely": decreasing,
        "mean_abs_ratio_early": round(early, 8),
        "mean_abs_ratio_late": round(late, 8),
        "note": (
            "|D_AB|/X im Mittel fallend (heuristisch o(X)-kompatibel)"
            if decreasing
            else "|D_AB|/X nicht klar fallend — o(X) offen"
        ),
    }


def analyze_delta_program(
    limit: int = 100_000,
    n_max: int = 8,
    m_max: int = 8,
    n_delta_bins: int = 8,
    n_perm: int = 499,
    seed: int = 42,
    asymptotic_checkpoints: int = 12,
) -> dict[str, Any]:
    """Defektprogramm §4.8.23: S/P-Split, δ-Buckets, V₄-Charaktere, D_AB(X)."""
    rng = np.random.default_rng(seed)
    primes = primes_above_three(limit)

    log_samples = [log_decomposition(n) for n in (72, 210, 2310, 30030)]
    for row in log_samples:
        row["S"] = round(row["S"], 6)
        row["P"] = round(row["P"], 6)
        row["log_N"] = round(row["log_N"], 6)
        row["residual"] = round(row["residual"], 12)

    prime_char_examples = [v4_character_vector_prime(p) for p in primes[:8]]
    m_examples = [M_representation(n) for n in (210, 30030)]

    deltas: list[float] = []
    prime_records: list[tuple[int, str, float]] = []
    for p in primes:
        lab = class_of(p)
        if lab is None:
            continue
        delta, _, _, _ = defect_delta(p, n_max, m_max)
        deltas.append(delta)
        prime_records.append((p, lab, delta))

    if not deltas:
        edges = [0.0, 1.0]
    else:
        qs = np.linspace(0.0, 1.0, n_delta_bins + 1)
        edges = [float(np.quantile(deltas, q)) for q in qs]
        edges[-1] = edges[-1] * (1.0 + 1e-12) + 1e-9

    delta_buckets: list[dict[str, Any]] = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        counts = np.zeros(4, dtype=np.int64)
        bin_deltas: list[float] = []
        for p, lab, delta in prime_records:
            if i < len(edges) - 2:
                in_bin = lo <= delta < hi
            else:
                in_bin = lo <= delta <= hi
            if not in_bin:
                continue
            counts[LABEL_TO_IDX[lab]] += 1
            bin_deltas.append(delta)
        n_total = int(counts.sum())
        if n_total == 0:
            continue
        obs = charge_observables(counts)
        row: dict[str, Any] = {
            "bin": i,
            "delta_lo": round(lo, 6),
            "delta_hi": round(hi, 6),
            "n_primes": n_total,
            "N_E": obs["N_E"],
            "N_A": obs["N_A"],
            "N_B": obs["N_B"],
            "N_C": obs["N_C"],
            "counts": counts_to_dict(counts),
            "charge": {k: round(v, 4) if isinstance(v, float) else v for k, v in obs.items()},
            "mean_delta": round(float(np.mean(bin_deltas)), 6),
            "fractions": {
                lab: round(int(counts[j]) / n_total, 4) for j, lab in enumerate(LABELS)
            },
        }
        if n_total >= 8:
            row["perm_p_chi2"] = round(perm_p_value_chi2(counts, n_perm, rng), 4)
        delta_buckets.append(row)

    global_counts = np.zeros(4, dtype=np.int64)
    for _, lab, _ in prime_records:
        global_counts[LABEL_TO_IDX[lab]] += 1
    n_all = int(global_counts.sum())
    global_obs = charge_observables(global_counts)
    global_chi2 = chi2_vs_uniform(global_counts)
    global_p = perm_p_value_chi2(global_counts, n_perm, rng)

    asymptotic_rows = asymptotic_D_AB_series(limit, asymptotic_checkpoints)
    o_x = assess_o_X(asymptotic_rows)

    return {
        "section": "4.8.23",
        "mode": "delta-program",
        "script": "eabc_shell_prime_distribution.py",
        "limit": limit,
        "n_max": n_max,
        "m_max": m_max,
        "n_primes": n_all,
        "log_split": {
            "formula": "log N = S(N) + P(N)",
            "S": "alpha*log(2) + beta*log(3)",
            "P": "sum e_i log(p_i)",
            "samples": log_samples,
        },
        "v4_characters": {
            "group": "V4 ≅ C2×C2 ≅ (Z/12Z)×",
            "table": V4_CHARACTERS,
            "prime_examples": prime_char_examples,
        },
        "M_representation_examples": m_examples,
        "delta_buckets": delta_buckets,
        "delta_bin_edges": [round(e, 6) for e in edges],
        "delta_stats": {
            "mean": round(float(np.mean(deltas)), 6),
            "median": round(float(np.median(deltas)), 6),
            "max": round(float(np.max(deltas)), 6),
        },
        "asymptotic_D_AB": {
            "series": asymptotic_rows,
            "o_X_assessment": o_x,
        },
        "test_priorities": {
            "A": "R_{n,m}-Geometrie und Schalen-Gitter",
            "B": "δ(p)-Wachstum, Beschränkung, Häufungen",
            "C": "N_A(δ), N_B(δ), … in δ-Buckets",
            "D": "D_AB, D_AC auf Schalenbändern",
            "E": "Asymptotik D_AB(X)=o(X)?",
        },
        "global": {
            "counts": counts_to_dict(global_counts),
            "fractions": {
                lab: round(int(global_counts[i]) / n_all, 4) for i, lab in enumerate(LABELS)
            },
            "charge": {
                k: round(v, 4) if isinstance(v, float) else v for k, v in global_obs.items()
            },
            "chi2_vs_uniform": round(global_chi2, 4),
            "perm_p_chi2": round(global_p, 4),
            "perm_p_D_AB": round(perm_p_value_D_abs(global_counts, "D_AB", n_perm, rng), 4),
        },
        "n_perm": n_perm,
        "seed": seed,
    }


def format_summary(result: dict[str, Any]) -> str:
    if result.get("mode") == "delta-program":
        g = result["global"]
        fr = g["fractions"]
        ch = g["charge"]
        ds = result["delta_stats"]
        o_x = result["asymptotic_D_AB"]["o_X_assessment"]
        base = (
            f"§4.8.23 delta-program: n={result['n_primes']} Primzahlen≤{result['limit']}, "
            f"δ-buckets={len(result['delta_buckets'])}; "
            f"EABC={fr['E']:.3f}/{fr['A']:.3f}/{fr['B']:.3f}/{fr['C']:.3f}; "
            f"D_AB={ch['D_AB']:.1f}, δ̄={ds['mean']:.2f}; "
            f"o(X)?={o_x.get('o_X_likely')} ({o_x.get('note', '')})"
        )
        return base

    g = result["global"]
    fr = g["fractions"]
    base = (
        f"n={result['n_primes']} Primzahlen>p≤{result['limit']}, "
        f"mode={result['shell_mode']}, buckets={result['n_buckets']}; "
        f"EABC={fr['E']:.3f}/{fr['A']:.3f}/{fr['B']:.3f}/{fr['C']:.3f}, "
        f"χ²={g['chi2_vs_uniform']}, perm_p={g['perm_p_chi2']}"
    )
    if result["shell_mode"] == "defect" and "charge" in g:
        ch = g["charge"]
        base += (
            f"; D_AB={ch['D_AB']:.1f}, D_AC={ch['D_AC']:.1f}, D_E={ch['D_E']:.1f}"
        )
        if "perm_p_D_AB" in g:
            base += (
                f" (perm_p D: {g.get('perm_p_D_AB')}/"
                f"{g.get('perm_p_D_AC')}/{g.get('perm_p_D_E')})"
            )
        if "delta_stats" in result and result["delta_stats"]["mean"] is not None:
            ds = result["delta_stats"]
            base += f"; δ̄={ds['mean']:.3f}"
    return base


def main() -> None:
    parser = argparse.ArgumentParser(
        description="EABC-Primverteilung auf Schalen R_{n,m} (§4.8.21–4.8.23)"
    )
    parser.add_argument(
        "--mode",
        choices=("shell", "delta-program"),
        default="shell",
        help="shell: §4.8.21–22; delta-program: §4.8.23 Defektprogramm",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=50_000,
        help="Obergrenze für Primzahl-Sieb (default: 50000)",
    )
    parser.add_argument(
        "--shell-mode",
        choices=("factor", "log-grid", "defect"),
        default="factor",
        help="Schalen-Zuordnung: factor, log-grid oder defect (§4.8.22)",
    )
    parser.add_argument(
        "--band-mode",
        choices=("nearest", "epsilon", "interval"),
        default="nearest",
        help="defect: nearest R, ε-Band oder Intervall R_k≤p<R_{k+1}",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=None,
        help="defect/epsilon: maximale |p−R_{n,m}| für Band-Zugehörigkeit",
    )
    parser.add_argument(
        "--null-model",
        choices=("perm", "uniform"),
        default="perm",
        help="Nullvergleich: Label-Permutation oder Gleichverteilung",
    )
    parser.add_argument("--n-max", type=int, default=8, help="Gitter max n (default: 8)")
    parser.add_argument("--m-max", type=int, default=8, help="Gitter max m (default: 8)")
    parser.add_argument(
        "--n-delta-bins", type=int, default=8, help="δ-Quantil-Buckets (delta-program)"
    )
    parser.add_argument(
        "--asymptotic-checkpoints",
        type=int,
        default=12,
        help="Checkpoints für D_AB(X) (delta-program)",
    )
    parser.add_argument("--n-perm", type=int, default=499, help="Permutationen (default: 499)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=str, default="", help="JSON-Ausgabedatei")
    args = parser.parse_args()

    if args.mode == "delta-program":
        result = analyze_delta_program(
            limit=args.limit,
            n_max=args.n_max,
            m_max=args.m_max,
            n_delta_bins=args.n_delta_bins,
            n_perm=args.n_perm,
            seed=args.seed,
            asymptotic_checkpoints=args.asymptotic_checkpoints,
        )
    else:
        result = analyze_distribution(
            limit=args.limit,
            shell_mode=args.shell_mode,
            n_max=args.n_max,
            m_max=args.m_max,
            n_perm=args.n_perm,
            seed=args.seed,
            epsilon=args.epsilon,
            band_mode=args.band_mode,
            null_model=args.null_model,
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
