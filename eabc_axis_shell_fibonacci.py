#!/usr/bin/env python3
"""
Achsenschalen und Fibonacci-Kompressibilität κ_F (§4.8.25)

Mehrgliedrige Schalen R² = Σ_{j=1}^k q_j^{α_j} mit Primachsen q_j ∈ {2,3,5,7,…}.
Fibonacci-Kompressibilität: s_F(R) = minimale Anzahl Fibonacci-Summanden, die R²
(exakt oder relativ) approximieren; κ_F(R) = k / s_F(R).

Observable:
  C_F(R) = κ_F(R) · (|D_AB| + |D_AC| + |D_BC|)

EABC-Klassen in ε-Umgebung des Schalenradius R; Korrelation κ_F vs. Asymmetrie
(Pearson/Spearman, Permutations-Null).

Referenz: collatz_eabc_zirkulationshypothese.md §4.8.20–4.8.25
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
    primes_above_three,
    sieve_primes,
)

DEFAULT_OUTPUT = Path("eabc_axis_shell_fibonacci.json")
PRIME_AXES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29)


def fibonacci_table(limit: int) -> list[int]:
    """Fibonacci-Zahlen F_1=1, F_2=1, … bis > limit (1-indiziert, ohne F_0=0)."""
    if limit < 1:
        return []
    fibs = [1, 1]
    while fibs[-1] <= limit:
        fibs.append(fibs[-1] + fibs[-2])
    return fibs


def min_fibonacci_summands(
    target: int,
    fibs: list[int],
    *,
    eps_rel: float = 0.0,
    max_summands: int | None = None,
) -> tuple[int, int | None]:
    """
    Minimale Anzahl s_F und beste Fibonacci-Summe für target = R².

    Liefert (s_F, approx_sum). Bei Nicht-Erreichbarkeit: (max_summands oder target, None).
    """
    if target < 1:
        return 1, 0
    tol = 0 if eps_rel <= 0.0 else max(1, int(math.ceil(eps_rel * target)))
    hi = target + tol
    cap = max_summands if max_summands is not None else min(target, 64)
    inf = 10**9
    dp = [inf] * (hi + 1)
    dp[0] = 0
    for f in fibs:
        if f > hi:
            break
        for s in range(hi, f - 1, -1):
            prev = dp[s - f]
            if prev < inf and prev + 1 <= cap and prev + 1 < dp[s]:
                dp[s] = prev + 1
    best_s = inf
    best_sum: int | None = None
    for t in range(target - tol, target + tol + 1):
        if 0 <= t <= hi and dp[t] < best_s:
            best_s = dp[t]
            best_sum = t
    if best_s < inf:
        return int(best_s), best_sum
    return cap, None


def kappa_F(k: int, s_F: int) -> float:
    """Kompressionsgrad κ_F = k / s_F; s_F ≥ 1."""
    s = max(1, s_F)
    return float(k) / float(s)


def axis_shell_r_squared(terms: list[tuple[int, int]]) -> int:
    return int(sum(int(q) ** int(e) for q, e in terms))


def generate_axis_shells(
    k: int,
    exp_max: int,
    *,
    axes: tuple[int, ...] = PRIME_AXES,
    min_exp: int = 1,
) -> list[dict[str, Any]]:
    """Schalenfamilie mit genau k aktiven Primachsen und Exponenten in [min_exp, exp_max]."""
    if k < 1 or k > len(axes):
        raise ValueError(f"k={k} außerhalb der Achsenliste (max {len(axes)})")
    chosen_axes = axes[:k]
    shells: list[dict[str, Any]] = []
    exp_range = range(min_exp, exp_max + 1)
    for exponents in product(exp_range, repeat=k):
        terms = [(chosen_axes[i], exponents[i]) for i in range(k)]
        r_sq = axis_shell_r_squared(terms)
        r = math.sqrt(r_sq)
        shells.append(
            {
                "k": k,
                "terms": [(int(q), int(e)) for q, e in terms],
                "r_squared": int(r_sq),
                "R": float(r),
                "axes_key": ",".join(f"{q}^{e}" for q, e in terms),
            }
        )
    return shells


def generate_twin_two_shells(exp_max: int, min_exp: int = 1) -> list[dict[str, Any]]:
    """Spezialfamilie R² = 2^{n1} + 2^{n2} (k=2, gleiche Achse)."""
    shells: list[dict[str, Any]] = []
    exp_range = range(min_exp, exp_max + 1)
    for n1, n2 in product(exp_range, repeat=2):
        terms = [(2, n1), (2, n2)]
        r_sq = axis_shell_r_squared(terms)
        shells.append(
            {
                "k": 2,
                "family": "twin-2",
                "terms": terms,
                "r_squared": int(r_sq),
                "R": float(math.sqrt(r_sq)),
                "axes_key": f"2^{n1}+2^{n2}",
            }
        )
    return shells


def charge_with_D_BC(counts: np.ndarray) -> dict[str, float]:
    obs = charge_observables(counts)
    n_b = int(counts[LABEL_TO_IDX["B"]])
    n_c = int(counts[LABEL_TO_IDX["C"]])
    obs["D_BC"] = float(n_b - n_c)
    obs["asymmetry_L1"] = abs(obs["D_AB"]) + abs(obs["D_AC"]) + abs(obs["D_BC"])
    return obs


def primes_in_radius_band(
    r: float,
    primes: list[int],
    epsilon: float,
) -> list[int]:
    lo = r - epsilon
    hi = r + epsilon
    return [p for p in primes if lo <= p <= hi]


def eabc_for_shell(
    shell: dict[str, Any],
    primes: list[int],
    epsilon: float,
) -> dict[str, Any]:
    """EABC-Verteilung in ε-Band um Schalenradius R."""
    r = float(shell["R"])
    band = primes_in_radius_band(r, primes, epsilon)
    counts = np.zeros(4, dtype=np.int64)
    for p in band:
        lab = class_of(p)
        if lab is None:
            continue
        counts[LABEL_TO_IDX[lab]] += 1
    n_total = int(counts.sum())
    charge = charge_with_D_BC(counts) if n_total > 0 else {
        "N_E": 0.0,
        "N_A": 0.0,
        "N_B": 0.0,
        "N_C": 0.0,
        "D_AB": 0.0,
        "D_AC": 0.0,
        "D_BC": 0.0,
        "D_E": 0.0,
        "asymmetry_L1": 0.0,
    }
    return {
        "epsilon": epsilon,
        "n_primes_band": n_total,
        "counts": {lab: int(counts[i]) for i, lab in enumerate(LABELS)},
        "charge": {k: round(v, 4) if isinstance(v, float) else v for k, v in charge.items()},
        "primes_sample": band[:8],
    }


def enrich_shell_fibonacci(
    shell: dict[str, Any],
    fibs: list[int],
    eps_rel: float,
) -> dict[str, Any]:
    k = int(shell["k"])
    r_sq = int(shell["r_squared"])
    s_F, fib_sum = min_fibonacci_summands(r_sq, fibs, eps_rel=eps_rel, max_summands=k)
    kap = kappa_F(k, s_F)
    out = dict(shell)
    out["s_F"] = int(s_F)
    out["fibonacci_sum"] = fib_sum
    out["kappa_F"] = round(kap, 6)
    out["compressible"] = bool(s_F < k)
    return out


def correlation_pair(x: np.ndarray, y: np.ndarray) -> dict[str, float | None]:
    if len(x) < 3:
        return {"pearson": None, "spearman": None, "n": float(len(x))}
    if np.std(x) == 0.0 or np.std(y) == 0.0:
        return {"pearson": None, "spearman": None, "n": float(len(x))}
    pearson = float(np.corrcoef(x, y)[0, 1])
    # Spearman via Rank
    rx = np.argsort(np.argsort(x)).astype(np.float64)
    ry = np.argsort(np.argsort(y)).astype(np.float64)
    spearman = float(np.corrcoef(rx, ry)[0, 1])
    return {"pearson": pearson, "spearman": spearman, "n": float(len(x))}


def perm_null_correlation(
    x: np.ndarray,
    y: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
    method: str = "pearson",
) -> dict[str, Any]:
    obs = correlation_pair(x, y)
    key = obs[method]
    if key is None or len(x) < 3:
        return {"observed": key, "null_mean": None, "null_std": None, "p_value": None, "n_perm": n_perm}
    vals: list[float] = []
    y_shuf = y.copy()
    for _ in range(n_perm):
        rng.shuffle(y_shuf)
        c = correlation_pair(x, y_shuf)[method]
        if c is not None:
            vals.append(abs(c))
    arr = np.array(vals, dtype=np.float64)
    p = float(np.mean(arr >= abs(key)))
    return {
        "observed": float(key),
        "null_mean": float(np.mean(arr)),
        "null_std": float(np.std(arr)),
        "p_value": round((p + 1.0 / (n_perm + 1)), 4),
        "n_perm": n_perm,
    }


def run_axis_shell_analysis(
    *,
    k_values: tuple[int, ...] = (2, 3, 4, 5),
    exp_max: int = 4,
    prime_limit: int = 10_000,
    epsilon: float = 2.0,
    eps_rel: float = 0.0,
    min_primes_band: int = 1,
    include_twin_two: bool = True,
    n_perm: int = 499,
    seed: int = 42,
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    primes = primes_above_three(prime_limit)

    shells_raw: list[dict[str, Any]] = []
    for k in k_values:
        shells_raw.extend(generate_axis_shells(k, exp_max))
    if include_twin_two and 2 in k_values:
        for sh in generate_twin_two_shells(exp_max):
            if not any(s["axes_key"] == sh["axes_key"] for s in shells_raw):
                shells_raw.append(sh)

    max_r_sq = max((int(s["r_squared"]) for s in shells_raw), default=1)
    fibs = fibonacci_table(max_r_sq * 2 + 1)

    rows: list[dict[str, Any]] = []
    for sh in shells_raw:
        enriched = enrich_shell_fibonacci(sh, fibs, eps_rel)
        eabc = eabc_for_shell(enriched, primes, epsilon)
        asym = float(eabc["charge"]["asymmetry_L1"])
        c_f = enriched["kappa_F"] * asym
        row = {
            **{k: enriched[k] for k in ("k", "terms", "r_squared", "R", "axes_key", "s_F", "kappa_F", "compressible")},
            "family": enriched.get("family", f"prime-{enriched['k']}"),
            "fibonacci_sum": enriched.get("fibonacci_sum"),
            "eabc": eabc,
            "C_F": round(c_f, 6),
            "asymmetry_L1": round(asym, 4),
        }
        rows.append(row)

    eligible = [r for r in rows if r["eabc"]["n_primes_band"] >= min_primes_band]
    kappa_arr = np.array([r["kappa_F"] for r in eligible], dtype=np.float64)
    asym_arr = np.array([r["asymmetry_L1"] for r in eligible], dtype=np.float64)
    c_f_arr = np.array([r["C_F"] for r in eligible], dtype=np.float64)

    corr_kappa_asym = correlation_pair(kappa_arr, asym_arr)
    perm_kappa_asym = perm_null_correlation(kappa_arr, asym_arr, n_perm, rng, "pearson")
    perm_spearman = perm_null_correlation(kappa_arr, asym_arr, n_perm, rng, "spearman")

    compressible = [r for r in rows if r["compressible"]]
    incompressible = [r for r in rows if not r["compressible"]]

    return {
        "section": "4.8.25",
        "script": "eabc_axis_shell_fibonacci.py",
        "k_values": list(k_values),
        "exp_max": exp_max,
        "prime_limit": prime_limit,
        "epsilon": epsilon,
        "eps_rel_fibonacci": eps_rel,
        "min_primes_band": min_primes_band,
        "n_shells_total": len(rows),
        "n_shells_eligible": len(eligible),
        "n_compressible": len(compressible),
        "correlation": {
            "kappa_F_vs_asymmetry_L1": corr_kappa_asym,
            "perm_null_pearson": perm_kappa_asym,
            "perm_null_spearman": perm_spearman,
            "hypothesis": "große κ_F korreliert mit starken EABC-Asymmetrien (explorativ)",
        },
        "summary_stats": {
            "kappa_F_mean": round(float(np.mean(kappa_arr)), 4) if len(kappa_arr) else None,
            "asymmetry_L1_mean": round(float(np.mean(asym_arr)), 4) if len(asym_arr) else None,
            "C_F_mean": round(float(np.mean(c_f_arr)), 4) if len(c_f_arr) else None,
        },
        "shells": rows,
        "n_perm": n_perm,
        "seed": seed,
    }


def format_summary(report: dict[str, Any]) -> str:
    corr = report["correlation"]["kappa_F_vs_asymmetry_L1"]
    perm = report["correlation"]["perm_null_pearson"]
    lines = [
        f"§{report['section']} Achsenschalen κ_F — {report['n_shells_total']} Schalen, "
        f"{report['n_shells_eligible']} mit ≥{report['min_primes_band']} Prim(en) im ε-Band",
        f"  κ_F vs. Asymmetrie: Pearson={corr.get('pearson')}, Spearman={corr.get('spearman')} (n={corr.get('n')})",
        f"  Perm-Null Pearson: p={perm.get('p_value')}, beobachtet={perm.get('observed')}",
        f"  Kompressibel (s_F<k): {report['n_compressible']}/{report['n_shells_total']}",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Achsenschalen κ_F und Observable C_F (§4.8.25)"
    )
    parser.add_argument("--smoke", action="store_true", help="kleines Gitter, schneller Lauf")
    parser.add_argument("--k-max", type=int, default=5, help="maximale Achsenzahl k (2..k-max)")
    parser.add_argument("--exp-max", type=int, default=4, help="maximaler Exponent je Achse")
    parser.add_argument("--prime-limit", type=int, default=10_000)
    parser.add_argument("--epsilon", type=float, default=2.0, help="ε-Band um R für EABC")
    parser.add_argument("--eps-rel", type=float, default=0.0, help="relative Fibonacci-Toleranz")
    parser.add_argument("--min-primes-band", type=int, default=1)
    parser.add_argument("--n-perm", type=int, default=499)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-twin-two", action="store_true")
    args = parser.parse_args()

    if args.smoke:
        k_values = (2, 3)
        exp_max = 2
        prime_limit = 2_000
        n_perm = 49
        min_primes_band = 1
    else:
        k_values = tuple(range(2, args.k_max + 1))
        exp_max = args.exp_max
        prime_limit = args.prime_limit
        n_perm = args.n_perm
        min_primes_band = args.min_primes_band

    report = run_axis_shell_analysis(
        k_values=k_values,
        exp_max=exp_max,
        prime_limit=prime_limit,
        epsilon=args.epsilon,
        eps_rel=args.eps_rel,
        min_primes_band=min_primes_band,
        include_twin_two=not args.no_twin_two,
        n_perm=n_perm,
        seed=args.seed,
    )
    print(format_summary(report))
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"JSON → {args.output}")


if __name__ == "__main__":
    main()
