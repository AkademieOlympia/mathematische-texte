from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isqrt
from statistics import mean
from typing import Dict, List, Set


# =========================================================
# Primzahlen, Zwillinge, Vierlinge
# =========================================================

def prime_sieve(limit: int) -> List[bool]:
    sieve = [True] * (limit + 1)
    if limit >= 0:
        sieve[0] = False
    if limit >= 1:
        sieve[1] = False

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = [False] * (((limit - start) // p) + 1)
    return sieve


def is_prime(n: int, prime_flags: List[bool]) -> bool:
    return 0 <= n < len(prime_flags) and prime_flags[n]


def is_twin_center(m: int, prime_flags: List[bool]) -> bool:
    return m - 1 >= 2 and m + 1 < len(prime_flags) and prime_flags[m - 1] and prime_flags[m + 1]


def standard_fourling_centers(limit: int, prime_flags: List[bool]) -> List[int]:
    centers = []
    for p in range(2, limit - 8 + 1):
        if prime_flags[p] and prime_flags[p + 2] and prime_flags[p + 6] and prime_flags[p + 8]:
            centers.append(p + 4)
    return centers


# =========================================================
# Resonanzräume
# =========================================================

def resonance_landscape(four_centers: List[int]) -> List[int]:
    vals = set()
    for c in four_centers:
        vals.add(c - 3)
        vals.add(c)
        vals.add(c + 3)
    return sorted(vals)


def resonance_flanks(four_centers: List[int]) -> List[int]:
    vals = set()
    for c in four_centers:
        vals.add(c - 3)
        vals.add(c + 3)
    return sorted(vals)


def nearest_value(x: int, values: List[int]) -> int:
    best = values[0]
    best_dist = abs(x - best)
    for v in values[1:]:
        d = abs(x - v)
        if d < best_dist or (d == best_dist and v < best):
            best = v
            best_dist = d
    return best


# =========================================================
# Relaxation
# =========================================================

@dataclass
class RelaxationResult:
    start: int
    status: str  # prime, twin_center, cycle, zero, max_steps
    terminal_value: int
    steps: int
    path: List[int]


def relax_difference(
    x0: int,
    anchorspace: List[int],
    prime_flags: List[bool],
    max_steps: int = 100,
) -> RelaxationResult:
    seen: Set[int] = set()
    path = [x0]
    x = x0

    for step in range(max_steps + 1):
        if x in seen:
            return RelaxationResult(x0, "cycle", x, step, path)
        seen.add(x)

        if is_prime(x, prime_flags):
            return RelaxationResult(x0, "prime", x, step, path)

        if is_twin_center(x, prime_flags):
            return RelaxationResult(x0, "twin_center", x, step, path)

        if x == 0:
            return RelaxationResult(x0, "zero", x, step, path)

        a = nearest_value(x, anchorspace)
        x = abs(x - a)
        path.append(x)

    return RelaxationResult(x0, "max_steps", x, max_steps, path)


def relax_modular(
    x0: int,
    anchorspace: List[int],
    prime_flags: List[bool],
    max_steps: int = 100,
) -> RelaxationResult:
    """
    Relaxation über Modulo:
    x_{n+1} = x_n mod nearest_anchor(x_n)
    """
    seen: Set[int] = set()
    path = [x0]
    x = x0

    for step in range(max_steps + 1):
        if x in seen:
            return RelaxationResult(x0, "cycle", x, step, path)
        seen.add(x)

        if is_prime(x, prime_flags):
            return RelaxationResult(x0, "prime", x, step, path)

        if is_twin_center(x, prime_flags):
            return RelaxationResult(x0, "twin_center", x, step, path)

        if x == 0:
            return RelaxationResult(x0, "zero", x, step, path)

        a = nearest_value(x, anchorspace)
        if a == 0:
            return RelaxationResult(x0, "zero", 0, step, path)

        x = x % a
        path.append(x)

    return RelaxationResult(x0, "max_steps", x, max_steps, path)


# =========================================================
# Batch-Auswertung
# =========================================================

@dataclass
class BatchStats:
    name: str
    n: int
    counts: Dict[str, int]
    mean_steps: float
    mean_steps_prime: float
    mean_steps_twin: float
    top_terminal_primes: List[tuple]
    top_terminal_twins: List[tuple]
    top_terminal_all: List[tuple]


def batch_run(
    xs: List[int],
    mode: str,
    anchorspace: List[int],
    prime_flags: List[bool],
    max_steps: int = 100,
    name: str = "",
) -> BatchStats:
    counts = Counter()
    terminal_primes = Counter()
    terminal_twins = Counter()
    terminal_all = Counter()

    all_steps = []
    prime_steps = []
    twin_steps = []

    for x in xs:
        if mode == "difference":
            res = relax_difference(x, anchorspace, prime_flags, max_steps=max_steps)
        elif mode == "modular":
            res = relax_modular(x, anchorspace, prime_flags, max_steps=max_steps)
        else:
            raise ValueError(f"Unbekannter Modus: {mode}")

        counts[res.status] += 1
        all_steps.append(res.steps)
        terminal_all[res.terminal_value] += 1

    return counts


# =========================================================
# Hauptprogramm
# =========================================================

def main() -> None:
    limit = 100_000
    prime_flags = prime_sieve(limit + 100)

    four_centers = standard_fourling_centers(limit, prime_flags)
    landscape = resonance_landscape(four_centers)
    flanks = resonance_flanks(four_centers)

    print(f"Schranke: {limit}")
    print(f"Anzahl Standardvierlingszentren: {len(four_centers)}")
    print(f"Erste Vierlingszentren: {four_centers[:25]}")
    print(f"Erste Flankenpunkte: {flanks[:30]}")

    sample_xs = [14, 30, 42, 60, 72, 105, 195, 822, 1000, 5000, 12345, 54321, 99991]

    print_sample_paths(
        sample_xs,
        "difference",
        flanks,
        prime_flags,
        title="Differenzregel / Flanken",
    )

    print_sample_paths(
        sample_xs,
        "difference",
        landscape,
        prime_flags,
        title="Differenzregel / Landschaft",
    )

    print_sample_paths(
        sample_xs,
        "modular",
        landscape,
        prime_flags,
        title="Modulo-Regel / Landschaft",
    )

    xs = list(range(2, limit + 1))

    stats_diff_flanks = batch_run(
        xs,
        "difference",
        flanks,
        prime_flags,
        name="Differenzregel / Flanken",
    )

    stats_diff_land = batch_run(
        xs,
        "difference",
        landscape,
        prime_flags,
        name="Differenzregel / Landschaft",
    )

    stats_mod_land = batch_run(
        xs,
        "modular",
        landscape,
        prime_flags,
        name="Modulo-Regel / Landschaft",
    )

    print_batch_stats(stats_diff_flanks)
    print_batch_stats(stats_diff_land)
    print_batch_stats(stats_mod_land)


if __name__ == "__main__":
    main()
