#!/usr/bin/env python3
"""
Experiment (§20.5): Fibonacci-nahe Rekonfigurationspunkte in EABC-Bernoulli-Signaturen.

Scannt n in Fenstern um Fibonacci-Zahlen F_k und vergleicht Observablen (V_n, σ, χ, ι_chir,
|ΔV_n|) mit einer gleichgroßen Zufallsstichprobe aus [1, max_n].

Epistemik: negativer Befund falsifiziert die Fibonacci-Kopplung in dieser Form — kein Beweis
bei positivem Befund. Vgl. collatz_eabc_bernoulli_uebersetzung.md §20.5.

Ausführung:
    python3 collatz_eabc_fibonacci_reconfig_test.py
    python3 collatz_eabc_fibonacci_reconfig_test.py --max-n 500 --window 2 --output report.json
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from collatz_eabc_bernoulli_sensor import (
    BernoulliRow,
    bernoulli_row,
    _sieve_primes,
)
from eabc_from_lean import is_prime_quadruplet

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_fibonacci_reconfig.json"


def fibonacci_up_to(limit: int) -> list[int]:
    """Fibonacci-Zahlen F_k mit 1 <= F_k <= limit."""
    if limit < 1:
        return []
    a, b = 1, 1
    out: list[int] = []
    while a <= limit:
        out.append(a)
        a, b = b, a + b
    return out


def delta_v_norm(prev: BernoulliRow, curr: BernoulliRow) -> int:
    """L1-Norm der Änderung von V_n zwischen aufeinanderfolgenden n."""
    return sum(abs(a - b) for a, b in zip(prev.v.as_tuple(), curr.v.as_tuple()))


def i_chir_flip(prev: BernoulliRow, curr: BernoulliRow) -> bool:
    return prev.i_chir != curr.i_chir and prev.i_chir != 0 and curr.i_chir != 0


def quadruplet_count(sig: list[int]) -> int:
    return sum(1 for p in sig if is_prime_quadruplet(p))


@dataclass(frozen=True, slots=True)
class SamplePoint:
    n: int
    near_fibonacci: int | None
    sigma: int
    chi: int
    i_chir: int
    v_total: int
    delta_v: int
    i_chir_flip: bool
    quadruplet_count: int


def build_samples(
    rows: dict[int, BernoulliRow],
    indices: list[int],
    fib_set: set[int],
    window: int,
) -> list[SamplePoint]:
    out: list[SamplePoint] = []
    for n in sorted(indices):
        if n not in rows:
            continue
        r = rows[n]
        prev = rows.get(n - 1)
        dv = delta_v_norm(prev, r) if prev is not None else 0
        flip = i_chir_flip(prev, r) if prev is not None else False
        near_fib = None
        for f in fib_set:
            if abs(n - f) <= window:
                near_fib = f
                break
        out.append(
            SamplePoint(
                n=n,
                near_fibonacci=near_fib,
                sigma=r.sigma,
                chi=r.chi,
                i_chir=r.i_chir,
                v_total=r.v.total,
                delta_v=dv,
                i_chir_flip=flip,
                quadruplet_count=quadruplet_count(r.prime_sig),
            )
        )
    return out


def _mean_abs(xs: list[int]) -> float:
    return statistics.mean(abs(x) for x in xs) if xs else 0.0


def compare_groups(fib_pts: list[SamplePoint], ref_pts: list[SamplePoint]) -> dict[str, Any]:
    def stats(pts: list[SamplePoint]) -> dict[str, float]:
        return {
            "count": len(pts),
            "mean_abs_sigma": _mean_abs([p.sigma for p in pts]),
            "mean_abs_chi": _mean_abs([p.chi for p in pts]),
            "mean_delta_v": statistics.mean([p.delta_v for p in pts]) if pts else 0.0,
            "i_chir_flip_rate": sum(1 for p in pts if p.i_chir_flip) / len(pts) if pts else 0.0,
            "quadruplet_rate": statistics.mean([p.quadruplet_count for p in pts]) if pts else 0.0,
        }

    fib_s = stats(fib_pts)
    ref_s = stats(ref_pts)
    return {
        "fibonacci_window": fib_s,
        "reference_sample": ref_s,
        "ratio_delta_v": (fib_s["mean_delta_v"] / ref_s["mean_delta_v"])
        if ref_s["mean_delta_v"]
        else None,
        "ratio_flip_rate": (fib_s["i_chir_flip_rate"] / ref_s["i_chir_flip_rate"])
        if ref_s["i_chir_flip_rate"]
        else None,
    }


def run_scan(
    max_n: int = 200,
    window: int = 1,
    seed: int = 42,
    ref_size: int | None = None,
) -> dict[str, Any]:
    primes = _sieve_primes(2 * max_n + 1)
    rows = {n: bernoulli_row(n, primes) for n in range(1, max_n + 1)}

    fibs = fibonacci_up_to(max_n)
    fib_set = set(fibs)
    fib_indices = sorted({n for f in fibs for n in range(max(1, f - window), min(max_n, f + window) + 1)})

    k = ref_size if ref_size is not None else len(fib_indices)
    rng = random.Random(seed)
    pool = [n for n in range(2, max_n + 1) if n not in fib_indices]
    ref_indices = rng.sample(pool, min(k, len(pool)))

    fib_pts = build_samples(rows, fib_indices, fib_set, window)
    ref_pts = build_samples(rows, ref_indices, fib_set, window)

    comparison = compare_groups(fib_pts, ref_pts)
    return {
        "max_n": max_n,
        "window": window,
        "fibonacci_values": fibs,
        "fibonacci_indices": fib_indices,
        "reference_indices": ref_indices,
        "comparison": comparison,
        "fibonacci_samples": [asdict(p) for p in fib_pts],
        "reference_samples": [asdict(p) for p in ref_pts],
        "label": "Experiment (§20.5) — kein Theorem",
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fibonacci-nahe Rekonfigurationspunkte in V_n (§20.5)"
    )
    parser.add_argument("--max-n", type=int, default=200)
    parser.add_argument("--window", type=int, default=1, help="Fenster ±window um F_k")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--ref-size", type=int, default=None)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.max_n < 3:
        raise SystemExit("--max-n muss ≥ 3 sein")
    if args.window < 0:
        raise SystemExit("--window muss ≥ 0 sein")

    report = run_scan(args.max_n, args.window, args.seed, args.ref_size)
    args.output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    cmp_ = report["comparison"]
    print(
        f"Fibonacci-Fenster: {cmp_['fibonacci_window']['count']} Punkte, "
        f"mean |ΔV|={cmp_['fibonacci_window']['mean_delta_v']:.3f}; "
        f"Referenz: mean |ΔV|={cmp_['reference_sample']['mean_delta_v']:.3f}"
    )
    print(f"Geschrieben: {args.output}")


if __name__ == "__main__":
    main()
