#!/usr/bin/env python3
"""Stufe 1 Generalangriff: naive κ_K-Kodierung entlang odd-to-odd-Collatz-Bahnen.

Referenz: CollatzEabc.Kappa.lean, collatz_generalangriff_2026.md.
Kein Collatz-Beweis — nur Injektivitäts- und Dynamiktreue-Tests.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

from eabc_from_lean import EClass, class_of

LETTERS = "EABC"


def nu2(m: int) -> int:
    v = 0
    while m % 2 == 0 and m > 0:
        m //= 2
        v += 1
    return v


def collatz_u(n: int) -> int:
    """Odd-to-odd Collatz-Schritt U(n)."""
    assert n % 2 == 1 and n > 0
    m = 3 * n + 1
    return m // (2 ** nu2(m))


def iterate_u(n: int, k: int) -> int:
    cur = n
    for _ in range(k):
        cur = collatz_u(cur)
    return cur


def kappa_prefix(n: int, k: int) -> list[str | None]:
    """Erste k EABC-Buchstaben (mod 12); None wenn keine EABC-Klasse."""
    return [class_of(iterate_u(n, i)).value if class_of(iterate_u(n, i)) else None for i in range(k)]


def kappa_prefix_word(n: int, k: int) -> str | None:
    letters = kappa_prefix(n, k)
    if any(x is None for x in letters):
        return None
    return "".join(letters)  # type: ignore[arg-type]


@dataclass(frozen=True, slots=True)
class KappaSummary:
    limit: int
    k: int
    odd_starts: int
    defined_count: int
    undefined_step_rate: float
    collision_pairs: int
    collision_rate: float
    dynamics_failures: int
    dynamics_checked: int
    dynamics_holds: bool
    verdict: str


def check_dynamics(n: int, k: int) -> bool:
    """κ_K(U(n))_i = κ_K(n)_{i+1} für i < K-1; letzter Buchstabe = Klasse von iterate_u(n, K)."""
    if n % 2 == 0:
        return True
    u_n = collatz_u(n)
    left = kappa_prefix(u_n, k)
    right_prefix = kappa_prefix(n, k)
    for i in range(k):
        if i + 1 < k:
            expected = right_prefix[i + 1]
        else:
            cls = class_of(iterate_u(n, k))
            expected = cls.value if cls else None
        if left[i] != expected:
            return False
    return True


def run_test(limit: int, k: int) -> tuple[dict, KappaSummary]:
    odd_starts = list(range(1, limit + 1, 2))
    words: dict[int, str | None] = {}
    undefined_steps = 0
    total_steps = 0
    bucket: dict[str | None, list[int]] = defaultdict(list)

    for n in odd_starts:
        w = kappa_prefix_word(n, k)
        words[n] = w
        pref = kappa_prefix(n, k)
        for letter in pref:
            total_steps += 1
            if letter is None:
                undefined_steps += 1
        bucket[w].append(n)

    collision_pairs = sum(len(v) * (len(v) - 1) // 2 for v in bucket.values() if v)
    defined_count = sum(1 for w in words.values() if w is not None)

    dynamics_checked = 0
    dynamics_failures = 0
    for n in odd_starts:
        dynamics_checked += 1
        if not check_dynamics(n, k):
            dynamics_failures += 1

    summary = KappaSummary(
        limit=limit,
        k=k,
        odd_starts=len(odd_starts),
        defined_count=defined_count,
        undefined_step_rate=undefined_steps / total_steps if total_steps else 0.0,
        collision_pairs=collision_pairs,
        collision_rate=collision_pairs / (len(odd_starts) * (len(odd_starts) - 1) // 2)
        if len(odd_starts) > 1
        else 0.0,
        dynamics_failures=dynamics_failures,
        dynamics_checked=dynamics_checked,
        dynamics_holds=dynamics_failures == 0,
        verdict=_verdict(k, defined_count, len(odd_starts), collision_pairs, dynamics_failures),
    )

    detail = {
        "summary": asdict(summary),
        "top_collisions": _top_collisions(bucket, k),
        "sample_undefined": _sample_undefined(odd_starts, k),
    }
    return detail, summary


def _verdict(k: int, defined: int, total: int, collisions: int, dyn_fail: int) -> str:
    if dyn_fail > 0:
        return "DYNAMIK VERLETZT (unerwartet — Lean sagt Shift+Append)"
    if collisions == 0 and defined == total:
        return f"INJEKTIV auf odd n≤N für K={k} (nur numerisch)"
    parts = []
    if defined < total:
        parts.append(f"teilweise undefiniert ({defined}/{total} volle Wörter)")
    if collisions > 0:
        parts.append(f"{collisions} Kollisionspaare bei festem K={k}")
    parts.append("naive κ nicht injektiv / nicht treu")
    return "; ".join(parts)


def _top_collisions(bucket: dict[str | None, list[int]], k: int, top: int = 5) -> list[dict]:
    items = [(w, ns) for w, ns in bucket.items() if w is not None and len(ns) > 1]
    items.sort(key=lambda x: len(x[1]), reverse=True)
    out = []
    for w, ns in items[:top]:
        out.append({"word": w, "count": len(ns), "starts": ns[:8], "k": k})
    return out


def _sample_undefined(odd_starts: list[int], k: int, limit: int = 8) -> list[dict]:
    samples = []
    for n in odd_starts:
        pref = kappa_prefix(n, k)
        if any(x is None for x in pref):
            samples.append({"n": n, "mod12": n % 12, "prefix": pref})
        if len(samples) >= limit:
            break
    return samples


def main() -> None:
    parser = argparse.ArgumentParser(description="κ_K Injektivitäts- und Dynamiktest")
    parser.add_argument("--limit", type=int, default=100_000)
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("collatz_kappa_test.json"),
    )
    args = parser.parse_args()
    detail, summary = run_test(args.limit, args.k)
    args.output.write_text(json.dumps(detail, indent=2), encoding="utf-8")
    print(json.dumps(asdict(summary), indent=2))
    print(f"→ {args.output}")


if __name__ == "__main__":
    main()
