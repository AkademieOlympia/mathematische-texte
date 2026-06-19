#!/usr/bin/env python3
"""EABC occupancy tree — Scan → Blockzustände → Reduktionsbaum.

Kernidee (Cook–Mertz / Williams): lange Historie wird nicht materialisiert,
sondern blockweise auf monoidische Zustände
    Z = (O, T, n)
reduziert und assoziativ gemerged:
    (Z_1 ⊕ Z_2) ⊕ Z_3 = Z_1 ⊕ (Z_2 ⊕ Z_3).

M = 60060, Kanäle C_M (HL-zulässige Einheiten mod M).
Für N: O(N) besetzte Kanäle, T(c) Erstbesetzung, n(c) Ereigniszähler,
rho(N) = |O(N)| / |C_M|.  EABC-Faser π₁₂(c) ∈ {E, A, B, C} via c mod 12.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from collections import Counter, defaultdict
from math import gcd, log2
from pathlib import Path
from typing import Any

import numpy as np

M = 60_060
HL_MODULI = (2, 3, 5, 7, 11, 13)
OFFSETS = (0, 2, 6, 8)
EABC = {1: "E", 5: "A", 7: "B", 11: "C"}
EABC_FIBERS = ("E", "A", "B", "C")

OccupancyState = tuple[set[int], dict[int, int], dict[int, int]]


# ── Schicht 0: Kanaluniversum ────────────────────────────────────────────────


def admissible_quad_channels(M: int = M) -> list[int]:
    """C_M: c ∈ (ℤ/Mℤ)× mit c, c+2, c+6, c+8 jeweils coprim zu M."""
    return [
        c
        for c in range(M)
        if gcd(c, M) == 1 and all(gcd(c + d, M) == 1 for d in OFFSETS)
    ]


def channel_universe() -> frozenset[int]:
    """Frozenset-Variante von admissible_quad_channels (für schnelle Membership)."""
    return frozenset(admissible_quad_channels())


def eabc_class(c: int) -> str:
    return EABC[c % 12]


eabc_fiber = eabc_class  # Alias (bestehende Tests/Skripte)


# ── Schicht 1: Monoid-Merge ───────────────────────────────────────────────────


def identity_state() -> OccupancyState:
    """Neutrales Element Z_0 = (∅, ∅, 0) des Besetzungs-Monoids."""
    return set(), {}, Counter()


def state_equal(Z1: OccupancyState, Z2: OccupancyState) -> bool:
    """Gleichheit auf (O, T, n); n vergleicht als Multimenge (Counter)."""
    O1, T1, n1 = Z1
    O2, T2, n2 = Z2
    return O1 == O2 and T1 == T2 and Counter(n1) == Counter(n2)


def merge_state(Z1: OccupancyState, Z2: OccupancyState) -> OccupancyState:
    """Monoid-Merge: Z_1 ⊕ Z_2 = (O_1 ∪ O_2, min(T_1, T_2), n_1 + n_2)."""
    O1, T1, n1 = Z1
    O2, T2, n2 = Z2
    O = O1 | O2
    T = dict(T1)
    for c, t in T2.items():
        T[c] = min(T.get(c, t), t)
    n = dict(n1)
    for c, k in n2.items():
        n[c] = n.get(c, 0) + k
    return O, T, n


def empty_state() -> OccupancyState:
    return identity_state()


# ── Schicht 2: Scan (lokaler Block) ──────────────────────────────────────────


def is_quadruplet(p: int, base_primes: list[int] | None = None) -> bool:
    """Prüft (p, p+2, p+6, p+8) auf Primheit (Referenzimplementierung)."""
    if p < 5:
        return False
    for d in OFFSETS:
        n = p + d
        if n < 2:
            return False
        if n < 4:
            continue
        if n % 2 == 0:
            return False
        r = int(math.isqrt(n))
        if base_primes is not None:
            for q in base_primes:
                if q > r:
                    break
                if n % q == 0:
                    return False
        else:
            for q in range(3, r + 1, 2):
                if n % q == 0:
                    return False
    return True


def small_primes_upto(n: int) -> list[int]:
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            sieve[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    return [p for p in range(2, n + 1) if sieve[p]]


def odd_segment_prime_flags(L: int, R: int, base_primes: list[int]) -> bytearray:
    assert L % 2 == 1 and R % 2 == 1
    m = (R - L) // 2 + 1
    flags = bytearray(b"\x01") * m
    for q in base_primes:
        if q == 2:
            continue
        q2 = q * q
        if q2 > R:
            break
        start = max(q2, ((L + q - 1) // q) * q)
        if start % 2 == 0:
            start += q
        idx = (start - L) // 2
        if idx < m:
            flags[idx::q] = b"\x00" * (((m - 1 - idx) // q) + 1)
    return flags


def scan_block(
    a: int,
    b: int,
    channels: frozenset[int] | list[int] | set[int],
    base_primes: list[int] | None = None,
) -> OccupancyState:
    """Lokaler Blockzustand Z_i für Primzahlvierlinge p ∈ [a, b] auf Kanälen C_M."""
    channel_set = channels if isinstance(channels, frozenset) else frozenset(channels)
    if b < max(5, a):
        return empty_state()

    if base_primes is None:
        root = int(math.isqrt(b + max(OFFSETS))) + 1
        base_primes = small_primes_upto(root)

    L = max(5, a)
    if L % 2 == 0:
        L += 1
    H = b
    if H % 2 == 0:
        H -= 1

    R = H + max(OFFSETS)
    if R % 2 == 0:
        R += 1

    flags = odd_segment_prime_flags(L, R, base_primes)
    arr = np.frombuffer(flags, dtype=np.uint8)
    npos = (H - L) // 2 + 1
    mask = arr[0:npos] & arr[1 : npos + 1] & arr[3 : npos + 3] & arr[4 : npos + 4]
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return empty_state()

    p = L + 2 * idx.astype(np.int64)
    residues = p % M
    channel_mask = np.fromiter((r in channel_set for r in residues), dtype=bool, count=len(residues))
    if not channel_mask.any():
        return empty_state()

    p_hit = p[channel_mask]
    c_hit = residues[channel_mask]

    O: set[int] = set()
    T: dict[int, int] = {}
    n: dict[int, int] = {}
    for pi, ci in zip(p_hit.tolist(), c_hit.tolist(), strict=True):
        O.add(ci)
        if ci not in T:
            T[ci] = pi
        n[ci] = n.get(ci, 0) + 1
    return O, T, n


# ── Schicht 3: Reduktionsbaum ─────────────────────────────────────────────────


def reduce_tree(states: list[OccupancyState]) -> OccupancyState:
    """Pairwise assoziatives Merge (Cook–Mertz / Williams-Baum)."""
    if not states:
        return empty_state()
    layer = list(states)
    while len(layer) > 1:
        nxt: list[OccupancyState] = []
        for i in range(0, len(layer), 2):
            if i + 1 < len(layer):
                nxt.append(merge_state(layer[i], layer[i + 1]))
            else:
                nxt.append(layer[i])
        layer = nxt
    return layer[0]


def occupancy_tree(
    N: int,
    block_size: int = 1_000_000,
    M_mod: int = M,
) -> dict[str, Any]:
    """Scan → Blockzustände → Reduktionsbaum → report_state."""
    channels = admissible_quad_channels(M_mod)
    channel_set = frozenset(channels)
    root = int(math.isqrt(N + max(OFFSETS))) + 1
    base_primes = small_primes_upto(root)

    states: list[OccupancyState] = []
    for a in range(1, N + 1, block_size):
        b = min(N, a + block_size - 1)
        states.append(scan_block(a, b, channel_set, base_primes))

    Z = reduce_tree(states)
    return report_state(Z, channels)


# ── Bericht / Entropie ────────────────────────────────────────────────────────


def entropy_counts(counter: Counter[int] | dict[int, int]) -> float:
    """Shannon-Entropie in Bit (log₂)."""
    total = sum(counter.values())
    if total == 0:
        return 0.0
    return -sum((v / total) * log2(v / total) for v in counter.values())


def report_state(Z: OccupancyState, channels: list[int]) -> dict[str, Any]:
    O, T, n = Z
    rho = len(O) / len(channels) if channels else 0.0
    fiber_total: dict[str, int] = defaultdict(int)
    fiber_seen: dict[str, int] = defaultdict(int)
    for c in channels:
        fiber_total[eabc_class(c)] += 1
        if c in O:
            fiber_seen[eabc_class(c)] += 1
    latest = sorted(T.items(), key=lambda x: x[1], reverse=True)[:10]
    return {
        "channels_total": len(channels),
        "occupied": len(O),
        "rho": rho,
        "entropy_60060": entropy_counts(n),
        "fiber_total": dict(fiber_total),
        "fiber_seen": dict(fiber_seen),
        "latest_first_hits": latest,
    }


def shannon_entropy(counts: dict[str, int]) -> float:
    """Nats-Entropie (Legacy-Hilfsfunktion)."""
    total = sum(counts.values())
    if total == 0:
        return 0.0
    h = 0.0
    for k in counts.values():
        if k > 0:
            p = k / total
            h -= p * math.log(p)
    return h


def fiber_statistics(O: set[int], T: dict[int, int], n: dict[int, int]) -> dict[str, Any]:
    channels = {f: 0 for f in EABC_FIBERS}
    events = {f: 0 for f in EABC_FIBERS}
    first_occ: dict[str, list[int]] = {f: [] for f in EABC_FIBERS}

    for c in O:
        f = eabc_class(c)
        channels[f] += 1
        if c in T:
            first_occ[f].append(T[c])
        events[f] += n.get(c, 0)

    return {
        "channels_by_fiber": channels,
        "events_by_fiber": events,
        "channel_entropy": shannon_entropy(channels),
        "event_entropy": shannon_entropy(events),
        "mean_first_occurrence": {
            f: (sum(first_occ[f]) / len(first_occ[f]) if first_occ[f] else None)
            for f in EABC_FIBERS
        },
    }


def latest_occupation(T: dict[int, int]) -> dict[str, Any] | None:
    if not T:
        return None
    t_max = max(T.values())
    ch = sorted(c for c, t in T.items() if t == t_max)
    return {
        "T_max": t_max,
        "channels": ch,
        "fiber": eabc_class(ch[0]),
    }


def scan_occupancy(N: int, block_width: int = 5_000_000) -> dict[str, Any]:
    """Erweiterter Lauf mit Zeitmessung und Faserdiagnostik (CLI / Tests)."""
    channels = channel_universe()
    root = int(math.isqrt(N + max(OFFSETS))) + 1
    base_primes = small_primes_upto(root)

    block_states: list[OccupancyState] = []
    L = 5
    if L % 2 == 0:
        L += 1

    t0 = time.time()
    while L <= N:
        H = min(N, L + block_width - 1)
        if H % 2 == 0:
            H -= 1
        block_states.append(scan_block(L, H, channels, base_primes))
        L = H + 2

    O, T, n = reduce_tree(block_states)
    elapsed = time.time() - t0

    card_C = len(channels)
    card_O = len(O)
    rho = card_O / card_C if card_C else 0.0

    return {
        "N": N,
        "M": M,
        "card_C_M": card_C,
        "card_O_N": card_O,
        "rho_N": rho,
        "total_events": sum(n.values()),
        "num_blocks": len(block_states),
        "latest_T": latest_occupation(T),
        "fiber_stats": fiber_statistics(O, T, n),
        "report": report_state((O, T, n), sorted(channels)),
        "elapsed_sec": elapsed,
    }


def print_report(report: dict[str, Any]) -> None:
    print("=== EABC_OCCUPANCY_TREE ===")
    print(f"N = {report['N']:,}   M = {report['M']:,}")
    print(f"|C_M| = {report['card_C_M']}")
    print(f"|O(N)| = {report['card_O_N']}")
    print(f"rho(N) = {report['rho_N']:.6f}")
    print(f"total events = {report['total_events']}")
    print(f"blocks = {report['num_blocks']}   elapsed = {report['elapsed_sec']:.2f}s")

    latest = report["latest_T"]
    if latest:
        ch = latest["channels"]
        ch_str = str(ch[0]) if len(ch) == 1 else str(ch)
        print(
            f"späteste T(c) = {latest['T_max']:,}  "
            f"Kanal {ch_str}  Faser {latest['fiber']}"
        )

    fs = report["fiber_stats"]
    print("\nEABC-Faserstatistik (besetzte Kanäle):")
    for f in EABC_FIBERS:
        print(f"  {f}: {fs['channels_by_fiber'][f]}")
    print(f"  H_channel = {fs['channel_entropy']:.6f} nats")


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC_OCCUPANCY_TREE")
    parser.add_argument("--N", type=int, default=10**8)
    parser.add_argument(
        "--block",
        type=int,
        default=1_000_000,
        help="Blockbreite für Baum-Reduktion",
    )
    parser.add_argument("--json", type=str, default=None, help="JSON-Ausgabepfad")
    parser.add_argument(
        "--simple",
        action="store_true",
        help="Nur occupancy_tree/report_state (Williams-API)",
    )
    args = parser.parse_args()

    if args.simple:
        result = occupancy_tree(args.N, block_size=args.block)
        for k, v in result.items():
            print(k, "=", v)
        return

    report = scan_occupancy(args.N, block_width=args.block)
    print_report(report)

    if args.json:
        out = Path(args.json)
        out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
        print(f"\nJSON: {out.resolve()}")


if __name__ == "__main__":
    main()
