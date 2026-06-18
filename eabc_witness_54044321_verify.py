#!/usr/bin/env python3
"""Verifikation EABC-Zeuge p=54_044_321 (ABCE, mod-60060-Kanal 50381).

Tao-Stil: **Beispiel/Zeuge** — kein Holonomie-Theorem.
Siehe collatz_eabc_zirkulationshypothese.md §4.4.
"""

from __future__ import annotations

import argparse

from eabc_hl_coefficient_hypotheses import admissible_mod

P_WITNESS = 54_044_321
MOD_REFINED = 60_060  # 2·3·5·7·11·13
HL_MODULI = (2, 3, 5, 7, 11, 13)
EXPECTED_RESIDUE = 50_381
EXPECTED_MOD12 = 5  # A-Startklasse → ABCE (γ⁺-Vierling)

EXPECTED_LATEST10 = [
    50_381,
    34_211,
    54_911,
    46_301,
    18_281,
    36_011,
    39_041,
    18_371,
    11_651,
    2_411,
]


def is_prime(n: int) -> bool:
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


def orientation(p: int) -> str:
    r = p % 12
    if r == 5:
        return "ABCE"
    if r == 11:
        return "CEAB"
    return f"other({r})"


def hl_admissible_mod60060() -> list[int]:
    """HL-zulässige Startrestklassen mod 60060 (p, p+2, p+6, p+8 alle ≠ 0 mod q)."""
    return [
        r
        for r in range(MOD_REFINED)
        if all(admissible_mod(q, r) for q in HL_MODULI)
    ]


def _prime_sieve(limit: int) -> bytearray:
    sieve = bytearray(b"\x01") * (limit + 9)
    sieve[0:2] = b"\x00\x00"
    for q in range(2, int(limit**0.5) + 1):
        if sieve[q]:
            sieve[q * q : limit + 9 : q] = b"\x00" * (((limit + 8 - q * q) // q) + 1)
    return sieve


def scan_latest_channels(limit: int, top_k: int = 10) -> list[tuple[int, int, str]]:
    """Erstes Auftreten je HL-Kanal mod 60060 (ABCE + CEAB) bis limit."""
    adm = set(hl_admissible_mod60060())
    sieve = _prime_sieve(limit + 8)
    first: dict[int, int] = {}
    orient: dict[int, str] = {}
    for p in range(5, limit + 1, 2):
        if not (
            sieve[p]
            and sieve[p + 2]
            and sieve[p + 6]
            and sieve[p + 8]
        ):
            continue
        r = p % MOD_REFINED
        if r not in adm or r in first:
            continue
        first[r] = p
        orient[r] = orientation(p)
    ranked = sorted(first.items(), key=lambda kv: kv[1], reverse=True)
    return [(r, fp, orient[r]) for r, fp in ranked[:top_k]]


def verify_witness() -> bool:
    p = P_WITNESS
    quad = (p, p + 2, p + 6, p + 8)
    ok = True

    print("=== Zeuge p=54_044_321 ===")
    for n in quad:
        prime = is_prime(n)
        print(f"  {n}: {'Primzahl' if prime else 'NICHT prim'}")
        ok &= prime

    hl_ok = all(admissible_mod(q, p) for q in HL_MODULI)
    print(f"  HL-zulässig mod {MOD_REFINED}: {hl_ok}")
    ok &= hl_ok

    print(f"  p mod 60060 = {p % MOD_REFINED} (erwartet {EXPECTED_RESIDUE})")
    ok &= p % MOD_REFINED == EXPECTED_RESIDUE

    print(f"  p mod 12 = {p % 12} (erwartet {EXPECTED_MOD12} → ABCE)")
    ok &= p % 12 == EXPECTED_MOD12

    print(f"  Orientierung: {orientation(p)}")
    ok &= orientation(p) == "ABCE"

    print(f"\n  Boxed: p={p}, p mod 60060={p % MOD_REFINED}, ABCE")
    return ok


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC-Zeuge 54044321 verifizieren")
    parser.add_argument(
        "--scan",
        action="store_true",
        help="10 späteste HL-Kanäle mod 60060 bis 10^8 (ABCE+CEAB)",
    )
    parser.add_argument("--limit", type=int, default=10**8)
    args = parser.parse_args()

    ok = verify_witness()

    if args.scan:
        print(f"\n=== Späteste HL-Kanäle mod {MOD_REFINED} bis {args.limit:,} ===")
        print(f"  HL-zulässige Kanäle: {len(hl_admissible_mod60060())}")

        latest = scan_latest_channels(args.limit, top_k=10)
        residues = [r for r, _, _ in latest]
        print("  Rang | Kanal | Orient. | erstes p")
        for i, (r, fp, o) in enumerate(latest, 1):
            mark = " ← Zeuge" if fp == P_WITNESS else ""
            print(f"  {i:2d}   | {r:5d} | {o:5s}   | {fp}{mark}")

        print(f"\n  Erwartet (10 späteste): {EXPECTED_LATEST10}")
        print(f"  Reproduziert:          {residues}")
        ok &= residues == EXPECTED_LATEST10

    print(f"\nVerifikation: {'OK' if ok else 'FEHLER'}")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
