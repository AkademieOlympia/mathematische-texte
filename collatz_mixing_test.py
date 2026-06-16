#!/usr/bin/env python3
"""
Numerische Tests für den Uniformitäts-Angriff (Angriff A + B).

Misst:
  - Wartezeit bis E-Typ-Block-Schritt entlang odd-to-odd-Trajektorien
  - Maximale Wartezeit vs. log2(n)
  - 2-adische Distanzen für LTE-Familien n = 2^{k+1}·3^r - 1 und n = 4·3^r - 1
  - Empirische mod-12-Übergangsmatrix und Spektrum
"""

from __future__ import annotations

import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Grundfunktionen
# ---------------------------------------------------------------------------

def nu2(m: int) -> int:
    if m == 0:
        return 0
    v = 0
    while m % 2 == 0:
        m //= 2
        v += 1
    return v


def U(n: int) -> int:
    """Odd-to-odd Collatz-Schritt."""
    assert n % 2 == 1 and n > 0
    return (3 * n + 1) // (2 ** nu2(3 * n + 1))


def eabc_class(n: int) -> str:
    r = n % 12
    return {1: "E", 5: "A", 7: "B", 11: "C"}.get(r, "?")


def block_step(n: int, max_halvings: int = 200) -> tuple[int, str, int]:
    """
    Ein Block-Schritt: von ungeradem n zu nächstem ungeradem m.
    Gibt (m, Startklasse, Anzahl odd-to-odd-Übergänge) zurück.
    """
    assert n % 2 == 1
    start_cls = eabc_class(n)
    steps = 0
    cur = n
    for _ in range(max_halvings):
        if cur % 2 == 0:
            cur //= 2
            continue
        if steps > 0:
            return cur, start_cls, steps
        nxt = U(cur)
        steps += 1
        cur = nxt
    raise RuntimeError(f"block_step: keine Konvergenz für n={n}")


def trajectory_blocks(n0: int, max_blocks: int = 500) -> list[tuple[int, str]]:
    """Liste von (Wert, EABC-Klasse) nach jedem Block."""
    out: list[tuple[int, str]] = []
    n = n0
    if n % 2 == 0:
        while n % 2 == 0 and n > 0:
            n //= 2
    if n <= 0:
        return out
    for _ in range(max_blocks):
        if n == 1:
            out.append((1, "E"))
            break
        cls = eabc_class(n)
        out.append((n, cls))
        n = U(n)
        if n <= 0:
            break
    return out


def wait_until_e_type(n0: int, max_blocks: int = 10_000) -> int | None:
    """Anzahl Block-Schritte bis erster E-Typ (exkl. Start, falls Start ≠ E)."""
    n = n0
    if n % 2 == 0:
        while n % 2 == 0 and n > 0:
            n //= 2
    if n <= 0:
        return None
    if eabc_class(n) == "E":
        return 0
    for k in range(1, max_blocks + 1):
        n = U(n)
        if n == 1:
            return k
        if eabc_class(n) == "E":
            return k
    return None


def dist2_adic(a: int, b: int) -> float:
    """2-adische Distanz |a-b|_2 = 2^{-nu2(a-b)} für a ≠ b."""
    if a == b:
        return 0.0
    return 2.0 ** (-nu2(abs(a - b)))


def lte_worst(k: int, r: int) -> int:
    return (2 ** (k + 1)) * (3 ** r) - 1


def lte_minimal(r: int) -> int:
    return 4 * (3 ** r) - 1


def successor_after_cA_block(n0: int, k: int) -> int:
    """Nachfolger nach C^k A-Block bei Start n0 mit nu2(n0+1) >= k+1."""
    n = n0
    for _ in range(k):
        n = U(n)
    return U(n)


def valuation_n_plus_one(n: int) -> int:
    return nu2(n + 1)


# ---------------------------------------------------------------------------
# Test A: Mischzeit / E-Typ-Wartezeit
# ---------------------------------------------------------------------------

def test_mixing_time(n_max: int = 1_000_000) -> dict:
    max_wait = 0
    max_wait_n = 0
    waits: list[int] = []
    by_start_class: dict[str, list[int]] = defaultdict(list)
    failures = 0

    for n in range(1, n_max + 1):
        w = wait_until_e_type(n)
        if w is None:
            failures += 1
            continue
        waits.append(w)
        cls = eabc_class(n if n % 2 == 1 else n // (2 ** nu2(n)))
        if n % 2 == 0:
            m = n
            while m % 2 == 0:
                m //= 2
            cls = eabc_class(m)
        by_start_class[cls].append(w)
        if w > max_wait:
            max_wait = w
            max_wait_n = n

    waits_sorted = sorted(waits)
    n_w = len(waits)
    pct = lambda p: waits_sorted[int(p * n_w)] if n_w else 0

    # Fit max_wait vs log2(n_max): check if max_wait ~ c * log2(n)
    log_bound = math.log2(n_max) if n_max > 1 else 1.0
    ratio = max_wait / log_bound if log_bound > 0 else float("inf")

    return {
        "n_max": n_max,
        "samples": n_w,
        "failures": failures,
        "max_wait": max_wait,
        "max_wait_n": max_wait_n,
        "max_wait_over_log2_n": ratio,
        "mean_wait": sum(waits) / n_w if n_w else 0,
        "median_wait": pct(0.5),
        "p90_wait": pct(0.9),
        "p99_wait": pct(0.99),
        "max_by_class": {c: max(v) if v else 0 for c, v in by_start_class.items()},
        "mean_by_class": {c: sum(v) / len(v) if v else 0 for c, v in by_start_class.items()},
    }


# ---------------------------------------------------------------------------
# Test B: 2-adische Distanzen LTE-Familien
# ---------------------------------------------------------------------------

def test_lte_2adic(k_max: int = 20, r_max: int = 15) -> dict:
    """Distanz von LTE-Starts zu E-Fixpunkt-Kandidaten und zu 1 (Zyklus)."""
    results = []
    min_dist_to_one = float("inf")
    min_pair = None

    for k in range(1, k_max + 1):
        for r in range(1, r_max + 1):
            n0 = lte_worst(k, r)
            n_succ = successor_after_cA_block(n0, k)
            v_succ = valuation_n_plus_one(n_succ)
            d_one = dist2_adic(n0, 1)
            d_succ_one = dist2_adic(n_succ, 1)
            # Uniformitäts-Schranke c·2^{-floor(log2 n)}
            logn = math.floor(math.log2(n0)) if n0 > 0 else 0
            bound = 2.0 ** (-logn) if logn > 0 else 1.0
            ratio = d_one / bound if bound > 0 else 0
            results.append({
                "k": k, "r": r, "n0": n0,
                "nu2_succ_plus1": v_succ,
                "dist2_to_1": d_one,
                "dist2_succ_to_1": d_succ_one,
                "uniform_bound": bound,
                "ratio_dist_over_bound": ratio,
                "eabc_succ": eabc_class(n_succ),
            })
            if d_one < min_dist_to_one:
                min_dist_to_one = d_one
                min_pair = (k, r, n0)

    # Spezialfall 4·3^r - 1
    special = []
    for r in range(1, r_max + 1):
        n0 = lte_minimal(r)
        n1 = U(n0)
        n2 = U(n1)
        special.append({
            "r": r, "n0": n0,
            "n2": n2,
            "nu2_n2_plus1": valuation_n_plus_one(n2),
            "eabc_n2": eabc_class(n2),
            "dist2_to_1": dist2_adic(n0, 1),
        })

    return {
        "lte_worst_samples": len(results),
        "min_dist2_to_1": min_dist_to_one,
        "min_dist_pair": min_pair,
        "worst_ratio_bound": max(r["ratio_dist_over_bound"] for r in results),
        "all_nu2_succ_le_2": all(r["nu2_succ_plus1"] <= 2 for r in results),
        "lte_worst_table": results[:12],
        "lte_minimal_table": special[:10],
    }


# ---------------------------------------------------------------------------
# Empirische mod-12-Matrix
# ---------------------------------------------------------------------------

def empirical_transition_matrix(n_max: int = 500_000) -> dict:
  counts: dict[str, Counter[str]] = {c: Counter() for c in "EABC"}
  for n in range(1, n_max + 1, 2):
      if n > n_max:
          break
      c0 = eabc_class(n)
      if c0 == "?":
          continue
      c1 = eabc_class(U(n))
      if c1 == "?":
          continue
      counts[c0][c1] += 1

  matrix = {}
  for s in "EABC":
      total = sum(counts[s].values())
      matrix[s] = {t: counts[s][t] / total if total else 0 for t in "EABC"}

  # Spektrum der 4x4-Matrix
  import numpy as np
  order = ["E", "A", "B", "C"]
  M = np.array([[matrix[s][t] for t in order] for s in order])
  eigvals = np.linalg.eigvals(M)
  moduli = sorted([abs(x) for x in eigvals], reverse=True)

  return {
      "n_odd_samples": sum(sum(counts[s].values()) for s in "EABC"),
      "matrix": matrix,
      "eigenvalue_moduli": [float(x) for x in moduli],
      "lambda2_abs": float(moduli[1]) if len(moduli) > 1 else 0.0,
      "min_entry": float(M.min()),
  }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    n_max = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
    print(f"=== Collatz Uniformitäts-Tests (n ≤ {n_max:,}) ===\n")

    print("Angriff A: E-Typ-Mischzeit ...")
    mix = test_mixing_time(n_max)
    print(json.dumps(mix, indent=2))

    print("\nAngriff B: LTE 2-adische Distanzen ...")
    adic = test_lte_2adic()
    print(json.dumps(adic, indent=2, default=str))

    print("\nEmpirische mod-12-Matrix (n ungerade ≤ 500k) ...")
    em = empirical_transition_matrix(min(n_max, 500_000))
    print(json.dumps(em, indent=2))

    out = Path(__file__).with_suffix(".json")
    out.write_text(json.dumps({"mixing": mix, "adic": adic, "empirical_matrix": em}, indent=2))
    print(f"\nErgebnisse gespeichert: {out}")


if __name__ == "__main__":
    main()
