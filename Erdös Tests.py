#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
erdos_binom_scan.py

Numerische Arbeitsumgebung für die Untersuchung der Aussage:

    Für 1 <= i < j <= n/2 existiert eine Primzahl p >= i mit
    p | gcd(binomial(n,i), binomial(n,j)).

Enthaltene Pakete:
1. Restfälle für i = 3, n <= 10
2. Bestimmung kleiner B(i)
3. Vollscan für 4 <= i <= 9
4. Tests des schwachen Escape-Lemmas
5. Tests des Kontroll-/Adjazenzarguments

Nur Standardbibliothek.
"""

from math import comb, gcd, isqrt
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple, Any
import json
import hashlib
import time


# =========================================================
# Grundfunktionen
# =========================================================

def vp(n: int, p: int) -> int:
    """p-adische Bewertung v_p(n)."""
    e = 0
    while n % p == 0 and n > 0:
        n //= p
        e += 1
    return e


def prime_factors(n: int) -> List[int]:
    """Gibt die verschiedenen Primfaktoren von n zurück."""
    if n < 2:
        return []
    out = []
    if n % 2 == 0:
        out.append(2)
        while n % 2 == 0:
            n //= 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 2
    if n > 1:
        out.append(n)
    return out


def factor_with_exponents(n: int) -> List[Tuple[int, int]]:
    """Primfaktorzerlegung mit Exponenten."""
    if n < 2:
        return []
    out = []
    e = 0
    while n % 2 == 0:
        n //= 2
        e += 1
    if e:
        out.append((2, e))
    d = 3
    while d * d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e:
            out.append((d, e))
        d += 2
    if n > 1:
        out.append((n, 1))
    return out


def primes_up_to(n: int) -> List[int]:
    """Sieb der Eratosthenes bis n."""
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:n + 1:p] = [False] * (((n - start) // p) + 1)
    return [k for k, ok in enumerate(sieve) if ok]


def is_S_smooth(x: int, S: List[int]) -> bool:
    """Test, ob x nur Primteiler aus S besitzt."""
    if x < 1:
        return False
    if x == 1:
        return True
    y = x
    for p in S:
        while y % p == 0:
            y //= p
    return y == 1


def witnesses(n: int, i: int, j: int) -> List[int]:
    """
    Liefert alle Primteiler p >= i von gcd(binomial(n,i), binomial(n,j)).
    """
    A = comb(n, i)
    B = comb(n, j)
    G = gcd(A, B)
    ps = prime_factors(G)
    return [p for p in ps if p >= i]


def witness_exists(n: int, i: int, j: int) -> Tuple[bool, List[int]]:
    ws = witnesses(n, i, j)
    return len(ws) > 0, ws


# =========================================================
# Paket A: Restfälle für i = 3
# =========================================================

def scan_i3_small(n_min: int = 6, n_max: int = 10) -> Dict[str, Any]:
    """
    Prüft alle Tripel (n,3,j) mit 4 <= j <= n//2 und n <= 10.
    """
    rows = []
    bad = []
    total = 0

    for n in range(n_min, n_max + 1):
        for j in range(4, n // 2 + 1):
            total += 1
            ws = witnesses(n, 3, j)
            row = {
                "n": n,
                "i": 3,
                "j": j,
                "witnesses": ws,
            }
            rows.append(row)
            if not ws:
                bad.append((n, 3, j))

    return {
        "package": "i3_small_scan",
        "n_min": n_min,
        "n_max": n_max,
        "tested": total,
        "bad_count": len(bad),
        "bad_examples": bad,
        "rows": rows,
    }


# =========================================================
# Paket B: B(i) für kleine i
# =========================================================

def compute_Bi(i: int, search_limit: int) -> int:
    """
    Berechnet B(i) = max{x : x und x-1 sind S_i-glatt}
    durch brute force bis search_limit.
    """
    S = primes_up_to(i)
    best = None
    for x in range(2, search_limit + 1):
        if is_S_smooth(x, S) and is_S_smooth(x - 1, S):
            best = x
    return best if best is not None else -1


def recompute_small_B() -> Dict[str, Any]:
    """
    Bestimmt die kleinen B(i) für i=4..9.
    Grenzen so gewählt, dass die bekannten Kandidaten erfasst werden.
    """
    limits = {
        4: 100,
        5: 500,
        6: 500,
        7: 5000,
        8: 5000,
        9: 5000,
    }
    out = {}
    for i in range(4, 10):
        Bi = compute_Bi(i, limits[i])
        out[i] = {
            "i": i,
            "search_limit": limits[i],
            "B_i": Bi,
            "N_i": Bi + i - 1 if Bi >= 0 else None,
            "S_i": primes_up_to(i),
        }
    return {
        "package": "recompute_small_B",
        "results": out,
    }


# =========================================================
# Paket C: Vollscan für 4 <= i <= 9
# =========================================================

DEFAULT_BOUNDS = {
    4: 12,
    5: 85,
    6: 86,
    7: 4381,
    8: 4382,
    9: 4383,
}

# Reduzierte Schranken für schnellen Testlauf (Vollscan dauert sonst Stunden)
QUICK_BOUNDS = {
    4: 12,
    5: 85,
    6: 86,
    7: 100,
    8: 100,
    9: 100,
}


def blocked_candidates_for_i(i: int, Nmax: int) -> List[Tuple[int, int, int]]:
    """
    Liefert alle blockierten Kandidaten (n,i,j) mit n <= Nmax.
    """
    bad = []
    for n in range(2 * i, Nmax + 1):
        for j in range(i + 1, n // 2 + 1):
            ok, _ = witness_exists(n, i, j)
            if not ok:
                bad.append((n, i, j))
    return bad


def full_small_scan(bounds: Dict[int, int] = None, verbose: bool = True) -> Dict[str, Any]:
    """
    Vollscan für i=4..9 mit vorgegebenen Schranken.
    """
    if bounds is None:
        bounds = DEFAULT_BOUNDS.copy()

    summary = {}
    all_bad = []
    grand_total = 0

    for i in sorted(bounds):
        N = bounds[i]
        total = 0
        bad = []
        if verbose:
            print(f"  full_small_scan: i={i}, N={N} ...", flush=True)

        for n in range(2 * i, N + 1):
            if verbose and (n - 2 * i) % 500 == 0 and n > 2 * i:
                print(f"    n={n}/{N}", flush=True)
            for j in range(i + 1, n // 2 + 1):
                total += 1
                ok, ws = witness_exists(n, i, j)
                if not ok:
                    bad.append((n, i, j))
                    all_bad.append((n, i, j))

        grand_total += total
        if verbose:
            print(f"    i={i}: {total} getestet, {len(bad)} bad", flush=True)
        summary[i] = {
            "i": i,
            "N_i": N,
            "tested": total,
            "bad_count": len(bad),
            "bad_examples": bad[:20],
        }

    return {
        "package": "full_small_scan",
        "bounds": bounds,
        "grand_total_tested": grand_total,
        "grand_bad_count": len(all_bad),
        "summary": summary,
        "all_bad_examples": all_bad[:100],
    }


# =========================================================
# Paket D: Schwaches Escape-Lemma testen
# =========================================================

def weak_escape_violations(n_max: int, i_max: int) -> Dict[str, Any]:
    """
    Sucht Verletzungen der schwachen Escape-Version:
    Wenn x = n-k einen Primteiler q > i mit
        q in (j-i, j], v_q(x)=1
    besitzt und x = q*m,
    dann sollen alle Primteiler r>i von m ebenfalls in (j-i, j] liegen.
    """
    violations = []
    tested_lonely_like = 0

    for i in range(3, i_max + 1):
        for n in range(2 * i, n_max + 1):
            for j in range(i + 1, n // 2 + 1):
                left = j - i
                right = j
                for k in range(i):
                    x = n - k
                    pfs = prime_factors(x)
                    for q in pfs:
                        if q > i and left < q <= right and vp(x, q) == 1:
                            tested_lonely_like += 1
                            m = x // q
                            big_rs = [r for r in prime_factors(m) if r > i]
                            bad_rs = [r for r in big_rs if not (left < r <= right)]
                            if bad_rs:
                                violations.append({
                                    "n": n,
                                    "i": i,
                                    "j": j,
                                    "k": k,
                                    "x": x,
                                    "q": q,
                                    "m": m,
                                    "big_rs": big_rs,
                                    "bad_rs": bad_rs,
                                })

    return {
        "package": "weak_escape_violations",
        "n_max": n_max,
        "i_max": i_max,
        "tested_lonely_like_cases": tested_lonely_like,
        "violation_count": len(violations),
        "violations": violations[:200],
    }


# =========================================================
# Paket E: Kontroll-/Adjazenzargument testen
# =========================================================

def controlled_positions(n: int, i: int) -> Tuple[List[Tuple[int, int, List[int]]], List[Tuple[int, int, List[int]]]]:
    """
    Kontrolliert/unkontrolliert rein nach der Blockdefinition:
    kontrolliert <=> n-t besitzt einen Primteiler > i.
    """
    controlled = []
    uncontrolled = []
    for t in range(i):
        x = n - t
        pf = prime_factors(x)
        if any(p > i for p in pf):
            controlled.append((t, x, pf))
        else:
            uncontrolled.append((t, x, pf))
    return controlled, uncontrolled


def adjacent_uncontrolled_pairs(n: int, i: int) -> List[Tuple[int, int, int]]:
    """
    Liefert Paare (t, n-t, n-(t+1)) benachbarter unkontrollierter Positionen.
    """
    controlled, uncontrolled = controlled_positions(n, i)
    unc_t = {t for t, _, _ in uncontrolled}
    pairs = []
    for t in range(i - 1):
        if t in unc_t and (t + 1) in unc_t:
            pairs.append((t, n - t, n - (t + 1)))
    return pairs


def adjacency_profile(n_max: int, i_min: int = 10, i_max: int = 30) -> Dict[str, Any]:
    """
    Liefert Profilinformationen über kontrollierte Positionen und benachbarte
    unkontrollierte Paare für viele (n,i).
    """
    rows = []
    no_adjacent = []

    for i in range(i_min, i_max + 1):
        for n in range(2 * i, n_max + 1):
            controlled, uncontrolled = controlled_positions(n, i)
            pairs = adjacent_uncontrolled_pairs(n, i)
            row = {
                "n": n,
                "i": i,
                "controlled_count": len(controlled),
                "uncontrolled_count": len(uncontrolled),
                "adjacent_pair_count": len(pairs),
                "adjacent_pairs_preview": pairs[:5],
            }
            rows.append(row)
            if not pairs:
                no_adjacent.append({
                    "n": n,
                    "i": i,
                    "controlled_count": len(controlled),
                    "uncontrolled_count": len(uncontrolled),
                })

    return {
        "package": "adjacency_profile",
        "n_max": n_max,
        "i_min": i_min,
        "i_max": i_max,
        "rows_count": len(rows),
        "no_adjacent_count": len(no_adjacent),
        "no_adjacent_examples": no_adjacent[:100],
        "rows_preview": rows[:100],
    }


# =========================================================
# Bericht / Export
# =========================================================

def sha256_of_object(obj: Any) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def write_json(filename: str, obj: Any) -> None:
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)


def compact_summary(report: Dict[str, Any]) -> Dict[str, Any]:
    """
    Erzeugt eine knappe Zusammenfassung für schnellen Blick.
    """
    summary = {
        "generated_at_unix": time.time(),
        "sections": {}
    }

    if "i3_small" in report:
        x = report["i3_small"]
        summary["sections"]["i3_small"] = {
            "tested": x["tested"],
            "bad_count": x["bad_count"],
        }

    if "small_B" in report:
        x = report["small_B"]["results"]
        summary["sections"]["small_B"] = {
            str(i): {
                "B_i": x[i]["B_i"],
                "N_i": x[i]["N_i"],
            } for i in x
        }

    if "small_scan" in report:
        x = report["small_scan"]
        summary["sections"]["small_scan"] = {
            "grand_total_tested": x["grand_total_tested"],
            "grand_bad_count": x["grand_bad_count"],
            "per_i": {
                str(i): {
                    "tested": x["summary"][i]["tested"],
                    "bad_count": x["summary"][i]["bad_count"],
                } for i in x["summary"]
            }
        }

    if "weak_escape" in report:
        x = report["weak_escape"]
        summary["sections"]["weak_escape"] = {
            "tested_lonely_like_cases": x["tested_lonely_like_cases"],
            "violation_count": x["violation_count"],
        }

    if "adjacency" in report:
        x = report["adjacency"]
        summary["sections"]["adjacency"] = {
            "rows_count": x["rows_count"],
            "no_adjacent_count": x["no_adjacent_count"],
        }

    summary["sha256"] = sha256_of_object(summary)
    return summary


# =========================================================
# Hauptlauf
# =========================================================

def run_all(
    run_i3_small: bool = True,
    run_small_B: bool = True,
    run_small_scan: bool = True,
    run_weak_escape: bool = True,
    run_adjacency: bool = True,
    quick: bool = True,
    weak_escape_n_max: int = 300,
    weak_escape_i_max: int = 20,
    adjacency_n_max: int = 250,
    adjacency_i_min: int = 10,
    adjacency_i_max: int = 25,
) -> Dict[str, Any]:
    report = {}
    bounds = QUICK_BOUNDS if quick else DEFAULT_BOUNDS
    if quick:
        print("Schnellmodus (quick=True): reduzierte Schranken für full_small_scan", flush=True)

    if run_i3_small:
        print("Paket A: i3_small ...", flush=True)
        report["i3_small"] = scan_i3_small()

    if run_small_B:
        print("Paket B: small_B ...", flush=True)
        report["small_B"] = recompute_small_B()

    if run_small_scan:
        print("Paket C: full_small_scan ...", flush=True)
        report["small_scan"] = full_small_scan(bounds=bounds)

    if run_weak_escape:
        print("Paket D: weak_escape ...", flush=True)
        report["weak_escape"] = weak_escape_violations(
            n_max=weak_escape_n_max,
            i_max=weak_escape_i_max,
        )

    if run_adjacency:
        print("Paket E: adjacency ...", flush=True)
        report["adjacency"] = adjacency_profile(
            n_max=adjacency_n_max,
            i_min=adjacency_i_min,
            i_max=adjacency_i_max,
        )

    report["summary"] = compact_summary(report)
    report["sha256"] = sha256_of_object(report)
    return report


if __name__ == "__main__":
    report = run_all(
        run_i3_small=True,
        run_small_B=True,
        run_small_scan=True,
        run_weak_escape=True,
        run_adjacency=True,
        quick=True,  # False für Vollscan (kann Stunden dauern)
        weak_escape_n_max=300,
        weak_escape_i_max=20,
        adjacency_n_max=250,
        adjacency_i_min=10,
        adjacency_i_max=25,
    )

    write_json("erdos_report_full.json", report)
    write_json("erdos_report_summary.json", report["summary"])

    print("Fertig.")
    print("Vollbericht:   erdos_report_full.json")
    print("Kurzbericht:   erdos_report_summary.json")
    print("SHA256 full:   ", report["sha256"])
    print("SHA256 summary:", report["summary"]["sha256"])
    print(json.dumps(report["summary"], ensure_ascii=False, indent=2))