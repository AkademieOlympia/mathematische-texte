#!/usr/bin/env python3
"""EABC-Vierlinge bis X: W_E, R_beta und mod-420-Diagnostik.

Stufen 0–3 (vgl. collatz_eabc_zirkulationshypothese.md §4.2):
  Stufe 0 — D=A−C, Q=A+C, W_E=D/Q, R_beta=D/Q^beta (Definition, keine Exponentenannahme)
  Stufe 1 — alpha_eff = d log|D|/d log Q (numerisches Signal, kein Satz)
  Stufe 2 — alpha_E = inf{beta : R_beta beschränkt} (Vermutung, erst wenn Daten es nahelegen)
  Stufe 3 — W_E Grenzwert ≠0? (Orientierung/Holonomie am Ende)

Kein Z_E = D/sqrt(Q) auf Definitionsebene — R_1_2 ist Heuristik/Diagnose-Alias.

CSV-Spalten R_beta:
  R_1_2  = R_{1/2}(X) = D(X)/sqrt(Q(X))  (Heuristik/Diagnose; Alias Z_E)
  R_2_3, R_3_4, R_9_10 -> Zwischennormierungen (beta = 2/3, 3/4, 0.9)
  R_1    = R_1(X) = D(X)/Q(X) = W_E(X)  (Stufe 3: Orientierung)
"""

import argparse
import csv
import math
import time
from pathlib import Path

import numpy as np

REGULAR_420 = [11, 101, 191, 221, 311, 401]
R_BETA_COLUMNS = ("R_1_2", "R_2_3", "R_3_4", "R_9_10", "R_1")


def small_primes_upto(n: int):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            sieve[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    return [p for p in range(2, n + 1) if sieve[p]]


def odd_segment_prime_flags(L: int, R: int, base_primes):
    """
    Prime flags for odd numbers in [L,R].
    L and R must be odd.
    index i represents n = L + 2*i.
    """
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


def default_checkpoints(X: int):
    cps = []
    for e in range(6, int(math.log10(X)) + 1):
        for m in [1, 2, 3, 5]:
            v = m * 10**e
            if v <= X:
                cps.append(v)
    if X not in cps:
        cps.append(X)
    return sorted(set(cps))


def emit_row(writer, X, total, abce, ceab, counts420, elapsed):
    denom = abce + ceab
    diff = abce - ceab
    R_1_2 = diff / (denom ** 0.5) if denom else 0.0
    R_2_3 = diff / (denom ** (2 / 3)) if denom else 0.0
    R_3_4 = diff / (denom ** 0.75) if denom else 0.0
    R_9_10 = diff / (denom ** 0.9) if denom else 0.0
    R_1 = diff / denom if denom else 0.0
    W_E = R_1
    Z_E = R_1_2
    sqrt_Q = math.sqrt(denom) if denom else 0.0
    absZ_E = abs(Z_E)
    log_Q = math.log(denom) if denom > 0 else ""
    log_absdiff = math.log(abs(diff)) if diff != 0 else ""
    alpha_eff = (log_absdiff / log_Q) if (log_Q and log_absdiff) else ""

    regular_counts = np.array([counts420[r] for r in REGULAR_420], dtype=np.float64)
    mean420 = regular_counts.mean() if regular_counts.sum() else 0.0
    W420 = ((regular_counts.max() - regular_counts.min()) / mean420) if mean420 else 0.0
    chi420 = 0.0
    if mean420:
        chi420 = float(((regular_counts - mean420) ** 2 / mean420).sum())

    row = {
        "X": X,
        "Q_total": total,
        "ABCE": abce,
        "CEAB": ceab,
        "diff": diff,
        "W_E": W_E,
        "Z_E": Z_E,
        "sqrt_Q": sqrt_Q,
        "absZ_E": absZ_E,
        "log_Q": log_Q,
        "log_absdiff": log_absdiff,
        "alpha_eff": alpha_eff,
        "R_1_2": R_1_2,
        "R_2_3": R_2_3,
        "R_3_4": R_3_4,
        "R_9_10": R_9_10,
        "R_1": R_1,
        "W420_range": W420,
        "chi2_420_df5": chi420,
        "mod420_5": counts420[5],
        "mod420_11": counts420[11],
        "mod420_101": counts420[101],
        "mod420_191": counts420[191],
        "mod420_221": counts420[221],
        "mod420_311": counts420[311],
        "mod420_401": counts420[401],
        "elapsed_sec": elapsed,
    }
    writer.writerow(row)
    print(
        f"X={X:>12}  Q={total:>8}  D={diff:+5}  "
        f"W_E={W_E:+.6e}  Z_E={Z_E:+.3f}  |Z_E|={absZ_E:.3f}  "
        f"W420={W420:.6e}  chi2_420={chi420:.3f}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--X", type=int, default=10**10)
    parser.add_argument(
        "--block",
        type=int,
        default=50_000_000,
        help="integer width per block; increase if RAM allows",
    )
    parser.add_argument("--out", type=str, default="eabc_quadruplets.csv")
    args = parser.parse_args()

    X = args.X
    block = args.block

    root = int(math.isqrt(X + 8)) + 1
    base_primes = small_primes_upto(root)

    checkpoints = default_checkpoints(X)
    cp_idx = 0

    total = 0
    abce = 0
    ceab = 0
    counts420 = np.zeros(420, dtype=np.int64)

    out_path = Path(args.out)
    start_time = time.time()

    fieldnames = [
        "X",
        "Q_total",
        "ABCE",
        "CEAB",
        "diff",
        "W_E",
        "Z_E",
        "sqrt_Q",
        "absZ_E",
        "log_Q",
        "log_absdiff",
        "alpha_eff",
        *R_BETA_COLUMNS,
        "W420_range",
        "chi2_420_df5",
        "mod420_5",
        "mod420_11",
        "mod420_101",
        "mod420_191",
        "mod420_221",
        "mod420_311",
        "mod420_401",
        "elapsed_sec",
    ]

    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        L = 5
        if L % 2 == 0:
            L += 1

        while L <= X:
            next_cp = checkpoints[cp_idx]
            H = min(X, next_cp, L + block - 1)
            if H % 2 == 0:
                H -= 1
            if H < L:
                cp_idx += 1
                continue

            R = H + 8
            if R % 2 == 0:
                R += 1

            flags = odd_segment_prime_flags(L, R, base_primes)
            arr = np.frombuffer(flags, dtype=np.uint8)
            npos = (H - L) // 2 + 1

            mask = arr[0:npos] & arr[1 : npos + 1] & arr[3 : npos + 3] & arr[4 : npos + 4]
            idx = np.flatnonzero(mask)
            if idx.size:
                p = L + 2 * idx.astype(np.int64)
                total += int(idx.size)

                residues12 = p % 12
                abce += int(np.count_nonzero(residues12 == 5))
                ceab += int(np.count_nonzero(residues12 == 11))

                residues420 = p % 420
                counts420 += np.bincount(residues420, minlength=420)

            if H >= next_cp or H + 1 >= next_cp:
                elapsed = time.time() - start_time
                emit_row(writer, next_cp, total, abce, ceab, counts420, elapsed)
                f.flush()
                cp_idx += 1
                if cp_idx >= len(checkpoints):
                    break

            L = H + 2

    print(f"\nCSV gespeichert unter: {out_path.resolve()}")


if __name__ == "__main__":
    main()
