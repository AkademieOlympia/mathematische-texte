#!/usr/bin/env python3
"""
prime_cluster_evolution_tests.py

Exploratory tests for an "evolutionary" hierarchy of minimal admissible
prime constellations k=2..7.

It compares:
- occurrence counts up to N
- allowed primorial channel counts
- survival proxies
- local EABC transport metrics
- true clusters vs false k-1 clusters by local zirkulation

Run:
    python prime_cluster_evolution_tests.py --N 10000000 --W 30
"""

import argparse
import math
import numpy as np

STATE = {1: 0, 5: 1, 7: 2, 11: 3}  # E,A,B,C

PATTERNS = {
    "twin_k2_0_2": [0, 2],
    "triplet_k3_0_2_6": [0, 2, 6],
    "quad_k4_0_2_6_8": [0, 2, 6, 8],
    "quint_k5_0_4_6_10_12": [0, 4, 6, 10, 12],
    "sext_k6_0_4_6_10_12_16": [0, 4, 6, 10, 12, 16],
    "sept_k7_0_2_6_8_12_18_20": [0, 2, 6, 8, 12, 18, 20],
}

SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19]


def sieve(n: int) -> np.ndarray:
    arr = np.ones(n + 1, dtype=bool)
    arr[:2] = False
    arr[4::2] = False
    arr[2] = True
    m = int(math.isqrt(n))
    for i in range(3, m + 1, 2):
        if arr[i]:
            arr[i*i::i] = False
    return arr


def starts_pattern(isprime: np.ndarray, offsets, limit: int) -> np.ndarray:
    ok = np.ones(limit + 1, dtype=bool)
    for o in offsets:
        ok &= isprime[o:o + limit + 1]
    return np.nonzero(ok)[0]


def primorial(primes):
    M = 1
    for p in primes:
        M *= p
    return M


def allowed_residue_count(M: int, offsets):
    count = 0
    for r in range(M):
        good = True
        for q in SMALL_PRIMES:
            if M % q == 0 and any((r + o) % q == 0 for o in offsets):
                good = False
                break
        if good:
            count += 1
    return count


def eabc_states(prime_arr):
    return np.array([STATE[int(p % 12)] for p in prime_arr if int(p % 12) in STATE], dtype=np.int8)


def local_metrics(starts, offsets, W, prime_arr, state_arr, max_samples=None, seed=1):
    starts = np.asarray(starts, dtype=np.int64)
    if max_samples is not None and len(starts) > max_samples:
        rng = np.random.default_rng(seed)
        starts = rng.choice(starts, max_samples, replace=False)

    maxoff = max(offsets)
    rows = []
    for p in starts:
        lo, hi = int(p) - W, int(p) + maxoff + W
        i = np.searchsorted(prime_arr, lo, "left")
        j = np.searchsorted(prime_arr, hi, "right")
        ss = state_arr[i:j]
        if len(ss) < 2:
            continue

        M = np.zeros((4, 4), dtype=float)
        for a, b in zip(ss[:-1], ss[1:]):
            M[a, b] += 1.0
        total = M.sum()
        if total == 0:
            continue

        fwd = (M[3, 0] + M[0, 1] + M[1, 2] + M[2, 3]) / total
        bwd = (M[0, 3] + M[1, 0] + M[2, 1] + M[3, 2]) / total
        zirc = fwd - bwd
        diag = np.trace(M) / total

        row_sums = M.sum(axis=1)
        ent = 0.0
        for r in range(4):
            if row_sums[r] > 0:
                probs = M[r] / row_sums[r]
                ent += (row_sums[r] / total) * (-sum(x * math.log(x) for x in probs if x > 0))

        with np.errstate(invalid="ignore", divide="ignore"):
            T = M / row_sums[:, None]
        T = np.nan_to_num(T, nan=0.0)
        eigs = sorted([abs(x) for x in np.linalg.eigvals(T)], reverse=True)
        lam2 = eigs[1] if len(eigs) > 1 else 0.0
        rows.append((zirc, ent, lam2, diag, fwd, bwd, total))

    return np.array(rows)


def auc_score(pos, neg):
    pos = np.asarray(pos)
    neg = np.asarray(neg)
    vals = np.concatenate([pos, neg])
    order = np.argsort(vals)
    ranks = np.empty(len(vals), dtype=float)
    i = 0
    while i < len(vals):
        j = i + 1
        while j < len(vals) and vals[order[j]] == vals[order[i]]:
            j += 1
        ranks[order[i:j]] = (i + j - 1) / 2 + 1
        i = j
    n1, n0 = len(pos), len(neg)
    return (ranks[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10_000_000)
    ap.add_argument("--W", type=int, default=30)
    ap.add_argument("--sample", type=int, default=20000)
    args = ap.parse_args()

    N, W = args.N, args.W
    print(f"Sieve up to {N:,} ...")
    isprime = sieve(N + 128)
    primes = np.nonzero(isprime[:N + 128])[0]
    prime_arr = primes[primes > 3]
    state_arr = eabc_states(prime_arr)

    starts_by_name = {}

    print("\n=== Cluster counts and Hardy-Littlewood-style density proxy ===")
    print("name,k,span,count,x_over_logxk,count_over_proxy")
    for name, offsets in PATTERNS.items():
        starts = starts_pattern(isprime, offsets, N)
        starts_by_name[name] = starts
        k = len(offsets)
        proxy = N / (math.log(N) ** k)
        ratio = len(starts) / proxy if proxy else float("nan")
        print(f"{name},{k},{max(offsets)},{len(starts)},{proxy:.2f},{ratio:.4f}")

    print("\n=== Primorial channel counts ===")
    for name, offsets in PATTERNS.items():
        print(f"\n{name}")
        for mprimes in ([2,3,5], [2,3,5,7], [2,3,5,7,11], [2,3,5,7,11,13]):
            M = primorial(mprimes)
            c = allowed_residue_count(M, offsets)
            print(f"  M={M:<6d} allowed_channels={c:<7d} fraction={c/M:.8f}")

    print(f"\n=== Local EABC transport metrics, W={W} ===")
    print("name,k,count,zirc,entropy,lambda2,diag")
    for name, offsets in PATTERNS.items():
        starts = starts_by_name[name]
        met = local_metrics(starts, offsets, W, prime_arr, state_arr, max_samples=args.sample)
        if len(met) == 0:
            continue
        m = met.mean(axis=0)
        print(f"{name},{len(offsets)},{len(starts)},{m[0]:.5f},{m[1]:.5f},{m[2]:.5f},{m[3]:.5f}")

    print("\n=== True clusters vs false k-1 clusters by zirc ===")
    rng = np.random.default_rng(123)
    for name, offsets in PATTERNS.items():
        k = len(offsets)
        true_starts = starts_by_name[name]
        if len(true_starts) < 10:
            print(f"{name}: skipped, too few true starts")
            continue

        # sample broad candidate starts to keep runtime bounded
        sample_size = min(1_000_000, N + 1)
        cands = rng.choice(np.arange(0, N + 1), size=sample_size, replace=False)
        cands.sort()
        pc = np.zeros(len(cands), dtype=np.int16)
        valid = np.ones(len(cands), dtype=bool)
        for o in offsets:
            idx = cands + o
            valid &= idx < len(isprime)
            pc += isprime[idx]
        cands = cands[valid]
        pc = pc[valid]
        false = cands[pc == k-1]

        if len(false) < 10:
            print(f"{name}: skipped, false starts={len(false)}")
            continue

        pos = local_metrics(true_starts, offsets, W, prime_arr, state_arr, max_samples=args.sample)[:, 0]
        neg = local_metrics(false, offsets, W, prime_arr, state_arr, max_samples=args.sample, seed=321)[:, 0]
        print(f"{name}: true={len(true_starts)} false(k-1)={len(false)} AUC_zirc={auc_score(pos, neg):.3f}")


if __name__ == "__main__":
    main()
