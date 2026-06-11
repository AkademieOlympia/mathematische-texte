#!/usr/bin/env python3
"""
cluster_evolution_graph.py

Builds an exploratory "evolution graph" for minimal admissible prime clusters.

Goal:
- Nodes = cluster patterns / EABC motifs
- Edges = extension / recombination / containment relations
- Measures:
    * occurrence counts up to N
    * pre-field EABC zirkulation before cluster start
    * around-field zirkulation around cluster
    * child/parent conditional rates
    * false-child comparison
    * graph export as CSV

Run:
    python cluster_evolution_graph.py --N 10000000 --W 30

Suggested later:
    python cluster_evolution_graph.py --N 100000000 --W 30
"""

import argparse
import csv
import math
from pathlib import Path
from collections import defaultdict
import numpy as np


STATE = {1: 0, 5: 1, 7: 2, 11: 3}
STATE_NAME = ["E", "A", "B", "C"]


# A small hierarchy of compact admissible prime constellations.
# You can add your own patterns here.
PATTERNS = {
    "twin_02": [0, 2],
    "cousin_04": [0, 4],
    "sexy_06": [0, 6],

    "triplet_026": [0, 2, 6],
    "triplet_046": [0, 4, 6],

    "quad_0268": [0, 2, 6, 8],

    "quint_0461012": [0, 4, 6, 10, 12],

    "sext_046101216": [0, 4, 6, 10, 12, 16],

    "sept_0268121820": [0, 2, 6, 8, 12, 18, 20],
}


def sieve(n: int) -> np.ndarray:
    """Fast enough Eratosthenes sieve for N up to around 1e8 on a normal desktop."""
    arr = np.ones(n + 1, dtype=bool)
    arr[:2] = False
    arr[4::2] = False
    arr[2] = True
    m = math.isqrt(n)
    for i in range(3, m + 1, 2):
        if arr[i]:
            arr[i*i::i] = False
    return arr


def starts_pattern(isprime: np.ndarray, offsets, limit: int) -> np.ndarray:
    ok = np.ones(limit + 1, dtype=bool)
    for o in offsets:
        ok &= isprime[o:o + limit + 1]
    return np.nonzero(ok)[0]


def eabc_motif_from_offsets(p: int, offsets) -> str:
    out = []
    for o in offsets:
        r = (p + o) % 12
        out.append(STATE_NAME[STATE[r]])
    return "".join(out)


def local_transport_metrics(prime_arr, state_arr, lo: int, hi: int):
    """
    Build local 4x4 EABC transition matrix from primes in [lo, hi].
    Returns zirkulation, entropy, lambda2, diag, total_transitions.
    """
    i = np.searchsorted(prime_arr, lo, side="left")
    j = np.searchsorted(prime_arr, hi, side="right")
    ss = state_arr[i:j]

    if len(ss) < 2:
        return None

    M = np.zeros((4, 4), dtype=float)
    for a, b in zip(ss[:-1], ss[1:]):
        M[a, b] += 1.0

    total = M.sum()
    if total <= 0:
        return None

    # Forward cycle from earlier EABC tests:
    # C -> E -> A -> B -> C
    fwd = (M[3, 0] + M[0, 1] + M[1, 2] + M[2, 3]) / total
    bwd = (M[0, 3] + M[1, 0] + M[2, 1] + M[3, 2]) / total
    zirc = fwd - bwd

    diag = np.trace(M) / total

    row_sums = M.sum(axis=1)
    entropy = 0.0
    for r in range(4):
        if row_sums[r] > 0:
            probs = M[r] / row_sums[r]
            entropy += (row_sums[r] / total) * (-sum(x * math.log(x) for x in probs if x > 0))

    with np.errstate(invalid="ignore", divide="ignore"):
        T = M / row_sums[:, None]
    T = np.nan_to_num(T, nan=0.0, posinf=0.0, neginf=0.0)

    eigs = sorted([abs(x) for x in np.linalg.eigvals(T)], reverse=True)
    lam2 = eigs[1] if len(eigs) > 1 else 0.0

    return {
        "zirc": zirc,
        "entropy": entropy,
        "lambda2": lam2,
        "diag": diag,
        "transitions": total,
        "fwd": fwd,
        "bwd": bwd,
    }


def mean_metrics_for_starts(starts, offsets, W, prime_arr, state_arr, mode="around", max_samples=None, seed=1):
    rng = np.random.default_rng(seed)
    starts = np.asarray(starts, dtype=np.int64)

    if max_samples is not None and len(starts) > max_samples:
        starts = rng.choice(starts, size=max_samples, replace=False)

    rows = []
    maxoff = max(offsets)

    for p in starts:
        if mode == "around":
            lo, hi = int(p) - W, int(p) + maxoff + W
        elif mode == "pre":
            lo, hi = int(p) - W, int(p) - 1
        elif mode == "post":
            lo, hi = int(p) + maxoff + 1, int(p) + maxoff + W
        else:
            raise ValueError("mode must be around, pre, or post")

        m = local_transport_metrics(prime_arr, state_arr, lo, hi)
        if m is not None:
            rows.append(m)

    if not rows:
        return {}

    keys = rows[0].keys()
    return {k: float(np.mean([r[k] for r in rows])) for k in keys}


def contains_pattern(parent_offsets, child_offsets):
    """
    Does parent contain child after translation?
    Returns list of shifts s such that {s+child_offsets} subset of parent_offsets.
    """
    P = set(parent_offsets)
    shifts = []
    for po in parent_offsets:
        for co in child_offsets:
            s = po - co
            if all((s + x) in P for x in child_offsets):
                shifts.append(s)
    return sorted(set(shifts))


def extension_edges(patterns):
    """
    Edges child -> parent if parent contains translated child
    and len(parent)=len(child)+1 or more.
    """
    edges = []
    names = list(patterns.keys())
    for a in names:
        for b in names:
            if a == b:
                continue
            ca, cb = patterns[a], patterns[b]
            if len(cb) <= len(ca):
                continue
            shifts = contains_pattern(cb, ca)
            if shifts:
                edges.append((a, b, shifts))
    return edges


def conditional_extension_rate(parent_starts, child_starts, shifts):
    """
    For child->parent containment. If parent contains child shifted by s,
    then a child occurrence at p corresponds to a parent at p-s.
    We measure fraction of child starts that extend to a parent.
    """
    parent_set = set(int(x) for x in parent_starts)
    child_starts = np.asarray(child_starts, dtype=np.int64)
    hits = 0
    for p in child_starts:
        if any((int(p) - s) in parent_set for s in shifts):
            hits += 1
    return hits / len(child_starts) if len(child_starts) else 0.0, hits


def false_kminus1_starts(isprime, offsets, N, max_count=200000, seed=123):
    """
    Candidate starts where exactly k-1 offsets are prime.
    This gives near-miss controls.
    """
    rng = np.random.default_rng(seed)
    cands = np.arange(0, N + 1, dtype=np.int64)
    if len(cands) > max_count:
        cands = rng.choice(cands, size=max_count, replace=False)
        cands.sort()

    pc = np.zeros(len(cands), dtype=np.int16)
    valid = np.ones(len(cands), dtype=bool)
    for o in offsets:
        idx = cands + o
        valid &= idx < len(isprime)
        pc += isprime[idx]
    cands = cands[valid]
    pc = pc[valid]
    return cands[pc == len(offsets)-1]


def auc_score(pos, neg):
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
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


def metric_values(starts, offsets, W, prime_arr, state_arr, key="zirc", mode="around", max_samples=None, seed=1):
    rng = np.random.default_rng(seed)
    starts = np.asarray(starts, dtype=np.int64)
    if max_samples is not None and len(starts) > max_samples:
        starts = rng.choice(starts, size=max_samples, replace=False)

    vals = []
    maxoff = max(offsets)
    for p in starts:
        if mode == "around":
            lo, hi = int(p) - W, int(p) + maxoff + W
        elif mode == "pre":
            lo, hi = int(p) - W, int(p) - 1
        else:
            lo, hi = int(p) + maxoff + 1, int(p) + maxoff + W

        m = local_transport_metrics(prime_arr, state_arr, lo, hi)
        if m is not None:
            vals.append(m[key])
    return np.asarray(vals, dtype=float)


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10_000_000)
    ap.add_argument("--W", type=int, default=30)
    ap.add_argument("--outdir", type=str, default="cluster_evolution_out")
    ap.add_argument("--sample", type=int, default=20000)
    args = ap.parse_args()

    N, W = args.N, args.W
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Sieve up to {N:,} ...")
    isprime = sieve(N + 256)

    primes = np.nonzero(isprime[:N + 256])[0]
    prime_arr = primes[primes > 3]
    state_arr = np.array([STATE[int(p % 12)] for p in prime_arr], dtype=np.int8)

    # Occurrences
    starts = {}
    for name, offsets in PATTERNS.items():
        starts[name] = starts_pattern(isprime, offsets, N)

    # Node table
    node_rows = []
    print("\n=== Nodes / cluster motifs ===")
    for name, offsets in PATTERNS.items():
        st = starts[name]
        around = mean_metrics_for_starts(st, offsets, W, prime_arr, state_arr, "around", args.sample)
        pre = mean_metrics_for_starts(st, offsets, W, prime_arr, state_arr, "pre", args.sample)
        post = mean_metrics_for_starts(st, offsets, W, prime_arr, state_arr, "post", args.sample)

        motif_examples = []
        for p in st[:5]:
            try:
                motif_examples.append(eabc_motif_from_offsets(int(p), offsets))
            except KeyError:
                motif_examples.append("?")
        motif_set = ",".join(sorted(set(motif_examples)))

        row = {
            "name": name,
            "k": len(offsets),
            "span": max(offsets),
            "count": len(st),
            "density": len(st) / N,
            "eabc_examples": motif_set,
            "around_zirc": around.get("zirc", float("nan")),
            "around_entropy": around.get("entropy", float("nan")),
            "around_lambda2": around.get("lambda2", float("nan")),
            "pre_zirc": pre.get("zirc", float("nan")),
            "post_zirc": post.get("zirc", float("nan")),
        }
        node_rows.append(row)

        print(
            f"{name:24s} k={row['k']} count={row['count']:8d} "
            f"around_zirc={row['around_zirc']:.4f} pre_zirc={row['pre_zirc']:.4f} "
            f"motifs={motif_set}"
        )

    write_csv(
        outdir / "nodes.csv",
        node_rows,
        ["name","k","span","count","density","eabc_examples",
         "around_zirc","around_entropy","around_lambda2","pre_zirc","post_zirc"]
    )

    # Extension graph
    print("\n=== Extension / containment edges ===")
    edge_rows = []
    for child, parent, shifts in extension_edges(PATTERNS):
        rate, hits = conditional_extension_rate(starts[parent], starts[child], shifts)
        row = {
            "child": child,
            "parent": parent,
            "shifts": ";".join(map(str, shifts)),
            "child_count": len(starts[child]),
            "parent_count": len(starts[parent]),
            "hits": hits,
            "conditional_rate": rate,
        }
        edge_rows.append(row)
        print(
            f"{child:24s} -> {parent:24s} "
            f"shifts={shifts} hits={hits:6d}/{len(starts[child]):6d} "
            f"rate={rate:.6f}"
        )

    write_csv(
        outdir / "edges.csv",
        edge_rows,
        ["child","parent","shifts","child_count","parent_count","hits","conditional_rate"]
    )

    # True-vs-false zirkulation controls
    print("\n=== True vs false near-miss controls ===")
    control_rows = []
    for name, offsets in PATTERNS.items():
        true_st = starts[name]
        if len(true_st) < 5:
            continue
        false_st = false_kminus1_starts(isprime, offsets, N, max_count=500000)
        if len(false_st) < 5:
            continue

        pos = metric_values(true_st, offsets, W, prime_arr, state_arr, "zirc", "around", args.sample)
        neg = metric_values(false_st, offsets, W, prime_arr, state_arr, "zirc", "around", args.sample)
        auc = auc_score(pos, neg)

        row = {
            "name": name,
            "true_count": len(true_st),
            "false_count_sampled": len(false_st),
            "auc_zirc": auc,
            "true_zirc_mean": float(np.mean(pos)),
            "false_zirc_mean": float(np.mean(neg)),
        }
        control_rows.append(row)

        print(
            f"{name:24s} true={len(true_st):8d} false={len(false_st):8d} "
            f"AUC_zirc={auc:.3f} mean_true={row['true_zirc_mean']:.4f} "
            f"mean_false={row['false_zirc_mean']:.4f}"
        )

    write_csv(
        outdir / "controls.csv",
        control_rows,
        ["name","true_count","false_count_sampled","auc_zirc","true_zirc_mean","false_zirc_mean"]
    )

    print(f"\nWrote CSV files to: {outdir.resolve()}")
    print("\nInterpretation guide:")
    print("  nodes.csv    = motif activity and local transport metrics")
    print("  edges.csv    = extension / containment graph")
    print("  controls.csv = true clusters vs near-miss false clusters")
    print("\nCore hypothesis to check:")
    print("  Higher-k clusters should become rarer, but their local around_zirc should rise.")
    print("  If true-vs-false AUC_zirc stays > 0.5, local loop coherence separates real clusters.")


if __name__ == "__main__":
    main()
