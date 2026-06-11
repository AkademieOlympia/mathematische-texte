import math
import time

import numpy as np
import pandas as pd

# ============================================================
# BM-60 Twin-Isotropie-Test -- schnelle Version
# ============================================================

N_LIST = np.array([10_000_000, 20_000_000, 50_000_000, 100_000_000, 200_000_000], dtype=np.int64)
MAX_N = int(np.max(N_LIST))
MOD = 60
CHANNELS = np.array([11, 17, 29, 41, 47, 59], dtype=np.int64)
W_LIST = np.array([1000], dtype=np.int64)

PERMS = 2000
RANDOM_SEED = 42
OUTPUT_PREFIX = "bm60_twin_n_scaling_w1000"

np.random.seed(RANDOM_SEED)


def sieve_primes(n):
    print(f"[1] Erzeuge Primzahlsieb bis {n:,} ...")
    t0 = time.time()

    is_prime = np.ones(n + 3, dtype=bool)
    is_prime[:2] = False

    for i in range(2, int(math.isqrt(n + 2)) + 1):
        if is_prime[i]:
            is_prime[i*i:n+3:i] = False

    print(f"    fertig in {time.time() - t0:.2f}s")
    return is_prime


def find_twins(is_prime, n, chunk_size=5_000_000):
    print("[2] Suche Primzahlzwillinge ...")
    t0 = time.time()

    chunks = []
    for lo in range(3, n - 1, chunk_size):
        hi = min(n - 1, lo + chunk_size)
        p = np.arange(lo, hi, dtype=np.int64)
        mask = is_prime[p] & is_prime[p + 2]
        chunks.append(p[mask])

    starts = np.concatenate(chunks) if chunks else np.array([], dtype=np.int64)
    starts = starts[starts > MOD]

    print(f"    Anzahl Zwillinge nach Anfangsfilter: {len(starts):,}")
    print(f"    fertig in {time.time() - t0:.2f}s")
    return starts


def exact_random_baseline(qs, n, w):
    """
    Exakte Fensterwahrscheinlichkeit:
    Wahrscheinlichkeit, dass [x, x+w] mindestens einen Treffer aus qs enthält.
    """
    qs = np.sort(qs.astype(np.int64))

    if len(qs) == 0:
        return 0.0

    max_start = n - w
    if max_start <= 0:
        return 0.0

    a = np.maximum(0, qs - w)
    b = np.minimum(max_start, qs)

    valid = a <= b
    a = a[valid]
    b = b[valid]

    if len(a) == 0:
        return 0.0

    order = np.argsort(a)
    a = a[order]
    b = b[order]

    total = 0
    cur_a = int(a[0])
    cur_b = int(b[0])

    for ai, bi in zip(a[1:], b[1:]):
        ai = int(ai)
        bi = int(bi)

        if ai <= cur_b:
            if bi > cur_b:
                cur_b = bi
        else:
            total += cur_b - cur_a
            cur_a, cur_b = ai, bi

    total += cur_b - cur_a

    return total / max_start


def cluster_probs_for_qs(qs, w_list):
    """
    Clusterwahrscheinlichkeit P(next gap < W) für alle W.
    """
    if len(qs) < 3:
        return np.full(len(w_list), np.nan)

    gaps = np.diff(qs)

    return np.array([np.mean(gaps < w) for w in w_list], dtype=float)


def stats_for_labels(starts, labels, channels, w_list, baselines):
    """
    Berechnet Delta je Kanal und W.
    baselines[channel_index, W_index]
    """
    deltas = np.zeros((len(channels), len(w_list)), dtype=float)
    clusters = np.zeros_like(deltas)
    counts = np.zeros(len(channels), dtype=int)

    for i, r in enumerate(channels):
        qs = starts[labels == r]
        counts[i] = len(qs)

        c = cluster_probs_for_qs(qs, w_list)
        clusters[i, :] = c
        deltas[i, :] = c - baselines[i, :]

    return counts, clusters, deltas


def spread_by_w(deltas, channels):
    """
    Smax je W: max(delta_r)-min(delta_r)
    """
    smax = np.nanmax(deltas, axis=0) - np.nanmin(deltas, axis=0)

    hot_idx = np.nanargmax(deltas, axis=0)
    cold_idx = np.nanargmin(deltas, axis=0)

    hot = channels[hot_idx]
    cold = channels[cold_idx]

    hot_delta = deltas[hot_idx, np.arange(deltas.shape[1])]
    cold_delta = deltas[cold_idx, np.arange(deltas.shape[1])]

    return smax, hot, cold, hot_delta, cold_delta


def run_for_n(starts_all, n):
    starts = starts_all[starts_all < n - 1]
    labels = starts % MOD

    print("\n============================================================")
    print(f"[N] N = {n:,}")
    print("============================================================")
    print(f"    Zwillinge nach Anfangsfilter: {len(starts):,}")

    print("\n[3] Häufigkeiten je Kanal")
    for r in CHANNELS:
        print(f"    r={r:2d}: {np.sum(labels == r):,}")

    print("\n[4] Berechne exakte Baselines je Kanal und W ...")
    t0 = time.time()

    baselines = np.zeros((len(CHANNELS), len(W_LIST)), dtype=float)

    for i, r in enumerate(CHANNELS):
        qs = starts[labels == r]
        for j, w in enumerate(W_LIST):
            baselines[i, j] = exact_random_baseline(qs, n, int(w))

    print(f"    fertig in {time.time() - t0:.2f}s")

    print("\n[5] Echte Kanalstatistik")
    counts, clusters, deltas = stats_for_labels(
        starts, labels, CHANNELS, W_LIST, baselines
    )

    s_real, hot, cold, hot_delta, cold_delta = spread_by_w(deltas, CHANNELS)

    stats_rows = []
    for i, r in enumerate(CHANNELS):
        for j, w in enumerate(W_LIST):
            stats_rows.append({
                "N": int(n),
                "channel": int(r),
                "W": int(w),
                "count": int(counts[i]),
                "cluster_prob": clusters[i, j],
                "random_prob": baselines[i, j],
                "delta": deltas[i, j],
            })

    stats_df = pd.DataFrame(stats_rows)
    print(stats_df.to_string(index=False))

    print("\n[6] Permutationstest ...")
    t0 = time.time()

    perm_s = np.zeros((PERMS, len(W_LIST)), dtype=float)

    for k in range(PERMS):
        if (k + 1) % 100 == 0:
            print(f"    Permutation {k+1}/{PERMS}")

        perm_labels = np.random.permutation(labels)

        _, _, d_perm = stats_for_labels(
            starts, perm_labels, CHANNELS, W_LIST, baselines
        )

        perm_s[k, :], _, _, _, _ = spread_by_w(d_perm, CHANNELS)

    print(f"    fertig in {time.time() - t0:.2f}s")

    summary_rows = []

    for j, w in enumerate(W_LIST):
        mean_perm = float(np.mean(perm_s[:, j]))
        std_perm = float(np.std(perm_s[:, j]))
        D = float(s_real[j] - mean_perm)
        z = D / std_perm if std_perm > 0 else np.nan
        p = float(np.mean(perm_s[:, j] >= s_real[j]))
        channel_deltas = {
            f"Delta_{int(channel)}": float(deltas[i, j])
            for i, channel in enumerate(CHANNELS)
        }

        summary_rows.append({
            "N": int(n),
            "W": int(w),
            "S_real": float(s_real[j]),
            "mean_perm": mean_perm,
            "std_perm": std_perm,
            "D_excess": D,
            "z": z,
            "p_value": p,
            "hot_channel": int(hot[j]),
            "hot_delta": float(hot_delta[j]),
            "cold_channel": int(cold[j]),
            "cold_delta": float(cold_delta[j]),
            **channel_deltas,
        })

    summary_df = pd.DataFrame(summary_rows)

    print("\n[7] Look-elsewhere Summary")
    print(summary_df.to_string(index=False))

    perm_df = pd.DataFrame(perm_s, columns=[f"W_{int(w)}" for w in W_LIST])
    perm_df.insert(0, "N", int(n))

    return stats_df, summary_df, perm_df


def main():
    print("============================================================")
    print("BM-60 Twin-Isotropie-Test -- N-Skalierung W=1000")
    print("============================================================")
    print(f"N_LIST = {[int(n) for n in N_LIST]}")
    print(f"MAX_N  = {MAX_N:,}")
    print(f"PERMS  = {PERMS}")
    print(f"W_LIST = {[int(w) for w in W_LIST]}")
    print("============================================================")

    is_prime = sieve_primes(MAX_N)
    starts_all = find_twins(is_prime, MAX_N)

    all_stats = []
    all_summaries = []
    all_perms = []

    for n in N_LIST:
        stats_df, summary_df, perm_df = run_for_n(starts_all, int(n))
        all_stats.append(stats_df)
        all_summaries.append(summary_df)
        all_perms.append(perm_df)

    stats_all_df = pd.concat(all_stats, ignore_index=True)
    summary_all_df = pd.concat(all_summaries, ignore_index=True)
    perm_all_df = pd.concat(all_perms, ignore_index=True)

    stats_all_df.to_csv(f"{OUTPUT_PREFIX}_channel_stats.csv", index=False)
    summary_all_df.to_csv(f"{OUTPUT_PREFIX}_look_elsewhere_summary.csv", index=False)
    perm_all_df.to_csv(f"{OUTPUT_PREFIX}_perm_smax.csv", index=False)

    print("\nCSV-Dateien geschrieben:")
    print(f"  {OUTPUT_PREFIX}_channel_stats.csv")
    print(f"  {OUTPUT_PREFIX}_look_elsewhere_summary.csv")
    print(f"  {OUTPUT_PREFIX}_perm_smax.csv")


if __name__ == "__main__":
    main()
