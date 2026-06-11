import math
import time

import numpy as np

N = 100_000_000
MOD = 60
W = 800
PERMS = 5000
CHANNELS = np.array([11, 17, 29, 41, 47, 59], dtype=np.int64)

np.random.seed(42)


def fast_sieve(n: int) -> np.ndarray:
    is_prime = np.ones(n + 3, dtype=bool)
    is_prime[:2] = False
    for i in range(2, int(math.isqrt(n + 2)) + 1):
        if is_prime[i]:
            is_prime[i * i : n + 3 : i] = False
    return is_prime


def get_twins(is_prime: np.ndarray, n: int) -> np.ndarray:
    p = np.where(is_prime[3 : n - 1] & is_prime[5 : n + 1])[0] + 3
    return p[p > MOD].astype(np.int64)


def exact_random_baseline(qs: np.ndarray, n: int, w: int) -> float:
    if len(qs) == 0:
        return 0.0

    max_start = n - w
    if max_start <= 0:
        return 0.0

    qs = np.sort(qs.astype(np.int64))
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
            cur_b = max(cur_b, bi)
        else:
            total += cur_b - cur_a
            cur_a, cur_b = ai, bi

    total += cur_b - cur_a
    return total / max_start


def cluster_prob(qs: np.ndarray, w: int) -> float:
    if len(qs) < 3:
        return float("nan")
    return float(np.mean(np.diff(qs) < w))


def delta_for_channel(starts, labels, r, baseline, w: int) -> float:
    qs = starts[labels == r]
    return cluster_prob(qs, w) - baseline


def main() -> None:
    t0 = time.time()
    print("Sieb und Zwillinge...")
    is_prime = fast_sieve(N)
    starts = get_twins(is_prime, N)
    labels = starts % MOD
    print("Zwillinge:", len(starts), "Zeit:", f"{time.time() - t0:.2f}s")
    _ = CHANNELS  # stellt sicher, dass der Kanal-Array geladen ist

    # Baselines für echte Kanäle einmal berechnen
    baseline_11 = exact_random_baseline(starts[labels == 11], N, W)
    baseline_47 = exact_random_baseline(starts[labels == 47], N, W)

    d11 = delta_for_channel(starts, labels, 11, baseline_11, W)
    d47 = delta_for_channel(starts, labels, 47, baseline_47, W)
    p_real = d47 - d11

    print("Realer Befund:")
    print("Delta_11 =", d11)
    print("Delta_47 =", d47)
    print("P_11_47 =", p_real)

    t1 = time.time()
    perm_vals = []
    for k in range(PERMS):
        if (k + 1) % 500 == 0:
            print("Permutation", k + 1, "/", PERMS)

        perm_labels = np.random.permutation(labels)
        # gleiche Baselines wie real: Test auf Labelstruktur
        pd11 = delta_for_channel(starts, perm_labels, 11, baseline_11, W)
        pd47 = delta_for_channel(starts, perm_labels, 47, baseline_47, W)

        perm_vals.append(pd47 - pd11)

    print("Permutationen fertig in", f"{time.time() - t1:.2f}s")

    perm_vals = np.array(perm_vals)

    mean_perm = float(np.mean(perm_vals))
    std_perm = float(np.std(perm_vals))
    z = (p_real - mean_perm) / std_perm if std_perm > 0 else float("nan")
    p_value = float(np.mean(perm_vals >= p_real))

    print()
    print("Permutationstest P_11_47 bei W=800")
    print("PERMS =", PERMS)
    print("mean_perm =", mean_perm)
    print("std_perm =", std_perm)
    print("P_real =", p_real)
    print("z =", z)
    print("p_value =", p_value)

    print()
    print("Quantile:")
    for q in [0.5, 0.9, 0.95, 0.975, 0.99]:
        print(q, float(np.quantile(perm_vals, q)))


if __name__ == "__main__":
    main()
