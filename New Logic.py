from math import isqrt
from statistics import mean, median

def sieve_primes(n):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            step = p
            start = p * p
            sieve[start:n+1:step] = b"\x00" * (((n - start) // step) + 1)
    return sieve

def is_prime_quad_start(p, prime):
    return prime[p] and prime[p+2] and prime[p+6] and prime[p+8]

def find_prime_quadruplets(limit):
    prime = sieve_primes(limit + 8)
    quads = []
    for p in range(5, limit + 1):
        if p % 12 in (5, 11):
            if is_prime_quad_start(p, prime):
                quads.append(p)
    return quads, prime

def small_sieve_candidate(p, small_primes):
    if p % 12 not in (5, 11):
        return False
    for q in small_primes:
        if q in (2, 3):
            continue
        residues = {(-d) % q for d in (0, 2, 6, 8)}
        if p % q in residues:
            return False
    return True

def next_compatible_candidates(start, stop, small_primes):
    return [
        p for p in range(start, stop + 1)
        if small_sieve_candidate(p, small_primes)
    ]

def eabc_class(p):
    r = p % 12
    if r == 1:
        return "E"
    if r == 5:
        return "A"
    if r == 7:
        return "B"
    if r == 11:
        return "C"
    return "X"

def transition_signature(p_now, p_next):
    gap = p_next - p_now
    return {
        "from": eabc_class(p_now),
        "to": eabc_class(p_next),
        "gap_mod_12": gap % 12,
        "gap_mod_5": gap % 5,
        "gap_mod_7": gap % 7,
        "gap_mod_11": gap % 11,
    }

def build_transition_memory(quads):
    memory = {}
    for i in range(len(quads) - 2):
        sig = transition_signature(quads[i], quads[i + 1])
        key = (
            sig["from"],
            sig["to"],
            sig["gap_mod_12"],
            sig["gap_mod_5"],
            sig["gap_mod_7"],
        )
        next_gap = quads[i + 2] - quads[i + 1]
        memory.setdefault(key, []).append(next_gap)
    return memory

def eabc_weighted_gap_prediction(current_index, quads, memory, fallback_window=20):
    if current_index < 1:
        return 1000
    sig = transition_signature(quads[current_index - 1], quads[current_index])
    key = (
        sig["from"],
        sig["to"],
        sig["gap_mod_12"],
        sig["gap_mod_5"],
        sig["gap_mod_7"],
    )
    if key in memory and len(memory[key]) >= 2:
        return int(median(memory[key]))
    gaps = [
        quads[i + 1] - quads[i]
        for i in range(max(0, current_index - fallback_window), current_index)
    ]
    return int(median(gaps)) if gaps else 1000

def quat_orientation(gap):
    r = gap % 12
    if r == 0:
        return "E"
    if r in (2, 10):
        return "i"
    if r in (4, 8):
        return "j"
    if r == 6:
        return "k"
    return "X"

def dynamical_signature(quads, n, window=10):
    if n < window + 2:
        return None
    gaps = [
        quads[i + 1] - quads[i]
        for i in range(n - window, n)
    ]
    acc = [
        gaps[i + 1] - gaps[i]
        for i in range(len(gaps) - 1)
    ]
    jerk = [
        acc[i + 1] - acc[i]
        for i in range(len(acc) - 1)
    ]
    avg_gap = mean(gaps)
    holonomy_energy = sum(
        (g - avg_gap) ** 2 for g in gaps
    )
    return {
        "gap_mean": avg_gap,
        "gap_std_energy": holonomy_energy,
        "last_gap": gaps[-1],
        "last_acc": acc[-1] if acc else 0,
        "last_jerk": jerk[-1] if jerk else 0,
        "orientation": quat_orientation(gaps[-1]),
    }

def weighted_target(current, signature):
    drift = (
        0.60 * signature["gap_mean"]
        + 0.25 * signature["last_gap"]
        + 0.10 * signature["last_acc"]
        + 0.05 * signature["last_jerk"]
    )
    return int(current + drift)

def baseline_target(quads, current_index, window):
    gaps = [
        quads[i + 1] - quads[i]
        for i in range(max(0, current_index - window), current_index)
    ]
    estimated_gap = int(median(gaps)) if gaps else 1000
    return quads[current_index] + estimated_gap

def pearson_corr(x_values, y_values):
    if not x_values or not y_values or len(x_values) != len(y_values):
        return None
    mx = mean(x_values)
    my = mean(y_values)
    cov = sum((x - mx) * (y - my) for x, y in zip(x_values, y_values))
    var_x = sum((x - mx) ** 2 for x in x_values)
    var_y = sum((y - my) ** 2 for y in y_values)
    if var_x == 0 or var_y == 0:
        return None
    return cov / ((var_x * var_y) ** 0.5)

def average_ranks(values):
    indexed = sorted(enumerate(values), key=lambda t: t[1])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i + 1
        while j < len(indexed) and indexed[j][1] == indexed[i][1]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            original_idx = indexed[k][0]
            ranks[original_idx] = avg_rank
        i = j
    return ranks

def spearman_corr(x_values, y_values):
    if not x_values or not y_values or len(x_values) != len(y_values):
        return None
    rx = average_ranks(x_values)
    ry = average_ranks(y_values)
    return pearson_corr(rx, ry)

def orientation_summary(results):
    groups = {}
    for r in results:
        orient = r.get("orientation")
        dist = r.get("candidate_distance")
        if orient is None or dist is None:
            continue
        groups.setdefault(orient, []).append(dist)
    summary = []
    for orient, dists in sorted(groups.items()):
        top10 = sum(d <= 10 for d in dists)
        summary.append({
            "orientation": orient,
            "count": len(dists),
            "median_dist": sorted(dists)[len(dists) // 2],
            "top10_rate": top10 / len(dists),
        })
    return summary

def holonomy_threshold_from_quads(quads, window, quantile=0.5):
    energies = []
    for i in range(window + 2, len(quads) - 1):
        sig = dynamical_signature(quads, i, window=window)
        if sig is not None:
            energies.append(sig["gap_std_energy"])
    if not energies:
        return None
    energies.sort()
    idx = int((len(energies) - 1) * quantile)
    return energies[idx]

def predict_next_quad_dynamical(
    current_index,
    quads,
    small_primes,
    window=10,
    search_radius=50000,
):
    current = quads[current_index]
    true_next = quads[current_index + 1]
    signature = dynamical_signature(quads, current_index, window=window)
    if signature is None:
        return None
    target = weighted_target(current, signature)

    candidates = next_compatible_candidates(
        max(5, target - search_radius),
        target + search_radius,
        small_primes
    )

    if not candidates:
        return None

    prediction = min(candidates, key=lambda x: abs(x - target))

    try:
        prediction_index = candidates.index(prediction)
        true_index = candidates.index(true_next)
        candidate_distance = abs(true_index - prediction_index)
    except ValueError:
        prediction_index = None
        true_index = None
        candidate_distance = None

    return {
        "current_quad_start": current,
        "true_next_quad_start": true_next,
        "target": target,
        "prediction": prediction,
        "prediction_error": prediction - true_next,
        "absolute_error": abs(prediction - true_next),
        "prediction_index": prediction_index,
        "true_index": true_index,
        "candidate_distance": candidate_distance,
        "candidate_count": len(candidates),
        "holonomy_energy": signature["gap_std_energy"],
        "gap_mean": signature["gap_mean"],
        "last_gap": signature["last_gap"],
        "last_acc": signature["last_acc"],
        "last_jerk": signature["last_jerk"],
        "orientation": signature["orientation"],
    }

def predict_next_quad_baseline(
    current_index,
    quads,
    small_primes,
    window=10,
    search_radius=50000,
):
    current = quads[current_index]
    true_next = quads[current_index + 1]
    target = baseline_target(quads, current_index, window)

    candidates = next_compatible_candidates(
        max(5, target - search_radius),
        target + search_radius,
        small_primes
    )
    if not candidates:
        return None

    prediction = min(candidates, key=lambda x: abs(x - target))
    try:
        prediction_index = candidates.index(prediction)
        true_index = candidates.index(true_next)
        candidate_distance = abs(true_index - prediction_index)
    except ValueError:
        prediction_index = None
        true_index = None
        candidate_distance = None

    return {
        "target": target,
        "prediction": prediction,
        "candidate_distance": candidate_distance,
        "prediction_index": prediction_index,
        "true_index": true_index,
    }

def predict_next_quad_hybrid(
    current_index,
    quads,
    small_primes,
    holonomy_threshold,
    window=10,
    search_radius=50000,
):
    signature = dynamical_signature(quads, current_index, window=window)
    if signature is None or holonomy_threshold is None:
        return predict_next_quad_baseline(
            current_index=current_index,
            quads=quads,
            small_primes=small_primes,
            window=window,
            search_radius=search_radius,
        )
    if signature["gap_std_energy"] < holonomy_threshold:
        return predict_next_quad_dynamical(
            current_index=current_index,
            quads=quads,
            small_primes=small_primes,
            window=window,
            search_radius=search_radius,
        )
    return predict_next_quad_baseline(
        current_index=current_index,
        quads=quads,
        small_primes=small_primes,
        window=window,
        search_radius=search_radius,
    )

def holonomy_bin_summary(results, bins=4):
    pairs = [
        (r["holonomy_energy"], r["candidate_distance"])
        for r in results
        if r["holonomy_energy"] is not None and r["candidate_distance"] is not None
    ]
    if not pairs:
        return []
    pairs.sort(key=lambda x: x[0])
    n = len(pairs)
    summaries = []
    for b in range(bins):
        start = (b * n) // bins
        end = ((b + 1) * n) // bins
        chunk = pairs[start:end]
        if not chunk:
            continue
        energies = [e for e, _ in chunk]
        distances = [d for _, d in chunk]
        top10_hits = sum(d <= 10 for d in distances)
        summaries.append({
            "bin": b + 1,
            "count": len(chunk),
            "energy_min": energies[0],
            "energy_max": energies[-1],
            "median_h": energies[len(energies) // 2],
            "median_dist": sorted(distances)[len(distances) // 2],
            "top10_rate": top10_hits / len(chunk),
        })
    return summaries

def run_dynamical_test(limit=2_000_000, windows=(5, 10, 20, 30, 50), show_bins=True):
    small_primes = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    quads, prime = find_prime_quadruplets(limit)
    print("Quaternionisch-dynamischer Test")
    print("W | Exakt% | Top-10% | MedianDist | MedianAbsErr | MedianH")
    print("-" * 76)
    all_results = {}
    for W in windows:
        results = []
        for i in range(W + 2, len(quads) - 1):
            r = predict_next_quad_dynamical(
                current_index=i,
                quads=quads,
                small_primes=small_primes,
                window=W
            )
            if r is not None:
                results.append(r)

        exact_hits = sum(
            r["prediction"] == r["true_next_quad_start"]
            for r in results
        )
        top10_hits = sum(
            r["candidate_distance"] is not None
            and r["candidate_distance"] <= 10
            for r in results
        )
        distances = sorted(
            r["candidate_distance"]
            for r in results
            if r["candidate_distance"] is not None
        )
        abs_errors = sorted(r["absolute_error"] for r in results)
        holonomies = sorted(r["holonomy_energy"] for r in results)
        median_dist = distances[len(distances) // 2] if distances else None
        median_abs_error = abs_errors[len(abs_errors) // 2] if abs_errors else None
        median_h = holonomies[len(holonomies) // 2] if holonomies else None

        print(
            f"{W:2d} | "
            f"{(exact_hits / len(results)) if results else 0:7.2%} | "
            f"{(top10_hits / len(results)) if results else 0:7.2%} | "
            f"{str(median_dist):>10} | "
            f"{str(median_abs_error):>12} | "
            f"{median_h:8.0f}"
        )

        if results:
            pairs = [
                (r["holonomy_energy"], r["candidate_distance"])
                for r in results
                if r["holonomy_energy"] is not None and r["candidate_distance"] is not None
            ]
            if pairs:
                hs = [h for h, _ in pairs]
                ds = [d for _, d in pairs]
                p_corr = pearson_corr(hs, ds)
                s_corr = spearman_corr(hs, ds)
                print(
                    "    Korrelation H~Distanz: "
                    f"Pearson={p_corr:.4f}  Spearman={s_corr:.4f}"
                )

            threshold = holonomy_threshold_from_quads(quads, window=W, quantile=0.5)
            hybrid_results = []
            for i in range(W + 2, len(quads) - 1):
                hr = predict_next_quad_hybrid(
                    current_index=i,
                    quads=quads,
                    small_primes=small_primes,
                    holonomy_threshold=threshold,
                    window=W,
                )
                if hr is not None:
                    hybrid_results.append(hr)
            hybrid_top10 = sum(
                r["candidate_distance"] is not None and r["candidate_distance"] <= 10
                for r in hybrid_results
            )
            hybrid_top10_rate = (hybrid_top10 / len(hybrid_results)) if hybrid_results else 0
            print(
                "    Hybrid (H<Median => dynamisch sonst baseline): "
                f"Top-10={hybrid_top10_rate:.2%}  Schwelle={threshold:.0f}"
            )

            orient_stats = orientation_summary(results)
            if orient_stats:
                print("    Orientierung (Dynamik):")
                print("    O | N   | MedianDist | Top-10%")
                for row in orient_stats:
                    print(
                        f"    {row['orientation']:>1} | "
                        f"{row['count']:>3} | "
                        f"{row['median_dist']:>10} | "
                        f"{row['top10_rate'] * 100:>6.2f}%"
                    )

        if show_bins and results:
            bins = holonomy_bin_summary(results, bins=4)
            if bins:
                print("    Holonomie-Bins (Q1 niedrig -> Q4 hoch):")
                print("    Bin | N   | MedianH   | MedianDist | Top-10%")
                for b in bins:
                    print(
                        f"    Q{b['bin']}  | "
                        f"{b['count']:>3} | "
                        f"{b['median_h']:>9.0f} | "
                        f"{b['median_dist']:>10} | "
                        f"{b['top10_rate'] * 100:>6.2f}%"
                    )
        all_results[W] = results
    return quads, all_results

if __name__ == "__main__":
    quads, dyn_results = run_dynamical_test(
        limit=2_000_000,
        windows=(5, 10, 20, 30, 50, 80),
        show_bins=True
    )