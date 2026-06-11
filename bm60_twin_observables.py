from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# BM60 Twin Observables
# ============================================================
#
# Observablen:
#
#   P_11_47(W) = Delta_47(W) - Delta_11(W)
#       Primäre Spiegelpolarität:
#       11 -> 47 unter r -> -r-2 mod 60.
#       Deutung: quadratische Störung 47/49 = 7^2.
#
#   P_17_41(W) = Delta_41(W) - Delta_17(W)
#       Vergleichs-Spiegelpaar.
#
#   S_seam(W) = Delta_59(W) - Delta_29(W)
#       Randnahtspannung:
#       29 und 59 sind selbstkonjugiert.
#       59 liegt an der Wheel-Naht (-1,+1), 29 im Inneren.
#
#   Edge_59(W) = Delta_59(W) - mean(Delta_rest)
#       Absolute Randdämpfung von 59 gegen alle anderen Kanäle.
#
# ============================================================


MOD = 60
CHANNELS = np.array([11, 17, 29, 41, 47, 59], dtype=np.int64)

DEFAULT_W_LIST = np.array(
    [300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000, 3000],
    dtype=np.int64,
)


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
    """
    Exakte Fensterwahrscheinlichkeit:
    Wahrscheinlichkeit, dass ein zufälliges Fenster [x, x+w]
    mindestens einen Treffer aus qs enthält.

    Berechnet als Länge der Vereinigung der Intervalle [q-w, q].
    """
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
            if bi > cur_b:
                cur_b = bi
        else:
            total += cur_b - cur_a
            cur_a, cur_b = ai, bi

    total += cur_b - cur_a

    return total / max_start


def cluster_prob(qs: np.ndarray, w: int) -> float:
    if len(qs) < 3:
        return float("nan")

    gaps = np.diff(qs)
    return float(np.mean(gaps < w))


def compute_baselines(
    starts: np.ndarray,
    labels: np.ndarray,
    n: int,
    w_list: np.ndarray,
) -> dict[tuple[int, int], float]:
    baselines: dict[tuple[int, int], float] = {}

    for r in CHANNELS:
        qs = starts[labels == r]
        for W in w_list:
            baselines[(int(r), int(W))] = exact_random_baseline(qs, n, int(W))

    return baselines


def compute_deltas(
    starts: np.ndarray,
    labels: np.ndarray,
    n: int,
    w_list: np.ndarray,
    baselines: dict[tuple[int, int], float],
) -> pd.DataFrame:
    rows = []

    for W in w_list:
        for r in CHANNELS:
            qs = starts[labels == r]

            c = cluster_prob(qs, int(W))
            rb = baselines[(int(r), int(W))]
            d = c - rb

            rows.append(
                {
                    "N": int(n),
                    "W": int(W),
                    "channel": int(r),
                    "count": int(len(qs)),
                    "cluster_prob": c,
                    "random_prob": rb,
                    "delta": d,
                }
            )

    return pd.DataFrame(rows)


def observables_from_delta_df(delta_df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for W, sub in delta_df.groupby("W"):
        d = {int(row.channel): float(row.delta) for row in sub.itertuples()}

        mean_rest_59 = np.mean([d[r] for r in CHANNELS if int(r) != 59])
        mean_all = np.mean([d[int(r)] for r in CHANNELS])

        rows.append(
            {
                "W": int(W),
                "Delta_11": d[11],
                "Delta_17": d[17],
                "Delta_29": d[29],
                "Delta_41": d[41],
                "Delta_47": d[47],
                "Delta_59": d[59],

                # Primäre Observablen
                "P_11_47": d[47] - d[11],
                "P_17_41": d[41] - d[17],
                "S_seam_29_59": d[59] - d[29],

                # Zusätzliche Randdiagnostik
                "Edge_59_vs_rest": d[59] - mean_rest_59,
                "Delta_59_minus_mean_all": d[59] - mean_all,

                # Dominanzen
                "Dom_11_47_vs_17_41": (d[47] - d[11]) - (d[41] - d[17]),
                "Dom_11_47_vs_seam_abs": (d[47] - d[11]) - abs(d[59] - d[29]),
            }
        )

    return pd.DataFrame(rows).sort_values("W")


def max_stat(curve_df: pd.DataFrame, column: str) -> dict:
    idx = curve_df[column].idxmax()
    row = curve_df.loc[idx]

    return {
        "observable": column,
        "M_real": float(row[column]),
        "W_peak": int(row["W"]),
    }


def min_stat(curve_df: pd.DataFrame, column: str) -> dict:
    idx = curve_df[column].idxmin()
    row = curve_df.loc[idx]

    return {
        "observable": column,
        "M_min_real": float(row[column]),
        "W_min": int(row["W"]),
    }


def compute_observable_curve_for_labels(
    starts: np.ndarray,
    labels: np.ndarray,
    n: int,
    w_list: np.ndarray,
    baselines: dict[tuple[int, int], float],
) -> pd.DataFrame:
    delta_df = compute_deltas(starts, labels, n, w_list, baselines)
    return observables_from_delta_df(delta_df)


def permutation_scan_tests(
    starts: np.ndarray,
    labels: np.ndarray,
    n: int,
    w_list: np.ndarray,
    baselines: dict[tuple[int, int], float],
    real_curve: pd.DataFrame,
    perms: int,
    seed: int,
) -> dict:
    rng = np.random.default_rng(seed)

    # Realstatistiken
    real_stats = {
        "P_11_47": max_stat(real_curve, "P_11_47"),
        "S_seam_positive": max_stat(real_curve, "S_seam_29_59"),
        "S_seam_negative": min_stat(real_curve, "S_seam_29_59"),
        "Edge_59_vs_rest_positive": max_stat(real_curve, "Edge_59_vs_rest"),
        "Edge_59_vs_rest_negative": min_stat(real_curve, "Edge_59_vs_rest"),
    }

    if perms <= 0:
        return {
            "permutations": 0,
            "real_stats": real_stats,
            "perm_tests": {},
        }

    perm_max_p1147 = []
    perm_max_seam = []
    perm_min_seam = []
    perm_absmax_seam = []
    perm_max_edge59 = []
    perm_min_edge59 = []
    perm_absmax_edge59 = []

    t0 = time.time()

    for k in range(perms):
        if (k + 1) % 500 == 0:
            print(f"    Permutation {k+1}/{perms}")

        perm_labels = rng.permutation(labels)
        curve = compute_observable_curve_for_labels(
            starts, perm_labels, n, w_list, baselines
        )

        p1147_vals = curve["P_11_47"].to_numpy()
        seam_vals = curve["S_seam_29_59"].to_numpy()
        edge_vals = curve["Edge_59_vs_rest"].to_numpy()

        perm_max_p1147.append(float(np.max(p1147_vals)))

        perm_max_seam.append(float(np.max(seam_vals)))
        perm_min_seam.append(float(np.min(seam_vals)))
        perm_absmax_seam.append(float(np.max(np.abs(seam_vals))))

        perm_max_edge59.append(float(np.max(edge_vals)))
        perm_min_edge59.append(float(np.min(edge_vals)))
        perm_absmax_edge59.append(float(np.max(np.abs(edge_vals))))

    elapsed = time.time() - t0
    print(f"    Permutationen fertig in {elapsed:.2f}s")

    perm_arrays = {
        "P_11_47_max": np.array(perm_max_p1147),
        "S_seam_max": np.array(perm_max_seam),
        "S_seam_min": np.array(perm_min_seam),
        "S_seam_absmax": np.array(perm_absmax_seam),
        "Edge_59_max": np.array(perm_max_edge59),
        "Edge_59_min": np.array(perm_min_edge59),
        "Edge_59_absmax": np.array(perm_absmax_edge59),
    }

    def one_sided_max_test(real_value: float, arr: np.ndarray) -> dict:
        mean = float(np.mean(arr))
        std = float(np.std(arr))
        z = (real_value - mean) / std if std > 0 else float("nan")
        p = float(np.mean(arr >= real_value))
        return {
            "real": real_value,
            "mean_perm": mean,
            "std_perm": std,
            "z": z,
            "p": p,
            "q50": float(np.quantile(arr, 0.50)),
            "q90": float(np.quantile(arr, 0.90)),
            "q95": float(np.quantile(arr, 0.95)),
            "q975": float(np.quantile(arr, 0.975)),
            "q99": float(np.quantile(arr, 0.99)),
        }

    def one_sided_min_test(real_value: float, arr_min: np.ndarray) -> dict:
        # Für negative Minima: extrem ist arr <= real_value
        mean = float(np.mean(arr_min))
        std = float(np.std(arr_min))
        z = (real_value - mean) / std if std > 0 else float("nan")
        p = float(np.mean(arr_min <= real_value))
        return {
            "real": real_value,
            "mean_perm": mean,
            "std_perm": std,
            "z": z,
            "p": p,
            "q01": float(np.quantile(arr_min, 0.01)),
            "q025": float(np.quantile(arr_min, 0.025)),
            "q05": float(np.quantile(arr_min, 0.05)),
            "q10": float(np.quantile(arr_min, 0.10)),
            "q50": float(np.quantile(arr_min, 0.50)),
        }

    def absmax_test(real_value_abs: float, arr_absmax: np.ndarray) -> dict:
        mean = float(np.mean(arr_absmax))
        std = float(np.std(arr_absmax))
        z = (real_value_abs - mean) / std if std > 0 else float("nan")
        p = float(np.mean(arr_absmax >= real_value_abs))
        return {
            "real_abs": real_value_abs,
            "mean_perm": mean,
            "std_perm": std,
            "z": z,
            "p": p,
            "q50": float(np.quantile(arr_absmax, 0.50)),
            "q90": float(np.quantile(arr_absmax, 0.90)),
            "q95": float(np.quantile(arr_absmax, 0.95)),
            "q975": float(np.quantile(arr_absmax, 0.975)),
            "q99": float(np.quantile(arr_absmax, 0.99)),
        }

    real_p1147 = real_stats["P_11_47"]["M_real"]

    real_seam_max = real_stats["S_seam_positive"]["M_real"]
    real_seam_min = real_stats["S_seam_negative"]["M_min_real"]
    real_seam_abs = max(abs(real_seam_max), abs(real_seam_min))

    real_edge_max = real_stats["Edge_59_vs_rest_positive"]["M_real"]
    real_edge_min = real_stats["Edge_59_vs_rest_negative"]["M_min_real"]
    real_edge_abs = max(abs(real_edge_max), abs(real_edge_min))

    perm_tests = {
        "P_11_47_max_scan": one_sided_max_test(real_p1147, perm_arrays["P_11_47_max"]),

        "S_seam_positive_scan": one_sided_max_test(real_seam_max, perm_arrays["S_seam_max"]),
        "S_seam_negative_scan": one_sided_min_test(real_seam_min, perm_arrays["S_seam_min"]),
        "S_seam_absmax_scan": absmax_test(real_seam_abs, perm_arrays["S_seam_absmax"]),

        "Edge_59_positive_scan": one_sided_max_test(real_edge_max, perm_arrays["Edge_59_max"]),
        "Edge_59_negative_scan": one_sided_min_test(real_edge_min, perm_arrays["Edge_59_min"]),
        "Edge_59_absmax_scan": absmax_test(real_edge_abs, perm_arrays["Edge_59_absmax"]),
    }

    return {
        "permutations": perms,
        "elapsed_perm_s": elapsed,
        "real_stats": real_stats,
        "perm_tests": perm_tests,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--n", type=int, default=200_000_000)
    parser.add_argument("--perms", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--out-prefix",
        type=str,
        default="bm60_twin_observables",
    )
    parser.add_argument(
        "--w-list",
        type=int,
        nargs="*",
        default=[int(x) for x in DEFAULT_W_LIST],
    )

    args = parser.parse_args()

    n = int(args.n)
    perms = int(args.perms)
    seed = int(args.seed)
    w_list = np.array(args.w_list, dtype=np.int64)

    print("============================================================")
    print("BM60 Twin Observables")
    print("============================================================")
    print("N      =", n)
    print("MOD    =", MOD)
    print("W_LIST =", [int(x) for x in w_list])
    print("PERMS  =", perms)
    print("seed   =", seed)
    print("============================================================")

    t0 = time.time()

    print("Sieb und Zwillinge...")
    is_prime = fast_sieve(n)
    starts = get_twins(is_prime, n)
    labels = starts % MOD
    print("Zwillinge:", len(starts), "Zeit:", f"{time.time() - t0:.2f}s")

    print("\nKanalhäufigkeiten:")
    counts = {int(r): int(np.sum(labels == r)) for r in CHANNELS}
    print(counts)

    print("\nBaselines...")
    t1 = time.time()
    baselines = compute_baselines(starts, labels, n, w_list)
    print("fertig in", f"{time.time() - t1:.2f}s")

    print("\nReale Delta-Kurven...")
    delta_df = compute_deltas(starts, labels, n, w_list, baselines)
    curve_df = observables_from_delta_df(delta_df)

    print(curve_df.to_string(index=False))

    print("\nPermutationstests...")
    scan_result = permutation_scan_tests(
        starts=starts,
        labels=labels,
        n=n,
        w_list=w_list,
        baselines=baselines,
        real_curve=curve_df,
        perms=perms,
        seed=seed,
    )

    out_prefix = Path(args.out_prefix)
    curve_path = out_prefix.with_name(f"{out_prefix.name}_curve_n{n}_p{perms}.csv")
    delta_path = out_prefix.with_name(f"{out_prefix.name}_deltas_n{n}_p{perms}.csv")
    summary_path = out_prefix.with_name(f"{out_prefix.name}_summary_n{n}_p{perms}.json")

    curve_df.to_csv(curve_path, index=False)
    delta_df.to_csv(delta_path, index=False)

    summary = {
        "script": "bm60_twin_observables.py",
        "N": n,
        "MOD": MOD,
        "W_LIST": [int(x) for x in w_list],
        "seed": seed,
        "counts": counts,
        "elapsed_total_s": time.time() - t0,
        **scan_result,
    }

    def _json_default(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)

    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=_json_default)

    print("\nExports:")
    print(" ", curve_path.resolve())
    print(" ", delta_path.resolve())
    print(" ", summary_path.resolve())

    print("\nKurzsummary:")
    for key, value in scan_result["real_stats"].items():
        print(key, ":", value)

    print("\nPermutationsergebnisse:")
    for key, value in scan_result["perm_tests"].items():
        print(key, ":", value)


if __name__ == "__main__":
    main()