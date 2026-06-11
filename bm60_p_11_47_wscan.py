import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np

from bm60_singular_series_twin import (
    bateman_horn_theoretical_pmf_mod60,
    compute_twin_prime_constant_C2,
    delta_normalized_resid_emp,
    hl_correction_vs_uniform,
    TWIN_MOD60,
)

MOD = 60
W_LIST = np.array([400, 600, 800, 1000, 1200, 1500, 2000, 3000], dtype=np.int64)
# Standard: 200M (Replikation des W-Scans gegenüber N=100M)
_DEFAULT_N = 200_000_000
_DEFAULT_PERMS = 5000
_DEFAULT_SEED = 42


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


def pol_11_47_for_W(starts: np.ndarray, labels: np.ndarray, W: int, baselines: dict) -> float:
    d11, d47 = deltas_11_47_for_w(starts, labels, W, baselines)
    return d47 - d11


def deltas_11_47_for_w(starts: np.ndarray, labels: np.ndarray, W: int, baselines: dict) -> tuple[float, float]:
    w = int(W)
    d11 = cluster_prob(starts[labels == 11], w) - baselines[(11, w)]
    d47 = cluster_prob(starts[labels == 47], w) - baselines[(47, w)]
    return d11, d47


def _counts_twin60(labels: np.ndarray) -> dict[int, int]:
    u, c = np.unique(labels, return_counts=True)
    return {int(a): int(b) for a, b in zip(u, c)}


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="BM60 Twin: W-Scan 11–47, scan-korrigierter Spiegeltest (Perm-Max).")
    p.add_argument(
        "-n",
        "--n",
        type=int,
        default=_DEFAULT_N,
        help=f"Obere Schranke fürs Sieb (Default: {_DEFAULT_N} = 200M).",
    )
    p.add_argument(
        "--perms",
        type=int,
        default=_DEFAULT_PERMS,
        help=f"Anzahl Permutationen (0 = nur W-Kurve, kein Look-elsewhere-Test, schnell; Default: {_DEFAULT_PERMS}).",
    )
    p.add_argument("--seed", type=int, default=_DEFAULT_SEED, help=f"Zufallssaat (Default: {_DEFAULT_SEED}).")
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Ordner für CSV/JSON (Default: Skriptverzeichnis).",
    )
    p.add_argument(
        "--no-export",
        action="store_true",
        help="Keine Ergebnisdateien schreiben.",
    )
    return p.parse_args()


def _export_results(
    out_dir: Path,
    n: int,
    perms: int,
    seed: int,
    elapsed_s: float,
    w_list: np.ndarray,
    real_curve: np.ndarray,
    real_curve_resid: np.ndarray,
    m_real: float,
    w_peak: int,
    m_resid: float,
    w_peak_resid: int,
    mean_perm: float,
    std_perm: float,
    z_scan: float,
    p_scan: float,
    mean_perm_resid: float,
    std_perm_resid: float,
    z_scan_resid: float,
    p_scan_resid: float,
    perm_max_arr: np.ndarray,
    perm_max_resid_arr: np.ndarray,
    extra: dict,
) -> None:
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = f"n{n}_p{perms}"

    curve_csv = out_dir / f"bm60_p_11_47_wscan_curve_{tag}.csv"
    with curve_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["N", "W", "Pol_11_47", "Pol_11_47_resid_emp_hl"])
        for W, pol, pr in zip(w_list, real_curve, real_curve_resid):
            w.writerow([n, int(W), float(pol), float(pr)])

    if len(perm_max_arr) == 0:
        qpm = None
        qpmr = None
    else:
        qpm = {str(q): float(np.quantile(perm_max_arr, q)) for q in (0.5, 0.9, 0.95, 0.975, 0.99)}
        qpmr = {str(q): float(np.nanquantile(perm_max_resid_arr, q)) for q in (0.5, 0.9, 0.95, 0.975, 0.99)}

    summary = {
        "script": "bm60_p_11_47_wscan.py",
        "N": n,
        "MOD": MOD,
        "W_LIST": [int(x) for x in w_list],
        "permutations": perms,
        "seed": seed,
        "elapsed_s": elapsed_s,
        "M_real": m_real,
        "W_peak": w_peak,
        "mean_perm_max": mean_perm,
        "std_perm_max": std_perm,
        "z_scan": z_scan,
        "p_scan": p_scan,
        "quantiles_perm_max": qpm,
        "M_resid_emp": m_resid,
        "W_peak_resid_emp": w_peak_resid,
        "mean_perm_max_resid_emp": mean_perm_resid,
        "std_perm_max_resid_emp": std_perm_resid,
        "z_scan_resid_emp": z_scan_resid,
        "p_scan_resid_emp": p_scan_resid,
        "quantiles_perm_max_resid_emp": qpmr,
        **extra,
    }
    json_path = out_dir / f"bm60_p_11_47_wscan_summary_{tag}.json"
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")


def main() -> None:
    args = _parse_args()
    n = int(args.n)
    perms = int(args.perms)
    if perms < 0:
        raise SystemExit("--perms muss >= 0 sein (0 = nur Kurve).")
    np.random.seed(int(args.seed))

    t_start = time.time()
    print("N =", n, "| MOD =", MOD, "| PERMS =", perms, "| seed =", int(args.seed))
    t0 = time.time()
    print("Sieb und Zwillinge...")
    is_prime = fast_sieve(n)
    starts = get_twins(is_prime, n)
    labels = starts % MOD
    print("Zwillinge:", len(starts), "Zeit:", f"{time.time() - t0:.2f}s")

    c2_info = compute_twin_prime_constant_C2()
    pmf_hl = bateman_horn_theoretical_pmf_mod60()
    count_by = _counts_twin60(labels)
    total_t = int(len(starts))
    kappa = hl_correction_vs_uniform(count_by, total_t)

    print()
    print("── Hardy–Littlewood (Bateman–Horn, Zwillingsmuster (0,2)) ──")
    print("C2 (partielles Produkt) ≈", c2_info.c2, f"(Primzahlen p≤{c2_info.prime_upper}, {c2_info.n_primes} Faktoren)")
    print("2*C2 ≈", c2_info.two_c2)
    print("Theoretische PMF fürs Muster p mod 60 ∈ {11,..,59} (Hauptterm):", pmf_hl[11], "pro Klasse (uniform 1/6)")
    print("Hinweis: Keine Differenz ändert d_47−d_11 multiplikativ nontrivial (nur gemeinsame Skalierung) —")
    print("  daher: P^theo_11/47 = (d_47−d_11)  für Linearkombination; empirische κ folgen unten.")
    n11 = int(count_by.get(11, 0))
    n47 = int(count_by.get(47, 0))
    print("Beobachtet n_11, n_47, T =", n11, n47, total_t, "| f_a = n_a/T:", f"{n11/total_t:.6f}", f"{n47/total_t:.6f}")
    print("κ_emp = f_a / (1/6) —", "11:", kappa[11], "47:", kappa[47], "(Abweichung von HL-Uniform)")

    print("Baselines pro W und Kanal 11/47 ...")
    t1 = time.time()
    baselines: dict[tuple[int, int], float] = {}
    for W in W_LIST:
        w = int(W)
        baselines[(11, w)] = exact_random_baseline(starts[labels == 11], n, w)
        baselines[(47, w)] = exact_random_baseline(starts[labels == 47], n, w)
    print("fertig in", f"{time.time() - t1:.2f}s")

    k11 = kappa[11]
    k47 = kappa[47]

    real_list: list[float] = []
    resid_list: list[float] = []
    for W in W_LIST:
        d11, d47 = deltas_11_47_for_w(starts, labels, int(W), baselines)
        real_list.append(d47 - d11)
        resid_list.append(delta_normalized_resid_emp(d11, d47, k11, k47))

    real_curve = np.array(real_list, dtype=float)
    real_curve_resid = np.array(resid_list, dtype=float)
    m_real = float(np.nanmax(real_curve))
    w_peak = int(W_LIST[int(np.nanargmax(real_curve))])
    m_resid = float(np.nanmax(real_curve_resid))
    w_peak_resid = int(W_LIST[int(np.nanargmax(real_curve_resid))])

    if perms == 0:
        print("Keine Permutationen (--perms 0), nur reale W-Kurve.")
        t2 = time.time()
        perm_max_arr = np.array([], dtype=float)
        perm_max_resid_arr = np.array([], dtype=float)
        p_scan = float("nan")
        z_scan = float("nan")
        mean_perm = float("nan")
        std_perm = float("nan")
        p_scan_resid = float("nan")
        z_scan_resid = float("nan")
        mean_perm_resid = float("nan")
        std_perm_resid = float("nan")
    else:
        print("Permutationen ...")
        t2 = time.time()
        perm_max: list[float] = []
        perm_max_resid: list[float] = []
        for k in range(perms):
            if (k + 1) % 500 == 0:
                print("Permutation", k + 1, "/", perms)

            perm_labels = np.random.permutation(labels)
            curve: list[float] = []
            curve_r: list[float] = []
            for W in W_LIST:
                d11, d47 = deltas_11_47_for_w(starts, perm_labels, int(W), baselines)
                curve.append(d47 - d11)
                curve_r.append(delta_normalized_resid_emp(d11, d47, k11, k47))
            a = np.array(curve, dtype=float)
            b = np.array(curve_r, dtype=float)
            perm_max.append(float(np.nanmax(a)))
            perm_max_resid.append(float(np.nanmax(b)))

        print("fertig in", f"{time.time() - t2:.2f}s")

        perm_max_arr = np.array(perm_max, dtype=float)
        perm_max_resid_arr = np.array(perm_max_resid, dtype=float)

        p_scan = float(np.mean(perm_max_arr >= m_real))
        std_perm = float(np.std(perm_max_arr))
        z_scan = (m_real - float(np.mean(perm_max_arr))) / std_perm if std_perm > 0 else float("nan")
        mean_perm = float(np.mean(perm_max_arr))

        p_scan_resid = float(np.mean(perm_max_resid_arr >= m_resid))
        std_perm_resid = float(np.std(perm_max_resid_arr))
        mean_perm_resid = float(np.mean(perm_max_resid_arr))
        z_scan_resid = (m_resid - mean_perm_resid) / std_perm_resid if std_perm_resid > 0 else float("nan")

    elapsed_s = time.time() - t_start

    print()
    print("Scan-korrigierter Spiegeltest")
    print("N =", n)
    print("W_LIST =", [int(w) for w in W_LIST])
    print("real_curve =", real_curve)
    print("M_real =", m_real)
    print("W_peak =", w_peak)
    print()
    print("P_11/47^resid(emp) = d_47/κ_47 - d_11/κ_11  (κ aus n_a/T vs 1/6)")
    print("real_curve_resid_emp =", real_curve_resid)
    print("M_resid_emp =", m_resid)
    print("W_peak_resid_emp =", w_peak_resid)
    print("mean_perm_max =", mean_perm)
    print("std_perm_max =", std_perm)
    print("z_scan =", z_scan)
    print("p_scan =", p_scan)
    print("mean_perm_max_resid_emp =", mean_perm_resid)
    print("std_perm_max_resid_emp =", std_perm_resid)
    print("z_scan_resid_emp =", z_scan_resid)
    print("p_scan_resid_emp =", p_scan_resid)
    print("Laufzeit (gesamt) =", f"{elapsed_s:.2f}s")

    if perms > 0:
        print()
        print("Quantile perm_max:")
        for q in [0.5, 0.9, 0.95, 0.975, 0.99]:
            print(q, float(np.quantile(perm_max_arr, q)))

        print()
        print("Quantile perm_max (resid emp):")
        for q in [0.5, 0.9, 0.95, 0.975, 0.99]:
            print(q, float(np.nanquantile(perm_max_resid_arr, q)))

    extra_export = {
        "twin_mod60": list(TWIN_MOD60),
        "C2_partial": c2_info.c2,
        "two_C2_partial": c2_info.two_c2,
        "C2_prime_upper": c2_info.prime_upper,
        "bateman_horn_pmf_mod60": {str(k): v for k, v in pmf_hl.items()},
        "counts_by_label_mod60": {str(k): v for k, v in count_by.items()},
        "kappa_emp_11_47": {"11": kappa[11], "47": kappa[47]},
    }

    if not args.no_export:
        _export_results(
            args.out_dir,
            n,
            perms,
            int(args.seed),
            elapsed_s,
            W_LIST,
            real_curve,
            real_curve_resid,
            m_real,
            w_peak,
            m_resid,
            w_peak_resid,
            mean_perm,
            std_perm,
            z_scan,
            p_scan,
            mean_perm_resid,
            std_perm_resid,
            z_scan_resid,
            p_scan_resid,
            perm_max_arr,
            perm_max_resid_arr,
            extra_export,
        )
        tag = f"n{n}_p{perms}"
        print()
        print("Export:", args.out_dir / f"bm60_p_11_47_wscan_curve_{tag}.csv")
        print("Export:", args.out_dir / f"bm60_p_11_47_wscan_summary_{tag}.json")


if __name__ == "__main__":
    main()
