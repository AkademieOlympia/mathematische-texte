"""
W-Scan mod 420: Kalt/Warm 311/107 (Spiegelpaar, entspricht 11/47 mod 60 im LCM-Raster 420).

Berechnet pro W: Delta-Clustering vs. exakte Random-Baseline, Polarität P = Delta_hot - Delta_cold.

Mit --perms > 0: scan-korrigierter Permutationstest (Look-elsewhere), analog BM60:
  M_real^420  = max_W P^420_{311/107}(W)     (feste Baselines, echte Label an Positionen)
  M_perm^420  = max_W P^420_{311/107,perm}(W) (Permutation der mod-420-Multimenge der Labels)
  p_scan^420  = Pr(M_perm^420 >= M_real^420) unter den Ziehungen.

Heuristik (Literatur-nah): p_scan^420 < 0,05 = „Vererbung“ bzw. Signal nicht durch Zufalls-Umlabelung
für die gleiche Räumlichkeit gedeckt; andernfalls bleibt ein positiver Hinweis.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd

_DEFAULT_N = 50_000_000
_DEFAULT_PERMS = 0
_DEFAULT_SEED = 42
MOD = 420
W_LIST = np.array([400, 600, 800, 1000, 1200, 1500, 2000, 3000], dtype=np.int64)

# 420-Spiegelpaar (CRT-Lift: 311≡11 (60), 107≡47 (60) u. a.)
COLD = 311
HOT = 107
_DEFAULT_CSV = "bm420_p_311_107_wscan.csv"


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


def _build_baselines(starts: np.ndarray, labels: np.ndarray, n: int) -> dict[tuple[int, int], float]:
    b: dict[tuple[int, int], float] = {}
    for W in W_LIST:
        w = int(W)
        b[(COLD, w)] = exact_random_baseline(starts[labels == COLD], n, w)
        b[(HOT, w)] = exact_random_baseline(starts[labels == HOT], n, w)
    return b


def p_311_107_for_w(starts: np.ndarray, labels: np.ndarray, W: int, baselines: dict, n: int) -> float:
    w = int(W)
    d_c = cluster_prob(starts[labels == COLD], w) - baselines[(COLD, w)]
    d_h = cluster_prob(starts[labels == HOT], w) - baselines[(HOT, w)]
    return d_h - d_c


def _export_perm_results(
    out_dir: Path,
    n: int,
    perms: int,
    seed: int,
    elapsed_s: float,
    w_list: np.ndarray,
    real_curve: np.ndarray,
    m_real: float,
    w_peak: int,
    mean_perm: float,
    std_perm: float,
    z_scan: float,
    p_scan: float,
    perm_max_arr: np.ndarray,
) -> None:
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = f"n{n}_p{perms}"
    curve_csv = out_dir / f"bm420_p_311_107_wscan_curve_{tag}.csv"
    with curve_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["N", "W", "P_311_107"])
        for W, pol in zip(w_list, real_curve):
            w.writerow([n, int(W), float(pol)])

    summary = {
        "script": "bm420_p_311_107_wscan.py",
        "MOD": MOD,
        "channels": {"cold": COLD, "hot": HOT},
        "N": n,
        "W_LIST": [int(x) for x in w_list],
        "permutations": perms,
        "seed": seed,
        "elapsed_s": elapsed_s,
        "M_real_420": m_real,
        "W_peak_420": w_peak,
        "M_real": m_real,
        "W_peak": w_peak,
        "mean_M_perm_420": mean_perm,
        "mean_perm_max": mean_perm,
        "std_M_perm_420": std_perm,
        "std_perm_max": std_perm,
        "z_scan_420": z_scan,
        "p_scan_420": p_scan,
        "z_scan": z_scan,
        "p_scan": p_scan,
        "vererbung_signifikant_05": bool(p_scan < 0.05),
        "quantiles_M_perm_420": {
            str(q): float(np.nanquantile(perm_max_arr, q)) for q in (0.5, 0.9, 0.95, 0.975, 0.99)
        },
    }
    (out_dir / f"bm420_p_311_107_wscan_summary_{tag}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mod-420 W-Scan 311/107; optional scan-korrigierter Permutationstest (Look-elsewhere-Max über W_LIST)."
    )
    parser.add_argument("-n", "--n", type=int, default=_DEFAULT_N, help="Oberes Sieblimit (Default: 50M).")
    parser.add_argument(
        "--perms",
        type=int,
        default=_DEFAULT_PERMS,
        help="Permutationen (0 = nur Kurve, Default; z. B. 5000 für vollen Test).",
    )
    parser.add_argument("--seed", type=int, default=_DEFAULT_SEED, help="Zufallssaat (Default: 42).")
    parser.add_argument(
        "-o",
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / _DEFAULT_CSV,
        help="CSV bei --perms 0 (Default: bm420_p_311_107_wscan.csv neben Skript).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Bei --perms > 0: Ordner für Kurve + Summary-JSON (mit n/p-Tag im Dateinamen).",
    )
    args = parser.parse_args()
    n = int(args.n)
    perms = int(args.perms)
    if perms < 0:
        raise SystemExit("--perms muss >= 0 sein.")
    np.random.seed(int(args.seed))

    t_start = time.time()
    print("N =", n, "MOD =", MOD, "COLD/HOT =", COLD, HOT, "| perms =", perms, "| seed =", int(args.seed))
    print("Sieb und Zwillinge...")
    t0 = time.time()
    is_prime = fast_sieve(n)
    starts = get_twins(is_prime, n)
    labels = starts % MOD
    print("Zwillinge:", len(starts), "Zeit:", f"{time.time() - t0:.2f}s")

    print("Baselines (311/107) pro W ...")
    t1 = time.time()
    baselines = _build_baselines(starts, labels, n)
    print("fertig in", f"{time.time() - t1:.2f}s")

    real_curve = np.array([p_311_107_for_w(starts, labels, int(W), baselines, n) for W in W_LIST], dtype=float)
    m_real = float(np.nanmax(real_curve))
    w_peak = int(W_LIST[int(np.nanargmax(real_curve))])

    if perms == 0:
        rows = []
        for W in W_LIST:
            qs_cold = starts[labels == COLD]
            qs_hot = starts[labels == HOT]
            wv = int(W)
            c_cold = cluster_prob(qs_cold, wv)
            c_hot = cluster_prob(qs_hot, wv)
            rb_cold = baselines[(COLD, wv)]
            rb_hot = baselines[(HOT, wv)]
            d_cold = c_cold - rb_cold
            d_hot = c_hot - rb_hot
            pol = d_hot - d_cold
            rows.append(
                {
                    "N": n,
                    "W": wv,
                    "cold": COLD,
                    "hot": HOT,
                    "count_cold": len(qs_cold),
                    "count_hot": len(qs_hot),
                    "Delta_cold": d_cold,
                    "Delta_hot": d_hot,
                    "P_311_107": pol,
                }
            )
        df = pd.DataFrame(rows)
        print()
        print(df.to_string(index=False))
        out = args.out.resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, index=False)
        print("geschrieben:", out)
        return

    print("Permutationen (Look-elsewhere-Max über W_LIST) ...")
    t2 = time.time()
    perm_max: list[float] = []
    for k in range(perms):
        if (k + 1) % 500 == 0:
            print("Permutation", k + 1, "/", perms)
        pl = np.random.permutation(labels)
        curve = np.array(
            [p_311_107_for_w(starts, pl, int(W), baselines, n) for W in W_LIST], dtype=float
        )
        perm_max.append(float(np.nanmax(curve)))

    print("fertig in", f"{time.time() - t2:.2f}s")
    perm_max_arr = np.array(perm_max, dtype=float)
    mean_perm = float(np.mean(perm_max_arr))
    std_perm = float(np.std(perm_max_arr))
    z_scan = (m_real - mean_perm) / std_perm if std_perm > 0 else float("nan")
    p_scan = float(np.mean(perm_max_arr >= m_real))
    elapsed_s = time.time() - t_start

    print()
    print("── BM420: M_real^420 vs. M_perm^420  (max_W P^420_{311/107}, Look-elsewhere) ──")
    print("W_LIST =", [int(w) for w in W_LIST])
    print("real_curve =", real_curve)
    print("M_real^420  =", m_real, "  (W_peak =", w_peak, ")")
    print("M_perm^420  = mean(max_W P perm) =", mean_perm, "| std(perm-Max) =", std_perm)
    print("z_scan^420  =", z_scan, "| p_scan^420 =", p_scan)
    if p_scan < 0.05:
        print("Hinweis: p_scan^420 < 0,05  →  Vererbung/Signal nach diesem Nnull formal signifikant.")
    else:
        print("Hinweis: p_scan^420 >= 0,05  →  bleibt ein positiver Hinweis, nicht formal signifikant.")
    print("Laufzeit (gesamt) =", f"{elapsed_s:.2f}s")
    print()
    print("Quantile M_perm^420 (over-W-Max, empir.):")
    for q in [0.5, 0.9, 0.95, 0.975, 0.99]:
        print(q, float(np.nanquantile(perm_max_arr, q)))

    out_dir = args.out_dir
    _export_perm_results(
        out_dir,
        n,
        perms,
        int(args.seed),
        elapsed_s,
        W_LIST,
        real_curve,
        m_real,
        w_peak,
        mean_perm,
        std_perm,
        z_scan,
        p_scan,
        perm_max_arr,
    )
    tag = f"n{n}_p{perms}"
    print()
    print("Export:", out_dir.resolve() / f"bm420_p_311_107_wscan_curve_{tag}.csv")
    print("Export:", out_dir.resolve() / f"bm420_p_311_107_wscan_summary_{tag}.json")


if __name__ == "__main__":
    main()
