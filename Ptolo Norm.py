import csv
import math
from bisect import bisect_left, bisect_right
from collections import defaultdict
from sympy import primerange

EABC_RESIDUES = {1, 5, 7, 11}


def _finite(xs):
    return [x for x in xs if isinstance(x, (int, float)) and not math.isnan(float(x))]


def _mean(xs):
    xs = _finite(xs)
    return sum(xs) / len(xs) if xs else float("nan")


def _stdev_sample(xs):
    """Stichproben-Standardabweichung (ddof=1), wie pandas Series.std()."""
    xs = _finite(xs)
    n = len(xs)
    if n < 2:
        return float("nan")
    m = sum(xs) / n
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1))


def _var_sample(xs):
    """Stichprobenvarianz (ddof=1), wie pandas Series.var()."""
    xs = _finite(xs)
    n = len(xs)
    if n < 2:
        return float("nan")
    m = sum(xs) / n
    return sum((x - m) ** 2 for x in xs) / (n - 1)


def _min_f(xs):
    xs = _finite(xs)
    return min(xs) if xs else float("nan")


def _max_f(xs):
    xs = _finite(xs)
    return max(xs) if xs else float("nan")


def summarize_by_radius(rows):
    """Gruppierung nach R; gleiche Kennzahlen wie die frühere pandas-agg."""
    by_r = defaultdict(list)
    for row in rows:
        by_r[row["R"]].append(row)

    out = []
    for R in sorted(by_r.keys()):
        grp = by_r[R]
        alphas = [r["alpha_arith"] for r in grp]
        Ts = [r["T_EABC"] for r in grp]
        out.append(
            {
                "R": R,
                "mean_alpha": _mean(alphas),
                "std_alpha": _stdev_sample(alphas),
                "var_alpha": _var_sample(alphas),
                "min_alpha": _min_f(alphas),
                "max_alpha": _max_f(alphas),
                "mean_T": _mean(Ts),
                "std_T": _stdev_sample(Ts),
            }
        )
    return out


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def print_table(rows, fieldnames, float_fmt="{:.12g}"):
    if not rows:
        print("(leer)")
        return

    def cell(v):
        if isinstance(v, float):
            if math.isnan(v):
                return "nan"
            return float_fmt.format(v)
        return str(v)

    str_rows = [[cell(row.get(k, "")) for k in fieldnames] for row in rows]
    widths = [max(len(fieldnames[i]), *(len(r[i]) for r in str_rows)) for i in range(len(fieldnames))]
    sep = "  "

    header = sep.join(fn.ljust(widths[i]) for i, fn in enumerate(fieldnames))
    print(header)
    print(sep.join("-" * w for w in widths))
    for sr in str_rows:
        print(sep.join(sr[i].ljust(widths[i]) for i in range(len(fieldnames))))


def ptolemaic_balance(a, b):
    """
    Ptolemäische Balance:
    pi * sqrt(3)/2 * logarithmischer Skalenfaktor.
    """
    if a < 3:
        a = 3
    if b <= a:
        return float("nan")

    L = math.log(math.log(b)) - math.log(math.log(a))

    return math.pi * (math.sqrt(3) / 2) * L


def rho_rec(m, R):
    """
    Harmonische Summe der primen EABC-Reste in H_R(m).
    """
    a = max(3, m - R)
    b = m + R

    s = 0.0
    count = 0

    for p in primerange(a, b + 1):
        if p % 12 in EABC_RESIDUES:
            s += 1.0 / p
            count += 1

    return s, count, a, b


def _build_eabc_prefix(m, max_R):
    """Einmaliger Primzahl-Scan pro Zentrum m bis Radius max_R."""
    a0 = max(3, m - max_R)
    b0 = m + max_R
    primes = []
    pref = [0.0]
    for p in primerange(a0, b0 + 1):
        if p % 12 in EABC_RESIDUES:
            p = int(p)
            primes.append(p)
            pref.append(pref[-1] + 1.0 / p)
    return {"m": m, "primes": primes, "pref": pref}


def _rho_rec_from_prefix(prefix_data, R):
    m = prefix_data["m"]
    a = max(3, m - R)
    b = m + R
    primes = prefix_data["primes"]
    pref = prefix_data["pref"]
    lo = bisect_left(primes, a)
    hi = bisect_right(primes, b)
    return pref[hi] - pref[lo], hi - lo, a, b


def alpha_arith(m, R, prefix_data=None):
    """
    Ptolemäisch normierte arithmetische Kopplung.
    """
    if prefix_data is not None and prefix_data.get("m") == m:
        rho, count, a, b = _rho_rec_from_prefix(prefix_data, R)
    else:
        rho, count, a, b = rho_rec(m, R)
    B = ptolemaic_balance(a, b)

    if B == 0 or math.isnan(B):
        alpha = float("nan")
    else:
        alpha = rho / B

    return {
        "m": m,
        "R": R,
        "a": a,
        "b": b,
        "N_primes_EABC": count,
        "rho_rec": rho,
        "B_ptol": B,
        "alpha_arith": alpha,
        "T_EABC": alpha * math.pi if not math.isnan(alpha) else float("nan"),
    }


def run_test(centers, radii):
    rows = []
    max_R = max(radii)
    for m in centers:
        prefix_data = _build_eabc_prefix(m, max_R)
        for R in radii:
            rows.append(alpha_arith(m, R, prefix_data=prefix_data))
    return rows


if __name__ == "__main__":

    base = 10_000_000

    centers = [
        base + 59,
        base + 101,
        base + 137,
        base + 233,
        base + 359,
        base + 719,
        base + 1439,
        base + 2879,
        base + 7919,
        base + 113160,
    ]

    radii = [
        1000,
        2000,
        5000,
        10000,
        20000,
        50000,
        100000,
        200000,
        500000,
        1000000,
    ]

    rows = run_test(centers, radii)

    raw_fields = ["m", "R", "a", "b", "N_primes_EABC", "rho_rec", "B_ptol", "alpha_arith", "T_EABC"]
    print_table(rows, raw_fields)

    summary = summarize_by_radius(rows)

    print("\n=== Radius-Zusammenfassung ===")
    summary_fields = [
        "R",
        "mean_alpha",
        "std_alpha",
        "var_alpha",
        "min_alpha",
        "max_alpha",
        "mean_T",
        "std_T",
    ]
    print_table(summary, summary_fields)

    write_csv("bamberger_ptolemaic_alpha_raw.csv", rows, raw_fields)
    write_csv("bamberger_ptolemaic_alpha_summary.csv", summary, summary_fields)
