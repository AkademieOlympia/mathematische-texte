import csv
import math
from bisect import bisect_left, bisect_right
from collections import defaultdict
from sympy import primerange

EABC_CLASSES = {
    1: "E",
    5: "A",
    7: "B",
    11: "C",
}

ORDER = ["E", "A", "B", "C"]
IDX = {c: i for i, c in enumerate(ORDER)}


def _finite(xs):
    return [x for x in xs if isinstance(x, (int, float)) and not math.isnan(float(x))]


def _mean(xs):
    xs = _finite(xs)
    return sum(xs) / len(xs) if xs else float("nan")


def _stdev_sample(xs):
    xs = _finite(xs)
    n = len(xs)
    if n < 2:
        return float("nan")
    m = sum(xs) / n
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1))


def _var_sample(xs):
    xs = _finite(xs)
    n = len(xs)
    if n < 2:
        return float("nan")
    m = sum(xs) / n
    return sum((x - m) ** 2 for x in xs) / (n - 1)


def summarize_by_radius(rows):
    """Gruppierung nach R (wie pandas groupby 'R')."""
    by_r = defaultdict(list)
    for row in rows:
        by_r[row["R"]].append(row)

    out = []
    for R in sorted(by_r.keys()):
        grp = by_r[R]
        out.append(
            {
                "R": R,
                "mean_H_iso": _mean([r["H_iso"] for r in grp]),
                "var_H_iso": _var_sample([r["H_iso"] for r in grp]),
                "mean_alpha_direct": _mean([r["alpha_direct"] for r in grp]),
                "std_alpha_direct": _stdev_sample([r["alpha_direct"] for r in grp]),
                "mean_alpha_matrix": _mean([r["alpha_matrix"] for r in grp]),
                "std_alpha_matrix": _stdev_sample([r["alpha_matrix"] for r in grp]),
                "mean_T_matrix": _mean([r["T_matrix"] for r in grp]),
                "std_T_matrix": _stdev_sample([r["T_matrix"] for r in grp]),
                "mean_matrix_factor": _mean([r["matrix_factor"] for r in grp]),
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
    print(sep.join(fn.ljust(widths[i]) for i, fn in enumerate(fieldnames)))
    print(sep.join("-" * w for w in widths))
    for sr in str_rows:
        print(sep.join(sr[i].ljust(widths[i]) for i in range(len(fieldnames))))


def eabc_prime_slots(m, R):
    a = max(3, m - R)
    b = m + R

    slots = []

    for p in primerange(a, b + 1):
        r = int(p % 12)
        if r in EABC_CLASSES:
            slots.append((int(p), EABC_CLASSES[r]))

    return slots, a, b


def transition_matrix(slots):
    P = [[0.0] * 4 for _ in range(4)]

    if len(slots) < 2:
        return P

    for (_, c1), (_, c2) in zip(slots[:-1], slots[1:]):
        i = IDX[c1]
        j = IDX[c2]
        P[i][j] += 1.0

    total = sum(sum(row) for row in P)

    if total > 0:
        inv = 1.0 / total
        P = [[x * inv for x in row] for row in P]

    return P


def _build_slot_prefix(m, max_R):
    """Einmaliger Aufbau für alle R eines festen Zentrums m."""
    a0 = max(3, m - max_R)
    b0 = m + max_R

    primes = []
    cls = []
    rho_pref = [0.0]

    for p in primerange(a0, b0 + 1):
        r = int(p % 12)
        if r in EABC_CLASSES:
            p = int(p)
            primes.append(p)
            cls.append(IDX[EABC_CLASSES[r]])
            rho_pref.append(rho_pref[-1] + 1.0 / p)

    n = len(primes)
    edge_pref = [[0] * 16]
    for k in range(1, n):
        row = edge_pref[-1].copy()
        row[cls[k - 1] * 4 + cls[k]] += 1
        edge_pref.append(row)

    return {
        "m": m,
        "primes": primes,
        "rho_pref": rho_pref,
        "edge_pref": edge_pref,
    }


def _window_stats_from_prefix(prefix_data, R):
    m = prefix_data["m"]
    a = max(3, m - R)
    b = m + R

    primes = prefix_data["primes"]
    rho_pref = prefix_data["rho_pref"]
    edge_pref = prefix_data["edge_pref"]

    lo = bisect_left(primes, a)
    hi = bisect_right(primes, b)
    n_slots = hi - lo
    rho = rho_pref[hi] - rho_pref[lo]

    P = [[0.0] * 4 for _ in range(4)]
    if n_slots >= 2:
        vec_hi = edge_pref[hi - 1]
        vec_lo = edge_pref[lo]
        counts = [vec_hi[i] - vec_lo[i] for i in range(16)]
        total = float(n_slots - 1)
        for idx, cnt in enumerate(counts):
            i, j = divmod(idx, 4)
            P[i][j] = cnt / total

    return P, rho, n_slots, a, b


def h_iso(P):
    target = 1.0 / 16.0
    return float(sum((P[i][j] - target) ** 2 for i in range(4) for j in range(4)))


def frob_norm(P):
    return float(math.sqrt(sum(P[i][j] ** 2 for i in range(4) for j in range(4))))


def rho_rec_from_slots(slots):
    return sum(1.0 / p for p, _ in slots)


def log_scale(a, b):
    a = max(3, a)
    if b <= a:
        return float("nan")

    return math.log(math.log(b)) - math.log(math.log(a))


def alpha_matrix_ptolemaic(m, R, prefix_data=None):
    if prefix_data is not None and prefix_data.get("m") == m:
        P, rho, n_slots, a, b = _window_stats_from_prefix(prefix_data, R)
    else:
        slots, a, b = eabc_prime_slots(m, R)
        P = transition_matrix(slots)
        rho = rho_rec_from_slots(slots)
        n_slots = len(slots)
    H = h_iso(P)
    L = log_scale(a, b)

    B = math.pi * (math.sqrt(3) / 2) * L

    frob = frob_norm(P)

    # Isotrope Frobeniusnorm ist 1/4.
    matrix_factor = 4.0 * frob

    if B == 0 or math.isnan(B):
        alpha_A = float("nan")
        alpha_M = float("nan")
    else:
        alpha_A = rho / B
        alpha_M = alpha_A * matrix_factor

    return {
        "m": m,
        "R": R,
        "a": a,
        "b": b,
        "N_slots": n_slots,
        "rho_rec": rho,
        "L": L,
        "B_ptol": B,
        "H_iso": H,
        "frob_P": frob,
        "matrix_factor": matrix_factor,
        "alpha_direct": alpha_A,
        "alpha_matrix": alpha_M,
        "T_direct": alpha_A * math.pi if not math.isnan(alpha_A) else float("nan"),
        "T_matrix": alpha_M * math.pi if not math.isnan(alpha_M) else float("nan"),
    }


def run(centers, radii):
    rows = []
    max_R = max(radii)

    for m in centers:
        prefix_data = _build_slot_prefix(m, max_R)
        for R in radii:
            rows.append(alpha_matrix_ptolemaic(m, R, prefix_data=prefix_data))

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

    rows = run(centers, radii)

    raw_fields = [
        "m",
        "R",
        "a",
        "b",
        "N_slots",
        "rho_rec",
        "L",
        "B_ptol",
        "H_iso",
        "frob_P",
        "matrix_factor",
        "alpha_direct",
        "alpha_matrix",
        "T_direct",
        "T_matrix",
    ]

    print("\n=== Rohdaten ===")
    print_table(rows, raw_fields)

    summary = summarize_by_radius(rows)

    summary_fields = [
        "R",
        "mean_H_iso",
        "var_H_iso",
        "mean_alpha_direct",
        "std_alpha_direct",
        "mean_alpha_matrix",
        "std_alpha_matrix",
        "mean_T_matrix",
        "std_T_matrix",
        "mean_matrix_factor",
    ]

    print("\n=== Zusammenfassung nach Radius ===")
    print_table(summary, summary_fields)

    write_csv("bamberger_ptolemaic_matrix_alpha_raw.csv", rows, raw_fields)
    write_csv("bamberger_ptolemaic_matrix_alpha_summary.csv", summary, summary_fields)
