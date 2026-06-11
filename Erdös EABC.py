import math
import random

# -------------------------------------------------
# Basisfunktionen
# -------------------------------------------------

def frac(x):
    return x - math.floor(x)

def family_mod12(p):
    r = p % 12
    if r == 1:
        return "E"
    elif r == 5:
        return "A"
    elif r == 7:
        return "B"
    elif r == 11:
        return "C"
    else:
        return None

def bamberg_phase(p, alpha, phi):
    fam = family_mod12(p)
    if fam is None:
        return None
    return frac(alpha * p + phi[fam])


def centered_mod1(x):
    y = x % 1.0
    if y >= 0.5:
        y -= 1.0
    return y


def residual_signature(quad, alpha, phi):
    p, q, r, s = quad

    thE = frac(alpha * p + phi["E"])
    thA = frac(alpha * q + phi["A"])
    thB = frac(alpha * r + phi["B"])
    thC = frac(alpha * s + phi["C"])

    predA = (2 * alpha + phi["A"] - phi["E"]) % 1.0
    predB = (6 * alpha + phi["B"] - phi["E"]) % 1.0
    predC = (8 * alpha + phi["C"] - phi["E"]) % 1.0

    obsA = (thA - thE) % 1.0
    obsB = (thB - thE) % 1.0
    obsC = (thC - thE) % 1.0

    return {
        "A_res": centered_mod1(obsA - predA),
        "B_res": centered_mod1(obsB - predB),
        "C_res": centered_mod1(obsC - predC),
    }


def summarize_residuals(quads, alpha, phi):
    valsA, valsB, valsC = [], [], []
    for quad in quads:
        res = residual_signature(quad, alpha, phi)
        valsA.append(res["A_res"])
        valsB.append(res["B_res"])
        valsC.append(res["C_res"])

    def mean_abs(xs):
        return sum(abs(x) for x in xs) / len(xs) if xs else None

    def variance(xs):
        if not xs:
            return None
        m = sum(xs) / len(xs)
        return sum((x - m) ** 2 for x in xs) / len(xs)

    return {
        "mean_abs_A": mean_abs(valsA),
        "mean_abs_B": mean_abs(valsB),
        "mean_abs_C": mean_abs(valsC),
        "var_A": variance(valsA),
        "var_B": variance(valsB),
        "var_C": variance(valsC),
    }


# -------------------------------------------------
# Kreisgeometrie
# -------------------------------------------------

def circular_gaps(phases):
    pts = sorted(phases)
    ext = pts + [pts[0] + 1.0]
    return [ext[i + 1] - ext[i] for i in range(len(pts))]

def circular_width(phases):
    # minimale Bogenlänge, die alle Punkte enthält
    gaps = circular_gaps(phases)
    return 1.0 - max(gaps)

def normalized_gap_signature(phases):
    # sortierte Gap-Signatur als Formvektor
    gaps = circular_gaps(phases)
    return sorted(gaps)

def l2_distance(v, w):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(v, w)))

# -------------------------------------------------
# Primzahlfunktionen
# -------------------------------------------------

def primes_up_to(N):
    sieve = [True] * (N + 1)
    sieve[0:2] = [False, False]
    for p in range(2, int(N**0.5) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:N + 1:p] = [False] * (((N - start) // p) + 1)
    return [p for p in range(N + 1) if sieve[p]]

def prime_quadruplets(N):
    plist = set(primes_up_to(N + 10))
    out = []
    for p in sorted(plist):
        if p + 8 <= N and all(q in plist for q in (p, p + 2, p + 6, p + 8)):
            out.append((p, p + 2, p + 6, p + 8))
    return out


# -------------------------------------------------
# EABC-Primlisten
# -------------------------------------------------

def build_family_prime_lists(N):
    fam_primes = {"E": [], "A": [], "B": [], "C": []}
    for p in primes_up_to(N):
        fam = family_mod12(p)
        if fam in fam_primes:
            fam_primes[fam].append(p)
    return fam_primes


def random_family_tuples(N, count=1000, seed=0):
    rng = random.Random(seed)
    fam_primes = build_family_prime_lists(N)
    out = []
    for _ in range(count):
        out.append((
            rng.choice(fam_primes["E"]),
            rng.choice(fam_primes["A"]),
            rng.choice(fam_primes["B"]),
            rng.choice(fam_primes["C"]),
        ))
    return out


# -------------------------------------------------
# Phasen für allgemeine EABC-Tupel
# -------------------------------------------------

def tuple_phases_eabc(tup, alpha, phi):
    pE, pA, pB, pC = tup
    return [
        frac(alpha * pE + phi["E"]),
        frac(alpha * pA + phi["A"]),
        frac(alpha * pB + phi["B"]),
        frac(alpha * pC + phi["C"]),
    ]


def summarize_general_tuples(tuples_list, alpha, phi):
    if not tuples_list:
        return None

    records = []
    for tup in tuples_list:
        phases = tuple_phases_eabc(tup, alpha, phi)
        records.append({
            "tuple": tup,
            "phases": phases,
            "width": circular_width(phases),
            "gaps": normalized_gap_signature(phases),
        })

    vals = sorted(rec["width"] for rec in records)
    m = len(vals)

    def avg_gap_signature(recs):
        acc = [0.0, 0.0, 0.0, 0.0]
        for rec in recs:
            for i, x in enumerate(rec["gaps"]):
                acc[i] += x
        return [x / len(recs) for x in acc]

    def gap_signature_variance(recs):
        mean_sig = avg_gap_signature(recs)
        var = 0.0
        for rec in recs:
            var += sum((rec["gaps"][i] - mean_sig[i]) ** 2 for i in range(4))
        return var / len(recs)

    return {
        "count": m,
        "mean_width": sum(vals) / m,
        "median_width": vals[m // 2] if m % 2 == 1 else 0.5 * (vals[m // 2 - 1] + vals[m // 2]),
        "fraction_below_0.10": sum(1 for x in vals if x < 0.10) / m,
        "fraction_below_0.20": sum(1 for x in vals if x < 0.20) / m,
        "fraction_below_0.30": sum(1 for x in vals if x < 0.30) / m,
        "avg_gap_signature": avg_gap_signature(records),
        "gap_signature_variance": gap_signature_variance(records),
    }


# -------------------------------------------------
# Vergleich: echte Vierlinge vs zufällige EABC-Primtupel
# -------------------------------------------------

def compare_real_vs_random_eabc(N, alpha, phi, seed=0):
    real_quads = prime_quadruplets(N)
    random_tups = random_family_tuples(N, count=len(real_quads), seed=seed)

    real_as_eabc = [(p, p + 2, p + 6, p + 8) for (p, _, _, _) in real_quads]

    return {
        "real_prime_quadruplets": summarize_general_tuples(real_as_eabc, alpha, phi),
        "random_EABC_prime_tuples": summarize_general_tuples(random_tups, alpha, phi),
    }


# -------------------------------------------------
# Hilfsfunktionen für Primzahl-Stichproben
# -------------------------------------------------

def nearest_prime_index(primes, x):
    lo, hi = 0, len(primes)
    while lo < hi:
        mid = (lo + hi) // 2
        if primes[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    if lo == 0:
        return 0
    if lo == len(primes):
        return len(primes) - 1
    if abs(primes[lo] - x) < abs(primes[lo - 1] - x):
        return lo
    return lo - 1


def random_prime_4tuples_global(N, count=100, seed=0):
    rng = random.Random(seed)
    plist = [p for p in primes_up_to(N) if p > 3]
    out = []
    for _ in range(count):
        tup = tuple(sorted(rng.sample(plist, 4)))
        out.append(tup)
    return out


def random_prime_4tuples_local(real_quads, N, window=5000, seed=0):
    """
    Für jeden echten Vierling werden 4 zufällige Primzahlen aus einem lokalen Fenster
    um die Größenordnung von p gezogen.
    """
    rng = random.Random(seed)
    plist = [p for p in primes_up_to(N) if p > 3]
    out = []

    for quad in real_quads:
        p = quad[0]
        candidates = [q for q in plist if abs(q - p) <= window]
        if len(candidates) < 4:
            idx = nearest_prime_index(plist, p)
            lo = max(0, idx - 100)
            hi = min(len(plist), idx + 101)
            candidates = plist[lo:hi]
        tup = tuple(sorted(rng.sample(candidates, 4)))
        out.append(tup)

    return out


# -------------------------------------------------
# Allgemeine 4-Tupel-Phasen
# -------------------------------------------------

def tuple_phases_plain(tup, alpha, phi_mode="bamberg", phi=None):
    """
    tup = (q1, q2, q3, q4) sortiert.
    phi_mode:
      - 'bamberg' : 0, 1/4, 1/2, 3/4 auf die sortierten Plätze
      - 'allzero' : alle 0
    """
    q1, q2, q3, q4 = tup

    if phi_mode == "bamberg":
        if phi is None:
            phi = [0.0, 0.25, 0.5, 0.75]
        return [
            frac(alpha * q1 + phi[0]),
            frac(alpha * q2 + phi[1]),
            frac(alpha * q3 + phi[2]),
            frac(alpha * q4 + phi[3]),
        ]
    if phi_mode == "allzero":
        return [
            frac(alpha * q1),
            frac(alpha * q2),
            frac(alpha * q3),
            frac(alpha * q4),
        ]
    raise ValueError(f"Unbekannter phi_mode: {phi_mode}")


def summarize_plain_tuples(tuples_list, alpha, phi_mode="bamberg", phi=None):
    if not tuples_list:
        return None

    records = []
    for tup in tuples_list:
        phases = tuple_phases_plain(tup, alpha, phi_mode=phi_mode, phi=phi)
        records.append({
            "tuple": tup,
            "phases": phases,
            "width": circular_width(phases),
            "gaps": normalized_gap_signature(phases),
        })

    vals = sorted(rec["width"] for rec in records)
    m = len(vals)

    def avg_gap_signature(recs):
        acc = [0.0, 0.0, 0.0, 0.0]
        for rec in recs:
            for i, x in enumerate(rec["gaps"]):
                acc[i] += x
        return [x / len(recs) for x in acc]

    def gap_signature_variance(recs):
        mean_sig = avg_gap_signature(recs)
        var = 0.0
        for rec in recs:
            var += sum((rec["gaps"][i] - mean_sig[i]) ** 2 for i in range(4))
        return var / len(recs)

    return {
        "count": m,
        "mean_width": sum(vals) / m,
        "median_width": vals[m // 2] if m % 2 == 1 else 0.5 * (vals[m // 2 - 1] + vals[m // 2]),
        "fraction_below_0.10": sum(1 for x in vals if x < 0.10) / m,
        "fraction_below_0.20": sum(1 for x in vals if x < 0.20) / m,
        "fraction_below_0.30": sum(1 for x in vals if x < 0.30) / m,
        "avg_gap_signature": avg_gap_signature(records),
        "gap_signature_variance": gap_signature_variance(records),
    }


# -------------------------------------------------
# Andere Primzahlmuster
# -------------------------------------------------

def prime_pattern_tuples(N, offsets):
    """
    offsets z.B. [0, 2, 6, 8]
    Liefert alle Tupel (p+o für o in offsets), so dass alle p+o prim und ≤ N sind.
    """
    plist = set(primes_up_to(N + max(offsets) + 10))
    out = []
    for p in sorted(plist):
        if p + max(offsets) <= N and all((p + o) in plist for o in offsets):
            out.append(tuple(p + o for o in offsets))
    return out


# -------------------------------------------------
# Hauptvergleich: Vierlinge vs. Kontrollen
# -------------------------------------------------

def compare_quadruplets_vs_controls(N, alpha, seed=0, local_window=5000):
    real_quads = prime_quadruplets(N)

    global_random = random_prime_4tuples_global(N, count=len(real_quads), seed=seed)
    local_random = random_prime_4tuples_local(real_quads, N, window=local_window, seed=seed)

    pattern_0246 = prime_pattern_tuples(N, [0, 2, 4, 6])
    pattern_04610 = prime_pattern_tuples(N, [0, 4, 6, 10])
    pattern_0268 = real_quads

    return {
        "real_0268": summarize_plain_tuples(pattern_0268, alpha, phi_mode="bamberg"),
        "random_global_4primes": summarize_plain_tuples(global_random, alpha, phi_mode="bamberg"),
        "random_local_4primes": summarize_plain_tuples(local_random, alpha, phi_mode="bamberg"),
        "pattern_0246": summarize_plain_tuples(pattern_0246, alpha, phi_mode="bamberg") if pattern_0246 else None,
        "pattern_04610": summarize_plain_tuples(pattern_04610, alpha, phi_mode="bamberg") if pattern_04610 else None,
    }


def print_comparison_block(title, result):
    print("\n" + "=" * 72)
    print(title)
    for key, val in result.items():
        print(f"\n[{key}]")
        print(val)


# -------------------------------------------------
# Familien-Discrepancy
# -------------------------------------------------

def discrepancy(points, bins=50):
    pts = sorted(points)
    n = len(pts)
    if n == 0:
        return 0.0
    best = 0.0
    grid = [i / bins for i in range(bins + 1)]
    for a in grid:
        for b in grid:
            if b < a:
                continue
            count = sum(1 for x in pts if a <= x <= b)
            expected = (b - a) * n
            best = max(best, abs(count - expected))
    return best / n

def family_discrepancies(N, alpha, phi):
    data = {fam: [] for fam in ["E", "A", "B", "C"]}
    all_data = []
    for p in primes_up_to(N):
        fam = family_mod12(p)
        if fam is None:
            continue
        x = bamberg_phase(p, alpha, phi)
        data[fam].append(x)
        all_data.append(x)
    out = {fam: discrepancy(data[fam]) for fam in data}
    out["ALL"] = discrepancy(all_data)
    return out

# -------------------------------------------------
# Vierlings-Analyse
# -------------------------------------------------

def quadruplet_phase_record(quad, alpha, phi, mode="bamberg", rng=None):
    """
    mode:
      - 'bamberg'   : E,A,B,C fest auf p,p+2,p+6,p+8
      - 'allzero'   : alle Phasenverschiebungen = 0
      - 'permuted'  : phi-Werte zufällig auf E,A,B,C permutiert
      - 'randomphi' : komplett neue Zufalls-phasen
    """
    p, q, r, s = quad

    if mode == "bamberg":
        phases = [
            frac(alpha * p + phi["E"]),
            frac(alpha * q + phi["A"]),
            frac(alpha * r + phi["B"]),
            frac(alpha * s + phi["C"]),
        ]

    elif mode == "allzero":
        phases = [
            frac(alpha * p),
            frac(alpha * q),
            frac(alpha * r),
            frac(alpha * s),
        ]

    elif mode == "permuted":
        if rng is None:
            rng = random.Random(0)
        vals = [phi["E"], phi["A"], phi["B"], phi["C"]]
        rng.shuffle(vals)
        phases = [
            frac(alpha * p + vals[0]),
            frac(alpha * q + vals[1]),
            frac(alpha * r + vals[2]),
            frac(alpha * s + vals[3]),
        ]

    elif mode == "randomphi":
        if rng is None:
            rng = random.Random(0)
        vals = [rng.random() for _ in range(4)]
        phases = [
            frac(alpha * p + vals[0]),
            frac(alpha * q + vals[1]),
            frac(alpha * r + vals[2]),
            frac(alpha * s + vals[3]),
        ]

    else:
        raise ValueError(f"Unbekannter mode: {mode}")

    width = circular_width(phases)
    gaps = normalized_gap_signature(phases)

    return {
        "quad": quad,
        "phases": phases,
        "width": width,
        "gaps": gaps,
    }

def quadruplet_phase_widths(N, alpha, phi, mode="bamberg", seed=0):
    rng = random.Random(seed)
    widths = []
    for quad in prime_quadruplets(N):
        rec = quadruplet_phase_record(quad, alpha, phi, mode=mode, rng=rng)
        widths.append((quad[0], rec["width"], rec["phases"], rec["gaps"]))
    return widths

def average_gap_signature(records):
    m = len(records)
    if m == 0:
        return None
    acc = [0.0, 0.0, 0.0, 0.0]
    for rec in records:
        for i, x in enumerate(rec["gaps"]):
            acc[i] += x
    return [x / m for x in acc]

def gap_signature_variance(records):
    """
    Misst, wie stabil die Vierlingsform ist.
    Kleine Varianz = ähnliche Kreisgestalt über viele Vierlinge.
    """
    m = len(records)
    if m == 0:
        return None
    mean_sig = average_gap_signature(records)
    var = 0.0
    for rec in records:
        var += sum((rec["gaps"][i] - mean_sig[i]) ** 2 for i in range(4))
    return var / m

def summarize_quadruplet_widths(N, alpha, phi, mode="bamberg", seed=0):
    rng = random.Random(seed)
    quads = prime_quadruplets(N)
    records = [quadruplet_phase_record(quad, alpha, phi, mode=mode, rng=rng) for quad in quads]

    if not records:
        return None

    vals = sorted(rec["width"] for rec in records)
    m = len(vals)
    mean = sum(vals) / m
    median = vals[m // 2] if m % 2 == 1 else 0.5 * (vals[m // 2 - 1] + vals[m // 2])

    below_010 = sum(1 for x in vals if x < 0.10) / m
    below_020 = sum(1 for x in vals if x < 0.20) / m
    below_030 = sum(1 for x in vals if x < 0.30) / m

    avg_sig = average_gap_signature(records)
    sig_var = gap_signature_variance(records)

    return {
        "count_quadruplets": m,
        "mean_width": mean,
        "median_width": median,
        "fraction_below_0.10": below_010,
        "fraction_below_0.20": below_020,
        "fraction_below_0.30": below_030,
        "avg_gap_signature": avg_sig,
        "gap_signature_variance": sig_var,
    }


def synthetic_quadruplets(N, count=None):
    quads = [(n, n + 2, n + 6, n + 8) for n in range(1, N - 8)]
    if count is not None:
        quads = quads[:count]
    return quads


def summarize_given_quads(quads, alpha, phi, mode="bamberg", seed=0):
    rng = random.Random(seed)
    records = [quadruplet_phase_record(quad, alpha, phi, mode=mode, rng=rng) for quad in quads]
    if not records:
        return None
    vals = sorted(rec["width"] for rec in records)
    m = len(vals)
    return {
        "count": m,
        "mean_width": sum(vals) / m,
        "median_width": vals[m // 2] if m % 2 == 1 else 0.5 * (vals[m // 2 - 1] + vals[m // 2]),
        "avg_gap_signature": average_gap_signature(records),
        "gap_signature_variance": gap_signature_variance(records),
    }


# -------------------------------------------------
# Alpha-Scan + Nullmodelle
# -------------------------------------------------

def scan_alphas(N, alphas, phi, seed=0):
    rows = []

    for name, alpha in alphas:
        disc = family_discrepancies(N, alpha, phi)

        bamberg = summarize_quadruplet_widths(N, alpha, phi, mode="bamberg", seed=seed)
        allzero = summarize_quadruplet_widths(N, alpha, phi, mode="allzero", seed=seed)
        permuted = summarize_quadruplet_widths(N, alpha, phi, mode="permuted", seed=seed)
        randomphi = summarize_quadruplet_widths(N, alpha, phi, mode="randomphi", seed=seed)

        row = {
            "alpha_name": name,
            "alpha_value": alpha,
            "disc_E": disc["E"],
            "disc_A": disc["A"],
            "disc_B": disc["B"],
            "disc_C": disc["C"],
            "disc_ALL": disc["ALL"],

            "bamberg_mean_width": bamberg["mean_width"] if bamberg else None,
            "bamberg_gap_var": bamberg["gap_signature_variance"] if bamberg else None,

            "allzero_mean_width": allzero["mean_width"] if allzero else None,
            "allzero_gap_var": allzero["gap_signature_variance"] if allzero else None,

            "permuted_mean_width": permuted["mean_width"] if permuted else None,
            "permuted_gap_var": permuted["gap_signature_variance"] if permuted else None,

            "randomphi_mean_width": randomphi["mean_width"] if randomphi else None,
            "randomphi_gap_var": randomphi["gap_signature_variance"] if randomphi else None,
        }
        rows.append(row)

    return rows

def print_scan_table(rows):
    print("\nAlpha-Scan:")
    for row in rows:
        print("-" * 70)
        print(f"alpha = {row['alpha_name']} = {row['alpha_value']}")
        print("Discrepancies:")
        print({
            "E": row["disc_E"],
            "A": row["disc_A"],
            "B": row["disc_B"],
            "C": row["disc_C"],
            "ALL": row["disc_ALL"],
        })
        print("Mean width:")
        print({
            "bamberg": row["bamberg_mean_width"],
            "allzero": row["allzero_mean_width"],
            "permuted": row["permuted_mean_width"],
            "randomphi": row["randomphi_mean_width"],
        })
        print("Gap-signature variance:")
        print({
            "bamberg": row["bamberg_gap_var"],
            "allzero": row["allzero_gap_var"],
            "permuted": row["permuted_gap_var"],
            "randomphi": row["randomphi_gap_var"],
        })

# -------------------------------------------------
# Hauptteil
# -------------------------------------------------

if __name__ == "__main__":
    N = 200000

    phi = {
        "E": 0.0,
        "A": 0.25,
        "B": 0.5,
        "C": 0.75,
    }

    alphas = [
        ("inv_phi", (math.sqrt(5) - 1) / 2),
        ("sqrt2-1", math.sqrt(2) - 1),
        ("sqrt3-1", math.sqrt(3) - 1),
        ("log2", math.log(2)),
        ("log3", math.log(3)),
        ("log(3/2)", math.log(3 / 2)),
        ("alpha_fs", 1 / 137.035999),
        ("inv_alpha_fs_mod1", 137.035999 % 1),
    ]

    alpha = (math.sqrt(5) - 1) / 2

    print("Discrepancies:")
    print(family_discrepancies(N, alpha, phi))

    test_alphas = [
        ("inv_phi", (math.sqrt(5) - 1) / 2),
        ("sqrt2-1", math.sqrt(2) - 1),
        ("log(3/2)", math.log(3 / 2)),
        ("alpha_fs", 1 / 137.035999),
    ]

    for name, alpha_cmp in test_alphas:
        print("\n" + "=" * 72)
        print(f"Vergleich für alpha = {name} = {alpha_cmp}")
        result = compare_real_vs_random_eabc(N, alpha_cmp, phi, seed=0)
        print(result["real_prime_quadruplets"])
        print(result["random_EABC_prime_tuples"])

    control_compare_alphas = [
        ("sqrt2-1", math.sqrt(2) - 1),
        ("log(3/2)", math.log(3 / 2)),
        ("inv_phi", (math.sqrt(5) - 1) / 2),
        ("alpha_fs", 1 / 137.035999),
    ]

    for name, alpha_ctrl in control_compare_alphas:
        ctrl_result = compare_quadruplets_vs_controls(
            N, alpha_ctrl, seed=0, local_window=5000
        )
        print_comparison_block(
            f"Vergleich echter Vierlinge gegen Kontrollen für alpha={name}={alpha_ctrl}",
            ctrl_result,
        )

    real_quads = prime_quadruplets(N)
    fake_quads = synthetic_quadruplets(N, count=len(real_quads))

    print("\nReal prime quadruplets:")
    print(summarize_given_quads(real_quads, alpha, phi, mode="bamberg", seed=0))

    print("\nSynthetic quadruplets:")
    print(summarize_given_quads(fake_quads, alpha, phi, mode="bamberg", seed=0))

    print("\nQuadruplet summary (bamberg):")
    print(summarize_quadruplet_widths(N, alpha, phi, mode="bamberg", seed=0))

    print("\nQuadruplet summary (allzero):")
    print(summarize_quadruplet_widths(N, alpha, phi, mode="allzero", seed=0))

    print("\nQuadruplet summary (permuted):")
    print(summarize_quadruplet_widths(N, alpha, phi, mode="permuted", seed=0))

    print("\nQuadruplet summary (randomphi):")
    print(summarize_quadruplet_widths(N, alpha, phi, mode="randomphi", seed=0))

    rows = scan_alphas(N, alphas, phi, seed=0)
    print_scan_table(rows)
