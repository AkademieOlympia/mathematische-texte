import csv
import math
import random
import sys
from pathlib import Path
from statistics import mean

# =========================================================
# Basis
# =========================================================

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
    return None


def is_finite_number(x):
    return x is not None and not math.isnan(x) and math.isfinite(x)


def canonical_eabc(tup):
    d = {}
    for q in tup:
        fam = family_mod12(q)
        if fam is None:
            return None
        d[fam] = q
    if set(d.keys()) != {"E", "A", "B", "C"}:
        return None
    return (d["E"], d["A"], d["B"], d["C"])


def eabc_compatible(tup):
    return canonical_eabc(tup) is not None


def stretch_k_eabc_feasible(k):
    """
    True genau dann, wenn es p gibt, sodass (p, p+2k, p+6k, p+8k) mod 12
    die vier Familien {1,5,7,11} (E/A/B/C) je einmal trifft.
    Für andere k ist stretched_quadruplets_eabc_only leer — oder ältere CSVs
    können Zeilen ohne Phasen enthalten (dann mean_width=nan).
    """
    target = frozenset({1, 5, 7, 11})
    offs = (0, 2 * k, 6 * k, 8 * k)
    for p12 in range(12):
        if {(p12 + d) % 12 for d in offs} == target:
            return True
    return False


def circular_gaps(phases):
    pts = sorted(phases)
    ext = pts + [pts[0] + 1.0]
    return [ext[i + 1] - ext[i] for i in range(len(pts))]


def circular_width(phases):
    return 1.0 - max(circular_gaps(phases))


def coherence(phases):
    z = sum(complex(math.cos(2 * math.pi * t), math.sin(2 * math.pi * t)) for t in phases) / len(phases)
    return abs(z)


def phase_vector_eabc(tup, alpha, phi):
    eabc = canonical_eabc(tup)
    if eabc is None:
        return None
    qE, qA, qB, qC = eabc
    return [
        frac(alpha * qE + phi["E"]),
        frac(alpha * qA + phi["A"]),
        frac(alpha * qB + phi["B"]),
        frac(alpha * qC + phi["C"]),
    ]


# =========================================================
# Primzahlen
# =========================================================

def primes_up_to(N):
    sieve = [True] * (N + 1)
    sieve[0:2] = [False, False]
    for p in range(2, int(N**0.5) + 1):
        if sieve[p]:
            sieve[p * p : N + 1 : p] = [False] * (((N - p * p) // p) + 1)
    return [p for p in range(N + 1) if sieve[p]]


def is_prime_small(n, prime_set=None):
    if n < 2:
        return False
    if prime_set is not None:
        return n in prime_set
    if n % 2 == 0:
        return n == 2
    r = int(n**0.5)
    f = 3
    while f <= r:
        if n % f == 0:
            return False
        f += 2
    return True


def prime_quadruplets(N):
    plist = primes_up_to(N + 10)
    pset = set(plist)
    out = []
    for p in plist:
        if p + 8 <= N and all(q in pset for q in (p, p + 2, p + 6, p + 8)):
            out.append((p, p + 2, p + 6, p + 8))
    return out


# =========================================================
# Feature-Hilfen
# =========================================================

def sieve_score_tuple(tup, small_primes):
    total = 0
    good = 0
    for r in small_primes:
        if r in (2, 3):
            continue
        total += 1
        if all(x % r != 0 for x in tup):
            good += 1
    return good / total if total else 0.0


def product_fourth_root(tup):
    P = 1
    for x in tup:
        P *= x
    return P ** 0.25


def midpoint(tup):
    return 0.5 * (min(tup) + max(tup))


def radius_error_stretched(tup):
    R = product_fourth_root(tup)
    m = midpoint(tup)
    return abs(R - m)


def tuple_is_real_prime_quadruplet(tup, prime_set=None):
    a = sorted(tup)
    if len(a) != 4:
        return False
    diffs = [a[i + 1] - a[i] for i in range(3)]
    if diffs != [2, 4, 2]:
        return False
    return all(is_prime_small(x, prime_set) for x in a)


# =========================================================
# Gedehnte Vierlinge
# =========================================================

def stretched_quadruplets_eabc_only(N, k_values):
    """
    Nur EABC-kompatible gedehnte Vierlinge:
    (p, p+2k, p+6k, p+8k)
    """
    out = []
    for k in k_values:
        max_p = N - 8 * k
        for p in range(5, max_p + 1, 2):
            tup = (p, p + 2 * k, p + 6 * k, p + 8 * k)
            if eabc_compatible(tup):
                out.append((k, tup))
    return out


def sample_random_4tuples_eabc_only(N, count, seed=0):
    rng = random.Random(seed)
    plist = [p for p in primes_up_to(N) if p > 3]
    out = []
    trials = 0
    while len(out) < count and trials < 50 * count:
        trials += 1
        tup = tuple(sorted(rng.sample(plist, 4)))
        if eabc_compatible(tup):
            out.append(tup)
    return out


# =========================================================
# Feature-Bau
# =========================================================

def build_feature_row(kind, tup, alpha, phi, small_primes, k=None, prime_set=None):
    phases = phase_vector_eabc(tup, alpha, phi)
    gaps = sorted(circular_gaps(phases)) if phases is not None else None

    row = {
        "kind": kind,
        "k": k if k is not None else "",
        "q1": tup[0],
        "q2": tup[1],
        "q3": tup[2],
        "q4": tup[3],
        "midpoint": midpoint(tup),
        "root4_product": product_fourth_root(tup),
        "radius_error": radius_error_stretched(tup),
        "sieve_score": sieve_score_tuple(tup, small_primes),
        "is_real_prime_quadruplet": 1 if tuple_is_real_prime_quadruplet(tup, prime_set=prime_set) else 0,
    }

    if phases is None:
        if kind == "stretched_quad":
            raise ValueError(
                "stretched_quad ohne EABC-Kanone — tup muss eabc_compatible sein: "
                f"tup={tup}, k={k}"
            )
        row.update({
            "phase_E": "",
            "phase_A": "",
            "phase_B": "",
            "phase_C": "",
            "width": "",
            "coherence": "",
            "gap1": "",
            "gap2": "",
            "gap3": "",
            "gap4": "",
        })
    else:
        row.update({
            "phase_E": phases[0],
            "phase_A": phases[1],
            "phase_B": phases[2],
            "phase_C": phases[3],
            "width": circular_width(phases),
            "coherence": coherence(phases),
            "gap1": gaps[0],
            "gap2": gaps[1],
            "gap3": gaps[2],
            "gap4": gaps[3],
        })

    return row


# =========================================================
# Datensatz bauen
# =========================================================

def _resolve_data_path(p: str | Path) -> Path:
    path = Path(p)
    return path if path.is_absolute() else Path(__file__).resolve().parent / path


def build_dataset(
    N,
    alpha,
    phi,
    outfile="vierlings_dataset_clean.csv",
    k_values=(1, 5, 7, 11, 13),
    max_real=None,
    max_stretched=None,
    max_random=None,
    sieve_bound=100,
    seed=0,
):
    rng = random.Random(seed)
    path = _resolve_data_path(outfile)

    skipped_k = [k for k in k_values if not stretch_k_eabc_feasible(k)]
    k_values = [k for k in k_values if stretch_k_eabc_feasible(k)]
    if skipped_k:
        print(
            "[Vierlings Sonde] k ohne gültige EABC-Mod12-Struktur übersprungen: "
            + ", ".join(str(k) for k in skipped_k),
            file=sys.stderr,
        )

    plist = primes_up_to(N + 10)
    pset = set(plist)
    small_prs = primes_up_to(sieve_bound)

    real_quads = [t for t in prime_quadruplets(N) if eabc_compatible(t)]
    if max_real is not None:
        real_quads = real_quads[:max_real]

    stretched = stretched_quadruplets_eabc_only(N, k_values)
    # echte Vierlinge als stretched k=1 sind okay; optional herausfiltern:
    # stretched = [(k, t) for (k, t) in stretched if not (k == 1 and tuple_is_real_prime_quadruplet(t, pset))]
    if max_stretched is not None and len(stretched) > max_stretched:
        stretched = rng.sample(stretched, max_stretched)

    count_rand = len(real_quads) if max_random is None else max_random
    random_tups = sample_random_4tuples_eabc_only(N, count=count_rand, seed=seed)

    rows = []

    for tup in real_quads:
        rows.append(build_feature_row(
            kind="real_quad",
            tup=tup,
            alpha=alpha,
            phi=phi,
            small_primes=small_prs,
            k=1,
            prime_set=pset,
        ))

    for k, tup in stretched:
        rows.append(build_feature_row(
            kind="stretched_quad",
            tup=tup,
            alpha=alpha,
            phi=phi,
            small_primes=small_prs,
            k=k,
            prime_set=pset,
        ))

    for tup in random_tups:
        rows.append(build_feature_row(
            kind="random_4tuple",
            tup=tup,
            alpha=alpha,
            phi=phi,
            small_primes=small_prs,
            k="",
            prime_set=pset,
        ))

    fieldnames = list(rows[0].keys()) if rows else []
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return {
        "outfile": str(path),
        "num_real": len(real_quads),
        "num_stretched": len(stretched),
        "num_random": len(random_tups),
        "total_rows": len(rows),
    }


# =========================================================
# Auswertung
# =========================================================

def load_dataset(filename="vierlings_dataset_clean.csv"):
    path = _resolve_data_path(filename)
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def to_float(x):
    if x == "" or x is None:
        return None
    y = float(x)
    if math.isnan(y) or not math.isfinite(y):
        return None
    return y


def to_int(x):
    if x == "" or x is None:
        return None
    return int(x)


def prepare_rows(rows):
    out = []
    for r in rows:
        out.append({
            "kind": r["kind"],
            "k": None if r.get("k", "") == "" else float(r["k"]),
            "midpoint": to_float(r.get("midpoint")),
            "root4_product": to_float(r.get("root4_product")),
            "radius_error": to_float(r.get("radius_error")),
            "sieve_score": to_float(r.get("sieve_score")),
            "width": to_float(r.get("width")),
            "coherence": to_float(r.get("coherence")),
            "gap1": to_float(r.get("gap1")),
            "gap2": to_float(r.get("gap2")),
            "gap3": to_float(r.get("gap3")),
            "gap4": to_float(r.get("gap4")),
            "is_real_prime_quadruplet": to_int(r.get("is_real_prime_quadruplet")),
        })
    return out


def safe_mean(vals):
    vals = [x for x in vals if is_finite_number(x)]
    return mean(vals) if vals else None


def effect_size(real_vals, other_vals):
    real_vals = [x for x in real_vals if is_finite_number(x)]
    other_vals = [x for x in other_vals if is_finite_number(x)]
    if len(real_vals) < 2 or len(other_vals) < 2:
        return None
    mu1 = mean(real_vals)
    mu2 = mean(other_vals)
    v1 = sum((x - mu1) ** 2 for x in real_vals) / (len(real_vals) - 1)
    v2 = sum((x - mu2) ** 2 for x in other_vals) / (len(other_vals) - 1)
    pooled = math.sqrt((v1 + v2) / 2.0)
    if pooled == 0:
        return float("inf") if mu1 != mu2 else 0.0
    return abs(mu1 - mu2) / pooled


def _fmt_mu(mu):
    return f"{mu:10.6f}" if mu is not None else "       —   "


def summarize_by_kind(rows, features):
    groups = {}
    for r in rows:
        groups.setdefault(r["kind"], []).append(r)

    print("\nZusammenfassung nach Klassen:")
    for kind, chunk in groups.items():
        print("-" * 70)
        print(kind, f"(n={len(chunk)})")
        for feat in features:
            vals = [r[feat] for r in chunk]
            mu = safe_mean(vals)
            if mu is not None:
                print(f"{feat:15s}: mean = {mu:.6f}")


def rank_features(rows, features, other_kind=None):
    real = [r for r in rows if r["kind"] == "real_quad"]
    if other_kind is None:
        other = [r for r in rows if r["kind"] != "real_quad"]
        title = "echte Vierlinge vs alle anderen"
    else:
        other = [r for r in rows if r["kind"] == other_kind]
        title = f"echte Vierlinge vs {other_kind}"

    scores = []
    for feat in features:
        real_vals = [r[feat] for r in real]
        other_vals = [r[feat] for r in other]
        score = effect_size(real_vals, other_vals)
        if score is not None:
            scores.append((feat, score, safe_mean(real_vals), safe_mean(other_vals)))
    scores.sort(key=lambda x: x[1], reverse=True)

    print(f"\nFeature-Ranking: {title}")
    print("-" * 70)
    for feat, score, mu_real, mu_other in scores:
        print(
            f"{feat:15s}  effect={score:8.4f}   real={_fmt_mu(mu_real)}   other={_fmt_mu(mu_other)}"
        )


# =========================================================
# Kombinierter Score
# =========================================================

def combined_score(r):
    """
    Startscore auf Basis Ihrer bisherigen Resultate:
    starkes Gewicht auf sieve_score,
    dann coherence, width, gap4.
    """
    if not all(is_finite_number(r[x]) for x in ("sieve_score", "coherence", "width", "gap4")):
        return None
    return (
        5.0 * r["sieve_score"]
        + 2.0 * r["coherence"]
        - 2.0 * r["width"]
        + 1.0 * r["gap4"]
    )


def rank_by_combined_score(rows, target_kind="real_quad", top_n=30):
    scored = []
    for r in rows:
        s = combined_score(r)
        if s is not None:
            rr = dict(r)
            rr["combined_score"] = s
            scored.append(rr)

    scored.sort(key=lambda x: x["combined_score"], reverse=True)

    print(f"\nTop {top_n} nach kombiniertem Score")
    print("-" * 70)
    hits = 0
    for i, r in enumerate(scored[:top_n], start=1):
        is_hit = (r["kind"] == target_kind)
        hits += int(is_hit)
        print(
            f"{i:2d}. kind={r['kind']:14s} score={r['combined_score']:.6f} "
            f"sieve={r['sieve_score']:.6f} coh={r['coherence']:.6f} "
            f"width={r['width']:.6f} gap4={r['gap4']:.6f} k={r['k']}"
        )
    print("-" * 70)
    print(f"Treffer {target_kind} unter Top {top_n}: {hits}")


def best_k_values(rows):
    stretched = [r for r in rows if r["kind"] == "stretched_quad" and r["k"] is not None]
    by_k = {}
    for r in stretched:
        by_k.setdefault(int(r["k"]), []).append(r)

    print("\nEABC-kompatible gedehnte Vierlinge nach k:")
    print("-" * 70)
    for k in sorted(by_k):
        chunk = by_k[k]
        n = len(chunk)
        widths = [r["width"] for r in chunk]
        cohs = [r["coherence"] for r in chunk]
        sieve = [r["sieve_score"] for r in chunk]
        nw = sum(1 for x in widths if is_finite_number(x))
        nc = sum(1 for x in cohs if is_finite_number(x))
        mw = safe_mean(widths)
        mc = safe_mean(cohs)
        ms = safe_mean(sieve)
        w_s = f"{mw:.6f}" if mw is not None else "—"
        c_s = f"{mc:.6f}" if mc is not None else "—"
        s_s = f"{ms:.6f}" if ms is not None else "—"
        phase_hint = ""
        if nw < n or nc < n:
            phase_hint = (
                f"   (Phasen nur {nw}/{n} width, {nc}/{n} coherence — ggf. alte CSV "
                "oder k ohne EABC-Struktur)"
            )
        print(
            f"k={k:2d}   n={n:5d}   mean_width={w_s}   mean_coherence={c_s}   "
            f"mean_sieve={s_s}{phase_hint}"
        )


# =========================================================
# Beispiel
# =========================================================

if __name__ == "__main__":
    phi = {"E": 0.0, "A": 0.25, "B": 0.5, "C": 0.75}
    alpha = math.sqrt(2) - 1
    # alternativ: alpha = math.log(3 / 2)

    stats = build_dataset(
        N=200000,
        alpha=alpha,
        phi=phi,
        outfile="vierlings_dataset_clean.csv",
        k_values=(1, 5, 7, 11, 13),
        max_real=None,
        max_stretched=20000,
        max_random=5000,
        sieve_bound=100,
        seed=0,
    )
    print(stats)

    rows = load_dataset("vierlings_dataset_clean.csv")
    rows = prepare_rows(rows)

    features = [
        "sieve_score",
        "radius_error",
        "width",
        "coherence",
        "gap1",
        "gap2",
        "gap3",
        "gap4",
    ]

    summarize_by_kind(rows, features)
    rank_features(rows, features, other_kind=None)
    rank_features(rows, features, other_kind="stretched_quad")
    best_k_values(rows)
    rank_by_combined_score(rows, target_kind="real_quad", top_n=30)
