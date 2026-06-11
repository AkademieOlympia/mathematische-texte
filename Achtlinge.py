#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===========================================================
EABC–Morley–Walter–Zeta Testsystem bis 10^9
===========================================================

Dieses Programm untersucht:

1. Primcluster
   - Zwillinge
   - Drillinge
   - Vierlinge
   - Fünflinge
   - Sechslinge
   - Achtlinge

2. EABC-Holonomie
   - ABCE
   - CEAB

3. Bernoulli-Potenziale

4. einfache Morley-/Walter-Proxys

5. Korrelationen mit Zeta-Nullstellen

-----------------------------------------------------------
BENÖTIGTE DATEI
-----------------------------------------------------------

Rieman100000.txt
oder
Rieman100000.tex

mit einer Nullstelle pro Zeile.

Beispiel:

14.134725141734693790
21.022039638771554993
25.010857580145688763
...

-----------------------------------------------------------
EMPFOHLENE AUSFÜHRUNG
-----------------------------------------------------------

python3 eabc_zeta_test.py

-----------------------------------------------------------
HINWEIS
-----------------------------------------------------------

Bis 10^9 wird segmented sieve empfohlen.
Das Programm unterstützt:

- klassisches sieve bis kleinere Grenzen
- segmented sieve für große Grenzen

===========================================================
"""

import math
import json
import numpy as np
from collections import defaultdict
from pathlib import Path
from datetime import datetime, timezone

# =========================================================
# EINSTELLUNGEN
# =========================================================

ANALYSIS_NMAX = 10**7
PRIME_CACHE_LIMIT = 10**8

USE_SEGMENTED = True

SEGMENT_SIZE = 10**6

ZETA_FILE_CANDIDATES = (
    "zeros6.npy",
)

MAX_ZEROS = 2000
SHUFFLE_RUNS = 200
SHUFFLE_SEED = 12345

PRIME_CACHE_FILE = "primes_cache_1e8.npy"

# =========================================================
# ZETA-NULLSTELLEN LADEN
# =========================================================

def resolve_zeta_file(candidates):

    script_dir = Path(__file__).resolve().parent
    checked = []
    seen = set()

    for name in candidates:

        by_cwd = Path(name)
        resolved_cwd = by_cwd.resolve()

        if resolved_cwd not in seen:
            checked.append(resolved_cwd)
            seen.add(resolved_cwd)

        if by_cwd.is_file():
            return str(by_cwd)

        by_script_dir = script_dir / name
        resolved_script = by_script_dir.resolve()

        if resolved_script not in seen:
            checked.append(resolved_script)
            seen.add(resolved_script)

        if by_script_dir.is_file():
            return str(by_script_dir)

    checked_lines = "\n".join(f"- {p}" for p in checked)
    raise FileNotFoundError(
        "Keine Zeta-Datei gefunden. Gepruefte Pfade:\n"
        f"{checked_lines}\n"
        "Lege eine der erwarteten Dateien ab oder passe ZETA_FILE_CANDIDATES an."
    )


def load_zeros(filename, max_zeros=2000):

    if Path(filename).suffix.lower() == ".npy":
        zeros = np.asarray(np.load(filename, allow_pickle=False), dtype=float).reshape(-1)
        return zeros[:max_zeros]

    zeros = []

    with open(filename, "r", encoding="utf-8") as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            try:
                val = float(line)
                zeros.append(val)

            except:
                continue

            if len(zeros) >= max_zeros:
                break

    return np.array(zeros)


# =========================================================
# BERNOUILLI B2
# =========================================================

def B2(x):

    return x*x - x + 1.0/6.0


# =========================================================
# EABC-KLASSE
# =========================================================

def eabc(n):

    r = n % 12

    if r == 1:
        return "E"

    elif r == 5:
        return "A"

    elif r == 7:
        return "B"

    elif r == 11:
        return "C"

    return None


# =========================================================
# HOLONOMIE VIERLING
# =========================================================

def holonomy(m):

    r = m % 12

    if r == 9:
        return +1     # ABCE

    elif r == 3:
        return -1     # CEAB

    return 0


# ============================================================
# EABC-HOLONOMIE-OPERATOR
# ============================================================

EABC = ["E", "A", "B", "C"]
IDX = {x: i for i, x in enumerate(EABC)}

# Zwei echte Gegenrichtungen
CYCLE_PLUS = ["A", "B", "C", "E"]    # ABCE (vorwaerts)
CYCLE_MINUS = ["A", "E", "C", "B"]   # AECB (rueckwaerts)


def cycle_matrix(cycle):
    """
    Erzeugt gerichtete Uebergangsmatrix eines EABC-Umlaufs.
    """
    M = np.zeros((4, 4), dtype=float)
    for i in range(len(cycle)):
        a = cycle[i]
        b = cycle[(i + 1) % len(cycle)]
        M[IDX[a], IDX[b]] = 1.0
    return M


P_plus = cycle_matrix(CYCLE_PLUS)
P_minus = cycle_matrix(CYCLE_MINUS)

# Chiraler Transportoperator
H_chiral = P_plus - P_minus
# Symmetrischer Bandoperator
H_band = P_plus + P_minus
# Diskreter Laplaceoperator auf dem Band
L_band = H_band.T @ H_band
# Chiraler Laplace-/Energieoperator
L_chiral = H_chiral.T @ H_chiral


def eabc_spectrum_report():
    print()
    print("="*60)
    print("EABC-HOLONOMIE-SPEKTRUM")
    print("="*60)
    for name, M in [
        ("P_plus ABCE", P_plus),
        ("P_minus AECB", P_minus),
        ("H_band = P_plus + P_minus", H_band),
        ("H_chiral = P_plus - P_minus", H_chiral),
        ("L_band", L_band),
        ("L_chiral", L_chiral),
    ]:
        eig = np.linalg.eigvals(M)
        print()
        print(name)
        print("Matrix:")
        print(M)
        print("Eigenwerte:")
        print(np.round(eig, 6))


def project_eabc_mode(Y_plus, Y_minus):
    """
    Projektionen der numerischen Daten auf Band- und Chiralmoden.
    """
    Y_plus = np.array(Y_plus, dtype=float)
    Y_minus = np.array(Y_minus, dtype=float)
    Y_band = Y_plus + Y_minus
    Y_chiral = Y_plus - Y_minus
    return Y_band, Y_chiral


# =========================================================
# SEGMENTED SIEVE
# =========================================================

def small_primes(limit):

    sieve = np.ones(limit + 1, dtype=np.bool_)

    sieve[:2] = False

    for p in range(2, int(limit**0.5) + 1):

        if sieve[p]:
            sieve[p*p::p] = False

    return np.nonzero(sieve)[0]


def segmented_sieve(nmax, segment_size=10**6):

    root = int(math.isqrt(nmax)) + 1

    base_primes = small_primes(root)

    yield 2

    low = 3

    while low <= nmax:

        high = min(low + segment_size - 1, nmax)

        size = high - low + 1

        sieve = np.ones(size, dtype=np.bool_)

        for p in base_primes:

            start = max(p*p, ((low + p - 1)//p)*p)

            for k in range(start, high+1, p):

                sieve[k-low] = False

        for i in range(size):

            if sieve[i]:

                val = low + i

                if val % 2 == 1:
                    yield val

        low = high + 1


# =========================================================
# HASHSET PRIMZAHLEN
# =========================================================

def resolve_prime_cache_file(filename):

    script_dir = Path(__file__).resolve().parent
    return script_dir / filename


def resolve_prime_cache_meta_file(cache_path):

    return cache_path.with_suffix(cache_path.suffix + ".meta.json")


def load_prime_cache_meta(meta_path):

    if not meta_path.is_file():
        return None

    try:
        with open(meta_path, "r", encoding="utf-8") as f:
            meta = json.load(f)
    except Exception:
        return None

    if not isinstance(meta, dict):
        return None

    return meta


def load_prime_cache(cache_path, limit):

    if not cache_path.is_file():
        return None

    arr = np.load(cache_path, allow_pickle=False)
    arr = np.asarray(arr, dtype=np.int64).reshape(-1)

    if arr.size == 0:
        return None

    meta_path = resolve_prime_cache_meta_file(cache_path)
    meta = load_prime_cache_meta(meta_path)

    if meta is not None:
        meta_limit = int(meta.get("target_limit", -1))
        meta_max_prime = int(meta.get("max_prime", -1))
        meta_count = int(meta.get("prime_count", -1))
        arr_max_prime = int(arr[-1])
        arr_count = int(arr.size)

        if meta_limit >= limit and meta_max_prime == arr_max_prime and meta_count == arr_count:
            arr = arr[arr <= limit]
            if arr.size == 0:
                return None
            print(f"Lade Primzahl-Cache (Meta bestaetigt): {cache_path}")
            return set(arr.tolist())

    cached_limit = int(arr[-1])
    if cached_limit < limit:
        # Der letzte Cache-Wert ist die groesste Primzahl <= Limit und kann
        # daher kleiner als das Limit sein, obwohl der Cache vollstaendig ist.
        tail_size = limit - cached_limit
        max_tail_check = 100000

        if tail_size <= max_tail_check:
            base_primes = arr[arr <= int(math.isqrt(limit))].tolist()

            def is_prime_tail(n):
                if n < 2:
                    return False
                if n % 2 == 0:
                    return n == 2
                r = int(math.isqrt(n))
                for p in base_primes:
                    if p > r:
                        break
                    if n % p == 0:
                        return False
                return True

            has_missing_prime = any(
                is_prime_tail(x) for x in range(cached_limit + 1, limit + 1)
            )

            if has_missing_prime:
                print(
                    f"Primzahl-Cache zu klein ({cached_limit} < {limit}), berechne neu..."
                )
                return None

            print(
                "Primzahl-Cache endet unter dem Limit, "
                "Restintervall ist aber primfrei. Nutze Cache."
            )
        else:
            print(
                f"Primzahl-Cache zu klein ({cached_limit} < {limit}), berechne neu..."
            )
            return None

    arr = arr[arr <= limit]

    if arr.size == 0:
        return None

    print(f"Lade Primzahl-Cache: {cache_path}")
    return set(arr.tolist())


def save_prime_cache(cache_path, primes, target_limit):

    arr = np.array(sorted(primes), dtype=np.uint32)
    np.save(cache_path, arr, allow_pickle=False)
    print(f"Primzahl-Cache gespeichert: {cache_path}")

    meta_path = resolve_prime_cache_meta_file(cache_path)
    meta = {
        "target_limit": int(target_limit),
        "max_prime": int(arr[-1]) if arr.size else 0,
        "prime_count": int(arr.size),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)
        f.write("\n")
    print(f"Primzahl-Cache-Meta gespeichert: {meta_path}")


def build_prime_set(limit, cache_filename=None):

    cache_path = None
    if cache_filename:
        cache_path = resolve_prime_cache_file(cache_filename)
        cached = load_prime_cache(cache_path, limit)
        if cached is not None:
            return cached

    primes = set()

    if USE_SEGMENTED:

        for p in segmented_sieve(limit, SEGMENT_SIZE):
            primes.add(p)

    else:

        sieve = np.ones(limit + 1, dtype=np.bool_)
        sieve[:2] = False

        for p in range(2, int(limit**0.5)+1):

            if sieve[p]:
                sieve[p*p::p] = False

        arr = np.nonzero(sieve)[0]

        primes = set(arr.tolist())

    if cache_path is not None:
        save_prime_cache(cache_path, primes, limit)

    return primes


# =========================================================
# CLUSTER-DEFINITIONEN
# =========================================================

CLUSTERS = {

    "twin_left":
        [-4,-2],

    "twin_mid":
        [-2,2],

    "twin_right":
        [2,4],

    "triple_left":
        [-4,-2,2],

    "triple_right":
        [-2,2,4],

    "quad":
        [-4,-2,2,4],
}

CLUSTERS_P = {
    "quint1":
        [0,2,6,8,12],
    "quint2":
        [0,4,6,10,12],
    "six":
        [0,4,6,10,12,16],
}


# =========================================================
# TESTFELDER
# =========================================================

def bernoulli_potential(m, qs=(5,7,11,13,17,19)):

    s = 0.0

    for q in qs:

        x = (m % q)/q

        s += B2(x)

    return s


def morley_proxy(points):

    d = []

    for i in range(len(points)):
        for j in range(i+1, len(points)):
            d.append(abs(points[j]-points[i]))

    d = np.array(d, dtype=float)

    return np.std(d)


def walter_proxy(points):

    diffs = np.diff(points)

    return np.std(diffs)


# =========================================================
# ZETA-SIGNAL
# =========================================================

def zeta_signal(x, zeros):

    lx = math.log(x)

    return np.sum(np.cos(zeros * lx))


def chi4(n):

    r = n % 4
    if r == 1:
        return 1.0
    if r == 3:
        return -1.0
    return 0.0


# =========================================================
# ANALYSE
# =========================================================

def analyze(primes, zeros, analysis_limit):

    counts = defaultdict(int)

    zeta_corr_data = defaultdict(list)

    values = []

    for m in range(15, analysis_limit, 15):

        h = holonomy(m)

        y = 0

        for eps in [-4,-2,2,4]:

            if (m+eps) in primes:
                y += 1

        zsig = zeta_signal(m, zeros)

        values.append((m, y, zsig, h))

        # ---------------------------------------------
        # Cluster testen
        # ---------------------------------------------

        for name, pattern in CLUSTERS.items():

            ok = True

            pts = []

            for d in pattern:

                x = m + d

                if x not in primes:
                    ok = False
                    break

                pts.append(x)

            if ok:

                counts[name] += 1

                B = bernoulli_potential(m)

                M = morley_proxy(pts)

                W = walter_proxy(pts)

                zeta_corr_data[name].append(
                    (B, M, W, h, zsig)
                )

    max_delta = max(max(pattern) for pattern in CLUSTERS_P.values())

    p_candidates = [p for p in primes if 5 <= p < analysis_limit]

    for p in p_candidates:

        if p + max_delta >= analysis_limit:
            continue

        for name, pattern in CLUSTERS_P.items():

            ok = True

            for d in pattern:
                if (p + d) not in primes:
                    ok = False
                    break

            if ok:
                counts[name] += 1

    return counts, zeta_corr_data, values


# =========================================================
# KORRELATION
# =========================================================

def corr(a, b):

    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)

    if len(a) < 3:
        return np.nan

    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan

    return np.corrcoef(a,b)[0,1]


def binned_corr(Y, Z, bin_size=10000):

    Y = np.array(Y, dtype=float)
    Z = np.array(Z, dtype=float)

    n = len(Y) // bin_size

    if n < 3:
        return np.nan

    Yb = Y[:n*bin_size].reshape(n, bin_size).mean(axis=1)
    Zb = Z[:n*bin_size].reshape(n, bin_size).mean(axis=1)

    return np.corrcoef(Yb, Zb)[0,1]


def autocorr_lags(x, max_lag=50):

    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 3:
        return []

    x = x - np.mean(x)
    var = np.var(x)
    if var == 0:
        return []

    lags = []
    upper = min(max_lag, n - 2)
    for lag in range(1, upper + 1):
        c = np.mean(x[:-lag] * x[lag:]) / var
        lags.append((lag, float(c)))
    return lags


def power_spectrum_top(x, top_k=5):

    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 8:
        return []

    x = x - np.mean(x)
    fftv = np.fft.rfft(x)
    power = (np.abs(fftv) ** 2) / n
    freqs = np.fft.rfftfreq(n, d=1.0)

    if power.size <= 1:
        return []

    # DC-Komponente ignorieren.
    idx = np.argsort(power[1:])[::-1] + 1
    idx = idx[:top_k]

    out = []
    for i in idx:
        f = float(freqs[i])
        p = float(power[i])
        period = np.inf if f == 0 else float(1.0 / f)
        out.append((f, period, p))
    return out


def zero_mode_test(field, zeros, m_values, max_zeros_list=(50, 100, 250, 500, 1000, 2000)):
    """
    Testet ein Feld gegen partielle Zeta-Nullstellendetektoren.
    field: reelles oder komplexes Feld auf m_values
    zeros: geladene Nullstellen gamma
    m_values: zugehoerige Mittelpunkte m
    """
    field = np.array(field)
    m_values = np.array(m_values, dtype=float)
    print()
    print("="*60)
    print("NULLSTELLENTEST")
    print("="*60)
    for K in max_zeros_list:
        gammas = zeros[:K]
        Z_cos = np.array([
            np.sum(np.cos(gammas * np.log(m)))
            for m in m_values
        ])
        Z_sin = np.array([
            np.sum(np.sin(gammas * np.log(m)))
            for m in m_values
        ])
        print()
        print(f"Erste {K} Nullstellen")
        if np.iscomplexobj(field):
            F_re = field.real
            F_im = field.imag
            F_abs = np.abs(field)
            F_phase = np.angle(field)
            print("Re(Psi)   vs Z_cos:", corr(F_re, Z_cos))
            print("Im(Psi)   vs Z_sin:", corr(F_im, Z_sin))
            print("|Psi|     vs Z_cos:", corr(F_abs, Z_cos))
            print("arg(Psi)  vs Z_sin:", corr(F_phase, Z_sin))
            for bs in [100, 500, 1000, 2000, 5000, 10000]:
                print(f"  Bin {bs:6d} |Psi| vs Z_cos:",
                      binned_corr(F_abs, Z_cos, bs))
                print(f"  Bin {bs:6d} arg(Psi) vs Z_sin:",
                      binned_corr(F_phase, Z_sin, bs))
        else:
            print("field vs Z_cos:", corr(field, Z_cos))
            print("field vs Z_sin:", corr(field, Z_sin))
            for bs in [100, 500, 1000, 2000, 5000, 10000]:
                print(f"  Bin {bs:6d} field vs Z_cos:",
                      binned_corr(field, Z_cos, bs))
                print(f"  Bin {bs:6d} field vs Z_sin:",
                      binned_corr(field, Z_sin, bs))


def compute_binned_series(Y, Z, bins):

    series = []
    for bs in bins:
        series.append((int(bs), binned_corr(Y, Z, int(bs))))
    return series


def print_binned_series(label, Y, Z, bins):

    print(label)
    series = compute_binned_series(Y, Z, bins)
    for bs, rho in series:
        print(f"  Bin {bs:7d}: {rho}")
    return series


def sign_of(x, eps=1e-12):

    if not np.isfinite(x):
        return 0
    if x > eps:
        return 1
    if x < -eps:
        return -1
    return 0


def sign_intervals(series):

    valid = [(bs, rho) for bs, rho in series if np.isfinite(rho)]
    if not valid:
        return []

    intervals = []
    start_bs = None
    end_bs = None
    current_sign = 0

    for bs, rho in valid:
        s = sign_of(rho)
        if s == 0:
            continue
        if current_sign == 0:
            start_bs = bs
            end_bs = bs
            current_sign = s
            continue
        if s == current_sign:
            end_bs = bs
            continue
        intervals.append((current_sign, start_bs, end_bs))
        start_bs = bs
        end_bs = bs
        current_sign = s

    if current_sign != 0:
        intervals.append((current_sign, start_bs, end_bs))

    return intervals


def sign_transitions(series):

    valid = [(bs, rho) for bs, rho in series if np.isfinite(rho)]
    transitions = []

    prev_sign = 0
    prev_bs = None

    for bs, rho in valid:
        s = sign_of(rho)
        if s == 0:
            continue
        if prev_sign != 0 and s != prev_sign:
            transitions.append((prev_bs, bs, prev_sign, s))
        prev_sign = s
        prev_bs = bs

    return transitions


def print_sign_report(label, series):

    finite = [(bs, rho) for bs, rho in series if np.isfinite(rho)]

    print(label)
    if not finite:
        print("  Keine gueltigen Multi-Scale-Werte.")
        return

    first_positive = next((bs for bs, rho in finite if rho > 0), None)
    first_negative = next((bs for bs, rho in finite if rho < 0), None)

    print("  Erstes positives Bin B*:", first_positive)
    print("  Erstes negatives Bin:", first_negative)

    intervals = sign_intervals(series)
    if intervals:
        for s, b0, b1 in intervals:
            tag = "+" if s > 0 else "-"
            print(f"  Intervall {tag}: [{b0}, {b1}]")
    else:
        print("  Keine stabilen Vorzeichenintervalle erkannt.")

    transitions = sign_transitions(series)
    if transitions:
        for b_prev, b_now, s_prev, s_now in transitions:
            a = "+" if s_prev > 0 else "-"
            b = "+" if s_now > 0 else "-"
            print(f"  Vorzeichenwechsel: {a} -> {b} zwischen {b_prev} und {b_now}")
    else:
        print("  Kein Vorzeichenwechsel gefunden.")


def summarize_null_distribution(obs, null_values):

    null = np.asarray(null_values, dtype=float)
    null = null[np.isfinite(null)]

    if not np.isfinite(obs) or null.size == 0:
        return {
            "n": int(null.size),
            "mean": np.nan,
            "std": np.nan,
            "q05": np.nan,
            "q95": np.nan,
            "p_two_sided": np.nan,
        }

    mean = float(np.mean(null))
    std = float(np.std(null))
    q05 = float(np.quantile(null, 0.05))
    q95 = float(np.quantile(null, 0.95))
    p_two = (1 + int(np.sum(np.abs(null) >= abs(obs)))) / (null.size + 1)

    return {
        "n": int(null.size),
        "mean": mean,
        "std": std,
        "q05": q05,
        "q95": q95,
        "p_two_sided": float(p_two),
    }


def shuffle_significance(Y_plus_full, Y_minus_full, Z, bins, runs=200, seed=12345):

    Y_plus_full = np.asarray(Y_plus_full, dtype=float)
    Y_minus_full = np.asarray(Y_minus_full, dtype=float)
    Z = np.asarray(Z, dtype=float)
    bins = [int(bs) for bs in bins]

    obs_band = {
        bs: binned_corr(Y_plus_full + Y_minus_full, Z, bs)
        for bs in bins
    }
    obs_chiral = {
        bs: binned_corr(Y_plus_full - Y_minus_full, Z, bs)
        for bs in bins
    }

    null_band = {bs: [] for bs in bins}
    null_chiral = {bs: [] for bs in bins}

    rng = np.random.default_rng(seed)

    for _ in range(runs):
        yp = np.array(Y_plus_full, copy=True)
        ym = np.array(Y_minus_full, copy=True)
        rng.shuffle(yp)
        rng.shuffle(ym)

        y_band = yp + ym
        y_chiral = yp - ym

        for bs in bins:
            null_band[bs].append(binned_corr(y_band, Z, bs))
            null_chiral[bs].append(binned_corr(y_chiral, Z, bs))

    stats_band = {
        bs: summarize_null_distribution(obs_band[bs], null_band[bs])
        for bs in bins
    }
    stats_chiral = {
        bs: summarize_null_distribution(obs_chiral[bs], null_chiral[bs])
        for bs in bins
    }

    return obs_band, obs_chiral, stats_band, stats_chiral


def shuffle_test(Y_plus, Y_minus, Z, bin_size=10000, n_shuffles=500, seed=42):

    rng = np.random.default_rng(seed)
    Y_plus = np.array(Y_plus, dtype=float)
    Y_minus = np.array(Y_minus, dtype=float)
    Z = np.array(Z, dtype=float)

    n = min(len(Y_plus), len(Y_minus), len(Z))
    if n < 3:
        return np.nan, np.nan, np.nan, np.nan

    Y_plus = Y_plus[:n]
    Y_minus = Y_minus[:n]
    Z = Z[:n]

    real_diff = Y_plus - Y_minus
    real_corr = binned_corr(real_diff, Z, bin_size)

    combined = np.concatenate([Y_plus, Y_minus])
    sims = []
    for _ in range(n_shuffles):
        rng.shuffle(combined)
        yp = combined[:n]
        ym = combined[n:]
        sims.append(binned_corr(yp - ym, Z, bin_size))

    sims = np.array([x for x in sims if not np.isnan(x)], dtype=float)
    if sims.size == 0 or not np.isfinite(real_corr):
        return real_corr, np.nan, np.nan, np.nan

    p_value = float(np.mean(np.abs(sims) >= abs(real_corr)))
    return real_corr, float(np.mean(sims)), float(np.std(sims)), p_value


def print_shuffle_report(channel_label, obs_map, stats_map, bins):

    print(channel_label)
    for bs in bins:
        obs = obs_map[bs]
        st = stats_map[bs]
        print(
            f"  Bin {bs:7d}: obs={obs}, null_mu={st['mean']}, "
            f"null_q05={st['q05']}, null_q95={st['q95']}, p2={st['p_two_sided']}, n={st['n']}"
        )


def print_shuffle_report_with_flags(channel_label, obs_map, stats_map, bins, alpha=0.05):

    print(channel_label)
    sig_bins = []

    for bs in bins:
        obs = obs_map[bs]
        st = stats_map[bs]
        p2 = st["p_two_sided"]
        is_sig = np.isfinite(p2) and p2 < alpha
        if is_sig:
            sig_bins.append(bs)
        mark = " *" if is_sig else ""
        print(
            f"  Bin {bs:7d}: obs={obs}, null_mu={st['mean']}, "
            f"null_q05={st['q05']}, null_q95={st['q95']}, p2={p2}, n={st['n']}{mark}"
        )

    if sig_bins:
        print(f"  Signifikant (p2 < {alpha}): {sig_bins}")
    else:
        print(f"  Keine signifikanten Skalen (p2 < {alpha}).")


# =========================================================
# MAIN
# =========================================================

def main():

    print("="*60)
    print("Lade Zeta-Nullstellen...")
    print("="*60)

    zeta_file = resolve_zeta_file(ZETA_FILE_CANDIDATES)
    print("Verwendete Datei:", zeta_file)

    zeros = load_zeros(zeta_file, MAX_ZEROS)

    print("Anzahl Nullstellen:", len(zeros))

    print("="*60)
    print("Baue Primzahlmenge...")
    print("="*60)
    print(f"Analyse bis: {ANALYSIS_NMAX}")
    print(f"Cache-Aufbau bis: {PRIME_CACHE_LIMIT}")

    primes = build_prime_set(PRIME_CACHE_LIMIT, PRIME_CACHE_FILE)

    print("Anzahl Primzahlen:", len(primes))

    print("="*60)
    print("Analysiere...")
    print("="*60)

    counts, data, values = analyze(primes, zeros, ANALYSIS_NMAX)

    print()
    print("="*60)
    print("CLUSTER")
    print("="*60)

    for k,v in counts.items():

        print(f"{k:12s}: {v}")

    print()
    print("="*60)
    print("KORRELATIONEN")
    print("="*60)

    for name, arr in data.items():

        if len(arr) < 5:
            continue

        B = [x[0] for x in arr]
        M = [x[1] for x in arr]
        W = [x[2] for x in arr]
        H = [x[3] for x in arr]
        Z = [x[4] for x in arr]

        print()
        print(name)

        print("Bernoulli vs Zeta:",
              corr(B,Z))

        print("Morley vs Zeta:",
              corr(M,Z))

        print("Walter vs Zeta:",
              corr(W,Z))

        print("Holonomie vs Zeta:",
              corr(H,Z))

    print()
    print("="*60)
    print("GESAMTFELD")
    print("="*60)

    M_VALUES = []
    Mvals = [v[0] for v in values]
    M_VALUES.extend(Mvals)
    Y = [v[1] for v in values]
    Z = [v[2] for v in values]
    H = [v[3] for v in values]

    fixed_bins = [1000, 5000, 10000, 50000, 100000]
    multiscale_bins = np.unique(np.logspace(2, 6, 40).astype(int))

    Y_plus = [y for y, h in zip(Y, H) if h == 1]
    Z_plus = [z for z, h in zip(Z, H) if h == 1]

    Y_minus = [y for y, h in zip(Y, H) if h == -1]
    Z_minus = [z for z, h in zip(Z, H) if h == -1]

    # Punktweise Kanaele auf derselben m-Achse:
    # banddichte-sensitiv (Summe) und Phasendifferenzkanal ABCE-CEAB.
    Y_plus_full = [y if h == 1 else 0.0 for y, h in zip(Y, H)]
    Y_minus_full = [y if h == -1 else 0.0 for y, h in zip(Y, H)]
    Y_band = [yp + ym for yp, ym in zip(Y_plus_full, Y_minus_full)]
    Y_chiral = [yp - ym for yp, ym in zip(Y_plus_full, Y_minus_full)]
    Z_chi4 = [chi4(m) * z for m, z in zip(Mvals, Z)]

    print("Primrand vs Zeta roh:",
          corr(Y,Z))

    for bs in fixed_bins:
        print(f"Primrand vs Zeta binned {bs}:",
              binned_corr(Y, Z, bs))

    print()
    print("Multi-Scale Primrand vs Zeta")
    series_all = print_binned_series("", Y, Z, multiscale_bins)
    print_sign_report("Vorzeichenanalyse Primrand", series_all)

    print()
    print("ABCE-Phasenkanal (Y_plus) vs Zeta roh:",
          corr(Y_plus, Z_plus))
    for bs in fixed_bins:
        print(f"ABCE-Phasenkanal vs Zeta binned {bs}:",
              binned_corr(Y_plus, Z_plus, bs))
    series_plus = print_binned_series("ABCE-Phasenkanal Multi-Scale", Y_plus, Z_plus, multiscale_bins)
    print_sign_report("Vorzeichenanalyse ABCE-Phase", series_plus)

    print()
    print("CEAB-Phasenkanal (Y_minus) vs Zeta roh:",
          corr(Y_minus, Z_minus))
    for bs in fixed_bins:
        print(f"CEAB-Phasenkanal vs Zeta binned {bs}:",
              binned_corr(Y_minus, Z_minus, bs))
    series_minus = print_binned_series("CEAB-Phasenkanal Multi-Scale", Y_minus, Z_minus, multiscale_bins)
    print_sign_report("Vorzeichenanalyse CEAB-Phase", series_minus)

    print()
    print("Bandkanal Y_band = Y_ABCE + Y_CEAB vs Zeta roh:",
          corr(Y_band, Z))
    print("Phasendifferenzkanal Y_phase = Y_ABCE - Y_CEAB vs Zeta roh:",
          corr(Y_chiral, Z))
    for bs in fixed_bins:
        print(f"Band-Summe vs Zeta binned {bs}:",
              binned_corr(Y_band, Z, bs))
        print(f"Phasendifferenz vs Zeta binned {bs}:",
              binned_corr(Y_chiral, Z, bs))

    series_band = print_binned_series("Band-Summe Multi-Scale", Y_band, Z, multiscale_bins)
    print_sign_report("Vorzeichenanalyse Band-Summe", series_band)

    series_chiral = print_binned_series("Phasendifferenz Multi-Scale", Y_chiral, Z, multiscale_bins)
    print_sign_report("Vorzeichenanalyse Phasendifferenz", series_chiral)

    print()
    print("Phasendifferenzkanal vs Chi4-Zeta-Detektor roh:",
          corr(Y_chiral, Z_chi4))
    for bs in fixed_bins:
        print(f"Phasendifferenzkanal vs Chi4-Zeta binned {bs}:",
              binned_corr(Y_chiral, Z_chi4, bs))
    series_chiral_chi4 = print_binned_series(
        "Phasendifferenzkanal vs Chi4-Zeta Multi-Scale",
        Y_chiral,
        Z_chi4,
        multiscale_bins
    )
    print_sign_report("Vorzeichenanalyse Phasendifferenz vs Chi4-Zeta", series_chiral_chi4)

    eabc_spectrum_report()

    Y_band, Y_phase = project_eabc_mode(Y_plus_full, Y_minus_full)
    print()
    print("="*60)
    print("EABC-HOLONOMIEKANAL")
    print("="*60)
    print("Bandenergie:", np.mean(Y_band**2))
    print("Phasendifferenzenergie:", np.mean(Y_phase**2))
    print("Band vs Zeta roh:", corr(Y_band, Z))
    print("Phasendifferenz vs Zeta roh:", corr(Y_phase, Z))
    for bs in [100, 500, 1000, 2000, 5000, 10000]:
        print(f"Bin {bs:6d} Band vs Zeta:",
              binned_corr(Y_band, Z, bs))
        print(f"Bin {bs:6d} Phasendifferenz vs Zeta:",
              binned_corr(Y_phase, Z, bs))

    print()
    print("="*60)
    print("KOMPLEXES PRIMVIERLINGSFELD")
    print("="*60)
    Psi = Y_band + 1j * Y_phase
    amplitude = np.abs(Psi)
    phase = np.angle(Psi)
    phi = np.unwrap(np.angle(Psi))
    dphi = np.diff(phi)

    print("Amplitude vs Zeta roh:", corr(amplitude, Z))
    print("Phase vs Zeta roh:", corr(phase, Z))
    print("dphi mean:", np.mean(dphi))
    print("dphi std:", np.std(dphi))
    print("dphi vs Zeta roh:", corr(dphi, np.asarray(Z)[1:]))

    ac = autocorr_lags(dphi, max_lag=40)
    print("dphi ACF (erste Lags):")
    if ac:
        for lag, val in ac[:10]:
            print(f"  lag {lag:2d}: {val:+.4f}")
    else:
        print("  keine stabile ACF berechenbar")

    spec = power_spectrum_top(dphi, top_k=6)
    print("dphi Spektrum (Top-Frequenzen):")
    if spec:
        for f, per, p in spec:
            print(f"  f={f:.6f}, Periode={per:.2f}, Power={p:.6f}")
    else:
        print("  kein stabiles Spektrum berechenbar")

    for bs in [100, 500, 1000, 2000, 5000, 10000]:
        print(f"Bin {bs:6d} Amplitude vs Zeta:",
              binned_corr(amplitude, Z, bs))
        print(f"Bin {bs:6d} Phase vs Zeta:",
              binned_corr(phase, Z, bs))
        print(f"Bin {bs:6d} dphi vs Zeta:",
              binned_corr(dphi, np.asarray(Z)[1:], bs))

    zero_mode_test(Y_band, zeros, M_VALUES)
    zero_mode_test(Y_phase, zeros, M_VALUES)
    zero_mode_test(Psi, zeros, M_VALUES)

    print()
    print("="*60)
    print("SIGNIFIKANZKONTROLLE (SHUFFLE)")
    print("="*60)
    print(f"Runs: {SHUFFLE_RUNS}, Seed: {SHUFFLE_SEED}")
    print("Shuffle: Y_ABCE und Y_CEAB werden getrennt entlang m permutiert (Trefferzahlen bleiben erhalten).")

    obs_band, obs_chiral, stats_band, stats_chiral = shuffle_significance(
        Y_plus_full,
        Y_minus_full,
        Z,
        fixed_bins,
        runs=SHUFFLE_RUNS,
        seed=SHUFFLE_SEED
    )

    print_shuffle_report("Band-Summe (Y_band) vs Zeta", obs_band, stats_band, fixed_bins)
    print_shuffle_report("Phasendifferenz (Y_phase) vs Zeta", obs_chiral, stats_chiral, fixed_bins)

    print()
    print("Multi-Scale Shuffle (Band-Summe / Phasendifferenz)")
    obs_band_ms, obs_chiral_ms, stats_band_ms, stats_chiral_ms = shuffle_significance(
        Y_plus_full,
        Y_minus_full,
        Z,
        multiscale_bins,
        runs=SHUFFLE_RUNS,
        seed=SHUFFLE_SEED
    )
    print_shuffle_report_with_flags(
        "Band-Summe Multi-Scale vs Zeta",
        obs_band_ms,
        stats_band_ms,
        multiscale_bins,
        alpha=0.05
    )
    print_shuffle_report_with_flags(
        "Phasendifferenz Multi-Scale vs Zeta",
        obs_chiral_ms,
        stats_chiral_ms,
        multiscale_bins,
        alpha=0.05
    )

    print()
    print("Shuffle-Test Phasendifferenz")
    for bs in [1000, 5000, 10000, 20000, 30000, 50000]:
        real, mu, sd, p = shuffle_test(
            Y_plus_full,
            Y_minus_full,
            Z,
            bs,
            n_shuffles=500,
            seed=42
        )
        print(
            f"Bin {bs:6d}: real={real:+.4f}, shuffle={mu:+.4f} ± {sd:.4f}, p={p:.4f}"
        )

    print()
    print("Shuffle-Test Phasendifferenz (Chi4-Zeta-Detektor)")
    for bs in [1000, 5000, 10000, 20000, 30000, 50000]:
        real, mu, sd, p = shuffle_test(
            Y_plus_full,
            Y_minus_full,
            Z_chi4,
            bs,
            n_shuffles=500,
            seed=42
        )
        print(
            f"Bin {bs:6d}: real={real:+.4f}, shuffle={mu:+.4f} ± {sd:.4f}, p={p:.4f}"
        )

    print()
    print("FERTIG.")
    print("="*60)


if __name__ == "__main__":

    main()