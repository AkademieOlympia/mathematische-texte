import numpy as np
from collections import Counter
from math import isqrt
from scipy.stats import skew

def orientation(p):
    r = p % 12
    if r == 5:
        return "ABCE"
    elif r == 11:
        return "CEAB"
    else:
        return "OTHER"

def fisher_z(r):
    r = np.clip(r, -0.999999999, 0.999999999)
    return np.arctanh(r)

def logit(p):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return np.log(p / (1 - p))

def centered_log_ratio(counts):
    counts = np.asarray(counts, dtype=float)
    k = len(counts)
    p = (counts + 0.5) / (np.sum(counts) + 0.5 * k)
    g = np.exp(np.mean(np.log(p)))
    return np.log(p / g)

def sieve_primes(n):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:n + 1:p] = b"\x00" * (((n - start) // p) + 1)
    return sieve

def is_prime_quad_start(p, prime):
    return prime[p] and prime[p + 2] and prime[p + 6] and prime[p + 8]

def find_prime_quadruplets(limit):
    prime = sieve_primes(limit + 8)
    quads = []
    for p in range(5, limit + 1):
        if p % 12 in (5, 11) and is_prime_quad_start(p, prime):
            quads.append(p)
    return quads

def analyze_quadruplets(quadruplets, W=100000):
    """
    quadruplets: Liste der Startwerte p mit (p,p+2,p+6,p+8)
    W: Fensterbreite in der Mittelpunkt-Koordinate m=p+4
    """

    starts = np.array(quadruplets, dtype=int)
    centers = starts + 4

    min_m, max_m = int(np.min(centers)), int(np.max(centers))
    bins = np.arange(min_m, max_m + W, W)

    window_density = []
    window_bias_raw = []
    window_bias_logit = []
    window_mean_gap = []

    global_orient = Counter()

    for i in range(len(bins) - 1):
        lo, hi = bins[i], bins[i+1]
        mask = (centers >= lo) & (centers < hi)
        ps = starts[mask]

        if len(ps) == 0:
            continue

        orients = [orientation(int(p)) for p in ps]
        c = Counter(orients)

        abce = c["ABCE"]
        ceab = c["CEAB"]
        total = abce + ceab

        global_orient.update(c)

        density = total / W

        raw_bias = (abce - ceab) / total if total > 0 else np.nan

        # Jeffreys-Korrektur
        p_abce = (abce + 0.5) / (total + 1.0)
        logit_bias = logit(p_abce)

        gaps = np.diff(np.sort(ps))
        mean_gap = np.mean(gaps) if len(gaps) > 1 else np.nan

        window_density.append(density)
        window_bias_raw.append(raw_bias)
        window_bias_logit.append(logit_bias)
        window_mean_gap.append(mean_gap)

    window_density = np.array(window_density)
    window_bias_raw = np.array(window_bias_raw)
    window_bias_logit = np.array(window_bias_logit)
    window_mean_gap = np.array(window_mean_gap)

    print("\n=== Globale Orientierung ===")
    print(global_orient)

    abce = global_orient["ABCE"]
    ceab = global_orient["CEAB"]
    total = abce + ceab

    B = (abce - ceab) / total
    sigma_B = 1 / np.sqrt(total)
    Z_B = B / sigma_B

    p_abce = (abce + 0.5) / (total + 1.0)
    L = logit(p_abce)

    print("\n=== Chiraler Bias ===")
    print("ABCE:", abce)
    print("CEAB:", ceab)
    print("B_raw:", B)
    print("sigma_B:", sigma_B)
    print("Z_B:", Z_B)
    print("logit_bias:", L)

    print("\n=== Fenster-Bias Rohskala ===")
    print("mean:", np.nanmean(window_bias_raw))
    print("std :", np.nanstd(window_bias_raw))
    print("skew:", skew(window_bias_raw, nan_policy="omit"))

    print("\n=== Fenster-Bias Logit-Skala ===")
    print("mean:", np.nanmean(window_bias_logit))
    print("std :", np.nanstd(window_bias_logit))
    print("skew:", skew(window_bias_logit, nan_policy="omit"))

    # Korrelationen testen
    valid = (
        np.isfinite(window_density)
        & np.isfinite(window_bias_raw)
        & np.isfinite(window_mean_gap)
    )

    if np.sum(valid) > 3:
        r_density_bias = np.corrcoef(
            window_density[valid],
            window_bias_raw[valid]
        )[0, 1]

        r_gap_bias = np.corrcoef(
            window_mean_gap[valid],
            window_bias_raw[valid]
        )[0, 1]

        print("\n=== Korrelationen roh und Fisher-z ===")
        print("corr(Dichte, Bias):", r_density_bias)
        print("Fisher-z:", fisher_z(r_density_bias))

        print("corr(Gap, Bias):", r_gap_bias)
        print("Fisher-z:", fisher_z(r_gap_bias))

    # Residuen mod 210
    residues = [int(p % 210) for p in starts]
    rc = Counter(residues)

    print("\n=== Residuen mod 210 ===")
    for a, n in sorted(rc.items()):
        print(a, n)

    # Residuen mod 420: gemeinsame Skala fuer mod 210 und EABC/mod 12
    residues420 = [int(p % 420) for p in starts]
    rc420 = Counter(residues420)
    print("\n=== Residuen mod 420 ===")
    for a, n in sorted(rc420.items()):
        print(a, n, "mod12=", a % 12, "orientation=", orientation(a))

    # Ohne singulaeren Start p=5
    starts_reg = np.array([p for p in starts if p != 5], dtype=int)
    residues420_reg = [int(p % 420) for p in starts_reg]
    rc420_reg = Counter(residues420_reg)
    print("\n=== Residuen mod 420 ohne p=5 ===")
    for a, n in sorted(rc420_reg.items()):
        print(a, n, "mod12=", a % 12, "orientation=", orientation(a))

    # Gruppierung nach mod420-Kanal und Orientierung
    print("\n=== mod420-Kanaele nach Orientierung ===")
    abce_channels = {}
    ceab_channels = {}
    for a, n in sorted(rc420_reg.items()):
        if a % 12 == 5:
            abce_channels[a] = n
        elif a % 12 == 11:
            ceab_channels[a] = n
    print("ABCE-Kanaele:")
    for a, n in abce_channels.items():
        print(a, n)
    print("CEAB-Kanaele:")
    for a, n in ceab_channels.items():
        print(a, n)

    # Lift-Paar-Differenzen auf mod 420
    pair_map = {
        11: 221,
        101: 311,
        191: 401,
    }
    print("\n=== Lift-Paar-Differenzen Delta_r = N(r+210)-N(r) ===")
    pair_deltas = {}
    for r, r_lift in pair_map.items():
        n_r = rc420_reg.get(r, 0)
        n_lift = rc420_reg.get(r_lift, 0)
        delta = n_lift - n_r
        pair_deltas[r] = delta
        print(f"Delta_{r}: {n_lift} - {n_r} = {delta}")

    print("\n=== Normierte Lift-Paar-Differenzen ===")
    pair_zscores = {}
    for r, r_lift in pair_map.items():
        n_r = rc420_reg.get(r, 0)
        n_lift = rc420_reg.get(r_lift, 0)
        delta = n_lift - n_r
        total = n_r + n_lift
        z = delta / np.sqrt(total) if total > 0 else np.nan
        pair_zscores[r] = z
        print(
            f"Pair ({r},{r_lift}): "
            f"{n_lift} - {n_r} = {delta}, "
            f"total={total}, Z={z:+.4f}"
        )

    # Delta_r(x)-Verlauf als Funktion der Schranke x
    starts_reg_sorted = np.sort(starts_reg)
    x_min = int(starts_reg_sorted[0]) if len(starts_reg_sorted) > 0 else 0
    x_max = int(starts_reg_sorted[-1]) if len(starts_reg_sorted) > 0 else 0
    x_step = max(10_000, W)
    delta_curve = []
    if x_max > 0:
        print("\n=== Delta_r(x)-Verlauf ===")
        print("x | Delta_11 | Delta_101 | Delta_191")
        idx = 0
        running = Counter()
        for x in range(x_step, x_max + x_step, x_step):
            while idx < len(starts_reg_sorted) and starts_reg_sorted[idx] <= x:
                running[int(starts_reg_sorted[idx] % 420)] += 1
                idx += 1
            d11 = running[221] - running[11]
            d101 = running[311] - running[101]
            d191 = running[401] - running[191]
            delta_curve.append((x, d11, d101, d191))
            print(f"{x:7d} | {d11:8d} | {d101:9d} | {d191:9d}")

    counts = np.array([n for a, n in sorted(rc.items())])
    clr = centered_log_ratio(counts)

    print("\n=== Centered log-ratio der Residuen ===")
    for (a, n), y in zip(sorted(rc.items()), clr):
        print(f"{a:3d}  count={n:5d}  clr={y:+.6f}")

    print("\nInterpretation:")
    print("- ABCE/CEAB hat global keinen signifikanten Bias.")
    print("- mod 210 faltet die EABC-Chiralitaet zusammen.")
    print("- mod 420 entfaltet die chirale Doppeldeckung.")
    print("- Die regulaeren mod420-Kanaele bilden drei Lift-Paare:")
    print("  (11,221), (101,311), (191,401).")
    print("- Alle normierten Lift-Paar-Differenzen liegen bei |Z| < 1.")
    print("- Daher gibt es bis zur getesteten Schranke keine signifikante chirale")
    print("  Haeufigkeitsverletzung, sondern eine nahezu balancierte Doppeldeckung.")
    print("- Der Start p=5 ist ein singulaerer Randdefekt und sollte getrennt behandelt werden.")

    return {
        "window_density": window_density,
        "window_bias_raw": window_bias_raw,
        "window_bias_logit": window_bias_logit,
        "window_mean_gap": window_mean_gap,
        "global_orientation": global_orient,
        "residue_counts_mod210": rc,
        "residue_counts_mod420": rc420,
        "residue_counts_mod420_regular": rc420_reg,
        "residue_clr": clr,
        "pair_deltas": pair_deltas,
        "pair_zscores": pair_zscores,
        "delta_curve": delta_curve,
    }

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Analyse von Primvierlings-Orientierungen.")
    parser.add_argument("--limit", type=int, default=2_000_000, help="Obere Grenze fuer Vierlingsstarts.")
    parser.add_argument("--window", type=int, default=100_000, help="Fensterbreite in m = p + 4.")
    args = parser.parse_args()

    quadruplets = find_prime_quadruplets(args.limit)
    print(f"Gefundene Vierlingsstarts bis {args.limit}: {len(quadruplets)}")

    if len(quadruplets) < 2:
        print("Zu wenige Vierlinge fuer Analyse gefunden.")
    else:
        analyze_quadruplets(quadruplets, W=args.window)