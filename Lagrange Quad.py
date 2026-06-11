import math
import statistics
from itertools import permutations


def four_square_reps(n, positive_only=False, ordered=False):
    """
    Erzeuge Vierquadrate-Darstellungen von n.

    Rückgabe:
        - ordered=False: ungeordnete Tupel (E, A, B, C) mit E >= A >= B >= C >= 0
        - ordered=True : alle verschiedenen Permutationen solcher Tupel

    positive_only=True erzwingt C >= 1 statt C >= 0.
    """
    reps = []
    min_c = 1 if positive_only else 0

    limit = int(math.isqrt(n))
    for E in range(limit, -1, -1):
        E2 = E * E
        if E2 > n:
            continue
        for A in range(min(E, int(math.isqrt(n - E2))), -1, -1):
            A2 = A * A
            if E2 + A2 > n:
                continue
            for B in range(min(A, int(math.isqrt(n - E2 - A2))), -1, -1):
                B2 = B * B
                if E2 + A2 + B2 > n:
                    continue
                rest = n - E2 - A2 - B2
                C = int(math.isqrt(rest))
                if C < min_c or C > B:
                    continue
                if C * C == rest:
                    reps.append((E, A, B, C))

    if not ordered:
        return reps

    all_ordered = set()
    for rep in reps:
        for p in set(permutations(rep)):
            all_ordered.add(p)
    return sorted(all_ordered, reverse=True)


def square_lengths(rep, use_squares=True):
    """
    Liefert die vier Längen für das Kreis-/Sehnenmodell.
    Standard: l_i = E^2, A^2, B^2, C^2
    Optional weichere Variante: l_i = E, A, B, C
    """
    if use_squares:
        return [x * x for x in rep]
    return list(rep)


def normalized_weights(rep, use_squares=True):
    lengths = square_lengths(rep, use_squares=use_squares)
    total = sum(lengths)
    if total == 0:
        return [0.0, 0.0, 0.0, 0.0]
    return [x / total for x in lengths]


def entropy(weights):
    s = 0.0
    for p in weights:
        if p > 0:
            s -= p * math.log(p)
    return s


def concentration(weights):
    return sum(p * p for p in weights)


def balance_index(rep, use_squares=True):
    lengths = square_lengths(rep, use_squares=use_squares)
    mx = max(lengths)
    mn = min(lengths)
    if mx == 0:
        return 1.0
    return mn / mx


def variance_from_uniform(weights):
    return sum((p - 0.25) ** 2 for p in weights) / 4.0


def closure_function(r, lengths):
    """
    f(r) = sum arcsin(l_i / (2r))
    Gesucht ist f(r) = pi
    """
    s = 0.0
    for l in lengths:
        x = l / (2.0 * r)
        # numerische Absicherung
        if x < 0:
            return float("nan")
        if x > 1:
            return float("nan")
        s += math.asin(x)
    return s


def has_chord_closure(rep, use_squares=True, tol=1e-12):
    """
    Existenztest für die Gleichung

        sum_i arcsin(l_i / (2r)) = pi

    mit l_i = E^2, A^2, B^2, C^2 (oder alternativ E,A,B,C)

    Da f(r) streng fallend ist, existiert genau dann eine Lösung,
    wenn f(r_min) >= pi (numerisch: >= pi - tol), wobei r_min = max(l_i)/2.
    """
    lengths = square_lengths(rep, use_squares=use_squares)
    mx = max(lengths)
    if mx == 0:
        return False

    r_min = mx / 2.0
    f_min = closure_function(r_min, lengths)
    return f_min >= math.pi - tol


def chord_closure_radius(rep, use_squares=True, tol=1e-12, max_iter=200):
    """
    Löse per Bisektion:

        sum_i arcsin(l_i / (2r)) = pi

    Rückgabe:
        Radius r oder None, falls keine Lösung existiert.
    """
    lengths = square_lengths(rep, use_squares=use_squares)
    mx = max(lengths)
    if mx == 0:
        return None

    r_lo = mx / 2.0
    f_lo = closure_function(r_lo, lengths)

    if math.isnan(f_lo) or f_lo < math.pi - tol:
        return None

    if abs(f_lo - math.pi) <= tol:
        return r_lo

    r_hi = max(sum(lengths), mx, 1.0)
    f_hi = closure_function(r_hi, lengths)

    # obere Grenze vergrößern, bis f_hi < pi
    grow_count = 0
    while (math.isnan(f_hi) or f_hi >= math.pi) and grow_count < 100:
        r_hi *= 2.0
        f_hi = closure_function(r_hi, lengths)
        grow_count += 1

    if math.isnan(f_hi) or f_hi >= math.pi:
        return None

    # Bisektion
    for _ in range(max_iter):
        r_mid = 0.5 * (r_lo + r_hi)
        f_mid = closure_function(r_mid, lengths)

        if math.isnan(f_mid):
            return None

        if abs(f_mid - math.pi) <= tol:
            return r_mid

        if f_mid > math.pi:
            r_lo = r_mid
        else:
            r_hi = r_mid

    return 0.5 * (r_lo + r_hi)


def analyze_rep(rep, use_squares=True):
    """
    Analysiere ein einzelnes Tupel.
    """
    lengths = square_lengths(rep, use_squares=use_squares)
    weights = normalized_weights(rep, use_squares=use_squares)

    fit = has_chord_closure(rep, use_squares=use_squares)
    kappa = chord_closure_radius(rep, use_squares=use_squares) if fit else None

    return {
        "rep": rep,
        "lengths": lengths,
        "weights": weights,
        "entropy": entropy(weights),
        "concentration": concentration(weights),
        "balance": balance_index(rep, use_squares=use_squares),
        "variance_uniform": variance_from_uniform(weights),
        "fit": fit,
        "kappa": kappa,
    }


def analyze_n(n, positive_only=False, ordered=False, use_squares=True):
    """
    Analysiere eine Zahl n.
    """
    reps = four_square_reps(n, positive_only=positive_only, ordered=ordered)
    rep_data = [analyze_rep(rep, use_squares=use_squares) for rep in reps]

    fit_data = [d for d in rep_data if d["fit"] and d["kappa"] is not None]
    kappas = [d["kappa"] for d in fit_data]
    entropies = [d["entropy"] for d in rep_data]

    result = {
        "n": n,
        "num_reps": len(reps),
        "num_fit": len(fit_data),
        "fit_ratio": (len(fit_data) / len(reps)) if reps else 0.0,
        "kappa_min": min(kappas) if kappas else None,
        "kappa_max": max(kappas) if kappas else None,
        "kappa_mean": statistics.mean(kappas) if kappas else None,
        "entropy_min": min(entropies) if entropies else None,
        "entropy_max": max(entropies) if entropies else None,
        "entropy_mean": statistics.mean(entropies) if entropies else None,
        "representations": rep_data,
    }
    return result


def analyze_range(N, positive_only=False, ordered=False, use_squares=True):
    """
    Analysiere alle n = 1,...,N
    """
    return [
        analyze_n(
            n,
            positive_only=positive_only,
            ordered=ordered,
            use_squares=use_squares,
        )
        for n in range(1, N + 1)
    ]


def print_summary_table(results, limit=None):
    """
    Kompakte Tabelle.
    """
    print(
        f"{'n':>4} {'#Q':>4} {'#fit':>5} {'rho':>8} "
        f"{'H_min':>10} {'H_max':>10} {'k_min':>12} {'k_max':>12}"
    )
    print("-" * 72)

    rows = results if limit is None else results[:limit]
    for r in rows:
        def fmt(x, digits=6):
            return "-" if x is None else f"{x:.{digits}f}"

        print(
            f"{r['n']:>4} "
            f"{r['num_reps']:>4} "
            f"{r['num_fit']:>5} "
            f"{r['fit_ratio']:>8.4f} "
            f"{fmt(r['entropy_min']):>10} "
            f"{fmt(r['entropy_max']):>10} "
            f"{fmt(r['kappa_min']):>12} "
            f"{fmt(r['kappa_max']):>12}"
        )


def print_representations_for_n(result, max_rows=20):
    """
    Detaillierte Ausgabe für eine einzelne Zahl n.
    """
    print(f"\n=== n = {result['n']} ===")
    print(f"Anzahl Darstellungen: {result['num_reps']}")
    print(f"Einrastfähig:         {result['num_fit']}")
    print(f"Einrastquote rho(n):  {result['fit_ratio']:.4f}")
    print()

    rows = result["representations"][:max_rows]
    for d in rows:
        rep = d["rep"]
        print(
            f"rep={rep} | "
            f"L={d['lengths']} | "
            f"H={d['entropy']:.6f} | "
            f"K={d['concentration']:.6f} | "
            f"bal={d['balance']:.6f} | "
            f"fit={d['fit']} | "
            f"kappa={'-' if d['kappa'] is None else f'{d['kappa']:.12f}'}"
        )


def first_n_primes(n):
    primes = []
    x = 2
    while len(primes) < n:
        is_prime = True
        r = int(math.isqrt(x))
        for p in primes:
            if p > r:
                break
            if x % p == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(x)
        x += 1 if x == 2 else 2
    return primes


def analyze_prime_list(primes, positive_only=False, ordered=False, use_squares=True):
    results = [
        analyze_n(
            p,
            positive_only=positive_only,
            ordered=ordered,
            use_squares=use_squares,
        )
        for p in primes
    ]
    return results


def print_prime_summary(results, max_rows=30):
    print(
        f"{'p':>6} {'#Q':>5} {'#fit':>6} {'rho':>8} "
        f"{'k_min':>12} {'k_max':>12}"
    )
    print("-" * 60)

    for r in results[:max_rows]:
        def fmt(x):
            return "-" if x is None else f"{x:.6f}"

        print(
            f"{r['n']:>6} "
            f"{r['num_reps']:>5} "
            f"{r['num_fit']:>6} "
            f"{r['fit_ratio']:>8.4f} "
            f"{fmt(r['kappa_min']):>12} "
            f"{fmt(r['kappa_max']):>12}"
        )


def aggregate_prime_results(results):
    total = len(results)
    with_fit = sum(1 for r in results if r["num_fit"] > 0)
    without_fit = [r["n"] for r in results if r["num_fit"] == 0]

    avg_reps = statistics.mean(r["num_reps"] for r in results)
    max_reps = max(r["num_reps"] for r in results)
    argmax_reps = [r["n"] for r in results if r["num_reps"] == max_reps]

    print("\n=== Gesamtübersicht erste Primzahlen ===")
    print(f"Anzahl Primzahlen:                {total}")
    print(f"Mit mindestens 1 fit-Darstellung: {with_fit}")
    print(f"Ohne fit-Darstellung:             {len(without_fit)}")
    print(f"Mittlere Anzahl Darstellungen:    {avg_reps:.3f}")
    print(f"Maximale Anzahl Darstellungen:    {max_reps}")
    print(f"Erreicht bei Primzahl(en):        {argmax_reps[:10]}")
    print()
    print("Primzahlen ohne fit-Darstellung:")
    print(without_fit)


def critical_closure_value(rep, use_squares=True):
    """
    f(r_min) am kritischen Radius r_min = max(l_i)/2.
    """
    lengths = square_lengths(rep, use_squares=use_squares)
    mx = max(lengths)
    if mx == 0:
        return None
    r_min = mx / 2.0
    return closure_function(r_min, lengths)


def closure_gap_at_critical(rep, use_squares=True):
    """gap = pi - f(r_min)."""
    val = critical_closure_value(rep, use_squares=use_squares)
    if val is None:
        return None
    return math.pi - val


def analyze_nonfit_prime_details(prime_results, use_squares=True):
    """
    Extrahiert alle Primzahlen ohne fit-Darstellung und analysiert
    ihre Repräsentationen im Detail.
    """
    nonfit_prime_results = [r for r in prime_results if r["num_fit"] == 0]
    detailed = []

    for r in nonfit_prime_results:
        reps = r["representations"]

        enriched = []
        for d in reps:
            crit = critical_closure_value(d["rep"], use_squares=use_squares)
            gap = closure_gap_at_critical(d["rep"], use_squares=use_squares)

            enriched.append({
                **d,
                "critical_value": crit,
                "critical_gap": gap,
            })

        # Sortierung: zuerst kleinste gap, also "am ehesten noch fit"
        enriched.sort(
            key=lambda x: (
                float("inf") if x["critical_gap"] is None else x["critical_gap"],
                x["rep"]
            )
        )

        detailed.append({
            "p": r["n"],
            "num_reps": r["num_reps"],
            "representations": enriched,
            "best_gap": enriched[0]["critical_gap"] if enriched else None,
            "best_rep": enriched[0]["rep"] if enriched else None,
        })

    return detailed


def print_nonfit_prime_summary(nonfit_details):
    """
    Kompakte Übersicht aller nicht-einrastenden Primzahlen.
    """
    print("\n=== Übersicht der nicht-einrastenden Primzahlen ===")
    print(f"{'p':>6} {'#Q':>5} {'beste rep':>18} {'best_gap':>14}")
    print("-" * 56)

    for item in nonfit_details:
        p = item["p"]
        num_reps = item["num_reps"]
        best_rep = item["best_rep"]
        best_gap = item["best_gap"]

        rep_str = "-" if best_rep is None else str(best_rep)
        gap_str = "-" if best_gap is None else f"{best_gap:.12f}"

        print(f"{p:>6} {num_reps:>5} {rep_str:>18} {gap_str:>14}")


def print_best_nonfit_candidates(nonfit_details, top_k=20):
    """
    Zeigt die non-fit-Primzahlen, sortiert nach kleinster gap.
    Das sind die 'fast einrastenden' Kandidaten.
    """
    items = sorted(
        nonfit_details,
        key=lambda x: float("inf") if x["best_gap"] is None else x["best_gap"]
    )

    print(f"\n=== Top {top_k} fast einrastende non-fit-Primzahlen ===")
    print(f"{'p':>6} {'beste rep':>18} {'best_gap':>14}")
    print("-" * 44)

    for item in items[:top_k]:
        p = item["p"]
        rep = item["best_rep"]
        gap = item["best_gap"]

        rep_str = "-" if rep is None else str(rep)
        gap_str = "-" if gap is None else f"{gap:.12f}"

        print(f"{p:>6} {rep_str:>18} {gap_str:>14}")


def print_nonfit_prime_details(nonfit_details, max_primes=None, max_reps_per_prime=None):
    """
    Detaillierte Ausgabe für nicht-einrastende Primzahlen.

    Für jede Darstellung werden ausgegeben:
    - rep
    - Längen
    - Entropie
    - Konzentration
    - Balance
    - kritischer Wert f(r_min)
    - gap = pi - f(r_min)
    """
    print("\n=== Detailanalyse der nicht-einrastenden Primzahlen ===")

    items = nonfit_details if max_primes is None else nonfit_details[:max_primes]

    for item in items:
        print(f"\n--- Primzahl p = {item['p']} | #Q = {item['num_reps']} ---")
        reps = item["representations"]
        reps = reps if max_reps_per_prime is None else reps[:max_reps_per_prime]

        for d in reps:
            rep = d["rep"]
            lengths = d["lengths"]
            H = d["entropy"]
            K = d["concentration"]
            bal = d["balance"]
            crit = d["critical_value"]
            gap = d["critical_gap"]

            crit_str = "-" if crit is None else f"{crit:.12f}"
            gap_str = "-" if gap is None else f"{gap:.12f}"

            print(
                f"rep={rep} | "
                f"L={lengths} | "
                f"H={H:.6f} | "
                f"K={K:.6f} | "
                f"bal={bal:.6f} | "
                f"f(r_min)={crit_str} | "
                f"pi-f(r_min)={gap_str}"
            )


def export_nonfit_prime_details_to_csv(nonfit_details, filename="nonfit_primes_details.csv"):
    """
    Exportiert alle Details der nicht-einrastenden Primzahlen in eine CSV-Datei.
    """
    import csv

    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "p",
            "num_reps_for_p",
            "rep_E",
            "rep_A",
            "rep_B",
            "rep_C",
            "lengths",
            "entropy",
            "concentration",
            "balance",
            "critical_value",
            "critical_gap",
        ])

        for item in nonfit_details:
            p = item["p"]
            num_reps = item["num_reps"]

            for d in item["representations"]:
                E, A, B, C = d["rep"]
                writer.writerow([
                    p,
                    num_reps,
                    E, A, B, C,
                    str(d["lengths"]),
                    d["entropy"],
                    d["concentration"],
                    d["balance"],
                    d["critical_value"],
                    d["critical_gap"],
                ])

    print(f"\nCSV exportiert: {filename}")


def mean_safe(values):
    return statistics.mean(values) if values else None


def stdev_safe(values):
    if len(values) >= 2:
        return statistics.pstdev(values)
    if len(values) == 1:
        return 0.0
    return None


def min_safe(values):
    return min(values) if values else None


def max_safe(values):
    return max(values) if values else None


def flatten_prime_representation_data(prime_results, use_squares=True):
    flat = []
    for r in prime_results:
        p = r["n"]
        for d in r["representations"]:
            rep = d["rep"]
            crit = critical_closure_value(rep, use_squares=use_squares)
            gap = closure_gap_at_critical(rep, use_squares=use_squares)
            flat.append({
                "p": p,
                "rep": rep,
                "lengths": d["lengths"],
                "entropy": d["entropy"],
                "concentration": d["concentration"],
                "balance": d["balance"],
                "variance_uniform": d["variance_uniform"],
                "fit": d["fit"],
                "kappa": d["kappa"],
                "critical_value": crit,
                "critical_gap": gap,
            })
    return flat


def describe_metric(rows, key):
    vals = [x[key] for x in rows if x[key] is not None]
    return {
        "n": len(vals),
        "mean": mean_safe(vals),
        "min": min_safe(vals),
        "max": max_safe(vals),
        "stdev": stdev_safe(vals),
    }


def compare_fit_vs_nonfit(prime_results, use_squares=True):
    flat = flatten_prime_representation_data(prime_results, use_squares=use_squares)
    fit_rows = [x for x in flat if x["fit"]]
    nonfit_rows = [x for x in flat if not x["fit"]]

    return {
        "prime_count": len(prime_results),
        "use_squares": use_squares,
        "fit_rows": fit_rows,
        "nonfit_rows": nonfit_rows,
        "fit_metrics": {
            "entropy": describe_metric(fit_rows, "entropy"),
            "concentration": describe_metric(fit_rows, "concentration"),
            "balance": describe_metric(fit_rows, "balance"),
            "kappa": describe_metric(fit_rows, "kappa"),
        },
        "nonfit_metrics": {
            "entropy": describe_metric(nonfit_rows, "entropy"),
            "concentration": describe_metric(nonfit_rows, "concentration"),
            "balance": describe_metric(nonfit_rows, "balance"),
            "critical_value": describe_metric(nonfit_rows, "critical_value"),
            "critical_gap": describe_metric(nonfit_rows, "critical_gap"),
        },
    }


def print_fit_vs_nonfit_report(stats_result):
    def fmt_metric(label, d):
        if d["n"] == 0:
            print(f"  {label}: (keine Daten)")
            return
        print(
            f"  {label}: n={d['n']}  mean={d['mean']:.6f}  min={d['min']:.6f}  "
            f"max={d['max']:.6f}  stdev={d['stdev']:.6f}"
        )

    print("\n=== Statistik fit vs. non-fit (alle Darstellungen) ===")
    print(f"Primzahlen in Stichprobe: {stats_result['prime_count']}")
    print(f"use_squares: {stats_result['use_squares']}")

    fit_rows = stats_result["fit_rows"]
    nonfit_rows = stats_result["nonfit_rows"]

    print("\n--- fit=True ---")
    print(f"  Anzahl Darstellungen: {len(fit_rows)}")
    for key, label in [
        ("entropy", "entropy"),
        ("concentration", "concentration"),
        ("balance", "balance"),
        ("kappa", "kappa"),
    ]:
        fmt_metric(label, stats_result["fit_metrics"][key])

    print("\n--- fit=False ---")
    print(f"  Anzahl Darstellungen: {len(nonfit_rows)}")
    for key, label in [
        ("entropy", "entropy"),
        ("concentration", "concentration"),
        ("balance", "balance"),
        ("critical_value", "critical_value"),
        ("critical_gap", "critical_gap"),
    ]:
        fmt_metric(label, stats_result["nonfit_metrics"][key])


def export_fit_vs_nonfit_summary_to_csv(stats_result, filename="fit_vs_nonfit_summary.csv"):
    import csv

    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["group", "metric", "n", "mean", "min", "max", "stdev"])

        for metric, d in stats_result["fit_metrics"].items():
            writer.writerow(["fit", metric, d["n"], d["mean"], d["min"], d["max"], d["stdev"]])

        for metric, d in stats_result["nonfit_metrics"].items():
            writer.writerow(["nonfit", metric, d["n"], d["mean"], d["min"], d["max"], d["stdev"]])

    print(f"\nCSV exportiert: {filename}")


def select_best_representation_for_prime(result, use_squares=True):
    reps = result["representations"]
    fit_rows = [d for d in reps if d["fit"] and d["kappa"] is not None]

    if fit_rows:
        best = sorted(
            fit_rows,
            key=lambda d: (
                d["kappa"],
                -d["entropy"],
                -d["balance"],
                d["rep"],
            )
        )[0]

        return {
            "p": result["n"],
            "status": "fit",
            "rep": best["rep"],
            "lengths": best["lengths"],
            "entropy": best["entropy"],
            "concentration": best["concentration"],
            "balance": best["balance"],
            "variance_uniform": best["variance_uniform"],
            "kappa": best["kappa"],
            "critical_gap": 0.0,
        }

    enriched = []
    for d in reps:
        gap = closure_gap_at_critical(d["rep"], use_squares=use_squares)
        enriched.append((d, gap))

    best_d, best_gap = sorted(
        enriched,
        key=lambda item: (
            float("inf") if item[1] is None else item[1],
            -item[0]["entropy"],
            -item[0]["balance"],
            item[0]["rep"],
        )
    )[0]

    return {
        "p": result["n"],
        "status": "nonfit",
        "rep": best_d["rep"],
        "lengths": best_d["lengths"],
        "entropy": best_d["entropy"],
        "concentration": best_d["concentration"],
        "balance": best_d["balance"],
        "variance_uniform": best_d["variance_uniform"],
        "kappa": None,
        "critical_gap": best_gap,
    }


def analyze_best_representations_for_primes(primes, positive_only=False, ordered=False, use_squares=True):
    results = analyze_prime_list(
        primes,
        positive_only=positive_only,
        ordered=ordered,
        use_squares=use_squares,
    )

    best_rows = [
        select_best_representation_for_prime(r, use_squares=use_squares)
        for r in results
    ]

    return {
        "mode": "hard" if use_squares else "soft",
        "prime_results": results,
        "best_rows": best_rows,
    }


def compare_best_hard_vs_soft(primes, positive_only=False, ordered=False):
    hard = analyze_best_representations_for_primes(
        primes,
        positive_only=positive_only,
        ordered=ordered,
        use_squares=True,
    )
    soft = analyze_best_representations_for_primes(
        primes,
        positive_only=positive_only,
        ordered=ordered,
        use_squares=False,
    )

    hard_map = {row["p"]: row for row in hard["best_rows"]}
    soft_map = {row["p"]: row for row in soft["best_rows"]}

    combined = []
    for p in primes:
        h = hard_map[p]
        s = soft_map[p]

        combined.append({
            "p": p,
            "hard_status": h["status"],
            "hard_rep": h["rep"],
            "hard_entropy": h["entropy"],
            "hard_concentration": h["concentration"],
            "hard_balance": h["balance"],
            "hard_kappa": h["kappa"],
            "hard_gap": h["critical_gap"],
            "soft_status": s["status"],
            "soft_rep": s["rep"],
            "soft_entropy": s["entropy"],
            "soft_concentration": s["concentration"],
            "soft_balance": s["balance"],
            "soft_kappa": s["kappa"],
            "soft_gap": s["critical_gap"],
        })

    return {
        "hard": hard,
        "soft": soft,
        "combined": combined,
    }


def print_best_representation_table(best_analysis, max_rows=40):
    rows = best_analysis["best_rows"][:max_rows]
    mode = best_analysis["mode"]

    print(f"\n=== Beste Darstellung pro Primzahl | Modus: {mode} ===")
    print(
        f"{'p':>6} {'status':>8} {'rep':>18} {'H':>10} {'K':>10} {'bal':>10} {'kappa':>12} {'gap':>12}"
    )
    print("-" * 96)

    def fmt(x, digits=6):
        return "-" if x is None else f"{x:.{digits}f}"

    for r in rows:
        print(
            f"{r['p']:>6} "
            f"{r['status']:>8} "
            f"{str(r['rep']):>18} "
            f"{fmt(r['entropy']):>10} "
            f"{fmt(r['concentration']):>10} "
            f"{fmt(r['balance']):>10} "
            f"{fmt(r['kappa']):>12} "
            f"{fmt(r['critical_gap']):>12}"
        )


def print_hard_vs_soft_comparison(compare_result, max_rows=60, only_differences=False):
    rows = compare_result["combined"]

    if only_differences:
        rows = [
            r for r in rows
            if (
                r["hard_status"] != r["soft_status"]
                or r["hard_rep"] != r["soft_rep"]
            )
        ]

    rows = rows[:max_rows]

    print("\n=== Harte vs. weiche Version: beste Darstellung pro Primzahl ===")
    print(
        f"{'p':>6} "
        f"{'hard':>8} {'hard_rep':>18} {'hard_gap':>12} "
        f"{'soft':>8} {'soft_rep':>18} {'soft_gap':>12}"
    )
    print("-" * 92)

    def fmt(x, digits=6):
        return "-" if x is None else f"{x:.{digits}f}"

    for r in rows:
        print(
            f"{r['p']:>6} "
            f"{r['hard_status']:>8} {str(r['hard_rep']):>18} {fmt(r['hard_gap']):>12} "
            f"{r['soft_status']:>8} {str(r['soft_rep']):>18} {fmt(r['soft_gap']):>12}"
        )


def summarize_hard_vs_soft(compare_result):
    rows = compare_result["combined"]

    hard_fit = sum(1 for r in rows if r["hard_status"] == "fit")
    soft_fit = sum(1 for r in rows if r["soft_status"] == "fit")

    improved = sum(
        1 for r in rows
        if r["hard_status"] == "nonfit" and r["soft_status"] == "fit"
    )

    worsened = sum(
        1 for r in rows
        if r["hard_status"] == "fit" and r["soft_status"] == "nonfit"
    )

    both_nonfit = [
        r for r in rows
        if r["hard_status"] == "nonfit" and r["soft_status"] == "nonfit"
    ]

    print("\n=== Zusammenfassung hart vs. weich ===")
    print(f"Primzahlen gesamt:                  {len(rows)}")
    print(f"Fit im harten Modell:              {hard_fit}")
    print(f"Fit im weichen Modell:             {soft_fit}")
    print(f"Verbesserung hard->soft:           {improved}")
    print(f"Verschlechterung hard->soft:       {worsened}")
    print(f"In beiden Modellen non-fit:        {len(both_nonfit)}")

    if both_nonfit:
        print("Beispiele (erste 20):")
        print([r["p"] for r in both_nonfit[:20]])


def export_best_representation_comparison_to_csv(compare_result, filename="best_representation_hard_vs_soft.csv"):
    import csv

    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "p",
            "hard_status", "hard_rep", "hard_entropy", "hard_concentration", "hard_balance", "hard_kappa", "hard_gap",
            "soft_status", "soft_rep", "soft_entropy", "soft_concentration", "soft_balance", "soft_kappa", "soft_gap",
        ])

        for r in compare_result["combined"]:
            writer.writerow([
                r["p"],
                r["hard_status"], str(r["hard_rep"]), r["hard_entropy"], r["hard_concentration"], r["hard_balance"], r["hard_kappa"], r["hard_gap"],
                r["soft_status"], str(r["soft_rep"]), r["soft_entropy"], r["soft_concentration"], r["soft_balance"], r["soft_kappa"], r["soft_gap"],
            ])

    print(f"\nCSV exportiert: {filename}")


if __name__ == "__main__":
    # Basislauf 1..50
    results = analyze_range(
        N=50,
        positive_only=False,
        ordered=False,
        use_squares=True,
    )
    print_summary_table(results, limit=50)
    print_representations_for_n(results[24])  # n = 25
    print_representations_for_n(results[48])  # n = 49

    # Erste 1000 Primzahlen
    primes_1000 = first_n_primes(1000)

    prime_results = analyze_prime_list(
        primes_1000,
        positive_only=False,
        ordered=False,
        use_squares=True,
    )

    print_prime_summary(prime_results, max_rows=40)
    aggregate_prime_results(prime_results)

    # Statistik fit vs. non-fit
    fit_nonfit_stats = compare_fit_vs_nonfit(
        prime_results,
        use_squares=True,
    )
    print_fit_vs_nonfit_report(fit_nonfit_stats)
    export_fit_vs_nonfit_summary_to_csv(
        fit_nonfit_stats,
        filename="fit_vs_nonfit_summary.csv",
    )

    # Non-fit-Details
    nonfit_details = analyze_nonfit_prime_details(
        prime_results,
        use_squares=True,
    )
    print_nonfit_prime_summary(nonfit_details)
    print_best_nonfit_candidates(nonfit_details, top_k=20)
    print_nonfit_prime_details(
        nonfit_details,
        max_primes=None,
        max_reps_per_prime=None,
    )
    export_nonfit_prime_details_to_csv(
        nonfit_details,
        filename="nonfit_primes_details.csv",
    )

    # Beste Darstellung pro Primzahl / hart vs. weich
    best_hard = analyze_best_representations_for_primes(
        primes_1000,
        positive_only=False,
        ordered=False,
        use_squares=True,
    )
    best_soft = analyze_best_representations_for_primes(
        primes_1000,
        positive_only=False,
        ordered=False,
        use_squares=False,
    )
    print_best_representation_table(best_hard, max_rows=40)
    print_best_representation_table(best_soft, max_rows=40)

    hard_soft_compare = compare_best_hard_vs_soft(
        primes_1000,
        positive_only=False,
        ordered=False,
    )
    summarize_hard_vs_soft(hard_soft_compare)
    print_hard_vs_soft_comparison(
        hard_soft_compare,
        max_rows=80,
        only_differences=True,
    )
    export_best_representation_comparison_to_csv(
        hard_soft_compare,
        filename="best_representation_hard_vs_soft.csv",
    )