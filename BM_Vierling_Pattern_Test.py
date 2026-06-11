#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BM_Vierling_Pattern_Test.py

Statistische Analyse einer Liste von Primzahlvierlingen:
    (p, p+2, p+6, p+8)

Erwartete CSV-Spalten u.a.:
    Index,
    E_Restklasse_1,
    A_Restklasse_5,
    B_Restklasse_7,
    C_Restklasse_11,
    Mittel,
    Norm_Mittel_unter_1,
    Dreiteilungen_wie_oft_v3,
    Potenz_3_teilt_Mittel

Die Datei darf auch die leeren Doppelkomma-Spalten aus deinem Export enthalten.

Aufruf:
    python BM_Vierling_Pattern_Test.py [pfad/zur.csv]

Ohne Argument wird ``primzahlvierlinge_1000.csv`` im Skript- oder Arbeitsverzeichnis
gesucht.
"""

import sys
import math
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from collections import Counter

try:
    from scipy.stats import chisquare, fisher_exact
    SCIPY_OK = True
except Exception:
    SCIPY_OK = False


# ------------------------------------------------------------
# Hilfsfunktionen
# ------------------------------------------------------------

def v_p(n: int, p: int) -> int:
    """p-adische Bewertung v_p(n)."""
    c = 0
    while n % p == 0 and n != 0:
        n //= p
        c += 1
    return c


def primes_upto(n: int):
    """Kleine Primzahlen bis n."""
    if n < 2:
        return []
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for k in range(2, int(n ** 0.5) + 1):
        if sieve[k]:
            sieve[k*k:n+1:k] = False
    return [i for i, b in enumerate(sieve) if b]


def pair_count_within(sorted_x, W: int) -> int:
    """
    Anzahl Paare (i,j), i<j, mit x_j - x_i <= W.
    Zwei-Zeiger-Verfahren.
    """
    x = np.asarray(sorted_x)
    n = len(x)
    j = 0
    count = 0
    for i in range(n):
        while j < n and x[j] - x[i] <= W:
            j += 1
        count += max(0, j - i - 1)
    return int(count)


def nearest_gap_flags(sorted_x, W: int):
    """
    Markiert, ob ein Zentrum einen Nachbarn im Abstand <= W hat.
    """
    x = np.asarray(sorted_x)
    n = len(x)
    flags = np.zeros(n, dtype=bool)
    if n <= 1:
        return flags
    gaps = np.diff(x)
    flags[:-1] |= gaps <= W
    flags[1:] |= gaps <= W
    return flags


def summarize_gaps(m):
    gaps = np.diff(np.sort(m))
    return {
        "n_gaps": len(gaps),
        "min": int(np.min(gaps)),
        "q05": float(np.quantile(gaps, 0.05)),
        "q25": float(np.quantile(gaps, 0.25)),
        "median": float(np.median(gaps)),
        "mean": float(np.mean(gaps)),
        "q75": float(np.quantile(gaps, 0.75)),
        "q95": float(np.quantile(gaps, 0.95)),
        "max": int(np.max(gaps)),
    }


def empirical_pvalue(obs, sims, side="greater"):
    sims = np.asarray(sims)
    if side == "greater":
        return (np.sum(sims >= obs) + 1) / (len(sims) + 1)
    if side == "less":
        return (np.sum(sims <= obs) + 1) / (len(sims) + 1)
    if side == "two-sided":
        center = np.mean(sims)
        return (np.sum(np.abs(sims - center) >= abs(obs - center)) + 1) / (len(sims) + 1)
    raise ValueError("side must be greater, less, or two-sided")


# ------------------------------------------------------------
# CSV robust einlesen
# ------------------------------------------------------------

def load_vierlinge_csv(path: str) -> pd.DataFrame:
    """
    Liest auch die CSV mit doppelten Kommas.
    Pandas erzeugt dann Unnamed-Spalten, die wir ignorieren.
    """
    df = pd.read_csv(path)

    # Nur echte Spalten behalten, keine leeren Unnamed-Spalten.
    keep = [c for c in df.columns if not str(c).startswith("Unnamed")]
    df = df[keep].copy()

    # Spaltennamen säubern
    df.columns = [str(c).strip() for c in df.columns]

    needed = ["Index", "Mittel"]
    for c in needed:
        if c not in df.columns:
            raise ValueError(f"Pflichtspalte fehlt: {c}")

    # EABC-Spalten automatisch finden
    eabc_cols = [
        c for c in df.columns
        if "Restklasse" in c and any(tag in c for tag in ["E_", "A_", "B_", "C_"])
    ]

    if len(eabc_cols) != 4:
        raise ValueError(
            f"Ich habe nicht genau 4 EABC-Spalten gefunden, sondern {len(eabc_cols)}: {eabc_cols}"
        )

    # Numerisch machen
    for c in ["Index", "Mittel"] + eabc_cols:
        df[c] = pd.to_numeric(df[c], errors="raise").astype(np.int64)

    if "Dreiteilungen_wie_oft_v3" in df.columns:
        df["Dreiteilungen_wie_oft_v3"] = pd.to_numeric(
            df["Dreiteilungen_wie_oft_v3"], errors="coerce"
        ).astype("Int64")

    if "Potenz_3_teilt_Mittel" in df.columns:
        df["Potenz_3_teilt_Mittel"] = pd.to_numeric(
            df["Potenz_3_teilt_Mittel"], errors="coerce"
        ).astype("Int64")

    # Aus den vier EABC-Spalten die natürliche Vierlingsordnung rekonstruieren
    vals = df[eabc_cols].to_numpy()
    vals_sorted = np.sort(vals, axis=1)

    df["p_start"] = vals_sorted[:, 0]
    df["q1"] = vals_sorted[:, 0]
    df["q2"] = vals_sorted[:, 1]
    df["q3"] = vals_sorted[:, 2]
    df["q4"] = vals_sorted[:, 3]

    df["v3_calc"] = df["Mittel"].map(lambda x: v_p(int(x), 3))
    df["pow3_calc"] = 3 ** df["v3_calc"]

    df["start_mod12"] = df["p_start"] % 12
    df["center_mod12"] = df["Mittel"] % 12

    df["chirality"] = np.where(
        df["start_mod12"] == 5,
        "ABCE",
        np.where(df["start_mod12"] == 11, "CEAB", "OTHER")
    )

    return df


# ------------------------------------------------------------
# Admissibles Zentrumsgitter
# ------------------------------------------------------------

def admissible_centers(max_m: int, sieve_bound: int = 97):
    """
    Erzeugt Zentren m <= max_m, die auf kleinen Primzahlen nicht sofort verboten sind.

    Bedingung:
        m-4, m-2, m+2, m+4
    dürfen durch keine kleine Primzahl l <= sieve_bound teilbar sein.

    Zusätzlich:
        m % 3 == 0
        m ist positiv
    """
    m = np.arange(3, max_m + 1, 3, dtype=np.int64)

    ok = np.ones(len(m), dtype=bool)
    offsets = np.array([-4, -2, 2, 4], dtype=np.int64)

    for ell in primes_upto(sieve_bound):
        if ell in (2, 3):
            continue

        # Verboten, wenn irgendeine der vier Zahlen 0 mod ell ist.
        bad = np.zeros(len(m), dtype=bool)
        for off in offsets:
            bad |= ((m + off) % ell == 0)
        ok &= ~bad

    return m[ok]


def weighted_sample_without_replacement(rng, population, n, weights=None):
    """
    Ziehung ohne Zurücklegen.
    Für große Populationen und n=1000 noch ausreichend schnell.
    """
    population = np.asarray(population)
    if weights is None:
        idx = rng.choice(len(population), size=n, replace=False)
    else:
        w = np.asarray(weights, dtype=float)
        w = w / w.sum()
        idx = rng.choice(len(population), size=n, replace=False, p=w)
    return np.sort(population[idx])


# ------------------------------------------------------------
# Hauptanalyse
# ------------------------------------------------------------

def _default_vierlinge_csv() -> Path:
    """Project CSV next to this script, else ./primzahlvierlinge_1000.csv."""
    parent = Path(__file__).resolve().parent
    for base in (parent, Path.cwd()):
        cand = base / "primzahlvierlinge_1000.csv"
        if cand.is_file():
            return cand
    return parent / "primzahlvierlinge_1000.csv"


def main():
    parser = argparse.ArgumentParser(
        description="BM-Statistik für Primzahlvierlinge (CSV mit Mittel / EABC)."
    )
    parser.add_argument(
        "csv",
        nargs="?",
        default=None,
        help="Pfad zur Vierlings-CSV (Vorgabe: primzahlvierlinge_1000.csv im Skript- oder Arbeitsverzeichnis)",
    )
    args = parser.parse_args()
    path_obj = Path(args.csv) if args.csv else _default_vierlinge_csv()
    if not path_obj.is_file():
        print(
            "Keine CSV gefunden. Nutzung:\n"
            "  python BM_Vierling_Pattern_Test.py [vierlinge.csv]\n"
            f"Erwartet u.a.: {_default_vierlinge_csv()}"
        )
        sys.exit(1)
    path = str(path_obj)

    # Parameter
    RNG_SEED = 12345
    N_MONTE_CARLO = 2000
    SIEVE_BOUND = 97
    WINDOWS = [30, 60, 90, 120, 300, 1000, 3000, 5000, 10000]

    rng = np.random.default_rng(RNG_SEED)

    df = load_vierlinge_csv(path)
    df = df.sort_values("Mittel").reset_index(drop=True)
    n_all = len(df)
    df2 = df[df["Mittel"] > 1000].copy().reset_index(drop=True)
    df = df2

    m = df["Mittel"].to_numpy(dtype=np.int64)
    n = len(m)
    max_m = int(m.max())

    print("\n============================================================")
    print("BM-Primzahlvierlingsanalyse")
    print("============================================================")
    print(f"Anzahl Vierlinge: {n} (Mittel > 1000, zuvor {n_all} in CSV)")
    print(f"kleinstes Zentrum: {m.min()}")
    print(f"größtes Zentrum:  {m.max()}")
    print(f"Siebgrenze Baseline: Primzahlen <= {SIEVE_BOUND}")
    print(f"Monte-Carlo-Läufe: {N_MONTE_CARLO}")

    # --------------------------------------------------------
    # 1. Strukturchecks
    # --------------------------------------------------------

    print("\n------------------------------------------------------------")
    print("1. Strukturchecks")
    print("------------------------------------------------------------")

    is_quad = (
        (df["q2"] - df["q1"] == 2)
        & (df["q3"] - df["q1"] == 6)
        & (df["q4"] - df["q1"] == 8)
    )

    center_ok = df["Mittel"] == df["p_start"] + 4
    div3_ok = df["Mittel"] % 3 == 0
    mod12_ok = df["start_mod12"].isin([5, 11])

    print(f"Echte Form (p,p+2,p+6,p+8): {is_quad.sum()} / {n}")
    print(f"Zentrum m = p + 4:             {center_ok.sum()} / {n}")
    print(f"Zentrum durch 3 teilbar:       {div3_ok.sum()} / {n}")
    print(f"Startklasse 5 oder 11 mod 12:  {mod12_ok.sum()} / {n}")

    if not is_quad.all():
        print("\nWARNUNG: Einige Zeilen sind keine echten Vierlinge:")
        print(df.loc[~is_quad, ["Index", "q1", "q2", "q3", "q4", "Mittel"]].head(20))

    # --------------------------------------------------------
    # 2. Chiralität
    # --------------------------------------------------------

    print("\n------------------------------------------------------------")
    print("2. EABC-Chiralität")
    print("------------------------------------------------------------")

    chir_counts = df["chirality"].value_counts()
    print(chir_counts.to_string())

    if set(chir_counts.index).issubset({"ABCE", "CEAB"}):
        a = int(chir_counts.get("ABCE", 0))
        c = int(chir_counts.get("CEAB", 0))
        print(f"\nABCE : CEAB = {a} : {c}")
        print(f"Differenz = {a - c}")

        # einfacher Binomial-z-Score gegen 50:50
        z = (a - n/2) / math.sqrt(n/4)
        print(f"z gegen 50:50 ≈ {z:.3f}")

    # Runs gleicher Chiralität
    ch = df["chirality"].to_numpy()
    runs = []
    current = ch[0]
    length = 1
    for x in ch[1:]:
        if x == current:
            length += 1
        else:
            runs.append((current, length))
            current = x
            length = 1
    runs.append((current, length))

    run_lengths = np.array([r[1] for r in runs])
    print(f"\nRuns gleicher Chiralität: {len(runs)}")
    print(f"maximale Run-Länge:       {run_lengths.max()}")
    print(f"mittlere Run-Länge:       {run_lengths.mean():.3f}")

    # --------------------------------------------------------
    # 3. 3-adische Verteilung
    # --------------------------------------------------------

    print("\n------------------------------------------------------------")
    print("3. 3-adische Verteilung der Zentren")
    print("------------------------------------------------------------")

    v3_counts = Counter(df["v3_calc"])
    for k in sorted(v3_counts):
        print(f"v3={k}: {v3_counts[k]}")

    if "Dreiteilungen_wie_oft_v3" in df.columns:
        mismatch = df["Dreiteilungen_wie_oft_v3"].astype(float) != df["v3_calc"].astype(float)
        print(f"\nAbweichungen zur CSV-v3-Spalte: {int(mismatch.sum())}")

    if "Potenz_3_teilt_Mittel" in df.columns:
        mismatch_pow = df["Potenz_3_teilt_Mittel"].astype(float) != df["pow3_calc"].astype(float)
        print(f"Abweichungen zur CSV-Potenz-Spalte: {int(mismatch_pow.sum())}")

    # Chi-Quadrat-Test gegen geometrische Verteilung, Binning v3>=6
    # P(v3=k | 3|m) = 2 / 3^k für k>=1
    bins = [1, 2, 3, 4, 5]
    obs = []
    exp = []
    for k in bins:
        obs.append(v3_counts.get(k, 0))
        exp.append(n * (2 / (3 ** k)))

    obs_ge6 = sum(v for k, v in v3_counts.items() if k >= 6)
    exp_ge6 = n * sum(2 / (3 ** k) for k in range(6, 50))

    obs.append(obs_ge6)
    exp.append(exp_ge6)

    obs = np.array(obs, dtype=float)
    exp = np.array(exp, dtype=float)

    # Erwartungen auf gleiche Summe normieren
    exp *= obs.sum() / exp.sum()

    print("\nBinning für Chi-Quadrat: v3=1,2,3,4,5,>=6")
    print("observed:", obs.astype(int).tolist())
    print("expected:", [round(x, 2) for x in exp])

    if SCIPY_OK:
        chi, pval = chisquare(obs, f_exp=exp)
        print(f"Chi² = {chi:.3f}, p ≈ {pval:.4g}")
    else:
        print("scipy nicht verfügbar: Chi²-p-Wert übersprungen.")

    # --------------------------------------------------------
    # 4. Abstände
    # --------------------------------------------------------

    print("\n------------------------------------------------------------")
    print("4. Abstandsstatistik der Zentren")
    print("------------------------------------------------------------")

    gap_summary = summarize_gaps(m)
    for k, v in gap_summary.items():
        print(f"{k:>8}: {v}")

    gaps = np.diff(m)
    print("\nKleinste 20 Zentrum-Abstände:")
    smallest_idx = np.argsort(gaps)[:20]
    for idx in smallest_idx:
        print(
            f"Index sortiert {idx+1:4d}->{idx+2:4d}: "
            f"m={m[idx]} -> {m[idx+1]}, Δ={gaps[idx]}"
        )

    # --------------------------------------------------------
    # 5. Clusterzählungen
    # --------------------------------------------------------

    print("\n------------------------------------------------------------")
    print("5. Clusterzählungen im echten Datensatz")
    print("------------------------------------------------------------")

    obs_cluster = {}
    for W in WINDOWS:
        pc = pair_count_within(m, W)
        frac = nearest_gap_flags(m, W).mean()
        obs_cluster[W] = (pc, frac)
        print(f"W={W:5d}: Paare <=W = {pc:6d}, Anteil mit Nachbar <=W = {frac:.4f}")

    # --------------------------------------------------------
    # 6. Admissible Baseline
    # --------------------------------------------------------

    print("\n------------------------------------------------------------")
    print("6. Monte-Carlo gegen admissibles Siebgitter")
    print("------------------------------------------------------------")

    adm = admissible_centers(max_m=max_m, sieve_bound=SIEVE_BOUND)
    print(f"Admissible Zentren bis {max_m}: {len(adm)}")
    print(f"Admissible Dichte bezogen auf Vielfache von 3: {len(adm) / (max_m // 3):.6f}")

    # Inhomogene Gewichtung: Primzahlvierlingsdichte ~ 1 / log(x)^4
    # Für kleine m log stabilisieren.
    weights = 1.0 / np.maximum(np.log(adm), 2.0) ** 4
    weights = weights / weights.sum()

    sim_results = {W: [] for W in WINDOWS}
    sim_min_gap = []
    sim_median_gap = []
    sim_max_run = []

    for r in range(N_MONTE_CARLO):
        sample = weighted_sample_without_replacement(rng, adm, n, weights=weights)
        sgaps = np.diff(sample)

        sim_min_gap.append(np.min(sgaps))
        sim_median_gap.append(np.median(sgaps))

        # Clusterpaare
        for W in WINDOWS:
            sim_results[W].append(pair_count_within(sample, W))

        # Chiralität aus Zentrum mod 12:
        # m ≡ 9 => p ≡ 5 => ABCE
        # m ≡ 3 => p ≡ 11 => CEAB
        sim_ch = np.where(sample % 12 == 9, "ABCE", "CEAB")
        run_max = 1
        cur = 1
        for i in range(1, len(sim_ch)):
            if sim_ch[i] == sim_ch[i-1]:
                cur += 1
                run_max = max(run_max, cur)
            else:
                cur = 1
        sim_max_run.append(run_max)

    print("\nMonte-Carlo p-Werte:")
    print("Interpretation: kleiner p-Wert bei 'greater' heißt: echter Wert ist größer als Baseline.")

    for W in WINDOWS:
        obs_pc = obs_cluster[W][0]
        sims = np.array(sim_results[W])
        p_greater = empirical_pvalue(obs_pc, sims, side="greater")
        p_less = empirical_pvalue(obs_pc, sims, side="less")
        print(
            f"W={W:5d}: obs={obs_pc:6d}, "
            f"sim_mean={sims.mean():8.2f}, sim_sd={sims.std(ddof=1):8.2f}, "
            f"p_greater={p_greater:.4f}, p_less={p_less:.4f}"
        )

    obs_min_gap = int(np.min(gaps))
    obs_med_gap = float(np.median(gaps))
    p_min_less = empirical_pvalue(obs_min_gap, sim_min_gap, side="less")
    p_med_less = empirical_pvalue(obs_med_gap, sim_median_gap, side="less")

    print("\nGap-Vergleich:")
    print(
        f"min_gap: obs={obs_min_gap}, "
        f"sim_mean={np.mean(sim_min_gap):.2f}, "
        f"p_less={p_min_less:.4f}"
    )
    print(
        f"median_gap: obs={obs_med_gap:.2f}, "
        f"sim_mean={np.mean(sim_median_gap):.2f}, "
        f"p_less={p_med_less:.4f}"
    )

    obs_max_run = int(run_lengths.max())
    p_run_greater = empirical_pvalue(obs_max_run, sim_max_run, side="greater")

    print("\nChiralitäts-Runs gegen Baseline:")
    print(
        f"max_run: obs={obs_max_run}, "
        f"sim_mean={np.mean(sim_max_run):.2f}, "
        f"p_greater={p_run_greater:.4f}"
    )

    # --------------------------------------------------------
    # 7. Korrelation v3 vs lokale Cluster
    # --------------------------------------------------------

    print("\n------------------------------------------------------------")
    print("7. Korrelation: hohe 3-Adik vs lokale Cluster")
    print("------------------------------------------------------------")

    for W in [300, 1000, 3000, 10000]:
        clustered = nearest_gap_flags(m, W)
        high_v3 = df["v3_calc"].to_numpy() >= 3

        a = int(np.sum(high_v3 & clustered))
        b = int(np.sum(high_v3 & ~clustered))
        c = int(np.sum(~high_v3 & clustered))
        d = int(np.sum(~high_v3 & ~clustered))

        print(f"\nW={W}")
        print("Tabelle [[v3>=3 & Cluster, v3>=3 & kein Cluster],")
        print("         [v3<3  & Cluster, v3<3  & kein Cluster]]")
        print([[a, b], [c, d]])

        if SCIPY_OK:
            odds, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            print(f"Fisher exact: odds={odds:.4f}, p≈{p:.4g}")
        else:
            print("scipy nicht verfügbar: Fisher-Test übersprungen.")

    # --------------------------------------------------------
    # 8. Fazit-Hinweise
    # --------------------------------------------------------

    print("\n============================================================")
    print("Kurzinterpretation")
    print("============================================================")
    print("""
1. EABC-Restklassen, Zentrum m=p+4 und 3|m sind strukturell erzwungen.
   Diese Muster sind wichtig für das BM-Modell, aber nicht als Statistik-Signal zu zählen.

2. Die echte Statistik liegt in:
   - den Zentrum-Abständen Δ_i,
   - Clusterpaaren innerhalb W,
   - Runs der Chiralität ABCE/CEAB,
   - und der Frage, ob hohe v3(m) lokal häufiger in Clustern auftaucht.

3. Ein p-Wert < 0.01 in mehreren Fenstern W, stabil gegen verschiedene Siebgrenzen,
   wäre ein ernstzunehmender Hinweis. Einzelne kleine p-Werte sind wegen Mehrfachtests
   vorsichtig zu behandeln.

4. Besonders interessant für BM:
   Falls hohe 3-Adik, also v3(m)>=3 oder >=4, systematisch mit lokalen Clustern
   korreliert, wäre das ein nichttriviales 3-adisches Resonanzsignal.
""")


if __name__ == "__main__":
    main()