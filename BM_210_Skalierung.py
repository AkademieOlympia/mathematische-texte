#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BM_210_Skalierung.py

Skaliert die mod-210 3+1-Analyse der Primzahlvierlinge über mehrere N.

Getestet wird:
    Q(p) = (p, p+2, p+6, p+8)

Modulo 210:
    {5} ∪ {11,101,191}

Nach Entfernung des Ankers 5:
    drei Hauptfamilien 11,101,191

Für jedes N werden berechnet:
- Anzahl Vierlinge
- Familienhäufigkeiten
- Gleichverteilung der drei Familien
- Tripelneutralität
- Monochromie
- zyklische Orientierung
- Permutations-p-Werte

Ausgabe:
    bm210_scaling_results.csv
"""

import math
import time
import numpy as np
import pandas as pd
from collections import Counter

try:
    from scipy.stats import chisquare
    SCIPY_OK = True
except Exception:
    SCIPY_OK = False


# ------------------------------------------------------------
# Parameter
# ------------------------------------------------------------

N_VALUES = [
    1_000_000,
    2_000_000,
    5_000_000,
    10_000_000,
    20_000_000,
    50_000_000,
]

N_SIM = 5000
SEED = 12345

FAMILIES = [11, 101, 191]
FAMILY_SET = set(FAMILIES)

CYCLE_PLUS = {
    (11, 101, 191),
    (101, 191, 11),
    (191, 11, 101),
}

CYCLE_MINUS = {
    (11, 191, 101),
    (191, 101, 11),
    (101, 11, 191),
}


# ------------------------------------------------------------
# Primzahlsieb
# ------------------------------------------------------------

def sieve_bool(n: int) -> np.ndarray:
    is_prime = np.ones(n + 1, dtype=bool)
    is_prime[:2] = False

    limit = int(n ** 0.5)
    for k in range(2, limit + 1):
        if is_prime[k]:
            is_prime[k*k:n+1:k] = False

    return is_prime


def find_prime_quadruplet_starts(N: int) -> np.ndarray:
    """
    Startwerte p mit p,p+2,p+6,p+8 <= N prim.
    """
    is_prime = sieve_bool(N)

    max_p = N - 8
    candidates = np.arange(2, max_p + 1, dtype=np.int64)

    mask = (
        is_prime[candidates]
        & is_prime[candidates + 2]
        & is_prime[candidates + 6]
        & is_prime[candidates + 8]
    )

    return candidates[mask]


# ------------------------------------------------------------
# Tripel-Statistik
# ------------------------------------------------------------

def count_triplet_statistics(families: np.ndarray):
    """
    Zählt neutrale, monochrome und orientierte Tripel.
    """
    n = len(families)
    n_triples = max(0, n - 2)

    neutral = 0
    mono = 0
    plus = 0
    minus = 0

    triple_counts = Counter()

    for i in range(n_triples):
        tri = tuple(int(x) for x in families[i:i+3])
        triple_counts[tri] += 1

        if set(tri) == FAMILY_SET:
            neutral += 1
        if len(set(tri)) == 1:
            mono += 1
        if tri in CYCLE_PLUS:
            plus += 1
        if tri in CYCLE_MINUS:
            minus += 1

    return {
        "n_triples": n_triples,
        "neutral": neutral,
        "mono": mono,
        "plus": plus,
        "minus": minus,
        "orient_diff": plus - minus,
        "triple_counts": triple_counts,
    }


def permutation_baseline(families: np.ndarray, n_sim: int, seed: int):
    """
    Permutationstest:
    Erhält exakt die Familienhäufigkeiten, zerstört aber Reihenfolge.
    """
    rng = np.random.default_rng(seed)

    obs = count_triplet_statistics(families)

    sim_neutral = np.zeros(n_sim, dtype=np.int64)
    sim_mono = np.zeros(n_sim, dtype=np.int64)
    sim_orient_diff = np.zeros(n_sim, dtype=np.int64)

    for s in range(n_sim):
        perm = rng.permutation(families)
        stat = count_triplet_statistics(perm)

        sim_neutral[s] = stat["neutral"]
        sim_mono[s] = stat["mono"]
        sim_orient_diff[s] = stat["orient_diff"]

    def p_greater(obs_value, sims):
        return (np.sum(sims >= obs_value) + 1) / (len(sims) + 1)

    def p_less(obs_value, sims):
        return (np.sum(sims <= obs_value) + 1) / (len(sims) + 1)

    def p_two_sided(obs_value, sims):
        center = np.mean(sims)
        return (np.sum(np.abs(sims - center) >= abs(obs_value - center)) + 1) / (len(sims) + 1)

    return {
        "obs": obs,

        "neutral_perm_mean": float(sim_neutral.mean()),
        "neutral_perm_sd": float(sim_neutral.std(ddof=1)),
        "neutral_p_greater": float(p_greater(obs["neutral"], sim_neutral)),
        "neutral_p_less": float(p_less(obs["neutral"], sim_neutral)),

        "mono_perm_mean": float(sim_mono.mean()),
        "mono_perm_sd": float(sim_mono.std(ddof=1)),
        "mono_p_greater": float(p_greater(obs["mono"], sim_mono)),
        "mono_p_less": float(p_less(obs["mono"], sim_mono)),

        "orient_perm_mean": float(sim_orient_diff.mean()),
        "orient_perm_sd": float(sim_orient_diff.std(ddof=1)),
        "orient_p_two_sided": float(p_two_sided(obs["orient_diff"], sim_orient_diff)),
    }


# ------------------------------------------------------------
# Übergangsmatrix
# ------------------------------------------------------------

def transition_matrix_3(families: np.ndarray):
    idx = {v: i for i, v in enumerate(FAMILIES)}
    M = np.zeros((3, 3), dtype=np.int64)

    for a, b in zip(families[:-1], families[1:]):
        a = int(a)
        b = int(b)
        if a in idx and b in idx:
            M[idx[a], idx[b]] += 1

    return M


# ------------------------------------------------------------
# Ein N auswerten
# ------------------------------------------------------------

def analyze_N(N: int, n_sim: int, seed: int):
    t0 = time.time()

    starts = find_prime_quadruplet_starts(N)
    residues = starts % 210

    counts210 = Counter(int(x) for x in residues)

    anchor_count = counts210.get(5, 0)

    mask_main = np.isin(residues, FAMILIES)
    families = residues[mask_main].astype(np.int64)

    fam_counts = Counter(int(x) for x in families)

    n_quads = len(starts)
    n_main = len(families)

    # Gleichverteilung der drei Familien
    fam_obs = np.array([
        fam_counts.get(11, 0),
        fam_counts.get(101, 0),
        fam_counts.get(191, 0),
    ], dtype=float)

    if SCIPY_OK and n_main > 0:
        fam_exp = np.array([n_main / 3] * 3)
        chi, p_family_equal = chisquare(fam_obs, f_exp=fam_exp)
    else:
        chi, p_family_equal = float("nan"), float("nan")

    # Übergangsmatrix
    M = transition_matrix_3(families)

    # Tripel + Permutation
    perm = permutation_baseline(families, n_sim=n_sim, seed=seed + N)
    obs = perm["obs"]

    n_triples = obs["n_triples"]

    elapsed = time.time() - t0

    row = {
        "N": N,
        "quadruplets": n_quads,
        "anchor_5": anchor_count,
        "family_total": n_main,

        "family_11": fam_counts.get(11, 0),
        "family_101": fam_counts.get(101, 0),
        "family_191": fam_counts.get(191, 0),

        "family_chi2": chi,
        "family_equal_p": p_family_equal,

        "trans_11_11": M[0, 0],
        "trans_11_101": M[0, 1],
        "trans_11_191": M[0, 2],
        "trans_101_11": M[1, 0],
        "trans_101_101": M[1, 1],
        "trans_101_191": M[1, 2],
        "trans_191_11": M[2, 0],
        "trans_191_101": M[2, 1],
        "trans_191_191": M[2, 2],

        "triples": n_triples,

        "neutral": obs["neutral"],
        "neutral_rate": obs["neutral"] / n_triples if n_triples > 0 else float("nan"),
        "neutral_perm_mean": perm["neutral_perm_mean"],
        "neutral_perm_sd": perm["neutral_perm_sd"],
        "neutral_p_greater": perm["neutral_p_greater"],
        "neutral_p_less": perm["neutral_p_less"],

        "mono": obs["mono"],
        "mono_rate": obs["mono"] / n_triples if n_triples > 0 else float("nan"),
        "mono_perm_mean": perm["mono_perm_mean"],
        "mono_perm_sd": perm["mono_perm_sd"],
        "mono_p_greater": perm["mono_p_greater"],
        "mono_p_less": perm["mono_p_less"],

        "plus": obs["plus"],
        "minus": obs["minus"],
        "orient_diff": obs["orient_diff"],
        "orient_perm_mean": perm["orient_perm_mean"],
        "orient_perm_sd": perm["orient_perm_sd"],
        "orient_p_two_sided": perm["orient_p_two_sided"],

        "elapsed_sec": elapsed,
    }

    return row


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    print("\n============================================================")
    print("BM 210 Skalierungstest")
    print("============================================================")
    print(f"N_VALUES = {N_VALUES}")
    print(f"N_SIM    = {N_SIM}")
    print(f"SEED     = {SEED}")

    rows = []

    for N in N_VALUES:
        print("\n------------------------------------------------------------")
        print(f"Analysiere N={N}")
        print("------------------------------------------------------------")

        row = analyze_N(N, n_sim=N_SIM, seed=SEED)
        rows.append(row)

        print(f"Vierlinge:      {row['quadruplets']}")
        print(f"Anker 5:        {row['anchor_5']}")
        print(
            f"Familien:       11={row['family_11']}, "
            f"101={row['family_101']}, "
            f"191={row['family_191']}"
        )
        print(f"family_equal_p: {row['family_equal_p']:.6f}")

        print(
            f"Neutral:        {row['neutral']} "
            f"rate={row['neutral_rate']:.6f}, "
            f"p_greater={row['neutral_p_greater']:.6f}, "
            f"p_less={row['neutral_p_less']:.6f}"
        )

        print(
            f"Mono:           {row['mono']} "
            f"rate={row['mono_rate']:.6f}, "
            f"p_greater={row['mono_p_greater']:.6f}, "
            f"p_less={row['mono_p_less']:.6f}"
        )

        print(
            f"Orient:         plus={row['plus']}, minus={row['minus']}, "
            f"diff={row['orient_diff']}, "
            f"p_two={row['orient_p_two_sided']:.6f}"
        )

        print(f"Laufzeit:       {row['elapsed_sec']:.2f} s")

        # Zwischenspeichern nach jedem N
        df = pd.DataFrame(rows)
        df.to_csv("bm210_scaling_results.csv", index=False)

    print("\n============================================================")
    print("Fertig")
    print("============================================================")
    print("CSV gespeichert: bm210_scaling_results.csv")

    df = pd.DataFrame(rows)

    print("\nKompaktübersicht:")
    cols = [
        "N",
        "quadruplets",
        "family_11",
        "family_101",
        "family_191",
        "family_equal_p",
        "neutral_rate",
        "neutral_p_greater",
        "mono_rate",
        "mono_p_less",
        "orient_diff",
        "orient_p_two_sided",
    ]
    print(df[cols].to_string(index=False))

    print("\nInterpretation:")
    print("""
1. Stabiler 1⊕3-Befund:
   Nur die Klassen 5,11,101,191 sollten auftreten.
   Die Klasse 5 bleibt der Anker.

2. QCD-Singulett-Analogon:
   neutral_p_greater < 0.01 über mehrere N wäre stark.
   Das hieße: vollständige Tripel {11,101,191} sind überproduziert.

3. Ausschluss-/Pauli-Analogon:
   mono_p_less < 0.01 über mehrere N wäre stark.
   Das hieße: FFF-Tripel werden systematisch vermieden.

4. Chiralitäts-Analogon:
   orient_p_two_sided < 0.01 mit stabilem Vorzeichen der orient_diff
   wäre ein Hinweis auf bevorzugte Drehrichtung.

5. Falls alles p≈0.1...0.9 bleibt:
   Die 1⊕3-Struktur ist real, aber die Reihenfolge verhält sich weitgehend
   permutationell zufällig.
""")


if __name__ == "__main__":
    main()