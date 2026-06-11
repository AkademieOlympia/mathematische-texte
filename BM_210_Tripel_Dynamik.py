#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BM_210_Tripel_Dynamik.py

Erzeugt Primzahlvierlinge Q(p)=(p,p+2,p+6,p+8) bis N
und testet die mod-210 3+1-Struktur:

    {5} ∪ {11,101,191}

Getestet werden:
- Häufigkeiten der mod-210-Familien
- Tripelneutralität der drei Hauptfamilien
- Monochrome Tripel
- zyklische Orientierung
- Übergangsmatrix 3×3
- Permutations-Baselines
"""

import sys
import math
import numpy as np
from collections import Counter

try:
    from scipy.stats import chisquare
    SCIPY_OK = True
except Exception:
    SCIPY_OK = False


def sieve_bool(n: int) -> np.ndarray:
    """Primzahlsieb bis n."""
    is_prime = np.ones(n + 1, dtype=bool)
    is_prime[:2] = False

    limit = int(n ** 0.5)
    for k in range(2, limit + 1):
        if is_prime[k]:
            is_prime[k*k:n+1:k] = False

    return is_prime


def find_prime_quadruplets(N: int):
    """
    Findet Startwerte p mit:
        p, p+2, p+6, p+8 <= N
    alle prim.
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


def transition_matrix(families, fam_values):
    idx = {v: i for i, v in enumerate(fam_values)}
    M = np.zeros((len(fam_values), len(fam_values)), dtype=int)

    for a, b in zip(families[:-1], families[1:]):
        if a in idx and b in idx:
            M[idx[a], idx[b]] += 1

    return M


def permutation_test_counts(families, allowed, n_sim=10000, seed=12345):
    rng = np.random.default_rng(seed)

    n = len(families)

    obs_neutral = 0
    obs_mono = 0
    obs_plus = 0
    obs_minus = 0

    cycle_plus = {
        (11, 101, 191),
        (101, 191, 11),
        (191, 11, 101),
    }

    cycle_minus = {
        (11, 191, 101),
        (191, 101, 11),
        (101, 11, 191),
    }

    triples = []
    for i in range(n - 2):
        tri = tuple(families[i:i+3])
        triples.append(tri)

        if set(tri) == allowed:
            obs_neutral += 1
        if len(set(tri)) == 1:
            obs_mono += 1
        if tri in cycle_plus:
            obs_plus += 1
        if tri in cycle_minus:
            obs_minus += 1

    sim_neutral = []
    sim_mono = []
    sim_plus_minus_diff = []

    for _ in range(n_sim):
        perm = rng.permutation(families)

        neu = 0
        mono = 0
        plus = 0
        minus = 0

        for i in range(n - 2):
            tri = tuple(perm[i:i+3])

            if set(tri) == allowed:
                neu += 1
            if len(set(tri)) == 1:
                mono += 1
            if tri in cycle_plus:
                plus += 1
            if tri in cycle_minus:
                minus += 1

        sim_neutral.append(neu)
        sim_mono.append(mono)
        sim_plus_minus_diff.append(plus - minus)

    sim_neutral = np.array(sim_neutral)
    sim_mono = np.array(sim_mono)
    sim_plus_minus_diff = np.array(sim_plus_minus_diff)

    def p_greater(obs, sims):
        return (np.sum(sims >= obs) + 1) / (len(sims) + 1)

    def p_less(obs, sims):
        return (np.sum(sims <= obs) + 1) / (len(sims) + 1)

    def p_twosided(obs, sims):
        center = np.mean(sims)
        return (np.sum(np.abs(sims - center) >= abs(obs - center)) + 1) / (len(sims) + 1)

    return {
        "triples": triples,
        "obs_neutral": obs_neutral,
        "obs_mono": obs_mono,
        "obs_plus": obs_plus,
        "obs_minus": obs_minus,
        "obs_orient_diff": obs_plus - obs_minus,

        "sim_neutral_mean": float(sim_neutral.mean()),
        "sim_neutral_sd": float(sim_neutral.std(ddof=1)),
        "p_neutral_greater": p_greater(obs_neutral, sim_neutral),
        "p_neutral_less": p_less(obs_neutral, sim_neutral),

        "sim_mono_mean": float(sim_mono.mean()),
        "sim_mono_sd": float(sim_mono.std(ddof=1)),
        "p_mono_greater": p_greater(obs_mono, sim_mono),
        "p_mono_less": p_less(obs_mono, sim_mono),

        "sim_orient_diff_mean": float(sim_plus_minus_diff.mean()),
        "sim_orient_diff_sd": float(sim_plus_minus_diff.std(ddof=1)),
        "p_orient_twosided": p_twosided(obs_plus - obs_minus, sim_plus_minus_diff),
    }


def main():
    if len(sys.argv) >= 2:
        N = int(sys.argv[1])
    else:
        N = 10_000_000

    if len(sys.argv) >= 3:
        N_SIM = int(sys.argv[2])
    else:
        N_SIM = 10000

    print("\n============================================================")
    print("BM mod-210 Tripel-Dynamik")
    print("============================================================")
    print(f"N = {N}")
    print(f"Permutationen = {N_SIM}")

    starts = find_prime_quadruplets(N)
    centers = starts + 4
    residues = starts % 210

    print("\n------------------------------------------------------------")
    print("1. Primzahlvierlinge")
    print("------------------------------------------------------------")
    print(f"Anzahl Vierlinge: {len(starts)}")
    print(f"Erster Start:     {starts[0] if len(starts) else None}")
    print(f"Letzter Start:    {starts[-1] if len(starts) else None}")

    unique = sorted(set(residues.tolist()))
    print(f"Unique residues mod 210: {unique}")

    counts210 = Counter(residues.tolist())

    print("\nmod-210-Häufigkeiten der relevanten Klassen:")
    for r in [5, 11, 101, 191]:
        print(f"{r:3d}: {counts210.get(r, 0)}")

    other = sum(c for r, c in counts210.items() if r not in {5, 11, 101, 191})
    print(f"andere Klassen: {other}")

    # Hauptfamilien nach Entfernung des Ankers
    main_mask = np.isin(residues, [11, 101, 191])
    families = residues[main_mask]

    print("\n------------------------------------------------------------")
    print("2. Drei Hauptfamilien")
    print("------------------------------------------------------------")
    print(f"Anzahl ohne Ankerklasse: {len(families)}")

    fam_counts = Counter(families.tolist())
    for r in [11, 101, 191]:
        print(f"{r:3d}: {fam_counts.get(r, 0)}")

    obs = np.array([fam_counts.get(11, 0), fam_counts.get(101, 0), fam_counts.get(191, 0)])
    exp = np.array([len(families) / 3] * 3)

    print("\nGleichverteilungstest der drei Hauptfamilien:")
    print("observed:", obs.tolist())
    print("expected:", [round(x, 3) for x in exp])

    if SCIPY_OK:
        chi, pval = chisquare(obs, f_exp=exp)
        print(f"chi²={chi:.4f}, p={pval:.6g}")
    else:
        print("scipy nicht verfügbar.")

    # Übergangsmatrix
    print("\n------------------------------------------------------------")
    print("3. Übergangsmatrix 3×3")
    print("------------------------------------------------------------")

    fam_values = [11, 101, 191]
    M = transition_matrix(families, fam_values)

    print("Zeilen/Spalten:", fam_values)
    print(M)
    print("Zeilensummen:", M.sum(axis=1).tolist())

    print("\nZeilenweise χ² gegen Gleichverteilung:")
    if SCIPY_OK:
        for i, r in enumerate(fam_values):
            row = M[i]
            nrow = row.sum()
            if nrow > 0:
                exp_row = np.array([nrow / 3] * 3)
                chi, pval = chisquare(row, f_exp=exp_row)
                print(f"Start {r:3d}: row={row.tolist()}, chi²={chi:.4f}, p={pval:.6g}")
    else:
        print("scipy nicht verfügbar.")

    # Tripeltests
    print("\n------------------------------------------------------------")
    print("4. Tripelneutralität / Monochromie / Orientierung")
    print("------------------------------------------------------------")

    allowed = {11, 101, 191}
    res = permutation_test_counts(
        families=families,
        allowed=allowed,
        n_sim=N_SIM,
        seed=12345,
    )

    n_triples = max(0, len(families) - 2)
    print(f"Gleitende Dreierfenster: {n_triples}")

    print("\nTripelneutralität:")
    print(f"beobachtet neutral: {res['obs_neutral']}")
    print(f"Rate:               {res['obs_neutral'] / n_triples:.6f}")
    print(f"Permutation mean:   {res['sim_neutral_mean']:.3f}")
    print(f"Permutation sd:     {res['sim_neutral_sd']:.3f}")
    print(f"p_greater:          {res['p_neutral_greater']:.6f}")
    print(f"p_less:             {res['p_neutral_less']:.6f}")

    print("\nMonochrome Tripel:")
    print(f"beobachtet mono:    {res['obs_mono']}")
    print(f"Rate:               {res['obs_mono'] / n_triples:.6f}")
    print(f"Permutation mean:   {res['sim_mono_mean']:.3f}")
    print(f"Permutation sd:     {res['sim_mono_sd']:.3f}")
    print(f"p_greater:          {res['p_mono_greater']:.6f}")
    print(f"p_less:             {res['p_mono_less']:.6f}")

    print("\nZyklische Orientierung neutraler Tripel:")
    print(f"plus:               {res['obs_plus']}")
    print(f"minus:              {res['obs_minus']}")
    print(f"Differenz:          {res['obs_orient_diff']}")
    print(f"Perm diff mean:     {res['sim_orient_diff_mean']:.3f}")
    print(f"Perm diff sd:       {res['sim_orient_diff_sd']:.3f}")
    print(f"p_two_sided:        {res['p_orient_twosided']:.6f}")

    # Top Tripel
    print("\n------------------------------------------------------------")
    print("5. Häufigste Tripel")
    print("------------------------------------------------------------")

    triple_counts = Counter(res["triples"])
    for tri, c in triple_counts.most_common(15):
        print(f"{tri}: {c}")

    # Export kleiner CSV
    out = f"bm210_triplet_summary_N{N}.csv"
    with open(out, "w", encoding="utf-8") as f:
        f.write("metric,value\n")
        f.write(f"N,{N}\n")
        f.write(f"quadruplets,{len(starts)}\n")
        f.write(f"family_11,{fam_counts.get(11, 0)}\n")
        f.write(f"family_101,{fam_counts.get(101, 0)}\n")
        f.write(f"family_191,{fam_counts.get(191, 0)}\n")
        f.write(f"triples,{n_triples}\n")
        f.write(f"neutral,{res['obs_neutral']}\n")
        f.write(f"neutral_rate,{res['obs_neutral'] / n_triples}\n")
        f.write(f"neutral_perm_mean,{res['sim_neutral_mean']}\n")
        f.write(f"neutral_p_greater,{res['p_neutral_greater']}\n")
        f.write(f"neutral_p_less,{res['p_neutral_less']}\n")
        f.write(f"mono,{res['obs_mono']}\n")
        f.write(f"mono_rate,{res['obs_mono'] / n_triples}\n")
        f.write(f"mono_perm_mean,{res['sim_mono_mean']}\n")
        f.write(f"mono_p_greater,{res['p_mono_greater']}\n")
        f.write(f"mono_p_less,{res['p_mono_less']}\n")
        f.write(f"plus,{res['obs_plus']}\n")
        f.write(f"minus,{res['obs_minus']}\n")
        f.write(f"orient_diff,{res['obs_orient_diff']}\n")
        f.write(f"orient_p_two_sided,{res['p_orient_twosided']}\n")

    print("\n------------------------------------------------------------")
    print("6. Export")
    print("------------------------------------------------------------")
    print(f"Summary gespeichert: {out}")

    print("\n============================================================")
    print("Kurzfazit-Hinweis")
    print("============================================================")
    print("""
Interpretation:
- p_neutral_greater klein: überzufällige RGB-/Singulettbildung.
- p_neutral_less klein: Unterdrückung vollständiger Dreiermischung.
- p_mono_less klein: Vermeidung monochromer FFF-Blöcke.
- p_orient_two_sided klein: bevorzugte chirale Orientierung.

Für BM/QCD wäre besonders interessant:
1. p_neutral_greater < 0.01
   oder
2. p_mono_less < 0.01 stabil über N hinweg
   oder
3. p_orient_two_sided < 0.01 mit stabiler Vorzeichenrichtung.
""")


if __name__ == "__main__":
    main()