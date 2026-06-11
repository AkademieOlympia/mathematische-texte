#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BM_210_Operator_Permutation.py

Permutationstest für den SU(3)-artigen Oktettoperator der drei
mod-210-Hauptfamilien der Primzahlvierlinge.

Grundstruktur:
    {5} ∪ {11,101,191}

Nach Entfernung der Ankerklasse 5 wird die Folge der Familien
    11,101,191
betrachtet.

Aus der Übergangsmatrix M wird ein Markov-Kern P gebaut.
Dann:
    Delta = P - J/3

Getestet wird:
    ||Delta||_F

gegen Permutationsbaseline:
    gleiche Familienzahlen, zufällige Reihenfolge.

Zusätzlich:
    ||A||/||Delta||
    J_CP-artiger Wert
    Eigenwerte von H = S + iA
    PMNS-Winkel, aber nur als Nebenprodukt
"""

import time
import math
import numpy as np
import pandas as pd
from collections import Counter


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
    Findet Startwerte p mit p,p+2,p+6,p+8 <= N prim.
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
# Operatorfunktionen
# ------------------------------------------------------------

def transition_matrix_3(seq: np.ndarray) -> np.ndarray:
    idx = {v: i for i, v in enumerate(FAMILIES)}
    M = np.zeros((3, 3), dtype=float)

    for a, b in zip(seq[:-1], seq[1:]):
        a = int(a)
        b = int(b)
        if a in idx and b in idx:
            M[idx[a], idx[b]] += 1.0

    return M


def row_normalize(M: np.ndarray) -> np.ndarray:
    P = M.astype(float).copy()
    row_sums = P.sum(axis=1)

    for i in range(3):
        if row_sums[i] > 0:
            P[i, :] /= row_sums[i]
        else:
            P[i, :] = 1.0 / 3.0

    return P


def phase_fix_columns(U: np.ndarray) -> np.ndarray:
    U = U.copy()
    for j in range(U.shape[1]):
        k = np.argmax(np.abs(U[:, j]))
        phase = np.angle(U[k, j])
        U[:, j] *= np.exp(-1j * phase)
    return U


def operator_observables_from_sequence(seq: np.ndarray):
    """
    Baut M, P, Delta, H aus Familienfolge.
    """
    M = transition_matrix_3(seq)
    P = row_normalize(M)

    J = np.ones((3, 3), dtype=float)
    P0 = J / 3.0

    Delta = P - P0
    S = 0.5 * (Delta + Delta.T)
    A = 0.5 * (Delta - Delta.T)

    H = S + 1j * A
    H = H - np.trace(H) / 3.0 * np.eye(3)

    norm_delta = np.linalg.norm(Delta, "fro")
    norm_S = np.linalg.norm(S, "fro")
    norm_A = np.linalg.norm(A, "fro")

    relative_octet = norm_delta  # weil ||J/3||_F = 1
    chiral_fraction = norm_A / norm_delta if norm_delta > 0 else np.nan

    # Eigenwerte H
    eigvals, eigvecs = np.linalg.eigh(H)
    idx = np.argsort(eigvals)
    E = eigvals[idx].real
    U = eigvecs[:, idx]
    U = phase_fix_columns(U)

    absU2 = np.abs(U) ** 2

    # PMNS-artige Winkel, nur diagnostisch
    s13_sq = absU2[0, 2]
    c13_sq = 1.0 - s13_sq

    if c13_sq > 0:
        s12_sq = absU2[0, 1] / c13_sq
        s23_sq = absU2[1, 2] / c13_sq
    else:
        s12_sq = np.nan
        s23_sq = np.nan

    def safe_angle(s2):
        if np.isnan(s2):
            return np.nan
        s2 = max(0.0, min(1.0, float(s2)))
        return math.degrees(math.asin(math.sqrt(s2)))

    theta12 = safe_angle(s12_sq)
    theta13 = safe_angle(s13_sq)
    theta23 = safe_angle(s23_sq)

    Jcp = np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0]))

    d21 = E[1] - E[0]
    d31 = E[2] - E[0]
    ratio_21_31 = abs(d21 / d31) if d31 != 0 else np.nan

    return {
        "M": M,
        "P": P,
        "Delta": Delta,
        "S": S,
        "A": A,
        "H": H,

        "norm_delta": norm_delta,
        "norm_S": norm_S,
        "norm_A": norm_A,
        "relative_octet": relative_octet,
        "chiral_fraction": chiral_fraction,

        "E1": E[0],
        "E2": E[1],
        "E3": E[2],
        "ratio_21_31": ratio_21_31,

        "theta12": theta12,
        "theta13": theta13,
        "theta23": theta23,
        "Jcp": Jcp,
    }


def empirical_pvalues(obs, sims):
    sims = np.asarray(sims)

    p_greater = (np.sum(sims >= obs) + 1) / (len(sims) + 1)
    p_less = (np.sum(sims <= obs) + 1) / (len(sims) + 1)

    center = np.mean(sims)
    p_two = (np.sum(np.abs(sims - center) >= abs(obs - center)) + 1) / (len(sims) + 1)

    z = (obs - np.mean(sims)) / np.std(sims, ddof=1) if np.std(sims, ddof=1) > 0 else np.nan

    return p_greater, p_less, p_two, z


# ------------------------------------------------------------
# Ein N analysieren
# ------------------------------------------------------------

def analyze_N(N: int, n_sim: int, seed: int):
    t0 = time.time()

    starts = find_prime_quadruplet_starts(N)
    residues = starts % 210

    counts = Counter(int(x) for x in residues)

    mask = np.isin(residues, FAMILIES)
    seq = residues[mask].astype(np.int64)

    obs = operator_observables_from_sequence(seq)

    rng = np.random.default_rng(seed + N)

    sim_delta = np.zeros(n_sim)
    sim_chiral_fraction = np.zeros(n_sim)
    sim_norm_A = np.zeros(n_sim)
    sim_Jcp = np.zeros(n_sim)
    sim_ratio = np.zeros(n_sim)

    for s in range(n_sim):
        perm = rng.permutation(seq)
        sim_obs = operator_observables_from_sequence(perm)

        sim_delta[s] = sim_obs["norm_delta"]
        sim_chiral_fraction[s] = sim_obs["chiral_fraction"]
        sim_norm_A[s] = sim_obs["norm_A"]
        sim_Jcp[s] = sim_obs["Jcp"]
        sim_ratio[s] = sim_obs["ratio_21_31"]

    p_delta_g, p_delta_l, p_delta_two, z_delta = empirical_pvalues(
        obs["norm_delta"], sim_delta
    )

    p_chi_g, p_chi_l, p_chi_two, z_chi = empirical_pvalues(
        obs["chiral_fraction"], sim_chiral_fraction
    )

    p_A_g, p_A_l, p_A_two, z_A = empirical_pvalues(
        obs["norm_A"], sim_norm_A
    )

    p_J_g, p_J_l, p_J_two, z_J = empirical_pvalues(
        abs(obs["Jcp"]), np.abs(sim_Jcp)
    )

    p_ratio_g, p_ratio_l, p_ratio_two, z_ratio = empirical_pvalues(
        obs["ratio_21_31"], sim_ratio
    )

    elapsed = time.time() - t0

    row = {
        "N": N,
        "quadruplets": len(starts),
        "anchor_5": counts.get(5, 0),
        "family_total": len(seq),
        "family_11": counts.get(11, 0),
        "family_101": counts.get(101, 0),
        "family_191": counts.get(191, 0),

        "obs_delta": obs["norm_delta"],
        "sim_delta_mean": float(sim_delta.mean()),
        "sim_delta_sd": float(sim_delta.std(ddof=1)),
        "delta_z": z_delta,
        "delta_p_greater": p_delta_g,
        "delta_p_less": p_delta_l,
        "delta_p_two": p_delta_two,

        "obs_norm_A": obs["norm_A"],
        "sim_norm_A_mean": float(sim_norm_A.mean()),
        "sim_norm_A_sd": float(sim_norm_A.std(ddof=1)),
        "norm_A_z": z_A,
        "norm_A_p_greater": p_A_g,
        "norm_A_p_less": p_A_l,
        "norm_A_p_two": p_A_two,

        "obs_chiral_fraction": obs["chiral_fraction"],
        "sim_chiral_fraction_mean": float(sim_chiral_fraction.mean()),
        "sim_chiral_fraction_sd": float(sim_chiral_fraction.std(ddof=1)),
        "chiral_fraction_z": z_chi,
        "chiral_fraction_p_greater": p_chi_g,
        "chiral_fraction_p_less": p_chi_l,
        "chiral_fraction_p_two": p_chi_two,

        "obs_Jcp": obs["Jcp"],
        "abs_Jcp": abs(obs["Jcp"]),
        "sim_abs_Jcp_mean": float(np.abs(sim_Jcp).mean()),
        "sim_abs_Jcp_sd": float(np.abs(sim_Jcp).std(ddof=1)),
        "abs_Jcp_z": z_J,
        "abs_Jcp_p_greater": p_J_g,
        "abs_Jcp_p_less": p_J_l,
        "abs_Jcp_p_two": p_J_two,

        "obs_ratio_21_31": obs["ratio_21_31"],
        "sim_ratio_mean": float(np.nanmean(sim_ratio)),
        "sim_ratio_sd": float(np.nanstd(sim_ratio, ddof=1)),
        "ratio_z": z_ratio,
        "ratio_p_greater": p_ratio_g,
        "ratio_p_less": p_ratio_l,
        "ratio_p_two": p_ratio_two,

        "theta12": obs["theta12"],
        "theta13": obs["theta13"],
        "theta23": obs["theta23"],

        "elapsed_sec": elapsed,
    }

    return row, obs, sim_delta


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    print("\n============================================================")
    print("BM mod-210 Operator-Permutationstest")
    print("============================================================")
    print("Familien:", FAMILIES)
    print("N_VALUES:", N_VALUES)
    print("N_SIM:", N_SIM)

    rows = []

    last = None

    for N in N_VALUES:
        print("\n------------------------------------------------------------")
        print(f"Analysiere N={N}")
        print("------------------------------------------------------------")

        row, obs, sim_delta = analyze_N(N, N_SIM, SEED)
        rows.append(row)
        last = (N, row, obs, sim_delta)

        print(f"Vierlinge:          {row['quadruplets']}")
        print(f"Familien total:     {row['family_total']}")
        print(
            f"Familien:           11={row['family_11']}, "
            f"101={row['family_101']}, "
            f"191={row['family_191']}"
        )

        print("\nOktettstärke ||Delta||_F:")
        print(f"  obs:              {row['obs_delta']:.6f}")
        print(f"  sim mean:         {row['sim_delta_mean']:.6f}")
        print(f"  sim sd:           {row['sim_delta_sd']:.6f}")
        print(f"  z:                {row['delta_z']:+.3f}")
        print(f"  p_greater:        {row['delta_p_greater']:.6f}")
        print(f"  p_less:           {row['delta_p_less']:.6f}")

        print("\nAntisymmetrischer Anteil ||A||:")
        print(f"  obs:              {row['obs_norm_A']:.6f}")
        print(f"  sim mean:         {row['sim_norm_A_mean']:.6f}")
        print(f"  z:                {row['norm_A_z']:+.3f}")
        print(f"  p_greater:        {row['norm_A_p_greater']:.6f}")

        print("\nChiralitätsfraktion ||A||/||Delta||:")
        print(f"  obs:              {row['obs_chiral_fraction']:.6f}")
        print(f"  sim mean:         {row['sim_chiral_fraction_mean']:.6f}")
        print(f"  z:                {row['chiral_fraction_z']:+.3f}")
        print(f"  p_greater:        {row['chiral_fraction_p_greater']:.6f}")

        print("\nJarlskog-artig |Jcp|:")
        print(f"  obs |J|:          {row['abs_Jcp']:.6f}")
        print(f"  sim mean |J|:     {row['sim_abs_Jcp_mean']:.6f}")
        print(f"  z:                {row['abs_Jcp_z']:+.3f}")
        print(f"  p_greater:        {row['abs_Jcp_p_greater']:.6f}")

        print(f"\nWinkel diagnostisch: θ12={row['theta12']:.2f}, θ13={row['theta13']:.2f}, θ23={row['theta23']:.2f}")
        print(f"Laufzeit: {row['elapsed_sec']:.2f} s")

        pd.DataFrame(rows).to_csv("bm210_operator_permutation_results.csv", index=False)

    print("\n============================================================")
    print("Fertig")
    print("============================================================")
    print("CSV gespeichert: bm210_operator_permutation_results.csv")

    df = pd.DataFrame(rows)

    cols = [
        "N",
        "quadruplets",
        "obs_delta",
        "sim_delta_mean",
        "delta_z",
        "delta_p_greater",
        "obs_norm_A",
        "sim_norm_A_mean",
        "norm_A_z",
        "norm_A_p_greater",
        "obs_chiral_fraction",
        "sim_chiral_fraction_mean",
        "chiral_fraction_z",
        "chiral_fraction_p_greater",
        "abs_Jcp",
        "sim_abs_Jcp_mean",
        "abs_Jcp_z",
        "abs_Jcp_p_greater",
        "theta13",
    ]

    print("\nKompaktübersicht:")
    print(df[cols].to_string(index=False))

    print("\nInterpretation:")
    print("""
1. delta_p_greater klein:
   Der echte Oktettoperator ist stärker als bei zufälliger Reihenfolge.
   Dann wäre der Oktettrest sequenziell real.

2. delta_p_greater unauffällig:
   Der Oktettrest ist kompatibel mit Permutationsrauschen.

3. norm_A_p_greater klein:
   Der antisymmetrische/chirale Anteil ist stärker als Zufall.

4. abs_Jcp_p_greater klein:
   Der CP-artige Jarlskog-Wert ist stärker als Zufall.

5. Wenn alle p-Werte unauffällig bleiben:
   Die 1⊕3-Familienstruktur ist real, aber der Neutrino-/SU(3)-Operator
   entsteht aus Reihenfolgenrauschen.
""")


if __name__ == "__main__":
    main()
    clear