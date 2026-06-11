#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BM_210_Neutrino_Scaling.py

Skalierungstest des BM mod-210 Neutrino-Analogons.

Für mehrere Schranken N:
1. Primzahlvierlinge Q(p)=(p,p+2,p+6,p+8) erzeugen.
2. Startklassen p mod 210 bestimmen.
3. Ankerklasse 5 entfernen.
4. Drei Hauptfamilien 11,101,191 als Flavorbasis benutzen.
5. Übergangsmatrix P_N bauen.
6. Oktettanteil Delta = P_N - J/3 bestimmen.
7. Hermiteschen Hamiltonian H = S + iA bilden.
8. PMNS-artige Winkel, Jarlskog-Indikator und Splitting-Verhältnis extrahieren.

Ausgabe:
    bm210_neutrino_scaling.csv
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
# Übergangsmatrix
# ------------------------------------------------------------

def transition_matrix_3(families: np.ndarray) -> np.ndarray:
    idx = {v: i for i, v in enumerate(FAMILIES)}
    M = np.zeros((3, 3), dtype=float)

    for a, b in zip(families[:-1], families[1:]):
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


# ------------------------------------------------------------
# Neutrino-Analogon
# ------------------------------------------------------------

def phase_fix_columns(U: np.ndarray) -> np.ndarray:
    """
    Dreht jede Eigenvektor-Spalte so, dass ihr größter Eintrag reell positiv ist.
    """
    U = U.copy()
    for j in range(U.shape[1]):
        k = np.argmax(np.abs(U[:, j]))
        phase = np.angle(U[k, j])
        U[:, j] *= np.exp(-1j * phase)
    return U


def neutrino_observables_from_P(P: np.ndarray):
    """
    Aus P wird:
        Delta = P - J/3
        H = Sym(Delta) + i Anti(Delta)
    gebildet.

    Danach Diagonalisierung und Extraktion:
        theta12, theta13, theta23
        J_CP
        Splitting-Verhältnis
        Oktettstärke
        Chiralitätsanteil
    """

    J = np.ones((3, 3), dtype=float)
    P_singlet = J / 3.0

    Delta = P - P_singlet

    S = 0.5 * (Delta + Delta.T)
    A = 0.5 * (Delta - Delta.T)

    H = S + 1j * A

    # Spur entfernen
    H = H - np.trace(H) / 3.0 * np.eye(3)

    herm_error = np.linalg.norm(H - np.conjugate(H.T), ord="fro")

    eigvals, eigvecs = np.linalg.eigh(H)

    idx = np.argsort(eigvals)
    E = eigvals[idx].real
    U = eigvecs[:, idx]
    U = phase_fix_columns(U)

    unit_error = np.linalg.norm(np.conjugate(U.T) @ U - np.eye(3), ord="fro")

    absU2 = np.abs(U) ** 2

    # PMNS-artige Winkel
    s13_sq = absU2[0, 2]
    c13_sq = 1.0 - s13_sq

    if c13_sq > 0:
        s12_sq = absU2[0, 1] / c13_sq
        s23_sq = absU2[1, 2] / c13_sq
    else:
        s12_sq = np.nan
        s23_sq = np.nan

    def safe_angle_deg(s2):
        if np.isnan(s2):
            return np.nan
        s2 = max(0.0, min(1.0, float(s2)))
        return math.degrees(math.asin(math.sqrt(s2)))

    theta13 = safe_angle_deg(s13_sq)
    theta12 = safe_angle_deg(s12_sq)
    theta23 = safe_angle_deg(s23_sq)

    # Jarlskog-artiger Indikator
    Jcp = np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0]))

    # Splittings
    d21 = E[1] - E[0]
    d31 = E[2] - E[0]
    d32 = E[2] - E[1]

    ratio_21_31 = abs(d21 / d31) if d31 != 0 else np.nan
    ratio_32_31 = abs(d32 / d31) if d31 != 0 else np.nan

    norm_singlet = np.linalg.norm(P_singlet, ord="fro")
    norm_delta = np.linalg.norm(Delta, ord="fro")
    norm_S = np.linalg.norm(S, ord="fro")
    norm_A = np.linalg.norm(A, ord="fro")

    relative_octet = norm_delta / norm_singlet if norm_singlet > 0 else np.nan
    chiral_fraction = norm_A / norm_delta if norm_delta > 0 else np.nan

    return {
        "P": P,
        "Delta": Delta,
        "S": S,
        "A": A,
        "H": H,
        "E": E,
        "U": U,
        "absU2": absU2,

        "theta12": theta12,
        "theta13": theta13,
        "theta23": theta23,

        "s12_sq": s12_sq,
        "s13_sq": s13_sq,
        "s23_sq": s23_sq,

        "Jcp": Jcp,

        "E1": E[0],
        "E2": E[1],
        "E3": E[2],

        "d21": d21,
        "d31": d31,
        "d32": d32,

        "ratio_21_31": ratio_21_31,
        "ratio_32_31": ratio_32_31,

        "norm_delta": norm_delta,
        "norm_S": norm_S,
        "norm_A": norm_A,
        "relative_octet": relative_octet,
        "chiral_fraction": chiral_fraction,

        "herm_error": herm_error,
        "unit_error": unit_error,
    }


# ------------------------------------------------------------
# Einzelnes N analysieren
# ------------------------------------------------------------

def analyze_N(N: int):
    t0 = time.time()

    starts = find_prime_quadruplet_starts(N)
    residues = starts % 210

    counts = Counter(int(x) for x in residues)

    main_mask = np.isin(residues, FAMILIES)
    family_sequence = residues[main_mask].astype(np.int64)

    M = transition_matrix_3(family_sequence)
    P = row_normalize(M)

    obs = neutrino_observables_from_P(P)

    elapsed = time.time() - t0

    absU2 = obs["absU2"]

    row = {
        "N": N,
        "quadruplets": len(starts),
        "anchor_5": counts.get(5, 0),
        "family_total": len(family_sequence),
        "family_11": counts.get(11, 0),
        "family_101": counts.get(101, 0),
        "family_191": counts.get(191, 0),

        "T_11_11": M[0, 0],
        "T_11_101": M[0, 1],
        "T_11_191": M[0, 2],
        "T_101_11": M[1, 0],
        "T_101_101": M[1, 1],
        "T_101_191": M[1, 2],
        "T_191_11": M[2, 0],
        "T_191_101": M[2, 1],
        "T_191_191": M[2, 2],

        "relative_octet": obs["relative_octet"],
        "norm_delta": obs["norm_delta"],
        "norm_S": obs["norm_S"],
        "norm_A": obs["norm_A"],
        "chiral_fraction": obs["chiral_fraction"],

        "theta12": obs["theta12"],
        "theta13": obs["theta13"],
        "theta23": obs["theta23"],

        "s12_sq": obs["s12_sq"],
        "s13_sq": obs["s13_sq"],
        "s23_sq": obs["s23_sq"],

        "Jcp": obs["Jcp"],

        "E1": obs["E1"],
        "E2": obs["E2"],
        "E3": obs["E3"],

        "d21": obs["d21"],
        "d31": obs["d31"],
        "d32": obs["d32"],
        "ratio_21_31": obs["ratio_21_31"],
        "ratio_32_31": obs["ratio_32_31"],

        "Ue1_sq": absU2[0, 0],
        "Ue2_sq": absU2[0, 1],
        "Ue3_sq": absU2[0, 2],
        "Umu1_sq": absU2[1, 0],
        "Umu2_sq": absU2[1, 1],
        "Umu3_sq": absU2[1, 2],
        "Utau1_sq": absU2[2, 0],
        "Utau2_sq": absU2[2, 1],
        "Utau3_sq": absU2[2, 2],

        "herm_error": obs["herm_error"],
        "unit_error": obs["unit_error"],
        "elapsed_sec": elapsed,
    }

    return row, M, P, obs


# ------------------------------------------------------------
# Ausgabehelfer
# ------------------------------------------------------------

def print_matrix(name, M, precision=5):
    print(f"\n{name}")
    with np.printoptions(precision=precision, suppress=True):
        print(M)


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    print("\n============================================================")
    print("BM mod-210 Neutrino-Skalierung")
    print("============================================================")
    print("Familien:", FAMILIES)
    print("N_VALUES:", N_VALUES)

    rows = []

    last = None

    for N in N_VALUES:
        print("\n------------------------------------------------------------")
        print(f"Analysiere N={N}")
        print("------------------------------------------------------------")

        row, M, P, obs = analyze_N(N)
        rows.append(row)
        last = (N, row, M, P, obs)

        print(f"Vierlinge:      {row['quadruplets']}")
        print(f"Familien total: {row['family_total']}")
        print(
            f"Familien:       11={row['family_11']}, "
            f"101={row['family_101']}, "
            f"191={row['family_191']}"
        )

        print(f"rel. Oktett:    {row['relative_octet']:.6f}")
        print(f"Chiralanteil:   {row['chiral_fraction']:.6f}")

        print(
            f"Winkel:         "
            f"θ12={row['theta12']:.3f}°, "
            f"θ13={row['theta13']:.3f}°, "
            f"θ23={row['theta23']:.3f}°"
        )

        print(f"J_CP:           {row['Jcp']:+.6f}")
        print(f"|Δ21/Δ31|:      {row['ratio_21_31']:.6f}")
        print(f"Hermitizität:   {row['herm_error']:.3e}")
        print(f"Unitarität:     {row['unit_error']:.3e}")
        print(f"Laufzeit:       {row['elapsed_sec']:.2f} s")

        pd.DataFrame(rows).to_csv("bm210_neutrino_scaling.csv", index=False)

    print("\n============================================================")
    print("Fertig")
    print("============================================================")
    print("CSV gespeichert: bm210_neutrino_scaling.csv")

    df = pd.DataFrame(rows)

    cols = [
        "N",
        "quadruplets",
        "relative_octet",
        "chiral_fraction",
        "theta12",
        "theta13",
        "theta23",
        "Jcp",
        "ratio_21_31",
        "Ue3_sq",
    ]

    print("\nKompaktübersicht:")
    print(df[cols].to_string(index=False))

    # Detailausgabe für größtes N
    N, row, M, P, obs = last

    print("\n============================================================")
    print(f"Detailausgabe für größtes N={N}")
    print("============================================================")

    print_matrix("Übergangsmatrix M", M, precision=0)
    print_matrix("Markov-Kern P", P, precision=6)
    print_matrix("Delta = P - J/3", obs["Delta"], precision=6)
    print_matrix("H = S + iA", obs["H"], precision=6)
    print_matrix("|U|^2", obs["absU2"], precision=6)

    print("\nInterpretation:")
    print("""
1. Stabiler kleiner theta13-Wert über N wäre der wichtigste BM-Neutrino-Hinweis.

2. Wenn theta12/theta23 stark schwanken, ist die PMNS-Analogie instabil.

3. Wenn |Δ21/Δ31| nicht deutlich kleiner wird, entsteht keine reale
   Neutrino-Massenhierarchie.

4. Wenn J_CP mit stabilem Vorzeichen bleibt, kann man von einem kleinen
   operatorischen CP-Analogon sprechen.

5. Wenn relative_octet gegen 0 fällt, verschwindet die Dynamik asymptotisch
   in reiner Gleichmischung.

6. Wenn relative_octet stabil > 0 bleibt, existiert ein kleiner, aber messbarer
   arithmetischer Oktettoperator.
""")


if __name__ == "__main__":
    main()