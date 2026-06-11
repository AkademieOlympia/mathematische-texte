#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BM_210_Neutrino_Test.py

Neutrino-Analogtest für den mod-210-Familienoperator der Primzahlvierlinge.

Drei Hauptfamilien:
    F1 = 11
    F2 = 101
    F3 = 191

Interpretation als Flavorbasis:
    |11>  ~ |nu_e>
    |101> ~ |nu_mu>
    |191> ~ |nu_tau>

Aus Übergangskern P:
    P = J/3 + Delta

Wir bauen einen hermiteschen Hamiltonian:
    H = Sym(Delta) + i Anti(Delta)

Dann:
    H |m_i> = E_i |m_i>

Die Eigenvektoren bilden eine PMNS-artige Mischungsmatrix U.
Aus U extrahieren wir:
    theta_12, theta_13, theta_23
    Jarlskog-artigen CP-Wert
    Oszillationswahrscheinlichkeiten P(alpha -> beta)
"""

import numpy as np
import math


# ------------------------------------------------------------
# 1. Dein empirischer Markov-Kern P aus N=10^7
# ------------------------------------------------------------

P = np.array([
    [0.317610, 0.352201, 0.330189],
    [0.386760, 0.320557, 0.292683],
    [0.359589, 0.287671, 0.352740],
], dtype=float)

families = ["11", "101", "191"]
flavors = [r"$\nu_e$", r"$\nu_\mu$", r"$\nu_\tau$"]


# ------------------------------------------------------------
# 2. SU(3)-Singulett/Oktett-Zerlegung
# ------------------------------------------------------------

J = np.ones((3, 3), dtype=float)
P_singlet = J / 3.0
Delta = P - P_singlet

S = 0.5 * (Delta + Delta.T)
A = 0.5 * (Delta - Delta.T)

# Hermitescher Neutrino-Analog-Hamiltonian
# A ist reell antisymmetrisch, also ist iA hermitesch.
H = S + 1j * A

# Optional: Spur entfernen
H = H - np.trace(H) / 3.0 * np.eye(3)


# ------------------------------------------------------------
# 3. Diagonalisierung: Massbasis-Analog
# ------------------------------------------------------------

eigvals, eigvecs = np.linalg.eigh(H)

# Sortieren nach Eigenwert
idx = np.argsort(eigvals)
eigvals = eigvals[idx]
eigvecs = eigvecs[:, idx]

# PMNS-artige Mischungsmatrix
# Konvention: flavor = U * mass
U = eigvecs

# Phasen bereinigen: Jede Spalte so drehen, dass größter Eintrag reell positiv ist.
for j in range(3):
    k = np.argmax(np.abs(U[:, j]))
    phase = np.angle(U[k, j])
    U[:, j] *= np.exp(-1j * phase)


# ------------------------------------------------------------
# 4. Mischungswinkel aus |U|
# Standard-PMNS-Konvention näherungsweise:
#
# s13^2 = |U_e3|^2
# s12^2 = |U_e2|^2 / (1 - |U_e3|^2)
# s23^2 = |U_mu3|^2 / (1 - |U_e3|^2)
# ------------------------------------------------------------

absU2 = np.abs(U) ** 2

s13_sq = absU2[0, 2]
c13_sq = 1.0 - s13_sq

s12_sq = absU2[0, 1] / c13_sq if c13_sq > 0 else np.nan
s23_sq = absU2[1, 2] / c13_sq if c13_sq > 0 else np.nan

theta13 = math.degrees(math.asin(math.sqrt(max(0.0, min(1.0, s13_sq)))))
theta12 = math.degrees(math.asin(math.sqrt(max(0.0, min(1.0, s12_sq)))))
theta23 = math.degrees(math.asin(math.sqrt(max(0.0, min(1.0, s23_sq)))))


# ------------------------------------------------------------
# 5. Jarlskog-artiger CP-Indikator
#
# J_CP = Im(U_e1 U_mu2 U_e2* U_mu1*)
# ------------------------------------------------------------

Jcp = np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0]))


# ------------------------------------------------------------
# 6. Massensplittings-Analog
#
# E_i sind dimensionslos. Für ein echtes Neutrino-Modell bräuchte man
# eine physikalische Skala. Hier betrachten wir nur Verhältnisse.
# ------------------------------------------------------------

E = eigvals.real

d21 = E[1] - E[0]
d31 = E[2] - E[0]
d32 = E[2] - E[1]

ratio_21_31 = abs(d21 / d31) if d31 != 0 else np.nan


# ------------------------------------------------------------
# 7. Oszillationswahrscheinlichkeiten
#
# Analog:
#     A(alpha -> beta; t) = sum_i U_beta,i exp(-i E_i t) U_alpha,i*
#     P = |A|^2
# ------------------------------------------------------------

def oscillation_matrix(t):
    phase = np.diag(np.exp(-1j * E * t))
    Aamp = U @ phase @ np.conjugate(U.T)
    Prob = np.abs(Aamp) ** 2
    return Prob


def print_matrix(name, M, precision=5):
    print(f"\n{name}")
    with np.printoptions(precision=precision, suppress=True):
        print(M)


# ------------------------------------------------------------
# 8. Ausgabe
# ------------------------------------------------------------

print("\n============================================================")
print("BM mod-210 Neutrino-Analogtest")
print("============================================================")

print("\nFlavorbasis:")
for f, nu in zip(families, flavors):
    print(f"  Familie {f}  ~  {nu}")

print_matrix("P Markov-Kern", P)
print_matrix("P_singlet = J/3", P_singlet)
print_matrix("Delta = P - J/3", Delta)
print_matrix("S = symmetrischer Anteil", S)
print_matrix("A = antisymmetrischer Anteil", A)
print_matrix("H = S + iA, spurfrei", H)

print("\nHermitizitätscheck:")
print("||H - H†||_F =", np.linalg.norm(H - np.conjugate(H.T), "fro"))

print("\nNormen:")
print(f"||J/3||_F        = {np.linalg.norm(P_singlet, 'fro'):.8f}")
print(f"||Delta||_F      = {np.linalg.norm(Delta, 'fro'):.8f}")
print(f"relative Oktettstärke = {np.linalg.norm(Delta, 'fro') / np.linalg.norm(P_singlet, 'fro'):.8f}")
print(f"||S||_F          = {np.linalg.norm(S, 'fro'):.8f}")
print(f"||A||_F          = {np.linalg.norm(A, 'fro'):.8f}")
print(f"Chiralitätsanteil ||A||/||Delta|| = {np.linalg.norm(A, 'fro') / np.linalg.norm(Delta, 'fro'):.8f}")

print("\nEigenwerte des H-Neutrino-Analogons:")
for i, x in enumerate(E, start=1):
    print(f"  E_{i} = {x:+.8f}")

print("\nMassensplitting-Analog:")
print(f"  Δ21 = {d21:+.8f}")
print(f"  Δ31 = {d31:+.8f}")
print(f"  Δ32 = {d32:+.8f}")
print(f"  |Δ21/Δ31| = {ratio_21_31:.8f}")

print_matrix("|U|^2, PMNS-artige Mischungsmatrix-Beträge", absU2)

print("\nUnitaritätscheck U:")
print(
    "||U†U - I||_F =",
    np.linalg.norm(np.conjugate(U.T) @ U - np.eye(3), "fro"),
)

print("\nPMNS-artige Mischungswinkel:")
print(f"  theta_12 ≈ {theta12:.3f}°")
print(f"  theta_13 ≈ {theta13:.3f}°")
print(f"  theta_23 ≈ {theta23:.3f}°")

print("\nJarlskog-artiger CP-Indikator:")
print(f"  J_CP ≈ {Jcp:+.8e}")

print("\nOszillationswahrscheinlichkeiten P(alpha -> beta; t):")

for t in [0, 5, 10, 20, 50, 100]:
    Prob = oscillation_matrix(t)
    print_matrix(f"t = {t}", Prob, precision=4)
    print("  Zeilensummen:", np.sum(Prob, axis=1))

print("\nInterpretation:")
print("""
1. Wenn |U|^2 nahe an der Einheitsmatrix liegt:
   kaum Flavor-Mischung.

2. Wenn |U|^2 breit verteilt ist:
   starke PMNS-artige Mischung.

3. Wenn theta_13 klein, theta_12/theta_23 groß:
   qualitativ neutrinoähnlich.

4. Wenn J_CP nahe 0:
   keine CP-artige Phase im arithmetischen Operator.

5. Wenn |Δ21/Δ31| deutlich kleiner als 1:
   hierarchische Massensplittings analog zu realen Neutrinos.

6. Da H aus dem mod-210-Operator stammt, sind alle Werte dimensionslos.
   Für echte Physik bräuchte man eine externe Skala.
""")