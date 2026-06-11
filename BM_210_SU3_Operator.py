#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BM_210_SU3_Operator.py

Analogon: drei mod-210-Familien (11, 101, 191) wie Farbtripel eines SU(3)-Kopplungsmusters.

Aus gleitenden Quadruplett-Starts berechnet:
- Roh-Übergangsmatrix 3×3 (nur aufeinanderfolgende Hauptfamilien)
- Zeilen-stochastischer Markov-Kern und Spektrum
- näherungsweise Stationärverteilung (Power-Iteration)
- Erwartungen der diagonalen SU(3)-Köpfchen λ₃, λ₈ (Standard-Gell-Mann auf dem
  empirischen Wahrscheinlichkeitsvektor in der Reihenfolge |11|, |101|, |191⟩)

Die Datei war mit Terminalausgabe von BM_210_Skalierung.py überschrieben worden;
dies ist die gültige Ersatz-Implementierung.
"""

from __future__ import annotations

import math
import sys
import numpy as np
from numpy.linalg import eigvals

MOD = 210
FAMILIES = np.array([11, 101, 191])
IDX = {v: i for i, v in enumerate(FAMILIES)}


def sieve_bool(n: int) -> np.ndarray:
    is_prime = np.ones(n + 1, dtype=bool)
    is_prime[:2] = False
    lim = int(n**0.5)
    for k in range(2, lim + 1):
        if is_prime[k]:
            is_prime[k * k : n + 1 : k] = False
    return is_prime


def find_prime_quadruplet_starts(N: int) -> np.ndarray:
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


def transition_matrix_3(families: np.ndarray) -> np.ndarray:
    """Zählt Übergänge i→j in Reihenfolge der Familienliste (nur Hauptfamilien)."""
    M = np.zeros((3, 3), dtype=np.int64)
    for a, b in zip(families[:-1], families[1:]):
        a_i, b_j = int(a), int(b)
        if a_i in IDX and b_j in IDX:
            M[IDX[a_i], IDX[b_j]] += 1
    return M


def row_stochastic(M: np.ndarray) -> np.ndarray:
    """Normiert Zeilen; leere Zeilen → gleichverteilt (Fallback)."""
    M = np.asarray(M, dtype=float)
    P = np.zeros((3, 3), dtype=float)
    for i in range(3):
        s = float(np.sum(M[i]))
        if s > 0.0:
            P[i] = M[i] / s
        else:
            P[i] = 1.0 / 3.0
    return P


def stationary_row(P: np.ndarray, iterations: int = 8000, tol: float = 1e-12) -> np.ndarray:
    """Annähernd stationäre Zeilenverteilung π mit π ≈ π P für zeilen-stochastisches P."""
    pi = np.array([1.0, 1.0, 1.0], dtype=float) / 3.0
    for _ in range(iterations):
        nxt = pi @ P
        if float(np.linalg.norm(nxt - pi, ord=1)) < tol:
            break
        pi = nxt
    pi = pi / float(np.sum(np.clip(pi, 0.0, None)) + 1e-15)
    return pi


def gellmann_expectations(p: np.ndarray) -> tuple[float, float]:
    """⟨λ₃⟩, ⟨λ₈⟩ für Diagonal-ρ = diag(p) in Basis (11,101,191). Standard-Norm SU(3)."""
    # λ₃ = diag(1, -1, 0)
    lam3_expect = float(p[0] - p[1])
    # λ₈ = (1/sqrt(3))*diag(1, 1, -2)
    inv_sqrt3 = 1.0 / math.sqrt(3.0)
    lam8_expect = float(inv_sqrt3 * (p[0] + p[1] - 2.0 * p[2]))
    return lam3_expect, lam8_expect


def analyze(N: int) -> None:
    starts = find_prime_quadruplet_starts(N)
    residues = (starts % MOD).astype(np.int64)
    mask = np.isin(residues, FAMILIES)
    families = residues[mask]

    print("\n" + "=" * 60)
    print("BM mod-210 – SU(3)-Operator-Analog")
    print("=" * 60)
    print(f"N = {N}")
    print(f"Quadruppelstarts: {len(starts)}")
    print(f"Sequenzlänge (nur 11|101|191): {len(families)}")
    if len(families) < 5:
        print("Zu wenige Hauptfamilien-Einträge für sinnvolle Statistik.")
        return

    # relative Häufigkeiten ohne Reihenfolge
    fc = np.array(
        [
            np.count_nonzero(families == FAMILIES[0]),
            np.count_nonzero(families == FAMILIES[1]),
            np.count_nonzero(families == FAMILIES[2]),
        ],
        dtype=float,
    )
    fc_p = fc / fc.sum()

    M = transition_matrix_3(families)
    print("\nRoh-Zählung Übergänge M[i→j], Zeilen/Spalten =", FAMILIES.tolist())
    print(M)
    print("Zeilensummen:", M.sum(axis=1).tolist())

    P = row_stochastic(M)
    lam = eigvals(P)
    order = np.argsort(-np.real(lam))
    lam = lam[order]
    pi = stationary_row(P)

    print("\nZeilen-Markov-Kern P (Leerzeile → Gleichverteilung):")
    print(np.round(P, 6))

    print("\nEigenwerte von P:")
    for k in range(3):
        z = lam[k]
        print(f"  λ_{k+1} = {z.real:.6f} {'+' if z.imag >= 0 else '-'} {abs(z.imag):.6f}j")

    print("\nnäherungsweise Stationärverteilung π (π ≈ π P):")
    for i in range(3):
        print(f"  P({int(FAMILIES[i])}) ≈ {pi[i]:.6f}")

    l3_hat, l8_hat = gellmann_expectations(fc_p)
    l3_sta, l8_sta = gellmann_expectations(pi)
    print("\nGell-Mann-diagonale Köpfchen auf empirischer Verteilung p (ohne Reihenfolge):")
    print(f"  ⟨λ₃⟩̂ = {l3_hat:+.6f}   ⟨λ₈⟩̂ = {l8_hat:+.6f}")
    print("\nDieselben auf π aus dem Markov-Kern:")
    print(f"  ⟨λ₃⟩π = {l3_sta:+.6f}   ⟨λ₈⟩π = {l8_sta:+.6f}")
    print(
        "(Vollständig gemischtes SU(3)-Singulett hätte λ₃ = λ₈ = 0; "
        "reine Gleichverteilung |pᵢ| = 1/3 ebenfalls.)"
    )

    J = np.ones((3, 3))
    P_uniform = J / 3.0

    Delta = P - P_uniform

    norm_uniform = float(np.linalg.norm(P_uniform, ord="fro"))
    norm_delta = float(np.linalg.norm(Delta, ord="fro"))
    relative_delta = norm_delta / norm_uniform

    print("\n‖P - J/3‖_F =", norm_delta)
    print("‖J/3‖_F     =", norm_uniform)
    print("relative Oktettstärke =", relative_delta)


def main() -> None:
    N = int(sys.argv[1]) if len(sys.argv) >= 2 else 10_000_000
    analyze(N)


if __name__ == "__main__":
    main()
