#!/usr/bin/env python3

import numpy as np
from collections import Counter
from math import isqrt
from numpy.linalg import eigvals

# ==========================================
# Parameter
# ==========================================

N = 100_000_000

HALF_CHANNELS = [11, 101, 191, 221, 311, 401]
INDEX = {r: i for i, r in enumerate(HALF_CHANNELS)}

# ==========================================
# Primsieb
# ==========================================

def prime_sieve(limit):
    sieve = np.ones(limit + 1, dtype=np.bool_)
    sieve[:2] = False

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p::p] = False

    return sieve


# ==========================================
# Vierlinge
# ==========================================

def quadruplets(limit):

    sieve = prime_sieve(limit + 8)

    starts = []

    for p in range(5, limit + 1):

        if (
            sieve[p]
            and sieve[p + 2]
            and sieve[p + 6]
            and sieve[p + 8]
        ):
            starts.append(p)

    return np.array(starts, dtype=np.int64)


# ==========================================
# Kanalbestimmung
# ==========================================

def channel_mod420(p):

    r = p % 420

    if r in INDEX:
        return r

    return None


# ==========================================
# Übergangsmatrix
# ==========================================

def build_transition_matrix(channels):

    M = np.zeros((6, 6), dtype=np.float64)

    for a, b in zip(channels[:-1], channels[1:]):

        ia = INDEX[a]
        ib = INDEX[b]

        M[ia, ib] += 1

    row_sums = M.sum(axis=1)

    for i in range(6):
        if row_sums[i] > 0:
            M[i] /= row_sums[i]

    return M


# ==========================================
# Dirac-Operator
# ==========================================

def dirac_operator(T):

    L = np.eye(6) - T

    D = np.block([
        [np.zeros((6, 6)), L],
        [L.T, np.zeros((6, 6))]
    ])

    return D


# ==========================================
# Fourier auf Ring
# ==========================================

def ring_fourier(counts):

    v = np.array(counts, dtype=float)

    v /= v.sum()

    return np.fft.fft(v)


# ==========================================
# Hauptprogramm
# ==========================================

def main():

    print("Suche Vierlinge ...")

    starts = quadruplets(N)

    print("Vierlinge:", len(starts))

    channels = []

    for p in starts:

        r = channel_mod420(p)

        if r is not None:
            channels.append(r)

    channels = np.array(channels)

    print("Kanäle:", len(channels))

    # ------------------------
    # Häufigkeiten
    # ------------------------

    cnt = Counter(channels)

    print("\nHäufigkeiten")

    for c in HALF_CHANNELS:
        print(c, cnt[c])

    # ------------------------
    # Übergangsmatrix
    # ------------------------

    T = build_transition_matrix(channels)

    print("\nTransition Matrix")
    np.set_printoptions(
        precision=4,
        suppress=True,
        linewidth=200
    )
    print(T)

    # ------------------------
    # Eigenwerte
    # ------------------------

    eigT = eigvals(T)

    eigT = eigT[np.argsort(-np.abs(eigT))]

    print("\nEigenwerte(T)")
    print(eigT)

    # Spektrallücke

    gap = 1.0 - np.abs(eigT[1])

    print("\nSpektrallücke")
    print(gap)

    # ------------------------
    # Dirac
    # ------------------------

    D = dirac_operator(T)

    eigD = np.sort(np.real(eigvals(D)))

    print("\nKleinste Dirac Eigenwerte")
    print(eigD[:20])

    # ------------------------
    # Fourier
    # ------------------------

    counts = [cnt[c] for c in HALF_CHANNELS]

    F = ring_fourier(counts)

    print("\nFourier Moden")

    for k, val in enumerate(F):

        print(
            k,
            abs(val),
            np.angle(val)
        )

    # ------------------------
    # Chiraler Modus
    # ------------------------

    chi = np.array([
        +1,
        -1,
        +1,
        -1,
        +1,
        -1
    ], dtype=float)

    chi /= np.linalg.norm(chi)

    chi2 = T @ chi

    overlap = np.dot(chi, chi2)

    print("\nChiraler Overlap")
    print(overlap)

    print("\nFertig.")


if __name__ == "__main__":
    main()