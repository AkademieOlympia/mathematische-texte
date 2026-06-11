import numpy as np
from collections import defaultdict


def build_transition_matrix(states_sequence):
    # eindeutige Zustände
    states = sorted(list(set(states_sequence)))
    index = {s: i for i, s in enumerate(states)}
    k = len(states)

    P = np.zeros((k, k))

    for i in range(len(states_sequence) - 1):
        a = states_sequence[i]
        b = states_sequence[i + 1]
        P[index[a], index[b]] += 1

    # Normierung
    for i in range(k):
        row_sum = P[i].sum()
        if row_sum > 0:
            P[i] /= row_sum

    return P, states


def spectral_analysis(P):
    eigvals, eigvecs = np.linalg.eig(P)

    idx = np.argsort(-np.abs(eigvals))
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    return eigvals, eigvecs


def stationary_distribution(P):
    w, v = np.linalg.eig(P.T)
    i = np.argmin(np.abs(w - 1))
    pi = np.real(v[:, i])

    if pi.sum() < 0:
        pi = -pi

    pi /= pi.sum()
    return pi


def run_markov_pipeline(states_sequence):
    P, states = build_transition_matrix(states_sequence)

    eigvals, eigvecs = spectral_analysis(P)
    pi = stationary_distribution(P)

    print("\n=== Zustaende ===")
    for s in states:
        print(s)

    print("\n=== Uebergangsmatrix ===")
    print(P)

    print("\n=== Fuehrende Eigenwerte ===")
    for i in range(min(10, len(eigvals))):
        print(i, eigvals[i])

    print("\n=== Stationaere Verteilung ===")
    for s, p in zip(states, pi):
        print(f"{str(s):10s}: {p:.4f}")

    return P, eigvals, eigvecs, pi


if __name__ == "__main__":
    states_sequence = ["A", "B", "A", "C", "B", "B", "A", "D", "A", "C"]
    run_markov_pipeline(states_sequence)