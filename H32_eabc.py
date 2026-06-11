# SageMath: EABC + Modulo-30-Hülle H32 (32 Slots), Übergangsmatrix, Spektraldaten
#
# Aufruf:
#   sage "H32_eabc.py"              # search_BSN32 bis 10^7
#   sage -python "H32_eabc.py"
# In Sage: load("H32_eabc.py"); test_H32(113160)  für Einzel-Demo

from collections import Counter

import numpy as np

try:
    from sage.all import factor, is_prime  # type: ignore[import-untyped]
except ImportError:
    raise SystemExit(
        "SageMath wird benötigt (sage.all.is_prime, factor).\n"
        'Bitte z.B. ausführen: sage "H32_eabc.py"'
    ) from None


# EABC-Klassen modulo 12
def eabc_class(n):
    r = int(n) % 12
    if r == 1:
        return "E"
    if r == 5:
        return "A"
    if r == 7:
        return "B"
    if r == 11:
        return "C"
    return None


# kanonische zweite reduzierte Modulo-30-Hülle: 32 Slots
H32_pos = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59]
H32 = sorted([-r for r in H32_pos] + H32_pos)


def occupied_prime_slots(M, H):
    return [r for r in H if is_prime(M + r)]


def missing_slots(M, H):
    return [r for r in H if not is_prime(M + r)]


def eabc_counts(M, R):
    return Counter(eabc_class(M + r) for r in R)


def transition_matrix(M, R):
    labels = ["E", "A", "B", "C"]
    idx = {x: i for i, x in enumerate(labels)}
    R_sorted = sorted(R)
    classes = [eabc_class(M + r) for r in R_sorted]
    K = np.zeros((4, 4))
    for a, b in zip(classes[:-1], classes[1:]):
        if a in idx and b in idx:
            K[idx[a], idx[b]] += 1
    return K, classes


def spectral_data(K):
    s = np.linalg.svd(K, compute_uv=False)
    mass = np.sum(s**2)
    lambdas = (s**2) / mass if mass > 0 else np.zeros_like(s)
    theta = 1 - lambdas[0] if mass > 0 else None
    entropy = -sum(float(l) * float(np.log(l)) for l in lambdas if l > 0)
    entropy_bits = entropy / np.log(2)
    entropy_norm = entropy / np.log(4)
    return s, lambdas, theta, entropy, entropy_bits, entropy_norm


def is_symmetric(R):
    S = set(R)
    return all(-r in S for r in S)


def is_eabc_balanced_12(M, P):
    if len(P) != 12:
        return False
    c = eabc_counts(M, P)
    return all(c[x] == 3 for x in ["E", "A", "B", "C"])


def search_BSN32(limit=10_000_000):
    hits = []
    for M in range(30, limit + 1, 30):
        if is_prime(M):
            continue

        P = occupied_prime_slots(M, H32)
        Q = missing_slots(M, H32)

        if len(P) != 12:
            continue

        if not is_symmetric(P):
            continue

        if not is_eabc_balanced_12(M, P):
            continue

        K, _ = transition_matrix(M, P)
        s, lambdas, theta, entropy, entropy_bits, entropy_norm = spectral_data(K)

        hits.append(
            {
                "M": M,
                "P": P,
                "Q": Q,
                "counts": eabc_counts(M, P),
                "singular": s,
                "lambda": lambdas,
                "theta_rank": theta,
                "S_ent": entropy,
                "S_ent_bits": entropy_bits,
                "S_ent_norm": entropy_norm,
                "theta_slot": len(Q) / len(H32),
                "factor_M": factor(M),
            }
        )

        print(
            "HIT",
            M,
            "theta_rank",
            theta,
            "S_ent",
            entropy,
            "P",
            P,
            "Q",
            Q,
            "sing",
            s,
        )

    return hits


def test_H32(M):
    P = occupied_prime_slots(M, H32)
    Q = missing_slots(M, H32)
    K, classes = transition_matrix(M, P)
    s, lambdas, theta, entropy, entropy_bits, entropy_norm = spectral_data(K)

    return {
        "M": M,
        "factor_M": factor(M),
        "H_size": len(H32),
        "P_size": len(P),
        "Q_size": len(Q),
        "P": P,
        "Q": Q,
        "counts": eabc_counts(M, P),
        "classes": classes,
        "K": K,
        "singular": s,
        "lambda": lambdas,
        "theta_rank": theta,
        "S_ent": entropy,
        "S_ent_bits": entropy_bits,
        "S_ent_norm": entropy_norm,
        "theta_slot": len(Q) / len(H32),
    }


if __name__ == "__main__":
    print("[H32_eabc.py] BSN32-Suche auf Vielfachen von 30 — nicht 13erKollaps (BM13-12-Krone).")
    hits32 = search_BSN32(10_000_000)
    print("Treffer:", len(hits32))
