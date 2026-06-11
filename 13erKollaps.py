import sympy as sp
import numpy as np
from collections import Counter

print(
    "[13erKollaps.py] BM13-12er-Krone + H32-Slot-Diagnose "
    "(|H32|, |P|, |Q|, theta_slot, P, Q — wie H32_eabc.py)."
)


def _is_prime_fast(n: int) -> bool:
    """Primzahltest: zuerst Sage (bei SageMath-Python), sonst gmpy2, sonst SymPy."""
    n = int(n)
    if n < 2:
        return False
    try:
        from sage.arith.misc import is_prime as sage_is_prime

        return bool(sage_is_prime(n))
    except ImportError:
        pass
    try:
        import gmpy2

        return bool(gmpy2.is_prime(n))
    except ImportError:
        pass
    return bool(sp.isprime(n))


def eabc_class(n):
    r = n % 12
    if r == 1:
        return "E"
    if r == 5:
        return "A"
    if r == 7:
        return "B"
    if r == 11:
        return "C"
    return None


# Modulo-30-Hülle H32 (32 Slots), konsistent mit H32_eabc.py
H32_POS = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59]
H32 = sorted([-r for r in H32_POS] + H32_POS)


def h32_occupied_slots(M):
    """Offsets r ∈ H32 mit M+r prim (gleiche Logik wie occupied_prime_slots)."""
    M = int(M)
    return [r for r in H32 if _is_prime_fast(M + r)]


def h32_missing_slots(M):
    M = int(M)
    return [r for r in H32 if not _is_prime_fast(M + r)]


def h32_diagnose(M):
    P = h32_occupied_slots(M)
    Q = h32_missing_slots(M)
    return {
        "H_size": len(H32),
        "P": P,
        "Q": Q,
        "P_size": len(P),
        "Q_size": len(Q),
        "theta_slot": len(Q) / len(H32),
    }


def nearest_primes_around(M, k=12, R=500):
    """Die k Primzahlen mit kleinstem Abstand zu M im Radius R (ohne volles Fenster-scannen)."""
    M = int(M)
    primes = []
    for d in range(R + 1):
        cand = [M] if d == 0 else [M - d, M + d]
        for n in cand:
            if n >= 2 and _is_prime_fast(n):
                primes.append(n)
        if len(primes) >= k:
            break
    uniq = sorted(set(primes), key=lambda p: abs(p - M))
    nearest = uniq[:k]
    return sorted(nearest)

def eabc_balance_score(primes):
    counts = Counter(eabc_class(p) for p in primes)
    target = len(primes) / 4
    defect = sum((counts[x] - target)**2 for x in ["E","A","B","C"])
    max_defect = len(primes) ** 2
    return 1 - defect / max_defect, counts


def is_perfect_eabc(counts, k=12):
    target = k // 4
    return all(counts[x] == target for x in ["E", "A", "B", "C"])


def symmetry_score(M, primes):
    r = sorted([p-M for p in primes])
    if len(r) % 2 != 0:
        return 0
    half = len(r)//2
    defect = sum(abs(r[i] + r[-1-i]) for i in range(half))
    R = max(abs(x) for x in r) if r else 1
    return np.exp(-defect / (R + 1e-9)), defect, r


def max_radius(rel):
    return max(abs(x) for x in rel)


def is_perfect_sym(rel):
    rel = sorted(rel)
    return all(rel[i] + rel[-1 - i] == 0 for i in range(len(rel) // 2))


def transition_matrix(primes):
    classes = [eabc_class(p) for p in sorted(primes)]
    labels = ["E","A","B","C"]
    idx = {x:i for i,x in enumerate(labels)}
    K = np.zeros((4,4))
    for a,b in zip(classes[:-1], classes[1:]):
        if a in idx and b in idx:
            K[idx[a], idx[b]] += 1
    return K

def rank_score(K):
    if np.linalg.norm(K) == 0:
        return 0, []
    s = np.linalg.svd(K, compute_uv=False)
    power = np.sum(s**2)
    if power == 0:
        return 0, s
    dominant = s[0]**2 / power
    return 1 - dominant, s


def theta_rank(singular):
    """Anteil der nicht-dominanten Singulärwert-Energie (1 − σ₀²/Σσᵢ²)."""
    s = np.array(singular, dtype=float)
    if s.size == 0:
        return 0.0
    denom = np.sum(s**2)
    if denom == 0.0:
        return 0.0
    return 1.0 - (s[0] ** 2 / denom)


def chi_geom(N):
    return 1.0 - N ** (-1 / 3)


def sperre_score(M, small_primes=(2, 3, 5, 7, 11, 13)):
    score = 0
    for q in small_primes:
        if M % q == 0:
            score += 1/q
    return score

def cluster_score(M, k=12, R=500):
    primes = nearest_primes_around(M, k=k, R=R)
    if len(primes) < k:
        return None
    s_sym, defect, rel = symmetry_score(M, primes)
    s_eabc, counts = eabc_balance_score(primes)
    K = transition_matrix(primes)
    s_rank, sing = rank_score(K)
    s_sperre = sperre_score(M)
    
    total = (
        0.30*s_sym +
        0.30*s_eabc +
        0.25*s_rank +
        0.15*s_sperre
    )
    return {
        "M": M,
        "score": total,
        "sym": s_sym,
        "sym_defect": defect,
        "eabc": s_eabc,
        "counts": counts,
        "rank": s_rank,
        "singular": sing,
        "sperre": s_sperre,
        "primes": primes,
        "rel": rel,
    }

candidates = []
for M in range(1000, 200000):
    if _is_prime_fast(M):
        continue
    res = cluster_score(M, k=12, R=1000)
    if res:
        candidates.append(res)

top = sorted(candidates, key=lambda x: x["score"], reverse=True)[:20]

chi13 = chi_geom(13)
theta_ideal = 4 / 7

for c in top:
    th = theta_rank(c["singular"])
    delta_47 = abs(th - theta_ideal)
    h32 = h32_diagnose(c["M"])
    print(
        c["M"],
        "score", round(c["score"], 6),
        "theta", round(th, 6),
        "delta13", round(abs(th - chi13), 6),
        "delta47", round(delta_47, 8),
        "Rmax", max_radius(c["rel"]),
        "perfectEABC", is_perfect_eabc(c["counts"]),
        "perfectSym", is_perfect_sym(c["rel"]),
        "|H32|",
        h32["H_size"],
        "|P|",
        h32["P_size"],
        "|Q|",
        h32["Q_size"],
        "theta_slot",
        h32["theta_slot"],
        "P",
        h32["P"],
        "Q",
        h32["Q"],
        "counts",
        c["counts"],
        "rel",
        c["rel"],
        "sing",
        np.array2string(c["singular"], precision=8),
    )

perfect = []
for c in candidates:
    theta = theta_rank(c["singular"])
    delta47 = abs(theta - theta_ideal)
    if (
        is_perfect_eabc(c["counts"])
        and is_perfect_sym(c["rel"])
        and delta47 < 1e-6
    ):
        perfect.append(c)

perfect = sorted(perfect, key=lambda x: (x["M"], max_radius(x["rel"])))

for c in perfect[:20]:
    theta = theta_rank(c["singular"])
    h32 = h32_diagnose(c["M"])
    print(
        c["M"],
        "Rmax",
        max_radius(c["rel"]),
        "theta",
        theta,
        "|H32|",
        h32["H_size"],
        "|P|",
        h32["P_size"],
        "|Q|",
        h32["Q_size"],
        "theta_slot",
        h32["theta_slot"],
        "P",
        h32["P"],
        "Q",
        h32["Q"],
        "rel",
        c["rel"],
        "sing",
        c["singular"],
    )