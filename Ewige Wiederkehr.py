import numpy as np
import random
import bisect
import math

N = 30_000_000
MOD = 420
W = 10000
TARGET_HOT = 401
TARGET_COLD = 311
PERMS = 2000

# ------------------------------------------------------------
# Sieb
# ------------------------------------------------------------
is_prime = np.ones(N + 20, dtype=bool)
is_prime[:2] = False

for i in range(2, int(N**0.5) + 1):
    if is_prime[i]:
        is_prime[i*i:N+20:i] = False

# ------------------------------------------------------------
# Vierlinge
# ------------------------------------------------------------
quadruplets = []
for p in range(5, N - 10):
    if is_prime[p] and is_prime[p+2] and is_prime[p+6] and is_prime[p+8]:
        quadruplets.append(p)

quadruplets = np.array(quadruplets)
labels = np.array([p % MOD for p in quadruplets])

print("N:", N)
print("Anzahl Vierlinge:", len(quadruplets))
print("MOD:", MOD)
print("W:", W)

# ------------------------------------------------------------
# Hilfsfunktionen
# ------------------------------------------------------------
def delta_for_class(qs, N, W, trials=10000):
    """
    qs = Positionen der Vierlinge in einer Klasse.
    Gibt Cluster, Random, Delta zurück.
    """
    qs = sorted(qs)
    if len(qs) < 3:
        return None, None, None

    gaps = [qs[i+1] - qs[i] for i in range(len(qs)-1)]

    cluster = sum(1 for g in gaps if g < W) / len(gaps)

    count = 0
    for _ in range(trials):
        n = random.randint(0, N - W)
        j = bisect.bisect_left(qs, n)
        if j < len(qs) and qs[j] <= n + W:
            count += 1

    rand = count / trials
    delta = cluster - rand

    return cluster, rand, delta


def score_for_labels(quadruplets, labels, hot=401, cold=311):
    """
    Berechnet S = Delta_hot - Delta_cold.
    """
    q_hot = quadruplets[labels == hot]
    q_cold = quadruplets[labels == cold]

    c_hot, r_hot, d_hot = delta_for_class(q_hot, N, W)
    c_cold, r_cold, d_cold = delta_for_class(q_cold, N, W)

    return {
        "hot_count": len(q_hot),
        "cold_count": len(q_cold),
        "hot_cluster": c_hot,
        "hot_random": r_hot,
        "hot_delta": d_hot,
        "cold_cluster": c_cold,
        "cold_random": r_cold,
        "cold_delta": d_cold,
        "S": d_hot - d_cold
    }

# ------------------------------------------------------------
# Echter Wert
# ------------------------------------------------------------
real = score_for_labels(quadruplets, labels, TARGET_HOT, TARGET_COLD)

print("\n--- Echter Befund ---")
print("Hot class:", TARGET_HOT)
print("Cold class:", TARGET_COLD)
print("Hot count:", real["hot_count"])
print("Cold count:", real["cold_count"])
print("Delta hot:", real["hot_delta"])
print("Delta cold:", real["cold_delta"])
print("S = Delta_hot - Delta_cold:", real["S"])

# ------------------------------------------------------------
# Permutationstest
# ------------------------------------------------------------
perm_scores = []

for k in range(PERMS):
    permuted_labels = np.random.permutation(labels)
    res = score_for_labels(quadruplets, permuted_labels, TARGET_HOT, TARGET_COLD)
    perm_scores.append(res["S"])

perm_scores = np.array(perm_scores)

p_value = np.mean(perm_scores >= real["S"])
mean_perm = np.mean(perm_scores)
std_perm = np.std(perm_scores)

z_perm = (real["S"] - mean_perm) / std_perm if std_perm > 0 else float("nan")

print("\n--- Permutationstest ---")
print("Permutationen:", PERMS)
print("Mittelwert S_perm:", mean_perm)
print("Std S_perm:", std_perm)
print("Reales S:", real["S"])
print("Permutation-z:", z_perm)
print("p-Wert:", p_value)

# Optional: Quantile
print("\nQuantile S_perm:")
for q in [0.5, 0.9, 0.95, 0.975, 0.99]:
    print(q, np.quantile(perm_scores, q))