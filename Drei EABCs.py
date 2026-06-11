import numpy as np

# =========================================
# 1. PRIME SIEVE (fast)
# =========================================
N = 10_000_000

is_prime = np.ones(N+10, dtype=bool)
is_prime[:2] = False

for i in range(2, int(N**0.5)+1):
    if is_prime[i]:
        is_prime[i*i:N+1:i] = False

# =========================================
# 2. EABC CLASS
# =========================================
def eabc_class(p):
    r = p % 12
    if r == 1: return "E"
    if r == 5: return "A"
    if r == 7: return "B"
    if r == 11: return "C"
    return None

# =========================================
# 3. PHASE (mod 36)
# =========================================
def phase(p):
    return (p % 36) // 12   # 0,1,2

# =========================================
# 4. FIND QUADRUPLETS
# =========================================
quadruplets = []
phases = {0: [], 1: [], 2: []}

for p in range(5, N-10):
    if not is_prime[p]:
        continue

    if is_prime[p+2] and is_prime[p+6] and is_prime[p+8]:
        quadruplets.append(p)
        ph = phase(p)
        phases[ph].append(p)

# =========================================
# 5. BASIC STATS
# =========================================
total = len(quadruplets)

print("Total quadruplets:", total)
print()

for ph in [0,1,2]:
    count = len(phases[ph])
    print(f"Phase {ph}: {count} ({count/total:.4f})")

# =========================================
# 6. GAP ANALYSIS
# =========================================
def gaps(arr):
    return np.diff(arr) if len(arr) > 1 else []

print("\nMean gaps:")
for ph in [0,1,2]:
    g = gaps(phases[ph])
    if len(g) > 0:
        print(f"Phase {ph}: mean gap = {np.mean(g):.2f}")

# =========================================
# 7. CLUSTER TEST (window density)
# =========================================
W = 10000
def cluster_score(arr):
    arr = np.array(arr)
    counts = []
    for i in range(0, N, W):
        counts.append(np.sum((arr >= i) & (arr < i+W)))
    return np.var(counts)

print("\nCluster variance:")
for ph in [0,1,2]:
    print(f"Phase {ph}: var = {cluster_score(phases[ph]):.4f}")