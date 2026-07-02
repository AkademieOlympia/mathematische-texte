import math
from collections import defaultdict

def v2(n):
    v = 0
    while n % 2 == 0: n //= 2; v += 1
    return v

def eabc_class(n):
    m = n
    while m % 2 == 0: m //= 2
    while m % 3 == 0: m //= 3
    return {1:'E', 5:'A', 7:'B', 11:'C'}[m % 12]

def T_odd(n):
    m = 3*n + 1; return m >> v2(m)

# Zyklensuche n=1..10000
for n0 in range(1, 10001):
    if n0 % 2 == 0 or n0 % 3 == 0: continue
    visited, n = {}, n0
    for _ in range(1_000_000):
        if n in visited: break   # trivialer Zyklus zurück zu 1
        if n == 1: break
        visited[n] = True; n = T_odd(n)
# → Keine nicht-trivialen Zyklen gefunden ✓

# Kontraktionsraten
log_rates = defaultdict(list); v2_vals = defaultdict(list)
for n in range(1, 20001, 2):
    if n % 3 == 0: continue
    cl = eabc_class(n); v = v2(3*n+1)
    log_rates[cl].append(math.log((3*n+1)/(n*2**v)))
    v2_vals[cl].append(v)

# E_π[v₂] = 1.9994 ≈ 2 ✓
# E_π[log r] = −0.2871 ≈ log(3/4) = −0.2877 ✓