# Referenzen (gut zum Nachschlagen):
# https://en.wikipedia.org/wiki/Regular_icosahedron
# https://en.wikipedia.org/wiki/Alternating_group_A5

import itertools, numpy as np, math, collections

phi = (1 + 5**0.5) / 2

# 12 Ikosaeder-Ecken in Standardlage
verts = []
for s1 in (-1,1):
    for s2 in (-1,1):
        verts.append((0, s1, s2*phi))
        verts.append((s1, s2*phi, 0))
        verts.append((s1*phi, 0, s2))

# Duplikate entfernen
V = []
for v in verts:
    if v not in V:
        V.append(v)
V = np.array(V, float)
assert len(V) == 12

def tet_volume(idx):
    a,b,c,d = V[list(idx)]
    M = np.vstack([b-a, c-a, d-a]).T
    return abs(np.linalg.det(M))/6.0

eps = 1e-9

# alle echten Tetraeder
tets = []
for comb in itertools.combinations(range(12), 4):
    if tet_volume(comb) > eps:
        tets.append(tuple(comb))

print("Anzahl echter Tetraeder:", len(tets))  # 420

# Disjunktheits-Grade
tset = [set(t) for t in tets]
deg = [0]*len(tets)
E = 0
for i in range(len(tets)):
    for j in range(i+1, len(tets)):
        if tset[i].isdisjoint(tset[j]):
            E += 1
            deg[i] += 1
            deg[j] += 1

print("Anzahl ungerichteter disjunkter Paare:", E)  # 12450
print("Gradverteilung:", collections.Counter(deg))  # 57/59/61

# markierter Eckpunkt
marked = 0
tets_with_mark = [t for t in tets if marked in t]
print("Tetraeder mit markiertem Eckpunkt:", len(tets_with_mark))  # 140

def count_second(t1):
    rest = [i for i in range(12) if i not in t1]
    c = 0
    for comb in itertools.combinations(rest, 4):
        if tet_volume(comb) > eps:
            c += 1
    return c

counts = [count_second(t) for t in tets_with_mark]
print("Restverteilung (57/59/61):", collections.Counter(counts))
print("Summe (Paare mit genau einem markierten):", sum(counts))  # 8300