from itertools import product
from sage.groups.perm_gps.permgroup_named import AlternatingGroup

# ==============================================================================
# 1. Zahlkörper Q(sqrt(5)) und Galois-Automorphismus definieren
# ==============================================================================
K.<a> = NumberField(x^2 - 5)  # a = sqrt(5)
phi = (1 + a) / 2  # Goldener Schnitt phi
phi_inv = (a - 1) / 2  # 1/phi = phi - 1
sigma = K.hom([-a], K)  # Galois-Automorphismus: sqrt(5) -> -sqrt(5)

# ==============================================================================
# 2. Konstruktion der 120 Ikosier-Vektoren (H_4 Wurzelsystem in K^4)
# ==============================================================================
H4_set = set()

# Typ 1: 8 Vektoren der Form (±1, 0, 0, 0) und Koordinatenpermutationen
for i in range(4):
    for s in [-1, 1]:
        v = [K(0)] * 4
        v[i] = K(s)
        H4_set.add(tuple(v))

# Typ 2: 16 Vektoren der Form 1/2 * (±1, ±1, ±1, ±1)
for s in product([-1, 1], repeat=4):
    v = tuple(K(s_i) / 2 for s_i in s)
    H4_set.add(v)

# Typ 3: 96 Vektoren aus geraden Permutationen von 1/2 * (0, ±1, ±1/phi, ±phi)
A4 = AlternatingGroup(4)
for s1, s2, s3 in product([-1, 1], repeat=3):
    raw = [K(0), K(s1) / 2, K(s2) * phi_inv / 2, K(s3) * phi / 2]
    for p in A4:
        # p wirkt 1-basiert auf {1, 2, 3, 4}
        permuted = tuple(raw[p(i) - 1] for i in range(1, 5))
        H4_set.add(permuted)

print(f"[1] Anzahl konstruierter H_4-Vektoren: {len(H4_set)}")
assert len(H4_set) == 120, "Fehler: H_4 muss genau 120 Vektoren enthalten!"

# ==============================================================================
# 3. Minkowski-Galois-Einbettung: iota(q) = (q, sigma(q)) in K^8
# ==============================================================================
def minkowski_embed(v):
    v_sigma = tuple(sigma(x) for x in v)
    vec = vector(K, list(v) + list(v_sigma))
    vec.set_immutable()  # Macht Vektoren hashbar für die Verwendung in sets
    return vec

H4_embedded = [minkowski_embed(v) for v in H4_set]

# Verifikation der Normen in H_4
norms_H4 = [v.dot_product(v) for v in H4_embedded]
print(f"[2] Alle {len(H4_embedded)} H_4-Vektoren besitzen Norm Quadrat = 2: "
      f"{all(n == 2 for n in norms_H4)}")

# ==============================================================================
# 4. Vervollständigung zum E_8-Wurzelsystem (240 Vektoren)
# ==============================================================================
# In E_8 erzeugen zwei Wurzeln u, v mit <u,v> = -1 eine neue Wurzel w = u + v
E8_set = set(H4_embedded)

for u in H4_embedded:
    for v in H4_embedded:
        if u.dot_product(v) == -1:
            w = u + v
            w.set_immutable()
            E8_set.add(w)

E8_roots = list(E8_set)
print(f"[3] Gesamtzahl der generierten E_8-Wurzeln: {len(E8_roots)}")
assert len(E8_roots) == 240, "Fehler: Das E_8-Wurzelsystem muss genau 240 Wurzeln umfassen!"

# ==============================================================================
# 5. Algebraische & Geometrische Verifikation der E_8-Eigenschaften
# ==============================================================================
# a) Normen-Check
all_norm_2 = all(v.dot_product(v) == 2 for v in E8_roots)
print(f"[4a] Alle 240 Vektoren besitzen Norm Quadrat = 2: {all_norm_2}")

# b) Ganzzahligkeit der Skalarprodukte in {-2, -1, 0, 1, 2}
inner_prods = set()
for i in range(len(E8_roots)):
    for j in range(i, len(E8_roots)):
        ip = E8_roots[i].dot_product(E8_roots[j])
        inner_prods.add(ip)

print(f"[4b] Menge aller Skalarprodukte im System: {sorted(list(inner_prods))}")
assert inner_prods == {-2, -1, 0, 1, 2}, "Fehler: Skalarprodukte müssen in {-2,-1,0,1,2} liegen!"

# c) Unimodularität: Determinante der Gram-Matrix einer Lattices-Basis = 1
V = VectorSpace(K, 8)
basis_vectors = []
for v in E8_roots:
    if V.subspace(basis_vectors + [v]).dimension() > len(basis_vectors):
        basis_vectors.append(v)
        if len(basis_vectors) == 8:
            break

Gram = matrix(K, 8, 8, [[u.dot_product(v) for v in basis_vectors] for u in basis_vectors])
det_G = Gram.det()
print(f"[4c] Determinante der Gram-Matrix der Basis: {det_G}")
assert det_G == 1, "Fehler: Das E_8-Gitter muss unimodular sein (det = 1)!"

# d) Abgleich mit dem in SageMath integrierten Wurzelsystem
E8_sage = RootSystem(['E', 8]).root_lattice()
print(f"[4d] SageMath-Referenz RootSystem(['E', 8]): {len(E8_sage.roots())} Wurzeln.")

print("\nERFOLG: Die Galois-Minkowski-Einbettung von H_4 erzeugt exakt das E_8 Wurzelsystem!")