#!/usr/bin/env python3
import itertools
import numpy as np

# -----------------------------
# Geometry: Cube vertices (enthält 2 reguläre Tetraeder; Ikosaeder enthält keine)
# -----------------------------
def cube_vertices():
    """Würfel (±1,±1,±1) – enthält genau 2 einbeschriebene reguläre Tetraeder (Stella Octangula)."""
    pts = []
    for x in (-1.0, 1.0):
        for y in (-1.0, 1.0):
            for z in (-1.0, 1.0):
                pts.append((x, y, z))
    return np.array(pts, dtype=float)

# -----------------------------
# Find regular tetrahedra among vertices
# -----------------------------
def is_regular_tetra(indices, V, tol=1e-7):
    # 6 pairwise distances
    pts = V[list(indices)]
    d = []
    for i in range(4):
        for j in range(i+1,4):
            d.append(np.linalg.norm(pts[i]-pts[j]))
    d = np.array(d)
    return np.max(d) - np.min(d) < tol

def find_inscribed_tetrahedra(V, tol=1e-7):
    tetras = []
    for comb in itertools.combinations(range(len(V)), 4):
        if is_regular_tetra(comb, V, tol=tol):
            tetras.append(tuple(comb))
    # Deduplicate by sorting indices (already sorted) and unique
    tetras = sorted(set(tetras))
    return tetras

# -----------------------------
# Surface lattice points in barycentric coords
# a+b+c+d = n, and at least one coordinate is 0
# -----------------------------
def bary_surface_points(n):
    pts = []
    # Enumerate a,b,c and compute d = n-a-b-c
    # Restrict to boundary: min(a,b,c,d)==0
    for a in range(n+1):
        for b in range(n+1-a):
            for c in range(n+1-a-b):
                d = n - a - b - c
                if min(a,b,c,d) == 0:
                    pts.append((a,b,c,d))
    return np.array(pts, dtype=np.int32)

def bary_to_xyz(bary, tet_vertices):
    # bary: (N,4) integer, sum=n
    # map to xyz by normalized barycentric weights
    n = np.sum(bary[0])
    w = bary.astype(float) / float(n)
    # tet_vertices: (4,3)
    return w @ tet_vertices

# -----------------------------
# Hashing points with quantization for robust set operations
# -----------------------------
def quantize_points(P, q=1e-6):
    # convert to integer grid keys
    return np.round(P / q).astype(np.int64)

def keys_set(Q):
    # Q: (N,3) int64
    return { (int(x),int(y),int(z)) for x,y,z in Q }

# -----------------------------
# Main test: for n=1..Nmax
# -----------------------------
def test_up_to(Nmax=1000, q=1e-6, tol_tetra=1e-7, choose_triplet=(0,1,0), verbose=True):
    V = cube_vertices()
    tetras = find_inscribed_tetrahedra(V, tol=tol_tetra)
    if verbose:
        print(f"Found {len(tetras)} inscribed regular tetrahedra among cube vertices.")
        for i,t in enumerate(tetras[:10]):
            print(" tetra", i, t)

    if len(tetras) < 2:
        raise RuntimeError("Not enough inscribed tetrahedra found. Try increasing tol_tetra.")

    # pick three tetrahedra for RGB
    iR, iG, iB = choose_triplet
    TR = np.array(V[list(tetras[iR])], dtype=float)
    TG = np.array(V[list(tetras[iG])], dtype=float)
    TB = np.array(V[list(tetras[iB])], dtype=float)

    # precompute target tetra surfaces (all tetrahedra) as sets of keys for each n on the fly
    results = []

    for n in range(1, Nmax+1):
        bary = bary_surface_points(n)

        PR = bary_to_xyz(bary, TR)
        PG = bary_to_xyz(bary, TG)
        PB = bary_to_xyz(bary, TB)

        KR = keys_set(quantize_points(PR, q=q))
        KG = keys_set(quantize_points(PG, q=q))
        KB = keys_set(quantize_points(PB, q=q))

        # Exclusive (one-color) points
        ER = KR - (KG | KB)
        EG = KG - (KR | KB)
        EB = KB - (KR | KG)

        # Build all tetra target surface-keysets to match against
        tetra_surfaces = []
        for t in tetras:
            T = np.array(V[list(t)], dtype=float)
            P = bary_to_xyz(bary, T)
            K = keys_set(quantize_points(P, q=q))
            tetra_surfaces.append(K)

        # Check if any exclusive set equals any tetra surface
        hit = []
        for color, Eset in (("R",ER), ("G",EG), ("B",EB)):
            if len(Eset) == 0:
                continue
            for j, Ktet in enumerate(tetra_surfaces):
                if Eset == Ktet:
                    hit.append((color, j))
                    break

        if hit:
            results.append((n, hit))
            if verbose:
                print(f"n={n}: monochrome tetra emerged: {hit}")

    return results

if __name__ == "__main__":
    # choose_triplet: (0,1,0) = R und B nutzen Tetra 0, G nutzt Tetra 1 (Würfel hat nur 2 Tetraeder)
    hits = test_up_to(Nmax=100, q=1e-6, tol_tetra=1e-7, choose_triplet=(0,1,0), verbose=True)
    print("\n=== SUMMARY ===")
    print(f"Total hits: {len(hits)}")
    if hits:
        ns = [n for n,_ in hits]
        print("Hit n values:", ns)