# e8_weyl_shape_rigidity_audit.sage
#
# Exact rational SageMath audit for Weyl(E8) rigidity on
# V_shape = Sym_0(8), dim = 35.
#
# Tests:
#   G1: simple reflections generate the E8 Cartan data exactly.
#   G2: induced conjugation action preserves Sym_0(8).
#   G3: Fix_W(Sym_0(8)) = {0}.
#   G4: dimension of W-invariant symmetric bilinear forms on Sym_0(8) is 1.
#
# Interpretation:
#   G3 => no nonzero Weyl-fixed linear shear direction.
#   G4 => no preferred quadratic shear direction: every invariant symmetric
#         bilinear form is proportional to the Frobenius metric.
#
# IMPORTANT:
#   This is a representation-theoretic statement at the unperturbed E8 point.
#   It does NOT imply that all higher-order Weyl-invariant anisotropies vanish,
#   and it does NOT establish Navier-Stokes regularity.

from sage.all import *

Q = QQ
n = 8

# -------------------------------------------------------------------------
# 1. Bourbaki-compatible E8 simple roots in R^8, all norm^2 = 2.
# -------------------------------------------------------------------------
half = Q(1)/2
roots = [
    vector(Q, [ half,-half,-half,-half,-half,-half,-half, half]),
    vector(Q, [ 1, 1, 0, 0, 0, 0, 0, 0]),
    vector(Q, [-1, 1, 0, 0, 0, 0, 0, 0]),
    vector(Q, [ 0,-1, 1, 0, 0, 0, 0, 0]),
    vector(Q, [ 0, 0,-1, 1, 0, 0, 0, 0]),
    vector(Q, [ 0, 0, 0,-1, 1, 0, 0, 0]),
    vector(Q, [ 0, 0, 0, 0,-1, 1, 0, 0]),
    vector(Q, [ 0, 0, 0, 0, 0,-1, 1, 0]),
]

Cartan = matrix(Q, [[a.dot_product(b) for b in roots] for a in roots])

expected_Cartan = matrix(Q, [
[ 2, 0,-1, 0, 0, 0, 0, 0],
[ 0, 2, 0,-1, 0, 0, 0, 0],
[-1, 0, 2,-1, 0, 0, 0, 0],
[ 0,-1,-1, 2,-1, 0, 0, 0],
[ 0, 0, 0,-1, 2,-1, 0, 0],
[ 0, 0, 0, 0,-1, 2,-1, 0],
[ 0, 0, 0, 0, 0,-1, 2,-1],
[ 0, 0, 0, 0, 0, 0,-1, 2],
])

assert Cartan == expected_Cartan

I8 = identity_matrix(Q, n)

def reflection(alpha):
    # norm^2(alpha)=2, hence s_alpha(v)=v-(v.alpha)alpha
    return I8 - alpha.column() * alpha.row()

S = [reflection(a) for a in roots]
for s in S:
    assert s.transpose()*s == I8
    assert s*s == I8
    assert s.det() == -1

# -------------------------------------------------------------------------
# 2. Rational basis of Sym_0(8), dimension 35.
#
# Off-diagonal B_ij = E_ij + E_ji, 28 vectors.
# Diagonal H_i = E_ii - E_{i+1,i+1}, 7 vectors.
# -------------------------------------------------------------------------
shape_basis = []
shape_names = []

for i in range(n):
    for j in range(i+1, n):
        M = matrix(Q, n, n, 0)
        M[i,j] = 1
        M[j,i] = 1
        shape_basis.append(M)
        shape_names.append("S_%d_%d" % (i+1,j+1))

for i in range(n-1):
    M = matrix(Q, n, n, 0)
    M[i,i] = 1
    M[i+1,i+1] = -1
    shape_basis.append(M)
    shape_names.append("H_%d" % (i+1))

m = len(shape_basis)
assert m == 35

def shape_coords(X):
    assert X == X.transpose()
    assert X.trace() == 0

    out = []

    # Off-diagonal coordinates
    for i in range(n):
        for j in range(i+1, n):
            out.append(X[i,j])

    # Diagonal coordinates:
    # sum c_i (E_ii-E_(i+1)(i+1)) = diag(d_1,...,d_8)
    # gives c_i = d_1+...+d_i.
    running = Q(0)
    for i in range(n-1):
        running += X[i,i]
        out.append(running)

    v = vector(Q, out)

    # exact reconstruction certificate
    Y = matrix(Q, n, n, 0)
    for c,B in zip(v,shape_basis):
        Y += c*B
    assert Y == X
    return v

def induced_shape_rep(s):
    # column j is coordinate vector of s B_j s^T
    cols = []
    for B in shape_basis:
        X = s*B*s.transpose()
        assert X == X.transpose()
        assert X.trace() == 0
        cols.append(shape_coords(X))
    return matrix(Q, m, m, lambda i,j: cols[j][i])

R = [induced_shape_rep(s) for s in S]

for r in R:
    assert r*r == identity_matrix(Q,m)

# -------------------------------------------------------------------------
# 3. G3: no Weyl-fixed shape vector.
# Simple generators suffice.
# -------------------------------------------------------------------------
Afix = block_matrix([[r - identity_matrix(Q,m)] for r in R], subdivide=False)
fixed_kernel = Afix.right_kernel()
d_fixed = fixed_kernel.dimension()

# -------------------------------------------------------------------------
# 4. G4: invariant symmetric bilinear forms.
#
# Unknown G is symmetric 35x35, so 630 variables.
# Impose r^T G r = G for all eight simple reflections.
# Since the shape basis is rational but not orthonormal, this directly tests
# bilinear-form invariance and does not assume Euclidean coordinates in R^35.
# -------------------------------------------------------------------------
pairs = [(i,j) for i in range(m) for j in range(i,m)]
Nvar = len(pairs)
assert Nvar == m*(m+1)//2 == 630

def sym_matrix_from_coeff_vector(c):
    G = matrix(Q,m,m,0)
    for a,(i,j) in enumerate(pairs):
        G[i,j] = c[a]
        G[j,i] = c[a]
    return G

# Build sparse linear constraint operator column-by-column.
# Each column corresponds to one symmetric basis matrix E^(a).
# Rows concatenate all entries of r^T G r - G over the 8 generators.
nrows = len(R)*m*m
Mcon = matrix(Q, nrows, Nvar, sparse=True)

for a,(i,j) in enumerate(pairs):
    G = matrix(Q,m,m,0)
    G[i,j] = 1
    G[j,i] = 1
    if i == j:
        G[i,i] = 1

    row0 = 0
    for r in R:
        D = r.transpose()*G*r - G
        for p in range(m):
            for q in range(m):
                val = D[p,q]
                if val:
                    Mcon[row0 + p*m + q, a] = val
        row0 += m*m

Kbil = Mcon.right_kernel()
d_inv_bilinear = Kbil.dimension()

# -------------------------------------------------------------------------
# 5. Frobenius Gram matrix must be one invariant form.
# -------------------------------------------------------------------------
GF = matrix(Q,m,m,lambda i,j:
            sum(shape_basis[i][p,q]*shape_basis[j][p,q]
                for p in range(n) for q in range(n)))

assert GF.det() != 0

for r in R:
    assert r.transpose()*GF*r == GF

# Verify uniqueness if kernel dimension is 1.
frobenius_unique = False
if d_inv_bilinear == 1:
    kb = Kbil.basis()[0]
    G0 = sym_matrix_from_coeff_vector(kb)

    # Find a nonzero entry and compare proportionality exactly.
    ratio = None
    for i in range(m):
        for j in range(m):
            if G0[i,j] != 0:
                ratio = GF[i,j] / G0[i,j]
                break
        if ratio is not None:
            break
    assert ratio is not None
    assert GF == ratio*G0
    frobenius_unique = True

print("="*88)
print("E8 WEYL SHAPE-RIGIDITY AUDIT -- EXACT QQ")
print("="*88)
print("shape dimension                         :", m)
print("number of simple Weyl generators        :", len(R))
print("Fix_W(Sym_0(8)) dimension               :", d_fixed)
print("invariant symmetric bilinear dim        :", d_inv_bilinear)
print("Frobenius metric invariant              :", True)
print("Frobenius unique up to scale            :", frobenius_unique)
print("-"*88)

if d_fixed == 0:
    print("G3  NO WEYL-FIXED LINEAR SHEAR DIRECTION       : EXACT PASS")
else:
    print("G3  NO WEYL-FIXED LINEAR SHEAR DIRECTION       : FAIL")

if d_inv_bilinear == 1 and frobenius_unique:
    print("G4  UNIQUE QUADRATIC WEYL METRIC UP TO SCALE   : EXACT PASS")
else:
    print("G4  UNIQUE QUADRATIC WEYL METRIC UP TO SCALE   : NOT ESTABLISHED")

print()
print("Scope:")
print("  - exact representation-theoretic statement at the E8 symmetric point")
print("  - excludes preferred directions only at linear/quadratic order")
print("  - does NOT exclude higher-order Weyl-invariant anisotropies")
print("  - does NOT imply dynamical PDE coercivity or Navier-Stokes regularity")
