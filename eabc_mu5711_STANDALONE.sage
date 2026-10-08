# eabc_mu5711_STANDALONE.sage
# SageMath 10.5+ / 10.9
#
# Self-contained production bundle:
#   Rehm/Coxeter fiber generator
#   corrected Wilson G2(2) automorphism audit
#   optimized (5,7,11) production orbit reduction
#
# No external .sage files are required.
#

# ===== MODULE 1: REHM FIBERS =====
# eabc_rehm_fiber_generator.sage
# SageMath 10.5+
#
# Exact generator for Hans Peter Rehm's sets
#
#   P(p) = { pi in C : n(pi)=p, pi > 0, pi == 1 (mod 2) },  p > 2,
#
# where C is Coxeter's maximal order of integral Cayley octaves.
#
# Source convention:
# Coxeter/Rehm Q-basis e0=1,e1,...,e7 with oriented products determined by
#     e_{rho(1)} = e_{rho(2)} e_{rho(4)}
# for powers rho of the 7-cycle (1 2 3 4 5 6 7).
#
# Coxeter Z-basis:
#   (1, e1, e2, e3, h, e1*h, e2*h, e3*h),
#   h = 1/2 (e1 + e2 + e3 - e4).
#
# IMPORTANT:
# Rehm fixes any total additive ordering on C; here we choose lexicographic
# order in this Coxeter Z-basis. Since pi == 1 mod 2 forces the first
# coefficient odd, "pi > 0" becomes c[0] > 0 for this lex choice.
#
# The generator enumerates the norm-p shell of the Coxeter lattice using
# the integral Gram matrix 2G (an E8 Gram matrix), then applies exactly
# Rehm's two normalizing filters:
#   (i) coefficient vector == (1,0,...,0) mod 2
#   (ii) lex-positive.
#
# It does NOT quotient by units. For odd p, Rehm's definition itself is the
# normalization. Proposition 6.4 proves |P(p)| = p^3 + 1.
#
# This is intentionally separate from the G2/deficit backend.

from sage.all import *

QQ8 = VectorSpace(QQ, 8)
ZZ8 = VectorSpace(ZZ, 8)

def qv(xs):
    return QQ8([QQ(x) for x in xs])

ONE = qv([1,0,0,0,0,0,0,0])
E = [qv([1 if j == i else 0 for j in range(8)]) for i in range(8)]

# Rehm/Coxeter Fano orientation:
# e_{rho(2)} e_{rho(4)} = e_{rho(1)}
FANO_REHM = [
    (2,4,1),
    (3,5,2),
    (4,6,3),
    (5,7,4),
    (6,1,5),
    (7,2,6),
    (1,3,7),
]

_MUL = {}
for a,b,c in FANO_REHM:
    _MUL[(a,b)] = ( 1,c)
    _MUL[(b,c)] = ( 1,a)
    _MUL[(c,a)] = ( 1,b)
    _MUL[(b,a)] = (-1,c)
    _MUL[(c,b)] = (-1,a)
    _MUL[(a,c)] = (-1,b)

def basis_mul(i,j):
    if i == 0:
        return E[j]
    if j == 0:
        return E[i]
    if i == j:
        return -ONE
    s,k = _MUL[(i,j)]
    return s*E[k]

def omul(x,y):
    x,y = qv(x),qv(y)
    z = QQ8.zero()
    for i in range(8):
        if x[i] == 0:
            continue
        for j in range(8):
            if y[j] == 0:
                continue
            z += x[i]*y[j]*basis_mul(i,j)
    return z

def norm2(x):
    x=qv(x)
    return sum(a*a for a in x)

def inner(x,y):
    x,y=qv(x),qv(y)
    return sum(x[i]*y[i] for i in range(8))

# Coxeter Z-basis
h = (E[1] + E[2] + E[3] - E[4]) / 2
COXETER_BASIS = [
    E[0], E[1], E[2], E[3],
    h,
    omul(E[1],h),
    omul(E[2],h),
    omul(E[3],h),
]

B = matrix(QQ, 8, 8, [x for b in COXETER_BASIS for x in b])
# Rows are basis vectors in the orthonormal e-coordinate system.
G = B * B.transpose()

# 2G is integral and unimodular; norm in Rehm's convention is
# n(sum c_i b_i) = c^T G c.
G2 = matrix(ZZ, 2*G)

def coxeter_coords_to_octonion(c):
    c = vector(ZZ,c)
    return qv(c * B)

def octonion_to_coxeter_coords(x):
    x = qv(x)
    c = vector(QQ,x) * B.inverse()
    if any(a.denominator() != 1 for a in c):
        raise ValueError("x is not in Coxeter's order C.")
    return vector(ZZ,[ZZ(a) for a in c])

def in_coxeter_order(x):
    try:
        octonion_to_coxeter_coords(x)
        return True
    except ValueError:
        return False

def rehm_mod2_normalized(c):
    """
    pi == 1 mod 2 in C iff its Coxeter Z-basis coefficients satisfy
    c == (1,0,...,0) mod 2.
    """
    c = vector(ZZ,c)
    return (c[0] % 2 == 1 and all(c[i] % 2 == 0 for i in range(1,8)))

def rehm_positive(c):
    """
    Lexicographic positive order in Coxeter Z-basis.
    For normalized vectors c[0] is odd, hence nonzero.
    """
    c = vector(ZZ,c)
    for a in c:
        if a != 0:
            return a > 0
    return False

def shell_coefficients_norm_p(p, up_to_sign=False):
    """
    Enumerate all Coxeter-basis coefficient vectors c with n(c)=p.

    We enumerate in the integral lattice with Gram matrix 2G.
    Its squared length is 2*n(c), hence target length 2p.
    """
    p = ZZ(p)
    if p <= 0:
        raise ValueError("p must be positive.")

    L = IntegralLattice(G2)
    # short_vectors(n) returns vectors of length < n.
    # Hence request up through length 2p.
    sv = L.short_vectors(2*p + 1, up_to_sign_flag=up_to_sign)
    return [vector(ZZ,c) for c in sv[2*p]]

def generate_rehm_P(p, audit_full_shell=True):
    """
    Generate Rehm's P(p) for an odd prime p.

    Exact definition (Rehm Prop. 5.6):
        P(p) = {pi in C | n(pi)=p, pi>0, pi == 1 mod 2}.

    Returns octonions in the orthonormal e0,...,e7 coordinates used
    by the geometric backend.
    """
    p = ZZ(p)
    if p <= 2 or not p.is_prime():
        raise ValueError("This generator implements Rehm P(p) for odd primes p > 2.")

    shell = shell_coefficients_norm_p(p, up_to_sign=False)

    if audit_full_shell:
        expected_shell = 240*(p**3 + 1)
        if len(shell) != expected_shell:
            raise AssertionError(
                "Full norm-p shell count failed: "
                f"{len(shell)} != 240*(p^3+1)={expected_shell}"
            )

    coeffs = [
        c for c in shell
        if rehm_mod2_normalized(c) and rehm_positive(c)
    ]

    expected = p**3 + 1
    if len(coeffs) != expected:
        raise AssertionError(
            f"Rehm P({p}) cardinality failed: {len(coeffs)} != {expected}"
        )

    reps = [coxeter_coords_to_octonion(c) for c in coeffs]

    # Exact norm audit
    for i,pi in enumerate(reps):
        if norm2(pi) != p:
            raise AssertionError(f"Norm failure in P({p}) at index {i}.")

    return reps

def audit_coxeter_basis():
    """
    Structural checks on the lattice convention.
    """
    assert G2.det() == 1, "2G should be an E8 unimodular Gram matrix."
    assert all(x in ZZ for x in G2.list())
    # Norm-one shell should have 240 vectors.
    shell1 = shell_coefficients_norm_p(1, up_to_sign=False)
    assert len(shell1) == 240
    return {
        "det_2G": G2.det(),
        "norm1_shell": len(shell1),
        "G": G,
        "2G": G2,
    }

def install_into_backend(backend, p):
    """
    Convenience adapter for eabc_s6_9_octavian_backend.sage FiberBackend.
    """
    reps = generate_rehm_P(p)
    backend.add_fiber(p, reps)
    return len(reps)


# ===== MODULE 2: CORRECTED G2(2) AUTOMORPHISMS =====
# SageMath 10.5+
# Corrected AUT-1 ... AUT-6 audit for Rehm/Coxeter order
#
# Requires:
#   eabc_rehm_fiber_generator.sage
#
# Key correction:
# Wilson's chosen integral-octonion order is shifted relative to the
# Rehm/Coxeter order used here.  If S : i_t -> i_{t+1}, the Wilson
# stabilizer generators are transferred by A -> S^(-1) A S.
#
# Row-vector convention throughout.

from sage.all import *
from collections import Counter, deque

# bundled: Rehm generator already included

W2R = {0:1, 1:3, 2:5, 3:7, 4:2, 5:4, 6:6}

def wi(t):
    return E[W2R[ZZ(t) % 7]]

def linmap_from_wilson_images(imgw):
    inv = {r:t for t,r in W2R.items()}
    images = [E[0]] + [imgw[inv[r]] for r in range(1,8)]
    return matrix(QQ, [list(v) for v in images])

def wilson_shift():
    return linmap_from_wilson_images({t:wi(t+1) for t in range(7)})

def wilson_b2():
    # Wilson monomial generator i_t -> i_{2t}
    return linmap_from_wilson_images({t:wi(2*t) for t in range(7)})

def wilson_c():
    # (i0,...,i6) -> (-i0,-i1,i6,i3,i5,i4,i2)
    d = {
        0:-wi(0), 1:-wi(1), 2:wi(6), 3:wi(3),
        4:wi(5), 5:wi(4), 6:wi(2)
    }
    return linmap_from_wilson_images(d)

def wilson_d():
    # (i0,...,i6) -> (i0,-i1,-i4,-i3,-i2,-i6,-i5)
    d = {
        0:wi(0), 1:-wi(1), 2:-wi(4), 3:-wi(3),
        4:-wi(2), 5:-wi(6), 6:-wi(5)
    }
    return linmap_from_wilson_images(d)

def wilson_epsilon():
    d = {}
    d[0] = wi(6); d[6] = wi(2); d[2] = wi(0)
    d[1] = ( wi(1)+wi(3)-wi(4)+wi(5) )/2
    d[3] = (-wi(1)+wi(3)+wi(4)+wi(5) )/2
    d[4] = ( wi(1)-wi(3)+wi(4)+wi(5) )/2
    d[5] = (-wi(1)-wi(3)-wi(4)+wi(5) )/2
    return linmap_from_wilson_images(d)

def orth_to_coxeter(A):
    M = B * A * B.inverse()
    assert all(q.denominator() == 1 for q in M.list())
    return matrix(ZZ, M)

def corrected_wilson_generators():
    S = wilson_shift()
    # Transfer Wilson's 0-order stabilizer to our Rehm/Coxeter order.
    Asrc = {
        "b2": wilson_b2(),
        "c": wilson_c(),
        "d": wilson_d(),
        "eps": wilson_epsilon(),
    }
    return {
        name: orth_to_coxeter(S.inverse() * A * S)
        for name,A in Asrc.items()
    }

def unitrow(i):
    return vector(ZZ, [1 if k == i else 0 for k in range(8)])

def audit_generator(M, name):
    # AUT-1
    assert unitrow(0) * M == unitrow(0)

    # AUT-2, row convention
    assert M * G2 * M.transpose() == G2

    # AUT-3: all 64 Coxeter-basis products
    images = [
        coxeter_coords_to_octonion(unitrow(i) * M)
        for i in range(8)
    ]
    for i in range(8):
        for j in range(8):
            lhs = omul(images[i], images[j])
            bij = omul(COXETER_BASIS[i], COXETER_BASIS[j])
            cij = octonion_to_coxeter_coords(bij)
            rhs = coxeter_coords_to_octonion(cij * M)
            assert lhs == rhs, (name, i, j)
    return True

def rehm_sign_normalize_coeff(c):
    c = vector(ZZ,c)
    # Since c == e0 mod 2, c[0] is odd and nonzero.
    return c if c[0] > 0 else -c

def induced_permutation(M, Pp):
    idx = {
        tuple(octonion_to_coxeter_coords(x)): i
        for i,x in enumerate(Pp)
    }
    img = []
    for x in Pp:
        c = octonion_to_coxeter_coords(x)
        d = rehm_sign_normalize_coeff(c * M)
        assert rehm_mod2_normalized(d)
        assert rehm_positive(d)
        key = tuple(d)
        assert key in idx
        img.append(idx[key])
    assert len(set(img)) == len(Pp)
    return Permutation([i+1 for i in img])  # Sage permutations are 1-based

def orbit_lengths_from_generator_perms(perms, n):
    # Use zero-based lookup tables for a simple deterministic BFS.
    tabs = []
    for s in perms:
        tabs.append([s(i+1)-1 for i in range(n)])

    seen = [False]*n
    lengths = []
    for seed in range(n):
        if seen[seed]:
            continue
        q = [seed]
        seen[seed] = True
        size = 0
        while q:
            x = q.pop()
            size += 1
            for tab in tabs:
                y = tab[x]
                if not seen[y]:
                    seen[y] = True
                    q.append(y)
        lengths.append(size)
    return sorted(lengths, reverse=True)

def run_aut1_6():
    gens = corrected_wilson_generators()

    for name,M in gens.items():
        audit_generator(M,name)
        print("PASS AUT-1/2/3:", name)

    G = MatrixGroup(list(gens.values()))
    assert G.order() == 12096
    print("PASS AUT-4: group order =", G.order())

    expected = {
        5:  {"orbits":[126],     "max_gen_order":6},
        7:  {"orbits":[288,56],  "max_gen_order":6},
        11: {"orbits":[756,576], "max_gen_order":6},
    }

    rows = []
    for p in [5,7,11]:
        Pp = generate_rehm_P(p)
        assert len(Pp) == p**3 + 1

        sigmas = {
            name: induced_permutation(M,Pp)
            for name,M in gens.items()
        }
        orders = {name:s.order() for name,s in sigmas.items()}
        maxord = max(orders.values())
        orbits = orbit_lengths_from_generator_perms(list(sigmas.values()),len(Pp))

        assert maxord == expected[p]["max_gen_order"]
        assert orbits == expected[p]["orbits"]
        assert all(12096 % n == 0 for n in orbits)

        # Certificate: image contains an element of order > 2.
        assert maxord > 2
        kernel = 1

        rows.append((p,orders,orbits,kernel))
        print("p =",p,"generator orders =",orders,
              "orbits =",orbits,"=> K_p = 1")

    print("PASS AUT-6f")
    return gens,G,rows


# ===== MODULE 3: MU_5711 PRODUCTION =====
# eabc_mu5711_production.sage
# SageMath 10.5+
#
# Production orbit reduction for
#     X = P(5) x P(7) x P(11)
# under the diagonal G2(2)-action.
#
# Core optimization:
# P(5) is one transitive G2(2)-orbit of size 126.
# Fix i0 in P(5), let H = Stab_G(i0), |H| = 12096/126 = 96.
#
# Then:
#   G \ (P5 x P7 x P11)  <->  H \ (P7 x P11)
#
# Hence we never traverse 57,734,208 triples.
# We traverse only 344*1332 = 458,208 pairs.
#
# Full G-orbit size corresponding to an H-orbit O_H is
#     |O_G| = 126 * |O_H|.
#
# Geometry is evaluated exactly once per H-orbit representative,
# using the Rehm/Coxeter multiplication convention itself.
#
# Requires:
#   eabc_rehm_fiber_generator.sage
#   eabc_aut1_6_corrected.sage
#
# Output:
#   mu_5711_exact.json
#   mu_5711_exact.txt
#   mu_5711_orbits.csv

from sage.all import *
from collections import Counter
import json, csv, os, time

BASE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "."

# bundled: AUT/Rehm modules already included

# ----------------------------------------------------------------------
# Exact Rehm-convention geometry
# ----------------------------------------------------------------------

def realpart_rehm(x):
    return QQ(x[0])

def imagpart_rehm(x):
    x = vector(QQ, x)
    y = vector(QQ, x)
    y[0] = 0
    return y

def associator_rehm(x,y,z):
    return omul(omul(x,y),z) - omul(x,omul(y,z))

def phi_rehm(u,v,w):
    return inner(omul(u,v), w)

def triple_geometry_rehm(pi,rho,sigma):
    """
    Exact G2 geometry in the SAME Fano/Rehm multiplication convention
    used by eabc_rehm_fiber_generator.sage.

    Returns:
      detG_norm
      phi_raw
      phi2_norm
      t
    with exact identity
      t = detG_norm - phi2_norm
    and theorem-level bound 0 <= t <= 1.
    """
    u = imagpart_rehm(pi)
    v = imagpart_rehm(rho)
    w = imagpart_rehm(sigma)

    su,sv,sw = norm2(u), norm2(v), norm2(w)
    if 0 in (su,sv,sw):
        raise ValueError("Zero imaginary projection in Rehm representative.")

    Gm = matrix(QQ, [
        [inner(u,u), inner(u,v), inner(u,w)],
        [inner(v,u), inner(v,v), inner(v,w)],
        [inner(w,u), inner(w,v), inner(w,w)],
    ])

    D = QQ(su*sv*sw)
    detG_norm = Gm.det()/D

    ph = phi_rehm(u,v,w)
    phi2_norm = ph*ph/D

    a = associator_rehm(u,v,w)
    t = norm2(a)/(4*D)

    assert t == detG_norm - phi2_norm
    assert 0 <= t <= 1

    return {
        "detG_norm": detG_norm,
        "phi_raw": ph,
        "phi2_norm": phi2_norm,
        "t": t,
    }

# ----------------------------------------------------------------------
# Build the three faithful permutation actions
# ----------------------------------------------------------------------

def build_fibers_and_actions():
    print("=== Generate Rehm fibers ===")
    P5  = generate_rehm_P(5)
    P7  = generate_rehm_P(7)
    P11 = generate_rehm_P(11)

    assert len(P5)  == 126
    assert len(P7)  == 344
    assert len(P11) == 1332

    print("|P5| =",len(P5),"|P7| =",len(P7),"|P11| =",len(P11))

    gensM = corrected_wilson_generators()

    # Repeat the hard generator audits in the production script.
    for name,M in gensM.items():
        audit_generator(M,name)

    Gmat = MatrixGroup(list(gensM.values()))
    assert Gmat.order() == 12096

    sig5  = {name: induced_permutation(M,P5)  for name,M in gensM.items()}
    sig7  = {name: induced_permutation(M,P7)  for name,M in gensM.items()}
    sig11 = {name: induced_permutation(M,P11) for name,M in gensM.items()}

    # Deterministic fidelity certificates.
    for p,sigs in [(5,sig5),(7,sig7),(11,sig11)]:
        orders = {name:s.order() for name,s in sigs.items()}
        assert max(orders.values()) > 2
        print("p =",p,"generator permutation orders =",orders)

    return P5,P7,P11,gensM,sig5,sig7,sig11

def combined_permutation_group(sig5,sig7,sig11):
    """
    Faithful action on the disjoint union
        P5 | P7 | P11
    of degree 126+344+1332 = 1802.

    The block offsets allow us to recover each fiber action from a
    single permutation group element.
    """
    n5,n7,n11 = 126,344,1332
    N = n5+n7+n11

    names = list(sig5.keys())
    assert set(names) == set(sig7.keys()) == set(sig11.keys())

    S = SymmetricGroup(N)
    cg = []

    for name in names:
        images = []

        # P5 block: 1 ... n5
        images += [sig5[name](i) for i in range(1,n5+1)]

        # P7 block
        images += [n5 + sig7[name](j) for j in range(1,n7+1)]

        # P11 block
        images += [n5+n7 + sig11[name](k) for k in range(1,n11+1)]

        g = S(images)
        cg.append(g)

    G = PermutationGroup(cg)
    assert G.order() == 12096

    return G

# ----------------------------------------------------------------------
# Stabilizer reduction:
# G\(P5 x P7 x P11) == H\(P7 x P11), |H|=96
# ----------------------------------------------------------------------

def stabilizer_reduction(G):
    # First point in the P5 block.
    H = G.stabilizer(1)
    assert H.order() == 96

    print("PASS: |G| =",G.order())
    print("PASS: P5 stabilizer |H| =",H.order())
    print("Reduced pair space = 344 * 1332 =",344*1332)

    return H

def h_generator_tables(H):
    """
    Convert H generators to zero-based integer lookup arrays on P7/P11.
    """
    n5,n7,n11 = 126,344,1332
    off7 = n5
    off11 = n5+n7

    tabs = []
    for h in H.gens():
        t7 = [h(off7+j+1) - off7 - 1 for j in range(n7)]
        t11 = [h(off11+k+1) - off11 - 1 for k in range(n11)]

        assert sorted(t7) == list(range(n7))
        assert sorted(t11) == list(range(n11))
        tabs.append((t7,t11))

    return tabs

def enumerate_pair_orbits(H, progress_every=50000):
    """
    BFS on only 458,208 pair-indices.

    Pair id:
        q = j*n11 + k.

    Returns list of:
        {
          rep_pair_id,
          j,
          k,
          H_orbit_size,
          G_orbit_size,
          stabilizer_size
        }
    """
    n7,n11 = 344,1332
    N = n7*n11
    tabs = h_generator_tables(H)

    visited = bytearray(N)
    orbits = []

    processed = 0
    t0 = time.time()

    for seed in range(N):
        if visited[seed]:
            continue

        # Since seeds are scanned increasingly, seed is the canonical
        # minimum pair-id in its orbit.
        visited[seed] = 1
        queue = [seed]
        head = 0
        hsize = 0

        while head < len(queue):
            q = queue[head]
            head += 1
            hsize += 1

            j,k = divmod(q,n11)

            for t7,t11 in tabs:
                yj = t7[j]
                yk = t11[k]
                y = yj*n11 + yk

                if not visited[y]:
                    visited[y] = 1
                    queue.append(y)

        processed += hsize

        # Full diagonal G-orbit:
        # |G orbit| = [G:H] * |H orbit| = 126*hsize.
        gsize = 126*hsize
        assert 12096 % gsize == 0
        stab = 12096 // gsize

        j0,k0 = divmod(seed,n11)

        orbits.append({
            "rep_pair_id": seed,
            "j": j0,
            "k": k0,
            "H_orbit_size": hsize,
            "G_orbit_size": gsize,
            "stabilizer_size": stab,
        })

        if progress_every and processed % progress_every < hsize:
            print(
                f"pairs {processed}/{N}; "
                f"triplet orbits={len(orbits)}; "
                f"elapsed={time.time()-t0:.1f}s"
            )

    assert processed == N

    # A1-A3 on full G orbits.
    total_full = sum(o["G_orbit_size"] for o in orbits)
    assert total_full == 57734208
    assert all(12096 % o["G_orbit_size"] == 0 for o in orbits)
    assert all(1 <= o["G_orbit_size"] <= 12096 for o in orbits)
    assert len(orbits) >= 4773

    return orbits

# ----------------------------------------------------------------------
# Independent Burnside cross-check on H acting on P7 x P11
# ----------------------------------------------------------------------

def burnside_H(H):
    """
    Exact Burnside count for H acting diagonally on P7 x P11.

      N = 1/|H| sum_h f7(h) f11(h)

    This is equivalent to the G-orbit count on P5 x P7 x P11.
    Only 96 group elements are involved.
    """
    n5,n7,n11 = 126,344,1332
    off7 = n5
    off11 = n5+n7

    S = ZZ(0)
    rows = []

    for h in H:
        f7 = sum(1 for j in range(n7)
                 if h(off7+j+1) == off7+j+1)
        f11 = sum(1 for k in range(n11)
                  if h(off11+k+1) == off11+k+1)

        term = f7*f11
        S += term
        rows.append((h.order(), f7, f11, term))

    assert S % H.order() == 0
    N = S // H.order()

    return ZZ(N), rows

def burnside_full_G(G):
    """
    Optional independent full-G Burnside count:
      1/12096 sum_g f5(g) f7(g) f11(g)

    This is slower than burnside_H but still only loops over 12,096
    permutations of degree 1802, not 57 million triples.
    """
    n5,n7,n11 = 126,344,1332
    off7=n5
    off11=n5+n7

    S = ZZ(0)

    for g in G:
        f5 = sum(1 for i in range(n5)
                 if g(i+1) == i+1)
        f7 = sum(1 for j in range(n7)
                 if g(off7+j+1) == off7+j+1)
        f11 = sum(1 for k in range(n11)
                  if g(off11+k+1) == off11+k+1)

        S += f5*f7*f11

    assert S % G.order() == 0
    return ZZ(S // G.order())

# ----------------------------------------------------------------------
# Geometry once per orbit and exact mu_{5,7,11}
# ----------------------------------------------------------------------

def attach_geometry(orbits,P5,P7,P11):
    """
    Every G-orbit has a representative with first coordinate P5[0].
    Compute exact Rehm-convention geometry once for that representative.
    """
    pi0 = P5[0]

    for n,o in enumerate(orbits,1):
        rho = P7[o["j"]]
        sigma = P11[o["k"]]

        geo = triple_geometry_rehm(pi0,rho,sigma)
        t = QQ(geo["t"])

        # A4
        assert 0 <= t <= 1

        o["t"] = t
        o["detG_norm"] = QQ(geo["detG_norm"])
        o["phi2_norm"] = QQ(geo["phi2_norm"])

        if n % 500 == 0:
            print("geometry",n,"/",len(orbits))

    return orbits

def aggregate_mu(orbits):
    hist = Counter()
    total = ZZ(0)

    for o in orbits:
        w = ZZ(o["G_orbit_size"])
        hist[QQ(o["t"])] += w
        total += w

    assert total == 57734208

    mu = {t: QQ(n)/total for t,n in sorted(hist.items())}

    # A5
    assert sum(mu.values()) == 1

    M1 = sum(t*prob for t,prob in mu.items())
    M2 = sum(t*t*prob for t,prob in mu.items())
    var = M2-M1*M1
    passoc = mu.get(QQ(0),QQ(0))

    return {
        "total": total,
        "histogram": hist,
        "mu": mu,
        "support_size": len(mu),
        "P_assoc": passoc,
        "M1": M1,
        "M2": M2,
        "variance": var,
    }

# ----------------------------------------------------------------------
# Serialization
# ----------------------------------------------------------------------

def save_results(orbits,stats,N_bfs,N_burnside,N_burnside_full=None):
    csv_path = os.path.join(BASE,"mu_5711_orbits.csv")
    json_path = os.path.join(BASE,"mu_5711_exact.json")
    txt_path = os.path.join(BASE,"mu_5711_exact.txt")

    with open(csv_path,"w",newline="") as f:
        w=csv.writer(f)
        w.writerow([
            "orbit_id","j_P7","k_P11",
            "H_orbit_size","G_orbit_size","stabilizer_size",
            "t","detG_norm","phi2_norm"
        ])
        for idx,o in enumerate(orbits):
            w.writerow([
                idx,o["j"],o["k"],
                o["H_orbit_size"],o["G_orbit_size"],o["stabilizer_size"],
                str(o["t"]),str(o["detG_norm"]),str(o["phi2_norm"])
            ])

    payload = {
        "total_triples": int(stats["total"]),
        "triplet_orbits_BFS": int(N_bfs),
        "triplet_orbits_Burnside_H": int(N_burnside),
        "triplet_orbits_Burnside_full_G":
            None if N_burnside_full is None else int(N_burnside_full),
        "support_size": int(stats["support_size"]),
        "P_assoc": str(stats["P_assoc"]),
        "M1": str(stats["M1"]),
        "M2": str(stats["M2"]),
        "variance": str(stats["variance"]),
        "mu": [
            {
                "t": str(t),
                "count": int(stats["histogram"][t]),
                "probability": str(prob)
            }
            for t,prob in stats["mu"].items()
        ]
    }

    with open(json_path,"w") as f:
        json.dump(payload,f,indent=2)

    with open(txt_path,"w") as f:
        f.write("EABC mu_{5,7,11} production audit\n")
        f.write("==================================\n")
        f.write(f"total_triples = {stats['total']}\n")
        f.write(f"N_orb_BFS = {N_bfs}\n")
        f.write(f"N_orb_Burnside_H = {N_burnside}\n")
        f.write(f"N_orb_Burnside_full_G = {N_burnside_full}\n")
        f.write(f"support_size = {stats['support_size']}\n")
        f.write(f"P_assoc = {stats['P_assoc']}\n")
        f.write(f"M1 = {stats['M1']}\n")
        f.write(f"M2 = {stats['M2']}\n")
        f.write(f"variance = {stats['variance']}\n")
        f.write("\nmu:\n")
        for t,prob in stats["mu"].items():
            f.write(
                f"t={t}\tcount={stats['histogram'][t]}"
                f"\tprob={prob}\n"
            )

    print("saved:",csv_path)
    print("saved:",json_path)
    print("saved:",txt_path)

# ----------------------------------------------------------------------
# Main production pipeline
# ----------------------------------------------------------------------

def run_mu5711_production(run_full_G_burnside=True):
    print("=== EABC OCT-G3e/f production ===")

    P5,P7,P11,gensM,sig5,sig7,sig11 = build_fibers_and_actions()

    G = combined_permutation_group(sig5,sig7,sig11)
    H = stabilizer_reduction(G)

    print("\n=== Direct orbit enumeration on H \\ (P7 x P11) ===")
    orbits = enumerate_pair_orbits(H)
    N_bfs = ZZ(len(orbits))
    print("N_orb BFS =",N_bfs)

    print("\n=== Burnside on H ===")
    N_burnside,burn_rows = burnside_H(H)
    print("N_orb Burnside(H) =",N_burnside)

    # A6
    assert N_bfs == N_burnside

    # A7 / integrality is already asserted inside Burnside.
    assert N_burnside >= 4773

    N_full = None
    if run_full_G_burnside:
        print("\n=== Independent full-G Burnside cross-check ===")
        N_full = burnside_full_G(G)
        print("N_orb Burnside(G) =",N_full)
        assert N_full == N_bfs

    print("\n=== Lazy exact geometry: once per orbit ===")
    attach_geometry(orbits,P5,P7,P11)

    stats = aggregate_mu(orbits)

    print("\n=== FINAL ===")
    print("triplet orbits =",N_bfs)
    print("support size =",stats["support_size"])
    print("P_assoc =",stats["P_assoc"])
    print("M1 =",stats["M1"])
    print("M2 =",stats["M2"])
    print("variance =",stats["variance"])

    save_results(orbits,stats,N_bfs,N_burnside,N_full)

    return {
        "G":G,
        "H":H,
        "orbits":orbits,
        "stats":stats,
        "N_bfs":N_bfs,
        "N_burnside":N_burnside,
        "N_burnside_full":N_full,
    }

if __name__ == "__main__":
    run_mu5711_production(run_full_G_burnside=True)
