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

if __name__ == "__main__":
    print(audit_coxeter_basis())
    for p in [3,5,7]:
        Xp = generate_rehm_P(p)
        print("p =", p, "|P(p)| =", len(Xp), "target =", p**3 + 1)
