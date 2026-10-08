# EABC / (5,5,5)-Vakuum — exakter Volltest mit echter Rehm-P(5)-Faser
# SageMath 10.5+
#
# Dieses Skript verwendet dieselbe Rehm/Coxeter-Fasererzeugung und
# dieselbe Fano-Multiplikationskonvention wie der bestehende EABC-Backend-Code.
# Es berechnet ALLE 126^3 = 2,000,376 geordneten Tripel exakt.
#
# Interpretation:
#   forced_assoc   = durch Alternativität erzwungener Sockel
#   assoc_distinct = zusätzliche Nullassoziatoren mit 3 verschiedenen Elementen
#   assoc_total    = gesamte Zahl assoziativer Tripel
#
# Wichtig: "Artin-Kondensation" wird hier nicht als Identität behauptet,
# sondern als universeller Alternativitäts-Sockel gezählt. Der Rest wird
# unmittelbar aus der konkreten P(5)-Faser bestimmt.

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



from collections import Counter
from itertools import product
import time

def associator_exact(x, y, z):
    return omul(omul(x, y), z) - omul(x, omul(y, z))

def exact_555_test(P5, keep_examples=10, progress_every=250000):
    n = len(P5)
    assert n == 126, f"Erwartet |P(5)|=126, erhalten {n}"
    assert len(set(tuple(x) for x in P5)) == n, "P(5) enthält doppelte Elemente."
    assert all(norm2(x) == 5 for x in P5), "Norm-5-Audit fehlgeschlagen."

    total = n**3
    forced = 3*n**2 - 2*n
    distinct_total = n*(n-1)*(n-2)

    assoc_total = ZZ(0)
    assoc_distinct = ZZ(0)
    nonassoc_distinct = ZZ(0)

    # Exakter Direktcheck der vier Alternativitätsklassen.
    forced_direct = Counter()
    examples_assoc = []
    examples_nonassoc = []

    t0 = time.time()
    done = 0

    for i, j, k in product(range(n), repeat=3):
        a, b, c = P5[i], P5[j], P5[k]
        zero = (associator_exact(a,b,c) == 0)

        if i == j == k:
            cls = "u=v=w"
        elif i == j:
            cls = "u=v"
        elif j == k:
            cls = "v=w"
        elif i == k:
            cls = "u=w"
        else:
            cls = "distinct"

        if cls != "distinct":
            # Wegen Alternativität muss der Assoziator EXAKT verschwinden.
            assert zero, ("Alternativitäts-Sockel verletzt", i, j, k, a, b, c)
            forced_direct[cls] += 1
        elif zero:
            assoc_distinct += 1
            if len(examples_assoc) < keep_examples:
                examples_assoc.append((i,j,k))
        else:
            nonassoc_distinct += 1
            if len(examples_nonassoc) < keep_examples:
                examples_nonassoc.append((i,j,k))

        if zero:
            assoc_total += 1

        done += 1
        if progress_every and done % progress_every == 0:
            print(f"  {done:,}/{total:,}  ({100*done/total:.1f} %)  {time.time()-t0:.1f}s")

    # Harte Identitäten.
    assert forced_direct["u=v=w"] == n
    assert forced_direct["u=v"] == n*(n-1)
    assert forced_direct["v=w"] == n*(n-1)
    assert forced_direct["u=w"] == n*(n-1)
    assert sum(forced_direct.values()) == forced
    assert assoc_total == forced + assoc_distinct
    assert assoc_distinct + nonassoc_distinct == distinct_total
    assert assoc_total <= total

    forced_ratio = QQ(forced)/QQ(total)
    distinct_assoc_ratio = QQ(assoc_distinct)/QQ(distinct_total)
    total_assoc_ratio = QQ(assoc_total)/QQ(total)

    print("\n" + "="*76)
    print("EABC (5,5,5)-VAKUUM — EXAKTER VOLLTEST")
    print("="*76)
    print(f"|P(5)|                              = {n}")
    print(f"Gesamtzahl Tripel                   = {total:,}")
    print()
    print("--- Universeller Alternativitäts-Sockel ---")
    print(f"u=v=w                               = {n:,}")
    print(f"genau zwei gleich                   = {3*n*(n-1):,}")
    print(f"erzwungen assoziativ                = {forced:,}")
    print(f"Sockelanteil                        = {float(forced_ratio):.8%}")
    print()
    print("--- Drei paarweise verschiedene Elemente ---")
    print(f"distinct-Tupel                      = {distinct_total:,}")
    print(f"davon zusätzlich assoziativ         = {assoc_distinct:,}")
    print(f"davon nichtassoziativ               = {nonassoc_distinct:,}")
    print(f"zusätzlicher assoc.-Anteil          = {float(distinct_assoc_ratio):.8%}")
    print()
    print("--- Gesamt ---")
    print(f"assoziative Tripel insgesamt        = {assoc_total:,}")
    print(f"Gesamt-Assoziativitätsanteil        = {float(total_assoc_ratio):.8%}")
    print(f"Überschuss über Sockel              = {assoc_distinct:,}")
    print()
    print("--- Direkt geprüfter Sockel ---")
    for key in ["u=v=w","u=v","v=w","u=w"]:
        print(f"{key:36s} = {forced_direct[key]:,}")
    print()
    print("--- Beispiel-Tripel: distinct + assoziativ ---")
    print(examples_assoc)
    print("--- Beispiel-Tripel: distinct + nichtassoziativ ---")
    print(examples_nonassoc)
    print("="*76)

    return {
        "n": n,
        "total": total,
        "forced_assoc": forced,
        "assoc_distinct": assoc_distinct,
        "assoc_total": assoc_total,
        "nonassoc_distinct": nonassoc_distinct,
        "forced_ratio": forced_ratio,
        "distinct_assoc_ratio": distinct_assoc_ratio,
        "total_assoc_ratio": total_assoc_ratio,
        "forced_direct": dict(forced_direct),
        "examples_assoc": examples_assoc,
        "examples_nonassoc": examples_nonassoc,
    }

if __name__ == "__main__":
    print("[1] Rehm/Coxeter-P(5) wird exakt erzeugt ...")
    P5 = generate_rehm_P(5)
    assert len(P5) == 126
    print("    PASS: |P(5)| = 126")
    print("    PASS: alle Normen = 5")
    print()
    print("[2] Starte vollständigen 126^3-Assoziatorlauf ...")
    result = exact_555_test(P5)
    print("\nPASS: (5,5,5)-Volltest abgeschlossen.")
