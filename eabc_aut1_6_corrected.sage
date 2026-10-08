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

load("eabc_rehm_fiber_generator.sage")

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

if __name__ == "__main__":
    run_aut1_6()
