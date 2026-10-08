# EABC m=1 Dimensionless Defect Audit
# Exact symbolic surrogate: K_{1,13}; replace A1 later by a genuine Hat patch adjacency.
import sympy as sp

n_leaves = 13
N = n_leaves + 1

# 1. Hat surrogate graph
A1 = sp.zeros(N)
for j in range(1, N):
    A1[0, j] = 1
    A1[j, 0] = 1

D1 = sp.diag(*[n_leaves] + [1]*n_leaves)
L1 = D1 - A1

lam = sp.symbols("lam")
print("charpoly(L1) =", sp.factor(L1.charpoly(lam).as_expr()))
print("Spec(L1) =", L1.eigenvals())

# 2. Local G(theta) deformation sector
kappa, kappa3, kappa4, g, mu = sp.symbols(
    "kappa kappa3 kappa4 g mu", real=True
)
a = kappa3 / sp.factorial(3)
b = kappa4 / sp.factorial(4)

Htheta = sp.Matrix([
    [0,     kappa, 0],
    [kappa, a,     kappa],
    [0,     kappa, b],
])

print("charpoly(Htheta) =",
      sp.factor(Htheta.charpoly(mu).as_expr()))

# 3. Fano quadratic contraction
# For fully antisymmetric octonion structure constants:
# sum_{c,d} C[a,c,d] C[b,c,d] = 6 delta_ab.
QF = 6 * sp.eye(7)
print("Spec(QF) =", QF.eigenvals())

# 4. Separable defect Hamiltonian
# H0 = kron(L1,I3,I7) + kron(I14,Htheta,I7) + g*kron(I14,I3,QF)
# Its exact spectrum is:
#     epsilon = lambda_L + mu_theta + 6*g
# with lambda_L in {0,1,14}.
#
# Multiplicities of lambda_L:
#     0  -> 1
#     1  -> 12
#     14 -> 1
# and every level is additionally 7-fold degenerate in the Fano factor.

# 5. Rank bookkeeping
rank_Pi4 = 4
d_theta = 3
rank_IF = N * d_theta * rank_Pi4
print("rank(I_N ⊗ I_3 ⊗ Pi4) =", rank_IF)  # 168, not 4

# To obtain a genuine 4D physical space use rank-1 projections
# on a chosen Hat mode and theta mode:
# P_phys = P_hat(rank 1) ⊗ P_theta(rank 1) ⊗ Pi4(rank 4).

# 6. Newton identities for a true 4x4 H_phys
s1, s2, s3, s4 = sp.symbols("s1 s2 s3 s4")
e1 = s1
e2 = (s1**2 - s2) / 2
e3 = (s1**3 - 3*s1*s2 + 2*s3) / 6
e4 = (s1**4 - 6*s1**2*s2 + 3*s2**2 + 8*s1*s3 - 6*s4) / 24

assert sp.simplify(s1 - e1) == 0
assert sp.simplify(s2 - s1*e1 + 2*e2) == 0
assert sp.simplify(s3 - s2*e1 + s1*e2 - 3*e3) == 0
assert sp.simplify(
    s4 - s3*e1 + s2*e2 - s1*e3 + 4*e4
) == 0

print("Newton identities: PASS")
