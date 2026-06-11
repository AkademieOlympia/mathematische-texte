import numpy as np
import math

# -------------------------
# Quaternionen
# -------------------------
def quat_mul(q1, q2):
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    return (
        a1*a2 - b1*b2 - c1*c2 - d1*d2,
        a1*b2 + b1*a2 + c1*d2 - d1*c2,
        a1*c2 - b1*d2 + c1*a2 + d1*b2,
        a1*d2 + b1*c2 - c1*b2 + d1*a2,
    )

def quat_conj(q):
    return (q[0], -q[1], -q[2], -q[3])

def quat_norm(q):
    return math.sqrt(sum(x*x for x in q))

def quat_inverse(q):
    n2 = sum(x*x for x in q)
    c = quat_conj(q)
    return tuple(x / n2 for x in c)

# -------------------------
# Übergang → Kraus-Operator
# -------------------------
def transition_operator(Q1, Q2):
    delta = quat_mul(quat_inverse(Q1), Q2)
    norm = quat_norm(delta)
    
    if norm == 0:
        return np.zeros((4,4))
    
    # normierter Operator
    vec = np.array(delta) / norm
    
    # als Matrix (outer product)
    return np.outer(vec, vec)

# -------------------------
# C*-Operator D
# -------------------------
def build_Cstar_operator(quads, N=300):
    N = min(N, len(quads))
    D = np.zeros((4,4))
    
    for i in range(N - 1):
        B = transition_operator(quads[i], quads[i+1])
        D += B.T @ B   # B†B
    
    return D

# -------------------------
# Spektrum
# -------------------------
def spectrum(D):
    return np.linalg.eigvalsh(D)