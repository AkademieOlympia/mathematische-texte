# Test Suite for the Albert Algebra H_3(O)
# Läuft mit CPython + NumPy (kein SageMath nötig).
import numpy as np

print("=== ALBERT ALGEBRA H_3(O) ALGEBRAIC GATE TEST ===")

# Fano-Tripel e_i e_j = e_k (Indizes 1..7), Standard e1 e2 = e3
_FANO_TRIPLES = (
    (1, 2, 3),
    (1, 4, 5),
    (1, 7, 6),
    (2, 4, 6),
    (2, 5, 7),
    (3, 4, 7),
    (3, 6, 5),
)
_FANO = np.zeros((8, 8), dtype=int)
for _i, _j, _k in _FANO_TRIPLES:
    _FANO[_i, _j] = _k
    _FANO[_j, _i] = -_k
    _FANO[_j, _k] = _i
    _FANO[_k, _j] = -_i
    _FANO[_k, _i] = _j
    _FANO[_i, _k] = -_j

print("[Init] Oktonionen-Basis geladen.")
print("[Init] Konstruiere Test-Matrix X in H_3(O)...")

# X = [[alpha_1, x_3, \bar{x}_2], [ \bar{x}_3, alpha_2, x_1], [x_2, \bar{x}_1, alpha_3]]


def create_sample_albert_element(seed=42):
    rng = np.random.default_rng(seed)
    alpha = rng.uniform(-2, 2, 3)
    x = [rng.uniform(-1, 1, 8) for _ in range(3)]
    return alpha, x


def octonion_norm_sq(o):
    o = np.asarray(o, dtype=float).reshape(8)
    return float(np.dot(o, o))


def octonion_conjugate(o):
    res = np.asarray(o, dtype=float).reshape(8).copy()
    res[1:] *= -1.0
    return res


def octonion_mul(o1, o2):
    a = np.asarray(o1, dtype=float).reshape(8)
    b = np.asarray(o2, dtype=float).reshape(8)
    res = np.zeros(8, dtype=float)
    res[0] = a[0] * b[0] - np.dot(a[1:], b[1:])
    res[1:] = a[0] * b[1:] + b[0] * a[1:]
    for i in range(1, 8):
        for j in range(1, 8):
            if i == j:
                continue
            target = int(_FANO[i, j])
            if target == 0:
                continue
            res[abs(target)] += (1.0 if target > 0 else -1.0) * a[i] * b[j]
    return res


def compute_cubic_norm(alpha, x):
    a1, a2, a3 = (float(v) for v in alpha)
    x1, x2, x3 = x
    n1 = a1 * a2 * a3
    n2 = (
        -a1 * octonion_norm_sq(x1)
        - a2 * octonion_norm_sq(x2)
        - a3 * octonion_norm_sq(x3)
    )
    triple_prod = octonion_mul(octonion_mul(x1, x2), x3)
    return n1 + n2 + 2.0 * float(triple_prod[0])


def quadratic_sigma(alpha, x):
    return (
        float(alpha[0] * alpha[1] + alpha[1] * alpha[2] + alpha[2] * alpha[0])
        - octonion_norm_sq(x[0])
        - octonion_norm_sq(x[1])
        - octonion_norm_sq(x[2])
    )


def albert_adjoint_matrix(X):
    # X^# = X∘X − T(X) X + σ(X) I  (kubische Jordan-Algebra)
    T = float(X[0, 0, 0] + X[1, 1, 0] + X[2, 2, 0])
    alpha = np.array([X[0, 0, 0], X[1, 1, 0], X[2, 2, 0]])
    x = [X[1, 2], X[2, 0], X[0, 1]]
    sigma = quadratic_sigma(alpha, x)
    return jordan(X, X) - T * X + sigma * ident_matrix()


def as_matrix(alpha, x):
    a1, a2, a3 = (float(v) for v in alpha)
    x1, x2, x3 = (np.asarray(xi, dtype=float).reshape(8) for xi in x)
    M = np.zeros((3, 3, 8))
    M[0, 0, 0] = a1
    M[1, 1, 0] = a2
    M[2, 2, 0] = a3
    M[1, 2] = x1
    M[2, 1] = octonion_conjugate(x1)
    M[2, 0] = x2
    M[0, 2] = octonion_conjugate(x2)
    M[0, 1] = x3
    M[1, 0] = octonion_conjugate(x3)
    return M


def mat_mul(X, Y):
    Z = np.zeros((3, 3, 8))
    for i in range(3):
        for j in range(3):
            acc = np.zeros(8)
            for k in range(3):
                acc = acc + octonion_mul(X[i, k], Y[k, j])
            Z[i, j] = acc
    return Z


def jordan(X, Y):
    return 0.5 * (mat_mul(X, Y) + mat_mul(Y, X))


def ident_matrix():
    I = np.zeros((3, 3, 8))
    I[0, 0, 0] = I[1, 1, 0] = I[2, 2, 0] = 1.0
    return I


def _pass(ok):
    return "PASS" if ok else "FAIL"


def test_algebraic_gates(alpha, x, tol=1e-10):
    dim_ok = 3 + 3 * 8 == 27

    N_X = compute_cubic_norm(alpha, x)
    norm_ok = np.isfinite(N_X)

    X = as_matrix(alpha, x)
    X_adj = albert_adjoint_matrix(X)
    I = ident_matrix()

    eps_1 = float(np.max(np.abs(jordan(X, X_adj) - N_X * I)))

    X_ss = albert_adjoint_matrix(X_adj)
    eps_2 = float(np.max(np.abs(X_ss - N_X * X)))

    T = float(np.sum(alpha))
    sigma = quadratic_sigma(alpha, x)
    X2 = jordan(X, X)
    X3 = jordan(X, X2)
    eps_ch = float(np.max(np.abs(X3 - T * X2 + sigma * X - N_X * I)))

    print("\n--- ERGEBNISSE DES ALGEBRAISCHEN GATES ---")
    print(f"Gate J1 (Dimension 27): {_pass(dim_ok)}")
    print(f"Gate J2 (Kubische Norm N(X) wohldefiniert): {_pass(norm_ok)}")
    print("Gate J3 (Aut(J)=F4, Str0(J)=E6(-26)): THEORETICAL")
    print(f"Residuum epsilon_1 (X o X^# == N(X)I): {eps_1:.3e} [{_pass(eps_1 < tol)}]")
    print(f"Residuum epsilon_2 ((X^#)^# == N(X)X): {eps_2:.3e} [{_pass(eps_2 < tol)}]")
    print(f"Gate J5 (Jordan-Cayley-Hamilton): {eps_ch:.3e} [{_pass(eps_ch < tol)}]")
    return eps_1 < tol and eps_2 < tol and eps_ch < tol and dim_ok and norm_ok


if __name__ == "__main__":
    alpha_vals, x_vals = create_sample_albert_element()
    print(f"Diagonale alpha: {alpha_vals}")
    N_X = compute_cubic_norm(alpha_vals, x_vals)
    print(f"Kubische Norm N(X) berechnet: {N_X:.6f}")
    ok = test_algebraic_gates(alpha_vals, x_vals)
    raise SystemExit(0 if ok else 1)
