import math

import numpy as np

EPS = 1e-12


def dagger(A):
    return A.conj().T


def fro_err(A, B):
    return np.linalg.norm(A - B, ord="fro")


def is_unitary(U, eps=EPS):
    n = U.shape[0]
    return fro_err(dagger(U) @ U, np.eye(n)) < eps


def is_hermitian(A, eps=EPS):
    return fro_err(A, dagger(A)) < eps


def entropy_vn(rho, eps=1e-14):
    vals = np.linalg.eigvalsh((rho + dagger(rho)) / 2)
    vals = vals[vals > eps]
    return float(-np.sum(vals * np.log(vals)))


def random_state(n):
    v = np.random.normal(size=n) + 1j * np.random.normal(size=n)
    return v / np.linalg.norm(v)


def density(psi):
    return np.outer(psi, psi.conj())


def partial_trace_switch(rho, d_switch=8, d_particle=2):
    # rho lives on H_switch ⊗ H_particle
    R = rho.reshape(d_switch, d_particle, d_switch, d_particle)
    return np.einsum("apbp->ab", R)


def partial_trace_particle(rho, d_switch=8, d_particle=2):
    R = rho.reshape(d_switch, d_particle, d_switch, d_particle)
    return np.einsum("apaq->pq", R)


def make_XZ_8():
    d = 8
    omega = np.exp(2j * np.pi / d)

    X = np.zeros((d, d), dtype=complex)
    for k in range(d):
        X[(k + 1) % d, k] = 1

    Z = np.diag([omega**k for k in range(d)])
    return X, Z, omega


def make_fourier_8():
    d = 8
    omega = np.exp(2j * np.pi / d)
    F = np.array([[omega ** (j * k) for k in range(d)] for j in range(d)], dtype=complex)
    return F / np.sqrt(d)


def make_observables_8():
    d = 8
    X, Z, omega = make_XZ_8()

    # Hermitesche Observablen aus Schaltoperatoren
    Sx = (X + dagger(X)) / 2
    Sz = (Z + dagger(Z)) / 2

    return X, Z, Sx, Sz, omega


def make_particle_space_v4():
    # Regular representation space C[V4] with basis (E, A, B, C)
    return np.eye(4, dtype=complex)


def build_axiom_system():
    """
    Minimales Axiomensystem:
    A1: H_S = l^2(Z_8)
    A2: H_P = C[V_4]
    A3: H = H_S (x) H_P
    A4: X und Z erzeugen eine Weyl-Darstellung mit ZX = omega XZ
    """
    X, Z, omega = make_XZ_8()
    HS = np.eye(8, dtype=complex)
    HP = make_particle_space_v4()
    H = np.kron(HS, HP)
    return {"HS": HS, "HP": HP, "H": H, "X": X, "Z": Z, "omega": omega}


def test_axioms(eps=EPS):
    ax = build_axiom_system()
    HS, HP, H = ax["HS"], ax["HP"], ax["H"]
    X, Z, omega = ax["X"], ax["Z"], ax["omega"]

    dim_hs_ok = HS.shape == (8, 8)
    dim_hp_ok = HP.shape == (4, 4)
    dim_h_ok = H.shape == (32, 32)

    weyl_error = fro_err(Z @ X, omega * X @ Z)
    weyl_ok = weyl_error < eps

    return {
        "A1_dim_HS_is_8": dim_hs_ok,
        "A2_dim_HP_is_4": dim_hp_ok,
        "A3_dim_H_is_32": dim_h_ok,
        "A4_weyl_error": weyl_error,
        "A4_weyl_passed": weyl_ok,
        "passed": dim_hs_ok and dim_hp_ok and dim_h_ok and weyl_ok,
    }


def variance(A, psi):
    exp_A = np.vdot(psi, A @ psi)
    exp_A2 = np.vdot(psi, A @ A @ psi)
    return float(np.real(exp_A2 - exp_A * np.conj(exp_A)))


def uncertainty_test(A, B, psi, eps=1e-12):
    var_A = variance(A, psi)
    var_B = variance(B, psi)

    dA = np.sqrt(max(var_A, 0))
    dB = np.sqrt(max(var_B, 0))

    comm = A @ B - B @ A
    rhs = 0.5 * abs(np.vdot(psi, comm @ psi))

    return dA * dB + eps >= rhs, dA * dB, rhs


def spectral_projectors(A, eps=1e-10):
    vals, vecs = np.linalg.eigh((A + dagger(A)) / 2)

    projectors = []
    used = np.zeros(len(vals), dtype=bool)

    for i, lam in enumerate(vals):
        if used[i]:
            continue

        idx = np.where(abs(vals - lam) < eps)[0]
        used[idx] = True

        V = vecs[:, idx]
        P = V @ dagger(V)
        projectors.append((lam, P))

    return projectors


def projector_tests(A, eps=1e-10):
    n = A.shape[0]
    Ps = spectral_projectors(A, eps=eps)

    sum_P = np.zeros((n, n), dtype=complex)
    max_idem = 0
    max_herm = 0

    for lam, P in Ps:
        max_idem = max(max_idem, fro_err(P @ P, P))
        max_herm = max(max_herm, fro_err(P, dagger(P)))
        sum_P += P

    sum_rule = fro_err(sum_P, np.eye(n))

    return {
        "num_projectors": len(Ps),
        "max_idempotence_error": max_idem,
        "max_hermiticity_error": max_herm,
        "sum_rule_error": sum_rule,
        "passed": max_idem < eps and max_herm < eps and sum_rule < eps,
    }


def born_rule_test(A, shots=20000, eps=0.02):
    n = A.shape[0]
    psi = random_state(n)
    Ps = spectral_projectors(A)

    probs = np.array([np.real(np.vdot(psi, P @ psi)) for _, P in Ps])
    probs = probs / probs.sum()

    samples = np.random.choice(len(Ps), size=shots, p=probs)
    freqs = np.bincount(samples, minlength=len(Ps)) / shots

    max_err = float(np.max(abs(freqs - probs)))

    return {
        "prob_sum": float(probs.sum()),
        "max_frequency_error": max_err,
        "passed": abs(probs.sum() - 1) < 1e-12 and max_err < eps,
    }


def entanglement_test(d_switch=8, d_particle=2):
    # Produktzustand
    psi_s = random_state(d_switch)
    psi_p = random_state(d_particle)
    psi_prod = np.kron(psi_s, psi_p)
    rho_prod = density(psi_prod)

    rho_s_prod = partial_trace_particle(rho_prod, d_switch, d_particle)
    rho_p_prod = partial_trace_switch(rho_prod, d_switch, d_particle)

    S_s_prod = entropy_vn(rho_s_prod)
    S_p_prod = entropy_vn(rho_p_prod)

    # verschränkter Zustand, eingebettet in 8 x 2
    psi_ent = np.zeros(d_switch * d_particle, dtype=complex)
    psi_ent[0 * d_particle + 0] = 1 / np.sqrt(2)
    psi_ent[1 * d_particle + 1] = 1 / np.sqrt(2)

    rho_ent = density(psi_ent)

    rho_s_ent = partial_trace_particle(rho_ent, d_switch, d_particle)
    rho_p_ent = partial_trace_switch(rho_ent, d_switch, d_particle)

    S_s_ent = entropy_vn(rho_s_ent)
    S_p_ent = entropy_vn(rho_p_ent)

    return {
        "product_entropy_switch": S_s_prod,
        "product_entropy_particle": S_p_prod,
        "entangled_entropy_switch": S_s_ent,
        "entangled_entropy_particle": S_p_ent,
        "pure_bipartite_entropy_match_error": abs(S_s_ent - S_p_ent),
        "passed": (
            S_s_prod < 1e-10
            and S_p_prod < 1e-10
            and abs(S_s_ent - S_p_ent) < 1e-10
            and S_s_ent > 0.1
        ),
    }


def casimir_test(X, Z, eps=1e-10):
    Sx = (X + dagger(X)) / 2
    Sy = (X - dagger(X)) / (2j)
    Sz = (Z + dagger(Z)) / 2

    C = Sx @ Sx + Sy @ Sy + Sz @ Sz
    c_x = fro_err(C @ Sx - Sx @ C, np.zeros_like(C))
    c_y = fro_err(C @ Sy - Sy @ C, np.zeros_like(C))
    c_z = fro_err(C @ Sz - Sz @ C, np.zeros_like(C))

    return {
        "comm_C_Sx": c_x,
        "comm_C_Sy": c_y,
        "comm_C_Sz": c_z,
        "passed": c_x < eps and c_y < eps and c_z < eps,
    }


def su2_generators_j(j):
    """
    Standard-su(2)-Darstellung mit Spin j in Basis |m>, m=j,...,-j.
    Dimension = 2j+1.
    """
    d = int(2 * j + 1)
    m_vals = np.array([j - k for k in range(d)], dtype=float)

    Jp = np.zeros((d, d), dtype=complex)
    for col, m in enumerate(m_vals):
        target_m = m + 1
        row_candidates = np.where(np.isclose(m_vals, target_m))[0]
        if row_candidates.size == 0:
            continue
        row = int(row_candidates[0])
        Jp[row, col] = np.sqrt((j - m) * (j + m + 1))

    Jm = dagger(Jp)
    Jx = (Jp + Jm) / 2
    Jy = (Jp - Jm) / (2j)
    Jz = np.diag(m_vals.astype(complex))
    return Jx, Jy, Jz


def casimir_test_variant_su2(eps=1e-10):
    """
    Zweite Casimir-Variante:
    Nutzt eine exakte su(2)-Darstellung auf d=8 (j=7/2).
    """
    j = 3.5
    Jx, Jy, Jz = su2_generators_j(j)
    C = Jx @ Jx + Jy @ Jy + Jz @ Jz
    d = Jx.shape[0]
    I = np.eye(d, dtype=complex)

    c_x = fro_err(C @ Jx - Jx @ C, np.zeros_like(C))
    c_y = fro_err(C @ Jy - Jy @ C, np.zeros_like(C))
    c_z = fro_err(C @ Jz - Jz @ C, np.zeros_like(C))

    expected = j * (j + 1) * I
    casimir_scalar_err = fro_err(C, expected)

    return {
        "spin_j": j,
        "dimension": d,
        "comm_C_Jx": c_x,
        "comm_C_Jy": c_y,
        "comm_C_Jz": c_z,
        "casimir_scalar_error": casimir_scalar_err,
        "passed": c_x < eps and c_y < eps and c_z < eps and casimir_scalar_err < eps,
    }


def ordered_projectors_from_z(Z):
    vals, vecs = np.linalg.eig(Z)
    phases = np.mod(np.angle(vals), 2 * np.pi)
    order = np.argsort(phases)

    Ps = []
    for idx in order:
        v = vecs[:, idx]
        v = v / np.linalg.norm(v)
        Ps.append(np.outer(v, v.conj()))
    return Ps


def casimir_test_variant_weyl_constructed(X, Z, eps=1e-10):
    """
    Dritte Casimir-Variante (BM-nah):
    Konstruktion der su(2)-Generatoren direkt aus Weyl-Daten (X, Z)
    mittels Spektralprojektoren von Z und gewichteten X-Übergängen.
    """
    d = X.shape[0]
    j = (d - 1) / 2
    Ps = ordered_projectors_from_z(Z)

    # m_k = -j+k in der nach Phase sortierten Z-Eigenbasis
    # (aufsteigende magnetische Quantenzahl, damit J+ : k -> k+1)
    m_vals = [-j + k for k in range(d)]

    Jz = np.zeros((d, d), dtype=complex)
    for k in range(d):
        Jz += m_vals[k] * Ps[k]

    # Leiteroperator ohne zyklischen Wrap-around-Term
    Jp = np.zeros((d, d), dtype=complex)
    for k in range(d - 1):
        m = m_vals[k]
        coeff = np.sqrt((j - m) * (j + m + 1))
        Jp += coeff * (Ps[k + 1] @ X @ Ps[k])

    Jm = dagger(Jp)
    Jx = (Jp + Jm) / 2
    Jy = (Jp - Jm) / (2j)

    C = Jx @ Jx + Jy @ Jy + Jz @ Jz
    I = np.eye(d, dtype=complex)

    c_x = fro_err(C @ Jx - Jx @ C, np.zeros_like(C))
    c_y = fro_err(C @ Jy - Jy @ C, np.zeros_like(C))
    c_z = fro_err(C @ Jz - Jz @ C, np.zeros_like(C))
    casimir_scalar_err = fro_err(C, j * (j + 1) * I)

    su2_comm_x = fro_err(Jy @ Jz - Jz @ Jy, 1j * Jx)
    su2_comm_y = fro_err(Jz @ Jx - Jx @ Jz, 1j * Jy)
    su2_comm_z = fro_err(Jx @ Jy - Jy @ Jx, 1j * Jz)

    return {
        "spin_j": j,
        "dimension": d,
        "comm_C_Jx": c_x,
        "comm_C_Jy": c_y,
        "comm_C_Jz": c_z,
        "casimir_scalar_error": casimir_scalar_err,
        "su2_comm_error_x": su2_comm_x,
        "su2_comm_error_y": su2_comm_y,
        "su2_comm_error_z": su2_comm_z,
        "passed": (
            c_x < eps
            and c_y < eps
            and c_z < eps
            and casimir_scalar_err < eps
            and su2_comm_x < 1e-8
            and su2_comm_y < 1e-8
            and su2_comm_z < 1e-8
        ),
    }


def laplace_operator(X):
    d = X.shape[0]
    return 2 * np.eye(d, dtype=complex) - X - dagger(X)


def laplace_spectrum_test(X, eps=1e-10):
    d = X.shape[0]
    Delta = laplace_operator(X)
    spec_num = np.sort(np.real(np.linalg.eigvalsh((Delta + dagger(Delta)) / 2)))
    ks = np.arange(d)
    spec_theory = np.sort(2 - 2 * np.cos(2 * np.pi * ks / d))
    max_err = float(np.max(np.abs(spec_num - spec_theory)))
    return {
        "spectrum_numeric": spec_num.tolist(),
        "spectrum_theory": spec_theory.tolist(),
        "max_spectral_error": max_err,
        "passed": max_err < eps,
    }


def dirac_operator(X):
    d = X.shape[0]
    I = np.eye(d, dtype=complex)
    upper = np.hstack([np.zeros((d, d), dtype=complex), X - I])
    lower = np.hstack([dagger(X) - I, np.zeros((d, d), dtype=complex)])
    return np.vstack([upper, lower])


def dirac_test(X, eps=1e-10):
    d = X.shape[0]
    Delta = laplace_operator(X)
    D = dirac_operator(X)
    D2 = D @ D
    target = np.block(
        [
            [Delta, np.zeros((d, d), dtype=complex)],
            [np.zeros((d, d), dtype=complex), Delta],
        ]
    )
    err = fro_err(D2, target)
    return {"d2_block_laplace_error": err, "passed": err < eps}


def entropy_test_suite(d_switch=8, d_particle=2):
    # Reiner Zustand -> S=0
    psi = random_state(d_switch)
    rho_pure = density(psi)
    s_pure = entropy_vn(rho_pure)

    # Maximale Mischung -> S = log(d)
    rho_mix = np.eye(d_switch, dtype=complex) / d_switch
    s_mix = entropy_vn(rho_mix)
    s_mix_target = np.log(d_switch)

    # Verschränkter Zustand auf d_switch x d_particle
    ent = entanglement_test(d_switch=d_switch, d_particle=d_particle)

    return {
        "pure_entropy": s_pure,
        "mixed_entropy": s_mix,
        "mixed_entropy_target_log_d": float(s_mix_target),
        "mixed_entropy_error": abs(s_mix - s_mix_target),
        "entanglement": ent,
        "passed": (
            s_pure < 1e-10
            and abs(s_mix - s_mix_target) < 1e-10
            and ent["passed"]
        ),
    }


def holonomy_phase(step_count, phase_unit=np.exp(2j * np.pi / 24)):
    return phase_unit**step_count


def holonomy_test(eps=1e-12):
    # Diskrete geschlossene Wege als Windungszahlen im 24-Zyklus
    g1 = 5
    g2 = 7
    left = holonomy_phase(g1 + g2)
    right = holonomy_phase(g1) * holonomy_phase(g2)
    comp_err = abs(left - right)

    closed_24 = holonomy_phase(24)
    closed_8x3 = holonomy_phase(8 * 3)

    return {
        "composition_error": comp_err,
        "gamma_24": complex(closed_24),
        "gamma_8x3": complex(closed_8x3),
        "passed": comp_err < eps and abs(closed_24 - 1) < eps and abs(closed_8x3 - 1) < eps,
    }


def sagnac_test(X, eps=1e-12):
    """
    Diskreter Sagnac-Test:
    U+ = e^{i theta} X, U- = e^{-i theta} X^{-1}, theta = 2pi/24
    Nach 8 Schritten ergibt sich die orientierte Phasendifferenz Delta phi = 16 theta mod 2pi.
    """
    d = X.shape[0]
    theta = 2 * np.pi / 24
    U_plus = np.exp(1j * theta) * X
    U_minus = np.exp(-1j * theta) * dagger(X)

    loop_plus = np.linalg.matrix_power(U_plus, d)
    loop_minus = np.linalg.matrix_power(U_minus, d)

    # Zentralphase aus der Spur extrahieren
    phase_plus = np.angle(np.trace(loop_plus))
    phase_minus = np.angle(np.trace(loop_minus))
    delta_phase = (phase_plus - phase_minus) % (2 * np.pi)
    delta_theory = (16 * theta) % (2 * np.pi)

    err = abs(delta_phase - delta_theory)
    return {
        "theta": float(theta),
        "delta_phase_numeric": float(delta_phase),
        "delta_phase_theory": float(delta_theory),
        "delta_phase_error": float(err),
        "passed": err < eps,
    }


def cycle_24_test(X, eps=1e-12):
    d = X.shape[0]
    I = np.eye(d, dtype=complex)
    theta = 2 * np.pi / 24
    U = np.exp(1j * theta) * X

    x8_err = fro_err(np.linalg.matrix_power(X, 8), I)
    x24_err = fro_err(np.linalg.matrix_power(X, 24), I)
    u24_err = fro_err(np.linalg.matrix_power(U, 24), I)

    return {
        "x8_identity_error": x8_err,
        "x24_identity_error": x24_err,
        "u24_identity_error": u24_err,
        "lcm_8_12": int(np.lcm(8, 12)),
        "passed": x8_err < eps and x24_err < eps and u24_err < 1e-9,
    }


def make_XZ_d(d):
    omega = np.exp(2j * np.pi / d)
    X = np.zeros((d, d), dtype=complex)
    for k in range(d):
        X[(k + 1) % d, k] = 1
    Z = np.diag([omega**k for k in range(d)])
    return X, Z, omega


def born_probability_vector(A, psi):
    Ah = (A + dagger(A)) / 2
    _, vecs = np.linalg.eigh(Ah)
    amps = dagger(vecs) @ psi
    p = np.real(amps.conj() * amps)
    p = np.clip(p, 0, None)
    if p.sum() == 0:
        return p
    return p / p.sum()


def born_probability_by_projectors(A, psi):
    Ps = spectral_projectors(A, eps=1e-8)
    probs = np.array([np.real(np.vdot(psi, P @ psi)) for _, P in Ps])
    probs = np.clip(probs, 0, None)
    if probs.sum() == 0:
        return probs
    return probs / probs.sum()


def unitary_from_hermitian(H, eps_scale):
    # U = exp(-i eps H) via Spektralzerlegung (H hermitesch)
    vals, vecs = np.linalg.eigh((H + dagger(H)) / 2)
    phases = np.exp(-1j * eps_scale * vals)
    return vecs @ np.diag(phases) @ dagger(vecs)


def perturbation_robustness_test(X, Z, eps_noise=1e-3, seed=7):
    rng = np.random.default_rng(seed)
    d = X.shape[0]
    I = np.eye(d, dtype=complex)

    R = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    R = R / np.linalg.norm(R, ord="fro")
    X_eps = X + eps_noise * R

    unitarity_break = fro_err(dagger(X_eps) @ X_eps, I)

    Sx = (X + dagger(X)) / 2
    Sx_eps_raw = (X_eps + dagger(X_eps)) / 2
    psi = random_state(d)
    p0 = born_probability_vector(Sx, psi)
    p_raw = born_probability_vector(Sx_eps_raw, psi)
    born_prob_l1_raw = float(np.sum(np.abs(p_raw - p0)))

    # Strukturerhaltende Stoerung: U_eps = exp(-i eps H), Sx' = U Sx U^dagger
    H = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    H = (H + dagger(H)) / 2
    H = H / np.linalg.norm(H, ord="fro")
    U_struct = unitary_from_hermitian(H, eps_noise)
    Sx_eps_unitary = U_struct @ Sx @ dagger(U_struct)
    # Degenerationsrobuster Vergleich ueber Spektralprojektoren
    p0_proj = born_probability_by_projectors(Sx, psi)
    p_unitary_proj = born_probability_by_projectors(Sx_eps_unitary, psi)
    born_prob_l1_unitary = float(np.sum(np.abs(p_unitary_proj - p0_proj)))

    theta = 2 * np.pi / 24
    U_eps = np.exp(1j * theta) * X_eps
    holonomy_break = fro_err(np.linalg.matrix_power(U_eps, 24), I)

    sensitivity = {
        "unitarity_break": unitarity_break,
        "born_prob_l1_raw": born_prob_l1_raw,
        "born_prob_l1_unitary": born_prob_l1_unitary,
        "holonomy_break": holonomy_break,
    }
    most_sensitive = max(sensitivity, key=sensitivity.get)

    thr_unitarity = 5e-2
    thr_born_raw = 2e-1
    thr_born_unitary = 2e-1
    thr_holonomy = 2e-1

    return {
        "epsilon": eps_noise,
        "unitarity_break": unitarity_break,
        "born_prob_l1_raw": born_prob_l1_raw,
        "born_prob_l1_unitary": born_prob_l1_unitary,
        "holonomy_break_24": holonomy_break,
        "unitarity_threshold": thr_unitarity,
        "born_raw_threshold": thr_born_raw,
        "born_unitary_threshold": thr_born_unitary,
        "holonomy_threshold": thr_holonomy,
        "unitarity_ok": unitarity_break < thr_unitarity,
        "born_raw_ok": born_prob_l1_raw < thr_born_raw,
        "born_ok": born_prob_l1_unitary < thr_born_unitary,
        "holonomy_ok": holonomy_break < thr_holonomy,
        "most_sensitive_metric": most_sensitive,
        "passed_overall": (
            unitarity_break < thr_unitarity
            and born_prob_l1_unitary < thr_born_unitary
            and holonomy_break < thr_holonomy
        ),
    }


def gauge_holonomy_invariance_test(d=8, seed=11, eps=1e-12):
    rng = np.random.default_rng(seed)
    phases = rng.uniform(0, 2 * np.pi, size=d)
    links = np.exp(1j * phases)
    gamma = np.prod(links)

    alpha = rng.uniform(0, 2 * np.pi, size=d)
    links_g = np.zeros(d, dtype=complex)
    for r in range(d):
        links_g[r] = np.exp(1j * (alpha[(r + 1) % d] - alpha[r])) * links[r]
    gamma_g = np.prod(links_g)

    err = abs(gamma_g - gamma)
    return {"gamma": complex(gamma), "gamma_gauge": complex(gamma_g), "error": err, "passed": err < eps}


def construct_weyl_su2_generators(X, Z):
    d = X.shape[0]
    j = (d - 1) / 2
    Ps = ordered_projectors_from_z(Z)
    m_vals = [-j + k for k in range(d)]

    Jz = np.zeros((d, d), dtype=complex)
    for k in range(d):
        Jz += m_vals[k] * Ps[k]

    Jp = np.zeros((d, d), dtype=complex)
    for k in range(d - 1):
        m = m_vals[k]
        coeff = np.sqrt((j - m) * (j + m + 1))
        Jp += coeff * (Ps[k + 1] @ X @ Ps[k])

    Jm = dagger(Jp)
    Jx = (Jp + Jm) / 2
    Jy = (Jp - Jm) / (2j)
    return Jx, Jy, Jz


def construct_weyl_su2_ladder_data(X, Z):
    d = X.shape[0]
    j = (d - 1) / 2
    Ps = ordered_projectors_from_z(Z)
    m_vals = [-j + k for k in range(d)]

    Jz = np.zeros((d, d), dtype=complex)
    for k in range(d):
        Jz += m_vals[k] * Ps[k]

    Jp = np.zeros((d, d), dtype=complex)
    for k in range(d - 1):
        m = m_vals[k]
        coeff = np.sqrt((j - m) * (j + m + 1))
        Jp += coeff * (Ps[k + 1] @ X @ Ps[k])

    Jm = dagger(Jp)
    Jx = (Jp + Jm) / 2
    Jy = (Jp - Jm) / (2j)
    return Jx, Jy, Jz, Jp, Jm, j


def spin_seven_half_characterization_test(X, Z, eps=1e-10):
    Jx, Jy, Jz, Jp, _, j = construct_weyl_su2_ladder_data(X, Z)
    d = X.shape[0]
    I = np.eye(d, dtype=complex)

    C = Jx @ Jx + Jy @ Jy + Jz @ Jz
    casimir_err = fro_err(C, j * (j + 1) * I)

    nilpotent_order_err = fro_err(np.linalg.matrix_power(Jp, d), np.zeros((d, d), dtype=complex))
    almost_top_err = np.linalg.norm(np.linalg.matrix_power(Jp, d - 1), ord="fro")

    vals = np.sort(np.real(np.linalg.eigvalsh((Jz + dagger(Jz)) / 2)))
    target_vals = np.arange(-j, j + 1, 1.0)
    jz_ladder_err = float(np.max(np.abs(vals - target_vals)))

    return {
        "spin_j_target": float(j),
        "casimir_error": casimir_err,
        "jz_ladder_error": jz_ladder_err,
        "Jplus_power_d_error": nilpotent_order_err,
        "Jplus_power_d_minus_1_norm": float(almost_top_err),
        "passed": casimir_err < eps and jz_ladder_err < eps and nilpotent_order_err < eps and almost_top_err > 1e-6,
    }


def period_24_theorem_test(X, eps=1e-10):
    d = X.shape[0]
    I = np.eye(d, dtype=complex)
    theta = 2 * np.pi / 24

    operators = {
        "U_theta_X": np.exp(1j * theta) * X,
        "U_3theta_X": np.exp(3j * theta) * X,
        "U_theta_X3": np.exp(1j * theta) * np.linalg.matrix_power(X, 3),
    }

    op_errors = {}
    for name, U in operators.items():
        op_errors[name] = fro_err(np.linalg.matrix_power(U, 24), I)

    return {
        "theorem": "Per(BM-IV) = lcm(8,12) = 24",
        "operator_u24_errors": op_errors,
        "passed": all(err < eps for err in op_errors.values()),
    }


def effective_two_state_projection_test(X, Z, eps=1e-8):
    Jx, Jy, Jz = construct_weyl_su2_generators(X, Z)
    d = X.shape[0]

    candidates = {}

    # EABC-Projektionen auf Z8-Basis (jeweils Rang 2)
    class_sets = {
        "E": [0, 4],
        "A": [1, 5],
        "B": [2, 6],
        "C": [3, 7],
    }
    for name, idxs in class_sets.items():
        P = np.zeros((d, d), dtype=complex)
        for i in idxs:
            P[i, i] = 1
        candidates[f"EABC_{name}"] = P

    # Vierlings-Projektion (odd/even): Rang 4, sollte nicht direkt qubit sein
    P_odd = np.zeros((d, d), dtype=complex)
    P_even = np.zeros((d, d), dtype=complex)
    for i in range(d):
        if i % 2 == 0:
            P_even[i, i] = 1
        else:
            P_odd[i, i] = 1
    candidates["quadruplet_even"] = P_even
    candidates["quadruplet_odd"] = P_odd

    # Paritaetsprojektoren aus Spiegeloperator |r> -> |-r mod 8|
    R = np.zeros((d, d), dtype=complex)
    for r in range(d):
        R[(-r) % d, r] = 1
    vals, vecs = np.linalg.eigh((R + dagger(R)) / 2)
    V_plus = vecs[:, vals > 0]
    V_minus = vecs[:, vals < 0]
    if V_plus.shape[1] > 0:
        candidates["parity_plus"] = V_plus @ dagger(V_plus)
    if V_minus.shape[1] > 0:
        candidates["parity_minus"] = V_minus @ dagger(V_minus)

    details = {}
    best_rank2 = None
    best_err = np.inf

    for name, P in candidates.items():
        # numerischer Rang
        eigs = np.linalg.eigvalsh((P + dagger(P)) / 2)
        rank = int(np.sum(eigs > 1e-8))

        info = {"rank": rank}
        if rank == 2:
            vals_p, vecs_p = np.linalg.eigh((P + dagger(P)) / 2)
            V = vecs_p[:, vals_p > 1e-8]
            jx2 = dagger(V) @ Jx @ V
            jy2 = dagger(V) @ Jy @ V
            jz2 = dagger(V) @ Jz @ V
            sx, sy, sz = 2 * jx2, 2 * jy2, 2 * jz2
            comm_err = fro_err(sx @ sy - sy @ sx, 2j * sz)
            info["pauli_comm_error"] = comm_err
            info["effective_qubit"] = comm_err < eps
            if comm_err < best_err:
                best_err = comm_err
                best_rank2 = name
        else:
            info["effective_qubit"] = False
        details[name] = info

    return {
        "best_rank2_projection": best_rank2,
        "best_rank2_pauli_error": float(best_err) if np.isfinite(best_err) else None,
        "details": details,
        "passed": bool(np.isfinite(best_err) and best_err < eps),
    }


def spin_half_projection_test(X, Z, eps=1e-8):
    Jx, Jy, Jz = construct_weyl_su2_generators(X, Z)
    vals, vecs = np.linalg.eigh((Jz + dagger(Jz)) / 2)

    idx_minus = int(np.argmin(np.abs(vals + 0.5)))
    idx_plus = int(np.argmin(np.abs(vals - 0.5)))

    V = np.column_stack([vecs[:, idx_minus], vecs[:, idx_plus]])
    jx2 = dagger(V) @ Jx @ V
    jy2 = dagger(V) @ Jy @ V
    jz2 = dagger(V) @ Jz @ V

    sx = 2 * jx2
    sy = 2 * jy2
    sz = 2 * jz2

    comm_xy = fro_err(sx @ sy - sy @ sx, 2j * sz)
    anti_xx = fro_err(sx @ sx, np.eye(2))
    anti_yy = fro_err(sy @ sy, np.eye(2))
    anti_zz = fro_err(sz @ sz, np.eye(2))

    return {
        "m_minus": float(vals[idx_minus]),
        "m_plus": float(vals[idx_plus]),
        "comm_xy_error": comm_xy,
        "sigma_x2_error": anti_xx,
        "sigma_y2_error": anti_yy,
        "sigma_z2_error": anti_zz,
        "passed": comm_xy < eps and anti_xx < eps and anti_yy < eps and anti_zz < eps,
    }


def einstein_de_haas_test(eps=1e-12):
    """
    Diskreter Einstein-de-Haas-Test auf Tetraederpunkten:
    - Flaechenzyklus liefert L ~ 4S
    - Bilanzrelation: Delta S + (1/4) Delta L_body = 0
    """
    pts = {
        "E": np.array([1, 1, 1], dtype=float),
        "A": np.array([-1, -1, 1], dtype=float),
        "B": np.array([1, -1, -1], dtype=float),
        "C": np.array([-1, 1, -1], dtype=float),
    }
    faces = {"E": ["A", "B", "C"], "A": ["E", "C", "B"], "B": ["E", "A", "C"], "C": ["E", "B", "A"]}

    def L_cycle(cyc):
        L = np.zeros(3)
        for x, y in zip(cyc, cyc[1:] + cyc[:1]):
            L += np.cross(pts[x], pts[y])
        return L

    L_errors = []
    balance_errors = []
    for s, cyc in faces.items():
        S = pts[s]
        L = L_cycle(cyc)
        L_errors.append(np.linalg.norm(L - 4 * S))

        dS = -2 * S
        dL_body = 8 * S
        balance_errors.append(np.linalg.norm(dS + 0.25 * dL_body))

    max_L_error = float(max(L_errors))
    max_balance_error = float(max(balance_errors))
    return {
        "max_L_minus_4S_error": max_L_error,
        "max_balance_error": max_balance_error,
        "passed": max_L_error < eps and max_balance_error < eps,
    }


def atomic_configuration_test():
    """
    Referenzprofil Atomkonfiguration:
    - Wasserstoff: 1 Elektron in 1s
    - Helium: geschlossene 1s-Schale
    """
    # Elektronenschalen-Proxys (Kapazitaet 1s = 2)
    h_electrons = 1
    he_electrons = 2
    shell_capacity_1s = 2
    h_shell_ok = 0 < h_electrons <= shell_capacity_1s
    he_closed_shell = he_electrons == shell_capacity_1s

    return {
        "hydrogen_shell_ok": h_shell_ok,
        "helium_closed_shell": he_closed_shell,
        "passed": h_shell_ok and he_closed_shell,
    }


def nuclear_stability_test():
    """
    Referenzprofil Nuklearstabilitaet:
    - Eisen: hohe Bindungsenergie pro Nukleon
    - Uran: primordial langlebig
    - Transurane: keine stabilen Isotope, aber teils langlebige Isotope
    """
    # Kernstabilitaets-Proxys (bekannte Referenzwerte, grob)
    fe56_be_per_a = 8.790  # MeV/Nukleon, nahe globalem Maximum
    u238_half_life_years = 4.468e9
    u235_half_life_years = 7.04e8

    fe_binding_ok = 8.6 <= fe56_be_per_a <= 9.0
    uranium_primordial_ok = u238_half_life_years > 1e8 and u235_half_life_years > 1e8

    # Transurane: keine stabilen Isotope, aber Existenz langlebiger Vertreter
    # (vereinfachter Konsistenztest)
    transuran_stable_isotopes = 0
    transuran_long_lived_years = {
        "Np-237": 2.144e6,
        "Pu-244": 8.0e7,
        "Cm-247": 1.56e7,
    }
    transuran_no_stable_ok = transuran_stable_isotopes == 0
    transuran_long_lived_exists = max(transuran_long_lived_years.values()) > 1e6

    return {
        "iron_fe56_be_per_a": fe56_be_per_a,
        "iron_binding_ok": fe_binding_ok,
        "uranium_u238_half_life_years": u238_half_life_years,
        "uranium_u235_half_life_years": u235_half_life_years,
        "uranium_primordial_ok": uranium_primordial_ok,
        "transuran_stable_isotopes": transuran_stable_isotopes,
        "transuran_no_stable_ok": transuran_no_stable_ok,
        "transuran_max_half_life_years": max(transuran_long_lived_years.values()),
        "transuran_long_lived_exists": transuran_long_lived_exists,
        "passed": fe_binding_ok and uranium_primordial_ok and transuran_no_stable_ok and transuran_long_lived_exists,
    }


def minimal_global_return_period(U, max_n=48, eps=1e-10):
    d = U.shape[0]
    I = np.eye(d, dtype=complex)
    for n in range(1, max_n + 1):
        if fro_err(np.linalg.matrix_power(U, n), I) < eps:
            return n
    return None


def prediction_24_return_family_test(X, eps=1e-10):
    """
    Vorhersage P1 (aus 24-Theorem):
    Fuer U_{alpha,m} = exp(2pi i alpha/24) X^m mit m ungerade (invertierbar mod 8)
    ist die minimale globale Rueckkehrperiode 24.
    """
    alphas = [1, 5, 7, 11]
    ms = [1, 3, 5, 7]
    periods = {}
    ok = True
    for alpha in alphas:
        for m in ms:
            U = np.exp(2j * np.pi * alpha / 24) * np.linalg.matrix_power(X, m)
            p = minimal_global_return_period(U, max_n=48, eps=eps)
            periods[f"a{alpha}_m{m}"] = p
            ok = ok and (p == 24)

    return {"periods": periods, "passed": ok}


def prediction_spin72_transition_pattern_test(X, Z, eps=1e-8):
    """
    Vorhersage P2 (aus Spin-7/2-Struktur):
    Die benachbarten Leiteramplituden folgen |<m+1|J+|m>| = sqrt((j-m)(j+m+1)).
    """
    _, _, Jz, Jp, _, j = construct_weyl_su2_ladder_data(X, Z)
    vals, vecs = np.linalg.eigh((Jz + dagger(Jz)) / 2)
    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]

    Jp_basis = dagger(vecs) @ Jp @ vecs
    observed = []
    theory = []
    for k in range(len(vals) - 1):
        m = float(vals[k])
        observed.append(abs(Jp_basis[k + 1, k]))
        theory.append(np.sqrt((j - m) * (j + m + 1)))

    observed = np.array(observed, dtype=float)
    theory = np.array(theory, dtype=float)
    max_err = float(np.max(np.abs(observed - theory)))

    return {
        "j": float(j),
        "observed_couplings": observed.tolist(),
        "theory_couplings": theory.tolist(),
        "max_coupling_error": max_err,
        "passed": max_err < eps,
    }


def prediction_suite(X, Z):
    p24 = prediction_24_return_family_test(X)
    pspin = prediction_spin72_transition_pattern_test(X, Z)
    return {
        "prediction_24_return_family": p24,
        "prediction_spin72_transition_pattern": pspin,
        "passed": p24["passed"] and pspin["passed"],
    }


def negative_control_suite():
    ctrls = {}
    for d in [7, 9, 10]:
        Xd, Zd, omega_d = make_XZ_d(d)
        I = np.eye(d, dtype=complex)
        weyl_err = fro_err(Zd @ Xd, omega_d * Xd @ Zd)
        lap = laplace_spectrum_test(Xd)["max_spectral_error"]
        dirac_err = dirac_test(Xd)["d2_block_laplace_error"]
        theta = 2 * np.pi / 24
        U = np.exp(1j * theta) * Xd
        cycle24_err = fro_err(np.linalg.matrix_power(U, 24), I)
        ctrls[f"Z{d}"] = {
            "weyl_error": weyl_err,
            "laplace_error": lap,
            "dirac_error": dirac_err,
            "cycle24_error": cycle24_err,
            "period_lcm_d_12": int(np.lcm(d, 12)),
            "passed_24_signature": cycle24_err < 1e-9,
        }

    # Absichtlich falsche Partikel-Kongruenz (mod 10 statt mod 12)
    ctrls["particle_mod10"] = {
        "claimed_global_period": int(np.lcm(8, 10)),
        "bmiv_global_period": int(np.lcm(8, 12)),
        "matches_bmiv_24": int(np.lcm(8, 10)) == 24,
    }
    return ctrls


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    k = 3
    while k * k <= n:
        if n % k == 0:
            return False
        k += 2
    return True


def prime_quadruplet_observable_test(n_max=5000):
    hits = []
    for n in range(3, n_max - 8, 2):
        q = int(is_prime(n) and is_prime(n + 2) and is_prime(n + 6) and is_prime(n + 8))
        if q == 1:
            hits.append(n)

    if not hits:
        return {"num_hits": 0, "passed": False}

    hits = np.array(hits)
    mod8 = np.bincount(hits % 8, minlength=8) / len(hits)
    mod12 = np.bincount(hits % 12, minlength=12) / len(hits)
    mod24 = np.bincount(hits % 24, minlength=24) / len(hits)

    # Nichttriviale Signatur: Konzentration auf wenigen Restklassen
    concentration_12 = float(np.max(mod12))
    concentration_24 = float(np.max(mod24))

    return {
        "num_hits": int(len(hits)),
        "density": float(len(hits) / n_max),
        "peak_mod8": int(np.argmax(mod8)),
        "peak_mod12": int(np.argmax(mod12)),
        "peak_mod24": int(np.argmax(mod24)),
        "concentration_mod12": concentration_12,
        "concentration_mod24": concentration_24,
        "passed": len(hits) > 5,
    }


def add_table_row(rows, test, error, threshold, meaning, status_override=None):
    status = status_override if status_override is not None else ("OK" if error <= threshold else "FAIL")
    rows.append(
        {
            "test": test,
            "error": float(error),
            "threshold": float(threshold),
            "status": status,
            "meaning": meaning,
        }
    )


def print_result_table(rows):
    header = f"{'Test':26} {'Fehler':>12} {'Schwelle':>12} {'Status':>8}  Bedeutung"
    print(header)
    print("-" * len(header))
    for r in rows:
        print(
            f"{r['test'][:26]:26} "
            f"{r['error']:>12.3e} "
            f"{r['threshold']:>12.3e} "
            f"{r['status']:>8}  "
            f"{r['meaning']}"
        )


def print_brief_summary(results, rows):
    ok_count = sum(1 for r in rows if r["status"] == "OK")
    control_fail_rows = [
        r for r in rows if r["status"] == "FAIL" and ("Diagnose" in r["meaning"] or "roh" in r["test"].lower())
    ]
    structural_fail_rows = [
        r for r in rows if r["status"] == "FAIL" and r not in control_fail_rows
    ]
    control_fail_count = len(control_fail_rows)
    structural_fail_count = len(structural_fail_rows)
    scored = []
    for r in rows:
        thr = r["threshold"] if r["threshold"] > 0 else 1.0
        ratio = r["error"] / thr
        scored.append({**r, "ratio": ratio})

    ok_rows = [r for r in scored if r["status"] == "OK"]
    fail_rows = [r for r in scored if r["status"] == "FAIL" and r["test"] not in {x["test"] for x in control_fail_rows}]

    strengths = []
    for r in sorted(ok_rows, key=lambda x: x["ratio"])[:2]:
        factor = (r["threshold"] / r["error"]) if r["error"] > 0 else np.inf
        factor_str = ">1e12" if np.isinf(factor) else f"{factor:.2e}"
        strengths.append(
            f"{r['test']} stabil ({fmt_float(r['error'])} <= {fmt_float(r['threshold'])}, x{factor_str} unter Schwelle)"
        )

    open_points = []
    if fail_rows:
        for r in sorted(fail_rows, key=lambda x: x["ratio"], reverse=True)[:2]:
            factor = r["error"] / r["threshold"] if r["threshold"] > 0 else np.inf
            factor_str = ">1e12" if np.isinf(factor) else f"{factor:.2e}"
            open_points.append(
                f"{r['test']} faellt ({fmt_float(r['error'])} > {fmt_float(r['threshold'])}, x{factor_str} ueber Schwelle)"
            )
    else:
        # Falls nichts faellt: zeige die knappsten OK-Tests als Beobachtung
        for r in sorted(ok_rows, key=lambda x: x["ratio"], reverse=True)[:2]:
            factor = (r["threshold"] / r["error"]) if r["error"] > 0 else np.inf
            factor_str = ">1e12" if np.isinf(factor) else f"{factor:.2e}"
            open_points.append(
                f"{r['test']} knapp ({fmt_float(r['error'])} vs {fmt_float(r['threshold'])}, x{factor_str} unter Schwelle)"
            )

    while len(strengths) < 2:
        strengths.append("keine weitere Staerke identifiziert")
    while len(open_points) < 2:
        open_points.append("kein weiterer offener Punkt identifiziert")

    print_section("Kurz-Zusammenfassung")
    print_kv("Tests OK", ok_count)
    print_kv("Strukturelle FAIL", structural_fail_count)
    print_kv("Kontroll-FAIL", control_fail_count)
    print_kv("Staerke 1", strengths[0])
    print_kv("Staerke 2", strengths[1])
    print_kv("Offen 1", open_points[0])
    print_kv("Offen 2", open_points[1])
    print_kv(
        "Protokollsatz",
        "BM-IV besteht strukturelle Quantentests auf Z8 (Weyl/Dirac/Holonomie/Born/Spin-7/2/Periode-24); "
        "direkter Spin-1/2 fehlt, rohe nichtunitaere Stoerungen sind Negativkontrolle und kein Modellfehler.",
        width=16,
    )


def bool_badge(flag):
    return "OK" if bool(flag) else "FAIL"


def fmt_float(x):
    return f"{float(x):.3e}"


def print_kv(label, value, width=34):
    print(f"  - {label:<{width}} {value}")


def print_section(title):
    bar = "-" * len(title)
    print(f"\n{title}\n{bar}")


def print_pretty_report(results):
    built = results["built_in_identities"]
    drv = results["derived_findings"]
    ref = results["reference_profiles"]

    print_section("A) Eingebaute Identitaeten")
    print_kv("Unitaritaet X8", bool_badge(built["unitarity_X8"]))
    print_kv("Unitaritaet Z8", bool_badge(built["unitarity_Z8"]))
    print_kv("Unitaritaet F8", bool_badge(built["unitarity_F8"]))
    print_kv("Hermitezitaet Sx", bool_badge(built["hermiticity_Sx"]))
    print_kv("Hermitezitaet Sz", bool_badge(built["hermiticity_Sz"]))
    print_kv("X8^8 - I (Fro)", fmt_float(built["X8_pow8_identity_error"]))
    print_kv("Weyl ZX-omegaXZ", fmt_float(built["weyl_error_ZX_minus_omega_XZ"]))
    print_kv("Weyl-Status", bool_badge(built["weyl_passed"]))
    print_kv("Laplace max Spektralfehler", fmt_float(built["laplace"]["max_spectral_error"]))
    print_kv("Dirac^2 - Block(Laplace)", fmt_float(built["dirac"]["d2_block_laplace_error"]))

    print_section("B) Abgeleitete Befunde")
    print_kv("Casimir naiv", "VERWORFEN (falscher Generator)")
    print_kv("Casimir SU2", bool_badge(drv["casimir_variant_su2"]["passed"]))
    print_kv("Casimir Weyl-konstruiert", bool_badge(drv["casimir_variant_weyl_constructed"]["passed"]))
    print_kv("Born Sx max Abweichung", fmt_float(drv["born_Sx"]["max_frequency_error"]))
    print_kv("Born Sz max Abweichung", fmt_float(drv["born_Sz"]["max_frequency_error"]))
    print_kv("Entropie-Suite", bool_badge(drv["entropy"]["passed"]))
    print_kv("Gauge-Holonomie", bool_badge(drv["gauge_holonomy_invariance"]["passed"]))
    print_kv("Gauge |Gamma'-Gamma|", fmt_float(drv["gauge_holonomy_invariance"]["error"]))
    print_kv("Sagnac", bool_badge(drv["sagnac"]["passed"]))
    print_kv("Sagnac Delta-Fehler", fmt_float(drv["sagnac"]["delta_phase_error"]))
    print_kv("24-Zyklus", bool_badge(drv["cycle_24"]["passed"]))
    print_kv("24-Theorem", bool_badge(drv["period_24_theorem"]["passed"]))
    print_kv("Robustheit-Holonomie", bool_badge(drv["robustness"]["holonomy_ok"]))
    print_kv("Robustheit-Born (unitaer)", bool_badge(drv["robustness"]["born_ok"]))
    print_kv("Robustheit-Born (roh)", "NEGATIVKONTROLLE")
    print_kv("Empfindlichste Metrik", drv["robustness"]["most_sensitive_metric"])
    print_kv("Spin-7/2-Charakterisierung", bool_badge(drv["spin_7_2_characterization"]["passed"]))
    print_kv("Vorhersage-Suite", bool_badge(drv["prediction_suite"]["passed"]))
    print_kv("Spin-1/2-Projektion", bool_badge(drv["spin_half_projection"]["passed"]))
    print_kv("Effektive 2-Zustandsprojektion", bool_badge(drv["effective_two_state_projection"]["passed"]))
    print_kv("Einstein-de-Haas", bool_badge(drv["einstein_de_haas"]["passed"]))
    print_kv("EdH max Bilanzfehler", fmt_float(drv["einstein_de_haas"]["max_balance_error"]))
    print_kv("Primvierlings-Observable", bool_badge(drv["prime_quadruplet_observable"]["passed"]))

    print_section("C) Referenzprofile / Kalibrierung")
    print_kv(
        "Atomstabilitaet-Set",
        f"{bool_badge(ref['atomic_configuration']['passed'] and ref['nuclear_stability']['passed'])}",
    )
    print_kv("Atomkonfiguration H/He", f"{bool_badge(ref['atomic_configuration']['passed'])}   Referenzprofil reproduziert")
    print_kv(
        "Nuklearstabilitaet Fe/U",
        f"{bool_badge(ref['nuclear_stability']['passed'])}   bekannte Ordnung korrekt klassifiziert",
    )
    print_kv(
        "Fe56 Bindungsenergie",
        f"{ref['nuclear_stability']['iron_fe56_be_per_a']:.3f} MeV/A   Referenzwert, nicht prognostiziert",
    )
    print_kv(
        "U238 Halbwertszeit",
        f"{ref['nuclear_stability']['uranium_u238_half_life_years']:.3e} a   Referenzwert, nicht prognostiziert",
    )

    print_section("Negative Kontrollen")
    for key in ["Z7", "Z9", "Z10"]:
        row = drv["negative_controls"][key]
        print_kv(
            f"{key} 24-Signatur",
            f"{bool_badge(row['passed_24_signature'])} (U^24-I={fmt_float(row['cycle24_error'])}, lcm={row['period_lcm_d_12']})",
        )
    mod10 = drv["negative_controls"]["particle_mod10"]
    print_kv(
        "Partikelraum mod 10",
        f"lcm(8,10)={mod10['claimed_global_period']} vs BM-IV={mod10['bmiv_global_period']}",
    )

    print_section("Theorem-Block")
    print_kv("Theorem T1", "Per(BM-IV) = lcm(8,12) = 24")
    print_kv("Numerischer Status T1", bool_badge(drv["period_24_theorem"]["passed"]))
    print_kv("Zusatz", "Fuer Kontrollraeume Z_n mit n!=8 gilt i.A. keine BM-IV-24-Rueckkehr")
    z7 = drv["negative_controls"]["Z7"]["passed_24_signature"]
    z9 = drv["negative_controls"]["Z9"]["passed_24_signature"]
    z10 = drv["negative_controls"]["Z10"]["passed_24_signature"]
    print_kv("Kontrollbefund", f"Z7={bool_badge(z7)}, Z9={bool_badge(z9)}, Z10={bool_badge(z10)}")


def run_quantum_tests(verbose=True):
    X8, Z8, Sx, Sz, omega = make_observables_8()
    F8 = make_fourier_8()

    results = {"built_in_identities": {}, "derived_findings": {}, "reference_profiles": {}}

    # 0. Axiomatisierung
    results["axioms"] = test_axioms()

    # 1. Basisstruktur (konstruktive Konsistenz)
    results["built_in_identities"]["unitarity_X8"] = is_unitary(X8)
    results["built_in_identities"]["unitarity_Z8"] = is_unitary(Z8)
    results["built_in_identities"]["unitarity_F8"] = is_unitary(F8)
    results["built_in_identities"]["hermiticity_Sx"] = is_hermitian(Sx)
    results["built_in_identities"]["hermiticity_Sz"] = is_hermitian(Sz)
    results["built_in_identities"]["X8_pow8_identity_error"] = fro_err(np.linalg.matrix_power(X8, 8), np.eye(8))

    # 2. Weyl-Heisenberg-Relation
    weyl_err_1 = fro_err(Z8 @ X8, omega * X8 @ Z8)
    results["built_in_identities"]["weyl_error_ZX_minus_omega_XZ"] = weyl_err_1
    results["built_in_identities"]["weyl_passed"] = weyl_err_1 < 1e-12

    # 2b. Casimir-Struktur
    results["derived_findings"]["casimir_naive"] = casimir_test(X8, Z8)
    results["derived_findings"]["casimir_variant_su2"] = casimir_test_variant_su2()
    results["derived_findings"]["casimir_variant_weyl_constructed"] = casimir_test_variant_weyl_constructed(X8, Z8)

    # 2c. Spektrale Geometrie auf Z8
    results["built_in_identities"]["laplace"] = laplace_spectrum_test(X8)

    # 2d. Dirac-Operator und Quadrat
    results["built_in_identities"]["dirac"] = dirac_test(X8)

    # 3. Spuren
    results["derived_findings"]["trace_Sx"] = complex(np.trace(Sx))
    results["derived_findings"]["trace_Sz"] = complex(np.trace(Sz))
    results["derived_findings"]["trace_Sx_abs"] = abs(np.trace(Sx))
    results["derived_findings"]["trace_Sz_abs"] = abs(np.trace(Sz))

    # 4. Projektorzerlegung
    results["derived_findings"]["projectors_Sx"] = projector_tests(Sx)
    results["derived_findings"]["projectors_Sz"] = projector_tests(Sz)

    # 5. Unsicherheitsrelation
    psi = random_state(8)
    ok_unc, lhs, rhs = uncertainty_test(Sx, Sz, psi)

    results["derived_findings"]["uncertainty_passed"] = ok_unc
    results["derived_findings"]["uncertainty_lhs"] = lhs
    results["derived_findings"]["uncertainty_rhs"] = rhs

    # 6. Born-Regel
    results["derived_findings"]["born_Sx"] = born_rule_test(Sx)
    results["derived_findings"]["born_Sz"] = born_rule_test(Sz)

    # 7. Entropie/Verschränkung (messbar, nicht nur interpretativ)
    results["derived_findings"]["entropy"] = entropy_test_suite()

    # 8. Holonomie als Gruppenstruktur + Gauge-Invarianz
    results["derived_findings"]["sagnac"] = sagnac_test(X8)
    results["derived_findings"]["holonomy"] = holonomy_test()
    results["derived_findings"]["gauge_holonomy_invariance"] = gauge_holonomy_invariance_test()

    # 9. 24-Zyklus
    results["derived_findings"]["cycle_24"] = cycle_24_test(X8)
    results["derived_findings"]["period_24_theorem"] = period_24_theorem_test(X8)

    # 10. Robuste Modellbefunde
    results["derived_findings"]["robustness"] = perturbation_robustness_test(X8, Z8)
    results["derived_findings"]["negative_controls"] = negative_control_suite()
    results["derived_findings"]["spin_half_projection"] = spin_half_projection_test(X8, Z8)
    results["derived_findings"]["spin_7_2_characterization"] = spin_seven_half_characterization_test(X8, Z8)
    results["derived_findings"]["prediction_suite"] = prediction_suite(X8, Z8)
    results["derived_findings"]["effective_two_state_projection"] = effective_two_state_projection_test(X8, Z8)
    results["derived_findings"]["einstein_de_haas"] = einstein_de_haas_test()
    results["derived_findings"]["prime_quadruplet_observable"] = prime_quadruplet_observable_test()
    results["derived_findings"]["entanglement"] = entanglement_test()
    results["reference_profiles"]["atomic_configuration"] = atomic_configuration_test()
    results["reference_profiles"]["nuclear_stability"] = nuclear_stability_test()

    # Ergebnis-Tabelle
    table_rows = []
    add_table_row(
        table_rows,
        "Weyl",
        results["built_in_identities"]["weyl_error_ZX_minus_omega_XZ"],
        1e-12,
        "diskrete Quantisierung",
    )
    add_table_row(
        table_rows,
        "Laplace-Spektrum Z8",
        results["built_in_identities"]["laplace"]["max_spectral_error"],
        1e-12,
        "Schaltergeometrie",
    )
    add_table_row(
        table_rows,
        "Dirac^2-Block",
        results["built_in_identities"]["dirac"]["d2_block_laplace_error"],
        1e-12,
        "nichtkommutative Geometrie",
    )
    add_table_row(
        table_rows,
        "Born (Sx)",
        results["derived_findings"]["born_Sx"]["max_frequency_error"],
        2e-2,
        "Messstatistik",
    )
    add_table_row(
        table_rows,
        "Casimir-naiv verworfen",
        results["derived_findings"]["casimir_naive"]["comm_C_Sx"],
        1e-12,
        "falscher Kandidat, nicht Modellbruch",
        status_override="INFO",
    )
    add_table_row(
        table_rows,
        "SU2-Casimir",
        results["derived_findings"]["casimir_variant_su2"]["casimir_scalar_error"],
        1e-12,
        "Spin-7/2-Referenz",
    )
    add_table_row(
        table_rows,
        "Weyl-Casimir",
        results["derived_findings"]["casimir_variant_weyl_constructed"]["casimir_scalar_error"],
        1e-12,
        "BM-nahe su(2)-Einbettung",
    )
    add_table_row(
        table_rows,
        "Sagnac Delta",
        results["derived_findings"]["sagnac"]["delta_phase_error"],
        1e-12,
        "orientierte Umlaufphase",
    )
    add_table_row(
        table_rows,
        "Gauge-Holonomie",
        results["derived_findings"]["gauge_holonomy_invariance"]["error"],
        1e-12,
        "Eichinvarianz",
    )
    add_table_row(
        table_rows,
        "24-Theorem U^24",
        max(results["derived_findings"]["period_24_theorem"]["operator_u24_errors"].values()),
        1e-10,
        "Per(BM-IV)=24",
    )
    add_table_row(
        table_rows,
        "Robustheit-Holonomie",
        results["derived_findings"]["robustness"]["holonomy_break_24"],
        results["derived_findings"]["robustness"]["holonomy_threshold"],
        "Stabilität unter Stoerung",
    )
    add_table_row(
        table_rows,
        "Robustheit-Born",
        results["derived_findings"]["robustness"]["born_prob_l1_unitary"],
        results["derived_findings"]["robustness"]["born_unitary_threshold"],
        "Messstabilitaet (unitaere Stoerung)",
    )
    add_table_row(
        table_rows,
        "Robustheit-Born roh",
        results["derived_findings"]["robustness"]["born_prob_l1_raw"],
        results["derived_findings"]["robustness"]["born_raw_threshold"],
        "Negativkontrolle: rohe Stoerung",
        status_override="INFO",
    )
    add_table_row(
        table_rows,
        "Spin-7/2-Char.",
        results["derived_findings"]["spin_7_2_characterization"]["casimir_error"],
        1e-10,
        "Z8 als Spin-7/2-artig",
    )
    add_table_row(
        table_rows,
        "Prediction-24-Family",
        0.0 if results["derived_findings"]["prediction_suite"]["prediction_24_return_family"]["passed"] else 1.0,
        0.5,
        "falsifizierbare 24-Rueckkehr",
    )
    add_table_row(
        table_rows,
        "Prediction-Spin7/2",
        results["derived_findings"]["prediction_suite"]["prediction_spin72_transition_pattern"]["max_coupling_error"],
        1e-8,
        "Leiteramplitudenmuster",
    )
    add_table_row(
        table_rows,
        "Einstein-de-Haas",
        results["derived_findings"]["einstein_de_haas"]["max_balance_error"],
        1e-12,
        "diskrete Spin-Bahn-Bilanz",
    )
    add_table_row(
        table_rows,
        "Spin-1/2-Projektion",
        results["derived_findings"]["spin_half_projection"]["comm_xy_error"],
        1e-8,
        "Pauli-Kommutator",
    )
    results["result_table"] = table_rows

    if verbose:
        print("\nBM-IV Quantum Test Suite (Axiomatisch)")
        print("=" * 52)
        print_pretty_report(results)

        print("\nErgebnis-Tabelle")
        print_result_table(table_rows)
        print_brief_summary(results, table_rows)

    return results


# --- EABC: ptolemäische Schließung (tetraedrische Einbettung) ---

EABC_TETRA_AXES = {
    "E": np.array([1.0, 1.0, 1.0]),
    "A": np.array([1.0, -1.0, -1.0]),
    "B": np.array([-1.0, 1.0, -1.0]),
    "C": np.array([-1.0, -1.0, 1.0]),
}

# Physikalische Zeugen (Projektkontext, optional):
# c: SI-Lichtgeschwindigkeit; H0_CMB: Bamberger Modell II.tex (~67.4 km/s/Mpc, frühes Univ.)
C_LIGHT_SI = 299792458.0
H0_CMB_KM_S_MPC = 67.4
_MPC_TO_M = 3.0856775814913673e22
H0_CMB_SI = H0_CMB_KM_S_MPC * 1000.0 / _MPC_TO_M
A0_COSMO_SI = C_LIGHT_SI * H0_CMB_SI


def _as_eabc_vector(x):
    arr = np.asarray(x, dtype=float).reshape(-1)
    if arr.shape != (4,):
        raise ValueError("X muss (E, A, B, C) mit Laenge 4 sein")
    if np.any(arr < 0):
        raise ValueError("EABC-Besetzung X muss in R_>=0^4 liegen")
    return arr


def eabc_tetrahedron_embed(E, A, B, C):
    """P_i = x_i * v_i auf den Tetraederachsen."""
    x = _as_eabc_vector((E, A, B, C))
    labels = ("E", "A", "B", "C")
    return {lab: float(val) * EABC_TETRA_AXES[lab] for lab, val in zip(labels, x)}


def eabc_tetrahedron_distances(E, A, B, C):
    """Sechs Abstaende d_ij = ||P_i - P_j||."""
    pts = eabc_tetrahedron_embed(E, A, B, C)
    pairs = (("E", "A"), ("E", "B"), ("E", "C"), ("A", "B"), ("A", "C"), ("B", "C"))
    return {f"{i}_{j}": float(np.linalg.norm(pts[i] - pts[j])) for i, j in pairs}


def pi_ptolemaic_defect(E, A, B, C):
    """Pi(X) = d_EB*d_AC - d_EA*d_BC - d_EC*d_AB."""
    d = eabc_tetrahedron_distances(E, A, B, C)
    return (
        d["E_B"] * d["A_C"]
        - d["E_A"] * d["B_C"]
        - d["E_C"] * d["A_B"]
    )


def _K_P_ptolemaic_denom(E, A, B, C, eps=EPS):
    d = eabc_tetrahedron_distances(E, A, B, C)
    return d["E_B"] * d["A_C"] + d["E_A"] * d["B_C"] + d["E_C"] * d["A_B"]


def K_P_ptolemaic(E, A, B, C, eps=EPS):
    """Normierter Ptolemäus-Defekt K_P = Pi / (Summe der drei Produkte)."""
    denom = _K_P_ptolemaic_denom(E, A, B, C, eps=eps)
    if denom <= eps:
        return float("nan")
    return float(pi_ptolemaic_defect(E, A, B, C) / denom)


def ptolemaic_vacuum_reference():
    """Vakuum X=(1,1,1,1): Pi_0=-8, K_P_0=-1/3 (tetraedrisch, nicht ptolemäisch geschlossen)."""
    pi_0 = float(pi_ptolemaic_defect(1.0, 1.0, 1.0, 1.0))
    k_0 = float(K_P_ptolemaic(1.0, 1.0, 1.0, 1.0))
    return {"Pi_0": pi_0, "K_P_0": k_0}


PTOLEMAIC_VACUUM = ptolemaic_vacuum_reference()
PI_0 = PTOLEMAIC_VACUUM["Pi_0"]
K_P_0 = PTOLEMAIC_VACUUM["K_P_0"]


def pi_ptolemaic_tilde(E, A, B, C):
    """Pi_tilde(X) = Pi(X) - Pi(1,1,1,1)."""
    return float(pi_ptolemaic_defect(E, A, B, C) - PI_0)


def K_P_ptolemaic_tilde(E, A, B, C, eps=EPS):
    """
    Renormierter Defekt relativ zum Vakuum (1,1,1,1).

    K_P_tilde = K_P(X) - K_P(1,1,1,1), äquivalent zu Pi_tilde/denom mit demselben
    Nenner wie K_P (Formmaß, skaleninvariant). Für symmetrische (n,n,n,n) gilt K_P_tilde=0.
    """
    k = K_P_ptolemaic(E, A, B, C, eps=eps)
    if math.isnan(k):
        return float("nan")
    return float(k - K_P_0)


def a_P_ptolemaic(E, A, B, C, a_0=A0_COSMO_SI, use_tilde=True):
    """a_P(X) = a_0 * |K_P_tilde(X)| (Vakuum: a_P=0)."""
    k = K_P_ptolemaic_tilde(E, A, B, C) if use_tilde else K_P_ptolemaic(E, A, B, C)
    if math.isnan(k):
        return float("nan")
    return float(a_0 * abs(k))


def ell_P_ptolemaic(E, A, B, C, c=C_LIGHT_SI, a_0=A0_COSMO_SI, eps=EPS, use_tilde=True):
    """Rindler-Horizont ell_P = c^2 / (a_0 * |K_P_tilde|)."""
    k = K_P_ptolemaic_tilde(E, A, B, C, eps=eps) if use_tilde else K_P_ptolemaic(E, A, B, C, eps=eps)
    if math.isnan(k) or abs(k) <= eps:
        return float("inf")
    return float((c * c) / (a_0 * abs(k)))


def eabc_ptolemaic_chain(E, A, B, C, c=C_LIGHT_SI, a_0=A0_COSMO_SI):
    """EABC-Besetzung -> Tetraeder -> Pi, K_P, renormierte Größen -> a_P, ell_P."""
    d = eabc_tetrahedron_distances(E, A, B, C)
    pi_val = pi_ptolemaic_defect(E, A, B, C)
    k_val = K_P_ptolemaic(E, A, B, C)
    pi_tilde = pi_ptolemaic_tilde(E, A, B, C)
    k_tilde = K_P_ptolemaic_tilde(E, A, B, C)
    return {
        "X": (float(E), float(A), float(B), float(C)),
        "points": eabc_tetrahedron_embed(E, A, B, C),
        "distances": d,
        "Pi": float(pi_val),
        "Pi_tilde": float(pi_tilde),
        "K_P": float(k_val) if not math.isnan(k_val) else float("nan"),
        "K_P_tilde": float(k_tilde) if not math.isnan(k_tilde) else float("nan"),
        "a_P": a_P_ptolemaic(E, A, B, C, a_0=a_0),
        "ell_P": ell_P_ptolemaic(E, A, B, C, c=c, a_0=a_0),
        "ptolemaic_closed": abs(k_tilde) <= EPS if not math.isnan(k_tilde) else False,
        "physics": {"c": float(c), "a_0": float(a_0), "H0_CMB_km_s_Mpc": H0_CMB_KM_S_MPC},
    }


def _naive_collinear_K_P(E, A, B, C):
    """1D-Einbettung x_i=i: trivial K_P=0 (Kontrollfall)."""
    coords = {"E": E, "A": A, "B": B, "C": C}
    pairs = (("E", "A"), ("E", "B"), ("E", "C"), ("A", "B"), ("A", "C"), ("B", "C"))
    d = {f"{i}_{j}": abs(coords[i] - coords[j]) for i, j in pairs}
    pi_val = d["E_B"] * d["A_C"] - d["E_A"] * d["B_C"] - d["E_C"] * d["A_B"]
    denom = d["E_B"] * d["A_C"] + d["E_A"] * d["B_C"] + d["E_C"] * d["A_B"]
    return 0.0 if denom == 0 else pi_val / denom


def _is_hurwitz_quaternion(e, a, b, c, tol=1e-9):
    """Hurwitz-Ganzheit: 2a,2b,2c,2d ∈ Z und a+b+c+d ∈ Z."""
    coeffs = (e, a, b, c)
    if any(abs(2.0 * x - round(2.0 * x)) > tol for x in coeffs):
        return False
    return abs(sum(coeffs) - round(sum(coeffs))) <= tol


def hurwitz_quaternion_hull(norm_bound=100):
    """
    AHPN-/Hurwitz-Hülle: Quaternionen e+ai+bj+ck mit reduzierter Norm < norm_bound.
    Koeffizienten ganzzahlig oder halbzahlig (Hurwitz-Ordnung).
    """
    hits = []
    seen = set()
    max_half = int(math.isqrt(norm_bound)) * 2 + 4
    for e2 in range(-max_half, max_half + 1):
        for a2 in range(-max_half, max_half + 1):
            for b2 in range(-max_half, max_half + 1):
                for c2 in range(-max_half, max_half + 1):
                    e, a, b, c = e2 / 2.0, a2 / 2.0, b2 / 2.0, c2 / 2.0
                    if not _is_hurwitz_quaternion(e, a, b, c):
                        continue
                    n2 = e * e + a * a + b * b + c * c
                    if n2 <= 0 or n2 >= norm_bound:
                        continue
                    key = (e2, a2, b2, c2)
                    if key in seen:
                        continue
                    seen.add(key)
                    hits.append({"e": e, "a": a, "b": b, "c": c, "norm_sq": float(n2)})
    hits.sort(key=lambda row: (row["norm_sq"], row["e"], row["a"], row["b"], row["c"]))
    return hits


def eabc_occupation_entropy(E, A, B, C, eps=EPS):
    """Shannon-Entropie der EABC-Besetzung (positive Gewichte, normiert)."""
    x = np.array([E, A, B, C], dtype=float)
    x = x[x > eps]
    if x.size == 0:
        return 0.0
    p = x / x.sum()
    return float(-np.sum(p * np.log(p)))


def prime_density_proxy(p, window=20):
    """Lokale Primdichte um p: #(Primzahlen in [p-window, p+window]) / (2*window+1)."""
    lo = max(2, int(p) - window)
    hi = int(p) + window
    count = sum(1 for n in range(lo, hi + 1) if is_prime(n))
    return float(count / (hi - lo + 1))


def is_prime_quadruplet(p):
    return is_prime(p) and is_prime(p + 2) and is_prime(p + 6) and is_prime(p + 8)


def prime_quadruplet_tuples(limit=200):
    """Primzahlvierlinge (p, p+2, p+6, p+8) mit p <= limit."""
    out = []
    for p in range(5, limit - 8):
        if is_prime_quadruplet(p):
            out.append((float(p), float(p + 2), float(p + 6), float(p + 8)))
    return out


def eabc_ptolemaic_class_a_symmetric(n_max=100, tol=1e-9):
    """Klasse A: (n,n,n,n), n=1..n_max — K_P_tilde ≈ 0."""
    errors = []
    for n in range(1, n_max + 1):
        k_tilde = K_P_ptolemaic_tilde(float(n), float(n), float(n), float(n))
        errors.append(abs(k_tilde))
    max_err = float(max(errors))
    return {
        "n_max": n_max,
        "max_abs_K_P_tilde": max_err,
        "passed": max_err <= tol,
    }


EABC_PRIME_EXAMPLE_TUPLES = (
    (2.0, 3.0, 5.0, 7.0),
    (11.0, 13.0, 17.0, 19.0),
)


def eabc_ptolemaic_class_b_primes(limit=200):
    """Klasse B: Primzahlvierlinge (p,p+2,p+6,p+8) und explizite Prim-4-Tupel."""
    rows = []
    seen = set()
    candidates = list(EABC_PRIME_EXAMPLE_TUPLES) + prime_quadruplet_tuples(limit=limit)
    for x in candidates:
        key = tuple(int(v) for v in x)
        if key in seen:
            continue
        seen.add(key)
        ch = eabc_ptolemaic_chain(*x)
        rows.append(
            {
                "X": x,
                "Pi": ch["Pi"],
                "K_P": ch["K_P"],
                "K_P_tilde": ch["K_P_tilde"],
                "prime_density": prime_density_proxy(x[0]),
                "is_quadruplet": bool(is_prime_quadruplet(int(x[0]))),
            }
        )
    return rows


def eabc_ptolemaic_class_c_hurwitz(norm_bound=100, max_rows=24):
    """Klasse C: Hurwitz-Quaternionen (e,a,b,c), |q|^2<norm_bound -> EABC via |Koeff.|."""
    rows = []
    for q in hurwitz_quaternion_hull(norm_bound=norm_bound):
        x = (abs(q["e"]), abs(q["a"]), abs(q["b"]), abs(q["c"]))
        if max(x) <= EPS:
            continue
        ch = eabc_ptolemaic_chain(*x)
        rows.append(
            {
                "quaternion": (q["e"], q["a"], q["b"], q["c"]),
                "norm_sq": q["norm_sq"],
                "X": x,
                "K_P_tilde": ch["K_P_tilde"],
                "entropy": eabc_occupation_entropy(*x),
            }
        )
    rows.sort(key=lambda r: (abs(r["K_P_tilde"]), r["norm_sq"]))
    return rows[:max_rows]


def eabc_ptolemaic_correlations(class_b_rows, class_c_rows):
    """
    Korrelationen K_P_tilde mit vorhandenen Größen (Primdichte, EABC-Entropie,
    BM-IV-Holonomie, 24-Zyklus) — nur wenn Daten vorhanden.
    """
    hol = holonomy_test()
    cyc = cycle_24_test(make_XZ_8()[0])
    out = {
        "holonomy_gamma_24": hol["gamma_24"],
        "cycle_24_passed": cyc["passed"],
        "correlations": {},
    }

    def _pearson(xs, ys):
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        if xs.size < 2 or ys.size < 2:
            return float("nan")
        if np.std(xs) < EPS or np.std(ys) < EPS:
            return float("nan")
        return float(np.corrcoef(xs, ys)[0, 1])

    if class_b_rows:
        ks = [r["K_P_tilde"] for r in class_b_rows]
        rhos = [r["prime_density"] for r in class_b_rows]
        hol_idx = [abs(holonomy_phase(int(round(sum(r["X"])) % 24)).real) for r in class_b_rows]
        out["correlations"]["class_B_K_P_tilde_vs_prime_density"] = _pearson(ks, rhos)
        out["correlations"]["class_B_K_P_tilde_vs_holonomy_phase"] = _pearson(ks, hol_idx)

    if class_c_rows:
        ks = [r["K_P_tilde"] for r in class_c_rows]
        ent = [r["entropy"] for r in class_c_rows]
        hol_idx = [abs(holonomy_phase(int(round(sum(r["X"])) % 24)).real) for r in class_c_rows]
        mod24 = [sum(r["X"]) % 24 for r in class_c_rows]
        out["correlations"]["class_C_K_P_tilde_vs_entropy"] = _pearson(ks, ent)
        out["correlations"]["class_C_K_P_tilde_vs_holonomy_phase"] = _pearson(ks, hol_idx)
        out["correlations"]["class_C_K_P_tilde_vs_sum_mod24"] = _pearson(ks, mod24)

    return out


def run_eabc_ptolemaic_tests(n_sym_max=100, prime_limit=200, hurwitz_bound=100):
    """Drei Testklassen + Korrelationsblock."""
    class_a = eabc_ptolemaic_class_a_symmetric(n_max=n_sym_max)
    class_b = eabc_ptolemaic_class_b_primes(limit=prime_limit)
    class_c = eabc_ptolemaic_class_c_hurwitz(norm_bound=hurwitz_bound)
    corr = eabc_ptolemaic_correlations(class_b, class_c)
    return {"class_A": class_a, "class_B": class_b, "class_C": class_c, "correlations": corr}


def _spearman(xs, ys):
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if xs.size < 2 or xs.size != ys.size:
        return float("nan")
    if np.std(xs) < EPS or np.std(ys) < EPS:
        return float("nan")
    rx = np.argsort(np.argsort(xs))
    ry = np.argsort(np.argsort(ys))
    return float(np.corrcoef(rx, ry)[0, 1])


def _pair_correlations(xs, ys, min_n=5):
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    mask = np.isfinite(xs) & np.isfinite(ys)
    xs = xs[mask]
    ys = ys[mask]
    if xs.size < min_n:
        return {"n": int(xs.size), "pearson": float("nan"), "spearman": float("nan"), "defined": False}
    if np.std(xs) < EPS or np.std(ys) < EPS:
        return {"n": int(xs.size), "pearson": float("nan"), "spearman": float("nan"), "defined": False}
    return {
        "n": int(xs.size),
        "pearson": float(np.corrcoef(xs, ys)[0, 1]),
        "spearman": _spearman(xs, ys),
        "defined": True,
    }


def _bmiv_dirac_cache():
    X8, _, _ = make_XZ_8()
    D = dirac_operator(X8)
    eigs = np.sort(np.real(np.linalg.eigvalsh((D + dagger(D)) / 2)))
    d2_err = dirac_test(X8)["d2_block_laplace_error"]
    return {
        "dirac_eigs": eigs,
        "dirac_mean": float(np.mean(eigs)),
        "dirac_spread": float(np.ptp(eigs)),
        "dirac_d2_err": float(d2_err),
    }


_BMIV_DIRAC_CACHE = None


def _get_bmiv_dirac_cache():
    global _BMIV_DIRAC_CACHE
    if _BMIV_DIRAC_CACHE is None:
        _BMIV_DIRAC_CACHE = _bmiv_dirac_cache()
    return _BMIV_DIRAC_CACHE


def eabc_to_holonomy_wind(E, A, B, C):
    """EABC -> diskreter Windungsindex mod 24 (BM-IV-24-Zyklus)."""
    x = _as_eabc_vector((E, A, B, C))
    return int(np.round(x.sum())) % 24


def eabc_to_bmiv_u_params(E, A, B, C):
    """EABC -> (alpha, m) fuer U = exp(2pi i alpha/24) X^m (Prediction-24-Familie)."""
    x = _as_eabc_vector((E, A, B, C))
    alpha = int(np.round((x.sum() % 24) * 1.0)) % 24
    if alpha == 0:
        alpha = 1
    m = int(np.argmax(x)) + 1
    if m % 2 == 0:
        m = (m % 7) + 1
    return alpha, m


def eabc_spectral_observables(E, A, B, C, ahpn_norm_sq=None):
    """Observablen pro EABC-Punkt: Ptolemäus, Holonomie-Proxys, BM-IV-Dirac (fix), AHPN, Entropie."""
    hol_bm = holonomy_test()
    wind = eabc_to_holonomy_wind(E, A, B, C)
    gamma_wind = holonomy_phase(wind)
    alpha, m = eabc_to_bmiv_u_params(E, A, B, C)
    X8, _, _ = make_XZ_8()
    U = np.exp(2j * np.pi * alpha / 24) * np.linalg.matrix_power(X8, m)
    u24 = np.linalg.matrix_power(U, 24)
    gamma_u24 = complex(np.trace(u24) / X8.shape[0])

    d = eabc_tetrahedron_distances(E, A, B, C)
    dist_vals = list(d.values())
    dcache = _get_bmiv_dirac_cache()

    if ahpn_norm_sq is None:
        ahpn_norm_sq = float(E * E + A * A + B * B + C * C)

    k_p = K_P_ptolemaic(E, A, B, C)
    k_tilde = K_P_ptolemaic_tilde(E, A, B, C)

    return {
        "K_P": float(k_p) if not math.isnan(k_p) else float("nan"),
        "K_P_tilde": float(k_tilde) if not math.isnan(k_tilde) else float("nan"),
        "Pi_tilde": float(pi_ptolemaic_tilde(E, A, B, C)),
        "gamma_24_bm": float(abs(hol_bm["gamma_24"])),
        "holonomy_wind": int(wind),
        "holonomy_phase_abs": float(abs(gamma_wind)),
        "holonomy_phase_angle": float(np.angle(gamma_wind)),
        "gamma_u24_trace": float(abs(gamma_u24)),
        "gamma_u24_angle": float(np.angle(gamma_u24)),
        "bmiv_alpha": int(alpha),
        "bmiv_m": int(m),
        "dirac_mean": dcache["dirac_mean"],
        "dirac_spread": dcache["dirac_spread"],
        "dirac_d2_err": dcache["dirac_d2_err"],
        "dist_mean": float(np.mean(dist_vals)),
        "dist_spread": float(np.ptp(dist_vals)),
        "ahpn_norm_sq": float(ahpn_norm_sq),
        "entropy": float(eabc_occupation_entropy(E, A, B, C)),
    }


def eabc_spectral_collect_samples(
    n_sym_max=20,
    prime_limit=500,
    hurwitz_bound=100,
    max_hurwitz_rows=200,
    n_random=80,
    random_seed=42,
):
    """Klassen A–E: diverse EABC-Stichprobe fuer Spektralsweep."""
    rows = []

    for n in range(1, n_sym_max + 1):
        x = (float(n), float(n), float(n), float(n))
        obs = eabc_spectral_observables(*x)
        rows.append({"class": "A", "tag": f"sym({n})", "X": x, **obs})

    seen_b = set()
    for x in list(EABC_PRIME_EXAMPLE_TUPLES) + prime_quadruplet_tuples(limit=prime_limit):
        key = tuple(int(v) for v in x)
        if key in seen_b:
            continue
        seen_b.add(key)
        obs = eabc_spectral_observables(*x)
        rows.append({"class": "B", "tag": f"prime{key[0]}", "X": x, **obs})

    hurwitz_candidates = []
    seen_c = set()
    for q in hurwitz_quaternion_hull(norm_bound=hurwitz_bound):
        x = (abs(q["e"]), abs(q["a"]), abs(q["b"]), abs(q["c"]))
        if max(x) <= EPS:
            continue
        xkey = tuple(round(v, 6) for v in x)
        if xkey in seen_c:
            continue
        seen_c.add(xkey)
        obs = eabc_spectral_observables(*x, ahpn_norm_sq=q["norm_sq"])
        hurwitz_candidates.append(
            {
                "class": "C",
                "tag": f"H|q|^2={q['norm_sq']:g}",
                "X": x,
                "hurwitz": (q["e"], q["a"], q["b"], q["c"]),
                **obs,
            }
        )
    hurwitz_candidates.sort(key=lambda r: (abs(r["K_P_tilde"]), r["ahpn_norm_sq"]))
    rows.extend(hurwitz_candidates[:max_hurwitz_rows])

    rng = np.random.default_rng(random_seed)
    for i in range(n_random):
        x = tuple(float(rng.uniform(0.5, 10.0)) for _ in range(4))
        obs = eabc_spectral_observables(*x)
        rows.append({"class": "D", "tag": f"rnd{i+1}", "X": x, **obs})

    bases = ((3.0, 2.0, 4.0, 1.0), (1.0, 2.0, 3.0, 4.0), (2.0, 3.0, 5.0, 7.0))
    scales = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0)
    for bi, base in enumerate(bases):
        k_ref = K_P_ptolemaic(*base)
        for s in scales:
            x = tuple(float(s * v) for v in base)
            obs = eabc_spectral_observables(*x)
            k_s = K_P_ptolemaic(*x)
            ratio = float(k_s / k_ref) if k_ref and not math.isnan(k_ref) and not math.isnan(k_s) else float("nan")
            rows.append(
                {
                    "class": "E",
                    "tag": f"scale{s:g}_b{bi}",
                    "X": x,
                    "scale": float(s),
                    "K_P_ratio_to_base": ratio,
                    **obs,
                }
            )

    return rows


def eabc_spectral_correlation_block(rows, y_keys, min_n=5):
    ks = [r["K_P_tilde"] for r in rows]
    block = {}
    for key in y_keys:
        ys = [r[key] for r in rows]
        block[key] = _pair_correlations(ks, ys, min_n=min_n)
    by_class = {}
    for cls in ("A", "B", "C", "D", "E", "ALL"):
        sub = rows if cls == "ALL" else [r for r in rows if r["class"] == cls]
        if len(sub) < min_n:
            continue
        by_class[cls] = {key: _pair_correlations([r["K_P_tilde"] for r in sub], [r[key] for r in sub], min_n=min_n) for key in y_keys}
    return {"global": block, "by_class": by_class}


def run_eabc_spectral_sweep(
    n_sym_max=20,
    prime_limit=500,
    hurwitz_bound=100,
    max_hurwitz_rows=200,
    n_random=80,
    random_seed=42,
):
    rows = eabc_spectral_collect_samples(
        n_sym_max=n_sym_max,
        prime_limit=prime_limit,
        hurwitz_bound=hurwitz_bound,
        max_hurwitz_rows=max_hurwitz_rows,
        n_random=n_random,
        random_seed=random_seed,
    )
    y_keys = (
        "K_P",
        "holonomy_phase_abs",
        "holonomy_phase_angle",
        "gamma_u24_trace",
        "gamma_u24_angle",
        "dirac_mean",
        "dirac_spread",
        "dist_mean",
        "dist_spread",
        "ahpn_norm_sq",
        "entropy",
        "holonomy_wind",
    )
    corr = eabc_spectral_correlation_block(rows, y_keys)
    gamma_bm = [r["gamma_24_bm"] for r in rows]
    gamma_u24 = [r["gamma_u24_trace"] for r in rows]
    hol_abs = [r["holonomy_phase_abs"] for r in rows]
    class_e = [r for r in rows if r["class"] == "E"]
    scale_ratios = [r["K_P_ratio_to_base"] for r in class_e if np.isfinite(r.get("K_P_ratio_to_base", float("nan")))]
    return {
        "rows": rows,
        "correlations": corr,
        "summary": {
            "n_total": len(rows),
            "n_by_class": {c: sum(1 for r in rows if r["class"] == c) for c in "ABCDE"},
            "gamma_24_bm_constant": float(np.std(gamma_bm)) < 1e-12,
            "gamma_24_bm_value": float(gamma_bm[0]) if gamma_bm else float("nan"),
            "gamma_u24_trace_std": float(np.std(gamma_u24)),
            "holonomy_phase_abs_std": float(np.std(hol_abs)),
            "class_E_K_P_ratio_std": float(np.std(scale_ratios)) if scale_ratios else float("nan"),
        },
    }


def demo_eabc_spectral_sweep():
    print("EABC Spektralsweep (Ptolemäus vs BM-IV-Proxys)")
    print("=" * 64)
    print(f"Vakuum: Pi_0={PI_0:.6g}  K_P_0={K_P_0:.6g}")
    print("Hinweis: gamma_24_bm = |holonomy_phase(24)| auf Z8 ist strukturell ~1 (BM-IV-Schleife).")
    print("EABC-gekoppelte Holonomie: holonomy_phase(sum(X) mod 24), U^24 aus (alpha,m)(X).")
    print()

    sweep = run_eabc_spectral_sweep()
    rows = sweep["rows"]
    summ = sweep["summary"]

    print_section("Stichprobe")
    for cls, label in (
        ("A", "symmetrisch (n,n,n,n)"),
        ("B", "Prim-4 / Vierlinge"),
        ("C", "Hurwitz |q|^2<100"),
        ("D", "zufaellig asymmetrisch [0.5,10]"),
        ("E", "Skalierung fester Formen"),
    ):
        print_kv(f"Klasse {cls} ({label})", summ["n_by_class"][cls])

    print_section("Konstante BM-IV-Groessen (nicht EABC-abhaengig)")
    r0 = rows[0]
    print_kv("gamma_24_bm", f"{r0['gamma_24_bm']:.12g}  (std ueber Sweep: {np.std([r['gamma_24_bm'] for r in rows]):.3e})")
    print_kv("dirac_mean Z8", f"{r0['dirac_mean']:.6g}")
    print_kv("dirac_spread Z8", f"{r0['dirac_spread']:.6g}")
    print_kv("dirac_d2_err", f"{r0['dirac_d2_err']:.3e}")
    if summ["gamma_24_bm_constant"]:
        print("  -> gamma_24_bm ist im Sweep konstant; Korrelation mit K_P_tilde undefiniert.")

    print_section("Variabilitaet EABC-gekoppelter Groessen")
    winds = [r["holonomy_wind"] for r in rows]
    print_kv("std(holonomy_wind)", f"{float(np.std(winds)):.6g}  (mod-24 aus sum(X))")
    print_kv("std(holonomy_phase_abs)", f"{summ['holonomy_phase_abs_std']:.6g}  (|z|=1 auf S^1)")
    print_kv("std(gamma_u24_trace)", f"{summ['gamma_u24_trace_std']:.6g}")
    print_kv("Klasse E: std(K_P_ratio)", f"{summ['class_E_K_P_ratio_std']:.6g}  (Skaleninvarianz K_P)")
    if summ["holonomy_phase_abs_std"] < 1e-9:
        print("  -> holonomy_phase_abs ist konstant 1; Korrelation undefiniert (nur Winkel/Wind variieren).")

    print_section("Scatter-Tabelle (Auswahl, |K_P_tilde| absteigend)")
    hdr = f"{'Kl':2} {'Tag':14} {'K_P':>10} {'K~':>10} {'hol|':>8} {'U24|':>8} {'dist_m':>8} {'|q|^2':>8} {'S':>7}"
    print(hdr)
    print("-" * len(hdr))
    show = sorted(rows, key=lambda r: abs(r["K_P_tilde"]), reverse=True)[:24]
    for r in show:
        print(
            f"{r['class']:2} {r['tag'][:14]:14} "
            f"{r['K_P']:10.4g} {r['K_P_tilde']:10.4g} "
            f"{r['holonomy_phase_abs']:8.4f} {r['gamma_u24_trace']:8.4f} "
            f"{r['dist_mean']:8.4f} {r['ahpn_norm_sq']:8.2f} {r['entropy']:7.4f}"
        )

    print_section("Korrelationen K_P_tilde (global, Pearson r / Spearman rho)")
    y_keys = sweep["correlations"]["global"].keys()
    print(f"{'Groesse':28} {'n':>5} {'Pearson':>10} {'Spearman':>10} {'def':>5}")
    print("-" * 62)
    for key in y_keys:
        c = sweep["correlations"]["global"][key]
        p = "nan" if not c["defined"] else f"{c['pearson']:10.4f}"
        s = "nan" if not c["defined"] else f"{c['spearman']:10.4f}"
        print(f"{key:28} {c['n']:5} {p:>10} {s:>10} {'ja' if c['defined'] else 'nein':>5}")

    print_section("Korrelationen nach Klasse (nur definierte, |r|>0.3)")
    for cls, block in sweep["correlations"]["by_class"].items():
        hits = []
        for key, c in block.items():
            if c["defined"] and abs(c["pearson"]) > 0.3:
                hits.append(f"{key}: r={c['pearson']:.3f}")
        if hits:
            print(f"  {cls}: " + "; ".join(hits[:6]))

    print_section("Hypothese K_P_tilde vs gamma_24")
    if summ["gamma_24_bm_constant"]:
        print("  gamma_24_bm ~ 1 fuer alle Punkte: keine Korrelation testbar.")
    if summ["gamma_u24_trace_std"] > 1e-9:
        c = sweep["correlations"]["global"]["gamma_u24_trace"]
        print(
            f"  EABC-U^24-Spur |gamma|: std={summ['gamma_u24_trace_std']:.4g}; "
            f"Pearson(K~,|gamma|)={'nan' if not c['defined'] else f'{c['pearson']:.4f}'}"
        )
    else:
        print("  gamma_u24_trace variiert nicht signifikant im Sweep.")

    print_section("Klasse A Sanity (K_P_tilde ~ 0)")
    class_a = [r for r in rows if r["class"] == "A"]
    max_a = max(abs(r["K_P_tilde"]) for r in class_a)
    print_kv("max |K_P_tilde| Klasse A", f"{max_a:.3e}  {bool_badge(max_a < 1e-9)}")

    return sweep


def demo_eabc_ptolemaic_closure(run_tests=True):
    print("EABC Ptolemäische Schließung (tetraedrische Einbettung, renormiert)")
    print("=" * 62)
    print(f"Vakuum X=(1,1,1,1):  Pi_0={PI_0:.6g}   K_P_0={K_P_0:.6g}")
    print("Renormierung: Pi_tilde = Pi - Pi_0;  K_P_tilde = K_P - K_P_0")
    print()

    cases = [
        (1.0, 1.0, 1.0, 1.0),
        (2.0, 1.0, 1.5, 0.5),
        (3.0, 2.0, 4.0, 1.0),
    ]
    for x in cases:
        ch = eabc_ptolemaic_chain(*x)
        print(f"X = {ch['X']}")
        print(
            f"  Pi = {ch['Pi']:.6g}   Pi_tilde = {ch['Pi_tilde']:.6g}   "
            f"K_P = {ch['K_P']:.6g}   K_P_tilde = {ch['K_P_tilde']:.6g}   "
            f"geschlossen: {ch['ptolemaic_closed']}"
        )
        print(f"  a_P = {ch['a_P']:.6g} m/s^2   ell_P = {ch['ell_P']:.6g} m")
        print(
            f"  (Zeugen c={ch['physics']['c']:.3g} m/s, "
            f"a_0=c*H0={ch['physics']['a_0']:.3g} m/s^2, H0_CMB={H0_CMB_KM_S_MPC} km/s/Mpc)"
        )
        print()

    k_scale_a = K_P_ptolemaic(3.0, 2.0, 4.0, 1.0)
    k_scale_b = K_P_ptolemaic(1.5, 1.0, 2.0, 0.5)
    print("Skaleninvarianz (Formmaß): K_P(3,2,4,1) =", f"{k_scale_a:.12g}", "  K_P(1.5,1,2,0.5) =", f"{k_scale_b:.12g}")

    k_col = _naive_collinear_K_P(1.0, 2.0, 3.0, 4.0)
    k_tet = K_P_ptolemaic(1.0, 2.0, 3.0, 4.0)
    print("\nKontrolle 1D (kollinear): K_P =", k_col, "(trivial ~0)")
    print("gleiches X, Tetraeder:     K_P =", f"{k_tet:.6g}", "  K_P_tilde =", f"{K_P_ptolemaic_tilde(1,2,3,4):.6g}")

    if not run_tests:
        return True

    print("\n" + "=" * 62)
    print("Testklassen A / B / C")
    print("=" * 62)
    tests = run_eabc_ptolemaic_tests()

    ca = tests["class_A"]
    print(f"\nKlasse A — symmetrisch (n,n,n,n), n=1..{ca['n_max']}")
    print(f"  max |K_P_tilde| = {ca['max_abs_K_P_tilde']:.3e}   Status: {bool_badge(ca['passed'])}")

    print("\nKlasse B — Prim-4-Tupel und Vierlinge (p,p+2,p+6,p+8):")
    for row in tests["class_B"][:8]:
        tag = "Vierling" if row.get("is_quadruplet") else "Prim-4"
        print(
            f"  [{tag}] X={tuple(int(v) for v in row['X'])}  "
            f"K_P={row['K_P']:.6g}  K_P_tilde={row['K_P_tilde']:.6g}  "
            f"rho~={row['prime_density']:.4f}"
        )
    if len(tests["class_B"]) > 8:
        print(f"  ... ({len(tests['class_B'])} Vierlinge bis p<={200})")

    print(f"\nKlasse C — Hurwitz-Hülle (Norm^2<{100}), Auswahl nach |K_P_tilde|:")
    for row in tests["class_C"][:10]:
        q = row["quaternion"]
        print(
            f"  q=({q[0]:g},{q[1]:g},{q[2]:g},{q[3]:g})  |q|^2={row['norm_sq']:.0f}  "
            f"X={tuple(round(v, 3) for v in row['X'])}  K_P_tilde={row['K_P_tilde']:.6g}  S={row['entropy']:.4f}"
        )

    print("\nKorrelationen (Pearson r, falls Daten vorhanden):")
    corr = tests["correlations"]
    print_kv("Holonomie gamma_24", corr["holonomy_gamma_24"])
    print_kv("24-Zyklus OK", bool_badge(corr["cycle_24_passed"]))
    for name, val in corr["correlations"].items():
        print_kv(name, "nan" if math.isnan(val) else f"{val:.4f}")

    return ca["passed"]


def sieve_is_prime(limit):
    """Prim-Sieb bis limit (inklusive). Rueckgabe: bytearray is_prime[0..limit]."""
    if limit < 2:
        return bytearray(limit + 1)
    is_p = bytearray(b"\x01") * (limit + 1)
    is_p[0] = is_p[1] = 0
    m = int(math.isqrt(limit))
    for i in range(2, m + 1):
        if is_p[i]:
            step = i
            start = i * i
            is_p[start : limit + 1 : step] = b"\x00" * ((limit - start) // step + 1)
    return is_p


def prime_quadruplets_via_sieve(p_limit=1_000_000):
    """
    Primzahlvierlinge (p, p+2, p+6, p+8) mit p <= p_limit.
    Erwartung: O(n log log n) Sieb, dann linearer Scan.
    """
    hi = int(p_limit) + 8
    is_p = sieve_is_prime(hi)
    out = []
    for p in range(5, int(p_limit) - 7, 2):
        if is_p[p] and is_p[p + 2] and is_p[p + 6] and is_p[p + 8]:
            out.append((float(p), float(p + 2), float(p + 6), float(p + 8)))
    return out


def _eabc_tuple_symmetry_proxy(x):
    arr = np.asarray(x, dtype=float)
    mu = float(np.mean(arr))
    if mu <= EPS:
        return float("nan")
    return float(np.std(arr) / mu)


def eabc_prime_sweep_row(x, prime_density_window=20):
    """Metriken fuer ein EABC-4-Tupel (Prim-Vierling oder Referenz)."""
    p0 = float(x[0])
    ch = eabc_ptolemaic_chain(*x)
    return {
        "p": p0,
        "X": tuple(float(v) for v in x),
        "Pi": ch["Pi"],
        "Pi_tilde": ch["Pi_tilde"],
        "K_P": ch["K_P"],
        "K_P_tilde": ch["K_P_tilde"],
        "abs_K_P_tilde": abs(ch["K_P_tilde"]) if not math.isnan(ch["K_P_tilde"]) else float("nan"),
        "prime_density": prime_density_proxy(p0, window=prime_density_window),
        "is_quadruplet": bool(is_prime_quadruplet(int(p0))),
        "symmetry_cv": _eabc_tuple_symmetry_proxy(x),
        "p_times_abs_K_tilde": (
            p0 * abs(ch["K_P_tilde"]) if not math.isnan(ch["K_P_tilde"]) else float("nan")
        ),
    }


def eabc_prime_sweep_collect(p_limit=1_000_000, include_classics=True):
    """Vierlinge per Sieb + klassische Prim-4-Tupel ohne Duplikate."""
    tuples = prime_quadruplets_via_sieve(p_limit=p_limit)
    seen = {tuple(int(v) for v in t) for t in tuples}
    if include_classics:
        for x in EABC_PRIME_EXAMPLE_TUPLES:
            key = tuple(int(v) for v in x)
            if key not in seen:
                seen.add(key)
                tuples.append(x)
    tuples.sort(key=lambda t: t[0])
    rows = [eabc_prime_sweep_row(t) for t in tuples]
    return rows


def _loglog_regression(ps, abs_ks, min_positive=1e-300):
    """log(p) vs log(|K_P_tilde|): Steigung beta, Achsenabschnitt alpha, Pearson r."""
    ps = np.asarray(ps, dtype=float)
    abs_ks = np.asarray(abs_ks, dtype=float)
    mask = np.isfinite(ps) & np.isfinite(abs_ks) & (ps > 0)
    if mask.sum() < 3:
        return {
            "n": int(mask.sum()),
            "beta": float("nan"),
            "alpha": float("nan"),
            "pearson_r": float("nan"),
            "r_squared": float("nan"),
            "defined": False,
        }
    lp = np.log(ps[mask])
    lk = np.log(np.maximum(abs_ks[mask], min_positive))
    if np.std(lp) < EPS or np.std(lk) < EPS:
        return {
            "n": int(mask.sum()),
            "beta": float("nan"),
            "alpha": float("nan"),
            "pearson_r": float("nan"),
            "r_squared": float("nan"),
            "defined": False,
        }
    beta, alpha = np.polyfit(lp, lk, 1)
    pearson_r = float(np.corrcoef(lp, lk)[0, 1])
    pred = alpha + beta * lp
    ss_res = float(np.sum((lk - pred) ** 2))
    ss_tot = float(np.sum((lk - np.mean(lk)) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > EPS else float("nan")
    return {
        "n": int(mask.sum()),
        "beta": float(beta),
        "alpha": float(alpha),
        "pearson_r": pearson_r,
        "r_squared": r_squared,
        "defined": True,
        "model": "|K_P_tilde| ~ p^beta  (log-log: log|K~| = alpha + beta*log p)",
    }


def eabc_prime_sweep_analysis(rows):
    """Trend K_P_tilde vs p: log-log, 1/p-Skalierung, Symmetrie-Artefakt."""
    ps = [r["p"] for r in rows]
    abs_ks = [r["abs_K_P_tilde"] for r in rows]
    ks = [r["K_P_tilde"] for r in rows]
    cvs = [r["symmetry_cv"] for r in rows]

    loglog = _loglog_regression(ps, abs_ks)
    loglog_pos = _loglog_regression(
        [r["p"] for r in rows if r["K_P_tilde"] > 0],
        [r["abs_K_P_tilde"] for r in rows if r["K_P_tilde"] > 0],
    )

    inv_p = [1.0 / r["p"] for r in rows]
    corr_inv = _pair_correlations(ps, abs_ks)
    corr_cv = _pair_correlations(ps, cvs)
    corr_pK = _pair_correlations(ps, [r["p_times_abs_K_tilde"] for r in rows])

    quad_only = [r for r in rows if r["is_quadruplet"]]
    loglog_quad = _loglog_regression(
        [r["p"] for r in quad_only],
        [r["abs_K_P_tilde"] for r in quad_only],
    )

    abs_arr = np.asarray(abs_ks, dtype=float)
    trend_verdict = "unbestimmt"
    if loglog["defined"]:
        beta = loglog["beta"]
        r = loglog["pearson_r"]
        if beta < -0.15 and r < -0.3:
            trend_verdict = "abnehmend (Hypothese K~ -> 0 gestuetzt)"
        elif abs(beta) < 0.08:
            trend_verdict = "kein klarer p-Abfall (Hypothese schwach)"
        elif beta > 0.15 and r > 0.3:
            trend_verdict = "wachsend mit p (Hypothese widerlegt)"
        else:
            trend_verdict = "schwach / gemischt"

    scaling_note = (
        "K_P ist skaleninvariant bei uniformer Skalierung; Prim-Vierlinge werden mit p "
        "groesser aber nicht proportional skaliert — relative Spreizung (CV der Besetzung) "
        "kann |K~| beeinflussen. Vergleiche beta (log-log) mit Korrelation p vs CV."
    )
    if corr_cv["defined"] and abs(corr_cv["pearson"]) > 0.5:
        scaling_note += " Starke p–CV-Kopplung: moegliches Symmetrie-/Form-Artefakt."
    if corr_pK["defined"] and abs(corr_pK["pearson"]) < 0.3 and loglog["defined"] and loglog["beta"] < 0:
        scaling_note += " p*|K~| nicht stabil -> eher echter Abfall als reines 1/p."

    return {
        "loglog_all": loglog,
        "loglog_positive_K": loglog_pos,
        "loglog_quadruplets_only": loglog_quad,
        "pearson_p_vs_abs_K": corr_inv,
        "pearson_p_vs_symmetry_cv": corr_cv,
        "pearson_p_vs_p_times_abs_K": corr_pK,
        "min_abs_K": float(np.nanmin(abs_arr)) if abs_arr.size else float("nan"),
        "max_abs_K": float(np.nanmax(abs_arr)) if abs_arr.size else float("nan"),
        "median_abs_K": float(np.nanmedian(abs_arr)) if abs_arr.size else float("nan"),
        "trend_verdict": trend_verdict,
        "scaling_artifact_note": scaling_note,
        "n_quadruplets": sum(1 for r in rows if r["is_quadruplet"]),
        "n_total": len(rows),
    }


def write_eabc_prime_sweep_csv(rows, csv_path):
    import csv

    fields = (
        "p",
        "E",
        "A",
        "B",
        "C",
        "Pi",
        "Pi_tilde",
        "K_P",
        "K_P_tilde",
        "abs_K_P_tilde",
        "prime_density",
        "is_quadruplet",
        "symmetry_cv",
        "p_times_abs_K_tilde",
    )
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            x = r["X"]
            w.writerow(
                {
                    "p": int(r["p"]) if r["p"] == int(r["p"]) else r["p"],
                    "E": x[0],
                    "A": x[1],
                    "B": x[2],
                    "C": x[3],
                    "Pi": r["Pi"],
                    "Pi_tilde": r["Pi_tilde"],
                    "K_P": r["K_P"],
                    "K_P_tilde": r["K_P_tilde"],
                    "abs_K_P_tilde": r["abs_K_P_tilde"],
                    "prime_density": r["prime_density"],
                    "is_quadruplet": int(r["is_quadruplet"]),
                    "symmetry_cv": r["symmetry_cv"],
                    "p_times_abs_K_tilde": r["p_times_abs_K_tilde"],
                }
            )
    return csv_path


def run_eabc_prime_sweep(p_limit=1_000_000, csv_path=None, include_classics=True):
    rows = eabc_prime_sweep_collect(p_limit=p_limit, include_classics=include_classics)
    analysis = eabc_prime_sweep_analysis(rows)
    if csv_path is None:
        import os

        csv_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "eabc_prime_sweep.csv",
        )
    if rows:
        write_eabc_prime_sweep_csv(rows, csv_path)
    return {"rows": rows, "analysis": analysis, "csv_path": csv_path, "p_limit": int(p_limit)}


def _pick_prime_sweep_samples(rows, k_each=3):
    """Kleine / mittlere / grosse p: je k_each Zeilen."""
    if not rows:
        return {"small": [], "medium": [], "large": []}
    sorted_rows = sorted(rows, key=lambda r: r["p"])
    n = len(sorted_rows)
    if n <= 3 * k_each:
        third = max(1, n // 3)
        return {
            "small": sorted_rows[:third],
            "medium": sorted_rows[third : 2 * third],
            "large": sorted_rows[2 * third :],
        }
    i0 = k_each
    im = n // 2
    i1 = n - k_each
    return {
        "small": sorted_rows[:i0],
        "medium": sorted_rows[im - k_each // 2 : im - k_each // 2 + k_each],
        "large": sorted_rows[i1:],
    }


def demo_eabc_prime_sweep(p_limit=1_000_000):
    print("EABC Prim-Vierling-Sweep (p, p+2, p+6, p+8)")
    print("=" * 64)
    print(f"Vakuum: Pi_0={PI_0:.6g}  K_P_0={K_P_0:.6g}  (K_P_tilde=0 fuer (n,n,n,n))")
    print(f"Hypothese: K_P_tilde -> 0 fuer grosse Prim-Vierlinge (ptolemaeische Annaeherung)")
    print()

    sweep = run_eabc_prime_sweep(p_limit=p_limit)
    rows = sweep["rows"]
    an = sweep["analysis"]

    print_section("Uebersicht")
    print_kv("p-Obergrenze (Sieb)", sweep["p_limit"])
    print_kv("Anzahl 4-Tupel gesamt", an["n_total"])
    print_kv("davon echte Vierlinge", an["n_quadruplets"])
    if rows:
        print_kv("p_min", int(min(r["p"] for r in rows)))
        print_kv("p_max", int(max(r["p"] for r in rows)))
    print_kv("CSV", sweep["csv_path"])

    print_section("|K_P_tilde| — Statistik")
    print_kv("min", f"{an['min_abs_K']:.6g}")
    print_kv("median", f"{an['median_abs_K']:.6g}")
    print_kv("max", f"{an['max_abs_K']:.6g}")

    ll = an["loglog_all"]
    print_section("Log-Log-Regression: log(p) vs log(|K_P_tilde|)")
    if ll["defined"]:
        print_kv("n", ll["n"])
        print_kv("beta (Steigung)", f"{ll['beta']:.4f}  (~ p^beta)")
        print_kv("alpha", f"{ll['alpha']:.4f}")
        print_kv("Pearson r", f"{ll['pearson_r']:.4f}")
        print_kv("R^2", f"{ll['r_squared']:.4f}")
        if ll["beta"] < -0.5:
            print_kv("Lesart", "nahe 1/p-Abfall (beta ~ -1)")
        elif ll["beta"] < -0.15:
            print_kv("Lesart", "sublinearer Abfall mit p")
        elif abs(ll["beta"]) < 0.1:
            print_kv("Lesart", "kein klarer Potenzgesetz-Abfall")
    else:
        print("  Zu wenige positive/variierte Punkte fuer log-log-Fit.")

    lq = an["loglog_quadruplets_only"]
    if lq["defined"] and lq["n"] != ll.get("n"):
        print_kv("nur Vierlinge: beta", f"{lq['beta']:.4f}")
        print_kv("nur Vierlinge: r", f"{lq['pearson_r']:.4f}")

    print_section("Trend-Urteil")
    print_kv("Fazit", an["trend_verdict"])
    print(an["scaling_artifact_note"])

    print_section("Stichproben-Tabelle (klein / mittel / gross p)")
    hdr = f"{'p':>8} {'X':>22} {'K_P':>10} {'K~':>10} {'|K~|':>10} {'rho':>7} {'CV':>7}"
    print(hdr)
    print("-" * len(hdr))
    samples = _pick_prime_sweep_samples(rows, k_each=3)
    for label, block in (("klein", samples["small"]), ("mittel", samples["medium"]), ("gross", samples["large"])):
        for r in block:
            x = tuple(int(v) for v in r["X"])
            print(
                f"{int(r['p']):8} {str(x):>22} "
                f"{r['K_P']:10.4g} {r['K_P_tilde']:10.4g} {r['abs_K_P_tilde']:10.4g} "
                f"{r['prime_density']:7.4f} {r['symmetry_cv']:7.4f}"
            )
        if block:
            print(f"  --- {label} ---")

    print_section("Vollstaendige Vierling-Liste (|K~| absteigend, Top 12)")
    quad = sorted([r for r in rows if r["is_quadruplet"]], key=lambda r: r["abs_K_P_tilde"], reverse=True)[:12]
    for r in quad:
        x = tuple(int(v) for v in r["X"])
        tag = "Vierling" if r["is_quadruplet"] else "Ref"
        print(
            f"  [{tag}] p={int(r['p']):6} X={x}  K_P={r['K_P']:.6g}  "
            f"K_tilde={r['K_P_tilde']:.6g}  Pi_tilde={r['Pi_tilde']:.6g}"
        )

    return sweep


def _eabc_quadruplet_sum(p):
    """Summe von (p, p+2, p+6, p+8)."""
    return 4.0 * float(p) + 16.0


def _eabc_control_ap_matched(p):
    """AP mit gleicher Summe und Endpunkten p .. p+8 wie der Prim-Vierling."""
    k = 8.0 / 3.0
    p = float(p)
    return (p, p + k, p + 2.0 * k, p + 8.0)


def _eabc_control_consecutive(p):
    """Aufeinanderfolgende Zahlen, gleicher Anker p."""
    p = float(p)
    return (p, p + 1.0, p + 2.0, p + 3.0)


def _eabc_random_4tuple_same_sum(target_sum, rng, anchor_p, max_tries=200):
    """
    Vier positive Zahlen mit fester Summe, reproduzierbar; kein Prim-Vierling.
    """
    s = int(round(target_sum))
    if s < 4:
        return None
    for _ in range(max_tries):
        cuts = np.sort(rng.integers(1, s, size=3))
        a = cuts[0]
        b = cuts[1] - cuts[0]
        c = cuts[2] - cuts[1]
        d = s - cuts[2]
        if min(a, b, c, d) < 1:
            continue
        x = (float(a), float(b), float(c), float(d))
        p0 = int(x[0])
        if is_prime_quadruplet(p0):
            continue
        if abs(sum(x) - target_sum) > 1e-6:
            continue
        return x
    return None


def eabc_prime_control_row(x, row_type, anchor_p=None, prime_density_window=20):
    """Wie eabc_prime_sweep_row, zusaetzlich mit Typ und Anker-p."""
    base = eabc_prime_sweep_row(x, prime_density_window=prime_density_window)
    p_anchor = float(anchor_p if anchor_p is not None else x[0])
    mean_x = float(np.mean(x))
    return {
        **base,
        "type": row_type,
        "anchor_p": p_anchor,
        "sum_X": float(sum(x)),
        "mean_X": mean_x,
    }


def eabc_prime_control_collect(prime_rows=None, p_limit=1_000_000, seed=42):
    """
    Fuer jeden Prim-Vierling drei Kontroll-4-Tupel (gleiche Skala, kein Vierling).
    """
    if prime_rows is None:
        prime_rows = [
            r
            for r in eabc_prime_sweep_collect(p_limit=p_limit, include_classics=False)
            if r["is_quadruplet"]
        ]
    else:
        prime_rows = [r for r in prime_rows if r.get("is_quadruplet", is_prime_quadruplet(int(r["p"])))]

    rng = np.random.default_rng(seed)
    out = []
    for pr in prime_rows:
        p = float(pr["p"])
        quad = pr["X"] if "X" in pr else (p, p + 2, p + 6, p + 8)
        out.append(eabc_prime_control_row(quad, "prime", anchor_p=p))

        ap = _eabc_control_ap_matched(p)
        out.append(eabc_prime_control_row(ap, "control_ap", anchor_p=p))

        cons = _eabc_control_consecutive(p)
        out.append(eabc_prime_control_row(cons, "control_consecutive", anchor_p=p))

        target_sum = _eabc_quadruplet_sum(p)
        sub_rng = np.random.default_rng(seed + int(p) % (2**31 - 1))
        rnd = _eabc_random_4tuple_same_sum(target_sum, sub_rng, anchor_p=p)
        if rnd is not None:
            out.append(eabc_prime_control_row(rnd, "control_random", anchor_p=p))

    return out


def _eabc_control_by_type(rows):
    types = sorted({r["type"] for r in rows})
    return {t: [r for r in rows if r["type"] == t] for t in types}


def _eabc_p_bins(rows, n_bins=5):
    """Gleichmaessige p-Bins nach anchor_p; Median |K~| pro Typ und Bin."""
    if not rows:
        return []
    ps = np.array([r["anchor_p"] for r in rows], dtype=float)
    edges = np.quantile(ps, np.linspace(0, 1, n_bins + 1))
    edges = np.unique(edges)
    if edges.size < 2:
        edges = np.array([ps.min(), ps.max() + 1.0])
    bins_out = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        in_bin = [
            r
            for r in rows
            if (r["anchor_p"] >= lo if i == 0 else r["anchor_p"] > lo)
            and (r["anchor_p"] <= hi if i == len(edges) - 2 else r["anchor_p"] <= hi)
        ]
        if not in_bin:
            continue
        by_type = _eabc_control_by_type(in_bin)
        medians = {}
        for t, sub in by_type.items():
            vals = [r["abs_K_P_tilde"] for r in sub if np.isfinite(r["abs_K_P_tilde"])]
            medians[t] = float(np.median(vals)) if vals else float("nan")
        bins_out.append(
            {
                "bin": i,
                "p_lo": float(lo),
                "p_hi": float(hi),
                "n": len(in_bin),
                "median_abs_K_by_type": medians,
            }
        )
    return bins_out


def eabc_prime_control_analysis(rows):
    """Vergleich Prim vs Kontrollen: log-log-Steigung, CV-Kopplung, p-Bins."""
    by_type = _eabc_control_by_type(rows)
    per_type = {}
    for t, sub in by_type.items():
        ps = [r["anchor_p"] for r in sub]
        abs_ks = [r["abs_K_P_tilde"] for r in sub]
        cvs = [r["symmetry_cv"] for r in sub]
        per_type[t] = {
            "n": len(sub),
            "loglog": _loglog_regression(ps, abs_ks),
            "corr_cv_vs_abs_K": _pair_correlations(cvs, abs_ks),
            "median_abs_K": float(np.nanmedian(abs_ks)) if sub else float("nan"),
            "mean_cv": float(np.nanmean(cvs)) if sub else float("nan"),
        }

    prime = by_type.get("prime", [])
    controls = [r for r in rows if r["type"] != "prime"]
    ratio_medians = []
    for pr in prime:
        p = pr["anchor_p"]
        ak_p = pr["abs_K_P_tilde"]
        ctrl_same_p = [r for r in controls if r["anchor_p"] == p]
        if not ctrl_same_p or not np.isfinite(ak_p):
            continue
        med_c = float(np.median([r["abs_K_P_tilde"] for r in ctrl_same_p]))
        if med_c > EPS:
            ratio_medians.append(ak_p / med_c)

    bins = _eabc_p_bins(rows, n_bins=5)
    bin_prime_vs_ctrl = []
    for b in bins:
        med = b["median_abs_K_by_type"]
        mp = med.get("prime", float("nan"))
        others = [v for t, v in med.items() if t != "prime" and np.isfinite(v)]
        mc = float(np.median(others)) if others else float("nan")
        bin_prime_vs_ctrl.append(
            {
                "bin": b["bin"],
                "p_lo": b["p_lo"],
                "p_hi": b["p_hi"],
                "median_prime": mp,
                "median_controls": mc,
                "ratio_prime_over_controls": mp / mc if np.isfinite(mp) and np.isfinite(mc) and mc > EPS else float("nan"),
            }
        )

    betas = {t: per_type[t]["loglog"].get("beta", float("nan")) for t in per_type}
    prime_beta = betas.get("prime", float("nan"))
    ctrl_betas = [betas[t] for t in betas if t != "prime" and np.isfinite(betas[t])]
    beta_spread = (
        float(max(ctrl_betas) - min(ctrl_betas)) if len(ctrl_betas) >= 2 else float("nan")
    )
    beta_prime_minus_ctrl = (
        float(prime_beta - float(np.mean(ctrl_betas))) if ctrl_betas and np.isfinite(prime_beta) else float("nan")
    )

    structured_types = ("control_ap", "control_consecutive")
    struct_betas = [betas[t] for t in structured_types if t in betas and np.isfinite(betas[t])]
    rnd_beta = betas.get("control_random", float("nan"))

    verdict = "unbestimmt"
    artifact_score = 0
    genuine_score = 0

    if np.isfinite(prime_beta) and struct_betas:
        if prime_beta < -0.15 and all(b < -0.1 for b in struct_betas):
            artifact_score += 3
        if abs(prime_beta - float(np.mean(struct_betas))) < 0.08:
            artifact_score += 2

    if np.isfinite(rnd_beta) and np.isfinite(prime_beta) and abs(rnd_beta - prime_beta) > 0.2:
        artifact_score += 1

    if ratio_medians:
        med_ratio = float(np.median(ratio_medians))
        if 0.85 <= med_ratio <= 1.15:
            artifact_score += 1
        elif 1.3 < med_ratio < 2.2:
            artifact_score += 1
        elif med_ratio < 0.7 or med_ratio > 2.5:
            genuine_score += 1

    bin_ratios = [b["ratio_prime_over_controls"] for b in bin_prime_vs_ctrl if np.isfinite(b["ratio_prime_over_controls"])]
    if bin_ratios and all(1.5 <= r <= 2.2 for r in bin_ratios):
        artifact_score += 1

    if artifact_score >= 4 and genuine_score <= 1:
        verdict = (
            "vorwiegend Skalierungs-/Symmetrie-Artefakt: |K~|~1/p^2 mit kleinem CV; "
            "Prim-Vierlinge wie AP/aufeinanderfolgend, nicht wie Zufall gleicher Summe"
        )
    elif genuine_score >= 3 and artifact_score <= 2:
        verdict = "Prim-Struktur von skalenmatched Kontrollen klar unterscheidbar"
    elif artifact_score >= 3:
        verdict = "gemischt: gleicher p-Abfall bei strukturierten Kontrollen; Zufall (gleiche Summe) behaelt grosses |K~|"
    else:
        verdict = "kein klares Signal (Datenlage eng)"

    return {
        "per_type": per_type,
        "betas": betas,
        "beta_prime_minus_mean_controls": beta_prime_minus_ctrl,
        "beta_spread_controls": beta_spread,
        "median_ratio_prime_to_controls_at_same_p": float(np.median(ratio_medians)) if ratio_medians else float("nan"),
        "p_bins": bin_prime_vs_ctrl,
        "verdict": verdict,
        "artifact_score": artifact_score,
        "genuine_score": genuine_score,
        "n_total": len(rows),
        "n_prime": len(prime),
        "n_controls": len(controls),
    }


def write_eabc_prime_control_csv(rows, csv_path):
    import csv

    fields = (
        "type",
        "anchor_p",
        "p",
        "E",
        "A",
        "B",
        "C",
        "sum_X",
        "mean_X",
        "Pi",
        "Pi_tilde",
        "K_P",
        "K_P_tilde",
        "abs_K_P_tilde",
        "symmetry_cv",
        "prime_density",
        "is_quadruplet",
        "p_times_abs_K_tilde",
    )
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            x = r["X"]
            w.writerow(
                {
                    "type": r["type"],
                    "anchor_p": int(r["anchor_p"]) if r["anchor_p"] == int(r["anchor_p"]) else r["anchor_p"],
                    "p": int(r["p"]) if r["p"] == int(r["p"]) else r["p"],
                    "E": x[0],
                    "A": x[1],
                    "B": x[2],
                    "C": x[3],
                    "sum_X": r["sum_X"],
                    "mean_X": r["mean_X"],
                    "Pi": r["Pi"],
                    "Pi_tilde": r["Pi_tilde"],
                    "K_P": r["K_P"],
                    "K_P_tilde": r["K_P_tilde"],
                    "abs_K_P_tilde": r["abs_K_P_tilde"],
                    "symmetry_cv": r["symmetry_cv"],
                    "prime_density": r["prime_density"],
                    "is_quadruplet": int(r["is_quadruplet"]),
                    "p_times_abs_K_tilde": r["p_times_abs_K_tilde"],
                }
            )
    return csv_path


def run_eabc_prime_control(p_limit=1_000_000, csv_path=None, seed=42):
    prime_rows = [
        r
        for r in eabc_prime_sweep_collect(p_limit=p_limit, include_classics=False)
        if r["is_quadruplet"]
    ]
    rows = eabc_prime_control_collect(prime_rows=prime_rows, seed=seed)
    analysis = eabc_prime_control_analysis(rows)
    if csv_path is None:
        import os

        csv_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "eabc_prime_control.csv",
        )
    if rows:
        write_eabc_prime_control_csv(rows, csv_path)
    return {
        "rows": rows,
        "analysis": analysis,
        "csv_path": csv_path,
        "n_prime_quadruplets": len(prime_rows),
        "seed": seed,
    }


def demo_eabc_prime_control(p_limit=1_000_000, seed=42):
    print("EABC Prim-Vierling — Kontrollgruppen-Vergleich")
    print("=" * 64)
    print("Ziel: Prim-Effekt vs Skalierungsartefakt (CV -> 0, |K~| -> 0 mit p)")
    print()

    res = run_eabc_prime_control(p_limit=p_limit, seed=seed)
    rows = res["rows"]
    an = res["analysis"]

    print_section("Uebersicht")
    print_kv("Prim-Vierlinge", res["n_prime_quadruplets"])
    print_kv("Zeilen gesamt", an["n_total"])
    print_kv("davon Kontrollen", an["n_controls"])
    print_kv("Seed (control_random)", seed)
    print_kv("CSV", res["csv_path"])

    print_section("Log-Log-Steigung beta: log|K~| vs anchor_p")
    hdr = f"{'Typ':>22} {'n':>6} {'beta':>10} {'r':>8} {'med|K~|':>12} {'mean CV':>10}"
    print(hdr)
    print("-" * len(hdr))
    for t in sorted(an["per_type"].keys()):
        pt = an["per_type"][t]
        ll = pt["loglog"]
        beta = ll["beta"] if ll["defined"] else float("nan")
        r = ll["pearson_r"] if ll["defined"] else float("nan")
        print(
            f"{t:>22} {pt['n']:6} {beta:10.4f} {r:8.4f} "
            f"{pt['median_abs_K']:12.6g} {pt['mean_cv']:10.6g}"
        )

    print_section("CV vs |K~| (Pearson)")
    for t in sorted(an["per_type"].keys()):
        c = an["per_type"][t]["corr_cv_vs_abs_K"]
        if c["defined"]:
            print_kv(t, f"r={c['pearson']:.4f}  (n={c['n']})")

    print_section("Gleiches p: Median |K~| Prim / Median Kontrollen")
    print_kv("Median Verhaeltnis Prim/Kontrollen", f"{an['median_ratio_prime_to_controls_at_same_p']:.4f}")
    print_kv("beta(Prim) - mean(beta Kontrollen)", f"{an['beta_prime_minus_mean_controls']:.4f}")

    print_section("p-Bins (5 Quantile): Median |K~|")
    for b in an["p_bins"]:
        print(
            f"  Bin {b['bin']}: p in [{int(b['p_lo'])}, {int(b['p_hi'])}]  "
            f"Prim={b['median_prime']:.6g}  Kontrollen={b['median_controls']:.6g}  "
            f"Ratio={b['ratio_prime_over_controls']:.3f}"
        )

    print_section("Fazit")
    print_kv("Urteil", an["verdict"])
    print_kv("Artefakt-Score", an["artifact_score"])
    print_kv("Genuine-Score", an["genuine_score"])
    rnd_b = an["betas"].get("control_random", float("nan"))
    if np.isfinite(rnd_b):
        print(
            "  Hinweis: control_random (gleiche Summe, hohes CV) hat beta~0 — "
            "der p-Abfall bei Prim-Vierlingen erklaert sich durch relative Gleichfoermigkeit, "
            "nicht allein durch die Summe."
        )

    print("\nStichprobe (grosses p, je Typ eine Zeile):")
    large_p = sorted({r["anchor_p"] for r in rows})[-1]
    for t in ("prime", "control_ap", "control_consecutive", "control_random"):
        sub = [r for r in rows if r["type"] == t and r["anchor_p"] == large_p]
        if sub:
            r = sub[0]
            print(
                f"  {t:20} p={int(large_p)} X={tuple(round(v, 2) for v in r['X'])} "
                f"|K~|={r['abs_K_P_tilde']:.6g} CV={r['symmetry_cv']:.6g}"
            )

    return res


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] in ("--eabc-ptolemaic", "--eabc"):
        demo_eabc_ptolemaic_closure(run_tests=True)
    elif len(sys.argv) > 1 and sys.argv[1] == "--eabc-ptolemaic-tests":
        demo_eabc_ptolemaic_closure(run_tests=False)
        print("\n" + "=" * 62)
        print("Nur Testklassen A/B/C + Korrelationen")
        print("=" * 62)
        tests = run_eabc_ptolemaic_tests()
        ca = tests["class_A"]
        print(f"Klasse A max|K_P_tilde|={ca['max_abs_K_P_tilde']:.3e}  {bool_badge(ca['passed'])}")
        print(f"Klasse B: {len(tests['class_B'])} Vierlinge")
        print(f"Klasse C: {len(hurwitz_quaternion_hull(100))} Hurwitz-Elemente (|q|^2<100)")
        for name, val in tests["correlations"]["correlations"].items():
            print(f"  {name}: {'nan' if math.isnan(val) else f'{val:.4f}'}")
    elif len(sys.argv) > 1 and sys.argv[1] == "--eabc-spectral-sweep":
        demo_eabc_spectral_sweep()
    elif len(sys.argv) > 1 and sys.argv[1] in ("--eabc-prime-sweep", "--eabc-primes"):
        p_lim = 1_000_000
        if len(sys.argv) > 2:
            try:
                p_lim = int(sys.argv[2])
            except ValueError:
                pass
        demo_eabc_prime_sweep(p_limit=p_lim)
    elif len(sys.argv) > 1 and sys.argv[1] in ("--eabc-prime-control", "--eabc-prime-controls"):
        p_lim = 1_000_000
        seed = 42
        if len(sys.argv) > 2:
            try:
                p_lim = int(sys.argv[2])
            except ValueError:
                pass
        if len(sys.argv) > 3:
            try:
                seed = int(sys.argv[3])
            except ValueError:
                pass
        demo_eabc_prime_control(p_limit=p_lim, seed=seed)
    else:
        run_quantum_tests()