
import numpy as np

def unfolded_gaps(zeros: np.ndarray) -> np.ndarray:
    z = np.asarray(zeros, dtype=float)
    dz = z[1:] - z[:-1]
    mean_gap = 2.0 * np.pi / np.log(z[:-1] / (2.0 * np.pi))
    return dz / mean_gap

def build_blocks(u: np.ndarray, N: int, offset: int) -> np.ndarray:
    start = offset
    needed = start + 4 * N
    if needed > len(u):
        raise ValueError("Nicht genug entfaltete Gaps.")
    return np.array([u[start + 4*k:start + 4*k + 4] for k in range(N)], dtype=float)

def toy_local_bias(N: int, offset: int, scale: float = 0.035) -> np.ndarray:
    idx = np.arange(N, dtype=float)
    raw = 0.5 + 0.5 * np.abs(np.sin((idx + 1 + offset) / (N + 1)))
    return scale * raw

def block_cv(q: np.ndarray) -> np.ndarray:
    mu = np.mean(q, axis=1)
    sig = np.std(q, axis=1)
    return sig / (np.abs(mu) + 1e-12)

def block_norm_drift(q: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(q, axis=1) + 1e-12
    drift = np.zeros(len(norms))
    if len(norms) > 1:
        drift[:-1] = np.abs(np.log(norms[1:]) - np.log(norms[:-1]))
        drift[-1] = drift[-2]
    return drift

def eabc_extended_terms(q: np.ndarray):
    eps = 1e-12
    N = len(q)
    E, A, C, G = [], [], [], []
    for k, block in enumerate(q):
        u1, u2, u3, u4 = block
        norm = np.linalg.norm(block) + eps
        E_k = (abs(u1-u3) + abs(u2-u4)) / norm
        A_k = abs((u1+u2) - (u3+u4)) / norm
        G_k = (abs(u1 - 2*u2 + u3) + abs(u2 - 2*u3 + u4)) / norm
        if N == 1:
            C_k = 0.0
        else:
            left = np.linalg.norm(block - q[k-1]) if k > 0 else 0.0
            right = np.linalg.norm(q[k+1] - block) if k < N-1 else 0.0
            C_k = (left + right) / (2.0 * norm)
        E.append(E_k); A.append(A_k); C.append(C_k); G.append(G_k)
    return np.array(E), np.array(A), np.array(C), np.array(G)

def low_observable_base(q: np.ndarray, bias: np.ndarray) -> np.ndarray:
    CV = block_cv(q); F = block_norm_drift(q)
    return 2.0*CV + 2.0*bias + 1.0*F

def low_temperature_from_Sbar(Sbar: float, gamma: float, T_star: float):
    lam_L = np.exp(-gamma * Sbar)
    return lam_L, T_star * lam_L

def storage_term(q: np.ndarray) -> np.ndarray:
    norms2 = np.sum(q*q, axis=1)
    gamma = np.zeros(len(q))
    for k in range(len(q)):
        left = np.linalg.norm(q[k] - q[k-1]) if k > 0 else 0.0
        right = np.linalg.norm(q[k+1] - q[k]) if k < len(q)-1 else 0.0
        gamma[k] = 1e-6 + left + right
    return norms2 / gamma

def coherence_term(q: np.ndarray) -> np.ndarray:
    c = []
    for block in q:
        mean = np.mean(block); std = np.std(block)
        phases = np.array([1.0*mean+0.5*std, 2.0*mean+1.0*std, 3.0*mean+1.5*std, 4.0*mean+2.0*std])
        amp = np.mean(np.exp(1j*phases))
        c.append(np.abs(amp)**2)
    return np.array(c)

def tetra_phase_value(block: np.ndarray) -> complex:
    x1, x2, x3, x4 = block
    faces = [(x2, x3, x4), (x1, x3, x4), (x1, x2, x4), (x1, x2, x3)]
    vals = []
    for a, b, c in faces:
        omega = ((a - b) * (b - c) * (c - a)) / ((a*a + b*b + c*c + 1e-12) ** 1.5)
        vals.append(np.exp(1j * omega))
    return np.mean(vals)

def holonomy_arrays(q: np.ndarray):
    W_T = np.array([tetra_phase_value(block) for block in q], dtype=complex)
    Phi_T = np.angle(W_T)
    W_O = np.zeros_like(W_T, dtype=complex)
    if len(W_T) == 1:
        W_O[0] = 1.0 + 0.0j
    else:
        for k in range(len(W_T)-1):
            W_O[k] = W_T[k] * np.conj(W_T[k+1])
        W_O[-1] = W_O[-2]
    return Phi_T, np.angle(W_O)

def holonomy_term(q: np.ndarray) -> np.ndarray:
    Phi_T, Phi_O = holonomy_arrays(q)
    return (1.0 - np.cos(Phi_T)) + (1.0 - np.cos(Phi_O))

def high_raw_x(q: np.ndarray, alpha_s: float = 1.0, alpha_h: float = 0.3) -> np.ndarray:
    return alpha_s * storage_term(q) * coherence_term(q) + alpha_h * holonomy_term(q)

def high_scaled_y(q: np.ndarray, kappa: float = 8.0, zeta: float = 8.0) -> np.ndarray:
    x = high_raw_x(q)
    x_tilde = np.log(1.0 + kappa*x) / kappa
    return x_tilde / (1.0 + zeta*x_tilde)

def hubble_proxy(lambda_L_values, target_ratio: float = 1.0837):
    lam = np.array(lambda_L_values, dtype=float)
    eta = (target_ratio - 1.0) / (np.mean(lam) + 1e-12)
    return float(eta), 1.0 + eta*lam

def sigma8_proxy(TL_mean, TL_std):
    return 1.0 - TL_std/(TL_mean + 1e-12)

def theta_star_proxy(T_H_mean, T_H_std, T_L_mean, T_L_std, a=0.14, b=0.18):
    return float(100.0 * (T_H_std/(T_H_mean+1e-12)) * (T_H_mean/(T_L_mean+1e-12))**(-a) * (T_L_std/(T_L_mean+1e-12))**(-b))

def omega_m_proxy(T_H_mean, T_H_std, T_L_mean, T_L_std, R_H_mean):
    return float(T_H_std/(T_H_mean + 1e-12) + T_L_std/(T_L_mean + 1e-12) + (R_H_mean - 1.0))

def base_cosmo_summary(zeros, gamma=8.7, kappa=8.0, zeta=8.0, bias_scale=0.035, T_star=380000.0/137.0, Ns=(8,12,16), offsets=(0,1,2,3,4)):
    u = unfolded_gaps(zeros)
    raw_rows = []
    for N in Ns:
        for offset in offsets:
            q = build_blocks(u, N=N, offset=offset)
            bias = toy_local_bias(N, offset, scale=bias_scale)
            Sbar = float(np.mean(low_observable_base(q, bias)))
            lam_L, TL = low_temperature_from_Sbar(Sbar, gamma, T_star)
            TH_raw = T_star * high_scaled_y(q, kappa=kappa, zeta=zeta)
            raw_rows.append({
                "N": N, "offset": offset, "q": q, "bias": bias,
                "lambda_L": lam_L, "TL": TL, "TH_raw_mean": float(np.mean(TH_raw))
            })
    high_scale = 3000.0 / (np.mean([r["TH_raw_mean"] for r in raw_rows]) + 1e-12)
    TH_vals = np.array([r["TH_raw_mean"] * high_scale for r in raw_rows], dtype=float)
    TL_vals = np.array([r["TL"] for r in raw_rows], dtype=float)
    lam_vals = np.array([r["lambda_L"] for r in raw_rows], dtype=float)
    eta, RH_vals = hubble_proxy(lam_vals)

    TH_mean = float(np.mean(TH_vals)); TH_std = float(np.std(TH_vals))
    TL_mean = float(np.mean(TL_vals)); TL_std = float(np.std(TL_vals))
    RH_mean = float(np.mean(RH_vals)); RH_std = float(np.std(RH_vals))
    sigma8 = sigma8_proxy(TL_mean, TL_std)
    theta = theta_star_proxy(TH_mean, TH_std, TL_mean, TL_std)
    omega_m = omega_m_proxy(TH_mean, TH_std, TL_mean, TL_std, RH_mean)

    return {
        "rows": raw_rows,
        "T_H_mean": TH_mean, "T_H_std": TH_std,
        "T_L_mean": TL_mean, "T_L_std": TL_std,
        "R_H_mean": RH_mean, "R_H_std": RH_std,
        "sigma8_BM": sigma8, "theta_star_BM": theta,
        "omega_m_BM": omega_m,
    }

def split_only_eabc(rows, omega_m_bm, lam_dark=1.0, lam_bary=1.0):
    omega_b_list, omega_c_list = [], []
    for r in rows:
        q, bias = r["q"], r["bias"]
        CV = block_cv(q); F = block_norm_drift(q)
        E, A, C, G = eabc_extended_terms(q)
        bary = np.exp(-lam_bary*(CV + F + E + 0.5*C))
        dark = 1.0 - np.exp(-lam_dark*(C + G + A + 0.5*bias))
        bbar = np.mean(bary); cbar = np.mean(dark)
        omega_b = omega_m_bm * bbar / (bbar + cbar + 1e-12)
        omega_c = omega_m_bm * cbar / (bbar + cbar + 1e-12)
        omega_b_list.append(omega_b); omega_c_list.append(omega_c)
    omb = float(np.mean(omega_b_list)); omc = float(np.mean(omega_c_list))
    return {
        "omega_b_BM_mean": omb,
        "omega_c_BM_mean": omc,
        "omega_b_BM_std": float(np.std(omega_b_list)),
        "omega_c_BM_std": float(np.std(omega_c_list)),
        "baryon_fraction": float(omb/(omb+omc+1e-12)),
    }
