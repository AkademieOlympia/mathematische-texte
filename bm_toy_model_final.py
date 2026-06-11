
import numpy as np

# Finale konsistente Toy-Version mit BM-Observablen:
# T_H, T_L, R_H, sigma8_BM, theta*_BM, Omega_m^BM, n_s^BM

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

def low_observable(q: np.ndarray, bias: np.ndarray) -> np.ndarray:
    cv = block_cv(q)
    F = block_norm_drift(q)
    return 2.0 * cv + 2.0 * bias + 1.0 * F

def low_temperature(q: np.ndarray, bias: np.ndarray, gamma: float, T_star: float):
    S = low_observable(q, bias)
    Sbar = float(np.mean(S))
    lam_L = np.exp(-gamma * Sbar)
    T_L = T_star * lam_L
    return lam_L, T_L, Sbar

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
        for k in range(len(W_T) - 1):
            W_O[k] = W_T[k] * np.conj(W_T[k + 1])
        W_O[-1] = W_O[-2]
    Phi_O = np.angle(W_O)
    return Phi_T, Phi_O

def storage_term(q: np.ndarray) -> np.ndarray:
    norms2 = np.sum(q * q, axis=1)
    gamma = np.zeros(len(q))
    for k in range(len(q)):
        left = np.linalg.norm(q[k] - q[k - 1]) if k > 0 else 0.0
        right = np.linalg.norm(q[k + 1] - q[k]) if k < len(q) - 1 else 0.0
        gamma[k] = 1e-6 + left + right
    return norms2 / gamma

def coherence_term(q: np.ndarray) -> np.ndarray:
    c = []
    for block in q:
        mean = np.mean(block)
        std = np.std(block)
        phases = np.array([
            1.0 * mean + 0.5 * std,
            2.0 * mean + 1.0 * std,
            3.0 * mean + 1.5 * std,
            4.0 * mean + 2.0 * std,
        ])
        amp = np.mean(np.exp(1j * phases))
        c.append(np.abs(amp) ** 2)
    return np.array(c)

def holonomy_term(q: np.ndarray) -> np.ndarray:
    Phi_T, Phi_O = holonomy_arrays(q)
    return (1.0 - np.cos(Phi_T)) + (1.0 - np.cos(Phi_O))

def high_raw_x(q: np.ndarray, alpha_s: float = 1.0, alpha_h: float = 0.3) -> np.ndarray:
    return alpha_s * storage_term(q) * coherence_term(q) + alpha_h * holonomy_term(q)

def high_scaled_y(q: np.ndarray, kappa: float = 8.0, zeta: float = 8.0) -> np.ndarray:
    x = high_raw_x(q)
    x_tilde = np.log(1.0 + kappa * x) / kappa
    return x_tilde / (1.0 + zeta * x_tilde)

def hubble_proxy(lambda_L_values, target_ratio: float = 1.0837):
    lam = np.array(lambda_L_values, dtype=float)
    eta = (target_ratio - 1.0) / (np.mean(lam) + 1e-12)
    R = 1.0 + eta * lam
    return float(eta), R

def sigma8_proxy(T_L_mean: float, T_L_std: float):
    return 1.0 - T_L_std / (T_L_mean + 1e-12)

def theta_star_proxy_new(T_H_mean: float, T_H_std: float, T_L_mean: float, T_L_std: float,
                         a: float = 0.14, b: float = 0.18) -> float:
    rel_high = T_H_std / (T_H_mean + 1e-12)
    ratio_hl = T_H_mean / (T_L_mean + 1e-12)
    rel_low = T_L_std / (T_L_mean + 1e-12)
    return float(100.0 * rel_high * ratio_hl**(-a) * rel_low**(-b))

def omega_m_proxy(T_H_mean: float, T_H_std: float, T_L_mean: float, T_L_std: float, R_H_mean: float) -> float:
    return float(T_H_std/(T_H_mean + 1e-12) + T_L_std/(T_L_mean + 1e-12) + (R_H_mean - 1.0))

def fit_beta(logN: np.ndarray, logr: np.ndarray) -> float:
    A = np.vstack([logN, np.ones_like(logN)]).T
    m, c = np.linalg.lstsq(A, logr, rcond=None)[0]
    return float(-m)

def ns_proxy_from_sector_tilts(N_vals, rL_vals, rH_vals) -> dict:
    logN = np.log(np.asarray(N_vals, dtype=float))
    beta_L = fit_beta(logN, np.log(np.asarray(rL_vals, dtype=float)))
    beta_H = fit_beta(logN, np.log(np.asarray(rH_vals, dtype=float)))
    delta_beta = abs(beta_L - beta_H)
    n_s = 1.0 - delta_beta / 3.0
    return {
        "beta_L": float(beta_L),
        "beta_H": float(beta_H),
        "delta_beta": float(delta_beta),
        "n_s_BM": float(n_s),
    }

def run_model(
    zeros: np.ndarray,
    Ns=(8, 12, 16),
    offsets=(0, 1, 2, 3, 4),
    gamma: float = 8.7,
    kappa: float = 8.0,
    zeta: float = 8.0,
    bias_scale: float = 0.035,
    T_star: float = 380000.0 / 137.0,
    target_high_mean: float = 3000.0,
    theta_a: float = 0.14,
    theta_b: float = 0.18,
):
    u = unfolded_gaps(zeros)
    raw_rows = []

    for N in Ns:
        for offset in offsets:
            q = build_blocks(u, N=N, offset=offset)
            bias = toy_local_bias(N, offset, scale=bias_scale)
            lam_L, T_L_raw, Sbar = low_temperature(q, bias=bias, gamma=gamma, T_star=T_star)
            y = high_scaled_y(q, kappa=kappa, zeta=zeta)
            T_H_raw_arr = T_star * y
            raw_rows.append({
                "N": N,
                "offset": offset,
                "lambda_L": float(lam_L),
                "T_L_raw": float(T_L_raw),
                "Sbar": float(Sbar),
                "T_H_raw_mean": float(np.mean(T_H_raw_arr)),
                "T_H_raw_std_local": float(np.std(T_H_raw_arr)),
            })

    global_high_mean = np.mean([r["T_H_raw_mean"] for r in raw_rows])
    high_scale = target_high_mean / (global_high_mean + 1e-12)

    rows = []
    lambda_L_list = []

    for r in raw_rows:
        T_H = r["T_H_raw_mean"] * high_scale
        T_H_std_local = r["T_H_raw_std_local"] * high_scale
        rows.append({
            "N": r["N"],
            "offset": r["offset"],
            "lambda_L": r["lambda_L"],
            "T_L": r["T_L_raw"],
            "T_H": T_H,
            "T_H_std_local": T_H_std_local,
            "Sbar": r["Sbar"],
        })
        lambda_L_list.append(r["lambda_L"])

    eta, R = hubble_proxy(lambda_L_list)
    for i, row in enumerate(rows):
        row["R_H"] = float(R[i])

    T_H_vals = np.array([r["T_H"] for r in rows], dtype=float)
    T_L_vals = np.array([r["T_L"] for r in rows], dtype=float)
    R_vals = np.array([r["R_H"] for r in rows], dtype=float)

    T_H_mean = float(np.mean(T_H_vals))
    T_H_std = float(np.std(T_H_vals))
    T_L_mean = float(np.mean(T_L_vals))
    T_L_std = float(np.std(T_L_vals))
    R_H_mean = float(np.mean(R_vals))
    R_H_std = float(np.std(R_vals))

    sigma8 = sigma8_proxy(T_L_mean, T_L_std)
    theta_star = theta_star_proxy_new(T_H_mean, T_H_std, T_L_mean, T_L_std, a=theta_a, b=theta_b)
    omega_m = omega_m_proxy(T_H_mean, T_H_std, T_L_mean, T_L_std, R_H_mean)

    # sektorweise Rauigkeiten nach N
    byN = {}
    for N in Ns:
        subset = [r for r in rows if r["N"] == N]
        THN = np.array([r["T_H"] for r in subset], dtype=float)
        TLN = np.array([r["T_L"] for r in subset], dtype=float)
        byN[N] = {
            "TH_mean": float(np.mean(THN)),
            "TH_std": float(np.std(THN)),
            "TL_mean": float(np.mean(TLN)),
            "TL_std": float(np.std(TLN)),
            "rH": float(np.std(THN)/(np.mean(THN) + 1e-12)),
            "rL": float(np.std(TLN)/(np.mean(TLN) + 1e-12)),
        }

    N_vals = list(Ns)
    rL_vals = [byN[N]["rL"] for N in N_vals]
    rH_vals = [byN[N]["rH"] for N in N_vals]
    ns_info = ns_proxy_from_sector_tilts(N_vals, rL_vals, rH_vals)

    return {
        "params": {
            "gamma": gamma,
            "kappa": kappa,
            "zeta": zeta,
            "bias_scale": bias_scale,
            "theta_a": theta_a,
            "theta_b": theta_b,
            "T_star": T_star,
            "low_observable": "S_L = 2*CV + 2*B + F",
            "theta_formula": "100*(sigma(TH)/mean(TH))*(TH/TL)^(-a)*(sigma(TL)/TL)^(-b)"
        },
        "high_scale": float(high_scale),
        "eta": float(eta),
        "T_H_mean": T_H_mean,
        "T_H_std": T_H_std,
        "T_L_mean": T_L_mean,
        "T_L_std": T_L_std,
        "R_H_mean": R_H_mean,
        "R_H_std": R_H_std,
        "theta_star_BM": float(theta_star),
        "sigma8_BM": float(sigma8),
        "omega_m_BM": float(omega_m),
        "n_s_BM": float(ns_info["n_s_BM"]),
        "beta_L": float(ns_info["beta_L"]),
        "beta_H": float(ns_info["beta_H"]),
        "delta_beta": float(ns_info["delta_beta"]),
        "byN": byN,
        "rows": rows,
    }

if __name__ == "__main__":
    from pathlib import Path
    p = Path("zeros6.npy")
    zeros = np.load(p)
    result = run_model(zeros=zeros)
    for k, v in result.items():
        if k not in {"rows", "byN"}:
            print(k, "=", v)
