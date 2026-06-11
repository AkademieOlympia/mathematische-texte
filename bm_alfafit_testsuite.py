"""
BM / AlfaFit testsuite (mesh-free core tests + hooks for DEC spectrum)
Generated: 2026-01-25T12:26:01.814154Z

Assumes AQFT.json format:
- target_alpha_inverse (float)
- parameters: {N, alpha_infinity, c, N_total_transitions, N_coherent_transitions}
"""

from __future__ import annotations
import json, math

# Exact SI defining constants (2019 revision)
H = 6.62607015e-34          # J s (exact)
E_CHARGE = 1.602176634e-19  # C (exact)
C = 299792458.0             # m/s (exact)

def alpha_inv_pred(N: int, alpha_infinity: float, c_param: float) -> float:
    return alpha_infinity + c_param / math.sqrt(N)

def alpha_from_inv(alpha_inv: float) -> float:
    return 1.0 / alpha_inv

def eps0_from_alpha(alpha: float) -> float:
    # alpha = e^2 / (2 eps0 h c)  => eps0 = e^2 / (2 alpha h c)
    return (E_CHARGE**2) / (2.0 * alpha * H * C)

def mu0_from_eps0(eps0: float) -> float:
    return 1.0 / (eps0 * C**2)

def d_alpha_inv_dN(N: int, c_param: float) -> float:
    return -c_param / (2.0 * (N ** 1.5))

def run_mesh_free_tests(aqft: dict) -> dict:
    target_inv = float(aqft["target_alpha_inverse"])
    p = aqft["parameters"]
    N = int(p["N"])
    alpha_inf = float(p["alpha_infinity"])
    c_param = float(p["c"])

    pred = alpha_inv_pred(N, alpha_inf, c_param)
    delta = pred - target_inv

    alpha_target = alpha_from_inv(target_inv)
    eps0 = eps0_from_alpha(alpha_target)
    mu0 = mu0_from_eps0(eps0)

    out = {
        "N": N,
        "alpha_infinity": alpha_inf,
        "c": c_param,
        "alpha_inv_pred": pred,
        "alpha_inv_target": target_inv,
        "delta_alpha_inv": delta,
        "finite_term_c_over_sqrtN": c_param / math.sqrt(N),
        "d_alpha_inv_dN": d_alpha_inv_dN(N, c_param),
        "alpha_target": alpha_target,
        "eps0_from_alpha": eps0,
        "mu0_from_alpha": mu0,
        "c_required_for_exact_fit": (target_inv - alpha_inf) * math.sqrt(N),
        "alpha_infinity_required_for_exact_fit": target_inv - c_param / math.sqrt(N),
    }

    N_total = p.get("N_total_transitions")
    N_coh = p.get("N_coherent_transitions")
    if N_total and N_coh is not None:
        f_coh = float(N_coh) / float(N_total)
        out["coherence_fraction"] = f_coh
        out["coh_minus_alpha"] = f_coh - alpha_target
        out["N_coherent_ideal_float"] = alpha_target * float(N_total)
        out["N_coherent_nearest_int"] = int(round(alpha_target * float(N_total)))
    return out

if __name__ == "__main__":
    with open("AQFT.json","r",encoding="utf-8") as f:
        aqft = json.load(f)
    res = run_mesh_free_tests(aqft)
    print(json.dumps(res, indent=2))
