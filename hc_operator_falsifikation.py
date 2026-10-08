"""
hc_operator_falsifikation.py
Falsifikations-Suite: Cov(ψ) vs. L_C4+V vs. D+V.

Governance (#Energiedoku):
  Modell-Definition ≠ physikalische Interpretation.
  λ_min und Trennschärfe zum Nullmodell sind dimensionslose
  Fluktuationskennzahlen — keine Naturkonstanten.
  H_c / V: Modellannahmen; hier nur spektrale Robustheit vs. Nullmodell.

Bezug: Bamberg_Modell_Kern.lean (C4, α_B), Tag 2026-07-17-STABLE-SPECTRAL.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

BASE = Path("/Users/thomashoffbauer/Desktop/Mathematische Texte")
PRIMES_1B = BASE / "primes_cache_1e9.npy"
OUT_JSON = BASE / "hc_operator_falsifikation_report.json"

RNG = np.random.default_rng(20260717)


# --- Quaternion / GNS D (wie test_spektral_kern) ------------------------------

def quat_mul(q1, q2):
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    return (
        a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
        a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
        a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
        a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2,
    )


def quat_inverse(q):
    n2 = sum(x * x for x in q)
    c = (q[0], -q[1], -q[2], -q[3])
    return tuple(x / n2 for x in c)


def transition_operator(Q1, Q2):
    inv = quat_inverse(Q1)
    delta = quat_mul(inv, Q2)
    norm = math.sqrt(sum(x * x for x in delta))
    if norm == 0:
        return np.zeros((4, 4), dtype=float)
    vec = np.asarray(delta, dtype=float) / norm
    return np.outer(vec, vec)


def eabc_quaternion(p: int):
    r = p % 12
    if r == 5:
        return (1.0, 0.2 * p, 0.5 * p, 0.3 * p)
    if r == 11:
        return (1.0, 0.4 * p, 0.4 * p, 0.2 * p)
    return (1.0, float(p), 0.0, 0.0)


def build_D(quad_starts: np.ndarray, n_trans: int = 400) -> np.ndarray:
    quads = [eabc_quaternion(int(p)) for p in quad_starts[: n_trans + 1]]
    D = np.zeros((4, 4), dtype=float)
    for i in range(min(n_trans, len(quads) - 1)):
        B = transition_operator(quads[i], quads[i + 1])
        D += B.T @ B
    return 0.5 * (D + D.T)


# --- Strukturmatrizen ---------------------------------------------------------

def L_C4() -> np.ndarray:
    """Graph-Laplacian des C4-Orbit-Zyklus 0→1→2→3→0."""
    A = np.zeros((4, 4), dtype=float)
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    for i, j in edges:
        A[i, j] = A[j, i] = 1.0
    deg = np.diag(A.sum(axis=1))
    return deg - A


def V_bamberg(coupling: float = 1.0) -> np.ndarray:
    """
    Explizite Bamberg-Kopplung 5↔11 auf den Vierlings-Slots.
    Slot-Indexierung: 0:p, 1:p+2, 2:p+6, 3:p+8.
    Für p≡5: Slot0≡5, Slot2≡11 — komplementäre Injektion.
    Für p≡11: Slot0≡11, Slot2≡5 — dieselbe Kante {0,2}.
    """
    V = np.zeros((4, 4), dtype=float)
    V[0, 2] = V[2, 0] = -coupling  # anziehende / koppelnde Off-Diagonale
    # Spurfreiheit der Störung
    V -= np.diag(np.full(4, np.trace(V) / 4.0))
    return V


def V_null_permuted(coupling: float = 1.0, rng: np.random.Generator = RNG) -> np.ndarray:
    """Nullmodell: gleiche Stärke, aber falsche Slot-Paarung (kein {0,2})."""
    V = np.zeros((4, 4), dtype=float)
    # typische Fehlpaarung z.B. {1,3} (p+2 ↔ p+8), zerstört 5↔11-Injektion
    pairs = [(0, 1), (1, 3), (0, 3), (1, 2)]
    i, j = pairs[int(rng.integers(0, len(pairs)))]
    V[i, j] = V[j, i] = -coupling
    V -= np.diag(np.full(4, np.trace(V) / 4.0))
    return V


# --- Daten: ψ aus echten Vierlingen -------------------------------------------

def load_quad_starts(path: Path = PRIMES_1B) -> np.ndarray:
    primes = np.asarray(np.load(path, mmap_mode="r"), dtype=np.int64)

    def has(x):
        idx = np.searchsorted(primes, x)
        ok = idx < len(primes)
        idx = np.minimum(idx, len(primes) - 1)
        return ok & (primes[idx] == x)

    p = primes
    mask = (p > 3) & has(p + 2) & has(p + 6) & has(p + 8)
    return p[mask]


def build_psi(quad_starts: np.ndarray) -> np.ndarray:
    """
    Residualfeld ψ: normierte Slot-Phasen relativ zum Bamberg-Orbit.
    Φ: 5→0, 11→0.5, sonst 1; Residual = Φ - Orbitmittel.
    """
    offs = np.array([0, 2, 6, 8])
    members = quad_starts[:, None] + offs[None, :]
    r = members % 12
    phi = np.full(r.shape, 1.0, dtype=float)
    phi[r == 5] = 0.0
    phi[r == 11] = 0.5
    mean = phi.mean(axis=0, keepdims=True)
    psi = phi - mean
    # Normierung pro Zeile (dimensionslose Residualgröße)
    norms = np.linalg.norm(psi, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-12)
    return psi / norms


# --- Spektral-Diagnostik ------------------------------------------------------

def lambda_min(H: np.ndarray) -> float:
    H = 0.5 * (H + H.T)
    ev = np.linalg.eigvalsh(H)
    return float(np.min(ev))


def bootstrap_lambda_min_on_D(
    quad_starts: np.ndarray,
    V: np.ndarray,
    n_iterations: int = 400,
    n_trans: int = 400,
    rng: np.random.Generator = RNG,
) -> np.ndarray:
    """Bootstrap: resample Vierlingsstarts → D* → λ_min(D*+V)."""
    out = []
    n = len(quad_starts)
    for _ in range(n_iterations):
        idx = rng.choice(n, size=n, replace=True)
        D = build_D(quad_starts[idx], n_trans=min(n_trans, n - 1))
        out.append(lambda_min(D + V))
    return np.asarray(out, dtype=float)


def null_distribution_V(
    base: np.ndarray,
    coupling: float,
    n_iterations: int = 400,
    rng: np.random.Generator = RNG,
) -> np.ndarray:
    """Null: base + V_null (falsche Kanten) — zerstört 5↔11-Vorschrift."""
    return np.asarray(
        [lambda_min(base + V_null_permuted(coupling, rng)) for _ in range(n_iterations)],
        dtype=float,
    )


def bamberg_mode() -> np.ndarray:
    """Normierter Mode der komplementären Injektion Slot0↔Slot2 (5↔11)."""
    u = np.array([1.0, 0.0, -1.0, 0.0], dtype=float)
    return u / np.linalg.norm(u)


def rayleigh(H: np.ndarray, u: np.ndarray) -> float:
    H = 0.5 * (H + H.T)
    return float(u @ H @ u)


def separation_z(obs: float, null_samples: np.ndarray) -> float:
    mu = float(np.mean(null_samples))
    sd = float(np.std(null_samples))
    if sd < 1e-15:
        return float("inf") if abs(obs - mu) > 1e-12 else 0.0
    return (obs - mu) / sd


def _jsonify(x):
    if isinstance(x, dict):
        return {k: _jsonify(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonify(v) for v in x]
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        v = float(x)
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    return x


def analyze_candidate(
    name: str,
    H: np.ndarray,
    null_H_builder,
    n_null: int = 500,
    rng: np.random.Generator = RNG,
) -> dict:
    """
    Primärstatistik: Rayleigh-Quotient der Bamberg-Mode u=(1,0,-1,0)/√2.
    Null: gleiche Basis, aber V mit falscher Kantenpaarung.
    """
    u = bamberg_mode()
    H = 0.5 * (H + H.T)
    r_obs = rayleigh(H, u)
    lm = lambda_min(H)
    null_r = np.asarray([rayleigh(null_H_builder(rng), u) for _ in range(n_null)], dtype=float)
    null_lm = np.asarray([lambda_min(null_H_builder(rng)) for _ in range(n_null)], dtype=float)
    z_r = separation_z(r_obs, null_r)
    z_lm = separation_z(lm, null_lm)
    return {
        "candidate": name,
        "lambda_min": lm,
        "rayleigh_bamberg_mode": r_obs,
        "null_rayleigh_mean": float(np.mean(null_r)),
        "null_rayleigh_std": float(np.std(null_r)),
        "null_lambda_min_mean": float(np.mean(null_lm)),
        "null_lambda_min_std": float(np.std(null_lm)),
        "separation_z_rayleigh": z_r,
        "separation_z_lambda_min": z_lm,
        "distinguishable": abs(z_r) >= 3.0,
        "eigenvalues": np.linalg.eigvalsh(H).tolist(),
    }


def main():
    print("H_c-Falsifikations-Suite: Cov(ψ) | L_C4+V | D+V")
    coupling = 1.0
    L = L_C4()
    V = V_bamberg(coupling)

    # --- synthetisches Cov-Null (Reproduktion der Falsifikation) ---
    N = 28388
    psi_synth = RNG.normal(0, 0.01, (N, 4))
    psi_synth[:, 1] += psi_synth[:, 2] * 0.5
    H_cov = np.cov(psi_synth.T)
    # Null: Zeilenpermutation (wie zuvor — zerstört kaum Cov-Struktur)
    psi_shuf = psi_synth.copy()
    RNG.shuffle(psi_shuf)
    cov_null = np.min(np.real(np.linalg.eigvals(np.cov(psi_shuf.T))))
    cov_obs = lambda_min(H_cov)

    print("\n=== Kandidat Cov(ψ) [Replikation] ===")
    print(f"  λ_min original: {cov_obs:.6f}")
    print(f"  λ_min null:     {cov_null:.6f}")
    print(f"  unterscheidbar: {abs(cov_obs - cov_null) > 1e-6}")

    # --- echte Vierlinge ---
    quads = load_quad_starts()
    print(f"\n  Vierlinge geladen: n={len(quads):,}")
    psi = build_psi(quads)
    D = build_D(quads, n_trans=min(800, len(quads) - 1))

    H_LC4V = L + V
    H_DV = D + V

    boot_DV = bootstrap_lambda_min_on_D(quads, V, n_iterations=300, n_trans=400)

    res_cov = {
        "candidate": "Cov(psi)_synth",
        "lambda_min": cov_obs,
        "null_lambda_min": cov_null,
        "rayleigh_note": "Cov allein ohne explizite Strukturmatrix",
        "distinguishable": abs(cov_obs - cov_null) > 1e-6,
    }
    res_L = analyze_candidate(
        "L_C4 + V",
        H_LC4V,
        lambda rng: L + V_null_permuted(coupling, rng),
    )
    res_D = analyze_candidate(
        "D + V",
        H_DV,
        lambda rng: D + V_null_permuted(coupling, rng),
    )
    res_D["bootstrap_lambda_min_mean"] = float(np.mean(boot_DV))
    res_D["bootstrap_lambda_min_std"] = float(np.std(boot_DV))

    # Referenz: V=0 (Struktur aus)
    u = bamberg_mode()
    res_L["rayleigh_without_V"] = rayleigh(L, u)
    res_D["rayleigh_without_V"] = rayleigh(D, u)

    print("\n=== Kandidat L_C4 + V ===")
    print(f"  λ_min:              {res_L['lambda_min']:.6f}")
    print(f"  Rayleigh Bamberg:   {res_L['rayleigh_bamberg_mode']:.6f}")
    print(f"  Rayleigh ohne V:    {res_L['rayleigh_without_V']:.6f}")
    print(f"  Null Rayleigh μ/σ:  {res_L['null_rayleigh_mean']:.6f} / {res_L['null_rayleigh_std']:.6f}")
    print(f"  separation z (R):   {res_L['separation_z_rayleigh']:.3f}")
    print(f"  separation z (λ):   {res_L['separation_z_lambda_min']:.3f}")
    print(f"  unterscheidbar:     {res_L['distinguishable']}")

    print("\n=== Kandidat D + V ===")
    print(f"  λ_min:              {res_D['lambda_min']:.6f}")
    print(f"  Rayleigh Bamberg:   {res_D['rayleigh_bamberg_mode']:.6f}")
    print(f"  Rayleigh ohne V:    {res_D['rayleigh_without_V']:.6f}")
    print(f"  Null Rayleigh μ/σ:  {res_D['null_rayleigh_mean']:.6f} / {res_D['null_rayleigh_std']:.6f}")
    print(f"  separation z (R):   {res_D['separation_z_rayleigh']:.3f}")
    print(f"  Bootstrap λ μ/σ:    {res_D['bootstrap_lambda_min_mean']:.6f} / {res_D['bootstrap_lambda_min_std']:.6f}")
    print(f"  unterscheidbar:     {res_D['distinguishable']}")

    report = {
        "governance": (
            "dimensionslose Fluktuationskennzahlen; H_c/V Modellannahmen; "
            "keine Naturkonstanten"
        ),
        "primary_statistic": "Rayleigh-Quotient der Bamberg-Mode u=(1,0,-1,0)/sqrt(2)",
        "n_quadruplets": int(len(quads)),
        "coupling": coupling,
        "L_C4": L.tolist(),
        "V_bamberg": V.tolist(),
        "candidates": {
            "Cov_psi": res_cov,
            "L_C4_plus_V": res_L,
            "D_plus_V": res_D,
        },
        "verdict": {
            "Cov_psi_falsified": not res_cov["distinguishable"],
            "L_C4_plus_V_structure_signal": res_L["distinguishable"],
            "D_plus_V_structure_signal": res_D["distinguishable"],
        },
    }
    OUT_JSON.write_text(json.dumps(_jsonify(report), indent=2), encoding="utf-8")
    print(f"\nReport: {OUT_JSON}")

    v = report["verdict"]
    print("\n=== Verdict ===")
    print(f"  Cov(ψ) falsifiziert (kein Signal):     {v['Cov_psi_falsified']}")
    print(f"  L_C4+V vom Nullmodell trennbar:        {v['L_C4_plus_V_structure_signal']}")
    print(f"  D+V vom Nullmodell trennbar:           {v['D_plus_V_structure_signal']}")


if __name__ == "__main__":
    main()
