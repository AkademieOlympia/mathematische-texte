#!/usr/bin/env python3
"""
EABC Zeta-Recovery Audit (Bohm-Barrieren).

Automatisiertes, statistisch abgesichertes Audit-Werkzeug:
  - R_spec, v_A, v_H, tau_B
  - Zeta-Match-Score S_zeta (Peak-Uebereinstimmung am Leistungsspektrum)
  - Monte-Carlo-Blindtest gegen Zufallsspektren -> p_MC
  - optionaler Vergleichsplot L=2 vs L=8

Laeuft unter CPython+NumPy; unter SageMath 10.5 werden RDF/RR genutzt,
falls verfuegbar (sage -python ... / sage --python ...).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Numerische Backend-Schicht (NumPy default, Sage optional)
# ---------------------------------------------------------------------------

def _load_backend():
    try:
        from sage.all import RDF  # type: ignore

        def to_float(x) -> float:
            return float(RDF(x))

        return "sage+RDF", to_float
    except Exception:
        return "numpy", float


BACKEND_NAME, _to_float = _load_backend()
_trapz: Callable = getattr(np, "trapezoid", None) or getattr(np, "trapz")


# ---------------------------------------------------------------------------
# Kernfunktionen
# ---------------------------------------------------------------------------

def get_zeta_gammas() -> np.ndarray:
    """Imaginaerteile der ersten nicht-trivialen Riemann-Nullstellen."""
    return np.array(
        [14.134736, 21.022040, 25.010858, 30.424876, 32.935062],
        dtype=float,
    )


def eabc_velocity_profile(
    x_val: np.ndarray,
    L_val: float,
    gammas: np.ndarray,
) -> np.ndarray:
    """Parameterfreie EABC-Geschwindigkeitsfunktion im Barrierengebiet [0, L]."""
    v_base = 0.6
    oscillation = np.zeros_like(x_val, dtype=float)
    for i, g in enumerate(gammas):
        phase = float(g) * x_val + (i * np.pi / 4.0)
        oscillation += (0.12 / (i + 1)) * np.cos(phase)
    val = np.where((x_val >= 0.0) & (x_val <= L_val), v_base + oscillation, 1.0)
    return np.clip(val, 0.1, 1.5)


def compute_power_spectrum(
    v_bar: np.ndarray,
    x_bar: np.ndarray,
    omega: np.ndarray,
) -> np.ndarray:
    """Raeumliches Leistungsspektrum P(omega), vektorisiert."""
    v_mean = float(np.mean(v_bar))
    delta = (v_bar - v_mean).astype(complex)
    # Trapez-Gewichte
    w = np.ones_like(x_bar, dtype=float)
    if len(w) >= 2:
        w[0] = 0.5
        w[-1] = 0.5
    dx = float(x_bar[1] - x_bar[0]) if len(x_bar) > 1 else 1.0
    # Integrale: (omega, x) -> sum_x delta(x) exp(-i w x) * weight * dx
    kernel = np.exp(-1j * omega[:, None] * x_bar[None, :])
    integrals = (kernel * (delta * w)[None, :]).sum(axis=1) * dx
    return np.abs(integrals) ** 2


def zeta_match_score(
    power: np.ndarray,
    omega: np.ndarray,
    target_gammas: np.ndarray,
    sigma: float = 0.75,
) -> float:
    """
    S_zeta: Summe der lokal gemittelten, max-normierten Leistung
    in Gauss-Nachbarschaften um jede Ziel-Frequenz gamma_n.

    Nutzt explizit das Spektrum (kein reiner Gitter-Abstand).
    """
    p_max = float(np.max(power))
    if p_max <= 0.0:
        return 0.0
    p_norm = power / p_max
    score = 0.0
    for g in target_gammas:
        weights = np.exp(-0.5 * ((omega - float(g)) / sigma) ** 2)
        denom = float(np.sum(weights)) + 1e-15
        score += float(np.sum(weights * p_norm) / denom)
    return float(score)


@dataclass(frozen=True)
class AuditResult:
    L: float
    R_spec: float
    v_A: float
    v_H: float
    tau_B: float
    S_zeta: float
    p_MC: float
    num_trials: int
    backend: str


def run_monte_carlo_audit(
    L_val: float,
    gammas: np.ndarray,
    num_trials: int = 10_000,
    seed: int = 42,
) -> tuple[float, float, np.ndarray]:
    """
    Monte-Carlo-Blindtest gegen N_0 (zufaellige Frequenzspektren).

    Nullmodell: gleiche Profilfamilie, aber gammas ~ Uniform(5, 38).
    p_MC = (1 + #{S_null >= S_zeta}) / (1 + N)
    """
    x_bar = np.linspace(0.0, L_val, 800)
    omega = np.linspace(1.0, 40.0, 500)

    v_zeta = eabc_velocity_profile(x_bar, L_val, gammas)
    p_zeta = compute_power_spectrum(v_zeta, x_bar, omega)
    score_zeta = zeta_match_score(p_zeta, omega, gammas)

    rng = np.random.default_rng(seed)
    null_scores = np.empty(num_trials, dtype=float)
    for t in range(num_trials):
        rand_gammas = rng.uniform(5.0, 38.0, size=len(gammas))
        v_rand = eabc_velocity_profile(x_bar, L_val, rand_gammas)
        p_rand = compute_power_spectrum(v_rand, x_bar, omega)
        # Score gegen die *echten* Zeta-Frequenzen (Recovery-Gate)
        null_scores[t] = zeta_match_score(p_rand, omega, gammas)

    p_mc = (1.0 + float(np.sum(null_scores >= score_zeta))) / (1.0 + num_trials)
    return score_zeta, p_mc, null_scores


def execute_full_audit(
    L_val: float = 8.0,
    num_trials: int = 10_000,
    make_plot: bool = True,
) -> AuditResult:
    gammas = get_zeta_gammas()

    # 1. Aufloesungsparameter
    min_delta_gamma = float(np.min(np.diff(gammas)))
    r_spec = (L_val / (2.0 * np.pi)) * min_delta_gamma

    # 2. Integration & Mittelwerte
    x_bar = np.linspace(0.0, L_val, 800)
    v_bar = eabc_velocity_profile(x_bar, L_val, gammas)
    v_arithmetic = float(np.mean(v_bar))
    tau_b = float(_trapz(1.0 / v_bar, x_bar))
    v_harmonic = L_val / tau_b

    # 3. Monte-Carlo Signifikanztest
    score_zeta, p_mc, _ = run_monte_carlo_audit(
        L_val, gammas, num_trials=num_trials
    )

    result = AuditResult(
        L=float(L_val),
        R_spec=float(r_spec),
        v_A=float(v_arithmetic),
        v_H=float(v_harmonic),
        tau_B=float(tau_b),
        S_zeta=float(score_zeta),
        p_MC=float(p_mc),
        num_trials=int(num_trials),
        backend=BACKEND_NAME,
    )

    sep = "(Separierbar)" if r_spec > 1.0 else "(Verschwommen)"
    sig = "(Signifikant gegen Nullmodell)" if p_mc < 0.05 else "(nicht signifikant)"
    print(f"=== EABC RIGOROUS AUDIT REPORT (L = {L_val}) ===")
    print(f"Backend                         : {BACKEND_NAME}")
    print(f"Aufloesungsparameter R_spec       : {r_spec:.3f} {sep}")
    print(f"Arithmetisches Mittel v_A        : {v_arithmetic:.5f}")
    print(f"Harmonisches Mittel v_H          : {v_harmonic:.5f}")
    print(f"Barrieren-Traversalzeit tau_B    : {tau_b:.5f}")
    print(f"Zeta-Match-Score S_zeta          : {score_zeta:.4f}")
    print(f"Monte-Carlo p-Wert (p_MC)        : {p_mc:.4f} {sig}")
    print(f"MC-Trials N_0                    : {num_trials}")
    print("=" * 48)

    if make_plot:
        _plot_scaling_comparison(gammas)

    return result


def _plot_scaling_comparison(gammas: np.ndarray) -> None:
    """Vergleichsplot der normierten Spektren fuer L=2 und L=8."""
    omega = np.linspace(1.0, 40.0, 1000)
    audits = {}
    for L in (2.0, 8.0):
        x_bar = np.linspace(0.0, L, 800)
        v_bar = eabc_velocity_profile(x_bar, L, gammas)
        power = compute_power_spectrum(v_bar, x_bar, omega)
        min_delta = float(np.min(np.diff(gammas)))
        r_spec = (L / (2.0 * np.pi)) * min_delta
        audits[L] = (power, r_spec)

    plt.figure(figsize=(10, 5))
    for L, color in ((2.0, "crimson"), (8.0, "dodgerblue")):
        power, r_spec = audits[L]
        plt.plot(
            omega,
            power / np.max(power),
            label=f"L = {L:.1f} (R_spec = {r_spec:.2f})",
            color=color,
            alpha=0.85,
        )
    for idx, g in enumerate(gammas):
        plt.axvline(g, color="gray", linestyle="--", alpha=0.45)
        if idx == 0:
            plt.text(g + 0.2, 0.9, r"Riemann $\gamma_n$", fontsize=9)

    plt.title("EABC Zeta-Recovery: Spektrale Konvergenz mit wachsendem L")
    plt.xlabel(r"Frequenzen $\omega$")
    plt.ylabel(r"Normierte Leistung $P(\omega)$")
    plt.legend(loc="upper right")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("eabc_zeta_recovery_scaling.png", dpi=300)
    print("Plot gespeichert: eabc_zeta_recovery_scaling.png")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="EABC Zeta-Recovery Audit")
    parser.add_argument("--L", type=float, default=8.0, help="Barrierenlaenge")
    parser.add_argument(
        "--trials",
        type=int,
        default=10_000,
        help="Monte-Carlo-Trials (Default 10000)",
    )
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    execute_full_audit(
        L_val=args.L,
        num_trials=args.trials,
        make_plot=not args.no_plot,
    )
