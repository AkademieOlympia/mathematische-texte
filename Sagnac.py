#!/usr/bin/env python3
"""
Numerischer Test des arithmetischen Sagnac-Terms im Primvierlings-Modell.

Modellkern:
    H0      = diag(E) + g * K0
    H_sag   = diag(E) + g * K_sag + beta * diag(sigma * Phi)

mit
    K0[m,n]      = exp(-|J_m - J_n| / tau)
    K_sag[m,n]   = K0[m,n] * exp(i * (sigma_m*Phi_m - sigma_n*Phi_n))
    J_n          = p_n^2 * (1/sqrt(p_n)) = p_n^(3/2)
    Phi_n        = alpha * A_n * Omega_n
    Omega_n      = 1/sqrt(p_n)
    A_n          = 8  (effektive Zellflaeche der Vierlingszelle)

Die Ausgabe zeigt:
  - Anzahl gefundener Vierlinge
  - Spektralverschiebung gegenueber H0
  - Gap-Ratio <r> (Poisson ~0.386, GOE ~0.536, GUE ~0.603)
  - Chirale Aufspaltung via Vergleich alpha vs -alpha
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import Iterable, List

import numpy as np


@dataclass
class PrimeQuadruplet:
    p: int
    sigma: int  # +1 fuer ABCE, -1 fuer CEAB


def sieve_primes(limit: int) -> np.ndarray:
    """Gibt ein Bool-Array is_prime[0..limit] zurueck."""
    if limit < 2:
        return np.zeros(limit + 1, dtype=bool)
    is_prime = np.ones(limit + 1, dtype=bool)
    is_prime[:2] = False
    max_i = int(limit ** 0.5)
    for i in range(2, max_i + 1):
        if is_prime[i]:
            is_prime[i * i : limit + 1 : i] = False
    return is_prime


def find_prime_quadruplets(max_p: int, max_count: int | None = None) -> List[PrimeQuadruplet]:
    """Findet Vierlinge (p, p+2, p+6, p+8) mit p <= max_p."""
    is_prime = sieve_primes(max_p + 8)
    quads: List[PrimeQuadruplet] = []

    for p in range(5, max_p + 1, 2):
        if not (is_prime[p] and is_prime[p + 2] and is_prime[p + 6] and is_prime[p + 8]):
            continue

        r = p % 12
        if r == 5:  # ABCE
            sigma = +1
        elif r == 11:  # CEAB
            sigma = -1
        else:
            # Sollte fuer p > 3 bei Vierlingen nicht auftreten.
            continue

        quads.append(PrimeQuadruplet(p=p, sigma=sigma))
        if max_count is not None and len(quads) >= max_count:
            break

    return quads


def gap_ratio(evals: np.ndarray) -> float:
    """Mittelwert der benachbarten Gap-Ratios."""
    e = np.sort(np.real(evals))
    if e.size < 3:
        return float("nan")
    spacings = np.diff(e)
    s1 = spacings[:-1]
    s2 = spacings[1:]
    denom = np.maximum(s1, s2)
    ok = denom > 0
    if not np.any(ok):
        return float("nan")
    r = np.minimum(s1[ok], s2[ok]) / denom[ok]
    return float(np.mean(r))


def build_matrices(
    quads: List[PrimeQuadruplet],
    alpha: float,
    tau: float,
    coupling_g: float,
    beta: float,
    area: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Baut H0 und H_sag fuer eine gegebene Vierlingsliste."""
    p = np.array([q.p for q in quads], dtype=float)
    sigma = np.array([q.sigma for q in quads], dtype=float)

    omega = 1.0 / np.sqrt(p)
    j = p * p * omega
    phi = alpha * area * omega

    # Normierte Diagonalenergie.
    e_diag = (j - np.mean(j)) / (np.std(j) + 1e-12)
    h_diag = np.diag(e_diag)

    # Amplitudenkern.
    delta_j = np.abs(j[:, None] - j[None, :])
    k0 = np.exp(-delta_j / tau)
    np.fill_diagonal(k0, 0.0)

    # Sagnac-Phase.
    phase = np.exp(1j * (sigma[:, None] * phi[:, None] - sigma[None, :] * phi[None, :]))
    ksag = k0 * phase

    h0 = h_diag + coupling_g * k0
    hsag = h_diag + coupling_g * ksag + beta * np.diag(sigma * phi)

    return h0, hsag


def spectral_shift(e0: np.ndarray, e1: np.ndarray) -> float:
    """Mittlere absolute Verschiebung sortierter Eigenwerte."""
    a = np.sort(np.real(e0))
    b = np.sort(np.real(e1))
    n = min(len(a), len(b))
    if n == 0:
        return float("nan")
    return float(np.mean(np.abs(a[:n] - b[:n])))


def chiral_split_metric(
    quads: List[PrimeQuadruplet],
    alpha: float,
    tau: float,
    coupling_g: float,
    beta: float,
    area: float,
) -> float:
    """Vergleicht Spektren fuer +alpha und -alpha als Aufspaltungs-Metrik."""
    _, h_plus = build_matrices(quads, +alpha, tau, coupling_g, beta, area)
    _, h_minus = build_matrices(quads, -alpha, tau, coupling_g, beta, area)

    e_plus = np.sort(np.real(np.linalg.eigvalsh(h_plus)))
    e_minus = np.sort(np.real(np.linalg.eigvalsh(h_minus)))
    n = min(len(e_plus), len(e_minus))
    if n == 0:
        return float("nan")
    return float(np.mean(np.abs(e_plus[:n] - e_minus[:n])) / 2.0)


def parse_alpha_list(raw: str) -> List[float]:
    vals = []
    for part in raw.split(","):
        part = part.strip()
        if not part:
            continue
        vals.append(float(part))
    if not vals:
        raise ValueError("Keine gueltigen Alpha-Werte gefunden.")
    return vals


def run_experiment(
    max_p: int,
    max_quads: int | None,
    alpha_values: Iterable[float],
    tau: float,
    coupling_g: float,
    beta: float,
    area: float,
) -> None:
    quads = find_prime_quadruplets(max_p=max_p, max_count=max_quads)
    if len(quads) < 8:
        raise RuntimeError(
            f"Zu wenige Vierlinge gefunden ({len(quads)}). "
            "Bitte max-p erhoehen oder max-quads lockern."
        )

    sigma_vals = np.array([q.sigma for q in quads], dtype=int)
    frac_plus = np.mean(sigma_vals > 0)
    frac_minus = np.mean(sigma_vals < 0)

    print("=== Primvierlingsdatensatz ===")
    print(f"max_p                : {max_p}")
    print(f"anzahl vierlinge N   : {len(quads)}")
    print(f"chir. anteil +1      : {frac_plus:.3f}")
    print(f"chir. anteil -1      : {frac_minus:.3f}")
    print()

    # Referenz ohne Sagnac-Term
    h0, _ = build_matrices(quads, alpha=0.0, tau=tau, coupling_g=coupling_g, beta=beta, area=area)
    e0 = np.linalg.eigvalsh(h0)
    r0 = gap_ratio(e0)

    print("=== Referenz H0 (alpha=0) ===")
    print(f"<r>                  : {r0:.5f}")
    print("Referenzwerte: Poisson~0.386, GOE~0.536, GUE~0.603")
    print()

    print("=== Sweep ueber alpha ===")
    print(
        "alpha      <r>(H_sag)   shift_vs_H0   chiral_split(+/-alpha)   max|Im(ev)|"
    )
    print("-" * 78)

    for alpha in alpha_values:
        _, hsag = build_matrices(
            quads,
            alpha=alpha,
            tau=tau,
            coupling_g=coupling_g,
            beta=beta,
            area=area,
        )
        es = np.linalg.eigvals(hsag)
        # Numerisch sollte hsag hermitesch sein -> Im-Teile nahe 0.
        max_im = float(np.max(np.abs(np.imag(es))))
        es_real = np.real(es)

        rs = gap_ratio(es_real)
        shift = spectral_shift(e0, es_real)
        split = chiral_split_metric(
            quads,
            alpha=alpha,
            tau=tau,
            coupling_g=coupling_g,
            beta=beta,
            area=area,
        )

        print(
            f"{alpha:7.4f}    {rs:10.5f}   {shift:11.6f}        "
            f"{split:11.6f}        {max_im:9.2e}"
        )

    print()
    print("Interpretation (kurz):")
    print("- shift_vs_H0 > 0: Sagnac-Term veraendert das Spektrum.")
    print("- chiral_split > 0: orientierungsabhaengige Aufspaltung sichtbar.")
    print("- Aenderung von <r>: Hinweis auf geaenderte Spektralstatistik.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Testet den arithmetischen Sagnac-Term im Primvierlingsmodell."
    )
    parser.add_argument(
        "--max-p",
        type=int,
        default=300000,
        help="Suche Vierlinge bis p <= max-p (Default: 300000).",
    )
    parser.add_argument(
        "--max-quads",
        type=int,
        default=300,
        help="Maximale Anzahl Vierlinge fuer den Hamiltonian (Default: 300).",
    )
    parser.add_argument(
        "--alphas",
        type=str,
        default="0.00,0.02,0.05,0.10,0.20,0.35",
        help="Kommagetrennte Alpha-Werte (Default: 0.00,0.02,0.05,0.10,0.20,0.35).",
    )
    parser.add_argument(
        "--tau",
        type=float,
        default=140.0,
        help="Abklinglaenge des Kopplungskerns exp(-|dJ|/tau).",
    )
    parser.add_argument(
        "--g",
        type=float,
        default=0.35,
        help="Staerke der Off-Diagonal-Kopplung.",
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=1.0,
        help="Gewicht des diagonalen Sagnac-Terms beta*diag(sigma*Phi).",
    )
    parser.add_argument(
        "--area",
        type=float,
        default=8.0,
        help="Effektive Zellflaeche A_n (Default: 8.0).",
    )
    args = parser.parse_args()

    alphas = parse_alpha_list(args.alphas)
    max_quads = None if args.max_quads <= 0 else args.max_quads

    run_experiment(
        max_p=args.max_p,
        max_quads=max_quads,
        alpha_values=alphas,
        tau=args.tau,
        coupling_g=args.g,
        beta=args.beta,
        area=args.area,
    )


if __name__ == "__main__":
    main()
