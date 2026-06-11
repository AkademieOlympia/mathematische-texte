#!/usr/bin/env python3
"""
EABC Tschebyscheff-/Chebyshev-Bias: Dirichlet-Charakter-Projektion mod 420.

Kanäle (Restklassen): 11, 101, 191, 221, 311, 401
Besetzungszahlen aus Stiefel.py (N=100M, 4767 Vierlingskanäle).
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
from sympy import factorint, gcd, isprime, primitive_root

MOD = 420
CHANNELS = [11, 101, 191, 221, 311, 401]
COUNTS = [765, 831, 786, 809, 809, 767]
N = sum(COUNTS)


def channel_mod420(p: int) -> int | None:
    r = p % MOD
    return r if r in CHANNELS else None


def _discrete_log(a: int, gen: int, order: int, mod: int) -> int:
    """a = gen^k (mod mod), k in {0,...,order-1}."""
    x = a % mod
    gk = 1
    for k in range(order):
        if gk == x:
            return k
        gk = (gk * gen) % mod
    raise ValueError(f"{a} ist keine Einheit mod {mod}")


def _char_prime_power(a: int, q: int, order: int, gen: int, exp: int) -> complex:
    if gcd(a, q) != 1:
        return 0.0
    k = _discrete_log(a, gen, order, q)
    return np.exp(2j * np.pi * exp * k / order)


def build_dirichlet_characters_mod_420() -> list[tuple[str, int, callable]]:
    """
    Alle nichttrivialen Dirichlet-Charaktere mod 420 via CRT, gefiltert auf primitive.

    420 = 4 * 3 * 5 * 7,  (Z/420Z)* ≅ C2 × C2 × C4 × C6,  φ(420) = 96.
    """
    specs = {
        4: (2, 3),
        3: (2, 2),
        5: (4, 2),
        7: (6, 3),
    }
    proper_divisors = [
        d for d in range(2, MOD) if MOD % d == 0 and d < MOD
    ]

    characters: list[tuple[str, int, callable]] = []

    for e4 in range(2):
        for e3 in range(2):
            for e5 in range(4):
                for e7 in range(6):
                    if e4 == e3 == e5 == e7 == 0:
                        continue

                    def make_chi(e4=e4, e3=e3, e5=e5, e7=e7):
                        def chi(a: int) -> complex:
                            if gcd(a, MOD) != 1:
                                return 0.0
                            return (
                                _char_prime_power(a % 4, 4, *specs[4], e4)
                                * _char_prime_power(a % 3, 3, *specs[3], e3)
                                * _char_prime_power(a % 5, 5, *specs[5], e5)
                                * _char_prime_power(a % 7, 7, *specs[7], e7)
                            )

                        return chi

                    chi = make_chi()

                    # Ordnung: kleinstes m, sodass chi(a)^m = 1 für alle Einheiten a
                    order = 1
                    units = [a for a in range(1, MOD) if gcd(a, MOD) == 1]
                    for m in range(1, 97):
                        if all(abs(chi(a) ** m - 1.0) < 1e-9 for a in units):
                            order = m
                            break

                    # Primitivität: nicht induziert von echtem Teiler d | 420
                    induced = False
                    for d in proper_divisors:
                        reps: dict[int, complex] = {}
                        ok = True
                        for a in units:
                            key = a % d
                            val = chi(a)
                            if key in reps:
                                if abs(reps[key] - val) > 1e-9:
                                    ok = False
                                    break
                            else:
                                reps[key] = val
                        if ok:
                            induced = True
                            break

                    if not induced:
                        name = f"χ({e4},{e3},{e5},{e7})"
                        characters.append((name, order, chi))

    return characters


def project_delta_onto_characters(
    delta: np.ndarray, characters: list[tuple[str, int, callable]]
) -> list[dict]:
    """
    Projektion δ_i auf Charakter-Basis:
      c_χ = (1/√6) Σ_i δ_i · χ̅(r_i)
    wobei r_i die 6 Kanalrestklassen sind.
    """
    results = []
    for name, order, chi in characters:
        vals = np.array([chi(r) for r in CHANNELS], dtype=complex)
        # Orthogonale Projektion auf die 6-Kanal-Teilmenge
        coeff = np.dot(delta, np.conj(vals)) / math.sqrt(len(CHANNELS))
        results.append({
            "name": name,
            "order": order,
            "coeff": coeff,
            "magnitude": abs(coeff),
            "phase": np.angle(coeff),
            "chi_values": vals,
        })
    return sorted(results, key=lambda x: -x["magnitude"])


def ring_fourier(counts: list[int]) -> np.ndarray:
    v = np.array(counts, dtype=float)
    v /= v.sum()
    return np.fft.fft(v)


def order3_characters(characters: list[tuple[str, int, callable]]) -> list:
    return [(n, o, c) for n, o, c in characters if o == 3]


def main():
    print("=" * 70)
    print("EABC Tschebyscheff-/Chebyshev-Bias: Dirichlet-Projektion mod 420")
    print("=" * 70)

    p = np.array(COUNTS, dtype=float) / N
    delta = p - 1.0 / 6.0

    print(f"\nN = {N}")
    print(f"Kanäle: {CHANNELS}")
    print(f"counts: {COUNTS}")
    print(f"p:      {[round(x, 6) for x in p]}")
    print(f"δ:      {[round(x, 6) for x in delta]}")

    # B_{101,11}
    idx_11 = CHANNELS.index(11)
    idx_101 = CHANNELS.index(101)
    B_101_11 = COUNTS[idx_101] - COUNTS[idx_11]
    B_norm = B_101_11 / math.sqrt(N)
    print(f"\nB_{{101,11}} = N_101 - N_11 = {COUNTS[idx_101]} - {COUNTS[idx_11]} = {B_101_11}")
    print(f"B_{{101,11}} / √N = {B_norm:.6f}")

    # χ²-Test auf Gleichverteilung
    expected = N / 6.0
    chi2 = sum((c - expected) ** 2 / expected for c in COUNTS)
    # df=5, kritisch 11.07 bei α=0.05
    print(f"\nχ²-Test (Gleichverteilung, df=5): χ² = {chi2:.4f}  (kritisch 11.07 @ α=0.05)")

    # LIL-Grenze (heuristisch für √N-Skalierung)
    lil_bound = math.sqrt(2 * N * math.log(max(math.log(N), 3)))
    print(f"LIL-artige √N-Grenze (2N log log N)^½ ≈ {lil_bound:.2f}")
    print(f"|B_101,11| / LIL-Grenze = {abs(B_101_11) / lil_bound:.4f}")

    # Ring-Fourier (C6-Indexing)
    F = ring_fourier(COUNTS)
    print("\n--- Ring-Fourier (6 Kanäle als C₆) ---")
    for k in range(6):
        print(f"  k={k}: |F_k| = {abs(F[k]):.8f},  arg = {np.angle(F[k]):.6f}")

    k2_share = abs(F[2]) ** 2 + abs(F[4]) ** 2
    total_power = sum(abs(F[k]) ** 2 for k in range(1, 6))
    print(f"\n  k=2+k=4 Anteil an nicht-DC-Varianz: {100 * k2_share / total_power:.2f}%")

    # Dirichlet-Charaktere
    print("\n--- Primitive Dirichlet-Charaktere mod 420 ---")
    characters = build_dirichlet_characters_mod_420()
    print(f"  Anzahl primitiver Charaktere: {len(characters)}")

    # Charakterwerte auf den 6 Kanälen
    projections = project_delta_onto_characters(delta, characters)

    print("\n--- Top-10 Projektionen |c_χ| ---")
    for i, r in enumerate(projections[:10]):
        print(
            f"  {i+1}. {r['name']:20s}  ord={r['order']}  "
            f"|c|={r['magnitude']:.8f}  arg={r['phase']:.4f}"
        )

    strongest = projections[0]
    print(f"\nStärkster Charakter: {strongest['name']}, Ordnung {strongest['order']}")
    print(f"  |c| = {strongest['magnitude']:.8f}")
    print(f"  χ-Werte auf Kanälen: {[f'{v:.3f}' for v in strongest['chi_values']]}")

    # Ordnung-3 Charaktere
    ord3 = order3_characters(characters)
    print(f"\n--- Charaktere der Ordnung 3: {len(ord3)} ---")
    ord3_proj = project_delta_onto_characters(delta, ord3)
    for i, r in enumerate(ord3_proj[:5]):
        print(
            f"  {i+1}. {r['name']:20s}  |c|={r['magnitude']:.8f}  arg={r['phase']:.4f}"
        )

    # k=2 Fourier ↔ Dirichlet Ordnung 3 Verbindung
    print("\n--- k=2-Fourier ↔ Dirichlet Ordnung 3 ---")
    # k=2 FFT-Modus auf C6: exp(-2πi·2·j/6) = exp(-2πi·j/3)
    omega3 = np.exp(-2j * np.pi / 3)
    fft_k2_vec = np.array([omega3 ** j for j in range(6)])
    fft_k2_coeff = np.dot(delta, np.conj(fft_k2_vec)) / math.sqrt(6)
    print(f"  FFT k=2 Koeffizient (δ-Projektion): |c| = {abs(fft_k2_coeff):.8f}, arg = {np.angle(fft_k2_coeff):.4f}")
    print(f"  FFT k=2 aus counts:               |F_2| = {abs(F[2]):.8f}, arg = {np.angle(F[2]):.4f}")

    if ord3_proj:
        best_ord3 = ord3_proj[0]
        # Korrelation der χ-Vektoren mit FFT k=2 Richtung
        chi_vec = best_ord3["chi_values"]
        # Normalisiere
        chi_n = chi_vec / np.linalg.norm(chi_vec)
        fft_n = fft_k2_vec / np.linalg.norm(fft_k2_vec)
        overlap = abs(np.vdot(chi_n, fft_n))
        print(f"\n  Bester Ord-3-Charakter: {best_ord3['name']}")
        print(f"  |c_χ(δ)| = {best_ord3['magnitude']:.8f}")
        print(f"  Überlapp |⟨χ_ord3, FFT_k2⟩| = {overlap:.6f}")

        # Auch k=4 (konjugiert)
        omega3_conj = np.exp(2j * np.pi / 3)
        fft_k4_vec = np.array([omega3_conj ** j for j in range(6)])
        fft_k4_coeff = np.dot(delta, np.conj(fft_k4_vec)) / math.sqrt(6)
        print(f"  FFT k=4 Koeffizient: |c| = {abs(fft_k4_coeff):.8f}")

        # Vergleich: Verhältnis Projektionsstärken
        ratio = best_ord3["magnitude"] / abs(fft_k2_coeff) if abs(fft_k2_coeff) > 1e-12 else float("inf")
        print(f"  Verhältnis |c_ord3| / |c_FFT_k2| = {ratio:.4f}")

    # Hadamard-ähnliche Zusammenfassung: δ in Charakter-Basis (nur Ordnung 3)
    print("\n--- δ-Zerlegung (nur Ordnung-3-Anteil) ---")
    if ord3_proj:
        recon = np.zeros(6, dtype=complex)
        for r in ord3_proj:
            recon += r["coeff"] * r["chi_values"]
        # Realteil
        print(f"  Rekonstruktion δ (Ord3-Anteil): {[round(x, 6) for x in np.real(recon)]}")
        print(f"  Original δ:                    {[round(x, 6) for x in delta]}")
        frac_var = np.sum(np.real(recon) ** 2) / np.sum(delta ** 2)
        print(f"  Anteil erklärte δ-Varianz (Ord3): {100 * frac_var:.2f}%")

    print("\n" + "=" * 70)
    print("Fertig.")
    print("=" * 70)


if __name__ == "__main__":
    main()
