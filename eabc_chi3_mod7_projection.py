#!/usr/bin/env python3
"""
EABC χ₃ mod 7 Projektionsanalyse
Vergleich des kubischen Dirichlet-Charakters χ₃ mod 7 mit FFT k=2 Modus
auf den 6 Restklassen-Kanälen mod 420.
"""

import numpy as np
from itertools import product
from fractions import Fraction

# ── Daten ──────────────────────────────────────────────────────────────────
r_vals = np.array([11, 101, 191, 221, 311, 401])
counts = np.array([765, 831, 786, 809, 809, 767], dtype=float)
N = counts.sum()
mean = N / 6
delta_count = counts - mean
delta_prob = counts / N - 1 / 6

f_vals = np.array([+1 if r in (101, 221, 311) else -1 for r in r_vals], dtype=float)

# FFT k=2 Referenzwerte (vom Nutzer)
F2_mag_given = 0.016495
F2_var_explained = 0.597  # k=2/k=4

# ── χ₃ mod 7 ───────────────────────────────────────────────────────────────
omega = np.exp(2j * np.pi / 3)
chi3_table = {1: 1.0, 2: omega**2, 3: omega, 4: 1.0, 5: omega, 6: omega**2}

r_mod7 = r_vals % 7
chi3_r = np.array([chi3_table[int(x)] for x in r_mod7], dtype=complex)
chi3_bar_r = np.conj(chi3_r)

print("=" * 70)
print("1. χ₃(r mod 7) für alle 6 Kanäle")
print("=" * 70)
print(f"{'Kanal':>6} {'r':>5} {'r mod 7':>8} {'χ₃':>20} {'|χ₃|':>8} {'arg(χ₃)°':>10}")
for i, r in enumerate(r_vals):
    c = chi3_r[i]
    print(f"{i+1:>6} {r:>5} {r_mod7[i]:>8} {c:>20.6f} {abs(c):>8.4f} {np.degrees(np.angle(c)):>10.2f}")

# C3-Gruppen
print("\nC3-Gruppen (r mod 7):")
groups = {}
for i, r in enumerate(r_vals):
    g = r_mod7[i]
    groups.setdefault(g, []).append(int(r))
for g in sorted(groups):
    print(f"  mod7={g}: {groups[g]}  →  χ₃={chi3_table[g]:.4f}")

# ── Orthonormalisierung auf 6-Punkt-Menge ──────────────────────────────────
def inner_product(a, b):
    """Standard inner product (1/6) Σ a_n conj(b_n)"""
    return np.sum(a * np.conj(b)) / 6

def norm_sq(v):
    return np.real(inner_product(v, v))

def project(v, basis):
    """Projektion auf orthonormalisierte Basisvektoren (Liste)."""
    coeffs = [inner_product(v, b) for b in basis]
    recon = sum(c * b for c, b in zip(coeffs, basis))
    return coeffs, recon

# Konstante + χ₃, χ̄₃ Basis (Gram-Schmidt falls nötig)
ones = np.ones(6, dtype=complex)
# Orthonormalisiere
e0 = ones / np.sqrt(norm_sq(ones))
# χ₃ orthogonal zu e0? Prüfen
chi3_orth = chi3_r - inner_product(chi3_r, e0) * e0
e1 = chi3_orth / np.sqrt(norm_sq(chi3_orth))
chi3bar_orth = chi3_bar_r - inner_product(chi3_bar_r, e0) * e0 - inner_product(chi3_bar_r, e1) * e1
e2 = chi3bar_orth / np.sqrt(norm_sq(chi3bar_orth))

print("\n" + "=" * 70)
print("Orthonormalisierungs-Check")
print("=" * 70)
for name, ei, ej in [("e0,e0", e0, e0), ("e1,e1", e1, e1), ("e2,e2", e2, e2),
                      ("e0,e1", e0, e1), ("e0,e2", e0, e2), ("e1,e2", e1, e2)]:
    print(f"  <{name}> = {np.real(inner_product(ei, ej)):.6f}")

# ── 2. Projektion von δ ───────────────────────────────────────────────────
print("\n" + "=" * 70)
print("2. Projektion von δ auf χ₃-Basis")
print("=" * 70)

# Direkte Projektion (nicht orthonormalisiert, wie in Aufgabe)
c_chi3_count = inner_product(delta_count.astype(complex), chi3_r)
c_chi3bar_count = inner_product(delta_count.astype(complex), chi3_bar_r)
c_chi3_prob = inner_product(delta_prob.astype(complex), chi3_r)
c_chi3bar_prob = inner_product(delta_prob.astype(complex), chi3_bar_r)

# Orthonormalisierte Projektion
coeffs_count, recon_count = project(delta_count.astype(complex), [e0, e1, e2])
coeffs_prob, recon_prob = project(delta_prob.astype(complex), [e0, e1, e2])

var_total_count = norm_sq(delta_count.astype(complex))
var_total_prob = norm_sq(delta_prob.astype(complex))

# Varianz erklärt durch χ₃ + χ̄₃ (ohne Konstante)
var_chi3_count = abs(c_chi3_count)**2 / norm_sq(chi3_r) * norm_sq(chi3_r)  # = |c|^2 wenn nicht ortho
# Korrekter: erklärte Varianz = |<δ, e1>|² + |<δ, e2>|²
var_expl_count_ortho = abs(coeffs_count[1])**2 + abs(coeffs_count[2])**2
var_expl_prob_ortho = abs(coeffs_prob[1])**2 + abs(coeffs_prob[2])**2

# Alternative: nur χ₃ Richtung (reell kombiniert)
# δ ≈ Re(c) Re(χ₃) - Im(c) Im(χ₃) für reelle δ

print("\n--- δ (Counts) ---")
print(f"  c_{'{χ₃}'}     = {c_chi3_count:.6f}")
print(f"  |c_{'{χ₃}'}|²   = {abs(c_chi3_count)**2:.8f}")
print(f"  c_{'{χ̄₃}'}    = {c_chi3bar_count:.6f}")
print(f"  |c_{'{χ̄₃}'}|²  = {abs(c_chi3bar_count)**2:.8f}")
print(f"  |c_{'{χ₃}'}|    = {abs(c_chi3_count):.8f}")
print(f"  arg(c_{'{χ₃}'}) = {np.degrees(np.angle(c_chi3_count)):.2f}°")
print(f"  Gesamtvarianz Var(δ) = {var_total_count:.4f}")
print(f"  Erklärte Varianz (ortho e1+e2) = {var_expl_count_ortho:.8f}")
print(f"  Anteil erklärt = {var_expl_count_ortho/var_total_count:.4f} ({100*var_expl_count_ortho/var_total_count:.2f}%)")

print("\n--- δ (Probability) ---")
print(f"  c_{'{χ₃}'}     = {c_chi3_prob:.8f}")
print(f"  |c_{'{χ₃}'}|²   = {abs(c_chi3_prob)**2:.10f}")
print(f"  |c_{'{χ₃}'}|    = {abs(c_chi3_prob):.10f}")
print(f"  arg(c_{'{χ₃}'}) = {np.degrees(np.angle(c_chi3_prob)):.2f}°")
print(f"  Gesamtvarianz Var(δ) = {var_total_prob:.10f}")
print(f"  Erklärte Varianz (ortho e1+e2) = {var_expl_prob_ortho:.10f}")
print(f"  Anteil erklärt = {var_expl_prob_ortho/var_total_prob:.4f} ({100*var_expl_prob_ortho/var_total_prob:.2f}%)")

# ── 3. FFT k=2 Vergleich ──────────────────────────────────────────────────
print("\n" + "=" * 70)
print("3. Vergleich mit FFT k=2")
print("=" * 70)

# DFT auf 6-Ring: F_k = (1/6) Σ_n δ_n exp(-2πi k n / 6)
def dft_mode(v, k, normalized=True):
    n = np.arange(6)
    s = np.sum(v * np.exp(-2j * np.pi * k * n / 6))
    return s / 6 if normalized else s

F2_count = dft_mode(delta_count, 2)
F2_prob = dft_mode(delta_prob, 2)
F2_f = dft_mode(f_vals, 2)

print(f"\nDFT k=2 von δ (counts):  F₂ = {F2_count:.8f}")
print(f"  |F₂| = {abs(F2_count):.8f}  (gegeben: {F2_mag_given})")
print(f"  arg(F₂) = {np.degrees(np.angle(F2_count)):.2f}°")

print(f"\nDFT k=2 von δ (prob):    F₂ = {F2_prob:.8f}")
print(f"  |F₂| = {abs(F2_prob):.8f}  (unnorm: {abs(dft_mode(delta_prob, 2, False)):.6f})")

print(f"\nDFT k=2 von f(r):        F₂ = {F2_f:.8f}")
print(f"  |F₂| = {abs(F2_f):.8f}")
print(f"  arg(F₂) = {np.degrees(np.angle(F2_f)):.2f}°")

# Vergleich |c_χ₃| vs |F₂|
print(f"\n|Vergleich Magnitude|")
print(f"  |c_{'{χ₃}'}| (count) = {abs(c_chi3_count):.8f}")
print(f"  |F₂| (count)         = {abs(F2_count):.8f}")
print(f"  Ratio |c|/|F₂|      = {abs(c_chi3_count)/abs(F2_count):.4f}")

print(f"\n|Vergleich Phase|")
print(f"  arg(c_{'{χ₃}'}) = {np.degrees(np.angle(c_chi3_count)):.2f}°")
print(f"  arg(F₂)        = {np.degrees(np.angle(F2_count)):.2f}°")
print(f"  Δphase         = {np.degrees(np.angle(c_chi3_count) - np.angle(F2_count)):.2f}°")

# Korrelation f(r) mit k=2 Modus
k2_mode = np.array([np.exp(2j * np.pi * 2 * n / 6) for n in range(6)])
# Reelle Korrelation
corr_f_k2 = np.corrcoef(f_vals, np.real(k2_mode))[0, 1]
corr_f_k2_im = np.corrcoef(f_vals, np.imag(k2_mode))[0, 1]
corr_f_chi3_re = np.corrcoef(f_vals, np.real(chi3_r))[0, 1]
corr_f_chi3_im = np.corrcoef(f_vals, np.imag(chi3_r))[0, 1]

print(f"\nKorrelation f(r) mit k=2-Fouriermodus:")
print(f"  corr(f, Re(e^{{2πi·2n/6}})) = {corr_f_k2:.6f}")
print(f"  corr(f, Im(e^{{2πi·2n/6}})) = {corr_f_k2_im:.6f}")

print(f"\nKorrelation f(r) mit χ₃:")
print(f"  corr(f, Re(χ₃)) = {corr_f_chi3_re:.6f}")
print(f"  corr(f, Im(χ₃)) = {corr_f_chi3_im:.6f}")

# FFT erklärte Varianz berechnen
var_total = norm_sq(delta_prob.astype(complex))
F_all = [dft_mode(delta_prob, k) for k in range(6)]
var_fft = sum(abs(Fk)**2 for Fk in F_all)
print(f"\nFFT Varianz-Check: Σ|F_k|² = {var_fft:.10f}, Var(δ) = {var_total:.10f}")
F2_var = abs(F_all[2])**2
F4_var = abs(F_all[4])**2
print(f"  |F₂|² = {F2_var:.10f}, |F₄|² = {F4_var:.10f}")
print(f"  |F₂|²/|F₄|² = {F2_var/F4_var:.4f} (gegeben ≈ {F2_var_explained})")

# ── 4. Dirichlet-Charaktere mod 420 ───────────────────────────────────────
print("\n" + "=" * 70)
print("4. Dirichlet-Charaktere mod 420 auf den 6 Punkten")
print("=" * 70)

# 420 = 4 * 3 * 5 * 7
# Primitive Charaktere: Produkt von primitiven Charakteren mod jeder Primteilmenge

def chi4(a):
    """Quadratischer Charakter mod 4: χ₄(a) = 0 if even, (-1)^((a-1)/2) if odd"""
    a = int(a) % 4
    if a == 0:
        return 0
    return 1 if a == 1 else -1

def chi3_mod7(a):
    return chi3_table.get(int(a) % 7, 0)

def chi5(a):
    """Quadratischer Charakter mod 5 (Legendre-Symbol)"""
    a = int(a) % 5
    if a == 0:
        return 0
    legendre5 = {1: 1, 2: -1, 3: -1, 4: 1}
    return legendre5.get(a, 0)

def chi3_mod3(a):
    """Kubischer Charakter mod 3: only 1,2 nonzero"""
    a = a % 3
    if a == 0:
        return 0
    return 1  # only one nontrivial but mod 3 group is order 2

def primitive_characters_mod420():
    """Enumerate primitive Dirichlet characters mod 420 as products."""
    chars = []
    # χ mod 4: primitive is χ₄ (index 1), trivial at 0
    chi4_opts = [lambda a, s=1: 1, chi4]  # trivial and primitive
    
    # χ mod 3: only trivial (group order 2)
    chi3_3_opts = [lambda a: 1]
    
    # χ mod 5: primitive quadratic χ₅
    chi5_opts = [lambda a, s=1: 1, chi5]
    
    # χ mod 7: primitive cubic χ₇^(3) and quadratic χ₇^(2)
    # Cubic: kernel of x -> x^3
    omega7 = np.exp(2j * np.pi / 3)
    chi7_cubic = chi3_mod7
    chi7_quad = lambda a: (1 if pow(a % 7, 3, 7) in (1, 6) else -1) if a % 7 != 0 else 0
    # Better quadratic via Legendre
    legendre7 = {}
    for a in range(1, 7):
        legendre7[a] = 1 if pow(a, 3, 7) == 1 else -1  # (a|7)
    def chi7_q(a):
        a = int(a) % 7
        return legendre7.get(a, 0)
    
    chi7_opts = [lambda a: 1, chi7_cubic, chi7_q]
    
    for c4, c3, c5, c7 in product(chi4_opts, chi3_3_opts, chi5_opts, chi7_opts):
        if c4 == (lambda a, s=1: 1) and c5 == (lambda a, s=1: 1) and c7 == (lambda a: 1):
            continue  # skip fully trivial
        def char(a, _c4=c4, _c3=c3, _c5=c5, _c7=c7):
            if np.gcd(int(a), 420) > 1:
                return 0
            return _c4(a) * _c3(a) * _c5(a) * _c7(a)
        chars.append(char)
    return chars

chars = primitive_characters_mod420()
print(f"Anzahl primitiver Charakter-Produkte: {len(chars)}")

# Projektionen von f und δ
results_f = []
results_d = []
for idx, chi in enumerate(chars):
    chi_vals = np.array([chi(r) for r in r_vals], dtype=complex)
    if np.all(chi_vals == 0):
        continue
    c_f = inner_product(f_vals.astype(complex), chi_vals)
    c_d = inner_product(delta_prob.astype(complex), chi_vals)
    results_f.append((abs(c_f), c_f, idx, chi_vals.copy()))
    results_d.append((abs(c_d), c_d, idx, chi_vals.copy()))

results_f.sort(reverse=True)
results_d.sort(reverse=True)

print("\nTop 10 Projektionen |<f, χ>|:")
for i, (mag, c, idx, cv) in enumerate(results_f[:10]):
    print(f"  #{i+1}: |c|={mag:.6f}, c={c:.6f}, χ-Werte={cv}")

print("\nTop 10 Projektionen |<δ, χ>|:")
for i, (mag, c, idx, cv) in enumerate(results_d[:10]):
    print(f"  #{i+1}: |c|={mag:.6f}, c={c:.6f}, χ-Werte={cv}")

# Speziell χ₃ mod 7 induziert (nur χ₇ cubic, Rest trivial)
chi_induced = np.array([chi3_mod7(r) for r in r_vals], dtype=complex)
c_f_ind = inner_product(f_vals.astype(complex), chi_induced)
c_d_ind = inner_product(delta_prob.astype(complex), chi_induced)
print(f"\nχ₃ mod 7 (induziert, Rest trivial):")
print(f"  <f, χ₃> = {c_f_ind:.6f}, |<f, χ₃>| = {abs(c_f_ind):.6f}")
print(f"  <δ, χ₃> = {c_d_ind:.8f}, |<δ, χ₃>| = {abs(c_d_ind):.8f}")

# ── 5. Fit f(r) ≈ a₀ + a₁χ₃ + a₂χ̄₃ ──────────────────────────────────────
print("\n" + "=" * 70)
print("5. Fit f(r) ≈ a₀ + a₁χ₃ + a₂χ̄₃")
print("=" * 70)

# Design matrix: [1, Re(χ₃), Im(χ₃)] since f is real
X = np.column_stack([
    np.ones(6),
    np.real(chi3_r),
    np.imag(chi3_r)
])
# Least squares
coeffs, residuals, rank, sv = np.linalg.lstsq(X, f_vals, rcond=None)
a0, a1_re, a1_im = coeffs
a1 = a1_re + 1j * a1_im  # a₁ on χ₃, but we have Re/Im decomposition

f_fit = X @ coeffs
ss_res = np.sum((f_vals - f_fit)**2)
ss_tot = np.sum((f_vals - f_vals.mean())**2)
R2 = 1 - ss_res / ss_tot

print(f"  a₀ = {a0:.6f}")
print(f"  Koeffizient Re(χ₃): {a1_re:.6f}")
print(f"  Koeffizient Im(χ₃): {a1_im:.6f}")
print(f"  a₁ (komplex, auf χ₃): {a1_re + 1j*a1_im:.6f}")
print(f"  f_fit = {f_fit}")
print(f"  f_obs = {f_vals}")
print(f"  Residuen = {f_vals - f_fit}")
print(f"  R² = {R2:.6f}")

# Alternative: a₀ + a₁χ₃ + a₂χ̄₃ with complex coeffs via projection
a0_proj = inner_product(f_vals.astype(complex), e0)
a1_proj = inner_product(f_vals.astype(complex), e1)
a2_proj = inner_product(f_vals.astype(complex), e2)
f_fit2 = a0_proj * e0 + a1_proj * e1 + a2_proj * e2
ss_res2 = np.sum(np.abs(f_vals - np.real(f_fit2))**2)
R2_2 = 1 - ss_res2 / ss_tot

print(f"\n  (Orthonormal) a₀={a0_proj:.6f}, a₁={a1_proj:.6f}, a₂={a2_proj:.6f}")
print(f"  R² (ortho) = {R2_2:.6f}")

# ── Zusammenfassung ─────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("ZUSAMMENFASSUNG / FAZIT")
print("=" * 70)
print(f"""
χ₃-Tabelle auf 6 Kanälen:
  r mod 7 ∈ {{4,3,2}} → χ₃ ∈ {{1, ω, ω²}}

|⟨δ, χ₃⟩| (prob)  = {abs(c_chi3_prob):.8f}
|F₂| (prob)        = {abs(F2_prob):.8f}
Verhältnis         = {abs(c_chi3_prob)/abs(F2_prob):.4f}

Erklärte Varianz δ durch χ₃+χ̄₃: {100*var_expl_prob_ortho/var_total_prob:.2f}%
Erklärte Varianz δ durch FFT k=2:  {100*F2_var/var_total:.2f}%
FFT k=2/k=4 Ratio:                  {F2_var/F4_var:.4f}

f(r) Fit R² auf χ₃-Basis: {R2:.4f}
corr(f, Re(χ₃)) = {corr_f_chi3_re:.4f}
""")
