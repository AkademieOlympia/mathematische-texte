#!/usr/bin/env sage
# -*- coding: utf-8 -*-
# Lift von Pi_Gamma nach QQ(zeta_5) und Frobenius-Test (Pfad A, Schritt 1).
# Start: ./run_lift_zeta5.sh  oder  sage rundweg_lift_zeta5.sage

from sage.all import Matrix, QQ, CyclotomicField, vector


def init_operators_qq():
    """R, K, Gamma wie in Rundweg.py über QQ."""
    R = Matrix(QQ, [[0, 0, 0, 1],
                    [1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0]])
    K = Matrix(QQ, [[0, 0, 0, 1],
                    [0, 0, 1, 0],
                    [0, 1, 0, 0],
                    [1, 0, 0, 0]])
    Gamma = R * K
    return R, K, Gamma


def berechne_pi_gamma(Gamma, ring):
    """Pi_Gamma = (I + Gamma + Gamma^2 + Gamma^3) / 4."""
    I = Matrix.identity(ring, 4)
    return (I + Gamma + Gamma**2 + Gamma**3) / 4


def lift_matrix(M_qq, K_field):
    return Matrix(K_field, M_qq)


def eigenspace_daten(Pi):
    """Liste (lambda, dim, basis) aus eigenspaces_right."""
    out = []
    for ev, raum in Pi.eigenspaces_right():
        out.append((ev, raum.dimension(), raum.basis()))
    return sorted(out, key=lambda t: str(t[0]))


def sigma_auf_matrix(sigma, M):
    return M.apply_map(sigma)


def sigma_auf_vektor(sigma, v):
    return vector(v.base_ring(), [sigma(c) for c in v])


# --- Körper und Operatoren ---
K_field = CyclotomicField(5)
zeta = K_field.gen()

R_qq, K_qq, Gamma_qq = init_operators_qq()
Pi_qq = berechne_pi_gamma(Gamma_qq, QQ)

Gamma_K = lift_matrix(Gamma_qq, K_field)
Pi_K = berechne_pi_gamma(Gamma_K, K_field)

print("=== LIFT QQ(zeta_5) — Rundweg Pi_Gamma (Schritt 1) ===")
print(f"K = QQ(zeta_5), zeta = {zeta}, Grad [K:Q] = {K_field.degree()}")
print()

daten_qq = eigenspace_daten(Pi_qq)
daten_K = eigenspace_daten(Pi_K)

print("=== EIGENRÄUME: QQ vs. QQ(zeta_5) ===")
for (lam_q, dim_q, bas_q), (lam_k, dim_k, bas_k) in zip(daten_qq, daten_K):
    print(f"Eigenwert λ = {lam_q}")
    print(f"  Dimension über Q: {dim_q}  |  über K: {dim_k}  |  gleich: {dim_q == dim_k}")
    raum_lam_K = next(raum for ev, raum in Pi_K.eigenspaces_right() if ev == lam_k)
    for i, bq in enumerate(bas_q, start=1):
        bq_K = vector(K_field, bq)
        print(f"  QQ-Basis {i}: {tuple(bq)}  -> liegt in K-Eigenraum λ: {bq_K in raum_lam_K}")
    for i, bk in enumerate(bas_k, start=1):
        print(f"  Sage K-Basis {i}: {tuple(bk)}")
    print()

print("-" * 50)
print("Pi_Gamma über K (Lift):")
print(Pi_K)
print("-" * 50)

# --- Frobenius / Galois bei p = 5 ---
G = K_field.galois_group()
sigma = G.gen()
p = 5

print("=== FROBENIUS / GALOIS (p = 5) ===")
print(f"|Gal(K/Q)| = {G.order()}")
print(f"sigma = Galois-Generator: {sigma}")
print(f"sigma(zeta) = {sigma(zeta)}")
print(f"zeta^{p} = zeta^5 = {zeta**p}  (nicht Automorphismus: zeta |-> 1)")
print()

v0_K = vector(K_field, [0, 1, 0, -1])
sigma_v0 = sigma_auf_vektor(sigma, v0_K)
print("Wirkung von sigma auf Basis des λ=0-Raums (0,1,0,-1):")
print(f"  v0        = {tuple(v0_K)}")
print(f"  sigma(v0) = {tuple(sigma_v0)}")
print(f"  Fixpunkt? {sigma_v0 == v0_K}")
print()

raum1_K = next(raum for ev, raum in Pi_K.eigenspaces_right() if ev == K_field(1))
print(f"im(Pi_Gamma) ≅ Eigenspace λ=1, Dimension = {raum1_K.dimension()}")
for i, b in enumerate(raum1_K.basis(), start=1):
    sb = sigma_auf_vektor(sigma, b)
    print(f"  Basis {i}: {tuple(b)}")
    print(f"    sigma(Basis {i}) = {tuple(sb)}  |  unverändert? {sb == b}")
print()

sigma_Pi = sigma_auf_matrix(sigma, Pi_K)
print("sigma auf Pi_Gamma (koeffizientweise):")
print(f"  sigma(Pi) == Pi ? {sigma_Pi == Pi_K}")
print()

frobenius_nichtrivial_auf_K = sigma(zeta) != zeta
frobenius_nichtrivial_auf_basen = any(
    sigma_auf_vektor(sigma, vector(K_field, b)) != vector(K_field, b)
    for _, _, bas in daten_K
    for b in bas
)
print("=== FAZIT ===")
print(f"Galois nichttrivial auf zeta: {frobenius_nichtrivial_auf_K}  ({zeta} -> {sigma(zeta)})")
print(f"Galois nichttrivial auf Eigenraum-Basen (Sage): {frobenius_nichtrivial_auf_basen}")
print("QQ-Eigenräume bleiben unter sigma punktweise fix (erwartet für Gal(K/Q)-Wirkung).")
