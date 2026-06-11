#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# EABC/Kwant-Alpha-Verbesserungstest
#   python3 "Kwant Alpha Verbesserung.py"

import numpy as np
import math

print("EABC/Kwant-Alpha-Verbesserungstest")
print("==================================")

# ------------------------------------------------------------
# 1. Grunddaten
# ------------------------------------------------------------

N = 42
M = 7

channels = np.array([11, 17, 29, 41, 47, 59], dtype=float)

spectrum = np.array([
    -1.84775907,
    -1.41421356,
    -0.76536686,
    0.0,
    0.76536686,
    1.41421356,
    1.84775907
], dtype=float)

alpha_inv_old = 137.035999084
alpha_inv_codata2022 = 137.035999177

# ------------------------------------------------------------
# 2. Spektrale Größen
# ------------------------------------------------------------

positive_E = spectrum[spectrum > 1e-12]
gap = positive_E[0]

dk = 2 * math.pi / len(spectrum)
c_eff = gap / dk
rho_eff = 1 / gap

p0 = channels[0]
p_mean = np.mean(channels)
p_std = np.std(channels)

disorder = p_std / p_mean
disorder2 = disorder**2

T = 1 / (1 + disorder2)

B = 3 * N + p0
D = disorder2 / M
R = (1 - T) / B

alpha_inv_model = B + D + R
alpha_model = 1 / alpha_inv_model

# ------------------------------------------------------------
# 3. Ausgabe der Hauptformel
# ------------------------------------------------------------

print()
print("Hauptformel")
print("-----------")
print(f"N                     = {N}")
print(f"M                     = {M}")
print(f"p0                    = {p0:.0f}")
print(f"Kanäle                = {channels.astype(int).tolist()}")
print()
print(f"B = 3N+p0             = {B:.12f}")
print(f"mean(p)               = {p_mean:.12f}")
print(f"std(p)                = {p_std:.12f}")
print(f"disorder=std/mean     = {disorder:.12f}")
print(f"disorder^2            = {disorder2:.12f}")
print(f"T=1/(1+disorder^2)    = {T:.12f}")
print()
print(f"D=disorder^2/M        = {D:.12f}")
print(f"R=(1-T)/B             = {R:.12f}")
print()
print(f"alpha_inv_model       = {alpha_inv_model:.12f}")
print(f"alpha_model           = {alpha_model:.15f}")

# ------------------------------------------------------------
# 4. Vergleich gegen Referenzen
# ------------------------------------------------------------


def compare(label, target_inv):
    target_alpha = 1 / target_inv
    abs_err_inv = abs(alpha_inv_model - target_inv)
    rel_err_alpha = abs(alpha_model - target_alpha) / target_alpha

    print()
    print(label)
    print("-" * len(label))
    print(f"target alpha^-1       = {target_inv:.12f}")
    print(f"target alpha          = {target_alpha:.15f}")
    print(f"abs error alpha^-1    = {abs_err_inv:.12e}")
    print(f"rel error alpha       = {rel_err_alpha:.12e}")


compare("Vergleich gegen bisherige Referenz", alpha_inv_old)
compare("Vergleich gegen CODATA 2022", alpha_inv_codata2022)

# ------------------------------------------------------------
# 5. Varianten-Test: Welche Restform ist am besten?
# ------------------------------------------------------------

print()
print("Varianten-Test")
print("--------------")

variants = {
    "B": B,
    "B + D": B + D,
    "B + D + (1-T)/B": B + D + (1 - T) / B,
    "B + D + 1/(B*M)": B + D + 1 / (B * M),
    "B + D + (1-T)/(B*T)": B + D + (1 - T) / (B * T),
    "B + D + (1-T)^2/B": B + D + (1 - T) ** 2 / B,
    "B + D + disorder^4/(M*pi)": B + D + disorder**4 / (M * math.pi),
    "B + D + 1/(B*M*pi)": B + D + 1 / (B * M * math.pi),
}

target = alpha_inv_codata2022
rows = []

for name, inv_val in variants.items():
    alpha_val = 1 / inv_val
    abs_err = abs(inv_val - target)
    rel_err = abs(alpha_val - 1 / target) / (1 / target)
    rows.append((abs_err, rel_err, name, inv_val, alpha_val))

rows.sort(key=lambda x: x[0])

for abs_err, rel_err, name, inv_val, alpha_val in rows:
    print(
        f"{name:32s} "
        f"alpha^-1={inv_val:.12f}  "
        f"abs_err={abs_err:.12e}  "
        f"rel_err={rel_err:.12e}"
    )

# ------------------------------------------------------------
# 6. Restanalyse
# ------------------------------------------------------------

print()
print("Restanalyse")
print("-----------")

rest_after_D = target - (B + D)
rest_after_DR = target - (B + D + R)

print(f"Ziel CODATA 2022       = {target:.12f}")
print(f"B + D                  = {B + D:.12f}")
print(f"Rest nach B+D          = {rest_after_D:.12e}")
print(f"B + D + R              = {B + D + R:.12f}")
print(f"Rest nach B+D+R        = {rest_after_DR:.12e}")

print()
print("Fertig.")
