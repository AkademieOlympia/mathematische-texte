#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Datei: Kwant Alpha Delta.py
#
#   python3 "Kwant Alpha Delta.py"
# oder: python3 kwant_alpha_delta.py

import numpy as np
import math

# -----------------------------
# Eingabedaten
# -----------------------------

N = 42
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

alpha_phys_inv = 137.035999084
alpha_phys = 1 / alpha_phys_inv

# -----------------------------
# Grundgrößen
# -----------------------------

positive_E = spectrum[spectrum > 1e-12]
gap = positive_E[0]

M = len(spectrum)
dk = 2 * math.pi / M
c_eff = gap / dk
rho_eff = 1 / gap

p0 = channels[0]
channel_mean = np.mean(channels)
channel_std = np.std(channels)
disorder = channel_std / channel_mean
T_eff = 1 / (1 + disorder**2)

base_inv = 3 * N + p0
target_delta = alpha_phys_inv - base_inv

# -----------------------------
# Kandidaten für Delta
# -----------------------------

delta_candidates = {}

delta_candidates["gap"] = gap
delta_candidates["gap^2"] = gap**2
delta_candidates["gap^3"] = gap**3

delta_candidates["1/rho"] = 1 / rho_eff
delta_candidates["1/(rho*T)"] = 1 / (rho_eff * T_eff)
delta_candidates["1/(rho*c)"] = 1 / (rho_eff * c_eff)
delta_candidates["1/(rho*c*T)"] = 1 / (rho_eff * c_eff * T_eff)

delta_candidates["T*gap"] = T_eff * gap
delta_candidates["T*gap^2"] = T_eff * gap**2
delta_candidates["gap^2/T"] = gap**2 / T_eff
delta_candidates["gap^3/T"] = gap**3 / T_eff

delta_candidates["disorder"] = disorder
delta_candidates["disorder^2"] = disorder**2
delta_candidates["disorder/M"] = disorder / M
delta_candidates["disorder^2/M"] = disorder**2 / M

delta_candidates["gap/M"] = gap / M
delta_candidates["gap^2/M"] = gap**2 / M
delta_candidates["gap^3/M"] = gap**3 / M

delta_candidates["T/M"] = T_eff / M
delta_candidates["(1-T)/M"] = (1 - T_eff) / M
delta_candidates["(1-T)^2"] = (1 - T_eff)**2

delta_candidates["1/(N)"] = 1 / N
delta_candidates["1/(N*pi)"] = 1 / (N * math.pi)
delta_candidates["1/(3N+p0)"] = 1 / base_inv
delta_candidates["1/(channels_mean)"] = 1 / channel_mean
delta_candidates["1/(channels_std)"] = 1 / channel_std

# -----------------------------
# Kombinierte Kandidaten
# -----------------------------

for name, val in list(delta_candidates.items()):
    delta_candidates[f"{name}/pi"] = val / math.pi
    delta_candidates[f"{name}/2pi"] = val / (2 * math.pi)
    delta_candidates[f"{name}/12"] = val / 12
    delta_candidates[f"{name}/42"] = val / 42
    delta_candidates[f"{name}/137"] = val / 137

# -----------------------------
# Auswertung
# -----------------------------

print("EABC/Kwant-Alpha-Delta-Test")
print("===========================")
print()
print(f"N                 = {N}")
print(f"channels          = {channels.astype(int).tolist()}")
print(f"spectrum length   = {M}")
print()
print(f"gap               = {gap:.12f}")
print(f"c_eff             = {c_eff:.12f}")
print(f"rho_eff           = {rho_eff:.12f}")
print(f"T_eff             = {T_eff:.12f}")
print(f"disorder          = {disorder:.12f}")
print()
print(f"alpha_phys^-1     = {alpha_phys_inv:.12f}")
print(f"base_inv=3N+p0    = {base_inv:.12f}")
print(f"target_delta      = {target_delta:.12f}")
print()

results = []

for name, delta in delta_candidates.items():
    alpha_inv_model = base_inv + delta
    alpha_model = 1 / alpha_inv_model

    abs_error_inv = abs(alpha_inv_model - alpha_phys_inv)
    rel_error_alpha = abs(alpha_model - alpha_phys) / alpha_phys

    results.append((abs_error_inv, rel_error_alpha, name, delta, alpha_inv_model, alpha_model))

results.sort(key=lambda x: x[0])

print("Beste Delta-Kandidaten:")
print("-----------------------")
for abs_error_inv, rel_error_alpha, name, delta, alpha_inv_model, alpha_model in results[:25]:
    print(
        f"{name:25s} "
        f"delta={delta:.12f}  "
        f"alpha_inv={alpha_inv_model:.12f}  "
        f"abs_err_inv={abs_error_inv:.12f}  "
        f"rel_err_alpha={rel_error_alpha:.12e}"
    )

print()
print("Direkter Treffervergleich:")
print("--------------------------")
print(f"1/137                 = {1/137:.15f}")
print(f"1/137.035999084       = {alpha_phys:.15f}")
print(f"Unterschied           = {abs(1/137-alpha_phys):.15e}")
