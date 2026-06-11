#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Kwant/EABC-Feinstruktur-Prototyp
#
# Ausführen (Leerzeichen im Namen → immer in Anführungszeichen):
#   cd "/Pfad/zu/Mathematische Texte"
#   python3 "Kwant Alpha.py"
#
# Alternative ohne Leerzeichen im Dateinamen:
#   python3 kwant_alpha_eabc.py

import numpy as np
import math
from itertools import combinations

# -----------------------------
# 1. Deine Daten
# -----------------------------

N = 42
L = 6          # Kanäle
M = 7          # Phasen

channels = np.array([11, 17, 29, 41, 47, 59])
phases = np.arange(M)

# Spektrum aus deinem Lauf:
spectrum = np.array([
    -1.84775907,
    -1.41421356,
    -0.76536686,
    0.0,
    0.76536686,
    1.41421356,
    1.84775907
])

alpha_phys = 1 / 137.035999084


# -----------------------------
# 2. Effektive Modenstruktur
# -----------------------------

# Wir interpretieren die 7 Phasen als k-Moden:
k = 2 * np.pi * phases / M

# Zentriere k um 0
k_centered = np.fft.fftshift(k - np.pi)
E_centered = np.sort(spectrum)

# Für eine einfache Dirac-artige Steigung:
# benutze die beiden kleinsten positiven Energien
positive_E = spectrum[spectrum > 1e-9]
c_eff = positive_E[0] / (2 * np.pi / M)

# -----------------------------
# 3. Kanalgewichte aus Primzahlen
# -----------------------------

# Normierte Primkanal-Gewichte
weights = channels / np.sum(channels)

# Varianz / Unordnung der Kanäle
channel_mean = np.mean(channels)
channel_std = np.std(channels)
disorder = channel_std / channel_mean

# -----------------------------
# 4. Einfache Transmission
# -----------------------------

# Heuristische Transmission:
# je geordneter die Kanäle, desto höher die Transmission.
T_eff = 1 / (1 + disorder**2)

# -----------------------------
# 5. Zustandsdichte nahe Null
# -----------------------------

gap = positive_E[0]
rho_eff = 1 / gap

# -----------------------------
# 6. Kandidaten für alpha_N
# -----------------------------

# Verschiedene plausible dimensionslose Kombinationen
alpha_candidates = {}

alpha_candidates["T / (rho*c)"] = T_eff / (rho_eff * c_eff)

alpha_candidates["gap / (c*T)"] = gap / (c_eff * T_eff)

alpha_candidates["T * gap / c"] = T_eff * gap / c_eff

alpha_candidates["1 / (N*pi)"] = 1 / (N * math.pi)

alpha_candidates["1 / (N*pi + channels_mean)"] = 1 / (N * math.pi + channel_mean)

alpha_candidates["1 / (N*pi + std)"] = 1 / (N * math.pi + channel_std)

alpha_candidates["1 / (12^2 - 7)"] = 1 / (12**2 - 7)

alpha_candidates["1 / (3*N + 11)"] = 1 / (3*N + 11)


# -----------------------------
# 7. Ausgabe
# -----------------------------

print("EABC/Kwant-Alpha-Prototyp")
print("=========================")
print(f"N = {N} = {L} x {M}")
print(f"Kanäle: {channels.tolist()}")
print(f"Phasen: {phases.tolist()}")
print()
print("Spektrum:")
for x in spectrum:
    print(f"{x: .8f}")

print()
print("Effektive Größen:")
print(f"gap       = {gap:.10f}")
print(f"c_eff     = {c_eff:.10f}")
print(f"rho_eff   = {rho_eff:.10f}")
print(f"T_eff     = {T_eff:.10f}")
print(f"disorder  = {disorder:.10f}")
print()

print("Physikalische Feinstrukturkonstante:")
print(f"alpha_phys     = {alpha_phys:.12f}")
print(f"1/alpha_phys   = {1/alpha_phys:.9f}")
print()

print("Alpha-Kandidaten:")
ranked = []
for name, val in alpha_candidates.items():
    inv = 1 / val
    error = abs(val - alpha_phys) / alpha_phys
    ranked.append((error, name, val, inv))
    print(f"{name:30s} alpha={val:.12f}   1/alpha={inv:.6f}   rel.error={error:.6f}")

ranked.sort(key=lambda x: x[0])
best_err, best_name, best_val, best_inv = ranked[0]

print()
print("─" * 61)
print("Zusammenfassung (nach kleinstem rel. Fehler |alpha - alpha_phys|/alpha_phys)")
print("─" * 61)
print(f"Bester Kandidat : {best_name}")
print(f"  alpha_modell  = {best_val:.12f}")
print(f"  1/alpha_modell= {best_inv:.9f}")
print(f"  rel. Fehler   = {best_err:.12e}")
print(f"Referenz alpha_phys = {alpha_phys:.12f}")
if best_err < 1e-3:
    print("Hinweis: rel. Fehler < 1e-3 — heuristische Übereinstimmung mit alpha_phys.")
else:
    print("Hinweis: Heuristik trifft alpha_phys nur grob oder gar nicht (siehe Tabelle).")

# 8. Delta-Test (base_inv = 137.0 → alpha_inv = base_inv + delta)
print()
print("─────────────────────────────────────────────────────────────")
print("Delta-Test für Feinstrukturkorrektur")
print("─────────────────────────────────────────────────────────────")

alpha_phys_inv = 137.035999084
base_inv = 137.0
target_delta = alpha_phys_inv - base_inv

delta_candidates = {
    "gap/M": gap / len(spectrum),
    "gap^2/M": gap**2 / len(spectrum),
    "gap^3/M": gap**3 / len(spectrum),
    "(1-T_eff)/M": (1 - T_eff) / len(spectrum),
    "(1-T_eff)^2": (1 - T_eff) ** 2,
    "disorder^2/M": disorder**2 / len(spectrum),
    "1/(N*pi)": 1 / (N * math.pi),
    "1/(3N+11)": 1 / (3 * N + 11),
}

results = []

for name, delta in delta_candidates.items():
    alpha_inv_model = base_inv + delta
    alpha_model = 1 / alpha_inv_model
    abs_error_delta = abs(delta - target_delta)
    rel_error_alpha = abs(alpha_model - alpha_phys) / alpha_phys

    results.append(
        (
            abs_error_delta,
            rel_error_alpha,
            name,
            delta,
            alpha_inv_model,
            alpha_model,
        )
    )

results.sort(key=lambda x: x[0])

print(f"Ziel-delta       = {target_delta:.12f}")
print()

for abs_error_delta, rel_error_alpha, name, delta, alpha_inv_model, alpha_model in results:
    print(
        f"{name:18s} "
        f"delta={delta:.12f}  "
        f"1/alpha={alpha_inv_model:.12f}  "
        f"delta_err={abs_error_delta:.12f}  "
        f"rel_alpha_err={rel_error_alpha:.12e}"
    )

# 9. Zweitkorrektur (delta_1 = beste 1. Ordnung, Rest via eps)
print()
print("─────────────────────────────────────────────────────────────")
print("Zweitkorrektur-Test für Delta")
print("─────────────────────────────────────────────────────────────")

delta_1 = disorder**2 / len(spectrum)
target_delta = 137.035999084 - 137.0
rest_delta = target_delta - delta_1

print(f"delta_1 = disorder^2/M = {delta_1:.12f}")
print(f"target_delta          = {target_delta:.12f}")
print(f"rest_delta            = {rest_delta:.12f}")
print()

M_spec = len(spectrum)

second_order_candidates = {
    "(1-T)^2/137": (1 - T_eff) ** 2 / 137,
    "(1-T)^2/M": (1 - T_eff) ** 2 / M_spec,
    "(1-T)^2/(M*pi)": (1 - T_eff) ** 2 / (M_spec * math.pi),
    "(1-T)/137": (1 - T_eff) / 137,
    "(1-T)/(M*137)": (1 - T_eff) / (M_spec * 137),
    "disorder^4": disorder**4,
    "disorder^4/M": disorder**4 / M_spec,
    "disorder^4/(M*pi)": disorder**4 / (M_spec * math.pi),
    "gap^4/M": gap**4 / M_spec,
    "gap^4/(M*pi)": gap**4 / (M_spec * math.pi),
    "1/(137*M)": 1 / (137 * M_spec),
    "1/(137*M*pi)": 1 / (137 * M_spec * math.pi),
    "1/(137*M^2)": 1 / (137 * M_spec**2),
    "1/(137^2)": 1 / (137**2),
    "1/(137^2*M)": 1 / (137**2 * M_spec),
}

results2 = []

for name, eps in second_order_candidates.items():
    delta_total = delta_1 + eps
    alpha_inv_model = 137 + delta_total
    alpha_model = 1 / alpha_inv_model

    abs_err_delta = abs(delta_total - target_delta)
    rel_err_alpha = abs(alpha_model - alpha_phys) / alpha_phys

    results2.append(
        (
            abs_err_delta,
            rel_err_alpha,
            name,
            eps,
            delta_total,
            alpha_inv_model,
        )
    )

results2.sort(key=lambda x: x[0])

print("Beste Zweitkorrekturen:")
print("-----------------------")

for abs_err_delta, rel_err_alpha, name, eps, delta_total, alpha_inv_model in results2:
    print(
        f"{name:22s} "
        f"eps={eps:.12f}  "
        f"delta_total={delta_total:.12f}  "
        f"1/alpha={alpha_inv_model:.12f}  "
        f"delta_err={abs_err_delta:.12e}  "
        f"rel_alpha_err={rel_err_alpha:.12e}"
    )

# 10. Gedämpfte Zweitkorrektur: alpha_inv ≈ B + D + λ R
print()
print("─────────────────────────────────────────────────────────────")
print("Lambda-Test für gedämpfte Zweitkorrektur")
print("─────────────────────────────────────────────────────────────")

alpha_inv_target = 137.035999084

M_lam = len(spectrum)
B = 3 * N + 11

D = disorder**2 / M_lam
T = T_eff
R = (1 - T) / B

base = B + D
needed_R_factor = (alpha_inv_target - base) / R

print(f"Ziel alpha^-1        = {alpha_inv_target:.12f}")
print(f"B                    = {B:.12f}")
print(f"D=disorder^2/M       = {D:.12f}")
print(f"R=(1-T)/B            = {R:.12f}")
print(f"B+D                  = {base:.12f}")
print(f"benötigtes lambda    = {needed_R_factor:.12f}")
print()

lambda_candidates = {
    "1": 1.0,
    "T": T,
    "sqrt(T)": math.sqrt(T),
    "T^2": T**2,
    "1/(1+D)": 1 / (1 + D),
    "1/(1+disorder^2)": 1 / (1 + disorder**2),
    "M/(M+1)": M_lam / (M_lam + 1),
    "(M-1)/M": (M_lam - 1) / M_lam,
    "B/(B+1)": B / (B + 1),
    "(B-1)/B": (B - 1) / B,
    "137/(137+1)": 137 / 138,
    "1-D": 1 - D,
    "1-disorder^2/M": 1 - disorder**2 / M_lam,
    "1/(1+1/M)": 1 / (1 + 1 / M_lam),
    "1/(1+1/B)": 1 / (1 + 1 / B),
    "T + D": T + D,
    "T + R": T + R,
    "sqrt(T)/(1+D)": math.sqrt(T) / (1 + D),
    "T/(1+D)": T / (1 + D),
    "M*T/(M+1)": M_lam * T / (M_lam + 1),
}

lambda_results = []

for name, lam in lambda_candidates.items():
    alpha_inv_model = B + D + lam * R
    alpha_model = 1 / alpha_inv_model

    abs_err_inv = abs(alpha_inv_model - alpha_inv_target)
    rel_err_alpha = abs(alpha_model - alpha_phys) / alpha_phys

    lambda_results.append(
        (abs_err_inv, rel_err_alpha, name, lam, alpha_inv_model)
    )

lambda_results.sort(key=lambda x: x[0])

print("Beste Lambda-Kandidaten:")
print("------------------------")

for abs_err_inv, rel_err_alpha, name, lam, alpha_inv_model in lambda_results:
    print(
        f"{name:22s} "
        f"lambda={lam:.12f}  "
        f"alpha^-1={alpha_inv_model:.12f}  "
        f"abs_err={abs_err_inv:.12e}  "
        f"rel_alpha_err={rel_err_alpha:.12e}"
    )


def is_prime(n):
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    r = int(math.sqrt(n))
    for k in range(3, r + 1, 2):
        if n % k == 0:
            return False
    return True


def generate_primes(limit):
    return [n for n in range(2, limit + 1) if is_prime(n)]


def alpha_inv_eabc(N, M, channels, p0_mode="first"):
    channels = np.array(channels, dtype=float)

    if p0_mode == "first":
        p0 = channels[0]
    elif p0_mode == "min":
        p0 = np.min(channels)
    elif p0_mode == "mod11":
        candidates = [p for p in channels if int(p) % 12 == 11]
        p0 = candidates[0] if candidates else channels[0]
    else:
        p0 = channels[0]

    mean_p = np.mean(channels)
    std_p = np.std(channels)

    disorder = std_p / mean_p
    disorder2 = disorder**2

    B = 3 * N + p0
    D = disorder2 / M
    T = 1 / (1 + disorder2)
    R = (1 - T) / B

    alpha_inv = B + D + math.sqrt(T) * R

    return {
        "alpha_inv": alpha_inv,
        "B": B,
        "D": D,
        "T": T,
        "R": R,
        "sqrtT": math.sqrt(T),
        "disorder": disorder,
        "mean": mean_p,
        "std": std_p,
        "p0": p0,
    }


# 11. Stabilitätstest über viele Primkanal-Gruppen
print()
print("─────────────────────────────────────────────────────────────")
print("Stabilitätstest über viele Primkanal-Gruppen")
print("─────────────────────────────────────────────────────────────")

alpha_inv_target = 137.035999084
group_size = 6
M_test = 7

N_test = 42

primes = generate_primes(500)

results_stab = []

for i in range(0, len(primes) - group_size + 1):
    group = primes[i:i + group_size]

    data = alpha_inv_eabc(N_test, M_test, group, p0_mode="first")
    alpha_inv_model = data["alpha_inv"]

    abs_err = abs(alpha_inv_model - alpha_inv_target)
    rel_err = abs((1 / alpha_inv_model) - (1 / alpha_inv_target)) / (
        1 / alpha_inv_target
    )

    row = {
        "i": i,
        "group": group,
        "alpha_inv": alpha_inv_model,
        "abs_err": abs_err,
        "rel_err": rel_err,
        **data,
    }
    results_stab.append(row)

results_stab.sort(key=lambda x: x["abs_err"])

print(f"Ziel alpha^-1 = {alpha_inv_target:.12f}")
print(f"N             = {N_test}")
print(f"M             = {M_test}")
print(f"Gruppengröße  = {group_size}")
print()

print("Beste Gruppen:")
print("--------------")

for r in results_stab[:30]:
    print(
        f"i={r['i']:3d} "
        f"group={r['group']} "
        f"p0={r['p0']:.0f} "
        f"B={r['B']:.6f} "
        f"D={r['D']:.12f} "
        f"sqrtT={r['sqrtT']:.12f} "
        f"alpha^-1={r['alpha_inv']:.12f} "
        f"abs_err={r['abs_err']:.12e} "
        f"rel_err={r['rel_err']:.12e}"
    )

your_group = [11, 17, 29, 41, 47, 59]
your_data = alpha_inv_eabc(42, 7, your_group)

print()
print("Deine Gruppe:")
print("-------------")
print(f"group         = {your_group}")
print(f"alpha^-1      = {your_data['alpha_inv']:.12f}")
print(
    f"abs_err       = {abs(your_data['alpha_inv'] - alpha_inv_target):.12e}"
)
print(f"disorder      = {your_data['disorder']:.12f}")
print(f"D             = {your_data['D']:.12f}")
print(f"T             = {your_data['T']:.12f}")
print(f"sqrtT         = {your_data['sqrtT']:.12f}")
print(f"R             = {your_data['R']:.12f}")

print()
print("Fertig Stabilitätstest.")

# 12. Kombinatorischer Test: alle 6er mit p0 = 11 aus Primzahlen <= 59
print()
print("─────────────────────────────────────────────────────────────")
print(
    "Kombinatorischer Test: alle 6er-Gruppen aus Primzahlen "
    "<= 59 mit p0=11"
)
print("─────────────────────────────────────────────────────────────")

prime_pool = [p for p in generate_primes(59) if p >= 11]

combo_results = []

for group in combinations(prime_pool, 6):
    if group[0] != 11:
        continue

    data = alpha_inv_eabc(42, 7, list(group))
    alpha_inv_model = data["alpha_inv"]

    abs_err = abs(alpha_inv_model - alpha_inv_target)
    rel_err = abs((1 / alpha_inv_model) - (1 / alpha_inv_target)) / (
        1 / alpha_inv_target
    )

    residues = [p % 12 for p in group]

    combo_results.append(
        {
            "group": list(group),
            "residues": residues,
            "alpha_inv": alpha_inv_model,
            "abs_err": abs_err,
            "rel_err": rel_err,
            **data,
        }
    )

combo_results.sort(key=lambda x: x["abs_err"])

print(f"Anzahl getesteter Gruppen: {len(combo_results)}")
print()
print("Beste 30 Kombinationen:")
print("-----------------------")

for r in combo_results[:30]:
    print(
        f"group={r['group']} "
        f"mod12={r['residues']} "
        f"D={r['D']:.12f} "
        f"sqrtT={r['sqrtT']:.12f} "
        f"alpha^-1={r['alpha_inv']:.12f} "
        f"abs_err={r['abs_err']:.12e}"
    )

print()
print("Rang deiner Gruppe:")
print("-------------------")

your_group = [11, 17, 29, 41, 47, 59]

for rank, r in enumerate(combo_results, start=1):
    if r["group"] == your_group:
        print(f"Rang deiner Gruppe: {rank} von {len(combo_results)}")
        print(f"group        = {r['group']}")
        print(f"mod12        = {r['residues']}")
        print(f"alpha^-1     = {r['alpha_inv']:.12f}")
        print(f"abs_err      = {r['abs_err']:.12e}")
        print(f"rel_err      = {r['rel_err']:.12e}")
        break

print()
print("Fertig kombinatorischer Test.")

# 13. Erweiterter Kombinatorik-Test: Prime-Pools bis verschiedene Grenzen
print()
print("─────────────────────────────────────────────────────────────")
print(
    "Erweiterter Kombinatorik-Test: Prime-Pools bis verschiedene Grenzen"
)
print("─────────────────────────────────────────────────────────────")

your_group = [11, 17, 29, 41, 47, 59]

limits = [59, 97, 137, 211, 307]

for limit in limits:
    prime_pool = [p for p in generate_primes(limit) if p >= 11]

    combo_results = []

    for group in combinations(prime_pool, 6):
        if group[0] != 11:
            continue

        data = alpha_inv_eabc(42, 7, list(group))
        alpha_inv_model = data["alpha_inv"]

        abs_err = abs(alpha_inv_model - alpha_inv_target)
        rel_err = abs((1 / alpha_inv_model) - (1 / alpha_inv_target)) / (
            1 / alpha_inv_target
        )

        combo_results.append(
            {
                "group": list(group),
                "residues": [p % 12 for p in group],
                "alpha_inv": alpha_inv_model,
                "abs_err": abs_err,
                "rel_err": rel_err,
                **data,
            }
        )

    combo_results.sort(key=lambda x: x["abs_err"])

    your_rank = None
    your_entry = None

    for rank, r in enumerate(combo_results, start=1):
        if r["group"] == your_group:
            your_rank = rank
            your_entry = r
            break

    best = combo_results[0]

    print()
    print(f"Prime-Pool <= {limit}")
    print("-" * 40)
    print(f"Anzahl Gruppen       = {len(combo_results)}")
    print(f"Beste Gruppe         = {best['group']}")
    print(f"Beste mod12          = {best['residues']}")
    print(f"Beste alpha^-1       = {best['alpha_inv']:.12f}")
    print(f"Beste abs_err        = {best['abs_err']:.12e}")
    print(f"Beste D              = {best['D']:.12f}")
    print(f"Beste sqrtT          = {best['sqrtT']:.12f}")

    if your_rank is not None:
        print()
        print(
            f"Rang deiner Gruppe   = {your_rank} "
            f"von {len(combo_results)}"
        )
        print(f"Deine alpha^-1       = {your_entry['alpha_inv']:.12f}")
        print(f"Deine abs_err        = {your_entry['abs_err']:.12e}")
    else:
        print()
        print("Deine Gruppe ist in diesem Pool nicht enthalten.")

    print()
    print("Top 10:")
    for rank10, r in enumerate(combo_results[:10], start=1):
        print(
            f"rank={rank10} group={r['group']} "
            f"mod12={r['residues']} "
            f"alpha^-1={r['alpha_inv']:.12f} "
            f"abs_err={r['abs_err']:.12e}"
        )

print()
print("Fertig erweiterter Kombinatorik-Test.")

# 14. Ziel-Disorder-Test: Formel nach disorder invertieren (Bisektion)
print()
print("─────────────────────────────────────────────────────────────")
print("Ziel-Disorder-Test: Löse Formel nach disorder auf")
print("─────────────────────────────────────────────────────────────")

B = 137.0
M = 7.0


def alpha_inv_from_disorder(disorder):
    x = disorder**2
    T = 1 / (1 + x)
    D = x / M
    R = (1 - T) / B
    return B + D + math.sqrt(T) * R


lo = 0.0
hi = 2.0

for _ in range(200):
    mid = (lo + hi) / 2
    val = alpha_inv_from_disorder(mid)

    if val < alpha_inv_target:
        lo = mid
    else:
        hi = mid

target_disorder = (lo + hi) / 2
target_x = target_disorder**2
target_T = 1 / (1 + target_x)
target_D = target_x / M
target_R = (1 - target_T) / B
target_alpha_inv_model = alpha_inv_from_disorder(target_disorder)

print(f"Ziel alpha^-1          = {alpha_inv_target:.12f}")
print(f"Ziel disorder          = {target_disorder:.15f}")
print(f"Ziel disorder^2        = {target_x:.15f}")
print(f"Ziel T                 = {target_T:.15f}")
print(f"Ziel sqrt(T)           = {math.sqrt(target_T):.15f}")
print(f"Ziel D                 = {target_D:.15f}")
print(f"Ziel R                 = {target_R:.15f}")
print(f"Modell alpha^-1        = {target_alpha_inv_model:.15f}")
print()

print("Vergleich ausgewählter Gruppen gegen Ziel-Disorder")
print("--------------------------------------------------")

selected_groups = [
    [11, 17, 29, 41, 47, 59],
    [11, 23, 47, 59, 71, 73],
    [11, 41, 73, 83, 97, 113],
    [11, 127, 167, 191, 229, 263],
]

for group in selected_groups:
    data = alpha_inv_eabc(42, 7, group)
    disorder_group = data["disorder"]
    disorder_err = abs(disorder_group - target_disorder)
    alpha_err = abs(data["alpha_inv"] - alpha_inv_target)

    print(
        f"group={group} "
        f"disorder={disorder_group:.15f} "
        f"disorder_err={disorder_err:.15e} "
        f"alpha^-1={data['alpha_inv']:.12f} "
        f"alpha_err={alpha_err:.12e}"
    )

print()
print("Fertig Ziel-Disorder-Test.")
print("Fertig (Skript ohne Abbruch durchgelaufen).")