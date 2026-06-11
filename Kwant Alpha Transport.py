import numpy as np
import math
import os
import sys

try:
    import kwant
except ImportError:
    print("Kwant konnte nicht geladen werden (Modul fehlt oder lädt nicht).")
    print()
    print(
        "Hinweis: Unter Python 3.13 schlägt oft bereits »pip install kwant« "
        "fehl — das Projekt wird dann nicht erfolgreich gebaut "
        "(NumPy-/ABI-Anbindung)."
    )
    print()
    print("Typisch hilfreich:")
    print(
        "  • eigene Umgebung mit Python ≤ 3.12, dort »pip install kwant«;"
    )
    print(
        "  • oder conda/mamba: »conda install -n DEIN_ENV -c conda-forge kwant«"
    )
    print()
    print(
        "Wenn conda »already installed« meldet, liegt kwant oft nur in conda — "
        "dieses Skript dann mit Condas Python starten, z.B.:"
    )
    print(
        '  "$(conda info --base)/bin/python" "Kwant Alpha Transport.py"'
    )
    print()
    print(
        'Oder (wenn conda installiert ist): '
        'conda run -n base python "' + sys.argv[0] + '"'
    )
    print()
    print("Mit pyenv z.B.: pyenv install 3.12.x && pyenv local 3.12.x")
    raise SystemExit(1)


def _env_truthy(name: str) -> bool:
    v = os.environ.get(name, "")
    return v.lower() in ("1", "true", "yes", "on")


# Schnellmodus (weniger Rasterpunkte): z.B.
#   KWANT_TRANSPORT_QUICK=1 conda run -n base python "Kwant Alpha Transport.py"
QUICK = _env_truthy("KWANT_TRANSPORT_QUICK")


# ============================================================
# 1. Grunddaten
# ============================================================

alpha_inv_target = 137.035999084

N = 42
M = 7

channels = np.array([11, 17, 29, 41, 47, 59], dtype=float)
# Eine Zeile pro Kanalpotential — muss zur Rasterbreite L passen.
L = len(channels)

B = 3 * N + channels[0]

mean_p = np.mean(channels)
std_p = np.std(channels)
disorder = std_p / mean_p
x = disorder**2

T_disorder = 1 / (1 + x)

D = x / M
R = (1 - T_disorder) / B

alpha_inv_disorder = B + D + math.sqrt(T_disorder) * R

# Feinsweep: hier nur Energiebereich und Punktzahl ändern (Rest des Skripts unverändert).
FEINSWEEP_E_LO = -0.1
FEINSWEEP_E_HI = 0.1
FEINSWEEP_N = 400 if QUICK else 5000

# Raster für Kwant-Schleifen (voller Lauf kann sehr lange dauern).
COARSE_POINTS = 51 if QUICK else 401
ROBUST_GROB_POINTS = 81 if QUICK else 801
ROBUST_FINE_POINTS = 800 if QUICK else 8000
LINE_G_STEPS = 9 if QUICK else 17
LINE_GROB_POINTS = 81 if QUICK else 501
LINE_FINE_POINTS = 400 if QUICK else 2500

if QUICK:
    print("KWANT_TRANSPORT_QUICK=1 — reduzierte Raster für schnellen Testlauf.")
    print()


# ============================================================
# 2. Kwant-System bauen
# ============================================================

lat = kwant.lattice.square(norbs=1)

def make_system(channels, t=1.0, onsite_strength=1.0):
    syst = kwant.Builder()

    mean_p = np.mean(channels)
    std_p = np.std(channels)

    # normierte Kanalpotentiale
    potentials = onsite_strength * (channels - mean_p) / std_p

    # Streubereich: x = 0..M-1, y = 0..L-1
    for x_pos in range(M):
        for y_pos in range(L):
            syst[lat(x_pos, y_pos)] = potentials[y_pos]

    # Hoppings im Rechteck
    for x_pos in range(M):
        for y_pos in range(L):
            if x_pos + 1 < M:
                syst[lat(x_pos, y_pos), lat(x_pos + 1, y_pos)] = -t
            if y_pos + 1 < L:
                syst[lat(x_pos, y_pos), lat(x_pos, y_pos + 1)] = -t

    # Linker Lead
    sym_left = kwant.TranslationalSymmetry((-1, 0))
    lead_left = kwant.Builder(sym_left)

    for y_pos in range(L):
        lead_left[lat(0, y_pos)] = 0.0

    for y_pos in range(L):
        lead_left[lat(0, y_pos), lat(1, y_pos)] = -t
        if y_pos + 1 < L:
            lead_left[lat(0, y_pos), lat(0, y_pos + 1)] = -t

    # Rechter Lead
    sym_right = kwant.TranslationalSymmetry((1, 0))
    lead_right = kwant.Builder(sym_right)

    for y_pos in range(L):
        lead_right[lat(0, y_pos)] = 0.0

    for y_pos in range(L):
        lead_right[lat(0, y_pos), lat(1, y_pos)] = -t
        if y_pos + 1 < L:
            lead_right[lat(0, y_pos), lat(0, y_pos + 1)] = -t

    syst.attach_lead(lead_left)
    syst.attach_lead(lead_right)

    return syst.finalized()


# ============================================================
# 3. Kwant-Transmission berechnen
# ============================================================

def kwant_transmission(fsyst, energy):
    try:
        smat = kwant.smatrix(fsyst, energy)
        return smat.transmission(1, 0)
    except Exception:
        return np.nan


fsyst = make_system(channels, t=1.0, onsite_strength=1.0)

energies_coarse = np.linspace(-3.5, 3.5, COARSE_POINTS)

T_values = []

for E in energies_coarse:
    T_values.append(kwant_transmission(fsyst, E))

T_values = np.array(T_values, dtype=float)

valid = np.isfinite(T_values)

# Normierte Transmission: Kwant liefert absolute Kanaltransmission.
# Für Vergleich mit T_disorder normieren wir auf maximale Transmission.
T_max = np.nanmax(T_values)

if T_max > 0:
    T_norm = T_values / T_max
else:
    T_norm = T_values * np.nan


def transmission(E):
    return kwant_transmission(fsyst, E)


energies_fein = np.linspace(FEINSWEEP_E_LO, FEINSWEEP_E_HI, FEINSWEEP_N)

best_fein = None
alpha_target = alpha_inv_target

for E in energies_fein:
    T_raw = transmission(E)
    if not np.isfinite(T_raw) or not (T_max > 0):
        continue
    Tn = T_raw / T_max
    alpha_inv = B + D + np.sqrt(Tn) * (1 - Tn) / B
    err = abs(alpha_inv - alpha_target)
    if best_fein is None or err < best_fein[0]:
        best_fein = (err, E, Tn, T_raw, alpha_inv)

print()
print(
    "(err, E, T_norm, T_raw, alpha_inv) Feinsweep; "
    f"E in [{FEINSWEEP_E_LO}, {FEINSWEEP_E_HI}], {FEINSWEEP_N} Punkte;"
    " T_norm=T_kwant/T_max (T_max aus grobem Raster)"
)
print(best_fein)


# ============================================================
# 4. Suche Energie, bei der Kwant-T am besten zur Modell-T passt
# ============================================================

diff = np.abs(T_norm - T_disorder)
diff[~np.isfinite(diff)] = np.inf

best_idx = int(np.argmin(diff))
best_E = energies_coarse[best_idx]
T_kwant_best = T_norm[best_idx]
T_kwant_raw = T_values[best_idx]

# Kwant-basierte Alphaformel
R_kwant = (1 - T_kwant_best) / B

alpha_inv_kwant = B + D + math.sqrt(T_kwant_best) * R_kwant

alpha_err_kwant = abs(alpha_inv_kwant - alpha_inv_target)
alpha_err_disorder = abs(alpha_inv_disorder - alpha_inv_target)


# ============================================================
# 5. Ausgabe
# ============================================================

print("EABC/Kwant echter Transport-Test")
print("================================")
print()
print(f"N                         = {N}")
print(f"L x M                     = {L} x {M}")
print(f"channels                  = {channels.astype(int).tolist()}")
print(f"B = 3N+p0                 = {B:.12f}")
print()
print("Primkanal-Disorder:")
print("-------------------")
print(f"mean                      = {mean_p:.12f}")
print(f"std                       = {std_p:.12f}")
print(f"disorder=std/mean         = {disorder:.15f}")
print(f"x=disorder^2              = {x:.15f}")
print(f"D=x/M                     = {D:.15f}")
print(f"T_disorder                = {T_disorder:.15f}")
print(f"sqrt(T_disorder)          = {math.sqrt(T_disorder):.15f}")
print()
print("Alpha aus Disorder-T:")
print("---------------------")
print(f"alpha_inv_disorder        = {alpha_inv_disorder:.15f}")
print(f"target alpha_inv          = {alpha_inv_target:.15f}")
print(f"abs error                 = {alpha_err_disorder:.15e}")
print()
print("Kwant-Transport:")
print("----------------")
print(f"T_max raw                 = {T_max:.12f}")
print(f"beste Energie E           = {best_E:.12f}")
print(f"T_kwant_raw               = {T_kwant_raw:.12f}")
print(f"T_kwant_norm              = {T_kwant_best:.15f}")
print(f"|T_kwant_norm-T_disorder| = {abs(T_kwant_best - T_disorder):.15e}")
print()
print("Alpha aus Kwant-T:")
print("------------------")
print(f"alpha_inv_kwant           = {alpha_inv_kwant:.15f}")
print(f"target alpha_inv          = {alpha_inv_target:.15f}")
print(f"abs error                 = {alpha_err_kwant:.15e}")

print()
print("Top 20 Energien nach Nähe zu T_disorder:")
print("----------------------------------------")

order = np.argsort(diff)

for idx in order[:20]:
    if not np.isfinite(T_norm[idx]):
        continue

    E = energies_coarse[idx]
    Tk = T_norm[idx]
    Rk = (1 - Tk) / B
    alpha_inv_k = B + D + math.sqrt(Tk) * Rk
    alpha_err_k = abs(alpha_inv_k - alpha_inv_target)

    print(
        f"E={E: .6f}  "
        f"T_norm={Tk:.12f}  "
        f"alpha^-1={alpha_inv_k:.12f}  "
        f"alpha_err={alpha_err_k:.12e}"
    )

# ============================================================
# 6. Robustheitstest: Variation der Onsite-Stärke
# ============================================================

print()
print("─────────────────────────────────────────────────────────────")
print("Robustheitstest: Variation der Onsite-Stärke")
print("─────────────────────────────────────────────────────────────")

onsite_values = [
    0.25,
    0.50,
    0.75,
    0.90,
    1.00,
    1.10,
    1.25,
    1.50,
    2.00,
]

robust_results = []

for onsite_strength in onsite_values:

    fsyst_test = make_system(channels, t=1.0, onsite_strength=onsite_strength)

    energies_grob = np.linspace(-3.5, 3.5, ROBUST_GROB_POINTS)
    T_grob = []

    for E in energies_grob:
        T_grob.append(kwant_transmission(fsyst_test, E))

    T_grob = np.array(T_grob, dtype=float)
    T_max_test = np.nanmax(T_grob)

    if not np.isfinite(T_max_test) or T_max_test <= 0:
        continue

    energies_fine = np.linspace(-0.12, 0.12, ROBUST_FINE_POINTS)

    best_on = None

    for E in energies_fine:
        T_raw = kwant_transmission(fsyst_test, E)

        if not np.isfinite(T_raw):
            continue

        T_norm = T_raw / T_max_test

        if T_norm <= 0:
            continue

        alpha_inv_r = B + D + math.sqrt(T_norm) * (1 - T_norm) / B
        err = abs(alpha_inv_r - alpha_inv_target)

        if best_on is None or err < best_on[0]:
            best_on = (err, E, T_norm, T_raw, alpha_inv_r, T_max_test)

    if best_on is not None:
        robust_results.append((onsite_strength, *best_on))

print()
print("Onsite-Robustheit:")
print("------------------")

for (
    onsite_strength,
    err,
    E,
    T_norm,
    T_raw,
    alpha_inv,
    T_max_test,
) in robust_results:
    print(
        f"onsite={onsite_strength:5.2f}  "
        f"E_alpha={E: .12f}  "
        f"T_norm={T_norm:.12f}  "
        f"alpha^-1={alpha_inv:.12f}  "
        f"err={err:.12e}  "
        f"T_max={T_max_test:.12f}"
    )

print()
print("Fertig Onsite-Robustheitstest.")

# ============================================================
# 7. 2D-Resonanzlinien-Test: g und E variieren
# ============================================================

print()
print("─────────────────────────────────────────────────────────────")
print("2D-Resonanzlinien-Test: g in [0.8,1.2], E in [-0.15,0.15]")
print("─────────────────────────────────────────────────────────────")

g_values = np.linspace(0.8, 1.2, LINE_G_STEPS)

line_results = []

for g in g_values:
    print(f"Teste g = {g:.3f}")

    fsyst_g = make_system(channels, t=1.0, onsite_strength=g)

    energies_grob = np.linspace(-3.5, 3.5, LINE_GROB_POINTS)
    T_grob = []

    for E in energies_grob:
        T_grob.append(kwant_transmission(fsyst_g, E))

    T_grob = np.array(T_grob, dtype=float)
    T_max_g = np.nanmax(T_grob)

    if not np.isfinite(T_max_g) or T_max_g <= 0:
        continue

    energies_fine = np.linspace(-0.15, 0.15, LINE_FINE_POINTS)

    best_rl = None

    for E in energies_fine:
        T_raw = kwant_transmission(fsyst_g, E)

        if not np.isfinite(T_raw):
            continue

        T_norm = T_raw / T_max_g

        if T_norm <= 0:
            continue

        alpha_inv_line = B + D + math.sqrt(T_norm) * (1 - T_norm) / B
        err = abs(alpha_inv_line - alpha_inv_target)

        if best_rl is None or err < best_rl[0]:
            best_rl = (
                err,
                E,
                T_norm,
                T_raw,
                alpha_inv_line,
                T_max_g,
            )

    if best_rl is not None:
        line_results.append((g, *best_rl))

print()
print("Resonanzlinie E_alpha(g):")
print("-------------------------")

for (
    g,
    err,
    E,
    T_norm,
    T_raw,
    alpha_inv,
    T_max_g,
) in line_results:
    print(
        f"g={g:.3f}  "
        f"E_alpha={E: .12f}  "
        f"T_norm={T_norm:.12f}  "
        f"alpha^-1={alpha_inv:.12f}  "
        f"err={err:.12e}"
    )

print()
print("Beste Punkte insgesamt:")
print("-----------------------")

line_results_sorted = sorted(line_results, key=lambda x: x[1])

for (
    g,
    err,
    E,
    T_norm,
    T_raw,
    alpha_inv,
    T_max_g,
) in line_results_sorted[:10]:
    print(
        f"g={g:.3f}  "
        f"E_alpha={E: .12f}  "
        f"T_norm={T_norm:.12f}  "
        f"alpha^-1={alpha_inv:.12f}  "
        f"err={err:.12e}  "
        f"T_max={T_max_g:.12f}"
    )

print()
print("Fertig 2D-Resonanzlinien-Test.")
print("Fertig.")