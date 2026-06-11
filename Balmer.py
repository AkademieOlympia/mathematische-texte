import numpy as np
import pandas as pd

# ============================================================
# 1. Daten
# ============================================================

zeros = np.load("zeros6.npy").astype(float).ravel()
zeros = np.sort(zeros[np.isfinite(zeros)])
z = zeros[:500000]

WINDOW_LEN = 60
STEP = 10
K = 6
N_MC = 200

def zscore(x):
    x = np.asarray(x, dtype=float)
    s = np.std(x)
    if s == 0:
        return None
    return (x - np.mean(x)) / s

# ============================================================
# 2. Wasserstoff-Referenzen
# ============================================================

# Balmer: n1 = 2, n2 = 3..8
n_balmer = np.arange(3, 9, dtype=float)
balmer = 1/4 - 1/(n_balmer**2)
balmer_z = zscore(balmer)

# Lyman: n1 = 1, n2 = 2..7
n_lyman = np.arange(2, 8, dtype=float)
lyman = 1 - 1/(n_lyman**2)
lyman_z = zscore(lyman)

# ============================================================
# 3. Modi
# ============================================================

def build_modes(window):
    window = np.asarray(window, dtype=float)
    gaps = np.diff(window)
    cumgaps = np.cumsum(gaps)
    return {
        "raw": window,
        "gap": gaps,
        "cumgap": cumgaps,
    }

def transform_log(x):
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0):
        return None
    return np.log(x)

def rms(a, b):
    return float(np.sqrt(np.mean((a - b)**2)))

# ============================================================
# 4. Bester Treffer im Spektrum
# ============================================================

def best_match_in_spectrum(zeros_1d):
    best_balmer = None
    best_lyman = None

    for start in range(0, len(zeros_1d) - WINDOW_LEN + 1, STEP):
        window = zeros_1d[start:start + WINDOW_LEN]
        modes = build_modes(window)

        for mode_name, arr in modes.items():
            y_all = transform_log(arr)
            if y_all is None or len(y_all) < K:
                continue

            for j in range(0, len(y_all) - K + 1):
                y = y_all[j:j+K]
                yz = zscore(y)
                if yz is None:
                    continue

                e_balmer = rms(yz, balmer_z)
                e_lyman = rms(yz, lyman_z)

                row_balmer = {
                    "window_start": start,
                    "segment_start": start + j,
                    "segment_end": start + j + K - 1,
                    "mode": mode_name,
                    "error": e_balmer,
                    "values": y.tolist(),
                }

                row_lyman = {
                    "window_start": start,
                    "segment_start": start + j,
                    "segment_end": start + j + K - 1,
                    "mode": mode_name,
                    "error": e_lyman,
                    "values": y.tolist(),
                }

                if best_balmer is None or e_balmer < best_balmer["error"]:
                    best_balmer = row_balmer

                if best_lyman is None or e_lyman < best_lyman["error"]:
                    best_lyman = row_lyman

    return best_balmer, best_lyman

# ============================================================
# 5. Gap-shuffle Kontrolle
# ============================================================

def gap_shuffle_control(zeros_1d, rng):
    gaps = np.diff(zeros_1d)
    gaps_perm = rng.permutation(gaps)
    out = np.empty_like(zeros_1d)
    out[0] = zeros_1d[0]
    out[1:] = out[0] + np.cumsum(gaps_perm)
    return np.sort(out)

# ============================================================
# 6. Echter Datensatz
# ============================================================

best_balmer, best_lyman = best_match_in_spectrum(z)

print("Bester Balmer-Treffer:")
print(best_balmer)

print("\nBester Lyman-Treffer:")
print(best_lyman)

# ============================================================
# 7. Monte-Carlo
# ============================================================

rng = np.random.default_rng(12345)
mc_balmer = []
mc_lyman = []

for _ in range(N_MC):
    zc = gap_shuffle_control(z, rng)
    b, l = best_match_in_spectrum(zc)
    mc_balmer.append(b["error"])
    mc_lyman.append(l["error"])

mc_balmer = np.array(mc_balmer)
mc_lyman = np.array(mc_lyman)

p_balmer = (np.sum(mc_balmer <= best_balmer["error"]) + 1) / (N_MC + 1)
p_lyman = (np.sum(mc_lyman <= best_lyman["error"]) + 1) / (N_MC + 1)

print("\nBalmer p-Wert:", p_balmer)
print("Lyman p-Wert:", p_lyman)

# ============================================================
# 8. Export
# ============================================================

summary = pd.DataFrame([{
    "balmer_error": best_balmer["error"],
    "balmer_mode": best_balmer["mode"],
    "balmer_segment_start": best_balmer["segment_start"],
    "balmer_p_value": p_balmer,
    "lyman_error": best_lyman["error"],
    "lyman_mode": best_lyman["mode"],
    "lyman_segment_start": best_lyman["segment_start"],
    "lyman_p_value": p_lyman,
}])

summary.to_csv("balmer_lyman_nullstellen_test_summary.csv", index=False)
pd.DataFrame({"mc_balmer_error": mc_balmer, "mc_lyman_error": mc_lyman}).to_csv(
    "balmer_lyman_nullstellen_test_mc.csv", index=False
)