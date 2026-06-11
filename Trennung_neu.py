import sympy as sp
import numpy as np
import pandas as pd
import itertools as it
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt

# =====================================================
# 1. Intervall
# =====================================================

A = 113000
B = 113300
WINDOW_SIZE = B - A

# Skalierungstests in gleich breiten Fenstern
SCALE_BASES = [10**6, 10**7, 10**8]

# Zusatztest: breitere Fenster
SCALE_WIDTHS = [WINDOW_SIZE, 3000, 10000]


# =====================================================
# 2. Deine Spektralsumme S(M) hier einsetzen
# =====================================================

def S(M):
    """
    Ersatzmodell fuer die Spektralsumme, damit das Skript lauffaehig ist.
    Wenn du deine echte Formel hast, diese Funktion direkt ersetzen.
    """
    f = sp.factorint(M)
    omega = len(f)
    Omega = sum(f.values())
    tau = sp.divisor_count(M)
    return float(np.log1p(M) - 0.75 * omega - 0.35 * Omega + 0.02 * tau)


def build_dataframe(interval_start, interval_end):
    rows = []

    for M in range(interval_start, interval_end + 1):
        try:
            s_val = S(M)
        except Exception as e:
            print(f"Fehler bei M={M}: {e}")
            continue

        is_prime = sp.isprime(M)
        f = sp.factorint(M)
        omega = len(f)
        Omega = sum(f.values())
        tau = sp.divisor_count(M)

        rows.append(
            {
                "M": M,
                "S": float(s_val),
                "prime": int(is_prime),
                "composite": int((M > 1) and not is_prime),
                "omega": omega,
                "Omega": Omega,
                "tau": int(tau),
                "logM": float(np.log1p(M)),
                "factorization": f,
            }
        )

    df_local = pd.DataFrame(rows)
    if df_local.empty:
        raise RuntimeError("Keine gueltigen Daten erzeugt. Pruefe S(M) und Intervall.")

    primes_local = df_local[df_local["prime"] == 1]
    comps_local = df_local[df_local["composite"] == 1]
    if len(primes_local) == 0 or len(comps_local) == 0:
        raise RuntimeError("Zu wenige Klassen fuer Statistik (Primzahlen oder Komposita fehlen).")

    return df_local


def linear_r2(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot == 0:
        return np.nan
    return 1 - ss_res / ss_tot


def fit_linear_projection(df_local, feature_cols):
    X = df_local[feature_cols].to_numpy(dtype=float)
    y = df_local["S"].to_numpy(dtype=float)
    X_design = np.column_stack([np.ones(len(X)), X])
    beta, _, _, _ = np.linalg.lstsq(X_design, y, rcond=None)
    y_hat = X_design @ beta
    residual = y - y_hat
    return beta, y_hat, residual


def safe_auc(y_true, score):
    return roc_auc_score(y_true, score)


def orientation_free_auc(y_true, score):
    auc_raw = safe_auc(y_true, score)
    return max(auc_raw, 1 - auc_raw)


def add_sieve_features(df_in):
    """Erzeuge lokale Sieb-Features auf einem DataFrame-Fenster."""
    df_out = df_in.sort_values("M").copy()
    M_arr = df_out["M"].to_numpy(dtype=int)
    prime_arr = df_out[df_out["prime"] == 1]["M"].to_numpy(dtype=int)

    if len(prime_arr) == 0:
        # Fallback: neutrale Features
        df_out["nearest_prime_dist"] = 0.0
        df_out["local_gap_proxy"] = 0.0
    else:
        pos = np.searchsorted(prime_arr, M_arr)
        prev_p = np.where(pos > 0, prime_arr[pos - 1], prime_arr[0])
        next_idx = np.clip(pos, 0, len(prime_arr) - 1)
        next_p = prime_arr[next_idx]

        dist_prev = np.abs(M_arr - prev_p).astype(float)
        dist_next = np.abs(next_p - M_arr).astype(float)
        df_out["nearest_prime_dist"] = np.minimum(dist_prev, dist_next)
        df_out["local_gap_proxy"] = dist_prev + dist_next

    mod60 = (M_arr % 60).astype(float)
    df_out["sin60"] = np.sin(2 * np.pi * mod60 / 60.0)
    df_out["cos60"] = np.cos(2 * np.pi * mod60 / 60.0)
    df_out["coprime60"] = np.array([1.0 if np.gcd(int(m), 60) == 1 else 0.0 for m in mod60], dtype=float)

    return df_out


# =====================================================
# 3. Daten erzeugen
# =====================================================

df = build_dataframe(A, B)
primes = df[df["prime"] == 1]["S"]
comps = df[df["composite"] == 1]["S"]

# =====================================================
# 4. Basisstatistik
# =====================================================

print("\n=== Basisstatistik ===")
print("Anzahl Primzahlen:", len(primes))
print("Anzahl Zusammengesetzte:", len(comps))

print("\nPrimzahlen:")
print(primes.describe())

print("\nZusammengesetzte:")
print(comps.describe())

mu_c = comps.mean()
sigma_c = comps.std(ddof=1)
mu_p = primes.mean()
sigma_p = primes.std(ddof=1)

print("\n=== Abstände ===")
print("Mittelwert zusammengesetzt:", mu_c)
print("Mittelwert prim:", mu_p)
print("Differenz comp - prime:", mu_c - mu_p)
print("Z-Abstand der Primzahlen relativ zur Komposit-Verteilung:", (mu_p - mu_c) / sigma_c)

# =====================================================
# 5. Signifikanztests
# =====================================================

print("\n=== Signifikanztests ===")
welch = stats.ttest_ind(comps, primes, equal_var=False)
mann = stats.mannwhitneyu(comps, primes, alternative="two-sided")
ks = stats.ks_2samp(comps, primes)
print("Welch-t-Test:", welch)
print("Mann-Whitney-U:", mann)
print("Kolmogorov-Smirnov:", ks)

pooled = np.sqrt(
    ((len(comps) - 1) * sigma_c**2 + (len(primes) - 1) * sigma_p**2)
    / (len(comps) + len(primes) - 2)
)
cohen_d = (mu_c - mu_p) / pooled
print("Cohen d:", cohen_d)

# =====================================================
# 6. Klassifikationstest: Kann S(M) Primzahlen erkennen?
# =====================================================

y_true = df["prime"].values
score = df["S"].values
auc = safe_auc(y_true, score)

print("\n=== ROC-Test ===")
print("AUC:", auc)

fpr, tpr, thresholds = roc_curve(y_true, score)
youden = tpr - fpr
best_idx = np.argmax(youden)
best_threshold = thresholds[best_idx]

print("Bester Schwellenwert für S:", best_threshold)
print("entspricht S >", best_threshold)
print("TPR:", tpr[best_idx])
print("FPR:", fpr[best_idx])

# =====================================================
# 6b. Strukturtest: Korrelationen
# =====================================================

print("\n=== Korrelationsmatrix S/omega/Omega/tau ===")
print(df[["S", "omega", "Omega", "tau"]].corr())

# =====================================================
# 6c. Ablation + Residualtest
# =====================================================

print("\n=== Ablation (lineare Projektion auf Strukturfeatures) ===")
ablations = [["omega"], ["Omega"], ["tau"], ["omega", "Omega"], ["omega", "Omega", "tau"]]

for cols in ablations:
    _, y_hat, resid = fit_linear_projection(df, cols)
    r2 = linear_r2(df["S"].to_numpy(dtype=float), y_hat)
    auc_resid = orientation_free_auc(y_true, resid)
    print(f"Features={cols} -> R^2={r2:.6f}, AUC(|Richtung frei|) auf Residuum={auc_resid:.6f}")

full_cols = ["omega", "Omega", "tau", "logM"]
_, y_hat_full, resid_full = fit_linear_projection(df, full_cols)
r2_full = linear_r2(df["S"].to_numpy(dtype=float), y_hat_full)
auc_resid_full_raw = safe_auc(y_true, resid_full)
auc_resid_full = orientation_free_auc(y_true, resid_full)

print("\nVollmodell Features", full_cols)
print("R^2 Vollmodell:", r2_full)
print("AUC Residuum (roh):", auc_resid_full_raw)
print("AUC Residuum (richtungfrei):", auc_resid_full)

# =====================================================
# 6d. Carmichael-Schnelltest
# =====================================================

carmichael = [561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841]
print("\n=== Carmichael-Test ===")
for n in carmichael:
    try:
        print(n, S(n))
    except Exception as e:
        print(n, f"Fehler: {e}")

# =====================================================
# 6e. Skalierungstest
# =====================================================

print("\n=== Skalierungstest (mehrere Fensterbreiten) ===")
scale_rows = []

for width in SCALE_WIDTHS:
    for base in SCALE_BASES:
        a_s = base
        b_s = base + width

        try:
            df_s = build_dataframe(a_s, b_s)
        except Exception as e:
            print(f"Fenster [{a_s}, {b_s}] fehlgeschlagen: {e}")
            continue

        y_s = df_s["prime"].to_numpy(dtype=int)
        auc_s = safe_auc(y_s, df_s["S"].to_numpy(dtype=float))
        _, y_hat_s, resid_s = fit_linear_projection(df_s, full_cols)
        r2_s = linear_r2(df_s["S"].to_numpy(dtype=float), y_hat_s)
        auc_resid_s_raw = safe_auc(y_s, resid_s)
        auc_resid_s = orientation_free_auc(y_s, resid_s)

        scale_rows.append(
            {
                "A": a_s,
                "B": b_s,
                "width": width,
                "n": len(df_s),
                "primes": int(df_s["prime"].sum()),
                "auc_S": auc_s,
                "r2_struct": r2_s,
                "auc_resid_raw": auc_resid_s_raw,
                "auc_resid_oriented": auc_resid_s,
            }
        )

scale_df = pd.DataFrame(scale_rows)
if not scale_df.empty:
    print(scale_df.to_string(index=False))
else:
    print("Keine Skalierungsdaten berechnet.")

# =====================================================
# 6f. Residualanalyse (Goldbereich)
# =====================================================

print("\n=== Residualanalyse ===")

df["R"] = resid_full

# A) Fourieranalyse des Residuums (nach M sortiert)
R_series = df.sort_values("M")["R"].to_numpy(dtype=float)
R_centered = R_series - R_series.mean()

fft_vals = np.fft.rfft(R_centered)
fft_freq = np.fft.rfftfreq(len(R_centered), d=1.0)
fft_amp = np.abs(fft_vals)

if len(fft_amp) > 1:
    peak_idx = int(np.argmax(fft_amp[1:]) + 1)
    print("Dominante FFT-Frequenz:", fft_freq[peak_idx])
    print("Dominante FFT-Amplitude:", fft_amp[peak_idx])
else:
    print("FFT zu kurz für Peak-Analyse.")

# B) Modulo-Klassen (12 und 30)
mod12_mean = df.groupby(df["M"] % 12)["R"].mean().sort_index()
mod30_mean = df.groupby(df["M"] % 30)["R"].mean().sort_index()

print("\nResiduum-Mittel nach mod 12:")
print(mod12_mean.to_string())

print("\nResiduum-Mittel nach mod 30 (nur nichtleere Klassen):")
print(mod30_mean.to_string())

# C) Primabstände vs. Residuum auf Primzahlen
prime_df = df[df["prime"] == 1].sort_values("M").copy()
prime_vals = prime_df["M"].to_numpy(dtype=int)
prime_res = prime_df["R"].to_numpy(dtype=float)

gap_corr = np.nan
if len(prime_vals) >= 3:
    gaps = np.diff(prime_vals)
    # Residuum des "rechten" Prims je Gap
    res_right = prime_res[1:]
    if np.std(gaps) > 0 and np.std(res_right) > 0:
        gap_corr = float(np.corrcoef(gaps, res_right)[0, 1])

print("\nKorrelation Primabstand vs Residuum(rechtes Prim):", gap_corr)

# D) Autokorrelation des Residuums
max_lag_res = min(80, len(R_centered) - 1)
auto = []
lags = []
for lag in range(1, max_lag_res + 1):
    x = R_centered[:-lag]
    y = R_centered[lag:]
    den = np.std(x) * np.std(y)
    if den == 0:
        auto.append(np.nan)
    else:
        auto.append(float(np.mean((x - x.mean()) * (y - y.mean())) / den))
    lags.append(lag)

auto = np.array(auto, dtype=float)
if np.isfinite(auto).any():
    idx_max = int(np.nanargmax(np.abs(auto)))
    print("Stärkste |Autokorrelation| bei Lag:", lags[idx_max], "Wert:", auto[idx_max])
else:
    print("Autokorrelation nicht bestimmbar.")

# =====================================================
# 6g. Signifikanztest des Residuums (Permutation/Bootstrap)
# =====================================================

print("\n=== Signifikanztest Residuum ===")

N_PERM = 400
N_BOOT = 400
rng = np.random.default_rng(42)


def max_abs_acf(arr, max_lag):
    arr = np.asarray(arr, dtype=float)
    arr_c = arr - arr.mean()
    vals = []
    for lag in range(1, max_lag + 1):
        x = arr_c[:-lag]
        y = arr_c[lag:]
        den = np.std(x) * np.std(y)
        if den == 0:
            vals.append(np.nan)
        else:
            vals.append(float(np.mean((x - x.mean()) * (y - y.mean())) / den))
    vals = np.array(vals, dtype=float)
    if not np.isfinite(vals).any():
        return np.nan
    return float(np.nanmax(np.abs(vals)))


def block_shuffle(arr, block_len, rng_local):
    """Shuffle contiguous blocks, preserve within-block structure."""
    arr = np.asarray(arr)
    n = len(arr)
    if n == 0:
        return arr.copy()
    blocks = [arr[i:i + block_len] for i in range(0, n, block_len)]
    order = rng_local.permutation(len(blocks))
    return np.concatenate([blocks[i] for i in order])


def moving_block_bootstrap(arr, block_len, rng_local):
    """Moving-block bootstrap sample with replacement."""
    arr = np.asarray(arr)
    n = len(arr)
    if n == 0:
        return arr.copy()
    starts = np.arange(0, max(1, n - block_len + 1))
    out = []
    while sum(len(x) for x in out) < n:
        s = int(rng_local.choice(starts))
        out.append(arr[s:s + block_len])
    return np.concatenate(out)[:n]


# Beobachtete Kennzahlen
obs_gap_corr = float(np.abs(gap_corr)) if np.isfinite(gap_corr) else np.nan
obs_acf_peak = max_abs_acf(R_centered, max_lag_res)

# Permutationstest Gap-Korrelation
perm_gap = []
if len(prime_vals) >= 3 and np.std(np.diff(prime_vals)) > 0 and np.std(prime_res[1:]) > 0:
    gaps = np.diff(prime_vals)
    res_right = prime_res[1:]
    for _ in range(N_PERM):
        rp = rng.permutation(res_right)
        c = np.corrcoef(gaps, rp)[0, 1]
        perm_gap.append(abs(float(c)))
perm_gap = np.array(perm_gap, dtype=float)

if len(perm_gap) > 0 and np.isfinite(obs_gap_corr):
    p_gap = (1 + np.sum(perm_gap >= obs_gap_corr)) / (len(perm_gap) + 1)
    print("Gap-Korrelation |obs|:", obs_gap_corr)
    print("Gap-Permutation p-Wert:", p_gap)
else:
    p_gap = np.nan
    print("Gap-Permutation nicht bestimmbar.")

# Permutationstest ACF-Peak
perm_acf = []
if len(R_centered) > 10 and np.isfinite(obs_acf_peak):
    for _ in range(N_PERM):
        rp = rng.permutation(R_centered)
        perm_acf.append(max_abs_acf(rp, max_lag_res))
perm_acf = np.array(perm_acf, dtype=float)

if len(perm_acf) > 0 and np.isfinite(obs_acf_peak):
    p_acf = (1 + np.sum(perm_acf >= obs_acf_peak)) / (len(perm_acf) + 1)
    print("ACF-Peak |obs|:", obs_acf_peak)
    print("ACF-Permutation p-Wert:", p_acf)
else:
    p_acf = np.nan
    print("ACF-Permutation nicht bestimmbar.")

# Konservativer Block-basierter Nulltest (Abhaengigkeiten erhalten)
block_len_gap = max(4, int(round((len(prime_res[1:]) if len(prime_res) > 1 else 1) ** (1 / 3))))
block_len_r = max(6, int(round(len(R_centered) ** (1 / 3))))

block_null_gap = []
if len(prime_vals) >= 3 and np.std(np.diff(prime_vals)) > 0 and np.std(prime_res[1:]) > 0:
    gaps = np.diff(prime_vals)
    res_right = prime_res[1:]
    for _ in range(N_PERM):
        rr = block_shuffle(res_right, block_len_gap, rng)
        c = np.corrcoef(gaps, rr)[0, 1]
        block_null_gap.append(abs(float(c)))
block_null_gap = np.array(block_null_gap, dtype=float)

if len(block_null_gap) > 0 and np.isfinite(obs_gap_corr):
    p_gap_block = (1 + np.sum(block_null_gap >= obs_gap_corr)) / (len(block_null_gap) + 1)
    print("Gap-Block-Permutation p-Wert:", p_gap_block, "(Blocklaenge:", block_len_gap, ")")
else:
    p_gap_block = np.nan
    print("Gap-Block-Permutation nicht bestimmbar.")

block_null_acf = []
if len(R_centered) > 10 and np.isfinite(obs_acf_peak):
    for _ in range(N_PERM):
        rr = block_shuffle(R_centered, block_len_r, rng)
        block_null_acf.append(max_abs_acf(rr, max_lag_res))
block_null_acf = np.array(block_null_acf, dtype=float)

if len(block_null_acf) > 0 and np.isfinite(obs_acf_peak):
    p_acf_block = (1 + np.sum(block_null_acf >= obs_acf_peak)) / (len(block_null_acf) + 1)
    print("ACF-Block-Permutation p-Wert:", p_acf_block, "(Blocklaenge:", block_len_r, ")")
else:
    p_acf_block = np.nan
    print("ACF-Block-Permutation nicht bestimmbar.")

# Bootstrap-Konfidenzintervalle
boot_gap = []
if len(prime_vals) >= 8:
    gaps = np.diff(prime_vals)
    res_right = prime_res[1:]
    n_pairs = len(gaps)
    for _ in range(N_BOOT):
        idx_b = rng.integers(0, n_pairs, n_pairs)
        gb = gaps[idx_b]
        rb = res_right[idx_b]
        if np.std(gb) > 0 and np.std(rb) > 0:
            boot_gap.append(np.corrcoef(gb, rb)[0, 1])
boot_gap = np.array(boot_gap, dtype=float)

if len(boot_gap) > 10:
    ci_gap = np.quantile(boot_gap, [0.025, 0.975])
    print("Gap-Korrelation 95%-Bootstrap-CI:", ci_gap)
else:
    ci_gap = np.array([np.nan, np.nan], dtype=float)
    print("Gap-Bootstrap-CI nicht bestimmbar.")

# Block-Bootstrap-CI fuer Gap-Korrelation
boot_gap_block = []
if len(prime_vals) >= 8:
    gaps = np.diff(prime_vals)
    res_right = prime_res[1:]
    n_pairs = len(gaps)
    pair_mat = np.column_stack([gaps, res_right])
    for _ in range(N_BOOT):
        sample = moving_block_bootstrap(pair_mat, block_len_gap, rng)
        gb = sample[:, 0]
        rb = sample[:, 1]
        if np.std(gb) > 0 and np.std(rb) > 0:
            boot_gap_block.append(np.corrcoef(gb, rb)[0, 1])
boot_gap_block = np.array(boot_gap_block, dtype=float)

if len(boot_gap_block) > 10:
    ci_gap_block = np.quantile(boot_gap_block, [0.025, 0.975])
    print("Gap-Korrelation 95%-Block-Bootstrap-CI:", ci_gap_block)
else:
    ci_gap_block = np.array([np.nan, np.nan], dtype=float)
    print("Gap-Block-Bootstrap-CI nicht bestimmbar.")

boot_acf = []
if len(R_centered) > 10:
    n_r = len(R_centered)
    for _ in range(N_BOOT):
        idx_b = rng.integers(0, n_r, n_r)
        rb = R_centered[idx_b]
        boot_acf.append(max_abs_acf(rb, max_lag_res))
boot_acf = np.array(boot_acf, dtype=float)

if len(boot_acf) > 10:
    ci_acf = np.quantile(boot_acf, [0.025, 0.975])
    print("ACF-Peak 95%-Bootstrap-CI:", ci_acf)
else:
    ci_acf = np.array([np.nan, np.nan], dtype=float)
    print("ACF-Bootstrap-CI nicht bestimmbar.")

# Block-Bootstrap-CI fuer ACF-Peak
boot_acf_block = []
if len(R_centered) > 10:
    for _ in range(N_BOOT):
        rb = moving_block_bootstrap(R_centered, block_len_r, rng)
        boot_acf_block.append(max_abs_acf(rb, max_lag_res))
boot_acf_block = np.array(boot_acf_block, dtype=float)

if len(boot_acf_block) > 10:
    ci_acf_block = np.quantile(boot_acf_block, [0.025, 0.975])
    print("ACF-Peak 95%-Block-Bootstrap-CI:", ci_acf_block)
else:
    ci_acf_block = np.array([np.nan, np.nan], dtype=float)
    print("ACF-Block-Bootstrap-CI nicht bestimmbar.")

# =====================================================
# 6h. Explizite Zerlegung: E_mult + E_sieve
# =====================================================

print("\n=== Zerlegung E_mult + E_sieve ===")

df_sorted = add_sieve_features(df)

mult_cols = ["omega", "Omega", "tau", "logM"]
sieve_cols = ["nearest_prime_dist", "local_gap_proxy", "sin60", "cos60", "coprime60"]

_, yhat_mult, resid_mult = fit_linear_projection(df_sorted, mult_cols)
r2_mult = linear_r2(df_sorted["S"].to_numpy(dtype=float), yhat_mult)
auc_resid_mult = orientation_free_auc(df_sorted["prime"].to_numpy(dtype=int), resid_mult)

_, yhat_both, resid_both = fit_linear_projection(df_sorted, mult_cols + sieve_cols)
r2_both = linear_r2(df_sorted["S"].to_numpy(dtype=float), yhat_both)
auc_resid_both = orientation_free_auc(df_sorted["prime"].to_numpy(dtype=int), resid_both)

print("R^2(E_mult):", r2_mult)
print("AUC Residuum nach E_mult:", auc_resid_mult)
print("R^2(E_mult + E_sieve):", r2_both)
print("AUC Residuum nach E_mult+E_sieve:", auc_resid_both)
print("Delta R^2:", r2_both - r2_mult)
print("Delta AUC Residuum:", auc_resid_both - auc_resid_mult)

# =====================================================
# 6i. Generalisierung: train/test über Fenster
# =====================================================

print("\n=== Generalisierungstest (Train/Test über Fenster) ===")

train_base = 10**6
test_bases = [10**7, 10**8]
general_rows = []


def train_linear_beta(df_train, cols):
    X = df_train[cols].to_numpy(dtype=float)
    y = df_train["S"].to_numpy(dtype=float)
    Xd = np.column_stack([np.ones(len(X)), X])
    beta, _, _, _ = np.linalg.lstsq(Xd, y, rcond=None)
    return beta


def predict_linear_beta(df_any, cols, beta):
    X = df_any[cols].to_numpy(dtype=float)
    Xd = np.column_stack([np.ones(len(X)), X])
    return Xd @ beta


for width in [WINDOW_SIZE, 3000, 10000]:
    # Trainfenster
    df_train = build_dataframe(train_base, train_base + width)
    df_train = add_sieve_features(df_train)

    beta_mult = train_linear_beta(df_train, mult_cols)
    beta_both = train_linear_beta(df_train, mult_cols + sieve_cols)

    # Test auf anderen Skalen
    for base_t in test_bases:
        df_test = build_dataframe(base_t, base_t + width)
        df_test = add_sieve_features(df_test)
        y_test = df_test["S"].to_numpy(dtype=float)
        labels_test = df_test["prime"].to_numpy(dtype=int)

        yhat_mult_t = predict_linear_beta(df_test, mult_cols, beta_mult)
        yhat_both_t = predict_linear_beta(df_test, mult_cols + sieve_cols, beta_both)

        resid_mult_t = y_test - yhat_mult_t
        resid_both_t = y_test - yhat_both_t

        auc_mult_t = orientation_free_auc(labels_test, resid_mult_t)
        auc_both_t = orientation_free_auc(labels_test, resid_both_t)

        r2_mult_t = linear_r2(y_test, yhat_mult_t)
        r2_both_t = linear_r2(y_test, yhat_both_t)

        general_rows.append(
            {
                "train_A": train_base,
                "test_A": base_t,
                "width": width,
                "n_test": len(df_test),
                "primes_test": int(labels_test.sum()),
                "r2_mult_test": r2_mult_t,
                "r2_mult_sieve_test": r2_both_t,
                "auc_resid_mult_test": auc_mult_t,
                "auc_resid_mult_sieve_test": auc_both_t,
            }
        )

general_df = pd.DataFrame(general_rows)
if not general_df.empty:
    print(general_df.to_string(index=False))
else:
    print("Keine Generalisierungsdaten berechnet.")

# =====================================================
# 6j. Strenger OOS-Test (train-fixierte Standardisierung)
# =====================================================

print("\n=== Strenger OOS-Test (fixierte Skalierung + feste Richtung) ===")


def fit_standardized_beta(df_train, cols):
    X = df_train[cols].to_numpy(dtype=float)
    y = df_train["S"].to_numpy(dtype=float)
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd_safe = np.where(sd == 0, 1.0, sd)
    Xz = (X - mu) / sd_safe
    Xd = np.column_stack([np.ones(len(Xz)), Xz])
    beta, _, _, _ = np.linalg.lstsq(Xd, y, rcond=None)
    return beta, mu, sd_safe


def predict_standardized_beta(df_any, cols, beta, mu, sd_safe):
    X = df_any[cols].to_numpy(dtype=float)
    Xz = (X - mu) / sd_safe
    Xd = np.column_stack([np.ones(len(Xz)), Xz])
    return Xd @ beta


def signed_auc(labels, residual, train_direction):
    # Keine nachträgliche Richtungswahl im Test.
    score = train_direction * residual
    return safe_auc(labels, score)


strict_rows = []

for width in [WINDOW_SIZE, 3000, 10000]:
    df_train = build_dataframe(train_base, train_base + width)
    df_train = add_sieve_features(df_train)

    # Robuste, skalenärmere Siebfeatures
    eps = 1e-9
    df_train["npd_scaled"] = df_train["nearest_prime_dist"] / (df_train["logM"] + eps)
    df_train["lgp_scaled"] = df_train["local_gap_proxy"] / (df_train["logM"] + eps)

    mult_cols_strict = ["omega", "Omega", "tau", "logM"]
    sieve_cols_strict = ["npd_scaled", "lgp_scaled", "sin60", "cos60", "coprime60"]

    # Train auf standardized Features
    beta_m, mu_m, sd_m = fit_standardized_beta(df_train, mult_cols_strict)
    beta_b, mu_b, sd_b = fit_standardized_beta(df_train, mult_cols_strict + sieve_cols_strict)

    y_tr = df_train["S"].to_numpy(dtype=float)
    lbl_tr = df_train["prime"].to_numpy(dtype=int)
    pred_m_tr = predict_standardized_beta(df_train, mult_cols_strict, beta_m, mu_m, sd_m)
    pred_b_tr = predict_standardized_beta(df_train, mult_cols_strict + sieve_cols_strict, beta_b, mu_b, sd_b)

    resid_m_tr = y_tr - pred_m_tr
    resid_b_tr = y_tr - pred_b_tr

    # Richtungsfixierung nur aus TRAIN
    dir_m = 1.0 if resid_m_tr[lbl_tr == 1].mean() >= resid_m_tr[lbl_tr == 0].mean() else -1.0
    dir_b = 1.0 if resid_b_tr[lbl_tr == 1].mean() >= resid_b_tr[lbl_tr == 0].mean() else -1.0

    for base_t in test_bases:
        df_test = build_dataframe(base_t, base_t + width)
        df_test = add_sieve_features(df_test)
        df_test["npd_scaled"] = df_test["nearest_prime_dist"] / (df_test["logM"] + eps)
        df_test["lgp_scaled"] = df_test["local_gap_proxy"] / (df_test["logM"] + eps)

        y_te = df_test["S"].to_numpy(dtype=float)
        lbl_te = df_test["prime"].to_numpy(dtype=int)

        pred_m_te = predict_standardized_beta(df_test, mult_cols_strict, beta_m, mu_m, sd_m)
        pred_b_te = predict_standardized_beta(df_test, mult_cols_strict + sieve_cols_strict, beta_b, mu_b, sd_b)

        resid_m_te = y_te - pred_m_te
        resid_b_te = y_te - pred_b_te

        auc_m_signed = signed_auc(lbl_te, resid_m_te, dir_m)
        auc_b_signed = signed_auc(lbl_te, resid_b_te, dir_b)

        strict_rows.append(
            {
                "train_A": train_base,
                "test_A": base_t,
                "width": width,
                "n_test": len(df_test),
                "primes_test": int(lbl_te.sum()),
                "auc_resid_mult_signed": auc_m_signed,
                "auc_resid_mult_sieve_signed": auc_b_signed,
                "delta_auc_signed": auc_b_signed - auc_m_signed,
            }
        )

strict_df = pd.DataFrame(strict_rows)
if not strict_df.empty:
    print(strict_df.to_string(index=False))
else:
    print("Keine strengen OOS-Daten berechnet.")

# =====================================================
# 6k. Analytische Approximation und Kompressionsanalyse
# =====================================================

print("\n=== Analytische Approximation von S(M) ===")

# Hauptansatz: lineares multiplikatives Modell
beta_mult_raw = train_linear_beta(df_sorted, mult_cols)
coef_mult = pd.DataFrame(
    {
        "term": ["intercept"] + mult_cols,
        "coef": beta_mult_raw,
    }
)
print("Lineare Approximation S ≈ b0 + b1*omega + b2*Omega + b3*tau + b4*logM")
print(coef_mult.to_string(index=False))

# Erweiterter Ansatz mit Siebfeatures
beta_both_raw = train_linear_beta(df_sorted, mult_cols + sieve_cols)
coef_both = pd.DataFrame(
    {
        "term": ["intercept"] + (mult_cols + sieve_cols),
        "coef": beta_both_raw,
    }
)
print("\nErweiterte Approximation S ≈ E_mult + E_sieve")
print(coef_both.to_string(index=False))

# Kompressionsmetriken (In-Sample als Referenz)
S_true = df_sorted["S"].to_numpy(dtype=float)
pred_mult_ref = predict_linear_beta(df_sorted, mult_cols, beta_mult_raw)
pred_both_ref = predict_linear_beta(df_sorted, mult_cols + sieve_cols, beta_both_raw)

rmse_base = float(np.sqrt(np.mean((S_true - S_true.mean()) ** 2)))
rmse_mult = float(np.sqrt(np.mean((S_true - pred_mult_ref) ** 2)))
rmse_both = float(np.sqrt(np.mean((S_true - pred_both_ref) ** 2)))

comp_df = pd.DataFrame(
    [
        {
            "model": "mean_baseline",
            "rmse": rmse_base,
            "compression_vs_base": 0.0,
        },
        {
            "model": "E_mult",
            "rmse": rmse_mult,
            "compression_vs_base": 1.0 - (rmse_mult / rmse_base if rmse_base > 0 else np.nan),
        },
        {
            "model": "E_mult_plus_E_sieve",
            "rmse": rmse_both,
            "compression_vs_base": 1.0 - (rmse_both / rmse_base if rmse_base > 0 else np.nan),
        },
    ]
)
print("\nKompressionsanalyse (In-Sample RMSE):")
print(comp_df.to_string(index=False))

# =====================================================
# 6l. Minimalmodell-Suche (strenges OOS)
# =====================================================

print("\n=== Minimalmodell-Suche (strenges OOS) ===")

candidate_cols = ["omega", "Omega", "tau", "logM", "npd_scaled", "lgp_scaled", "sin60", "cos60", "coprime60"]
search_widths = [3000, 10000]
test_bases_mm = [10**7, 10**8]


def prepare_window_df(base, width):
    df_w = build_dataframe(base, base + width)
    df_w = add_sieve_features(df_w)
    eps = 1e-9
    df_w["npd_scaled"] = df_w["nearest_prime_dist"] / (df_w["logM"] + eps)
    df_w["lgp_scaled"] = df_w["local_gap_proxy"] / (df_w["logM"] + eps)
    return df_w


minimal_rows = []

# Fenster einmal vorberechnen (wichtig fuer Laufzeit)
window_cache = {}
for width in search_widths:
    window_cache[("train", width)] = prepare_window_df(train_base, width)
    for bt in test_bases_mm:
        window_cache[("test", bt, width)] = prepare_window_df(bt, width)

# Suche über kleine Modelle bis Größe 4
for k_size in [1, 2, 3, 4]:
    for cols_tuple in it.combinations(candidate_cols, k_size):
        cols = list(cols_tuple)

        auc_scores = []
        rmse_scores = []

        for width in search_widths:
            df_tr = window_cache[("train", width)]
            beta_c, mu_feat, sd_feat = fit_standardized_beta(df_tr, cols)

            ytr = df_tr["S"].to_numpy(dtype=float)
            ltr = df_tr["prime"].to_numpy(dtype=int)
            pred_tr = predict_standardized_beta(df_tr, cols, beta_c, mu_feat, sd_feat)
            resid_tr = ytr - pred_tr
            direction = 1.0 if resid_tr[ltr == 1].mean() >= resid_tr[ltr == 0].mean() else -1.0

            for bt in test_bases_mm:
                df_te = window_cache[("test", bt, width)]
                yte = df_te["S"].to_numpy(dtype=float)
                lte = df_te["prime"].to_numpy(dtype=int)
                pred_te = predict_standardized_beta(df_te, cols, beta_c, mu_feat, sd_feat)
                resid_te = yte - pred_te

                auc_scores.append(signed_auc(lte, resid_te, direction))
                rmse_scores.append(float(np.sqrt(np.mean((yte - pred_te) ** 2))))

        if len(auc_scores) == 0:
            continue

        minimal_rows.append(
            {
                "features": ",".join(cols),
                "k": k_size,
                "auc_oos_mean": float(np.mean(auc_scores)),
                "auc_oos_min": float(np.min(auc_scores)),
                "rmse_oos_mean": float(np.mean(rmse_scores)),
            }
        )

minimal_df = pd.DataFrame(minimal_rows)
if not minimal_df.empty:
    # hohe AUC, niedrige Komplexität, niedrige RMSE
    minimal_df = minimal_df.sort_values(
        by=["auc_oos_mean", "auc_oos_min", "k", "rmse_oos_mean"],
        ascending=[False, False, True, True],
    )
    print(minimal_df.head(12).to_string(index=False))

    print("\n=== Sparsity-Tradeoff (Score = AUC_mean - λ*k) ===")
    lambda_grid = [0.0, 0.001, 0.0025, 0.005, 0.01, 0.02]
    sparse_rows = []

    for lam in lambda_grid:
        tmp = minimal_df.copy()
        tmp["score_sparse"] = tmp["auc_oos_mean"] - lam * tmp["k"]
        tmp = tmp.sort_values(
            by=["score_sparse", "auc_oos_min", "rmse_oos_mean"],
            ascending=[False, False, True],
        )
        best = tmp.iloc[0]
        sparse_rows.append(
            {
                "lambda": lam,
                "features": best["features"],
                "k": int(best["k"]),
                "auc_oos_mean": float(best["auc_oos_mean"]),
                "auc_oos_min": float(best["auc_oos_min"]),
                "rmse_oos_mean": float(best["rmse_oos_mean"]),
                "score_sparse": float(best["score_sparse"]),
            }
        )

    sparse_df = pd.DataFrame(sparse_rows)
    print(sparse_df.to_string(index=False))
else:
    print("Keine Minimalmodelle berechnet.")
    sparse_df = pd.DataFrame()

# =====================================================
# 7. Plot
# =====================================================

plt.figure(figsize=(14, 7))
plt.hist(comps, bins=40, alpha=0.75, label="Zusammengesetzte Zahlen")
plt.hist(primes, bins=20, alpha=0.75, label="Primzahlen")
plt.axvline(mu_c, linestyle="--", label=f"Mittelwert Komposit = {mu_c:.2f}")
plt.axvline(mu_p, linestyle="--", label=f"Mittelwert Prim = {mu_p:.2f}")
plt.xlabel("S(M)")
plt.ylabel("Häufigkeit")
plt.title(f"Verteilung von S(M) im Intervall [{A}, {B}]")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("S_histogramm.png", dpi=200, bbox_inches="tight")
plt.close()

# Residual-Plot 1: Verlauf nach M
plt.figure(figsize=(14, 5))
plt.plot(df.sort_values("M")["M"], R_series, linewidth=1.2)
plt.title("Residuum R(M) über M")
plt.xlabel("M")
plt.ylabel("R(M)")
plt.grid(True, alpha=0.3)
plt.savefig("S_residuum_verlauf.png", dpi=200, bbox_inches="tight")
plt.close()

# Residual-Plot 2: FFT-Amplituden
plt.figure(figsize=(12, 5))
plt.plot(fft_freq, fft_amp, linewidth=1.2)
plt.title("FFT-Amplitudenspektrum des Residuums")
plt.xlabel("Frequenz")
plt.ylabel("|FFT|")
plt.grid(True, alpha=0.3)
plt.savefig("S_residuum_fft.png", dpi=200, bbox_inches="tight")
plt.close()

# Residual-Plot 3: Autokorrelation
plt.figure(figsize=(12, 5))
plt.plot(lags, auto, marker="o", markersize=3, linewidth=1.0)
plt.title("Autokorrelation des Residuums")
plt.xlabel("Lag")
plt.ylabel("ACF")
plt.grid(True, alpha=0.3)
plt.savefig("S_residuum_autokorr.png", dpi=200, bbox_inches="tight")
plt.close()

# =====================================================
# 8. Ausgabe speichern
# =====================================================

df.to_csv("S_test_intervall.csv", index=False)
print("\nCSV gespeichert als S_test_intervall.csv")

if not scale_df.empty:
    scale_df.to_csv("S_test_skalen.csv", index=False)
    print("CSV gespeichert als S_test_skalen.csv")

# Zusatz-CSV fuer Signifikanz/Zerlegung
summary_df = pd.DataFrame(
    [
        {
            "obs_gap_abs_corr": obs_gap_corr,
            "p_gap_perm": p_gap,
            "p_gap_block_perm": p_gap_block,
            "obs_acf_peak": obs_acf_peak,
            "p_acf_perm": p_acf,
            "p_acf_block_perm": p_acf_block,
            "gap_ci_low": ci_gap[0],
            "gap_ci_high": ci_gap[1],
            "gap_block_ci_low": ci_gap_block[0],
            "gap_block_ci_high": ci_gap_block[1],
            "acf_ci_low": ci_acf[0],
            "acf_ci_high": ci_acf[1],
            "acf_block_ci_low": ci_acf_block[0],
            "acf_block_ci_high": ci_acf_block[1],
            "block_len_gap": block_len_gap,
            "block_len_resid": block_len_r,
            "r2_mult": r2_mult,
            "r2_mult_sieve": r2_both,
            "auc_resid_mult": auc_resid_mult,
            "auc_resid_mult_sieve": auc_resid_both,
        }
    ]
)
summary_df.to_csv("S_residuum_signifikanz.csv", index=False)
print("CSV gespeichert als S_residuum_signifikanz.csv")

if not general_df.empty:
    general_df.to_csv("S_generalisierung.csv", index=False)
    print("CSV gespeichert als S_generalisierung.csv")

if not strict_df.empty:
    strict_df.to_csv("S_generalisierung_strikt.csv", index=False)
    print("CSV gespeichert als S_generalisierung_strikt.csv")

coef_mult.to_csv("S_approx_koeff_mult.csv", index=False)
coef_both.to_csv("S_approx_koeff_mult_sieve.csv", index=False)
comp_df.to_csv("S_kompression.csv", index=False)
print("CSV gespeichert als S_approx_koeff_mult.csv")
print("CSV gespeichert als S_approx_koeff_mult_sieve.csv")
print("CSV gespeichert als S_kompression.csv")

if not minimal_df.empty:
    minimal_df.to_csv("S_minimalmodell_oos.csv", index=False)
    print("CSV gespeichert als S_minimalmodell_oos.csv")

if not sparse_df.empty:
    sparse_df.to_csv("S_minimalmodell_sparse.csv", index=False)
    print("CSV gespeichert als S_minimalmodell_sparse.csv")
