import sympy as sp
import numpy as np
import pandas as pd
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

    # Einfache, stabile Surrogat-Summe als Startpunkt.
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

    # Intercept
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

z_gap = (mu_p - mu_c) / sigma_c
print("Z-Abstand der Primzahlen relativ zur Komposit-Verteilung:", z_gap)

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

# Effektgröße Cohen's d
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
score = df["S"].values  # aktuell: groessere S-Werte sprechen fuer Primzahlen

auc = safe_auc(y_true, score)

print("\n=== ROC-Test ===")
print("AUC:", auc)

fpr, tpr, thresholds = roc_curve(y_true, score)

# bester Youden-Index
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

ablations = [
    ["omega"],
    ["Omega"],
    ["tau"],
    ["omega", "Omega"],
    ["omega", "Omega", "tau"],
]

for cols in ablations:
    beta, y_hat, resid = fit_linear_projection(df, cols)
    r2 = linear_r2(df["S"].to_numpy(dtype=float), y_hat)
    auc_resid = orientation_free_auc(y_true, resid)
    print(f"Features={cols} -> R^2={r2:.6f}, AUC(|Richtung frei|) auf Residuum={auc_resid:.6f}")

# Vollmodell (mit logM) fuer robustere Entkopplung
full_cols = ["omega", "Omega", "tau", "logM"]
beta_full, y_hat_full, resid_full = fit_linear_projection(df, full_cols)
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

# =====================================================
# 8. Ausgabe speichern
# =====================================================

df.to_csv("S_test_intervall.csv", index=False)
print("\nCSV gespeichert als S_test_intervall.csv")

if not scale_df.empty:
    scale_df.to_csv("S_test_skalen.csv", index=False)
    print("CSV gespeichert als S_test_skalen.csv")
import sympy as sp
import numpy as np
import pandas as pd
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

    # Einfache, stabile Surrogat-Summe als Startpunkt.
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

    # Intercept
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

z_gap = (mu_p - mu_c) / sigma_c
print("Z-Abstand der Primzahlen relativ zur Komposit-Verteilung:", z_gap)

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

# Effektgröße Cohen's d
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
score = df["S"].values  # aktuell: groessere S-Werte sprechen fuer Primzahlen

auc = safe_auc(y_true, score)

print("\n=== ROC-Test ===")
print("AUC:", auc)

fpr, tpr, thresholds = roc_curve(y_true, score)

# bester Youden-Index
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

ablations = [
    ["omega"],
    ["Omega"],
    ["tau"],
    ["omega", "Omega"],
    ["omega", "Omega", "tau"],
]

for cols in ablations:
    beta, y_hat, resid = fit_linear_projection(df, cols)
    r2 = linear_r2(df["S"].to_numpy(dtype=float), y_hat)
    auc_resid = orientation_free_auc(y_true, resid)
    print(f"Features={cols} -> R^2={r2:.6f}, AUC(|Richtung frei|) auf Residuum={auc_resid:.6f}")

# Vollmodell (mit logM) fuer robustere Entkopplung
full_cols = ["omega", "Omega", "tau", "logM"]
beta_full, y_hat_full, resid_full = fit_linear_projection(df, full_cols)
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

print("\n=== Skalierungstest (gleich breite Fenster) ===")

scale_rows = []

for base in SCALE_BASES:
    a_s = base
    b_s = base + WINDOW_SIZE

    try:
        df_s = build_dataframe(a_s, b_s)
    except Exception as e:
        print(f"Fenster [{a_s}, {b_s}] fehlgeschlagen: {e}")
        continue

    y_s = df_s["prime"].to_numpy(dtype=int)
    auc_s = safe_auc(y_s, df_s["S"].to_numpy(dtype=float))

    _, y_hat_s, resid_s = fit_linear_projection(df_s, full_cols)
    r2_s = linear_r2(df_s["S"].to_numpy(dtype=float), y_hat_s)
    auc_resid_s = orientation_free_auc(y_s, resid_s)

    scale_rows.append(
        {
            "A": a_s,
            "B": b_s,
            "n": len(df_s),
            "primes": int(df_s["prime"].sum()),
            "auc_S": auc_s,
            "r2_struct": r2_s,
            "auc_resid": auc_resid_s,
        }
    )

scale_df = pd.DataFrame(scale_rows)

if not scale_df.empty:
    print(scale_df.to_string(index=False))
else:
    print("Keine Skalierungsdaten berechnet.")

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

# =====================================================
# 8. Ausgabe speichern
# =====================================================

df.to_csv("S_test_intervall.csv", index=False)
print("\nCSV gespeichert als S_test_intervall.csv")

if not scale_df.empty:
    scale_df.to_csv("S_test_skalen.csv", index=False)
    print("CSV gespeichert als S_test_skalen.csv")
import sympy as sp
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt

# =====================================================
# 1. Intervall
# =====================================================

A = 113000
B = 113300

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

    # Einfache, stabile Surrogat-Summe als Startpunkt.
    return float(np.log1p(M) - 0.75 * omega - 0.35 * Omega + 0.02 * tau)


# =====================================================
# 3. Daten erzeugen
# =====================================================

rows = []

for M in range(A, B + 1):
    try:
        s = S(M)
    except Exception as e:
        print(f"Fehler bei M={M}: {e}")
        continue

    is_prime = sp.isprime(M)

    rows.append({
        "M": M,
        "S": float(s),
        "prime": int(is_prime),
        "composite": int((M > 1) and not is_prime),
        "omega": len(sp.factorint(M)),
        "Omega": sum(sp.factorint(M).values()),
        "factorization": sp.factorint(M)
    })

df = pd.DataFrame(rows)

if df.empty:
    raise RuntimeError("Keine gueltigen Daten erzeugt. Pruefe S(M) und Intervall.")

primes = df[df["prime"] == 1]["S"]
comps = df[df["composite"] == 1]["S"]

if len(primes) == 0 or len(comps) == 0:
    raise RuntimeError("Zu wenige Klassen fuer Statistik (Primzahlen oder Komposita fehlen).")

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

z_gap = (mu_p - mu_c) / sigma_c
print("Z-Abstand der Primzahlen relativ zur Komposit-Verteilung:", z_gap)

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

# Effektgröße Cohen's d
pooled = np.sqrt(((len(comps)-1)*sigma_c**2 + (len(primes)-1)*sigma_p**2) / (len(comps)+len(primes)-2))
cohen_d = (mu_c - mu_p) / pooled

print("Cohen d:", cohen_d)

# =====================================================
# 6. Klassifikationstest: Kann S(M) Primzahlen erkennen?
# =====================================================

# Falls kleinere S-Werte eher Primzahlen bedeuten:
y_true = df["prime"].values
# In den aktuellen Daten gilt S(prim) > S(komposit),
# daher muss der Score direkt S sein.
score = df["S"].values

auc = roc_auc_score(y_true, score)

print("\n=== ROC-Test ===")
print("AUC:", auc)

fpr, tpr, thresholds = roc_curve(y_true, score)

# bester Youden-Index
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

print("\n=== Korrelationsmatrix S/omega/Omega ===")
print(df[["S", "omega", "Omega"]].corr())

# =====================================================
# 6c. Carmichael-Schnelltest
# =====================================================

carmichael = [
    561, 1105, 1729, 2465, 2821,
    6601, 8911, 10585, 15841
]

print("\n=== Carmichael-Test ===")
for n in carmichael:
    try:
        print(n, S(n))
    except Exception as e:
        print(n, f"Fehler: {e}")

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

# =====================================================
# 8. Ausgabe speichern
# =====================================================

df.to_csv("S_test_intervall.csv", index=False)
print("\nCSV gespeichert als S_test_intervall.csv")