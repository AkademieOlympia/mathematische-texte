import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp, mannwhitneyu

DATA_FILE = "vierlings_dataset.csv"
P_START_COLUMN = "q1"
PROJECTION_COLUMN = "radius_error"

df = pd.read_csv(DATA_FILE)

df["p_start"] = df[P_START_COLUMN]
df["projection_distance"] = df[PROJECTION_COLUMN]
df["p_mod60"] = df["p_start"] % 60

dist_11 = df[df["p_mod60"] == 11]["projection_distance"].dropna().to_numpy()
dist_41 = df[df["p_mod60"] == 41]["projection_distance"].dropna().to_numpy()

def ecdf(x):
    x = np.sort(x)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y

x11, y11 = ecdf(dist_11)
x41, y41 = ecdf(dist_41)

plt.figure(figsize=(10, 6))
plt.step(x11, y11, where="post", label="Kanal 11 mod 60")
plt.step(x41, y41, where="post", label="Kanal 41 mod 60")

plt.title("CDF der Projektdistanz $d_{proj}$ nach 60-smooth-Kanälen")
plt.xlabel("Projektdistanz $d_{proj}$")
plt.ylabel("Kumulative Wahrscheinlichkeit")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("projection_distance_cdf_11_vs_41.png", dpi=200)
plt.show()

ks = ks_2samp(dist_11, dist_41)
mw = mannwhitneyu(dist_11, dist_41, alternative="two-sided")

median_11 = np.median(dist_11)
median_41 = np.median(dist_41)

mean_11 = np.mean(dist_11)
mean_41 = np.mean(dist_41)

std_pooled = np.sqrt((np.var(dist_11, ddof=1) + np.var(dist_41, ddof=1)) / 2)
cohens_d = (mean_41 - mean_11) / std_pooled

print("Datenquelle:", DATA_FILE)
print("p_start-Spalte:", P_START_COLUMN)
print("Projektdistanz-Spalte:", PROJECTION_COLUMN)

print("\nKanal 11:")
print("N =", len(dist_11), "Median =", median_11, "Mean =", mean_11)

print("\nKanal 41:")
print("N =", len(dist_41), "Median =", median_41, "Mean =", mean_41)

print("\nDifferenzen:")
print("Median 41 - Median 11 =", median_41 - median_11)
print("Mean 41 - Mean 11 =", mean_41 - mean_11)
print("Cohen d =", cohens_d)

print("\nKolmogorov-Smirnov-Test:")
print("KS statistic =", ks.statistic)
print("p-value =", ks.pvalue)

print("\nMann-Whitney-U-Test:")
print("U statistic =", mw.statistic)
print("p-value =", mw.pvalue)

print("\n--- Kompakte Testausgabe ---")
print("KS,pKS,MW,pMW,Cohen-d")
print(f"{ks.statistic:.6g},{ks.pvalue:.6g},{mw.statistic:.6g},{mw.pvalue:.6g},{cohens_d:.6g}")

print("\nLaTeX:")
print(
    rf"KS={ks.statistic:.6g},\quad p_{{\mathrm{{KS}}}}={ks.pvalue:.6g},"
    rf"\quad MW={mw.statistic:.6g},\quad p_{{\mathrm{{MW}}}}={mw.pvalue:.6g},"
    rf"\quad \mathrm{{Cohen\text-}}d={cohens_d:.6g}."
)

if ks.pvalue > 0.05 and mw.pvalue > 0.05 and abs(cohens_d) < 0.2:
    print("\nIsotropie-Befund: statistisch sauber dokumentiert.")
else:
    print("\nIsotropie-Befund: nicht eindeutig; Tests oder Effektgröße sprechen gegen eine einfache Isotropie-Deutung.")