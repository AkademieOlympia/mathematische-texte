import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


# CSV laden
df = pd.read_csv("plot")   # <- anpassen!

# Zeiten in µs
aqg = df["aqg_us"].values
classic = df["classic_us"].values
gnfs = df["gnfs_us"].values

# =========================
# 1. Ratio Histogram
# =========================
ratio = aqg / classic
log_ratio = np.log10(ratio)

plt.figure()
plt.hist(log_ratio, bins=50)
plt.title("Log10 Ratio AQG / Classical")
plt.xlabel("log10(AQG / Classical)")
plt.ylabel("Frequency")
plt.axvline(0)  # Gleich schnell
plt.show()

# =========================
# 2. Scatter Plot
# =========================
plt.figure()
plt.scatter(classic, aqg, alpha=0.5)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Classical (µs)")
plt.ylabel("AQG (µs)")
plt.title("AQG vs Classical Runtime")

# Diagonale (gleich schnell)
x = np.linspace(min(classic), max(classic), 100)
plt.plot(x, x)

plt.show()

# =========================
# 3. Hybrid Simulation
# =========================
hybrid = np.minimum(aqg, classic)

print("Median AQG:", np.median(aqg))
print("Median Classical:", np.median(classic))
print("Median Hybrid:", np.median(hybrid))

print("Mean AQG:", np.mean(aqg))
print("Mean Classical:", np.mean(classic))
print("Mean Hybrid:", np.mean(hybrid))

# Speedup
speedup = classic / hybrid

plt.figure()
plt.hist(speedup, bins=50)
plt.title("Hybrid Speedup (Classical / Hybrid)")
plt.xlabel("Speedup factor")
plt.ylabel("Frequency")
plt.show()

# =========================
# 4. KMeans-Clustering
# =========================
# Nur gueltige positive Werte verwenden (log10 braucht > 0).
mask = (classic > 0) & (aqg > 0)
X = np.column_stack((np.log10(classic[mask]), np.log10(aqg[mask])))

kmeans = KMeans(n_clusters=2, random_state=42, n_init=10).fit(X)
labels = kmeans.labels_

plt.figure()
plt.scatter(X[:, 0], X[:, 1], c=labels, s=14, alpha=0.6)
plt.xlabel("log10(Classical)")
plt.ylabel("log10(AQG)")
plt.title("KMeans Cluster (n=2)")
plt.show()

print("Cluster-Zentren (log10-space):")
print(kmeans.cluster_centers_)

# Cluster-Interpretation: Welcher Cluster ist AQG-vorteilhaft?
ratio_masked = aqg[mask] / classic[mask]
for c in range(2):
    c_ratio = ratio_masked[labels == c]
    if c_ratio.size == 0:
        print(f"Cluster {c}: leer")
        continue
    aqg_better_share = np.mean(c_ratio < 1.0) * 100.0
    print(
        f"Cluster {c}: n={c_ratio.size}, "
        f"mean(aqg/classic)={np.mean(c_ratio):.4f}, "
        f"median={np.median(c_ratio):.4f}, "
        f"AQG_besser={aqg_better_share:.1f}%"
    )

# =========================
# 5. Regime-Bucket ueber log_ratio
# =========================
valid = (df["aqg_us"] > 0) & (df["classic_us"] > 0)
df_valid = df.loc[valid].copy()
df_valid["log_ratio"] = np.log10(df_valid["aqg_us"] / df_valid["classic_us"])

df_valid["regime"] = pd.cut(
    df_valid["log_ratio"],
    bins=[-np.inf, 0, 1, np.inf],
    labels=["AQG_fast", "neutral", "AQG_slow"],
)

print("\nMittelwert aqg_factor_count pro Regime:")
print(df_valid.groupby("regime")["aqg_factor_count"].mean())