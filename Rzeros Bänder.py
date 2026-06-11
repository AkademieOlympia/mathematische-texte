import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Fenster-/Blockstabilität der 3-Zonen-Struktur
zeros = np.load("zeros6.npy").astype(float).ravel()
zeros = np.sort(zeros[np.isfinite(zeros)])

BLOCK_SIZE = 10000
N_BLOCKS = 4
WINDOW_LEN = 80
STEP = 40
N_MC = 20
NBINS = 80

def build_cumgap(window):
    return np.cumsum(np.diff(np.asarray(window, dtype=float)))

def transform_log(x):
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0):
        return None
    return np.log(x)

def density_from_windows(zeros_1d, bins):
    hists = []
    for start in range(0, len(zeros_1d) - WINDOW_LEN + 1, STEP):
        y = transform_log(build_cumgap(zeros_1d[start:start + WINDOW_LEN]))
        if y is None or len(y) < 10:
            continue
        hist, _ = np.histogram(y, bins=bins, density=True)
        hists.append(hist)
    return np.vstack(hists).mean(axis=0)

def gap_shuffle_control(zeros_1d, rng):
    gaps = np.diff(zeros_1d)
    gaps_perm = rng.permutation(gaps)
    out = np.empty_like(zeros_1d)
    out[0] = zeros_1d[0]
    out[1:] = out[0] + np.cumsum(gaps_perm)
    return np.sort(out)

def gauss(x, A, mu, sigma):
    sigma = np.maximum(np.abs(sigma), 1e-6)
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def model(x, A1, mu1, s1, A2, mu2, s2, A3, mu3, s3, c):
    return gauss(x, A1, mu1, s1) - gauss(x, A2, mu2, s2) + gauss(x, A3, mu3, s3) + c

def fit_three_zone(block, rng):
    # fixed bins from block
    all_y = []
    for start in range(0, len(block) - WINDOW_LEN + 1, STEP):
        y = transform_log(build_cumgap(block[start:start + WINDOW_LEN]))
        if y is not None and len(y) >= 10:
            all_y.append(y)
    all_y = np.concatenate(all_y)
    bins = np.linspace(float(np.min(all_y)), float(np.max(all_y)), NBINS + 1)
    x = 0.5 * (bins[:-1] + bins[1:])

    obs_density = density_from_windows(block, bins)
    ctrls = []
    for _ in range(N_MC):
        ctrls.append(density_from_windows(gap_shuffle_control(block, rng), bins))
    ctrls = np.vstack(ctrls)
    ctrl_mean = ctrls.mean(axis=0)
    obs_resid = obs_density - ctrl_mean

    kernel = np.array([1, 2, 3, 2, 1], dtype=float)
    kernel /= kernel.sum()
    obs_resid_smooth = np.convolve(obs_resid, kernel, mode='same')

    mask = (x >= 1.0) & (x <= 5.3)
    xf = x[mask]
    yf = obs_resid_smooth[mask]

    p0 = [0.010, 3.17, 0.90, 0.110, 4.20, 0.18, 0.040, 4.62, 0.22, 0.0]
    lower = [0.0, 1.0, 0.05, 0.0, 3.7, 0.05, 0.0, 4.3, 0.05, -0.05]
    upper = [0.2, 4.2, 2.5, 0.4, 4.5, 0.6, 0.2, 5.3, 1.2, 0.05]

    pars, _ = curve_fit(model, xf, yf, p0=p0, bounds=(lower, upper), maxfev=50000)
    fit = model(xf, *pars)
    ss_res = float(np.sum((yf - fit)**2))
    ss_tot = float(np.sum((yf - np.mean(yf))**2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    rmse = float(np.sqrt(np.mean((yf - fit)**2)))
    return x, obs_resid_smooth, xf, yf, fit, pars, rmse, r2

rows = []
profiles = []
rng = np.random.default_rng(12345)

for b in range(N_BLOCKS):
    start = b * BLOCK_SIZE
    end = start + BLOCK_SIZE
    block = zeros[start:end]
    x, obs_resid_smooth, xf, yf, fit, pars, rmse, r2 = fit_three_zone(block, rng)
    rows.append({
        "block_id": b + 1,
        "zero_start_index": start,
        "zero_end_index": end - 1,
        "rmse": rmse,
        "r2": r2,
        "A_left": pars[0],
        "mu_left": pars[1],
        "sigma_left": pars[2],
        "A_mid": pars[3],
        "mu_mid": pars[4],
        "sigma_mid": pars[5],
        "A_right": pars[6],
        "mu_right": pars[7],
        "sigma_right": pars[8],
        "offset": pars[9],
        "left_area": float(pars[0] * abs(pars[2]) * np.sqrt(2*np.pi)),
        "mid_area": float(pars[3] * abs(pars[5]) * np.sqrt(2*np.pi)),
        "right_area": float(pars[6] * abs(pars[8]) * np.sqrt(2*np.pi)),
    })
    for xi, yi in zip(xf, yf):
        profiles.append({"block_id": b + 1, "x": xi, "residual_smooth": yi})

summary_df = pd.DataFrame(rows)
profiles_df = pd.DataFrame(profiles)

summary_path = "three_zone_block_stability_summary.csv"
profiles_path = "three_zone_block_stability_profiles.csv"
summary_df.to_csv(summary_path, index=False)
profiles_df.to_csv(profiles_path, index=False)

# plot key parameters across blocks
plot1 = "three_zone_block_stability_mus.png"
plt.figure(figsize=(8,5))
plt.plot(summary_df["block_id"], summary_df["mu_left"], marker="o", label="mu_left")
plt.plot(summary_df["block_id"], summary_df["mu_mid"], marker="o", label="mu_mid")
plt.plot(summary_df["block_id"], summary_df["mu_right"], marker="o", label="mu_right")
plt.xlabel("Block")
plt.ylabel("mu")
plt.title("Stabilität der Zonenlagen über Blöcke")
plt.legend()
plt.tight_layout()
plt.savefig(plot1, dpi=180)
plt.close()

plot2 = "three_zone_block_stability_sigmas.png"
plt.figure(figsize=(8,5))
plt.plot(summary_df["block_id"], summary_df["sigma_left"], marker="o", label="sigma_left")
plt.plot(summary_df["block_id"], summary_df["sigma_mid"], marker="o", label="sigma_mid")
plt.plot(summary_df["block_id"], summary_df["sigma_right"], marker="o", label="sigma_right")
plt.xlabel("Block")
plt.ylabel("sigma")
plt.title("Stabilität der Zonenbreiten über Blöcke")
plt.legend()
plt.tight_layout()
plt.savefig(plot2, dpi=180)
plt.close()

plot3 = "three_zone_block_stability_areas.png"
plt.figure(figsize=(8,5))
plt.plot(summary_df["block_id"], summary_df["left_area"], marker="o", label="left_area")
plt.plot(summary_df["block_id"], summary_df["mid_area"], marker="o", label="mid_area")
plt.plot(summary_df["block_id"], summary_df["right_area"], marker="o", label="right_area")
plt.xlabel("Block")
plt.ylabel("approximative Fläche")
plt.title("Stabilität der Zonenflächen über Blöcke")
plt.legend()
plt.tight_layout()
plt.savefig(plot3, dpi=180)
plt.close()

plot4 = "three_zone_block_stability_fitquality.png"
plt.figure(figsize=(8,5))
plt.plot(summary_df["block_id"], summary_df["r2"], marker="o", label="R²")
plt.plot(summary_df["block_id"], summary_df["rmse"], marker="o", label="RMSE")
plt.xlabel("Block")
plt.ylabel("Fitgüte")
plt.title("Fitgüte des 3-Zonen-Modells über Blöcke")
plt.legend()
plt.tight_layout()
plt.savefig(plot4, dpi=180)
plt.close()

print("Zusammenfassung:")
print(summary_df.to_string(index=False))
print("\nDateien:")
for p in [summary_path, profiles_path, plot1, plot2, plot3, plot4]:
    print(p)