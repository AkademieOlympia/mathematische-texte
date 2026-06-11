import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sympy import sieve

# --- KONFIGURATION ---
LIMIT = 500000
STEP = 1000
LABELS = ["E", "A", "DA", "C", "DB", "B"] # Hexagonale Anordnung
FIGSIZE = (9, 6.5)
EXPORT_BASENAME = "hexgon_ptolemaeus"
EXPORT_DPI = 300
SHOW_PLOT = False
TEST_TRIALS = 5000
TEST_RANDOM_SEED = 42
EXPORT_TEST_RESULTS = True
TEST_SUMMARY_NAME = "hexgon_quintuplet_test_summary.csv"
TEST_MINIMA_NAME = "hexgon_quintuplet_minima.csv"
MIN_PROMINENCE_QUANTILE = 0.75

# Geometrie-Setup: Hexagon-Vektoren
ANGLES = np.linspace(0, 2*np.pi, 6, endpoint=False)
VECTORS = {LABELS[i]: np.array([np.cos(ANGLES[i]), np.sin(ANGLES[i])]) for i in range(6)}

def is_sum_of_two_squares(n, p34_list):
    """Prüft die DA/DB-Zugehörigkeit (Zwei-Quadrate-Satz)."""
    if n == 1: return True
    temp = n
    for p in p34_list:
        if p*p > temp: break
        if temp % p == 0:
            c = 0
            while temp % p == 0:
                c += 1
                temp //= p
            if c % 2 != 0: return False
    return not (temp > 1 and temp % 4 == 3)

def calculate_simulation():
    p34 = [p for p in sieve.primerange(3, LIMIT + 1) if p % 4 == 3]
    prime_set = set(sieve.primerange(1, LIMIT + 1))
    
    stats = {l: {'count': 0, 'log_sum': 0.0} for l in LABELS}
    results = []
    
    for n in range(1, LIMIT + 1):
        ln = np.log1p(n)
        # Klassifizierung nach Ihrer Geometrie
        if n == 1: cat = 'E'
        elif n in prime_set:
            if n == 2: cat = 'C'
            elif n % 4 == 1: cat = 'A'
            else: cat = 'B'
        else:
            cat = 'DA' if is_sum_of_two_squares(n, p34) else 'DB'
        
        stats[cat]['count'] += 1
        stats[cat]['log_sum'] += ln
        
        if n % STEP == 0:
            # 1. Ptolemäische Spannung berechnen
            centers = {}
            for l in ["A", "DA", "B", "DB"]:
                r = stats[l]['log_sum'] / max(1, stats[l]['count'])
                centers[l] = r * VECTORS[l]
            
            dist = lambda p1, p2: np.linalg.norm(p1 - p2)
            a, b = dist(centers['A'], centers['DA']), dist(centers['DA'], centers['B'])
            c, d = dist(centers['B'], centers['DB']), dist(centers['DB'], centers['A'])
            e, f = dist(centers['A'], centers['B']), dist(centers['DA'], centers['DB'])
            
            # Spannung Delta_P = |ef - (ac + bd)|
            tension = np.abs(e * f - (a * c + b * d)) / (e * f + 1e-9)
            
            # 2. Entropie berechnen
            cnts = [stats[l]['count'] for l in LABELS]
            total = sum(cnts)
            probs = [c/total for c in cnts if c > 0]
            entropy = -sum(p * np.log(p) for p in probs)
            
            results.append({
                'N': n,
                'Ptolemy_Tension': tension,
                'Entropy': entropy,
                'Bary_X': np.average([VECTORS[l][0] for l in LABELS], weights=np.log1p(cnts)),
                'Bary_Y': np.average([VECTORS[l][1] for l in LABELS], weights=np.log1p(cnts))
            })
            
    return pd.DataFrame(results)


def find_local_minima(n_values, values, prominence_quantile=MIN_PROMINENCE_QUANTILE):
    n_values = np.asarray(n_values)
    values = np.asarray(values, dtype=float)
    if values.size < 3:
        return pd.DataFrame(columns=["N", "Ptolemy_Tension", "Prominence"])

    minima_mask = (values[1:-1] < values[:-2]) & (values[1:-1] <= values[2:])
    minima_positions = np.flatnonzero(minima_mask) + 1
    if minima_positions.size == 0:
        return pd.DataFrame(columns=["N", "Ptolemy_Tension", "Prominence"])

    left_drop = values[minima_positions - 1] - values[minima_positions]
    right_drop = values[minima_positions + 1] - values[minima_positions]
    prominence = np.minimum(left_drop, right_drop)
    threshold = np.quantile(prominence, prominence_quantile)
    keep_mask = prominence >= threshold
    return pd.DataFrame(
        {
            "N": n_values[minima_positions][keep_mask],
            "Ptolemy_Tension": values[minima_positions][keep_mask],
            "Prominence": prominence[keep_mask],
        }
    )


def nearest_reference_distances(target_values, reference_values):
    target_values = np.asarray(target_values, dtype=np.int64)
    reference_values = np.asarray(reference_values, dtype=np.int64)
    if target_values.size == 0 or reference_values.size == 0:
        empty = np.full(target_values.shape, np.nan, dtype=float)
        return empty, empty

    positions = np.searchsorted(reference_values, target_values)
    positions = np.clip(positions, 0, reference_values.size - 1)
    left_positions = np.clip(positions - 1, 0, reference_values.size - 1)

    right_values = reference_values[positions]
    left_values = reference_values[left_positions]
    use_left = np.abs(target_values - left_values) <= np.abs(right_values - target_values)
    nearest = np.where(use_left, left_values, right_values)
    distances = np.abs(target_values - nearest)
    return nearest.astype(float), distances.astype(float)


def shift_reference_grid(reference_values, shift_steps, n_min, n_max, step):
    reference_values = np.asarray(reference_values, dtype=np.int64)
    if reference_values.size == 0:
        return reference_values

    grid_size = int((n_max - n_min) // step) + 1
    normalized = (reference_values - n_min) // step
    shifted = (normalized + shift_steps) % grid_size
    return np.sort(n_min + shifted * step)


def run_quintuplet_proximity_test(n_grid, tension_values, quintuplet_grid, trials=TEST_TRIALS, seed=TEST_RANDOM_SEED):
    n_grid = np.asarray(n_grid, dtype=np.int64)
    quintuplet_grid = np.asarray(quintuplet_grid, dtype=np.int64)
    minima_df = find_local_minima(n_grid, tension_values)
    minima = minima_df["N"].to_numpy(dtype=np.int64)
    nearest_quintuplets, observed_distances = nearest_reference_distances(minima, quintuplet_grid)

    minima_df["Nearest_Quintuplet"] = nearest_quintuplets
    minima_df["Distance_To_Quintuplet"] = observed_distances

    if minima.size == 0 or quintuplet_grid.size == 0:
        summary = pd.DataFrame(
            [
                {
                    "minima_count": int(minima.size),
                    "prominence_quantile": float(MIN_PROMINENCE_QUANTILE),
                    "quintuplet_count": int(quintuplet_grid.size),
                    "observed_mean_distance": np.nan,
                    "random_null_mean_distance": np.nan,
                    "random_null_std_distance": np.nan,
                    "random_z_score": np.nan,
                    "random_p_value_left": np.nan,
                    "shift_null_mean_distance": np.nan,
                    "shift_null_std_distance": np.nan,
                    "shift_z_score": np.nan,
                    "shift_p_value_left": np.nan,
                    "trials": int(trials),
                    "step": int(STEP),
                }
            ]
        )
        return summary, minima_df

    observed_mean = float(np.mean(observed_distances))
    rng = np.random.default_rng(seed)
    random_null_means = np.empty(trials, dtype=float)
    shift_null_means = np.empty(trials, dtype=float)
    n_min = int(n_grid[0])
    n_max = int(n_grid[-1])
    grid_size = n_grid.size
    for idx in range(trials):
        sample = np.sort(rng.choice(n_grid, size=minima.size, replace=False))
        _, null_distances = nearest_reference_distances(sample, quintuplet_grid)
        random_null_means[idx] = np.mean(null_distances)

        shift_steps = int(rng.integers(1, grid_size))
        shifted_quintuplets = shift_reference_grid(quintuplet_grid, shift_steps, n_min, n_max, STEP)
        _, shifted_distances = nearest_reference_distances(minima, shifted_quintuplets)
        shift_null_means[idx] = np.mean(shifted_distances)

    random_null_mean = float(np.mean(random_null_means))
    random_null_std = float(np.std(random_null_means, ddof=1))
    random_z_score = (observed_mean - random_null_mean) / random_null_std if random_null_std > 0 else np.nan
    random_p_value_left = float((1 + np.count_nonzero(random_null_means <= observed_mean)) / (trials + 1))

    shift_null_mean = float(np.mean(shift_null_means))
    shift_null_std = float(np.std(shift_null_means, ddof=1))
    shift_z_score = (observed_mean - shift_null_mean) / shift_null_std if shift_null_std > 0 else np.nan
    shift_p_value_left = float((1 + np.count_nonzero(shift_null_means <= observed_mean)) / (trials + 1))

    summary = pd.DataFrame(
        [
            {
                "minima_count": int(minima.size),
                "prominence_quantile": float(MIN_PROMINENCE_QUANTILE),
                "quintuplet_count": int(quintuplet_grid.size),
                "observed_mean_distance": observed_mean,
                "random_null_mean_distance": random_null_mean,
                "random_null_std_distance": random_null_std,
                "random_z_score": random_z_score,
                "random_p_value_left": random_p_value_left,
                "shift_null_mean_distance": shift_null_mean,
                "shift_null_std_distance": shift_null_std,
                "shift_z_score": shift_z_score,
                "shift_p_value_left": shift_p_value_left,
                "trials": int(trials),
                "step": int(STEP),
            }
        ]
    )
    return summary, minima_df

# --- EXECUTION ---
df = calculate_simulation()

# Fünflinge für Plotting
ps = list(sieve.primerange(1, LIMIT))
pset = set(ps)
qs = [p for p in ps if all((p+d) in pset for d in [2,6,8,12]) or all((p+d) in pset for d in [4,6,10,12])]
sampled_quintuplets = np.unique(((np.array(qs) + STEP // 2) // STEP) * STEP)
sampled_quintuplets = sampled_quintuplets[(sampled_quintuplets >= STEP) & (sampled_quintuplets <= LIMIT)]
quintuplet_rows = df[df["N"].isin(sampled_quintuplets)]
test_summary_df, test_minima_df = run_quintuplet_proximity_test(
    df["N"].to_numpy(),
    df["Ptolemy_Tension"].to_numpy(),
    sampled_quintuplets,
)
if EXPORT_TEST_RESULTS:
    test_summary_df.to_csv(Path(__file__).with_name(TEST_SUMMARY_NAME), index=False)
    test_minima_df.to_csv(Path(__file__).with_name(TEST_MINIMA_NAME), index=False)

# Visualisierung
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=FIGSIZE, constrained_layout=True)

# Plot 1: Ptolemäische Spannung & Fünflinge
ax1.plot(df['N'], df['Ptolemy_Tension'], color='darkblue', label=r'Ptolemäische Spannung $\Delta_P$')
ax1.scatter(
    quintuplet_rows["N"],
    quintuplet_rows["Ptolemy_Tension"],
    color="orange",
    s=16,
    alpha=0.65,
    label="Fünflingsmarken",
    zorder=3,
)
ax1.set_title("Harmonische Stabilität: Ptolemäische Spannung vs. Fünfling-Singularitäten")
ax1.set_ylabel("Relative Spannung")
ax1.legend()

# Plot 2: Baryzentrischer Jitter (Spektralanalyse)
jitter = np.diff(df['Ptolemy_Tension'])
fft = np.abs(np.fft.rfft(jitter - np.mean(jitter)))
freqs = np.fft.rfftfreq(len(jitter), d=1) * (2*np.pi * (LIMIT/STEP))

ax2.plot(freqs, fft, color='purple')
ax2.set_xlim(0, 60)
ax2.set_title("Spektrale Signatur (Riemann-Resonanz)")
ax2.set_xlabel("Frequenz (Gamma)")

output_base = Path(__file__).with_name(EXPORT_BASENAME)
fig.savefig(output_base.with_suffix(".png"), dpi=EXPORT_DPI, bbox_inches="tight", pad_inches=0.05)
fig.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.05)
if SHOW_PLOT:
    plt.show()
else:
    plt.close(fig)