import numpy as np
from time import perf_counter

# ---- Load (fast) ----
zeros = np.load("zeros6.npy")  # path: /mnt/data/zeros6.npy if you run in this sandbox


# ---- Range queries (O(log N) per query) ----
def zeros_in_window(t: float, delta: float):
    """
    Return a view (slice) of zeros within [t-delta, t+delta].
    """
    lo = np.searchsorted(zeros, t - delta, side="left")
    hi = np.searchsorted(zeros, t + delta, side="right")
    return zeros[lo:hi]

def window_indices(t: float, delta: float):
    """
    Return (lo, hi) indices in zeros for [t-delta, t+delta].
    """
    lo = np.searchsorted(zeros, t - delta, side="left")
    hi = np.searchsorted(zeros, t + delta, side="right")
    return lo, hi


# ---- Nearest zero (O(log N)) ----
def nearest_zero(t: float):
    i = np.searchsorted(zeros, t, side="left")
    if i <= 0:
        return zeros[0]
    if i >= len(zeros):
        return zeros[-1]
    left, right = zeros[i - 1], zeros[i]
    return left if (t - left) <= (right - t) else right

def nearest_dist(t: float):
    z = nearest_zero(t)
    return abs(t - z), z


# ============================================================
# C1) Peak -> nearest zero distance statistics
# ============================================================
def c1_peak_to_zero(peaks_t):
    """
    peaks_t: iterable of peak locations t*
    Returns arrays: dists, nearest_zeros
    """
    peaks_t = np.asarray(peaks_t, dtype=np.float64)
    dists = np.empty_like(peaks_t)
    nz = np.empty_like(peaks_t)
    for j, t in enumerate(peaks_t):
        d, z = nearest_dist(float(t))
        dists[j] = d
        nz[j] = z
    return dists, nz


# ============================================================
# C2) "Hit mass" inside union of zero windows
#     Hit(delta) = integral_{union_n [t_n-d, t_n+d]} F(t) dt / integral F(t) dt
#
# Works for any discrete grid t_grid + values F(t_grid).
# We approximate integrals by trapezoid rule.
# ============================================================
def c2_hit_mass(t_grid, F, delta):
    """
    t_grid: sorted 1D grid of t values (float)
    F: same length values (nonnegative recommended)
    delta: window half-width around each zero

    Returns: hit_fraction in [0,1]
    """
    t_grid = np.asarray(t_grid, dtype=np.float64)
    F = np.asarray(F, dtype=np.float64)
    assert t_grid.ndim == 1 and F.ndim == 1 and t_grid.size == F.size
    assert np.all(np.diff(t_grid) > 0), "t_grid must be strictly increasing"

    # Total mass (trapezoid)
    total = np.trapezoid(F, t_grid)
    if total <= 0:
        return 0.0

    # Mark grid points that fall inside ANY [zero-delta, zero+delta]
    # Efficient approach:
    # For each grid point t, check if there exists a zero in [t-delta, t+delta]
    # That is equivalent to: count zeros in that interval > 0
    # Use searchsorted twice per t => O(M log N), but M is usually small (few thousand).
    inside = np.zeros(t_grid.size, dtype=bool)
    for i, t in enumerate(t_grid):
        lo = np.searchsorted(zeros, t - delta, side="left")
        hi = np.searchsorted(zeros, t + delta, side="right")
        inside[i] = (hi > lo)

    # Integrate only where inside==True.
    # To do trapezoid properly with a mask, integrate over segments where both endpoints are inside.
    mask = inside.astype(np.float64)
    hit = np.trapezoid(F * mask, t_grid)

    return float(hit / total)


# ============================================================
# Benchmark helpers
# ============================================================
def benchmark():
    # Benchmark: nearest_zero
    rng = np.random.default_rng(0)
    ts = rng.uniform(zeros[0], zeros[-1], size=200_000)

    t0 = perf_counter()
    s = 0.0
    for t in ts:
        s += nearest_zero(float(t))
    t1 = perf_counter()
    print(f"nearest_zero: {len(ts)/(t1-t0):,.0f} queries/sec (dummy sum={s:.3e})")

    # Benchmark: window query
    deltas = rng.uniform(0.01, 0.5, size=50_000)
    t0 = perf_counter()
    cnt = 0
    for t, d in zip(ts[:50_000], deltas):
        lo, hi = window_indices(float(t), float(d))
        cnt += (hi - lo)
    t1 = perf_counter()
    print(f"window_indices: {50_000/(t1-t0):,.0f} queries/sec (total hits={cnt})")


# ============================================================
# Example usage
# ============================================================
if __name__ == "__main__":
    # 1) C1 example: distances for some candidate peaks around 137
    peaks = [136.8, 137.0, 137.2, 138.1]
    peak_dists, nearest_zeros = c1_peak_to_zero(peaks)
    for pt, dist, nearest_z in zip(peaks, peak_dists, nearest_zeros):
        print(f"t*={pt:9.3f}  nearest zero={nearest_z:12.6f}  |Δ|={dist:9.6f}")

    # 2) C2 example: given a grid and a toy signal F(t)
    # (Replace signal with your real resonance curve, e.g. R_BC,k(t) or |ΔG|(t).)
    time_grid = np.linspace(131.0, 143.0, 2401)  # dt=0.005
    # toy signal: bump near 137
    signal = np.exp(-0.5*((time_grid-137.0)/0.12)**2)

    for win_delta in [0.05, 0.1, 0.2, 0.5]:
        hit_frac = c2_hit_mass(time_grid, signal, win_delta)
        print(f"C2 Hit(delta={win_delta}): {hit_frac:.4f}")

    # 3) Benchmarks
    benchmark()