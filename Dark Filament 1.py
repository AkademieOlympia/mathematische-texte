import os, math
import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# Config
# -------------------------
N_PART = 120_000          # DM "particles" on the sphere
N_MODES = 80              # random potential modes (controls coherence scale)
A_STEPS = 220             # number of frames
A_MAX = 1.25              # end "time" parameter (bigger => stronger nonlinearity)
ADHESION = True           # turn on simple stickiness
NU = 0.60                 # stickiness strength (only if ADHESION=True)

NBIN_LON = 720            # density map resolution (lambda bins)
NBIN_LAT = 360            # density map resolution (phi bins)

OUTDIR = "frames_dm"
DPI = 140

# -------------------------
# Utilities
# -------------------------
def fibonacci_sphere(n):
    """Approximately uniform points on S^2."""
    i = np.arange(n)
    phi = (1 + 5**0.5) / 2
    theta = 2*np.pi*i/phi
    z = 1 - 2*(i + 0.5)/n
    r = np.sqrt(np.maximum(0.0, 1 - z*z))
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    return np.stack([x,y,z], axis=1)

def normalize(v):
    n = np.linalg.norm(v, axis=1, keepdims=True)
    n = np.where(n < 1e-12, 1.0, n)
    return v / n

def random_potential_setup(n_modes, seed=7):
    rng = np.random.default_rng(seed)
    # random directions k on S^2 and random phases
    k = normalize(rng.normal(size=(n_modes, 3)))
    # wave numbers: small integers => large-scale structure on sky
    # tweak range for coarser/finer filament network
    kscale = rng.integers(low=2, high=16, size=n_modes).astype(np.float64)
    phase = rng.uniform(0, 2*np.pi, size=n_modes)
    amp = rng.normal(0, 1, size=n_modes) / (kscale**1.2)  # red-ish spectrum
    return k, kscale, phase, amp

def grad_sphere_potential(nvec, k, kscale, phase, amp):
    """
    Potential Phi(n) = sum amp_i cos( kscale_i * (k_i · n) + phase_i )
    Gradient on the sphere: grad_s = (I - n n^T) grad_3D Phi
    where grad_3D cos(u) = -sin(u) * d(u)/dn = -sin(u) * kscale * k_i
    """
    # u shape: (N, M)
    dot = nvec @ k.T                       # (N, M)
    u = dot * kscale[None, :] + phase[None, :]
    su = np.sin(u)                         # (N, M)

    # grad_3D Phi = sum_i [ -amp_i * sin(u_i) * kscale_i * k_i ]
    # compute as matrix product: su * (amp*kscale) then times k
    w = (-su) * (amp * kscale)[None, :]    # (N, M)
    g3 = w @ k                             # (N, 3)

    # project to tangent plane: g_s = g3 - (g3·n) n
    proj = np.sum(g3 * nvec, axis=1, keepdims=True)
    gs = g3 - proj * nvec
    return gs  # (N,3)

def lonlat_from_vec(v):
    x,y,z = v[:,0], v[:,1], v[:,2]
    lon = np.arctan2(y, x)                 # [-pi, pi]
    lat = np.arcsin(np.clip(z, -1, 1))     # [-pi/2, pi/2]
    return lon, lat

def density_map(lon, lat, nlon, nlat):
    # shift lon to [0, 2pi) for binning
    lon2 = (lon + np.pi) % (2*np.pi)
    # lat already [-pi/2, pi/2]
    H, _, _ = np.histogram2d(
        lat, lon2,
        bins=[nlat, nlon],
        range=[[-np.pi/2, np.pi/2], [0, 2*np.pi]]
    )
    # mild log scaling for visibility
    return np.log1p(H)

# -------------------------
# Main
# -------------------------
def main():
    os.makedirs(OUTDIR, exist_ok=True)

    # Initial particle directions on S^2
    q = fibonacci_sphere(N_PART)

    # Random potential modes
    k, kscale, phase, amp = random_potential_setup(N_MODES)

    # precompute a baseline displacement direction at initial positions
    s0 = -grad_sphere_potential(q, k, kscale, phase, amp)

    # normalize displacement magnitude somewhat
    s0_norm = np.linalg.norm(s0, axis=1, keepdims=True)
    s0 = s0 / np.maximum(s0_norm, 1e-9)

    for t in range(A_STEPS):
        # "time" parameter a in [0, A_MAX]
        a = A_MAX * (t / (A_STEPS - 1))

        # growth factor proxy D(a): smooth nonlinear ramp
        D = a

        # Zel'dovich step on sphere: n = normalize(q + D * s0)
        n = normalize(q + D * s0)

        if ADHESION:
            # Simple "stickiness": compress motion where displacement is strong,
            # mimicking caustic adhesion qualitatively (proxy for Burgers/adhesion).
            # We use an additional displacement computed at current n:
            s1 = -grad_sphere_potential(n, k, kscale, phase, amp)
            s1n = np.linalg.norm(s1, axis=1, keepdims=True)
            s1u = s1 / np.maximum(s1n, 1e-9)

            # compression factor reduces step in high-gradient regions
            comp = 1.0 / (1.0 + NU * s1n)
            n = normalize(n + 0.35 * D * comp * s1u)

        # build density map on lon/lat grid
        lon, lat = lonlat_from_vec(n)
        Z = density_map(lon, lat, NBIN_LON, NBIN_LAT)

        # render Mollweide using matplotlib's built-in projection
        fig = plt.figure(figsize=(11.5, 6.5), dpi=DPI)
        ax = fig.add_subplot(111, projection='mollweide')
        ax.grid(True, linewidth=0.6, alpha=0.5)

        # For pcolormesh: lon in [-pi,pi], lat in [-pi/2, pi/2]
        lon_edges = np.linspace(-np.pi, np.pi, NBIN_LON + 1)
        lat_edges = np.linspace(-np.pi/2, np.pi/2, NBIN_LAT + 1)

        # Z is indexed [lat, lon] but our histogram used lon in [0,2pi),
        # so we roll half-way to align with [-pi,pi]
        Z2 = np.roll(Z, shift=NBIN_LON//2, axis=1)

        ax.pcolormesh(lon_edges, lat_edges, Z2, shading='auto')

        ax.set_title(f"DM-Filamente (Zel’dovich{' + Adhesion' if ADHESION else ''}) | a={a:.3f}", pad=18)

        outpath = os.path.join(OUTDIR, f"frame_{t:04d}.png")
        fig.savefig(outpath, dpi=DPI)
        plt.close(fig)

    print(f"Frames written to {OUTDIR}/")
    print("Encode example:")
    print("  ffmpeg -framerate 30 -i frames_dm/frame_%04d.png -c:v libx264 -pix_fmt yuv420p dm_filaments.mp4")

if __name__ == "__main__":
    main()