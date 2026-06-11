"""
Mollweide-Visualisierung: Primzahlen als Punkte auf der Sphäre (Hopf-Abbildung)
nach Klassen mod 12 (E/A/B/C/P), projiziert in die Mollweide-Ebene.

Benötigt optional: zeros_gamma.txt (eine Zeile pro Gamma = Im(ρ) der Riemann-
Zeta-Nullstellen). Fehlt die Datei, werden synthetische Gammas verwendet.
Zum Erzeugen von zeros_gamma.txt siehe z.B. Zeta_Funktion_33_Nullpunkt.ipynb
(Imaginärteile der Nullstellen zeilenweise in eine Textdatei schreiben).
"""
import os
import math
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------
# Settings (tweakable)
# ---------------------------
T = 1000                 # number of primes / frames base
FRAME_STEP = 1           # set >1 to speed up (e.g. 2,5,10)
TRAIL_W = 150            # last W points as "trail"
ALPHA_DECAY = 250.0      # larger = slower fading of older points
BETA_GAP = 0.25          # weight of gap-signal in 4D vector
ALPHA_GAP = 0.8          # scaling inside tanh()
R = 1.0                  # Mollweide radius scale
OUTDIR = "frames"        # where PNG frames are saved
DPI = 140

# Marker shapes for overlay (no need for strong colors)
MARKERS = {"E": "o", "A": "s", "B": "^", "C": "D", "P": "x"}  # P = pre-class (2,3)

# ---------------------------
# Helpers
# ---------------------------
def sieve_primes(n):
    """Return first n primes (simple sieve-ish incremental)."""
    primes = []
    candidate = 2
    while len(primes) < n:
        is_prime = True
        r = int(math.isqrt(candidate))
        for p in primes:
            if p > r:
                break
            if candidate % p == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(candidate)
        candidate += 1 if candidate == 2 else 2  # after 2, skip evens
    return np.array(primes, dtype=np.int64)

def class_mod12(p):
    """E/A/B/C by residue mod 12 (for p>3). 2,3 -> 'P'."""
    if p in (2, 3):
        return "P"
    r = p % 12
    if r == 1:
        return "E"
    if r == 5:
        return "A"
    if r == 7:
        return "B"
    if r == 11:
        return "C"
    # Should not happen for primes > 3
    return "P"

def wrap_to_pi(x):
    """Wrap angle to [-pi, pi]."""
    return (x + math.pi) % (2 * math.pi) - math.pi

def normalize01(x, xmin, xmax):
    if xmax <= xmin:
        return 0.0
    return (x - xmin) / (xmax - xmin)

def hopf_map_to_s2(u):
    """Hopf map: u=(x1,x2,x3,x4) on S^3 -> (X,Y,Z) on S^2."""
    x1, x2, x3, x4 = u
    X = 2.0 * (x1 * x3 + x2 * x4)
    Y = 2.0 * (x2 * x3 - x1 * x4)
    Z = (x1 * x1 + x2 * x2) - (x3 * x3 + x4 * x4)
    return np.array([X, Y, Z], dtype=np.float64)

def mollweide_forward(lam, phi, R=1.0, eps=1e-12):
    """
    Mollweide projection from (lambda, phi) to (x,y).
    Needs solving: 2θ + sin(2θ) = π sin(phi).
    """
    # Clamp latitude to avoid numerical issues at poles
    phi = float(np.clip(phi, -math.pi/2 + 1e-9, math.pi/2 - 1e-9))
    lam = float(np.clip(lam, -math.pi, math.pi))

    # Newton solve for theta
    target = math.pi * math.sin(phi)
    theta = phi  # good initial guess
    for _ in range(10):
        f = 2.0 * theta + math.sin(2.0 * theta) - target
        fp = 2.0 + 2.0 * math.cos(2.0 * theta)
        if abs(fp) < eps:
            break
        d = f / fp
        theta -= d
        if abs(d) < 1e-12:
            break

    x = R * (2.0 * math.sqrt(2.0) / math.pi) * lam * math.cos(theta)
    y = R * math.sqrt(2.0) * math.sin(theta)
    return x, y

def ellipse_boundary(n=600, R=1.0):
    """Mollweide boundary ellipse: x in [-2sqrt2 R, +2sqrt2 R], y in [-sqrt2 R, +sqrt2 R]."""
    t = np.linspace(0, 2*np.pi, n)
    a = 2.0 * math.sqrt(2.0) * R
    b = math.sqrt(2.0) * R
    return a*np.cos(t), b*np.sin(t)

def load_gammas(path, T):
    """
    Load gamma values from file (one per line).
    If missing, create a synthetic fallback (demo only).
    """
    if os.path.exists(path):
        vals = []
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                vals.append(float(s))
                if len(vals) >= T:
                    break
        if len(vals) < T:
            raise RuntimeError(f"File '{path}' has only {len(vals)} gamma values, need >= {T}.")
        return np.array(vals[:T], dtype=np.float64)

    # Fallback: synthetic monotonically increasing "gamma-like" sequence (NOT real zeta zeros)
    # This is only to let you render the pipeline end-to-end.
    print("Hinweis: 'zeros_gamma.txt' nicht gefunden – verwende synthetische Gamma-Werte (keine echten Zeta-Nullstellen).")
    n = np.arange(1, T+1, dtype=np.float64)
    gam = 2*np.pi * n / np.log(n + 10.0)
    return gam

# ---------------------------
# Main pipeline
# ---------------------------
def build_points(primes, gammas):
    """
    For each t, build a 4D vector v_t, normalize to S^3, Hopf-map to S^2,
    turn into (lambda, phi), then Mollweide (x,y).
    Returns arrays: x[t], y[t], class[t], plus diagnostics.
    """
    T = len(primes)
    classes = [class_mod12(int(p)) for p in primes]

    logp = np.log(primes.astype(np.float64))
    logp_max = float(np.max(logp))
    A = logp / logp_max  # amplitude in (0,1]

    gaps = np.zeros(T, dtype=np.float64)
    gaps[1:] = primes[1:] - primes[:-1]

    # Gap signal (bounded)
    G = np.tanh(ALPHA_GAP * (gaps / np.maximum(logp, 1e-9)))

    # Phase coupling theta_t = gamma_t * log p_t (wrapped)
    theta = np.array([wrap_to_pi(float(gammas[t] * logp[t])) for t in range(T)], dtype=np.float64)

    # Build v_t by class; keep it smooth and bounded
    V = np.zeros((T, 4), dtype=np.float64)
    for t in range(T):
        p = int(primes[t])
        K = classes[t]
        ct = math.cos(theta[t])
        st = math.sin(theta[t])
        at = float(A[t])
        gt = float(G[t])

        if K == "E":
            V[t] = [at*ct, at*st, BETA_GAP*gt, 0.0]
        elif K == "A":
            V[t] = [0.0, at*ct, at*st, BETA_GAP*gt]
        elif K == "B":
            V[t] = [BETA_GAP*gt, 0.0, at*ct, at*st]
        elif K == "C":
            V[t] = [at*st, BETA_GAP*gt, 0.0, at*ct]
        else:
            # p in {2,3} (or fallback): put somewhere unobtrusive
            V[t] = [0.0, 0.0, 0.0, 1.0]

    # Normalize to S^3
    norms = np.linalg.norm(V, axis=1)
    norms = np.where(norms < 1e-12, 1.0, norms)
    U = V / norms[:, None]

    # Hopf map to S^2 -> (X,Y,Z) then (lambda, phi)
    lam = np.zeros(T, dtype=np.float64)
    phi = np.zeros(T, dtype=np.float64)
    for t in range(T):
        X, Y, Z = hopf_map_to_s2(U[t])
        lam[t] = math.atan2(Y, X)
        # clip Z for numerical safety
        Zc = float(np.clip(Z, -1.0, 1.0))
        phi[t] = math.asin(Zc)

    # Mollweide projection
    x = np.zeros(T, dtype=np.float64)
    y = np.zeros(T, dtype=np.float64)
    for t in range(T):
        x[t], y[t] = mollweide_forward(lam[t], phi[t], R=R)

    diag = {
        "logp": logp,
        "A": A,
        "gaps": gaps,
        "G": G,
        "theta": theta,
        "lam": lam,
        "phi": phi,
        "V": V,
        "U": U,
        "norms": norms
    }
    return x, y, classes, diag

def render_frames(x, y, classes, primes, diag, outdir=None):
    outdir = outdir if outdir is not None else OUTDIR
    os.makedirs(outdir, exist_ok=True)
    xb, yb = ellipse_boundary(R=R)

    T = len(primes)
    # pre-split indices by class (for quicker scatter)
    idx_by_class = {k: np.array([i for i,c in enumerate(classes) if c == k], dtype=int)
                    for k in set(classes)}

    # Determine static view limits
    xlim = (-2.1*math.sqrt(2)*R, 2.1*math.sqrt(2)*R)
    ylim = (-1.1*math.sqrt(2)*R, 1.1*math.sqrt(2)*R)

    frame_indices = list(range(0, T, FRAME_STEP))
    n_frames = len(frame_indices)
    for fi, t in enumerate(frame_indices):
        if (fi + 1) % 50 == 0 or fi == 0 or fi == n_frames - 1:
            print(f"  Frame {fi+1}/{n_frames} (t={t+1}) …")
        fig = plt.figure(figsize=(10.5, 6.2), dpi=DPI)
        ax = fig.add_axes([0.05, 0.08, 0.72, 0.86])  # main
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal", adjustable="box")
        ax.axis("off")

        # boundary ellipse
        ax.plot(xb, yb, linewidth=1.2)

        # trail window
        start = max(0, t - TRAIL_W)
        trail_idx = np.arange(start, t+1)

        # fade by age
        ages = (t - trail_idx).astype(np.float64)
        alphas = np.exp(-ages / ALPHA_DECAY)
        alphas = np.clip(alphas, 0.05, 1.0)

        # Plot each class as overlay with different marker.
        # We rely on marker shape + transparency; matplotlib's default colors are fine.
        for K in ["E","A","B","C","P"]:
            if K not in idx_by_class:
                continue
            idx = idx_by_class[K]
            # intersect with trail indices
            mask = (idx >= start) & (idx <= t)
            idx2 = idx[mask]
            if idx2.size == 0:
                continue
            # alpha per point based on age
            a = np.array([alphas[i - start] for i in idx2], dtype=np.float64)

            # size slightly increases with amplitude (optional, subtle)
            sizes = 18.0 + 35.0 * diag["A"][idx2]

            ax.scatter(
                x[idx2], y[idx2],
                s=sizes,
                alpha=a,
                marker=MARKERS.get(K, "o"),
                linewidths=0.6
            )

        # Highlight current point
        ax.scatter([x[t]], [y[t]], s=120, marker="*", linewidths=0.8)

        # Sidebar text (diagnostics)
        axr = fig.add_axes([0.79, 0.08, 0.20, 0.86])
        axr.axis("off")

        p = int(primes[t])
        K = classes[t]
        gap = int(diag["gaps"][t]) if t > 0 else 0
        theta = float(diag["theta"][t])
        lam = float(diag["lam"][t])
        phi = float(diag["phi"][t])

        txt = []
        txt.append(f"Frame: {fi+1}/{len(frame_indices)}")
        txt.append(f"t = {t+1}")
        txt.append(f"Prime p_t = {p}")
        txt.append(f"Klasse = {K}")
        txt.append("")
        txt.append(f"Δp = {gap}")
        txt.append(f"A = log(p)/log(p1000) = {diag['A'][t]:.4f}")
        txt.append(f"G = tanh(...) = {diag['G'][t]:.4f}")
        txt.append("")
        txt.append("Phase / Kugel:")
        txt.append(f"θ = wrap(γ·log p) = {theta:.4f}")
        txt.append(f"λ = atan2(Y,X)    = {lam:.4f}")
        txt.append(f"φ = asin(Z)       = {phi:.4f}")
        txt.append("")
        txt.append("Klassen (mod 12):")
        txt.append("E≡1 ○  A≡5 □  B≡7 △  C≡11 ◇  P=2,3 ×")

        axr.text(0.0, 1.0, "\n".join(txt), va="top", fontsize=10)

        # Save frame
        outpath = os.path.join(outdir, f"frame_{fi:04d}.png")
        fig.savefig(outpath, dpi=DPI)
        plt.close(fig)

def parse_args():
    import argparse
    p = argparse.ArgumentParser(
        description="Mollweide-Visualisierung: Primzahlen nach Klassen mod 12 (Hopf → S² → Mollweide)."
    )
    p.add_argument(
        "-n", "--num-primes",
        type=int,
        default=T,
        metavar="N",
        help=f"Anzahl Primzahlen / Frames (Standard: {T})",
    )
    p.add_argument(
        "-o", "--outdir",
        type=str,
        default=OUTDIR,
        metavar="DIR",
        help=f"Ordner für PNG-Frames (Standard: {OUTDIR})",
    )
    p.add_argument(
        "-g", "--gamma-file",
        type=str,
        default="zeros_gamma.txt",
        metavar="FILE",
        help="Datei mit Gamma-Werten (eine Zeile pro Im(ρ)); fehlt sie, werden synthetische Werte verwendet",
    )
    return p.parse_args()


def main():
    args = parse_args()
    n = args.num_primes
    primes = sieve_primes(n)
    gammas = load_gammas(args.gamma_file, n)

    x, y, classes, diag = build_points(primes, gammas)
    render_frames(x, y, classes, primes, diag, outdir=args.outdir)

    print(f"Done. Frames saved in: {args.outdir}/")
    print("To encode MP4 with ffmpeg (example):")
    print(f"  ffmpeg -framerate 30 -i {args.outdir}/frame_%04d.png -c:v libx264 -pix_fmt yuv420p eabc_overlay.mp4")


if __name__ == "__main__":
    main()