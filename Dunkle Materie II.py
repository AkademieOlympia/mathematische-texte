import os, math
import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# Config
# -------------------------
T_FRAMES = 220
OUTDIR = "frames_bm_dm"
DPI = 140

N_SEEDS = 6000            # number of BM "nodes" (points)
KNN = 6                   # neighbors for filament graph
EDGE_THR = 0.78           # draw edges if DM strength above this quantile
FADE = 220.0

# -------------------------
# Core geometry: S3->S2 (Hopf) + Mollweide
# -------------------------
def wrap_to_pi(x):
    return (x + math.pi) % (2*math.pi) - math.pi

def hopf_map(u):
    x1,x2,x3,x4 = u
    X = 2.0*(x1*x3 + x2*x4)
    Y = 2.0*(x2*x3 - x1*x4)
    Z = (x1*x1 + x2*x2) - (x3*x3 + x4*x4)
    return np.array([X,Y,Z], dtype=np.float64)

def mollweide_forward(lam, phi, R=1.0, eps=1e-12):
    phi = float(np.clip(phi, -math.pi/2 + 1e-9, math.pi/2 - 1e-9))
    lam = float(np.clip(lam, -math.pi, math.pi))
    target = math.pi * math.sin(phi)
    theta = phi
    for _ in range(10):
        f = 2.0*theta + math.sin(2.0*theta) - target
        fp = 2.0 + 2.0*math.cos(2.0*theta)
        if abs(fp) < eps:
            break
        d = f/fp
        theta -= d
        if abs(d) < 1e-12:
            break
    x = R * (2.0*math.sqrt(2.0)/math.pi) * lam * math.cos(theta)
    y = R * math.sqrt(2.0) * math.sin(theta)
    return x,y

# -------------------------
# BM plug-ins (YOU customize)
# -------------------------
def init_seeds(n, seed=4):
    """
    Default: random 4D seeds on R^4 with mild structure.
    Replace this with your EABC/prime/zero-based seed generator if desired.
    """
    rng = np.random.default_rng(seed)
    V = rng.normal(size=(n,4))
    # encourage "EABC-ish" anisotropy a bit
    V[:,0] *= 1.2
    V[:,1] *= 0.9
    V[:,2] *= 1.0
    V[:,3] *= 0.8
    return V

def evolve(V, t, params):
    """
    Toy evolution: phase drift + weak coupling.
    Replace with BM dynamics: pressure term, Okto/Tri weighting, etc.
    """
    a = params["a_max"] * (t/(params["T"]-1))
    # rotate in (e,a) and (b,c) planes with different angular speeds
    w1 = 0.55 + 0.25*math.sin(0.03*t)
    w2 = 0.33 + 0.20*math.cos(0.02*t)

    E,A,B,C = V[:,0], V[:,1], V[:,2], V[:,3]
    c1,s1 = math.cos(w1*a), math.sin(w1*a)
    c2,s2 = math.cos(w2*a), math.sin(w2*a)

    E2 = c1*E - s1*A
    A2 = s1*E + c1*A
    B2 = c2*B - s2*C
    C2 = s2*B + c2*C

    V2 = np.stack([E2,A2,B2,C2], axis=1)

    # mild nonlinear "contraction" to help form bundles
    V2 *= (1.0 / (1.0 + 0.08*np.linalg.norm(V2, axis=1, keepdims=True)))
    return V2

def dm_strength(V, t, params):
    """
    THIS is where your BM 'dark matter' definition goes.

    Output must be in [0,1] (or we normalize it).
    Example toy DM: prefer large 'hidden' energy in (b,c) relative to (e,a).
    """
    ea = np.linalg.norm(V[:,0:2], axis=1)
    bc = np.linalg.norm(V[:,2:4], axis=1)
    raw = bc / np.maximum(ea + bc, 1e-9)
    # squashing to [0,1]
    return np.clip(raw, 0.0, 1.0)

# -------------------------
# Filament extraction: KNN edges
# -------------------------
def knn_edges(XY, k):
    """
    Build undirected edges using brute-force kNN (ok for a few thousand points).
    Returns list of (i,j).
    """
    N = XY.shape[0]
    # pairwise distances
    D = np.sum((XY[:,None,:] - XY[None,:,:])**2, axis=2)
    np.fill_diagonal(D, np.inf)
    nn = np.argpartition(D, kth=k, axis=1)[:, :k]
    edges = set()
    for i in range(N):
        for j in nn[i]:
            a,b = (i,int(j)) if i < j else (int(j),i)
            edges.add((a,b))
    return list(edges)

# -------------------------
# Render
# -------------------------
def render():
    os.makedirs(OUTDIR, exist_ok=True)

    params = {"T": T_FRAMES, "a_max": 1.35}

    V0 = init_seeds(N_SEEDS)
    xb = np.linspace(0, 2*np.pi, 700)
    a = 2.0*math.sqrt(2.0)
    b = math.sqrt(2.0)
    ex = a*np.cos(xb)
    ey = b*np.sin(xb)

    for t in range(T_FRAMES):
        Vt = evolve(V0, t, params)
        # normalize to S^3
        nrm = np.linalg.norm(Vt, axis=1, keepdims=True)
        U = Vt / np.maximum(nrm, 1e-12)

        # Hopf -> (lambda,phi)
        XYZ = np.array([hopf_map(U[i]) for i in range(U.shape[0])])
        lam = np.arctan2(XYZ[:,1], XYZ[:,0])
        phi = np.arcsin(np.clip(XYZ[:,2], -1.0, 1.0))

        # Mollweide -> (x,y)
        XY = np.zeros((U.shape[0],2), dtype=np.float64)
        for i in range(U.shape[0]):
            XY[i,0], XY[i,1] = mollweide_forward(lam[i], phi[i], R=1.0)

        # DM strength + threshold
        Dm = dm_strength(Vt, t, params)
        thr = np.quantile(Dm, EDGE_THR)
        keep = np.where(Dm >= thr)[0]

        # edges on "high DM" subset
        XYk = XY[keep]
        edges = knn_edges(XYk, KNN)

        # plot
        fig = plt.figure(figsize=(11.5, 6.5), dpi=DPI)
        ax = fig.add_axes([0.05, 0.08, 0.92, 0.86])
        ax.set_aspect("equal", adjustable="box")
        ax.axis("off")
        ax.plot(ex, ey, linewidth=1.2)

        # draw edges
        for (i,j) in edges:
            p = XYk[i]; q = XYk[j]
            ax.plot([p[0], q[0]], [p[1], q[1]], linewidth=0.6, alpha=0.30)

        # draw nodes (faint background + strong DM nodes)
        ax.scatter(XY[:,0], XY[:,1], s=6, alpha=0.08)
        ax.scatter(XYk[:,0], XYk[:,1], s=9, alpha=0.40)

        ax.text(-2.75, 1.25, f"BM-DM Threads | frame {t+1}/{T_FRAMES}", fontsize=12)

        outpath = os.path.join(OUTDIR, f"frame_{t:04d}.png")
        fig.savefig(outpath, dpi=DPI)
        plt.close(fig)

    print(f"Frames written to {OUTDIR}/")
    print("Encode example:")
    print("  ffmpeg -framerate 30 -i frames_bm_dm/frame_%04d.png -c:v libx264 -pix_fmt yuv420p bm_dm_threads.mp4")

if __name__ == "__main__":
    render()