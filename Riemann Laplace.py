import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# ---------- Riemann zeros ----------
zeros = np.load("zeros6.npy").astype(float)  # your file
# pick a slice to avoid startup transients
Z = zeros[:200000]

def N_riemann(T):
    # Riemann-von Mangoldt counting function (smooth approx)
    # N(T) ~ (T/2π) log(T/2π) - (T/2π) + 7/8
    T = np.asarray(T)
    return (T/(2*np.pi)) * (np.log(T/(2*np.pi)) - 1.0) + 7.0/8.0

tZ = N_riemann(Z)
sZ = np.diff(tZ)
sZ = sZ / np.mean(sZ)  # normalize mean spacing to 1

# ---------- Tetrahedron graph Laplacian ----------
def tetra_points(n):
    pts = []
    idx = {}
    for i in range(n+1):
        for j in range(n+1-i):
            for k in range(n+1-i-j):
                idx[(i,j,k)] = len(pts)
                pts.append((i,j,k))
    return pts, idx

def laplacian_tetra(n):
    pts, idx = tetra_points(n)
    N = len(pts)
    dirs = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
    rows, cols, data = [], [], []
    deg = np.zeros(N)

    for a,(i,j,k) in enumerate(pts):
        for di,dj,dk in dirs:
            p = (i+di, j+dj, k+dk)
            b = idx.get(p, None)
            if b is not None:
                rows.append(a); cols.append(b); data.append(-1.0)
                deg[a] += 1.0

    rows.extend(range(N)); cols.extend(range(N)); data.extend(deg.tolist())
    L = sp.csr_matrix((data,(rows,cols)), shape=(N,N))
    return L

def unfolded_spacings_from_eigs(eigs, poly_deg=5):
    eigs = np.sort(eigs)
    # drop the 0 eigenvalue
    eigs = eigs[eigs > 1e-9]
    k = np.arange(1, len(eigs)+1)
    # fit smooth counting function k ≈ P(λ)
    P = np.polyfit(eigs, k, deg=poly_deg)
    Ns = np.polyval(P, eigs)
    s = np.diff(Ns)
    s = s / np.mean(s)
    return s

# choose n so matrix is manageable; increase for more eigenvalues
n = 18
L = laplacian_tetra(n)

# compute smallest eigenvalues (excluding zero) using shift-invert near 0
k_eigs = 800
vals = spla.eigsh(L, k=k_eigs, sigma=1e-8, which="LM",
                 return_eigenvectors=False, tol=1e-7, maxiter=5000)
vals = np.sort(vals)

sL = unfolded_spacings_from_eigs(vals, poly_deg=5)

# ---------- Plot comparison ----------
def wigner_gue_pdf(s):
    return (32.0/np.pi**2) * s**2 * np.exp(-4.0*s**2/np.pi)

bins = np.linspace(0, 3, 80)
centers = 0.5*(bins[1:]+bins[:-1])

hZ, _ = np.histogram(sZ[(sZ>=0)&(sZ<=3)], bins=bins, density=True)
hL, _ = np.histogram(sL[(sL>=0)&(sL<=3)], bins=bins, density=True)
wg = wigner_gue_pdf(centers)

plt.figure()
plt.plot(centers, hZ, label="Riemann zeros (unfolded)")
plt.plot(centers, hL, label=f"Tetra Laplacian n={n} (unfolded)")
plt.plot(centers, wg, label="GUE Wigner surmise")
plt.xlabel("s")
plt.ylabel("density")
plt.legend()
plt.title("Nearest-neighbor spacing comparison")
plt.show()

# quick diagnostics
print("mean spacings:", np.mean(sZ), np.mean(sL))
print("var spacings :", np.var(sZ), np.var(sL))