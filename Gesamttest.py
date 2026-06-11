import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange

# -----------------------------
# 0) Z2xZ2-Ladungen via mod 12 (für p>3 exakt)
# -----------------------------
def chi4(p):
    return 1 if (p % 4) == 1 else -1

def chi3(p):
    r = p % 12
    return 1 if r in (1,7) else -1

# -----------------------------
# 1) Build lab (Chebyshev-theta Bias + Gauge)
# -----------------------------
def build_lab(limit=200000, alpha=np.pi/60, beta=np.pi/60, window=5000, Ls=(120,240)):
    primes = [p for p in primerange(1, limit) if p > 3]
    r = np.array([p % 12 for p in primes], dtype=int)
    logp = np.log(np.array(primes, dtype=float))

    isE = (r == 1).astype(int)
    isA = (r == 5).astype(int)
    isB = (r == 7).astype(int)
    isC = (r == 11).astype(int)

    thetaE = np.cumsum(logp * isE)
    thetaA = np.cumsum(logp * isA)
    thetaB = np.cumsum(logp * isB)
    thetaC = np.cumsum(logp * isC)

    I = thetaE
    V = (thetaA + thetaB + thetaC)/3.0 - thetaE

    A = np.empty(len(primes), dtype=float)
    for i,p in enumerate(primes):
        A[i] = alpha * 0.5*(1 - chi4(p)) + beta * 0.5*(1 - chi3(p))

    theta = np.cumsum(A)
    W = np.exp(1j * theta)

    ang = np.unwrap(np.angle(W))
    winding = ang / (2*np.pi)

    loops = {}
    for L in Ls:
        if len(A) > L:
            loop_phase = np.convolve(A, np.ones(L), mode='valid')
            loops[L] = np.exp(1j * loop_phase)
        else:
            loops[L] = np.array([], dtype=complex)

    if len(primes) > window:
        k = np.ones(window)/window
        wind_s = np.convolve(winding, k, mode='valid')
        V_s = np.convolve(V, k, mode='valid')
        dw = np.gradient(wind_s)
        dV = np.gradient(V_s)
        sigma_top = np.where(np.abs(dV) > 1e-12, dw/dV, 0.0)
        x_axis = np.array(primes[window-1:])
    else:
        x_axis = np.array([])
        sigma_top = np.array([])

    return {
        "primes": np.array(primes),
        "r": r,
        "I": I,
        "V": V,
        "A": A,
        "W": W,
        "winding": winding,
        "loops": loops,
        "x_axis": x_axis,
        "sigma_top": sigma_top
    }

# -----------------------------
# 2) Robust helpers: plateaus + block permutation
# -----------------------------
def plateau_segments(y, smooth=2000, slope_thr=2e-4, min_len=2000):
    y = np.asarray(y, float)
    if len(y) < max(smooth, min_len) + 10:
        return [], np.array([]), np.array([])

    k = np.ones(smooth)/smooth
    ys = np.convolve(y, k, mode='valid')
    slope = np.gradient(ys)

    mask = np.abs(slope) < slope_thr
    segs = []
    i = 0
    while i < len(mask):
        if mask[i]:
            j = i
            while j < len(mask) and mask[j]:
                j += 1
            if (j - i) >= min_len:
                segs.append((i, j))
            i = j
        else:
            i += 1
    return segs, ys, slope

def plateau_score(segs, ys):
    if not segs or len(ys)==0:
        return 0.0
    score = 0.0
    for a,b in segs:
        v = np.var(ys[a:b])
        score += (b-a) * (1.0/(1.0+v))
    return float(score)

def block_permute(x, block=2000, rng=None):
    x = np.asarray(x)
    n = len(x)
    if rng is None:
        rng = np.random.default_rng()
    idx = np.arange(n)
    blocks = [idx[i:i+block] for i in range(0, n, block)]
    rng.shuffle(blocks)
    perm = np.concatenate(blocks)
    return x[perm]

def pvalue_plateau_block(y, n_perm=120, block=2000, smooth=2000, slope_thr=2e-4, min_len=2000, seed=0):
    segs, ys, _ = plateau_segments(y, smooth=smooth, slope_thr=slope_thr, min_len=min_len)
    obs = plateau_score(segs, ys)

    rng = np.random.default_rng(seed)
    sims = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        yp = block_permute(y, block=block, rng=rng)
        segs_p, ys_p, _ = plateau_segments(yp, smooth=smooth, slope_thr=slope_thr, min_len=min_len)
        sims[i] = plateau_score(segs_p, ys_p)

    p = (np.sum(sims >= obs) + 1) / (n_perm + 1)
    return obs, p

# -----------------------------
# 3) Resonanz p-value (120/240) via block permutation of A
# -----------------------------
def resonance_stat(loop_complex):
    if len(loop_complex) == 0:
        return 0.0
    return float(np.abs(np.mean(loop_complex)))

def pvalue_resonance_block(A, L, n_perm=180, block=2000, seed=1):
    A = np.asarray(A, float)
    if len(A) <= L + 10:
        return 0.0, 1.0

    loop_phase = np.convolve(A, np.ones(L), mode='valid')
    Wobs = np.exp(1j * loop_phase)
    obs = resonance_stat(Wobs)

    rng = np.random.default_rng(seed)
    sims = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        Ap = block_permute(A, block=block, rng=rng)
        lp = np.convolve(Ap, np.ones(L), mode='valid')
        sims[i] = resonance_stat(np.exp(1j * lp))

    p = (np.sum(sims >= obs) + 1) / (n_perm + 1)
    return obs, p

# -----------------------------
# 4) Nichtlokalität: E_{k+1} aus globaler Phase vs lokale Z2-bits
# -----------------------------
def logistic_fit_params(X, y, l2=1e-2, iters=1200, lr=0.06):
    X = np.asarray(X, float); y = np.asarray(y, float)
    n,d = X.shape
    w = np.zeros(d); b = 0.0
    for _ in range(iters):
        z = X @ w + b
        p = 1/(1+np.exp(-z))
        gw = (X.T @ (p-y))/n + l2*w
        gb = np.mean(p-y)
        w -= lr*gw; b -= lr*gb
    return w,b

def logistic_predict(X, w, b):
    z = np.asarray(X,float) @ w + b
    return 1/(1+np.exp(-z))

def logloss(y, p, eps=1e-12):
    y = np.asarray(y)
    p = np.clip(np.asarray(p), eps, 1-eps)
    return float(-np.mean(y*np.log(p) + (1-y)*np.log(1-p)))

def nonlocality_test_block(primes, r, W, train_frac=0.7, n_perm=160, block=2000, seed=2):
    isE = (r==1).astype(int)
    y = isE[1:]
    phi = np.angle(W[:-1])

    Xg = np.column_stack([np.cos(phi), np.sin(phi)])
    chi4_bit = np.array([(0 if chi4(p)==1 else 1) for p in primes[:-1]], dtype=float)
    chi3_bit = np.array([(0 if chi3(p)==1 else 1) for p in primes[:-1]], dtype=float)
    Xl = np.column_stack([chi4_bit, chi3_bit])

    n = len(y)
    ntr = int(train_frac*n)

    wl, bl = logistic_fit_params(Xl[:ntr], y[:ntr])
    wg, bg = logistic_fit_params(Xg[:ntr], y[:ntr])

    pl = logistic_predict(Xl[ntr:], wl, bl)
    pg = logistic_predict(Xg[ntr:], wg, bg)

    ll_l = logloss(y[ntr:], pl)
    ll_g = logloss(y[ntr:], pg)
    delta = ll_l - ll_g  # >0 => global besser

    rng = np.random.default_rng(seed)
    sims = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        phip = block_permute(phi, block=block, rng=rng)
        Xgp = np.column_stack([np.cos(phip), np.sin(phip)])
        wg_p, bg_p = logistic_fit_params(Xgp[:ntr], y[:ntr])
        pg_p = logistic_predict(Xgp[ntr:], wg_p, bg_p)
        ll_g_p = logloss(y[ntr:], pg_p)
        sims[i] = ll_l - ll_g_p

    p = (np.sum(sims >= delta) + 1) / (n_perm + 1)
    return delta, p, ll_l, ll_g

# -----------------------------
# 5) FDR (Benjamini–Hochberg)
# -----------------------------
def bh_fdr(pvals):
    pvals = np.asarray(pvals, float)
    m = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q_full = np.empty_like(q)
    q_full[order] = np.clip(q, 0, 1)
    return q_full

# -----------------------------
# 6) The scan
# -----------------------------
def scan_alpha_beta(
    limit=200000,
    alpha_grid=None,
    beta_grid=None,
    perms_fast=True
):
    if alpha_grid is None:
        # typischer Bereich: "ein paar Grad" pro Schritt; hier in Radiant
        alpha_grid = np.linspace(0, np.pi/20, 13)   # 0 ... 9°
    if beta_grid is None:
        beta_grid  = np.linspace(0, np.pi/20, 13)

    # Permutationsbudget: "fast" vs "robust"
    if perms_fast:
        nP_plateau = 120
        nP_res     = 180
        nP_nonloc  = 160
    else:
        nP_plateau = 300
        nP_res     = 500
        nP_nonloc  = 400

    # fixe Settings
    block = 2000
    smooth = 2000
    slope_thr = 2e-4
    min_len = 2000
    window = 5000

    results = []  # (alpha,beta, p_plateau,p120,p240, p_nonloc, delta, obs_plateau, stat120,stat240, ll_l,ll_g)

    for ia,alpha in enumerate(alpha_grid):
        for ib,beta in enumerate(beta_grid):
            lab = build_lab(limit=limit, alpha=alpha, beta=beta, window=window, Ls=(120,240))

            obs_pl, p_pl = pvalue_plateau_block(
                lab["winding"], n_perm=nP_plateau, block=block,
                smooth=smooth, slope_thr=slope_thr, min_len=min_len,
                seed=1000 + ia*100 + ib
            )

            stat120, p120 = pvalue_resonance_block(
                lab["A"], L=120, n_perm=nP_res, block=block, seed=2000 + ia*100 + ib
            )
            stat240, p240 = pvalue_resonance_block(
                lab["A"], L=240, n_perm=nP_res, block=block, seed=3000 + ia*100 + ib
            )

            delta, p_nonloc, ll_l, ll_g = nonlocality_test_block(
                lab["primes"], lab["r"], lab["W"],
                n_perm=nP_nonloc, block=block, seed=4000 + ia*100 + ib
            )

            results.append((alpha, beta, p_pl, p120, p240, p_nonloc, delta, obs_pl, stat120, stat240, ll_l, ll_g))

            print(f"alpha={alpha:.4f} beta={beta:.4f} | "
                  f"pP={p_pl:.3f} p120={p120:.3f} p240={p240:.3f} pNL={p_nonloc:.3f} Δ={delta:+.4e}")

    results = np.array(results, dtype=float)
    # columns
    # 0 alpha,1 beta,2 p_pl,3 p120,4 p240,5 p_nonloc,6 delta,7 obs_pl,8 stat120,9 stat240,10 ll_l,11 ll_g

    # FDR pro Test separat (konservativ, aber sauber)
    q_pl    = bh_fdr(results[:,2])
    q_120   = bh_fdr(results[:,3])
    q_240   = bh_fdr(results[:,4])
    q_nonlc = bh_fdr(results[:,5])

    # "Gesamt"-Kriterium: max(qs) klein UND delta>0
    q_max = np.maximum.reduce([q_pl, q_120, q_240, q_nonlc])
    good = (q_max < 0.10) & (results[:,6] > 0)  # FDR 10% als Start

    # Top nach q_max, dann nach delta
    order = np.lexsort((-results[:,6], q_max))
    top = results[order][:20]
    top_q = q_max[order][:20]

    print("\n=== TOP-20 (nach q_max dann Δ) ===")
    for i in range(len(top)):
        a,b,ppl,p120,p240,pnl,delta,obspl,st120,st240,ll_l,ll_g = top[i]
        print(f"{i+1:02d}) alpha={a:.5f} beta={b:.5f} | "
              f"qmax={top_q[i]:.3f} | "
              f"pP={ppl:.3f} p120={p120:.3f} p240={p240:.3f} pNL={pnl:.3f} | Δ={delta:+.3e}")

    # Heatmaps (p oder q)
    na = len(alpha_grid); nb = len(beta_grid)
    def gridify(col):
        return results[:,col].reshape(na, nb)

    P_pl  = gridify(2)
    P_120 = gridify(3)
    P_240 = gridify(4)
    P_nl  = gridify(5)
    Dlt   = gridify(6)

    Qmax  = q_max.reshape(na, nb)

    fig, axs = plt.subplots(2, 3, figsize=(14, 8))
    im0 = axs[0,0].imshow(P_pl,  origin='lower', aspect='auto')
    axs[0,0].set_title("p(Plateau)")
    plt.colorbar(im0, ax=axs[0,0], fraction=0.046)

    im1 = axs[0,1].imshow(P_120, origin='lower', aspect='auto')
    axs[0,1].set_title("p(Resonanz 120)")
    plt.colorbar(im1, ax=axs[0,1], fraction=0.046)

    im2 = axs[0,2].imshow(P_240, origin='lower', aspect='auto')
    axs[0,2].set_title("p(Resonanz 240)")
    plt.colorbar(im2, ax=axs[0,2], fraction=0.046)

    im3 = axs[1,0].imshow(P_nl,  origin='lower', aspect='auto')
    axs[1,0].set_title("p(Nichtlokalität)")
    plt.colorbar(im3, ax=axs[1,0], fraction=0.046)

    im4 = axs[1,1].imshow(Dlt,   origin='lower', aspect='auto')
    axs[1,1].set_title("Δlogloss (local - global)")
    plt.colorbar(im4, ax=axs[1,1], fraction=0.046)

    im5 = axs[1,2].imshow(Qmax,  origin='lower', aspect='auto')
    axs[1,2].set_title("q_max (FDR) über alle 4 Tests")
    plt.colorbar(im5, ax=axs[1,2], fraction=0.046)

    for ax in axs.ravel():
        ax.set_xlabel("beta index")
        ax.set_ylabel("alpha index")

    plt.tight_layout()
    plt.show()

    return {
        "results": results,
        "alpha_grid": alpha_grid,
        "beta_grid": beta_grid,
        "q_pl": q_pl, "q_120": q_120, "q_240": q_240, "q_nonloc": q_nonlc,
        "q_max": q_max,
        "good_mask": good
    }

# -----------------------------
# RUN SCAN
# -----------------------------
scan_out = scan_alpha_beta(
    limit=200000,
    alpha_grid=np.
    
    linspace(0, np.pi/20, 13),
    beta_grid=np.linspace(0, np.pi/20, 13),
    perms_fast=True
)
