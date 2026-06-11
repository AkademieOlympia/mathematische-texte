import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange

# -----------------------------
# 0) Z2xZ2-Ladungen via mod 12 (für p>3 exakt)
# -----------------------------
def chi4(p):
    # (-1/p): +1 für p ≡ 1 mod 4, -1 für p ≡ 3 mod 4
    return 1 if (p % 4) == 1 else -1

def chi3(p):
    # (-3/p) über mod 12 Klassen (für p>3):
    r = p % 12
    return 1 if r in (1,7) else -1

# -----------------------------
# 1) Build lab v2:
#    - Bias V als Chebyshev-theta
#    - Strom I als theta_E
#    - Gauge-Connection A(p) = (alpha*q4 + beta*q3) * w(p)
#      mit w(p)=log p (default), damit "Flux" kompatibel zu theta
# -----------------------------
def build_lab_v2(limit=200000, alpha=0.02, beta=0.02, window=5000, Ls=(120,240),
                weight="log", gamma=1.0):
    primes = [p for p in primerange(1, limit) if p > 3]
    r = np.array([p % 12 for p in primes], dtype=int)
    logp = np.log(np.array(primes, dtype=float))

    isE = (r == 1).astype(int)
    isA = (r == 5).astype(int)
    isB = (r == 7).astype(int)
    isC = (r == 11).astype(int)

    # Chebyshev theta je Klasse
    thetaE = np.cumsum(logp * isE)
    thetaA = np.cumsum(logp * isA)
    thetaB = np.cumsum(logp * isB)
    thetaC = np.cumsum(logp * isC)

    # Observablen
    I = thetaE
    V = (thetaA + thetaB + thetaC)/3.0 - thetaE

    # Gewichte für Connection
    if weight == "log":
        w = logp
    elif weight == "sqrtlog":
        w = np.sqrt(logp)
    elif weight == "ones":
        w = np.ones_like(logp)
    else:
        raise ValueError("weight must be one of: log, sqrtlog, ones")

    # Gauge-Connection A(p)
    # q4,q3 ∈ {0,1} aktivieren Phase je Z2-Komponente
    A = np.empty(len(primes), dtype=float)
    for i,p in enumerate(primes):
        q4 = 0.5*(1 - chi4(p))
        q3 = 0.5*(1 - chi3(p))
        A[i] = (alpha*q4 + beta*q3) * (w[i]**gamma)

    # Wilson line / Holonomie
    theta = np.cumsum(A)
    W = np.exp(1j * theta)

    # Winding
    ang = np.unwrap(np.angle(W))
    winding = ang / (2*np.pi)

    # Wilson loops über L
    loops = {}
    for L in Ls:
        if len(A) > L:
            loop_phase = np.convolve(A, np.ones(L), mode='valid')
            loops[L] = np.exp(1j * loop_phase)
        else:
            loops[L] = np.array([], dtype=complex)

    # sigma_top ~ d(winding)/dV geglättet
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
        "logp": logp,
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
# 2) Block-Permutation
# -----------------------------
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

# -----------------------------
# 3) Plateau-Test (relaxed) auf winding
# -----------------------------
def plateau_segments(y, smooth=1500, slope_thr=1e-3, min_len=1200):
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

def pvalue_plateau_block(y, n_perm=160, block=2000, seed=0,
                         smooth=1500, slope_thr=1e-3, min_len=1200):
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
# 4) Quantisierungs-Test auf sigma_top:
#    Idee: Plateaus zeigen sich als "stabile Levels" -> Peaks im Histogramm.
#    Statistik: max_peak_height (normalized) gegenüber Block-Permutation von sigma_top
# -----------------------------
def quantization_stat_sigma(sigma, bins=80):
    sigma = np.asarray(sigma, float)
    if len(sigma) < 1000:
        return 0.0
    # robust trimming
    lo, hi = np.quantile(sigma, [0.01, 0.99])
    s = sigma[(sigma >= lo) & (sigma <= hi)]
    if len(s) < 1000:
        return 0.0
    hist, _ = np.histogram(s, bins=bins, density=True)
    # "Peakiness" = maximaler Dichte-Peak
    return float(np.max(hist))

def pvalue_quantization_sigma_block(sigma, n_perm=200, block=2000, seed=3, bins=80):
    obs = quantization_stat_sigma(sigma, bins=bins)
    rng = np.random.default_rng(seed)
    sims = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        sp = block_permute(sigma, block=block, rng=rng)
        sims[i] = quantization_stat_sigma(sp, bins=bins)
    p = (np.sum(sims >= obs) + 1) / (n_perm + 1)
    return obs, p

# -----------------------------
# 5) Resonanz-Test (120/240): circular mean magnitude
# -----------------------------
def resonance_stat(loop_complex):
    if len(loop_complex) == 0:
        return 0.0
    return float(np.abs(np.mean(loop_complex)))

def pvalue_resonance_block(A, L, n_perm=240, block=2000, seed=1):
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
# 6) Nichtlokalitäts-Test (AB-Information):
#    Ziel: E_{k+1} vorhersagen
#    Local features: (chi4,chi3) bits
#    Global features: (cos phi, sin phi)
# -----------------------------
def logistic_fit_params(X, y, l2=1e-2, iters=1400, lr=0.06):
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

def nonlocality_test_block(primes, r, W, train_frac=0.7, n_perm=220, block=2000, seed=2):
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
# 7) FDR BH
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
# 8) Scan v2:
#    Tests:
#      p_plateau_wind
#      p_quant_sigma
#      p_res_120
#      p_res_240
#      p_nonlocal + delta
#    Gesamt: q_max über 5 p-Werte; zusätzlich delta>0 als AB-Kriterium
# -----------------------------
def scan_v2(
    limit=200000,
    alpha_grid=None,
    beta_grid=None,
    weight="log",
    gamma=1.0,
    perms_fast=True
):
    if alpha_grid is None:
        alpha_grid = np.linspace(0.002, 0.12, 13)  # absichtlich >0
    if beta_grid is None:
        beta_grid  = np.linspace(0.002, 0.12, 13)

    if perms_fast:
        nP_pl   = 160
        nP_q    = 200
        nP_res  = 240
        nP_nl   = 220
    else:
        nP_pl   = 350
        nP_q    = 450
        nP_res  = 600
        nP_nl   = 500

    block = 2000
    window = 5000

    results = []
    # cols:
    # alpha,beta, p_pl, p_q, p120, p240, p_nl, delta, obs_pl, obs_q, stat120, stat240, ll_l,ll_g
    for ia,a in enumerate(alpha_grid):
        for ib,b in enumerate(beta_grid):
            lab = build_lab_v2(limit=limit, alpha=a, beta=b, window=window, Ls=(120,240),
                               weight=weight, gamma=gamma)

            obs_pl, p_pl = pvalue_plateau_block(
                lab["winding"], n_perm=nP_pl, block=block, seed=1000+ia*100+ib,
                smooth=1500, slope_thr=1e-3, min_len=1200
            )

            obs_q, p_q = pvalue_quantization_sigma_block(
                lab["sigma_top"], n_perm=nP_q, block=block, seed=2000+ia*100+ib, bins=80
            )

            stat120, p120 = pvalue_resonance_block(
                lab["A"], L=120, n_perm=nP_res, block=block, seed=3000+ia*100+ib
            )
            stat240, p240 = pvalue_resonance_block(
                lab["A"], L=240, n_perm=nP_res, block=block, seed=4000+ia*100+ib
            )

            delta, p_nl, ll_l, ll_g = nonlocality_test_block(
                lab["primes"], lab["r"], lab["W"],
                n_perm=nP_nl, block=block, seed=5000+ia*100+ib
            )

            results.append((a,b, p_pl,p_q,p120,p240,p_nl, delta, obs_pl,obs_q, stat120,stat240, ll_l,ll_g))
            print(f"a={a:.4f} b={b:.4f} | pP={p_pl:.3f} pQ={p_q:.3f} p120={p120:.3f} p240={p240:.3f} pNL={p_nl:.3f} Δ={delta:+.3e}")

    results = np.array(results, dtype=float)

    # FDR pro p-Spalte
    q_pl  = bh_fdr(results[:,2])
    q_q   = bh_fdr(results[:,3])
    q_120 = bh_fdr(results[:,4])
    q_240 = bh_fdr(results[:,5])
    q_nl  = bh_fdr(results[:,6])

    q_max = np.maximum.reduce([q_pl,q_q,q_120,q_240,q_nl])

    # "physikalischer Treffer" = delta>0 UND q_max klein
    good = (results[:,7] > 0) & (q_max < 0.10)

    order = np.lexsort((-results[:,7], q_max))
    top = results[order][:25]
    top_q = q_max[order][:25]

    print("\n=== TOP-25 (nach q_max dann Δ) ===")
    for i,row in enumerate(top):
        a,b,ppl,pq,p120,p240,pnl,delta,obspl,obsq,st120,st240,ll_l,ll_g = row
        print(f"{i+1:02d}) a={a:.5f} b={b:.5f} | qmax={top_q[i]:.3f} | "
              f"pP={ppl:.3f} pQ={pq:.3f} p120={p120:.3f} p240={p240:.3f} pNL={pnl:.3f} | Δ={delta:+.3e}")

    # Heatmap q_max & delta
    na, nb = len(alpha_grid), len(beta_grid)
    Qmax = q_max.reshape(na, nb)
    Dlt  = results[:,7].reshape(na, nb)

    fig, axs = plt.subplots(1,2, figsize=(12,4))
    im0 = axs[0].imshow(Qmax, origin="lower", aspect="auto")
    axs[0].set_title("q_max (FDR) über 5 Tests")
    axs[0].set_xlabel("beta index"); axs[0].set_ylabel("alpha index")
    plt.colorbar(im0, ax=axs[0], fraction=0.046)

    im1 = axs[1].imshow(Dlt, origin="lower", aspect="auto")
    axs[1].set_title("Δlogloss (local - global)")
    axs[1].set_xlabel("beta index"); axs[1].set_ylabel("alpha index")
    plt.colorbar(im1, ax=axs[1], fraction=0.046)

    plt.tight_layout()
    plt.show()

    print(f"\n#Good (Δ>0 & qmax<0.10): {np.sum(good)} von {len(good)}")

    return {
        "results": results,
        "alpha_grid": alpha_grid,
        "beta_grid": beta_grid,
        "q_pl": q_pl, "q_q": q_q, "q_120": q_120, "q_240": q_240, "q_nl": q_nl,
        "q_max": q_max,
        "good_mask": good
    }

# -----------------------------
# RUN: Scan v2 (empfohlen)
# -----------------------------
out = scan_v2(
    limit=200000,
    alpha_grid=np.linspace(0.002, 0.12, 13),
    beta_grid=np.linspace(0.002, 0.12, 13),
    weight="log",
    gamma=1.0,
    perms_fast=True
)
