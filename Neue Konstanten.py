import itertools, numpy as np, math

phi = (1 + 5**0.5)/2

# 12 vertices in standard position
verts = []
for s1 in (-1,1):
    for s2 in (-1,1):
        verts.append((0, s1, s2*phi))
        verts.append((s1, s2*phi, 0))
        verts.append((s1*phi, 0, s2))
V = []
for v in verts:
    if v not in V:
        V.append(v)
V = np.array(V, float)
assert len(V)==12

def angle(u,v):
    cuv = float(np.dot(u,v) / (np.linalg.norm(u)*np.linalg.norm(v)))
    cuv = max(-1.0, min(1.0, cuv))
    return math.acos(cuv)

def tet_volume_idx(idx):
    A,B,C,D = V[list(idx)]
    return abs(np.linalg.det(np.vstack([B-A, C-A, D-A]).T))/6.0

def tet_deficits_idx(idx):
    A,B,C,D = V[list(idx)]
    P = [A,B,C,D]
    deficits = []
    for i in range(4):
        a = P[i]
        others = [P[j] for j in range(4) if j!=i]
        v0,v1,v2 = [o-a for o in others]
        th01 = angle(v0,v1)
        th02 = angle(v0,v2)
        th12 = angle(v1,v2)
        delta = 2*math.pi - (th01+th02+th12)
        deficits.append(delta)
    return np.array(deficits)

def D_curv_idx(idx):
    d = tet_deficits_idx(idx)
    return float(np.sum((d-math.pi)**2)), np.sort(d)  # (energy, signature)

# enumerate all nondegenerate tets
eps = 1e-9
tets = []
for comb in itertools.combinations(range(12), 4):
    if tet_volume_idx(comb) > eps:
        tets.append(tuple(comb))
assert len(tets) == 420

tset = [set(t) for t in tets]

# precompute geometric invariants
vol = np.array([tet_volume_idx(t) for t in tets])
Dval = np.zeros(len(tets))
Dsig = np.zeros((len(tets),4))
for i,t in enumerate(tets):
    dv, ds = D_curv_idx(t)
    Dval[i] = dv
    Dsig[i,:] = ds

Vref = float(np.max(vol))  # or np.mean(vol)

def P_press(Vt, p):
    return p*(Vref - Vt)**2

def E_geom(i, p, alpha=1.0, beta=1.0):
    return alpha*Dval[i] + beta*P_press(vol[i], p)

def pair_energy(i, j, p,
                alpha=1.0, beta=1.0,
                lamV=0.2, lamD=0.2):
    # E_pair(T1,T2;p) = E_geom(T2;p) + mismatch penalties w.r.t T1
    Eg2 = E_geom(j, p, alpha, beta)
    Vm = (vol[j] - vol[i])**2
    Dm = float(np.sum((Dsig[j]-Dsig[i])**2))
    return Eg2 + lamV*beta*Vm + lamD*alpha*Dm

def soft_count_omega(i, p, disjoint_js, q=0.6, tau_scale=0.02,
                    alpha=1.0, beta=1.0, lamV=0.2, lamD=0.2):
    Es = np.array([pair_energy(i,j,p,alpha,beta,lamV,lamD) for j in disjoint_js], float)
    # adaptive cutoff as quantile
    Ecut = float(np.quantile(Es, q))
    # adaptive tau
    s = float(np.std(Es)) + 1e-12
    tau = tau_scale * s
    w = 1.0/(1.0 + np.exp((Es - Ecut)/tau))
    return float(np.sum(w)), Ecut, tau, Es

def F_state(i, p, disjoint_js, Theta,
           alpha=1.0, beta=1.0, lamV=0.2, lamD=0.2, q=0.6, tau_scale=0.02):
    Omega, Ecut, tau, _ = soft_count_omega(i,p,disjoint_js,q,tau_scale,alpha,beta,lamV,lamD)
    # avoid log(0)
    Omega = max(Omega, 1e-12)
    return E_geom(i,p,alpha,beta) - Theta*math.log(Omega), Omega, Ecut, tau

# build disjoint lists for each tetrahedron
disjoint = []
for i in range(len(tets)):
    dis = [j for j in range(len(tets)) if i!=j and tset[i].isdisjoint(tset[j])]
    disjoint.append(dis)
# check degrees are 57/59/61
deg = [len(disjoint[i]) for i in range(len(tets))]
print("degree distribution:", {k:deg.count(k) for k in sorted(set(deg))})

# Example run: pick a tetrahedron index i0
i0 = 0
ps = np.linspace(0.0, 1.0, 200)

# Theta = eta*h_H ; here we set h_H=1 for dimensionless demo
eta = 0.01
h_H = 1.0
Theta = eta*h_H

Fs = []
Omegas = []
for p in ps:
    Fv, Om, Ecut, tau = F_state(i0, p, disjoint[i0], Theta,
                                alpha=1.0, beta=1.0, lamV=0.2, lamD=0.2, q=0.6, tau_scale=0.02)
    Fs.append(Fv)
    Omegas.append(Om)

Fs = np.array(Fs)
Omegas = np.array(Omegas)

# numerical derivative as schuett-stress
sigma = np.diff(Fs)/np.diff(ps)

print("Omega range:", float(np.min(Omegas)), float(np.max(Omegas)))
print("sigma range:", float(np.min(sigma)), float(np.max(sigma)))