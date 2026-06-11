import numpy as np
import kwant
from itertools import product

# --- EABC Klassifikation ---
def eabc_class(p):
    r = p % 12
    if r == 1: return "E"
    if r == 5: return "A"
    if r == 7: return "B"
    if r == 11: return "C"

def sigma(P):
    s = set()
    for p in P:
        c = eabc_class(p)
        if c is not None and c != "E":
            s.add(c)
    return sorted(s)


def participation_ratio(psi):
    """Teilnahmezahl PR = 1 / sum_i p_i^2 mit p_i = |ψ_i|^2 / sum_j |ψ_j|^2.

    Bei bereits normiertem ψ (sum |ψ_i|^2 = 1) entspricht das der Kurzform
    return 1.0 / np.sum(np.abs(psi) ** 4).
    """
    p = np.abs(psi) ** 2
    s = np.sum(p)
    if not np.isfinite(s) or s <= 0:
        return np.nan
    p = p / s
    return 1.0 / np.sum(p ** 2)


# --- Systembau ---
def make_system(L, M, P, t=1.0):
    lat = kwant.lattice.square(norbs=1)
    sys = kwant.Builder()

    for x, y in product(range(L), range(M)):
        eps = np.random.choice(P)
        sys[lat(x, y)] = eps

        if x > 0:
            sys[lat(x, y), lat(x-1, y)] = -t
        if y > 0:
            sys[lat(x, y), lat(x, y-1)] = -t

    # Leads (einfach gehalten)
    lead = kwant.Builder(kwant.TranslationalSymmetry((-1, 0)))
    for y in range(M):
        lead[lat(0, y)] = 0
    lead[lat.neighbors()] = -t

    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys.finalized()

# --- Transmission / Leitwert (Landauer: G = G_0 * T in Einheiten von 2 e^2/h)
def conductance(sys, energy):
    smatrix = kwant.smatrix(sys, energy)
    T = smatrix.transmission(1, 0)
    return T  # proportional zu G


transmission = conductance


def log_transport(sys, energies):
    """Mittelwert von ln(T + eps) über Energien (eps vermeidet ln 0)."""
    Ts = np.array([conductance(sys, E) for E in energies])
    return np.mean(np.log(Ts + 1e-12))


# --- Scan ---
def scan_system(P, energies=np.linspace(-0.1,0.1,200)):
    sys = make_system(6,7,P)
    Tvals = np.array([conductance(sys, E) for E in energies])
    Tnorm = Tvals / np.max(Tvals)
    return np.mean(Tnorm), np.std(Tnorm)


def disorder_scan(W_values):
    """Störstärke |W|: Pool P ~ U(-|W|, |W|), dann make_system wählt pro Site aus P.

    Rückgabe: Liste von (W_signiert, mean(T), std(T)) über Energien -0.2..0.2 (200 Punkte).
    """
    energies = np.linspace(-0.2, 0.2, 200)
    out = []
    for W in W_values:
        W = float(W)
        if not np.isfinite(W):
            continue
        a = abs(W)
        P = np.random.uniform(-a, a, size=100)
        sys = make_system(6, 7, P)
        Ts = [conductance(sys, E) for E in energies]
        out.append((W, float(np.mean(Ts)), float(np.std(Ts))))
    return out


def size_scaling(sizes, W):
    """Mittlere Transmission vs. Kantenlänge L (System L×L), festes Störintervall |W|.

    Pool P ~ U(-|W|, |W|), make_system(L, L, P); Energien -0.2 … 0.2 (100 Punkte).
    Rückgabe: Liste (L, mean(T)).
    """
    energies = np.linspace(-0.2, 0.2, 100)
    a = abs(float(W))
    if not np.isfinite(a):
        return []
    out = []
    for L in sizes:
        L = int(L)
        if L < 1:
            continue
        P = np.random.uniform(-a, a, size=100)
        sys = make_system(L, L, P)
        Ts = [conductance(sys, E) for E in energies]
        out.append((L, float(np.mean(Ts))))
    return out


# --- Testfälle ---
test_sets = {
    "AC": [11,17,29,41,47,59],
    "ABC": [5,7,11,13,17,19],
    "AB": [5,7,17,19,29],
    "A": [5,17,29,41]
}

results = {}

for name, P in test_sets.items():
    meanT, stdT = scan_system(P)
    results[name] = {
        "Sigma": sigma(P),
        "T_mean": meanT,
        "T_std": stdT
    }

# Gleichverteilte Kanalpotentiale auf [10, 60] (wie Streuwerte für random.choice(P))
np.random.seed(12345)
P_random = np.random.uniform(10, 60, size=6)
meanT, stdT = scan_system(P_random)
results["uniform_10_60"] = {
    "Sigma": sigma(P_random),
    "P": P_random.tolist(),
    "T_mean": meanT,
    "T_std": stdT,
}

np.random.seed(12346)
P = np.random.uniform(-200, 200, size=200)
meanT, stdT = scan_system(P)
results["uniform_-200_200_x200"] = {
    "Sigma": sigma(P),
    "n_P": len(P),
    "P_preview": P[:10].tolist(),
    "T_mean": meanT,
    "T_std": stdT,
}

sizes = [6, 10, 20, 30]
np.random.seed(12347)
W_scaling = 50.0
results["size_scaling"] = {
    "W": W_scaling,
    "sizes": sizes,
    "L_vs_meanT": size_scaling(sizes, W_scaling),
}

for k, v in results.items():
    print(k, v)