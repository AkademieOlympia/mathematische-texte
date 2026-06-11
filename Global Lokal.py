import numpy as np
from itertools import product
from collections import Counter

# ============================================================
# 1. Klein-Vierergruppe V4 = {E, A, B, C}
# ============================================================

SIGMAS = ["E", "A", "B", "C"]

# Multiplikationstabelle
mul_table = {
    ("E", "E"): "E", ("E", "A"): "A", ("E", "B"): "B", ("E", "C"): "C",
    ("A", "E"): "A", ("A", "A"): "E", ("A", "B"): "C", ("A", "C"): "B",
    ("B", "E"): "B", ("B", "A"): "C", ("B", "B"): "E", ("B", "C"): "A",
    ("C", "E"): "C", ("C", "A"): "B", ("C", "B"): "A", ("C", "C"): "E",
}

def mul(x: str, y: str) -> str:
    return mul_table[(x, y)]

def total_signature(state):
    """Produkt sigma1*sigma2*sigma3*sigma4."""
    s = "E"
    for x in state:
        s = mul(s, x)
    return s

# ============================================================
# 2. Neutraler Vierlingsraum X_E
# ============================================================

all_states = list(product(SIGMAS, repeat=4))
neutral_states = [s for s in all_states if total_signature(s) == "E"]

assert len(neutral_states) == 64, f"Erwartet 64 Zustände, erhalten: {len(neutral_states)}"

state_to_idx = {s: i for i, s in enumerate(neutral_states)}
idx_to_state = {i: s for s, i in state_to_idx.items()}

n = len(neutral_states)

# ============================================================
# 3. Elementare Operationen
# ============================================================

def shift(state):
    """Zyklischer Shift T(s1,s2,s3,s4) = (s2,s3,s4,s1)."""
    return (state[1], state[2], state[3], state[0])

def shift_inv(state):
    """Inverse Verschiebung T^{-1}(s1,s2,s3,s4) = (s4,s1,s2,s3)."""
    return (state[3], state[0], state[1], state[2])

def pair_flip(state, i, j, g):
    """
    Multipliziere Einträge an Position i und j mit g.
    i, j in {0,1,2,3}, i<j, g in {A,B,C}
    """
    s = list(state)
    s[i] = mul(g, s[i])
    s[j] = mul(g, s[j])
    return tuple(s)

def count_nontrivial(state):
    """Anzahl nichttrivialer Einträge in {A,B,C}."""
    return sum(1 for x in state if x != "E")

# ============================================================
# 4. Operatoren als Matrizen
# ============================================================

def permutation_operator(op):
    """
    Baut Matrix U für eine zustandserhaltende Permutation/Abbildung op
    auf dem neutralen Raum.
    Konvention: U |state> = |op(state)>
    """
    U = np.zeros((n, n), dtype=np.float64)
    for s in neutral_states:
        t = op(s)
        if t not in state_to_idx:
            raise ValueError(f"Bildzustand {t} liegt nicht im neutralen Sektor.")
        U[state_to_idx[t], state_to_idx[s]] = 1.0
    return U

# Shiftoperatoren
T = permutation_operator(shift)
T_inv = permutation_operator(shift_inv)

# Paarflip-Summe
pair_indices = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
nontrivial_group_elements = ["A", "B", "C"]

F_sum = np.zeros((n, n), dtype=np.float64)

for i, j in pair_indices:
    for g in nontrivial_group_elements:
        F_ij_g = permutation_operator(lambda s, i=i, j=j, g=g: pair_flip(s, i, j, g))
        F_sum += F_ij_g

# ============================================================
# 5. Potentialterm V
# ============================================================

def build_potential(mu=1.0, alpha=0.0):
    """
    Einfaches diagonales Potential:
        v(state) = mu * N_nt(state) + alpha * I_ABC_E(state)
    wobei I_ABC_E = 1 genau dann, wenn die vier Symbole {E,A,B,C}
    jeweils genau einmal vorkommen.
    """
    V = np.zeros((n, n), dtype=np.float64)
    for s in neutral_states:
        cnt = Counter(s)
        indicator_abce = int(
            cnt["E"] == 1 and cnt["A"] == 1 and cnt["B"] == 1 and cnt["C"] == 1
        )
        v = mu * count_nontrivial(s) + alpha * indicator_abce
        V[state_to_idx[s], state_to_idx[s]] = v
    return V

# ============================================================
# 6. Hamiltonoperator
# ============================================================

def build_hamiltonian(t=1.0, lam=0.2, mu=0.5, alpha=0.0):
    V = build_potential(mu=mu, alpha=alpha)
    H = -t * (T + T_inv) - lam * F_sum + V
    return H

# ============================================================
# 7. Diagonalisierung und Analyse
# ============================================================

def diagonalize_hamiltonian(t=1.0, lam=0.2, mu=0.5, alpha=0.0, top_k=12):
    H = build_hamiltonian(t=t, lam=lam, mu=mu, alpha=alpha)
    evals, evecs = np.linalg.eigh(H)

    print("=" * 70)
    print("BM-Vierlingsmodell: Spektrum im neutralen Sektor")
    print("=" * 70)
    print(f"Dimension des Raums: {n}")
    print(f"Parameter: t={t}, lambda={lam}, mu={mu}, alpha={alpha}")
    print()

    print(f"Die {top_k} kleinsten Eigenwerte:")
    for k in range(min(top_k, len(evals))):
        print(f"{k:2d}: {evals[k]: .8f}")
    print()

    # einfache Entartungsanalyse (gerundet)
    rounded = np.round(evals, 10)
    degeneracies = Counter(rounded)
    print("Häufigkeiten gerundeter Eigenwerte (erste 20):")
    for val, deg in sorted(degeneracies.items())[:20]:
        print(f"{val: .10f}  ->  {deg}")
    print()

    return H, evals, evecs

# ============================================================
# 8. T-Eigenmoden eines gegebenen Orbits anschauen
# ============================================================

def orbit_under_shift(state):
    """Erzeuge den T-Orbit eines Zustands."""
    orbit = []
    seen = set()
    s = state
    while s not in seen:
        seen.add(s)
        orbit.append(s)
        s = shift(s)
    return orbit

def restricted_shift_block(state):
    """
    Shiftblock auf dem Orbit eines Zustands.
    Nützlich zur Kontrolle des 4x4-Blocks für z.B. (A,B,C,E).
    """
    orbit = orbit_under_shift(state)
    m = len(orbit)
    idx = {s: i for i, s in enumerate(orbit)}
    M = np.zeros((m, m), dtype=np.float64)
    for s in orbit:
        M[idx[shift(s)], idx[s]] = 1.0
    return orbit, M


def print_dominant_components(evals, evecs, k=0, top_n=12):
    """
    Zeige die größten Basisanteile des k-ten Eigenzustands.
    """
    if not (0 <= k < evecs.shape[1]):
        raise IndexError(f"k={k} liegt außerhalb des Eigenzustandsbereichs")

    vec = evecs[:, k]
    probs = np.abs(vec) ** 2
    order = np.argsort(probs)[::-1]
    top_n = min(top_n, len(order))

    print("=" * 70)
    print(f"Dominante Komponenten des Eigenzustands k={k}, Eigenwert={evals[k]:.8f}")
    print("=" * 70)

    total_shown = 0.0
    for r in range(top_n):
        idx = order[r]
        state = idx_to_state[idx]
        weight = probs[idx]
        total_shown += weight
        print(f"{r + 1:2d}: {state}   |amp|^2 = {weight:.8f}")

    print(f"\nSumme der gezeigten Gewichte: {total_shown:.8f}")
    print(f"Restgewicht: {1.0 - total_shown:.8f}")
    print()


def classify_state_type(state):
    """
    Grobe Typklassifikation neutraler Vierlinge.
    """
    cnt = Counter(state)
    pattern = tuple(sorted(cnt.values(), reverse=True))

    if cnt["E"] == 4:
        return "EEEE"
    if cnt["E"] == 1 and cnt["A"] == 1 and cnt["B"] == 1 and cnt["C"] == 1:
        return "ABCE"
    if pattern == (2, 2):
        return "AABB-Typ"
    if pattern == (2, 1, 1):
        return "AABC-Typ"
    return f"sonst {pattern}"


def analyze_eigenstate_types(evals, evecs, k=0, top_n=20):
    if not (0 <= k < evecs.shape[1]):
        raise IndexError(f"k={k} liegt außerhalb des Eigenzustandsbereichs")

    vec = evecs[:, k]
    probs = np.abs(vec) ** 2
    order = np.argsort(probs)[::-1]
    top_n = min(top_n, len(order))

    print("=" * 70)
    print(f"Typanalyse des Eigenzustands k={k}, Eigenwert={evals[k]:.8f}")
    print("=" * 70)

    type_weights = Counter()
    for r in range(top_n):
        idx = order[r]
        state = idx_to_state[idx]
        weight = probs[idx]
        typ = classify_state_type(state)
        type_weights[typ] += weight

    for typ, weight in type_weights.most_common():
        print(f"{typ:12s} -> {weight:.8f}")
    print()


def type_weights_of_eigenstate(evecs, k=0):
    if not (0 <= k < evecs.shape[1]):
        raise IndexError(f"k={k} liegt außerhalb des Eigenzustandsbereichs")

    vec = evecs[:, k]
    probs = np.abs(vec) ** 2
    weights = Counter()

    for idx, prob in enumerate(probs):
        state = idx_to_state[idx]
        typ = classify_state_type(state)
        weights[typ] += prob

    return dict(weights)


def scan_alpha(alpha_values, t=1.0, lam=0.15, mu=0.4):
    print("=" * 90)
    print("Alpha-Scan: Typgewichte des Grundzustands")
    print("=" * 90)
    print(f"{'alpha':>8s} | {'E0':>12s} | {'AABB':>10s} | {'ABCE':>10s} | {'EEEE':>10s}")
    print("-" * 90)

    for alpha in alpha_values:
        H = build_hamiltonian(t=t, lam=lam, mu=mu, alpha=alpha)
        evals, evecs = np.linalg.eigh(H)
        weights = type_weights_of_eigenstate(evecs, k=0)

        aabb = weights.get("AABB-Typ", 0.0)
        abce = weights.get("ABCE", 0.0)
        eeee = weights.get("EEEE", 0.0)

        print(f"{alpha:8.3f} | {evals[0]:12.8f} | {aabb:10.6f} | {abce:10.6f} | {eeee:10.6f}")


def scan_alpha_fine(alpha_values, t=1.0, lam=0.15, mu=0.4):
    print("=" * 100)
    print("Feinscan Alpha: Typgewichte des Grundzustands")
    print("=" * 100)
    print(f"{'alpha':>8s} | {'E0':>12s} | {'AABB':>10s} | {'ABCE':>10s} | {'EEEE':>10s} | {'ABCE-AABB':>12s}")
    print("-" * 100)

    for alpha in alpha_values:
        H = build_hamiltonian(t=t, lam=lam, mu=mu, alpha=alpha)
        evals, evecs = np.linalg.eigh(H)
        weights = type_weights_of_eigenstate(evecs, k=0)

        aabb = weights.get("AABB-Typ", 0.0)
        abce = weights.get("ABCE", 0.0)
        eeee = weights.get("EEEE", 0.0)

        print(
            f"{alpha:8.3f} | {evals[0]:12.8f} | {aabb:10.6f} | "
            f"{abce:10.6f} | {eeee:10.6f} | {abce - aabb:12.6f}"
        )

# ============================================================
# 9. Beispielausführung
# ============================================================

if __name__ == "__main__":
    # Gesamtspektrum
    H, evals, evecs = diagonalize_hamiltonian(
        t=1.0,
        lam=0.15,
        mu=0.4,
        alpha=-0.3,   # bevorzugt ABCE-Typen leicht energetisch
        top_k=16
    )

    # 4er-Shiftorbit von (A,B,C,E)
    seed = ("A", "B", "C", "E")
    orbit, T_block = restricted_shift_block(seed)

    print("=" * 70)
    print("Shift-Orbit von (A,B,C,E)")
    print("=" * 70)
    for i, s in enumerate(orbit):
        print(f"{i}: {s}")

    print("\nShift-Matrix auf dem Orbit:")
    print(T_block)

    # Transport-Hamiltonian auf dem Orbit: H_trans = -t(T + T^{-1})
    t_val = 1.0
    H_trans_block = -t_val * (T_block + T_block.T)
    block_evals, _ = np.linalg.eigh(H_trans_block)

    print("\nH_trans-Block:")
    print(H_trans_block)
    print("\nEigenwerte des 4x4-Transportblocks:")
    print(block_evals)

    print_dominant_components(evals, evecs, k=0, top_n=16)
    analyze_eigenstate_types(evals, evecs, k=0, top_n=30)

    print_dominant_components(evals, evecs, k=1, top_n=16)
    analyze_eigenstate_types(evals, evecs, k=1, top_n=30)

    scan_alpha(np.linspace(-1.0, 0.5, 7), t=1.0, lam=0.15, mu=0.4)
    scan_alpha_fine(np.linspace(-0.8, 0.1, 10), t=1.0, lam=0.15, mu=0.4)