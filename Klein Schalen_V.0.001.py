import numpy as np
import pandas as pd
from dataclasses import dataclass
from scipy.sparse import dok_matrix
from scipy.sparse.linalg import eigsh

# -----------------------------------------
# Parameter
# -----------------------------------------

SHELL_MAX = 1
INNER_MAX = 4
N_EIGS = 120

params = {
    "lambda_2": 1.0,
    "lambda_3": 1.2,
    "lambda_5": 5.8,
    "lambda_E": 3.7,
    "lambda_A": 4.0,
    "lambda_B": 4.3,
    "lambda_C": 4.6,
    "mu": 1.5,
    "rho": 1.1,
    "g2": 0.35,
    "g3": 0.40,
    "g5": 0.75,
    "g_int": 0.28,
    "g_coll": 0.95,
}

# -----------------------------------------
# Zustand
# -----------------------------------------

@dataclass(frozen=True)
class State:
    nu2: int
    nu3: int
    sigma5: int
    eps5: int
    wE: int
    wA: int
    wB: int
    wC: int

    @property
    def inner_weight(self) -> int:
        return self.wE + self.wA + self.wB + self.wC

    def label(self) -> str:
        return (
            f"|{self.nu2},{self.nu3},{self.sigma5};"
            f"{self.eps5};{self.wE},{self.wA},{self.wB},{self.wC}>"
        )

def balance(st: State) -> int:
    return (st.wE - st.wC) ** 2 + (st.wA - st.wB) ** 2

def diag_energy(st: State) -> float:
    return (
        params["lambda_2"] * st.nu2
        + params["lambda_3"] * st.nu3
        + params["lambda_5"] * st.sigma5
        + params["lambda_E"] * st.wE
        + params["lambda_A"] * st.wA
        + params["lambda_B"] * st.wB
        + params["lambda_C"] * st.wC
        + params["mu"] * balance(st)
        + params["rho"] * st.eps5
    )

# -----------------------------------------
# Basis
# -----------------------------------------

def build_basis(shell_max: int, inner_max: int):
    basis = []
    for nu2 in range(shell_max + 1):
        for nu3 in range(shell_max + 1 - nu2):
            for sigma5 in range(shell_max + 1 - nu2 - nu3):
                for eps5 in (0, 1):
                    for wE in range(inner_max + 1):
                        for wA in range(inner_max + 1 - wE):
                            for wB in range(inner_max + 1 - wE - wA):
                                for wC in range(inner_max + 1 - wE - wA - wB):
                                    if eps5 == 1 and wA == 0:
                                        continue
                                    basis.append(State(nu2, nu3, sigma5, eps5, wE, wA, wB, wC))
    return basis

basis = build_basis(SHELL_MAX, INNER_MAX)
idx = {st: i for i, st in enumerate(basis)}
N = len(basis)

print("Dimension:", N)

# -----------------------------------------
# Sparse-Operator
# -----------------------------------------

mat = dok_matrix((N, N), dtype=np.float64)

for i, st in enumerate(basis):
    mat[i, i] = diag_energy(st)

def add_edge(a: State, b: State, w: float):
    i = idx.get(a)
    j = idx.get(b)
    if i is None or j is None or i == j:
        return
    mat[i, j] = mat.get((i, j), 0.0) + w
    mat[j, i] = mat.get((j, i), 0.0) + w

for st in basis:
    # äußere Schalen
    if st.nu2 + 1 + st.nu3 + st.sigma5 <= SHELL_MAX:
        add_edge(
            st,
            State(st.nu2 + 1, st.nu3, st.sigma5, st.eps5, st.wE, st.wA, st.wB, st.wC),
            params["g2"],
        )

    if st.nu3 + 1 + st.nu2 + st.sigma5 <= SHELL_MAX:
        add_edge(
            st,
            State(st.nu2, st.nu3 + 1, st.sigma5, st.eps5, st.wE, st.wA, st.wB, st.wC),
            params["g3"],
        )

    if st.sigma5 + 1 + st.nu2 + st.nu3 <= SHELL_MAX:
        add_edge(
            st,
            State(st.nu2, st.nu3, st.sigma5 + 1, st.eps5, st.wE, st.wA, st.wB, st.wC),
            params["g5"],
        )

    # innere Ringkopplung
    if st.wE >= 1:
        add_edge(
            st,
            State(st.nu2, st.nu3, st.sigma5, st.eps5, st.wE - 1, st.wA + 1, st.wB, st.wC),
            params["g_int"],
        )
        add_edge(
            st,
            State(st.nu2, st.nu3, st.sigma5, st.eps5, st.wE - 1, st.wA, st.wB, st.wC + 1),
            params["g_int"],
        )

    if st.wA >= 1:
        add_edge(
            st,
            State(st.nu2, st.nu3, st.sigma5, st.eps5, st.wE, st.wA - 1, st.wB + 1, st.wC),
            params["g_int"],
        )

    if st.wB >= 1:
        add_edge(
            st,
            State(st.nu2, st.nu3, st.sigma5, st.eps5, st.wE, st.wA, st.wB - 1, st.wC + 1),
            params["g_int"],
        )

    # 25-Kollaps A+A -> 25
    if st.wA >= 2 and (st.nu2 + st.nu3 + st.sigma5 + 1 <= SHELL_MAX):
        add_edge(
            st,
            State(st.nu2, st.nu3, st.sigma5 + 1, st.eps5, st.wE, st.wA - 2, st.wB, st.wC),
            params["g_coll"],
        )

op = mat.tocsr()

print("Nichtnull-Einträge:", op.nnz)

# -----------------------------------------
# Unterste Eigenwerte
# -----------------------------------------

vals, vecs = eigsh(op, k=min(N_EIGS, N - 2), which="SA")
order = np.argsort(vals)
vals = vals[order]
vecs = vecs[:, order]

eig_df = pd.DataFrame({
    "level": np.arange(1, len(vals) + 1),
    "eigenvalue": vals,
})
eig_df.to_csv("bm_dirac_sparse_inner4_eigenwerte.csv", index=False)

print("\nErste Eigenwerte:")
print(eig_df)

# -----------------------------------------
# Dominante Komponenten
# -----------------------------------------

rows = []
for lev in range(min(8, vecs.shape[1])):
    probs = vecs[:, lev] ** 2
    top_idx = np.argsort(probs)[-10:][::-1]
    for rank, j in enumerate(top_idx, start=1):
        st = basis[j]
        rows.append({
            "level": lev + 1,
            "rank_in_level": rank,
            "state": st.label(),
            "weight_sq": float(probs[j]),
            "diag_energy": diag_energy(st),
            "inner_weight": st.inner_weight,
            "balance": balance(st),
            "is_vollschale_1111": (
                st.wE == 1 and st.wA == 1 and st.wB == 1 and st.wC == 1
            ),
        })

dom_df = pd.DataFrame(rows)
dom_df.to_csv("bm_dirac_sparse_inner4_dominante_komponenten.csv", index=False)

print("\nDominante Komponenten der ersten Levels:")
print(dom_df)

# -----------------------------------------
# Vollschalen-Tabelle
# -----------------------------------------

fullshell = []
for st in basis:
    if st.wE == 1 and st.wA == 1 and st.wB == 1 and st.wC == 1:
        fullshell.append({
            "state": st.label(),
            "diag_energy": diag_energy(st),
            "nu2": st.nu2,
            "nu3": st.nu3,
            "sigma5": st.sigma5,
            "eps5": st.eps5,
        })

full_df = pd.DataFrame(fullshell).sort_values(["diag_energy", "state"])
full_df.to_csv("bm_dirac_sparse_inner4_vollschalen.csv", index=False)

print("\nVollschalen:")
print(full_df)

# -----------------------------------------
# Einfache Plots
# -----------------------------------------

import matplotlib.pyplot as plt

plt.figure(figsize=(10, 4))
plt.plot(eig_df["level"], eig_df["eigenvalue"], marker="o")
plt.xlabel("Levelindex")
plt.ylabel("Eigenwert")
plt.title("Sparse BM-Dirac mit INNER_MAX=4")
plt.tight_layout()
plt.savefig("bm_dirac_sparse_inner4_eigenwerte.png", dpi=180)
plt.close()

gaps = np.diff(vals)
plt.figure(figsize=(10, 4))
plt.plot(np.arange(1, len(gaps) + 1), gaps, marker="o")
plt.xlabel("Gapindex")
plt.ylabel("Levelabstand")
plt.title("Levelabstände, sparse BM-Dirac mit INNER_MAX=4")
plt.tight_layout()
plt.savefig("bm_dirac_sparse_inner4_gaps.png", dpi=180)
plt.close()