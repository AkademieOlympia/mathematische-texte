import numpy as np
import pandas as pd

from KleinSchalen import State, build_basis, build_operator

shell_max = 1
inner_max = 4
params = {
    "lambda_2": 1.0, "lambda_3": 1.2, "lambda_5": 5.8,
    "lambda_E": 3.7, "lambda_A": 4.0, "lambda_B": 4.3, "lambda_C": 4.6,
    "mu": 1.5, "rho": 1.1,
    "g2": 0.35, "g3": 0.40, "g5": 0.75, "g_int": 0.28, "g_coll": 0.95,
}
basis = build_basis(shell_max=shell_max, inner_max=inner_max)
idx = {st: i for i, st in enumerate(basis)}
op = build_operator(basis, params, shell_max=shell_max)

target_states = [
    State(0,0,0,0,1,1,1,1),  # F
    State(0,0,0,0,0,2,2,0),  # A
    State(0,0,0,0,2,0,0,2),  # E
]
labels = ["F", "A", "E"]

sub_idx = [idx[s] for s in target_states]
all_idx = set(range(op.shape[0]))
rest_idx = sorted(all_idx - set(sub_idx))

E0 = 16.6

# volle Matrix als csr oder dense-Zugriffe
# op muss scipy sparse csr_matrix sein
Heff2 = np.zeros((3, 3), dtype=float)

for a, ia in enumerate(sub_idx):
    for b, ib in enumerate(sub_idx):
        val = 0.0
        for mu in rest_idx:
            Emu = op[mu, mu]
            v1 = op[ia, mu]
            v2 = op[mu, ib]
            if abs(v1) > 1e-14 and abs(v2) > 1e-14 and abs(E0 - Emu) > 1e-12:
                val += (v1 * v2) / (E0 - Emu)
        Heff2[a, b] = val

df_eff2 = pd.DataFrame(Heff2, index=labels, columns=labels)
print("\nZweite-Ordnungs-Korrekturmatrix:")
print(df_eff2)

Heff_total = 16.6 * np.eye(3) + Heff2
df_total = pd.DataFrame(Heff_total, index=labels, columns=labels)
print("\nEffektive Gesamtmatrix:")
print(df_total)

vals2, vecs2 = np.linalg.eigh(Heff_total)
order = np.argsort(vals2)
vals2 = vals2[order]
vecs2 = vecs2[:, order]

print("\nEffektive Eigenwerte:")
for i, val in enumerate(vals2, start=1):
    print(i, val)

print("\nEffektive Eigenvektoren in Basis [F,A,E]:")
for i in range(3):
    print(f"Mode {i+1}: {vecs2[:, i]}")