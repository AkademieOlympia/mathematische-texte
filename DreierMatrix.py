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

# Zielzustände im balancierten Vierer-Subraum
target_states = [
    State(0, 0, 0, 0, 1, 1, 1, 1),  # F
    State(0, 0, 0, 0, 0, 2, 2, 0),  # A
    State(0, 0, 0, 0, 2, 0, 0, 2),  # E
]

labels = ["F=(1,1,1,1)", "A=(0,2,2,0)", "E=(2,0,0,2)"]

sub_idx = [idx[st] for st in target_states]

# Projektionsmatrix
D_eff = op[sub_idx, :][:, sub_idx].toarray()

eff_df = pd.DataFrame(D_eff, index=labels, columns=labels)
print("\nEffektive 3x3-Matrix im balancierten Vierer-Subraum:")
print(eff_df)

# Eigenwerte/Eigenvektoren des reduzierten Operators
vals_eff, vecs_eff = np.linalg.eigh(D_eff)
order = np.argsort(vals_eff)
vals_eff = vals_eff[order]
vecs_eff = vecs_eff[:, order]

print("\nEigenwerte von D_eff:")
for i, val in enumerate(vals_eff, start=1):
    print(f"{i}: {val:.6f}")

print("\nEigenvektoren von D_eff in der Basis [F, A, E]:")
for i in range(len(vals_eff)):
    print(f"Level {i+1}: {vecs_eff[:, i]}")