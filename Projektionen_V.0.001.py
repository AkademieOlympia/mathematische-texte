import numpy as np
import pandas as pd

# Basis und Eigenvektoren aus KleinSchalen erzeugen
from KleinSchalen import build_basis, build_operator, lowest_eigensystem, dominant_components

shell_max = 1
inner_max = 4
params = {
    "lambda_2": 1.0, "lambda_3": 1.2, "lambda_5": 5.8,
    "lambda_E": 3.7, "lambda_A": 4.0, "lambda_B": 4.3, "lambda_C": 4.6,
    "mu": 1.5, "rho": 1.1,
    "g2": 0.35, "g3": 0.40, "g5": 0.75, "g_int": 0.28, "g_coll": 0.95,
}
basis = build_basis(shell_max=shell_max, inner_max=inner_max)
op = build_operator(basis, params, shell_max=shell_max)
vals, vecs = lowest_eigensystem(op, k=120)

# Erwartet:
# - basis: Liste von State-Objekten
# - vecs: Eigenvektormatrix, Spalten = Eigenzustände

def sector_weights(basis, vecs, n_levels=20):
    rows = []

    shell_idx = []
    one_idx = []
    mix_idx = []
    full_idx = []
    sigma5_idx = []
    aa_idx = []

    for i, st in enumerate(basis):
        inner_weight = st.wE + st.wA + st.wB + st.wC

        if inner_weight == 0:
            shell_idx.append(i)
        elif inner_weight == 1:
            one_idx.append(i)
        elif (st.wE, st.wA, st.wB, st.wC) == (1, 1, 1, 1):
            full_idx.append(i)
        else:
            mix_idx.append(i)

        if st.sigma5 >= 1:
            sigma5_idx.append(i)
        if st.wA >= 2:
            aa_idx.append(i)

    for lev in range(min(n_levels, vecs.shape[1])):
        probs = vecs[:, lev] ** 2

        rows.append({
            "level": lev + 1,
            "W_shell": float(np.sum(probs[shell_idx])),
            "W_one": float(np.sum(probs[one_idx])),
            "W_mix": float(np.sum(probs[mix_idx])),
            "W_full": float(np.sum(probs[full_idx])),
            "W_25": float(np.sum(probs[sigma5_idx])),
            "W_AA": float(np.sum(probs[aa_idx])),
        })

    df = pd.DataFrame(rows)
    df["dominant_sector"] = df[["W_shell", "W_one", "W_mix", "W_full"]].idxmax(axis=1)
    return df
proj_df = sector_weights(basis, vecs, n_levels=120)
proj_df["eigenvalue"] = vals[:len(proj_df)]
proj_df.to_csv("bm_dirac_sparse_inner4_sector_weights.csv", index=False)

print("\nSektorprojektionen:")
print(proj_df)

full_010 = proj_df[proj_df["W_full"] > 0.10][["level", "eigenvalue", "W_full", "dominant_sector"]]
full_025 = proj_df[proj_df["W_full"] > 0.25][["level", "eigenvalue", "W_full", "dominant_sector"]]
full_050 = proj_df[proj_df["W_full"] > 0.50][["level", "eigenvalue", "W_full", "dominant_sector"]]

print("\nErstes Level mit W_full > 0.10:")
print(full_010.head(1).to_string(index=False))

print("\nErstes Level mit W_full > 0.25:")
print(full_025.head(1).to_string(index=False))

print("\nErstes Level mit W_full > 0.50:")
print(full_050.head(1).to_string(index=False))

dom_df = dominant_components(basis, vecs, n_levels=120, n_components=8)
for lev in range(min(120, vecs.shape[1])):
    print(f"\nDominante Komponenten von Level {lev + 1}:")
    sub = dom_df[dom_df["level"] == lev + 1].sort_values("rank_in_level")
    print(sub.to_string(index=False))

import matplotlib.pyplot as plt

plt.figure(figsize=(12, 6))
plt.plot(proj_df["level"], proj_df["W_shell"], marker="o", label="shell")
plt.plot(proj_df["level"], proj_df["W_one"], marker="o", label="one-channel")
plt.plot(proj_df["level"], proj_df["W_mix"], marker="o", label="mix")
plt.plot(proj_df["level"], proj_df["W_full"], marker="o", label="full-shell")
plt.xlabel("Level")
plt.ylabel("Projektionsgewicht")
plt.title("Sektorprojektionen des sparse BM-Dirac (INNER_MAX=4)")
plt.legend()
plt.tight_layout()
plt.savefig("bm_dirac_sparse_inner4_sector_weights.png", dpi=180)
plt.close()

plt.figure(figsize=(12, 5))
plt.plot(proj_df["level"], proj_df["W_25"], marker="o", label="25-shell")
plt.plot(proj_df["level"], proj_df["W_AA"], marker="o", label="A+A sector")
plt.xlabel("Level")
plt.ylabel("Projektionsgewicht")
plt.title("25-Schale vs. A+A-Sektor")
plt.legend()
plt.tight_layout()
plt.savefig("bm_dirac_sparse_inner4_25_vs_AA.png", dpi=180)
plt.close()