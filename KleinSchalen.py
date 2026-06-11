from dataclasses import dataclass
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import dok_matrix, csr_matrix
from scipy.sparse.linalg import eigsh


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


def diag_energy(st: State, p: Dict[str, float]) -> float:
    return (
        p["lambda_2"] * st.nu2
        + p["lambda_3"] * st.nu3
        + p["lambda_5"] * st.sigma5
        + p["lambda_E"] * st.wE
        + p["lambda_A"] * st.wA
        + p["lambda_B"] * st.wB
        + p["lambda_C"] * st.wC
        + p["mu"] * balance(st)
        + p["rho"] * st.eps5
    )


def build_basis(shell_max: int, inner_max: int) -> List[State]:
    basis: List[State] = []
    for nu2 in range(shell_max + 1):
        for nu3 in range(shell_max + 1 - nu2):
            for sigma5 in range(shell_max + 1 - nu2 - nu3):
                for eps5 in (0, 1):
                    for wE in range(inner_max + 1):
                        for wA in range(inner_max + 1 - wE):
                            for wB in range(inner_max + 1 - wE - wA):
                                for wC in range(inner_max + 1 - wE - wA - wB):
                                    # Wenn eps5 = 1, soll im A-Kanal mindestens eine Einheit sitzen.
                                    if eps5 == 1 and wA == 0:
                                        continue
                                    basis.append(
                                        State(nu2, nu3, sigma5, eps5, wE, wA, wB, wC)
                                    )
    return basis


def add_edge(
    mat: dok_matrix,
    idx: Dict[State, int],
    st1: State,
    st2: State,
    weight: float,
) -> None:
    i = idx.get(st1)
    j = idx.get(st2)
    if i is None or j is None or i == j:
        return
    mat[i, j] = mat.get((i, j), 0.0) + weight
    mat[j, i] = mat.get((j, i), 0.0) + weight


def build_operator(
    basis: List[State],
    p: Dict[str, float],
    shell_max: int,
) -> csr_matrix:
    idx = {st: i for i, st in enumerate(basis)}
    n = len(basis)
    mat = dok_matrix((n, n), dtype=np.float64)

    # Diagonale
    for i, st in enumerate(basis):
        mat[i, i] = diag_energy(st, p)

    # Off-diagonal
    for st in basis:
        # äußere Schalen
        if st.nu2 + 1 + st.nu3 + st.sigma5 <= shell_max:
            add_edge(
                mat, idx, st,
                State(st.nu2 + 1, st.nu3, st.sigma5, st.eps5, st.wE, st.wA, st.wB, st.wC),
                p["g2"],
            )

        if st.nu3 + 1 + st.nu2 + st.sigma5 <= shell_max:
            add_edge(
                mat, idx, st,
                State(st.nu2, st.nu3 + 1, st.sigma5, st.eps5, st.wE, st.wA, st.wB, st.wC),
                p["g3"],
            )

        if st.sigma5 + 1 + st.nu2 + st.nu3 <= shell_max:
            add_edge(
                mat, idx, st,
                State(st.nu2, st.nu3, st.sigma5 + 1, st.eps5, st.wE, st.wA, st.wB, st.wC),
                p["g5"],
            )

        # innere Ringkopplung E -> A und E -> C
        if st.wE >= 1:
            add_edge(
                mat, idx, st,
                State(st.nu2, st.nu3, st.sigma5, st.eps5, st.wE - 1, st.wA + 1, st.wB, st.wC),
                p["g_int"],
            )
            add_edge(
                mat, idx, st,
                State(st.nu2, st.nu3, st.sigma5, st.eps5, st.wE - 1, st.wA, st.wB, st.wC + 1),
                p["g_int"],
            )

        # A -> B
        if st.wA >= 1:
            add_edge(
                mat, idx, st,
                State(st.nu2, st.nu3, st.sigma5, st.eps5, st.wE, st.wA - 1, st.wB + 1, st.wC),
                p["g_int"],
            )

        # B -> C
        if st.wB >= 1:
            add_edge(
                mat, idx, st,
                State(st.nu2, st.nu3, st.sigma5, st.eps5, st.wE, st.wA, st.wB - 1, st.wC + 1),
                p["g_int"],
            )

        # 25-Kollaps: A + A -> 25-Schale
        if st.wA >= 2 and (st.nu2 + st.nu3 + st.sigma5 + 1 <= shell_max):
            add_edge(
                mat, idx, st,
                State(st.nu2, st.nu3, st.sigma5 + 1, st.eps5, st.wE, st.wA - 2, st.wB, st.wC),
                p["g_coll"],
            )

    return mat.tocsr()


def lowest_eigensystem(
    op: csr_matrix,
    k: int = 20,
    sigma: float | None = None,
) -> Tuple[np.ndarray, np.ndarray]:
    n = op.shape[0]
    k = min(k, n - 2)
    # which='SA' = kleinste algebraische Eigenwerte
    vals, vecs = eigsh(op, k=k, which="SA", sigma=sigma)
    order = np.argsort(vals)
    return vals[order], vecs[:, order]


def dominant_components(
    basis: List[State],
    eigvecs: np.ndarray,
    n_levels: int = 10,
    n_components: int = 8,
) -> pd.DataFrame:
    rows = []
    for lev in range(min(n_levels, eigvecs.shape[1])):
        probs = eigvecs[:, lev] ** 2
        top_idx = np.argsort(probs)[-n_components:][::-1]
        for rank, idx0 in enumerate(top_idx, start=1):
            st = basis[idx0]
            rows.append(
                {
                    "level": lev + 1,
                    "rank_in_level": rank,
                    "state": st.label(),
                    "weight_sq": float(probs[idx0]),
                    "nu2": st.nu2,
                    "nu3": st.nu3,
                    "sigma5": st.sigma5,
                    "eps5": st.eps5,
                    "wE": st.wE,
                    "wA": st.wA,
                    "wB": st.wB,
                    "wC": st.wC,
                    "inner_weight": st.inner_weight,
                    "balance": balance(st),
                    "is_vollschale_1111": (
                        st.wE == 1 and st.wA == 1 and st.wB == 1 and st.wC == 1
                    ),
                }
            )
    return pd.DataFrame(rows)


def basis_table(basis: List[State], p: Dict[str, float]) -> pd.DataFrame:
    rows = []
    for st in basis:
        rows.append(
            {
                "state": st.label(),
                "nu2": st.nu2,
                "nu3": st.nu3,
                "sigma5": st.sigma5,
                "eps5": st.eps5,
                "wE": st.wE,
                "wA": st.wA,
                "wB": st.wB,
                "wC": st.wC,
                "inner_weight": st.inner_weight,
                "balance": balance(st),
                "energy_diag": diag_energy(st, p),
                "is_vollschale_1111": (
                    st.wE == 1 and st.wA == 1 and st.wB == 1 and st.wC == 1
                ),
                "is_closed_equal": (
                    st.wE == st.wA == st.wB == st.wC
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(["energy_diag", "state"]).reset_index(drop=True)


def main() -> None:
    # Sehr kompakt und lokal stabil
    shell_max = 1
    inner_max = 4

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

    basis = build_basis(shell_max=shell_max, inner_max=inner_max)
    op = build_operator(basis, params, shell_max=shell_max)

    print(f"Dimension des Zustandsraums: {len(basis)}")
    print(f"Nichtnull-Einträge: {op.nnz}")

    eigvals, eigvecs = lowest_eigensystem(op, k=20)

    eig_df = pd.DataFrame(
        {"level": np.arange(1, len(eigvals) + 1), "eigenvalue": eigvals}
    )
    dom_df = dominant_components(basis, eigvecs, n_levels=10, n_components=8)
    bas_df = basis_table(basis, params)

    eig_df.to_csv("bm_dirac_sparse_eigenwerte.csv", index=False)
    dom_df.to_csv("bm_dirac_sparse_dominante_komponenten.csv", index=False)
    bas_df.to_csv("bm_dirac_sparse_basis.csv", index=False)

    print("\nErste 20 Eigenwerte:")
    print(eig_df.to_string(index=False))

    print("\nDominante Komponenten der ersten 6 Levels:")
    print(dom_df[dom_df["level"] <= 6].to_string(index=False))

    # Plot Eigenwerte
    plt.figure(figsize=(10, 4))
    plt.plot(eig_df["level"], eig_df["eigenvalue"], marker="o")
    plt.xlabel("Levelindex")
    plt.ylabel("Eigenwert")
    plt.title("Unterste Eigenwerte des sparse BM-Dirac-Operators")
    plt.tight_layout()
    plt.savefig("bm_dirac_sparse_eigenwerte.png", dpi=180)
    plt.close()

    # Plot Gaps
    gaps = np.diff(eigvals)
    plt.figure(figsize=(10, 4))
    plt.plot(np.arange(1, len(gaps) + 1), gaps, marker="o")
    plt.xlabel("Gapindex")
    plt.ylabel("Levelabstand")
    plt.title("Levelabstände des sparse BM-Dirac-Operators")
    plt.tight_layout()
    plt.savefig("bm_dirac_sparse_gaps.png", dpi=180)
    plt.close()

    # Dominante Komponenten der ersten 4 Levels
    _, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes = axes.ravel()
    for lev in range(4):
        sub = dom_df[dom_df["level"] == lev + 1].sort_values("rank_in_level")
        axes[lev].bar(range(len(sub)), sub["weight_sq"])
        axes[lev].set_xticks(range(len(sub)))
        axes[lev].set_xticklabels(sub["state"], rotation=70, fontsize=7)
        axes[lev].set_ylabel("Gewicht²")
        axes[lev].set_title(f"Level {lev+1}")
    plt.tight_layout()
    plt.savefig("bm_dirac_sparse_dominante_levels.png", dpi=180)
    plt.close()

    # Spezielle Tabellen
    vacuum_shell = bas_df[bas_df["inner_weight"] == 0].copy()
    sigma5_states = bas_df[bas_df["sigma5"] >= 1].copy()
    aa_states = bas_df[bas_df["wA"] >= 2].copy()
    voll = bas_df[bas_df["is_vollschale_1111"]].copy()

    vacuum_shell.to_csv("bm_dirac_sparse_vakuum_schale.csv", index=False)
    sigma5_states.to_csv("bm_dirac_sparse_sigma5_zustaende.csv", index=False)
    aa_states.to_csv("bm_dirac_sparse_AA_zustaende.csv", index=False)
    voll.to_csv("bm_dirac_sparse_vollschalen.csv", index=False)

    print("\nNiedrigste reine Schalenzustände:")
    print(vacuum_shell.head(10).to_string(index=False))

    print("\nNiedrigste sigma5-Zustände:")
    print(sigma5_states.head(10).to_string(index=False))

    print("\nNiedrigste A+A-Zustände:")
    print(aa_states.head(10).to_string(index=False))

    if len(voll) > 0:
        print("\nVollschalen (1,1,1,1):")
        print(voll.to_string(index=False))
    else:
        print("\nIn dieser Minimalbasis gibt es keine Vollschale (1,1,1,1),")
        print("weil INNER_MAX = 2 ist. Für Vollschalen brauchst du mindestens INNER_MAX = 4.")


if __name__ == "__main__":
    main()