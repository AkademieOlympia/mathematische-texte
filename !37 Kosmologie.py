import numpy as np
import matplotlib.pyplot as plt

def kappa_scan(limit=100000,
               window=5000,
               start_synth=137,
               L_list=(120,137,240,274),
               kappas=None,
               n_perm=120,
               block=200,
               only_ABC=True,
               mode="bias_coupled"):
    """
    Scan über kappa und sammle pvals/stats für mehrere L.
    """
    if kappas is None:
        # dicht um 0..2, plus ein paar größere Werte
        kappas = np.concatenate([
            np.linspace(0.0, 0.6, 13),
            np.linspace(0.7, 2.0, 14),
            np.array([2.5, 3.0, 4.0])
        ])

    rows = []
    for kappa in kappas:
        _, _, _, tests, meta = quantum_gravity_primes(
            limit=limit,
            window=window,
            start_synth=start_synth,
            L_list=L_list,
            n_perm=n_perm,
            block=block,
            mode=mode,
            kappa=float(kappa),
            only_ABC=only_ABC
        )
        row = {"kappa": float(kappa)}
        for L in L_list:
            row[f"p{L}"] = tests["pvals"][L]
            row[f"s{L}"] = tests["stats"][L]
        rows.append(row)

        print("kappa={:.3f} | ".format(kappa) + " ".join([f"p{L}={row[f'p{L}']:.4f}" for L in L_list]))

    return rows

def pick_best(rows, L=137):
    # kleinster p-Wert
    key = f"p{L}"
    best = min(rows, key=lambda r: r[key])
    return best

def plot_scan(rows, L_list=(120,137,240,274)):
    kappas = np.array([r["kappa"] for r in rows])

    # Plot p-values
    plt.figure(figsize=(12,6))
    for L in L_list:
        p = np.array([r[f"p{L}"] for r in rows])
        plt.plot(kappas, p, marker="o", linewidth=1.2, label=f"p{L}")
    plt.ylim(-0.02, 1.02)
    plt.title("Resonanz-Signifikanz vs. κ  (Block-Permutation p-Werte)")
    plt.xlabel("κ")
    plt.ylabel("p-Wert")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Plot stats
    plt.figure(figsize=(12,6))
    for L in L_list:
        s = np.array([r[f"s{L}"] for r in rows])
        plt.plot(kappas, s, marker="o", linewidth=1.2, label=f"stat{L}")
    plt.title("Resonanz-Stärke vs. κ  (|E[e^{iΦ_L}]|)")
    plt.xlabel("κ")
    plt.ylabel("stat")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.show()

# -------------------------
# RUN
# -------------------------
rows = kappa_scan(
    limit=100000,
    window=5000,
    start_synth=137,
    L_list=(120,137,240,274),
    kappas=np.linspace(0.0, 2.0, 21),  # z.B. 0.0,0.1,...,2.0
    n_perm=120,
    block=200,
    only_ABC=True,
    mode="bias_coupled"
)

best137 = pick_best(rows, L=137)
print("\nBEST (nach p137):")
print(best137)

plot_scan(rows, L_list=(120,137,240,274))
