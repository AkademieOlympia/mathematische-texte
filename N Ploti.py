import math
import pandas as pd
import matplotlib.pyplot as plt

CSV_PATH = "arith_struktur_v3_bis_100000.csv"

# ------------------------------------------------------------
# Laden
# ------------------------------------------------------------

df = pd.read_csv(CSV_PATH)

# Sicherheitscheck
required_cols = [
    "n", "sB", "mB", "sigma", "epsilon",
    "kappa_5", "kappa_7", "kappa_11", "kappa_13",
    "fam_support", "rho_E", "rho_A", "rho_B", "rho_C",
    "C_fam", "K_bal", "K_mass", "K_small"
]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Fehlende Spalten in CSV: {missing}")


# ------------------------------------------------------------
# Hilfsfunktionen
# ------------------------------------------------------------

def structural_distance_row(r1, r2, w_eps=1.0, w_rho=1.0, w_kappa=1.0):
    eps_diff = r1["epsilon"] - r2["epsilon"]

    rho_diff_sq = (
        (r1["rho_E"] - r2["rho_E"]) ** 2
        + (r1["rho_A"] - r2["rho_A"]) ** 2
        + (r1["rho_B"] - r2["rho_B"]) ** 2
        + (r1["rho_C"] - r2["rho_C"]) ** 2
    )

    kap_diff_sq = (
        (r1["kappa_5"] - r2["kappa_5"]) ** 2
        + (r1["kappa_7"] - r2["kappa_7"]) ** 2
        + (r1["kappa_11"] - r2["kappa_11"]) ** 2
        + (r1["kappa_13"] - r2["kappa_13"]) ** 2
    )

    return math.sqrt(
        w_eps * eps_diff**2 + w_rho * rho_diff_sq + w_kappa * kap_diff_sq
    )


# ------------------------------------------------------------
# Plot 1: Scatter epsilon vs C_fam
# ------------------------------------------------------------

def plot_scatter():
    sizes = 10 + 120 * (df["K_mass"] / df["K_mass"].max())

    plt.figure(figsize=(10, 7))
    sc = plt.scatter(
        df["epsilon"],
        df["C_fam"],
        c=df["kappa_13"],
        s=sizes,
        alpha=0.6
    )
    plt.xlabel("epsilon")
    plt.ylabel("C_fam")
    plt.title("Strukturraum: epsilon vs. C_fam (Farbe = kappa_13, Größe = K_mass)")
    cbar = plt.colorbar(sc)
    cbar.set_label("kappa_13")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("plot_scatter_epsilon_Cfam.png", dpi=200)
    plt.close()


# ------------------------------------------------------------
# Plot 2: Faserplots
# ------------------------------------------------------------

def plot_fiber(kernel_value):
    fib = df[df["mB"] == kernel_value].copy().sort_values(by="n")
    if fib.empty:
        print(f"Keine Daten für Faser mB={kernel_value}")
        return

    plt.figure(figsize=(10, 7))
    plt.plot(fib["sigma"], fib["K_bal"], marker="o", label="K_bal")
    plt.plot(fib["sigma"], fib["K_mass"], marker="s", label="K_mass")
    plt.plot(fib["sigma"], fib["K_small"], marker="^", label="K_small")

    for _, row in fib.iterrows():
        if row["n"] in [kernel_value, 420, 693, 945, 143]:
            plt.annotate(
                str(int(row["n"])),
                (row["sigma"], row["K_bal"]),
                fontsize=8,
                xytext=(4, 4),
                textcoords="offset points"
            )

    plt.xlabel("sigma")
    plt.ylabel("Krümmungsmaß")
    plt.title(f"Restkern-Faser mB={kernel_value}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"plot_fiber_mB_{kernel_value}.png", dpi=200)
    plt.close()


# ------------------------------------------------------------
# Plot 3: Distanz-Heatmap für eine Faser
# ------------------------------------------------------------

def plot_fiber_heatmap(kernel_value, max_points=20):
    fib = df[df["mB"] == kernel_value].copy().sort_values(by="n")
    if fib.empty:
        print(f"Keine Daten für Heatmap mB={kernel_value}")
        return

    fib = fib.head(max_points).reset_index(drop=True)
    labels = fib["n"].astype(int).astype(str).tolist()

    m = len(fib)
    dist = [[0.0 for _ in range(m)] for _ in range(m)]

    for i in range(m):
        for j in range(m):
            dist[i][j] = structural_distance_row(fib.iloc[i], fib.iloc[j])

    plt.figure(figsize=(9, 7))
    plt.imshow(dist, aspect="auto")
    plt.colorbar(label="d_str3")
    plt.xticks(range(m), labels, rotation=90)
    plt.yticks(range(m), labels)
    plt.title(f"Strukturdistanz-Heatmap für Faser mB={kernel_value}")
    plt.tight_layout()
    plt.savefig(f"plot_heatmap_mB_{kernel_value}.png", dpi=200)
    plt.close()


# ------------------------------------------------------------
# Plot 4: Optionaler Kappa-Pfad entlang einer Faser
# ------------------------------------------------------------

def plot_fiber_kappas(kernel_value):
    fib = df[df["mB"] == kernel_value].copy().sort_values(by="n")
    if fib.empty:
        print(f"Keine Daten für Kappa-Faser mB={kernel_value}")
        return

    plt.figure(figsize=(10, 7))
    plt.plot(fib["sigma"], fib["kappa_5"], marker="o", label="kappa_5")
    plt.plot(fib["sigma"], fib["kappa_7"], marker="s", label="kappa_7")
    plt.plot(fib["sigma"], fib["kappa_11"], marker="^", label="kappa_11")
    plt.plot(fib["sigma"], fib["kappa_13"], marker="d", label="kappa_13")

    plt.xlabel("sigma")
    plt.ylabel("kappa")
    plt.title(f"Mehrskalen-Kappa entlang der Faser mB={kernel_value}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"plot_kappas_mB_{kernel_value}.png", dpi=200)
    plt.close()


# ------------------------------------------------------------
# Hauptlauf
# ------------------------------------------------------------

if __name__ == "__main__":
    plot_scatter()

    for kernel in [35, 77, 143]:
        plot_fiber(kernel)
        plot_fiber_heatmap(kernel, max_points=20)
        plot_fiber_kappas(kernel)

    print("Plots geschrieben:")
    print(" - plot_scatter_epsilon_Cfam.png")
    print(" - plot_fiber_mB_35.png")
    print(" - plot_fiber_mB_77.png")
    print(" - plot_fiber_mB_143.png")
    print(" - plot_heatmap_mB_35.png")
    print(" - plot_heatmap_mB_77.png")
    print(" - plot_heatmap_mB_143.png")
    print(" - plot_kappas_mB_35.png")
    print(" - plot_kappas_mB_77.png")
    print(" - plot_kappas_mB_143.png")