import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


CSV_PATH = Path("trajectories.csv")


def load_trajectories(path: Path):
    trajectories = defaultdict(
        lambda: {
            "step": [],
            "offset": [],
            "energy": [],
            "gradient": [],
            "target_lock_q_step": [],
            "target_lock_q_offset": [],
        }
    )

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            probe_id = int(row["ProbeID"])
            step = int(row["Step"])
            offset = float(row["Offset"])
            trajectories[probe_id]["step"].append(step)
            trajectories[probe_id]["offset"].append(offset)
            trajectories[probe_id]["energy"].append(float(row["Energy"]))
            trajectories[probe_id]["gradient"].append(float(row["Gradient"]))

            if row.get("TargetLockQ", "0") == "1" or row.get("Status") == "target_lock_q":
                trajectories[probe_id]["target_lock_q_step"].append(step)
                trajectories[probe_id]["target_lock_q_offset"].append(offset)

    return dict(sorted(trajectories.items()))


def save_plot(trajectories, value_key: str, ylabel: str, filename: str, log_y: bool = False):
    fig, ax = plt.subplots(figsize=(10, 6))

    for probe_id, values in trajectories.items():
        ax.plot(values["step"], values[value_key], label=f"Sonde {probe_id}", linewidth=1.2)
        if value_key == "offset" and values["target_lock_q_step"]:
            ax.scatter(
                values["target_lock_q_step"],
                values["target_lock_q_offset"],
                marker="o",
                s=30,
                label=f"Sonde {probe_id} target_lock_q",
            )

    ax.set_title(f"{ylabel} pro Schritt")
    ax.set_xlabel("Schritt")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)

    if log_y:
        ax.set_yscale("log")

    if len(trajectories) <= 16:
        ax.legend(loc="best", fontsize=8, ncol=2)

    fig.tight_layout()
    fig.savefig(filename, dpi=200)
    plt.close(fig)


def save_gradient_magnitude_plot(trajectories, filename: str):
    fig, ax = plt.subplots(figsize=(10, 6))

    for probe_id, values in trajectories.items():
        magnitudes = [abs(value) for value in values["gradient"]]
        ax.plot(values["step"], magnitudes, label=f"Sonde {probe_id}", linewidth=1.2)

    ax.set_title("Betrag des Gradienten pro Schritt")
    ax.set_xlabel("Schritt")
    ax.set_ylabel("|Gradient|")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)

    if len(trajectories) <= 16:
        ax.legend(loc="best", fontsize=8, ncol=2)

    fig.tight_layout()
    fig.savefig(filename, dpi=200)
    plt.close(fig)


def extract_root_cluster(trajectories, cluster_size: int = 5):
    ordered = sorted(
        trajectories.items(),
        key=lambda item: item[1]["offset"][0]
    )
    cluster = dict(ordered[-cluster_size:])
    return dict(sorted(cluster.items()))


def save_root_cluster_plot(trajectories, filename: str):
    root_cluster = extract_root_cluster(trajectories)

    fig, ax = plt.subplots(figsize=(10, 6))
    for probe_id, values in root_cluster.items():
        ax.plot(values["step"], values["offset"], label=f"Sonde {probe_id}", linewidth=1.4)

    ax.set_title("Root-Cluster: Offset pro Schritt")
    ax.set_xlabel("Schritt")
    ax.set_ylabel("Offset")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)

    fig.tight_layout()
    fig.savefig(filename, dpi=200)
    plt.close(fig)


def main():
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"{CSV_PATH} wurde nicht gefunden.")

    trajectories = load_trajectories(CSV_PATH)
    if not trajectories:
        raise RuntimeError("Keine Trajektorien in der CSV gefunden.")

    save_plot(trajectories, "offset", "Offset", "trajectories_offset_vs_step.png")
    save_plot(trajectories, "energy", "Energie", "trajectories_energy_vs_step.png", log_y=True)
    save_plot(trajectories, "gradient", "Gradient", "trajectories_gradient_vs_step.png")
    save_gradient_magnitude_plot(trajectories, "trajectories_gradient_abs_vs_step.png")
    save_root_cluster_plot(trajectories, "trajectories_root_cluster_offset_vs_step.png")

    print("[PLOT] Diagramme gespeichert:")
    print("[PLOT] trajectories_offset_vs_step.png")
    print("[PLOT] trajectories_energy_vs_step.png")
    print("[PLOT] trajectories_gradient_vs_step.png")
    print("[PLOT] trajectories_gradient_abs_vs_step.png")
    print("[PLOT] trajectories_root_cluster_offset_vs_step.png")


if __name__ == "__main__":
    main()
