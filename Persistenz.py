import os
import sys

import numpy as np
import matplotlib

# Standard: kein blockierendes Fenster — Plot als PNG.
# Interaktiv: PERSISTENZ_SHOW=1 oder HEUREKA_SHOW=1
_SHOW = os.environ.get("PERSISTENZ_SHOW", os.environ.get("HEUREKA_SHOW", "0")) == "1"
if not _SHOW:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

def field_amplitude(x, scale=1.0):
    target = -1.0 * np.exp(-((x - 80) ** 2) / (10 * scale))
    trap = -0.8 * np.exp(-((x - 40) ** 2) / (1 * scale))  # Sehr schmal
    noise = 0.1 * np.sin(x * 5)  # Die "Maische"
    return target + trap + noise


# Topologische Güte: Integration über verschiedene Skalen (Radien)
def topo_persistence(x_pos):
    # Wir messen die "Lebensdauer" des Merkmals über 3 Skalen
    scales = [0.5, 1.0, 2.0]
    amps = [field_amplitude(x_pos, s) for s in scales]
    # Ein Merkmal ist persistent, wenn die Amplitude bei Vergrößerung der
    # Skale erhalten bleibt (Integrale Stabilität).
    return np.mean(amps) * (1.0 / (1.0 + np.std(amps)))


def run_demo() -> None:
    # Wir simulieren das Feld: Ein tiefes Tal (Ziel) und ein flacheres,
    # aber steiles Tal (die Falle).
    x = np.linspace(0, 100, 1000)

    pos_classic = 30.0
    pos_topo = 30.0
    path_classic = [pos_classic]
    path_topo = [pos_topo]

    for _ in range(50):
        grad_c = (field_amplitude(pos_classic + 0.1) - field_amplitude(pos_classic - 0.1)) / 0.2
        pos_classic -= 2.0 * grad_c
        path_classic.append(pos_classic)

        grad_t = (topo_persistence(pos_topo + 0.1) - topo_persistence(pos_topo - 0.1)) / 0.2
        pos_topo -= 30.0 * grad_t
        path_topo.append(pos_topo)

    plt.figure(figsize=(14, 7))
    plt.plot(x, field_amplitude(x), color="gray", label="Arithmetisches Feld (Maische)")
    plt.plot(path_classic, [field_amplitude(p) for p in path_classic], "ro-", label="Klassische Sonde (Bleibt hängen)")
    plt.plot(path_topo, [field_amplitude(p) for p in path_topo], "go-", label="Topologischer Autopilot (Überfliegt Falle)")
    plt.axvline(x=40, color="orange", linestyle="--", label="Falle (Geringe Persistenz)")
    plt.axvline(x=80, color="green", linestyle="--", label="Ziel: Faktor p (Hohe Persistenz)")
    plt.title("Topologisches Überfliegen: Persistenz vs. Amplitude")
    plt.xlabel("x")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.tight_layout()

    _out = os.path.join(os.path.dirname(os.path.abspath(__file__)) or ".", "persistenz_plot.png")
    plt.savefig(_out, dpi=150, bbox_inches="tight")
    print(f"[Persistenz] Plot gespeichert: {_out}", file=sys.stderr, flush=True)
    if _SHOW:
        plt.show()
    plt.close()


if __name__ == "__main__":
    run_demo()
