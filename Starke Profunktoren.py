"""
Signalling-Flux entlang Riemann-Nullstellen (ohne Sage: nutzt ``HigherOrder``).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
import HigherOrder as ho
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 1. Daten laden (Pfad neben dieser Datei)
zeros = np.load(_ROOT / "zeros6.npy")
sample_size = 10000
test_zeros = zeros[:sample_size]


def _prime_op_choi(p: int, phase: "ho.Quaternion") -> ho.QMatrix:
    """
    2×2-Operator [[p, phase], [phase*, 1]] in H, wie im ursprünglichen Sage-Skript.
    """
    return ho.QMatrix(
        [
            [ho.quat(p, 0, 0, 0), phase],
            [phase.conjugate(), ho.quat(1, 0, 0, 0)],
        ]
    )


def get_signalling_flux_at_t(t: float, p1: int = 17, p2: int = 19) -> float:
    """
    Signalling-Flux zwischen p1, p2 an der Gitterposition t (Höhe der Nullstelle).
    Par-Produkt (A*⊗B*)*, Block-Differenz A−D, normiert mit ||phase||.
    """
    u = float(t)
    phase = ho.quat(0.5, u, u**2, math.sin(u))
    op1 = _prime_op_choi(p1, phase)
    op2 = _prime_op_choi(p2, phase)
    par_res = ho.apply_par_product(op1, op2)
    a_block = par_res.submatrix(0, 0, 2, 2)
    d_block = par_res.submatrix(2, 2, 2, 2)
    num = (a_block - d_block).frobenius_norm()
    ph = math.sqrt(float(phase.reduced_norm()))
    if ph < 1e-30:
        return 0.0
    return float(num / ph)


# 2. Großskalige Berechnung
flux_values = [get_signalling_flux_at_t(float(t)) for t in test_zeros]

# 3. Phasenübergang (Rolling Varianz)
df = pd.DataFrame({"t": test_zeros, "flux": flux_values})
df["flux_rolling_std"] = df["flux"].rolling(window=100).std()

# 4. Visualisierung
fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12, 7))
ax0.plot(df["t"], df["flux"], color="blue", alpha=0.6, label="Signalling Flux")
ax0.set_title(f"Kausale Dynamik über {sample_size} Riemann-Nullstellen")
ax0.set_ylabel("Flux (Verschränkungs-Energie)")
ax0.grid(True)
ax0.legend()

ax1.plot(
    df["t"],
    df["flux_rolling_std"],
    color="red",
    label="Stabilität (Rolling Std)",
)
ax1.set_ylabel("Varianz (Phasensignal)")
ax1.set_xlabel("Imaginärteil t (Frequenz)")
ax1.grid(True)
ax1.legend()

fig.tight_layout()
out_png = _ROOT / "riemann_phase_transition.png"
fig.savefig(out_png, dpi=150)
plt.close(fig)
print(
    f"Analyse abgeschlossen. Ergebnisse in '{out_png.name}' gespeichert "
    f"({out_png})"
)
