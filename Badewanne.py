# Beispiel: N = 10007 * 10009 = 100160063
# (Repariert: Imports, zeros, schnelle Summation, kein blockierendes plt.show standardmäßig)
import os
import sys

import numpy as np
import matplotlib

_SHOW = os.environ.get("BADEWANNE_SHOW", os.environ.get("HEUREKA_SHOW", "0")) == "1"
if not _SHOW:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ZEROS_PATH = os.path.join(os.path.dirname(__file__), "zeros6.npy")
if os.path.isfile(_ZEROS_PATH):
    zeros = np.load(_ZEROS_PATH).astype(float).ravel()
    zeros = np.sort(zeros[np.isfinite(zeros)])
else:
    try:
        from mpmath import zetazero

        _n_fallback = int(os.environ.get("BADEWANNE_MP_ZEROS", "2000"))
        print(
            f"[Badewanne] {_ZEROS_PATH} fehlt — berechne {_n_fallback} γ mit mpmath …",
            file=sys.stderr,
        )
        zeros = np.array(
            [float(zetazero(k + 1).imag) for k in range(_n_fallback)],
            dtype=float,
        )
    except ImportError:
        print("[Badewanne] Weder zeros6.npy noch mpmath — nutze 40 feste γ.", file=sys.stderr)
        zeros = np.array(
            [
                14.134725,
                21.022040,
                25.010858,
                30.424876,
                32.935062,
                37.586178,
                40.918719,
                43.327073,
                48.005151,
                49.773832,
                52.970321,
                56.446247,
                59.347044,
                60.831779,
                65.112544,
                67.079811,
                69.546402,
                72.067158,
                75.704691,
                77.144840,
                79.337375,
                82.910380,
                84.735493,
                87.425275,
                88.809111,
                92.491899,
                94.651344,
                95.870634,
                98.831194,
                101.317851,
                103.725538,
                105.446623,
                107.168611,
                111.029535,
                111.874659,
                114.320220,
                116.226680,
                118.790782,
                121.370125,
                122.946829,
            ],
            dtype=float,
        )

N_test = 100160063
NUM_ZEROS = int(os.environ.get("BADEWANNE_NUM_ZEROS", "500000"))
NUM_ZEROS = min(NUM_ZEROS, len(zeros))
gammas_use = zeros[:NUM_ZEROS]

x_scan = np.linspace(2, int(np.sqrt(N_test)) + 100, 5000)


def spectral_probe(x, N, gammas, chunk_size=4000):
    x = np.asarray(x, dtype=float)
    gammas = np.asarray(gammas, dtype=float)
    ln_N = np.log(N)
    log_x = np.log(x)
    correlation = np.zeros_like(x, dtype=float)
    n = len(gammas)
    report_every = max(1, n // 10)

    for start in range(0, n, chunk_size):
        g = gammas[start : start + chunk_size][:, np.newaxis]
        correlation += np.sum(np.cos(g * log_x) * np.cos(g * ln_N), axis=0)
        if (start + chunk_size) % report_every < chunk_size or start + chunk_size >= n:
            pct = min(100, int(100 * min(start + chunk_size, n) / n))
            print(f"[Badewanne] {pct}% ({min(start + chunk_size, n)}/{n} γ)", flush=True)

    return correlation / (np.sqrt(x) * np.log(N))


if __name__ == "__main__":
    print(
        f"[Badewanne] N={N_test}, {len(gammas_use)} Nullstellen, {len(x_scan)} x-Punkte …",
        flush=True,
    )
    probe_results = spectral_probe(x_scan, N_test, gammas_use)

    out_png = os.path.join(os.path.dirname(__file__) or ".", "badewanne_spectrum.png")
    plt.figure(figsize=(10, 5))
    plt.plot(x_scan, probe_results)
    plt.axvline(x=10007, color="r", label="Tatsächlicher Faktor p")
    plt.title(f"Spektrale Sonde für N={N_test}")
    plt.xlabel("x")
    plt.ylabel("Korrelation")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"[Badewanne] Fertig: {out_png}", flush=True)
    if _SHOW:
        plt.show()
    plt.close()
