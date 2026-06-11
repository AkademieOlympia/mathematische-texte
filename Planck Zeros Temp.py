import numpy as np
import math
from numpy.polynomial import Polynomial

# Konstanten (wie in deinem Skript)
c_val = 299792458.0
kB = 1.380649e-23
e = 1.602176634e-19
R_K = 25812.80745
h = (e**2) * R_K

def planck_u_nu(nu, T):
    nu = np.asarray(nu, dtype=np.float64)
    x = (h * nu) / (kB * T)
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        factor = (8*np.pi*h*nu**3)/(c_val**3)
        return factor/(np.exp(x)-1.0)

def build_RH_template(z_grid, zeros, weights=None, chunk=500):
    z = np.asarray(z_grid, dtype=np.float64)
    if weights is None:
        weights = np.ones_like(zeros)
    weights = np.asarray(weights, dtype=np.float64)
    M = np.zeros_like(z)
    for i in range(0, len(zeros), chunk):
        g = zeros[i:i+chunk][:, None]
        w = weights[i:i+chunk][:, None]
        M += (w * np.cos(g * z[None, :])).sum(axis=0)
    return M

def normalize01(a, eps=1e-12):
    a = np.asarray(a, dtype=np.float64)
    return (a - a.min())/(np.ptp(a) + eps)

def estimate_temperature_from_zeros(nu, S_obs, zeros, T_grid,
                                   z_grid_points=6000,
                                   use_weights=True):
    """
    Gibt T_hat und scores(T) zurück.
    nu: Frequenzen [Hz]
    S_obs: Spektrum (proportional zu u_nu ist ok)
    zeros: array der Imaginärteile gamma (z.B. erste 10k)
    T_grid: Kandidaten-Temperaturen [K]
    """

    nu = np.asarray(nu, dtype=np.float64)
    S_obs = np.asarray(S_obs, dtype=np.float64)

    w = np.log(nu)
    const = math.log(h/kB)

    # Dämpfung hoher gamma (robuster als "alle gleich")
    if use_weights:
        weights = 1.0/np.sqrt(0.25 + zeros**2)
    else:
        weights = None

    # z-Bereich so wählen, dass alle T-Shifts abgedeckt sind
    z_min = (w.min() + const - math.log(np.max(T_grid))) - 0.2
    z_max = (w.max() + const - math.log(np.min(T_grid))) + 0.2
    z_grid = np.linspace(z_min, z_max, z_grid_points)

    # RH-Template M(z) einmalig vorrechnen
    M_grid = build_RH_template(z_grid, zeros, weights=weights, chunk=500)

    scores = []
    for T in T_grid:
        u = planck_u_nu(nu, T)
        R = (S_obs / (u + 1e-300)) - 1.0

        # Glatten Baseline-Anteil entfernen (u_true/u - 1 bei falschem T)
        if len(R) > 10:
            poly = Polynomial.fit(np.linspace(0, 1, len(R)), R, deg=5)
            R = R - poly(np.linspace(0, 1, len(R)))

        z = w + const - math.log(T)
        Mt = np.interp(z, z_grid, M_grid)

        Rt = R - R.mean()
        Mt0 = Mt - Mt.mean()
        denom = (np.linalg.norm(Rt) * np.linalg.norm(Mt0) + 1e-12)
        scores.append(float(np.dot(Rt, Mt0) / denom))

    scores = np.array(scores)
    T_hat = float(T_grid[np.argmax(scores)])
    return T_hat, scores


if __name__ == "__main__":
    import os
    import matplotlib.pyplot as plt

    # Riemann-Nullstellen: gamma_zeros_100k.txt, gamma_zeros_100.txt oder zeros_gamma.txt
    script_dir = os.path.dirname(os.path.abspath(__file__)) or "."
    max_zeros_use = 10_000  # bis zehntausendste Nullstelle
    zeros_path = None
    for name in ("gamma_zeros_100k.txt", "gamma_zeros_100.txt", "zeros_gamma.txt"):
        p = os.path.join(script_dir, name)
        if os.path.isfile(p):
            zeros_path = p
            break

    if zeros_path:
        zeros_all = np.loadtxt(zeros_path, dtype=np.float64)
        if zeros_all.ndim > 1:
            zeros_all = zeros_all.ravel()
        n_total = len(zeros_all)
        zeros = zeros_all[:max_zeros_use]
        print(f"{os.path.basename(zeros_path)}: {n_total:,} Nullstellen, davon {len(zeros):,} verwendet (bis 10.000.).")
        print("Erste 5 γ:", zeros[:5].tolist())
    else:
        # Fallback: 30 bekannte + 970 approx = 1000 (bessere Schätzung als nur 30)
        known = np.array([
            14.1347, 21.0220, 25.0109, 30.4249, 32.9351, 37.5862, 40.9187, 43.3271,
            48.0052, 49.7738, 52.9703, 56.4462, 59.3470, 60.8318, 65.1125, 67.0798,
            69.5464, 72.0672, 75.7047, 77.1448, 79.3374, 82.9104, 84.7355, 87.4253,
            88.8091, 92.4919, 94.6513, 95.8706, 98.8312, 101.3179
        ], dtype=np.float64)
        n = np.arange(31, max_zeros_use + 1, dtype=np.float64)
        approx = (2 * np.pi * n) / np.log(n)
        zeros = np.concatenate([known, approx])
        print(f"Keine Nullstellen-Datei – verwende {len(zeros)} approximate γ-Werte (bis Index 10.000).")

    nu = np.linspace(1e12, 800e12, 2000)
    w = np.log(nu)
    const = math.log(h / kB)
    T_grid = np.linspace(100, 15_000, 300)

    # Kurven für T_true = 1000 K und 5000 K extrahieren
    T_targets = [1000, 5000]
    curves = {}

    for T_true in T_targets:
        S_obs = planck_u_nu(nu, T_true)
        z_test = w + const - math.log(T_true)
        M_test = build_RH_template(z_test, zeros, weights=1.0 / np.sqrt(0.25 + zeros**2), chunk=500)
        M_norm = normalize01(M_test)
        S_obs = S_obs * (1.0 + 0.03 * (M_norm - 0.5))

        T_hat, scores = estimate_temperature_from_zeros(nu, S_obs, zeros, T_grid,
                                                         z_grid_points=4000, use_weights=True)
        curves[T_true] = (T_grid.copy(), scores.copy(), T_hat)

        csv_path = os.path.join(script_dir, f"planck_korrelation_T{T_true}K.csv")
        np.savetxt(csv_path, np.column_stack([T_grid, scores]),
                   header="T_K,Korrelation", delimiter=",", comments="")
        print(f"T_true={T_true} K: T_hat={T_hat:.1f} K, Korr_max={scores.max():.4f}")
        print(f"  → Kurve gespeichert: {csv_path}")

    # Plot: beide Kurven
    fig, ax = plt.subplots(figsize=(10, 5))
    for T_true in T_targets:
        T_grid_cur, scores_cur, T_hat = curves[T_true]
        ax.plot(T_grid_cur, scores_cur, lw=1.2, label=f"T_true = {T_true} K")
        ax.axvline(T_true, color=ax.get_lines()[-1].get_color(), ls="--", alpha=0.6)
        ax.axvline(T_hat, color=ax.get_lines()[-1].get_color(), ls=":", alpha=0.8)

    ax.set_xlabel("T (K)")
    ax.set_ylabel("Korrelation mit RH-Template")
    ax.set_title("Temperatur-Schätzung: Korrelationskurven für T_true = 1000 K und 5000 K (10.000 Nullstellen)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    grafik_dir = os.path.join(script_dir, "Grafik")
    os.makedirs(grafik_dir, exist_ok=True)

    out_path = os.path.join(grafik_dir, "planck_temp_curves_1000_5000.png")
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close()
    out_path_abs = os.path.abspath(out_path)
    if os.path.isfile(out_path_abs):
        print(f"Plot gespeichert: {out_path_abs}")
    else:
        print(f"FEHLER: Datei nicht erstellt: {out_path_abs}")

    # Planck-Strahlungskurven für verschiedene Temperaturen (u_nu vs ν)
    T_planck = [1000, 2500, 5000, 7500]  # K
    nu_plot = np.linspace(0.1e12, 1200e12, 2000)

    fig2, (ax2a, ax2b) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    # Oben: absolute Kurven (log-Skala, damit alle sichtbar)
    for T in T_planck:
        u = planck_u_nu(nu_plot, T)
        ax2a.plot(nu_plot / 1e12, u, lw=1.5, label=f"T = {T} K")
    ax2a.set_yscale("log")
    ax2a.set_ylabel("Energiedichte uν (J·s/m³)")
    ax2a.set_title("Planck-Strahlungskurven bei verschiedenen Temperaturen (absolut)")
    ax2a.legend()
    ax2a.grid(True, alpha=0.3, which="both")
    ax2a.set_xlim(0, 1200)
    ax2a.set_ylim(bottom=1e-20)

    # Unten: normalisiert (max=1) – Formänderung bei T
    for T in T_planck:
        u = planck_u_nu(nu_plot, T)
        u_norm = u / (u.max() + 1e-300)
        ax2b.plot(nu_plot / 1e12, u_norm, lw=1.5, label=f"T = {T} K")
    ax2b.set_xlabel("Frequenz ν (THz)")
    ax2b.set_ylabel("uν / max(uν)")
    ax2b.set_title("Planck-Strahlungskurven (normalisiert 0–1)")
    ax2b.legend()
    ax2b.grid(True, alpha=0.3)
    ax2b.set_xlim(0, 1200)
    ax2b.set_ylim(0, 1.05)
    plt.tight_layout()
    out_planck = os.path.join(grafik_dir, "planck_strahlungskurven.png")
    plt.savefig(out_planck, dpi=120, bbox_inches="tight")
    plt.close()
    out_planck_abs = os.path.abspath(out_planck)
    if os.path.isfile(out_planck_abs):
        print(f"Planck-Strahlungskurven gespeichert: {out_planck_abs}")
    else:
        print(f"FEHLER: Datei nicht erstellt: {out_planck_abs}")

    # Auch im Hauptordner ablegen (Fallback)
    for fname in ("planck_temp_curves_1000_5000.png", "planck_strahlungskurven.png"):
        src = os.path.join(grafik_dir, fname)
        dst = os.path.join(script_dir, fname)
        if os.path.isfile(src):
            import shutil
            shutil.copy2(src, dst)
            print(f"Kopie im Hauptordner: {os.path.abspath(dst)}")