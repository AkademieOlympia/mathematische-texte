import numpy as np
from scipy.optimize import curve_fit

# Konstanten (SI)
h = 6.62607015e-34
k_B = 1.380649e-23
c = 299792458

def planck_law(nu, T):
    """Klassische Planck-Kurve."""
    return (8 * np.pi * h * nu**3 / c**3) / (np.exp(h * nu / (k_B * T)) - 1)

def modulation_template(z, weights, frequencies):
    """Erzeugt das oszillatorische M(z) aus einer Frequenzliste (Nullstellen oder ln p)."""
    return sum(w * np.cos(f * z) for w, f in zip(weights, frequencies))

# Beispielhafter Workflow für die Residuen-Analyse
def analyze_residuals(nu_data, energy_data, T_guess, zeros_list):
    # 1. Fit der Temperatur T (Baseline)
    popt, _ = curve_fit(planck_law, nu_data, energy_data, p0=[T_guess])
    T_fit = popt[0]
    
    # 2. Extraktion der Residuen
    planck_base = planck_law(nu_data, T_fit)
    residuals = (energy_data - planck_base) / planck_base
    
    # 3. Transformation in den z-Raum (log-Frequenz)
    z_data = np.log(nu_data)
    
    return z_data, residuals, T_fit

# Hinweis: zeros_list würde hier die imaginärteile der Riemann-Nullstellen enthalten.


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import os

    # Beispiel: Synthetische Planck-Daten bei T = 1000 K
    T_true = 1000.0
    nu_data = np.linspace(1e12, 800e12, 2000)
    energy_data = planck_law(nu_data, T_true)

    # Optional: kleines Rauschen (1 %) für realistischere Residuen
    rng = np.random.default_rng(42)
    energy_data *= 1.0 + 0.01 * rng.standard_normal(len(energy_data))

    zeros_list = []  # Platzhalter (für Erweiterung mit Riemann-Nullstellen)

    z_data, residuals, T_fit = analyze_residuals(
        nu_data, energy_data, T_guess=T_true * 1.05, zeros_list=zeros_list
    )

    print(f"T_true = {T_true} K")
    print(f"T_fit  = {T_fit:.2f} K")
    print(f"Residuen: min={residuals.min():.6f}, max={residuals.max():.6f}, std={residuals.std():.6f}")

    # Plot: Residuen vs. z (log-Frequenz)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(z_data, residuals, "b.", ms=0.8, alpha=0.7)
    ax.axhline(0, color="gray", ls="--", alpha=0.5)
    ax.set_xlabel("z = ln(ν)")
    ax.set_ylabel("Relative Residuen (S/S_Planck − 1)")
    ax.set_title("Planck-Residuen nach Temperatur-Fit (T_true = 1000 K)")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    grafik_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Grafik")
    os.makedirs(grafik_dir, exist_ok=True)
    out_path = os.path.join(grafik_dir, "planck_residuen.png")
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Plot gespeichert: {out_path}")