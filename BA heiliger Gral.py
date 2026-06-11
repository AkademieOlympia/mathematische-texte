import os

import numpy as np

_ROOT = os.path.dirname(os.path.abspath(__file__))


def _resolve_zeros_path(zeros_path: str) -> str:
    if os.path.isfile(zeros_path):
        return zeros_path
    alt = os.path.join(_ROOT, os.path.basename(zeros_path))
    if os.path.isfile(alt):
        return alt
    raise FileNotFoundError(
        f"Nullstellen-Datei nicht gefunden: {zeros_path!r} (auch nicht {_ROOT})"
    )


# =========================================================================
# DAS BAMBERGER SIEGEL: SPEKTRALE TIEFENBOHRUNG (VERSION 2.0)
# =========================================================================


class BambergerSiegel:
    def __init__(self, zeros_path, N):
        print("[*] Initialisiere Siegel-Engine...")
        path = _resolve_zeros_path(zeros_path)
        self.zeros = np.load(path).astype(float).ravel()
        self.zeros = self.zeros[np.isfinite(self.zeros)]
        print(f"[*] Hardware-Check: {self.zeros.size:,} Nullstellen bereit.")

        self.N = float(N)
        self.ln_N = np.log(self.N)

    def drill(self, x_pos):
        """Führt den spektralen Zoom durch: 100k vs volle Liste (bis 2M)."""
        x_pos = float(x_pos)
        ln_x = np.log(x_pos)
        sqrt_x = np.sqrt(x_pos)

        n_low = min(100_000, self.zeros.size)
        n_high = self.zeros.size
        z_low = self.zeros[:n_low]
        z_high = self.zeros[:n_high]

        c_low = float(np.sum(np.cos(z_low * self.ln_N) * np.cos(z_low * ln_x)) / sqrt_x)
        c_high = float(np.sum(np.cos(z_high * self.ln_N) * np.cos(z_high * ln_x)) / sqrt_x)

        eps = 1e-15
        if abs(c_low) > eps:
            gain = abs(c_high / c_low)
        else:
            gain = 0.0 if abs(c_high) <= eps else float("inf")
        return c_low, c_high, gain

    def execute_scan(self, kandidaten):
        print(f"\nTarget N: {int(self.N)}")
        print(f"{'Kandidat x':<15} | {'100k-Resonanz':<15} | {'2M-Resonanz':<15} | {'Gain'}")
        print("-" * 75)
        
        for x in kandidaten:
            low, high, gain = self.drill(x)
            status = "!!! MONADE !!!" if gain > 10.0 else "Maische"
            gstr = f"{gain:>6.2f}x" if np.isfinite(gain) else "   infx"
            print(f"{x:<15.2f} | {low:<15.6f} | {high:<15.6f} | {gstr}  -> {status}")

# --- DURCHFÜHRUNG ---
if __name__ == "__main__":
    # Wir setzen N auf unser neues Ziel
    N_NEU = 31337 * 31357  # 982634309
    
    # Wir testen die Umgebung: Ein Geist (31300), das Ziel (31337) und ein Pseudo-Peak (31360)
    test_kandidaten = [31300.0, 31337.0, 31360.0]
    
    siegel = BambergerSiegel("zeros6.npy", N_NEU)
    siegel.execute_scan(test_kandidaten)