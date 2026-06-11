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


class HeiligerGralSiegel:
    def __init__(self, zeros_path, N):
        print(f"[*] Kalibrierung auf N = {N}...")
        path = _resolve_zeros_path(zeros_path)
        self.zeros = np.load(path).astype(float).ravel()
        self.zeros = self.zeros[np.isfinite(self.zeros)]
        self.N = float(N)
        self.ln_N = np.log(self.N)

    def messung(self, x, pwr):
        x = float(x)
        ln_x = np.log(x)
        n = max(0, min(int(pwr), self.zeros.size))
        if n == 0:
            return 0.0
        gamma_chunk = self.zeros[:n]
        # C(x) Berechnung
        interferenz = np.sum(np.cos(gamma_chunk * self.ln_N) * np.cos(gamma_chunk * ln_x))
        return float(interferenz / np.sqrt(x))

    def check_truth(self, kandidaten):
        print(f"\n{'Kandidat x':<12} | {'Gain':<8} | {'Phase-Lock':<12} | {'Status'}")
        print("-" * 60)
        
        for x in kandidaten:
            c_low = self.messung(x, 100_000)
            c_high = self.messung(x, 2_000_000)
            eps = 1e-15
            gain = abs(c_high / c_low) if abs(c_low) > eps else (0.0 if abs(c_high) <= eps else float("inf"))
            # Phase-Lock: Haben beide Messungen das gleiche Vorzeichen?
            phase_stable = np.sign(c_low) == np.sign(c_high)
            
            status = "REALER FAKTOR" if (gain > 15 and phase_stable) else "MAISCHE-GEIST"

            lock_str = "STABIL" if phase_stable else "FLACKERT"
            gain_str = f"{gain:>7.2f}x" if np.isfinite(gain) else "    infx"
            print(f"{x:<12.2f} | {gain_str:<8} | {lock_str:<12} | {status}")

if __name__ == "__main__":
    # N = 31337 * 31357 (982634309)
    N_TRUE = 31337 * 31357

    kandidaten = [31300.0, 31337.0, 31360.0]

    gral = HeiligerGralSiegel("zeros6.npy", N_TRUE)
    gral.check_truth(kandidaten)