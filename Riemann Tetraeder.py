import importlib.util
import os

import numpy as np

_ROOT = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "gral_miller_rabin", os.path.join(_ROOT, "Gral Miller Rabin.py")
)
_gral = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_gral)
SuperResAutopilot_64 = _gral.SuperResAutopilot_64


class TetrahedralResonanceEngine(SuperResAutopilot_64):
    def measure_tetra_resonance(self, a_val, b_val, c_val):
        """
        Berechnet die Resonanz basierend auf der tetraedrischen
        Verschränkung der drei Achsen.
        """
        prod = float(a_val) * float(b_val) * float(c_val)
        if prod <= 0.0 or not np.isfinite(prod):
            return 0.0

        volume_anchor = np.log(prod) / 3.0

        n = min(2_000_000, self.zeros.size)
        if n == 0:
            return 0.0
        gamma = self.zeros[:n]
        resonance = np.sum(np.cos(gamma * self.ln_N) * np.cos(gamma * volume_anchor))

        return float(resonance)

    def explore_geometry(self, center_x, search_steps=50):
        print(f"[*] Sondiere tetraedrischen Raum um {center_x}")
        print(f"{'Achse A':<15} | {'Achse B':<15} | {'Resonanz-Stabilität'}")
        print("-" * 55)

        best_stability = -np.inf
        cx = float(center_x)

        for i in range(-search_steps, search_steps):
            a_axis = cx + i
            if a_axis <= 1.0 or not np.isfinite(a_axis):
                continue
            b_axis = float(self.N) / a_axis
            if b_axis <= 0.0 or not np.isfinite(b_axis):
                continue
            c_axis = float(np.sqrt(self.N))

            res = self.measure_tetra_resonance(a_axis, b_axis, c_axis)

            if res > best_stability:
                best_stability = res

            if i % 10 == 0:
                print(f"{a_axis:<15.0f} | {b_axis:<15.2f} | {res:.6f}")

            ia = int(np.round(a_axis))
            if ia > 1 and self.N % ia == 0:
                print("\n" + "!" * 60)
                print("GEOMETRISCHER LOCK! Tetraeder-Achse A ist ein Faktor.")
                print("!" * 60)
                return ia

        return None


if __name__ == "__main__":
    N_demo = 982_634_309
    sqrt_n = float(np.sqrt(N_demo))
    eng = TetrahedralResonanceEngine("zeros6.npy", N_demo)
    eng.explore_geometry(sqrt_n, search_steps=30)
