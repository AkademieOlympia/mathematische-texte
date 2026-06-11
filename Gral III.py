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


def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(np.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True


class BambergerGralFinal:
    def __init__(self, zeros_path, N):
        path = _resolve_zeros_path(zeros_path)
        self.zeros = np.load(path).astype(float).ravel()
        self.zeros = self.zeros[np.isfinite(self.zeros)]
        self.N = float(N)
        self.ln_N = np.log(self.N)

    def scan_region(self, start, end):
        print(f"[*] Scanne Bereich {start} bis {end} für N={self.N}")
        print(f"{'x':<10} | {'Gain':<8} | {'Phase':<8} | {'Prim?':<6} | {'STATUS'}")
        print("-" * 60)

        for x in range(start, end + 1):
            c_low = self.calculate_c(x, 100_000)
            c_high = self.calculate_c(x, 2_000_000)

            eps = 1e-15
            if abs(c_low) > eps:
                gain = abs(c_high / c_low)
            else:
                gain = 0.0 if abs(c_high) <= eps else float("inf")

            phase_stable = np.sign(c_low) == np.sign(c_high)
            prim = is_prime(x)

            if phase_stable and gain > 10:
                status = "!!! ZIEL GEFUNDEN !!!" if prim else "Harmonischer Geist"
                prim_str = "JA" if prim else "NEIN"
                gain_str = f"{gain:>7.2f}x" if np.isfinite(gain) else "    infx"
                print(f"{x:<10} | {gain_str:<8} | STABIL   | {prim_str:<6} | {status}")

    def calculate_c(self, x, pwr):
        x = float(x)
        if x <= 0.0 or not np.isfinite(x):
            return 0.0
        ln_x = np.log(x)
        n = max(0, min(int(pwr), self.zeros.size))
        if n == 0:
            return 0.0
        gamma = self.zeros[:n]
        return float(np.sum(np.cos(gamma * self.ln_N) * np.cos(gamma * ln_x)) / np.sqrt(x))


if __name__ == "__main__":
    N = 31337 * 31357
    gral = BambergerGralFinal("zeros6.npy", N)
    gral.scan_region(31300, 31400)
