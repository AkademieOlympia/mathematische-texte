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


class MatchMakerSiegel:
    def __init__(self, zeros_path, N):
        path = _resolve_zeros_path(zeros_path)
        self.zeros = np.load(path).astype(float).ravel()
        self.zeros = self.zeros[np.isfinite(self.zeros)]
        self.N = int(N)
        self.ln_N = np.log(float(self.N))
        self.kandidaten = []

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

    def scan_and_lock(self, start, end):
        print(f"[*] Starte Finalen Match-Maker Scan für N = {self.N}")
        print(f"[*] Bereich: {start} bis {end}")
        print("-" * 65)

        for x in range(start, end + 1):
            c_low = self.calculate_c(x, 100_000)
            c_high = self.calculate_c(x, 2_000_000)

            eps = 1e-15
            if abs(c_low) > eps:
                gain = abs(c_high / c_low)
            else:
                gain = 0.0 if abs(c_high) <= eps else float("inf")

            phase_stable = np.sign(c_low) == np.sign(c_high)

            if phase_stable and gain > 10.0 and is_prime(x):
                gstr = f"{gain:.2f}x" if np.isfinite(gain) else "inf"
                print(f"[!] POTENZIELLER FAKTOR GEFUNDEN: {x} (Gain: {gstr})")
                self.kandidaten.append(x)

        print("-" * 65)
        print(f"[*] Prüfe Verschränkung für {len(self.kandidaten)} Kandidaten...")

        for i in range(len(self.kandidaten)):
            for j in range(i, len(self.kandidaten)):
                p, q = self.kandidaten[i], self.kandidaten[j]
                if p * q == self.N:
                    print("\n" + "=" * 65)
                    print("BINGO! RSA-SCHLÜSSEL GEKNACKT!")
                    print(f"Faktor p: {p}")
                    print(f"Faktor q: {q}")
                    print(f"Verifikation: {p} * {q} = {p * q}")
                    print("=" * 65)
                    return p, q

        print("[?] Kein perfektes Match gefunden. Scan-Bereich oder N prüfen.")
        return None


if __name__ == "__main__":
    # DAS ABSOLUT KORREKTE N (Faktoren-Produkt)
    N_TRUE = 982634309  # Hier lag der 1000er-Hund begraben

    engine = MatchMakerSiegel("zeros6.npy", N_TRUE)
    engine.scan_and_lock(31300, 31400)
