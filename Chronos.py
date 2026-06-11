import os

from mpmath import mp
import numpy as np

# Höchste Präzision für den 256-Bit Raum
mp.dps = 80


class ChronosEngine_256:
    def __init__(self, N, zeros, pilot_zeros: int = 100_000):
        self.N = mp.mpf(N)
        self.zeros = np.asarray(zeros, dtype=np.float64).ravel()
        self.zeros = self.zeros[np.isfinite(self.zeros)]
        self._pilot_n = min(int(pilot_zeros), self.zeros.size)
        self._ln_N = float(mp.log(self.N))
        self.tau = mp.mpc(0, 1)

    def tetra_resonance(self, x):
        """Summe cos(γ ln x) cos(γ ln N) — vektorisiert (numpy), sonst Minuten pro Aufruf."""
        ln_x = float(mp.log(mp.mpf(x)))
        n = self._pilot_n
        if n == 0:
            return 0.0
        g = self.zeros[:n]
        return float(np.sum(np.cos(g * ln_x) * np.cos(g * self._ln_N)))

    def get_phase_velocity(self, k):
        x0 = mp.mpf(210 * k + 13)
        x1 = mp.mpf(210 * (k + 1) + 13)
        res0 = self.tetra_resonance(x0)
        res1 = self.tetra_resonance(x1)
        velocity = res1 - res0
        return float(velocity), res0

    def final_lock(self, k):
        # Wir prüfen die 8 erlaubten Restklassen im 210er Gitter für diesen Block
        allowed = [r for r in range(210) if np.gcd(r, 210) == 1]
        print(f"[*] Initialisiere Deep-Drill in Block k={k}...", flush=True)

        n_int = int(self.N)
        for r in allowed:
            p_cand = int(210 * k + r)
            if p_cand > 1 and n_int % p_cand == 0:
                q_cand = n_int // p_cand
                print("\n" + "=" * 60, flush=True)
                print("!!! CHRONOS-LOCK ERFOLGREICH !!!", flush=True)
                print(f"Faktor p: {p_cand}", flush=True)
                print(f"Faktor q: {q_cand}", flush=True)
                print("Verifikation: p * q == N", flush=True)
                print("=" * 60, flush=True)
                return p_cand
        return None

    def run_geodesic_jump(self, k_start):
        k = int(k_start)
        print(
            "[*] Chronos-Engine aktiv. Navigiere auf Geodäten im 256-Bit Raum…",
            flush=True,
        )

        for step in range(1000):
            v_phi, resonance = self.get_phase_velocity(k)

            jump = int(10000 / (1 + abs(v_phi * 10**6)))
            jump = max(1, jump)

            if step % 10 == 0:
                print(
                    f"k: {k:<20} | Resonanz: {resonance:.6f} | Jump: {jump}",
                    flush=True,
                )

            if abs(resonance) > 0.1:
                print(f"\n[!!!] SINGULARITÄT ERREICHT BEI k = {k}", flush=True)
                return self.final_lock(k)

            k += jump

        print("[*] Geodäte: 1000 Schritte ohne Schwellen-Treffer.", flush=True)
        return None


def _load_zeros():
    root = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(root, "zeros6.npy")
    if os.path.isfile(path):
        z = np.load(path).astype(float).ravel()
        return z[np.isfinite(z)]
    return np.linspace(14.0, 5000.0, 5000, dtype=np.float64)


if __name__ == "__main__":
    N_demo = 982_634_309
    zeros = _load_zeros()
    eng = ChronosEngine_256(N_demo, zeros, pilot_zeros=100_000)
    k0 = max(0, int(float(np.sqrt(N_demo)) // 210))
    print(f"[*] Start k0={k0}, N={N_demo}, γ-Pilot={eng._pilot_n}\n", flush=True)
    eng.run_geodesic_jump(k0)
