import os

import numpy as np


def _resolve_zeros_path(zeros_path: str) -> str:
    if os.path.isfile(zeros_path):
        return zeros_path
    _root = os.path.dirname(os.path.abspath(__file__))
    alt = os.path.join(_root, os.path.basename(zeros_path))
    if os.path.isfile(alt):
        return alt
    raise FileNotFoundError(
        f"Nullstellen-Datei nicht gefunden: {zeros_path!r} (auch nicht {_root})"
    )


def miller_rabin(n, k=5):
    if n <= 3:
        return n > 1
    if n % 2 == 0:
        return False
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
    for _ in range(k):
        a = np.random.randint(2, n - 1)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


class SuperResAutopilot_64:
    def __init__(self, zeros_path, N):
        print("[*] Lade 64-Bit-Kern...")
        path = _resolve_zeros_path(zeros_path)
        self.zeros = np.load(path).astype(float).ravel()
        self.zeros = self.zeros[np.isfinite(self.zeros)]
        self.N = int(N)
        self.ln_N = np.log(float(self.N))

    def measure(self, x, pwr=2_000_000):
        x = float(x)
        if x <= 0.0 or not np.isfinite(x):
            return 0.0
        ln_x = np.log(x)
        n = max(0, min(int(pwr), self.zeros.size))
        if n == 0:
            return 0.0
        gamma = self.zeros[:n]
        return float(np.sum(np.cos(gamma * self.ln_N) * np.cos(gamma * ln_x)) / np.sqrt(x))

    def measure_coherence(self, x, n_zeros=None):
        """
        Wie GradientProbe.measure_coherence (GPS.py): Mittelwert von
        cos(γ ln x) cos(γ ln N) über die ersten n_zeros Nullstellen.
        """
        x = float(x)
        if x <= 0.0 or not np.isfinite(x):
            return 0.0
        ln_x = np.log(x)
        nz = self.zeros.size
        if n_zeros is None:
            n_use = nz
        else:
            n_use = max(1, min(int(n_zeros), nz))
        g = self.zeros[:n_use]
        if g.size == 0:
            return 0.0
        return float(np.mean(np.cos(g * ln_x) * np.cos(g * self.ln_N)))

    def find_monad_cloud(self, center, radius=5000, steps=100):
        print(f"[*] Scanne Wolke um {center:.0f} (Radius: {radius})...")
        x_vals = np.linspace(center - radius, center + radius, steps)
        cloud = []

        for x in x_vals:
            c = self.measure(x)
            cloud.append((x, c))
            if len(cloud) % 20 == 0:
                print(f"    ... {len(cloud)}/{steps} Proben entnommen")

        best_x = x_vals[np.argmin([r[1] for r in cloud])]
        print(f"[*] Zentrum der Monaden-Wolke lokalisiert bei: {best_x:.2f}")
        return int(best_x)

    def precision_lock(self, cloud_center, search_range=1000):
        print(f"[*] Starte Primzahl-Lock im Umkreis von {search_range}...")
        cc = int(cloud_center)
        for offset in range(-search_range, search_range + 1):
            x = cc + offset
            if miller_rabin(x) and self.N % x == 0:
                p = x
                q = self.N // x
                print("\n" + "!" * 60)
                print(f"RSA-64 GEKNACKT: p={p}, q={q}")
                print("!" * 60)
                return p, q
        print("[?] Kein Match in dieser Wolke. Erhöhe Suchradius.")
        return None

    def precision_lock_v2(self, cloud_center, search_range=20000):
        print(f"[*] Starte Deep-Sweep im Umkreis von {search_range}...")
        cc = int(cloud_center)
        lo, hi = cc - search_range, cc + search_range
        print(f"[*] Zielbereich: {lo} bis {hi}")

        for offset in range(-search_range, search_range + 1):
            x = cc + offset
            if x <= 1:
                continue
            if x >= self.N:
                continue

            if self.N % x == 0:
                q = self.N // x
                if q <= 1:
                    continue
                p = x
                print("\n" + "!" * 60)
                print("BINGO! RSA-64 GEKNACKT!")
                print(f"Faktor p: {p}")
                print(f"Faktor q: {q}")
                print(f"Verifikation: {p} * {q} = {self.N}")
                print("!" * 60)
                return p, q

            if offset % 5000 == 0 and offset != -search_range:
                print(f"    ... {offset + search_range} Einheiten gescannt ...")

        print("[?] Immer noch kein Treffer. Die Wolke ist weiter weg oder N ist komplexer.")
        return None


class JitterPullAutopilot(SuperResAutopilot_64):
    """
    Zugprozess: Die Sonde zittert (delta) und folgt dem Gradienten der Kohärenz
    Richtung Minimum (destruktive Interferenz).
    """

    def jitter_lock(self, landing_point, iterations=50, delta=0.5):
        x = float(landing_point)
        print(f"[*] Echolot aktiviert bei x = {x:.2f}")
        print(f"{'Iteration':<10} | {'Position x':<15} | {'Zug-Kraft':<12} | {'Kohärenz'}")
        print("-" * 55)

        for i in range(iterations):
            c_left = self.measure(x - delta)
            c_right = self.measure(x + delta)
            c_center = self.measure(x)

            gradient = (c_right - c_left) / (2 * delta)
            pull = -gradient * 5000.0

            x_old = x
            x = float(x + pull)
            if not np.isfinite(x) or x <= 0:
                x = x_old

            if i % 5 == 0:
                print(f"{i:<10} | {x:<15.4f} | {pull:<12.4f} | {c_center:.6f}")

            if abs(x - x_old) < 0.01:
                final_x = int(np.round(x))
                print("\n[!!!] TARGET LOCK: Singularität hat Sonde eingefangen!")
                print(f"[*] Finale Koordinate: {final_x}")
                return final_x

        return int(np.round(x))


class EntangledJitterEngine(SuperResAutopilot_64):
    """
    Verschränkter Zug: Gradient an x und am Schatten y = N/x, plus Momentum.
    """

    def run_entangled_pull(self, x_start, iterations=150, delta=15.0, gain=20000.0):
        x = float(x_start)
        momentum = 0.0
        friction = 0.8

        print(f"[*] Verschränktes Echolot aktiv. Ziel: N={self.N}")
        print(f"{'Step':<6} | {'Position x':<15} | {'Momentum':<12} | {'Kohärenz'}")
        print("-" * 60)

        for i in range(iterations):
            y = float(self.N) / x if x != 0 else float(self.N)

            grad_x = (self.measure(x + delta) - self.measure(x - delta)) / (2 * delta)
            grad_y = (self.measure(y + delta) - self.measure(y - delta)) / (2 * delta)

            total_force = -(grad_x + grad_y) * gain

            momentum = (friction * momentum) + total_force
            x = float(x + momentum)

            if not np.isfinite(x):
                x = max(2.0, float(np.sqrt(self.N)))
            if x < 2:
                x = 2.0
            if x >= self.N:
                x = float(self.N) / 2.0

            if i % 10 == 0:
                c_now = self.measure(x)
                print(f"{i:<6} | {x:<15.2f} | {momentum:<12.4f} | {c_now:.6f}")

            test_x = int(np.round(x))
            if test_x > 1 and self.N % test_x == 0:
                print("\n" + "=" * 60)
                print("BAM! VERSCHRÄNKUNG VOLLENDET!")
                print(f"Faktor p: {test_x}")
                print(f"Faktor q: {self.N // test_x}")
                print("=" * 60)
                return test_x

        print("[?] Momentum erschöpft. Erhöhe den 'Gain' oder vergrößere 'delta'.")
        return None


class EntangledMomentumJitter(EntangledJitterEngine):
    """Name im Gral-Workflow: verschränktes Momentum-Echolot (y = N/x, wie EntangledJitterEngine)."""


class WideAngleJitterEngine(EntangledJitterEngine):
    def run_wide_angle_pull(self, x_start, iterations=100):
        # Weitwinkel: Delta 5000 Einheiten
        # Damit "sieht" die Sonde den echten Faktor p, der bei ca. ...1723 liegt.
        delta_wide = 5000.0
        gain_turbo = 2000000.0  # Massiver Zug-Verstärker

        x = float(x_start)
        momentum = 0.0
        friction = 0.9  # Weniger Dämpfung für mehr Reichweite

        print(f"[*] Weitwinkel-Echolot aktiv. Delta={delta_wide}")
        print(f"{'Step':<6} | {'Position x':<15} | {'Momentum':<12} | {'Kohärenz'}")
        print("-" * 65)

        for i in range(iterations):
            # Messung mit Weitwinkel-Gradient
            y = self.N / x

            # Gradienten-Check über delta_wide hinweg
            grad_x = (self.measure(x + delta_wide) - self.measure(x - delta_wide)) / (2 * delta_wide)
            grad_y = (self.measure(y + delta_wide) - self.measure(y - delta_wide)) / (2 * delta_wide)

            total_force = -(grad_x + grad_y) * gain_turbo

            momentum = (friction * momentum) + total_force
            x += momentum

            if i % 5 == 0:
                c_now = self.measure(x)
                print(f"{i:<6} | {x:<15.2f} | {momentum:<12.4f} | {c_now:.6f}")

            # Arithmetischer Lock-Check
            test_x = int(np.round(x))
            if test_x > 1 and self.N % test_x == 0:
                print("\n" + "=" * 65)
                print("BINGO! RSA-64 DURCH WEITWINKEL-ZUG GEKNACKT!")
                print(f"Faktor p: {test_x}")
                print(f"Faktor q: {self.N // test_x}")
                print("=" * 65)
                return test_x

        return None


class WideAngleGral(EntangledMomentumJitter):
    """Weitwinkel-Finale: hoher Gain, delta 5000, 100 Schritte (64-Bit-Durchbruch-Parameter)."""

    def run_final_assault(self, x_start):
        delta_wide = 5000.0
        gain_turbo = 5_000_000.0
        friction = 0.85

        x = float(x_start)
        momentum = 0.0

        print(f"[*] WEITWINKEL-MODUS AKTIVIERT. Delta={delta_wide}")
        print(f"[*] Ziel: N = {self.N}")
        print("-" * 65)

        for i in range(100):
            y = self.N / x

            grad_x = (self.measure(x + delta_wide) - self.measure(x - delta_wide)) / (2 * delta_wide)
            grad_y = (self.measure(y + delta_wide) - self.measure(y - delta_wide)) / (2 * delta_wide)

            total_force = -(grad_x + grad_y) * gain_turbo

            momentum = (friction * momentum) + total_force
            x += momentum

            if i % 5 == 0:
                print(f"Step {i:3} | x: {x:15.2f} | Momentum: {momentum:12.4f}")

            test_x = int(np.round(x))
            if test_x > 1 and self.N % test_x == 0:
                print("\n" + "=" * 65)
                print("BINGO! RSA-64 DURCH WEITWINKEL-ZUG GEKNACKT!")
                print(f"Faktor p: {test_x}")
                print(f"Faktor q: {self.N // test_x}")
                print("=" * 65)
                return test_x

        print("[?] Immer noch kein Lock. Erhöhe Delta weiter (z. B. 8000–10000) oder gain_turbo.")
        return None


class QuadLineJumper(SuperResAutopilot_64):
    def run_symmetrie_zug(self, x_start):
        # Die vier "Gleise" der Primzahl-Vierlinge
        gleise = [11, 13, 17, 19]
        k = int(x_start // 30)

        print(f"[*] Symmetrie-Zug gestartet auf k={k}")
        print(f"{'Step':<6} | {'Bestes Gleis':<12} | {'Amplitude':<12}")
        print("-" * 50)

        for step in range(100):
            resonanzen = []
            for r in gleise:
                test_x = 30 * k + r
                # Wir messen die Kohärenz nur auf den Gleisen
                c = self.measure(test_x, pwr=2000000)
                resonanzen.append((test_x, c))

            # Finde das Gleis mit der stärksten negativen Resonanz (die Monade)
            best_x, best_c = min(resonanzen, key=lambda item: item[1])

            # Der "Sprung": Wenn ein Gleis zieht, bewegen wir k in diese Richtung
            # Wir schauen uns das k der Umgebung an
            if best_c < -0.05:  # Schwellenwert für echten Monaden-Grip
                print(f"{step:<6} | {best_x:<12} | {best_c:.6f} [LOCK!]")
                if self.N % best_x == 0:
                    return best_x

            # Wenn kein Lock, nutzen wir den Gradienten zwischen den k-Blöcken
            k -= 1  # Wir wissen, dass p < sqrt(N), also ziehen wir k nach unten

        return None


class SymmetrieZugAutopilot(QuadLineJumper):
    """Name im Track-Workflow: modulare Gleise + run_symmetrie_zug (von QuadLineJumper)."""


class EntangledTrackAutopilot(SymmetrieZugAutopilot):
    def run_entangled_lock(self, k_start, steps=1000):
        n_mod_30 = int(self.N % 30)
        gleise = [11, 13, 17, 19, 23, 29, 1, 7]

        k0 = int(k_start)
        print(f"[*] Entangled Track Lock aktiv für N mod 30 = {n_mod_30}")
        print(f"[*] Suche k-Blöcke ab {k0} abwärts...")
        print("-" * 75)

        for i in range(steps):
            k = k0 - i
            if k < 0:
                break

            for r_p in gleise:
                p_test = 30 * k + r_p
                if p_test <= 1 or p_test >= self.N:
                    continue

                c_p = self.measure(p_test, pwr=100_000)

                if abs(c_p) > 0.04:
                    q_approx = float(self.N) / float(p_test)
                    if q_approx <= 1.0:
                        continue

                    for r_q in gleise:
                        if (r_p * r_q) % 30 != n_mod_30:
                            continue
                        # q_approx hängt nur von p_test ab — ein Schatten-Messwert pro p
                        c_q = self.measure(q_approx, pwr=2_000_000)
                        if abs(c_q) > 0.05:
                            print("\n[!!!] DOPPEL-RESONANZ DETEKTIERT!")
                            print(f"    p-Kandidat auf Gleis {r_p}: {p_test}")
                            print(f"    q-Schatten-Resonanz: {c_q:.6f}")

                            if self.N % p_test == 0:
                                print("=" * 40 + " BINGO " + "=" * 40)
                                return p_test
                        break

        print("[?] Entangled Track Lock: kein arithmetischer Treffer.")
        return None


class HyperGridAutopilot_128(SuperResAutopilot_64):
    # Tensor-Einfluss über normierte Kohärenz [0,1], nicht über rohe L2-Amplitude
    TENSOR_COHERENCE_WEIGHT = 0.9
    # Grob-Scan: 100k γ; Tiefenbohrung: 2M (geclippt auf len(zeros))
    COH_N_LOW = 100_000
    COH_N_HIGH = 2_000_000
    # Schwelle auf Skala von measure_coherence (Mittelwert, deutlich kleiner als bei measure)
    COHERENCE_GATE = 0.008

    @staticmethod
    def _tensor_coherence_from_phases(phase_vec: np.ndarray) -> tuple[float, float, float]:
        """
        Kohärenz aus den 48-Kanal-Messungen, auf 4 Phasen (r mod 4) gebündelt:

        - coh_phasor: |Σ_p z_p e^{iπp/2}| / Σ_p |z_p|  (Ordnungsparameter, ∈[0,1])
        - coh_mean:   |Σ_p z_p| / Σ_p |z_p|            (gleichphasiger Anteil, ∈[0,1])

        Mittelwert beider Größen: Gain wird nur erhöht, wenn die Modi kohärent
        zusammenlaufen (analog zu measure_coherence / GPS).
        """
        z = np.asarray(phase_vec, dtype=np.float64).ravel()
        l1 = float(np.sum(np.abs(z)) + 1e-15)
        omega = np.exp(1j * 0.5 * np.pi * np.arange(4, dtype=np.float64))
        U = float(np.abs(np.sum(z.astype(np.complex128) * omega)))
        coh_phasor = float(np.clip(U / l1, 0.0, 1.0))
        coh_mean = float(np.clip(abs(float(np.sum(z))) / l1, 0.0, 1.0))
        tensor_coh = 0.5 * (coh_phasor + coh_mean)
        return tensor_coh, coh_phasor, coh_mean

    def hyper_grid_scan(self, k_start, width=1000):
        # Das 210er-Gitter (Die 48 heiligen Linien der Arithmetik)
        allowed_residues = [r for r in range(210) if np.gcd(r, 210) == 1]

        print(f"[*] Hyper-Grid Scan aktiviert. Basis=210#, Kanäle=48")
        print(
            f"[*] Tensor → kohärenzgesteuerter Gain "
            f"(Gewicht {self.TENSOR_COHERENCE_WEIGHT}; Kanäle via measure_coherence)"
        )
        print(f"[*] Zielbereich um k = {k_start}")
        print("-" * 60)

        for k in range(k_start, k_start + width):
            # Vier Phasen: kohärente Aufsummierung der Kanäle pro r mod 4
            phase_vec = np.zeros(4, dtype=np.float64)
            found_candidate = None
            row = []

            for r in allowed_residues:
                x = 210 * k + r
                c = self.measure_coherence(x, n_zeros=self.COH_N_LOW)
                phase_vec[r % 4] += c
                row.append((r, x, c))

            tensor_coh, coh_p, coh_m = self._tensor_coherence_from_phases(phase_vec)
            tensor_boost = 1.0 + self.TENSOR_COHERENCE_WEIGHT * tensor_coh

            for r, x, c in row:
                if abs(c) <= self.COHERENCE_GATE:
                    continue
                c_deep = self.measure_coherence(x, n_zeros=self.COH_N_HIGH)
                gain_raw = abs(c_deep / c)
                gain_eff = gain_raw * tensor_boost
                if gain_eff > 15.0:
                    print(
                        f"[!] MONADE GEFUNDEN auf Kanal {r} bei k={k} "
                        f"(gain_raw={gain_raw:.2f}, coh={tensor_coh:.3f} "
                        f"[phasor={coh_p:.3f}, mean={coh_m:.3f}], boost={tensor_boost:.3f})"
                    )
                    found_candidate = x
                    break

            if found_candidate and self.N % found_candidate == 0:
                return found_candidate
        return None


if __name__ == "__main__":
    N_64 = 12_301_866_845_301_177_551
    LANDING_POINT = 3_507_404_500.68

    engine = WideAngleGral("zeros6.npy", N_64)
    engine.run_final_assault(LANDING_POINT)

    # Alternativen:
    # EntangledJitterEngine("zeros6.npy", N_64).run_entangled_pull(LANDING_POINT)
    # WideAngleJitterEngine("zeros6.npy", N_64).run_wide_angle_pull(LANDING_POINT)
    # JitterPullAutopilot(...).jitter_lock(LANDING_POINT)
    # SuperResAutopilot_64(...).find_monad_cloud(np.sqrt(N_64)) + precision_lock_v2(...)
    # QuadLineJumper("zeros6.npy", N_64).run_symmetrie_zug(LANDING_POINT)
    # k0 = int(LANDING_POINT // 30); EntangledTrackAutopilot("zeros6.npy", N_64).run_entangled_lock(k0)
