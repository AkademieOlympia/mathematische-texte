"""
Autopilot-Kette: GradientProbe / TopologischerAutopilot (GPS.py)
-> OctonionicAutopilot -> TargetLockAutopilot -> RadarGuidedAutopilot
-> DampedRadarAutopilot (Radar + Richtungs-Dämpfung)
-> BreakoutAutopilot (dynamisches Radar + Kick bei Oszillation)
-> QuantumSwarmAutopilot (Mehrsonden-Schwarm + „Tunnel“ zum Besten)
-> PolarizedSwarmAutopilot (Tabu-Zonen / „Arithmetic Ghosts“)
-> MonadicZoomAutopilot (Kohärenz vs. wachsende Nullstellen-Anzahl).
"""
from __future__ import annotations

import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.abspath(__file__))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from GPS import TopologischerAutopilot  # noqa: E402


class OctonionicAutopilot(TopologischerAutopilot):
    def __init__(self, zeros, N):
        super().__init__(zeros, N)

    def berechne_oktonionischen_boost(self, x):
        """
        Misst die Nicht-Assoziativität am Punkt x.
        Ein hoher Wert bedeutet 'instabiles Gelände' -> Schneller fliegen!
        """
        x = float(x)
        if x <= 0.0 or not np.isfinite(x):
            return 1.0

        n = len(self.zeros)
        n = n - (n % 8)
        if n < 8:
            return 1.0

        parts = np.array_split(self.zeros[:n], 8)
        ln_x = np.log(x)
        oct_state = [float(np.sum(np.cos(sec * ln_x))) for sec in parts]
        associator_tension = float(np.std(oct_state))
        return 1.0 + associator_tension * 10.0

    def warp_flug(self, start_x, steps=15):
        x = float(start_x)
        pfad = [x]

        print(f"[*] Warp-Antrieb aktiviert bei x={x:.2f}")
        for i in range(steps):
            eps = 0.1
            g1 = self.berechne_topologische_guete(x - eps)
            g2 = self.berechne_topologische_guete(x + eps)
            topo_grad = (g2 - g1) / (2 * eps)

            boost = self.berechne_oktonionischen_boost(x)
            schrittweite = 50.0 * boost
            x = x - schrittweite * topo_grad
            pfad.append(x)

            if i % 3 == 0:
                print(f"    Sonden-Log {i}: Pos={x:.2f} | Boost={boost:.2f}x")

        return pfad


class TargetLockAutopilot(OctonionicAutopilot):
    def __init__(self, zeros, N, target_threshold=0.85):
        super().__init__(zeros, N)
        self.lock_threshold = target_threshold  # Schwelle für H2-Stabilität
        self.locked = False

    def navigations_logik(self, x):
        """Entscheidet zwischen Warp-Flug und Ziel-Anflug."""
        guete = self.berechne_topologische_guete(x)

        if guete > self.lock_threshold:
            if not self.locked:
                print(f"!!! TARGET LOCK bei x={x:.2f} (Guete: {guete:.4f}) !!!")
                self.locked = True

            boost = 1.0
            lernrate = 5.0  # Sanfte Landung
        else:
            boost = self.berechne_oktonionischen_boost(x)
            lernrate = 50.0  # Aggressives Überfliegen der Maische

        return boost, lernrate, guete

    def navigations_logik_v2(self, x, letzte_schritte):
        """
        Wie navigations_logik, aber mit Oszillations-Dämpfung (Zick-Zack in letzte_schritte)
        und hart begrenztem Boost im Warp-Modus.
        """
        guete = self.berechne_topologische_guete(x)

        if (
            len(letzte_schritte) > 2
            and np.sign(letzte_schritte[-1]) != np.sign(letzte_schritte[-2])
        ):
            damping = 0.1
        else:
            damping = 1.0

        if guete > self.lock_threshold:
            if not self.locked:
                print(f"!!! TARGET LOCK bei x={x:.2f} (Guete: {guete:.4f}) !!!")
            self.locked = True
            boost = 1.0
            lr = 1.0
        else:
            boost = self.berechne_oktonionischen_boost(x)
            boost = min(boost, 50.0)
            lr = 5.0 * damping

        return boost, lr, guete

    def präzisions_flug(self, start_x, steps=20, navigation_version="v1"):
        x = float(start_x)
        pfad = [x]
        self.locked = False
        letzte_schritte = []

        print(f"[*] Sonde gestartet. Suche stabile H2-Signatur für N={self.N}...")
        for i in range(steps):
            if navigation_version == "v2":
                boost, lr, current_g = self.navigations_logik_v2(x, letzte_schritte)
            else:
                boost, lr, current_g = self.navigations_logik(x)

            eps = 0.1
            g1 = self.berechne_topologische_guete(x - eps)
            g2 = self.berechne_topologische_guete(x + eps)
            topo_grad = (g2 - g1) / (2 * eps)

            schritt = (lr * boost) * topo_grad
            prev_x = x
            x = x - schritt
            letzte_schritte.append(x - prev_x)
            if len(letzte_schritte) > 64:
                letzte_schritte = letzte_schritte[-64:]
            pfad.append(x)

            if i % 2 == 0:
                mode = "LOCK" if self.locked else "WARP"
                print(f"    [{mode}] Schritt {i}: x={x:.2f} | Guete={current_g:.4f} | Boost={boost:.2f}x")

            if self.locked and abs(schritt) < 1e-5:
                print(f"[*] Landung erfolgreich bei x={x:.6f}")
                break

        return pfad


class RadarGuidedAutopilot(TargetLockAutopilot):
    def __init__(self, zeros, N, radar_range=10.0, target_threshold=0.85):
        super().__init__(zeros, N, target_threshold=target_threshold)
        self.radar_range = float(radar_range)

    def surround_scan(self, x):
        """Führt eine 360-Grad-Peilung im arithmetischen Raum durch."""
        x = float(x)
        p_links = x - self.radar_range
        p_rechts = x + self.radar_range

        g_links = self.berechne_topologische_guete(p_links)
        g_hier = self.berechne_topologische_guete(x)
        g_rechts = self.berechne_topologische_guete(p_rechts)

        if g_links > g_rechts and g_links > g_hier:
            return -1.0
        if g_rechts > g_links and g_rechts > g_hier:
            return 1.0

        s = float(np.sign(g_rechts - g_links))
        if s == 0.0:
            eps = 0.1
            g1 = self.berechne_topologische_guete(x - eps)
            g2 = self.berechne_topologische_guete(x + eps)
            s = float(np.sign(g2 - g1))
            if s == 0.0:
                s = 1.0
        return s

    def flug_mit_radar(self, start_x, steps=20, navigation_version="v1"):
        x = float(start_x)
        pfad = [x]
        self.locked = False
        letzte_schritte = []

        print(f"[*] Radar-Scan aktiviert bei x={x:.2f}...")
        for i in range(steps):
            richtung = self.surround_scan(x)
            if navigation_version == "v2":
                boost, lr, current_g = self.navigations_logik_v2(x, letzte_schritte)
            else:
                boost, lr, current_g = self.navigations_logik(x)

            schrittweite = richtung * (lr * boost) * 0.05
            x = x + schrittweite
            letzte_schritte.append(schrittweite)
            if len(letzte_schritte) > 64:
                letzte_schritte = letzte_schritte[-64:]
            pfad.append(x)

            mode = "LOCK" if self.locked else "WARP"
            pfeil = "->" if richtung > 0 else "<-" if richtung < 0 else "o"
            print(
                f"    [{mode}] Radar-Fix: x={x:.2f} | Richtung={pfeil} | Guete={current_g:.4f}"
            )

            if self.locked and abs(schrittweite) < 1e-4:
                print("[*] Punktlandung erfolgreich!")
                break

        return pfad


class MonadicZoomAutopilot(RadarGuidedAutopilot):
    """
    Spektrale „Tiefenbohrung“: Kohärenz an festem x bei wachsender Anzahl
    genutzter Nullstellen (Hardware-Power), um stabile Resonanz (Monade) von Rauschen zu trennen.
    """

    def __init__(self, zeros, N):
        super().__init__(zeros, N)
        self.max_hardware_power = len(self.zeros)

    def spektrale_tiefenbohrung(self, x_position, schritte=10):
        x_position = float(x_position)
        print(f"\n[*] Starte Tiefenbohrung an Position x = {x_position:.4f}")
        print(f"    {'Hardware-Power':<20} | {'Kohärenz':<15} | {'Stabilität'}")
        print("-" * 55)

        verlauf: list[float] = []
        pmax = self.max_hardware_power
        if pmax <= 1:
            power_stufen = np.array([1], dtype=int)
        elif pmax < 100_000:
            power_stufen = np.linspace(1, pmax, max(1, schritte)).astype(int)
        else:
            power_stufen = np.linspace(100_000, pmax, max(1, schritte)).astype(int)
        power_stufen = np.clip(np.unique(power_stufen), 1, pmax)

        for pwr in power_stufen:
            k = self.measure_coherence(x_position, n_zeros=int(pwr))
            verlauf.append(k)
            stabilität = float(np.std(verlauf[-3:])) if len(verlauf) >= 3 else 0.0
            print(f"    {int(pwr):<20,} | {k:<15.6f} | {stabilität:.6f}")

        if not verlauf:
            print("[?] Keine Messstufen — Abbruch.")
            return False, 0.0

        v0, v1 = verlauf[0], verlauf[-1]
        if abs(v0) < 1e-15:
            tendenz = 0.0
        else:
            tendenz = v1 / v0

        if tendenz > 1.5 and abs(v1) > 0.05:
            print(f"[!!!] MONADE BESTÄTIGT: Resonanz verstärkt sich um Faktor {tendenz:.2f}")
            return True, v1

        print("[?] GEISTER-PEAK: Resonanz ist instabil oder zu schwach.")
        return False, v1

    def scan_and_zoom(self, kandidaten_liste):
        ergebnisse = []
        for x in kandidaten_liste:
            is_monad, final_power = self.spektrale_tiefenbohrung(float(x))
            if is_monad:
                ergebnisse.append((x, final_power))
        return ergebnisse


def monadic_deep_drill(hw, N, test_points, n_low=100_000, n_high=2_000_000):
    """
    Beschießt Punkte x mit zwei Spektralauflösungen (Standard: 100k vs. 2M γ).

    `hw`: Objekt mit Attribut ``all_zeros`` (z. B. ``RiemannHardware`` aus Holonomy Okto)
    oder direkt ein 1D-``numpy``-Array von Nullstellen.
    """
    z = np.asarray(getattr(hw, "all_zeros", hw), dtype=float).ravel()
    n_z = z.size
    if n_z == 0:
        raise ValueError("Keine Nullstellen (all_zeros leer).")

    n_a = max(0, min(int(n_low), n_z))
    n_b = max(0, min(int(n_high), n_z))
    if n_a == 0:
        raise ValueError("n_low nach Clipping 0 — zu wenige Nullstellen.")

    ln_N = float(np.log(float(N)))
    results: dict[float, tuple[float, float]] = {}

    print(f"{'Punkt x':<15} | {f'{n_a:,} Zeros':<12} | {f'{n_b:,} Zeros':<12} | {'Trend'}")
    print("-" * 60)

    for x in test_points:
        x = float(x)
        if x <= 0.0 or not np.isfinite(x):
            continue
        ln_x = np.log(x)
        ga, gb = z[:n_a], z[:n_b]
        k_low = float(
            np.sum(np.cos(ga * ln_N) * np.cos(ga * ln_x)) / np.sqrt(x)
        )
        k_high = float(
            np.sum(np.cos(gb * ln_N) * np.cos(gb * ln_x)) / np.sqrt(x)
        )

        trend = (
            "AUFSTEIGEND (MONADE!)"
            if abs(k_high) > abs(k_low) * 1.5
            else "ZERFALL (MAISCHE)"
        )
        print(f"{x:<15.2f} | {k_low:<12.4f} | {k_high:<12.4f} | {trend}")
        results[x] = (k_low, k_high)

    return results


class QuantumSwarmAutopilot(RadarGuidedAutopilot):
    """Viele Sonden um ein Zentrum; lokaler Radar-Schritt, dann Zug zum besten Signal (0.7/0.3-Mix)."""

    def __init__(self, zeros, N, n_probes=10):
        super().__init__(zeros, N)
        self.n_probes = int(n_probes)
        self.positions: np.ndarray = np.array([], dtype=float)

    def schwarm_flug(self, center_x, radius=100.0, steps=15):
        center_x = float(center_x)
        if self.n_probes < 1:
            raise ValueError("n_probes muss mindestens 1 sein")

        self.positions = np.linspace(
            center_x - radius, center_x + radius, self.n_probes, dtype=float
        )
        print(f"[*] Schwarm-Start: {self.n_probes} Monaden ausgesetzt um x={center_x}")

        if steps <= 0:
            gu0 = np.array(
                [self.berechne_topologische_guete(float(p)) for p in self.positions],
                dtype=float,
            )
            k = int(np.argmax(gu0))
            return float(self.positions[k])

        guete_werte: list[float] = []
        for i in range(steps):
            guete_werte = []
            for j in range(self.n_probes):
                pj = float(self.positions[j])
                g = self.berechne_topologische_guete(pj)
                guete_werte.append(g)

                richtung = self.surround_scan(pj)
                self.positions[j] = pj + richtung * 5.0

            best_idx = int(np.argmax(guete_werte))
            best_x = float(self.positions[best_idx])
            max_g = guete_werte[best_idx]

            for j in range(self.n_probes):
                self.positions[j] = 0.7 * self.positions[j] + 0.3 * best_x

            print(f"    Schritt {i}: Bestes Signal bei x={best_x:.2f} | Güte={max_g:.4f}")

            if max_g > self.lock_threshold:
                print(f"[!!!] TARGET LOCK DURCH SCHWARM-KOHÄRENZ bei x={best_x:.4f}")
                return best_x

        return float(self.positions[int(np.argmax(guete_werte))])


class PolarizedSwarmAutopilot(QuantumSwarmAutopilot):
    """Schwarm wie QuantumSwarm, aber Güte wird nahe konfigurierbarer Tabu-Zonen abgeschwächt."""

    def __init__(self, zeros, N):
        super().__init__(zeros, N)
        self.tabu_zones = [15360.0, 15250.0]

    def berechne_topologische_guete(self, x):
        x = float(x)
        g = float(super().berechne_topologische_guete(x))

        for zone in self.tabu_zones:
            dist = abs(x - float(zone))
            if dist < 50.0:
                g *= dist / 50.0

        return g


class DampedRadarAutopilot(RadarGuidedAutopilot):
    """
    Radar geführter Flug mit Anti-Oszillation: Dämpfung bei Richtungswechsel
    (surround_scan). Die Methode navigations_logik_v2(x, current_direction) ist
    bewusst anders signiert als TargetLockAutopilot.navigations_logik_v2(x, letzte_schritte);
    daher sind präzisions_flug(..., v2) und flug_mit_radar(..., v2) hier deaktiviert —
    stattdessen präzisions_anflug() nutzen.
    """

    def __init__(
        self,
        zeros,
        N,
        max_step=20.0,
        radar_range=10.0,
        target_threshold=0.85,
    ):
        super().__init__(zeros, N, radar_range=radar_range, target_threshold=target_threshold)
        self.max_step = float(max_step)
        self.last_direction = 0.0

    def navigations_logik_v2(self, x, current_direction):
        guete = self.berechne_topologische_guete(x)

        damping = 1.0
        if current_direction != self.last_direction and self.last_direction != 0:
            damping = 0.05
            print("    [DAMPING] Richtungswechsel erkannt! Energie drosseln.")

        if guete > self.lock_threshold:
            self.locked = True
            boost = 1.0
            lr = 1.0
        else:
            boost = self.berechne_oktonionischen_boost(x)
            boost = min(boost, 10.0)
            lr = self.max_step * damping

        self.last_direction = float(current_direction)
        return boost, lr, guete

    def präzisions_flug(self, start_x, steps=20, navigation_version="v1"):
        if navigation_version == "v2":
            raise ValueError(
                "DampedRadarAutopilot: v2 erwartet hier current_direction, nicht letzte_schritte — "
                "nutze präzisions_anflug() oder einen RadarGuidedAutopilot für listenbasiertes v2."
            )
        return super().präzisions_flug(start_x, steps=steps, navigation_version=navigation_version)

    def flug_mit_radar(self, start_x, steps=20, navigation_version="v1"):
        if navigation_version == "v2":
            raise ValueError(
                "DampedRadarAutopilot: v2 hier nicht unterstützt — nutze präzisions_anflug()."
            )
        return super().flug_mit_radar(start_x, steps=steps, navigation_version=navigation_version)

    def präzisions_anflug(self, start_x, steps=30):
        x = float(start_x)
        pfad = [x]
        self.locked = False
        self.last_direction = 0.0

        print(f"[*] Dämpfungs-Autopilot aktiv. Ziel-Suche für N={self.N}...")

        for i in range(steps):
            richtung = self.surround_scan(x)
            boost, lr, current_g = self.navigations_logik_v2(x, richtung)

            schritt = richtung * lr * (boost / 10.0)
            x = x + schritt
            pfad.append(x)

            mode = "LOCK" if self.locked else "WARP"
            if i % 2 == 0:
                print(
                    f"    [{mode}] Schritt {i}: x={x:.2f} | Güte={current_g:.4f} | Schrittw={schritt:.2f}"
                )

            if self.locked and abs(schritt) < 1e-4:
                print(f"[*] Punktlandung bei x={x:.6f}!")
                break

        return pfad


class BreakoutAutopilot(DampedRadarAutopilot):
    """Erweitert DampedRadar: bei Richtungswechsel Radar-Reichweite aufblasen, nach vielen Flips Kick."""

    def __init__(self, zeros, N):
        super().__init__(zeros, N)
        self.oscillation_count = 0
        self.dynamic_radar = 10.0

    def navigations_logik_v3(self, x, richtung):
        guete = self.berechne_topologische_guete(x)

        if richtung != self.last_direction:
            self.oscillation_count += 1
            self.dynamic_radar = min(self.dynamic_radar * 1.5, 200.0)
        else:
            self.oscillation_count = max(0, self.oscillation_count - 1)
            self.dynamic_radar = 10.0

        kick = 50.0 if self.oscillation_count > 4 else 1.0
        boost, lr, g = self.navigations_logik_v2(x, richtung)
        return boost * kick, lr, g

    def flug_mit_ausbruch(self, start_x, steps=30):
        x = float(start_x)
        pfad = [x]
        self.locked = False
        self.last_direction = 0.0
        self.oscillation_count = 0
        self.dynamic_radar = 10.0

        print(f"[*] Breakout-Modus: dynamisches Radar + Kick für N={self.N}...")

        for i in range(steps):
            saved_rr = self.radar_range
            self.radar_range = float(self.dynamic_radar)
            try:
                richtung = RadarGuidedAutopilot.surround_scan(self, x)
            finally:
                self.radar_range = saved_rr

            boost, lr, current_g = self.navigations_logik_v3(x, richtung)

            schritt = richtung * lr * (boost / 10.0)
            x = x + schritt
            pfad.append(x)

            mode = "LOCK" if self.locked else "WARP"
            if i % 2 == 0:
                print(
                    f"    [{mode}] Schritt {i}: x={x:.2f} | Güte={current_g:.4f} | "
                    f"Schrittw={schritt:.2f} | Radar={self.dynamic_radar:.1f} | osc={self.oscillation_count}"
                )

            if self.locked and abs(schritt) < 1e-4:
                print(f"[*] Punktlandung bei x={x:.6f}!")
                break

        return pfad


def _load_zeros():
    path = os.path.join(_ROOT, "zeros6.npy")
    if os.path.isfile(path):
        z = np.load(path).astype(float).ravel()
        z = z[np.isfinite(z)]
        print(f"[Autopilot] zeros6.npy: {z.size:,} Nullstellen geladen (volle Liste für Güte/Kohärenz).", file=sys.stderr)
        return z
    print("[Autopilot] zeros6.npy fehlt — nutze 40 feste γ.", file=sys.stderr)
    return np.array(
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
        ],
        dtype=float,
    )


if __name__ == "__main__":
    zeros = _load_zeros()
    N_demo = 243_382_003
    print(f"[*] Arbeitsspeicher-Nullstellen: {zeros.size:,} (alle werden von GradientProbe/Autopilot genutzt)\n")

    print("--- OctonionicAutopilot.warp_flug ---")
    autopilot = OctonionicAutopilot(zeros, N_demo)
    pfad = autopilot.warp_flug(start_x=15_300.0)
    print(f"[*] Pfadlänge {len(pfad)} | Ende x={pfad[-1]:.4f}\n")

    print("--- TargetLockAutopilot.präzisions_flug ---")
    target = TargetLockAutopilot(zeros, N_demo, target_threshold=0.85)
    pfad_lock = target.präzisions_flug(start_x=15_300.0, steps=25)
    print(f"[*] Pfadlänge {len(pfad_lock)} | Ende x={pfad_lock[-1]:.4f} | locked={target.locked}\n")

    print("--- RadarGuidedAutopilot.flug_mit_radar ---")
    radar = RadarGuidedAutopilot(zeros, N_demo, radar_range=10.0, target_threshold=0.85)
    pfad_radar = radar.flug_mit_radar(start_x=15_300.0, steps=20)
    print(f"[*] Pfadlänge {len(pfad_radar)} | Ende x={pfad_radar[-1]:.4f} | locked={radar.locked}\n")

    print("--- MonadicZoomAutopilot.spektrale_tiefenbohrung ---")
    zoom = MonadicZoomAutopilot(zeros, N_demo)
    monad_ok, k_end = zoom.spektrale_tiefenbohrung(15_300.0, schritte=5)
    print(f"[*] Monade={monad_ok} | letzte Kohärenz={k_end:.6f}\n")

    print("--- monadic_deep_drill ---")
    from types import SimpleNamespace

    test_zentren = [15360.79, 15269.05, 15401.00]
    monadic_deep_drill(SimpleNamespace(all_zeros=zeros), N_demo, test_zentren)
    print()

    print("--- QuantumSwarmAutopilot.schwarm_flug ---")
    swarm = QuantumSwarmAutopilot(zeros, N_demo, n_probes=10)
    x_swarm = swarm.schwarm_flug(center_x=15_300.0, radius=100.0, steps=15)
    print(f"[*] Schwarm-Ergebnis x={x_swarm:.4f}\n")

    print("--- PolarizedSwarmAutopilot.schwarm_flug ---")
    polar = PolarizedSwarmAutopilot(zeros, N_demo)
    x_pol = polar.schwarm_flug(center_x=15_300.0, radius=100.0, steps=15)
    print(f"[*] Polarisierter Schwarm x={x_pol:.4f}\n")

    print("--- DampedRadarAutopilot.präzisions_anflug ---")
    damped = DampedRadarAutopilot(zeros, N_demo, max_step=20.0, radar_range=10.0, target_threshold=0.85)
    pfad_d = damped.präzisions_anflug(start_x=15_300.0, steps=30)
    print(f"[*] Pfadlänge {len(pfad_d)} | Ende x={pfad_d[-1]:.4f} | locked={damped.locked}\n")

    print("--- BreakoutAutopilot.flug_mit_ausbruch ---")
    breakout = BreakoutAutopilot(zeros, N_demo)
    pfad_bo = breakout.flug_mit_ausbruch(start_x=15_300.0, steps=30)
    print(f"[*] Pfadlänge {len(pfad_bo)} | Ende x={pfad_bo[-1]:.4f} | locked={breakout.locked}")
