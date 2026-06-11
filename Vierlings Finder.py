import numpy as np

def is_prime(n):
    if n < 2: return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0: return False
    return True

class AQGVierlingFinder:
    def __init__(self):
        self.phi = (1 + 5**0.5) / 2
        self.karl_factor = 1.5      # Skalierung 45°-Transformation
        self.anchor = 5005          # 5-7-11-13 Resonanz
        self.grid_const = 32        # 32-Glattheit

    def check_quadruplet_resonance(self, p_start):
        """
        Prüft, ob ein Primzahl-Vierling {p, p+2, p+6, p+8}
        mit dem 5005-Anker harmonisiert.
        """
        quad = [p_start, p_start + 2, p_start + 6, p_start + 8]
        if all(is_prime(p) for p in quad):
            # Summe der Resonanzen im Gitter
            resonance_sum = sum(quad) % self.grid_const
            # Prüfung der Nähe zum 5005-Anker (Mod-Resonanz)
            anchor_match = (sum(quad) * self.anchor) % self.grid_const
            return True, quad, anchor_match
        return False, None, None

    def search_channels(self, limit=10000):
        channels = []
        print(f"Starte Suche nach supraleitenden Kanälen im 32-Gitter...")
        
        for p in range(2, limit):
            is_res, quad, match = self.check_quadruplet_resonance(p)
            if is_res:
                # Hier greift die Morley-Logik: 
                # Nur "glatte" Ergebnisse werden als Kanal akzeptiert
                if match < 1.0: # Schwellenwert für arithmetische Glattheit
                    channels.append(quad)
                    print(f"Kanal gefunden bei Vierling {quad} (Resonanz-Güte: {match})")
        
        return channels


class AQGRouter:
    def __init__(self, channels):
        self.channels = channels  # Liste der gefundenen Vierlinge
        self.karl_factor = 1.5
        self.grid_const = 32

    def calculate_hop_resonance(self, q1, q2):
        """
        Prüft, ob zwei Vierlinge (q1, q2) über einen
        Bose-Einstein-Kanal verbunden sind.
        """
        # Distanz zwischen den Zentren der Vierlinge
        dist = abs(np.mean(q1) - np.mean(q2))

        # Stabilität im 32-Gitter nach Karl-Skalierung
        bridge_stability = (dist * self.karl_factor) % self.grid_const
        return bridge_stability < 1.0  # Schwellenwert für Kopplung

    def build_routing_table(self):
        routing_table = {}
        for i, q1 in enumerate(self.channels):
            connections = []
            for j, q2 in enumerate(self.channels):
                if i != j and self.calculate_hop_resonance(q1, q2):
                    connections.append(j)
            routing_table[i] = connections
        return routing_table


def aqg_compute_step(integer_input):
    # 1. Energetisierung: Transformation in den Riemann-Raum (Dodekaeder)
    energy_state = np.exp(1j * (integer_input / 32) * np.pi)

    # 2. Supraleitende Suche: Das Signal "rutscht" in den naechsten Vierling
    # Bei komplexem Zustand verwenden wir den Betrag als skalare Resonanzstaerke.
    resonance = np.mod(np.abs(energy_state * 5005), 1.0)

    # Supraleitungs-Effekt: nahe 0 => instantane Stabilisierung (Tunneling)
    if resonance < 0.05:
        stable_output = np.round(integer_input * 1.5)  # Karl-Sprung
    else:
        stable_output = integer_input  # Kein Tunneling moeglich

    # 3. Rekristallisation: Zurueck in die Ganzzahl-Welt
    return int(stable_output % 32)


if __name__ == "__main__":
    # --- Ausführung ---
    explorer = AQGVierlingFinder()
    found_channels = explorer.search_channels(5000)
    print(f"\nSuche beendet. {len(found_channels)} supraleitende BEK-Kanaele identifiziert.")

    # --- Integration Router ---
    router = AQGRouter(found_channels)
    table = router.build_routing_table()

    print("Routing-Tabelle erstellt.")
    for node, targets in list(table.items())[:5]:  # Zeige erste 5 Knoten
        print(f"Vierling-Knoten {node} {found_channels[node]} koppelt mit: {targets}")

print(f"\nSuche beendet. {len(found_channels)} supraleitende BEK-Kanäle identifiziert.")