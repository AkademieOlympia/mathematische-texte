from __future__ import annotations

import numpy as np


class GradientProbe:
    """
    Basis: Riemann-Nullstellen (gamma) und Ziel-N.
    Kohärenz ~ Mittelwert cos(γ log x) cos(γ log N) (wie Spektral-Sonden).
    Ohne ``n_zeros`` werden alle geladenen γ genutzt (kein stiller 5k-Deckel).
    """

    def __init__(self, zeros, N):
        self.zeros = np.asarray(zeros, dtype=float).ravel()
        self.N = float(N)
        if self.N <= 0:
            raise ValueError("N muss positiv sein.")
        self._ln_N = np.log(self.N)

    def measure_coherence(self, x: float, n_zeros: int | None = None) -> float:
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
        return float(np.mean(np.cos(g * ln_x) * np.cos(g * self._ln_N)))


class TopologischerAutopilot(GradientProbe):
    def __init__(self, zeros, N):
        super().__init__(zeros, N)
        self.persistence_threshold = 0.5 # Nur "langlebige" Strukturen zählen

    def berechne_topologische_guete(self, x, radien=[0.1, 0.5, 1.0, 2.0]):
        """
        Simuliert die Persistenz: Misst die Kohärenz über verschiedene 
        Skalierungen (Radien) der Robinson-Monade.
        """
        kohärenz_werte = []
        for r in radien:
            # Wir messen die "Dichte" in einer Umgebung r um x
            k = self.measure_coherence(x + r) + self.measure_coherence(x - r)
            kohärenz_werte.append(k)
        
        # Ein Merkmal ist "persistent", wenn es über alle Radien 
        # hinweg stabil (korreliert) bleibt.
        stabilität = np.std(kohärenz_werte) 
        return np.mean(kohärenz_werte) / (1 + stabilität)

    def smarter_flug(self, start_x, steps=20):
        """
        Navigiert die Sonde unter Umgehung des Rauschens.
        """
        x = start_x
        pfad = [x]
        
        for _ in range(steps):
            # Wir messen den Gradienten der TOPOLOGISCHEN GÜTE, nicht der rohen Amplitude
            eps = 0.1
            g1 = self.berechne_topologische_guete(x - eps)
            g2 = self.berechne_topologische_guete(x + eps)
            
            topo_grad = (g2 - g1) / (2 * eps)
            
            # Die Sonde "gleitet" zum stabilsten topologischen Punkt
            x = x - 50.0 * topo_grad 
            pfad.append(x)
            
        return pfad