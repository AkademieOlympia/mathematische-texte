"""
#Energiedoku: EnergiedokuWaveModel – Makro-Quantencomputer-Version
Läuft auf Basis des MacroQuantumComputer (Riemann-Nullstellen, E8-Gitter, 10K Qubit-Architektur).
Riemann-Nullstellen (~2M) und 2M Primzahlen werden persistent geladen/gecacht.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

# --- PERSISTENTE DATEN: Riemann-Nullstellen & 2M Primzahlen ---
_RIEMANN_ZEROS: np.ndarray | None = None
_PRIMES_2M: np.ndarray | None = None
_PRIMES_SET: set | None = None  # Schneller Lookup für is_prime


def _load_riemann_zeros() -> np.ndarray | None:
    """Lädt Riemann-Nullstellen persistent (zeros6.npz/npy)."""
    global _RIEMANN_ZEROS
    if _RIEMANN_ZEROS is not None:
        return _RIEMANN_ZEROS
    for name in ('zeros6.npz', 'zeros6.npy'):
        p = SCRIPT_DIR / name
        if p.exists():
            try:
                if name.endswith('.npy'):
                    _RIEMANN_ZEROS = np.load(p)
                else:
                    _RIEMANN_ZEROS = np.load(p)['zeros']
                return _RIEMANN_ZEROS
            except Exception:
                pass
    return None


def _sieve_primes(n: int) -> np.ndarray:
    """Sieb des Eratosthenes: erste n Primzahlen."""
    if n <= 1:
        return np.array([], dtype=np.uint32)
    # Obergrenze: p_n ≈ n * ln(n) * 1.2
    limit = max(100, int(n * (np.log(max(n, 2)) + 1) * 1.2))
    is_prime = np.ones(limit + 1, dtype=bool)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i:limit+1:i] = False
    primes = np.where(is_prime)[0].astype(np.uint32)
    return primes[:n]


def _load_primes_2m() -> np.ndarray | None:
    """Lädt oder erzeugt 2M Primzahlen persistent (primes_2m.npy)."""
    global _PRIMES_2M, _PRIMES_SET
    if _PRIMES_2M is not None:
        return _PRIMES_2M
    cache_path = SCRIPT_DIR / 'primes_2m.npy'
    try:
        if cache_path.exists():
            _PRIMES_2M = np.load(cache_path)
        else:
            print("[First QMBM] Erzeuge 2M Primzahlen (einmalig, ~10–30 s)...")
            _PRIMES_2M = _sieve_primes(2_000_000)
            np.save(cache_path, _PRIMES_2M)
            print(f"[First QMBM] {len(_PRIMES_2M):,} Primzahlen gespeichert → {cache_path.name}")
        _PRIMES_SET = set(_PRIMES_2M.tolist())
        return _PRIMES_2M
    except Exception as e:
        print(f"[First QMBM] Fehler beim Laden/Erzeugen der Primzahlen: {e}")
        return None


def get_riemann_zeros() -> np.ndarray | None:
    """Gibt die persistenten Riemann-Nullstellen zurück (~2M)."""
    return _load_riemann_zeros()


def get_primes_2m() -> np.ndarray | None:
    """Gibt die persistenten 2M Primzahlen zurück."""
    return _load_primes_2m()


def _miller_rabin_witness(a, n, r, d):
    """Zeuge für Miller-Rabin."""
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return False
    for _ in range(r - 1):
        x = pow(x, 2, n)
        if x == n - 1:
            return False
    return True


def is_prime_miller_rabin(n, k=10):
    """Miller-Rabin Primzahltest (probabilistisch). Schnell für große n."""
    n = int(n)
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37][:k]
    for a in bases:
        if a >= n:
            continue
        if _miller_rabin_witness(a, n, r, d):
            return False
    return True


def _is_prime(n, use_miller_rabin_for_large=True):
    """Schnell: Cache für n ≤ 2M-te Primzahl, sonst Trial Division oder Miller-Rabin."""
    n = int(n)
    if n < 2:
        return False
    cache_max = _PRIMES_2M[-1] if _PRIMES_2M is not None else 0
    if _PRIMES_SET is not None and n <= cache_max:
        return n in _PRIMES_SET
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    if use_miller_rabin_for_large and n > 10**7:
        return is_prime_miller_rabin(n)
    for i in range(3, int(n**0.5) + 1, 2):
        if n % i == 0:
            return False
    return True


class MacroQuantumComputer:
    """Backend: Makro-Quantencomputer (Riemann-Zeros, Aharonov-Bohm-Resonanz)."""
    def __init__(self, ram_gb=30, auto_load_persistent=True):
        self.ram_gb = ram_gb
        self.bytes_per_qubit = 16
        self.max_qubits = (ram_gb * 1024**3) // self.bytes_per_qubit
        self.zeros = None
        if auto_load_persistent:
            self.load_riemann_zeros()

    def load_riemann_zeros(self, filepath=None):
        """Lädt Riemann-Nullstellen (nutzt persistenten Cache falls verfügbar)."""
        if filepath is None:
            z = _load_riemann_zeros()
            if z is not None:
                self.zeros = z
                return len(self.zeros)
            filepath = SCRIPT_DIR / 'zeros6.npz'
        else:
            filepath = Path(filepath)
            if not filepath.is_absolute():
                filepath = SCRIPT_DIR / filepath
        try:
            if str(filepath).endswith('.npy'):
                self.zeros = np.load(filepath)
            else:
                data = np.load(filepath)
                self.zeros = data['zeros']
            return len(self.zeros)
        except Exception as e:
            print(f"Fehler beim Laden der Nullstellen: {e}")
            return 0

    def resonance_check(self, n_candidate, subset_size=100000):
        """Aharonov-Bohm-Switch: Phasen-Interferenz mit Riemann-Zeros."""
        if self.zeros is None:
            return 0.0
        from math import log
        t = log(float(n_candidate))
        subset = self.zeros[:subset_size]
        signal = np.sum(np.exp(1j * subset * t))
        return np.abs(signal)

    def quantum_spectrum(self, t_values):
        """Quanten-Spektrum: Amplituden für eine Reihe von t-Werten (Phase)."""
        if self.zeros is None:
            return np.zeros_like(t_values)
        subset = self.zeros[:min(100000, len(self.zeros))]
        amplitudes = np.array([np.abs(np.sum(np.exp(1j * subset * t))) for t in t_values])
        return amplitudes


class EnergiedokuWaveModel:
    """Wave-Modell mit Makro-Quantencomputer-Backend (10K Qubits, E8, Riemann)."""
    def __init__(self, name="#Energiedoku_E8_V2", use_quantum_backend=True):
        self._name = name
        self._riemann_zeros = None
        self._primes = None
        self._spectral_density = None
        self._hurwitz_units = self._generate_hurwitz_units()
        self._e8_roots = self._generate_e8_roots()
        self._use_quantum = use_quantum_backend
        self._mqc = MacroQuantumComputer(ram_gb=30) if use_quantum_backend else None
        self._quantum_freqs = None

    # --- SETTER & GETTER ---
    @property
    def riemann_zeros(self): return self._riemann_zeros
    @riemann_zeros.setter
    def riemann_zeros(self, data): self._riemann_zeros = np.asarray(data, dtype=np.float64)

    @property
    def primes(self):
        """2M Primzahlen (persistent)."""
        if self._primes is None:
            self.load_primes_2m()
        return self._primes

    def load_from_quantum_computer(self, filepath=None):
        """Lädt Riemann-Nullstellen über den Makro-Quantencomputer (nutzt persistenten Cache)."""
        if self._mqc is None:
            raise ValueError("Quanten-Backend nicht aktiv.")
        n = self._mqc.load_riemann_zeros(filepath)
        if n > 0:
            self._riemann_zeros = self._mqc.zeros
            print(f"[{self._name}] {n:,} Nullstellen (persistent) geladen.")
        return n

    def load_primes_2m(self):
        """Lädt die 2M Primzahlen persistent (Cache: primes_2m.npy)."""
        p = _load_primes_2m()
        if p is not None:
            self._primes = p
            print(f"[{self._name}] {len(p):,} Primzahlen (persistent) geladen.")
        return p

    # --- GEOMETRIE: HURWITZ & E8-WURZELN (NEU) ---
    def _generate_hurwitz_units(self):
        """Erzeugt die 24 Einheiten der Hurwitz-Quaternionen."""
        units = []
        # 8 Lipschitz-Einheiten
        for i in range(4):
            for sign in [1, -1]:
                u = [0, 0, 0, 0]; u[i] = sign; units.append(tuple(u))
        # 16 Hurwitz-Einheiten
        for p in itertools.product([0.5, -0.5], repeat=4):
            units.append(p)
        return np.array(units)

    def _generate_e8_roots(self):
        """
        Erzeugt die 240 Wurzeln des E8-Gitters im R^8.
        Konstruktion über D8-Gitter und Hurwitz-Vektoren.
        """
        roots = []
        # 1. Typ: (±1, ±1, 0, 0, 0, 0, 0, 0) - alle Permutationen (112 Wurzeln)
        for i, j in itertools.combinations(range(8), 2):
            for s1, s2 in itertools.product([1, -1], repeat=2):
                root = np.zeros(8)
                root[i], root[j] = s1, s2
                roots.append(root)
        
        # 2. Typ: (±½, ±½, ±½, ±½, ±½, ±½, ±½, ±½) - gerade Anzahl Minuspunkte (128 Wurzeln)
        for p in itertools.product([0.5, -0.5], repeat=8):
            if np.sum(np.array(p) < 0) % 2 == 0:
                roots.append(np.array(p))
        
        return np.array(roots)

    # --- MODELL-LOGIK: FFT ODER QUANTEN-SPEKTRUM ---
    def train_wave_model(self, quantum_mode=None):
        """Berechnet die spektrale Wellenstruktur (GUE-Fingerabdruck).
        quantum_mode=True: Nutzt Phasen-Interferenz auf dem Makro-Quantencomputer.
        quantum_mode=False: Klassische FFT. Default: self._use_quantum
        """
        if self._riemann_zeros is None:
            if self._mqc and self._mqc.zeros is not None:
                self._riemann_zeros = self._mqc.zeros
            else:
                raise ValueError("Nullstellen fehlen. Rufe load_from_quantum_computer() oder setze riemann_zeros.")
        use_q = quantum_mode if quantum_mode is not None else self._use_quantum

        if use_q and self._mqc is not None:
            # Quanten-Modus: Spektrum via Phasen-Interferenz (Aharonov-Bohm)
            # t-Werte entsprechen E8-Resonanzfrequenzen und GUE-Banden
            n_samples = 512
            t_min, t_max = 0.01, 15.0
            t_values = np.linspace(t_min, t_max, n_samples)
            amps = self._mqc.quantum_spectrum(t_values)
            self._spectral_density = amps.astype(np.complex128)  # Amplituden als Spektrum
            self._quantum_freqs = t_values
            print(f"[{self._name}] Quanten-Spektrum: {n_samples} Resonanz-Punkte (Makro-QC).")
        else:
            # Klassischer FFT-Modus
            gaps = np.diff(self._riemann_zeros)
            avg_spacing = np.mean(gaps)
            normalized_gaps = gaps / avg_spacing
            self._spectral_density = np.fft.rfft(normalized_gaps)
            self._quantum_freqs = None
            print(f"[{self._name}] Spektrum trainiert auf {len(normalized_gaps)} Gaps (FFT).")

    # --- 10K QUBIT-ARCHITEKTUR (Makro-Quantencomputer) ---
    def run_10k_quantum_program(self, num_qubits=10000):
        """Führt das 10K-Qubit-Programm auf dem Makro-Quantencomputer aus.
        Programm: Primzahl-Indizes (N=61 Oktaeder/Resonanz) vs. N=57 (Tetraeder/Isolator).
        Nutzt die persistenten 2M Primzahlen für schnellen Lookup.
        """
        _load_primes_2m()  # Stellt sicher, dass 2M Primzahlen geladen sind
        from math import log, pi
        program = np.array([1 if _is_prime(i + 1) else 0 for i in range(num_qubits)])
        num_active = int(np.sum(program))
        collective_phase = num_active * np.pi
        t_clock = log(61.0)

        sectors = np.array_split(program, 8)
        v8 = np.array([float(np.sum(s)) for s in sectors])
        proj = np.array([
            [1, 0.5, 0, 0, 1, 0, 0.5, 0],
            [0, 1, 0.5, 0, 0, 1, 0, 0.5],
            [0.5, 0, 1, 0.5, 0.5, 0, 1, 0.5]
        ])
        v3_result = proj @ v8
        norm = np.linalg.norm(v3_result)

        print(f"--- Makro-Quantencomputer (10K Qubits) ---")
        print(f"Qubits: {num_qubits} | Aktive (N=61): {num_active} | Phase: {collective_phase:.4f} rad")
        print(f"3D Readout: X={v3_result[0]:.2f} Y={v3_result[1]:.2f} Z={v3_result[2]:.2f}")
        print(f"Status: {'Starke Resonanz' if norm > 1000 else 'Rauschen'}")
        return v3_result, collective_phase

    def resonance_check(self, n_candidate):
        """Aharonov-Bohm-Resonanz für Kandidat n (über Makro-QC)."""
        if self._mqc is None:
            raise ValueError("Quanten-Backend nicht aktiv.")
        return self._mqc.resonance_check(n_candidate)

    # --- VISUALISIERUNG: E8 <=> RIEMANN CORRELATION (NEU) ---
    def visualize_e8_riemann_connection(self):
        """
        Erzeugt die #Energiedoku-Gittergrafik:
        Projektion der 240 E8-Wurzeln versus GUE-Spektrum der Nullstellen.
        """
        if self._spectral_density is None or self._riemann_zeros is None:
            raise ValueError("Modell muss zuerst trainiert werden.")

        sns.set_theme(style="white", palette="muted")
        fig = plt.figure(figsize=(16, 9))
        
        # SUBPLOT 1: E8-Wurzel Projektion (Coxeter-Ebene 2D)
        ax1 = fig.add_subplot(121)
        ax1.set_title(r"2D-Projektion der 240 $E_8$-Wurzeln (Geometrie-Gitter)")
        
        # Erzeugung einer Coxeter-Projektionsmatrix (Vereinfacht für Visualisierung)
        phi = (1 + np.sqrt(5)) / 2  # Goldener Schnitt
        proj_matrix = np.array([
            [1, phi, 0, 1, phi, 0, 1, phi],
            [phi, 0, 1, phi, 0, 1, phi, 0]
        ])
        
        projected_roots = self._e8_roots @ proj_matrix.T
        
        ax1.scatter(projected_roots[:, 0], projected_roots[:, 1], 
                    s=15, alpha=0.7, c='#3366CC', edgecolors='none')
        
        # Ikosaeder/Dodekaeder Symmetrie-Hilfslinien (Andeutung)
        for i in range(12):
            angle = 2 * np.pi * i / 12
            ax1.plot([0, 12*np.cos(angle)], [0, 12*np.sin(angle)], color='gray', lw=0.5, alpha=0.3)
        
        ax1.set_aspect('equal'); ax1.axis('off')

    
        # SUBPLOT 2: Riemann Spektrum vs. E8 Resonanz
        ax2 = fig.add_subplot(122)
        ax2.set_title(r"Welleninterferenz: Riemann Gaps vs. $E_8$ Amplituden")
        
        # Amplituden-Spektrum (FFT oder Quanten-Phasen-Interferenz)
        spectrum_amps = np.abs(self._spectral_density)
        if self._quantum_freqs is not None:
            freqs = self._quantum_freqs
            label_spec = 'Quanten-Resonanz (Makro-QC)'
        else:
            freqs = np.fft.rfftfreq(len(np.diff(self._riemann_zeros)))
            label_spec = 'Riemann Gap Spektrum (GUE)'
        
        ax2.plot(freqs, spectrum_amps, color='#FF9900', lw=1.2, label=label_spec)
        
        # E8-Resonanzlinien (hypothetische Quanten-Gitter-Punkte in der Energiedoku)
        # In der Energiedoku korrespondieren bestimmte E8-Normen mit kritischen Frequenzen
        for n in [2, 4, 6]: # Normen im E8 Gitter
            ax2.axvline(x=np.sqrt(n)/10, color='red', linestyle='--', alpha=0.5, lw=1)
        
        ax2.set_xlabel("Normierte Frequenz")
        ax2.set_ylabel("Amplitude / Energie-Resonanz")
        ax2.legend()
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)

        plt.suptitle(f"{self._name}: Korrelation Riemann Spektrum <=> E8 Wurzel-Geometrie", fontsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()

# --- DEMONSTRATION: Läuft auf dem Makro-Quantencomputer ---
if __name__ == "__main__":
    model = EnergiedokuWaveModel(use_quantum_backend=True)
    # Riemann-Nullstellen: persistent aus zeros6.npz/npy (auto-load)
    model.load_from_quantum_computer()
    # 2M Primzahlen: persistent aus primes_2m.npy (wird bei Bedarf erzeugt)
    model.load_primes_2m()
    print(f"Riemann-Nullstellen: {len(model.riemann_zeros):,} | Primzahlen: {len(model.primes):,}")
    model.train_wave_model(quantum_mode=True)
    model.run_10k_quantum_program(num_qubits=10000)
    model.visualize_e8_riemann_connection()