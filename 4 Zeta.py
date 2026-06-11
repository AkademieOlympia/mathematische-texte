import numpy as np
import matplotlib.pyplot as plt
import warnings

# ============================================================
# KONFIGURATION
# ============================================================
N = 200_000         # Bereichsgrenze
APPLY_LAMB_SHIFT = True

# Lamb-Shift Parameter
alpha = 1.0 / 137.035999
# Wir testen "was passiert bei 137" -> Wir nehmen alpha^1 (starke Kopplung)
# statt alpha^5, damit der Effekt deutlich ist.
coupling = alpha ** 1 

VISIBILITY_BOOST = 1.0 

spin_factors = {1: 0.5, 5: 1.0, 7: 1.0, 11: 1.5}

def safe_show():
    """
    Versucht den Plot zu speichern, statt ihn anzuzeigen, um Blockaden zu vermeiden.
    """
    try:
        output_file = "4_Zeta_137.png"
        print(f"Speichere Plot als '{output_file}' (statt Anzeige)...")
        plt.savefig(output_file, dpi=300)
        print("Plot gespeichert.")
    except KeyboardInterrupt:
        print("\nPlot durch Benutzer abgebrochen.")
    except Exception as e:
        print(f"\nFehler beim Speichern: {e}")

# ============================================================
# DATENGENERIERUNG (Optimiert & Vektorisiert)
# ============================================================
def prime_sieve(n):
    """Schnelles Sieb des Eratosthenes mit NumPy."""
    sieve = np.ones(n + 1, dtype=bool)
    sieve[0:2] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return np.nonzero(sieve)[0]

def calculate_lamb_shift_vectorized(primes, residues):
    """Vektorisierte Berechnung des Lamb-Shifts."""
    if not APPLY_LAMB_SHIFT:
        return np.zeros_like(primes, dtype=float)

    # Basis-Shift aus alpha/log(p)
    # Vermeide log(1) oder log(0), primes starten bei 5
    with np.errstate(divide='ignore'):
        base_shift = (coupling * VISIBILITY_BOOST) / np.log(primes)
    
    # Spin-Faktoren anwenden
    # E=1, A=5, B=7, C=11
    spin_mult = np.zeros_like(primes, dtype=float)
    spin_mult[residues == 1]  = spin_factors[1]
    spin_mult[residues == 5]  = spin_factors[5]
    spin_mult[residues == 7]  = spin_factors[7]
    spin_mult[residues == 11] = spin_factors[11]
    
    return base_shift * spin_mult

print(f"Generiere Primzahlen bis {N}...")
primes = prime_sieve(N)

# Nur Primzahlen > 3 (Filter für mod 12 Familien)
mask_p = primes > 3
active_primes = primes[mask_p]

print("Sortiere Primzahlen in Familien (Vektorisiert)...")
residues = active_primes % 12

# Arrays für Signale (0 bis N)
# Wir nutzen float für Amplitudenmodulation
signals = {
    "E (1 mod 12)": np.zeros(N + 1, dtype=float),
    "A (5 mod 12)": np.zeros(N + 1, dtype=float),
    "B (7 mod 12)": np.zeros(N + 1, dtype=float),
    "C (11 mod 12)": np.zeros(N + 1, dtype=float)
}

# Indizes für die Familien
mask_E = (residues == 1)
mask_A = (residues == 5)
mask_B = (residues == 7)
mask_C = (residues == 11)

# Amplituden berechnen
shifts = calculate_lamb_shift_vectorized(active_primes, residues)
amps = 1.0 - shifts

# Debug für p=137
idx_137 = np.where(active_primes == 137)[0]
val_137 = None
if len(idx_137) > 0:
    i = idx_137[0]
    val_137 = (active_primes[i], shifts[i], amps[i])
    print(f"--- Spezial-Check p=137 ---")
    print(f"Familie: {residues[i]} (5=A)")
    print(f"Shift: {shifts[i]:.6e}")
    print(f"Amp: {amps[i]:.6f}")
    print("---------------------------")

# Signale füllen
signals["E (1 mod 12)"][active_primes[mask_E]] = amps[mask_E]
signals["A (5 mod 12)"][active_primes[mask_A]] = amps[mask_A]
signals["B (7 mod 12)"][active_primes[mask_B]] = amps[mask_B]
signals["C (11 mod 12)"][active_primes[mask_C]] = amps[mask_C]

print("Berechne Fourier-Transformationen (Die 4 Zeta-Spektren)...")
ffts = {}
for label, sig in signals.items():
    # Wir nehmen den Betrag der FFT (Spektrum)
    # Nur positive Frequenzen bis N/2
    spec = np.abs(np.fft.rfft(sig))
    ffts[label] = spec

# ============================================================
# PLOTTING
# ============================================================
print("Erstelle Plot...")
try:
    fig, axes = plt.subplots(4, 1, figsize=(12, 16), sharex=True)
    
    # Farben für die Familien (E, A, B, C -> Gold, Rot, Grün, Blau)
    colors = {
        "E (1 mod 12)": "gold",
        "A (5 mod 12)": "crimson", 
        "B (7 mod 12)": "forestgreen",
        "C (11 mod 12)": "royalblue"
    }

    # Frequenzachse für rfft (0 bis 0.5)
    freqs = np.fft.rfftfreq(N + 1)
    
    # Plot für jede Familie
    for ax, (label, spec) in zip(axes, ffts.items()):
        color = colors.get(label, "black")
        
        # Plot spectrum
        ax.plot(freqs[1:], spec[1:], color=color, lw=0.5, alpha=0.8)
        
        # Titel & Labels
        ax.set_ylabel(f"Amplitude (korrigiert)", fontsize=10)
        ax.set_title(f"Familie {label}", fontsize=12, fontweight="bold", color=color)
        ax.grid(True, alpha=0.3)
        
        # Markiere 137 im entsprechenden Plot (A-Familie)
        if "A (5 mod 12)" in label and val_137:
            # Das ist ein Spektrum-Plot, p=137 ist hier keine Frequenz,
            # sondern ein Beitrag zum Signal.
            # Aber wir können einen Text einfügen, was passiert ist.
            shift_text = f"p=137 gedämpft um {val_137[1]:.1e}"
            ax.text(0.02, 0.9, shift_text, transform=ax.transAxes, 
                   fontsize=10, color='red', bbox=dict(facecolor='white', alpha=0.8))

    # Beschriftung und Layout
    title_text = f"Spektrale Zerlegung der Primzahlen (N={N}) in 4 Familien\n(Die '4 Zeta'-Hypothese)"
    if APPLY_LAMB_SHIFT:
        title_text += f" + Lamb-Shift (alpha^1, p=137 Effekt)"

    # Gesamt-Titel über alle Subplots
    fig.suptitle(title_text, fontsize=16)
    plt.xlabel("Frequenz (normiert)", fontsize=12)
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])

    safe_show()
    
except Exception as e:
    print(f"Fehler beim Plotting: {e}")
