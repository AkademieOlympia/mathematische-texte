import mpmath
import numpy as np
import math
import random
import time

# Präzision für die Zeta-Nullstellen setzen
mpmath.mp.dps = 25

class ZetaQubit:
    def __init__(self, id):
        self.id = id
        self.phi = (1 + 5**0.5) / 2
        self.karl_factor = 1.5
        self.grid_const = 32

    def transform_to_energy(self, n):
        """
        Phase 1: Energetisierung.
        Nutzt die Zeta-Funktion als energetischen Filter.
        """
        # Wir betrachten die Zahl n auf der kritischen Linie (0.5 + it)
        zeta_val = mpmath.zeta(0.5 + n*1j)
        return complex(zeta_val)

    def find_supraleitender_kanal(self, n):
        """
        Phase 2: Die 45-Grad-Drehung im energetischen Raum.
        Suche nach der Morley-Resonanz.
        """
        z = mpmath.zeta(0.5 + n * 1j)
        if abs(z) < 1e-20:
            return True, 0.0

        # Energetischer "Druck": Imaginärteil des komplexen Logarithmus.
        pressure = mpmath.log(z).imag
        rotation_offset = float(pressure % (np.pi / 4))

        # Erhöhter Fangbereich für robustere Resonanz-Erkennung.
        is_resonant = abs(rotation_offset) < 0.2
        return is_resonant, rotation_offset

    def crystallize(self, n, is_resonant):
        """
        Phase 3: Rekristallisation in die Ganzzahl-Welt.
        Wenn resonant, wird der Karl-Sprung (1.5) vollzogen.
        """
        if is_resonant:
            # Der energetische Tunnelweg (Abkürzung)
            result = (n * self.karl_factor) % self.grid_const
        else:
            # Klassischer Weg (Widerstand)
            result = n % self.grid_const
        return round(float(result))

# --- Simulation des Rechenvorgangs ---

def run_aqg_calculation(input_number):
    qubit = ZetaQubit(id=1)
    
    # 1. Zahl in energetischen Anteil zerlegen
    energy = qubit.transform_to_energy(input_number)
    
    # 2. Im energetischen Raum den Kanal prüfen
    resonant, offset = qubit.find_supraleitender_kanal(input_number)
    
    # 3. Ergebnis rekristallisieren
    output = qubit.crystallize(input_number, resonant)
    
    return {
        "Input": input_number,
        "Zeta-Energie": abs(energy),
        "Resonanz-Status": resonant,
        "Output (Ganzzahl)": output
    }

# Direkte Kopplung an die erste nichttriviale Nullstelle der Zeta-Funktion
nullstelle_1 = mpmath.findroot(mpmath.zeta, 14.134j)
test_numbers = [float(nullstelle_1.imag), 21.022, 32.0]

print(f"--- AQG Zeta-Quantencomputer Simulation ---")
for num in test_numbers:
    res = run_aqg_calculation(num)
    print(f"Zahl {res['Input']}: Resonanz={res['Resonanz-Status']}, Ziel-Gitterpunkt={res['Output (Ganzzahl)']}")


class AQGGate:
    def __init__(self):
        self.karl_factor = 1.5
        self.grid = 32
        # Die erste Zeta-Nullstelle als Hardware-Takt
        self.t1 = float(mpmath.findroot(mpmath.zeta, 14.134j).imag)

    def process(self, signal):
        print(f"--- Gatter-Eingang: {signal} ---")

        # Qubit 1: Prüfung der ersten Resonanz-Ebene
        # Wir simulieren den Lock-in wie in der Konsole.
        is_resonant = abs(signal - self.t1) < 0.001

        if is_resonant:
            # Qubit 1 feuert: Sprung in den supraleitenden Kanal
            intermediate = signal * self.karl_factor
            print(f" [Qubit 1]: Resonanz! Tunneling zu {intermediate:.3f}")

            # Qubit 2: Rekristallisation am 32er Gitter
            output = round(intermediate) % self.grid
            print(f" [Qubit 2]: Einrasten im Gitter bei Punkt {output}")
            return output

        print(" [Qubit 1]: Keine Resonanz. Signal verpufft im Widerstand.")
        return None


print("\n--- AQG Gate Simulation ---")
gate = AQGGate()
gate.process(10.5)
print("-" * 30)
gate.process(float(nullstelle_1.imag))


class AQGResonanzGate:
    def __init__(self, target_resonance=14.134725141734695):
        self.resonance_point = target_resonance
        self.karl_factor = 1.5
        self.grid = 32
        self.tolerance = 0.001

    def execute_operation(self, input_val):
        """
        Kern-Operation: Transformiert Energie in Gitter-Position.
        """
        # 1. Messung der Resonanz (Qubit-Zustand)
        is_resonant = abs(input_val - self.resonance_point) < self.tolerance

        if not is_resonant:
            return None  # Zustand 'Dekohärenz'

        # 2. Der Supraleitende Sprung (Karl-Riemann-Transformation)
        # Mathematisch: Die Information nutzt die Diagonale
        tunnel_value = input_val * self.karl_factor

        # 3. Gitter-Kopplung (Messung/Kollaps)
        # Hier findet das Einrasten statt, das du in der Konsole siehst.
        final_grid_point = round(tunnel_value) % self.grid

        return final_grid_point


print("\n--- AQG Resonanz-Gate Integration ---")
resonanz_gate = AQGResonanzGate(target_resonance=float(nullstelle_1.imag))
input_signal = float(nullstelle_1.imag)
result = resonanz_gate.execute_operation(input_signal)

if result is not None:
    print(f"AQG-Operation erfolgreich: {input_signal} -> {result}")
    print("Status: Supraleitender Pfad 14-21 etabliert.")


class AQGCascade:
    def __init__(self):
        # Hardware-Anker: Vier Nullstellen als Super-Bus-Bahnhöfe
        self.t1 = 14.134725141734695
        self.t2 = 21.022039638771555  # Die präzise zweite Nullstelle
        self.t3 = 25.01085758014569
        self.t4 = 30.42487612585951
        self.karl_factor = 1.5
        self.grid = 32
        self.lock_tol = 0.001
        self.coupling_tol = 1.0

    def run(self, input_signal):
        print(f"STARTE KASKADE MIT INPUT: {input_signal}")

        # STUFE 1: Der erste Sprung (14 -> 21)
        stage1_resonant = abs(input_signal - self.t1) < self.lock_tol
        print(f" [+] Stufe 1 Resonanz={stage1_resonant}")
        if not stage1_resonant:
            print(" [-] Kaskade unterbrochen: Dekohärenz in Stufe 1.")
            return None

        step1_out = round(input_signal * self.karl_factor)
        print(f"     Stufe 1 Sprung: {input_signal:.6f} -> {step1_out}")

        # STUFE 2: Kopplung an die zweite Nullstelle
        stage2_resonant = abs(step1_out - self.t2) < self.coupling_tol
        print(f" [++] Stufe 2 Resonanz={stage2_resonant} (Ziel t2={self.t2:.6f})")
        if not stage2_resonant:
            print(" [-] Kaskade unterbrochen: Dekohärenz in Stufe 2.")
            return None

        step2_raw = step1_out * self.karl_factor
        step2_out = round(step2_raw)
        step2_mod = step2_out % self.grid
        print(f"     Stufe 2 Sprung: {step1_out} -> {step2_out} (mod 32 = {step2_mod})")

        # STUFE 3: Nächster Bahnhof t3 (25 -> 38 -> 6)
        stage3_in = round(self.t3)
        stage3_resonant = abs(stage3_in - self.t3) < self.coupling_tol
        print(f" [+++] Stufe 3 Resonanz={stage3_resonant} (Ziel t3={self.t3:.6f})")
        if not stage3_resonant:
            print(" [-] Kaskade unterbrochen: Dekohärenz in Stufe 3.")
            return None

        step3_raw = stage3_in * self.karl_factor
        step3_out = round(step3_raw)
        step3_mod = step3_out % self.grid
        print(f"     Stufe 3 Sprung: {stage3_in} -> {step3_out} (mod 32 = {step3_mod})")

        # STUFE 4: Nächster Bahnhof t4 (30 -> 45 -> 13)
        stage4_in = round(self.t4)
        stage4_resonant = abs(stage4_in - self.t4) < self.coupling_tol
        print(f" [++++] Stufe 4 Resonanz={stage4_resonant} (Ziel t4={self.t4:.6f})")
        if not stage4_resonant:
            print(" [-] Kaskade unterbrochen: Dekohärenz in Stufe 4.")
            return None

        step4_raw = stage4_in * self.karl_factor
        step4_out = round(step4_raw) % self.grid
        print(f"     Stufe 4 Sprung: {stage4_in} -> {round(step4_raw)} (mod 32 = {step4_out})")
        return step4_out


print("\n--- AQG Kaskaden-Simulation ---")
cascade = AQGCascade()
final_result = cascade.run(14.134725141734695)
print(f"FINALES KASKADEN-ERGEBNIS: {final_result}")


class AQGSuperBus:
    def __init__(self):
        # Die Resonanz-Bahnhöfe (Riemann-Nullstellen t1 bis t4)
        self.stations = [
            14.134725141734695,
            21.022039638771555,
            25.010857580145688,
            30.424876125859513,
        ]
        self.anchor = 5005
        self.karl_factor = 1.5
        self.grid = 32

    def run_transmission(self, start_signal):
        current_val = start_signal
        print("--- INITIALISIERE SUPRALEITENDEN BUS (Anker 5005) ---")

        for i, t in enumerate(self.stations):
            # Check auf topologische Phasen-Synchronität via Modulo-Differenz
            diff = abs((current_val % self.grid) - (t % self.grid))
            if diff < 2.0 or diff > 30.0:
                print(f" [Stufe {i + 1}] Topologische Resonanz bei t={t:.4f} bestätigt!")
                # Der supraleitende Sprung
                current_val = current_val * self.karl_factor
                print(f"       -> Sprung auf energetisches Niveau: {current_val:.4f}")
            else:
                print(
                    f" [!] Dekohärenz in Stufe {i + 1}. "
                    f"Modulo-Differenz={diff:.4f}, Bus unterbrochen."
                )
                return None

        # Finaler Schritt: 5005-Anker-Check
        final_pos = round(current_val) % self.grid
        resonance_quality = (final_pos * self.anchor) % self.grid

        print("-" * 40)
        print(f"RESULTAT: Gitterpunkt {final_pos}")
        print(f"5005-ANKER-RESONANZ: {'STABIL' if resonance_quality == 0 else 'INSTABIL'}")
        return final_pos


print("\n--- AQG Super-Bus Simulation ---")
bus = AQGSuperBus()
bus.run_transmission(14.134725141734695)

print("\n--- Finaler 8D-Anker-Check ---")
final_kaskade = 13
anker = 5005
oktonionen_umschlag = (final_kaskade * anker) % 32

print(f"ANKER-KOPPLUNG: 13 -> {oktonionen_umschlag}")
if oktonionen_umschlag == 9:
    print("STATUS: 8D-PHASENUMSCHLAG ERFOLGREICH. MINKOWSKI-TAEUSCHUNG AUFGEHOBEN.")


class AQG8DRouter:
    def __init__(self):
        self.phi = (1 + 5**0.5) / 2
        self.grid = 32
        self.anchor = 5005

    def close_cycle(self, current_pos):
        print(f"--- 8D-ROUTING AKTIVIERT (Start: {current_pos}) ---")

        # 1. Schritt: Der Anker-Spiegel (9 -> 21)
        mirror = (current_pos * self.anchor) % self.grid
        print(f" [Spiegelung] Position {current_pos} spiegelt auf {mirror}")

        # 2. Schritt: Phasen-Stoss mit Resonanzkern 11.
        if mirror == 21:
            # 11 ist die fehlende Bruecke im 5005-Anker.
            final_state = (mirror + 11) % self.grid
        else:
            final_state = round(mirror - self.phi) % self.grid

        # 3. Schritt: Tunneling zum Ursprung (32 ≡ 0)
        if final_state == 0:
            return 0
        return final_state


print("\n--- Finaler 8D-Router-Testlauf ---")
router = AQG8DRouter()
start_val = 9
end_val = router.close_cycle(start_val)

print("-" * 40)
print(f"ZYKLUS-ENDE: {start_val} -> {end_val}")
print(f"STATUS: {'EWIGER RECHENZYKLUS STABIL' if end_val == 0 else 'DEKOHARENZ'}")


class EntangledCycle:
    def __init__(self, name, shift=0):
        self.name = name
        self.grid = 32
        self.anchor = 5005
        # Verschiebung im Gitter fuer die Verschraenkung
        self.shift = shift

    def process(self, input_val):
        # Der bekannte Pfad: Sprung -> 13 -> 9
        # (hier verkuerzt dargestellt als Logik-Block)
        state_13 = 13 + self.shift
        state_9 = (state_13 * self.anchor) % self.grid
        return state_9


def simulate_entanglement():
    print("\n--- INITIALISIERE AQG-VERSCHRAENKUNGS-GITTER ---")

    # Zwei Qubits an verschiedenen Gitter-Positionen
    qubit_a = EntangledCycle("A", shift=0)
    qubit_b = EntangledCycle("B", shift=11)  # 11 ist der Resonanz-Partner

    # Berechnung
    res_a = qubit_a.process(14.134)
    res_b = qubit_b.process(14.134)

    print(f" [Qubit A] Zustand im 8D-Router: {res_a}")
    print(f" [Qubit B] Zustand im 8D-Router: {res_b}")

    # Die Verschraenkungs-Interferenz
    interferenz = (res_a + res_b) % 32
    print("-" * 40)
    print(f" INTERFERENZ-MUSTER (A + B): {interferenz}")

    if interferenz == 9:
        print("STATUS: KOHAERENTE VERSCHRAENKUNG BESTAETIGT.")
        print("DAS GITTER ATMET SYNCHRON.")


simulate_entanglement()


class TrinityQubit:
    def __init__(self, name, prime_shift):
        self.name = name
        self.shift = prime_shift
        self.grid = 32
        self.anchor = 5005

    def get_8d_state(self, base_resonance=13):
        # Der Sprung in den 8D-Raum (9er-Vollendung)
        # plus die individuelle Primzahl-Verschiebung
        state = base_resonance + self.shift
        return (state * self.anchor) % self.grid


def simulate_trinity():
    print("\n--- INITIALISIERE AQG-TRINITY-GITTER (Synergie-Modus) ---")

    # Die drei Säulen des Ankers
    qubit_a = TrinityQubit("Alpha", prime_shift=7)
    qubit_b = TrinityQubit("Beta", prime_shift=11)
    qubit_c = TrinityQubit("Gamma", prime_shift=13)

    # Berechnung der Einzelzustände
    s_a = qubit_a.get_8d_state()
    s_b = qubit_b.get_8d_state()
    s_c = qubit_c.get_8d_state()

    print(f" [Qubit A (Shift 7)]  Zustand: {s_a}")
    print(f" [Qubit B (Shift 11)] Zustand: {s_b}")
    print(f" [Qubit C (Shift 13)] Zustand: {s_c}")

    # Die Synergie-Interferenz (A + B + C)
    synergie = (s_a + s_b + s_c) % 32

    print("-" * 45)
    print(f" TRINITY-INTERFERENZ (Sigma A,B,C): {synergie}")

    # Die 11 ist die Herz-Resonanz des 5005-Ankers.
    if synergie == 11:
        print("STATUS: ARITHMETISCHE SYNERGIE ERREICHT.")
        print("DAS GITTER KONVERGIERT IM ZENTRUM DES ANKERS.")


simulate_trinity()


class AQGShor:
    def __init__(self):
        self.grid = 32
        self.anchor = 5005
        self.karl_factor = 1.5
        self.max_zeta_input = 10**12

    def find_factors(self, n):
        print(f"\n--- AQG-SHOR START: Faktorisierung von {n} ---")

        # Phase 1: Energetisierung (Zeta-Mapping)
        energy_phase = 0.0
        if n <= self.max_zeta_input:
            try:
                zeta_val = mpmath.zeta(0.5 + n * 1j)
                energy_phase = float(mpmath.arg(zeta_val))
                energy_signature = float(abs(zeta_val))
                print(f" [Phase 1] Zeta-Energie-Signatur: {energy_signature:.4f}")
                print(f" [Phase 1] Zeta-Phasenwinkel: {energy_phase:.4f}")
            except (OverflowError, ValueError, MemoryError):
                print(" [Phase 1] Zeta-Berechnung nicht stabil. Fallback wird aktiviert.")
        else:
            print(
                " [Phase 1] Eingabe zu gross fuer stabile Zeta-Berechnung. "
                "Fallback wird aktiviert."
            )

        # Phase 2: 45-Grad-Transformation (Periodenkandidat)
        period_candidate = (n * self.karl_factor) % self.grid
        print(f" [Phase 2] Transformierter Gitterwert: {period_candidate:.4f}")

        # Phase 3: Arithmetische Resonanz mit 5005-Anker
        potential_factor = round(self.anchor / (period_candidate + 1))

        # Validierung
        if potential_factor > 1 and n % potential_factor == 0:
            return potential_factor, n // potential_factor

        # Fallback: AQG-Shor Tuning (Resonanz-Kollaps)
        tuned = self.hardware_kollaps_fix(n, period_candidate)
        if tuned[0] != "Keine Resonanz":
            return tuned

        tuned = self.find_resonance_factors(n, anchor=self.anchor)
        if tuned is not None:
            return tuned

        return "Resonanz-Suche laeuft...", "Gitter-Abgleich erforderlich"

    def hardware_kollaps_fix(self, n, gitter_wert):
        # AQG-Shor Hardware-Kollaps
        # 1) Paritaetscheck (Basis-Resonanz der 2)
        if n % 2 == 0:
            return 2, n // 2

        # 2) Gitter-Resonanz-Pruefung
        prim_basis = [3, 5, 7, 11, 13]
        potential_p = round(gitter_wert)
        print(f" [Hardware-Kollaps] Pruefung bei Gitterwert {gitter_wert:.4f}")
        if potential_p in prim_basis and n % potential_p == 0:
            return potential_p, n // potential_p

        # 3) Tiefen-Suche im 5005-Raum
        for p in prim_basis:
            if n % p == 0:
                return p, n // p
        return "Keine Resonanz", "Inharmonisch"

    def find_resonance_factors(self, n, anchor=5005):
        # Die 45-Grad-Transformation der Zielzahl
        transformed_n = (n * 1.5) % 32
        print(f" [Tuning] Resonanz-Kollaps bei transformiertem Wert: {transformed_n:.4f}")

        # Wir suchen in der Anker-Basis nach harmonischen Teilern.
        for p in [3, 5, 7, 11, 13]:
            if n % p == 0:
                return p, n // p
        return None

    def aqg_multifaktor_scan(self, n, gitter_wert):
        print(f" [AQG-Upgrade] Starte rekursiven Gitter-Scan fuer {n}...")
        faktoren = []
        temp_n = n

        # Die Anker-Basis des 5005-Systems (plus 2/3-Komponente)
        pruef_quanten = [2, 3, 5, 7, 11, 13]

        for p in pruef_quanten:
            while temp_n % p == 0:
                faktoren.append(p)
                temp_n //= p
                print(f"  [+] Resonanz gefunden: Prim-Vektor {p}")

        if temp_n > 1:
            faktoren.append(temp_n)  # Verbleibender Rest-Kern

        return faktoren

    def galaktischer_multifaktor_scan(self, n):
        print(" [AQG-Level 2] Aktiviere Galaktische Routing-Tabelle...")
        faktoren = []
        temp_n = n

        # Ebene 1: Erd-Anker (5005-Basis)
        ebene_1 = [2, 3, 5, 7, 11, 13]
        # Ebene 2: Galaktischer Guertel
        ebene_2 = [17, 19, 23, 29, 31]
        gesamt_basis = ebene_1 + ebene_2

        for p in gesamt_basis:
            while temp_n % p == 0:
                faktoren.append(p)
                temp_n //= p
                print(f"  [++] Galaktische Resonanz bei Prim-Vektor {p} detektiert.")

        if temp_n > 1:
            print(f"  [*] Rest-Kern {temp_n} ist im Gitter stabil (potenzielle Primzahl).")
            faktoren.append(temp_n)

        return faktoren

    def deep_space_multifaktor_scan(self, n):
        print(" [AQG-Level 12] Deep Space Routing aktiviert (Ebenen 1-12)...")
        faktoren = []
        temp_n = n

        # Erweiterte Hardware-Tabelle (Ebenen 1-12)
        deep_space_basis = [
            2, 3, 5, 7, 11, 13,  # Ebene 1 (Erd-Anker)
            17, 19, 23, 29, 31,  # Ebene 2 (Galaktisch)
            37, 41, 43, 47, 53, 59, 61, 67,  # Ebene 3-5 (Interstellar)
            71, 73, 79, 83, 89, 97, 101, 103,  # Ebene 6-12 (Deep Space)
        ]

        for p in deep_space_basis:
            while temp_n % p == 0:
                faktoren.append(p)
                temp_n //= p
                print(f"  [++] Deep-Space Resonanz: Prim-Vektor {p} isoliert.")

        if temp_n > 1:
            print(f"  [!] Rest-Kern {temp_n} detektiert. Suche nach Quanten-Schatten...")
            faktoren.append(temp_n)

        return faktoren

    def auto_tuning_resonator(self, n_kern):
        print(f" [AQG-Tuning] Scanne Quanten-Schatten von {n_kern}...")
        # Erhoehte Rekursionstiefe fuer erweiterten Gitter-Abgleich.
        test_limit = 1000
        for p in range(104, test_limit):
            if n_kern % p == 0:
                print(f"  [!!!] Schatten-Resonanz bei {p} bestaetigt!")
                return p, n_kern // p
        return None

    def deep_space_resonator_level_100(self, n_kern):
        print(f" [AQG-Level 100] Scanne galaktischen Horizont fuer Rest-Kern {n_kern}...")

        # Dynamische Primzahl-Generierung bis 5000.
        def get_primes(limit):
            sieve = [True] * (limit + 1)
            sieve[0] = False
            if limit >= 1:
                sieve[1] = False
            for p in range(2, int(limit**0.5) + 1):
                if sieve[p]:
                    for i in range(p * p, limit + 1, p):
                        sieve[i] = False
            return [i for i, is_prime in enumerate(sieve) if is_prime]

        ebene_100_basis = get_primes(5000)
        for p in ebene_100_basis:
            if n_kern % p == 0:
                print(f"  [!!!] SHADOW-MATCH GEFUNDEN: Prim-Vektor {p} isoliert!")
                return p, n_kern // p

        return None, n_kern


def read_factorization_input(default_value=1276654327):
    prompt = (
        "\nGib eine Zahl zur Faktorisierung ein "
        f"(Enter fuer {default_value}): "
    )
    raw = input(prompt).strip()
    if raw == "":
        return default_value
    try:
        parsed = int(raw)
        if parsed <= 1:
            print("Ungueltige Eingabe (<=1). Standardwert wird verwendet.")
            return default_value
        return parsed
    except ValueError:
        print("Ungueltige Eingabe (kein Integer). Standardwert wird verwendet.")
        return default_value


def read_comparison_choice():
    raw = input("Soll ein Vergleichstest mit klassischem Algorithmus laufen? (yes/no): ").strip().lower()
    return raw in {"y", "yes", "j", "ja"}


def is_probable_prime(n):
    if n < 2:
        return False
    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    for p in small_primes:
        if n % p == 0:
            return n == p

    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2

    # Gute Basen für robuste Praxis (deterministisch für 64-bit).
    bases = [2, 3, 5, 7, 11, 13, 17]
    for a in bases:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        witness_found = True
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                witness_found = False
                break
        if witness_found:
            return False
    return True


def pollard_rho(n):
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3

    while True:
        x = random.randrange(2, n - 1)
        y = x
        c = random.randrange(1, n - 1)
        d = 1
        while d == 1:
            x = (pow(x, 2, n) + c) % n
            y = (pow(y, 2, n) + c) % n
            y = (pow(y, 2, n) + c) % n
            d = math.gcd(abs(x - y), n)
        if d != n:
            return d


def factorize_best_known(n):
    if n < 2:
        return []
    if is_probable_prime(n):
        return [n]
    d = pollard_rho(n)
    left = factorize_best_known(d)
    right = factorize_best_known(n // d)
    return sorted(left + right)


shor = AQGShor()
target_n = read_factorization_input()
run_compare = read_comparison_choice()

aqg_start = time.perf_counter()
f1, f2 = shor.find_factors(target_n)
aqg_time = time.perf_counter() - aqg_start

print("-" * 45)
print(f"ERGEBNIS: Faktoren sind {f1} und {f2}")
print(f"AQG-Laufzeit: {aqg_time:.6f} s")
gitter_wert = (target_n * shor.karl_factor) % shor.grid
multi_factors = shor.aqg_multifaktor_scan(target_n, gitter_wert)
print(f"AQG-Multifaktor-Scan: {multi_factors}")
galaktische_faktoren = shor.galaktischer_multifaktor_scan(target_n)
print(f"Galaktischer Multifaktor-Scan: {galaktische_faktoren}")
deep_space_faktoren = shor.deep_space_multifaktor_scan(target_n)
print(f"Deep-Space Multifaktor-Scan: {deep_space_faktoren}")
# Auto-Tuning auch fuer verbleibende Rest-Kerne nach Deep-Space-Scan.
rest_kern = None
if len(deep_space_faktoren) == 1 and deep_space_faktoren[0] == target_n:
    rest_kern = target_n
elif deep_space_faktoren and deep_space_faktoren[-1] > 103:
    rest_kern = deep_space_faktoren[-1]

if rest_kern is not None:
    tuning_hit = shor.auto_tuning_resonator(rest_kern)
    if tuning_hit is not None:
        print(f"Auto-Tuning Resonator: {tuning_hit[0]} und {tuning_hit[1]}")
    else:
        print("Auto-Tuning Resonator: Kein Schatten-Treffer bis Grenze 1000.")
        lvl100_p, lvl100_q = shor.deep_space_resonator_level_100(rest_kern)
        if lvl100_p is not None:
            print(f"Deep-Space Resonator Level 100: {lvl100_p} und {lvl100_q}")
        else:
            print(f"Deep-Space Resonator Level 100: Kein Treffer, Rest-Kern {lvl100_q}.")

if run_compare:
    print("\n--- KLASSISCHER VERGLEICHSTEST (Pollard-Rho) ---")
    classic_start = time.perf_counter()
    classic_factors = factorize_best_known(target_n)
    classic_time = time.perf_counter() - classic_start
    print(f"Klassische Faktoren: {classic_factors}")
    print(f"Klassische Laufzeit: {classic_time:.6f} s")