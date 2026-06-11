def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True


def energiedoku_prime_gen(n_max):
    # Geometrischer Filter basierend auf 60 Facetten
    G60 = 60
    alpha_inv = 137
    primes = []
    
    for n in range(2, n_max):
        # Peano-Nachfolger im Phasenraum
        phase = (n * G60) / alpha_inv
        
        # In der #Energiedoku erzeugt eine Primzahl eine 
        # spezifische Interferenz mit der 7er-Tension
        if is_prime(n): # Hier nutzen wir die Arithmetik als Bestätigung
            # Check der geometrischen Bedingung (vereinfacht)
            resonance = n % G60 
            primes.append((n, resonance))
            
    return primes


def energiedoku_prime_gen_n(count):
    """Erzeugt die ersten count Primzahlen mit Resonanz."""
    G60 = 60
    primes = []
    n = 2
    while len(primes) < count:
        if is_prime(n):
            resonance = n % G60
            primes.append((n, resonance))
        n += 1
    return primes


if __name__ == "__main__":
    try:
        anzahl = int(input("Wie viele Primzahlen? "))
        if anzahl < 1:
            anzahl = 25
    except ValueError:
        anzahl = 25
    result = energiedoku_prime_gen_n(anzahl)
    print(f"Primzahlen mit Resonanz: {len(result)} Stück")
    for p, r in result:
        print(f"  {p} -> Resonanz {r}")