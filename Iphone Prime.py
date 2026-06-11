import json
import os
import random

# Hochpräziser Miller-Rabin Test für industrielle Sicherheit
def is_prime_miller_rabin(n, k=40):
    if n <= 3: return n == 2 or n == 3
    if n % 2 == 0: return False
    r, d = 0, n - 1
    while d % 2 == 0:
        r, d = r + 1, d // 2
    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1: break
        else: return False
    return True

class MorleyGenerator:
    def __init__(self, state_file="morley_state.json"):
        self.state_file = state_file
        self.load_state()

    def load_state(self):
        if os.path.exists(self.state_file):
            with open(self.state_file, 'r') as f:
                self.state = json.load(f)
        else:
            # Initialer Warp-Punkt (Beispiel n=1200 / 10^90)
            self.state = {"n_layer": 1200, "last_p": 0, "count": 0}

    def save_state(self):
        with open(self.state_file, 'w') as f:
            json.dump(self.state, f)

    def get_next_prime(self, start_p=None):
        p = start_p if start_p else self.state["last_p"]
        if p % 2 == 0: p += 1
        
        while True:
            p += 2
            if is_prime_miller_rabin(p):
                self.state["last_p"] = p
                self.state["count"] += 1
                self.save_state()
                return p

if __name__ == "__main__":
    # Anwendung
    gen = MorleyGenerator()
    # Optionaler Startpunkt aus Umgebung, sonst aus dem Zustand/fallback.
    env_start = os.getenv("MORLEY_START_P")
    if env_start is not None and env_start.isdigit():
        start_p = int(env_start)
    else:
        start_p = gen.state["last_p"] if gen.state["last_p"] > 0 else 3445331583263037223893915737290561921038

    next_secure_prime = gen.get_next_prime(start_p=start_p)
    print(f"Generierte Primzahl #{gen.state['count']}: {next_secure_prime}")