import numpy as np
import pandas as pd

# Define first 10 Carmichael numbers (Composite Almost Primes)
carmichael_numbers = [561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841, 29341]

# Load Riemann zeros
zeros = np.load('zeros6.npy')
gamma = zeros[:100000] # Use a significant subset for high resolution

def zeta_signal(x, gamma_list):
    ln_x = np.log(x)
    return np.sum(np.cos(gamma_list * ln_x))

def find_one_representation(n):
    limit = int(np.sqrt(n)) + 1
    for a in range(limit):
        for b in range(limit):
            if a*a + b*b > n: break
            for c in range(limit):
                rem = n - (a*a + b*b + c*c)
                if rem < 0: break
                d = int(np.sqrt(rem))
                if a*a + b*b + c*c + d*d == n:
                    return (a, b, c, d)
    return (0,0,0,0)

results = []
for n in carmichael_numbers:
    sig = zeta_signal(n, gamma)
    coords = find_one_representation(n)
    results.append({
        'n': n,
        'Zeta_Resonance': sig,
        'a': coords[0],
        'b': coords[1],
        'c': coords[2],
        'd': coords[3]
    })

df_carmichael = pd.DataFrame(results)

# Compare with some true primes of similar magnitude from the Vierlinge file
# (Already have some data from previous turns, but let's take a quick look)
# Prime 821 had resonance ~ -1481 (with 50k zeros). 
# Let's re-calculate a prime resonance for comparison with same 100k zeros.
prime_ref = 1481
prime_sig = zeta_signal(prime_ref, gamma)

df_carmichael.to_csv('carmichael_zeta_analysis.csv', index=False)
print(df_carmichael)
print(f"\nReference Prime {prime_ref} Resonance: {prime_sig}")