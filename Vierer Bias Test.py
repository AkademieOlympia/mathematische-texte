import numpy as np
import math
from scipy.stats import binomtest

pr = np.load("eabc_primes_2M.npy")

# Spalte 0: Primzahl, Spalte 1: Familienlabel
primes = pr[:, 0]
labels = pr[:, 1]

# E=0, A=1, B=2, C=3
ABCE = (1, 2, 3, 0)
CEAB = (3, 0, 1, 2)

count_abce = 0
count_ceab = 0
events = []

for i in range(len(labels) - 3):
    quad = tuple(labels[i:i+4])

    if quad == ABCE:
        count_abce += 1
        events.append((i, int(primes[i+3]), "A-B-C-E", count_abce, count_ceab, count_abce - count_ceab))

    elif quad == CEAB:
        count_ceab += 1
        events.append((i, int(primes[i+3]), "C-E-A-B", count_abce, count_ceab, count_abce - count_ceab))

total = count_abce + count_ceab
delta = count_abce - count_ceab
rel_asym = delta / total
log_bias = math.log(count_abce / count_ceab)

print("A-B-C-E:", count_abce)
print("C-E-A-B:", count_ceab)
print("Gesamt:", total)
print("Delta:", delta)
print("Relative Asymmetrie:", rel_asym)
print("Log-Bias:", log_bias)

test = binomtest(count_abce, total, 0.5)
print("p-Wert:", test.pvalue)
print("95%-Konfidenzintervall:", test.proportion_ci())