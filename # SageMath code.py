# SageMath code
# Docs (als Referenz):
# https://doc.sagemath.org/html/en/reference/calculus/sage/calculus/functions.html
# https://doc.sagemath.org/html/en/reference/functions/sage/functions/zeta.html
# https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/quatalg/quaternion_algebra.html

from sage.all import *
import numpy as np

############################
# 1) Quaternionen-Basis & E/A/B/C Map (mod 12)
############################

# Hamilton-Quaternionen über RR:
H = QuaternionAlgebra(RR, -1, -1)
i, j, k = H.gens()
one = H(1)

def class_basis_mod12(p):
    """
    Mappe Primzahlen >3 nach mod12-Klassen:
    1 -> 1, 5 -> i, 7 -> j, 11 -> k
    (p=2,3 behandeln wir separat)
    """
    r = Integer(p) % 12
    if r == 1:
        return one
    elif r == 5:
        return i
    elif r == 7:
        return j
    elif r == 11:
        return k
    else:
        # Für p=2,3 oder unerwartete Fälle:
        return one

def prime_quaternion_field(s, P=5000, use_log=True, include_2_3=False):
    """
    Q(s) = sum_{p<=P} w(p,s) * u(p)
    mit w(p,s) = (log p) * p^{-s} oder p^{-s}
    """
    sigma = RR(s.real())
    t = RR(s.imag())
    Qs = H(0)

    for p in prime_range(2, P+1):
        if (p in [2,3]) and (not include_2_3):
            continue
        u = class_basis_mod12(p)
        # Gewicht: log(p)*p^{-s} (komplex), als Paar (Re,Im) auf Quaternionen-Realteil abbilden?
        # Wir nehmen Betrag/Phase sauber: w = log(p)*exp(-(sigma+it)*log p)
        lp = RR(log(p))
        if use_log:
            amp = lp
        else:
            amp = RR(1)

        # complex weight:
        # p^{-s} = exp(-s*log p) = exp(-sigma*log p) * exp(-i*t*log p)
        mag = exp(-sigma*lp)
        phase = -t*lp
        w_re = amp * mag * cos(phase)
        w_im = amp * mag * sin(phase)

        # Einfache, transparente Einbettung:
        # Realteil der Gewichtung skaliert u(p)
        # Imaginärteil skaliert ebenfalls u(p) aber auf eine "orthogonale" Achse – hier: multipliziere mit i als Phasen-Träger
        # (Das ist eine Design-Entscheidung; du kannst die Einbettung später an dein Bamberg-Schema anpassen.)
        Qs += (w_re * u) + (w_im * (i*u))

    return Qs

def qnorm(Q):
    # Quaternionen-Norm (Euklidische Norm der Komponenten)
    return sqrt(Q.reduced_norm())

############################
# 2) Xi(s)-Landschaft (completed zeta)
############################

# Xi-Funktion: xi(s) = 1/2 * s*(s-1) * pi^{-s/2} * Gamma(s/2) * zeta(s)
def xi(s):
    return (RR(0.5) * s*(s-1) * (pi**(-s/2)) * gamma(s/2) * zeta(s))

############################
# 3) Scan/Heatmap Helfer
############################

def scan_grid(func, sigmas, ts):
    """
    func(sigma, t) -> float
    """
    A = np.zeros((len(ts), len(sigmas)), dtype=float)
    for it, t in enumerate(ts):
        for isg, sg in enumerate(sigmas):
            A[it, isg] = float(func(sg, t))
    return A

############################
# 4) Experimente
############################

# Parameter (zuerst klein lassen)
P = 5000
sigmas = np.linspace(0.2, 0.8, 121)
ts     = np.linspace(0, 60, 241)

# (a) Xi-Landschaft: log|xi(s)|
def xi_logabs(sg, t):
    s = CC(sg, t)
    val = xi(s)
    return float(log(abs(val) + 1e-30))

Xi = scan_grid(xi_logabs, sigmas, ts)

# (b) Quaternionenfeld-Norm: log ||Q(s)||
def Q_logabs(sg, t):
    s = CC(sg, t)
    Qs = prime_quaternion_field(s, P=P, use_log=True, include_2_3=False)
    return float(log(qnorm(Qs) + 1e-30))

Qn = scan_grid(Q_logabs, sigmas, ts)

############################
# 5) Visualisierung
############################

import matplotlib.pyplot as plt

# Xi heatmap
plt.figure()
plt.imshow(Xi, aspect='auto', origin='lower',
           extent=[sigmas[0], sigmas[-1], ts[0], ts[-1]])
plt.xlabel('sigma')
plt.ylabel('t')
plt.title('log |xi(sigma + i t)|')
plt.colorbar()
plt.show()

# Quaternion norm heatmap
plt.figure()
plt.imshow(Qn, aspect='auto', origin='lower',
           extent=[sigmas[0], sigmas[-1], ts[0], ts[-1]])
plt.xlabel('sigma')
plt.ylabel('t')
plt.title('log ||Q(sigma + i t)||  (P=%d)' % P)
plt.colorbar()
plt.show()

############################
# 6) 1D-Schnitt: Fix t und vergleiche Sigma-Profil
############################

t0 = 14.134725  # Nähe der ersten nichttrivialen Nullstelle (nur als Orientierung)
vals_sigma = []
for sg in sigmas:
    vals_sigma.append(Q_logabs(sg, t0))

plt.figure()
plt.plot(sigmas, vals_sigma)
plt.xlabel('sigma')
plt.ylabel('log ||Q(sigma + i t0)||')
plt.title('Sigma-Profil der Quaternionen-"Energie" bei t0=%.6f' % t0)
plt.show()