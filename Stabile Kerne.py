import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig, eigvals

# ============================================================
# REPRODUZIERBARKEIT
# ============================================================

np.random.seed(42)

# ============================================================
# PARAMETER
# ============================================================

N = 10_000_000

MODULUS = 60

# primtragende Klassen mod 60
states = [
    1,7,11,13,
    17,19,23,29,
    31,37,41,43,
    47,49,53,59
]

idx = {
    s:i for i,s in enumerate(states)
}

# ============================================================
# PRIMZAHLEN
# ============================================================

def primes_upto(n):

    sieve = np.ones(n + 1, dtype=bool)

    sieve[:2] = False

    for p in range(
        2,
        int(np.sqrt(n)) + 1
    ):

        if sieve[p]:
            sieve[p*p:n+1:p] = False

    return np.nonzero(sieve)[0]


print("Erzeuge Primzahlen ...")

primes = primes_upto(N)

print("Primzahlen:", len(primes))

# ============================================================
# ZUSTANDSSEQUENZ
# ============================================================

seq = []

for p in primes:

    r = p % MODULUS

    if r in states:
        seq.append(r)

seq = np.array(seq)

print("Zustände:", len(seq))

# ============================================================
# ÜBERGANGSMATRIX
# ============================================================

def transition_matrix(sequence):

    n = len(states)

    M = np.zeros((n,n))

    for a,b in zip(sequence[:-1], sequence[1:]):

        i = idx[a]
        j = idx[b]

        M[i,j] += 1

    for i in range(n):

        s = M[i].sum()

        if s > 0:
            M[i] /= s

    return M


M = transition_matrix(seq)

# ============================================================
# EIGENWERTE
# ============================================================

evals = eigvals(M)

evals_sorted = sorted(
    evals,
    key=lambda z: abs(z),
    reverse=True
)

print("\n=== Dominante Eigenwerte ===")

for z in evals_sorted[:10]:
    print(z)

# ============================================================
# STATIONÄRE VERTEILUNG
# ============================================================

vals, vecs = eig(M.T)

# Eigenwert nahe 1
k = np.argmin(
    np.abs(vals - 1)
)

pi = np.real(
    vecs[:,k]
)

pi = pi / np.sum(pi)

print("\n=== Stationäre Verteilung ===")

for s,p in zip(states, pi):

    print(
        f"{s:2d} : {p:.6f}"
    )

# ============================================================
# DETAILED BALANCE
# ============================================================

balance = np.zeros_like(M)

for i in range(len(states)):

    for j in range(len(states)):

        balance[i,j] = (
            pi[i] * M[i,j]
            -
            pi[j] * M[j,i]
        )

# ============================================================
# GERICHTETE STRÖME
# ============================================================

J = M - M.T

# antisymmetrische Stärke
flow_strength = np.sum(
    np.abs(J)
)

print("\n=== Flussstärke ===")

print(flow_strength)

# ============================================================
# DETAILLIERTE BALANCE
# ============================================================

balance_violation = np.sum(
    np.abs(balance)
)

print("\n=== Verletzung Detailed Balance ===")

print(balance_violation)

# ============================================================
# NETTO-FLÜSSE
# ============================================================

flows = []

for i in range(len(states)):

    for j in range(i+1, len(states)):

        f = J[i,j]

        flows.append(
            (
                abs(f),
                states[i],
                states[j],
                f
            )
        )

flows.sort(reverse=True)

print("\n=== Größte gerichtete Flüsse ===")

for x in flows[:20]:

    mag,a,b,f = x

    print(
        f"{a:2d} -> {b:2d} : {f:+.6f}"
    )

# ============================================================
# SPEKTRALRADIUS
# ============================================================

absvals = np.abs(evals_sorted)

rho = absvals[0]

lambda2 = absvals[1]

tau = -1 / np.log(lambda2)

print("\n=== Dynamik ===")

print("rho(M)     =", rho)
print("|lambda2|  =", lambda2)
print("tau        =", tau)

# ============================================================
# ENTROPIE
# ============================================================

eps = 1e-15

H = -np.sum(
    M * np.log(M + eps)
)

print("\n=== Entropie ===")

print(H)

# ============================================================
# FOURIER DER EIGENWERTE
# ============================================================

angles = np.angle(evals)

radii = np.abs(evals)

# ============================================================
# PLOT 1
# ============================================================

plt.figure(figsize=(8,8))

plt.scatter(
    np.real(evals),
    np.imag(evals),
    s=80
)

circle = plt.Circle(
    (0,0),
    1,
    fill=False
)

plt.gca().add_artist(circle)

plt.axhline(0,color='black')
plt.axvline(0,color='black')

plt.title(
    "Eigenwerte im Komplexen"
)

plt.grid()

plt.axis("equal")

plt.savefig(
    "mod60_eigenkreis.png",
    dpi=300
)

# ============================================================
# PLOT 2
# ============================================================

plt.figure(figsize=(10,8))

plt.imshow(
    J,
    cmap="seismic"
)

plt.colorbar()

plt.xticks(
    range(len(states)),
    states,
    rotation=90
)

plt.yticks(
    range(len(states)),
    states
)

plt.title(
    "Gerichtete Flüsse J_ij"
)

plt.savefig(
    "mod60_fluesse.png",
    dpi=300
)

# ============================================================
# PLOT 3
# ============================================================

plt.figure(figsize=(10,8))

plt.imshow(
    balance,
    cmap="coolwarm"
)

plt.colorbar()

plt.xticks(
    range(len(states)),
    states,
    rotation=90
)

plt.yticks(
    range(len(states)),
    states
)

plt.title(
    "Detailed-Balance-Verletzung"
)

plt.savefig(
    "mod60_balance.png",
    dpi=300
)

# ============================================================
# PLOT 4
# ============================================================

plt.figure(figsize=(12,5))

plt.bar(
    range(len(states)),
    pi
)

plt.xticks(
    range(len(states)),
    states
)

plt.title(
    "Stationäre Verteilung"
)

plt.grid()

plt.savefig(
    "mod60_stationaer.png",
    dpi=300
)

# ============================================================
# PLOT 5
# ============================================================

plt.figure(figsize=(12,5))

plt.bar(
    range(len(absvals)),
    absvals
)

plt.title(
    "Eigenwertbeträge"
)

plt.grid()

plt.savefig(
    "mod60_spektrum.png",
    dpi=300
)

print("\nPlots gespeichert.")
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigvals
from scipy.fft import fft, fftfreq

# ============================================================
# REPRODUZIERBARKEIT
# ============================================================

np.random.seed(42)

# ============================================================
# PARAMETER
# ============================================================

N = 10_000_000

WINDOW = 5000
STEP = 1000

SIMULATIONS = 25

MAX_LAG = 250

# ============================================================
# MODULORAUM
# ============================================================

MODULUS = 60

# primtragende Klassen mod 60
states = [
    1, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41,
    43, 47, 49, 53, 59
]

idx = {
    s:i for i,s in enumerate(states)
}

# ============================================================
# PRIMZAHLEN
# ============================================================

def primes_upto(n):

    sieve = np.ones(n + 1, dtype=bool)

    sieve[:2] = False

    for p in range(
        2,
        int(np.sqrt(n)) + 1
    ):

        if sieve[p]:
            sieve[p*p:n+1:p] = False

    return np.nonzero(sieve)[0]


print("Erzeuge Primzahlen ...")

primes = primes_upto(N)

print("Primzahlen:", len(primes))

# ============================================================
# ZUSTANDSSEQUENZ
# ============================================================

seq = []

for p in primes:

    r = p % MODULUS

    if r in states:
        seq.append(r)

seq = np.array(seq)

print("Zustände:", len(seq))

# ============================================================
# ÜBERGANGSMATRIX
# ============================================================

def transition_matrix(sequence):

    n = len(states)

    M = np.zeros((n,n))

    for a,b in zip(sequence[:-1], sequence[1:]):

        i = idx[a]
        j = idx[b]

        M[i,j] += 1

    for i in range(n):

        s = M[i].sum()

        if s > 0:
            M[i] /= s

    return M


M_real = transition_matrix(seq)

# ============================================================
# EIGENWERTE
# ============================================================

eig_real = eigvals(M_real)

# nach Betrag sortieren
eig_sorted = sorted(
    eig_real,
    key=lambda z: abs(z),
    reverse=True
)

print("\n=== Dominante Eigenwerte ===")

for x in eig_sorted[:10]:
    print(x)

# ============================================================
# SPEKTRALRADIUS
# ============================================================

absvals = np.abs(eig_sorted)

rho = absvals[0]

lambda2 = absvals[1]

print("\n=== Spektralgrößen ===")

print("rho(M) =", rho)
print("|lambda2| =", lambda2)

# ============================================================
# RELAXATIONSZEIT
# ============================================================

tau = -1.0 / np.log(lambda2)

print("Relaxationszeit tau =", tau)

# ============================================================
# ENTROPIE
# ============================================================

def entropy_matrix(M):

    eps = 1e-15

    return -np.sum(
        M * np.log(M + eps)
    )


H_real = entropy_matrix(M_real)

print("\n=== Entropie ===")
print(H_real)

# ============================================================
# RESONANZFUNKTION
# ============================================================

uniform = 1 / len(states)

def resonance(window):

    M = transition_matrix(window)

    return np.sum(
        (M - uniform)**2
    )

# ============================================================
# RESONANZFELD
# ============================================================

R_real = []

print("\nBerechne Resonanzfeld ...")

for start in range(
    0,
    len(seq)-WINDOW,
    STEP
):

    w = seq[start:start+WINDOW]

    R_real.append(
        resonance(w)
    )

R_real = np.array(R_real)

print("\n=== Resonanz ===")

print("min  :", R_real.min())
print("max  :", R_real.max())
print("mean :", R_real.mean())
print("std  :", R_real.std())

# ============================================================
# AUTOKORRELATION
# ============================================================

num_seq = np.array([
    idx[x]
    for x in seq
])

num_seq = (
    num_seq - np.mean(num_seq)
) / np.std(num_seq)

corr = []

print("\nBerechne Autokorrelation ...")

for lag in range(1, MAX_LAG):

    c = np.mean(
        num_seq[:-lag] *
        num_seq[lag:]
    )

    corr.append(c)

corr = np.array(corr)

print("\n=== Autokorrelation ===")

print("max :", corr.max())
print("min :", corr.min())

# ============================================================
# FOURIER
# ============================================================

Y = fft(R_real)

freq = fftfreq(len(R_real))

# ============================================================
# HL-SURROGATE
# ============================================================

print("\nErzeuge HL-Surrogate ...")

gaps = np.diff(primes)

def hl_surrogate(start, gaps, limit):

    shuffled = np.random.permutation(gaps)

    out = [start]

    s = start

    for g in shuffled:

        s += g

        if s > limit:
            break

        out.append(s)

    return np.array(out)

# ============================================================
# SPEICHER
# ============================================================

hl_res = []
hl_ent = []
hl_corr = []
hl_lambda2 = []

# ============================================================
# SIMULATIONEN
# ============================================================

for sim in range(SIMULATIONS):

    print("Simulation", sim+1)

    hp = hl_surrogate(
        11,
        gaps,
        N
    )

    hseq = []

    for x in hp:

        r = x % MODULUS

        if r in states:
            hseq.append(r)

    hseq = np.array(hseq)

    if len(hseq) < WINDOW:
        continue

    # Matrix
    Mh = transition_matrix(hseq)

    evh = eigvals(Mh)

    evh = sorted(
        evh,
        key=lambda z: abs(z),
        reverse=True
    )

    hl_lambda2.append(
        abs(evh[1])
    )

    # Entropie
    hl_ent.append(
        entropy_matrix(Mh)
    )

    # Resonanz
    Rh = []

    for start in range(
        0,
        len(hseq)-WINDOW,
        STEP
    ):

        w = hseq[start:start+WINDOW]

        Rh.append(
            resonance(w)
        )

    Rh = np.array(Rh)

    hl_res.append(
        Rh.mean()
    )

    # Autokorrelation
    hnum = np.array([
        idx[x]
        for x in hseq
    ])

    hnum = (
        hnum - np.mean(hnum)
    ) / np.std(hnum)

    hc = []

    for lag in range(1, MAX_LAG):

        c = np.mean(
            hnum[:-lag] *
            hnum[lag:]
        )

        hc.append(c)

    hc = np.array(hc)

    hl_corr.append(
        np.max(np.abs(hc))
    )

hl_res = np.array(hl_res)
hl_ent = np.array(hl_ent)
hl_corr = np.array(hl_corr)
hl_lambda2 = np.array(hl_lambda2)

# ============================================================
# Z-SCORES
# ============================================================

z_res = (
    R_real.mean()
    -
    hl_res.mean()
) / hl_res.std()

z_ent = (
    H_real
    -
    hl_ent.mean()
) / hl_ent.std()

z_corr = (
    np.max(np.abs(corr))
    -
    hl_corr.mean()
) / hl_corr.std()

z_lambda = (
    lambda2
    -
    hl_lambda2.mean()
) / hl_lambda2.std()

print("\n=== HL-VERGLEICH ===")

print("HL Resonanz :", hl_res.mean())
print("HL Entropie :", hl_ent.mean())
print("HL Corr     :", hl_corr.mean())
print("HL lambda2  :", hl_lambda2.mean())

print("\n=== Z-SCORES ===")

print("Resonanz :", z_res)
print("Entropie :", z_ent)
print("Korrelation :", z_corr)
print("lambda2 :", z_lambda)

# ============================================================
# PLOT 1
# ============================================================

plt.figure(figsize=(12,5))

plt.plot(R_real)

plt.title(
    f"Resonanzfeld mod {MODULUS}"
)

plt.grid()

plt.savefig(
    f"mod{MODULUS}_resonanz.png",
    dpi=300
)

# ============================================================
# PLOT 2
# ============================================================

plt.figure(figsize=(12,5))

plt.plot(
    range(1,MAX_LAG),
    corr
)

plt.title(
    f"Autokorrelation mod {MODULUS}"
)

plt.xlabel("Lag")
plt.ylabel("C(k)")

plt.grid()

plt.savefig(
    f"mod{MODULUS}_autokorrelation.png",
    dpi=300
)

# ============================================================
# PLOT 3
# ============================================================

plt.figure(figsize=(8,5))

plt.hist(
    hl_res,
    bins=10,
    alpha=0.7
)

plt.axvline(
    R_real.mean(),
    linewidth=3
)

plt.title(
    f"HL Vergleich mod {MODULUS}"
)

plt.grid()

plt.savefig(
    f"mod{MODULUS}_hl_resonanz.png",
    dpi=300
)

# ============================================================
# PLOT 4
# ============================================================

plt.figure(figsize=(7,7))

plt.imshow(
    M_real,
    cmap="viridis"
)

plt.colorbar()

plt.xticks(
    range(len(states)),
    states,
    rotation=90
)

plt.yticks(
    range(len(states)),
    states
)

plt.title(
    f"Übergangsmatrix mod {MODULUS}"
)

plt.savefig(
    f"mod{MODULUS}_matrix.png",
    dpi=300
)

# ============================================================
# PLOT 5
# ============================================================

vals_real = np.abs(eig_real)

plt.figure(figsize=(10,5))

plt.scatter(
    range(len(vals_real)),
    sorted(vals_real, reverse=True),
    s=100
)

plt.title(
    f"Eigenwertspektrum mod {MODULUS}"
)

plt.grid()

plt.savefig(
    f"mod{MODULUS}_eigenwerte.png",
    dpi=300
)

print("\nPlots gespeichert.")
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigvals
from scipy.fft import fft, fftfreq
from collections import defaultdict

# ============================================================
# REPRODUZIERBARKEIT
# ============================================================

np.random.seed(42)

# ============================================================
# PARAMETER
# ============================================================

N = 10_000_000

WINDOW = 5000
STEP = 1000

SIMULATIONS = 25

MAX_LAG = 200

# ============================================================
# PRIMZAHLEN
# ============================================================

def primes_upto(n):

    sieve = np.ones(n + 1, dtype=bool)

    sieve[:2] = False

    for p in range(
        2,
        int(np.sqrt(n)) + 1
    ):

        if sieve[p]:
            sieve[p*p:n+1:p] = False

    return np.nonzero(sieve)[0]


print("Erzeuge Primzahlen ...")

primes = primes_upto(N)

print("Primzahlen:", len(primes))

# ============================================================
# MOD-36-ZUSTÄNDE
# ============================================================

states36 = [
    1,5,7,11,
    13,17,19,23,
    25,29,31,35
]

idx36 = {
    s:i for i,s in enumerate(states36)
}

# ============================================================
# SEQUENZ
# ============================================================

seq36 = []

for p in primes:

    r = p % 36

    if r in states36:
        seq36.append(r)

seq36 = np.array(seq36)

print("Zustände:", len(seq36))

# ============================================================
# ÜBERGANGSMATRIX
# ============================================================

def transition_matrix(seq):

    n = len(states36)

    M = np.zeros((n,n))

    for a,b in zip(seq[:-1], seq[1:]):

        i = idx36[a]
        j = idx36[b]

        M[i,j] += 1

    for i in range(n):

        s = M[i].sum()

        if s > 0:
            M[i] /= s

    return M


M_real = transition_matrix(seq36)

# ============================================================
# EIGENWERTE
# ============================================================

eig_real = eigvals(M_real)

print("\n=== Eigenwerte ===")

for x in eig_real:
    print(x)

# ============================================================
# ENTROPIE
# ============================================================

def entropy_matrix(M):

    eps = 1e-15

    return -np.sum(
        M * np.log(M + eps)
    )


H_real = entropy_matrix(M_real)

print("\n=== Entropie ===")
print(H_real)

# ============================================================
# RESONANZFUNKTION
# ============================================================

uniform = 1 / len(states36)

def resonance(window):

    M = transition_matrix(window)

    return np.sum(
        (M - uniform)**2
    )

# ============================================================
# RESONANZFELD
# ============================================================

R_real = []

print("\nBerechne Resonanzfeld ...")

for start in range(
    0,
    len(seq36)-WINDOW,
    STEP
):

    w = seq36[start:start+WINDOW]

    R_real.append(
        resonance(w)
    )

R_real = np.array(R_real)

print("\n=== Resonanz ===")

print("min  :", R_real.min())
print("max  :", R_real.max())
print("mean :", R_real.mean())
print("std  :", R_real.std())

# ============================================================
# AUTOKORRELATION
# ============================================================

# numerische Projektion
num_seq = np.array([
    idx36[x]
    for x in seq36
])

num_seq = (
    num_seq - np.mean(num_seq)
) / np.std(num_seq)

corr = []

print("\nBerechne Autokorrelation ...")

for lag in range(1, MAX_LAG):

    c = np.mean(
        num_seq[:-lag] *
        num_seq[lag:]
    )

    corr.append(c)

corr = np.array(corr)

print("\nAutokorrelation:")
print("max :", corr.max())
print("min :", corr.min())

# ============================================================
# FOURIER
# ============================================================

Y = fft(R_real)

freq = fftfreq(len(R_real))

# ============================================================
# HARDY-LITTLEWOOD-SURROGATE
# ============================================================

print("\nErzeuge HL-Surrogate ...")

gaps = np.diff(primes)

def hl_surrogate(start, gaps, limit):

    shuffled = np.random.permutation(gaps)

    out = [start]

    s = start

    for g in shuffled:

        s += g

        if s > limit:
            break

        out.append(s)

    return np.array(out)

# ============================================================
# HL-ERGEBNISSE
# ============================================================

hl_means = []

hl_entropy = []

hl_corrmax = []

hl_eigs = []

for sim in range(SIMULATIONS):

    print("Simulation", sim+1)

    hp = hl_surrogate(
        11,
        gaps,
        N
    )

    hseq = []

    for x in hp:

        r = x % 36

        if r in states36:
            hseq.append(r)

    hseq = np.array(hseq)

    if len(hseq) < WINDOW:
        continue

    # Matrix
    Mh = transition_matrix(hseq)

    evh = eigvals(Mh)

    hl_eigs.append(evh)

    # Entropie
    hl_entropy.append(
        entropy_matrix(Mh)
    )

    # Resonanz
    Rh = []

    for start in range(
        0,
        len(hseq)-WINDOW,
        STEP
    ):

        w = hseq[start:start+WINDOW]

        Rh.append(
            resonance(w)
        )

    Rh = np.array(Rh)

    hl_means.append(
        Rh.mean()
    )

    # Autokorrelation
    hnum = np.array([
        idx36[x]
        for x in hseq
    ])

    hnum = (
        hnum - np.mean(hnum)
    ) / np.std(hnum)

    hc = []

    for lag in range(1, MAX_LAG):

        c = np.mean(
            hnum[:-lag] *
            hnum[lag:]
        )

        hc.append(c)

    hc = np.array(hc)

    hl_corrmax.append(
        np.max(np.abs(hc))
    )

hl_means = np.array(hl_means)
hl_entropy = np.array(hl_entropy)
hl_corrmax = np.array(hl_corrmax)

# ============================================================
# Z-SCORES
# ============================================================

z_res = (
    R_real.mean()
    -
    hl_means.mean()
) / hl_means.std()

z_ent = (
    H_real
    -
    hl_entropy.mean()
) / hl_entropy.std()

z_corr = (
    np.max(np.abs(corr))
    -
    hl_corrmax.mean()
) / hl_corrmax.std()

print("\n=== HL-VERGLEICH ===")

print("HL Resonanz mean :", hl_means.mean())
print("HL Entropie mean :", hl_entropy.mean())
print("HL Corr mean     :", hl_corrmax.mean())

print("\n=== Z-SCORES ===")

print("Resonanz :", z_res)
print("Entropie :", z_ent)
print("Korrelation :", z_corr)

# ============================================================
# PLOT 1
# ============================================================

plt.figure(figsize=(12,5))

plt.plot(R_real)

plt.title(
    "mod36 Resonanzfeld"
)

plt.grid()

plt.savefig(
    "mod36_resonanz.png",
    dpi=300
)

# ============================================================
# PLOT 2
# ============================================================

plt.figure(figsize=(12,5))

plt.plot(
    range(1,MAX_LAG),
    corr
)

plt.title(
    "Autokorrelation"
)

plt.xlabel("Lag")
plt.ylabel("C(k)")

plt.grid()

plt.savefig(
    "mod36_autokorrelation.png",
    dpi=300
)

# ============================================================
# PLOT 3
# ============================================================

plt.figure(figsize=(8,5))

plt.hist(
    hl_means,
    bins=10,
    alpha=0.7
)

plt.axvline(
    R_real.mean(),
    linewidth=3
)

plt.title(
    "Resonanz vs HL"
)

plt.grid()

plt.savefig(
    "mod36_hl_resonanz.png",
    dpi=300
)

# ============================================================
# PLOT 4
# ============================================================

plt.figure(figsize=(6,6))

plt.imshow(
    M_real,
    cmap="viridis"
)

plt.colorbar()

plt.xticks(
    range(len(states36)),
    states36,
    rotation=90
)

plt.yticks(
    range(len(states36)),
    states36
)

plt.title(
    "Übergangsmatrix mod36"
)

plt.savefig(
    "mod36_matrix.png",
    dpi=300
)

# ============================================================
# PLOT 5
# ============================================================

vals_real = np.abs(eig_real)

vals_hl = np.abs(
    hl_eigs[0]
)

x = np.arange(len(vals_real))

plt.figure(figsize=(10,5))

plt.scatter(
    x,
    vals_real,
    s=120,
    label="Real"
)

plt.scatter(
    x,
    vals_hl,
    s=120,
    label="HL"
)

plt.legend()

plt.title(
    "Eigenwertvergleich"
)

plt.grid()

plt.savefig(
    "mod36_eigenwerte.png",
    dpi=300
)

print("\nPlots gespeichert.")
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigvals
from collections import defaultdict

# ============================================================
# REPRODUZIERBARKEIT
# ============================================================

np.random.seed(42)

# ============================================================
# PARAMETER
# ============================================================

N = 10_000_000

WINDOW = 4000
STEP = 1000

SIMULATIONS = 25

# ============================================================
# MOD-210-ZUSTÄNDE
# ============================================================

states210 = [5, 11, 101, 191]

idx210 = {
    s:i for i,s in enumerate(states210)
}


# ============================================================
# PRIMZAHLEN
# ============================================================

def primes_upto(n):

    sieve = np.ones(n + 1, dtype=bool)

    sieve[:2] = False

    for p in range(
        2,
        int(np.sqrt(n)) + 1
    ):

        if sieve[p]:
            sieve[p*p:n+1:p] = False

    return np.nonzero(sieve)[0]


print("Erzeuge Primzahlen ...")

primes = primes_upto(N)

prime_set = set(primes)

print("Primzahlen:", len(primes))

# ============================================================
# VIERLINGE
# ============================================================

quadruplets = []

for p in primes:

    if (
        p+2 in prime_set and
        p+6 in prime_set and
        p+8 in prime_set
    ):

        quadruplets.append(
            (p,p+2,p+6,p+8)
        )

print("\nVierlinge:", len(quadruplets))

# ============================================================
# MOD-210-SEQUENZ
# ============================================================

seq210 = []

for q in quadruplets:

    r = q[0] % 210

    if r in states210:
        seq210.append(r)

seq210 = np.array(seq210)

# ============================================================
# HÄUFIGKEITEN
# ============================================================

counts = defaultdict(int)

for x in seq210:
    counts[x] += 1

print("\n=== mod 210 Häufigkeiten ===")

for k in states210:
    print(k, counts[k])

# ============================================================
# ÜBERGANGSMATRIX
# ============================================================

def transition_matrix_210(seq):

    M = np.zeros((4,4))

    for a,b in zip(seq[:-1], seq[1:]):

        i = idx210[a]
        j = idx210[b]

        M[i,j] += 1

    for i in range(4):

        s = M[i].sum()

        if s > 0:
            M[i] /= s

    return M

M_real = transition_matrix_210(seq210)

print("\n=== Übergangsmatrix mod 210 ===")

print(M_real)

# ============================================================
# EIGENWERTE
# ============================================================

eig_real = eigvals(M_real)

print("\n=== Eigenwerte ===")

for x in eig_real:
    print(x)

# ============================================================
# RESONANZFUNKTION
# ============================================================

def resonance_210(window):

    M = transition_matrix_210(window)

    return np.sum(
        (M - 0.25)**2
    )


def resonance_field_210(seq, window, step):

    n = len(seq)

    if n < 2:
        return np.array([])

    # Bei kurzer Sequenz das Fenster automatisch verkleinern.
    eff_window = min(window, n)

    values = []

    for start in range(0, n - eff_window + 1, step):

        w = seq[start:start+eff_window]

        values.append(
            resonance_210(w)
        )

    return np.array(values)


def safe_z_score(observed, samples):

    if len(samples) == 0:
        return np.nan

    sigma = samples.std()

    if sigma == 0:
        return np.nan

    return (observed - samples.mean()) / sigma

# ============================================================
# RESONANZFELD
# ============================================================

R_real = resonance_field_210(seq210, WINDOW, STEP)

print("\n=== Resonanzfeld ===")

if len(R_real) == 0:

    print("Keine Resonanzfenster berechenbar (Sequenz zu kurz).")

else:

    print("min  :", R_real.min())
    print("max  :", R_real.max())
    print("mean :", R_real.mean())
    print("std  :", R_real.std())

# ============================================================
# PERMUTATIONSTEST
# ============================================================

perm_means = []

perm_eigs = []

print("\nPermutationstests ...")

for sim in range(SIMULATIONS):

    perm = np.random.permutation(seq210)

    Mp = transition_matrix_210(perm)

    evp = eigvals(Mp)

    perm_eigs.append(evp)

    Rp = resonance_field_210(perm, WINDOW, STEP)

    if len(Rp) > 0:
        perm_means.append(
            Rp.mean()
        )

perm_means = np.array(perm_means)

print("\n=== Permutation ===")

if len(perm_means) == 0:
    print("Keine gültigen Permutationsfenster.")
else:
    print("mean :", perm_means.mean())
    print("std  :", perm_means.std())

z_perm = safe_z_score(
    R_real.mean() if len(R_real) > 0 else np.nan,
    perm_means
)

print("\nZ gegen Permutation:", z_perm)

# ============================================================
# HARDY-LITTLEWOOD-SURROGAT
# ============================================================

print("\nErzeuge HL-Surrogate ...")

gaps = np.diff(primes)

def hl_surrogate(start, gaps, limit):

    shuffled = np.random.permutation(gaps)

    out = [start]

    s = start

    for g in shuffled:

        s += g

        if s > limit:
            break

        out.append(s)

    return np.array(out)

hl_means = []

hl_eigs = []

for sim in range(SIMULATIONS):

    hp = hl_surrogate(
        11,
        gaps,
        N
    )

    hquads = []

    hset = set(hp)

    for p in hp:

        if (
            p+2 in hset and
            p+6 in hset and
            p+8 in hset
        ):

            hquads.append(
                p % 210
            )

    hquads = np.array([
        x for x in hquads
        if x in states210
    ])

    if len(hquads) < 2:
        continue

    Mh = transition_matrix_210(hquads)

    evh = eigvals(Mh)

    hl_eigs.append(evh)

    Rh = resonance_field_210(hquads, WINDOW, STEP)

    if len(Rh) == 0:
        continue

    hl_means.append(
        Rh.mean()
    )

hl_means = np.array(hl_means)

print("\n=== HL-Surrogat ===")

if len(hl_means) == 0:
    print("Keine gültigen HL-Fenster.")
else:
    print("mean :", hl_means.mean())
    print("std  :", hl_means.std())

z_hl = safe_z_score(
    R_real.mean() if len(R_real) > 0 else np.nan,
    hl_means
)

print("\nZ gegen HL:", z_hl)

# ============================================================
# PLOT 1
# ============================================================

plt.figure(figsize=(12,5))

plt.plot(R_real)

plt.title(
    "mod-210 Resonanzfeld"
)

plt.grid()

plt.savefig(
    "mod210_resonanzfeld.png",
    dpi=300
)

# ============================================================
# PLOT 2
# ============================================================

plt.figure(figsize=(8,5))

if len(perm_means) > 0:
    plt.hist(
        perm_means,
        bins=10,
        alpha=0.7,
        label="Permutation"
    )

if len(R_real) > 0:
    plt.axvline(
        R_real.mean(),
        linewidth=3,
        label="Original"
    )

if len(perm_means) > 0 or len(R_real) > 0:
    plt.legend()

plt.title(
    "Permutationstest mod210"
)

plt.grid()

plt.savefig(
    "mod210_permutation.png",
    dpi=300
)

# ============================================================
# PLOT 3
# ============================================================

plt.figure(figsize=(8,5))

if len(hl_means) > 0:
    plt.hist(
        hl_means,
        bins=10,
        alpha=0.7,
        label="HL"
    )

if len(R_real) > 0:
    plt.axvline(
        R_real.mean(),
        linewidth=3,
        label="Original"
    )

if len(hl_means) > 0 or len(R_real) > 0:
    plt.legend()

plt.title(
    "HL-Test mod210"
)

plt.grid()

plt.savefig(
    "mod210_hl.png",
    dpi=300
)

# ============================================================
# PLOT 4
# ============================================================

plt.figure(figsize=(6,6))

plt.imshow(
    M_real,
    cmap="viridis"
)

plt.colorbar()

plt.xticks(
    range(4),
    states210
)

plt.yticks(
    range(4),
    states210
)

plt.title(
    "Übergangsmatrix mod210"
)

plt.savefig(
    "mod210_matrix.png",
    dpi=300
)

# ============================================================
# PLOT 5
# ============================================================

vals_real = np.abs(eig_real)

vals_hl = np.abs(hl_eigs[0]) if len(hl_eigs) > 0 else np.full(4, np.nan)

x = np.arange(4)

plt.figure(figsize=(8,5))

plt.scatter(
    x,
    vals_real,
    s=120,
    label="Real"
)

plt.scatter(
    x,
    vals_hl,
    s=120,
    label="HL"
)

plt.legend()

plt.title(
    "Eigenwertvergleich mod210"
)

plt.grid()

plt.savefig(
    "mod210_eigenwerte.png",
    dpi=300
)

print("\nPlots gespeichert.")