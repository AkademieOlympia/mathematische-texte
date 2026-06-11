import numpy as np
from collections import Counter, defaultdict
import matplotlib.pyplot as plt

try:
    from scipy.fft import fft
except Exception:
    from numpy.fft import fft


STABLE_420 = [11, 101, 191, 221, 311, 401]


def sieve(n):
    isprime = np.ones(n + 10, dtype=bool)
    isprime[:2] = False
    for i in range(2, int(np.sqrt(n)) + 1):
        if isprime[i]:
            isprime[i * i:n + 10:i] = False
    return isprime


def generate_quadruplets(n):
    isprime = sieve(n + 10)
    primes = np.nonzero(isprime[:n + 1])[0]
    mask = (
        isprime[primes]
        & isprime[primes + 2]
        & isprime[primes + 6]
        & isprime[primes + 8]
    )
    return primes[mask]


def orientation(p):
    if p % 12 == 5:
        return "ABCE"
    if p % 12 == 11:
        return "CEAB"
    return "OTHER"


def analyze_mod420(quads):
    c420 = Counter(int(p % 420) for p in quads)
    orient = defaultdict(Counter)
    for p in quads:
        orient[int(p % 420)][orientation(int(p))] += 1
    return c420, orient


def compute_masses(quads):
    ori = Counter(orientation(int(p)) for p in quads)
    q = len(quads)

    if q == 0:
        return {
            "Q": 0,
            "ABCE": 0,
            "CEAB": 0,
            "m_nu": np.nan,
            "m_s": np.nan,
            "R": np.nan,
            "delta420": np.nan,
        }

    m_nu = 0.45 * abs(ori["ABCE"] - ori["CEAB"]) / q

    c420, _ = analyze_mod420(quads)
    vals = np.array([c420[r] for r in STABLE_420], dtype=float)
    delta420 = (vals.max() - vals.min()) / vals.mean() if vals.mean() > 0 else np.nan
    m_s = 0.45 * delta420 if np.isfinite(delta420) else np.nan
    r_ratio = np.inf if m_nu == 0 else m_s / m_nu

    return {
        "Q": q,
        "ABCE": ori["ABCE"],
        "CEAB": ori["CEAB"],
        "m_nu": m_nu,
        "m_s": m_s,
        "R": r_ratio,
        "delta420": delta420,
    }


def spectral_test(quads):
    seq = np.array([int(p % 420) for p in quads if int(p % 420) in STABLE_420], dtype=int)
    mapping = {11: 0, 101: 1, 191: 2, 221: 3, 311: 4, 401: 5}
    x = np.array([mapping[s] for s in seq], dtype=float)
    if len(x) == 0:
        return np.array([])
    fvals = np.abs(fft(x))
    return fvals[: len(fvals) // 2]


def rolling_r(quads, w=100):
    if len(quads) < w:
        return np.array([])
    vals = []
    for i in range(0, len(quads) - w + 1):
        block = quads[i:i + w]
        res = compute_masses(block)
        vals.append(res["R"])
    return np.array(vals, dtype=float)


if __name__ == "__main__":
    tests = [10**6, 3 * 10**6, 10**7, 3 * 10**7, 10**8]
    all_results = []

    for n in tests:
        print("\n===================================")
        print("N =", n)
        print("===================================")

        quads = generate_quadruplets(n)
        res = compute_masses(quads)

        print("Vierlinge:", res["Q"])
        print("Orientierung:", f"ABCE={res['ABCE']}", f"CEAB={res['CEAB']}")

        c420, orient = analyze_mod420(quads)
        print("\nmod420:")
        for r in STABLE_420:
            print(r, c420[r])

        print("\nmod420 mit Orientierung:")
        for r in STABLE_420:
            print(r, dict(orient[r]))

        print("\nm_nu =", res["m_nu"])
        print("m_s   =", res["m_s"])
        print("R     =", res["R"])
        print("|R-126| =", abs(res["R"] - 126))
        print("|R-137.036| =", abs(res["R"] - 137.036))
        all_results.append(res)

    print("\n===================================")
    print("Spektralanalyse")
    print("===================================")
    quads = generate_quadruplets(10**8)
    fvals = spectral_test(quads)
    if len(fvals) > 0:
        dominant = np.argsort(fvals)[-10:][::-1]
        print("Dominante Frequenzen:")
        for d in dominant:
            print(int(d), float(fvals[d]))
    else:
        print("Keine spektralen Daten vorhanden.")

    print("\n===================================")
    print("Rolling Window")
    print("===================================")
    rvals = rolling_r(quads, w=100)
    finite = rvals[np.isfinite(rvals)]
    if len(finite) > 0:
        print("R_min =", np.min(finite))
        print("R_max =", np.max(finite))
        print("R_mean =", np.mean(finite))
        print("R_std =", np.std(finite))
        print("inf-Anteil =", np.mean(~np.isfinite(rvals)))
    else:
        print("Keine endlichen R-Werte im Rolling-Test.")

    plt.figure(figsize=(10, 5))
    plt.plot(rvals)
    plt.axhline(126, linestyle="--")
    plt.axhline(137.036, linestyle=":")
    plt.title("Rolling R(N)")
    plt.xlabel("Window")
    plt.ylabel("R")
    plt.grid()
    plt.tight_layout()
    plot_path = "rolling_R_1e8.png"
    plt.savefig(plot_path, dpi=150)
    print(f"\nPlot gespeichert: {plot_path}")