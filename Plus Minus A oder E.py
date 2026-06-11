import sympy as sp
import numpy as np
from collections import Counter, defaultdict
from scipy.fft import rfft
from sklearn.linear_model import LinearRegression

INTERPRETATION = {
    "pair-AE": "sehr hohe FFT, niedrige Persistenz",
    "pair-AC": "hohe FFT, tendenziell unruhig",
    "multi": "unregelmaessig, mittlere FFT, niedrige Persistenz",
    "pair-BC": "nur mittlere FFT, entsprechend nur mittlere Modulation",
    "A-pure": "erhoehte reine Modulation, aber unter den Paar-Klassen",
    "C-pure": "niedrige FFT, hohe Stabilitaet",
    "E-pure": "niedrige FFT, hohe Stabilitaet",
    "B-pure": "sehr niedrige FFT, aber extrem hohe Stabilitaet",
}


def family(p):
    r = p % 12
    if r == 1:
        return "E"
    if r == 5:
        return "A"
    if r == 7:
        return "B"
    if r == 11:
        return "C"
    return None


def kernel_235(n):
    for p in [2, 3, 5]:
        while n % p == 0:
            n //= p
    return n


def classify_kernel(k):
    fac = sp.factorint(k)
    fams = []

    for p in fac:
        f = family(p)
        if f is not None:
            fams.append(f)

    if not fams:
        return "smooth"

    unique = sorted(set(fams))

    if len(unique) == 1:
        return f"{unique[0]}-pure"
    if len(unique) == 2:
        return "pair-" + "".join(unique)
    return "multi"


def logw(m):
    return np.log(float(m)) - np.log(float(sp.nextprime(int(m))))


def run_statistics(N=100000):
    counts = Counter()
    groups = defaultdict(set)

    for n in range(2, N + 1):
        k = kernel_235(n)
        cls = classify_kernel(k)

        counts[cls] += 1

        if cls != "smooth":
            groups[cls].add(k)

    groups = {cls: np.array(sorted(vals), dtype=float) for cls, vals in groups.items()}
    return counts, groups


def analyze_group(ms):
    ms = np.array(ms, dtype=float)

    if len(ms) < 2:
        return {
            "count": len(ms),
            "var": 0.0,
            "fft_peak": 0.0,
            "lag1": 0.0,
            "r2": 1.0,
        }

    ys = np.array([logw(m) for m in ms], dtype=float)
    X = np.log(ms).reshape(-1, 1)

    model = LinearRegression().fit(X, ys)
    y_pred = model.predict(X)
    res = ys - y_pred

    fft_vals = np.abs(rfft(res))
    if len(fft_vals) <= 1:
        peak_ratio = 0.0
    else:
        fft_sum = np.sum(fft_vals)
        peak_ratio = 0.0 if fft_sum == 0 else np.max(fft_vals[1:]) / fft_sum

    if len(res) < 2 or np.allclose(res, res[0]):
        lag1 = 0.0
    else:
        lag1 = np.corrcoef(res[:-1], res[1:])[0, 1]
        if np.isnan(lag1):
            lag1 = 0.0

    return {
        "count": len(ms),
        "var": float(np.var(res)),
        "fft_peak": float(peak_ratio),
        "lag1": float(lag1),
        "r2": float(model.score(X, ys)),
    }


def compute_fft_by_class(N=100000, min_group_size=50, verbose=True):
    _, groups = run_statistics(N)
    results = {}

    for key, ms in sorted(groups.items()):
        if len(ms) < min_group_size:
            continue
        results[key] = analyze_group(ms)

    if verbose:
        print(f"\n=== FFT nach Klassen fuer N = {N} ===")
        sorted_results = sorted(
            results.items(),
            key=lambda item: (-item[1]["fft_peak"], item[0])
        )
        for key, value in sorted_results:
            print(
                f"{key:10s} | count={value['count']:5d} | "
                f"fft={value['fft_peak']:.4f} | lag1={value['lag1']:.4f} | "
                f"var={value['var']:.6f} | r2={value['r2']:.4f}"
            )

    return results


def explain_model(model, feature_names=None):
    if feature_names is None:
        feature_names = ["E", "A", "B", "C", "AE", "AC", "BC"]

    coefs = np.asarray(model.coef_, dtype=float)

    print("\n=== Regressionsmodell ===")
    print(f"Intercept: {float(model.intercept_):.6f}")

    for name, value in sorted(
        zip(feature_names, coefs),
        key=lambda item: -abs(item[1])
    ):
        print(f"{name:>6s}: {value:+.6f}")


def full_analysis(N=100000, min_group_size=50):
    counts, groups = run_statistics(N)
    # Empirische Modulationshierarchie nach den bisher beobachteten FFT-Werten:
    # pair-AE > pair-AC > multi ~ pair-BC > A-pure > C-pure > E-pure > B-pure
    priority = {
        "pair-AE": 8,
        "pair-AC": 7,
        "multi": 6,
        "pair-BC": 5,
        "A-pure": 4,
        "C-pure": 3,
        "E-pure": 2,
        "B-pure": 1,
    }

    print("\n=== Haeufigkeiten ===")
    total = sum(counts.values())
    for key, value in sorted(counts.items()):
        print(f"{key:10s}: {value:8d} ({value / total:.3f})")

    print("\n=== Analyse ===")
    results = {}

    for key, ms in sorted(groups.items()):
        if len(ms) < min_group_size:
            continue
        results[key] = analyze_group(ms)

    print("\n--- nach qualitativer Hierarchie ---")
    sorted_results = sorted(
        results.items(),
        key=lambda item: (-priority.get(item[0], 0), -item[1]["fft_peak"])
    )

    for key, value in sorted_results:
        note = INTERPRETATION.get(key, "")
        print(
            f"{key:10s} | count={value['count']:5d} | "
            f"fft={value['fft_peak']:.4f} | var={value['var']:.6f} | "
            f"lag1={value['lag1']:.4f} | r2={value['r2']:.4f}"
        )
        if note:
            print(f"  -> {note}")

    return results


if __name__ == "__main__":
    full_analysis()