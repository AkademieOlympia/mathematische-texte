import math
from itertools import product

def normalize(q):
    s = sum(q)
    if s <= 0:
        raise ValueError("Summe der Komponenten muss > 0 sein.")
    return tuple(x / s for x in q)

def entropy(ps):
    out = 0.0
    for p in ps:
        if p > 0:
            out -= p * math.log(p)
    return out

def complexity(
    q,
    a=(1.0, 2.0, 3.5),      # Schalentiefen-Gewichte
    lambdas=(0.8, 1.2, 1.6),# Kopplungsgewichte lambda_12, lambda_13, lambda_23
    mu=1.0,                 # Gewicht Mischungsbeitrag
    nu=1.0                  # Gewicht Kopplungsbeitrag
):
    p1, p2, p3 = normalize(q)

    # Tiefe
    c_depth = a[0]*p1 + a[1]*p2 + a[2]*p3

    # Mischung
    c_mix = entropy((p1, p2, p3))

    # Kopplungen
    lam12, lam13, lam23 = lambdas
    c_coupl = lam12*p1*p2 + lam13*p1*p3 + lam23*p2*p3

    # Gesamt
    c_total = c_depth + mu*c_mix + nu*c_coupl

    return {
        "q": q,
        "p1": p1,
        "p2": p2,
        "p3": p3,
        "C_depth": c_depth,
        "C_mix": c_mix,
        "C_coupl": c_coupl,
        "C_total": c_total,
    }

def generate_states(max_q=10):
    states = []
    for q in product(range(max_q + 1), repeat=3):
        if sum(q) == 0:
            continue
        states.append(q)
    return states

def classify_state(rec):
    p1, p2, p3 = rec["p1"], rec["p2"], rec["p3"]
    support = sum(1 for x in (p1, p2, p3) if x > 1e-12)

    if support == 1:
        if p1 > 0:
            return "rein erste Schale"
        if p2 > 0:
            return "rein zweite Schale"
        return "rein dritte Schale"

    if support == 2:
        if p3 == 0:
            return "erste-zweite Mischung"
        if p2 == 0:
            return "erste-dritte Mischung"
        return "zweite-dritte Mischung"

    return "Dreischalen-Zustand"

def print_top_complex_states(max_q=8, top_k=20):
    states = generate_states(max_q=max_q)
    records = [complexity(q) for q in states]
    records.sort(key=lambda r: r["C_total"], reverse=True)

    print(f"Top {top_k} komplexeste Zustände bis q_i <= {max_q}:\n")
    header = (
        f"{'q':>12} | {'Typ':<24} | {'p=(p1,p2,p3)':<28} | "
        f"{'C_depth':>8} | {'C_mix':>8} | {'C_coupl':>8} | {'C_total':>8}"
    )
    print(header)
    print("-" * len(header))

    for rec in records[:top_k]:
        q = rec["q"]
        typ = classify_state(rec)
        pstr = f"({rec['p1']:.3f}, {rec['p2']:.3f}, {rec['p3']:.3f})"
        print(
            f"{str(q):>12} | {typ:<24} | {pstr:<28} | "
            f"{rec['C_depth']:8.3f} | {rec['C_mix']:8.3f} | "
            f"{rec['C_coupl']:8.3f} | {rec['C_total']:8.3f}"
        )

def compare_selected_states():
    examples = [
        (10, 0, 0),   # rein erste Schale
        (10, 2, 0),   # etwas zweite Schale
        (10, 0, 2),   # etwas dritte Schale
        (10, 2, 1),   # leichte Dreischalenmischung
        (8, 4, 2),    # stärker gemischt
        (5, 5, 5),    # maximal gemischt
        (2, 5, 8),    # tiefenlastig
    ]

    print("\nVergleich ausgewählter Zustände:\n")
    header = (
        f"{'q':>12} | {'Typ':<24} | {'p=(p1,p2,p3)':<28} | "
        f"{'C_depth':>8} | {'C_mix':>8} | {'C_coupl':>8} | {'C_total':>8}"
    )
    print(header)
    print("-" * len(header))

    for q in examples:
        rec = complexity(q)
        typ = classify_state(rec)
        pstr = f"({rec['p1']:.3f}, {rec['p2']:.3f}, {rec['p3']:.3f})"
        print(
            f"{str(q):>12} | {typ:<24} | {pstr:<28} | "
            f"{rec['C_depth']:8.3f} | {rec['C_mix']:8.3f} | "
            f"{rec['C_coupl']:8.3f} | {rec['C_total']:8.3f}"
        )

if __name__ == "__main__":
    print_top_complex_states(max_q=8, top_k=20)
    compare_selected_states()