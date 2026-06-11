"""
#Energiedoku: Fokus 137 – Quaternionen-Analyse ohne Sage.
"""

def quat_mul(q1, q2):
    """Quaternionen-Multiplikation: q1 * q2."""
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    return (
        a1*a2 - b1*b2 - c1*c2 - d1*d2,
        a1*b2 + b1*a2 + c1*d2 - d1*c2,
        a1*c2 - b1*d2 + c1*a2 + d1*b2,
        a1*d2 + b1*c2 - c1*b2 + d1*a2,
    )

def quat_conj(q):
    """Konjugiertes Quaternion: (a, -b, -c, -d)."""
    return (q[0], -q[1], -q[2], -q[3])

def quat_norm(q):
    """Norm |q|^2 = a^2 + b^2 + c^2 + d^2."""
    return sum(x*x for x in q)

def quat_real(q):
    """Realteil des Quaternions."""
    return q[0]

def analyze_energy_137():
    # 137 als Quaternion (Norm p=137): 10^2 + 6^2 + 1^2 + 0^2 = 137
    q_137 = (10, 6, 1, 0)

    # Benachbarte Primzahlen mit Repräsentanten (Summe von 4 Quadraten)
    # 131 = 9^2 + 7^2 + 1^2 + 0^2
    # 139 = 11^2 + 3^2 + 3^2 + 0^2
    neighbors = [
        (131, (9, 7, 1, 0)),
        (139, (11, 3, 3, 0)),
    ]

    n137 = quat_norm(q_137)
    print(f"--- #Energiedoku: Fokus 137 ---")
    print(f"Zustand q_137: {q_137} (Norm: {n137})")

    for p, q_p in neighbors:
        n_p = quat_norm(q_p)
        # Barandes-Korrelation: Re(conj(q_137) * q_p) / (|q_137| * |q_p|)
        conj_q137 = quat_conj(q_137)
        prod = quat_mul(conj_q137, q_p)
        correlation = quat_real(prod) / (n137**0.5 * n_p**0.5)
        print(f"Übergangschance 137 -> {p}: {correlation:.4f}")

if __name__ == "__main__":
    analyze_energy_137()
