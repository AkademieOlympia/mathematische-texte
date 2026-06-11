import matplotlib.pyplot as plt


def is_prime_small(n: int) -> bool:
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    r = int(n**0.5)
    f = 3
    while f <= r:
        if n % f == 0:
            return False
        f += 2
    return True


def energiedoku_pilot_wave_model():
    # 1. Raum: Hamilton-Quaternionen (Schnitt reell = 1, Koeffizienten bei i, j)
    #    q = 1 + x*i + y*j  →  Norm N(q) = 1 + x² + y²
    # 2. „Quaternionen-Prim“ hier: N(q) ist eine rationale Primzahl (wie im Sage-Skript)

    points_to_check = []
    for x in range(-3, 4):
        for y in range(-3, 4):
            norm_val = float(1 + x * x + y * y)
            is_p = is_prime_small(int(norm_val))
            points_to_check.append((x, y, norm_val, is_p))

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")

    X = [p[0] for p in points_to_check]
    Y = [p[1] for p in points_to_check]
    Z = [p[2] for p in points_to_check]
    colors = ["red" if p[3] else "blue" for p in points_to_check]

    ax.scatter(X, Y, Z, c=colors, s=100, label="Blau: Komposit / Rot: Prim-Ziele")
    ax.set_title("#Energiedoku: Quaternionen Pilot-Wellen-Feld (Schnitt)")
    ax.set_xlabel("Im(i)")
    ax.set_ylabel("Im(j)")
    ax.set_zlabel("Energie-Potential (Norm)")

    print(
        "Modell generiert: Die roten Punkte sind die 'stabilen Zustände' (Primzahlen),"
    )
    print("zu denen die Pilotwelle die Information im Raum lenkt.")
    n_rot = sum(1 for p in points_to_check if p[3])
    n_blau = len(points_to_check) - n_rot
    print(f"Punkte: {len(points_to_check)} (rot=Prim-Norm: {n_rot}, blau: {n_blau})")
    plt.show()


if __name__ == "__main__":
    energiedoku_pilot_wave_model()
