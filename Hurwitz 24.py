from math import isqrt
from fractions import Fraction

# ============================================================
# Grundformatierung
# ============================================================

def format_scalar(x):
    if isinstance(x, Fraction):
        if x.denominator == 1:
            return str(x.numerator)
        return f"{x.numerator}/{x.denominator}"
    return str(x)

def format_quaternion(q):
    labels = ["", "i", "j", "k"]
    parts = []
    for x, lab in zip(q, labels):
        if x == 0:
            continue

        if not isinstance(x, Fraction):
            x = Fraction(x, 1)

        if lab == "":
            parts.append(format_scalar(x))
        else:
            if x == 1:
                parts.append(lab)
            elif x == -1:
                parts.append(f"-{lab}")
            else:
                parts.append(f"{format_scalar(x)}{lab}")

    if not parts:
        return "0"

    s = parts[0]
    for part in parts[1:]:
        if part.startswith("-"):
            s += " - " + part[1:]
        else:
            s += " + " + part
    return s


# ============================================================
# Quaternionen-Arithmetik
# q = (a,b,c,d) entspricht a + bi + cj + dk
# ============================================================

def Q(a, b, c, d):
    return (Fraction(a), Fraction(b), Fraction(c), Fraction(d))

def q_add(x, y):
    return tuple(xi + yi for xi, yi in zip(x, y))

def q_neg(x):
    return tuple(-xi for xi in x)

def q_conj(x):
    a, b, c, d = x
    return (a, -b, -c, -d)

def q_norm(x):
    a, b, c, d = x
    return a*a + b*b + c*c + d*d

def q_mul(x, y):
    a, b, c, d = x
    e, f, g, h = y
    return (
        a*e - b*f - c*g - d*h,
        a*f + b*e + c*h - d*g,
        a*g - b*h + c*e + d*f,
        a*h + b*g - c*f + d*e,
    )

def q_eq(x, y):
    return all(xi == yi for xi, yi in zip(x, y))


# ============================================================
# Erzeugung von X_p
# ============================================================

def generate_integer_solutions(p):
    sols = set()
    B = isqrt(p)

    for a in range(-B, B + 1):
        a2 = a * a
        if a2 > p:
            continue
        for b in range(-B, B + 1):
            ab2 = a2 + b * b
            if ab2 > p:
                continue
            for c in range(-B, B + 1):
                abc2 = ab2 + c * c
                if abc2 > p:
                    continue
                d2 = p - abc2
                d = isqrt(d2)
                if d * d == d2:
                    for dd in {d, -d}:
                        sols.add(Q(a, b, c, dd))

    return sorted(sols)

def generate_half_integer_solutions(p):
    sols = set()
    target = 4 * p
    B = isqrt(target)
    odd_vals = [m for m in range(-B, B + 1) if m % 2 != 0]

    for m1 in odd_vals:
        s1 = m1 * m1
        if s1 > target:
            continue
        for m2 in odd_vals:
            s2 = s1 + m2 * m2
            if s2 > target:
                continue
            for m3 in odd_vals:
                s3 = s2 + m3 * m3
                if s3 > target:
                    continue
                rem = target - s3
                m4 = isqrt(rem)
                if m4 * m4 == rem and m4 % 2 == 1:
                    for mm4 in {m4, -m4}:
                        sols.add(Q(Fraction(m1, 2), Fraction(m2, 2), Fraction(m3, 2), Fraction(mm4, 2)))

    return sorted(sols)

def generate_Xp(p):
    int_part = generate_integer_solutions(p)
    half_part = generate_half_integer_solutions(p)
    return sorted(set(int_part + half_part))


# ============================================================
# Die 24 Hurwitz-Einheiten
# Norm 1 in O_H
# ============================================================

def hurwitz_units():
    units = set()

    # ±1, ±i, ±j, ±k
    basic = [
        Q(1,0,0,0), Q(-1,0,0,0),
        Q(0,1,0,0), Q(0,-1,0,0),
        Q(0,0,1,0), Q(0,0,-1,0),
        Q(0,0,0,1), Q(0,0,0,-1),
    ]
    units.update(basic)

    # (±1 ± i ± j ± k)/2
    for s1 in [-1, 1]:
        for s2 in [-1, 1]:
            for s3 in [-1, 1]:
                for s4 in [-1, 1]:
                    units.add(Q(Fraction(s1,2), Fraction(s2,2), Fraction(s3,2), Fraction(s4,2)))

    units = sorted(units)

    # Sicherheitscheck: Norm = 1
    bad = [u for u in units if q_norm(u) != 1]
    if bad:
        raise ValueError(f"Es gibt Einheiten mit Norm != 1: {bad}")

    if len(units) != 24:
        raise ValueError(f"Es wurden {len(units)} statt 24 Hurwitz-Einheiten erzeugt.")

    return units


# ============================================================
# Orbit-Berechnung
# ============================================================

def left_orbit(x, units):
    return {q_mul(u, x) for u in units}

def right_orbit(x, units):
    return {q_mul(x, u) for u in units}

def double_orbit(x, units):
    return {q_mul(u1, q_mul(x, u2)) for u1 in units for u2 in units}

def orbit_partition(X, orbit_func, units):
    Xset = set(X)
    seen = set()
    orbits = []

    for x in X:
        if x in seen:
            continue
        orb = orbit_func(x, units)
        orb = orb & Xset
        orbits.append(orb)
        seen |= orb

    return orbits


# ============================================================
# Bericht
# ============================================================

def print_orbit_report(p, mode="left", max_orbits=10, max_elements_per_orbit=8):
    X = generate_Xp(p)
    units = hurwitz_units()

    if mode == "left":
        orbit_func = left_orbit
        mode_name = "Linksorbits"
    elif mode == "right":
        orbit_func = right_orbit
        mode_name = "Rechtsorbits"
    elif mode == "double":
        orbit_func = double_orbit
        mode_name = "Doppelorbits"
    else:
        raise ValueError("mode muss 'left', 'right' oder 'double' sein")

    orbits = orbit_partition(X, orbit_func, units)
    orbit_sizes = sorted([len(O) for O in orbits])

    print("=" * 80)
    print(f"Orbitbericht für p={p} ({mode_name})")
    print("=" * 80)
    print(f"|X_{p}| = {len(X)}")
    print(f"Anzahl Hurwitz-Einheiten = {len(units)}")
    print(f"Anzahl Orbits = {len(orbits)}")
    print(f"Orbitgrößen (sortiert) = {orbit_sizes}")
    print()

    all_size_24 = all(s == 24 for s in orbit_sizes)
    print(f"Alle Orbitgrößen = 24? {all_size_24}")
    if p > 2:
        print(f"Erwartung p+1 = {p+1}, tatsächlich Anzahl Orbits = {len(orbits)}")
        print(f"Anzahl Orbits = p+1? {len(orbits) == p+1}")
    print()

    for idx, O in enumerate(orbits[:max_orbits], 1):
        O_sorted = sorted(O)
        rep = O_sorted[0]
        print(f"Orbit {idx}: Größe {len(O)}")
        print(f"  Repräsentant: {format_quaternion(rep)}")
        for y in O_sorted[:max_elements_per_orbit]:
            print(f"    - {format_quaternion(y)}")
        if len(O_sorted) > max_elements_per_orbit:
            print(f"    ... ({len(O_sorted) - max_elements_per_orbit} weitere)")
        print()


# ============================================================
# Zusätzliche Tests: Stabilität und freie Wirkung
# ============================================================

def stabilizer_left(x, units):
    return {u for u in units if q_mul(u, x) == x}

def test_free_left_action(p):
    X = generate_Xp(p)
    units = hurwitz_units()

    bad = []
    for x in X:
        stab = stabilizer_left(x, units)
        if len(stab) != 1:
            bad.append((x, stab))

    print("=" * 80)
    print(f"Test auf freie Linkswirkung für p={p}")
    print("=" * 80)
    if not bad:
        print("Die Linkswirkung der 24 Hurwitz-Einheiten auf X_p ist frei.")
        print("Also hat jedes Element trivialen Stabilisator.")
    else:
        print(f"Es gibt {len(bad)} Elemente mit nichttrivialem Stabilisator.")
        for x, stab in bad[:10]:
            print(f"x = {format_quaternion(x)}")
            print("Stabilisator:")
            for u in sorted(stab):
                print(f"  {format_quaternion(u)}")
            print()


# ============================================================
# Hauptprogramm
# ============================================================

if __name__ == "__main__":
    for p in [2, 3, 5, 7]:
        print_orbit_report(p, mode="left", max_orbits=6, max_elements_per_orbit=10)
        test_free_left_action(p)
        print()