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
# ============================================================

def Q(a, b, c, d):
    return (Fraction(a), Fraction(b), Fraction(c), Fraction(d))

def q_mul(x, y):
    a, b, c, d = x
    e, f, g, h = y
    return (
        a*e - b*f - c*g - d*h,
        a*f + b*e + c*h - d*g,
        a*g - b*h + c*e + d*f,
        a*h + b*g - c*f + d*e,
    )

def q_norm(x):
    a, b, c, d = x
    return a*a + b*b + c*c + d*d


# ============================================================
# X_p
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
    return sorted(set(generate_integer_solutions(p) + generate_half_integer_solutions(p)))


# ============================================================
# Hurwitz-Einheiten
# ============================================================

def hurwitz_units():
    units = set()

    # 8 "offensichtliche" Einheiten
    units.update([
        Q(1,0,0,0), Q(-1,0,0,0),
        Q(0,1,0,0), Q(0,-1,0,0),
        Q(0,0,1,0), Q(0,0,-1,0),
        Q(0,0,0,1), Q(0,0,0,-1),
    ])

    # 16 halbzahlig-hurwitzsche Einheiten
    for s1 in [-1, 1]:
        for s2 in [-1, 1]:
            for s3 in [-1, 1]:
                for s4 in [-1, 1]:
                    units.add(Q(Fraction(s1,2), Fraction(s2,2), Fraction(s3,2), Fraction(s4,2)))

    units = sorted(units)

    if len(units) != 24:
        raise ValueError(f"Es wurden {len(units)} statt 24 Einheiten erzeugt.")

    bad = [u for u in units if q_norm(u) != 1]
    if bad:
        raise ValueError(f"Einheiten mit Norm != 1 gefunden: {bad}")

    return units


# ============================================================
# Orbits
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
        orb = orbit_func(x, units) & Xset
        orbits.append(frozenset(orb))
        seen |= orb

    return orbits


# ============================================================
# Freie Wirkung testen
# ============================================================

def left_stabilizer(x, units):
    return {u for u in units if q_mul(u, x) == x}

def right_stabilizer(x, units):
    return {u for u in units if q_mul(x, u) == x}

def test_free_action(p, side="left"):
    X = generate_Xp(p)
    units = hurwitz_units()

    bad = []
    for x in X:
        stab = left_stabilizer(x, units) if side == "left" else right_stabilizer(x, units)
        if len(stab) != 1:
            bad.append((x, stab))

    print("=" * 80)
    print(f"Test auf freie {('Links' if side=='left' else 'Rechts')}wirkung für p={p}")
    print("=" * 80)
    if not bad:
        print("Die Wirkung ist frei.")
    else:
        print(f"Nichttriviale Stabilisatoren bei {len(bad)} Elementen.")
        for x, stab in bad[:10]:
            print(f"x = {format_quaternion(x)}")
            print("Stabilisator:")
            for u in sorted(stab):
                print(f"  {format_quaternion(u)}")
            print()


# ============================================================
# Vergleich Links- vs. Rechtsorbits
# ============================================================

def canonicalize_partition(orbits):
    return {frozenset(O) for O in orbits}

def compare_left_right_orbits(p, max_show=5):
    X = generate_Xp(p)
    units = hurwitz_units()

    left_orbits = orbit_partition(X, left_orbit, units)
    right_orbits = orbit_partition(X, right_orbit, units)

    left_set = canonicalize_partition(left_orbits)
    right_set = canonicalize_partition(right_orbits)

    print("=" * 80)
    print(f"Vergleich Links- und Rechtsorbits für p={p}")
    print("=" * 80)
    print(f"Anzahl Linksorbits : {len(left_orbits)}")
    print(f"Anzahl Rechtsorbits: {len(right_orbits)}")
    print(f"Partitionen identisch? {left_set == right_set}")
    print(f"Linksortgrößen : {sorted(len(O) for O in left_orbits)}")
    print(f"Rechtsortgrößen: {sorted(len(O) for O in right_orbits)}")
    print()

    if left_set != right_set:
        print("Beispielhafte Linksorbits:")
        for i, O in enumerate(left_orbits[:max_show], 1):
            rep = sorted(O)[0]
            print(f"L{i}: Größe {len(O)}, Rep {format_quaternion(rep)}")
        print()
        print("Beispielhafte Rechtsorbits:")
        for i, O in enumerate(right_orbits[:max_show], 1):
            rep = sorted(O)[0]
            print(f"R{i}: Größe {len(O)}, Rep {format_quaternion(rep)}")
        print()


# ============================================================
# Doppelorbits
# ============================================================

def double_orbit_partition(X, units):
    Xset = set(X)
    seen = set()
    orbits = []

    for x in X:
        if x in seen:
            continue
        orb = double_orbit(x, units) & Xset
        orbits.append(frozenset(orb))
        seen |= orb

    return orbits

def report_double_orbits(p, max_show=8):
    X = generate_Xp(p)
    units = hurwitz_units()
    dorbits = double_orbit_partition(X, units)

    print("=" * 80)
    print(f"Doppelorbits für p={p}")
    print("=" * 80)
    print(f"|X_{p}| = {len(X)}")
    print(f"Anzahl Doppelorbits = {len(dorbits)}")
    print(f"Doppelorbitgrößen = {sorted(len(O) for O in dorbits)}")
    print()

    for i, O in enumerate(dorbits[:max_show], 1):
        rep = sorted(O)[0]
        print(f"D{i}: Größe {len(O)}, Rep {format_quaternion(rep)}")


# ============================================================
# Gesamtbericht
# ============================================================

def full_report(p):
    X = generate_Xp(p)
    units = hurwitz_units()

    left_orbits = orbit_partition(X, left_orbit, units)
    right_orbits = orbit_partition(X, right_orbit, units)

    print("=" * 80)
    print(f"Gesamtbericht p={p}")
    print("=" * 80)
    print(f"|X_{p}| = {len(X)}")
    print(f"Linksortgrößen : {sorted(len(O) for O in left_orbits)}")
    print(f"Rechtsortgrößen: {sorted(len(O) for O in right_orbits)}")
    print(f"Anzahl Linksorbits : {len(left_orbits)}")
    print(f"Anzahl Rechtsorbits: {len(right_orbits)}")
    print(f"Links/Rechts identisch? {canonicalize_partition(left_orbits) == canonicalize_partition(right_orbits)}")
    print()


if __name__ == "__main__":
    for p in [2, 3, 5, 7]:
        full_report(p)
        test_free_action(p, side="left")
        test_free_action(p, side="right")
        compare_left_right_orbits(p)
        report_double_orbits(p)
        print()