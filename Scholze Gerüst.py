from math import isqrt, gcd
from fractions import Fraction

# ------------------------------------------------------------
# Grundfunktionen
# ------------------------------------------------------------

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

        val = format_scalar(x)

        if lab == "":
            parts.append(val)
        else:
            if x == 1:
                parts.append(lab)
            elif x == -1:
                parts.append(f"-{lab}")
            else:
                parts.append(f"{val}{lab}")

    if not parts:
        return "0"

    s = parts[0]
    for part in parts[1:]:
        if part.startswith("-"):
            s += " - " + part[1:]
        else:
            s += " + " + part
    return s


# ------------------------------------------------------------
# Erzeugung von X_p
# ------------------------------------------------------------

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
                        sols.add((a, b, c, dd))

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
                        sols.add((
                            Fraction(m1, 2),
                            Fraction(m2, 2),
                            Fraction(m3, 2),
                            Fraction(mm4, 2),
                        ))

    return sorted(sols)


def generate_Xp(p):
    int_part = generate_integer_solutions(p)
    half_part = generate_half_integer_solutions(p)
    return int_part + half_part


# ------------------------------------------------------------
# Modulare Reduktion
# ------------------------------------------------------------

def inv_mod(a, n):
    a %= n
    if gcd(a, n) != 1:
        raise ValueError(f"{a} ist nicht invertierbar modulo {n}")
    # erweiterter euklidischer Algorithmus
    t, new_t = 0, 1
    r, new_r = n, a
    while new_r != 0:
        q = r // new_r
        t, new_t = new_t, t - q * new_t
        r, new_r = new_r, r - q * new_r
    if t < 0:
        t += n
    return t


def reduce_fraction_mod_N(x, N):
    """
    Reduziert x modulo N.
    Für Fraction(num, den) wird num * den^{-1} mod N verwendet.
    Das setzt voraus, dass den invertierbar modulo N ist.
    """
    if isinstance(x, Fraction):
        num = x.numerator % N
        den = x.denominator % N
        den_inv = inv_mod(den, N)
        return (num * den_inv) % N
    return x % N


def reduce_quaternion_mod_N(q, N):
    return tuple(reduce_fraction_mod_N(x, N) for x in q)


def negate_mod_N(qmod, N):
    return tuple((-x) % N for x in qmod)


# ------------------------------------------------------------
# Kandidatenmenge X_p mod N
# ------------------------------------------------------------

def image_Xp_mod_N(p, N):
    Xp = generate_Xp(p)
    image = {reduce_quaternion_mod_N(q, N) for q in Xp}
    return sorted(image)


def split_preimages_by_residue(p, N):
    """
    Liefert ein Dictionary:
    residue mod N -> Liste von Urbildern in X_p
    """
    Xp = generate_Xp(p)
    buckets = {}
    for q in Xp:
        r = reduce_quaternion_mod_N(q, N)
        buckets.setdefault(r, []).append(q)
    return buckets


# ------------------------------------------------------------
# Einfache Tests für eine Kandidaten-Schrittmengenstruktur
# ------------------------------------------------------------

def symmetry_report_mod_N(p, N):
    residues = image_Xp_mod_N(p, N)
    residue_set = set(residues)

    missing_negatives = []
    for r in residues:
        if negate_mod_N(r, N) not in residue_set:
            missing_negatives.append(r)

    return {
        "p": p,
        "N": N,
        "cardinality": len(residues),
        "symmetric_under_negation": len(missing_negatives) == 0,
        "missing_negatives": missing_negatives,
    }


def print_mod_N_report(p, N, max_entries=50):
    print("=" * 80)
    print(f"Bericht für p={p}, N={N}")
    print("=" * 80)

    if gcd(2, N) != 1:
        print("WARNUNG: N ist gerade. Halbzahlige Reduktion modulo N ist in diesem Skript")
        print("nicht sauber modelliert, weil 2 modulo N nicht invertierbar ist.")
        print("Bitte zunächst nur ungerade N verwenden.")
        print()
        return

    Xp = generate_Xp(p)
    residues = image_Xp_mod_N(p, N)
    buckets = split_preimages_by_residue(p, N)
    sym = symmetry_report_mod_N(p, N)

    print(f"|X_{p}| = {len(Xp)}")
    print(f"|X_{p} mod {N}| = {len(residues)}")
    print(f"Symmetrisch unter q -> -q mod {N}: {sym['symmetric_under_negation']}")
    if not sym["symmetric_under_negation"]:
        print(f"Fehlende Gegenelemente: {len(sym['missing_negatives'])}")
    print()

    print("Erste Residuen modulo N:")
    for idx, r in enumerate(residues[:max_entries], 1):
        print(f"{idx:3d}: {r}")

    if len(residues) > max_entries:
        print(f"... ({len(residues) - max_entries} weitere)")
    print()

    print("Beispielhafte Urbilder pro Restklasse:")
    for idx, r in enumerate(residues[:min(15, len(residues))], 1):
        pre = buckets[r]
        print(f"{idx:3d}: Residuum {r}")
        for q in pre[:4]:
            print(f"     - {q}   ->   {format_quaternion(q)}")
        if len(pre) > 4:
            print(f"     ... ({len(pre) - 4} weitere Urbilder)")
    print()


# ------------------------------------------------------------
# Demo
# ------------------------------------------------------------

if __name__ == "__main__":
    test_cases = [
        (2, 3),
        (3, 5),
        (5, 3),
        (5, 7),
        (7, 5),
    ]

    for p, N in test_cases:
        print_mod_N_report(p, N, max_entries=30)