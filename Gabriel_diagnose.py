import math
from pathlib import Path

import mpmath as mp


mp.mp.dps = 30
LI_2 = mp.li(2)


def prime_pi_values(bound):
    """Berechnet pi(x) fuer alle x < bound per Sieb des Eratosthenes."""
    if bound <= 2:
        return [0] * max(bound, 0)

    is_prime = [True] * bound
    is_prime[0] = False
    is_prime[1] = False

    for candidate in range(2, int(math.sqrt(bound - 1)) + 1):
        if is_prime[candidate]:
            start = candidate * candidate
            is_prime[start:bound:candidate] = [False] * len(is_prime[start:bound:candidate])

    values = [0] * bound
    count = 0
    for x in range(bound):
        if is_prime[x]:
            count += 1
        values[x] = count

    return values


def offset_li(x):
    """SageMath-Li(x)-nahe offset-logarithmische Integralfunktion."""
    return mp.li(x) - LI_2


def von_koch_bound(x):
    return (1 / (8 * mp.pi)) * mp.sqrt(x) * mp.log(x)


def analyze_compression_point(bound, start=1000):
    prime_counts = prime_pi_values(bound)
    max_ratio = mp.mpf("0")
    critical_x = 0
    critical_deviation = mp.mpf("0")
    critical_bound = mp.mpf("0")

    for x in range(max(3, start), bound):
        li_x = offset_li(x)
        deviation = abs(prime_counts[x] - li_x)
        bound_x = von_koch_bound(x)
        ratio = deviation / bound_x

        if ratio > max_ratio:
            max_ratio = ratio
            critical_x = x
            critical_deviation = deviation
            critical_bound = bound_x

    return {
        "critical_x": critical_x,
        "max_ratio": max_ratio,
        "deviation": critical_deviation,
        "bound": critical_bound,
    }


def zeta_wave_resonance(x_val):
    gammas = [mp.mpf("14.134725"), mp.mpf("21.022040"), mp.mpf("25.010858")]
    oscillation = mp.mpf("0")

    for gamma in gammas:
        oscillation += 2 * mp.sin(gamma * mp.log(x_val)) / (gamma * mp.sqrt(x_val))

    return oscillation


def check_far_field(x_val):
    prime_counts = prime_pi_values(x_val + 1)
    li_x = offset_li(x_val)
    deviation = prime_counts[x_val] - li_x
    bound_x = von_koch_bound(x_val)
    ratio = abs(deviation) / bound_x

    return {
        "x": x_val,
        "pi_x": prime_counts[x_val],
        "li_x": li_x,
        "deviation": deviation,
        "bound": bound_x,
        "ratio": ratio,
    }


def primes_in_window(start, end):
    """Listet alle Primzahlen im geschlossenen Intervall [start, end]."""
    if end < start:
        return []

    is_prime = [True] * (end + 1)
    if end >= 0:
        is_prime[0] = False
    if end >= 1:
        is_prime[1] = False

    for candidate in range(2, int(math.sqrt(end)) + 1):
        if is_prime[candidate]:
            for multiple in range(candidate * candidate, end + 1, candidate):
                is_prime[multiple] = False

    return [value for value in range(max(2, start), end + 1) if is_prime[value]]


def analyze_critical_window(start, end):
    """Analysiert Primzahlen und Luecken im kritischen Eichbereich."""
    primes = primes_in_window(start, end)
    gaps = []
    max_gap = 0
    max_gap_pair = None

    for previous, current in zip(primes, primes[1:]):
        gap = current - previous
        gaps.append((previous, current, gap))

        if gap > max_gap:
            max_gap = gap
            max_gap_pair = (previous, current)

    return {
        "start": start,
        "end": end,
        "primes": primes,
        "gaps": gaps,
        "max_gap": max_gap,
        "max_gap_pair": max_gap_pair,
    }


def quaternion_energy_density(bound, local_primes=None):
    """Berechnet die gewichtete Hurwitz-Massendichte fuer Primzahl-Knoten."""
    if local_primes is None:
        local_primes = [1409, 1423, 1427, 1429]

    primes = primes_in_window(2, bound)
    total_continuous_volume = offset_li(bound)
    total_quaternion_states = sum(24 * prime for prime in primes)
    normalization_rho = total_quaternion_states / total_continuous_volume
    local_impulses = []

    for prime in local_primes:
        local_mass = 24 * prime
        share = local_mass / total_quaternion_states * 100
        local_impulses.append((prime, local_mass, share))

    return {
        "bound": bound,
        "prime_count": len(primes),
        "total_continuous_volume": total_continuous_volume,
        "total_quaternion_states": total_quaternion_states,
        "normalization_rho": normalization_rho,
        "local_impulses": local_impulses,
    }


def fmt(value, digits=4):
    return f"{float(value):.{digits}f}"


def build_report():
    compression = analyze_compression_point(3000, start=1000)
    resonance_1500 = zeta_wave_resonance(1500)
    resonance_8000 = zeta_wave_resonance(8000)
    far_field = check_far_field(100000)
    critical_window = analyze_critical_window(1400, 1450)
    energy_density = quaternion_energy_density(2000)
    gap_lines = [
        f"Lücke zwischen {previous} und {current}: {gap}"
        for previous, current, gap in critical_window["gaps"]
    ]
    local_impulse_lines = [
        f"{prime:<5} | {mass:<23} | {fmt(share, 4)} %"
        for prime, mass, share in energy_density["local_impulses"]
    ]

    lines = [
        "--- Testlauf 1: Kritischer Eichpunkt ---",
        "Suchfenster: 1000 <= x < 3000",
        f"Maximale Annäherung an die Wand bei x = {compression['critical_x']}",
        f"Auslastung der Rigiditäts-Schranke: {fmt(compression['max_ratio'] * 100, 2)}%",
        f"Abweichung am Eichpunkt: {fmt(compression['deviation'], 4)}",
        f"Schrankenwert am Eichpunkt: ±{fmt(compression['bound'], 4)}",
        "",
        "--- Testlauf 2: Resonanz-Vergleich ---",
        f"Theoretische Resonanz-Amplitude bei x=1500: {fmt(resonance_1500, 4)}",
        f"Theoretische Resonanz-Amplitude bei x=8000: {fmt(resonance_8000, 4)}",
        "",
        "--- Testlauf 3: Fernfeld-Stichprobe ---",
        f"Bei x = {far_field['x']}: pi(x) = {far_field['pi_x']}, Li(x) = {fmt(far_field['li_x'], 2)}",
        f"Abweichung = {fmt(far_field['deviation'], 2)}, Schranke = ±{fmt(far_field['bound'], 2)}",
        f"Auslastung im Fernfeld: {fmt(far_field['ratio'] * 100, 2)}%",
        "",
        "--- Testlauf 4: Kritisches Fenster ---",
        f"Analyse des kritischen Fensters [{critical_window['start']}, {critical_window['end']}]",
        f"Primzahlen in diesem Bereich ({len(critical_window['primes'])} Stück):",
        str(critical_window["primes"]),
        "Lokale Lückenstruktur (P_n - P_{n-1}):",
        *gap_lines,
        (
            "Maximale arithmetische Spannung im Fenster: "
            f"Lücke von {critical_window['max_gap']} zwischen "
            f"{critical_window['max_gap_pair'][0]} und {critical_window['max_gap_pair'][1]}"
        ),
        "",
        "--- Testlauf 5: Quaternionische Feld-Schnittstelle ---",
        f"Analyse bis x = {energy_density['bound']}",
        f"Anzahl Primzahl-Knoten (Sphären): {energy_density['prime_count']}",
        f"Kontinuierliches Hintergrundpotential Li(x): {fmt(energy_density['total_continuous_volume'], 4)}",
        f"Gesamtzahl der Hurwitz-Zustände (Masse): {energy_density['total_quaternion_states']}",
        f"Globale Energiedichte (Knotenmasse / Kontinuum): {fmt(energy_density['normalization_rho'], 2)}",
        "",
        "Lokale gravitative Auswirkung im kritischen Fenster:",
        "p     | Hurwitz-Zustände (24*p) | Lokaler Impuls-Anteil (%)",
        "-" * 55,
        *local_impulse_lines,
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    report = build_report()
    print(report)

    output_path = Path("bamberger_diagnose_protokoll.txt")
    output_path.write_text(report + "\n", encoding="utf-8")
    print(f"\nProtokoll gespeichert als '{output_path}'.")
