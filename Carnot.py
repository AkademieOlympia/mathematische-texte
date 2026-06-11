# pyright: reportMissingImports=false
# SageMath 10.5: Carnot-Ptolemaeus-Validierung im eabc-Modell
from functools import lru_cache
import os
import subprocess
import sys

SAGE_BIN = "/Applications/SageMath-10-8.app/Contents/Frameworks/Sage.framework/Versions/Current/local/bin/sage"

try:
    from sage.all import (
        acos,
        graphics_array,
        I,
        QQ,
        QuaternionAlgebra,
        RR,
        four_squares,
        is_prime,
        log,
        matrix,
        matrix_plot,
        next_prime,
        pi,
        sqrt,
        zeta,
    )
except ModuleNotFoundError:
    if __name__ == "__main__" and os.environ.get("CARNOT_REEXEC_UNDER_SAGE") != "1":
        env = os.environ.copy()
        env["CARNOT_REEXEC_UNDER_SAGE"] = "1"
        completed = subprocess.run([SAGE_BIN, os.path.abspath(__file__), *sys.argv[1:]], env=env)
        raise SystemExit(completed.returncode)
    raise RuntimeError(
        "Dieses Skript benoetigt SageMath. Starte es mit 'sage Carnot.py' "
        "oder direkt mit dem SageMath-10.8-Binaer."
    )

# Quaternionen-Algebra definieren
Q = QuaternionAlgebra(QQ, -1, -1)
QR = QuaternionAlgebra(RR, -1, -1)

# Vier symmetrische Testzustände des Carnot-Ptolemaeus-Zirkels
ZUSTAENDE = (
    Q([1, 0, 0, 0]),   # A
    Q([0, 1, 0, 0]),   # B
    Q([-1, 0, 0, 0]),  # C
    Q([0, -1, 0, 0]),  # D
)

def ptolemy_check(A, B, C, D):
    """Prüft die Sehnenviereck-Bedingung für vier quaternionische Zustände."""
    # Seitenlängen (Distanzen zwischen den Zuständen)
    a = sqrt((B - A).reduced_norm()) # Seite AB
    b = sqrt((C - B).reduced_norm()) # Seite BC
    c = sqrt((D - C).reduced_norm()) # Seite CD
    d = sqrt((A - D).reduced_norm()) # Seite DA
    
    # Diagonalen
    e = sqrt((C - A).reduced_norm()) # Diagonale AC
    f = sqrt((D - B).reduced_norm()) # Diagonale BD
    
    ptolemy_sum = a*c + b*d
    diagonal_prod = e*f
    
    diff = abs(ptolemy_sum - diagonal_prod)
    return ptolemy_sum, diagonal_prod, diff

def berechne_system_resonanz(rotation_offset, zustaende=ZUSTAENDE):
    """Summiert |zeta(1/2 + i*t)| ueber alle gegebenen Zustaende."""
    gesamt_resonanz = RR(0)
    for q in zustaende:
        # Fuer Quaternionen-Elemente ist die reduzierte Norm die stabile Groesse.
        norm_wert = RR(q.reduced_norm()) * RR(rotation_offset)
        s = RR(1) / 2 + I * norm_wert
        gesamt_resonanz += abs(zeta(s))

    return gesamt_resonanz


def scan_energieniveaus(start=14.0, ende=30.0, schritte=100, zustaende=ZUSTAENDE):
    """Sucht im vorgegebenen Intervall das Resonanzminimum."""
    if schritte < 2:
        raise ValueError("schritte muss mindestens 2 sein.")

    delta = (RR(ende) - RR(start)) / (schritte - 1)
    energieniveaus = [RR(start) + n * delta for n in range(schritte)]
    ergebnisse = [(e, berechne_system_resonanz(e, zustaende)) for e in energieniveaus]
    optimal_e, min_res = min(ergebnisse, key=lambda x: x[1])
    return optimal_e, min_res, ergebnisse


def find_prime_clusters(limit_triplets=100, limit_quads=50):
    """Findet Primzahldrillinge und Primzahlvierlinge in Standardmustern."""
    triplets = []
    quads = []
    p = 2

    while len(triplets) < limit_triplets or len(quads) < limit_quads:
        p = next_prime(p)

        # Klassisches Vierlingsmuster: (p, p+2, p+6, p+8)
        if is_prime(p + 2) and is_prime(p + 6) and is_prime(p + 8):
            if len(quads) < limit_quads:
                quads.append((p, p + 2, p + 6, p + 8))

        # Zulaessiges Drillingsmuster: (p, p+2, p+6)
        if is_prime(p + 2) and is_prime(p + 6):
            if len(triplets) < limit_triplets:
                triplets.append((p, p + 2, p + 6))

    return triplets, quads


@lru_cache(maxsize=None)
def quat_representation(n):
    """Gibt eine Darstellung von n als Summe von vier Quadraten zurueck."""
    return tuple(int(x) for x in four_squares(int(n)))


def integer_to_quaternion(n):
    """Bettet eine Zahl ueber ihre Vier-Quadrate-Darstellung in Q ein."""
    return Q(quat_representation(n))


def get_imaginary_vector(n):
    """Extrahiert den imaginaeren Anteil der Vier-Quadrate-Darstellung."""
    _, b, c, d = quat_representation(n)
    return Q([0, b, c, d])


def quaternion_imaginary_part(q):
    """Projiziert eine Quaternion auf ihren rein imaginaeren Anteil."""
    _, b, c, d = q.coefficient_tuple()
    return q.parent()([0, b, c, d])


def cluster_to_quaternions(cluster):
    """Uebersetzt ein Primzahl-Cluster in Quaternionen-Zustaende."""
    return tuple(integer_to_quaternion(n) for n in cluster)


def imag_vector_dot(u, v):
    """Skalarprodukt zweier rein imaginaerer Quaternionen."""
    _, ux, uy, uz = u.coefficient_tuple()
    _, vx, vy, vz = v.coefficient_tuple()
    return ux * vx + uy * vy + uz * vz


def cross_product_area(u, v):
    """Betrag des Kreuzprodukts der imaginaeren Anteile."""
    area_sq = u.reduced_norm() * v.reduced_norm() - imag_vector_dot(u, v) ** 2
    if area_sq < 0:
        area_sq = 0
    return sqrt(area_sq)


def calculate_cycle_work(quadruplet):
    """Summiert die Kreuzprodukt-Flaechen eines geschlossenen Viererzyklus."""
    vectors = [get_imaginary_vector(p) for p in quadruplet]
    work = RR(0)
    for index, current_vector in enumerate(vectors):
        next_vector = vectors[(index + 1) % len(vectors)]
        work += cross_product_area(current_vector, next_vector)
    return work


def calculate_cycle_work_for_states(states):
    """Summiert die Kreuzprodukt-Flaechen eines Zyklus aus Quaternionen-Zustaenden."""
    vectors = [quaternion_imaginary_part(state) for state in states]
    work = RR(0)
    for index, current_vector in enumerate(vectors):
        next_vector = vectors[(index + 1) % len(vectors)]
        work += cross_product_area(current_vector, next_vector)
    return work


def calculate_phase_shift(quadruplet_quats):
    """Berechnet die kumulative Phasenabweichung relativ zu 2*pi."""
    total_phase = RR(0)
    vectors = [quaternion_imaginary_part(q) for q in quadruplet_quats]

    for index, v1 in enumerate(vectors):
        v2 = vectors[(index + 1) % len(vectors)]
        norm_product = v1.reduced_norm() * v2.reduced_norm()
        if norm_product == 0:
            continue

        cos_theta = RR(imag_vector_dot(v1, v2) / sqrt(norm_product))
        cos_theta = min(RR(1), max(RR(-1), cos_theta))
        total_phase += acos(cos_theta)

    return total_phase - 2 * pi


def zeta_resonance(q, rotation_offset=1.0):
    """Misst die Zeta-Resonanz einer Quaternion ueber ihre reduzierte Norm."""
    t = RR(q.reduced_norm()) * RR(rotation_offset)
    s = RR(1) / 2 + I * t
    return abs(zeta(s))


def kanonische_energie_minimierung(starr_8, frei_4, rotation_offset=1.0, epsilon=1e-12):
    """
    Bewertet eine 8+4-Konfiguration.

    `starr_8` definiert den Carnot-Rahmen, `frei_4` die beweglichen Resonanzzustaende.
    """
    if len(starr_8) != 8:
        raise ValueError("starr_8 muss genau 8 Quaternionen enthalten.")
    if len(frei_4) != 4:
        raise ValueError("frei_4 muss genau 4 Quaternionen enthalten.")

    rahmen_energie = calculate_cycle_work_for_states(starr_8)
    resonanz = sum(zeta_resonance(q, rotation_offset) for q in frei_4)
    resonanz = max(RR(epsilon), RR(resonanz))
    return rahmen_energie / resonanz


def analysiere_8plus4_konfiguration(starr_8, frei_4, rotation_offset=1.0, epsilon=1e-12):
    """Gibt alle Kennzahlen der 8+4-Energiedoku-Konfiguration zurueck."""
    rahmen_energie = calculate_cycle_work_for_states(starr_8)
    resonanz = sum(zeta_resonance(q, rotation_offset) for q in frei_4)
    resonanz = max(RR(epsilon), RR(resonanz))
    return {
        "rahmen_energie": rahmen_energie,
        "resonanz": resonanz,
        "effizienz_index": kanonische_energie_minimierung(
            starr_8,
            frei_4,
            rotation_offset=rotation_offset,
            epsilon=epsilon,
        ),
    }


def validate_alpha_coupling(index, mu_target):
    """
    Berechnet eine heuristische Alpha^-1-Vorhersage aus dem eabc-Defekt.
    """
    defekt = mu_target - index
    if defekt == 0:
        raise ValueError("index und mu_target duerfen nicht identisch sein.")
    return (index / defekt) * (RR(8) / RR(12))


def local_eabc_index_for_triplet(triplet):
    """Ein lokaler eabc-Index aus den imaginaeren Vier-Quadrate-Anteilen."""
    q_sum = sum(get_imaginary_vector(x).reduced_norm() for x in triplet)
    return RR(q_sum) / RR(12)


def batch_process_kanon(limit=1000, observable=None, mu_target=1836.15267):
    """
    Fuehrt einen statistischen Kanon-Check fuer viele Primzahldrillinge durch.
    """
    if limit < 1:
        raise ValueError("limit muss mindestens 1 sein.")

    observable = RR(observable if observable is not None else "0.13387")
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    triplets, _ = find_prime_clusters(limit_triplets=limit, limit_quads=0)
    beta_distribution = []

    for triplet in triplets:
        xi_local = local_eabc_index_for_triplet(triplet)
        beta_local = (RR(mu_target) - xi_local) / observable
        beta_distribution.append(float(beta_local))

    mean_beta = sum(beta_distribution) / len(beta_distribution)
    variance = sum((value - mean_beta) ** 2 for value in beta_distribution) / len(beta_distribution)
    std_beta = variance ** 0.5

    return {
        "count": len(beta_distribution),
        "mean_beta": mean_beta,
        "std_beta": std_beta,
        "min_beta": min(beta_distribution),
        "max_beta": max(beta_distribution),
        "beta_distribution": beta_distribution,
        "clusters": triplets,
    }


def analyze_outlier(distribution, clusters, target_beta=10.9446):
    """Findet den Drilling mit der groessten Abweichung vom Ziel-Beta."""
    if not distribution or not clusters or len(distribution) != len(clusters):
        raise ValueError("distribution und clusters muessen gleich lang und nicht leer sein.")

    diffs = [abs(beta - target_beta) for beta in distribution]
    outlier_idx = diffs.index(max(diffs))
    return {
        "cluster": clusters[outlier_idx],
        "beta": distribution[outlier_idx],
        "delta": diffs[outlier_idx],
        "index": outlier_idx,
    }


def calculate_stress_matrix(cluster):
    """Berechnet eine 3x3-Stress-Matrix fuer einen Primzahldrilling."""
    pts = [QR(state.coefficient_tuple()) for state in cluster_to_quaternions(cluster)]
    stress_matrix = matrix(RR, 3, 3)

    for row in range(3):
        for col in range(3):
            product_imag = quaternion_imaginary_part(pts[row] * pts[col])
            stress_matrix[row, col] = sqrt(product_imag.reduced_norm())

    return stress_matrix


def plot_stress_comparison(kanon_matrix, outlier_matrix, output_path):
    """Speichert zwei Stress-Matrizen als nebeneinander liegende Heatmaps."""
    p1 = matrix_plot(kanon_matrix, cmap="viridis")
    p2 = matrix_plot(outlier_matrix, cmap="inferno")
    combined = graphics_array([p1, p2])
    combined.save(
        output_path,
        figsize=[12, 5],
        axes=False,
        title="Stress-Matrix Vergleich: Kanon vs. Ausreisser",
    )
    return output_path


def ptolemy_defect(points):
    """Gibt nur den Ptolemaeus-Defekt eines Quaternionen-Vierlings zurueck."""
    return ptolemy_check(*points)[2]


def simulate_resilience(
    quadruplet,
    iterations=50,
    step_size=0.1,
    perturbation=None,
    gradient_eps=0.05,
):
    """Simuliert die Rueckkehr zum Ptolemaeus-Gleichgewicht nach einer Stoerung."""
    states = [QR(state.coefficient_tuple()) for state in cluster_to_quaternions(quadruplet)]
    perturbation = perturbation or QR([0, RR(1) / 2, RR(1) / 2, RR(1) / 2])
    states[2] += perturbation

    direction_x = QR([0, 1, 0, 0])
    direction_y = QR([0, 0, 1, 0])
    direction_z = QR([0, 0, 0, 1])
    current_pts = list(states)
    current_step = RR(step_size)
    eps = RR(gradient_eps)
    history = []

    for _ in range(iterations):
        err = ptolemy_defect(current_pts)
        history.append(float(err))

        gradient_components = []
        for direction in (direction_x, direction_y, direction_z):
            plus_pts = list(current_pts)
            minus_pts = list(current_pts)
            delta = QR([0, *(eps * coeff for coeff in direction.coefficient_tuple()[1:])])
            plus_pts[2] += delta
            minus_pts[2] -= delta
            plus_err = ptolemy_defect(plus_pts)
            minus_err = ptolemy_defect(minus_pts)
            gradient_components.append((plus_err - minus_err) / (2 * eps))

        candidate_pts = list(current_pts)
        candidate_pts[2] -= QR([0, *(current_step * component for component in gradient_components)])
        candidate_err = ptolemy_defect(candidate_pts)

        if candidate_err <= err:
            current_pts = candidate_pts
            current_step = min(RR(step_size), current_step * RR(1.05))
        else:
            current_step *= RR(0.5)

    final_error = ptolemy_defect(current_pts)
    return {
        "history": history,
        "start_error": history[0],
        "end_error": float(final_error),
        "final_states": tuple(current_pts),
    }


def find_virtual_complement(
    triplet,
    iterations=80,
    step_size=0.15,
    gradient_eps=0.05,
):
    """
    Findet zu einem Primzahldrilling ein virtuelles viertes Quaternion D,
    das den Ptolemaeus-Defekt moeglichst klein macht.
    """
    A, B, C = [QR(state.coefficient_tuple()) for state in cluster_to_quaternions(triplet)]
    current_D = A - B + C

    directions = (
        QR([1, 0, 0, 0]),
        QR([0, 1, 0, 0]),
        QR([0, 0, 1, 0]),
        QR([0, 0, 0, 1]),
    )
    current_step = RR(step_size)
    eps = RR(gradient_eps)
    history = []

    for _ in range(iterations):
        current_error = ptolemy_defect((A, B, C, current_D))
        history.append(float(current_error))

        gradient_components = []
        for direction in directions:
            plus_D = current_D + eps * direction
            minus_D = current_D - eps * direction
            plus_error = ptolemy_defect((A, B, C, plus_D))
            minus_error = ptolemy_defect((A, B, C, minus_D))
            gradient_components.append((plus_error - minus_error) / (2 * eps))

        candidate_D = current_D - QR([current_step * component for component in gradient_components])
        candidate_error = ptolemy_defect((A, B, C, candidate_D))

        if candidate_error <= current_error:
            current_D = candidate_D
            current_step = min(RR(step_size), current_step * RR(1.05))
        else:
            current_step *= RR(0.5)

    final_error = ptolemy_defect((A, B, C, current_D))
    return {
        "triplet": triplet,
        "A": A,
        "B": B,
        "C": C,
        "D_virt": current_D,
        "virtual_energy": current_D.reduced_norm(),
        "start_error": history[0],
        "end_error": float(final_error),
        "history": history,
    }


def analysiere_vierling(vierling, start=14.0, ende=30.0, schritte=100):
    """Berechnet Ptolemaeus-Defekt und Zeta-Minimum fuer einen Vierling."""
    zustaende = cluster_to_quaternions(vierling)
    p_sum, d_prod, diff = ptolemy_check(*zustaende)
    optimal_e, min_res, _ = scan_energieniveaus(start, ende, schritte, zustaende)
    return {
        "vierling": vierling,
        "zustaende": zustaende,
        "ptolemy_sum": p_sum,
        "diagonal_prod": d_prod,
        "ptolemy_diff": diff,
        "optimal_e": optimal_e,
        "min_res": min_res,
        "cycle_work": calculate_cycle_work(vierling),
    }


def ranglisten_fuer_vierlinge(vierlinge, top_n=3):
    """Erzeugt Ranglisten nach Ptolemaeus-Defekt und Zeta-Resonanz."""
    analysen = [analysiere_vierling(vierling) for vierling in vierlinge]
    nach_ptolemaeus = sorted(analysen, key=lambda item: (item["ptolemy_diff"], item["vierling"]))
    nach_resonanz = sorted(analysen, key=lambda item: (item["min_res"], item["vierling"]))
    return nach_ptolemaeus[:top_n], nach_resonanz[:top_n]


def drucke_rangliste(titel, eintraege, schluessel):
    """Gibt eine kompakte Rangliste fuer Vierlinge aus."""
    print(titel)
    for index, eintrag in enumerate(eintraege, start=1):
        print(
            f"{index}. {eintrag['vierling']} | "
            f"{schluessel}: {eintrag[schluessel].n()} | "
            f"E*: {eintrag['optimal_e'].n()}"
        )


def main():
    A, B, C, D = ZUSTAENDE
    p_sum, d_prod, error = ptolemy_check(A, B, C, D)
    optimal_e, min_res, _ = scan_energieniveaus()
    triplets, quads = find_prime_clusters()
    erster_drilling = triplets[0]
    erster_vierling = quads[0]
    beste_ptolemaeus, beste_resonanz = ranglisten_fuer_vierlinge(quads, top_n=3)
    fav_ptolemy = beste_ptolemaeus[0]
    fav_zeta = beste_resonanz[0]
    resilienz_ptolemy = simulate_resilience(fav_ptolemy["vierling"])
    resilienz_zeta = simulate_resilience(fav_zeta["vierling"])
    starr_8 = beste_ptolemaeus[0]["zustaende"] + beste_ptolemaeus[1]["zustaende"]
    frei_4 = fav_zeta["zustaende"]
    konfiguration_8plus4 = analysiere_8plus4_konfiguration(starr_8, frei_4)
    virtuelles_komplement = find_virtual_complement(triplets[-1])
    natuerlicher_vierling_phase = calculate_phase_shift(fav_ptolemy["zustaende"])
    geheilter_drilling_phase = calculate_phase_shift(
        tuple(QR(state.coefficient_tuple()) for state in cluster_to_quaternions(triplets[-1]))
        + (virtuelles_komplement["D_virt"],)
    )

    print(f"Ptolemaeus-Summe (ac + bd): {p_sum.n()}")
    print(f"Diagonal-Produkt (ef):     {d_prod.n()}")
    print(f"Abweichung (Entropie-Mass): {error.n()}")
    print(f"Optimales Energieniveau fuer den Zirkel: {optimal_e.n()}")
    print(f"Minimale Zeta-Amplitude:               {min_res.n()}")
    print()
    print("### Kanon-Konfiguration fuer #Energiedoku ###")
    print(f"Erster Drillings-Kanon: {erster_drilling} -> eabc-Resonanz aktiv")
    print(f"100. Drillings-Kanon:   {triplets[-1]}")
    print(
        f"Vier-Quadrate-Darstellung von {erster_drilling[0]}: "
        f"{quat_representation(erster_drilling[0])}"
    )
    print("-" * 40)
    print(f"Erster Vierlings-Kanon: {erster_vierling} -> Ptolemaeus-stabil")
    print(f"50. Vierlings-Kanon:    {quads[-1]}")
    print(
        f"Vier-Quadrate-Darstellung von {erster_vierling[0]}: "
        f"{quat_representation(erster_vierling[0])}"
    )
    print()
    drucke_rangliste("Top-3 Vierlinge nach kleinstem Ptolemaeus-Defekt:", beste_ptolemaeus, "ptolemy_diff")
    print()
    drucke_rangliste("Top-3 Vierlinge nach staerkster Zeta-Resonanz:", beste_resonanz, "min_res")
    print()
    print(f"Arbeitsflaeche (Ptolemaeus-Favorit): {fav_ptolemy['cycle_work'].n()}")
    print(f"Arbeitsflaeche (Zeta-Favorit):     {fav_zeta['cycle_work'].n()}")
    print()
    print(
        "Resilienz-Check (Ptolemaeus-Favorit): "
        f"Start-Fehler {resilienz_ptolemy['start_error']:.4f}, "
        f"End-Fehler {resilienz_ptolemy['end_error']:.4f}"
    )
    print(
        "Resilienz-Check (Zeta-Favorit):      "
        f"Start-Fehler {resilienz_zeta['start_error']:.4f}, "
        f"End-Fehler {resilienz_zeta['end_error']:.4f}"
    )
    print()
    print("8+4 Energiedoku-Konfiguration:")
    print(f"Rahmen-Energie: {konfiguration_8plus4['rahmen_energie'].n()}")
    print(f"Resonanzsumme:  {konfiguration_8plus4['resonanz'].n()}")
    print(f"Effizienz-Index: {konfiguration_8plus4['effizienz_index'].n()}")
    print()
    print(f"Virtuelles Komplement fuer Drilling {virtuelles_komplement['triplet']}:")
    print(f"Virtuelle Quaternione: {virtuelles_komplement['D_virt']}")
    print(
        "Benötigtes Energieniveau (virtuelle Norm): "
        f"{float(virtuelles_komplement['virtual_energy']):.2f}"
    )
    print(
        "Ptolemaeus-Defekt: "
        f"{virtuelles_komplement['start_error']:.4f} -> "
        f"{virtuelles_komplement['end_error']:.4f}"
    )
    print()
    print(f"Phasenverschiebung Natürlicher Vierling: {float(natuerlicher_vierling_phase):.6f} rad")
    print(f"Phasenverschiebung Geheilter Drilling:   {float(geheilter_drilling_phase):.6f} rad")
    print()
    phys_mu = RR("1836.15267")
    eabc_index = konfiguration_8plus4["effizienz_index"]
    diff = phys_mu - eabc_index
    phase_corr = abs(geheilter_drilling_phase) / (2 * pi)
    zeta_corr = RR(1) / konfiguration_8plus4["resonanz"]
    multiplikativ_index = eabc_index * (1 + phase_corr / 100)
    additiv_phase_index = eabc_index + phase_corr

    if phase_corr == 0:
        alpha_phase = RR(0)
    else:
        alpha_phase = diff / phase_corr
    kalibrierter_phase_index = eabc_index + alpha_phase * phase_corr

    combined_corr = phase_corr * zeta_corr
    if combined_corr == 0:
        beta_combined = RR(0)
    else:
        beta_combined = diff / combined_corr
    kalibrierter_combined_index = eabc_index + beta_combined * combined_corr

    print("Korrekturmodelle fuer die eabc-Effizienz:")
    print(f"Original Index: {eabc_index.n()}")
    print(f"Zielwert (Proton/Elektron): {phys_mu.n()}")
    print(f"Urspruengliche Differenz: {diff.n()}")
    print(f"Multiplikativer Phasenterm: {multiplikativ_index.n()}")
    print(f"Restdifferenz multiplikativ: {(phys_mu - multiplikativ_index).n()}")
    print(f"Additiver Phasenterm (unkalibriert): {additiv_phase_index.n()}")
    print(f"Restdifferenz additiv: {(phys_mu - additiv_phase_index).n()}")
    print(f"Kalibrierter Phasenkoeffizient alpha: {alpha_phase.n()}")
    print(f"Kalibriertes Phasenmodell: {kalibrierter_phase_index.n()}")
    print(f"Kombinierte Observable Phase/Zeta: {combined_corr.n()}")
    print(f"Kalibrierter Kombinationskoeffizient beta: {beta_combined.n()}")
    print(f"Kalibriertes Kombinationsmodell: {kalibrierter_combined_index.n()}")
    print()
    batch_stats = batch_process_kanon(limit=1000, observable=combined_corr, mu_target=phys_mu)
    print("### STATISTISCHER KANON-CHECK ###")
    print(f"Anzahl Cluster: {batch_stats['count']}")
    print(f"Mittelwert Beta: {batch_stats['mean_beta']:.4f}")
    print(f"Standardabweichung: {batch_stats['std_beta']:.4f}")
    print(f"Minimum Beta: {batch_stats['min_beta']:.4f}")
    print(f"Maximum Beta: {batch_stats['max_beta']:.4f}")
    print()
    outlier_info = analyze_outlier(batch_stats["beta_distribution"], batch_stats["clusters"])
    kanon_cluster = fav_ptolemy["vierling"][:3]
    kanon_stress_matrix = calculate_stress_matrix(kanon_cluster)
    stress_matrix = calculate_stress_matrix(outlier_info["cluster"])
    stress_plot_path = plot_stress_comparison(
        kanon_stress_matrix,
        stress_matrix,
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "stress_comparison.png"),
    )
    print("### KOLLAPS-ANALYSE: DER INSTABILSTE DRILLING ###")
    print(f"Idealer Kanon (Stress-Referenz): {kanon_cluster}")
    print(f"Cluster: {outlier_info['cluster']}")
    print(
        f"Abweichendes Beta: {outlier_info['beta']:.4f} "
        f"(Delta: {outlier_info['delta']:.4f})"
    )
    print(f"Determinante der Stress-Matrix: {stress_matrix.det().n()}")
    print(f"Heatmap gespeichert unter: {stress_plot_path}")
    print()
    alpha_inv = RR("137.035999")
    kopplung = diff / alpha_inv
    print("Alpha-Kopplung im eabc-Modell:")
    print(f"Inverses Feinstrukturmass: {alpha_inv.n()}")
    print(f"Benoetigter Korrekturfaktor pro Alpha-Einheit: {kopplung.n()}")
    print()
    alpha_inv_pred = (eabc_index / diff) * (log(phys_mu) / phys_mu) * RR(137)
    alpha_inv_geo = (eabc_index / RR(12)) - (diff * pi)
    print("Alpha-Vorhersage im eabc-Modell:")
    print(f"Bamberger Index (Struktur): {eabc_index.n()}")
    print(f"Vorhergesagtes Alpha^-1 (kalibriert): {alpha_inv_pred.n()}")
    print(f"Vorhergesagtes Alpha^-1 (geometrisch): {alpha_inv_geo.n()}")
    print(f"Referenzwert (CODATA): {alpha_inv.n()}")
    print(f"Abweichung geometrisch: {(alpha_inv - alpha_inv_geo).n()}")
    print()
    alpha_calc = validate_alpha_coupling(eabc_index, phys_mu)
    print("Bamberger Eich-Modul:")
    print(f"Theoretisches Alpha^-1: {alpha_calc.n()}")
    print(f"Abweichung zum CODATA-Wert: {(alpha_inv - alpha_calc).n()}")


if __name__ == "__main__":
    main()