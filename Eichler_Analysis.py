from __future__ import annotations

from collections import Counter
from itertools import permutations

import matplotlib.pyplot as plt

from Eichler import (
    SIGMAS,
    State,
    local_shell_energy_complex,
    neighbor_states,
    print_run,
    schalenbrenner_run,
    schalenbrenner_run_ext,
)

DEFAULT_COMPLEXITY_T_VALUES = tuple(0.5 + 0.25 * i for i in range(13))


def analyze_temperatures(
    temperatures: list[float] | tuple[float, ...] = (0.01, 0.1, 0.35, 0.7, 1.5),
    shells: int = 20,
) -> list[dict[str, object]]:
    print(f"\nThermodynamische Analyse ({shells} Schalen)")
    print(f"{'Temp':>6} | {'Avg Energy':>10} | {'Eichler-Ratio':>15} | Final State")
    print("-" * 70)
    results: list[dict[str, object]] = []

    for t in temperatures:
        history = schalenbrenner_run(temperature=t, shells=shells)
        avg_energy = sum(r.energy for r in history) / len(history)

        # Eichler-Ratio: perfekte Permutation von (E, A, B, C).
        hits = sum(1 for r in history if all(Counter(r.state)[s] == 1 for s in SIGMAS))
        ratio = hits / len(history)

        final_state = history[-1].state
        print(f"{t:6.2f} | {avg_energy:10.3f} | {ratio:15.2%} | {final_state}")
        results.append(
            {
                "temperature": t,
                "avg_energy": avg_energy,
                "eichler_ratio": ratio,
                "final_state": final_state,
            }
        )

    return results


def plot_temperature_analysis(
    results: list[dict[str, object]],
    output_path: str = "Eichler_Temperatur_Analyse.png",
) -> None:
    temperatures = [float(entry["temperature"]) for entry in results]
    avg_energies = [float(entry["avg_energy"]) for entry in results]
    ratios = [100.0 * float(entry["eichler_ratio"]) for entry in results]

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    axes[0].plot(temperatures, avg_energies, marker="o", linewidth=2, color="tab:blue")
    axes[0].set_ylabel("Mittlere Energie")
    axes[0].set_title("Eichler-Schalenbrenner: Temperaturanalyse")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(temperatures, ratios, marker="s", linewidth=2, color="tab:orange")
    axes[1].set_xlabel("Temperatur")
    axes[1].set_ylabel("Eichler-Ratio [%]")
    axes[1].set_ylim(0, 105)
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    print(f"\nPlot gespeichert: {output_path}")


def find_critical_point(t_min: float = 0.1, t_max: float = 2.5, steps: int = 25) -> None:
    if steps <= 1:
        raise ValueError("steps muss groesser als 1 sein.")

    print("\nBestimmung der kritischen Masse (Phasenuebergang)")
    print(f"{'Temp T':>8} | {'Var(|axis|)':>12} | {'Integritaet':>12} | {'Zustand'}")
    print("-" * 65)

    for i in range(steps):
        t = t_min + (t_max - t_min) * (i / (steps - 1))
        history = schalenbrenner_run(temperature=t, shells=100, seed=42)

        norms = [r.axis_norm for r in history]
        avg_norm = sum(norms) / len(norms)
        variance = sum((n - avg_norm) ** 2 for n in norms) / len(norms)

        # Integritaet = Anteil der Schalen mit allen vier Basiselementen.
        integrity = sum(1 for r in history if len(set(r.state)) == 4) / len(history)

        if variance < 1e-5:
            status = "FEST (Eichler-Kristall)"
        elif variance < 0.05:
            status = "FLUESSIG (Uebergang)"
        else:
            status = "GASFOERMIG (Chaos)"

        print(f"{t:8.2f} | {variance:12.5f} | {integrity:12.1%} | {status}")


def analyze_specific_heat(
    temperatures: list[float] | tuple[float, ...] = (0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 15.0),
) -> None:
    print("\nThermische Analyse: Spezifische Waerme & Phasenuebergang")
    print(f"{'Temp':>6} | {'Avg Energy':>10} | {'Spec Heat (Cv)':>15} | {'Eichler-Ratio'}")
    print("-" * 70)

    for t in temperatures:
        history = schalenbrenner_run(temperature=t, shells=100, seed=42)

        energies = [r.energy for r in history]
        avg_e = sum(energies) / len(energies)

        # Varianz der Energie fuer die effektive Waermekapazitaet.
        var_e = sum((e - avg_e) ** 2 for e in energies) / len(energies)
        cv = var_e / (t ** 2)

        hits = sum(1 for r in history if len(set(r.state)) == 4)
        ratio = hits / len(history)

        print(f"{t:6.2f} | {avg_e:10.3f} | {cv:15.5f} | {ratio:14.2%}")


def analyze_complexity_impact(
    complexities: list[int] | tuple[int, ...] = (4, 8, 12),
    temperatures: list[float] | tuple[float, ...] = DEFAULT_COMPLEXITY_T_VALUES,
) -> None:
    print("\nKomplexitaets-Check: Verschiebung der kritischen Masse")
    print(f"{'Compl.':>7} | {'Temp':>6} | {'Spec Heat':>12} | {'Eichler-Ratio'}")
    print("-" * 60)

    for c in complexities:
        for t in temperatures:
            history = schalenbrenner_run_ext(temperature=t, shells=100, complexity=c)

            energies = [r.energy for r in history]
            avg_e = sum(energies) / len(energies)
            var_e = sum((e - avg_e) ** 2 for e in energies) / len(energies)
            cv = var_e / (t ** 2)

            ratio = sum(1 for r in history if len(set(r.state)) == 4) / len(history)
            print(f"{c:7d} | {t:6.1f} | {cv:12.5f} | {ratio:14.2%}")


def collect_transition_rows(
    complexity: int,
    t_values: list[float] | tuple[float, ...],
) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []

    for t in t_values:
        history = schalenbrenner_run_ext(temperature=t, shells=100, complexity=complexity)
        energies = [r.energy for r in history]
        avg_e = sum(energies) / len(energies)
        var_e = sum((e - avg_e) ** 2 for e in energies) / len(energies)
        cv = var_e / (t ** 2)
        ratio = sum(1 for r in history if len(set(r.state)) == 4) / len(history)
        rows.append(
            {
                "temperature": t,
                "avg_energy": avg_e,
                "var_energy": var_e,
                "cv": cv,
                "ratio": ratio,
            }
        )

    return rows


def estimate_critical_temperature_by_complexity(
    complexities: list[int] | tuple[int, ...] = (4, 8, 12),
    t_values: list[float] | tuple[float, ...] = DEFAULT_COMPLEXITY_T_VALUES,
) -> None:
    print("\nGeschaetzte kritische Temperatur pro Komplexitaet")
    print(f"{'Compl.':>7} | {'T_c via max Cv':>13} | {'Ratio < 90% ab':>15}")
    print("-" * 60)

    for c in complexities:
        rows = collect_transition_rows(c, t_values)
        tc_cv = max(rows, key=lambda row: row["cv"])["temperature"]
        tc_ratio = next((row["temperature"] for row in rows if row["ratio"] < 0.90), None)
        tc_ratio_text = f"{tc_ratio:.2f}" if tc_ratio is not None else "nicht erreicht"
        print(f"{c:7d} | {tc_cv:13.2f} | {tc_ratio_text:>15}")


def analyze_transition_indicators_by_complexity(
    complexities: list[int] | tuple[int, ...] = (4, 8, 12),
    t_values: list[float] | tuple[float, ...] = DEFAULT_COMPLEXITY_T_VALUES,
) -> None:
    print("\nRobuste Uebergangsindikatoren pro Komplexitaet")
    print(f"{'Compl.':>7} | {'Peak Var(E)':>11} | {'Staerkster Ratio-Fall':>22} | {'Delta Ratio':>11}")
    print("-" * 78)

    for c in complexities:
        rows = collect_transition_rows(c, t_values)
        peak_var_row = max(rows, key=lambda row: row["var_energy"])

        ratio_drop_index = None
        largest_drop = float("-inf")
        for i in range(len(rows) - 1):
            drop = rows[i]["ratio"] - rows[i + 1]["ratio"]
            if drop > largest_drop:
                largest_drop = drop
                ratio_drop_index = i

        if ratio_drop_index is None:
            drop_window = "nicht bestimmbar"
            drop_text = "0.00%"
        else:
            t_left = rows[ratio_drop_index]["temperature"]
            t_right = rows[ratio_drop_index + 1]["temperature"]
            drop_window = f"{t_left:.2f}->{t_right:.2f}"
            drop_text = f"{largest_drop:.2%}"

        print(
            f"{c:7d} | {peak_var_row['temperature']:11.2f} | "
            f"{drop_window:>22} | {drop_text:>11}"
        )


def analyze_local_transition_window(
    complexities: list[int] | tuple[int, ...] = (4, 8, 12),
    t_min: float = 2.0,
    t_max: float = 2.8,
    step: float = 0.05,
) -> None:
    if step <= 0:
        raise ValueError("step muss positiv sein.")
    if t_max < t_min:
        raise ValueError("t_max muss groesser oder gleich t_min sein.")

    count = int(round((t_max - t_min) / step)) + 1
    t_values = tuple(round(t_min + step * i, 10) for i in range(count))

    print("\nLokaler Feinscan des Uebergangsfensters")
    print(
        f"Fenster: T in [{t_min:.2f}, {t_max:.2f}] mit Schrittweite {step:.2f}"
    )
    print(f"{'Compl.':>7} | {'Peak Var(E)':>11} | {'Peak Cv':>8} | {'Ratio < 90% ab':>15} | {'Staerkster Fall':>18}")
    print("-" * 86)

    for c in complexities:
        rows = collect_transition_rows(c, t_values)
        peak_var_row = max(rows, key=lambda row: row["var_energy"])
        peak_cv_row = max(rows, key=lambda row: row["cv"])
        first_ratio_drop = next((row["temperature"] for row in rows if row["ratio"] < 0.90), None)

        ratio_drop_index = None
        largest_drop = float("-inf")
        for i in range(len(rows) - 1):
            drop = rows[i]["ratio"] - rows[i + 1]["ratio"]
            if drop > largest_drop:
                largest_drop = drop
                ratio_drop_index = i

        if ratio_drop_index is None:
            drop_window = "nicht bestimmbar"
        else:
            t_left = rows[ratio_drop_index]["temperature"]
            t_right = rows[ratio_drop_index + 1]["temperature"]
            drop_window = f"{t_left:.2f}->{t_right:.2f}"

        ratio_text = f"{first_ratio_drop:.2f}" if first_ratio_drop is not None else "nicht erreicht"
        print(
            f"{c:7d} | {peak_var_row['temperature']:11.2f} | {peak_cv_row['temperature']:8.2f} | "
            f"{ratio_text:>15} | {drop_window:>18}"
        )


def scan_energy_gap(state: State, complexity: int = 4) -> None:
    print(f"\nTQFT-Analyse: Vermessung der Energy Gap (Komplexitaet {complexity})")
    print(f"{'Zustands-Typ':<24} | {'Energie E':>10} | {'Abstand zur Basis (Gap)'}")
    print("-" * 72)

    base_energy = local_shell_energy_complex(state, shell=0, complexity=complexity)
    print(f"{'Grundzustand (ABCE)':<24} | {base_energy:10.3f} | {'REF (0.000)'}")

    all_energies = []
    for candidate in neighbor_states(state):
        if candidate == state:
            continue
        energy = local_shell_energy_complex(candidate, shell=0, complexity=complexity)
        all_energies.append((energy, candidate))

    all_energies.sort(key=lambda item: item[0])

    for energy, candidate in all_energies[:5]:
        gap = energy - base_energy
        if len(set(candidate)) == 4:
            type_str = "Symmetrie-Shift"
        else:
            type_str = "Defekt (Symmetriebruch)"
        print(f"{type_str:<24} | {energy:10.3f} | {gap:+10.3f}")

    first_positive_gap = next(
        ((energy, candidate) for energy, candidate in all_energies if energy - base_energy > 1e-9),
        None,
    )
    if first_positive_gap is None:
        print("\nErster echter Gap: nicht gefunden (lokal voll entartet).")
    else:
        energy, candidate = first_positive_gap
        gap = energy - base_energy
        if len(set(candidate)) == 4:
            type_str = "Symmetrie-Shift"
        else:
            type_str = "Defekt (Symmetriebruch)"
        print("\nErster echter positiver Gap:")
        print(f"{type_str:<24} | {energy:10.3f} | {gap:+10.3f} | Zustand {candidate}")


def is_even_permutation(base: State, permuted: State) -> bool:
    position = {value: idx for idx, value in enumerate(base)}
    values = [position[value] for value in permuted]

    inversions = 0
    for i in range(len(values)):
        for j in range(i + 1, len(values)):
            if values[i] > values[j]:
                inversions += 1

    return inversions % 2 == 0


def compute_long_run_stats(temperature: float, shells: int, complexity: int) -> dict[str, float]:
    history = schalenbrenner_run_ext(
        temperature=temperature,
        shells=shells,
        complexity=complexity,
        seed=42,
    )

    changes = 0
    state_persistence: list[int] = []
    current_streak = 1

    for i in range(1, len(history)):
        if history[i].state != history[i - 1].state:
            changes += 1
            state_persistence.append(current_streak)
            current_streak = 1
        else:
            current_streak += 1

    state_persistence.append(current_streak)
    avg_persistence = sum(state_persistence) / len(state_persistence) if state_persistence else 0.0

    base: State = ("E", "A", "B", "C")
    even_perms = {
        perm for perm in permutations(base) if is_even_permutation(base, perm)  # type: ignore[arg-type]
    }

    even_count = 0
    odd_count = 0
    balanced_count = 0
    for record in history:
        if len(set(record.state)) != 4:
            continue
        balanced_count += 1
        if record.state in even_perms:
            even_count += 1
        else:
            odd_count += 1

    bias_ratio = even_count / odd_count if odd_count > 0 else float("inf")
    return {
        "changes": float(changes),
        "avg_persistence": avg_persistence,
        "even_count": float(even_count),
        "odd_count": float(odd_count),
        "bias_ratio": bias_ratio,
        "eichler_ratio": balanced_count / len(history) if history else 0.0,
    }


def test_chebyshev_stabilization(shells: int = 2000, complexity: int = 12) -> None:
    print(f"\nLangzeit-Test: Dynamische Stabilisierung (Schalen: {shells})")
    print("Suche nach Oszillationen im Bereich der kritischen Masse (T=2.05)")
    print("-" * 75)

    stats = compute_long_run_stats(temperature=2.05, shells=shells, complexity=complexity)
    print(f"Anzahl der Zustandswechsel: {int(stats['changes'])}")
    print(f"Mittlere Verweildauer (Persistence): {stats['avg_persistence']:.2f} Schalen")
    print(f"Gerade Permutationen (4n+1 Typ): {int(stats['even_count'])}")
    print(f"Ungerade Permutationen (4n+3 Typ): {int(stats['odd_count'])}")
    print(
        f"Bias-Verhaeltnis: {stats['bias_ratio']:.4f}"
        if stats["odd_count"] > 0
        else "Bias-Verhaeltnis: inf"
    )


def compare_long_run_temperatures(
    temperatures: list[float] | tuple[float, ...] = (1.5, 2.05, 2.5),
    shells: int = 2000,
    complexity: int = 12,
) -> None:
    print(f"\nLangzeitvergleich ueber den Uebergangsbereich ({shells} Schalen, complexity={complexity})")
    print(f"{'Temp':>5} | {'Wechsel':>8} | {'Persistenz':>10} | {'Even':>6} | {'Odd':>6} | {'Bias':>8} | {'Eichler-Ratio':>13}")
    print("-" * 80)

    for temperature in temperatures:
        stats = compute_long_run_stats(temperature=temperature, shells=shells, complexity=complexity)
        bias_text = f"{stats['bias_ratio']:.4f}" if stats["odd_count"] > 0 else "inf"
        print(
            f"{temperature:5.2f} | {int(stats['changes']):8d} | {stats['avg_persistence']:10.2f} | "
            f"{int(stats['even_count']):6d} | {int(stats['odd_count']):6d} | {bias_text:>8} | "
            f"{stats['eichler_ratio']:13.2%}"
        )


if __name__ == "__main__":
    make_plot = False
    run_specific_heat = True
    run_complexity_impact = True
    run_transition_indicators = True
    run_local_transition_window = True
    run_energy_gap = False
    run_chebyshev_stabilization = False
    run_long_run_temperature_compare = True
    specific_heat_temperatures = [2.0, 4.0, 6.0, 10.0, 20.0]
    gap_state: State = ("E", "A", "B", "C")
    gap_complexity = 4

    # 1. Den detaillierten Standard-Run wie bisher
    run = schalenbrenner_run(temperature=0.35)
    print_run(run)

    # 2. Die neue System-Analyse
    results = analyze_temperatures()

    # Optionaler Plot fuer die thermodynamische Uebersicht.
    if make_plot:
        plot_temperature_analysis(results)

    # Optionale Zusatzanalyse mit echtem Schalter.
    # find_critical_point()
    if run_specific_heat:
        analyze_specific_heat(temperatures=specific_heat_temperatures)
    if run_complexity_impact:
        analyze_complexity_impact()
        estimate_critical_temperature_by_complexity()
    if run_transition_indicators:
        analyze_transition_indicators_by_complexity()
    if run_local_transition_window:
        analyze_local_transition_window()
    if run_energy_gap:
        scan_energy_gap(gap_state, complexity=gap_complexity)
    if run_chebyshev_stabilization:
        test_chebyshev_stabilization()
    if run_long_run_temperature_compare:
        compare_long_run_temperatures()
