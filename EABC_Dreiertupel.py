import math
import sys
from itertools import product

try:
    import tkinter as tk
    from tkinter import simpledialog, messagebox
except Exception:  # pragma: no cover - Fallback fuer Umgebungen ohne Tk
    tk = None
    simpledialog = None
    messagebox = None


FAMILY_BY_RESIDUE = {
    1: "E",
    5: "A",
    7: "B",
    11: "C",
}

PATTERN_BY_SIZE = {
    2: ("E", "A"),
    3: ("E", "A", "B"),
    4: ("E", "A", "B", "C"),
    5: ("E", "A", "B", "C", "E"),
}


def ask_limit(default=200):
    """Liest N sichtbar ueber die Konsole, Dialog nur optional per --dialog."""
    prefer_dialog = "--dialog" in sys.argv

    if prefer_dialog and tk is not None and simpledialog is not None:
        root = tk.Tk()
        root.withdraw()
        try:
            value = simpledialog.askinteger(
                "EABC-Dreiertupel",
                "Obergrenze N eingeben:",
                initialvalue=default,
                minvalue=5,
            )
        finally:
            root.destroy()
        if value is None:
            raise SystemExit("Abgebrochen.")
        return value

    print("Obergrenze N eingeben.", flush=True)
    raw = input(f"N [{default}]: ").strip()
    if not raw:
        return default
    return int(raw)


def should_show_plot():
    """Erlaubt einen reinen Konsolenlauf per --no-plot."""
    return "--no-plot" not in sys.argv


def parse_preview(default=20):
    """Liest optional --preview N oder --preview=N aus den Argumenten."""
    args = sys.argv[1:]
    for index, arg in enumerate(args):
        if arg.startswith("--preview="):
            value = arg.split("=", 1)[1].strip()
            if not value:
                raise SystemExit("Bitte einen Wert fuer --preview angeben.")
            preview = int(value)
            if preview < 1:
                raise SystemExit("--preview muss >= 1 sein.")
            return preview

        if arg == "--preview":
            if index + 1 >= len(args):
                raise SystemExit("Bitte einen Wert nach --preview angeben.")
            preview = int(args[index + 1])
            if preview < 1:
                raise SystemExit("--preview muss >= 1 sein.")
            return preview

    return default


def sieve(limit):
    """Berechnet alle Primzahlen <= limit."""
    if limit < 2:
        return []

    prime_mask = [True] * (limit + 1)
    prime_mask[0] = prime_mask[1] = False
    for p in range(2, math.isqrt(limit) + 1):
        if prime_mask[p]:
            for multiple in range(p * p, limit + 1, p):
                prime_mask[multiple] = False
    return [n for n, is_prime in enumerate(prime_mask) if is_prime]


def classify_primes(primes):
    """Ordnet Primzahlen den EABC-Klassen via mod 12 zu."""
    families = {"E": [], "A": [], "B": [], "C": []}
    for p in primes:
        if p <= 3:
            continue
        residue = p % 12
        family = FAMILY_BY_RESIDUE.get(residue)
        if family is not None:
            families[family].append(p)
    return families


def evenly_sample(values, target_count):
    """Waehlt gleichmaessig ueber die Liste verteilte Werte aus."""
    if len(values) <= target_count:
        return values
    if target_count <= 1:
        return [values[-1]]

    step = (len(values) - 1) / (target_count - 1)
    indices = []
    seen = set()
    for i in range(target_count):
        index = round(i * step)
        if index not in seen:
            indices.append(index)
            seen.add(index)
    return [values[index] for index in indices]


def prepare_pattern_families(families, pattern, max_combinations):
    """Reduziert grosse Suchraeume gleichmaessig auf eine handhabbare Stichprobe."""
    selected = {family: list(families[family]) for family in pattern}
    estimated = 1
    for family in pattern:
        estimated *= len(selected[family])

    if estimated <= max_combinations:
        return selected, False, estimated

    unique_families = tuple(dict.fromkeys(pattern))
    per_family_target = max(2, int(max_combinations ** (1 / len(pattern))))

    for family in unique_families:
        selected[family] = evenly_sample(selected[family], per_family_target)

    estimated = 1
    for family in pattern:
        estimated *= len(selected[family])
    return selected, True, estimated


def fermat_two_squares_prime(prime):
    """Findet fuer p == 1 mod 4 eine Darstellung p = a^2 + b^2."""
    if prime <= 1 or prime % 4 != 1:
        return None

    for a in range(math.isqrt(prime), 0, -1):
        remainder = prime - a * a
        b = math.isqrt(remainder)
        if b * b == remainder:
            return tuple(sorted((a, b), reverse=True))
    return None


def annotate_three_group(item):
    """Reduziert Emissionen auf drei Endgruppen nach aktiven A/B/C-Familien."""
    families = tuple(item["pattern"])
    values = item["values"]
    active_abc = tuple(family for family in "ABC" if family in families)
    a_values = [value for family, value in zip(families, values) if family == "A"]
    a_value = a_values[0] if a_values else None

    annotated = dict(item)
    annotated["final_group"] = len(active_abc)
    annotated["final_group_name"] = "".join(active_abc) if active_abc else "leer"
    annotated["distinct_prime_count"] = len(set(values))
    annotated["a_value"] = a_value
    annotated["a_fermat"] = (
        fermat_two_squares_prime(a_value) if a_value is not None else None
    )
    return annotated


def build_emissions(families, max_combinations=250_000):
    """Erzeugt Emissionen mit 2 bis 5 Zahlen gemaess den EABC-Mustern."""
    emissions = []
    used_sampling = False

    for size, pattern in PATTERN_BY_SIZE.items():
        print(f"Berechne Muster {''.join(pattern)}...", flush=True)
        prepared, sampled, estimated = prepare_pattern_families(
            families, pattern, max_combinations
        )
        used_sampling = used_sampling or sampled
        value_lists = [prepared[family] for family in pattern]

        for values in product(*value_lists):
            e_value = math.prod(values)
            if cube_root_digit_group(e_value) != size:
                continue

            emissions.append(
                {
                    "size": size,
                    "pattern": "".join(pattern),
                    "values": values,
                    "e_value": e_value,
                    "root_floor": cube_root_floor(e_value),
                    "inverse_e": inverse_value(e_value),
                    "inverse_cuberoot": inverse_cube_root(e_value),
                    "sampled": sampled,
                    "estimated": estimated,
                }
            )

        found = sum(1 for item in emissions if item["size"] == size)
        info = " (mit Stichprobe)" if sampled else ""
        print(
            f"  -> {found} Emissionen fuer Gruppe {size}{info}",
            flush=True,
        )

    emissions = [annotate_three_group(item) for item in emissions]
    emissions.sort(key=lambda item: (item["final_group"], -item["e_value"]))
    return emissions, used_sampling


def cube_root(value):
    """Berechnet die reelle dritte Wurzel."""
    return value ** (1.0 / 3.0)


def cube_root_floor(value):
    """Ganzzahliger Anteil der reellen dritten Wurzel."""
    return int(cube_root(value))


def cube_root_digit_group(value):
    """Gruppiert nach der Stellenzahl von int(cbrt(E))."""
    return len(str(cube_root_floor(value)))


def inverse_cube_root(value):
    """Berechnet das Inverse der reellen dritten Wurzel."""
    return 1.0 / cube_root(value)


def inverse_value(value):
    """Berechnet das Inverse eines positiven Wertes."""
    return 1.0 / value


def print_summary(limit, families, emissions, used_sampling, preview=20):
    print(f"\nPrimzahlen bis N = {limit}")
    print(f"E: {len(families['E'])} Elemente")
    print(f"A: {len(families['A'])} Elemente")
    print(f"B: {len(families['B'])} Elemente")
    print(f"C: {len(families['C'])} Elemente")
    print(f"\nAnzahl der Emissionen: {len(emissions)}")
    print("\nEmissionen nach Gruppengroesse:")
    for size in range(2, 6):
        count = sum(1 for item in emissions if item["size"] == size)
        pattern = "".join(PATTERN_BY_SIZE[size])
        print(f"{size} Zahlen ({pattern}): {count}")

    print("\nReduzierte Drei-Gruppen-Liste:")
    for final_group in range(1, 4):
        count = sum(1 for item in emissions if item["final_group"] == final_group)
        name = {1: "A", 2: "AB", 3: "ABC"}[final_group]
        print(f"{final_group} aktive ABC-Familien ({name}): {count}")

    if used_sampling:
        print(
            "\nHinweis: Fuer grosse Suchraeume wurde gleichmaessig gestichprobt, "
            "damit die Emissionen berechenbar bleiben."
        )

    print("\nAuslese nach den drei Endgruppen:")
    for final_group in range(1, 4):
        name = {1: "A", 2: "AB", 3: "ABC"}[final_group]
        group_items = [item for item in emissions if item["final_group"] == final_group]
        if not group_items:
            continue

        print(f"\nGruppe {final_group} ({name}) - groesste Emissionen:")
        for index, item in enumerate(group_items[:preview], start=1):
            values_text = ", ".join(str(value) for value in item["values"])
            fermat_text = ""
            if item["a_value"] is not None and item["a_fermat"] is not None:
                a_left, a_right = item["a_fermat"]
                fermat_text = (
                    f", A={item['a_value']}={a_left}^2+{a_right}^2"
                )
            print(
                f"{index:2d}. {item['pattern']}=({values_text}) -> "
                f"E={item['e_value']}, int(cbrt(E))={item['root_floor']}, "
                f"Originalgruppe={item['size']}, Endgruppe={item['final_group_name']}, "
                f"verschiedene Primzahlen={item['distinct_prime_count']}, "
                f"1/cbrt(E)={item['inverse_cuberoot']:.6f}{fermat_text}"
            )


def plot_e_axis(emissions, max_points=300):
    """Traegt die inverse Darstellung der groessten Emissionen ein."""
    if not emissions:
        return

    import matplotlib.pyplot as plt

    shown = emissions[:max_points]
    inverse_e_values = [item["inverse_e"] for item in shown]
    inverse_cube_roots = [item["inverse_cuberoot"] for item in shown]
    final_groups = [item["final_group"] for item in shown]
    x_positions = range(1, len(shown) + 1)
    group_colors = {
        1: "forestgreen",
        2: "darkorange",
        3: "darkred",
    }

    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.plot(x_positions, inverse_e_values, linewidth=0.8, color="steelblue", alpha=0.65)
    for final_group in range(1, 4):
        xs = [x for x, group in zip(x_positions, final_groups) if group == final_group]
        ys = [y for y, group in zip(inverse_e_values, final_groups) if group == final_group]
        if xs:
            ax1.scatter(
                xs,
                ys,
                s=18,
                color=group_colors[final_group],
                label=f"Endgruppe {final_group}",
            )
    ax1.set_title("Inverse Darstellung der EABC-Emissionen")
    ax1.set_xlabel("Index der nach Endgruppen sortierten Emissionen")
    ax1.set_ylabel("1/E", color="darkblue")
    ax1.tick_params(axis="y", labelcolor="darkblue")
    ax1.legend(loc="upper right")

    ax2 = ax1.twinx()
    ax2.plot(x_positions, inverse_cube_roots, linewidth=1.0, color="darkred", alpha=0.8)
    ax2.set_ylabel("1/cbrt(E)", color="darkred")
    ax2.tick_params(axis="y", labelcolor="darkred")

    plt.tight_layout()
    plt.show()


def main():
    preview = parse_preview()
    limit = ask_limit()
    print(f"Starte Berechnung bis N={limit}...", flush=True)
    primes = sieve(limit)
    families = classify_primes(primes)
    emissions, used_sampling = build_emissions(families)
    print_summary(limit, families, emissions, used_sampling, preview=preview)
    if should_show_plot():
        plot_e_axis(emissions)

    if tk is not None and messagebox is not None:
        root = tk.Tk()
        root.withdraw()
        try:
            messagebox.showinfo(
                "EABC-Dreiertupel",
                f"Fertig. {len(emissions)} Emissionen wurden erzeugt.",
            )
        finally:
            root.destroy()


if __name__ == "__main__":
    main()
