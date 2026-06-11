# =================================================================
# ARITHMETISCHER HURWITZ-EXPLORER (ORIGINAL HURWITZ)
# Implementierung: Textsuche + originale Hurwitz-Register/Orbits
# =================================================================

import os
import math
import hashlib
import shutil
import subprocess
from pathlib import Path
from fractions import Fraction
from dataclasses import dataclass

DEFAULT_ZETTEL_DIRS = (
    Path.home() / "Obsidian" / "Zettel",
    Path.home() / "Documents" / "Obsidian" / "Zettel",
    Path.home() / "Library" / "Mobile Documents" / "iCloud~md~obsidian" / "Documents" / "Zettel",
    Path.home() / "Library" / "Mobile Documents" / "iCloud~md~obsidian" / "Documents" / "Obsidian" / "Zettel",
)


def resolve_default_zettel_dir():
    """Nimmt den ersten existierenden Standardpfad fuer die Obsidian-Zettel."""
    for candidate in DEFAULT_ZETTEL_DIRS:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "Kein Standard-Zettelverzeichnis gefunden. Geprueft wurden: "
        + ", ".join(str(path) for path in DEFAULT_ZETTEL_DIRS)
    )


def file_contains_text(file_path, search_term):
    """Prueft robust, ob ein Text in einer Datei vorkommt."""
    try:
        content = Path(file_path).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return False
    return search_term.casefold() in content.casefold()


def first_match_excerpt(file_path, search_term, max_len=140):
    """Erste passende Textstelle aus einer Datei."""
    try:
        lines = Path(file_path).read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError:
        return ""

    needle = search_term.casefold()
    for line in lines:
        if needle in line.casefold():
            text = line.strip()
            return text[: max_len - 3] + "..." if len(text) > max_len else text
    return ""


def rank_text_matches(target_dir, all_files, search_term):
    """
    Ranked Trefferliste fuer Inhaltsuche.

    Sortierung:
    1. mehr Vorkommen
    2. Suchtext auch im Dateinamen
    3. fruehere erste Trefferposition
    4. kuerzerer Dateiname
    5. alphabetischer Dateiname
    """
    needle = search_term.casefold()
    matches = []

    for idx, file_name in enumerate(all_files):
        file_path = Path(target_dir) / file_name
        try:
            content = file_path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue

        haystack = content.casefold()
        occurrence_count = haystack.count(needle)
        if occurrence_count == 0:
            continue

        first_position = haystack.find(needle)
        file_name_bonus = 1 if needle in file_name.casefold() else 0
        excerpt = first_match_excerpt(file_path, search_term)
        matches.append(
            SearchMatch(
                index=idx,
                file_name=file_name,
                occurrence_count=occurrence_count,
                first_position=first_position,
                file_name_bonus=file_name_bonus,
                excerpt=excerpt,
            )
        )

    matches.sort(
        key=lambda match: (
            -match.occurrence_count,
            -match.file_name_bonus,
            match.first_position,
            len(match.file_name),
            match.file_name.casefold(),
        )
    )
    return matches


def match_label(match, rank):
    """Kurzes Label für die interaktive Trefferauswahl."""
    safe_name = sanitize_dialog_text(match.file_name)
    return f"{rank}. {safe_name} | Treffer={match.occurrence_count}"


def sanitize_dialog_text(text):
    """Entfernt problematische Zeichen für osascript/Dialogtexte."""
    cleaned = text.replace("\x00", " ")
    cleaned = cleaned.replace("\r", " ").replace("\n", " ").replace("\t", " ")
    cleaned = "".join(ch if ch.isprintable() else " " for ch in cleaned)
    return " ".join(cleaned.split())


def choose_match_interactively(ranked_matches):
    """Trefferauswahl per Klick auf macOS, sonst per Konsolen-Fallback."""
    if len(ranked_matches) <= 1:
        return ranked_matches[0]

    labels = [match_label(match, i) for i, match in enumerate(ranked_matches, start=1)]
    if shutil.which("osascript"):
        escaped = ', '.join(
            '"' + sanitize_dialog_text(label).replace("\\", "\\\\").replace('"', '\\"') + '"'
            for label in labels
        )
        script = (
            f'set theChoices to {{{escaped}}}\n'
            'set selectedItem to choose from list theChoices with prompt "Treffer auswählen" '
            'default items {item 1 of theChoices}\n'
            'if selectedItem is false then\n'
            '    return ""\n'
            'end if\n'
            'return item 1 of selectedItem\n'
        )
        result = subprocess.run(
            ["osascript", "-e", script],
            capture_output=True,
            text=True,
            check=False,
        )
        selected_label = result.stdout.strip()
        if selected_label:
            return ranked_matches[labels.index(selected_label)]

    print("\nMehrere Treffer gefunden. Bitte Auswahl eingeben:")
    for label in labels:
        print(f"  {label}")
    choice = input("Nummer des gewünschten Treffers [1]: ").strip()
    if not choice:
        return ranked_matches[0]
    try:
        index = max(1, min(len(ranked_matches), int(choice))) - 1
    except ValueError:
        index = 0
    return ranked_matches[index]


def show_selected_match_dialog(match):
    """Zeigt nach der Auswahl die zugehörige Textstelle in einem Dialog."""
    if not shutil.which("osascript"):
        return

    title = sanitize_dialog_text(match.file_name)
    excerpt = sanitize_dialog_text(match.excerpt or "Keine Textstelle verfügbar.")
    message = (
        f"Datei: {title}\n"
        f"Trefferzahl: {match.occurrence_count}\n"
        f"Erste Position: {match.first_position}\n\n"
        f"Textstelle:\n{excerpt}"
    )
    escaped_message = message.replace("\\", "\\\\").replace('"', '\\"')
    script = (
        f'display dialog "{escaped_message}" '
        f'with title "Ausgewählter Treffer" buttons {{"OK"}} default button "OK"'
    )
    subprocess.run(
        ["osascript", "-e", script],
        capture_output=True,
        text=True,
        check=False,
    )


def open_note_in_obsidian(file_path):
    """Öffnet den ausgewählten Zettel direkt in Obsidian."""
    file_path = str(Path(file_path).expanduser())
    if shutil.which("open"):
        result = subprocess.run(
            ["open", "-a", "Obsidian", file_path],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode == 0:
            return True

        subprocess.run(
            ["open", file_path],
            capture_output=True,
            text=True,
            check=False,
        )
        return False
    return False


@dataclass(frozen=True)
class HurwitzQubit:
    """Ein Bitpaar als hurwitzsch interpretierte Qubit-Zelle."""
    pair_index: int
    bit_pair: str
    basis_unit: tuple


@dataclass(frozen=True)
class HurwitzRegister:
    """Instanziierter Hurwitz-Registerzustand für einen Dateikandidaten."""
    index: int
    file_name: str
    bitstring: str
    qubits: tuple
    lipschitz_state: tuple
    hurwitz_phase: tuple
    hurwitz_state: tuple
    norm: Fraction
    orbit_rep: tuple
    signature: tuple


@dataclass(frozen=True)
class HurwitzOracle:
    """Instanziiertes Oracle für eine Hurwitz-Signaturklasse."""
    target_register: HurwitzRegister
    marked_indices: tuple

    def describe(self):
        """Kurze textuelle Beschreibung des Oracles."""
        return {
            "target_file": self.target_register.file_name,
            "target_bits": self.target_register.bitstring,
            "signature": self.target_register.signature,
            "marked_indices": list(self.marked_indices),
        }


@dataclass(frozen=True)
class SearchMatch:
    """Gerankter Treffer in einem Zettel."""
    index: int
    file_name: str
    occurrence_count: int
    first_position: int
    file_name_bonus: int
    excerpt: str


@dataclass(frozen=True)
class HurwitzExplorer:
    """Vollständiger Explorer-Zustand mit Hurwitz-Registern und Oracle."""
    target_dir: Path
    search_term: str
    all_files: tuple
    ranked_matches: tuple
    matched_indices: tuple
    target_idx: int
    n_qubits: int
    target_register: HurwitzRegister
    hurwitz_oracle: HurwitzOracle

    @classmethod
    def from_search_term(cls, target_dir, search_term):
        """Erzeugt einen Explorer aus Verzeichnis und Suchbegriff."""
        target_dir = Path(target_dir).expanduser()
        if not target_dir.exists():
            raise FileNotFoundError(f"Pfad existiert nicht: {target_dir}")

        all_files = tuple(
            file_name
            for file_name in os.listdir(target_dir)
            if os.path.isfile(target_dir / file_name)
        )
        if not all_files:
            raise ValueError("Ordner ist leer.")

        ranked_matches = tuple(rank_text_matches(target_dir, all_files, search_term))
        if not ranked_matches:
            raise ValueError("Suchtext in keinem Zettel gefunden.")

        matched_indices = tuple(match.index for match in ranked_matches)
        target_idx = matched_indices[0]

        n_qubits = int(math.ceil(math.log2(len(all_files))))
        target_register = instantiate_hurwitz_register(
            target_idx, n_qubits, all_files[target_idx]
        )
        hurwitz_oracle = instantiate_hurwitz_oracle(
            all_files,
            n_qubits,
            target_idx,
            marked_indices=matched_indices,
        )
        return cls(
            target_dir=target_dir,
            search_term=search_term,
            all_files=all_files,
            ranked_matches=ranked_matches,
            matched_indices=matched_indices,
            target_idx=target_idx,
            n_qubits=n_qubits,
            target_register=target_register,
            hurwitz_oracle=hurwitz_oracle,
        )


def format_scalar(x):
    """Formatiert rationale Quaternionen-Koeffizienten lesbar."""
    if isinstance(x, Fraction):
        if x.denominator == 1:
            return str(x.numerator)
        return f"{x.numerator}/{x.denominator}"
    return str(x)


def format_quaternion(q):
    """Gibt q = a + bi + cj + dk als lesbaren String aus."""
    labels = ["", "i", "j", "k"]
    parts = []
    for value, label in zip(q, labels):
        if value == 0:
            continue
        if label == "":
            parts.append(format_scalar(value))
        elif value == 1:
            parts.append(label)
        elif value == -1:
            parts.append(f"-{label}")
        else:
            parts.append(f"{format_scalar(value)}{label}")

    if not parts:
        return "0"

    result = parts[0]
    for part in parts[1:]:
        if part.startswith("-"):
            result += " - " + part[1:]
        else:
            result += " + " + part
    return result


def Q(a, b, c, d):
    """Quaternion a + bi + cj + dk mit rationalen Koeffizienten."""
    return (Fraction(a), Fraction(b), Fraction(c), Fraction(d))


def q_mul(x, y):
    """Hamilton-Produkt zweier Quaternionen."""
    a, b, c, d = x
    e, f, g, h = y
    return (
        a * e - b * f - c * g - d * h,
        a * f + b * e + c * h - d * g,
        a * g - b * h + c * e + d * f,
        a * h + b * g - c * f + d * e,
    )


def q_norm(x):
    """Quaternionen-Norm a^2 + b^2 + c^2 + d^2."""
    a, b, c, d = x
    return a * a + b * b + c * c + d * d


def hurwitz_units():
    """Die 24 Hurwitz-Einheiten der Norm 1."""
    units = {
        Q(1, 0, 0, 0), Q(-1, 0, 0, 0),
        Q(0, 1, 0, 0), Q(0, -1, 0, 0),
        Q(0, 0, 1, 0), Q(0, 0, -1, 0),
        Q(0, 0, 0, 1), Q(0, 0, 0, -1),
    }

    for s1 in (-1, 1):
        for s2 in (-1, 1):
            for s3 in (-1, 1):
                for s4 in (-1, 1):
                    units.add(
                        Q(
                            Fraction(s1, 2),
                            Fraction(s2, 2),
                            Fraction(s3, 2),
                            Fraction(s4, 2),
                        )
                    )

    return sorted(units)


HURWITZ_UNITS = hurwitz_units()
HALF_UNITS = [u for u in HURWITZ_UNITS if any(value.denominator == 2 for value in u)]
BITPAIR_TO_UNIT = {
    "00": Q(1, 0, 0, 0),
    "01": Q(0, 1, 0, 0),
    "10": Q(0, 0, 1, 0),
    "11": Q(0, 0, 0, 1),
}


def choose_hurwitz_phase(file_name):
    """Ordnet einer Datei deterministisch eine halbzahlig-hurwitzsche Phase zu."""
    digest = hashlib.sha256(file_name.encode("utf-8")).digest()
    return HALF_UNITS[int.from_bytes(digest[:2], "big") % len(HALF_UNITS)]


def build_hurwitz_state(index, n_qubits, file_name):
    """
    Konstruiert eine Hurwitz-artige Repräsentation eines Qubit-Basiszustands.

    Die Bitpaare werden auf 1, i, j, k abgebildet; eine zusätzliche
    halbzahlig-hurwitzsche Einheit codiert eine dateispezifische Phase.
    """
    bitstring = format(index, f"0{n_qubits}b")
    padded_bits = bitstring if len(bitstring) % 2 == 0 else "0" + bitstring
    bit_pairs = [padded_bits[i:i + 2] for i in range(0, len(padded_bits), 2)]

    word_units = [BITPAIR_TO_UNIT[pair] for pair in bit_pairs]
    lipschitz_state = Q(1, 0, 0, 0)
    for unit in word_units:
        lipschitz_state = q_mul(lipschitz_state, unit)

    hurwitz_phase = choose_hurwitz_phase(file_name)
    hurwitz_state = q_mul(lipschitz_state, hurwitz_phase)

    return {
        "bitstring": bitstring,
        "bit_pairs": bit_pairs,
        "word_units": word_units,
        "lipschitz_state": lipschitz_state,
        "hurwitz_phase": hurwitz_phase,
        "hurwitz_state": hurwitz_state,
        "norm": q_norm(hurwitz_state),
    }


def q_key(q):
    """Hashbare Normalform einer rationalen Quaternion."""
    return tuple((value.numerator, value.denominator) for value in q)


def left_orbit(x):
    """Linksorbit eines Quaternionenzustands unter den 24 Hurwitz-Einheiten."""
    return {q_mul(unit, x) for unit in HURWITZ_UNITS}


def instantiate_hurwitz_register(index, n_qubits, file_name):
    """
    Erzeugt einen vollständigen Hurwitz-Registerzustand für einen Dateikandidaten.
    """
    info = build_hurwitz_state(index, n_qubits, file_name)
    orbit_rep = min(q_key(candidate) for candidate in left_orbit(info["hurwitz_state"]))
    signature = (orbit_rep, q_key(info["hurwitz_phase"]))
    qubits = tuple(
        HurwitzQubit(pair_index=i, bit_pair=pair, basis_unit=unit)
        for i, (pair, unit) in enumerate(zip(info["bit_pairs"], info["word_units"]))
    )
    return HurwitzRegister(
        index=index,
        file_name=file_name,
        bitstring=info["bitstring"],
        qubits=qubits,
        lipschitz_state=info["lipschitz_state"],
        hurwitz_phase=info["hurwitz_phase"],
        hurwitz_state=info["hurwitz_state"],
        norm=info["norm"],
        orbit_rep=orbit_rep,
        signature=signature,
    )


def hurwitz_signature(index, n_qubits, file_name):
    """Liefert die Hurwitz-Signatur eines instanziierten Registers."""
    register = instantiate_hurwitz_register(index, n_qubits, file_name)
    return {
        "state": {
            "bitstring": register.bitstring,
            "lipschitz_state": register.lipschitz_state,
            "hurwitz_phase": register.hurwitz_phase,
            "hurwitz_state": register.hurwitz_state,
            "norm": register.norm,
        },
        "orbit_rep": register.orbit_rep,
        "phase_key": q_key(register.hurwitz_phase),
        "signature": register.signature,
        "register": register,
    }


def matching_hurwitz_indices(all_files, n_qubits, target_idx):
    """Alle Indizes mit derselben Hurwitz-Signatur wie der Zieltreffer."""
    target_sig = hurwitz_signature(target_idx, n_qubits, all_files[target_idx])["signature"]
    matches = []
    for idx, file_name in enumerate(all_files):
        sig = hurwitz_signature(idx, n_qubits, file_name)["signature"]
        if sig == target_sig:
            matches.append(idx)
    return matches


def instantiate_hurwitz_oracle(all_files, n_qubits, target_idx, marked_indices=None):
    """Instanziiert ein Oracle fuer die Hurwitz-Klasse des Zieltreffers."""
    target_register = instantiate_hurwitz_register(target_idx, n_qubits, all_files[target_idx])
    if marked_indices is None:
        marked_indices = tuple(matching_hurwitz_indices(all_files, n_qubits, target_idx))
    return HurwitzOracle(
        target_register=target_register,
        marked_indices=tuple(marked_indices),
    )


def print_hurwitz_preview(all_files, n_qubits, target_idx, max_rows=8):
    """Zeigt eine kleine Hurwitz-Vorschau der ersten Kandidaten und des Ziels."""
    preview_count = min(len(all_files), max_rows)
    print("\n[HURWITZ] Quaternionische Kodierung der Kandidaten:")
    for idx in range(preview_count):
        register = instantiate_hurwitz_register(idx, n_qubits, all_files[idx])
        marker = " <== Treffer" if idx == target_idx else ""
        print(
            f"  {idx:>4} | {register.bitstring} | "
            f"{format_quaternion(register.hurwitz_state)} | "
            f"OrbitRep={register.orbit_rep}{marker}"
        )

    if target_idx >= preview_count:
        register = instantiate_hurwitz_register(target_idx, n_qubits, all_files[target_idx])
        print("  ...")
        print(
            f"  {target_idx:>4} | {register.bitstring} | "
            f"{format_quaternion(register.hurwitz_state)} | "
            f"OrbitRep={register.orbit_rep} <== Treffer"
        )

    target_register = instantiate_hurwitz_register(target_idx, n_qubits, all_files[target_idx])
    matches = matching_hurwitz_indices(all_files, n_qubits, target_idx)
    print("\n[HURWITZ] Zielzustand:")
    print(f"  Datei: {target_register.file_name}")
    print(f"  Bits: {target_register.bitstring}")
    print(f"  Bitpaare: {[qubit.bit_pair for qubit in target_register.qubits]}")
    print(f"  Instanziierte Hurwitz-Qubits: {len(target_register.qubits)}")
    print(f"  Lipschitz-Teil: {format_quaternion(target_register.lipschitz_state)}")
    print(f"  Hurwitz-Phase: {format_quaternion(target_register.hurwitz_phase)}")
    print(f"  Gesamtzustand: {format_quaternion(target_register.hurwitz_state)}")
    print(f"  Norm: {target_register.norm}")
    print(f"  Orbit-Repräsentant: {target_register.orbit_rep}")
    print(f"  Treffer via Hurwitz-Signatur: {matches[:8]}{' ...' if len(matches) > 8 else ''}")


# --- HAUPTPROGRAMM ---

def run_explorer():
    print("=== ARITHMETISCHER HURWITZ-EXPLORER ===")
    search_term = input("Suchtext: ").strip()
    if not search_term:
        print("Fehler: Kein Suchtext eingegeben.")
        return

    try:
        target_dir = resolve_default_zettel_dir()
        print(f"Standard-Zettelverzeichnis: {target_dir}")
        explorer = HurwitzExplorer.from_search_term(target_dir, search_term)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Fehler: {exc}")
        return

    print(f"Gefundene Elemente im Suchraum: {len(explorer.all_files)}")
    print(f"Hurwitz-Registerbreite: {explorer.n_qubits} Bit")
    print(f"Hurwitz-Zielzustand: {format_quaternion(explorer.target_register.hurwitz_state)}")

    print(f"\nTreffer: {len(explorer.matched_indices)} Datei(en)")
    for pos, match in enumerate(explorer.ranked_matches[:10], start=1):
        print(
            f"  {pos}. {match.file_name} "
            f"(Treffer={match.occurrence_count}, Pos={match.first_position})"
        )
        if match.excerpt:
            print(f"     {match.excerpt}")
    if len(explorer.ranked_matches) > 10:
        print(f"  ... und {len(explorer.ranked_matches) - 10} weitere")

    selected_match = choose_match_interactively(explorer.ranked_matches)
    show_selected_match_dialog(selected_match)
    selected_path = explorer.target_dir / explorer.all_files[selected_match.index]
    opened_in_obsidian = open_note_in_obsidian(selected_path)
    selected_register = instantiate_hurwitz_register(
        selected_match.index,
        explorer.n_qubits,
        explorer.all_files[selected_match.index],
    )
    selected_oracle = instantiate_hurwitz_oracle(
        explorer.all_files,
        explorer.n_qubits,
        selected_match.index,
        marked_indices=explorer.matched_indices,
    )

    print("\n[ORIGINAL HURWITZ]")
    print(f"  Ziel-Datei: {selected_oracle.target_register.file_name}")
    print(
        f"  Hurwitz-markierte Indizes: "
        f"{list(selected_oracle.marked_indices[:12])}"
        f"{' ...' if len(selected_oracle.marked_indices) > 12 else ''}"
    )
    print(f"  Orbit-Repräsentant: {selected_register.orbit_rep}")

    print("\nAusgewählter Treffer:")
    print(f"  Datei: {selected_register.file_name}")
    print(f"  Geöffnet in Obsidian: {'ja' if opened_in_obsidian else 'Fallback/Dateiöffnung'}")
    print(f"  Trefferzahl: {selected_match.occurrence_count}")
    print(f"  Erste Position: {selected_match.first_position}")
    if selected_match.excerpt:
        print(f"  Textstelle: {selected_match.excerpt}")

# Starten
if __name__ == "__main__":
    run_explorer()