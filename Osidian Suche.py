#!/usr/bin/env python3
"""
Obsidian-Suche: Durchsucht die gesamte Textmenge im Obsidian-Vault (Zettel auf iCloud).
Eingabe des Suchbegriffs per Kommandozeile oder interaktiv.

Verwendung:
  python3 "Osidian Suche.py" Suchbegriff
  python3 "Osidian Suche.py" -b Suchbegriff                # Mit Benchmark (Spotlight vs. Grover/Quanten)
  python3 "Osidian Suche.py" -s ergebnis.txt Suchbegriff   # Ergebnisse speichern
  python3 "Osidian Suche.py"                                # Interaktive Eingabe

  Benchmark separat: python3 "Osidian Suche Benchmark.py" Suchbegriff
"""

import os
import sys
import argparse
from pathlib import Path

# Standardpfad: Obsidian Zettel-Vault auf iCloud
DEFAULT_VAULT = Path.home() / "Library/Mobile Documents/iCloud~md~obsidian/Documents/Zettel"


def find_md_files(vault_path: Path) -> list[Path]:
    """Sammelt alle .md-Dateien im Vault (inkl. Unterordner)."""
    if not vault_path.exists():
        return []
    files = []
    for root, _, fnames in os.walk(vault_path):
        for f in fnames:
            if f.endswith(".md") and not f.startswith("."):
                files.append(Path(root) / f)
    return files


def search_in_file(filepath: Path, term: str, case_sensitive: bool = False) -> list[tuple[int, str]]:
    """Sucht nach term in einer Datei. Gibt (Zeile, Inhalt) zurück."""
    term_lower = term if case_sensitive else term.lower()
    matches = []
    try:
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            for i, line in enumerate(f, 1):
                check = line if case_sensitive else line.lower()
                if term_lower in check:
                    matches.append((i, line.rstrip()))
    except (OSError, UnicodeDecodeError):
        pass
    return matches


def search_vault(vault_path: Path, term: str, case_sensitive: bool = False) -> list[tuple[Path, list]]:
    """Durchsucht den gesamten Vault. Gibt [(Datei, [(Zeile, Inhalt), ...]), ...] zurück."""
    files = find_md_files(vault_path)
    results = []
    for fp in files:
        matches = search_in_file(fp, term, case_sensitive)
        if matches:
            rel = fp.relative_to(vault_path) if vault_path in fp.parents or fp.parent == vault_path else fp.name
            results.append((rel, matches))
    return results


def main():
    parser = argparse.ArgumentParser(description="Obsidian-Vault durchsuchen (ganze Textmenge)")
    parser.add_argument("suchbegriff", nargs="?", help="Suchbegriff (optional, sonst interaktiv)")
    parser.add_argument("-v", "--vault", default=str(DEFAULT_VAULT), help="Pfad zum Obsidian-Vault")
    parser.add_argument("-s", "--save", metavar="DATEI", help="Ergebnisse in Datei speichern")
    parser.add_argument("-b", "--benchmark", action="store_true", help="Benchmark: Vergleich mit Spotlight (Apple-Suche)")
    parser.add_argument("-i", "--case-sensitive", action="store_true", help="Groß-/Kleinschreibung beachten")
    args = parser.parse_args()

    vault = Path(args.vault).expanduser()
    if not vault.exists():
        print(f"Vault nicht gefunden: {vault}")
        sys.exit(1)

    term = args.suchbegriff
    if not term:
        term = input("Suchbegriff eingeben: ").strip()
    if not term:
        print("Kein Suchbegriff angegeben.")
        sys.exit(0)

    print(f"Durchsuche Vault: {vault}")
    print(f"Suchbegriff: {term}\n")

    results = search_vault(vault, term, args.case_sensitive)

    if args.benchmark:
        try:
            from importlib.util import spec_from_file_location, module_from_spec
            bench_path = Path(__file__).parent / "Osidian Suche Benchmark.py"
            spec = spec_from_file_location("benchmark", bench_path)
            bench_mod = module_from_spec(spec)
            spec.loader.exec_module(bench_mod)
            bench = bench_mod.run_benchmark(vault, term, args.case_sensitive)
            bench_mod.print_benchmark_report(bench, term)
        except Exception as e:
            print(f"Benchmark fehlgeschlagen: {e}")
        print()

    lines_out = []
    lines_out.append(f"=== Obsidian-Suche: \"{term}\" ===\n")
    lines_out.append(f"Vault: {vault}\n")
    lines_out.append(f"Treffer: {len(results)} Datei(en)\n")

    for rel_path, matches in results:
        lines_out.append(f"\n--- {rel_path} ({len(matches)} Treffer) ---")
        for line_no, content in matches[:10]:  # max 10 Zeilen pro Datei
            snippet = content[:120] + "..." if len(content) > 120 else content
            lines_out.append(f"  Z{line_no}: {snippet}")
        if len(matches) > 10:
            lines_out.append(f"  ... und {len(matches) - 10} weitere")
        lines_out.append("")

    text = "\n".join(lines_out)
    print(text)

    if args.save:
        with open(args.save, "w", encoding="utf-8") as f:
            f.write(text)
        print(f"Ergebnisse gespeichert: {args.save}")


if __name__ == "__main__":
    main()
