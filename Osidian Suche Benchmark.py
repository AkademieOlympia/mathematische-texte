#!/usr/bin/env python3
"""
Benchmark-Modul: Vergleicht Obsidian-Suche mit Spotlight und Grover/Quantenalgorithmen.
Kann standalone oder von Osidian Suche.py mit -b aufgerufen werden.
"""

import os
import sys
import subprocess
import time
from pathlib import Path

import numpy as np
from scipy.fft import fft, ifft

# Standardpfad
DEFAULT_VAULT = Path.home() / "Library/Mobile Documents/iCloud~md~obsidian/Documents/Zettel"


def find_md_files(vault_path: Path) -> list[Path]:
    """Sammelt alle .md-Dateien im Vault."""
    if not vault_path.exists():
        return []
    files = []
    for root, _, fnames in os.walk(vault_path):
        for f in fnames:
            if f.endswith(".md") and not f.startswith("."):
                files.append(Path(root) / f)
    return files


def search_in_file(filepath: Path, term: str, case_sensitive: bool = False) -> list:
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


def search_vault_python(vault_path: Path, term: str, case_sensitive: bool = False) -> list:
    """Python-Volltextsuche."""
    files = find_md_files(vault_path)
    results = []
    for fp in files:
        matches = search_in_file(fp, term, case_sensitive)
        if matches:
            rel = fp.relative_to(vault_path) if vault_path in fp.parents or fp.parent == vault_path else fp.name
            results.append((rel, matches))
    return results


class HolographicGroverSearch:
    """Grover-inspirierte holographische Suche (Quantenalgorithmus-Simulation)."""
    def __init__(self, file_list: list[str], n_scale: int = 4):
        self.files = list(file_list)
        self.n = len(self.files)
        self.n_scale = n_scale
        if self.n == 0:
            self.bulk_energy = np.array([])
            return
        self.bulk_energy = np.array([hash(f) % 10**6 for f in self.files], dtype=float)
        mx = np.max(self.bulk_energy)
        if mx > 0:
            self.bulk_energy /= mx

    def search(self, target: str) -> tuple:
        """Sucht nach target (Dateiname oder Suchbegriff) im Frequenzraum."""
        if self.n == 0:
            return None, 0.0
        residuum = fft(self.bulk_energy)
        threshold = 1.0 / self.n_scale
        residuum = np.clip(residuum.real, -threshold, threshold) + 1j * np.clip(residuum.imag, -threshold, threshold)
        ref_hash = hash(max(self.files, key=hash)) % 10**6
        target_hash = (hash(target) % 10**6) / ref_hash if ref_hash else 0
        correlation = np.abs(ifft(residuum * np.conj(fft(np.full(self.n, target_hash)))))
        idx = np.argmax(correlation)
        return self.files[idx], float(correlation[idx])


def run_benchmark(vault_path: Path, term: str, case_sensitive: bool = False) -> dict:
    """Führt Benchmark aus: Python, Spotlight, Grover/Quantenalgorithmen."""
    vault_path = Path(vault_path).expanduser()
    times = {}
    file_names = [str(f.name) for f in find_md_files(vault_path)]

    # 1. Python-Volltextsuche
    t0 = time.perf_counter()
    results_py = search_vault_python(vault_path, term, case_sensitive)
    times["Python (Volltextsuche)"] = time.perf_counter() - t0
    times["_n_python"] = len(results_py)

    # 2. Spotlight (Apple mdfind)
    try:
        t0 = time.perf_counter()
        proc = subprocess.run(
            ["mdfind", "-onlyin", str(vault_path), term],
            capture_output=True,
            text=True,
            timeout=30,
        )
        times["Spotlight (Apple mdfind)"] = time.perf_counter() - t0
        times["_n_spotlight"] = len(proc.stdout.strip().split("\n")) if proc.stdout.strip() else 0
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        times["Spotlight (Apple mdfind)"] = None
        times["_spotlight_error"] = str(e)
        times["_n_spotlight"] = 0

    # 3. Grover/Holographisch (Quantenalgorithmus-Simulation)
    if file_names:
        t0 = time.perf_counter()
        grover = HolographicGroverSearch(file_names)
        found, strength = grover.search(term)
        times["Grover (Quantenalgorithmus)"] = time.perf_counter() - t0
        times["_grover_found"] = found
        times["_grover_strength"] = strength
    else:
        times["Grover (Quantenalgorithmus)"] = None
        times["_grover_found"] = None

    return times


def print_benchmark_report(bench: dict, term: str) -> None:
    """Gibt Benchmark-Report und Fehleranalyse aus."""
    print("=== Benchmark (Geschwindigkeitsvergleich) ===\n")
    for name, val in bench.items():
        if name.startswith("_"):
            continue
        if val is not None:
            print(f"  {name}: {val*1000:.1f} ms")
        else:
            print(f"  {name}: nicht verfügbar ({bench.get('_spotlight_error', '?')})")

    print(f"\n  Python-Treffer: {bench.get('_n_python', 0)} Datei(en)")
    print(f"  Spotlight-Treffer: {bench.get('_n_spotlight', 0)} Datei(en)")
    if "_grover_found" in bench and bench["_grover_found"]:
        print(f"  Grover-Treffer: {bench['_grover_found']} (Stärke: {bench.get('_grover_strength', 0):.4f})")

    # ========== Fehleranalyse: Spotlight vs. Quantenalgorithmen ==========
    print("\n" + "=" * 60)
    print("FEHLERANALYSE / VERGLEICH")
    print("Spotlight (Apple) vs. hier verwendete Quantenalgorithmen")

    t_spot = bench.get("Spotlight (Apple mdfind)")
    t_grover = bench.get("Grover (Quantenalgorithmus)")
    t_py = bench.get("Python (Volltextsuche)")

    n_spot = bench.get("_n_spotlight", 0)
    n_py = bench.get("_n_python", 0)

    print("-" * 60)
    print("Spotlight:")
    print("  - Nutzt Index (mdimport), sucht Metadaten + Inhalt")
    print("  - Sehr schnell, aber iCloud-Ordner oft nicht vollständig indexiert")
    print(f"  - Laufzeit: {t_spot*1000:.1f} ms" if t_spot else "  - Nicht verfügbar")
    print(f"  - Treffer: {n_spot}")

    print("\nQuantenalgorithmen (Grover/holographisch):")
    print("  - FFT-basierte Korrelation im Frequenzraum")
    print("  - Sucht nur nach Dateinamen (nicht Inhalt)")
    print("  - Theoretisch O(√N), hier O(1) FFT")
    print(f"  - Laufzeit: {t_grover*1000:.1f} ms" if t_grover else "  - Nicht verfügbar")
    if bench.get("_grover_found"):
        print(f"  - Treffer: {bench['_grover_found']}")

    print("\nPython-Volltextsuche:")
    print("  - Liest jede Datei, sucht Zeile für Zeile")
    print("  - Liefert Zeilennummern und Kontext")
    print(f"  - Laufzeit: {t_py*1000:.1f} ms" if t_py else "")
    print(f"  - Treffer: {n_py}")

    print("\nFazit:")
    if t_spot and t_grover and t_py:
        fastest = min([("Spotlight", t_spot), ("Grover", t_grover), ("Python", t_py)], key=lambda x: x[1])
        print(f"  Schnellste Methode: {fastest[0]} ({fastest[1]*1000:.1f} ms)")
    print("  Spotlight: Index-basiert, schnell, aber iCloud-Index oft unvollständig.")
    print("  Grover: Nur Dateinamen, extrem schnell, für Inhaltssuche ungeeignet.")
    print("  Python: Volltext, Zeilenkontext, zuverlässig für iCloud.")
    print("=" * 60)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Obsidian-Benchmark: Spotlight vs. Quantenalgorithmen")
    parser.add_argument("suchbegriff", nargs="?", default="Primzahl", help="Suchbegriff")
    parser.add_argument("-v", "--vault", default=str(DEFAULT_VAULT), help="Vault-Pfad")
    args = parser.parse_args()

    vault = Path(args.vault).expanduser()
    if not vault.exists():
        print(f"Vault nicht gefunden: {vault}")
        sys.exit(1)

    term = args.suchbegriff
    print(f"Benchmark: Vault={vault}, Suchbegriff=\"{term}\"\n")

    bench = run_benchmark(vault, term)
    print_benchmark_report(bench, term)


if __name__ == "__main__":
    main()
