#!/usr/bin/env python3
"""
Grover II: Erweiterte Grover-Version für Inhaltssuche (Content Search).
Sucht im Dateiinhalt, nicht nur in Dateinamen. Nutzt Grover-inspirierte
Amplitudenverstärkung im Frequenzraum.
"""

import numpy as np
from scipy.fft import fft, ifft
from pathlib import Path


class GroverContentSearch:
    """
    Grover-inspirierte Inhaltssuche.
    - Liest Dateiinhalte
    - Berechnet Relevanz-Scores (Suchbegriff-Vorkommen)
    - Amplitude-Amplification im Frequenzraum (Grover-Iteration)
    - Liefert sortierte Treffer mit Zeilenkontext
    """

    def __init__(self, file_paths: list[Path], vault_root: Path | None = None):
        self.file_paths = list(file_paths)
        self.vault_root = vault_root
        self.n = len(self.file_paths)
        self._content_cache: dict[int, str] = {}

    def _read_content(self, idx: int) -> str:
        """Liest und cached den Dateiinhalt."""
        if idx in self._content_cache:
            return self._content_cache[idx]
        try:
            with open(self.file_paths[idx], "r", encoding="utf-8", errors="replace") as f:
                content = f.read()
        except (OSError, UnicodeDecodeError):
            content = ""
        self._content_cache[idx] = content
        return content

    def _relevance_scores(self, term: str, case_sensitive: bool = False) -> np.ndarray:
        """Berechnet Relevanz pro Datei: Vorkommen des Suchbegriffs (normalisiert 0..1)."""
        term_lower = term if case_sensitive else term.lower()
        scores = np.zeros(self.n, dtype=float)
        for i in range(self.n):
            content = self._read_content(i)
            check = content if case_sensitive else content.lower()
            count = check.count(term_lower)
            # Normierung: 1 + log(1+count) für abnehmenden Zusatznutzen
            scores[i] = np.log1p(count) if count > 0 else 0.0
        mx = np.max(scores)
        if mx > 0:
            scores /= mx
        return scores

    def _grover_amplify(self, scores: np.ndarray, iterations: int = 1) -> np.ndarray:
        """
        Grover-Amplitudenverstärkung: hebt Treffer hervor.
        Oracle: Markiert Treffer durch Phasenflip.
        Diffusion: Inversion um den Mittelwert.
        """
        if self.n == 0 or iterations < 1:
            return scores
        # Oracle: Treffer (score > 0) bekommen negatives Vorzeichen
        amplitudes = np.where(scores > 0, -scores, scores)
        for _ in range(iterations - 1):
            mean = np.mean(amplitudes)
            amplitudes = 2 * mean - amplitudes
        mean = np.mean(amplitudes)
        amplified = 2 * mean - amplitudes
        return np.abs(amplified)

    def _fft_amplify(self, scores: np.ndarray) -> np.ndarray:
        """FFT-basierte Verstärkung im Frequenzraum (holographische Variante)."""
        if self.n == 0:
            return scores
        spec = fft(scores)
        # Betonung der dominierenden Frequenzen
        spec *= (1.0 + 0.5 * np.abs(spec) / (np.max(np.abs(spec)) + 1e-12))
        return np.abs(ifft(spec).real)

    def search(
        self,
        term: str,
        case_sensitive: bool = False,
        top_k: int = 20,
        grover_iterations: int = 2,
        use_fft: bool = True,
    ) -> list[tuple[Path, float, list[tuple[int, str]]]]:
        """
        Sucht nach term im Dateiinhalt.
        Gibt [(rel_path, full_path, Score, [(Zeile, Inhalt), ...]), ...] zurück, sortiert nach Score.
        """
        scores = self._relevance_scores(term, case_sensitive)
        if grover_iterations > 0:
            scores = self._grover_amplify(scores, grover_iterations)
        if use_fft:
            scores = self._fft_amplify(scores)

        # Top-K Indizes
        order = np.argsort(scores)[::-1]
        results = []
        term_lower = term if case_sensitive else term.lower()
        for idx in order[:top_k]:
            if scores[idx] <= 0:
                break
            fp = self.file_paths[idx]
            content = self._read_content(idx)
            lines = content.split("\n")
            matches = []
            for i, line in enumerate(lines, 1):
                check = line if case_sensitive else line.lower()
                if term_lower in check:
                    matches.append((i, line.strip()[:150]))
            if matches:
                try:
                    rel = fp.relative_to(self.vault_root) if self.vault_root else fp.name
                except ValueError:
                    rel = fp.name
                results.append((Path(rel), fp, float(scores[idx]), matches[:5]))
        return results


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


# Import für find_md_files
import os
import subprocess
import sys
import threading

DEFAULT_VAULT = Path.home() / "Library/Mobile Documents/iCloud~md~obsidian/Documents/Zettel"


def open_file(path: Path) -> None:
    """Öffnet die Datei (Obsidian, falls .md)."""
    try:
        if sys.platform == "darwin":
            subprocess.run(["open", str(path)], check=False)
        elif sys.platform == "win32":
            os.startfile(str(path))
        else:
            subprocess.run(["xdg-open", str(path)], check=False)
    except Exception:
        pass


def run_interactive():
    """Interaktive Terminal-Oberfläche (Fallback wenn tkinter fehlt)."""
    vault = Path(DEFAULT_VAULT).expanduser()
    if not vault.exists():
        print(f"Vault nicht gefunden: {vault}")
        return
    files = find_md_files(vault)
    grover = GroverContentSearch(files, vault)
    print(f"Grover II – Interaktive Suche ({len(files)} Notizen)")
    print("Suchbegriff eingeben, Enter. Zahl 1–N für Öffnen. Leer = Ende.\n")
    while True:
        try:
            term = input("Suchbegriff: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nTschüss.")
            break
        if not term:
            break
        results = grover.search(term, top_k=25)
        if not results:
            print("Keine Treffer.\n")
            continue
        print(f"\n{len(results)} Treffer:\n")
        for i, (rel, full_path, score, matches) in enumerate(results, 1):
            print(f"  [{i}] {rel} (Score: {score:.3f})")
            for ln, snip in matches[:2]:
                print(f"      Z{ln}: {(snip[:80] + '...') if len(snip) > 80 else snip}")
            print()
        try:
            choice = input("Zahl zum Öffnen (1–{}), Enter = weiter: ".format(len(results))).strip()
            if choice:
                idx = int(choice)
                if 1 <= idx <= len(results):
                    open_file(results[idx - 1][1])
                    print("→ Geöffnet.\n")
        except (ValueError, KeyboardInterrupt):
            pass
        print()


def run_gui():
    """Dialog-orientierte Oberfläche: Suche mit klickbaren Ergebnissen."""
    try:
        import tkinter as tk
        from tkinter import ttk, scrolledtext, messagebox
    except ImportError:
        print("tkinter nicht verfügbar (pyenv-Python oft ohne Tcl/Tk).")
        print("→ Starte interaktive Terminal-Oberfläche.\n")
        run_interactive()
        return

    vault = Path(DEFAULT_VAULT).expanduser()
    if not vault.exists():
        messagebox.showerror("Fehler", f"Vault nicht gefunden:\n{vault}")
        return

    files = find_md_files(vault)
    grover = GroverContentSearch(files, vault)

    root = tk.Tk()
    root.title("Grover II – Obsidian Suche")
    root.geometry("700x600")
    root.minsize(500, 400)

    # Suchzeile
    f_search = ttk.Frame(root, padding=8)
    f_search.pack(fill=tk.X)
    ttk.Label(f_search, text="Suchbegriff:").pack(side=tk.LEFT, padx=(0, 0))
    entry = ttk.Entry(f_search, width=40)
    entry.pack(side=tk.LEFT, padx=4, fill=tk.X, expand=True)
    entry.focus()
    entry.bind("<Return>", lambda e: do_search())

    status = ttk.Label(f_search, text=f"  {len(files)} Notizen geladen")
    status.pack(side=tk.LEFT)

    def do_search():
        term = entry.get().strip()
        if not term:
            return
        status.config(text=" Suche… ")
        root.update()

        def search_thread():
            res = grover.search(term, top_k=25)
            root.after(0, lambda: show_results(res))

        threading.Thread(target=search_thread, daemon=True).start()

    def show_results(results):
        result_area.delete("1.0", tk.END)
        if not results:
            result_area.insert(tk.END, "Keine Treffer.")
            status.config(text=" Keine Treffer ")
            return
        status.config(text=f" {len(results)} Treffer ")
        for rel, full_path, score, matches in results:
            result_area.insert(tk.END, f"\n{rel} (Score: {score:.3f})\n", "header")
            for ln, snip in matches:
                result_area.insert(tk.END, f"  Z{ln}: {snip}\n", "snippet")
            result_area.insert(tk.END, "  ", "link")
            tag = f"open_{id(full_path)}"
            result_area.insert(tk.END, "▶ In Obsidian öffnen\n", tag)
            result_area.tag_config(tag, foreground="blue", underline=True)
            result_area.tag_bind(tag, "<Button-1>", lambda e, p=full_path: open_file(p))
            result_area.tag_bind(tag, "<Enter>", lambda e, t=tag: result_area.tag_config(t, foreground="darkblue"))
            result_area.tag_bind(tag, "<Leave>", lambda e, t=tag: result_area.tag_config(t, foreground="blue"))
        result_area.see("1.0")

    ttk.Button(f_search, text="Suchen", command=do_search).pack(side=tk.LEFT, padx=4)

    # Ergebnisbereich
    result_area = scrolledtext.ScrolledText(root, wrap=tk.WORD, font=("Helvetica", 11), padx=8, pady=8)
    result_area.pack(fill=tk.BOTH, expand=True)
    result_area.tag_config("header", font=("Helvetica", 11, "bold"))
    result_area.tag_config("snippet", foreground="gray30")

    result_area.insert(tk.END, "Suchbegriff eingeben und Enter drücken oder Suchen klicken.\n")
    result_area.insert(tk.END, "Klicken Sie auf „▶ In Obsidian öffnen“ um die Notiz zu öffnen.","")

    root.mainloop()


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Grover II – Inhaltssuche (Grover-inspiriert)")
    parser.add_argument("suchbegriff", nargs="?", help="Suchbegriff (ohne = GUI)")
    parser.add_argument("-v", "--vault", default=str(DEFAULT_VAULT), help="Vault-Pfad")
    parser.add_argument("-k", "--top", type=int, default=10, help="Anzahl Top-Treffer")
    parser.add_argument("-i", "--case-sensitive", action="store_true")
    parser.add_argument("-g", "--gui", action="store_true", help="GUI starten")
    args = parser.parse_args()

    if args.gui or (not args.suchbegriff and len(sys.argv) == 1):
        try:
            run_gui()
        except Exception as e:
            print(f"GUI fehlgeschlagen: {e}")
            print("→ Starte interaktive Terminal-Oberfläche.\n")
            run_interactive()
        return

    vault = Path(args.vault).expanduser()
    if not vault.exists():
        print(f"Vault nicht gefunden: {vault}")
        return

    files = find_md_files(vault)
    print(f"Grover II – Inhaltssuche in {len(files)} Dateien")
    print(f"Suchbegriff: {args.suchbegriff}\n")

    grover = GroverContentSearch(files, vault)
    results = grover.search(args.suchbegriff, case_sensitive=args.case_sensitive, top_k=args.top)

    for rel, full_path, score, matches in results:
        print(f"--- {rel} (Score: {score:.4f}) ---")
        for line_no, snippet in matches:
            print(f"  Z{line_no}: {snippet}")
        print()


if __name__ == "__main__":
    main()
