#!/opt/miniconda3/bin/python3
"""
Analytische-Qubit-Suche als Mac-Anwendung.

Liest TXT-, PDF-, CSV- und TEX-Dateien aus den Benutzerordnern und von lokalen
Laufwerken einmal in eine Korpusdatei. Ist die Datei schon da, bleibt sie
unverändert. Nur „Korpus erneuern“ schreibt sie neu.
Die Suche läuft über diesen Korpus, mit Dialogeingabe der Suchbegriffe.
"""

import json
import math
import os
import queue
import subprocess
import sys
import threading
import time

# pyenv-Python 3.13 hat hier kein Tkinter. Miniconda hat Tkinter, NumPy und pypdf.
_TK_PYTHON = "/opt/miniconda3/bin/python3"
try:
    import tkinter as tk
    from tkinter import messagebox, simpledialog, ttk
    import numpy as np
except ModuleNotFoundError:
    if os.path.isfile(_TK_PYTHON) and os.path.realpath(sys.executable) != os.path.realpath(_TK_PYTHON):
        os.execv(_TK_PYTHON, [_TK_PYTHON, os.path.abspath(__file__), *sys.argv[1:]])
    sys.stderr.write(
        "Tkinter oder NumPy fehlt in diesem Python.\n"
        f"Start: {_TK_PYTHON} {os.path.abspath(__file__)}\n"
    )
    raise SystemExit(1)

try:
    from pypdf import PdfReader
except ImportError:
    PdfReader = None

CORPUS_PATH = os.path.expanduser(
    "~/Library/Application Support/AnalytischeQubitSuche/korpus.jsonl"
)
HEADER_LEN = 480
MAX_BYTES = 30 * 1024 * 1024
MAX_PDF_BYTES = 250 * 1024 * 1024
MAX_CHARS = 400_000
MAX_PDF_PAGES = 80
EXTENSIONS = {".txt", ".pdf", ".csv", ".tex"}
DEFAULT_QUERY = "analytische qubit, analytisches qubit, analytic qubit"
ROOT_LABEL = "Benutzer und lokale Laufwerke"


def as_utf8(value):
    """Entfernt einzelne UTF-16-Surrogate, die UTF-8 und Tk nicht schreiben können."""
    if value is None:
        return ""
    if not isinstance(value, str):
        value = str(value)
    return value.encode("utf-8", "replace").decode("utf-8")


SKIP_DIR_NAMES = {
    "Library", "System", "Applications", "node_modules", "venv",
    "site-packages", "__pycache__", "build", "Caches", "Cache",
    "DerivedData", "CoreSimulator", "Trash", "Pods",
    "texlive", "texmf-dist", "Cellar", "Homebrew", "miniconda3",
    "anaconda3", "micromamba", ".git",
}
SKIP_SUFFIXES = (
    ".app", ".photoslibrary", ".musiclibrary", ".tvlibrary",
    ".bundle", ".framework", ".plugin", ".kext", ".sparsebundle",
)
PRUNE_PREFIXES = (
    "/System", "/private", "/dev", "/cores", "/Library", "/Applications",
    "/usr", "/opt", "/bin", "/sbin", "/Network", "/net", "/var",
)
_NETWORK_TYPES = {"smbfs", "nfs", "afpfs", "webdav", "cifs", "ftpfs", "osxfuse", "macfuse", "fuse"}
_NETWORK_MOUNTS = None


def seiberg_witten_frequency(u=2.5, lambda_scale=1.0):
    """AGM-Periode wie im Processing-Qubit der D9-Engine."""
    k2 = 2.0 * (lambda_scale ** 2) / (u + (lambda_scale ** 2))
    k_val = math.sqrt(max(0.0, min(1.0, k2)))
    a, b = 1.0, math.sqrt(1.0 - k_val * k_val)
    for _ in range(10):
        a, b = 0.5 * (a + b), math.sqrt(a * b)
    elliptic_k = math.pi / (2.0 * a)
    return math.sqrt(2.0) / (math.pi * math.sqrt(u + (lambda_scale ** 2))) * elliptic_k


def lock_analytic_qubit(score):
    """Score auf E8-Wurzel, V4-Sektor und Krylov-Tiefe legen."""
    weight = int(score)
    return weight % 240, weight % 4, weight % 64


def analytic_qubit_search(candidates, keywords_lower):
    """Grover-Schritt über die Treffer des Korpus. Gibt Treffer und Kennzahlen zurück."""
    n = len(candidates)
    if n == 0:
        return [], {"N": 0, "M": 0, "sw": 0.0, "k_classic": 0.0, "k": 0}

    scores = np.zeros(n, dtype=np.float64)
    for i, (path, text) in enumerate(candidates):
        score = sum(text.count(kw) for kw in keywords_lower)
        fname = os.path.basename(path).lower()
        if any(kw in fname for kw in keywords_lower):
            score += 5
        scores[i] = score

    marked = np.where(scores > 0)[0]
    m = len(marked)
    if m == 0:
        return [], {"N": n, "M": 0, "sw": 0.0, "k_classic": 0.0, "k": 0}

    sw_freq = seiberg_witten_frequency()
    k_classic = (math.pi / 4.0) * math.sqrt(n / m)
    k_qubit = max(1, int(math.floor(k_classic / (1.0 + 0.2 * sw_freq))))

    state = np.full(n, 1.0 / math.sqrt(n), dtype=np.float64)
    for _ in range(k_qubit):
        state[marked] *= -1.0
        state = 2.0 * np.mean(state) - state

    probabilities = state ** 2
    order = np.lexsort((-scores, -probabilities))
    results = []
    for idx in order:
        if scores[idx] <= 0:
            continue
        root, sector, depth = lock_analytic_qubit(scores[idx])
        results.append((candidates[idx][0], float(scores[idx]), float(probabilities[idx]), root, sector, depth))
    info = {"N": n, "M": m, "sw": sw_freq, "k_classic": k_classic, "k": k_qubit}
    return results, info


def format_header(files, created, root):
    payload = json.dumps(
        {"kind": "meta", "version": 2, "files": files, "created": created, "root": root},
        ensure_ascii=False,
    )
    if len(payload) > HEADER_LEN:
        payload = json.dumps({"kind": "meta", "version": 2, "files": files, "created": created, "root": "/"})
    return payload.ljust(HEADER_LEN) + "\n"


def read_meta(path):
    try:
        with open(path, "r", encoding="utf-8") as handle:
            obj = json.loads(handle.readline())
    except (OSError, json.JSONDecodeError, UnicodeError):
        return None
    if isinstance(obj, dict) and obj.get("kind") == "meta":
        return obj
    return None


def split_terms(raw):
    raw = (raw or "").strip()
    if not raw:
        return []
    if "," in raw:
        return [part.strip().lower() for part in raw.split(",") if part.strip()]
    return [raw.lower()]


def _network_mount_points():
    global _NETWORK_MOUNTS
    if _NETWORK_MOUNTS is not None:
        return _NETWORK_MOUNTS
    found = set()
    try:
        out = subprocess.check_output(["mount"], text=True, errors="replace")
    except (OSError, subprocess.CalledProcessError):
        _NETWORK_MOUNTS = found
        return found
    for line in out.splitlines():
        if " on " not in line or "(" not in line:
            continue
        path = line.split(" on ", 1)[1].rsplit(" (", 1)[0]
        fstype = line.rsplit("(", 1)[-1].split(",", 1)[0].strip().lower()
        if fstype in _NETWORK_TYPES or fstype.startswith("nfs"):
            found.add(os.path.realpath(path))
    _NETWORK_MOUNTS = found
    return found


def _on_network(path):
    real = os.path.realpath(path)
    for mount in _network_mount_points():
        if real == mount or real.startswith(mount + os.sep):
            return True
    return False


def _pruned(path):
    for prefix in PRUNE_PREFIXES:
        if path == prefix or path.startswith(prefix + os.sep):
            return True
    return _on_network(path)


def _skip_dir(name, path):
    if name.startswith(".") or name in SKIP_DIR_NAMES:
        return True
    lower = name.lower()
    if lower.endswith(SKIP_SUFFIXES):
        return True
    return _pruned(path)


def scan_roots():
    """Benutzerordner zuerst, danach lokale Laufwerke. Kein Netz, kein System."""
    roots = []
    if os.path.isdir("/Users"):
        roots.append("/Users")
    volumes = "/Volumes"
    if os.path.isdir(volumes):
        try:
            names = sorted(os.listdir(volumes))
        except OSError:
            names = []
        for name in names:
            if name.startswith(".") or name.startswith("com.apple"):
                continue
            path = os.path.join(volumes, name)
            if not os.path.isdir(path) or _skip_dir(name, path):
                continue
            try:
                if os.path.samefile(path, "/"):
                    continue
            except OSError:
                pass
            roots.append(path)
    return roots or ["/Users"]


def iter_source_files(roots, cancel_event=None):
    """TXT, PDF, CSV und TEX, ohne System-, Cache- und Programmordner."""
    if isinstance(roots, (str, os.PathLike)):
        roots = [os.fspath(roots)]
    for root in roots:
        if cancel_event is not None and cancel_event.is_set():
            return
        if not os.path.isdir(root) or _pruned(root):
            continue
        for current, dirs, files in os.walk(root, onerror=lambda _err: None, followlinks=False):
            if cancel_event is not None and cancel_event.is_set():
                return
            if _pruned(current):
                dirs.clear()
                continue
            dirs[:] = [name for name in dirs if not _skip_dir(name, os.path.join(current, name))]
            for name in files:
                if name.startswith("."):
                    continue
                if os.path.splitext(name)[1].lower() not in EXTENSIONS:
                    continue
                yield os.path.join(current, name)


def _read_plain(path):
    try:
        with open(path, "rb") as handle:
            raw = handle.read(MAX_BYTES)
    except OSError:
        return ""
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("cp1252", "replace")
    return text[:MAX_CHARS]


def extract_text(path):
    try:
        size = os.path.getsize(path)
    except OSError:
        return ""
    if size <= 0:
        return ""
    ext = os.path.splitext(path)[1].lower()
    try:
        if ext == ".pdf":
            if size > MAX_PDF_BYTES:
                return ""
            text = _read_pdf(path)
        else:
            if size > MAX_BYTES:
                return ""
            text = _read_plain(path)
    except Exception:
        return ""
    return as_utf8(text).lower()[:MAX_CHARS]


def _read_pdf(path):
    if PdfReader is None:
        return ""
    try:
        reader = PdfReader(path)
        parts = []
        total = 0
        for index, page in enumerate(reader.pages):
            if index >= MAX_PDF_PAGES or total >= MAX_CHARS:
                break
            piece = page.extract_text() or ""
            parts.append(piece)
            total += len(piece)
        return "".join(parts)
    except Exception:
        return ""


def _write_json_line(handle, obj):
    handle.write(as_utf8(json.dumps(obj, ensure_ascii=False)) + "\n")


def build_corpus(dest, roots=None, cancel_event=None, progress=None):
    """Schreibt den Korpus nach dest. Die bestehende Datei wird erst am Ende ersetzt."""
    if roots is None:
        roots = scan_roots()
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".tmp"
    created = time.strftime("%Y-%m-%d %H:%M:%S")
    seen = 0
    indexed = 0
    cancelled = False
    with open(tmp, "w", encoding="utf-8", errors="replace") as out:
        out.write(format_header(0, created, ROOT_LABEL))
        for path in iter_source_files(roots, cancel_event):
            seen += 1
            if cancel_event is not None and cancel_event.is_set():
                cancelled = True
                break
            try:
                text = extract_text(path)
                if text.strip():
                    _write_json_line(out, {"path": as_utf8(path), "text": text})
                    indexed += 1
                    if indexed % 25 == 0:
                        out.flush()
            except Exception:
                continue
            if progress and seen % 20 == 0:
                progress(seen, indexed, as_utf8(path))
        if not cancelled:
            out.seek(0)
            out.write(format_header(indexed, created, ROOT_LABEL))
            out.flush()
    if cancelled:
        try:
            os.remove(tmp)
        except OSError:
            pass
        return None
    os.replace(tmp, dest)
    if progress:
        progress(seen, indexed, "")
    return {"files": indexed, "seen": seen, "created": created, "root": ROOT_LABEL}


def search_corpus(corpus_path, terms, cancel_event=None, progress=None):
    keywords = [term.lower() for term in terms if term.strip()]
    candidates = []
    checked = 0
    with open(corpus_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if cancel_event is not None and cancel_event.is_set():
                return [], {"N": 0, "M": 0, "sw": 0.0, "k_classic": 0.0, "k": 0, "cancelled": True}
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except (json.JSONDecodeError, UnicodeError):
                continue
            if not isinstance(record, dict) or record.get("kind") == "meta":
                continue
            text = as_utf8(record.get("text") or "")
            path = as_utf8(record.get("path") or "")
            checked += 1
            name = os.path.basename(path).lower()
            if any(kw in name or kw in text for kw in keywords):
                candidates.append((path, text))
            if progress and checked % 200 == 0:
                progress(checked, len(candidates))
    results, info = analytic_qubit_search(candidates, keywords)
    info["checked"] = checked
    return results, info


class GroverApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Analytische Qubits — Grover-Suche")
        self.geometry("1080x680")
        self.minsize(780, 460)
        self.events = queue.Queue()
        self.cancel_event = threading.Event()
        self.busy = False
        self._build_widgets()
        self.after(100, self._poll)
        self.after(200, self._startup)
        self.after(250, self._bring_front)
        if sys.platform == "darwin":
            try:
                self.createcommand("::tk::mac::ReopenApplication", self._bring_front)
            except tk.TclError:
                pass

    def _bring_front(self):
        self.deiconify()
        self.lift()
        try:
            self.attributes("-topmost", True)
            self.after(300, lambda: self.attributes("-topmost", False))
            self.focus_force()
        except tk.TclError:
            pass

    def _build_widgets(self):
        bar = ttk.Frame(self, padding=12)
        bar.pack(fill="x")
        self.status = tk.StringVar(value="Starte …")
        ttk.Label(bar, textvariable=self.status, wraplength=860).pack(fill="x")
        self.detail = tk.StringVar(value="")
        ttk.Label(bar, textvariable=self.detail, wraplength=860).pack(fill="x", pady=(4, 0))

        buttons = ttk.Frame(self, padding=(12, 0, 12, 8))
        buttons.pack(fill="x")
        self.search_button = ttk.Button(buttons, text="Suchbegriffe …", command=self.on_search)
        self.search_button.pack(side="left")
        self.refresh_button = ttk.Button(buttons, text="Korpus erneuern", command=self.on_refresh)
        self.refresh_button.pack(side="left", padx=(8, 0))
        self.stop_button = ttk.Button(buttons, text="Abbrechen", command=self.on_cancel, state="disabled")
        self.stop_button.pack(side="left", padx=(8, 0))
        ttk.Label(buttons, text=CORPUS_PATH).pack(side="left", padx=(16, 0))

        table_frame = ttk.Frame(self, padding=(12, 0, 12, 12))
        table_frame.pack(fill="both", expand=True)
        columns = ("rank", "score", "prob", "root", "sector", "depth", "path")
        self.table = ttk.Treeview(table_frame, columns=columns, show="headings")
        headings = {
            "rank": "Rang", "score": "Score", "prob": "P %",
            "root": "Wurzel", "sector": "Sektor", "depth": "Tiefe", "path": "Datei",
        }
        widths = {"rank": 60, "score": 70, "prob": 70, "root": 70, "sector": 70, "depth": 70, "path": 640}
        for key in columns:
            self.table.heading(key, text=headings[key])
            self.table.column(key, width=widths[key], anchor="w")
        scroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.table.yview)
        self.table.configure(yscrollcommand=scroll.set)
        self.table.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")
        self.table.bind("<Double-1>", self._open_selected)

    def _startup(self):
        if os.path.isfile(CORPUS_PATH) and os.path.getsize(CORPUS_PATH) > 0:
            meta = read_meta(CORPUS_PATH) or {}
            count = meta.get("files")
            created = meta.get("created") or "unbekannt"
            count_text = f"{count} Dateien" if isinstance(count, int) else "Dateien"
            self.status.set(f"Korpus vorhanden ({count_text}, {created}). Er wird nicht verändert.")
            self.detail.set("Suche liest nur diese Datei. „Korpus erneuern“ schreibt sie neu.")
            return
        roots = scan_roots()
        self.status.set("Kein Korpus. Benutzerordner und lokale Laufwerke werden eingelesen.")
        self.detail.set(
            "Quellen: " + ", ".join(roots) + ". System, Programme, TeX Live und Caches bleiben ausgenommen."
        )
        self._start_build()

    def _set_busy(self, busy):
        self.busy = busy
        state = "disabled" if busy else "normal"
        self.search_button.configure(state=state)
        self.refresh_button.configure(state=state)
        self.stop_button.configure(state="normal" if busy else "disabled")

    def _start_build(self):
        self.cancel_event.clear()
        self._set_busy(True)
        self.detail.set("Einlesen läuft …")

        def progress(seen, indexed, path):
            self.events.put(("progress", f"Gesehen {seen}, im Korpus {indexed}. {path}"))

        def work():
            try:
                meta = build_corpus(CORPUS_PATH, cancel_event=self.cancel_event, progress=progress)
            except Exception as exc:
                self.events.put(("error", as_utf8(exc)))
                return
            self.events.put(("built", meta))

        threading.Thread(target=work, daemon=True).start()

    def on_refresh(self):
        if self.busy:
            return
        if os.path.isfile(CORPUS_PATH):
            replace = messagebox.askyesno(
                "Korpus erneuern",
                "Die vorhandene Korpusdatei wird erst ersetzt, wenn das Einlesen durch ist.",
                parent=self,
            )
            if not replace:
                return
        self.status.set("Korpus wird neu eingelesen …")
        self._start_build()

    def on_cancel(self):
        self.cancel_event.set()
        self.detail.set("Abbruch angefordert. Die bisherige Korpusdatei bleibt stehen.")

    def on_search(self):
        if self.busy:
            return
        if not os.path.isfile(CORPUS_PATH):
            messagebox.showinfo("Kein Korpus", "Zuerst muss der Korpus eingelesen sein.", parent=self)
            return
        raw = simpledialog.askstring(
            "Suchbegriffe",
            "Begriffe durch Komma trennen.\nOhne Komma gilt die ganze Eingabe als eine Phrase.",
            initialvalue=DEFAULT_QUERY,
            parent=self,
        )
        if raw is None:
            return
        terms = split_terms(raw)
        if not terms:
            messagebox.showinfo("Suchbegriffe", "Bitte mindestens einen Begriff eingeben.", parent=self)
            return
        self._start_search(terms)

    def _start_search(self, terms):
        self.cancel_event.clear()
        self._set_busy(True)
        self.status.set("Suche im Korpus: " + ", ".join(terms))
        for item in self.table.get_children():
            self.table.delete(item)

        def progress(checked, hits):
            self.events.put(("progress", f"Korpuszeilen {checked}, Treffer bisher {hits}."))

        def work():
            try:
                results, info = search_corpus(
                    CORPUS_PATH, terms, cancel_event=self.cancel_event, progress=progress,
                )
            except Exception as exc:
                self.events.put(("error", as_utf8(exc)))
                return
            self.events.put(("results", (results, info, terms)))

        threading.Thread(target=work, daemon=True).start()

    def _poll(self):
        try:
            while True:
                kind, payload = self.events.get_nowait()
                if kind == "progress":
                    self.detail.set(as_utf8(payload))
                elif kind == "error":
                    self._set_busy(False)
                    self.status.set("Fehler.")
                    messagebox.showerror("Fehler", as_utf8(payload), parent=self)
                elif kind == "built":
                    self._set_busy(False)
                    if payload is None:
                        self.status.set("Einlesen abgebrochen. Vorhandener Korpus, falls es ihn gab, ist unverändert.")
                    else:
                        self.status.set(
                            f"Korpus geschrieben: {payload['files']} Dateien, {payload['created']}."
                        )
                        self.detail.set(CORPUS_PATH)
                elif kind == "results":
                    self._show_results(*payload)
        except queue.Empty:
            pass
        self.after(100, self._poll)

    def _show_results(self, results, info, terms):
        self._set_busy(False)
        if info.get("cancelled"):
            self.status.set("Suche abgebrochen.")
            return
        shown = results[:200]
        for rank, (path, score, prob, root, sector, depth) in enumerate(shown, 1):
            self.table.insert(
                "", "end",
                values=(rank, f"{score:.0f}", f"{prob * 100:.2f}", root, sector, depth, path),
            )
        self.status.set(
            f"{info['M']} Treffer in {info.get('checked', 0)} Korpusdateien. "
            f"Angezeigt: {len(shown)}. ω = {info['sw']:.4f}, Iterationen {info['k']}."
        )
        self.detail.set("Begriffe: " + ", ".join(terms) + ". Doppelklick öffnet die Datei.")

    def _open_selected(self, _event):
        selected = self.table.selection()
        if not selected:
            return
        path = self.table.item(selected[0], "values")[-1]
        if path and os.path.exists(path):
            subprocess.run(["open", path], check=False)


def main():
    log_dir = os.path.dirname(CORPUS_PATH)
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, "app.log")
    log_handle = open(log_path, "a", encoding="utf-8", errors="replace")
    sys.stdout = log_handle
    sys.stderr = log_handle

    def _excepthook(exc_type, exc, tb):
        import traceback
        traceback.print_exception(exc_type, exc, tb, file=log_handle)
        log_handle.flush()

    sys.excepthook = _excepthook
    log_handle.write("\n--- start %s ---\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
    log_handle.flush()
    app = GroverApp()
    app.mainloop()


if __name__ == "__main__":
    main()
