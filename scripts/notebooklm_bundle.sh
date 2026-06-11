#!/usr/bin/env bash
# Erzeugt eine priorisierte NotebookLM-Quellensammlung (nicht alle 3700+ Dateien).
set -euo pipefail
cd "$(dirname "$0")/.."

OUT="_notebooklm/sources"
MAX_BYTES=$((50 * 1024 * 1024))
MAX_FILES=45

mkdir -p "$OUT"
rm -f "$OUT"/*

python3 - "$OUT" "$MAX_BYTES" "$MAX_FILES" <<'PY'
import os
import shutil
import subprocess
import sys
from pathlib import Path

out = Path(sys.argv[1])
max_bytes = int(sys.argv[2])
max_files = int(sys.argv[3])

skip_prefixes = (
    ".agents/", ".claude/", "node_modules/", "external/",
    "ptolemaeus-lean/.lake/", ".venv/", ".gitnexus/",
)
skip_names = {"AGENTS.md", "CLAUDE.md"}

def score(path: str) -> int:
    p = path.lower()
    s = 0
    if "eabc" in p:
        s += 120
    if path.count("/") == 0:
        s += 60
    if p.endswith(".md"):
        s += 40
    if p.endswith(".pdf"):
        s += 30
    if "readme" in p:
        s += 25
    return s

def safe_name(path: str) -> str:
    base = path.replace("/", "__")
    return base[:180] if len(base) > 180 else base

tracked = subprocess.check_output(["git", "ls-files"], text=True).splitlines()
candidates = []
for path in tracked:
    if any(path.startswith(pref) for pref in skip_prefixes):
        continue
    if Path(path).name in skip_names:
        continue
    # NotebookLM-Upload: PDF und Markdown zuverlässig; .tex liefert oft HTTP 400
    if not any(path.lower().endswith(ext) for ext in (".md", ".pdf")):
        continue
    if not os.path.isfile(path):
        continue
    size = os.path.getsize(path)
    if size > max_bytes:
        continue
    candidates.append((score(path), size, path))

candidates.sort(key=lambda t: (-t[0], t[1]))

manifest_lines = [
    "# Mathematische Texte — Inventar",
    "",
    f"Git-tracked gesamt: {len(tracked)} Dateien",
    f"In dieses Bundle aufgenommen: bis zu {max_files} priorisierte Quellen",
    "",
    "## Aufgenommene Quellen",
    "",
]

chosen = []
used_names = set()
for sc, size, path in candidates:
    if len(chosen) >= max_files:
        break
    dest_name = safe_name(path)
    if dest_name in used_names:
        dest_name = f"{len(chosen):03d}_{dest_name}"
    used_names.add(dest_name)
    shutil.copy2(path, out / dest_name)
    chosen.append((path, dest_name, size, sc))
    manifest_lines.append(f"- `{path}` → `{dest_name}` ({size // 1024} KiB, score {sc})")

manifest_lines += [
    "",
    "## Nicht aufgenommen (Auszug)",
    "",
    "Große PDFs, Skills, Build-Artefakte und Dateien jenseits des Limits.",
    f"Volles Repo: https://github.com/AkademieOlympia/mathematische-texte",
    f"Lean-Kern: https://github.com/AkademieOlympia/eabc-renorm",
]

overview = out / "00_overview.md"
overview.write_text(
    "# Mathematische Texte — NotebookLM-Bundle\n\n"
    "Priorisierte Auswahl (EABC-PDFs und Markdown; TeX nur im Git-Repo).\n"
    "Vollständiges Inventar: `MANIFEST.md`.\n",
    encoding="utf-8",
)

manifest = out / "MANIFEST.md"
manifest.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

print(f"==> {len(chosen)} Dateien nach {out}")
for path, dest, size, sc in chosen[:12]:
    print(f"  + {dest}  ({path})")
if len(chosen) > 12:
    print(f"  … und {len(chosen) - 12} weitere")
PY
