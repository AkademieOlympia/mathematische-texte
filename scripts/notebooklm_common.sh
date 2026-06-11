#!/usr/bin/env bash
# Gemeinsame Hilfen für NotebookLM-Skripte (Mathematische Texte).
set -euo pipefail

NOTEBOOK_TITLE="Mathematische Texte"
REPO_URL="https://github.com/AkademieOlympia/mathematische-texte"
RELATED_URL="https://github.com/AkademieOlympia/eabc-renorm"

resolve_nlm() {
  local here
  here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
  NLM="$here/.venv-notebooklm/bin/notebooklm"
  if [[ ! -x "$NLM" ]] && [[ -x "$HOME/Projects/eabc-renorm/.venv-notebooklm/bin/notebooklm" ]]; then
    NLM="$HOME/Projects/eabc-renorm/.venv-notebooklm/bin/notebooklm"
  fi
  export NLM
}

require_nlm() {
  resolve_nlm
  if [[ ! -x "$NLM" ]]; then
    echo "NotebookLM-CLI fehlt. Einrichten mit:" >&2
    echo "  bash scripts/setup_notebooklm.sh" >&2
    echo "  # oder die venv aus ~/Projects/eabc-renorm nutzen" >&2
    exit 1
  fi
  if ! "$NLM" auth check --test --json 2>/dev/null | grep -q '"status": "ok"'; then
    echo "NotebookLM nicht angemeldet. Bitte einmalig:" >&2
    echo "  $NLM login" >&2
    exit 2
  fi
}

resolve_notebook_id() {
  local title=$1
  "$NLM" list --json | python3 -c "
import json, sys
title = sys.argv[1]
for nb in json.load(sys.stdin).get('notebooks', []):
    if nb.get('title') == title:
        print(nb['id'])
        break
" "$title"
}
