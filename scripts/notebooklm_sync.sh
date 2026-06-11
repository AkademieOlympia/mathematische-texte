#!/usr/bin/env bash
# Verbindet Mathematische Texte mit einem NotebookLM-Notebook.
set -euo pipefail
cd "$(dirname "$0")/.."
# shellcheck source=notebooklm_common.sh
source "$(dirname "$0")/notebooklm_common.sh"

require_nlm
bash "$(dirname "$0")/notebooklm_bundle.sh"

NOTEBOOK_ID=$(resolve_notebook_id "$NOTEBOOK_TITLE")

if [[ -z "${NOTEBOOK_ID:-}" ]]; then
  echo "==> Erstelle Notebook: $NOTEBOOK_TITLE"
  NOTEBOOK_ID=$("$NLM" create "$NOTEBOOK_TITLE" --use --json | python3 -c "
import json, sys
data = json.load(sys.stdin)
print(data.get('notebook', data)['id'])
")
else
  echo "==> Nutze bestehendes Notebook: $NOTEBOOK_TITLE ($NOTEBOOK_ID)"
  "$NLM" use "$NOTEBOOK_ID"
fi

echo "==> Lade Quellen hoch"
"$NLM" source add "$REPO_URL" || true
"$NLM" source add "$RELATED_URL" || true
for f in _notebooklm/sources/*; do
  [[ -f "$f" ]] || continue
  echo "  + $(basename "$f")"
  "$NLM" source add "$f" || echo "    (übersprungen)" >&2
done

echo
echo "==> NotebookLM verbunden"
echo "    Titel:   $NOTEBOOK_TITLE"
echo "    ID:      $NOTEBOOK_ID"
echo "    Öffnen:  https://notebooklm.google.com/"
