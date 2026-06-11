#!/usr/bin/env bash
# Startet rundweg_schwung_resonanz.sage mit SageMath (nicht pyenv/python3).
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SAGE_BIN=""
if command -v sage >/dev/null 2>&1; then
  SAGE_BIN="$(command -v sage)"
elif [[ -x "${HOME}/.antigravity/antigravity/bin/sage" ]]; then
  SAGE_BIN="${HOME}/.antigravity/antigravity/bin/sage"
fi

if [[ -z "${SAGE_BIN}" ]]; then
  echo "Fehler: 'sage' nicht im PATH und nicht unter ~/.antigravity/antigravity/bin/sage." >&2
  echo "SageMath installieren, dann erneut versuchen. Siehe RUN.md in diesem Ordner." >&2
  exit 1
fi

exec "${SAGE_BIN}" "$DIR/rundweg_schwung_resonanz.sage" "$@"
