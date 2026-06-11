#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

if [[ -x "$HOME/Projects/eabc-renorm/.venv-notebooklm/bin/notebooklm" ]]; then
  echo "Nutze bestehende venv aus eabc-renorm"
  exit 0
fi

if [[ ! -d .venv-notebooklm ]]; then
  python3 -m venv .venv-notebooklm
fi
.venv-notebooklm/bin/pip install -q "notebooklm-py[browser]"
.venv-notebooklm/bin/playwright install chromium
echo "Fertig. Als Nächstes: .venv-notebooklm/bin/notebooklm login"
