#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

if ! command -v lake >/dev/null 2>&1; then
  if [[ -x "${HOME}/.elan/bin/lake" ]]; then
    export PATH="${HOME}/.elan/bin:${PATH}"
  else
    echo "lake nicht gefunden. Elan liegt nicht unter ${HOME}/.elan/bin." >&2
    exit 127
  fi
fi

echo "Projekt: $(pwd)"
echo "Toolchain: $(tr -d '[:space:]' < lean-toolchain)"
echo "Baue lokale Beweise: EabcToposProofs"
lake build EabcToposProofs
echo "Lean-Validierung erfolgreich."
