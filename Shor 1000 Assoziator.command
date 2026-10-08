#!/bin/bash
# Startet die Assoziator-Engine im Terminal. Doppelklick in Finder reicht.
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
SRC="$DIR/Shor 1000 Assoziator.cpp"
BIN="$DIR/shor_1000_assoziator"

if [[ ! -x "$BIN" || "$SRC" -nt "$BIN" ]]; then
  echo "Kompiliere (c++17 -O2) …"
  c++ -std=c++17 -O2 -Wall -Wextra -o "$BIN" "$SRC"
fi

exec "$BIN" "$@"
