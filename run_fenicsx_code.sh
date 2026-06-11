#!/bin/bash
# Führt FEniCSx-Code aus und schreibt Ergebnis in JSON-Datei

HELPER_SCRIPT="/Users/thomashoffbauer/Desktop/Mathematische Texte/fenicsx_helper.py"
CODE_FILE="$1"
RESULT_FILE="$2"

# Führe Helper-Skript aus
cat "$CODE_FILE" | /opt/miniconda3/envs/fenicsx-env/bin/python "$HELPER_SCRIPT" > "$RESULT_FILE" 2>&1
