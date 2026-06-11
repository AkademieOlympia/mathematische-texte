#!/bin/bash
# Wrapper-Skript zum Ausführen von FEniCSx-Skripten mit der richtigen Conda-Umgebung

# Verwende den vollständigen Python-Pfad der Conda-Umgebung
PYTHON="/opt/miniconda3/envs/fenicsx-env/bin/python"

# Führe das Skript aus
"$PYTHON" "$@"
