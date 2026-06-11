#!/bin/bash
# Hofbi-Autonomous-Batch

TARGETS=("861163313" "1234567891" "1073741823")
ZETA_TERMS=50000
SONDEN=14

for N in "${TARGETS[@]}"; do
    echo "-----------------------------------------------"
    echo "[BATCH] Starte Ionisation fuer N=$N"
    
    # Der Scanner kalibriert das Sonden-Gitter jetzt selbst um 1/sqrt(N).
    ./rydberg_scanner_m5 "$N" "auto" "$ZETA_TERMS" "$SONDEN"
    
    # Hier koennte eine Logik folgen, die trajectories.csv auswertet
    # und bei Erfolg zum naechsten Ziel springt.
done