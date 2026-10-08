#!/usr/bin/env bash
# ==========================================================================
# D9/E8 MASTER ENGINE v28.0 - AUTOMATIC MAC DEPLOYMENT & SETUP SCRIPT
# ==========================================================================
set -e

echo "=========================================================================="
echo " 🚀 D9/E8 MASTER ENGINE v28.0: MAC ENVIRONMENT SETUP & DEPLOYMENT"
echo "=========================================================================="
echo ""

TARGET_DIR="${HOME}/Desktop/Mathematische Texte"
echo "--- [1/4] Prüfe Zielordner: $TARGET_DIR ---"
mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR"
echo "✓ Ordner ist bereit."
echo ""

echo "--- [2/4] Einrichten der Python-Umgebung (venv) ---"
if [ ! -d "venv" ]; then
    echo "Erstelle venv-Umgebung..."
    python3 -m venv venv
fi

source venv/bin/activate
echo "✓ Virtual Environment aktiviert: $(which python)"
echo ""

echo "--- [3/4] Installiere Abhängigkeiten (Streamlit, Numba, Seaborn etc.) ---"
pip install --quiet --upgrade pip
pip install --quiet streamlit numba numpy pandas matplotlib seaborn

echo "✓ Alle Pakete erfolgreich installiert!"
echo ""

echo "--- [4/4] Status & Startanweisung ---"
if [ -f "streamlit_d9_physics_web_ui-v10.py" ]; then
    echo "🚀 Starte Streamlit Dashboard v10.0..."
    streamlit run streamlit_d9_physics_web_ui-v10.py
else
    echo "⚠️ 'streamlit_d9_physics_web_ui-v10.py' wurde noch nicht in diesem Ordner abgelegt."
    echo "Bitte lade die Datei aus dem Studio-Panel herunter und speichere sie in:"
    echo "  $TARGET_DIR/"
    echo ""
    echo "Danach einfach erneut ausführen mit:"
    echo "  ./deploy_mac_workspace.sh"
fi
