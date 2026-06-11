import hashlib
import json
from datetime import datetime

def generate_research_fingerprint(files):
    """
    Erstellt einen kryptographischen Fingerabdruck für die #Energiedoku.
    Nutzt SHA-256, um die Unveränderlichkeit der Forschungsdaten zu garantieren.
    """
    combined_hash = hashlib.sha256()
    file_manifest = {}

    for file_path in files:
        try:
            with open(file_path, "rb") as f:
                file_content = f.read()
                file_hash = hashlib.sha256(file_content).hexdigest()
                file_manifest[file_path] = file_hash
                combined_hash.update(file_content)
        except FileNotFoundError:
            print(f"Warnung: Datei {file_path} nicht gefunden.")

    root_hash = combined_hash.hexdigest()
    
    # Metadaten für die zensurresistente Dokumentation
    metadata = {
        "title": "#Energiedoku: Quaternion Prime Model",
        "author": "Bamberg Research",
        "timestamp": datetime.now().isoformat(),
        "root_hash": root_hash,
        "manifest": file_manifest,
        "note": "Kryptographischer Beweis zur Riemannschen Vermutung."
    }
    
    return metadata

# Beispielhafte Anwendung für deine Dateien in Bamberg:
# Liste hier deine TeX-Files, Sage-Skripte oder Daten-Tabellen auf
research_files = ["riemann_quaternion_model.sage", "thesis_draft.pdf"]
fingerprint = generate_research_fingerprint(research_files)

print(json.dumps(fingerprint, indent=4))

# Speichern als 'proof_of_existence.json' für den IPFS-Upload
with open("proof_of_existence.json", "w") as f:
    json.dump(fingerprint, f, indent=4)