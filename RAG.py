from pathlib import Path
from .config import RAW_DATA_DIR
from .rag import LocalRAG

def ingest_all_sources():
    """Liest alle validen Rohdaten aus Physik-, Mathematik- und Wikipedia-Ordnern ein."""
    rag = LocalRAG()
    
    subdirs = {
        "physics": RAW_DATA_DIR / "physics",
        "mathematics": RAW_DATA_DIR / "mathematics",
        "wikipedia": RAW_DATA_DIR / "wikipedia"
    }
    
    # Erlaubte Text-Dateiendungen
    VALID_EXTENSIONS = {".txt", ".md", ".tex", ".json", ".sage", ".py"}
    
    for source_type, path in subdirs.items():
        if not path.exists():
            continue
        for file_path in path.glob("**/*.*"):
            # Ignoriere macOS Meta-Dateien (.DS_Store) und nicht-Textformate
            if file_path.is_file() and not file_path.name.startswith(".") and file_path.suffix.lower() in VALID_EXTENSIONS:
                print(f"Ingestiere [{source_type}]: {file_path.name}")
                rag.ingest_file(str(file_path), source_type=source_type)

if __name__ == "__main__":
    ingest_all_sources()