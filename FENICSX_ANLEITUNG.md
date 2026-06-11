# FEniCSx Nutzung - Anleitung

## Problem
Der SageMath-Kernel startet nicht (Singular fehlt), daher können wir FEniCSx nicht direkt im Jupyter-Notebook mit SageMath verwenden.

## Lösung: Standalone-Skripte

### Option 1: Standalone-Skript ausführen (empfohlen)

```bash
# Direkt ausführen:
./fenicsx_standalone.py

# Oder mit dem Wrapper:
./run_fenicsx.sh fenicsx_standalone.py

# Oder mit vollständigem Pfad:
/opt/miniconda3/envs/fenicsx-env/bin/python fenicsx_standalone.py
```

Das Skript führt alle Beispiele aus:
- Gitter erstellen
- Funktionenraum definieren
- Poisson-Gleichung lösen
- Visualisierung (falls pyvista verfügbar)

### Option 2: Interaktive Python-Shell

```bash
# Starten Sie eine interaktive Python-Shell mit FEniCSx:
python -i fenicsx_interaktiv.py

# Dann können Sie interaktiv arbeiten:
>>> domain = mesh.create_unit_square(MPI.COMM_WORLD, 8, 8)
>>> V = functionspace(domain, ("Lagrange", 1))
>>> # ... weitere Berechnungen
```

### Option 3: Eigene Skripte erstellen

Erstellen Sie Ihre eigenen Python-Skripte basierend auf `fenicsx_example.py`:

```bash
# Kopieren Sie das Beispiel:
cp fenicsx_example.py mein_problem.py

# Bearbeiten Sie es:
# nano mein_problem.py  # oder Ihr bevorzugter Editor

# Führen Sie es aus:
./run_fenicsx.sh mein_problem.py
```

## Verfügbare Dateien

- `fenicsx_standalone.py` - Vollständiges Beispiel-Skript
- `fenicsx_example.py` - Einfaches Beispiel
- `fenicsx_interaktiv.py` - Für interaktive Nutzung
- `run_fenicsx.sh` - Wrapper-Skript zum Ausführen

## Wichtige Hinweise

1. **Python-Version**: FEniCSx benötigt Python 3.12 (nicht 3.13 wie SageMath)
2. **MPI-Fehler**: Der MPI-Finalisierungsfehler am Ende ist normal und kann ignoriert werden
3. **Visualisierung**: pyvista sollte automatisch funktionieren, falls installiert

## Beispiel: Eigenes Problem lösen

```python
#!/opt/miniconda3/envs/fenicsx-env/bin/python
from mpi4py import MPI
from dolfinx import mesh, fem, default_scalar_type
from dolfinx.fem import functionspace
from dolfinx.fem.petsc import LinearProblem
from ufl import TrialFunction, TestFunction, dx, grad, inner
import numpy as np

# Ihr Code hier...
domain = mesh.create_unit_square(MPI.COMM_WORLD, 20, 20)
# ... weitere Berechnungen
```

## Hilfe

- FEniCSx Dokumentation: https://fenicsproject.org/
- Beispiele: https://github.com/FEniCS/dolfinx/tree/main/python/demo
