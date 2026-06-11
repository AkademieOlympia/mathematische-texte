# Anleitung: FEniCSx Kernel in Jupyter verwenden

## Problem
Das Notebook verwendet den SageMath-Kernel statt des FEniCSx-Kernels.

## Lösung: Kernel im Jupyter-Interface wechseln

### Methode 1: Über das Kernel-Menü (empfohlen)

1. **Öffnen Sie das Notebook in Jupyter Lab**
2. **Klicken Sie rechts oben auf den Kernel-Namen** (zeigt aktuell "SageMath 10.8" oder ähnlich)
3. **Wählen Sie "Python 3.12 (fenicsx-env)" aus der Liste**
4. **Warten Sie, bis der Kernel neu gestartet wurde** (Sie sehen "Kernel starting..." oder ähnlich)
5. **Führen Sie die Zellen erneut aus**

### Methode 2: Über das Menü

1. **Klicken Sie auf "Kernel" im Menü**
2. **Wählen Sie "Change Kernel"**
3. **Wählen Sie "Python 3.12 (fenicsx-env)"**
4. **Warten Sie auf den Neustart**
5. **Führen Sie die Zellen erneut aus**

### Methode 3: Kernel neu starten

1. **Kernel → Restart Kernel**
2. **Dann: Kernel → Change Kernel → "Python 3.12 (fenicsx-env)"**
3. **Zellen erneut ausführen**

## Überprüfung

Nach dem Kernel-Wechsel sollte die erste Zelle erfolgreich ausführen:

```python
from mpi4py import MPI
from dolfinx import mesh, fem
import numpy as np

print("FEniCSx erfolgreich importiert!")
print(f"MPI Größe: {MPI.COMM_WORLD.size}")
```

**Erwartete Ausgabe:**
```
FEniCSx erfolgreich importiert!
MPI Größe: 1
```

**Falls Sie immer noch Fehler sehen:**
- Stellen Sie sicher, dass der Kernel wirklich "Python 3.12 (fenicsx-env)" ist
- Starten Sie Jupyter Lab neu
- Prüfen Sie, ob der Kernel in der Liste erscheint: `jupyter kernelspec list`

## Alternative: Direktes Python-Skript

Falls Jupyter weiterhin Probleme macht, können Sie auch direkt Python-Skripte verwenden:

```bash
conda activate fenicsx-env
python fenicsx_example.py
```
