#!/usr/bin/env python3
"""
Test-Skript für FEniCSx - läuft direkt mit der Conda-Umgebung
"""

import sys
print(f"Python-Pfad: {sys.executable}")
print(f"Python-Version: {sys.version}")

try:
    from mpi4py import MPI
    print(f"✓ mpi4py erfolgreich importiert")
    print(f"MPI Größe: {MPI.COMM_WORLD.size}")
except Exception as e:
    print(f"✗ Fehler beim Import von mpi4py: {e}")
    sys.exit(1)

try:
    from dolfinx import mesh, fem
    print(f"✓ dolfinx erfolgreich importiert")
    print(f"Version: {dolfinx.__version__}")
except Exception as e:
    print(f"✗ Fehler beim Import von dolfinx: {e}")
    sys.exit(1)

print("\n✓ Alle Imports erfolgreich! FEniCSx ist bereit.")
