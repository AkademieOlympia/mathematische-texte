#!/opt/miniconda3/envs/fenicsx-env/bin/python
"""
FEniCSx Interaktives Skript - für interaktive Nutzung
"""

from mpi4py import MPI
from dolfinx import mesh, fem, default_scalar_type, plot
from dolfinx.fem import functionspace
from dolfinx.fem.petsc import LinearProblem
from ufl import TrialFunction, TestFunction, dx, grad, inner
import numpy as np

print("FEniCSx interaktive Umgebung geladen!")
print(f"MPI Größe: {MPI.COMM_WORLD.size}")
print("\nVerfügbare Variablen:")
print("  - MPI: MPI-Interface")
print("  - mesh: Mesh-Funktionen")
print("  - fem: Finite Element Method Funktionen")
print("  - np: NumPy")
print("\nBeispiel:")
print("  domain = mesh.create_unit_square(MPI.COMM_WORLD, 8, 8)")
print("  V = functionspace(domain, ('Lagrange', 1))")
print("\nStarten Sie mit: python -i fenicsx_interaktiv.py")
