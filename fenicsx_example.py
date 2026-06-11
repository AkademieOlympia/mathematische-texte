#!/opt/miniconda3/envs/fenicsx-env/bin/python
"""
Einfaches FEniCSx Beispiel - Poisson-Gleichung
"""

from mpi4py import MPI
from dolfinx import mesh, fem, plot
import numpy as np

# Erstelle ein einfaches Gitter (Einheitsquadrat)
domain = mesh.create_unit_square(MPI.COMM_WORLD, 8, 8, mesh.CellType.triangle)

# Definiere den Funktionenraum (lineare Lagrange-Elemente)
from dolfinx.fem import functionspace
from dolfinx.fem.petsc import LinearProblem
from ufl import TrialFunction, TestFunction, dx, grad, inner

V = functionspace(domain, ("Lagrange", 1))

# Definiere die Randbedingungen
# Dirichlet-Randbedingung: u = 0 auf dem gesamten Rand
boundary_facets = mesh.locate_entities_boundary(
    domain, domain.topology.dim - 1, lambda x: np.full(x.shape[1], True)
)
boundary_dofs = fem.locate_dofs_topological(V, domain.topology.dim - 1, boundary_facets)
bc = fem.dirichletbc(0.0, boundary_dofs, V)

# Definiere die Bilinearform und Linearform
u = TrialFunction(V)
v = TestFunction(V)

# Poisson-Gleichung: -Δu = f mit f = 1
from dolfinx import default_scalar_type
f = fem.Constant(domain, default_scalar_type(1.0))
a = inner(grad(u), grad(v)) * dx
L = f * v * dx

# Löse das System
problem = LinearProblem(a, L, bcs=[bc])
uh = problem.solve()

# Ausgabe
print(f"Lösung berechnet! Maximaler Wert: {uh.x.array.max():.4f}")
print(f"Minimaler Wert: {uh.x.array.min():.4f}")

# Optional: Visualisierung mit pyvista
try:
    import pyvista
    pyvista.set_jupyter_backend("static")
    
    # Erstelle Plotter
    plotter = pyvista.Plotter()
    topology, cell_types, geometry = plot.vtk_mesh(domain, domain.topology.dim)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    
    # Füge Lösung hinzu
    grid.point_data["u"] = uh.x.array
    grid.set_active_scalars("u")
    
    plotter.add_mesh(grid, show_edges=True)
    plotter.show_axes()
    plotter.show()
    print("\nVisualisierung geöffnet!")
except Exception as e:
    print(f"Visualisierung nicht verfügbar: {e}")
