#!/opt/miniconda3/envs/fenicsx-env/bin/python
"""
FEniCSx Standalone-Skript - funktioniert ohne Jupyter
Führt alle Beispiele aus dem Notebook aus
"""

from mpi4py import MPI
from dolfinx import mesh, fem, default_scalar_type, plot
from dolfinx.fem import functionspace
from dolfinx.fem.petsc import LinearProblem
from ufl import TrialFunction, TestFunction, dx, grad, inner
import numpy as np

print("=" * 60)
print("FEniCSx Standalone Beispiel")
print("=" * 60)

# 1. Import testen
print("\n1. Import-Test:")
print(f"   Python: {MPI.COMM_WORLD.size} MPI-Prozess(e)")
print("   ✓ Alle Module erfolgreich importiert")

# 2. Gitter erstellen
print("\n2. Gitter erstellen:")
domain = mesh.create_unit_square(MPI.COMM_WORLD, 10, 10, mesh.CellType.triangle)
print(f"   ✓ Gitter erstellt: {domain.geometry.x.shape[0]} Knoten")

# 3. Funktionenraum
print("\n3. Funktionenraum definieren:")
V = functionspace(domain, ("Lagrange", 1))
print(f"   ✓ Funktionenraum: {V.dofmap.index_map.size_global} Freiheitsgrade")

# 4. Poisson-Gleichung lösen
print("\n4. Poisson-Gleichung lösen:")
# Randbedingungen
boundary_facets = mesh.locate_entities_boundary(
    domain, domain.topology.dim - 1, lambda x: np.full(x.shape[1], True)
)
boundary_dofs = fem.locate_dofs_topological(V, domain.topology.dim - 1, boundary_facets)
bc = fem.dirichletbc(default_scalar_type(0.0), boundary_dofs, V)

# Bilinearform und Linearform
u = TrialFunction(V)
v = TestFunction(V)
f = fem.Constant(domain, default_scalar_type(1.0))
a = inner(grad(u), grad(v)) * dx
L = f * v * dx

# Löse
problem = LinearProblem(a, L, bcs=[bc])
uh = problem.solve()

print(f"   ✓ Lösung berechnet!")
print(f"   Max: {uh.x.array.max():.4f}, Min: {uh.x.array.min():.4f}")

# 5. Optional: Visualisierung
print("\n5. Visualisierung:")
try:
    import pyvista
    pyvista.set_jupyter_backend("static")
    
    topology, cell_types, geometry = plot.vtk_mesh(domain, domain.topology.dim)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    grid.point_data["u"] = uh.x.array
    grid.set_active_scalars("u")
    
    plotter = pyvista.Plotter()
    plotter.add_mesh(grid, show_edges=True)
    plotter.show_axes()
    plotter.show()
    print("   ✓ Visualisierung geöffnet!")
except Exception as e:
    print(f"   ⚠ Visualisierung nicht verfügbar: {e}")

print("\n" + "=" * 60)
print("✓ Alle Berechnungen erfolgreich abgeschlossen!")
print("=" * 60)
