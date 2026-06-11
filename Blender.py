"""
E8-Root-Visualisierung für Blender.
Läuft nur innerhalb von Blender. Ausführen:
  blender --background --python Blender.py
oder in Blender: Scripting-Workspace -> Skript laden -> Ausführen
"""

try:
    import bpy
except ModuleNotFoundError:
    print("Fehler: Dieses Skript muss in Blender ausgeführt werden.")
    print("  Option 1: blender --background --python Blender.py")
    print("  Option 2: Blender öffnen -> Scripting -> Skript laden -> Run")
    raise SystemExit(1)

import numpy as np
from itertools import product, combinations

def create_material(name, color, emit=0):
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    bsdf = nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = (*color, 1.0)
        bsdf.inputs["Emission Strength"].default_value = emit
    return mat

def get_e8_roots():
    roots = []
    # Bosonen (112)
    for i, j in combinations(range(8), 2):
        for s1, s2 in product([-1, 1], repeat=2):
            root = np.zeros(8)
            root[i], root[j] = s1, s2
            roots.append((root, "Boson"))
    # Fermionen inkl. Neutrinos (128)
    np.random.seed(42)
    for signs in product([-0.5, 0.5], repeat=8):
        if sum(1 for s in signs if s < 0) % 2 == 0:
            type_label = "Neutrino" if np.random.rand() < 0.1 else "Fermion"
            roots.append((np.array(signs), type_label))
    return roots

def project_8d_to_3d(point_8d):
    """8D -> 3D Projektion."""
    p = np.asarray(point_8d)
    return tuple((p[:3] + p[3:6] * 0.5).tolist())

# Materialien vorbereiten
mats = {
    "Boson": create_material("Mat_Boson", (1.0, 0.8, 0.1), emit=5.0),
    "Fermion": create_material("Mat_Fermion", (0.1, 0.2, 0.8)),
    "Neutrino": create_material("Mat_Neutrino", (0.0, 1.0, 1.0), emit=10.0),
}

# Szene leeren (optional)
if bpy.context.scene.objects:
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()

# Instanziierung in der Szene
roots_data = get_e8_roots()
for i, (root, label) in enumerate(roots_data):
    pos = project_8d_to_3d(root)
    bpy.ops.mesh.primitive_uv_sphere_add(radius=0.1, location=pos)
    obj = bpy.context.active_object
    obj.name = f"{label}_{i}"
    obj.data.materials.append(mats[label])

    if label == "Boson":
        obj.scale = (0.5, 0.5, 0.5)

print(f"E8-Visualisierung: {len(roots_data)} Wurzeln erstellt.")
