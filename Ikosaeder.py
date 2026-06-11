import bpy
import math

# ----------------------------
# Szene leeren
# ----------------------------
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# ----------------------------
# Parameter
# ----------------------------
phi = (1 + 5 ** 0.5) / 2   # Goldener Schnitt
cube_size = 2.0
rect_scale = 1.0

# ----------------------------
# Würfel
# ----------------------------
bpy.ops.mesh.primitive_cube_add(size=cube_size, location=(0,0,0))
cube = bpy.context.object
cube.name = "Outer_Cube"
cube.display_type = 'WIRE'

# ----------------------------
# Ikosaeder erzeugen
# ----------------------------
verts = [
    (-1,  phi, 0), (1,  phi, 0), (-1, -phi, 0), (1, -phi, 0),
    (0, -1,  phi), (0, 1,  phi), (0, -1, -phi), (0, 1, -phi),
    (phi, 0, -1), (phi, 0, 1), (-phi, 0, -1), (-phi, 0, 1)
]

faces = [
    (0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),
    (1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
    (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),
    (4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1)
]

mesh = bpy.data.meshes.new("Icosahedron")
mesh.from_pydata(verts, [], faces)
mesh.update()

ico_obj = bpy.data.objects.new("Icosahedron", mesh)
bpy.context.collection.objects.link(ico_obj)

# Skalieren, damit es in den Würfel passt
scale_factor = cube_size / (2 * phi)
ico_obj.scale = (scale_factor, scale_factor, scale_factor)

# ----------------------------
# Goldene Rechtecke
# ----------------------------
def create_golden_rectangle(axis):
    width = 2
    height = 2 * phi

    rect_verts = [
        (-width/2, -height/2, 0),
        ( width/2, -height/2, 0),
        ( width/2,  height/2, 0),
        (-width/2,  height/2, 0),
    ]

    rect_mesh = bpy.data.meshes.new("GoldenRectangle")
    rect_mesh.from_pydata(rect_verts, [], [(0,1,2,3)])
    rect_mesh.update()

    rect_obj = bpy.data.objects.new("GoldenRectangle", rect_mesh)
    bpy.context.collection.objects.link(rect_obj)

    if axis == 'XY':
        rect_obj.rotation_euler = (0,0,0)
    elif axis == 'YZ':
        rect_obj.rotation_euler = (0, math.pi/2, 0)
    elif axis == 'ZX':
        rect_obj.rotation_euler = (math.pi/2, 0, 0)

    rect_obj.display_type = 'WIRE'
    rect_obj.scale = (rect_scale, rect_scale, rect_scale)

create_golden_rectangle('XY')
create_golden_rectangle('YZ')
create_golden_rectangle('ZX')

# ----------------------------
# Material für Ikosaeder
# ----------------------------
mat = bpy.data.materials.new(name="IcoMaterial")
mat.use_nodes = True
bsdf = mat.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.2, 0.6, 1.0, 1)
ico_obj.data.materials.append(mat)

print("Modell erfolgreich erstellt.")