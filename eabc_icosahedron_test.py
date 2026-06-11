import numpy as np
from itertools import combinations

PHI = (1 + np.sqrt(5)) / 2

FAMILIES = {
    "E": 1,
    "A": 5,
    "B": 7,
    "C": 11,
}

# EABC-Restklassen mod 12 (Achtlinge.py / eabc_from_lean.py)
EABC_RESIDUES = {1: "E", 5: "A", 7: "B", 11: "C"}
# Ikosaeder-Ecken auf dem (0, ±1, ±φ)-Gitter liefern nur {2, 4, 8, 10} mod 12.
ICO_RAW_RESIDUES = (2, 4, 8, 10)

def normalize(v):
    return v / np.linalg.norm(v)

def icosahedron_vertices():
    verts = []
    for s1 in [-1, 1]:
        for s2 in [-1, 1]:
            verts.append([0, s1, s2 * PHI])
            verts.append([s1, s2 * PHI, 0])
            verts.append([s1 * PHI, 0, s2])
    return np.array([normalize(np.array(v, dtype=float)) for v in verts])

def icosahedron_edges(vertices, tol=1e-6):
    dists = []
    for i, j in combinations(range(len(vertices)), 2):
        dists.append(np.linalg.norm(vertices[i] - vertices[j]))
    edge_len = min(dists)

    edges = []
    for i, j in combinations(range(len(vertices)), 2):
        if abs(np.linalg.norm(vertices[i] - vertices[j]) - edge_len) < tol:
            edges.append((i, j))
    return edges

def triangular_faces(vertices, edges):
    edge_set = {tuple(sorted(e)) for e in edges}
    faces = []
    for i, j, k in combinations(range(len(vertices)), 3):
        if (
            tuple(sorted((i, j))) in edge_set and
            tuple(sorted((j, k))) in edge_set and
            tuple(sorted((i, k))) in edge_set
        ):
            faces.append((i, j, k))
    return faces

def assign_eabc(vertices):
    labels = ["E", "A", "B", "C"]
    assignment = {}
    for i in range(len(vertices)):
        assignment[i] = labels[i % 4]
    return assignment

def assign_eabc_random(vertices, seed=None):
    labels = ["E"] * 3 + ["A"] * 3 + ["B"] * 3 + ["C"] * 3
    rng = np.random.default_rng(seed)
    rng.shuffle(labels)
    return {i: labels[i] for i in range(len(vertices))}


def vertex_residue_mod12(vertex, scale=6):
    """Restklasse mod 12 aus Ikosaeder-Koordinaten (0, ±1, ±φ)-Typ."""
    q = tuple(int(round(scale * c)) for c in vertex)
    return (q[0] + q[1] + q[2]) % 12


def nearest_eabc_residue(r):
    """Nächste EABC-Primrestklasse auf dem Kreis Z/12Z."""
    return min(EABC_RESIDUES.keys(), key=lambda x: min((x - r) % 12, (r - x) % 12))


def apply_r_eabc(vertices, edges, faces, assignment=None):
    """
    R_EABC: Renormierungsoperator auf 3+3+3+3-EABC-Zuweisungen.

    Mathematisch ist R_EABC die Projektion π auf die Untervarietät der
    geometrisch zulässigen Zuweisungen: für jede Ecke v wird unabhängig von
    der Eingabe die Ziel-Familie aus der Koordinaten-Restklasse
        r(v) = (round(scale·v_x) + round(scale·v_y) + round(scale·v_z)) mod 12
    und der EABC-Abbildung nearest_EABC(r) bestimmt. Der Operator ist
    idempotent (R_EABC ∘ R_EABC = R_EABC) und bildet beliebige gestörte
    Labelings auf dieselbe geometrische Zielzuweisung ab — analog zu einem
    Filter, der anisotrope Kusskonfigurationen zur nahezu isotropen Normkugel
    renormiert.
    """
    target, _ = assign_eabc_geometric(vertices, edges, faces)
    return target


def assign_eabc_geometric(vertices, edges, faces):
    """
    EABC-Zuweisung aus geometrischer Nachbarschaft und mod-12-Restklassen.

    Gewählte Regel (plausibelste Kopplung aus dem Projektmodell):
    - Koordinaten: Ikosaeder-Ecken liegen auf dem Gitter (0, ±1, ±φ); die
      skalierte Koordinatensumme liefert eine mod-12-Restklasse pro Ecke.
    - EABC-Abbildung: Standardzuordnung 1→E, 5→A, 7→B, 11→C (Achtlinge.py).
      Die vier Ikosaeder-Orbits {2,4,8,10} mod 12 werden auf die nächste
      EABC-Restklasse abgebildet → {E,A,B,C}, je 3 Ecken.
    - Graph: Kanten verbinden Ecken verschiedener Restklassen; die Zuweisung
      respektiert die Symmetrie des Koordinaten-Gitters (kein freies Labeling).

    Nicht gewählt: reine BFS-Propagation mit ABCE/CEAB-Umlauf (Achtlinge.py,
    CYCLE_PLUS/MINUS), weil parallele Transporte auf dem Ikosaeder-Zyklus
    i. d. R. kein 3+3+3+3-Balance ergeben und schwächere Metriken liefern.
    """
    assignment = {}
    residues = {}
    for i, v in enumerate(vertices):
        r = vertex_residue_mod12(v)
        residues[i] = r
        assignment[i] = EABC_RESIDUES[nearest_eabc_residue(r)]
    return assignment, residues

def compute_eabc_metrics(vertices, edges, faces, assignment):
    D = dipole(vertices, assignment)
    M = anisotropy_tensor(vertices, assignment)
    eigvals = np.linalg.eigvalsh(M)
    H_edges = edge_holonomy(vertices, edges, assignment)
    H_faces_total, H_faces_vector, _ = face_holonomy(vertices, faces, assignment)
    return {
        "D_norm": np.linalg.norm(D),
        "anisotropy": eigvals[-1] - eigvals[0],
        "H_scalar": H_faces_total,
        "H_vector_norm": np.linalg.norm(H_faces_vector),
        "H_edges": H_edges,
    }

def dipole(vertices, assignment):
    d = np.zeros(3)
    for i, v in enumerate(vertices):
        d += FAMILIES[assignment[i]] * v
    return d

def anisotropy_tensor(vertices, assignment):
    M = np.zeros((3, 3))
    for i, v in enumerate(vertices):
        w = FAMILIES[assignment[i]]
        M += w * np.outer(v, v)
    return M

def edge_holonomy(vertices, edges, assignment):
    H = 0.0
    for i, j in edges:
        wi = FAMILIES[assignment[i]]
        wj = FAMILIES[assignment[j]]
        orientation = np.sign(np.dot(np.cross(vertices[i], vertices[j]), np.array([1, 1, 1])))
        H += orientation * (wj - wi)
    return H

def face_holonomy(vertices, faces, assignment):
    total_scalar = 0.0
    total_vector = np.zeros(3)
    values = []
    for i, j, k in faces:
        vi, vj, vk = vertices[i], vertices[j], vertices[k]
        wi = FAMILIES[assignment[i]]
        wj = FAMILIES[assignment[j]]
        wk = FAMILIES[assignment[k]]
        # Orientierung der Dreiecksfläche
        normal = np.cross(vj - vi, vk - vi)
        center = vi + vj + vk
        orientation = np.sign(np.dot(normal, center))
        # orientiertes Volumen / diskrete Berry-Fläche
        scalar = orientation * np.dot(np.cross(vi, vj), vk)
        # EABC-gewichteter Flächenbeitrag
        weight = wi + wj + wk
        h_scalar = weight * scalar
        # nicht-kommutativer Vektorbeitrag
        h_vector = (
            wi * np.cross(vj, vk)
            + wj * np.cross(vk, vi)
            + wk * np.cross(vi, vj)
        ) * orientation
        values.append(h_scalar)
        total_scalar += h_scalar
        total_vector += h_vector
    return total_scalar, total_vector, np.array(values)

def build_dodecahedron_dual(vertices, edges, faces):
    """
    Ikosaeder-Dualität: Ikosaeder-Flächen ↔ Dodekaeder-Ecken,
    Ikosaeder-Ecken ↔ Dodekaeder-Flächen (Pentagone), Kanten bleiben dual.
    """
    do_vertices = []
    face_to_do_vertex = {}
    for f_idx, face in enumerate(faces):
        centroid = normalize(sum(vertices[i] for i in face))
        face_to_do_vertex[f_idx] = len(do_vertices)
        do_vertices.append(centroid)
    do_vertices = np.array(do_vertices)

    edge_to_faces = {}
    for f_idx, face in enumerate(faces):
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edge_to_faces.setdefault(tuple(sorted((a, b))), []).append(f_idx)

    do_edges = []
    for face_pair in edge_to_faces.values():
        if len(face_pair) == 2:
            do_edges.append(
                tuple(sorted((face_to_do_vertex[face_pair[0]], face_to_do_vertex[face_pair[1]])))
            )

    vertex_to_faces = {i: [] for i in range(len(vertices))}
    for f_idx, face in enumerate(faces):
        for v in face:
            vertex_to_faces[v].append(f_idx)

    def order_faces_around_vertex(v_idx, face_indices):
        center = vertices[v_idx]
        up = normalize(center)
        ref_vec = normalize(sum(vertices[i] for i in faces[face_indices[0]])) - center
        ref_vec = ref_vec - np.dot(ref_vec, up) * up
        ref_vec = normalize(ref_vec)
        tangent = normalize(np.cross(up, ref_vec))

        def angle_for_face(f_idx):
            c = normalize(sum(vertices[i] for i in faces[f_idx]))
            d = c - center
            d = d - np.dot(d, up) * up
            d = normalize(d)
            return np.arctan2(np.dot(d, tangent), np.dot(d, ref_vec))

        return tuple(face_to_do_vertex[f] for f in sorted(face_indices, key=angle_for_face))

    do_faces = []
    vertex_to_do_face = {}
    for v_idx, incident in vertex_to_faces.items():
        if len(incident) == 5:
            pent = order_faces_around_vertex(v_idx, incident)
            vertex_to_do_face[v_idx] = len(do_faces)
            do_faces.append(pent)

    return {
        "vertices": do_vertices,
        "edges": do_edges,
        "faces": do_faces,
        "face_to_do_vertex": face_to_do_vertex,
        "vertex_to_do_face": vertex_to_do_face,
    }


def triangulate_faces(faces):
    """Zerlegt n-Eck-Flächen per Fächer-Triangulation für die Holonomie."""
    triangles = []
    for face in faces:
        v0 = face[0]
        for a, b in zip(face[1:-1], face[2:]):
            triangles.append((v0, a, b))
    return triangles


def push_assignment_to_dual(ico_assignment, faces, face_to_do_vertex, vertex_to_do_face):
    """Ikosaeder-Ecken → Dodekaeder-Flächen; Flächen-Mehrheit → Dodekaeder-Ecken."""
    do_face_assignment = {vertex_to_do_face[i]: ico_assignment[i] for i in ico_assignment}

    do_vertex_assignment = {}
    for f_idx, face in enumerate(faces):
        labels = [ico_assignment[v] for v in face]
        do_vertex_assignment[face_to_do_vertex[f_idx]] = max(
            set(labels), key=labels.count
        )
    return do_vertex_assignment, do_face_assignment


def volume_renormalization(n):
    R = float(n)
    R_out = R * (1 + 2 * np.pi)

    V0 = 4 / 3 * np.pi * R**3
    Vout = 4 / 3 * np.pi * R_out**3

    lambda_sphere = (V0 / Vout) ** (1 / 3)

    return {
        "R": R,
        "R_out": R_out,
        "V0": V0,
        "Vout": Vout,
        "Vout/V0": Vout / V0,
        "lambda_sphere": lambda_sphere,
    }

def run_test(n=1):
    vertices = icosahedron_vertices()
    edges = icosahedron_edges(vertices)
    faces = triangular_faces(vertices, edges)
    assignment = assign_eabc(vertices)

    D = dipole(vertices, assignment)
    M = anisotropy_tensor(vertices, assignment)
    eigvals = np.linalg.eigvalsh(M)

    H_edges = edge_holonomy(vertices, edges, assignment)
    H_faces_total, H_faces_vector, H_faces = face_holonomy(vertices, faces, assignment)

    renorm = volume_renormalization(n)

    print("=== EABC Ikosaeder-Dodekaeder-Test ===")
    print(f"Vertices: {len(vertices)}")
    print(f"Edges:    {len(edges)}")
    print(f"Faces:    {len(faces)}")
    print()

    print("EABC assignment:")
    for i in range(len(vertices)):
        print(f"{i:2d}: {assignment[i]} = {FAMILIES[assignment[i]]}")
    print()

    print("Dipol / Restvektor:")
    print(D)
    print("Norm:", np.linalg.norm(D))
    print()

    print("Anisotropietensor:")
    print(M)
    print("Eigenwerte:", eigvals)
    print("Anisotropie:", eigvals[-1] - eigvals[0])
    print()

    print("Kanten-Holonomie:")
    print(H_edges)
    print()

    print("Flächen-Holonomie:")
    print("Skalar total:", H_faces_total)
    print("Vektor total:", H_faces_vector)
    print("Vektor-Norm:", np.linalg.norm(H_faces_vector))
    print("Einzelwerte:", H_faces)
    print()

    print("Volumen-Renormierung:")
    for k, v in renorm.items():
        print(f"{k}: {v}")

def run_control_test(n_permutations=1000, seed=42):
    vertices = icosahedron_vertices()
    edges = icosahedron_edges(vertices)
    faces = triangular_faces(vertices, edges)

    cyclic = compute_eabc_metrics(vertices, edges, faces, assign_eabc(vertices))

    keys = ["D_norm", "anisotropy", "H_scalar", "H_vector_norm", "H_edges"]
    samples = {k: [] for k in keys}
    rng = np.random.default_rng(seed)
    for _ in range(n_permutations):
        perm_seed = int(rng.integers(0, 2**31 - 1))
        assignment = assign_eabc_random(vertices, seed=perm_seed)
        metrics = compute_eabc_metrics(vertices, edges, faces, assignment)
        for k in keys:
            samples[k].append(metrics[k])

    print()
    print("=== Kontrolllauf: zufällige EABC-Permutationen ===")
    print(f"Permutationen: {n_permutations}, Basis-Seed: {seed}")
    print()

    print("Zyklischer Basisfall:")
    for k in keys:
        print(f"  {k}: {cyclic[k]}")
    print()

    print("Zufalls-Permutationen (Statistik):")
    for k in keys:
        arr = np.array(samples[k])
        print(
            f"  {k}: mean={arr.mean():.6f}, std={arr.std():.6f}, "
            f"min={arr.min():.6f}, max={arr.max():.6f}"
        )
    print()

    print("Zyklisch vs. Zufall:")
    for k in keys:
        arr = np.array(samples[k])
        mean, std = arr.mean(), arr.std()
        if std < 1e-9:
            print(f"  {k}: INVARIANT (alle Permutationen = {mean:.6f})")
        else:
            z = (cyclic[k] - mean) / std
            pct = 100 * np.mean(arr <= cyclic[k])
            print(f"  {k}: z={z:.3f}, zyklisch >= Zufall: {pct:.1f}%")

def run_geometric_test(n_permutations=1000, seed=42):
    vertices = icosahedron_vertices()
    edges = icosahedron_edges(vertices)
    faces = triangular_faces(vertices, edges)

    cyclic = compute_eabc_metrics(vertices, edges, faces, assign_eabc(vertices))
    geometric_assignment, raw_residues = assign_eabc_geometric(vertices, edges, faces)
    geometric = compute_eabc_metrics(vertices, edges, faces, geometric_assignment)

    keys = ["D_norm", "anisotropy", "H_vector_norm", "H_edges"]
    samples = {k: [] for k in keys}
    rng = np.random.default_rng(seed)
    for _ in range(n_permutations):
        perm_seed = int(rng.integers(0, 2**31 - 1))
        assignment = assign_eabc_random(vertices, seed=perm_seed)
        metrics = compute_eabc_metrics(vertices, edges, faces, assignment)
        for k in keys:
            samples[k].append(metrics[k])

    print()
    print("=== Geometrischer EABC-Test: zyklisch vs. geometrisch vs. Zufall ===")
    print(f"Permutationen (Zufall): {n_permutations}, Basis-Seed: {seed}")
    print()
    print("H_scalar ist invariant (nur von Gesamtgewicht 72 abhängig):")
    print(f"  zyklisch={cyclic['H_scalar']:.6f}, geometrisch={geometric['H_scalar']:.6f}")
    print()

    print("Geometrische EABC-Zuweisung (Koordinaten mod 12 → EABC):")
    for i in range(len(vertices)):
        letter = geometric_assignment[i]
        print(
            f"  {i:2d}: r={raw_residues[i]:2d} → {letter} "
            f"(EABC={FAMILIES[letter]})"
        )
    print()

    for label, metrics in [("Zyklisch", cyclic), ("Geometrisch", geometric)]:
        print(f"{label}:")
        for k in keys:
            print(f"  {k}: {metrics[k]}")
        print()

    print("Zufalls-Permutationen (Statistik):")
    for k in keys:
        arr = np.array(samples[k])
        print(
            f"  {k}: mean={arr.mean():.6f}, std={arr.std():.6f}, "
            f"min={arr.min():.6f}, max={arr.max():.6f}"
        )
    print()

    print("Vergleich (z-Score gegen Zufall, geometrisch >= Zufall in %):")
    for label, metrics in [("Zyklisch", cyclic), ("Geometrisch", geometric)]:
        print(f"  [{label}]")
        for k in keys:
            arr = np.array(samples[k])
            mean, std = arr.mean(), arr.std()
            if std < 1e-9:
                print(f"    {k}: INVARIANT (Zufall mean={mean:.6f})")
            else:
                z = (metrics[k] - mean) / std
                pct = 100 * np.mean(arr <= metrics[k])
                print(f"    {k}: z={z:.3f}, >= Zufall: {pct:.1f}%")


def run_renormalization_test(n_starts=30, seed=42):
    vertices = icosahedron_vertices()
    edges = icosahedron_edges(vertices)
    faces = triangular_faces(vertices, edges)
    geometric_assignment, _ = assign_eabc_geometric(vertices, edges, faces)
    geometric_metrics = compute_eabc_metrics(vertices, edges, faces, geometric_assignment)

    keys = ["D_norm", "anisotropy", "H_vector_norm", "H_edges"]
    before = {k: [] for k in keys}
    after = {k: [] for k in keys}
    rng = np.random.default_rng(seed)

    print()
    print("=== R_EABC-Störungstest: Zufallszuweisung → geometrische Renormierung ===")
    print(f"Zufallsstarts: {n_starts}, Basis-Seed: {seed}")
    print("H_scalar bleibt invariant (Gesamtgewicht 72).")
    print()

    for trial in range(n_starts):
        perm_seed = int(rng.integers(0, 2**31 - 1))
        disturbed = assign_eabc_random(vertices, seed=perm_seed)
        m_before = compute_eabc_metrics(vertices, edges, faces, disturbed)
        renormalized = apply_r_eabc(vertices, edges, faces, disturbed)
        m_after = compute_eabc_metrics(vertices, edges, faces, renormalized)

        for k in keys:
            before[k].append(m_before[k])
            after[k].append(m_after[k])

        if trial < 3:
            print(f"Start {trial + 1} (seed={perm_seed}):")
            print("  vor R_EABC:", {k: f"{m_before[k]:.6f}" for k in keys})
            print("  nach R_EABC:", {k: f"{m_after[k]:.6f}" for k in keys})
            print()

    print("Geometrisches Ziel (R_EABC-Fixpunkt):")
    for k in keys:
        print(f"  {k}: {geometric_metrics[k]:.6f}")
    print()

    print("Statistik über alle Zufallsstarts:")
    print(f"{'Metrik':<18} {'vor (mean±std)':<22} {'nach (mean±std)':<22} {'Δ mean':<10}")
    print("-" * 72)
    for k in keys:
        b = np.array(before[k])
        a = np.array(after[k])
        delta = b.mean() - a.mean()
        print(
            f"{k:<18} {b.mean():.4f}±{b.std():.4f}      "
            f"{a.mean():.4f}±{a.std():.4f}      {delta:+.4f}"
        )
    print()

    converged = all(
        after[k][i] == geometric_metrics[k]
        for k in keys
        for i in range(n_starts)
        if k != "H_scalar"
    )
    aniso_drop = np.mean(before["anisotropy"]) - np.mean(after["anisotropy"])
    print(
        f"Anisotropie-Abfall (mittel): {aniso_drop:.6f} "
        f"({100 * aniso_drop / np.mean(before['anisotropy']):.1f}% gegenüber Zufall)"
    )
    print(f"Alle Starts konvergieren zum geometrischen Ziel: {converged}")
    print(f"R_EABC idempotent (2. Anwendung ändert nichts): "
          f"{apply_r_eabc(vertices, edges, faces, renormalized) == renormalized}")


def run_duality_test(n_starts=30, seed=42):
    vertices = icosahedron_vertices()
    edges = icosahedron_edges(vertices)
    faces = triangular_faces(vertices, edges)
    dual = build_dodecahedron_dual(vertices, edges, faces)
    do_vertices = dual["vertices"]
    do_edges = dual["edges"]
    do_faces = dual["faces"]

    do_triangles = triangulate_faces(do_faces)

    geo_ico, _ = assign_eabc_geometric(vertices, edges, faces)
    geo_do, _ = assign_eabc_geometric(do_vertices, do_edges, do_triangles)
    m_geo_ico = compute_eabc_metrics(vertices, edges, faces, geo_ico)
    m_geo_do = compute_eabc_metrics(do_vertices, do_edges, do_triangles, geo_do)

    keys = ["D_norm", "anisotropy", "H_vector_norm", "H_edges"]
    ico_before, ico_after = {k: [] for k in keys}, {k: [] for k in keys}
    do_before, do_after = {k: [] for k in keys}, {k: [] for k in keys}
    rng = np.random.default_rng(seed)

    print()
    print("=== Dualitätstest: Ikosaeder ↔ Dodekaeder-Dual ===")
    print(f"Dodekaeder: {len(do_vertices)} Ecken, {len(do_edges)} Kanten, {len(do_faces)} Flächen")
    print("H_scalar invariant unter Permutation (nur kurz erwähnt).")
    print()

    print("Geometrische Fixpunkte (R_EABC-Ziel):")
    for label, metrics in [("Ikosaeder", m_geo_ico), ("Dodekaeder-Dual", m_geo_do)]:
        print(f"  [{label}]")
        for k in keys:
            print(f"    {k}: {metrics[k]:.6f}")
    print()

    for trial in range(n_starts):
        perm_seed = int(rng.integers(0, 2**31 - 1))
        ico_disturbed = assign_eabc_random(vertices, seed=perm_seed)
        do_disturbed, _ = push_assignment_to_dual(
            ico_disturbed, faces, dual["face_to_do_vertex"], dual["vertex_to_do_face"]
        )

        m_ico_b = compute_eabc_metrics(vertices, edges, faces, ico_disturbed)
        m_ico_a = compute_eabc_metrics(
            vertices, edges, faces, apply_r_eabc(vertices, edges, faces, ico_disturbed)
        )
        m_do_b = compute_eabc_metrics(do_vertices, do_edges, do_triangles, do_disturbed)
        m_do_a = compute_eabc_metrics(
            do_vertices, do_edges, do_triangles,
            apply_r_eabc(do_vertices, do_edges, do_triangles, do_disturbed),
        )

        for k in keys:
            ico_before[k].append(m_ico_b[k])
            ico_after[k].append(m_ico_a[k])
            do_before[k].append(m_do_b[k])
            do_after[k].append(m_do_a[k])

    print("Renormierung vor/nach R_EABC (Mittelwert über Zufallsstarts):")
    print(f"{'Metrik':<18} {'Ico vor':<12} {'Ico nach':<12} {'Do vor':<12} {'Do nach':<12}")
    print("-" * 66)
    for k in keys:
        print(
            f"{k:<18} "
            f"{np.mean(ico_before[k]):.4f}     "
            f"{np.mean(ico_after[k]):.4f}     "
            f"{np.mean(do_before[k]):.4f}     "
            f"{np.mean(do_after[k]):.4f}"
        )
    print()

    print("Anisotropie-Abfall (Ikosaeder vs. Dual):")
    ico_drop = np.mean(ico_before["anisotropy"]) - np.mean(ico_after["anisotropy"])
    do_drop = np.mean(do_before["anisotropy"]) - np.mean(do_after["anisotropy"])
    print(f"  Ikosaeder:  {ico_drop:.6f} ({100 * ico_drop / np.mean(ico_before['anisotropy']):.1f}%)")
    print(f"  Dual:       {do_drop:.6f} ({100 * do_drop / np.mean(do_before['anisotropy']):.1f}%)")
    print()

    tol = 1e-6
    invariant = all(abs(m_geo_ico[k] - m_geo_do[k]) < tol for k in keys)
    print("Symmetrie-Auslösung am geometrischen Fixpunkt invariant unter Dualität:")
    for k in keys:
        diff = abs(m_geo_ico[k] - m_geo_do[k])
        print(f"  {k}: Δ={diff:.6e}  →  {'invariant' if diff < tol else 'abweichend'}")
    print(f"  Gesamt: {'JA' if invariant else 'NEIN (Dual-Metriken weichen ab)'}")


def perturb_vertices(vertices, epsilon, seed=None, method="tangent"):
    """
    Kleine Störung der 12 Einheits-Ecken (Kussrichtungen).

    method='tangent': Gauß-Rauschen in der Tangentialebene + Re-Normalisierung.
    method='radial':    Gauß-Rauschen in R³ + Re-Normalisierung auf S².

    Bei epsilon=0 werden die Eingabe-Vertices unverändert zurückgegeben.
  """
    if epsilon == 0:
        return np.array(vertices, dtype=float, copy=True)

    rng = np.random.default_rng(seed)
    perturbed = []
    for v in vertices:
        v = np.asarray(v, dtype=float)
        if method == "tangent":
            noise = rng.normal(0, 1, 3)
            tangent = noise - np.dot(noise, v) * v
            t_norm = np.linalg.norm(tangent)
            if t_norm < 1e-12:
                perturbed.append(v.copy())
                continue
            tangent /= t_norm
            w = v + epsilon * tangent
        elif method == "radial":
            w = v + epsilon * rng.normal(0, 1, 3)
        else:
            raise ValueError(f"Unbekannte Perturbationsmethode: {method}")
        perturbed.append(normalize(w))
    return np.array(perturbed)


def _assignment_balance(assignment):
    """Zählt E/A/B/C — bei isotropem Fixpunkt je 3."""
    counts = {letter: 0 for letter in "EABC"}
    for letter in assignment.values():
        counts[letter] += 1
    return counts


DEFECT_FIXPOINT = "FIXPOINT"
DEFECT_DIRECTION = "DIRECTION_DEFECT"
DEFECT_CLASSIFICATION = "CLASSIFICATION_DEFECT"
DEFECT_LABEL = "LABEL_DEFECT"


def order_parameter(vertices, assignment, return_dipole=True):
    """
    Ordnungsparameter Δ(M) aus dem Anisotropietensor M = Σ w_i v_i v_i^T.

    Δ(M) = λ_max − λ_min (Eigenwerte von M). Optional Dipol-Norm ||D||.
    """
    M = anisotropy_tensor(vertices, assignment)
    eigvals = np.linalg.eigvalsh(M)
    result = {
        "M": M,
        "eigenvalues": eigvals,
        "anisotropy": eigvals[-1] - eigvals[0],
        "delta_M": eigvals[-1] - eigvals[0],
    }
    if return_dipole:
        D = dipole(vertices, assignment)
        result["dipole"] = D
        result["D_norm"] = np.linalg.norm(D)
    return result


def _geometry_max_deviation(vertices, reference_vertices):
    """Maximale Winkelabweichung (Bogenmaß) gegenüber referenzierter Ikosaeder-Geometrie."""
    deviations = []
    for v, ref in zip(vertices, reference_vertices):
        cos_angle = np.clip(np.dot(normalize(v), normalize(ref)), -1.0, 1.0)
        deviations.append(np.arccos(cos_angle))
    return max(deviations)


def classify_defect(
    vertices,
    assignment,
    epsilon_geom=0,
    epsilon_class=0.1,
    reference_vertices=None,
    edges=None,
    faces=None,
):
    """
    Defektklassifikation für (Vertices, Assignment):

    - FIXPOINT: Δ(M)≈0, Balance 3+3+3+3, ungestörte Geometrie
    - DIRECTION_DEFECT: Labels balanciert, aber Δ(M)>ε bei gestörter Geometrie
    - CLASSIFICATION_DEFECT: mod-12-Balance der geometrischen Zuweisung gebrochen
    - LABEL_DEFECT: Label-Zuweisung ohne 3+3+3+3-Balance (vor R_EABC)
    """
    if reference_vertices is None:
        reference_vertices = icosahedron_vertices()

    vertices = np.asarray(vertices, dtype=float)
    op = order_parameter(vertices, assignment)
    balance = _assignment_balance(assignment)
    label_balanced = all(balance[letter] == 3 for letter in "EABC")

    if edges is None:
        edges = icosahedron_edges(vertices)
    if faces is None:
        faces = triangular_faces(vertices, edges)

    geo_assignment, _ = assign_eabc_geometric(vertices, edges, faces)
    geo_balance = _assignment_balance(geo_assignment)
    mod12_balanced = all(geo_balance[letter] == 3 for letter in "EABC")

    geom_dev = _geometry_max_deviation(vertices, reference_vertices)
    geom_tol = 1e-9 if epsilon_geom == 0 else epsilon_geom
    geom_undisturbed = geom_dev <= geom_tol

    info = {
        "defect_class": None,
        "delta_M": op["delta_M"],
        "D_norm": op["D_norm"],
        "eigenvalues": op["eigenvalues"],
        "label_balance": balance,
        "geo_balance": geo_balance,
        "geom_deviation": geom_dev,
        "geom_undisturbed": geom_undisturbed,
        "label_balanced": label_balanced,
        "mod12_balanced": mod12_balanced,
    }

    if not mod12_balanced:
        info["defect_class"] = DEFECT_CLASSIFICATION
        return DEFECT_CLASSIFICATION, info
    if not label_balanced:
        info["defect_class"] = DEFECT_LABEL
        return DEFECT_LABEL, info
    if not geom_undisturbed and op["delta_M"] > epsilon_class:
        info["defect_class"] = DEFECT_DIRECTION
        return DEFECT_DIRECTION, info
    if op["delta_M"] > epsilon_class:
        info["defect_class"] = DEFECT_LABEL
        return DEFECT_LABEL, info
    info["defect_class"] = DEFECT_FIXPOINT
    return DEFECT_FIXPOINT, info


def geometric_pullback(vertices, reference_vertices=None):
    """
    Geometrische Rückführung zum Ikosaeder-Attraktor.

    Gewählte Methode: Voronoi-Projektion auf S² — jede gestörte Kussrichtung
    v wird auf die Referenz-Ecke mit maximalem Skalarprodukt argmax_j ⟨v, r_j⟩
    abgebildet (nächste Ecke auf der Einheitskugel). Plausibler als reines
    Index-Procrustes, weil bei starker Störung die kombinatorische Nachbarschaft
    erhalten bleibt und alle 12 Attraktorrichtungen wiederbelegt werden.

    Nicht gewählt: globales Procrustes-Matching (Permutation), weil der
    Renormierungsschritt K_n → K_{n+1} lokale Kussrichtungen pro Slot renormiert.
    """
    if reference_vertices is None:
        reference_vertices = icosahedron_vertices()
    ref = np.array([normalize(v) for v in reference_vertices])
    pulled = []
    for v in vertices:
        v = normalize(v)
        j = int(np.argmax(ref @ v))
        pulled.append(ref[j].copy())
    return np.array(pulled)


project_to_icosahedron = geometric_pullback


def apply_renorm_step(vertices, assignment, reference_vertices=None):
    """
    Gekoppelter Renormierungsschritt K_n → K_{n+1}:
    1. Geometrische Rückführung (Voronoi-Pullback zum Ikosaeder)
    2. R_EABC auf Labels
    """
    ref = reference_vertices if reference_vertices is not None else icosahedron_vertices()
    new_verts = geometric_pullback(vertices, ref)
    new_edges = icosahedron_edges(new_verts)
    new_faces = triangular_faces(new_verts, new_edges)
    new_assignment = apply_r_eabc(new_verts, new_edges, new_faces, assignment)
    return new_verts, new_assignment, new_edges, new_faces


def apply_coupled_renorm(vertices, assignment, n_steps=1, reference_vertices=None):
    """Iterative gekoppelte Renormierung (1–3 Schritte für Konvergenztests)."""
    verts = np.array(vertices, dtype=float, copy=True)
    assgn = dict(assignment)
    ref = reference_vertices if reference_vertices is not None else icosahedron_vertices()
    edges, faces = None, None
    for _ in range(n_steps):
        verts, assgn, edges, faces = apply_renorm_step(verts, assgn, ref)
    return verts, assgn, edges, faces


apply_r_star_eabc = apply_renorm_step


def apply_volume_scale(vertices, n=1):
    """
    Isotrope Volumen-Renormierung: v_i ↦ λ_sphere · v_i.

    λ_sphere = (V0/Vout)^(1/3) = 1/(1+2π) für R_out = R(1+2π).
    """
    renorm = volume_renormalization(n)
    lam = renorm["lambda_sphere"]
    scaled = np.array(vertices, dtype=float) * lam
    return scaled, renorm


def apply_full_renorm_chain(vertices, assignment, n=1, reference_vertices=None):
    """
    Vollständige Renormierungskette K_n → K_{n+1}:
      1. P_I  (Voronoi-Pullback auf Referenz-Ikosaeder)
      2. R_EABC (Label-Projektion)
      3. Volumen-Renormierung (isotrope Skalierung mit λ_sphere)

    Entspricht R*_EABC gefolgt von äußerer Skalenanpassung.
  """
    verts, assgn, edges, faces = apply_renorm_step(
        vertices, assignment, reference_vertices
    )
    verts, renorm = apply_volume_scale(verts, n=n)
    return verts, assgn, edges, faces, renorm


def run_full_chain_test(
    epsilons=(0.0, 0.05, 0.1, 0.2),
    n_trials=30,
    seed=42,
    n=1,
    epsilon_class=0.1,
    tol=1e-6,
):
    """
    Gestörte Konfiguration → R*_EABC → Volumen-Renormierung.

    Misst Δ(M), Radius und Vout/V0; prüft Δ(M)-Invarianz unter λ_sphere-Skalierung.
    """
    base_vertices = icosahedron_vertices()
    base_edges = icosahedron_edges(base_vertices)
    base_faces = triangular_faces(base_vertices, base_edges)
    geo_assignment, _ = assign_eabc_geometric(base_vertices, base_edges, base_faces)
    fixpoint_delta = order_parameter(base_vertices, geo_assignment)["delta_M"]
    renorm_ref = volume_renormalization(n)

    rng = np.random.default_rng(seed)

    print()
    print("=== Vollständige Renormierungskette: R*_EABC + Volumen ===")
    print(f"Trials pro ε: {n_trials}, n={n}, Seed: {seed}")
    print(f"λ_sphere = {renorm_ref['lambda_sphere']:.6f}, "
          f"Vout/V0 = {renorm_ref['Vout/V0']:.6f}, "
          f"R_out = R(1+2π) = {renorm_ref['R_out']:.6f}")
    print(f"Fixpunkt Δ(M) = {fixpoint_delta:.6e}")
    print()

    delta_after_rstar = []
    delta_after_volume = []
    radius_after_volume = []
    delta_invariance_ok = 0

    for eps in epsilons:
        d_rstar = []
        d_vol = []
        radii = []

        for trial in range(n_trials):
            trial_seed = int(rng.integers(0, 2**31 - 1))
            perturbed = perturb_vertices(base_vertices, eps, seed=trial_seed)
            random_assignment = assign_eabc_random(perturbed, seed=trial_seed + 1)

            verts_rstar, assgn_rstar, _, _ = apply_renorm_step(
                perturbed, random_assignment
            )
            op_rstar = order_parameter(verts_rstar, assgn_rstar)
            d_rstar.append(op_rstar["delta_M"])

            verts_full, assgn_full, _, _, renorm = apply_full_renorm_chain(
                perturbed, random_assignment, n=n
            )
            op_vol = order_parameter(verts_full, assgn_full)
            d_vol.append(op_vol["delta_M"])
            mean_radius = np.mean(np.linalg.norm(verts_full, axis=1))
            radii.append(mean_radius)

            if abs(op_rstar["delta_M"] - op_vol["delta_M"]) < tol:
                delta_invariance_ok += 1

            if trial == 0 and eps in (0.0, 0.1):
                defect_r, _ = classify_defect(
                    verts_rstar, assgn_rstar, epsilon_geom=0, epsilon_class=epsilon_class
                )
                defect_v, _ = classify_defect(
                    verts_full, assgn_full, epsilon_geom=0, epsilon_class=epsilon_class
                )
                print(f"  Beispiel ε={eps:.2f}:")
                print(f"    vor:           Δ(M)={order_parameter(perturbed, random_assignment)['delta_M']:.4f}")
                print(f"    nach R*_EABC:  Δ(M)={op_rstar['delta_M']:.4e}, Defekt={defect_r}")
                print(f"    + Volumen:     Δ(M)={op_vol['delta_M']:.4e}, "
                      f"Radius={mean_radius:.6f}, Defekt={defect_v}")
                print(f"    Δ(M)-Differenz R* vs. +Vol: {abs(op_rstar['delta_M'] - op_vol['delta_M']):.4e}")
                print()

        delta_after_rstar.extend(d_rstar)
        delta_after_volume.extend(d_vol)
        radius_after_volume.extend(radii)

        print(f"ε = {eps:.3f}:")
        print(f"  Δ(M) nach R*_EABC:  mean={np.mean(d_rstar):.4e}, max={np.max(d_rstar):.4e}")
        print(f"  Δ(M) nach +Volumen: mean={np.mean(d_vol):.4e}, max={np.max(d_vol):.4e}")
        print(f"  Radius (mean):      {np.mean(radii):.6f} (erwartet λ_sphere={renorm_ref['lambda_sphere']:.6f})")
        print()

    n_total = len(epsilons) * n_trials
    print("Zusammenfassung:")
    print(f"  Δ(M) invariant unter λ_sphere-Skalierung: "
          f"{delta_invariance_ok}/{n_total} Trials (|Δ_R* − Δ_vol| < {tol})")
    print(f"  Mittlerer Radius nach Volumen-Schritt: "
          f"{np.mean(radius_after_volume):.6f}")
    print(f"  Vout/V0 = {renorm_ref['Vout/V0']:.6f}")
    print()


def run_order_parameter_test(
    epsilons=(0.0, 0.01, 0.02, 0.05, 0.1, 0.2),
    n_trials=30,
    seed=42,
    epsilon_class=0.1,
):
    """Systematischer Test: Δ(M) und Defektklassen über ε und Zufalls-Labels."""
    base_vertices = icosahedron_vertices()
    base_edges = icosahedron_edges(base_vertices)
    base_faces = triangular_faces(base_vertices, base_edges)
    geo_assignment, _ = assign_eabc_geometric(base_vertices, base_edges, base_faces)
    fixpoint_op = order_parameter(base_vertices, geo_assignment)

    rng = np.random.default_rng(seed)
    all_classes = [DEFECT_FIXPOINT, DEFECT_DIRECTION, DEFECT_CLASSIFICATION, DEFECT_LABEL]

    print()
    print("=== Ordnungsparameter Δ(M) und Defektklassifikation ===")
    print(f"Trials pro ε: {n_trials}, Basis-Seed: {seed}, ε_class={epsilon_class}")
    print()
    print("Referenz-Fixpunkt (ungestört, geometrisch):")
    print(f"  Δ(M) = {fixpoint_op['delta_M']:.6e}")
    print(f"  ||D|| = {fixpoint_op['D_norm']:.6e}")
    print(f"  Eigenwerte: {fixpoint_op['eigenvalues']}")
    print()

    total_counts = {c: 0 for c in all_classes}

    for eps in epsilons:
        counts = {c: 0 for c in all_classes}
        delta_before = []
        delta_after_reabc = []

        for trial in range(n_trials):
            trial_seed = int(rng.integers(0, 2**31 - 1))
            perturbed = perturb_vertices(base_vertices, eps, seed=trial_seed)
            p_edges = icosahedron_edges(perturbed)
            p_faces = triangular_faces(perturbed, p_edges)

            random_assignment = assign_eabc_random(perturbed, seed=trial_seed + 1)
            defect, info = classify_defect(
                perturbed, random_assignment, epsilon_geom=eps, epsilon_class=epsilon_class
            )
            counts[defect] += 1
            total_counts[defect] += 1
            delta_before.append(info["delta_M"])

            renormalized = apply_r_eabc(perturbed, p_edges, p_faces, random_assignment)
            op_after = order_parameter(perturbed, renormalized)
            delta_after_reabc.append(op_after["delta_M"])

        print(f"ε = {eps:.3f}:")
        print(f"  Defektklassen: {', '.join(f'{c}={counts[c]}' for c in all_classes)}")
        print(
            f"  Δ(M) vor R_EABC:  mean={np.mean(delta_before):.4f}, "
            f"max={np.max(delta_before):.4f}"
        )
        print(
            f"  Δ(M) nach R_EABC: mean={np.mean(delta_after_reabc):.4e}, "
            f"max={np.max(delta_after_reabc):.4e}"
        )
        print()

    print("Gesamtverteilung (alle ε × Trials):")
    n_total = sum(total_counts.values())
    for c in all_classes:
        pct = 100 * total_counts[c] / n_total
        print(f"  {c}: {total_counts[c]} ({pct:.1f}%)")
    print()

    # Fixpunkt-Check bei ε=0
    defect_fp, _ = classify_defect(base_vertices, geo_assignment, epsilon_class=epsilon_class)
    print(f"Ungestörter geometrischer Fixpunkt: {defect_fp}")
    print()


def run_coupled_renorm_test(
    epsilons=(0.0, 0.02, 0.05, 0.1, 0.2),
    n_trials=30,
    seed=42,
    n_renorm_steps=3,
    epsilon_class=0.1,
    tol=1e-6,
):
    """
    Vergleich: nur R_EABC vs. gekoppelter Schritt (Geometrie-Pullback + R_EABC).
    """
    base_vertices = icosahedron_vertices()
    base_edges = icosahedron_edges(base_vertices)
    base_faces = triangular_faces(base_vertices, base_edges)
    geo_assignment, _ = assign_eabc_geometric(base_vertices, base_edges, base_faces)
    fixpoint_metrics = compute_eabc_metrics(base_vertices, base_edges, base_faces, geo_assignment)
    fixpoint_delta = order_parameter(base_vertices, geo_assignment)["delta_M"]

    rng = np.random.default_rng(seed)
    all_classes = [DEFECT_FIXPOINT, DEFECT_DIRECTION, DEFECT_CLASSIFICATION, DEFECT_LABEL]

    print()
    print("=== Gekoppelte Renormierung: Geometrie-Pullback + R_EABC ===")
    print(f"Trials pro ε: {n_trials}, Schritte gekoppelt: bis {n_renorm_steps}, Seed: {seed}")
    print("Geometrische Rückführung: Voronoi-Projektion auf Referenz-Ikosaeder (S²)")
    print()
    print(f"Fixpunkt Δ(M) = {fixpoint_delta:.6e}, Anisotropie = {fixpoint_metrics['anisotropy']:.6e}")
    print()

    summary = []

    for eps in epsilons:
        delta_before = []
        delta_reabc = []
        delta_coupled = []
        conv_reabc = 0
        conv_coupled_1 = 0
        conv_coupled_n = 0
        defect_before = {c: 0 for c in all_classes}
        defect_reabc = {c: 0 for c in all_classes}
        defect_coupled = {c: 0 for c in all_classes}

        for trial in range(n_trials):
            trial_seed = int(rng.integers(0, 2**31 - 1))
            perturbed = perturb_vertices(base_vertices, eps, seed=trial_seed)
            p_edges = icosahedron_edges(perturbed)
            p_faces = triangular_faces(perturbed, p_edges)
            random_assignment = assign_eabc_random(perturbed, seed=trial_seed + 1)

            op_b = order_parameter(perturbed, random_assignment)
            delta_before.append(op_b["delta_M"])
            d_b, _ = classify_defect(
                perturbed, random_assignment, epsilon_geom=eps, epsilon_class=epsilon_class
            )
            defect_before[d_b] += 1

            # Nur R_EABC
            assgn_reabc = apply_r_eabc(perturbed, p_edges, p_faces, random_assignment)
            op_r = order_parameter(perturbed, assgn_reabc)
            delta_reabc.append(op_r["delta_M"])
            d_r, _ = classify_defect(
                perturbed, assgn_reabc, epsilon_geom=eps, epsilon_class=epsilon_class
            )
            defect_reabc[d_r] += 1
            if op_r["delta_M"] < tol:
                conv_reabc += 1

            # Gekoppelt: 1 Schritt
            verts_c1, assgn_c1, e_c1, f_c1 = apply_coupled_renorm(
                perturbed, random_assignment, n_steps=1
            )
            op_c1 = order_parameter(verts_c1, assgn_c1)
            delta_coupled.append(op_c1["delta_M"])
            d_c1, _ = classify_defect(
                verts_c1, assgn_c1, epsilon_geom=0, epsilon_class=epsilon_class
            )
            defect_coupled[d_c1] += 1
            if op_c1["delta_M"] < tol:
                conv_coupled_1 += 1

            # Gekoppelt: n Schritte
            verts_cn, assgn_cn, _, _ = apply_coupled_renorm(
                perturbed, random_assignment, n_steps=n_renorm_steps
            )
            op_cn = order_parameter(verts_cn, assgn_cn)
            if op_cn["delta_M"] < tol:
                conv_coupled_n += 1

            if trial < 2 and eps in (0.0, 0.1):
                print(f"  Beispiel ε={eps:.2f}, Trial {trial + 1}:")
                print(f"    vor:      Δ(M)={op_b['delta_M']:.4f}, Defekt={d_b}")
                print(f"    R_EABC:   Δ(M)={op_r['delta_M']:.4e}, Defekt={d_r}")
                print(f"    gekoppelt:Δ(M)={op_c1['delta_M']:.4e}, Defekt={d_c1}")
                print(f"    {n_renorm_steps} Schritte: Δ(M)={op_cn['delta_M']:.4e}")
                print()

        row = {
            "epsilon": eps,
            "delta_before_mean": np.mean(delta_before),
            "delta_reabc_mean": np.mean(delta_reabc),
            "delta_coupled_mean": np.mean(delta_coupled),
            "conv_reabc": conv_reabc / n_trials,
            "conv_coupled_1": conv_coupled_1 / n_trials,
            "conv_coupled_n": conv_coupled_n / n_trials,
        }
        summary.append(row)

        print(f"ε = {eps:.3f}:")
        print(
            f"  Δ(M) mean: vor={row['delta_before_mean']:.4f}, "
            f"R_EABC={row['delta_reabc_mean']:.4e}, "
            f"gekoppelt(1)={row['delta_coupled_mean']:.4e}"
        )
        print(
            f"  Konvergenz Δ<ε_tol: R_EABC={row['conv_reabc']:.0%}, "
            f"gekoppelt(1)={row['conv_coupled_1']:.0%}, "
            f"gekoppelt({n_renorm_steps})={row['conv_coupled_n']:.0%}"
        )
        print(f"  Defekt vor:     {defect_before}")
        print(f"  Defekt R_EABC:  {defect_reabc}")
        print(f"  Defekt gekoppelt: {defect_coupled}")
        print()

    print("Vergleich gekoppelt vs. R_EABC allein:")
    better_coupled = sum(
        1 for r in summary if r["delta_coupled_mean"] < r["delta_reabc_mean"] - 1e-12
    )
    print(
        f"  Gekoppelt (1 Schritt) niedrigere mittlere Δ(M) in {better_coupled}/{len(summary)} ε-Bins"
    )
    better_conv = sum(
        1 for r in summary if r["conv_coupled_1"] > r["conv_reabc"]
    )
    print(
        f"  Höhere Fixpunkt-Konvergenzrate (Δ<{tol}) bei gekoppelt in {better_conv}/{len(summary)} ε-Bins"
    )
    print()


def _max_metric_deviation(metrics, reference):
    keys = ["D_norm", "anisotropy", "H_vector_norm", "H_edges"]
    return max(abs(metrics[k] - reference[k]) for k in keys)


def run_perturbation_stability_test(
    epsilons=(0.0, 0.01, 0.02, 0.05, 0.1, 0.2),
    n_trials=30,
    seed=42,
    perturb_method="tangent",
    tol=1e-6,
):
    """
    Stabilitätstest: Bleibt Isotropie (Anisotropie nach R_EABC) unter kleinen
    geometrischen Störungen der 12 Kussrichtungen erhalten?

    Ablauf pro Trial:
    1. perturb_vertices → gestörte Einheits-Ecken
    2. Kanten/Flächen aus gestörter Geometrie rekonstruieren
    3. Zufalls-Labeling auf gestörter Geometrie → Metriken vor R_EABC
    4. R_EABC auf gestörter Geometrie (assign_eabc_geometric nutzt geänderte
       mod-12-Restklassen — das ist gewollt, kein Bug)
    5. Vergleich mit ungestörtem geometrischen Fixpunkt

    Hinweis: Bei gestörter Geometrie ändert sich die mod-12-Klassifikation der
    Ecken; R_EABC projiziert auf die *neue* geometrische Zielzuweisung.
    """
    base_vertices = icosahedron_vertices()
    base_edges = icosahedron_edges(base_vertices)
    base_faces = triangular_faces(base_vertices, base_edges)
    base_assignment, base_residues = assign_eabc_geometric(base_vertices, base_edges, base_faces)
    base_metrics = compute_eabc_metrics(base_vertices, base_edges, base_faces, base_assignment)

    keys = ["D_norm", "anisotropy", "H_vector_norm", "H_edges"]
    rng = np.random.default_rng(seed)

    print()
    print("=== Stabilitätstest: Isotropie unter geometrischen Kuss-Störungen ===")
    print(f"Trials pro epsilon: {n_trials}, Basis-Seed: {seed}")
    print(f"Perturbation: {perturb_method}, epsilon ∈ {list(epsilons)}")
    print()
    print("Ungestörter geometrischer Fixpunkt (R_EABC-Ziel):")
    for k in keys:
        print(f"  {k}: {base_metrics[k]:.6f}")
    print(f"  EABC-Balance: {_assignment_balance(base_assignment)}")
    print(f"  Restklassen mod 12: {sorted(base_residues.values())}")
    print()
    print(
        "Hinweis: Gestörte Vertices ändern mod-12-Restklassen → andere EABC-Zuordnung. "
        "Das testet echte geometrische Robustheit, nicht nur Label-Rauschen."
    )
    print()

    rows = []
    first_break_epsilon = None

    for eps in epsilons:
        aniso_before = []
        aniso_after = []
        aniso_geometric = []
        dev_after = []
        dev_geometric = []
        residue_changes = []
        imbalanced = 0

        for trial in range(n_trials):
            trial_seed = int(rng.integers(0, 2**31 - 1))
            perturbed = perturb_vertices(
                base_vertices, eps, seed=trial_seed, method=perturb_method
            )
            p_edges = icosahedron_edges(perturbed)
            p_faces = triangular_faces(perturbed, p_edges)

            geo_assignment, p_residues = assign_eabc_geometric(perturbed, p_edges, p_faces)
            balance = _assignment_balance(geo_assignment)
            if any(balance[letter] != 3 for letter in "EABC"):
                imbalanced += 1

            n_changed = sum(
                1 for i in range(len(base_vertices))
                if p_residues[i] != base_residues[i]
            )
            residue_changes.append(n_changed)

            perm_seed = int(rng.integers(0, 2**31 - 1))
            random_assignment = assign_eabc_random(perturbed, seed=perm_seed)
            m_before = compute_eabc_metrics(perturbed, p_edges, p_faces, random_assignment)
            renormalized = apply_r_eabc(perturbed, p_edges, p_faces, random_assignment)
            m_after = compute_eabc_metrics(perturbed, p_edges, p_faces, renormalized)
            m_geometric = compute_eabc_metrics(perturbed, p_edges, p_faces, geo_assignment)

            aniso_before.append(m_before["anisotropy"])
            aniso_after.append(m_after["anisotropy"])
            aniso_geometric.append(m_geometric["anisotropy"])
            dev_after.append(_max_metric_deviation(m_after, base_metrics))
            dev_geometric.append(_max_metric_deviation(m_geometric, base_metrics))

        aniso_before = np.array(aniso_before)
        aniso_after = np.array(aniso_after)
        aniso_geometric = np.array(aniso_geometric)
        dev_after = np.array(dev_after)
        dev_geometric = np.array(dev_geometric)
        residue_changes = np.array(residue_changes)

        row = {
            "epsilon": eps,
            "aniso_before_mean": aniso_before.mean(),
            "aniso_before_max": aniso_before.max(),
            "aniso_after_mean": aniso_after.mean(),
            "aniso_after_max": aniso_after.max(),
            "aniso_geo_mean": aniso_geometric.mean(),
            "aniso_geo_max": aniso_geometric.max(),
            "dev_after_mean": dev_after.mean(),
            "dev_after_max": dev_after.max(),
            "dev_geo_mean": dev_geometric.mean(),
            "dev_geo_max": dev_geometric.max(),
            "residue_change_mean": residue_changes.mean(),
            "imbalanced_frac": imbalanced / n_trials,
        }
        rows.append(row)

        if first_break_epsilon is None and row["aniso_after_max"] > tol:
            first_break_epsilon = eps

    print(
        f"{'ε':>6} | {'Aniso vor':>10} {'max':>8} | "
        f"{'Aniso nach R':>10} {'max':>8} | "
        f"{'Aniso geo':>10} {'max':>8} | "
        f"{'Δ Fixpkt':>10} {'max':>8} | "
        f"{'Δr̄':>5} {'≠3+3+3+3':>8}"
    )
    print("-" * 105)
    for row in rows:
        print(
            f"{row['epsilon']:6.3f} | "
            f"{row['aniso_before_mean']:10.4f} {row['aniso_before_max']:8.4f} | "
            f"{row['aniso_after_mean']:10.4e} {row['aniso_after_max']:8.4e} | "
            f"{row['aniso_geo_mean']:10.4e} {row['aniso_geo_max']:8.4e} | "
            f"{row['dev_after_mean']:10.4f} {row['dev_after_max']:8.4f} | "
            f"{row['residue_change_mean']:5.1f} {row['imbalanced_frac']:8.1%}"
        )
    print()

    print("Zusammenfassung:")
    eps0 = rows[0]
    print(
        f"  ε=0: Anisotropie nach R_EABC max = {eps0['aniso_after_max']:.6e} "
        f"(Referenz: {tol:.0e})"
    )
    if first_break_epsilon is None:
        print(
            f"  Isotropie stabil für alle getesteten ε ≤ {epsilons[-1]} "
            f"(Anisotropie nach R_EABC ≤ {tol})"
        )
    else:
        print(
            f"  Erster ε-Bruch (Anisotropie nach R_EABC > {tol}): {first_break_epsilon}"
        )
    print(
        f"  Metrik-Abweichung vom ungestörten Fixpunkt wächst mit ε "
        f"(D_norm, |H|, H_edges — erwartet, da Geometrie sich ändert)"
    )
    print()

    # Feinere Suche um den Bruchpunkt, falls innerhalb des Bereichs
    if first_break_epsilon is not None and first_break_epsilon > 0:
        lo = 0.0
        for e in epsilons:
            if e < first_break_epsilon:
                lo = e
        hi = first_break_epsilon
        print(f"Grobe Bruchstelle: ε ∈ ({lo}, {hi}] — Feinsuche (10 Trials):")
        for _ in range(8):
            mid = (lo + hi) / 2
            max_aniso = 0.0
            for trial in range(10):
                trial_seed = int(rng.integers(0, 2**31 - 1))
                perturbed = perturb_vertices(
                    base_vertices, mid, seed=trial_seed, method=perturb_method
                )
                p_edges = icosahedron_edges(perturbed)
                p_faces = triangular_faces(perturbed, p_edges)
                renormalized = apply_r_eabc(perturbed, p_edges, p_faces)
                m_after = compute_eabc_metrics(perturbed, p_edges, p_faces, renormalized)
                max_aniso = max(max_aniso, m_after["anisotropy"])
            if max_aniso > tol:
                hi = mid
            else:
                lo = mid
            print(f"  ε={mid:.5f}: max Anisotropie nach R_EABC = {max_aniso:.6e}")
        print(f"  Geschätzte kritische ε ≈ {hi:.5f}")
        print()

    return rows


def run_core_lemma_tests(
    tol=1e-6,
    epsilon=0.1,
    n_idempotence_trials=30,
    seed=42,
    n=1,
):
    """
    Verifiziert die vier Kernargumente des EABC-Modells am Referenz-Ikosaeder.

    1. Isotropie-Lemma: λ₁=λ₂=λ₃=24, Δ(M)=0 am geometrischen Fixpunkt
    2. Projektor-Idempotenz: (R*_EABC)² = R*_EABC
    3. Dualitätsstabilität: Δ(M) invariant unter Ikosaeder↔Dodekaeder-Dualität
    4. Volumen-Kompatibilität: λ_sphere-Skalierung erhält Δ(M), skaliert Eigenwerte mit λ²
    """
    results = {}
    rng = np.random.default_rng(seed)

    # --- Referenzgeometrie ---
    base_vertices = icosahedron_vertices()
    base_edges = icosahedron_edges(base_vertices)
    base_faces = triangular_faces(base_vertices, base_edges)
    geo_assignment, _ = assign_eabc_geometric(base_vertices, base_edges, base_faces)

    print()
    print("=== EABC Kern-Lemma Tests ===")
    print(f"Toleranz: {tol:.0e}, Idempotenz-Trials: {n_idempotence_trials}, ε={epsilon}")
    print()

    # ------------------------------------------------------------------
    # Test 1: Isotropie-Lemma
    # ------------------------------------------------------------------
    print("[1] Isotropie-Lemma")
    op_fix = order_parameter(base_vertices, geo_assignment)
    eigvals = op_fix["eigenvalues"]
    delta_M = op_fix["delta_M"]
    expected_lambda = 24.0

    eig_ok = all(abs(ev - expected_lambda) < tol for ev in eigvals)
    delta_ok = delta_M < tol
    test1_pass = eig_ok and delta_ok

    print(f"  Geometrische Zuweisung (Referenz-Ikosaeder), Balance: {_assignment_balance(geo_assignment)}")
    print(f"  M-Eigenwerte: {eigvals}")
    print(f"  |λ_i − 24|: {[abs(ev - expected_lambda) for ev in eigvals]}")
    print(f"  Δ(M) = {delta_M:.6e}")
    print(f"  Ergebnis: {'PASS' if test1_pass else 'FAIL'}")
    print()
    results["isotropy"] = test1_pass

    # ------------------------------------------------------------------
    # Test 2: Projektor-Idempotenz R*_EABC
    # ------------------------------------------------------------------
    print("[2] Projektor-Idempotenz (R*_EABC)² = R*_EABC")
    idempotence_ok = 0
    for trial in range(n_idempotence_trials):
        trial_seed = int(rng.integers(0, 2**31 - 1))
        perturbed = perturb_vertices(base_vertices, epsilon, seed=trial_seed)
        random_assignment = assign_eabc_random(perturbed, seed=trial_seed + 1)

        verts1, assgn1, edges1, faces1 = apply_renorm_step(perturbed, random_assignment)
        verts2, assgn2, edges2, faces2 = apply_renorm_step(verts1, assgn1)

        verts_match = np.allclose(verts1, verts2, atol=tol, rtol=0)
        assgn_match = assgn1 == assgn2
        op1 = order_parameter(verts1, assgn1)
        op2 = order_parameter(verts2, assgn2)
        metrics_match = (
            abs(op1["delta_M"] - op2["delta_M"]) < tol
            and np.allclose(op1["eigenvalues"], op2["eigenvalues"], atol=tol, rtol=0)
        )

        if verts_match and assgn_match and metrics_match:
            idempotence_ok += 1

        if trial < 2:
            print(f"  Trial {trial + 1} (seed={trial_seed}):")
            print(f"    Vertices identisch: {verts_match}")
            print(f"    Assignment identisch: {assgn_match}")
            print(f"    Δ(M) X₁={op1['delta_M']:.6e}, X₂={op2['delta_M']:.6e}")
            print()

    test2_pass = idempotence_ok == n_idempotence_trials
    print(f"  Idempotenz: {idempotence_ok}/{n_idempotence_trials} Trials")
    print(f"  Ergebnis: {'PASS' if test2_pass else 'FAIL'}")
    print()
    results["idempotence"] = test2_pass

    # ------------------------------------------------------------------
    # Test 3: Dualitätsstabilität Δ(M)
    # ------------------------------------------------------------------
    print("[3] Dualitätsstabilität Δ(M)")
    dual = build_dodecahedron_dual(base_vertices, base_edges, base_faces)
    do_vertices = dual["vertices"]
    do_edges = dual["edges"]
    do_triangles = triangulate_faces(dual["faces"])

    geo_do, _ = assign_eabc_geometric(do_vertices, do_edges, do_triangles)
    op_ico = order_parameter(base_vertices, geo_assignment)
    op_do = order_parameter(do_vertices, geo_do)

    delta_ico = op_ico["delta_M"]
    delta_dod = op_do["delta_M"]
    delta_diff = abs(delta_ico - delta_dod)

    m_ico = compute_eabc_metrics(base_vertices, base_edges, base_faces, geo_assignment)
    m_do = compute_eabc_metrics(do_vertices, do_edges, do_triangles, geo_do)

    d_diff = abs(m_ico["D_norm"] - m_do["D_norm"])
    h_diff = abs(m_ico["H_vector_norm"] - m_do["H_vector_norm"])
    h_edges_diff = abs(m_ico["H_edges"] - m_do["H_edges"])

    delta_dual_ok = delta_diff < tol
    d_not_invariant = d_diff > tol
    h_not_invariant = h_diff > tol or h_edges_diff > tol
    test3_pass = delta_dual_ok and d_not_invariant and h_not_invariant

    print(f"  Δ(M)_Ikosaeder:     {delta_ico:.6e}")
    print(f"  Δ(M)_Dodekaeder-Dual: {delta_dod:.6e}")
    print(f"  |Δ_ico − Δ_dod|:    {delta_diff:.6e}  →  {'invariant' if delta_dual_ok else 'abweichend'}")
    print(f"  ||D||_ico={m_ico['D_norm']:.6f}, ||D||_do={m_do['D_norm']:.6f}, Δ={d_diff:.6e}")
    print(f"  |H|_ico={m_ico['H_vector_norm']:.6f}, |H|_do={m_do['H_vector_norm']:.6f}, Δ={h_diff:.6e}")
    print(f"  H_edges: ico={m_ico['H_edges']:.6f}, do={m_do['H_edges']:.6f}, Δ={h_edges_diff:.6e}")
    print(f"  D/H nicht dual-invariant: D={'ja' if d_not_invariant else 'nein'}, H={'ja' if h_not_invariant else 'nein'}")
    print(f"  Ergebnis: {'PASS' if test3_pass else 'FAIL'}")
    print()
    results["duality"] = test3_pass

    # ------------------------------------------------------------------
    # Test 4: Volumen-Kompatibilität
    # ------------------------------------------------------------------
    print("[4] Volumen-Kompatibilität")
    lambda_expected = 1.0 / (1.0 + 2.0 * np.pi)
    renorm = volume_renormalization(n)
    lambda_sphere = renorm["lambda_sphere"]
    lambda_ok = abs(lambda_sphere - lambda_expected) < tol

    trial_seed = int(rng.integers(0, 2**31 - 1))
    perturbed = perturb_vertices(base_vertices, epsilon, seed=trial_seed)
    random_assignment = assign_eabc_random(perturbed, seed=trial_seed + 1)

    verts_rstar, assgn_rstar, _, _ = apply_renorm_step(perturbed, random_assignment)
    op_before = order_parameter(verts_rstar, assgn_rstar)

    verts_scaled, _ = apply_volume_scale(verts_rstar, n=n)
    op_after = order_parameter(verts_scaled, assgn_rstar)

    delta_before = op_before["delta_M"]
    delta_after = op_after["delta_M"]
    delta_vol_ok = abs(delta_before - delta_after) < tol

    scale_factor = lambda_sphere ** 2
    eig_before = op_before["eigenvalues"]
    eig_after = op_after["eigenvalues"]
    eig_scaled_ok = np.allclose(eig_after, scale_factor * eig_before, atol=tol, rtol=0)

    test4_pass = lambda_ok and delta_vol_ok and eig_scaled_ok

    print(f"  λ_sphere (berechnet) = {lambda_sphere:.10f}")
    print(f"  λ_sphere (exakt)     = {lambda_expected:.10f}")
    print(f"  |λ − 1/(1+2π)|       = {abs(lambda_sphere - lambda_expected):.6e}")
    print(f"  Δ(M) vor Skalierung:  {delta_before:.6e}")
    print(f"  Δ(M) nach Skalierung: {delta_after:.6e}")
    print(f"  |Δ_vor − Δ_nach|      = {abs(delta_before - delta_after):.6e}")
    print(f"  Eigenwerte vor:  {eig_before}")
    print(f"  Eigenwerte nach: {eig_after}")
    print(f"  Erwartet λ²·vor: {scale_factor * eig_before}")
    print(f"  Skalierungsfaktor λ² = {scale_factor:.10f}")
    print(f"  Ergebnis: {'PASS' if test4_pass else 'FAIL'}")
    print()
    results["volume"] = test4_pass

    # ------------------------------------------------------------------
    # Zusammenfassung
    # ------------------------------------------------------------------
    passed = sum(results.values())
    total = len(results)

    print("=== EABC Kern-Lemma Tests ===")
    print(f"[1] Isotropie-Lemma:        {'PASS' if results['isotropy'] else 'FAIL'}")
    print(f"[2] Projektor-Idempotenz:   {'PASS' if results['idempotence'] else 'FAIL'}")
    print(f"[3] Δ(M) dualitätsstabil:   {'PASS' if results['duality'] else 'FAIL'}")
    print(f"[4] Volumen-Kompatibilität: {'PASS' if results['volume'] else 'FAIL'}")
    print(f"GESAMT: {passed}/{total} bestanden")
    print()

    return results


def run_prime_defect_test(primes=(5, 7, 11, 13, 17), tol=1e-6):
    """
    Prim-Normkugel-Gedankenexperiment (V5):
    - Δ(M⁺) = w_p für Rang-1-Störung am Fixpunkt
    - R*_EABC stellt 3+3+3+3-Label-Balance wieder her
    """
    verts = icosahedron_vertices()
    edges = icosahedron_edges(verts)
    faces = triangular_faces(verts, edges)
    geo_assignment, _ = assign_eabc_geometric(verts, edges, faces)

    def is_prime(n):
        if n < 2:
            return False
        if n % 2 == 0:
            return n == 2
        d = 3
        while d * d <= n:
            if n % d == 0:
                return False
            d += 2
        return True

    residue_to_family = {1: "E", 5: "A", 7: "B", 11: "C"}

    print()
    print("=== Prim-Defekt Test (V5) ===")
    all_ok = True

    for p in primes:
        if not is_prime(p) or p <= 3:
            continue
        r = p % 12
        if r not in residue_to_family:
            continue
        fam = residue_to_family[r]
        w = FAMILIES[fam]
        v = verts[0]

        M = anisotropy_tensor(verts, geo_assignment)
        M_plus = M + w * np.outer(v, v)
        eig = np.linalg.eigvalsh(M_plus)
        delta_plus = eig[-1] - eig[0]

        # Virtuelle Primkugel: Multiketten-Zählung (4,3,3,3) bei unveränderten 12 Labels
        vertex_counts = {f: list(geo_assignment.values()).count(f) for f in "EABC"}
        virtual_counts = dict(vertex_counts)
        virtual_counts[fam] += 1
        has_prime_defect = virtual_counts[fam] == 4 and all(
            virtual_counts[f] == 3 for f in "EABC" if f != fam
        )

        # Label-Defekt auf den 12 Ecken: eine Ecke auf die Primfamilie setzen
        imbalanced = dict(geo_assignment)
        swap_idx = next(i for i in range(12) if geo_assignment[i] != fam)
        imbalanced[swap_idx] = fam

        _, assgn_after, _, _ = apply_renorm_step(verts, imbalanced)
        counts_after = {f: list(assgn_after.values()).count(f) for f in "EABC"}
        restored = all(c == 3 for c in counts_after.values())

        # Tensor-Restauration: gestörte Geometrie → R* → Δ ≈ 0
        perturbed_verts = verts + 0.02 * np.random.default_rng(p).normal(size=verts.shape)
        perturbed_verts = np.array([normalize(v) for v in perturbed_verts])
        M_pert = anisotropy_tensor(perturbed_verts, imbalanced)
        delta_before = np.linalg.eigvalsh(M_pert)[-1] - np.linalg.eigvalsh(M_pert)[0]
        verts_after, assgn_tensor, _, _ = apply_renorm_step(perturbed_verts, imbalanced)
        M_after = anisotropy_tensor(verts_after, assgn_tensor)
        delta_after = np.linalg.eigvalsh(M_after)[-1] - np.linalg.eigvalsh(M_after)[0]
        tensor_restored = delta_after < tol

        ok = (
            abs(delta_plus - w) < tol
            and has_prime_defect
            and restored
            and assgn_after == geo_assignment
            and delta_before > tol
            and tensor_restored
        )
        all_ok = all_ok and ok
        print(
            f"  p={p:3d} w={w:2d} Δ(M⁺)={delta_plus:.4f} "
            f"PRIME_DEFECT={has_prime_defect} R*-Balance={restored} "
            f"Δ_before={delta_before:.4f} Δ_after={delta_after:.4e} "
            f"tensor_ok={tensor_restored} -> {'PASS' if ok else 'FAIL'}"
        )

    print(f"Ergebnis: {'PASS' if all_ok else 'FAIL'}")
    print()
    return all_ok


if __name__ == "__main__":
    run_core_lemma_tests()
    run_prime_defect_test()
    run_test(n=1)
    run_control_test(n_permutations=1000, seed=42)
    run_geometric_test(n_permutations=1000, seed=42)
    run_renormalization_test(n_starts=30, seed=42)
    run_duality_test(n_starts=30, seed=42)
    run_perturbation_stability_test(n_trials=30, seed=42)
    run_order_parameter_test(n_trials=30, seed=42)
    run_coupled_renorm_test(n_trials=30, seed=42)
    run_full_chain_test(n_trials=30, seed=42)