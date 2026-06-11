# bm_fib_ico_generator.py
# Generates a Fibonacci-discretized icosahedron shell (and its dodecahedral dual)
# in the style discussed for the Bamberg-Modell.
#
# Output: BM_FIB_ICO_<n>.json

import json, math, itertools, datetime, os

def fib_upto(n: int):
    F = [0, 1]
    for _ in range(2, n+1):
        F.append(F[-1] + F[-2])
    return F

def fibonacci_icosahedron_vertices(n: int, R: float):
    F = fib_upto(2*n+2)
    Fn, Fnp1, F2n1 = F[n], F[n+1], F[2*n+1]
    s = R / math.sqrt(F2n1)
    verts = []
    for a, b in itertools.product([1, -1], repeat=2):
        verts.append((0.0, s*a*Fn, s*b*Fnp1))
    for a, b in itertools.product([1, -1], repeat=2):
        verts.append((s*a*Fn, s*b*Fnp1, 0.0))
    for a, b in itertools.product([1, -1], repeat=2):
        verts.append((s*a*Fnp1, 0.0, s*b*Fn))
    return verts, {"Fn": Fn, "Fnp1": Fnp1, "F2n1": F2n1, "scale": s}

def dist2(p, q):
    return (p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2

def dot(a,b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def norm(a):
    return math.sqrt(dot(a,a))

def sub(a,b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def unit(a, eps=1e-15):
    na = norm(a)
    if na < eps:
        return (0.0,0.0,0.0)
    return (a[0]/na, a[1]/na, a[2]/na)

def edges_by_min_distance(verts, tol=1e-10):
    pairs = []
    for i, j in itertools.combinations(range(len(verts)), 2):
        d2 = dist2(verts[i], verts[j])
        pairs.append((d2, i, j))
    pairs.sort(key=lambda x: x[0])
    d2_min = pairs[0][0]
    thr = max(tol, tol*d2_min)
    edges = [(i, j) for (d2, i, j) in pairs if abs(d2 - d2_min) <= thr]
    return d2_min, edges

def faces_from_edges(num_verts, edges):
    adj = [set() for _ in range(num_verts)]
    for i,j in edges:
        adj[i].add(j); adj[j].add(i)
    faces = []
    for i in range(num_verts):
        for j in adj[i]:
            if j <= i:
                continue
            common = adj[i].intersection(adj[j])
            for k in common:
                if k <= j:
                    continue
                faces.append((i,j,k))
    return faces, adj

def face_center(verts, tri):
    a,b,c = (verts[tri[0]], verts[tri[1]], verts[tri[2]])
    return ((a[0]+b[0]+c[0])/3.0, (a[1]+b[1]+c[1])/3.0, (a[2]+b[2]+c[2])/3.0)

def normalize_to_radius(p, R):
    u = unit(p)
    return (R*u[0], R*u[1], R*u[2])

def angle_at_vertex(u, v, w):
    a = sub(u, v)
    b = sub(w, v)
    na, nb = norm(a), norm(b)
    if na == 0 or nb == 0:
        return None
    c = dot(a,b)/(na*nb)
    c = max(-1.0, min(1.0, c))
    return math.acos(c)

def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    r = int(math.isqrt(n))
    f = 3
    while f <= r:
        if n % f == 0:
            return False
        f += 2
    return True

def primes_from(start, count):
    out = []
    x = start
    while len(out) < count:
        if is_prime(x):
            out.append(x)
        x += 1
    return out

def mod12_class(p):
    r = p % 12
    if r == 1: return "e(1)"
    if r == 5: return "a(5)"
    if r == 7: return "b(7)"
    if r == 11: return "c(11)"
    return f"other({r})"

def build_bm_fib_ico(n=21, p_anchor=7919):
    R = math.log(p_anchor)
    verts, meta = fibonacci_icosahedron_vertices(n, R)
    d2_min, edges = edges_by_min_distance(verts)
    faces, adj = faces_from_edges(len(verts), edges)
    dual = [normalize_to_radius(face_center(verts, f), R) for f in faces]

    radii = [norm(v) for v in verts]
    edge_lengths = [math.sqrt(dist2(verts[i], verts[j])) for i,j in edges]

    vertex_angles = {}
    for v in range(len(verts)):
        neigh = sorted(list(adj[v]))
        angs = []
        for u,w in itertools.combinations(neigh, 2):
            ang = angle_at_vertex(verts[u], verts[v], verts[w])
            if ang is not None:
                angs.append(ang)
        vertex_angles[str(v)] = {
            "degree": len(neigh),
            "angles_rad": angs,
            "angles_rad_mean": sum(angs)/len(angs) if angs else None
        }

    labels = primes_from(p_anchor, 12)
    label_info = []
    for idx, p in enumerate(labels):
        label_info.append({
            "vertex_index": idx,
            "prime": p,
            "prime_log": math.log(p),
            "mod12": p % 12,
            "mod12_class": mod12_class(p)
        })

    obj = {
        "type": "BM_FIB_ICO",
        "version": "v1",
        "created_utc": datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z",
        "params": {
            "level_n": n,
            "p_anchor": p_anchor,
            "R": R,
            "phi_approx": f"F_{n+1}/F_{n}",
            "Fn": meta["Fn"],
            "Fnp1": meta["Fnp1"],
            "F2n1": meta["F2n1"],
            "scale_s": meta["scale"]
        },
        "counts": {"V": len(verts), "E": len(edges), "F": len(faces), "dual_V": len(dual)},
        "vertices": [{"i": i, "x": v[0], "y": v[1], "z": v[2], "r": radii[i]} for i,v in enumerate(verts)],
        "edges": [{"i": i, "u": u, "v": v, "length": math.sqrt(dist2(verts[u], verts[v]))} for i,(u,v) in enumerate(edges)],
        "faces": [{"i": i, "v": list(f)} for i,f in enumerate(faces)],
        "dual_vertices": [{"i": i, "x": d[0], "y": d[1], "z": d[2], "r": norm(d)} for i,d in enumerate(dual)],
        "stats": {
            "radii_min": min(radii),
            "radii_max": max(radii),
            "edge_length_min": min(edge_lengths),
            "edge_length_max": max(edge_lengths),
            "edge_length_mean": sum(edge_lengths)/len(edge_lengths),
            "d2_min": d2_min,
        },
        "vertex_angle_stats": vertex_angles,
        "labels_default": {
            "scheme": "BM1000_like_sequential_primes_by_index",
            "labels": label_info
        }
    }
    return obj

if __name__ == "__main__":
    n = int(os.environ.get("BM_FIB_N", "21"))
    p_anchor = int(os.environ.get("BM_FIB_P", "7919"))
    obj = build_bm_fib_ico(n=n, p_anchor=p_anchor)
    out = f"BM_FIB_ICO_{n}.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
    print("wrote", out)
