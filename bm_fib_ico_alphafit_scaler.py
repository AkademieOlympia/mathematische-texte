#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bm_fib_ico_alphafit_scaler.py

Purpose
-------
Scale BM Fibonacci-Ikosaeder label-optimizations to large prime pools
(e.g. AlphaFit "850000" = first 850,000 primes or any prime list you provide).

Core idea
---------
The geometry (Ikosaeder graph + convex-hull faces + edge weights) depends only on the
level n (Fibonacci approximation of φ) and NOT on which primes you assign. So we:

1) precompute the icosahedron graph + robust hull faces for a chosen n
2) precompute per-edge weights (radial/dual coupling, dihedral coupling)
3) slide a window of 12 primes across a large list (stride configurable)
4) optimize a chosen objective per window using simulated annealing (fast on 12 nodes)
5) write compact results (CSV) + optional detailed stream (JSONL) for later analysis

Input primes file formats
-------------------------
- .txt / .csv: one integer prime per line (commas allowed; first integer per line used)
- .json: either [p1,p2,...] or {"primes":[...]}.

Outputs
-------
CSV rows include:
- anchor_index, anchor_prime (= primes[i]), objective, energies, same_class_edges, runtime_ms
- optional: labeling (12 primes in vertex order) as a semicolon-separated string

JSONL (optional) stores a full record per anchor (incl. perm + labeling list).

Example
-------
python bm_fib_ico_alphafit_scaler.py \
  --primes-file example_primes_200.txt \
  --n-level 34 \
  --objective radial_smooth \
  --start-index 0 --end-index 50 --stride 1 \
  --iters 120000 --T0 1.0 --Tend 1e-5 \
  --seeds 1 2 3 \
  --out-csv alphafit_run_results.csv \
  --out-jsonl alphafit_run_results.jsonl

Notes on scaling to 850000
--------------------------
- This script is designed to stream results: it does NOT keep all labelings in RAM unless you want it.
- For 850k primes you likely want a stride > 1 (e.g. 10, 100) at first.
- You can resume by setting --start-index to the next anchor.
"""

import argparse, json, math, random, itertools, csv, time, os
from typing import List, Tuple, Dict

# ----------------- vector utils -----------------
def dot(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def sub(a,b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def norm(a): return math.sqrt(dot(a,a))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def unit(a, eps=1e-15):
    n = norm(a)
    if n < eps: return (0.0,0.0,0.0)
    return (a[0]/n,a[1]/n,a[2]/n)

def angle(u,v):
    nu, nv = norm(u), norm(v)
    if nu == 0 or nv == 0: return None
    c = dot(u,v)/(nu*nv)
    c = max(-1.0, min(1.0, c))
    return math.acos(c)

# ----------------- Fibonacci icosahedron -----------------
def fib_upto(n: int) -> List[int]:
    F = [0, 1]
    for _ in range(2, n+1):
        F.append(F[-1] + F[-2])
    return F

def fibonacci_icosahedron_vertices_unit(n: int) -> List[Tuple[float,float,float]]:
    """
    Build 12 vertices on a sphere using Fibonacci approximation.
    Uses (0, ±F_n, ±F_{n+1}), (±F_n, ±F_{n+1}, 0), (±F_{n+1}, 0, ±F_n)
    then normalizes to unit radius via sqrt(F_{2n+1}) identity.
    """
    F = fib_upto(2*n+2)
    Fn, Fnp1, F2n1 = F[n], F[n+1], F[2*n+1]  # F_n^2 + F_{n+1}^2 = F_{2n+1}
    s = 1.0 / math.sqrt(F2n1)
    verts = []
    for a, b in itertools.product([1, -1], repeat=2):
        verts.append((0.0, s*a*Fn, s*b*Fnp1))
    for a, b in itertools.product([1, -1], repeat=2):
        verts.append((s*a*Fn, s*b*Fnp1, 0.0))
    for a, b in itertools.product([1, -1], repeat=2):
        verts.append((s*a*Fnp1, 0.0, s*b*Fn))
    return verts

def dist2(p,q):
    return (p[0]-q[0])**2+(p[1]-q[1])**2+(p[2]-q[2])**2

def edges_by_min_distance(verts: List[Tuple[float,float,float]], tol=1e-10) -> List[Tuple[int,int]]:
    pairs=[]
    for i,j in itertools.combinations(range(len(verts)),2):
        pairs.append((dist2(verts[i],verts[j]), i, j))
    pairs.sort(key=lambda x:x[0])
    d2_min=pairs[0][0]
    thr=max(tol, tol*d2_min)
    edges=[(i,j) for (d2,i,j) in pairs if abs(d2-d2_min)<=thr]
    return edges

# ----------------- robust convex hull faces -----------------
def convex_hull_faces(verts, eps=1e-11):
    # oriented faces with outward normals (dot(n, center)>0).
    nV=len(verts)
    faces_set=set()
    oriented=[]
    for i,j,k in itertools.combinations(range(nV),3):
        A,B,C = verts[i], verts[j], verts[k]
        n = cross(sub(B,A), sub(C,A))
        nn = norm(n)
        if nn < eps:
            continue
        pos=neg=0
        for l in range(nV):
            if l in (i,j,k): 
                continue
            s = dot(n, sub(verts[l], A))
            if s > eps: pos += 1
            elif s < -eps: neg += 1
            if pos and neg:
                break
        if pos and neg:
            continue
        center=((A[0]+B[0]+C[0])/3.0, (A[1]+B[1]+C[1])/3.0, (A[2]+B[2]+C[2])/3.0)
        if dot(n, center) < 0:
            j,k = k,j
        tri_sorted=tuple(sorted((i,j,k)))
        if tri_sorted not in faces_set:
            faces_set.add(tri_sorted)
            oriented.append((i,j,k))
    return oriented

def edge_to_faces(faces):
    e2f={}
    for fi,(i,j,k) in enumerate(faces):
        for u,v in [(i,j),(j,k),(k,i)]:
            if u>v: u,v=v,u
            e2f.setdefault((u,v), []).append(fi)
    return e2f

def face_normals_and_radials(verts, faces):
    normals=[]; radials=[]
    for (i,j,k) in faces:
        A,B,C=verts[i],verts[j],verts[k]
        n=unit(cross(sub(B,A), sub(C,A)))
        center=((A[0]+B[0]+C[0])/3.0, (A[1]+B[1]+C[1])/3.0, (A[2]+B[2]+C[2])/3.0)
        r=unit(center)
        normals.append(n); radials.append(r)
    return normals, radials

def compute_edge_weights(verts, edges, faces, normals, radials):
    e2f=edge_to_faces(faces)
    w_rad={}
    w_dih={}
    for (u,v) in edges:
        key=(u,v) if u<v else (v,u)
        fids=e2f.get(key,[])
        if len(fids)!=2:
            continue
        f1,f2=fids
        d=unit(sub(verts[v], verts[u]))
        r1,r2=radials[f1], radials[f2]
        cr=(abs(dot(d,r1))+abs(dot(d,r2)))/2.0
        w_rad[(u,v)]=cr; w_rad[(v,u)]=cr
        dih=angle(normals[f1], normals[f2])
        w_dih[(u,v)]=dih; w_dih[(v,u)]=dih
    # normalize to mean 1
    if w_rad:
        m=sum(w_rad.values())/len(w_rad)
        if m!=0:
            for k in list(w_rad.keys()):
                w_rad[k]/=m
    if w_dih:
        m=sum(w_dih.values())/len(w_dih)
        if m!=0:
            for k in list(w_dih.keys()):
                w_dih[k]/=m
    return w_rad, w_dih

# ----------------- prime classes -----------------
def mod12_class(p:int)->str:
    r=p%12
    if r==1: return "e"
    if r==5: return "a"
    if r==7: return "b"
    if r==11:return "c"
    return "other"

# ----------------- annealing optimizer on 12 nodes -----------------
def anneal(V, energy_fn, iters, T0, Tend, seed):
    rnd=random.Random(seed)
    perm=list(range(V))
    rnd.shuffle(perm)
    e=energy_fn(perm)
    best=perm[:]; best_e=e
    for t in range(1, iters+1):
        frac=t/iters
        T=T0*((Tend/T0)**frac)
        i=rnd.randrange(V)
        j=rnd.randrange(V-1)
        if j>=i: j+=1
        perm[i],perm[j]=perm[j],perm[i]
        e2=energy_fn(perm)
        de=e2-e
        if de<=0 or rnd.random()<math.exp(-de/max(T,1e-12)):
            e=e2
            if e<best_e:
                best_e=e; best=perm[:]
        else:
            perm[i],perm[j]=perm[j],perm[i]
    return best, best_e

# ----------------- objectives -----------------
def build_objective(edges, w_rad, w_dih, primes_window, objective, alpha=1.0, beta=5.0):
    """
    primes_window: list of 12 primes (candidate labels). We permute assignment to vertices.
    objective:
      - "smooth": sum |ln p_u - ln p_v|
      - "radial_smooth": sum w_rad(u,v) * |ln p_u - ln p_v|
      - "dihedral_smooth": sum w_dih(u,v) * |ln p_u - ln p_v|
      - "combo_radial": radial_smooth + beta * same_class_edges
    """
    logs=[math.log(p) for p in primes_window]
    cls=[mod12_class(p) for p in primes_window]

    def smooth_term(perm, wmap=None):
        s=0.0
        for u,v in edges:
            w=1.0 if wmap is None else wmap.get((u,v), 1.0)
            s += w * abs(logs[perm[u]] - logs[perm[v]])
        return s

    def same_class_edges(perm):
        s=0.0
        for u,v in edges:
            cu,cv = cls[perm[u]], cls[perm[v]]
            if cu in ("e","a","b","c") and cu==cv:
                s += 1.0
        return s

    if objective=="smooth":
        return lambda perm: smooth_term(perm, None)
    if objective=="radial_smooth":
        return lambda perm: smooth_term(perm, w_rad)
    if objective=="dihedral_smooth":
        return lambda perm: smooth_term(perm, w_dih)
    if objective=="combo_radial":
        return lambda perm: alpha*smooth_term(perm, w_rad) + beta*same_class_edges(perm)
    raise ValueError(f"unknown objective: {objective}")

# ----------------- IO -----------------
def load_primes(path:str)->List[int]:
    path_lower=path.lower()
    if path_lower.endswith(".json"):
        with open(path,"r",encoding="utf-8") as f:
            obj=json.load(f)
        if isinstance(obj,list):
            return [int(x) for x in obj]
        if isinstance(obj,dict) and "primes" in obj:
            return [int(x) for x in obj["primes"]]
        raise ValueError("JSON must be a list [p1,p2,...] or {'primes':[...]}")

    primes=[]
    with open(path,"r",encoding="utf-8") as f:
        for line in f:
            line=line.strip()
            if not line:
                continue
            # allow CSV-like: take first integer token
            tok=line.replace(","," ").split()[0]
            primes.append(int(tok))
    return primes

def write_csv_header(path, fieldnames):
    exists=os.path.exists(path) and os.path.getsize(path)>0
    if not exists:
        with open(path,"w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()

def append_csv(path, fieldnames, row):
    with open(path,"a",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f, fieldnames=fieldnames)
        w.writerow(row)

def append_jsonl(path, record):
    with open(path,"a",encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False)+"\n")

# ----------------- main -----------------
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--primes-file", required=True, help="txt/csv/json prime list")
    ap.add_argument("--n-level", type=int, default=34, help="Fibonacci level n (e.g. 13,21,34)")
    ap.add_argument("--objective", default="radial_smooth",
                    choices=["smooth","radial_smooth","dihedral_smooth","combo_radial"])
    ap.add_argument("--window", type=int, default=12, help="must be 12 for icosahedron")
    ap.add_argument("--start-index", type=int, default=0)
    ap.add_argument("--end-index", type=int, default=-1, help="exclusive; -1 = max")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--iters", type=int, default=120000)
    ap.add_argument("--T0", type=float, default=1.0)
    ap.add_argument("--Tend", type=float, default=1e-5)
    ap.add_argument("--seeds", type=int, nargs="+", default=[1,2,3], help="run multiple seeds, keep best")
    ap.add_argument("--out-csv", default="alphafit_run_results.csv")
    ap.add_argument("--out-jsonl", default="", help="optional jsonl stream output")
    ap.add_argument("--include-labeling", action="store_true", help="store labeling as semicolon string in CSV")
    args=ap.parse_args()

    if args.window != 12:
        raise ValueError("This scaler currently targets the icosahedron (window=12).")

    primes=load_primes(args.primes_file)
    N=len(primes)
    end = args.end_index if args.end_index != -1 else (N-args.window+1)
    end = min(end, N-args.window+1)
    if args.start_index < 0 or args.start_index >= end:
        raise ValueError("start-index out of range.")
    if args.stride <= 0:
        raise ValueError("stride must be >= 1.")

    # precompute geometry for n-level on unit sphere
    verts=fibonacci_icosahedron_vertices_unit(args.n_level)
    edges=edges_by_min_distance(verts)
    faces=convex_hull_faces(verts)
    normals, radials = face_normals_and_radials(verts, faces)
    w_rad, w_dih = compute_edge_weights(verts, edges, faces, normals, radials)

    # output setup
    fieldnames=[
        "anchor_index","anchor_prime","n_level","objective","iters","T0","Tend","seeds",
        "best_energy","same_class_edges","smooth_energy","runtime_ms"
    ]
    if args.include_labeling:
        fieldnames.append("labeling_vertex_order")

    write_csv_header(args.out_csv, fieldnames)

    # run
    for i in range(args.start_index, end, args.stride):
        window=primes[i:i+args.window]
        # build objective once per window
        energy_fn = build_objective(edges, w_rad, w_dih, window, args.objective)

        # track best across seeds
        best_perm=None
        best_energy=None
        best_ms=None
        t0=time.time()
        for s in args.seeds:
            perm, e = anneal(12, energy_fn, args.iters, args.T0, args.Tend, seed=s)
            if best_energy is None or e < best_energy:
                best_energy=e
                best_perm=perm
        runtime_ms=int((time.time()-t0)*1000)

        # compute reporting metrics
        # smooth (unweighted) + same-class count (on best labeling)
        logs=[math.log(p) for p in window]
        cls=[mod12_class(p) for p in window]

        def smooth_unweighted(perm):
            s=0.0
            for u,v in edges:
                s += abs(logs[perm[u]] - logs[perm[v]])
            return s

        def same_class(perm):
            c=0
            for u,v in edges:
                cu,cv = cls[perm[u]], cls[perm[v]]
                if cu in ("e","a","b","c") and cu==cv:
                    c+=1
            return c

        smoothE = smooth_unweighted(best_perm)
        sameC = same_class(best_perm)

        row={
            "anchor_index": i,
            "anchor_prime": window[0],
            "n_level": args.n_level,
            "objective": args.objective,
            "iters": args.iters,
            "T0": args.T0,
            "Tend": args.Tend,
            "seeds": ",".join(str(x) for x in args.seeds),
            "best_energy": best_energy,
            "same_class_edges": sameC,
            "smooth_energy": smoothE,
            "runtime_ms": runtime_ms
        }

        if args.include_labeling:
            labeling=[window[best_perm[v]] for v in range(12)]
            row["labeling_vertex_order"]=";".join(str(p) for p in labeling)

        append_csv(args.out_csv, fieldnames, row)

        if args.out_jsonl:
            record={
                **row,
                "window_primes": window,
                "perm_vertex_to_window_index": best_perm,
                "labeling_vertex_order": [window[best_perm[v]] for v in range(12)]
            }
            append_jsonl(args.out_jsonl, record)

        # minimal progress indicator without flooding
        if (i-args.start_index) % (args.stride*50) == 0:
            print(f"[{i}/{end}) best_energy={best_energy:.6g} same={sameC} runtime_ms={runtime_ms}")

if __name__=="__main__":
    main()
