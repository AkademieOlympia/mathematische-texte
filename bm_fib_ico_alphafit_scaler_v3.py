#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bm_fib_ico_alphafit_scaler.py

AlphaFit-skalierbarer Label-Optimizer für das BM Fibonacci-Ikosaeder.

NEU (v3):
- Primzahlen werden via bm_data_primes.load_primes geladen

Vorher (v2):
- --resume  : setzt start-index automatisch auf (letzter anchor_index + stride)
- --jobs N  : Multiprocessing über Anker-Fenster (parallel)
- batchweise Ausgabe (CSV/JSONL) -> streambar, RAM-schonend

Siehe README_BM_FIB_ICO_AlphaFitScaler.md
"""

import argparse, json, math, random, itertools, csv, time, os, sys
from typing import List, Tuple, Dict, Optional
from multiprocessing import Pool, cpu_count

# Data loaders (separate modules)
from bm_data_primes import load_primes

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
            tok=line.replace(","," ").split()[0]
            primes.append(int(tok))
    return primes

def write_csv_header_if_needed(path, fieldnames):
    exists=os.path.exists(path) and os.path.getsize(path)>0
    if not exists:
        with open(path,"w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()

def append_csv_rows(path, fieldnames, rows):
    with open(path,"a",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f, fieldnames=fieldnames)
        for r in rows:
            w.writerow(r)

def append_jsonl_records(path, records):
    if not path:
        return
    with open(path,"a",encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False)+"\n")

def tail_last_anchor_index(csv_path: str) -> Optional[int]:
    """
    Return last anchor_index in an existing CSV (or None).
    Fast-ish tail read without pandas.
    """
    if not (os.path.exists(csv_path) and os.path.getsize(csv_path)>0):
        return None
    with open(csv_path, "rb") as f:
        f.seek(0, os.SEEK_END)
        end = f.tell()
        # read last ~8KB
        size = min(8192, end)
        f.seek(end - size)
        data = f.read().decode("utf-8", errors="ignore")
    lines = [ln for ln in data.splitlines() if ln.strip()]
    if len(lines) <= 1:
        return None
    # last non-header line
    last = lines[-1]
    # If the tail chunk ends in the middle of a line, scan from end backwards to find a valid CSV line
    # We'll parse by splitting comma and taking first field.
    for ln in reversed(lines):
        if ln.lower().startswith("anchor_index"):
            continue
        try:
            first = ln.split(",", 1)[0].strip()
            return int(first)
        except Exception:
            continue
    return None

# ----------------- Multiprocessing globals -----------------
G_EDGES=None
G_W_RAD=None
G_W_DIH=None

def _init_worker(edges, w_rad, w_dih):
    global G_EDGES, G_W_RAD, G_W_DIH
    G_EDGES=edges
    G_W_RAD=w_rad
    G_W_DIH=w_dih

def _worker_task(args):
    """
    args: (anchor_index, window_primes, n_level, objective, iters, T0, Tend, seeds, include_labeling, beta)
    """
    (i, window, n_level, objective, iters, T0, Tend, seeds, include_labeling, beta, out_jsonl) = args
    t0=time.time()

    energy_fn = build_objective(G_EDGES, G_W_RAD, G_W_DIH, window, objective, alpha=1.0, beta=beta)

    best_perm=None
    best_energy=None
    for s in seeds:
        perm, e = anneal(12, energy_fn, iters, T0, Tend, seed=int(s))
        if best_energy is None or e < best_energy:
            best_energy=e
            best_perm=perm

    # reporting metrics: smooth (unweighted) + same-class count (on best labeling)
    logs=[math.log(p) for p in window]
    cls=[mod12_class(p) for p in window]

    def smooth_unweighted(perm):
        s=0.0
        for u,v in G_EDGES:
            s += abs(logs[perm[u]] - logs[perm[v]])
        return s

    def same_class(perm):
        c=0
        for u,v in G_EDGES:
            cu,cv = cls[perm[u]], cls[perm[v]]
            if cu in ("e","a","b","c") and cu==cv:
                c+=1
        return c

    smoothE = smooth_unweighted(best_perm)
    sameC = same_class(best_perm)
    runtime_ms=int((time.time()-t0)*1000)

    row={
        "anchor_index": i,
        "anchor_prime": window[0],
        "n_level": n_level,
        "objective": objective,
        "iters": iters,
        "T0": T0,
        "Tend": Tend,
        "seeds": ",".join(str(x) for x in seeds),
        "best_energy": best_energy,
        "same_class_edges": sameC,
        "smooth_energy": smoothE,
        "runtime_ms": runtime_ms
    }

    record=None
    if out_jsonl:
        record={**row,
                "window_primes": window,
                "perm_vertex_to_window_index": best_perm,
                "labeling_vertex_order": [window[best_perm[v]] for v in range(12)]}
    if include_labeling:
        labeling=[window[best_perm[v]] for v in range(12)]
        row["labeling_vertex_order"]=";".join(str(p) for p in labeling)

    return row, record

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
    ap.add_argument("--beta", type=float, default=5.0, help="beta for combo_radial")
    ap.add_argument("--out-csv", default="alphafit_run_results.csv")
    ap.add_argument("--out-jsonl", default="", help="optional jsonl stream output")
    ap.add_argument("--include-labeling", action="store_true", help="store labeling as semicolon string in CSV")
    ap.add_argument("--resume", action="store_true", help="auto-continue from existing CSV (last anchor_index + stride)")
    ap.add_argument("--jobs", type=int, default=1, help="parallel jobs (multiprocessing). 1 = no parallelism")
    ap.add_argument("--batch-size", type=int, default=200, help="how many anchors to compute before writing to disk")
    ap.add_argument("--progress-every", type=int, default=2000, help="print progress every N anchors")
    args=ap.parse_args()

    if args.window != 12:
        raise ValueError("This scaler currently targets the icosahedron (window=12).")
    if args.stride <= 0:
        raise ValueError("stride must be >= 1.")
    if args.jobs < 1:
        raise ValueError("jobs must be >= 1.")
    if args.batch_size < 1:
        raise ValueError("batch-size must be >= 1.")

    primes=load_primes(args.primes_file)
    N=len(primes)
    end = args.end_index if args.end_index != -1 else (N-args.window+1)
    end = min(end, N-args.window+1)

    # resume logic
    if args.resume:
        last = tail_last_anchor_index(args.out_csv)
        if last is not None:
            args.start_index = last + args.stride
            print(f"[resume] last_anchor_index={last} -> start_index={args.start_index}")

    if args.start_index < 0 or args.start_index >= end:
        raise ValueError("start-index out of range after resume/end constraints.")

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

    write_csv_header_if_needed(args.out_csv, fieldnames)

    # prepare anchor indices
    anchors=list(range(args.start_index, end, args.stride))
    total=len(anchors)
    print(f"[run] primes={N}, anchors={total}, n={args.n_level}, objective={args.objective}, jobs={args.jobs}")

    # worker pool (optional)
    pool=None
    if args.jobs > 1:
        jobs = min(args.jobs, cpu_count())
        pool = Pool(processes=jobs, initializer=_init_worker, initargs=(edges, w_rad, w_dih))
    else:
        _init_worker(edges, w_rad, w_dih)

    processed=0
    t_global=time.time()
    try:
        # batch over anchors to keep memory stable
        for b0 in range(0, total, args.batch_size):
            b_anchors = anchors[b0:b0+args.batch_size]
            tasks=[]
            for i in b_anchors:
                window=primes[i:i+args.window]
                tasks.append((i, window, args.n_level, args.objective, args.iters, args.T0, args.Tend,
                              args.seeds, args.include_labeling, args.beta, args.out_jsonl))

            if pool:
                results = pool.map(_worker_task, tasks, chunksize=max(1, len(tasks)//(args.jobs*4)))
            else:
                results = list(map(_worker_task, tasks))

            rows=[]
            records=[]
            for row, rec in results:
                rows.append(row)
                if rec is not None:
                    records.append(rec)

            # ensure stable order by anchor_index
            rows.sort(key=lambda r: r["anchor_index"])
            records.sort(key=lambda r: r["anchor_index"]) if records else None

            append_csv_rows(args.out_csv, fieldnames, rows)
            append_jsonl_records(args.out_jsonl, records)

            processed += len(rows)
            if processed % args.progress_every == 0 or processed == total:
                dt = time.time() - t_global
                rate = processed/dt if dt>0 else 0.0
                print(f"[progress] {processed}/{total} anchors ({rate:.2f}/s) -> {args.out_csv}")

    finally:
        if pool:
            pool.close()
            pool.join()

    dt=time.time()-t_global
    print(f"[done] processed={processed}, seconds={dt:.2f}, out={args.out_csv}")

if __name__=="__main__":
    main()
