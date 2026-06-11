
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import math
import os
import random
import time
import heapq
from dataclasses import dataclass
from typing import Tuple, Iterator, Optional, Callable, List, Iterable, Dict

# ============================================================
# 0) Utilities
# ============================================================

def isqrt(n: int) -> int:
    return int(math.isqrt(n))

def gcd(a: int, b: int) -> int:
    return math.gcd(a, b)

def round_nearest_rational(num: int, den: int) -> int:
    # Round num/den to nearest integer, ties away from zero
    assert den > 0
    if num >= 0:
        q, r = divmod(num, den)
        if 2*r > den:  return q + 1
        if 2*r < den:  return q
        return q + 1
    else:
        return -round_nearest_rational(-num, den)

# ============================================================
# 1) Jacobi + signature class (A/B/E/C-like)
# ============================================================

def jacobi(a: int, n: int) -> int:
    a %= n
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be positive odd")
    result = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            r = n % 8
            if r in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0

def v4_signature(n: int) -> Tuple[int,int]:
    # (chi4, chi3) = ((-1/n), (-3/n))
    return (jacobi(-1, n), jacobi(-3, n))

def abc_like_class(n: int) -> str:
    # Convention:
    # C-like: (-1,-1), E-like: (1,1), (1,-1)->A, (-1,1)->B
    chi4, chi3 = v4_signature(n)
    if (chi4, chi3) == (-1, -1): return "C"
    if (chi4, chi3) == ( 1,  1): return "E"
    if (chi4, chi3) == ( 1, -1): return "A"
    if (chi4, chi3) == (-1,  1): return "B"
    return "U"

# ============================================================
# 2) Directions + caching
# ============================================================

_DIRS_CACHE: Dict[Tuple[int,int], List[Tuple[int,int]]] = {}
_DIRS_S1_CACHE: Dict[Tuple[str,int,int], List[Tuple[float,float]]] = {}

def surrogate_dirs(m: int, seed: int = 0) -> List[Tuple[int,int]]:
    # Fast integer directions, reduced by gcd, dedup up to sign
    rng = random.Random(seed)
    dirs: List[Tuple[int,int]] = []
    seen = set()
    while len(dirs) < m:
        x = rng.randint(-200, 200)
        y = rng.randint(-200, 200)
        if x == 0 and y == 0:
            continue
        g = math.gcd(abs(x), abs(y))
        x //= g; y //= g
        k = (x, y)
        kn = (-x, -y)
        if k in seen or kn in seen:
            continue
        seen.add(k)
        dirs.append(k)
    return dirs

def normalize2(x: float, y: float) -> Tuple[float,float]:
    r = math.hypot(x, y)
    if r == 0.0:
        return (0.0, 0.0)
    return (x/r, y/r)

SQRT3 = 1.7320508075688772
SQRT3_OVER_2 = 0.8660254037844386

def dirvec_gaussian(dx: int, dy: int) -> Tuple[float,float]:
    return normalize2(float(dx), float(dy))

def dirvec_eisenstein(da: int, db: int) -> Tuple[float,float]:
    x = float(da) - 0.5*float(db)
    y = SQRT3_OVER_2 * float(db)
    return normalize2(x, y)

def get_surrogate_dirs_cached(m: int, seed: int) -> List[Tuple[int,int]]:
    key = (m, seed)
    hit = _DIRS_CACHE.get(key)
    if hit is not None:
        return hit
    dirs = surrogate_dirs(m, seed=seed)
    _DIRS_CACHE[key] = dirs
    return dirs

def get_dirs_s1_cached(ring: str, m: int, seed: int) -> List[Tuple[float,float]]:
    key = (ring, m, seed)
    hit = _DIRS_S1_CACHE.get(key)
    if hit is not None:
        return hit
    raw = get_surrogate_dirs_cached(m, seed)
    if ring == "gaussian":
        s1 = [dirvec_gaussian(dx, dy) for dx, dy in raw]
    else:
        s1 = [dirvec_eisenstein(da, db) for da, db in raw]
    _DIRS_S1_CACHE[key] = s1
    return s1

# ============================================================
# 3) Phase / cheap scoring
# ============================================================

def phi_class_offset(cls: str) -> float:
    if cls == "E": return 0.0
    if cls == "A": return math.pi/2
    if cls == "B": return math.pi
    if cls == "C": return 3*math.pi/2
    return 0.0

def phi_of_n(n: int) -> float:
    cls = abc_like_class(n)
    base = phi_class_offset(cls)
    M = 2147483647
    x = math.sqrt(n) / M
    frac = x - math.floor(x)
    return base + 2*math.pi*frac

def target_vec(n: int) -> Tuple[float,float]:
    phi = phi_of_n(n)
    return (math.cos(phi), math.sin(phi))

def alignment_topk(dirs_s1: List[Tuple[float,float]], t: Tuple[float,float], topk: int = 8) -> float:
    dots = [rx*t[0] + ry*t[1] for rx,ry in dirs_s1 if (rx != 0.0 or ry != 0.0)]
    if not dots:
        return -1.0
    dots.sort(reverse=True)
    k = min(topk, len(dots))
    return sum(dots[:k]) / k

def ring_alignment(n: int, m: int, seed: int, ring: str, topk: int = 8) -> float:
    t = target_vec(n)
    dirs_s1 = get_dirs_s1_cached(ring, m, seed)
    return alignment_topk(dirs_s1, t, topk=topk)

def split_budget(B: int, Ag: float, Ae: float) -> Tuple[int,int]:
    Ag = max(0.0, Ag)
    Ae = max(0.0, Ae)
    s = Ag + Ae
    if s <= 1e-12:
        Bg = B//2
        return (Bg, B-Bg)
    Bg = int(round(B * (Ag / s)))
    Bg = max(0, min(B, Bg))
    return (Bg, B-Bg)

# Two cheap scorers. Both return "higher is better".
def score_phase_gaussian(g: Tuple[int,int], n: int) -> float:
    M = 97
    cls = abc_like_class(n)
    targets = {"E": (0,0), "A": (M//3, M//3), "B": (2*M//3, M//3), "C": (M//2, 2*M//3)}
    tx, ty = targets.get(cls, (0,0))
    x = g[0] % M
    y = g[1] % M
    dx = min((x-tx) % M, (tx-x) % M)
    dy = min((y-ty) % M, (ty-y) % M)
    return -(dx*dx + dy*dy)

def score_phase_eisenstein(g: Tuple[int,int], n: int) -> float:
    M = 97
    cls = abc_like_class(n)
    targets = {"E": (0,0), "A": (M//3, M//3), "B": (2*M//3, M//3), "C": (M//2, 2*M//3)}
    ta, tb = targets.get(cls, (0,0))
    a = g[0] % M
    b = g[1] % M
    da = min((a-ta) % M, (ta-a) % M)
    db = min((b-tb) % M, (tb-b) % M)
    return -(da*da + db*db)

def score_l1_gaussian(g: Tuple[int,int], n: int) -> float:
    base = isqrt(n)
    return -(abs(g[0]-base) + abs(g[1]))

def score_l1_eisenstein(g: Tuple[int,int], n: int) -> float:
    base = isqrt(n)
    return -(abs(g[0]-base) + abs(g[1]))

# ============================================================
# 4) Tuple-based Z[i]
# ============================================================

GI = Tuple[int,int]

def gi_mul(x: GI, y: GI) -> GI:
    a,b = x
    c,d = y
    return (a*c - b*d, a*d + b*c)

def gi_conj(x: GI) -> GI:
    return (x[0], -x[1])

def gi_sub(x: GI, y: GI) -> GI:
    return (x[0] - y[0], x[1] - y[1])

def gi_norm(x: GI) -> int:
    a,b = x
    return a*a + b*b

def gi_is_zero(x: GI) -> bool:
    return x[0] == 0 and x[1] == 0

def gi_divmod(alpha: GI, beta: GI) -> Tuple[GI, GI]:
    if gi_is_zero(beta):
        raise ZeroDivisionError("division by zero in Z[i]")
    N = gi_norm(beta)
    z = gi_mul(alpha, gi_conj(beta))
    qa = round_nearest_rational(z[0], N)
    qb = round_nearest_rational(z[1], N)
    q = (qa, qb)
    r = gi_sub(alpha, gi_mul(q, beta))
    return q, r

def gi_gcd(x: GI, y: GI) -> GI:
    a, b = x, y
    while not gi_is_zero(b):
        _, r = gi_divmod(a, b)
        a, b = b, r
    return a

# ============================================================
# 5) Tuple-based Z[ω]
# ============================================================

EI = Tuple[int,int]

def ei_mul(x: EI, y: EI) -> EI:
    a,b = x
    c,d = y
    return (a*c - b*d, a*d + b*c - b*d)

def ei_conj(x: EI) -> EI:
    a,b = x
    return (a - b, -b)

def ei_sub(x: EI, y: EI) -> EI:
    return (x[0] - y[0], x[1] - y[1])

def ei_norm(x: EI) -> int:
    a,b = x
    return a*a - a*b + b*b

def ei_is_zero(x: EI) -> bool:
    return x[0] == 0 and x[1] == 0

def ei_to_complex(x: EI) -> complex:
    a,b = x
    return complex(a - b/2.0, SQRT3_OVER_2*b)

def nearest_eisenstein_from_complex(z: complex) -> EI:
    # cube-coord rounding
    x, y = z.real, z.imag
    v = (2.0 / SQRT3) * y
    u = x + 0.5*v
    w = -u - v

    ru, rv, rw = round(u), round(v), round(w)
    du, dv, dw = abs(ru - u), abs(rv - v), abs(rw - w)

    if du > dv and du > dw:
        ru = -rv - rw
    elif dv > dw:
        rv = -ru - rw
    else:
        rw = -ru - rv

    return (int(ru), int(rv))

class EiStats:
    __slots__ = ("div_calls", "neighbor_calls")
    def __init__(self):
        self.div_calls = 0
        self.neighbor_calls = 0

EI_STATS = EiStats()

def ei_divmod_adaptive(alpha: EI, beta: EI,
                       neighbor_safety: bool = True,
                       adaptive: bool = True,
                       trigger_ratio: float = 0.92) -> Tuple[EI, EI]:
    EI_STATS.div_calls += 1

    if ei_is_zero(beta):
        raise ZeroDivisionError("division by zero in Z[ω]")

    Nbeta = ei_norm(beta)
    t = ei_mul(alpha, ei_conj(beta))
    z = ei_to_complex(t) / float(Nbeta)

    q0 = nearest_eisenstein_from_complex(z)
    r0 = ei_sub(alpha, ei_mul(q0, beta))

    if not neighbor_safety:
        return q0, r0

    if adaptive and ei_norm(r0) <= int(trigger_ratio * Nbeta):
        return q0, r0

    EI_STATS.neighbor_calls += 1

    neighbors = [(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)]
    qa0, qb0 = q0

    best_q = q0
    best_r = r0
    best_n = ei_norm(r0)

    for da,db in neighbors:
        q = (qa0 + da, qb0 + db)
        r = ei_sub(alpha, ei_mul(q, beta))
        nr = ei_norm(r)
        if nr < best_n:
            best_n = nr
            best_q = q
            best_r = r

    return best_q, best_r

def ei_gcd(x: EI, y: EI, neighbor_safety: bool = True) -> EI:
    a, b = x, y
    while not ei_is_zero(b):
        _, r = ei_divmod_adaptive(a, b, neighbor_safety=neighbor_safety, adaptive=True, trigger_ratio=0.92)
        a, b = b, r
    return a

# ============================================================
# 6) Candidate streaming + Top-B heap
# ============================================================

@dataclass(frozen=True)
class CandidateSpec:
    T: int = 40
    S: int = 40
    base_samples: int = 8
    dir_count: int = 120
    lambdas: Tuple[int, ...] = (1,2,3,5,8)
    dir_seed_offset: int = 12345

def candidates_gaussian_stream(n: int, spec: CandidateSpec, seed: int) -> Iterator[GI]:
    rng = random.Random(seed)
    base = isqrt(n)
    pairs = [(rng.randint(-spec.T, spec.T), rng.randint(-spec.S, spec.S))
             for _ in range(spec.base_samples)]
    dirs = get_surrogate_dirs_cached(spec.dir_count, seed + spec.dir_seed_offset)

    for t,s in pairs:
        g0 = (base + t, s)
        yield g0
        a0, b0 = g0
        for dx,dy in dirs:
            for lam in spec.lambdas:
                yield (a0 + lam*dx, b0 + lam*dy)

def candidates_eisenstein_stream(n: int, spec: CandidateSpec, seed: int) -> Iterator[EI]:
    rng = random.Random(seed)
    base = isqrt(n)
    pairs = [(rng.randint(-spec.T, spec.T), rng.randint(-spec.S, spec.S))
             for _ in range(spec.base_samples)]
    dirs = get_surrogate_dirs_cached(spec.dir_count, seed + spec.dir_seed_offset + 777)

    for t,s in pairs:
        g0 = (base + t, s)
        yield g0
        a0, b0 = g0
        for da,db in dirs:
            for lam in spec.lambdas:
                yield (a0 + lam*da, b0 + lam*db)

def top_b_candidates_stream(candidates: Iterable[Tuple[int,int]],
                            score_fn: Callable[[Tuple[int,int]], float],
                            B: int) -> List[Tuple[int,int]]:
    if B <= 0:
        return []
    heap: List[Tuple[float, int, Tuple[int,int]]] = []
    counter = 0
    for cand in candidates:
        s = score_fn(cand)
        item = (-s, counter, cand)  # store -score => worst on top
        counter += 1
        if len(heap) < B:
            heapq.heappush(heap, item)
        else:
            if item > heap[0]:
                heapq.heapreplace(heap, item)
    return [x[2] for x in sorted(heap, reverse=True)]

# ============================================================
# 7) Prefilter + Shadow
# ============================================================

@dataclass
class PrefilterResult:
    success: bool
    factor: Optional[int]
    rank_hit: Optional[int]
    tested: int
    t_prefilter: float
    ring_used: str

def try_factor_gaussian(n: int, g: GI) -> Optional[int]:
    d = gi_gcd((n, 0), g)
    f = gcd(gi_norm(d), n)
    return f if 1 < f < n else None

def try_factor_eisenstein(n: int, g: EI, neighbor_safety: bool = True) -> Optional[int]:
    d = ei_gcd((n, 0), g, neighbor_safety=neighbor_safety)
    f = gcd(ei_norm(d), n)
    return f if 1 < f < n else None

def prefilter_gaussian_fast(n: int, spec: CandidateSpec, B: int, score_mode: str, seed: int) -> PrefilterResult:
    t0 = time.time()
    score_fn_local = (lambda g: score_phase_gaussian(g, n)) if score_mode == "phase" else (lambda g: score_l1_gaussian(g, n))
    ordered = top_b_candidates_stream(candidates_gaussian_stream(n, spec, seed=seed), score_fn_local, B)

    tested = 0
    for idx, g in enumerate(ordered, start=1):
        tested += 1
        f = try_factor_gaussian(n, g)
        if f is not None:
            return PrefilterResult(True, f, idx, tested, time.time()-t0, "gaussian")
    return PrefilterResult(False, None, None, tested, time.time()-t0, "gaussian")

def prefilter_eisenstein_fast(n: int, spec: CandidateSpec, B: int, score_mode: str, seed: int, neighbor_safety: bool = True) -> PrefilterResult:
    t0 = time.time()
    score_fn_local = (lambda g: score_phase_eisenstein(g, n)) if score_mode == "phase" else (lambda g: score_l1_eisenstein(g, n))
    ordered = top_b_candidates_stream(candidates_eisenstein_stream(n, spec, seed=seed), score_fn_local, B)

    tested = 0
    for idx, g in enumerate(ordered, start=1):
        tested += 1
        f = try_factor_eisenstein(n, g, neighbor_safety=neighbor_safety)
        if f is not None:
            return PrefilterResult(True, f, idx, tested, time.time()-t0, "eisenstein")
    return PrefilterResult(False, None, None, tested, time.time()-t0, "eisenstein")

def prefilter_pipeline(n: int, spec: CandidateSpec, B: int, score_mode: str, seed: int, neighbor_safety: bool = True) -> PrefilterResult:
    cls = abc_like_class(n)

    # C-like Shadow-Shift
    if cls == "C":
        for k in (5,7):
            N = n * k
            res = prefilter_pipeline(N, spec, B, score_mode, seed=seed+1000+k, neighbor_safety=neighbor_safety)
            if res.success and res.factor:
                f = gcd(res.factor, n)
                if 1 < f < n:
                    return PrefilterResult(True, f, res.rank_hit, res.tested, res.t_prefilter, "shadow")
        return PrefilterResult(False, None, None, 0, 0.0, "shadow")

    if cls == "A":
        return prefilter_gaussian_fast(n, spec, B, score_mode, seed=seed)
    if cls == "B":
        return prefilter_eisenstein_fast(n, spec, B, score_mode, seed=seed, neighbor_safety=neighbor_safety)

    if cls == "E":
        m_align = min(spec.dir_count, 120)
        seed_align = seed + 2222
        Ag = ring_alignment(n, m_align, seed_align, "gaussian", topk=8)
        Ae = ring_alignment(n, m_align, seed_align, "eisenstein", topk=8)
        Bg, Be = split_budget(B, Ag, Ae)

        if Ag >= Ae:
            r1 = prefilter_gaussian_fast(n, spec, Bg, score_mode, seed=seed+1) if Bg > 0 else PrefilterResult(False, None, None, 0, 0.0, "none")
            if r1.success: return r1
            r2 = prefilter_eisenstein_fast(n, spec, Be, score_mode, seed=seed+2, neighbor_safety=neighbor_safety) if Be > 0 else PrefilterResult(False, None, None, 0, 0.0, "none")
        else:
            r1 = prefilter_eisenstein_fast(n, spec, Be, score_mode, seed=seed+1, neighbor_safety=neighbor_safety) if Be > 0 else PrefilterResult(False, None, None, 0, 0.0, "none")
            if r1.success: return r1
            r2 = prefilter_gaussian_fast(n, spec, Bg, score_mode, seed=seed+2) if Bg > 0 else PrefilterResult(False, None, None, 0, 0.0, "none")

        if r2.success:
            r2.tested += r1.tested
            r2.t_prefilter += r1.t_prefilter
            return r2
        return PrefilterResult(False, None, None, r1.tested + r2.tested, r1.t_prefilter + r2.t_prefilter, "E-both")

    return prefilter_gaussian_fast(n, spec, B, score_mode, seed=seed)

# ============================================================
# 8) Primality test switch: Miller-Rabin ON/OFF
# ============================================================

def trial_division_is_prime(n: int, limit: int = 200000) -> bool:
    # VERY slow for large n; intended only when MR is disabled and n is small.
    if n < 2:
        return False
    small = [2,3,5,7,11,13,17,19,23,29,31,37]
    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False
    r = isqrt(n)
    if r > limit:
        # avoid getting stuck forever; caller can decide what to do
        return False
    f = 41
    while f <= r:
        if n % f == 0:
            return False
        f += 2
    return True

def miller_rabin_is_prime(n: int, rounds: int = 12, rng: Optional[random.Random] = None) -> bool:
    if n < 2:
        return False
    small = [2,3,5,7,11,13,17,19,23,29,31,37]
    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False

    # n-1 = d*2^s
    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2

    if rng is None:
        rng = random.Random(0)

    def witness(a: int) -> bool:
        x = pow(a, d, n)
        if x == 1 or x == n-1:
            return False
        for _ in range(s-1):
            x = (x*x) % n
            if x == n-1:
                return False
        return True

    for _ in range(rounds):
        a = rng.randrange(2, n-1)
        if witness(a):
            return False
    return True

def is_prime(n: int, use_mr: bool, rng: random.Random) -> bool:
    if use_mr:
        return miller_rabin_is_prime(n, rounds=12, rng=rng)
    return trial_division_is_prime(n, limit=200000)

# ============================================================
# 9) Pollard-Rho fallback (optional) that uses is_prime()
# ============================================================

def pollard_rho(n: int, rng: random.Random) -> int:
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    while True:
        c = rng.randrange(1, n-1)
        x = rng.randrange(0, n-1)
        y = x
        d = 1
        f = lambda v: (pow(v, 2, n) + c) % n
        while d == 1:
            x = f(x)
            y = f(f(y))
            d = math.gcd(abs(x-y), n)
        if d != n:
            return d

def factorize(n: int, use_mr: bool, rng: random.Random) -> List[int]:
    # returns prime factors (not necessarily sorted)
    if n == 1:
        return []
    if is_prime(n, use_mr, rng):
        return [n]
    d = pollard_rho(n, rng)
    return factorize(d, use_mr, rng) + factorize(n//d, use_mr, rng)

# ============================================================
# 10) Semiprime generation for benchmarks
# ============================================================

def random_probable_prime(bits: int, rng: random.Random, use_mr: bool) -> int:
    # When MR is off, this becomes "find a small-ish prime" and may fail for larger bits.
    # We'll just keep trying and rely on is_prime with chosen backend.
    while True:
        x = rng.getrandbits(bits) | 1 | (1 << (bits-1))
        if is_prime(x, use_mr, rng):
            return x

def make_semiprime(bits: int, rng: random.Random, use_mr: bool) -> Tuple[int,int,int]:
    p = random_probable_prime(bits//2, rng, use_mr)
    q = random_probable_prime(bits//2, rng, use_mr)
    return p*q, p, q

# ============================================================
# 11) Runner (single mode + compare mode)
# ============================================================

def run_once(n: int, spec: CandidateSpec, B: int, score_mode: str, seed: int,
             neighbor_safety: bool, use_fallback: bool, use_mr: bool) -> Dict:
    rng = random.Random(seed)
    cls = abc_like_class(n)

    EI_STATS.div_calls = 0
    EI_STATS.neighbor_calls = 0

    t0 = time.time()
    res = prefilter_pipeline(n, spec, B, score_mode, seed=rng.randrange(1,10**9), neighbor_safety=neighbor_safety)
    t_pref = time.time() - t0

    out = {
        "n": n,
        "cls": cls,
        "prefilter_ok": res.success,
        "prefilter_factor": res.factor,
        "prefilter_ring": res.ring_used,
        "prefilter_tested": res.tested,
        "t_prefilter": t_pref,
        "ei_div_calls": EI_STATS.div_calls,
        "ei_neighbor_calls": EI_STATS.neighbor_calls,
    }

    if res.success:
        out["final_factors"] = sorted([res.factor, n//res.factor])
        out["t_total"] = t_pref
        out["fallback_used"] = False
        return out

    if not use_fallback:
        out["final_factors"] = None
        out["t_total"] = t_pref
        out["fallback_used"] = False
        return out

    # fallback factorization
    t1 = time.time()
    facs = sorted(factorize(n, use_mr=use_mr, rng=rng))
    t_fallback = time.time() - t1

    out["final_factors"] = facs
    out["fallback_used"] = True
    out["t_fallback"] = t_fallback
    out["t_total"] = t_pref + t_fallback
    return out

def format_row(cols: List[str], widths: List[int]) -> str:
    return "  ".join(c.ljust(w) for c,w in zip(cols, widths))

def print_table(rows: List[List[str]]):
    widths = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    for idx, r in enumerate(rows):
        print(format_row(r, widths))
        if idx == 0:
            print(format_row(["-"*w for w in widths], widths))

def main():
    ap = argparse.ArgumentParser(description="Fast ring-prefilter benchmark with MR on/off compare.")
    ap.add_argument("--trials", type=int, default=20, help="number of semiprimes")
    ap.add_argument("--bits", type=int, default=64, help="bit-size of semiprimes (try 64/96/128)")
    ap.add_argument("--B", type=int, default=400, help="budget for prefilter")
    ap.add_argument("--score", choices=["phase","l1"], default="phase", help="cheap candidate score")
    ap.add_argument("--neighbor-safety", action="store_true", help="enable Eisenstein neighbor safety (adaptive)")
    ap.add_argument("--fallback", action="store_true", help="enable Pollard-Rho fallback factorization")
    ap.add_argument("--mr", choices=["on","off"], default="on", help="primality backend for prime gen and fallback")
    ap.add_argument("--compare", action="store_true", help="compare same n-set with MR on vs MR off")
    ap.add_argument("--seed", type=int, default=1234567, help="master seed")
    args = ap.parse_args()

    spec = CandidateSpec()

    # Build a fixed dataset of semiprimes (generated with MR on for reliability)
    rng = random.Random(args.seed)
    dataset: List[int] = []
    ps: List[int] = []
    qs: List[int] = []

    # For fairness: create dataset with MR=ON always (ensures actual primes for p,q),
    # then run algorithms with MR toggled on/off on the *same* n.
    for i in range(args.trials):
        n, p, q = make_semiprime(args.bits, rng, use_mr=True)
        dataset.append(n); ps.append(p); qs.append(q)

    def run_suite(use_mr: bool) -> Dict:
        total_t = 0.0
        pref_succ = 0
        fallback_used = 0
        rows = [["i","cls","pref_ok","ring","tested","t_pref","fallback","t_total","factors"]]
        for i, n in enumerate(dataset, start=1):
            res = run_once(
                n=n, spec=spec, B=args.B, score_mode=args.score, seed=args.seed + i*99991,
                neighbor_safety=args.neighbor_safety,
                use_fallback=args.fallback,
                use_mr=use_mr
            )
            total_t += res["t_total"]
            if res["prefilter_ok"]:
                pref_succ += 1
            if res.get("fallback_used"):
                fallback_used += 1
            facs = res["final_factors"]
            facs_s = "-" if facs is None else ("*".join(str(x) for x in facs) if len(facs) <= 4 else f"{facs[0]}*...*{facs[-1]}")
            rows.append([
                str(i),
                res["cls"],
                "Y" if res["prefilter_ok"] else "N",
                res["prefilter_ring"],
                str(res["prefilter_tested"]),
                f"{res['t_prefilter']:.4f}",
                "Y" if res.get("fallback_used") else "N",
                f"{res['t_total']:.4f}",
                facs_s
            ])
        return {
            "rows": rows,
            "total_t": total_t,
            "mean_t": total_t / len(dataset),
            "pref_succ": pref_succ,
            "fallback_used": fallback_used
        }

    if args.compare:
        suite_on = run_suite(use_mr=True)
        suite_off = run_suite(use_mr=False)

        print("\n=== Dataset (generated with MR=ON) ===")
        print(f"trials={args.trials}  bits={args.bits}  B={args.B}  score={args.score}  neighbor_safety={args.neighbor_safety}  fallback={args.fallback}")

        print("\n=== Run A: MR=ON ===")
        print_table(suite_on["rows"][: min(len(suite_on["rows"]), 1 + min(12, args.trials))])
        print(f"\nSummary MR=ON: prefilter_success={suite_on['pref_succ']}/{args.trials}  fallback_used={suite_on['fallback_used']}/{args.trials}  mean_total_time={suite_on['mean_t']:.4f}s")

        print("\n=== Run B: MR=OFF (trial division prime test) ===")
        print_table(suite_off["rows"][: min(len(suite_off["rows"]), 1 + min(12, args.trials))])
        print(f"\nSummary MR=OFF: prefilter_success={suite_off['pref_succ']}/{args.trials}  fallback_used={suite_off['fallback_used']}/{args.trials}  mean_total_time={suite_off['mean_t']:.4f}s")

        print("\nNote: With MR=OFF, fallback factorization/primality checks can become very slow or may refuse to label large numbers as prime.")
        return

    use_mr = (args.mr == "on")
    print("\n=== Single run mode ===")
    print(f"trials={args.trials}  bits={args.bits}  B={args.B}  score={args.score}  neighbor_safety={args.neighbor_safety}  fallback={args.fallback}  mr={args.mr}")

    suite = run_suite(use_mr=use_mr)
    print_table(suite["rows"])
    print(f"\nSummary: prefilter_success={suite['pref_succ']}/{args.trials}  fallback_used={suite['fallback_used']}/{args.trials}  mean_total_time={suite['mean_t']:.4f}s")

if __name__ == "__main__":
    main()
