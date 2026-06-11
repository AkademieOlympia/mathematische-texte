from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import List, Tuple, Optional, Dict, Iterable
import random
import time
import math
import csv
from math import isqrt, sqrt

# Optional plots
try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False


# ============================================================
# 0) Utilities: rounding, gcd, simple Miller-Rabin, primes
# ============================================================

def round_nearest_rational(num: int, den: int) -> int:
    """
    Round num/den to nearest integer, ties away from zero.
    den > 0.
    """
    assert den > 0
    if num >= 0:
        q, r = divmod(num, den)
        if 2*r > den:  return q + 1
        if 2*r < den:  return q
        return q + 1
    else:
        return -round_nearest_rational(-num, den)

def igcd(a: int, b: int) -> int:
    return math.gcd(a, b)

def is_probable_prime(n: int, k: int = 16, rng: Optional[random.Random] = None) -> bool:
    """
    Probabilistic Miller-Rabin.
    """
    if n < 2:
        return False
    small_primes = [2,3,5,7,11,13,17,19,23,29,31,37]
    for p in small_primes:
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
        rng = random.Random()

    def witness(a: int) -> bool:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            return False
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                return False
        return True  # composite witness

    for _ in range(k):
        a = rng.randrange(2, n-1)
        if witness(a):
            return False
    return True

def random_prime(bits: int, rng: Optional[random.Random] = None) -> int:
    if rng is None:
        rng = random.Random()
    while True:
        x = rng.getrandbits(bits) | 1 | (1 << (bits-1))
        if is_probable_prime(x, k=16, rng=rng):
            return x

def jacobi(a: int, n: int) -> int:
    """
    Jacobi symbol (a/n) for odd n > 0.
    """
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

def sigma_class(n: int) -> str:
    """
    Placeholder E/A/B/C classification via two characters.
    You can replace mapping to your exact scheme.
    """
    j1 = jacobi(-1, n)
    j3 = jacobi(-3, n)
    mapping = {
        ( 1,  1): "E",
        ( 1, -1): "A",
        (-1,  1): "B",
        (-1, -1): "C",
    }
    return mapping.get((j1, j3), "U")


# ============================================================
# 1) Gaussian integers Z[i]: GI + divmod + gcd
# ============================================================

@dataclass(frozen=True)
class GI:
    a: int
    b: int

    def __add__(self, other: GI) -> GI:
        return GI(self.a + other.a, self.b + other.b)

    def __sub__(self, other: GI) -> GI:
        return GI(self.a - other.a, self.b - other.b)

    def __mul__(self, other: GI) -> GI:
        return GI(self.a*other.a - self.b*other.b, self.a*other.b + self.b*other.a)

    def conj(self) -> GI:
        return GI(self.a, -self.b)

    def norm(self) -> int:
        return self.a*self.a + self.b*self.b

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

def gi_divmod(alpha: GI, beta: GI) -> Tuple[GI, GI]:
    if beta.is_zero():
        raise ZeroDivisionError("division by zero in Z[i]")
    N = beta.norm()
    z = alpha * beta.conj()
    qa = round_nearest_rational(z.a, N)
    qb = round_nearest_rational(z.b, N)
    q = GI(qa, qb)
    r = alpha - (q * beta)
    return q, r

def gi_normalize(z: GI) -> GI:
    """
    Normalize up to units {±1, ±i}.
    """
    if z.is_zero():
        return z
    candidates = [
        z,
        GI(-z.a, -z.b),
        GI(-z.b, z.a),
        GI(z.b, -z.a),
    ]
    for c in candidates:
        if c.a > 0 or (c.a == 0 and c.b > 0):
            return c
    return candidates[0]

def gi_gcd(x: GI, y: GI) -> GI:
    a, b = x, y
    while not b.is_zero():
        _, r = gi_divmod(a, b)
        a, b = b, r
    return gi_normalize(a)


# ============================================================
# 2) Eisenstein integers Z[ω]: EI + divmod + gcd
# ============================================================

@dataclass(frozen=True)
class EI:
    a: int
    b: int  # a + b*ω

    def __add__(self, other: EI) -> EI:
        return EI(self.a + other.a, self.b + other.b)

    def __sub__(self, other: EI) -> EI:
        return EI(self.a - other.a, self.b - other.b)

    def __mul__(self, other: EI) -> EI:
        # (a+bω)(c+dω) = (ac-bd) + (ad+bc-bd)ω
        a, b = self.a, self.b
        c, d = other.a, other.b
        return EI(a*c - b*d, a*d + b*c - b*d)

    def conj(self) -> EI:
        # conj(a+bω) = (a-b) + (-b)ω
        return EI(self.a - self.b, -self.b)

    def norm(self) -> int:
        # a^2 - ab + b^2
        a, b = self.a, self.b
        return a*a - a*b + b*b

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def to_complex(self) -> complex:
        # ω = -1/2 + i*sqrt(3)/2
        return complex(self.a - self.b/2.0, (sqrt(3)/2.0)*self.b)

def nearest_eisenstein_from_complex(z: complex) -> EI:
    """
    Nearest Eisenstein integer to complex z using cube-coordinate rounding.
    """
    x, y = z.real, z.imag
    v = (2.0 / sqrt(3.0)) * y
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

    return EI(int(ru), int(rv))

def ei_divmod(alpha: EI, beta: EI, neighbor_safety: bool = True) -> Tuple[EI, EI]:
    if beta.is_zero():
        raise ZeroDivisionError("division by zero in Z[ω]")

    N = beta.norm()
    t = alpha * beta.conj()         # Eisenstein integer
    z = t.to_complex() / float(N)   # approximate complex ratio

    q0 = nearest_eisenstein_from_complex(z)

    if not neighbor_safety:
        r = alpha - (q0 * beta)
        return q0, r

    # Safety-net: check q0 and 6 neighbors, pick minimal remainder norm
    neighbors = [
        EI(0,0),
        EI(1,0), EI(-1,0),
        EI(0,1), EI(0,-1),
        EI(1,1), EI(-1,-1),  # ±(1+ω)
    ]

    best_q = q0
    best_r = alpha - (q0 * beta)
    best_norm = best_r.norm()

    for delta in neighbors:
        q = EI(q0.a + delta.a, q0.b + delta.b)
        r = alpha - (q * beta)
        nr = r.norm()
        if nr < best_norm:
            best_norm = nr
            best_q = q
            best_r = r

    return best_q, best_r

def ei_normalize(z: EI) -> EI:
    """
    Normalize up to units {±1, ±ω, ±ω^2}.
    """
    if z.is_zero():
        return z

    def mul_omega(t: EI) -> EI:
        # (a+bω)*ω = (-b) + (a-b)ω
        return EI(-t.b, t.a - t.b)

    def mul_minus_one(t: EI) -> EI:
        return EI(-t.a, -t.b)

    candidates = []
    t = z
    for _ in range(3):
        candidates.append(t)
        candidates.append(mul_minus_one(t))
        t = mul_omega(t)

    for c in candidates:
        if c.a > 0 or (c.a == 0 and c.b > 0):
            return c
    return candidates[0]

def ei_gcd(x: EI, y: EI, neighbor_safety: bool = True) -> EI:
    a, b = x, y
    while not b.is_zero():
        _, r = ei_divmod(a, b, neighbor_safety=neighbor_safety)
        a, b = b, r
    return ei_normalize(a)


# ============================================================
# 3) E8 direction set (real 240 roots) + projection -> 2D dirs
# ============================================================

def e8_roots() -> List[Tuple[float, ...]]:
    """
    Generate the 240 roots of E8:
      - 112 roots: permutations of (±1, ±1, 0,0,0,0,0,0)
      - 128 roots: (±1/2,...,±1/2) with even number of minus signs
    Returned as float tuples length 8.
    """
    roots = []

    # 112: choose i<j positions for ±1,±1
    for i in range(8):
        for j in range(i+1, 8):
            for si in (-1.0, 1.0):
                for sj in (-1.0, 1.0):
                    v = [0.0]*8
                    v[i] = si
                    v[j] = sj
                    roots.append(tuple(v))

    # 128: half-integer with even number of minus signs
    # Each coord is ±1/2
    for mask in range(1<<8):
        # count minus bits
        minus = 0
        v = []
        for k in range(8):
            if (mask >> k) & 1:
                v.append(-0.5)
                minus += 1
            else:
                v.append(0.5)
        if minus % 2 == 0:
            roots.append(tuple(v))

    assert len(roots) == 240
    return roots

def project_dirs_from_e8(m: int, seed: int = 0) -> List[Tuple[int,int]]:
    """
    Project E8 roots to 2D directions, then quantize to integer direction vectors.
    For reproducibility we pick a seeded random 2D projection.
    """
    rng = random.Random(seed)
    roots = e8_roots()

    # Random projection vectors u,v in R^8
    u = [rng.uniform(-1,1) for _ in range(8)]
    v = [rng.uniform(-1,1) for _ in range(8)]

    def dot(a,b):
        return sum(x*y for x,y in zip(a,b))

    # Compute 2D projected points for all roots
    pts = []
    for r in roots:
        x = dot(r, u)
        y = dot(r, v)
        pts.append((x,y))

    # Convert to integer directions (dx,dy) with gcd reduction
    dirs = []
    for x,y in pts:
        # scale to integers
        scale = 10_000
        dx = int(round(x*scale))
        dy = int(round(y*scale))
        if dx == 0 and dy == 0:
            continue
        g = math.gcd(abs(dx), abs(dy))
        dx //= g
        dy //= g
        # avoid duplicates up to sign
        if (dx,dy) in dirs or (-dx,-dy) in dirs:
            continue
        dirs.append((dx,dy))
        if len(dirs) >= m:
            break

    # If we didn't get enough (unlikely), pad with surrogate
    if len(dirs) < m:
        dirs += surrogate_dirs(m-len(dirs), seed=seed+1)
    return dirs[:m]

def surrogate_dirs(m: int, seed: int = 0) -> List[Tuple[int,int]]:
    """
    Simple surrogate integer direction set in Z^2.
    """
    rng = random.Random(seed)
    dirs = []
    while len(dirs) < m:
        x = rng.randint(-200,200)
        y = rng.randint(-200,200)
        if x == 0 and y == 0:
            continue
        g = math.gcd(abs(x), abs(y))
        x //= g; y //= g
        if (x,y) in dirs or (-x,-y) in dirs:
            continue
        dirs.append((x,y))
    return dirs

def get_directions(kind: str, m: int, seed: int) -> List[Tuple[int,int]]:
    if kind == "surrogate":
        return surrogate_dirs(m, seed=seed)
    if kind == "e8":
        return project_dirs_from_e8(m, seed=seed)
    raise ValueError(f"unknown direction kind: {kind}")


# ============================================================
# 4) Candidate generation (shared API) for GI and EI
# ============================================================

@dataclass(frozen=True)
class CandidateSpec:
    # Base neighborhood sampling
    T: int = 50            # real jitter range
    S: int = 50            # imag/omega jitter range
    base_samples: int = 25 # number of (t,s) pairs to try

    # Direction offsets
    dir_kind: str = "e8"   # "e8" or "surrogate"
    dir_count: int = 240
    lambdas: Tuple[int,...] = (1,2,3,5,8,13,21)  # scaling list

    # Optional sigma-dependent target / hooks
    use_sigma: bool = True

def base_pairs(spec: CandidateSpec, rng: random.Random) -> List[Tuple[int,int]]:
    pairs = []
    for _ in range(spec.base_samples):
        t = rng.randint(-spec.T, spec.T)
        s = rng.randint(-spec.S, spec.S)
        pairs.append((t,s))
    return pairs

def candidates_gaussian(n: int, spec: CandidateSpec, seed: int) -> List[GI]:
    rng = random.Random(seed)
    base = isqrt(n)
    dirs = get_directions(spec.dir_kind, spec.dir_count, seed=seed+1234)
    pairs = base_pairs(spec, rng)

    out: List[GI] = []
    for (t,s) in pairs:
        g0 = GI(base + t, s)
        out.append(g0)
        for (dx,dy) in dirs:
            for lam in spec.lambdas:
                out.append(GI(g0.a + lam*dx, g0.b + lam*dy))
    return out

def candidates_eisenstein(n: int, spec: CandidateSpec, seed: int) -> List[EI]:
    rng = random.Random(seed)
    base = isqrt(n)
    dirs = get_directions(spec.dir_kind, spec.dir_count, seed=seed+4321)
    pairs = base_pairs(spec, rng)

    out: List[EI] = []
    for (t,s) in pairs:
        g0 = EI(base + t, s)  # base + s*ω
        out.append(g0)
        for (dx,dy) in dirs:
            for lam in spec.lambdas:
                # interpret (dx,dy) as axial increments a + b*ω
                out.append(EI(g0.a + lam*dx, g0.b + lam*dy))
    return out


# ============================================================
# 5) Score functions (plug-in point for your Torus/Interferenzscore)
# ============================================================

def score_l1_gaussian(g: GI, n: int, sig: str) -> float:
    # small |Re-base| and |Im| preferred
    base = isqrt(n)
    return abs(g.a - base) + abs(g.b)

def score_l1_eisenstein(g: EI, n: int, sig: str) -> float:
    base = isqrt(n)
    return abs(g.a - base) + abs(g.b)

def score_phase_placeholder_gaussian(g: GI, n: int, sig: str) -> float:
    """
    Placeholder for a 'torus phase' style score:
    map (a,b) modulo M and prefer some sigma-dependent target.
    """
    M = 97
    # sigma-dependent "target" (completely heuristic placeholder)
    targets = {"E": (0,0), "A": (M//3, M//3), "B": (2*M//3, M//3), "C": (M//2, 2*M//3)}
    tx, ty = targets.get(sig, (0,0))
    x = g.a % M
    y = g.b % M
    # squared distance on torus (wrap-around)
    dx = min((x-tx) % M, (tx-x) % M)
    dy = min((y-ty) % M, (ty-y) % M)
    return dx*dx + dy*dy

def score_phase_placeholder_eisenstein(g: EI, n: int, sig: str) -> float:
    M = 97
    targets = {"E": (0,0), "A": (M//3, M//3), "B": (2*M//3, M//3), "C": (M//2, 2*M//3)}
    ta, tb = targets.get(sig, (0,0))
    a = g.a % M
    b = g.b % M
    da = min((a-ta) % M, (ta-a) % M)
    db = min((b-tb) % M, (tb-b) % M)
    return da*da + db*db

def get_score_fn(ring: str, name: str):
    if ring == "gaussian":
        if name == "l1":
            return score_l1_gaussian
        if name == "phase":
            return score_phase_placeholder_gaussian
    if ring == "eisenstein":
        if name == "l1":
            return score_l1_eisenstein
        if name == "phase":
            return score_phase_placeholder_eisenstein
    raise ValueError(f"unknown score fn: ring={ring}, name={name}")


# ============================================================
# 6) Prefilter: ring-gcd -> int gcd -> factor
# ============================================================

def factor_from_ring_norm(n: int, d_norm: int) -> Optional[int]:
    f = igcd(d_norm, n)
    if f != 1 and f != n:
        return f
    return None

def try_factor_gaussian(n: int, g: GI) -> Optional[int]:
    d = gi_gcd(GI(n, 0), g)
    return factor_from_ring_norm(n, d.norm())

def try_factor_eisenstein(n: int, g: EI, neighbor_safety: bool = True) -> Optional[int]:
    d = ei_gcd(EI(n, 0), g, neighbor_safety=neighbor_safety)
    return factor_from_ring_norm(n, d.norm())

@dataclass
class PrefilterResult:
    success: bool
    factor: Optional[int]
    rank_hit: Optional[int]
    tested: int
    t_prefilter: float
    ring_used: str


# ============================================================
# 7) Baseline fallback: Pollard-rho (simple)
# ============================================================

def pollard_rho(n: int, rng: random.Random, max_iters: int = 100_000) -> Optional[int]:
    if n % 2 == 0:
        return 2
    if is_probable_prime(n, k=16, rng=rng):
        return n

    while True:
        c = rng.randrange(1, n-1)
        x = rng.randrange(2, n-1)
        y = x
        d = 1

        f = lambda v: (pow(v, 2, n) + c) % n

        it = 0
        while d == 1 and it < max_iters:
            x = f(x)
            y = f(f(y))
            d = igcd(abs(x - y), n)
            it += 1

        if d == n or d == 1:
            # restart
            continue
        return d


# ============================================================
# 8) Runner: dataset generation, ablations, logging, plots
# ============================================================

@dataclass
class TrialRow:
    bits: int
    n: int
    sigma: str
    ring: str
    use_sigma: bool
    dir_kind: str
    dir_count: int
    lambdas: str
    base_samples: int
    T: int
    S: int
    score: str
    use_score: bool
    B: int
    neighbor_safety: bool
    shadow: bool

    success_prefilter: bool
    factor_prefilter: Optional[int]
    rank_hit: Optional[int]
    tested: int
    t_prefilter: float

    success_pipeline: bool
    factor_pipeline: Optional[int]
    t_total: float
    fallback_used: bool

def choose_ring(sig: str, use_sigma: bool) -> str:
    """
    Simple policy: if using sigma, pick gaussian for {E,B}, eisenstein for {A,C}.
    Otherwise default gaussian.
    Adjust to match your scheme.
    """
    if not use_sigma:
        return "gaussian"
    if sig in ("E", "B"):
        return "gaussian"
    if sig in ("A", "C"):
        return "eisenstein"
    return "gaussian"

def prefilter(n: int, spec: CandidateSpec, B: int, ring: str, score_name: str,
              use_score: bool, neighbor_safety: bool, seed: int) -> PrefilterResult:
    t0 = time.time()
    rng = random.Random(seed)
    sig = sigma_class(n)

    # generate candidates
    if ring == "gaussian":
        cands = candidates_gaussian(n, spec, seed=seed)
        score_fn = get_score_fn("gaussian", score_name)
        if use_score:
            scored = [(score_fn(g, n, sig), g) for g in cands]
            scored.sort(key=lambda x: x[0])
            ordered = [g for _, g in scored]
        else:
            rng.shuffle(cands)
            ordered = cands
        ordered = ordered[:min(B, len(ordered))]

        tested = 0
        for idx, g in enumerate(ordered, start=1):
            tested += 1
            f = try_factor_gaussian(n, g)
            if f is not None:
                t1 = time.time()
                return PrefilterResult(True, f, idx, tested, t1-t0, ring)

    elif ring == "eisenstein":
        cands = candidates_eisenstein(n, spec, seed=seed)
        score_fn = get_score_fn("eisenstein", score_name)
        if use_score:
            scored = [(score_fn(g, n, sig), g) for g in cands]
            scored.sort(key=lambda x: x[0])
            ordered = [g for _, g in scored]
        else:
            rng.shuffle(cands)
            ordered = cands
        ordered = ordered[:min(B, len(ordered))]

        tested = 0
        for idx, g in enumerate(ordered, start=1):
            tested += 1
            f = try_factor_eisenstein(n, g, neighbor_safety=neighbor_safety)
            if f is not None:
                t1 = time.time()
                return PrefilterResult(True, f, idx, tested, t1-t0, ring)
    else:
        raise ValueError("ring must be gaussian or eisenstein")

    t1 = time.time()
    return PrefilterResult(False, None, None, len(ordered), t1-t0, ring)

def gen_dataset(bits: int, count: int, seed: int) -> List[Tuple[int,int,int,str]]:
    """
    Returns list of (n,p,q,sigma). p,q included for truth but you can omit in logs.
    """
    rng = random.Random(seed)
    out = []
    for _ in range(count):
        p = random_prime(bits//2, rng=rng)
        q = random_prime(bits//2, rng=rng)
        n = p*q
        sig = sigma_class(n)
        out.append((n,p,q,sig))
    return out

def run_experiments(
    bits_list: Iterable[int] = (64, 96, 128),
    samples_per_bits: int = 200,
    seed: int = 0,
    B_list: Iterable[int] = (50, 200, 1000),
    score_list: Iterable[str] = ("l1", "phase"),
    ablations: Optional[List[Dict]] = None,
    output_csv: str = "results.csv",
    make_plots: bool = True
) -> List[TrialRow]:

    if ablations is None:
        # Minimal but telling ablation set:
        ablations = [
            dict(use_sigma=True,  use_score=True,  dir_kind="e8",       neighbor_safety=True,  shadow=False),
            dict(use_sigma=True,  use_score=False, dir_kind="e8",       neighbor_safety=True,  shadow=False),
            dict(use_sigma=True,  use_score=True,  dir_kind="surrogate",neighbor_safety=True,  shadow=False),
            dict(use_sigma=False, use_score=True,  dir_kind="e8",       neighbor_safety=True,  shadow=False),
        ]

    rows: List[TrialRow] = []
    rng = random.Random(seed)

    # CSV header
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(TrialRow.__annotations__.keys()))
        writer.writeheader()

        for bits in bits_list:
            dataset = gen_dataset(bits, samples_per_bits, seed=rng.randrange(1,10**9))

            for (n,p,q,sig) in dataset:
                for score_name in score_list:
                    for B in B_list:
                        for a in ablations:
                            # Spec per trial (you can also ablate dir_count, lambdas etc.)
                            spec = CandidateSpec(
                                T=50, S=50, base_samples=10,
                                dir_kind=a["dir_kind"], dir_count=240,
                                lambdas=(1,2,3,5,8,13,21),
                                use_sigma=a["use_sigma"]
                            )

                            ring = choose_ring(sig, use_sigma=a["use_sigma"])
                            t_all0 = time.time()

                            # Shadow hook: currently off in algebra (keep for later extension)
                            shadow = a.get("shadow", False)
                            n_eff = n  # if shadow: n*k (implement if desired)

                            pre = prefilter(
                                n=n_eff, spec=spec, B=B,
                                ring=ring, score_name=score_name,
                                use_score=a["use_score"],
                                neighbor_safety=a["neighbor_safety"],
                                seed=rng.randrange(1,10**9)
                            )

                            fallback_used = False
                            factor_pipe = pre.factor
                            success_pipe = pre.success

                            if not pre.success:
                                # fallback Pollard-rho (simple)
                                fallback_used = True
                                f_rho = pollard_rho(n_eff, rng=rng, max_iters=50_000)
                                if f_rho is not None and f_rho not in (1, n_eff):
                                    factor_pipe = f_rho
                                    success_pipe = True
                                else:
                                    factor_pipe = None
                                    success_pipe = False

                            t_all1 = time.time()

                            row = TrialRow(
                                bits=bits, n=n, sigma=sig,
                                ring=ring, use_sigma=a["use_sigma"],
                                dir_kind=a["dir_kind"], dir_count=spec.dir_count,
                                lambdas=",".join(map(str, spec.lambdas)),
                                base_samples=spec.base_samples, T=spec.T, S=spec.S,
                                score=score_name, use_score=a["use_score"],
                                B=B, neighbor_safety=a["neighbor_safety"],
                                shadow=shadow,

                                success_prefilter=pre.success,
                                factor_prefilter=pre.factor,
                                rank_hit=pre.rank_hit,
                                tested=pre.tested,
                                t_prefilter=pre.t_prefilter,

                                success_pipeline=success_pipe,
                                factor_pipeline=factor_pipe,
                                t_total=t_all1 - t_all0,
                                fallback_used=fallback_used
                            )

                            rows.append(row)
                            writer.writerow(asdict(row))

    if make_plots and HAS_MPL:
        plot_summary(rows, prefix="plot_")
    return rows


# ============================================================
# 9) Plotting: success_rate vs B, expected time vs B, hit-rank CDF
# ============================================================

def plot_summary(rows: List[TrialRow], prefix: str = "plot_") -> None:
    # Grouping keys
    def key_base(r: TrialRow):
        return (r.bits, r.score, r.use_sigma, r.use_score, r.dir_kind, r.ring)

    # Build sets
    bits_set = sorted(set(r.bits for r in rows))
    score_set = sorted(set(r.score for r in rows))
    B_set = sorted(set(r.B for r in rows))

    # 1) Success rate vs B (prefilter)
    for bits in bits_set:
        for score in score_set:
            plt.figure()
            for (use_sigma, use_score, dir_kind, ring) in sorted(set((r.use_sigma, r.use_score, r.dir_kind, r.ring) for r in rows if r.bits==bits and r.score==score)):
                ys = []
                for B in B_set:
                    sub = [r for r in rows if r.bits==bits and r.score==score and r.use_sigma==use_sigma and r.use_score==use_score and r.dir_kind==dir_kind and r.ring==ring and r.B==B]
                    if not sub:
                        ys.append(0.0)
                        continue
                    ys.append(sum(1 for r in sub if r.success_prefilter)/len(sub))
                plt.plot(B_set, ys, marker="o", label=f"{ring} | sigma={use_sigma} score={use_score} dirs={dir_kind}")
            plt.xlabel("Budget B")
            plt.ylabel("Prefilter success rate")
            plt.title(f"Prefilter success vs B | bits={bits} | score={score}")
            plt.legend()
            plt.tight_layout()
            plt.savefig(f"{prefix}success_prefilter_bits{bits}_{score}.png")
            plt.close()

    # 2) Expected total time vs B (pipeline)
    for bits in bits_set:
        for score in score_set:
            plt.figure()
            for (use_sigma, use_score, dir_kind, ring) in sorted(set((r.use_sigma, r.use_score, r.dir_kind, r.ring) for r in rows if r.bits==bits and r.score==score)):
                ys = []
                for B in B_set:
                    sub = [r for r in rows if r.bits==bits and r.score==score and r.use_sigma==use_sigma and r.use_score==use_score and r.dir_kind==dir_kind and r.ring==ring and r.B==B]
                    if not sub:
                        ys.append(0.0)
                        continue
                    ys.append(sum(r.t_total for r in sub)/len(sub))
                plt.plot(B_set, ys, marker="o", label=f"{ring} | sigma={use_sigma} score={use_score} dirs={dir_kind}")
            plt.xlabel("Budget B")
            plt.ylabel("Mean pipeline time (s)")
            plt.title(f"Mean pipeline time vs B | bits={bits} | score={score}")
            plt.legend()
            plt.tight_layout()
            plt.savefig(f"{prefix}time_pipeline_bits{bits}_{score}.png")
            plt.close()

    # 3) Hit-rank CDF for successful prefilter cases (one plot per bits/score)
    for bits in bits_set:
        for score in score_set:
            plt.figure()
            for (use_sigma, use_score, dir_kind, ring) in sorted(set((r.use_sigma, r.use_score, r.dir_kind, r.ring) for r in rows if r.bits==bits and r.score==score)):
                hits = [r.rank_hit for r in rows if r.bits==bits and r.score==score and r.use_sigma==use_sigma and r.use_score==use_score and r.dir_kind==dir_kind and r.ring==ring and r.success_prefilter and r.rank_hit is not None]
                if len(hits) < 10:
                    continue
                hits_sorted = sorted(hits)
                xs = hits_sorted
                ys = [(i+1)/len(hits_sorted) for i in range(len(hits_sorted))]
                plt.plot(xs, ys, label=f"{ring} | sigma={use_sigma} score={use_score} dirs={dir_kind}")
            plt.xlabel("Hit rank (smaller is better)")
            plt.ylabel("CDF")
            plt.title(f"Hit-rank CDF | bits={bits} | score={score}")
            plt.legend()
            plt.tight_layout()
            plt.savefig(f"{prefix}rank_cdf_bits{bits}_{score}.png")
            plt.close()


# ============================================================
# 10) Main entry (adjust parameters here)
# ============================================================

if __name__ == "__main__":
    rows = run_experiments(
        bits_list=(64, 96, 128),
        samples_per_bits=100,     # raise to 500/1000 later
        seed=0,
        B_list=(50, 200, 1000),
        score_list=("l1", "phase"),
        output_csv="results.csv",
        make_plots=True
    )
    print(f"Done. Wrote {len(rows)} rows to results.csv")
    if HAS_MPL:
        print("Plots written as plot_*.png")
    else:
        print("Matplotlib not available; skipped plots.")