import numpy as np
import sympy as sp
import random

# ---------- Quaternion ----------
def quat_mult(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array(
        [
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
        ]
    )


def quat_inv(q):
    conj = np.array([q[0], -q[1], -q[2], -q[3]])
    return conj / np.dot(q, q)


# ---------- Rotation ----------
def make_u(theta, axis):
    axis = np.array(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    c = np.cos(theta / 2)
    s = np.sin(theta / 2)
    return np.array([c, axis[0] * s, axis[1] * s, axis[2] * s])


# ---------- Projektion ----------
def nearest_prime_unbiased(x):
    x0 = int(abs(x))
    if x0 < 3:
        return 2
    try:
        p1 = sp.prevprime(x0)
    except:
        p1 = 2
    p2 = sp.nextprime(x0)
    return p1 if abs(x - p1) < abs(x - p2) else p2


def project(Qv):
    return [nearest_prime_unbiased(x) for x in Qv]


def random_local():
    base = random.randint(10_000, 200_000)
    offsets = sorted(np.random.randint(0, 20, size=4))
    return [nearest_prime_unbiased(base + o) for o in offsets]


# ---------- Twins ----------
def twin_count(vec):
    vec = sorted(vec)
    count = 0
    for i in range(len(vec) - 1):
        if vec[i + 1] - vec[i] == 2 and sp.isprime(vec[i]) and sp.isprime(vec[i + 1]):
            count += 1
    return count


def gap_pattern(vec):
    v = sorted(vec)
    return [v[i + 1] - v[i] for i in range(3)]


def dist_to_242(g):
    return (g[0] - 2) ** 2 + (g[1] - 4) ** 2 + (g[2] - 2) ** 2


# ---------- Vierlinge ----------
def quadruplets(N):
    res = []
    for p in range(5, N):
        if all(sp.isprime(p + k) for k in [0, 2, 6, 8]):
            res.append(p)
    return res


def Q(p):
    return np.array([p, p + 2, p + 6, p + 8], dtype=float)


# ---------- Scan ----------
quad = quadruplets(200000)

n_pairs = len(quad) - 1
random_hits = 0
for _ in range(n_pairs):
    Qv = np.array([random.randint(5, 200_000) for _ in range(4)], dtype=float)
    if twin_count(project(Qv)) > 0:
        random_hits += 1
twin_rate_random = random_hits / len(quad)
print("twin_rate_random (baseline):", twin_rate_random)

axes = [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 1, 0),
    (1, 0, 1),
    (0, 1, 1),
    (1, 1, 1),
]

thetas = np.linspace(0.1, 2 * np.pi, 12)

results = []

for axis in axes:
    for theta in thetas:
        if min(abs(theta), abs(theta - 2 * np.pi)) < 0.2:
            continue
        u = make_u(theta, axis)
        u_inv = quat_inv(u)

        twin_hits = 0
        acc_d242 = 0.0
        acc_d242_r = 0.0

        for i in range(n_pairs):
            q1 = Q(quad[i])
            q2 = Q(quad[i + 1])
            mid = (q1 + q2) / 2

            Qr = quat_mult(quat_mult(u, q1 - mid), u_inv) + mid
            Qp = project(Qr)

            gaps = gap_pattern(Qp)
            d242 = dist_to_242(gaps)

            Qrnd = [
                nearest_prime_unbiased(random.randint(10_000, 200_000))
                for _ in range(4)
            ]
            gaps_r = gap_pattern(Qrnd)
            d242_r = dist_to_242(gaps_r)

            acc_d242 += d242
            acc_d242_r += d242_r

            if twin_count(Qp) > 0:
                twin_hits += 1

        mean_d242 = acc_d242 / n_pairs
        mean_d242_r = acc_d242_r / n_pairs
        twin_rate_BM = twin_hits / len(quad)
        score = twin_rate_BM - twin_rate_random
        print(
            "Δ(2-4-2 distance) =",
            acc_d242 / len(quad) - acc_d242_r / len(quad),
        )
        results.append(
            (axis, theta, twin_rate_BM, score, mean_d242, mean_d242_r)
        )
        print(
            axis,
            theta,
            twin_rate_BM,
            "score",
            score,
            "mean_d242",
            mean_d242,
            "mean_d242_r",
            mean_d242_r,
        )

# ---------- Top Ergebnisse (nach score) ----------
results.sort(key=lambda x: x[3], reverse=True)

print("\nTop 5 (score):")
for r in results[:5]:
    print(r)
