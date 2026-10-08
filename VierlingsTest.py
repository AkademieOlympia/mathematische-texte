"""
VierlingsTest — G(2,6,8) mit demselben Hurwitz-/Gap-Mechanismus wie DrillingsTest.

Klasse: q_n = (n, n+2, n+6, n+8)

Zwei Gap-Operatoren:
  scalar : alte skalare Toy-Kopplung (Typ II, B∼√5/n) — nur Referenz
  Rg4    : echte 4×4-Gap-Matrix R_g^{(4)} wie die 3×3-Form bei Drillingen

Dynamik:
  F_α(q) = Ad_u(q) · (1 + α B(q)),  B = R / √N
  u = (1+i+j+k)/2,  Ad_u(q) = u q u^{-1}

Governance:
  Gap-Topologie ≠ Primzahlspezifik. Vergleich M_B^{(3)} vs M_B^{(4)}.
"""

from __future__ import annotations

import math

import numpy as np

# Hurwitz-Einheit (wie DrillingsTest)
u = np.array([0.5, 0.5, 0.5, 0.5], dtype=float)


def q_four(n: float) -> np.ndarray:
    return np.array([n, n + 2, n + 6, n + 8], dtype=float)


def q_triple(n: float) -> np.ndarray:
    """Drilling G(2,4,6): (n,n+2,n+6,6) — Referenz."""
    return np.array([n, n + 2, n + 6, 6.0], dtype=float)


def qmul(q, r):
    a, b, c, d = q
    e, f, g, h = r
    return np.array(
        [
            a * e - b * f - c * g - d * h,
            a * f + b * e + c * h - d * g,
            a * g - b * h + c * e + d * f,
            a * h + b * g - c * f + d * e,
        ],
        dtype=float,
    )


def qconj(q):
    a, b, c, d = q
    return np.array([a, -b, -c, -d], dtype=float)


def qnorm(q):
    return float(np.dot(q, q))


def Ad_u(q):
    return qmul(qmul(u, q), qconj(u))


def normalize(q, eps: float = 1e-30):
    return q / math.sqrt(max(qnorm(q), eps))


# ---------------------------------------------------------------------------
# Gap-Operatoren
# ---------------------------------------------------------------------------


def R_g_scalar_toy(q: np.ndarray) -> float:
    """
    Alte skalare Toy-Kopplung (VierlingsTest v1):
      gaps=(b−a,c−b,d−c,d−a) → w → |w|.
    Auf G(2,6,8): |w|=2√5 konstant ⇒ Typ II.
    """
    a, b, c, d = (float(x) for x in q)
    gaps = np.array([b - a, c - b, d - c, d - a], dtype=float)
    w = np.array(
        [
            gaps[0] - gaps[2],
            gaps[1] - gaps[3],
            gaps[0] + gaps[1] - gaps[3],
        ],
        dtype=float,
    )
    return float(np.linalg.norm(w))


def R_g3_quat(q: np.ndarray, eps: float = 1e-30) -> np.ndarray:
    """
    Drilling-Operator (DrillingsTest): 3×3-Gap-Matrix auf (q0,q1,q2),
    Ergebnis als reiner Vektoranteil (0,w).
    """
    p0, p1, p2, _ = (float(x) for x in q)
    g01 = abs(p0 - p1)
    g02 = abs(p0 - p2)
    g12 = abs(p1 - p2)
    i01 = 1.0 / max(g01, eps)
    i02 = 1.0 / max(g02, eps)
    i12 = 1.0 / max(g12, eps)
    M = np.array(
        [[0.0, i01, i02], [i01, 0.0, i12], [i02, i12, 0.0]], dtype=float
    )
    w = M @ np.array([p0, p1, p2], dtype=float)
    return np.array([0.0, w[0], w[1], w[2]], dtype=float)


def R_g4_quat(q: np.ndarray, eps: float = 1e-30) -> np.ndarray:
    """
    Echte 4×4-Gap-Matrix auf allen Komponenten:
      M_ii=0, M_ij=1/|q_i−q_j|,  R = M q ∈ H ≅ R⁴.
    """
    x = np.asarray(q, dtype=float)
    M = np.zeros((4, 4), dtype=float)
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            M[i, j] = 1.0 / max(abs(x[i] - x[j]), eps)
    return M @ x


def drive_B(q: np.ndarray, mode: str) -> np.ndarray | float:
    nrm = math.sqrt(max(qnorm(q), 1e-30))
    if mode == "scalar":
        return R_g_scalar_toy(q) / nrm
    if mode == "Rg3":
        return R_g3_quat(q) / nrm
    if mode == "Rg4":
        return R_g4_quat(q) / nrm
    raise ValueError(mode)


def B_norm(q: np.ndarray, mode: str) -> float:
    B = drive_B(q, mode)
    if isinstance(B, float):
        return abs(B)
    return float(np.linalg.norm(B))


# ---------------------------------------------------------------------------
# Dynamik + exakte E_g-Kettenregel (wie DrillingsTest)
# ---------------------------------------------------------------------------


def hurwitz_step(q, alpha: float, mode: str) -> np.ndarray:
    rot = Ad_u(q)
    B = drive_B(q, mode)
    if isinstance(B, float):
        # skalare Toy: multipliziere mit (1+α B) als reeller Faktor
        return rot * (1.0 + alpha * B)
    return qmul(rot, np.array([1.0, 0.0, 0.0, 0.0]) + alpha * B)


def Eg_exact_norm(q: np.ndarray, mode: str) -> float:
    """
    |d/dα normalize(F_α³(q))|_{α=0}| via Kettenregel.
    Für mode=scalar: degenerierte skalare Summe (nur Referenz).
    """
    if mode == "scalar":
        # v1-Formel: E = V0+V1+V2 mit V=α·|w|/|q| — hier |E/α|
        q0, q1, q2 = q, Ad_u(q), Ad_u(Ad_u(q))
        return abs(
            B_norm(q0, "scalar") + B_norm(q1, "scalar") + B_norm(q2, "scalar")
        )

    q0 = q
    q1 = Ad_u(q0)
    q2 = Ad_u(q1)
    B0 = drive_B(q0, mode)
    B1 = drive_B(q1, mode)
    B2 = drive_B(q2, mode)
    dq1 = qmul(q1, B0)
    dq2 = qmul(q2, B1) + Ad_u(dq1)
    dq3 = qmul(Ad_u(q2), B2) + Ad_u(dq2)
    nrm = math.sqrt(max(qnorm(q0), 1e-30))
    x = q0 / nrm
    v = dq3 / nrm
    v = v - float(np.dot(v, x)) * x
    return float(np.linalg.norm(v))


# ---------------------------------------------------------------------------
# Asymptotik / Klassifikation
# ---------------------------------------------------------------------------


def classify_growth(mode: str, seed_fn, ns: list[int]) -> dict:
    """Schätzt deg: n·B und B für große n."""
    rows = []
    for n in ns:
        q = seed_fn(float(n))
        b = B_norm(q, mode)
        rows.append({"n": n, "B": b, "nB": n * b})
    B_last = rows[-1]["B"]
    nB_last = rows[-1]["nB"]
    # Heuristik
    if B_last > 0.05:
        typ = "I_resonant"
    elif nB_last > 0.1 and B_last < 0.05:
        typ = "II_critical"
    else:
        typ = "III_subcritical_or_small"
    return {"type": typ, "rows": rows, "B_infty_est": B_last, "nB_est": nB_last}


def asymptotic_Rg4_MB() -> float:
    """
    Analytische Asymptotik für R_g4 auf G(2,6,8):
      q = n·1 + δ, δ=(0,2,6,8)
      M konstant (Gaps n-unabhängig), M1 := M·1 ≠ 0
      R = M q = n M1 + Mδ  ⇒  |R|/|q| → |M1|/2
    """
    delta = np.array([0.0, 2.0, 6.0, 8.0])
    # baue M auf δ-Differenzen (= q-Differenzen)
    M = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            if i != j:
                M[i, j] = 1.0 / abs(delta[i] - delta[j])
    ones = np.ones(4)
    M1 = M @ ones
    return float(np.linalg.norm(M1) / 2.0)  # |q|∼2n


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------


def report_comparison() -> None:
    ns = [10, 100, 1000, 10000, 100000, 1000000]
    print("=" * 64)
    print("VIERLING vs DRILLING — echter Hurwitz-Gap-Stack")
    print("=" * 64)
    print(
        "claim-boundary: classifies gap growth degree; "
        "not a prime-quadruplet selector."
    )

    # --- Toy scalar (Typ II Referenz) ---
    print("\n[A] SCALAR TOY on G(2,6,8)  (alte Normierung)")
    print("    erwartet: B∼√5/n → 0  (Typ II)")
    for n in ns:
        b = B_norm(q_four(n), "scalar")
        print(f"    n={n:8d}  B={b:.6e}  nB={n*b:.6f}")
    print(f"    lim nB = √5 ≈ {math.sqrt(5):.6f}")

    # --- Rg4 echte Matrix ---
    MB4 = asymptotic_Rg4_MB()
    print("\n[B] R_g^{(4)} on G(2,6,8)  (Hurwitz-kompatibel)")
    print(f"    analytisch: M_B^(4) = |M·1|/2 = {MB4:.12f}")
    print("    n            B              nB           |B−M_B|")
    for n in ns:
        b = B_norm(q_four(n), "Rg4")
        print(
            f"    {n:8d}  {b:.10f}  {n*b:12.6f}  {abs(b-MB4):.3e}"
        )
    cls4 = classify_growth("Rg4", q_four, ns)
    print(f"    Klassifikation: {cls4['type']}")

    # --- Rg3 auf Drilling-Referenz ---
    MB3 = math.sqrt(85.0 / 216.0)
    print("\n[C] R_g^{(3)} on G(2,4,6)  (Drilling-Referenz)")
    print(f"    analytisch: M_B^(3) = √(85/216) = {MB3:.12f}")
    print("    n            B              |B−M_B|")
    for n in ns:
        b = B_norm(q_triple(n), "Rg3")
        print(f"    {n:8d}  {b:.10f}  {abs(b-MB3):.3e}")

    # --- Rg3 auf Vierling (nur erste 3 Komponenten) ---
    print("\n[D] R_g^{(3)} applied to Vierling (q0,q1,q2) only")
    print("    n            B")
    for n in ns:
        b = B_norm(q_four(n), "Rg3")
        print(f"    {n:8d}  {b:.10f}")

    # --- E_g Vergleich ---
    print("\n[E] |E_g| exact (chain rule) — Period-3")
    print(
        "    n         |E| Rg4     |E| Rg3-drill   "
        "|E|/M_B4    |E|/M_B3"
    )
    for n in [10, 100, 1000, 10000, 100000, 1000000]:
        e4 = Eg_exact_norm(q_four(n), "Rg4")
        e3 = Eg_exact_norm(q_triple(n), "Rg3")
        print(
            f"    {n:8d}  {e4:.8f}  {e3:.8f}  "
            f"{e4/MB4:.6f}  {e3/MB3:.6f}"
        )

    print("\n[F] VERDICT")
    b4_big = B_norm(q_four(1_000_000), "Rg4")
    if b4_big > 0.1:
        print(
            f"    G(2,6,8) mit R_g^(4) ist Typ I (resonant): "
            f"B→M_B^(4)≈{MB4:.6f}>0."
        )
        print(
            "    Die Degeneration B∼1/n war ein Artefakt der skalaren "
            "Toy-Kopplung (Differenzprojektion), nicht der Gap-Klasse."
        )
    else:
        print(
            "    G(2,6,8) bleibt auch unter R_g^(4) kritisch/subkritisch."
        )
    print(
        f"    M_B^(3)=√(85/216)≈{MB3:.6f}  vs  "
        f"M_B^(4)=|M1|/2≈{MB4:.6f}"
    )
    print(
        "    Primvierlinge tragen keine zusätzliche Information jenseits "
        "der Gap-Klasse G(2,6,8)."
    )


def report_prime_quads(limit: int = 5000) -> None:
    print("\n[G] Primvierlinge unter R_g^(4) (nur Klasse, kein Selektor)")
    print("    p           B_Rg4      B_scalar")
    count = 0
    for n in range(5, limit, 2):
        if all(
            _is_prime(n + o) for o in (0, 2, 6, 8)
        ):
            b4 = B_norm(q_four(n), "Rg4")
            bs = B_norm(q_four(n), "scalar")
            print(f"    {n:6d}  {b4:.8f}  {bs:.8e}")
            count += 1
            if count >= 12:
                break


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True


if __name__ == "__main__":
    report_comparison()
    report_prime_quads()
