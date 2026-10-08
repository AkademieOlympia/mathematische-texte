"""
GapResonanceSymbolic — algebraische Triplet-Invarianten.

Für δ=(0,a,b) mit 0<a<b:
  M aus Gap-Differenzen,
  ρ² = |M1|²/3,
  T = tr(PMP),
  D = λ₁λ₂  (Produkt der nichttrivialen Eigenwerte von PMP).

Spiegelinvarianten:
  u = b
  v = a(b−a)

Geschlossene Formen (exakt):
  ρ² = 2(u⁴ + v²) / (3 u² v²)
  T  = −2(u² + v) / (3 u v)
  D  = −(u⁴ − 6 u² v + v²) / (3 u² v²)

Damit R(a,b) = F(u,v). Auf dem Kegel {u>0, 0<v≤u²/4} ist F injektiv:
  ρ²/T² = (3/2)·g(w), w=v/u², g(w)=(1+w²)/(1+w)² streng fallend auf (0,1],
  danach u = −2(1+w)/(3 T w) eindeutig. Klassifikation bis auf Spiegelung.

Claim-boundary: Gap-Geometrie; keine Primzahlspezifik.
Lean: KeplerHurwitz/GapResonanceTriplet.lean
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


a, b = sp.symbols("a b", positive=True)
u, v = sp.symbols("u v", positive=True)


def build_M(a_sym, b_sym):
    """Gap-Matrix für δ=(0,a,b), 0<a<b."""
    return sp.Matrix(
        [
            [0, 1 / a_sym, 1 / b_sym],
            [1 / a_sym, 0, 1 / (b_sym - a_sym)],
            [1 / b_sym, 1 / (b_sym - a_sym), 0],
        ]
    )


def invariants_ab():
    M = build_M(a, b)
    ones = sp.Matrix([1, 1, 1])
    M1 = sp.simplify(M * ones)
    rho2 = sp.factor(sp.simplify(M1.dot(M1) / 3))
    P = sp.eye(3) - ones * ones.T / 3
    PMP = sp.simplify(P * M * P)
    T = sp.factor(sp.simplify(PMP.trace()))
    # λ1 λ2 = Summe der 2×2-Hauptminoren (da det(PMP)=0, Eigenwerte 0,λ1,λ2)
    D = sp.factor(
        sp.simplify(
            PMP[0, 0] * PMP[1, 1]
            - PMP[0, 1] * PMP[1, 0]
            + PMP[0, 0] * PMP[2, 2]
            - PMP[0, 2] * PMP[2, 0]
            + PMP[1, 1] * PMP[2, 2]
            - PMP[1, 2] * PMP[2, 1]
        )
    )
    return rho2, T, D


def reduce_to_uv(expr):
    """Schreibe spiegelinvarianten Ausdruck als F(u,v), u=b, v=a(b-a)."""
    num, den = sp.fraction(sp.together(expr))
    x = sp.symbols("x")

    def reduce_poly(poly):
        p_expr = sp.expand(poly.subs(b, u)).subs(a, x)
        r = sp.rem(sp.Poly(p_expr, x), sp.Poly(x**2 - u * x + v, x))
        r_expr = sp.expand(r.as_expr())
        poly_r = sp.Poly(r_expr, x)
        Acoef = sp.factor(sp.simplify(poly_r.nth(1)))
        Bcoef = sp.factor(sp.simplify(poly_r.nth(0)))
        return Acoef, Bcoef

    An, Bn = reduce_poly(sp.expand(num))
    Ad, Bd = reduce_poly(sp.expand(den))
    if An != 0 or Ad != 0:
        raise ValueError(f"not mirror-invariant reduction: A_num={An}, A_den={Ad}")
    return sp.factor(sp.simplify(Bn / Bd))


def main() -> int:
    print("=" * 72)
    print("SYMBOLISCHE TRIPLET-INVARIANTEN  —  R = F(u,v)")
    print("=" * 72)

    rho2, T, D = invariants_ab()
    print("\n[1] Ausdrücke in (a,b)")
    print(f"  ρ² = {rho2}")
    print(f"  T  = {T}")
    print(f"  D  = {D}")

    for name, expr in (("ρ²", rho2), ("T", T), ("D", D)):
        mir = sp.simplify(expr.subs(a, b - a))
        ok = sp.simplify(expr - mir) == 0
        print(f"  {name} spiegelinvariant? {ok}")

    print("\n[2] Reduktion auf u=b, v=a(b−a)")
    rho2_uv = reduce_to_uv(rho2)
    T_uv = reduce_to_uv(T)
    D_uv = reduce_to_uv(D)
    print(f"  ρ² = {rho2_uv}")
    print(f"  T  = {T_uv}")
    print(f"  D  = {D_uv}")

    # Verifikation Rücksubstitution
    print("\n[3] Rücksubstitution")
    back = {
        "ρ²": sp.simplify(rho2 - rho2_uv.subs({u: b, v: a * (b - a)})) == 0,
        "T": sp.simplify(T - T_uv.subs({u: b, v: a * (b - a)})) == 0,
        "D": sp.simplify(D - D_uv.subs({u: b, v: a * (b - a)})) == 0,
    }
    for k, ok in back.items():
        print(f"  {k} recovers? {ok}")

    # Geschlossene Ziel-Formeln (kanonisch)
    rho2_closed = 2 * (u**4 + v**2) / (3 * u**2 * v**2)
    T_closed = -2 * (u**2 + v) / (3 * u * v)
    D_closed = -(u**4 - 6 * u**2 * v + v**2) / (3 * u**2 * v**2)
    print("\n[4] Kanonische geschlossene Formen")
    print(f"  ρ² = 2(u^4 + v^2) / (3 u^2 v^2)")
    print(f"  T  = -2(u^2 + v) / (3 u v)")
    print(f"  D  = -(u^4 - 6 u^2 v + v^2) / (3 u^2 v^2)")
    print(
        f"  match reduce? "
        f"ρ²={sp.simplify(rho2_uv - rho2_closed)==0}, "
        f"T={sp.simplify(T_uv - T_closed)==0}, "
        f"D={sp.simplify(D_uv - D_closed)==0}"
    )

    # Numerische Stichprobe: gleiche (u,v) ⇒ gleiche R; Spiegelung teilt (u,v)
    print("\n[5] Stichprobe Spiegelklassen")
    samples = [(1, 2), (2, 6), (1, 5), (3, 7), (5, 12)]
    for aa, bb in samples:
        uu, vv = bb, aa * (bb - aa)
        aa2, bb2 = bb - aa, bb  # Spiegel
        vv2 = aa2 * (bb2 - aa2)
        r1 = float(rho2_closed.subs({u: uu, v: vv}))
        r2 = float(rho2_closed.subs({u: bb2, v: vv2}))
        print(
            f"  (0,{aa},{bb}) ~ (0,{aa2},{bb2}):  "
            f"u,v=({uu},{vv})  ρ²={r1:.10f}  mirror ρ²={r2:.10f}  "
            f"v-match={vv==vv2}"
        )

    # Invertierung: aus T,D → u,v?
    print("\n[6] Invertierung (T,D) → (u,v)")
    t_sym, d_sym = sp.symbols("t d", real=True)
    # v = -2 u^2 / (3 t u + 2)  from T
    Vsol = sp.simplify(sp.solve(sp.Eq(t_sym, T_closed), v)[0])
    print(f"  aus T: v = {Vsol}")
    # plug into D
    Deq = sp.Eq(d_sym, sp.simplify(D_closed.subs(v, Vsol)))
    Deq_n = sp.numer(sp.together(Deq.lhs - Deq.rhs))
    Deq_n = sp.factor(sp.simplify(sp.expand(Deq_n)))
    print(f"  D-Gleichung (Zähler in u): {Deq_n}")
    # Degree in u
    poly = sp.Poly(sp.expand(Deq_n), u)
    print(f"  Grad in u: {poly.degree()}")

    # Zeige: ρ²,T,D sind abhängig — eine Relation
    # Eliminiere v between T and rho2
    rel = sp.resultant(
        sp.numer(sp.together(T_closed - t_sym)),
        sp.numer(sp.together(rho2_closed - sp.symbols("r2"))),
        v,
    )
    print(f"  resultant_v(T−t, ρ²−r2) Grad_u={sp.degree(sp.expand(rel), u)}")

    # Injektivität: w = v/u², ρ²/T² = (3/2) g(w), g streng fallend auf (0,1]
    print("\n[7] Injektivität von F (algebraisch)")
    w = sp.symbols("w", positive=True)
    g = (1 + w**2) / (1 + w) ** 2
    # ρ²/T² = (3/2) g(w)
    ratio = sp.simplify(rho2_closed / T_closed**2)
    ratio_w = sp.simplify(ratio.subs(v, w * u**2))
    target = sp.Rational(3, 2) * g
    print(f"  rho2/T^2 |_{{v=w*u^2}} = {ratio_w}")
    print(f"  (3/2)g(w)              = {sp.simplify(target)}")
    print(f"  Identitaet? {sp.simplify(ratio_w - target) == 0}")
    # g(w)−g(w') Zähler = 2(w−w')(ww'−1) ⇒ streng fallend auf (0,1]
    wp = sp.symbols("w_p", positive=True)
    gdiff = sp.factor(
        sp.expand((1 + w**2) * (1 + wp) ** 2 - (1 + wp**2) * (1 + w) ** 2)
    )
    print(f"  g-Zählerdifferenz = {gdiff}")
    gdiff_ok = sp.simplify(gdiff - 2 * (w - wp) * (w * wp - 1)) == 0
    print(f"  = 2(w−w')(ww'−1)? {gdiff_ok}")
    # Rekonstruktion u aus T,w
    u_rec = sp.simplify(-2 * (1 + w) / (3 * T_closed.subs(v, w * u**2) * w))
    print(f"  u-Rekonstruktion = {u_rec}  (soll = u: {sp.simplify(u_rec - u) == 0})")

    print("\n[8] Hauptsatz-Status")
    print(
        "  BEWIESEN: R(a,b)=F(u,v), u=b, v=a(b−a).\n"
        "  BEWIESEN: F injektiv auf {u>0, 0<v≤u²/4}.\n"
        "  Damit R(a,b)=R(c,d) ⇔ (c,d)=(a,b) oder (c,d)=(b−a,b).\n"
        "  Lean: KeplerHurwitz/GapResonanceTriplet.lean"
    )

    # Raster als Sanity-Check
    print("\n[9] Raster-Sanity")

    def Fpack(uu, vv):
        return (
            float(rho2_closed.subs({u: uu, v: vv})),
            float(T_closed.subs({u: uu, v: vv})),
            float(D_closed.subs({u: uu, v: vv})),
        )

    pts = []
    for bb in range(2, 81):
        for aa in range(1, bb):
            if aa > bb - aa:
                continue
            uu, vv = bb, aa * (bb - aa)
            pts.append(((aa, bb), Fpack(uu, vv)))
    buckets = {}
    collisions = 0
    for ab, f in pts:
        key = tuple(round(x, 10) for x in f)
        if key in buckets and buckets[key] != ab:
            collisions += 1
        else:
            buckets[key] = ab
    print(
        f"  kanonische Punkte: {len(pts)}  Signaturen: {len(buckets)}  "
        f"Kollisionen@1e-10: {collisions}"
    )

    injectivity_ok = (
        sp.simplify(ratio_w - target) == 0
        and gdiff_ok
        and sp.simplify(u_rec - u) == 0
    )

    out = Path(__file__).resolve().parent / "gap_resonance_symbolic.json"
    payload = {
        "claim": (
            "Triplet resonance signature is F(u,v) with u=b, v=a(b-a); "
            "classifies up to reflection; not about primes"
        ),
        "rho2_uv": str(rho2_closed),
        "T_uv": str(T_closed),
        "D_uv": str(D_closed),
        "u": "b",
        "v": "a*(b-a)",
        "mirror": "(a,b) -> (b-a, b)",
        "proved": [
            "R(a,b) depends only on (u,v)",
            "F injective on {u>0, 0<v<=u^2/4} via rho2/T^2=(3/2)g(w)",
            "classification up to reflection a <-> b-a",
        ],
        "open": [],
        "lean": "KeplerHurwitz/GapResonanceTriplet.lean",
        "raster_collisions": collisions,
        "raster_points": len(pts),
        "injectivity_symbolic": injectivity_ok,
    }
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")

    all_ok = all(back.values()) and collisions == 0 and injectivity_ok
    print(f"\nVERDICT: {'PASS' if all_ok else 'PARTIAL'}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
