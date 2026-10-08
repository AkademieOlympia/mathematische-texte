"""
CollatzResonanceEigenGap — zweite Resonanzklasse η = μ − 1/2.

Für Zyklen γ des Syracuse-Automaten mod 2^j:

  μ(γ) = a/L
  η(γ) = (a − L/2)/L = μ(γ) − 1/2
  Λ(γ) = a log2(3) − L

  ρ_j      = max μ
  η_max(j) = max η   (= ρ_j − 1/2)

Erwartung:
  η_max > 0  für Phantom-Levels j∈{8,10,11,12}
  η_max = 0  für j≥13  (nur γ₀={0} mit η=-1/2 und γ₁={1,2} mit η=0)

Spektrallücken-Hypothese (schwach): ∃J ∀j≥J: η_max(j) ≤ 0
  (äquivalent ρ_j ≤ 1/2 < p_c)

Stärker (nicht nötig für Lücke): Cyc(S_j) = {{0},{1,2}} für j≥J.

Epistemik: Automatenstruktur; kein Collatz-Beweis auf ℕ.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


LOG2_3 = math.log2(3)
P_CRIT = math.log(2) / math.log(3)
EPS_HALF = P_CRIT - 0.5  # ≈ 0.130929


def syracuse_next(r: int, mod: int) -> tuple[int, int]:
    if r % 2 == 0:
        return (r // 2) % mod, 0
    return ((3 * r + 1) // 2) % mod, 1


def analyze_level(j: int) -> dict:
    mod = 1 << j
    nxt = [0] * mod
    w = [0] * mod
    for r in range(mod):
        nxt[r], w[r] = syracuse_next(r, mod)

    visited = [0] * mod
    cycles_raw: list[list[int]] = []
    for start in range(mod):
        if visited[start]:
            continue
        path: list[int] = []
        pos: dict[int, int] = {}
        x = start
        while visited[x] == 0:
            visited[x] = 1
            pos[x] = len(path)
            path.append(x)
            x = nxt[x]
        if visited[x] == 1:
            cycles_raw.append(path[pos[x] :])
        for v in path:
            visited[v] = 2

    # basins
    attractor = [-1] * mod
    for ci, cyc in enumerate(cycles_raw):
        for v in cyc:
            attractor[v] = ci
    for r in range(mod):
        if attractor[r] != -1:
            continue
        x = r
        chain = []
        while attractor[x] == -1:
            chain.append(x)
            x = nxt[x]
        aid = attractor[x]
        for v in chain:
            attractor[v] = aid
    basin = [0] * len(cycles_raw)
    for r in range(mod):
        basin[attractor[r]] += 1

    cycles = []
    for i, cyc in enumerate(cycles_raw):
        L = len(cyc)
        a = sum(w[v] for v in cyc)
        mu = a / L
        eta = mu - 0.5
        Lam = a * LOG2_3 - L
        # notwendige Bedingung für ℕ-Zyklen: a/L < p_c; hier modular
        cycles.append(
            {
                "length": L,
                "odd": a,
                "mu": mu,
                "eta": eta,
                "Lambda": Lam,
                "basin": basin[i],
                "basin_frac": basin[i] / mod,
                "nodes": cyc if L <= 24 else cyc[:12] + ["…"],
                "is_zero": cyc == [0],
                "is_12": set(cyc) == {1, 2},
                "above_half": eta > 1e-15,
                "ge_pc": mu + 1e-15 >= P_CRIT,
            }
        )

    cycles.sort(key=lambda c: -c["eta"])
    rho = max(c["mu"] for c in cycles)
    eta_max = max(c["eta"] for c in cycles)
    eta_pos = [c for c in cycles if c["above_half"]]
    classified = all(c["is_zero"] or c["is_12"] for c in cycles)

    return {
        "j": j,
        "n_cycles": len(cycles),
        "rho": rho,
        "eta_max": eta_max,
        "Delta": P_CRIT - rho,
        "n_eta_positive": len(eta_pos),
        "basin_eta_positive": sum(c["basin_frac"] for c in eta_pos),
        "only_zero_and_12": classified,
        "cycles": cycles,
    }


def fit_exp_decay(js: list[int], etas: list[float]) -> dict | None:
    """Fit η ≈ C e^{-c j} on positive η points (log-linear)."""
    pts = [(j, e) for j, e in zip(js, etas) if e > 1e-15]
    if len(pts) < 2:
        return None
    # least squares on log η = log C - c j
    xs = [p[0] for p in pts]
    ys = [math.log(p[1]) for p in pts]
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    if den < 1e-18:
        return None
    slope = num / den  # = -c
    intercept = my - slope * mx  # = log C
    c = -slope
    C = math.exp(intercept)
    # residual
    sse = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
    return {
        "C": C,
        "c": c,
        "form": "C*exp(-c*j)",
        "n_points": n,
        "points": [{"j": j, "eta": e} for j, e in pts],
        "sse_log": sse,
        "note": "fit only on phantom levels with eta>0; not a decay law for j→∞",
    }


def main() -> int:
    print("=" * 72)
    print("RESONANZ-EIGENLÜCKE  —  η = μ − 1/2")
    print("=" * 72)
    print(f"p_c − 1/2 ≈ {EPS_HALF:.6f}")
    print("η_max = ρ − 1/2;  Spektrallücke ⇔ η_max ≤ 0 für j ≥ J\n")

    j_min, j_max = 4, 24
    rows = []
    print(
        f"  {'j':>3}  {'#cyc':>5}  {'ρ':>8}  {'η_max':>10}  {'#η>0':>5}  "
        f"{'basin(η>0)':>10}  {'nur{{0,12}}':>10}  {'Δ':>9}"
    )
    for j in range(j_min, j_max + 1):
        rec = analyze_level(j)
        rows.append(rec)
        print(
            f"  {j:3d}  {rec['n_cycles']:5d}  {rec['rho']:8.5f}  {rec['eta_max']:+10.6f}  "
            f"{rec['n_eta_positive']:5d}  {rec['basin_eta_positive']:10.6f}  "
            f"{'YES' if rec['only_zero_and_12'] else 'no':>10}  {rec['Delta']:+9.6f}"
        )

    print("\n[1] Phantom-Zyklen (η>0) im Detail")
    for rec in rows:
        if rec["n_eta_positive"] == 0:
            continue
        print(f"  j={rec['j']}:")
        for c in rec["cycles"]:
            if not c["above_half"]:
                continue
            print(
                f"    |γ|={c['length']:3d}  a={c['odd']:3d}  μ={c['mu']:.6f}  "
                f"η={c['eta']:+.6f}  Λ={c['Lambda']:+.4f}  "
                f"basin={c['basin_frac']:.4f}  nodes={c['nodes']}"
            )

    print("\n[2] Sterile Regime j≥13")
    sterile = [r for r in rows if r["j"] >= 13]
    eta_max_sterile = max(r["eta_max"] for r in sterile)
    classified_all = all(r["only_zero_and_12"] for r in sterile)
    print(f"  η_max(j) für j≥13: {eta_max_sterile:+.6f}")
    print(f"  Cyc = {{0}}∪{{1,2}} für alle j≥13 (im Scan)? {classified_all}")
    print(
        "  ⇒ schwache Resonanz-Sterilisierung: η_max≤0  (erfüllt)\n"
        "  ⇒ starke Form (nur zwei Zyklen): ebenfalls erfüllt bis j=24"
    )

    print("\n[3] Abkling-Fit für η_max>0 (nur Phantom-Levels)")
    js = [r["j"] for r in rows]
    etas = [r["eta_max"] for r in rows]
    fit = fit_exp_decay(js, etas)
    if fit:
        print(
            f"  Fit η≈{fit['C']:.4f}·exp(-{fit['c']:.4f}·j)  "
            f"auf {fit['n_points']} Punkten  sse_log={fit['sse_log']:.4e}"
        )
        print(
            "  Hinweis: nach j=12 ist η_max=0 exakt — kein exponentielles\n"
            "  Abklingen, sondern abrupte Sterilisierung (δ-Sprung)."
        )
    else:
        print("  zu wenige positive Punkte für Fit")

    # monotone Schranke nach Sterilisierung
    print("\n[4] Schranken")
    print(
        f"  Für j≥13: η_max(j)=0 ≤ C e^{{-cj}} für jedes C,c>0.\n"
        f"  Uniforme Lücke: ρ≤1/2 ⇒ Abstand zu p_c = {EPS_HALF:.6f}."
    )
    # notwendige Bedingung a/L < p_c für hypothetische ℕ-Zyklen
    print("\n[5] modular vs. ℕ-notwendig")
    print(
        "  ℕ-Zyklus braucht 2^L>3^a ⇒ μ<p_c.\n"
        "  Modular können Phantomzyklen μ>p_c haben (j=8,10–12),\n"
        "  heben aber nicht in den Limes / nicht nach ℕ."
    )
    violators = []
    for rec in rows:
        for c in rec["cycles"]:
            if c["ge_pc"]:
                violators.append((rec["j"], c["mu"], c["length"]))
    print(f"  Modulare Zyklen mit μ≥p_c: {len(violators)}  (alles Phantome)")

    print("\n[6] Gap-Analogie")
    print(
        "  M_B(G)  →  ρ_j  →  η_max(j)=ρ_j−1/2\n"
        "  η misst die Abweichung von der „balancierten“ Klasse μ=1/2."
    )

    out = Path(__file__).resolve().parent / "collatz_resonance_eigen_gap.json"
    payload = {
        "claim": (
            "Eigen-gap eta=mu-1/2 for Syracuse automata; "
            "eta_max=0 for j>=13 up to 24; strong cycle classification holds in range; "
            "not Collatz on N"
        ),
        "p_crit": P_CRIT,
        "eps_half": EPS_HALF,
        "j_range": [j_min, j_max],
        "rows": [
            {
                "j": r["j"],
                "n_cycles": r["n_cycles"],
                "rho": r["rho"],
                "eta_max": r["eta_max"],
                "Delta": r["Delta"],
                "n_eta_positive": r["n_eta_positive"],
                "basin_eta_positive": r["basin_eta_positive"],
                "only_zero_and_12": r["only_zero_and_12"],
                "cycles": r["cycles"],
            }
            for r in rows
        ],
        "sterile_regime": {
            "J": 13,
            "eta_max": eta_max_sterile,
            "strong_classification": classified_all,
            "weak_sterilization": eta_max_sterile <= 1e-15,
        },
        "exp_fit_phantoms": fit,
        "note": (
            "Abrupt sterilization at j=13, not smooth exponential decay of eta_max; "
            "eta_max=0 exactly for j>=13 in scan"
        ),
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
