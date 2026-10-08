"""
CollatzSurvivorEntropy — Survivor-Menge S_L und Λ-Spektrum.

Präzisierung:
  Finite Checks  ∀ n < 2^L₀ : HitTime(n) ≤ C log₂ n
  sind KEINE Zylinder-Zertifikate für alle L.

Survivor-Menge (relativ zu Schranke C·L):
  S_L = { ω ∈ {0,1}^L :
          ∃ n ∈ Z_ω (Repr. n = r_ω), HitTime(n) > C·L }

Bzw. strenger (log-n-Schranke) mit Repräsentant r:
  HitTime(r) > C·log₂(max(r,2)).

Messgrößen:
  ρ(L) = |S_L| / 2^L          (Überlebensdichte)
  h(L) = (1/L) log₂ |S_L|     (Entropierate)
  Λ_max(L) = max_{ω∈S_L} Λ_ω  (Drift der Survivor)

Vermutung (Spektrallücke):
  limsup h(L) < 1
  und/oder  Λ_max(L) ≤ −ε L  auf S_L (für große L)

Epistemik: Entropie-/Drift-Diagnostik; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


def syracuse_step(n: int) -> int:
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def hit_time(n: int, max_steps: int) -> int | None:
    if n in (1, 2):
        return 0
    x = n
    for k in range(1, max_steps + 1):
        x = syracuse_step(x)
        if x in (1, 2):
            return k
    return None


def word_and_affine_from_residue(r: int, L: int):
    """Wort + (a, Λ) entlang L Syracuse-Schritten vom Residuum r mod 2^L."""
    bits = []
    a = 0
    x = r
    mod = 1 << L
    for _ in range(L):
        if x % 2 == 0:
            bits.append(0)
            x = (x // 2) % mod
        else:
            bits.append(1)
            a += 1
            x = ((3 * x + 1) // 2) % mod
    Lam = a * math.log2(3) - L
    return "".join(map(str, bits)), a, Lam


def scan_survivors(L: int, C: float, mode: str = "CL") -> dict:
    """
    mode='CL':   Survivor wenn HitTime(r) > C*L
    mode='Clog': Survivor wenn HitTime(r) > C*log2(r)  (r≥2)
    Repräsentant: r = 0..2^L-1, r≥2 (bzw. r+2^L falls r<2).
    """
    mod = 1 << L
    survivors = []
    Lam_all = []
    budget = int(max(C * L + 50, C * L * 2, 200))

    for r in range(mod):
        n = r if r > 1 else r + mod  # positiver Repräsentant > 1
        word, a, Lam = word_and_affine_from_residue(r, L)
        Lam_all.append(Lam)
        ht = hit_time(n, max_steps=budget)
        if mode == "CL":
            thresh = C * L
            is_surv = ht is None or ht > thresh
        else:
            thresh = C * math.log2(n) + 2
            is_surv = ht is None or ht > thresh
        if is_surv:
            survivors.append(
                {
                    "r": r,
                    "n": n,
                    "word": word,
                    "a": a,
                    "Lambda": Lam,
                    "ht": ht,
                    "thresh": thresh,
                }
            )

    nS = len(survivors)
    rho = nS / mod
    h = (math.log2(nS) / L) if nS > 0 else float("-inf")
    Lams = [s["Lambda"] for s in survivors]
    return {
        "L": L,
        "C": C,
        "mode": mode,
        "n_cyl": mod,
        "n_S": nS,
        "rho": rho,
        "h": h,
        "Lambda_max_S": max(Lams) if Lams else None,
        "Lambda_mean_S": sum(Lams) / len(Lams) if Lams else None,
        "Lambda_min_S": min(Lams) if Lams else None,
        "Lambda_max_all": max(Lam_all),
        "Lambda_mean_all": sum(Lam_all) / len(Lam_all),
        "frac_S_Lambda_neg": (
            sum(1 for x in Lams if x < 0) / len(Lams) if Lams else None
        ),
        "worst": sorted(survivors, key=lambda s: -s["Lambda"])[:8],
        "budget": budget,
    }


def growth_theta(rows: list[dict], r: int = 2) -> list[dict]:
    """Schätze |S_{L+r}| / (|S_L| 2^r)."""
    by_L = {row["L"]: row for row in rows}
    out = []
    for L, row in sorted(by_L.items()):
        L2 = L + r
        if L2 not in by_L or row["n_S"] == 0:
            continue
        ratio = by_L[L2]["n_S"] / (row["n_S"] * (1 << r))
        out.append(
            {
                "L": L,
                "L2": L2,
                "theta_hat": ratio,
                "n_S_L": row["n_S"],
                "n_S_L2": by_L[L2]["n_S"],
            }
        )
    return out


def main() -> int:
    C = 16.0
    levels = [8, 10, 12, 14, 16]  # 2^16 = 65536 machbar
    print("=" * 72)
    print("SURVIVOR-ENTROPIE  —  S_L, ρ(L), h(L), Λ_max(S_L)")
    print("=" * 72)
    print("Finite Checks ∀n<2^L₀ ≠ Zertifikat ∀L.  Frage: Entropielücke?\n")

    rows_CL, rows_Clog = [], []
    for L in levels:
        for mode, bucket in (("CL", rows_CL), ("Clog", rows_Clog)):
            rec = scan_survivors(L, C=C, mode=mode)
            bucket.append(rec)
            print(
                f"  [{mode:4s}] L={L:2d}  |S|={rec['n_S']:6d}  "
                f"ρ={rec['rho']:.4f}  h={rec['h']:.4f}  "
                f"Λmax_S={rec['Lambda_max_S']}  "
                f"Λmean_S={rec['Lambda_mean_S']}  "
                f"fracΛ<0={rec['frac_S_Lambda_neg']}"
            )

    print("\n[Wachstum θ̂ = |S_{L+2}| / (|S_L|·4)]  mode=Clog")
    thetas = growth_theta(rows_Clog, r=2)
    for t in thetas:
        print(
            f"  L={t['L']}→{t['L2']}:  θ̂={t['theta_hat']:.4f}  "
            f"|S| {t['n_S_L']}→{t['n_S_L2']}"
        )

    print("\n[Λ auf Survivors vs alle Zylinder] mode=Clog")
    for rec in rows_Clog:
        print(
            f"  L={rec['L']:2d}  Λmax_S={rec['Lambda_max_S']}  "
            f"Λmax_all={rec['Lambda_max_all']:.3f}  "
            f"Λmean_all={rec['Lambda_mean_all']:.3f}"
        )

    # Spektrallücken-Diagnose
    print("\n[Diagnose]")
    hs = [r["h"] for r in rows_Clog if r["n_S"] > 0]
    rhos = [r["rho"] for r in rows_Clog]
    print(f"  h(L) Clog: {['%.3f'%x for x in hs]}")
    print(f"  ρ(L) Clog: {['%.4f'%x for x in rhos]}")
    if all(r["n_S"] == 0 for r in rows_Clog):
        print("  S_L leer für Clog → keine Survivor unter C=16 auf diesen L")
        print("  (konsistent mit M(L)/L≤15.4, aber NICHT ∀L)")
    else:
        print(
            "  Gesucht: limsup h(L)<1 und/oder Λmax(S_L)≤−εL\n"
            "  (echte Entropielücke / Driftlücke auf Nicht-Hit-Ästen)"
        )

    # CL-Mode: oft viele Survivors weil C*L << typische HitTime für große n
    print("\n[Hinweis mode=CL]")
    print(
        "  Schranke C·L ist für große Repräsentanten r≈2^L oft zu klein\n"
        "  (HitTime ~ c log r ~ c L, aber Konstante). Clog ist die\n"
        "  zur LogHitBound passende Survivor-Definition."
    )

    print("\n[Zwischensatz-Kandidat]")
    print(
        "  Jeder unendliche Nicht-Hit-Zylinderpfad erfüllt\n"
        "  liminf Λ_L < 0  (oder stärker: Λ_L ≤ −εL).\n"
        "  Dann |A|≤2^{−εL} ⇒ R fällt in den Hit-Zylinder.\n"
        "  Analog Gap: spec(PMP) ↔ Λ auf Survivor-Zylindern."
    )

    out = Path(__file__).resolve().parent / "collatz_survivor_entropy.json"
    payload = {
        "claim": (
            "Survivor entropy S_L; finite L0 checks do not imply forall L; "
            "seek entropy gap limsup h<1 or Lambda gap on survivors; "
            "not a Collatz proof"
        ),
        "C": C,
        "levels": levels,
        "CL": rows_CL,
        "Clog": rows_Clog,
        "theta_Clog": thetas,
        "candidate_theorem": (
            "Every infinite non-hit cylinder path has liminf Lambda_L < 0 "
            "(or Lambda_L <= -eps L), implying radius fall into hit cylinders"
        ),
        "gap_analogy": "Gap spec(PMP) <-> Collatz Lambda on survivor cylinders",
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
