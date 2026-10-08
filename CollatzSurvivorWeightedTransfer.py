"""
CollatzSurvivorWeightedTransfer — gewichtete Survivor-Transferlücke.

Hierarchie (nicht identisch):
  Entropielücke:  ρ(M_unweighted) < 2
  Driftlücke:     ρ(M_weighted) < 1   mit w = 2^{ΔΛ}

  Entropie ⇏ Drift; Drift ⇏ Entropie.
  Gesucht: gewichtete Survivor-Transferlücke
    sup_{ω∈Survivor} limsup Λ_L(ω)/L < 0
  + |x^*| Kontrolle.

S_L(C): Clog-Survivor (HitTime > C log₂ n + 2), konsistent mit Atlas.
ΔΛ beim Bit b:  -1 (b=0)  bzw.  log₂3 − 1 (b=1)
w = 2^{ΔΛ}:     1/2      bzw.  3/2

Epistemik: Spektraldiagnostik; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


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


def word_lambda(word: str) -> float:
    a = word.count("1")
    return a * math.log2(3) - len(word)


def enumerate_survivors_clog(L: int, C: float) -> list[dict]:
    mod = 1 << L
    out = []
    for r in range(mod):
        n = r if r > 1 else r + mod
        budget = int(max(C * math.log2(n) * 3 + 50, 500))
        ht = hit_time(n, budget)
        thresh = C * math.log2(n) + 2
        if not (ht is None or ht > thresh):
            continue
        x, bits, a = r, [], 0
        for _ in range(L):
            if x % 2 == 0:
                bits.append("0")
                x = (x // 2) % mod
            else:
                bits.append("1")
                a += 1
                x = ((3 * x + 1) // 2) % mod
        word = "".join(bits)
        out.append({"r": r, "n": n, "word": word, "a": a, "Lambda": word_lambda(word)})
    return out


def enumerate_nohit_L(L: int) -> list[dict]:
    """S_L ohne C: Orbit trifft {1,2} nicht in den ersten L Schritten."""
    mod = 1 << L
    out = []
    for r in range(mod):
        n = r if r > 1 else r + mod
        if hit_time(n, L) is not None:
            continue
        x, bits, a = r, [], 0
        for _ in range(L):
            if x % 2 == 0:
                bits.append("0")
                x = (x // 2) % mod
            else:
                bits.append("1")
                a += 1
                x = ((3 * x + 1) // 2) % mod
        word = "".join(bits)
        out.append({"word": word, "Lambda": word_lambda(word), "a": a, "n": n})
    return out


W_EVEN = 0.5
W_ODD = 1.5  # 2^{log2(3)-1} = 3/2
D_EVEN = -1.0
D_ODD = math.log2(3) - 1.0


def shift_matrices(survivors: list[dict]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    words = [s["word"] for s in survivors]
    idx = {w: i for i, w in enumerate(words)}
    n = len(words)
    M = np.zeros((n, n), dtype=float)
    Mw = np.zeros((n, n), dtype=float)
    for u in words:
        for b, wgt, _d in (("0", W_EVEN, D_EVEN), ("1", W_ODD, D_ODD)):
            v = u[1:] + b
            if v in idx:
                i, j = idx[u], idx[v]
                M[i, j] = 1.0
                Mw[i, j] = wgt
    return M, Mw, words


def spectral_radius(M: np.ndarray) -> float:
    if M.size == 0 or M.shape[0] == 0:
        return 0.0
    return float(np.max(np.abs(np.linalg.eigvals(M))))


def extension_weighted_growth(S_L: list[dict], S_Lp1: list[dict]) -> dict:
    """Rechteckige Extension: mittlere/maximale Gewichts-Verstärkung."""
    by = {s["word"]: s for s in S_Lp1}
    rows = []
    for s in S_L:
        w = s["word"]
        child_w = []
        for b, wgt in (("0", W_EVEN), ("1", W_ODD)):
            if w + b in by:
                child_w.append(wgt)
        rows.append(
            {
                "n_children": len(child_w),
                "sum_w": float(sum(child_w)),
                "max_w": float(max(child_w)) if child_w else 0.0,
            }
        )
    if not rows:
        return {"n": 0}
    sums = np.array([r["sum_w"] for r in rows])
    return {
        "n_nodes": len(rows),
        "mean_branch_weight": float(sums.mean()),  # ≈ lokal ρ-artig
        "max_branch_weight": float(sums.max()),
        "frac_branch_weight_lt_1": float(np.mean(sums < 1.0)),
        "frac_branch_weight_lt_2": float(np.mean(sums < 2.0)),
        "mean_n_children": float(np.mean([r["n_children"] for r in rows])),
    }


def lambda_rate_stats(survivors: list[dict]) -> dict:
    if not survivors:
        return {"n": 0}
    L = len(survivors[0]["word"])
    rates = np.array([s["Lambda"] / L for s in survivors])
    return {
        "L": L,
        "n": len(rates),
        "max_Lambda_over_L": float(rates.max()),
        "mean_Lambda_over_L": float(rates.mean()),
        "p90_Lambda_over_L": float(np.percentile(rates, 90)),
        "frac_Lambda_over_L_neg": float(np.mean(rates < 0)),
    }


def main() -> int:
    Cs = [8.0, 10.0, 12.0]
    levels = [8, 10, 12, 14]
    print("=" * 72)
    print("GEWICHTETE SURVIVOR-TRANSFERLÜCKE")
    print("  Entropie: ρ(M)<2     Drift: ρ(M_w)<1 mit w=2^{ΔΛ}")
    print("  Entropie ⇏ Drift (getrennte Spektren)")
    print("=" * 72)

    payload: dict = {
        "claim": (
            "Weighted survivor transfer gap is the precise next goal; "
            "entropy gap and drift gap are related but not equivalent; "
            "not a Collatz proof"
        ),
        "target": (
            "sup_{survivor} limsup Lambda_L/L < 0  via  rho(M_weighted)<1, "
            "plus x^* control"
        ),
        "weights": {"even": W_EVEN, "odd": W_ODD, "dLambda_even": D_EVEN, "dLambda_odd": D_ODD},
        "hierarchy_note": (
            "rho(M)<2 thins paths; rho(M_w)<1 forces negative mean resonance; "
            "neither implies the other alone"
        ),
        "clog": {},
        "nohit_L": {},
    }

    for C in Cs:
        print(f"\n### Clog-Survivor C={C}")
        rows = []
        by_L = {L: enumerate_survivors_clog(L, C) for L in levels}
        for L in levels:
            S = by_L[L]
            M, Mw, _ = shift_matrices(S)
            rho = spectral_radius(M)
            rhow = spectral_radius(Mw)
            rates = lambda_rate_stats(S)
            rec = {
                "L": L,
                "n_S": len(S),
                "rho_unweighted": rho,
                "rho_unweighted_over_2": rho / 2.0,
                "rho_weighted": rhow,
                "entropy_gap": rho < 2.0 - 1e-12,
                "weighted_gap": rhow < 1.0 - 1e-12,
                **{k: rates[k] for k in rates if k != "n"},
            }
            # extension growth if L+1 available conceptually — use L to L for shift only
            rows.append(rec)
            print(
                f"  L={L:2d} |S|={len(S):5d}  ρ(M)={rho:7.4f}  ρ/2={rho/2:6.4f}  "
                f"ρ(M_w)={rhow:7.4f}  maxΛ/L={rates.get('max_Lambda_over_L')}"
            )
        # extension branch weights L→L+1
        ext = []
        for L in levels:
            S0 = by_L[L]
            S1 = enumerate_survivors_clog(L + 1, C)
            eg = extension_weighted_growth(S0, S1)
            eg["L"] = L
            eg["n_S_Lp1"] = len(S1)
            ext.append(eg)
            print(
                f"  ext {L}→{L+1}: mean_branch_w={eg.get('mean_branch_weight')}  "
                f"frac_wsum<1={eg.get('frac_branch_weight_lt_1')}  "
                f"mean_kids={eg.get('mean_n_children')}"
            )
        payload["clog"][str(C)] = {"shift": rows, "extension_weights": ext}

    print("\n### NoHit-by-L (C-frei; zum Vergleich — oft keine Entropielücke)")
    nh_rows = []
    for L in (8, 10, 12):
        S = enumerate_nohit_L(L)
        M, Mw, _ = shift_matrices(S)
        rho, rhow = spectral_radius(M), spectral_radius(Mw)
        rates = lambda_rate_stats(S)
        nh_rows.append(
            {
                "L": L,
                "n_S": len(S),
                "rho_unweighted": rho,
                "rho_unweighted_over_2": rho / 2.0,
                "rho_weighted": rhow,
                "max_Lambda_over_L": rates.get("max_Lambda_over_L"),
            }
        )
        print(
            f"  L={L} |NoHit|={len(S)}  ρ(M)={rho:.4f}  ρ/2={rho/2:.4f}  "
            f"ρ(M_w)={rhow:.4f}  maxΛ/L={rates.get('max_Lambda_over_L')}"
        )
    payload["nohit_L"] = nh_rows

    print("\n[Konsolidierte Lesart]")
    print(
        "  1) ρ(M)<2  = entropische Verdünnung (kombinatorisch)\n"
        "  2) ρ(M_w)<1 = keine nichtnegativ-resonante Eigenmode (dynamisch)\n"
        "  3) Zielsatz: limsup Λ_L/L < 0 uniform auf Survivor + |x*| control\n"
        "  4) Endliche S_L=∅ für C=16 ist hart endlich — nicht die Spektrallücke"
    )

    out = Path(
        "/Users/thomashoffbauer/Projects/Kepler-Hurrwitz/docs/exports/"
        "collatz_survivor_weighted_transfer.json"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    desk = Path(__file__).resolve().parent / "collatz_survivor_weighted_transfer.json"
    try:
        desk.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    except OSError:
        pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
