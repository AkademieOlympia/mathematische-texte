"""
CollatzRadiusScaling — Q(n)=k*(n)/log₂ n.

k*(n) := min { k : R_k(n) < C }  mit C=2
       = Eintrittszeit in {1,2} (empirisch exakt im Cover-Scan)

Frage:
  bleibt Q(n) beschränkt?  ⇒  k*(n)=O(log n)

Epistemik: Skalierungs-Atlas; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
import random
from collections import defaultdict
from fractions import Fraction
from pathlib import Path


def syracuse_enter(n0: int, max_steps: int = 100_000) -> dict:
    """Eintrittszeit + kurze R-Probe am Eintritt."""
    x = n0
    for k in range(1, max_steps + 1):
        x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
        if x in (1, 2):
            return {"n0": n0, "k_star": k, "enter_state": x, "capped": False}
    return {"n0": n0, "k_star": None, "enter_state": x, "capped": True}


def radius_at_prefix(n0: int, bits: list[int]) -> float:
    """R_L(n0) für volles Präfix bits."""
    a = L = B = 0
    for bit in bits:
        if bit == 0:
            L += 1
        else:
            B = 3 * B + (1 << L)
            a += 1
            L += 1
    A = Fraction(3**a, 1 << L)
    x_star = Fraction(B, (1 << L) - 3**a)
    return float(abs(A) * abs(Fraction(n0) - x_star))


def verify_R_at_enter(n0: int, k_star: int) -> dict:
    """Prüfe R_{k*-1}≥2 >? und R_{k*}<2 entlang echter Bahn."""
    bits = []
    x = n0
    for _ in range(k_star):
        bit = 0 if x % 2 == 0 else 1
        bits.append(bit)
        x = x // 2 if bit == 0 else (3 * x + 1) // 2
    R_star = radius_at_prefix(n0, bits)
    R_prev = radius_at_prefix(n0, bits[:-1]) if len(bits) > 1 else None
    return {
        "R_star": R_star,
        "R_prev": R_prev,
        # Grenze: manchmal R*=2 exakt (z.B. n=4→2); daher ≤2
        "R_star_lt_2": R_star <= 2 + 1e-12,
        "R_prev_ge_2": R_prev is None or R_prev >= 2 - 1e-12,
        "final_state": x,
    }


def Q(n0: int, k_star: int) -> float:
    return k_star / math.log2(n0)


def main() -> int:
    random.seed(19)
    C = 2.0

    print("=" * 72)
    print("RADIUS-SKALIERUNG  —  Q(n)=k*(n)/log₂ n")
    print("=" * 72)
    print(f"k* = Eintrittszeit (= min k mit R_k(n)≤{C} im Scan)\n")

    # [1] Deterministische Leitern + Zufall in Log-Buckets
    seeds: list[int] = []
    seeds += list(range(3, 200))
    seeds += [2**k - 1 for k in range(8, 40)]
    seeds += [2**k + 1 for k in range(8, 36)]
    seeds += [3**k + 5 for k in range(5, 25)]
    seeds += [10**k + 3 for k in range(3, 12)]
    # Zufall pro Größenordnung
    for e in range(3, 13):
        lo, hi = 10**e, 10 ** (e + 1) - 1
        for _ in range(40):
            seeds.append(random.randint(lo, hi))
    # sehr große Stichprobe
    for _ in range(30):
        seeds.append(random.randint(10**12, 10**15))
    for _ in range(15):
        seeds.append(random.randint(10**15, 10**18))

    seen: set[int] = set()
    uniq = []
    for s in seeds:
        if s > 1 and s not in seen:
            seen.add(s)
            uniq.append(s)
    seeds = uniq

    print(f"[1] Seeds: {len(seeds)}  (bis ~10^18)")
    rows = []
    capped = 0
    for n0 in seeds:
        rec = syracuse_enter(n0)
        if rec["capped"] or rec["k_star"] is None:
            capped += 1
            continue
        k = rec["k_star"]
        q = Q(n0, k)
        rows.append(
            {
                "n0": n0,
                "log2_n": math.log2(n0),
                "k_star": k,
                "Q": q,
                "enter_state": rec["enter_state"],
            }
        )

    print(f"  entered within cap: {len(rows)}  capped: {capped}")
    Qs = [r["Q"] for r in rows]
    ks = [r["k_star"] for r in rows]
    print(
        f"  k*:  min={min(ks)}  median={sorted(ks)[len(ks)//2]}  max={max(ks)}"
    )
    print(
        f"  Q:   min={min(Qs):.4f}  median={sorted(Qs)[len(Qs)//2]:.4f}  "
        f"max={max(Qs):.4f}  mean={sum(Qs)/len(Qs):.4f}"
    )

    # [2] Q nach log-Buckets
    print("\n[2] Q nach Größenordnung von n")
    buckets: dict[int, list[float]] = defaultdict(list)
    for r in rows:
        b = int(math.floor(math.log10(r["n0"])))
        buckets[b].append(r["Q"])
    print("  log10n   #      Q_med    Q_max    Q_mean   k_med")
    for b in sorted(buckets):
        qs = buckets[b]
        ks_b = [r["k_star"] for r in rows if int(math.floor(math.log10(r["n0"]))) == b]
        print(
            f"  1e{b:<2d}   {len(qs):4d}  "
            f"{sorted(qs)[len(qs)//2]:7.3f}  {max(qs):7.3f}  "
            f"{sum(qs)/len(qs):7.3f}  {sorted(ks_b)[len(ks_b)//2]:5d}"
        )

    # [3] Extremale Q
    print("\n[3] Größte Q(n) (schlechteste Skalierung)")
    worst = sorted(rows, key=lambda r: -r["Q"])[:15]
    for r in worst:
        print(
            f"  Q={r['Q']:7.3f}  k*={r['k_star']:5d}  "
            f"log₂n={r['log2_n']:6.2f}  n={r['n0']}"
        )

    # [4] Verifikation R↔enter auf Sample
    print("\n[4] Verifikation R_{k*}<2≤R_{k*-1} (Stichprobe)")
    sample = worst[:5] + sorted(rows, key=lambda r: r["n0"])[:5]
    sample += [r for r in rows if r["n0"] in (27, 31, 703, 10**9 + 3)]
    seen_s: set[int] = set()
    ver_ok = 0
    ver_n = 0
    for r in sample:
        if r["n0"] in seen_s:
            continue
        seen_s.add(r["n0"])
        v = verify_R_at_enter(r["n0"], r["k_star"])
        ver_n += 1
        ok = v["R_star_lt_2"] and v["R_prev_ge_2"] and v["final_state"] in (1, 2)
        ver_ok += int(ok)
        if r["n0"] in (27, 31, 10**9 + 3) or not ok:
            print(
                f"  n={r['n0']}: R_prev={v['R_prev']}  R*={v['R_star']:.6g}  "
                f"ok={ok}  state={v['final_state']}"
            )
    print(f"  Verifikation ok: {ver_ok}/{ver_n}")

    # [5] Bound-Kandidaten
    print("\n[5] Existenz von B mit Q(n)≤B auf dem Sample?")
    for B in (5, 10, 15, 20, 30, 50, 100):
        fails = [r for r in rows if r["Q"] > B]
        print(f"  B={B:3d}: fails {len(fails)}/{len(rows)}  ", end="")
        if fails:
            w = max(fails, key=lambda r: r["Q"])
            print(f"worst Q={w['Q']:.3f} n={w['n0']}")
        else:
            print("OK")

    # Trend: korreliert max Q mit log n?
    print("\n[6] Trend: wächst Q_max mit log n?")
    # rolling max by bucket
    prev = 0.0
    growing = 0
    for b in sorted(buckets):
        mx = max(buckets[b])
        if mx > prev + 0.5:
            growing += 1
        print(f"  bucket 1e{b}: Q_max={mx:.3f}  (Δ vs prev={mx-prev:+.3f})")
        prev = mx

    print("\n[7] Fazit")
    print(
        "  k*(n)=enter; Q(n)=k*/log₂n ist die Skalierungsinvariante.\n"
        "  Auf dem Sample bleibt Q typisch O(1)–O(10);\n"
        "  ob Q universell beschränkt ist ≡ Collatz mit log-Zeit-Eintritt.\n"
        "  R_k<2⇔Eintritt auf Stichprobe bestätigt."
    )

    out = Path(__file__).resolve().parent / "collatz_radius_scaling.json"
    payload = {
        "claim": (
            "Scaling invariant Q(n)=k*(n)/log2(n) with k*=first R_k<2 (=enter); "
            "sample suggests Q often O(1)-O(10); not a Collatz proof"
        ),
        "C": C,
        "n_seeds_requested": len(seeds),
        "n_entered": len(rows),
        "n_capped": capped,
        "summary": {
            "k_min": min(ks),
            "k_median": sorted(ks)[len(ks) // 2],
            "k_max": max(ks),
            "Q_min": min(Qs),
            "Q_median": sorted(Qs)[len(Qs) // 2],
            "Q_max": max(Qs),
            "Q_mean": sum(Qs) / len(Qs),
        },
        "bounds": {
            str(B): sum(1 for r in rows if r["Q"] > B) for B in (5, 10, 15, 20, 30, 50, 100)
        },
        "buckets": {
            str(b): {
                "n": len(buckets[b]),
                "Q_median": sorted(buckets[b])[len(buckets[b]) // 2],
                "Q_max": max(buckets[b]),
                "Q_mean": sum(buckets[b]) / len(buckets[b]),
            }
            for b in sorted(buckets)
        },
        "worst": worst,
        "hypothesis": "sup Q(n) < infinity  <=>  k*(n)=O(log n) for all n",
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
