"""
CollatzScaleGenerator — Skalierungs-Audit um die Basis (C,L)=(8,8).

Modi:
  scale_L  — festes C=8, L ∈ {8,9,10,11}   (primär: Spektrale Stabilität)
  scale_C  — festes L=8, C ∈ {8,9,10,11,12} (Phasenraum / Azyklizität)
  grid     — volles Gitter L×C

Pro Zelle:
  |V_L|, ρ(M_w), Komposita M3∘M2 und M3∘M5 (falls m verfügbar),
  erste Blocklänge m∈{2..6} mit ρ(M_{w,m})<1.

Epistemik: Kartierung; kein Collatz-Beweis; Lean-Induktion separat.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from CollatzBlockWeightedTransfer import (
    build_base_Mw,
    block_matrix_walk_filter,
    enumerate_survivor_words,
    get_block_matrix,
    nilpotency_index,
    rational_subinvariant_from_perron,
    spectral,
)


def audit_cell(C: float, L: int, m_scan: list[int], composite_pairs: list[tuple[int, int]]) -> dict:
    nodes = enumerate_survivor_words(L, C)
    n = len(nodes)
    cell: dict = {
        "C": C,
        "L": L,
        "n_V": n,
        "empty": n == 0,
    }
    if n == 0:
        cell["status"] = "empty_survivor_set"
        return cell

    cache: dict[int, np.ndarray] = {}
    Mw = get_block_matrix(nodes, C, L, 1, cache)
    sp0 = spectral(Mw)
    cell["rho_Mw"] = sp0["rho"]
    cell["nnz_Mw"] = sp0["nnz"]
    cell["Mw_nilpotent"] = nilpotency_index(Mw) is not None

    # block scan
    blocks = {}
    first_contractive_m = None
    for m in m_scan:
        try:
            Mb = get_block_matrix(nodes, C, L, m, cache)
        except Exception as e:
            blocks[str(m)] = {"error": str(e)}
            continue
        sp = spectral(Mb)
        nil_k = nilpotency_index(Mb)
        rec = {
            "rho": sp["rho"],
            "rho_lt_1": sp["rho"] < 1.0 - 1e-12,
            "nilpotency_index": nil_k,
            "nnz": sp["nnz"],
        }
        blocks[str(m)] = rec
        if first_contractive_m is None and rec["rho_lt_1"]:
            first_contractive_m = m
    cell["blocks"] = blocks
    cell["m_min_rho_lt_1"] = first_contractive_m

    # composites
    comps = {}
    for mB, mA in composite_pairs:
        key = f"M{mB}_M{mA}"
        if mB not in cache and mB != 1:
            try:
                get_block_matrix(nodes, C, L, mB, cache)
            except Exception as e:
                comps[key] = {"error": str(e)}
                continue
        if mA not in cache and mA != 1:
            try:
                get_block_matrix(nodes, C, L, mA, cache)
            except Exception as e:
                comps[key] = {"error": str(e)}
                continue
        if mB not in cache or mA not in cache:
            # ensure loaded
            get_block_matrix(nodes, C, L, mB, cache)
            get_block_matrix(nodes, C, L, mA, cache)
        Mcomp = cache[mB] @ cache[mA]
        sp = spectral(Mcomp)
        nil_k = nilpotency_index(Mcomp)
        cert = None
        theta = None
        if sp["rho"] < 1e-12:
            theta = 0.5
            cert = {"regime": "nilpotent", "theta": 0.5, "theta_lt_1": True}
        else:
            cert = rational_subinvariant_from_perron(
                Mcomp, np.array(sp["v"], float), sp["rho"], margin=1.01
            )
            if cert and cert.get("theta_lt_1"):
                theta = cert["theta"]
        comps[key] = {
            "mB": mB,
            "mA": mA,
            "rho": sp["rho"],
            "rho_lt_1": sp["rho"] < 1.0 - 1e-12,
            "nilpotency_index": nil_k,
            "theta_contractive": theta,
            "subinvariant": (
                {k: cert[k] for k in cert if k != "v_int"} if cert else None
            ),
        }
    cell["composites"] = comps

    # summary flags
    cell["status"] = "ok"
    cell["has_composite_gap"] = any(
        c.get("rho_lt_1") for c in comps.values() if isinstance(c, dict) and "rho" in c
    )
    cell["M3_M2_rho"] = comps.get("M3_M2", {}).get("rho")
    cell["M3_M5_rho"] = comps.get("M3_M5", {}).get("rho")
    cell["M3_M5_theta"] = comps.get("M3_M5", {}).get("theta_contractive")
    return cell


def run_scale(
    mode: str,
    L_range: list[int],
    C_range: list[float],
    out_dir: Path,
) -> dict:
    m_scan = [2, 3, 4, 5, 6]
    composite_pairs = [(3, 2), (3, 5)]

    if mode == "scale_L":
        pairs = [(C_range[0], L) for L in L_range]
        title = f"scale_L C={C_range[0]} L={L_range}"
    elif mode == "scale_C":
        pairs = [(C, L_range[0]) for C in C_range]
        title = f"scale_C L={L_range[0]} C={C_range}"
    else:
        pairs = [(C, L) for L in L_range for C in C_range]
        title = f"grid L={L_range} C={C_range}"

    print("=" * 72)
    print(f"SCALE AUDIT — {title}")
    print("=" * 72)

    cells = []
    print(
        f"{'C':>5} {'L':>3} {'|V|':>5} {'ρ(Mw)':>9} {'m*':>3} "
        f"{'ρ32':>9} {'ρ35':>9} {'θ35':>8} {'gap?':>5}"
    )
    for C, L in pairs:
        cell = audit_cell(C, L, m_scan, composite_pairs)
        cells.append(cell)
        if cell["empty"]:
            print(f"{C:5.1f} {L:3d} {'∅':>5}")
            continue
        th35 = cell.get("M3_M5_theta")
        th35_s = f"{th35:.4f}" if th35 is not None else "—"
        r32 = cell["M3_M2_rho"] if cell["M3_M2_rho"] is not None else float("nan")
        r35 = cell["M3_M5_rho"] if cell["M3_M5_rho"] is not None else float("nan")
        print(
            f"{C:5.1f} {L:3d} {cell['n_V']:5d} {cell['rho_Mw']:9.4f} "
            f"{str(cell['m_min_rho_lt_1'] or '—'):>3} "
            f"{r32:9.4f} {r35:9.4f} {th35_s:>8} "
            f"{str(cell['has_composite_gap']):>5}"
        )

    # analysis hooks
    nonempty = [c for c in cells if not c.get("empty")]
    analysis = {
        "n_cells": len(cells),
        "n_nonempty": len(nonempty),
        "n_with_composite_gap": sum(1 for c in nonempty if c.get("has_composite_gap")),
        "M3_M2_stays_lt_1": all(
            (c.get("M3_M2_rho") is not None and c["M3_M2_rho"] < 1)
            for c in nonempty
            if c.get("M3_M2_rho") is not None
        ),
        "M3_M5_stays_lt_1": all(
            (c.get("M3_M5_rho") is not None and c["M3_M5_rho"] < 1)
            for c in nonempty
            if c.get("M3_M5_rho") is not None
        ),
        "theta_M3_M5_by_cell": [
            {"C": c["C"], "L": c["L"], "theta": c.get("M3_M5_theta"), "rho": c.get("M3_M5_rho")}
            for c in nonempty
        ],
    }

    payload = {
        "claim": (
            "Scaling audit around (C,L)=(8,8); composite gap stability vs C/L; "
            "not a Collatz proof; not a Lean induction"
        ),
        "mode": mode,
        "L_range": L_range,
        "C_range": C_range,
        "m_scan": m_scan,
        "composite_pairs": [{"mB": a, "mA": b} for a, b in composite_pairs],
        "cells": cells,
        "analysis": analysis,
    }

    stem = {
        "scale_L": "h7_scale_L_C8",
        "scale_C": "h7_scale_C_L8",
        "grid": "h7_scale_grid",
    }[mode]
    out = out_dir / f"{stem}.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    print("analysis:", json.dumps(analysis, indent=2))
    return payload


def main() -> int:
    ap = argparse.ArgumentParser(description="Collatz scaling audit around (8,8)")
    ap.add_argument(
        "--mode",
        choices=["scale_L", "scale_C", "grid"],
        default="scale_L",
        help="Primary: scale_L (C=8, vary L). Secondary: scale_C. Full: grid.",
    )
    ap.add_argument("--L", type=str, default="8,9,10,11")
    ap.add_argument("--C", type=str, default="8,9,10,11,12")
    ap.add_argument(
        "--out",
        type=str,
        default=str(Path("/Users/thomashoffbauer/Projects/Kepler-Hurrwitz/docs/exports")),
    )
    args = ap.parse_args()
    L_range = [int(x) for x in args.L.split(",") if x.strip()]
    C_range = [float(x) for x in args.C.split(",") if x.strip()]
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.mode == "scale_L":
        # festes C=8
        C_range = [8.0]
    elif args.mode == "scale_C":
        L_range = [8]

    run_scale(args.mode, L_range, C_range, out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
