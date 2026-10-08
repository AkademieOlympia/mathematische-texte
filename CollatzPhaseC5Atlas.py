"""
CollatzPhaseC5Atlas — Phase C.5: (C,L)-Raster für M₃M₂ / M₃M₅.

Default: vollständiger Durchlauf (keine interaktive Pause).
Bei jedem neuen IsResonant-Punkt (≠ bekannte Chirurgie-Insel):
  * SINGULARITY-WARNUNG in Log + JSON
  * Auto-Pfad-Dekonstruktion der Selbstschleifen
Optional: --pause-on-resonant für manuelles Anhalten.

Golden anchors (Lean):
  (8,8)  CompositePhaseCAuditReport / CompositeGapAudit
  (8,11) SurgicalAudit_L11 (chirurgisch; Baseline resonant)

Epistemik: numerischer Atlas; kein Collatz-Beweis; keine Induktion.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

from CollatzBlockWeightedTransfer import (
    enumerate_survivor_words,
    get_block_matrix,
    nilpotency_index,
    spectral,
)
from CollatzL11BreakAnalysis import extract_path_details

# Lean-verankerte Zellen
GOLDEN_CORE = {(8, 8)}
KNOWN_SURGERY = {(8, 11)}  # Baseline resonant; ForbiddenSet Not-Patch
PAIRS = [(3, 2), (3, 5)]


def classify_composite(rho: float, nil_k: int | None) -> str:
    if rho >= 1.0 - 1e-12:
        return "IsResonant"
    if nil_k is not None or rho < 1e-12:
        return "IsNilpotent"
    return "IsContractive"  # 0 < ρ < 1


def audit_cell(C: float, L: int) -> dict:
    t0 = time.time()
    nodes = enumerate_survivor_words(L, C)
    cell: dict = {
        "C": int(C) if float(C).is_integer() else C,
        "L": L,
        "n_V": len(nodes),
        "IsEmpty": len(nodes) == 0,
        "IsNilpotent": False,
        "IsResonant": False,
        "IsContractive": False,
        "composites": {},
        "elapsed_s": 0.0,
    }
    Ci, Li = int(round(C)), L
    cell["is_golden_core"] = (Ci, Li) in GOLDEN_CORE
    cell["is_known_surgery"] = (Ci, Li) in KNOWN_SURGERY

    if not nodes:
        cell["status"] = "IsEmpty"
        cell["elapsed_s"] = time.time() - t0
        return cell

    cache: dict[int, np.ndarray] = {}
    any_res = False
    all_nil = True
    any_contractive = False

    for mB, mA in PAIRS:
        key = f"M{mB}_M{mA}"
        MB = get_block_matrix(nodes, C, L, mB, cache)
        MA = get_block_matrix(nodes, C, L, mA, cache)
        M = MB @ MA
        sp = spectral(M)
        nil_k = nilpotency_index(M)
        tag = classify_composite(sp["rho"], nil_k)
        diag = np.diag(M)
        loops = [
            {"word": nodes[i], "weight": float(diag[i])}
            for i in np.where(diag > 1e-12)[0]
        ]
        cell["composites"][key] = {
            "rho": sp["rho"],
            "nilpotency_index": nil_k,
            "nnz": int(np.count_nonzero(M)),
            "tag": tag,
            "self_loops": loops,
        }
        if tag == "IsResonant":
            any_res = True
            all_nil = False
        elif tag == "IsContractive":
            all_nil = False
            any_contractive = True
        elif tag != "IsNilpotent":
            all_nil = False

    cell["IsResonant"] = any_res
    cell["IsNilpotent"] = all_nil and not any_res
    cell["IsContractive"] = (not any_res) and any_contractive and not all_nil
    if any_res:
        cell["status"] = "IsResonant"
    elif cell["IsEmpty"]:
        cell["status"] = "IsEmpty"
    elif cell["IsNilpotent"]:
        cell["status"] = "IsNilpotent"
    else:
        cell["status"] = "IsContractive"

    cell["elapsed_s"] = round(time.time() - t0, 3)
    return cell


def deconstruct_cell_loops(C: float, L: int, cell: dict) -> dict:
    """Pfad-Dekonstruktion aller Selbstschleifen in resonanten Kompositen."""
    nodes = enumerate_survivor_words(L, C)
    out: dict = {"C": C, "L": L, "composites": {}}
    for key, comp in cell["composites"].items():
        if comp.get("tag") != "IsResonant":
            continue
        mB, mA = (3, 2) if key == "M3_M2" else (3, 5)
        words = [lp["word"] for lp in comp.get("self_loops", [])]
        if not words:
            # fallback: top Perron support not available; skip
            continue
        block = {}
        for w in words:
            block[w] = extract_path_details(w, nodes, C, L, mA, mB)
        out["composites"][key] = block
    return out


def run_atlas(
    C_min: int,
    C_max: int,
    L_min: int,
    L_max: int,
    out_dir: Path,
    pause_on_resonant: bool,
    deconstruct: bool,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    cells: list[dict] = []
    warnings: list[dict] = []
    deconstructions: list[dict] = []

    print("=" * 72)
    print(
        f"PHASE C.5 ATLAS  C={C_min}..{C_max}  L={L_min}..{L_max}  "
        f"pause={pause_on_resonant}  deconstruct={deconstruct}"
    )
    print("=" * 72)
    print(
        f"{'C':>3} {'L':>3} {'|V|':>6} {'status':>14} "
        f"{'ρ32':>8} {'ρ35':>8} {'note'}"
    )

    for L in range(L_min, L_max + 1):
        for C in range(C_min, C_max + 1):
            cell = audit_cell(float(C), L)
            cells.append(cell)
            c32 = cell.get("composites", {}).get("M3_M2", {})
            c35 = cell.get("composites", {}).get("M3_M5", {})
            r32 = c32.get("rho")
            r35 = c35.get("rho")
            note = ""
            if cell.get("is_golden_core"):
                note = "GOLDEN(8,8)"
            elif cell.get("is_known_surgery"):
                note = "KNOWN_SURGERY"

            if cell["IsResonant"]:
                is_new = (C, L) not in KNOWN_SURGERY
                warn = {
                    "C": C,
                    "L": L,
                    "n_V": cell["n_V"],
                    "rho_M3_M2": r32,
                    "rho_M3_M5": r35,
                    "self_loops_M3_M2": c32.get("self_loops", []),
                    "self_loops_M3_M5": c35.get("self_loops", []),
                    "known_surgery": not is_new,
                    "message": (
                        "KNOWN surgical island (8,11)"
                        if not is_new
                        else "NEW resonant cell — candidate window-lifting"
                    ),
                }
                warnings.append(warn)
                note = ("⚠ NEW " if is_new else "⚠ ") + note
                print(
                    f"{C:3d} {L:3d} {cell['n_V']:6d} {cell['status']:>14} "
                    f"{(r32 if r32 is not None else float('nan')):8.4f} "
                    f"{(r35 if r35 is not None else float('nan')):8.4f} {note}"
                )
                print(f"    SINGULARITY: {warn['message']}", flush=True)

                if deconstruct:
                    try:
                        dec = deconstruct_cell_loops(float(C), L, cell)
                        deconstructions.append(dec)
                        dec_path = (
                            out_dir
                            / f"h7_phaseC5_deconstruct_C{C}_L{L}.json"
                        )
                        dec_path.write_text(
                            json.dumps(dec, indent=2) + "\n", encoding="utf-8"
                        )
                        print(f"    wrote {dec_path.name}", flush=True)
                    except Exception as e:
                        warn["deconstruct_error"] = str(e)
                        print(f"    deconstruct FAILED: {e}", flush=True)

                if pause_on_resonant and is_new:
                    print(
                        "    PAUSE (--pause-on-resonant). Enter to continue, q to abort:",
                        flush=True,
                    )
                    line = sys.stdin.readline()
                    if line.strip().lower() in {"q", "quit", "abort"}:
                        print("aborted by user")
                        break
            else:
                print(
                    f"{C:3d} {L:3d} {cell['n_V']:6d} {cell['status']:>14} "
                    f"{(r32 if r32 is not None else float('nan')):8.4f} "
                    f"{(r35 if r35 is not None else float('nan')):8.4f} {note}",
                    flush=True,
                )
        else:
            continue
        break  # aborted inner

    # summary stats
    n_res = sum(1 for c in cells if c.get("IsResonant"))
    n_nil = sum(1 for c in cells if c.get("IsNilpotent"))
    n_empty = sum(1 for c in cells if c.get("IsEmpty"))
    n_con = sum(1 for c in cells if c.get("IsContractive"))
    new_res = [w for w in warnings if not w.get("known_surgery")]

    payload = {
        "claim": (
            "Phase C.5 (C,L) atlas for M3M2/M3M5; singularity warnings + "
            "optional deconstruct; not a Collatz proof; not induction"
        ),
        "golden_lean_anchors": {
            "core_C8_L8": "CompositePhaseCAuditReport / CompositeGapAudit",
            "surgery_C8_L11": "SurgicalAudit_L11",
        },
        "grid": {
            "C": [C_min, C_max],
            "L": [L_min, L_max],
            "pairs": ["M3_M2", "M3_M5"],
        },
        "mode": {
            "pause_on_resonant": pause_on_resonant,
            "auto_deconstruct": deconstruct,
        },
        "counts": {
            "cells": len(cells),
            "IsResonant": n_res,
            "IsNilpotent": n_nil,
            "IsContractive": n_con,
            "IsEmpty": n_empty,
            "new_resonant_warnings": len(new_res),
        },
        "singularity_warnings": warnings,
        "new_resonant_cells": new_res,
        "cells": cells,
    }

    out = out_dir / (
        f"h7_phaseC5_atlas_C{C_min}-{C_max}_L{L_min}-{L_max}.json"
    )
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    # compact matrix summary
    matrix = []
    for L in range(L_min, L_max + 1):
        row = []
        for C in range(C_min, C_max + 1):
            hit = next(
                (c for c in cells if c["C"] == C and c["L"] == L), None
            )
            row.append(hit["status"] if hit else "?")
        matrix.append({"L": L, "by_C": row})
    summary = {
        "C_range": list(range(C_min, C_max + 1)),
        "status_matrix_rows_by_L": matrix,
        "counts": payload["counts"],
        "new_resonant_cells": new_res,
    }
    sum_path = out_dir / (
        f"h7_phaseC5_atlas_summary_C{C_min}-{C_max}_L{L_min}-{L_max}.json"
    )
    sum_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    print(f"wrote {sum_path}")
    print(
        f"counts: res={n_res} nil={n_nil} con={n_con} empty={n_empty} "
        f"new_res={len(new_res)}"
    )
    return payload


def main() -> int:
    ap = argparse.ArgumentParser(description="Phase C.5 (C,L) atlas")
    ap.add_argument("--C-min", type=int, default=8)
    ap.add_argument("--C-max", type=int, default=16)
    ap.add_argument("--L-min", type=int, default=8)
    ap.add_argument("--L-max", type=int, default=20)
    ap.add_argument(
        "--pause-on-resonant",
        action="store_true",
        help="Interactive pause on NEW resonant cells (default: full run)",
    )
    ap.add_argument(
        "--no-deconstruct",
        action="store_true",
        help="Skip auto path deconstruction on resonant cells",
    )
    ap.add_argument(
        "--out",
        type=str,
        default=str(
            Path("/Users/thomashoffbauer/Projects/Kepler-Hurrwitz/docs/exports")
        ),
    )
    args = ap.parse_args()
    run_atlas(
        args.C_min,
        args.C_max,
        args.L_min,
        args.L_max,
        Path(args.out),
        pause_on_resonant=args.pause_on_resonant,
        deconstruct=not args.no_deconstruct,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
