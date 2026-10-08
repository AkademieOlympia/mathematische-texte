"""
PhaseD_Embedding — Goldene Drift-Quelle (C,L)=(8,8): Perron-Sektor → Primvierling-Label.

Claim-Rand:
* Primvierling-Koordinaten = topologisches Labeling der 2 Perron-Knoten, kein 41D-Basisumbau.
* Die 2×2-Sektor-Matrix wird als Element von Mat₂(ℝ) ≅ Spann{I,X,Y,Z} gelesen (Hurwitz-nah).
* KEIN σ-Tuning; 1D-Normgleichheit ||q||=ρ ist tautologisch und kein Validierungskriterium.
* Nicht: „alle Clog-Komposite sind Hurwitz-Operatoren“.
* Kein Collatz-Beweis.

Datenquelle: live M₃M₅ am aktiven Support; rechter Perron-Support = Lean vCert-Wörter.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from CollatzBlockWeightedTransfer import (
    enumerate_survivor_words,
    get_block_matrix,
    nilpotency_index,
    spectral,
    word_of_residue,
)

# Repo Primvierling helpers (optional CEAB orbit on labels)
_REPO = Path("/Users/thomashoffbauer/Projects/Kepler-Hurrwitz")
if str(_REPO / "src") not in sys.path:
    sys.path.insert(0, str(_REPO / "src"))
from kepler_hurwitz.primvierling import (  # noqa: E402
    quat_norm,
    symmetry_shift_ceab,
)

GOLDEN = {"C": 8.0, "L": 8}
REF_RHO_M3M5 = 0.94921875  # Lean/numeric freeze at (8,8)
LEAN_PERRON_WORDS = ("01111100", "11111100")


def active_kernel(M: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    act = np.where((M.sum(axis=0) + M.sum(axis=1)) > 1e-15)[0]
    if act.size == 0:
        return M.copy(), act
    return M[np.ix_(act, act)].copy(), act


def frobenius(A: np.ndarray, B: np.ndarray) -> float:
    return float(np.sum(A * B))


def orthonormal_basis(n: int) -> dict[str, np.ndarray]:
    """
    Feste 4D-Unterraum-Basis in Mat_n(ℝ) (Frobenius-orthonormal):
      E0 ~ I
      E1 skew (tridiag ±1)
      E2 skew (next superdiag)
      E3 symmetric traceless (diag alternating)
    """
    E0 = np.eye(n) / np.sqrt(n)

    S1 = np.zeros((n, n))
    for i in range(n - 1):
        S1[i, i + 1] = 1.0
        S1[i + 1, i] = -1.0
    n1 = np.linalg.norm(S1, "fro")
    E1 = S1 / n1 if n1 > 0 else S1

    S2 = np.zeros((n, n))
    for i in range(n - 2):
        S2[i, i + 2] = 1.0
        S2[i + 2, i] = -1.0
    # orthogonalize vs E1
    S2 = S2 - frobenius(S2, E1) * E1
    n2 = np.linalg.norm(S2, "fro")
    E2 = S2 / n2 if n2 > 0 else S2

    D = np.diag([1.0 if i % 2 == 0 else -1.0 for i in range(n)])
    D = D - (np.trace(D) / n) * np.eye(n)  # traceless
    for E in (E0, E1, E2):
        D = D - frobenius(D, E) * E
    n3 = np.linalg.norm(D, "fro")
    E3 = D / n3 if n3 > 0 else D

    return {"scalar": E0, "i": E1, "j": E2, "k": E3}


class HurwitzEmbedding:
    """M → q = a + bi + cj + dk via Frobenius projection (no clog/σ)."""

    def __init__(self, matrix_T: np.ndarray):
        self.T = np.asarray(matrix_T, dtype=float)
        if self.T.ndim != 2 or self.T.shape[0] != self.T.shape[1]:
            raise ValueError("matrix_T must be square")
        self.basis = orthonormal_basis(self.T.shape[0])

    def project_to_hq(self) -> dict[str, float]:
        return {name: frobenius(self.T, E) for name, E in self.basis.items()}

    @staticmethod
    def quaternion_norm(coeffs: dict[str, float]) -> float:
        return float(
            np.sqrt(
                coeffs["scalar"] ** 2
                + coeffs["i"] ** 2
                + coeffs["j"] ** 2
                + coeffs["k"] ** 2
            )
        )

    def reconstruct(self, coeffs: dict[str, float] | None = None) -> np.ndarray:
        c = coeffs or self.project_to_hq()
        return sum(c[name] * self.basis[name] for name in ("scalar", "i", "j", "k"))

    def residual_frobenius(self) -> float:
        R = self.T - self.reconstruct()
        return float(np.linalg.norm(R, "fro"))


def build_golden_composites(C: float = 8.0, L: int = 8) -> dict:
    nodes = enumerate_survivor_words(L, C)
    cache: dict[int, np.ndarray] = {}
    for m in (2, 3, 5):
        get_block_matrix(nodes, C, L, m, cache)
    out = {"C": C, "L": L, "n_V": len(nodes), "nodes": nodes, "composites": {}}
    for key, mB, mA in (("M3_M2", 3, 2), ("M3_M5", 3, 5)):
        M = cache[mB] @ cache[mA]
        sp = spectral(M)
        K, act = active_kernel(M)
        out["composites"][key] = {
            "rho": sp["rho"],
            "nilpotency_index": nilpotency_index(M),
            "n_full": int(M.shape[0]),
            "n_active": int(K.shape[0]),
            "active_indices": act.tolist(),
            "active_words": [nodes[i] for i in act],
            "kernel": K.tolist(),
        }
    return out


def embed_report(kernel: np.ndarray, rho_ref: float, label: str) -> dict:
    emb = HurwitzEmbedding(kernel)
    coeffs = emb.project_to_hq()
    qn = emb.quaternion_norm(coeffs)
    return {
        "label": label,
        "n": int(kernel.shape[0]),
        "rho_ref": rho_ref,
        "coeffs": coeffs,
        "quaternion_norm": qn,
        "abs_diff_to_rho": abs(qn - rho_ref),
        "rel_diff_to_rho": abs(qn - rho_ref) / rho_ref if rho_ref > 1e-15 else None,
        "residual_frobenius": emb.residual_frobenius(),
        "match_within_1e-3": abs(qn - rho_ref) < 1e-3,
        "claim": (
            "||q||≈ρ would be a coincidence for this fixed 4D subspace; "
            "mismatch does not license σ-calibration"
        ),
    }


def extract_perron_sector(kernel: np.ndarray, words: list[str]) -> dict:
    """
    Perron-Sektor der aktiven Kernel-Matrix.

    Claim-Rand: 1D-Reduktion auf E_ρ liefert trivial [ρ] — ||q||=ρ wäre tautologisch.
    Der informative Anteil ist Support/Submatrix-Struktur, nicht die Norm-Gleichung.
    """
    vals, vecs = np.linalg.eig(kernel)
    k = int(np.argmax(np.abs(vals)))
    rho = float(np.real(vals[k]))
    vr = np.real(vecs[:, k])
    if vr.min() < 0:
        vr = -vr
    vr = np.maximum(vr, 0.0)
    if vr.max() > 0:
        vr = vr / vr.max()

    valsL, vecsL = np.linalg.eig(kernel.T)
    kL = int(np.argmax(np.abs(valsL)))
    vl = np.real(vecsL[:, kL])
    if vl.min() < 0:
        vl = -vl
    vl = np.maximum(vl, 0.0)
    if vl.max() > 0:
        vl = vl / vl.max()

    r_idx = np.where(vr > 1e-8)[0]
    l_idx = np.where(vl > 1e-8)[0]
    union = sorted(set(r_idx.tolist()) | set(l_idx.tolist()))
    den = float(vl @ vr)
    Pi = np.outer(vr, vl) / den if abs(den) > 1e-15 else np.zeros_like(kernel)
    abs_vals = sorted((float(np.abs(v)) for v in vals), reverse=True)

    return {
        "rho": rho,
        "n_nonzero_eigenvalues_eps": int(sum(1 for a in abs_vals if a > 1e-10)),
        "top_abs_eigenvalues": abs_vals[:8],
        "right_support": [
            {"word": words[i], "v": float(vr[i]), "index": int(i)} for i in r_idx
        ],
        "left_support": [
            {"word": words[i], "v": float(vl[i]), "index": int(i)} for i in l_idx
        ],
        "right_support_matrix": kernel[np.ix_(r_idx, r_idx)].tolist(),
        "union_words": [words[i] for i in union],
        "union_matrix": kernel[np.ix_(union, union)].tolist(),
        "projector_residual": float(np.linalg.norm(kernel @ Pi - rho * Pi)),
        "tautology_warning": (
            "Compressing to 1D Perron ray yields [ρ]; matching ||q|| to ρ is tautological. "
            "Use support structure, not norm equality, as the Phase-D signal."
        ),
        "matches_lean_M3M5_support_words": sorted(words[i] for i in r_idx)
        == list(LEAN_PERRON_WORDS),
    }


def minimal_residue_for_word(word: str, L: int) -> int:
    """Kleinstes Residuum r∈[0,2^L), das das Survivor-Wort erzeugt."""
    for r in range(1 << L):
        if word_of_residue(r, L) == word:
            return r
    raise ValueError(f"no residue for word {word}")


def word_to_primvierling_label(word: str, L: int = 8) -> dict:
    """
    Topologisches Label π=(a,b,c,e) für ein Survivor-Wort.

    Regel (deterministisch, kein Primquadrupel-Anspruch):
      * 8 Bits → vier 2-Bit-Blöcke → (a,b,c,e) ∈ {0,1,2,3}^4
      * zusätzlich: kanonisches Residuum r und quat_norm des Labels
      * CEAB-Orbit auf dem Label (Repo-Symmetrie), rein als Label-Invariante
    """
    if len(word) != L or any(c not in "01" for c in word):
        raise ValueError(f"expected length-{L} bitstring, got {word!r}")
    bits = [int(c) for c in word]
    coords = []
    for k in range(4):
        coords.append(bits[2 * k] * 2 + bits[2 * k + 1])
    a, b, c, e = coords
    pi = (a, b, c, e)
    r = minimal_residue_for_word(word, L)
    return {
        "word": word,
        "residue_min": r,
        "pi": {"a": a, "b": b, "c": c, "e": e, "tuple": list(pi)},
        "quat_norm_label": quat_norm(pi),
        "ceab_shift": list(symmetry_shift_ceab(pi)),
        "ceab_orbit": [list(pi), list(symmetry_shift_ceab(pi))],
        "claim": (
            "Topological Primvierling label from bit pairs; "
            "not a prime quadruplet; not a Collatz proof"
        ),
    }


def mat2_to_hurwitz_coeffs(M: np.ndarray) -> dict:
    """
    Zerlegung M ∈ Mat₂(ℝ) in Spann{I, X, Y, Z} (reelle Pauli-ähnliche Basis):
      I=[[1,0],[0,1]], X=[[0,1],[1,0]], Y=[[0,-1],[1,0]], Z=[[1,0],[0,-1]]
    Dann q ~ (a,b,c,d) mit M = a I + b X + c Y + d Z.
    """
    M = np.asarray(M, dtype=float)
    if M.shape != (2, 2):
        raise ValueError("expected 2×2 sector matrix")
    I = np.array([[1.0, 0.0], [0.0, 1.0]])
    X = np.array([[0.0, 1.0], [1.0, 0.0]])
    Y = np.array([[0.0, -1.0], [1.0, 0.0]])
    Z = np.array([[1.0, 0.0], [0.0, -1.0]])
    # orthonormal up to scale: ⟨A,B⟩_F / ⟨B,B⟩_F
    def coeff(B: np.ndarray) -> float:
        return frobenius(M, B) / frobenius(B, B)

    a, b, c, d = coeff(I), coeff(X), coeff(Y), coeff(Z)
    recon = a * I + b * X + c * Y + d * Z
    qnorm = float(np.sqrt(a * a + b * b + c * c + d * d))
    return {
        "a_scalar": a,
        "b_i": b,
        "c_j": c,
        "d_k": d,
        "quaternion_euclid_norm": qnorm,
        "spectral_radius_2x2": float(np.max(np.abs(np.linalg.eigvals(M)))),
        "frobenius_norm_2x2": float(np.linalg.norm(M, "fro")),
        "reconstruction_error": float(np.linalg.norm(M - recon, "fro")),
        "basis": "I,X,Y,Z on Mat2(R)",
    }


def map_perron_sector_to_primvierling(sector: dict, L: int = 8) -> dict:
    """
    Exaktes Mapping {v1,v2} → {π1,π2} + Hurwitz-Lesart der 2×2-Sektor-Matrix.
    """
    right = sector["right_support"]
    words = [row["word"] for row in right]
    if sorted(words) != sorted(LEAN_PERRON_WORDS):
        raise ValueError(f"unexpected Perron words {words}")
    # stabile Reihenfolge: Lean-Freeze-Reihenfolge
    ordered = list(LEAN_PERRON_WORDS)
    labels = {w: word_to_primvierling_label(w, L) for w in ordered}
    # 2×2 in derselben Reihenfolge
    idx = {row["word"]: i for i, row in enumerate(right)}
    Mraw = np.array(sector["right_support_matrix"], float)
    perm = [idx[w] for w in ordered]
    M = Mraw[np.ix_(perm, perm)]
    hq = mat2_to_hurwitz_coeffs(M)
    return {
        "claim": (
            "Golden drift source: Perron 2-node sector labeled by Primvierling bit-pair "
            "coords; Mat2≅span{I,X,Y,Z}; not a global H_q theorem"
        ),
        "nodes": ordered,
        "labels": labels,
        "sector_matrix_ordered": M.tolist(),
        "hurwitz_mat2": hq,
        "rho_sector": sector["rho"],
        "note": (
            "quaternion_euclid_norm is the Euclidean norm of (a,b,c,d), "
            "not required to equal ρ; spectral_radius_2x2 should equal ρ"
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--mode",
        choices=("frobenius", "perron", "labeling", "both"),
        default="labeling",
        help="labeling = Perron sector + Primvierling map (default Phase-D handoff)",
    )
    ap.add_argument(
        "--out",
        type=str,
        default=str(
            Path("/Users/thomashoffbauer/Projects/Kepler-Hurrwitz/docs/exports")
        ),
    )
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print(f"PHASE D EMBEDDING  golden cell (8,8)  mode={args.mode}")
    print("=" * 72)

    data = build_golden_composites(**GOLDEN)
    kernel_path = out_dir / "h7_golden_cell_kernel_M3M5_C8L8.json"
    # export M3M5 active kernel (real data; not fictional 12x12)
    export_kernel = {
        "claim": (
            "Active support kernel of M3M5 at (C,L)=(8,8); "
            "not a 12x12 fiction; not a Collatz proof"
        ),
        "C": data["C"],
        "L": data["L"],
        "n_V": data["n_V"],
        "composite": "M3_M5",
        "rho": data["composites"]["M3_M5"]["rho"],
        "n_active": data["composites"]["M3_M5"]["n_active"],
        "active_words": data["composites"]["M3_M5"]["active_words"],
        "matrix": data["composites"]["M3_M5"]["kernel"],
    }
    kernel_path.write_text(json.dumps(export_kernel, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {kernel_path}  n_active={export_kernel['n_active']}")

    reports = []
    perron = None
    if args.mode in ("frobenius", "both"):
        for key in ("M3_M5", "M3_M2"):
            comp = data["composites"][key]
            K = np.array(comp["kernel"], float)
            rho = float(comp["rho"])
            if K.size == 0:
                continue
            rep = embed_report(K, rho if rho > 1e-15 else 0.0, key)
            rep["nilpotency_index"] = comp["nilpotency_index"]
            reports.append(rep)
            print(
                f"\n{key} Frobenius: n={rep['n']} ρ={rho:.6f}  "
                f"||q||={rep['quaternion_norm']:.6f}  |Δ|={rep['abs_diff_to_rho']:.6f}"
            )

    labeling = None
    if args.mode in ("perron", "labeling", "both"):
        comp = data["composites"]["M3_M5"]
        K = np.array(comp["kernel"], float)
        words = comp["active_words"]
        perron = extract_perron_sector(K, words)
        print("\nM3_M5 Perron sector:")
        print(f"  ρ={perron['rho']:.6f}  nonzero_eigs≈{perron['n_nonzero_eigenvalues_eps']}")
        print(f"  right support: {[w['word'] for w in perron['right_support']]}")
        print(f"  left support:  {[w['word'] for w in perron['left_support']]}")
        print(f"  right 2×2:\n{np.array(perron['right_support_matrix'])}")
        print(f"  Lean support match: {perron['matches_lean_M3M5_support_words']}")
        print(f"  NOTE: {perron['tautology_warning']}")
        perron_path = out_dir / "h7_phaseD_perron_sector_M3M5_C8L8.json"
        perron_path.write_text(
            json.dumps(
                {
                    "claim": (
                        "Perron sector of active M3M5 at (8,8); "
                        "1D ρ-match is tautological; not a Collatz proof"
                    ),
                    "golden": GOLDEN,
                    "sector": perron,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        print(f"wrote {perron_path}")

    if args.mode in ("labeling", "both"):
        assert perron is not None
        labeling = map_perron_sector_to_primvierling(perron, L=int(GOLDEN["L"]))
        print("\nPrimvierling labeling (Perron → π):")
        for w in labeling["nodes"]:
            lab = labeling["labels"][w]
            print(
                f"  {w}  →  π={tuple(lab['pi']['tuple'])}  "
                f"r={lab['residue_min']}  ||π||²={lab['quat_norm_label']}  "
                f"CEAB→{tuple(lab['ceab_shift'])}"
            )
        hq = labeling["hurwitz_mat2"]
        print(
            f"  Mat2 Hurwitz coeffs: a={hq['a_scalar']:.6f} b={hq['b_i']:.6f} "
            f"c={hq['c_j']:.6f} d={hq['d_k']:.6f}"
        )
        print(
            f"  ρ(2×2)={hq['spectral_radius_2x2']:.6f}  "
            f"||(a,b,c,d)||={hq['quaternion_euclid_norm']:.6f}  "
            f"recon_err={hq['reconstruction_error']:.2e}"
        )
        lab_path = out_dir / "h7_phaseD_primvierling_labeling_C8L8.json"
        lab_path.write_text(
            json.dumps({"golden": GOLDEN, **labeling}, indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"wrote {lab_path}")

    payload = {
        "claim": (
            "Phase D golden drift source: Perron 2-node sector + Primvierling labeling; "
            "not a Collatz proof; not induction"
        ),
        "golden": GOLDEN,
        "reference_rho_M3M5": REF_RHO_M3M5,
        "reject_sigma_calibration": True,
        "mode": args.mode,
        "embeddings": reports,
        "perron_sector_M3M5": perron,
        "primvierling_labeling": labeling,
        "kernel_export": str(kernel_path),
    }
    out = out_dir / "h7_phaseD_embedding_C8L8.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
