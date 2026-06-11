# bm_qhe_model_v2.py
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, List, Set

import json
import numpy as np
import math

try:
    from scipy.sparse import coo_matrix, csr_matrix, diags
    from scipy.sparse.linalg import eigsh
except ImportError as e:
    raise ImportError("Dieses Modul benötigt scipy (sparse + eigsh).") from e

MODEL_VERSION = "0.4.0"

def _as_path(p: str | Path) -> Path:
    return p if isinstance(p, Path) else Path(p)

def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def _json_dump(obj: Dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")

def _json_load(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))

def _invalidate(d: Dict[str, Any], keys: List[str]) -> None:
    for k in keys:
        d.pop(k, None)

@dataclass
class SIConstants:
    """SI: h, e, c sind exakt; m_e ist empfohlen (CODATA)."""
    h: float = 6.62607015e-34
    c: float = 299792458.0
    e: float = 1.602176634e-19
    eps0: float = 8.8541878128e-12
    m_e: float = 9.1093837139e-31

    @property
    def hbar(self) -> float:
        return self.h / (2.0 * np.pi)

    @property
    def compton_frequency(self) -> float:
        return self.m_e * self.c**2 / self.h

    @property
    def compton_omega(self) -> float:
        return 2.0 * np.pi * self.compton_frequency

    @property
    def compton_tau_bar(self) -> float:
        return 1.0 / self.compton_omega

    @property
    def flux_quantum_h_over_e(self) -> float:
        return self.h / self.e

    @property
    def von_klitzing(self) -> float:
        return self.h / (self.e**2)

    @property
    def conductance_quantum(self) -> float:
        return (self.e**2) / self.h


@dataclass
class BambergQHEModel:
    """Bamberg-Kugel + Walter/Rest + U(1)-Gauge (QHE/Morley), Modus (C)."""

    vertices: np.ndarray
    faces: np.ndarray
    face_channel: np.ndarray
    radius_m: float = 1.0

    # Kepler-Normierung (optional): sqrt(12) für Flächen, sqrt(18) für Volumina
    use_kepler: bool = False
    kepler_area_root: float = 12.0**0.5
    kepler_volume_root: float = 18.0**0.5

    si: SIConstants = field(default_factory=SIConstants)
    tau0_s: Optional[float] = None

    N_phi_global: int = 0
    delta_N_phi_walter: int = 0
    theta_morley: float = 0.0

    w_walter: float = 0.1
    w_rest: float = 0.9

    eig_k: int = 80
    store_eigenvectors: bool = True

    _cache: Dict[str, Any] = field(default_factory=dict, init=False, repr=False)

    @property
    def nV(self) -> int:
        return int(self.vertices.shape[0])

    @property
    def nF(self) -> int:
        return int(self.faces.shape[0])

    @property
    def tau(self) -> float:
        return self.si.compton_tau_bar if self.tau0_s is None else float(self.tau0_s)

    @tau.setter
    def tau(self, v: float) -> None:
        self.tau0_s = float(v)
        _invalidate(self._cache, ["energy_scale"])

    @property
    def weights(self) -> Tuple[float, float]:
        return (float(self.w_walter), float(self.w_rest))

    @weights.setter
    def weights(self, wr: Tuple[float, float]) -> None:
        wW, wR = float(wr[0]), float(wr[1])
        s = wW + wR
        if s <= 0:
            raise ValueError("w_walter + w_rest muss > 0 sein.")
        self.w_walter = wW / s
        self.w_rest = wR / s
        _invalidate(self._cache, ["L_mix", "spectrum_A", "edge_phases"])

    @property
    def flux(self) -> Tuple[int, int]:
        return (int(self.N_phi_global), int(self.delta_N_phi_walter))

    @flux.setter
    def flux(self, v: Tuple[int, int]) -> None:
        self.N_phi_global = int(v[0])
        self.delta_N_phi_walter = int(v[1])
        _invalidate(self._cache, ["phase_faces", "edge_phases", "L_A", "qhe_summary"])

    def edges(self) -> np.ndarray:
        if "edges" in self._cache:
            return self._cache["edges"]
        edge_set: Set[Tuple[int, int]] = set()
        for a, b, c in self.faces.astype(int):
            edge_set.add(tuple(sorted((a, b))))
            edge_set.add(tuple(sorted((b, c))))
            edge_set.add(tuple(sorted((c, a))))
        E = np.array(sorted(edge_set), dtype=int)
        self._cache["edges"] = E
        return E

    def edge_index(self) -> Dict[Tuple[int, int], int]:
        if "edge_index" in self._cache:
            return self._cache["edge_index"]
        E = self.edges()
        idx = {(int(u), int(v)): i for i, (u, v) in enumerate(E)}
        self._cache["edge_index"] = idx
        return idx

    def walter_boundary_vertices(self) -> np.ndarray:
        if "boundary_vertices" in self._cache:
            return self._cache["boundary_vertices"]
        F = self.faces.astype(int)
        ch = self.face_channel.astype(int)
        nV = self.nV
        seenW = np.zeros(nV, dtype=bool)
        seenR = np.zeros(nV, dtype=bool)
        for (a, b, c), isW in zip(F, ch):
            if isW == 1:
                seenW[a] = seenW[b] = seenW[c] = True
            else:
                seenR[a] = seenR[b] = seenR[c] = True
        boundary = np.where(seenW & seenR)[0].astype(int)
        self._cache["boundary_vertices"] = boundary
        return boundary

    def phase_faces(self) -> np.ndarray:
        if "phase_faces" in self._cache:
            return self._cache["phase_faces"]
        maskW = (self.face_channel.astype(int) == 1)
        nW = int(np.sum(maskW))
        nF = self.nF

        phi_base = 2.0 * np.pi * int(self.N_phi_global) / max(nF, 1)
        phi_extra = 2.0 * np.pi * int(self.delta_N_phi_walter) / max(nW, 1) if nW > 0 else 0.0

        phi = np.full(nF, phi_base, dtype=float)
        phi[maskW] += phi_extra

        self._cache["phase_faces"] = phi
        return phi

    @staticmethod
    def _cotangent(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
        ba = b - a
        ca = c - a
        cross = np.linalg.norm(np.cross(ba, ca))
        dot = float(np.dot(ba, ca))
        if cross <= 0:
            return 0.0
        return dot / cross

    def _assemble_sparse_L_M(self, face_mask: Optional[np.ndarray] = None,
                             phases: Optional[Dict[Tuple[int, int], complex]] = None):
        nV = self.nV
        if face_mask is None:
            face_mask = np.ones(self.nF, dtype=bool)
        face_mask = face_mask.astype(bool)

        I: List[int] = []
        J: List[int] = []
        data: List[complex] = []
        deg = np.zeros(nV, dtype=complex if phases is not None else float)
        mdiag = np.zeros(nV, dtype=float)

        def add_edge(u: int, v: int, w: float) -> None:
            uu, vv = (u, v) if u < v else (v, u)
            phase = 1.0 + 0j
            if phases is not None and (uu, vv) in phases:
                phase = phases[(uu, vv)]
            I.extend([u, v])
            J.extend([v, u])
            data.extend([-w * phase, -w * np.conjugate(phase)])
            deg[u] += w
            deg[v] += w

        V = self.vertices
        F = self.faces.astype(int)
        for idx, (i, j, k) in enumerate(F):
            if not face_mask[idx]:
                continue
            vi, vj, vk = V[i], V[j], V[k]
            area = 0.5 * np.linalg.norm(np.cross(vj - vi, vk - vi))
            mdiag[i] += area / 3.0
            mdiag[j] += area / 3.0
            mdiag[k] += area / 3.0

            cot_i = self._cotangent(vi, vj, vk)
            cot_j = self._cotangent(vj, vk, vi)
            cot_k = self._cotangent(vk, vi, vj)
            add_edge(j, k, 0.5 * cot_i)
            add_edge(k, i, 0.5 * cot_j)
            add_edge(i, j, 0.5 * cot_k)

        for u in range(nV):
            I.append(u)
            J.append(u)
            data.append(deg[u])

        L = coo_matrix((np.array(data), (np.array(I), np.array(J))), shape=(nV, nV)).tocsr()
        return L, mdiag

    def L_full(self):
        if "L_full" in self._cache:
            return self._cache["L_full"]
        L, m = self._assemble_sparse_L_M()
        self._cache["L_full"] = (L, m)
        return L, m

    def solve_edge_thetas_from_face_flux(self, ridge: float = 1e-6) -> Dict[str, Any]:
        cache_key = f"theta_e_r{ridge:.0e}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        E = self.edges()
        nE = E.shape[0]
        nF = self.nF
        eidx = self.edge_index()

        F = self.faces.astype(int)
        b = self.phase_faces().astype(float)
        b = (b + np.pi) % (2 * np.pi) - np.pi

        rows, cols, vals = [], [], []
        for fi, (a, bV, c) in enumerate(F):
            for (u, v) in [(a, bV), (bV, c), (c, a)]:
                uu, vv = (u, v) if u < v else (v, u)
                j = eidx[(uu, vv)]
                sgn = +1.0 if (u < v) else -1.0
                rows.append(fi)
                cols.append(j)
                vals.append(sgn)

        A = coo_matrix((np.array(vals), (np.array(rows), np.array(cols))), shape=(nF, nE)).tocsr()
        AtA = (A.T @ A).tocsr()
        rhs = A.T @ b
        AtA_reg = AtA + ridge * diags(np.ones(nE), 0, format="csr")

        try:
            from scipy.sparse.linalg import spsolve
            theta = spsolve(AtA_reg, rhs)
        except Exception:
            from scipy.sparse.linalg import cg
            theta, info = cg(AtA_reg, rhs, tol=1e-10, maxiter=2000)
            if info != 0:
                raise RuntimeError(f"cg failed with info={info}")

        theta = np.array(theta, dtype=float)
        r = (A @ theta) - b
        r = (r + np.pi) % (2 * np.pi) - np.pi

        out = {
            "theta_e": theta,
            "residual_face": r,
            "residual_rmse": float(np.sqrt(np.mean(r**2))),
            "residual_maxabs": float(np.max(np.abs(r))),
            "ridge": float(ridge),
            "nE": int(nE),
            "nF": int(nF),
        }
        self._cache[cache_key] = out
        self._cache["theta_e"] = theta
        return out

    def solve_edge_thetas_unwrapped(self, ridge: float = 1e-6, max_iter: int = 8, k_cap: int = 5) -> Dict[str, Any]:
        """
        Verbesserte Face-Flux-Constraint-Lösung mit ganzzahligem 2π-Lift:

            sum(face edges) sgn*theta_e = phi_face + 2π k_face,   k_face ∈ ℤ.

        Algorithmus (iterativ):
          1) solve normal equations für theta_e (wie in solve_edge_thetas_from_face_flux)
          2) residual r_face = A theta - phi_face  (in R)
          3) update k_face := clip(round(r_face/(2π)), [-k_cap, k_cap])
          4) solve erneut für b := phi_face + 2π k_face

        Rückgabe:
          - theta_e
          - k_face (Lift)
          - residual statistics (wrapped + unwrapped)
        """
        cache_key = f"theta_e_unwrapped_r{ridge:.0e}_it{max_iter}_cap{k_cap}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        # Build incidence A once
        E = self.edges()
        nE = E.shape[0]
        nF = self.nF
        eidx = self.edge_index()
        F = self.faces.astype(int)

        rows, cols, vals = [], [], []
        for fi, (a, bV, c) in enumerate(F):
            for (u, v) in [(a, bV), (bV, c), (c, a)]:
                uu, vv = (u, v) if u < v else (v, u)
                j = eidx[(uu, vv)]
                sgn = +1.0 if (u < v) else -1.0
                rows.append(fi); cols.append(j); vals.append(sgn)

        from scipy.sparse import coo_matrix
        from scipy.sparse.linalg import spsolve
        from scipy.sparse import diags

        A = coo_matrix((np.array(vals), (np.array(rows), np.array(cols))), shape=(nF, nE)).tocsr()
        AtA = (A.T @ A).tocsr()
        AtA_reg = AtA + ridge * diags(np.ones(nE), 0, format="csr")

        phi = self.phase_faces().astype(float)  # in [0, 2π) typically small
        # Keep phi unwrapped as given; we'll add 2π*k.
        k_face = np.zeros(nF, dtype=int)
        theta = np.zeros(nE, dtype=float)

        for it in range(int(max_iter)):
            b = phi + (2.0 * np.pi) * k_face
            rhs = A.T @ b
            theta = spsolve(AtA_reg, rhs)
            theta = np.array(theta, dtype=float)

            r = (A @ theta) - phi  # unwrapped residual against original phi
            k_new = np.rint(r / (2.0 * np.pi)).astype(int)
            if k_cap is not None:
                k_new = np.clip(k_new, -int(k_cap), int(k_cap))

            if np.array_equal(k_new, k_face):
                # converged
                break
            k_face = k_new

        # Final residuals
        b_final = phi + (2.0 * np.pi) * k_face
        r_unwrapped = (A @ theta) - b_final
        r_wrapped = (r_unwrapped + np.pi) % (2.0 * np.pi) - np.pi

        out = {
            "theta_e": theta,
            "k_face": k_face,
            "iterations": int(it + 1),
            "residual_rmse_unwrapped": float(np.sqrt(np.mean(r_unwrapped ** 2))),
            "residual_maxabs_unwrapped": float(np.max(np.abs(r_unwrapped))),
            "residual_rmse_wrapped": float(np.sqrt(np.mean(r_wrapped ** 2))),
            "residual_maxabs_wrapped": float(np.max(np.abs(r_wrapped))),
            "ridge": float(ridge),
            "nE": int(nE),
            "nF": int(nF),
        }
        # Store a canonical theta for later use
        self._cache[cache_key] = out
        self._cache["theta_e"] = theta
        self._cache["k_face"] = k_face
        return out

    def chern_from_connection(self) -> Dict[str, Any]:
        """
        'Harter' Chern-Output direkt aus der diskreten Connection (Face-Flux):
          C = (1/2π) * sum_faces phi_face.

        Da die Face-Flux-Konstruktion hier aus ganzzahligen Flussquanten besteht,
        ist C (bis auf numerische Rundung) ein Integer.

        Hinweis: Das ist der Chern der *vorgegebenen* U(1)-Connection (Monopol-Fluss),
        nicht der aus einem Projektor/gefüllten Band extrahierte Chern.
        """
        phi = np.array(self.phase_faces(), dtype=float)
        total = float(np.sum(phi))
        C = total / (2.0 * np.pi)
        C_round = int(np.rint(C))
        return {
            "C_raw": float(C),
            "C_round": int(C_round),
            "total_flux_rad": float(total),
            "total_flux_quanta": float(C),
            "N_phi_global": int(self.N_phi_global),
            "delta_N_phi_walter": int(self.delta_N_phi_walter),
            "expected_total_quanta": int(self.N_phi_global) + int(self.delta_N_phi_walter),
        }

    def hall_quantities_from_nu(self, nu: int) -> Dict[str, Any]:
        """
        Quantum-Hall Größen in SI aus einem (hypothetisch) integer Füllfaktor ν.
        """
        nu = int(nu)
        if nu == 0:
            raise ValueError("nu darf nicht 0 sein.")
        return {
            "nu": nu,
            "sigma_xy_S": float(nu * self.si.conductance_quantum),      # Siemens
            "R_H_ohm": float(self.si.von_klitzing / nu),                # Ohm
            "R_K_ohm": float(self.si.von_klitzing),                     # Ohm
            "Phi0_Wb": float(self.si.flux_quantum_h_over_e),            # Weber
        }


    def edge_phases(self, ridge: float = 1e-6) -> Dict[Tuple[int, int], complex]:
        cache_key = f"edge_phases_r{ridge:.0e}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        sol = self.solve_edge_thetas_from_face_flux(ridge=ridge)
        theta = sol["theta_e"]
        E = self.edges()

        phases: Dict[Tuple[int, int], complex] = {}
        for k, (u, v) in enumerate(E):
            phases[(int(u), int(v))] = np.exp(1j * theta[k])

        if self.theta_morley and abs(self.theta_morley) > 1e-15:
            maskW = (self.face_channel.astype(int) == 1)
            walter_edges: Set[Tuple[int, int]] = set()
            for (a, b, c), isW in zip(self.faces.astype(int), maskW):
                if not isW:
                    continue
                for (u, v) in [(a, b), (b, c), (c, a)]:
                    uu, vv = (u, v) if u < v else (v, u)
                    walter_edges.add((uu, vv))
            tex = np.exp(1j * float(self.theta_morley))
            for ekey in walter_edges:
                phases[ekey] *= tex

        self._cache[cache_key] = phases
        return phases

    def L_magnetic(self, ridge: float = 1e-6):
        cache_key = f"L_A_r{ridge:.0e}"
        if cache_key in self._cache:
            return self._cache[cache_key]
        phases = self.edge_phases(ridge=ridge)
        L, m = self._assemble_sparse_L_M(phases=phases)
        self._cache[cache_key] = (L, m)
        return L, m

    def spectrum_magnetic(self, ridge: float = 1e-6):
        """
        Spektrum des magnetischen Operators (generalisiertes Problem über M):
          L x = λ M x
        Implementiert über symmetrisiertes A = M^{-1/2} L M^{-1/2}.
        """
        cache_key = f"spectrum_A_r{ridge:.0e}_k{int(self.eig_k)}_vecs{int(self.store_eigenvectors)}"
        if cache_key in self._cache:
            return self._cache[cache_key]
        L, m = self.L_magnetic(ridge=ridge)
        vals, vecs = self._smallest_k_eigs(L, m, k=int(self.eig_k), want_vecs=bool(self.store_eigenvectors))
        out = {"eigs": vals, "vecs": vecs}
        self._cache[cache_key] = out
        return out


    @staticmethod
    def _smallest_k_eigs(L: csr_matrix, mdiag: np.ndarray, k: int, want_vecs: bool):
        """Kleinste k Eigenwerte des symmetrisierten Problems A = M^{-1/2} L M^{-1/2}.

        Robustheit:
          - k wird automatisch auf <= N-2 gekappt (ARPACK-Bedingung).
          - Für sehr kleine N oder große k fällt die Routine auf dense eigh zurück.
        """
        mdiag = np.maximum(mdiag, 1e-18)
        # Build A (symmetric/Hermitian)
        Minv_sqrt = diags(1.0 / np.sqrt(mdiag), 0, format="csr")
        A = Minv_sqrt @ (L @ Minv_sqrt)
        N = int(A.shape[0])

        # Cap k for ARPACK; if N is tiny use dense
        k = int(k)
        if N <= 200 or k >= (N - 1):
            k_eff = max(1, min(k, N))
            Ad = A.toarray()
            Ad = 0.5 * (Ad + Ad.conj().T)
            w, vecs = np.linalg.eigh(Ad)
            w = np.real(w)
            w[w < 0] = 0.0
            idx = np.argsort(w)
            w = w[idx][:k_eff]
            if want_vecs:
                vecs = vecs[:, idx][:, :k_eff]
                return w, vecs
            return w, None

        k_eff = max(1, min(k, N - 2))
        if want_vecs:
            vals, vecs = eigsh(A, k=k_eff, which="SM", tol=1e-8)
            vals = np.real(vals)
            vals = np.maximum(vals, 0.0)
            idx = np.argsort(vals)
            return vals[idx], vecs[:, idx]
        vals = eigsh(A, k=k_eff, which="SM", return_eigenvectors=False, tol=1e-8)
        vals = np.real(vals)
        vals = np.maximum(vals, 0.0)
        vals.sort()
        return vals, None

    def estimate_landau_degeneracy(self, eigs: np.ndarray, gap_factor: float = 8.0, max_levels: int = 4) -> Dict[str, Any]:
        eigs = np.array(eigs, dtype=float)
        if eigs.size < 5:
            return {"clusters": [(0, int(eigs.size) - 1)], "sizes": [int(eigs.size)], "gaps": []}
        gaps = np.diff(eigs)
        base = float(np.median(gaps[:min(20, len(gaps))]))
        if base <= 0:
            base = float(np.mean(gaps) + 1e-12)
        thr = gap_factor * base

        clusters = []
        start = 0
        for i, g in enumerate(gaps):
            if g > thr:
                clusters.append((start, i))
                start = i + 1
                if len(clusters) >= max_levels:
                    break
        clusters.append((start, min(start + 200, len(eigs) - 1)))
        sizes = [j - i + 1 for (i, j) in clusters]
        return {
            "gap_base": base,
            "gap_threshold": thr,
            "clusters": clusters,
            "sizes": sizes,
        }

    def estimate_qhe_nu_and_chern(self, ridge: float = 1e-6) -> Dict[str, Any]:
        sp = self.spectrum_magnetic(ridge=ridge)
        eigs = sp["eigs"]
        cl = self.estimate_landau_degeneracy(eigs)
        expected_g0 = int(self.N_phi_global) + 1 if self.N_phi_global >= 0 else None
        g0 = int(cl["sizes"][0]) if cl["sizes"] else None

        nu_guess = None
        if expected_g0 is not None and g0 is not None and abs(g0 - expected_g0) <= max(1, expected_g0 // 10):
            nu_guess = 1

        return {
            "N_phi_global": int(self.N_phi_global),
            "delta_N_phi_walter": int(self.delta_N_phi_walter),
            "expected_lowest_LL_degeneracy": expected_g0,
            "observed_lowest_cluster_size": g0,
            "cluster_sizes": cl["sizes"],
            "cluster_ranges": cl["clusters"],
            "nu_guess_if_lowest_filled": nu_guess,
            "sigma_xy_if_nu_guess_SI": (nu_guess * self.si.conductance_quantum) if nu_guess is not None else None,
            "R_H_if_nu_guess_ohm": (self.si.von_klitzing / nu_guess) if (nu_guess not in (None, 0)) else None,
        }

    def edge_localization_scores(self, eigvecs: np.ndarray) -> np.ndarray:
        boundary = self.walter_boundary_vertices()
        if boundary.size == 0:
            return np.zeros(eigvecs.shape[1], dtype=float)
        V = eigvecs
        norms = np.sqrt(np.sum(np.abs(V) ** 2, axis=0))
        norms = np.where(norms > 0, norms, 1.0)
        Vn = V / norms
        probs = np.sum(np.abs(Vn[boundary, :]) ** 2, axis=0)
        return np.array(probs, dtype=float)

    def detect_edge_modes(self, ridge: float = 1e-6, top_k: int = 10, threshold: float = 0.25) -> Dict[str, Any]:
        sp = self.spectrum_magnetic(ridge=ridge)
        if sp["vecs"] is None:
            return {"error": "Keine Eigenvektoren gespeichert. Setze store_eigenvectors=True."}
        eigs = sp["eigs"]
        vecs = sp["vecs"]
        scores = self.edge_localization_scores(vecs)
        idx = np.argsort(scores)[::-1]
        top = idx[:min(top_k, idx.size)]
        edge_candidates = [int(i) for i in idx if float(scores[i]) >= threshold]
        return {
            "threshold": float(threshold),
            "top_k": int(top_k),
            "top_indices": top.tolist(),
            "top_scores": [float(scores[i]) for i in top],
            "top_eigs": [float(eigs[i]) for i in top],
            "edge_candidate_indices": edge_candidates[:50],
            "edge_candidate_count": int(len(edge_candidates)),
        }

    def energy_scale(self) -> float:
        if "energy_scale" in self._cache:
            return float(self._cache["energy_scale"])
        hbar = self.si.hbar
        scale = (hbar ** 2) / (2.0 * self.si.m_e * (self.radius_m ** 2))
        self._cache["energy_scale"] = float(scale)
        return float(scale)

    # --------------------
    # Geometry: area/volume + Kepler invariants
    # --------------------
    def polyhedron_area(self) -> float:
        """Gesamtfläche der triangulierten Oberfläche (in m^2, wenn vertices in m skaliert sind)."""
        V = self.vertices
        F = self.faces.astype(int)
        area = 0.0
        for (i,j,k) in F:
            vi, vj, vk = V[i], V[j], V[k]
            area += 0.5 * float(np.linalg.norm(np.cross(vj-vi, vk-vi)))
        return float(area)

    def polyhedron_volume(self) -> float:
        """
        Volumen des (orientierten) Polyeders über Divergenz-/Tetraederformel.
        Annahme: Oberfläche ist geschlossen und die Vertices sind um den Ursprung zentriert.
        Für Dreiecksface (i,j,k): V += dot(vi, cross(vj, vk))/6.
        """
        Vv = self.vertices
        F = self.faces.astype(int)
        vol = 0.0
        for (i,j,k) in F:
            vi, vj, vk = Vv[i], Vv[j], Vv[k]
            vol += float(np.dot(vi, np.cross(vj, vk))) / 6.0
        return float(abs(vol))

    def kepler_invariants(self) -> Dict[str, Any]:
        """
        Kepler-Schalter: Flächen durch sqrt(12), Volumen durch sqrt(18) normieren.
        Liefert außerdem Vergleich zu Kugel-Referenzen (4πR^2, 4/3πR^3).
        """
        A = self.polyhedron_area()
        V = self.polyhedron_volume()
        A_kepler = A / float(self.kepler_area_root)
        V_kepler = V / float(self.kepler_volume_root)

        R = float(self.radius_m)
        A_sphere = 4.0 * np.pi * (R**2)
        V_sphere = (4.0/3.0) * np.pi * (R**3)

        return {
            "use_kepler": bool(self.use_kepler),
            "kepler_area_root": float(self.kepler_area_root),
            "kepler_volume_root": float(self.kepler_volume_root),
            "A_poly": float(A),
            "V_poly": float(V),
            "A_kepler": float(A_kepler),
            "V_kepler": float(V_kepler),
            "A_sphere": float(A_sphere),
            "V_sphere": float(V_sphere),
            "A_poly_over_A_sphere": float(A / A_sphere) if A_sphere > 0 else None,
            "V_poly_over_V_sphere": float(V / V_sphere) if V_sphere > 0 else None,
            "A_kepler_over_A_sphere": float(A_kepler / A_sphere) if A_sphere > 0 else None,
            "V_kepler_over_V_sphere": float(V_kepler / V_sphere) if V_sphere > 0 else None,
        }

    # --------------------
    # QHE: Projektor-basierter Chern Marker (Bianco–Resta Stil)
    # --------------------
    def chern_projector_marker(self, ridge: float = 1e-6, n_occ: Optional[int] = None,
                               coords: str = "xy", use_dense: bool = True) -> Dict[str, Any]:
        """
        Real-space Chern Marker (heuristisch, aber 'projektor-basiert'):

          C ≈ -2π i Tr( P [X,P][Y,P] )

        wobei X,Y diagonale Matrizen der Koordinaten sind (Embedding).
        Für Kugel ist das keine perfekte torus-topologie, aber als Diagnose für 'topologische'
        Struktur / Robustheit nützlich.

        Parameter:
          n_occ: Anzahl besetzter Zustände (Eigenvektoren) aus dem magnetischen Spektrum.
                 Default: Größe des niedrigsten Landau-Clusters (falls detektierbar).
          coords: "xy", "yz", "zx" (Embedding-Achsen) oder "pca" (beste Ebene).
          use_dense: Für kleine nV sinnvoll (nV<=500). Für große Meshes später optimieren.
        """
        sp = self.spectrum_magnetic(ridge=ridge)
        eigs = sp["eigs"]
        vecs = sp["vecs"]
        if vecs is None:
            return {"error": "Eigenvektoren fehlen. Setze store_eigenvectors=True und berechne spectrum_magnetic."}

        # decide n_occ
        if n_occ is None:
            cl = self.estimate_landau_degeneracy(eigs)
            n_occ = int(cl["sizes"][0]) if cl.get("sizes") else min(5, vecs.shape[1])

        n_occ = max(1, min(int(n_occ), vecs.shape[1]))
        Psi = np.array(vecs[:, :n_occ], dtype=complex)

        # Coordinates
        Xv = self.vertices.astype(float) * float(self.radius_m)
        if coords == "xy":
            x = Xv[:,0]; y = Xv[:,1]
        elif coords == "yz":
            x = Xv[:,1]; y = Xv[:,2]
        elif coords == "zx":
            x = Xv[:,2]; y = Xv[:,0]
        elif coords == "pca":
            Xc = Xv - Xv.mean(axis=0, keepdims=True)
            U, S, VT = np.linalg.svd(Xc, full_matrices=False)
            p1, p2 = VT[0], VT[1]
            x = Xc @ p1
            y = Xc @ p2
        else:
            raise ValueError("coords must be one of: xy,yz,zx,pca")

        x = np.array(x, dtype=float)
        y = np.array(y, dtype=float)

        if use_dense:
            P = Psi @ Psi.conj().T
            dx = (x[:,None] - x[None,:])
            dy = (y[:,None] - y[None,:])
            XP = dx * P
            YP = dy * P
            term = P @ XP @ YP
            C = -2.0*np.pi * np.imag(np.trace(term))
            return {
                "coords": coords,
                "n_occ": int(n_occ),
                "C_projector": float(C),
                "note": "Für Kugel-Embedding als Diagnose; nicht identisch mit Torus-Chern.",
            }
        else:
            P = Psi @ Psi.conj().T
            dx = (x[:,None] - x[None,:])
            dy = (y[:,None] - y[None,:])
            XP = dx * P
            YP = dy * P
            term = P @ XP @ YP
            C = -2.0*np.pi * np.imag(np.trace(term))
            return {"coords": coords, "n_occ": int(n_occ), "C_projector": float(C)}

    # --------------------
    # Spectral flow test: vary Walter-excess or global flux
    # --------------------
    def spectral_flow(self, ridge: float = 1e-6, vary: str = "dNw",
                      values: Optional[List[int]] = None, k: int = 20) -> Dict[str, Any]:
        """
        Spektralfluss-Test:
          vary="dNw": variiert delta_N_phi_walter
          vary="Nphi": variiert N_phi_global

        Gibt für jeden Schritt die ersten k Eigenwerte zurück plus einfache Gap-Diagnostik.
        """
        if values is None:
            values = list(range(0, 6))
        values = [int(v) for v in values]
        base_Nphi = int(self.N_phi_global)
        base_dNw = int(self.delta_N_phi_walter)

        rows = []
        for v in values:
            if vary == "dNw":
                self.flux = (base_Nphi, v)
            elif vary == "Nphi":
                self.flux = (v, base_dNw)
            else:
                raise ValueError("vary must be 'dNw' or 'Nphi'")

            self.store_eigenvectors = False
            self.eig_k = int(k)
            sp = self.spectrum_magnetic(ridge=ridge)
            eigs = np.array(sp["eigs"], dtype=float)
            gap1 = float(eigs[1]-eigs[0]) if eigs.size >= 2 else None
            gap10 = float(eigs[min(10, eigs.size-1)] - eigs[min(9, eigs.size-1)]) if eigs.size >= 11 else None

            rows.append({
                "value": int(v),
                "first_eigs": [float(x) for x in eigs[:min(k, eigs.size)]],
                "gap_01": gap1,
                "gap_9_10": gap10
            })

        self.flux = (base_Nphi, base_dNw)
        return {"vary": vary, "values": values, "ridge": float(ridge), "k": int(k), "rows": rows}

    # --------------------
    # Angular momentum diagnostics (dimensions / l-estimates)
    # --------------------
    def angular_momentum_diagnostics(self, ridge: float = 1e-6, k: int = 30) -> Dict[str, Any]:
        """
        Für eine Kugel gilt im Kontinuum: Laplace-Beltrami Eigenwerte
           λ_l = l(l+1)/R^2
        und Drehimpuls:
           L^2 = ħ^2 l(l+1)
        Daraus folgt:
           L2_est = ħ^2 * λ * R^2
           l_est = (-1 + sqrt(1+4 λ R^2))/2

        Wir nutzen das als Dimensions-/Strukturtest für dein Modell.
        """
        self.eig_k = int(k)
        self.store_eigenvectors = False
        sp = self.spectrum_magnetic(ridge=ridge)
        eigs = np.array(sp["eigs"], dtype=float)

        R = float(self.radius_m)
        hbar = float(self.si.hbar)

        l_est = []
        L2_est = []
        for lam in eigs[:k]:
            x = float(lam) * (R**2)
            l = 0.5 * (-1.0 + math.sqrt(max(0.0, 1.0 + 4.0*x)))
            l_est.append(float(l))
            L2_est.append(float((hbar**2) * x))

        l_round = [int(round(v)) for v in l_est]
        from collections import Counter
        cnt = Counter(l_round)

        return {
            "ridge": float(ridge),
            "k": int(k),
            "R_m": float(R),
            "hbar_Js": float(hbar),
            "first_eigs": [float(x) for x in eigs[:min(12, eigs.size)]],
            "l_est_first": [float(x) for x in l_est[:min(12, len(l_est))]],
            "l_round_counts": dict(cnt),
            "L2_est_first_SI": [float(x) for x in L2_est[:min(12, len(L2_est))]],
            "note": "Wenn Cluster wie (l=0:1, l=1:3, l=2:5, ...) auftauchen, ist die Kugelstruktur gut getroffen."
        }


    def to_dict_meta(self) -> Dict[str, Any]:
        return {
            "model_version": MODEL_VERSION,
            "geometry": {"nV": self.nV, "nF": self.nF, "radius_m": float(self.radius_m), "use_kepler": bool(self.use_kepler), "kepler_area_root": float(self.kepler_area_root), "kepler_volume_root": float(self.kepler_volume_root)},
            "params": {
                "tau0_s": None if self.tau0_s is None else float(self.tau0_s),
                "N_phi_global": int(self.N_phi_global),
                "delta_N_phi_walter": int(self.delta_N_phi_walter),
                "theta_morley": float(self.theta_morley),
                "w_walter": float(self.w_walter),
                "w_rest": float(self.w_rest),
                "eig_k": int(self.eig_k),
                "store_eigenvectors": bool(self.store_eigenvectors),
            },
            "si": {
                "h": float(self.si.h),
                "c": float(self.si.c),
                "e": float(self.si.e),
                "eps0": float(self.si.eps0),
                "m_e": float(self.si.m_e),
            },
            "cache_keys": sorted(list(self._cache.keys())),
        }

    def save(self, folder: str | Path, include_cache: bool = True) -> Path:
        folder = _as_path(folder)
        _ensure_dir(folder)
        _json_dump(self.to_dict_meta(), folder / "meta.json")

        arrays = {
            "vertices": self.vertices.astype(np.float64),
            "faces": self.faces.astype(np.int64),
            "face_channel": self.face_channel.astype(np.int8),
        }
        if include_cache:
            if "phase_faces" in self._cache:
                arrays["phase_faces"] = np.array(self._cache["phase_faces"], dtype=np.float64)
            if "theta_e" in self._cache:
                arrays["theta_e"] = np.array(self._cache["theta_e"], dtype=np.float64)
            if "edges" in self._cache:
                arrays["edges"] = np.array(self._cache["edges"], dtype=np.int64)
            # store one magnetic spectrum if present
            for k in list(self._cache.keys()):
                if k.startswith("spectrum_A_r"):
                    arrays["eigs_A"] = np.array(self._cache[k]["eigs"], dtype=np.float64)
                    if self._cache[k]["vecs"] is not None:
                        arrays["vecs_A"] = np.array(self._cache[k]["vecs"])
                    arrays["spectrum_A_key"] = np.array([k], dtype=object)
                    break

        np.savez_compressed(folder / "data.npz", **arrays)
        return folder

    @classmethod
    def load(cls, folder: str | Path) -> "BambergQHEModel":
        folder = _as_path(folder)
        meta = _json_load(folder / "meta.json")
        data = np.load(folder / "data.npz", allow_pickle=True)

        si = SIConstants(
            h=float(meta["si"]["h"]),
            c=float(meta["si"]["c"]),
            e=float(meta["si"]["e"]),
            eps0=float(meta["si"]["eps0"]),
            m_e=float(meta["si"]["m_e"]),
        )

        obj = cls(
            vertices=data["vertices"],
            faces=data["faces"],
            face_channel=data["face_channel"],
            radius_m=float(meta["geometry"]["radius_m"]),
            use_kepler=bool(meta["geometry"].get("use_kepler", False)),
            kepler_area_root=float(meta["geometry"].get("kepler_area_root", 12.0**0.5)),
            kepler_volume_root=float(meta["geometry"].get("kepler_volume_root", 18.0**0.5)),
            si=si,
            tau0_s=meta["params"]["tau0_s"],
            N_phi_global=int(meta["params"]["N_phi_global"]),
            delta_N_phi_walter=int(meta["params"]["delta_N_phi_walter"]),
            theta_morley=float(meta["params"]["theta_morley"]),
            w_walter=float(meta["params"]["w_walter"]),
            w_rest=float(meta["params"]["w_rest"]),
            eig_k=int(meta["params"]["eig_k"]),
            store_eigenvectors=bool(meta["params"]["store_eigenvectors"]),
        )

        if "phase_faces" in data.files:
            obj._cache["phase_faces"] = data["phase_faces"]
        if "theta_e" in data.files:
            obj._cache["theta_e"] = data["theta_e"]
        if "edges" in data.files:
            obj._cache["edges"] = data["edges"]
        if "eigs_A" in data.files:
            key = str(data["spectrum_A_key"][0]) if "spectrum_A_key" in data.files else "spectrum_A_r1e-06"
            obj._cache[key] = {"eigs": data["eigs_A"], "vecs": data["vecs_A"] if "vecs_A" in data.files else None}
        return obj


def _load_mesh_json(mesh_json: Path):
    d = json.loads(mesh_json.read_text(encoding="utf-8"))
    V = np.array(d["vertices"], dtype=float)
    F = np.array(d["faces"], dtype=int)
    ch = np.array(d.get("face_channel", np.zeros(len(F), dtype=int)), dtype=int)
    r = float(d.get("radius", 1.0))
    return V, F, ch, r


def main():
    import argparse
    ap = argparse.ArgumentParser(prog="bm_qhe_model_v2", description="BambergQHEModel v0.2")
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_init = sub.add_parser("init")
    ap_init.add_argument("--mesh", required=True)
    ap_init.add_argument("--out", required=True)
    ap_init.add_argument("--R", type=float, default=1.0)

    ap_run = sub.add_parser("run")
    ap_run.add_argument("--state", required=True)
    ap_run.add_argument("--Nphi", type=int, default=0)
    ap_run.add_argument("--dNw", type=int, default=0)
    ap_run.add_argument("--theta", type=float, default=0.0)
    ap_run.add_argument("--ridge", type=float, default=1e-6)
    ap_run.add_argument("--k", type=int, default=80)
    ap_run.add_argument("--save", action="store_true")

    ap_diag = sub.add_parser("diag")
    ap_diag.add_argument("--state", required=True)
    ap_diag.add_argument("--ridge", type=float, default=1e-6)
    ap_diag.add_argument("--edge-threshold", type=float, default=0.25)
    ap_diag.add_argument("--topk", type=int, default=10)

    ap_unwrap = sub.add_parser("unwrap")
    ap_unwrap.add_argument("--state", required=True)
    ap_unwrap.add_argument("--ridge", type=float, default=1e-6)
    ap_unwrap.add_argument("--max-iter", type=int, default=8)
    ap_unwrap.add_argument("--k-cap", type=int, default=5)

    ap_chern = sub.add_parser("chern")
    ap_chern.add_argument("--state", required=True)
    ap_pchern = sub.add_parser("projector-chern")
    ap_pchern.add_argument("--state", required=True)
    ap_pchern.add_argument("--ridge", type=float, default=1e-6)
    ap_pchern.add_argument("--n-occ", type=int, default=0)
    ap_pchern.add_argument("--coords", type=str, default="xy")

    ap_flow = sub.add_parser("flow")
    ap_flow.add_argument("--state", required=True)
    ap_flow.add_argument("--ridge", type=float, default=1e-6)
    ap_flow.add_argument("--vary", type=str, default="dNw")
    ap_flow.add_argument("--values", type=str, default="0,1,2,3,4,5")
    ap_flow.add_argument("--k", type=int, default=20)

    ap_ldiag = sub.add_parser("ldiag")
    ap_ldiag.add_argument("--state", required=True)
    ap_ldiag.add_argument("--ridge", type=float, default=1e-6)
    ap_ldiag.add_argument("--k", type=int, default=30)

    ap_kepler = sub.add_parser("kepler")
    ap_kepler.add_argument("--state", required=True)


    args = ap.parse_args()

    if args.cmd == "init":
        V, F, ch, r = _load_mesh_json(_as_path(args.mesh))
        m = BambergQHEModel(vertices=V, faces=F, face_channel=ch, radius_m=float(args.R))
        m.save(args.out, include_cache=True)
        print("[init] saved:", args.out)

    elif args.cmd == "run":
        m = BambergQHEModel.load(args.state)
        m.store_eigenvectors = True
        m.eig_k = int(args.k)
        m.theta_morley = float(args.theta)
        m.flux = (int(args.Nphi), int(args.dNw))

        sol = m.solve_edge_thetas_from_face_flux(ridge=float(args.ridge))
        sp = m.spectrum_magnetic(ridge=float(args.ridge))

        print("[run] residual rmse:", sol["residual_rmse"], "maxabs:", sol["residual_maxabs"])
        print("[run] first 12 eigs:", sp["eigs"][:12].tolist())

        if args.save:
            m.save(args.state, include_cache=True)
            print("[run] saved updated state")

    elif args.cmd == "diag":
        m = BambergQHEModel.load(args.state)
        sol = m.solve_edge_thetas_from_face_flux(ridge=float(args.ridge))
        qhe = m.estimate_qhe_nu_and_chern(ridge=float(args.ridge))
        edge = m.detect_edge_modes(ridge=float(args.ridge), top_k=int(args.topk), threshold=float(args.edge_threshold))
        print("[diag] residual rmse:", sol["residual_rmse"], "maxabs:", sol["residual_maxabs"])
        print("[diag] qhe:", json.dumps(qhe, indent=2))
        print("[diag] edge:", json.dumps(edge, indent=2))


    elif args.cmd == "unwrap":
        m = BambergQHEModel.load(args.state)
        out = m.solve_edge_thetas_unwrapped(ridge=float(args.ridge), max_iter=int(args.max_iter), k_cap=int(args.k_cap))
        print("[unwrap]", json.dumps(out, indent=2))

    elif args.cmd == "chern":
        m = BambergQHEModel.load(args.state)
        ch = m.chern_from_connection()
        print("[chern]", json.dumps(ch, indent=2))

    elif args.cmd == "projector-chern":
        m = BambergQHEModel.load(args.state)
        n_occ = int(args.n_occ) if int(args.n_occ) > 0 else None
        m.store_eigenvectors = True
        out = m.chern_projector_marker(ridge=float(args.ridge), n_occ=n_occ, coords=str(args.coords), use_dense=True)
        print("[projector-chern]", json.dumps(out, indent=2))

    elif args.cmd == "flow":
        m = BambergQHEModel.load(args.state)
        vals = [int(x.strip()) for x in str(args.values).split(",") if x.strip() != ""]
        out = m.spectral_flow(ridge=float(args.ridge), vary=str(args.vary), values=vals, k=int(args.k))
        print("[flow]", json.dumps(out, indent=2))

    elif args.cmd == "ldiag":
        m = BambergQHEModel.load(args.state)
        out = m.angular_momentum_diagnostics(ridge=float(args.ridge), k=int(args.k))
        print("[ldiag]", json.dumps(out, indent=2))

    elif args.cmd == "kepler":
        m = BambergQHEModel.load(args.state)
        out = m.kepler_invariants()
        print("[kepler]", json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
