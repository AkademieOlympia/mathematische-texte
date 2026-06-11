# bm_qhe_model.py
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, List

import json
import numpy as np

try:
    from scipy.sparse import csr_matrix, coo_matrix, diags
    from scipy.sparse.linalg import eigsh
except ImportError as e:
    raise ImportError("Dieses Modul benötigt scipy (sparse + eigsh).") from e


# ----------------------------
# Utilities
# ----------------------------

MODEL_VERSION = "0.1.0"

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


# ----------------------------
# Physical constants helper
# ----------------------------

@dataclass
class SIConstants:
    """
    SI: h, e, c sind (seit 2019) exakt definiert.
    m_e ist empfohlen (CODATA), kann aber als Parameter gesetzt werden.
    """
    h: float = 6.62607015e-34          # J s (exakt)
    c: float = 299792458.0             # m/s (exakt)
    e: float = 1.602176634e-19         # C (exakt)
    eps0: float = 8.8541878128e-12     # F/m (aus mu0/c^2; hier als Zahl)
    m_e: float = 9.1093837139e-31      # kg (empfohlen)
    
    @property
    def hbar(self) -> float:
        return self.h / (2.0 * np.pi)

    @property
    def compton_frequency(self) -> float:
        # f_C = m c^2 / h
        return self.m_e * self.c**2 / self.h

    @property
    def compton_omega(self) -> float:
        return 2.0 * np.pi * self.compton_frequency

    @property
    def compton_period(self) -> float:
        return 1.0 / self.compton_frequency

    @property
    def compton_tau_bar(self) -> float:
        # reduced Compton time taū = 1/omega = ħ/(m c^2)
        return 1.0 / self.compton_omega

    @property
    def flux_quantum_h_over_e(self) -> float:
        # Φ0 = h/e  (QHE-Flussquantum)
        return self.h / self.e

    @property
    def von_klitzing(self) -> float:
        # R_K = h / e^2
        return self.h / (self.e**2)

    @property
    def conductance_quantum(self) -> float:
        # e^2/h
        return (self.e**2) / self.h


# ----------------------------
# Main model class
# ----------------------------

@dataclass
class BambergQHEModel:
    """
    Bamberg-Kugel + zweikanalige Walter/Rest-Zerlegung + U(1)-Gauge (Morley/QHE)
    
    Modus (C):
      - globaler Grundfluss N_phi
      - Walter-Exzessfluss delta_N_phi_W (zusätzliche Flussquanten "im Walter-Kanal")
      
    Caches:
      - Mesh derived: edges, face areas, etc.
      - Operatoren: L, M, magnetischer L_A
      - Spektren: eigenvalues, ggf. eigenvectors (optional)
    """
    # ----- geometry -----
    vertices: np.ndarray                 # (nV,3)
    faces: np.ndarray                    # (nF,3) int
    face_channel: np.ndarray             # (nF,) 1=Walter, 0=Rest
    radius_m: float = 1.0                # Radius in Meter (oder BM-Skala -> später mappen)
    
    # ----- physical params -----
    si: SIConstants = field(default_factory=SIConstants)

    # internal time: Peano tick Δn=1 -> Δt
    tau0_s: Optional[float] = None       # default: taū_C (reduced Compton time)
    
    # QHE / Gauge:
    N_phi_global: int = 0                # globale Flussquanten (Haldane sphere)
    delta_N_phi_walter: int = 0          # zusätzlicher Fluss im Walter-Kanal (droplet)
    theta_morley: float = 2.0*np.pi/3.0  # optional: Morley-Phase als Grundphase
    
    # Zehntel-Kanalgewicht (operatorisch)
    w_walter: float = 0.1
    w_rest: float = 0.9
    
    # solver params
    eig_k: int = 40
    store_eigenvectors: bool = False
    
    # internal caches
    _cache: Dict[str, Any] = field(default_factory=dict, init=False, repr=False)

    # --------------------
    # Getters / setters (Properties)
    # --------------------
    @property
    def nV(self) -> int:
        return int(self.vertices.shape[0])

    @property
    def nF(self) -> int:
        return int(self.faces.shape[0])

    @property
    def tau(self) -> float:
        # default to reduced Compton time
        if self.tau0_s is None:
            return self.si.compton_tau_bar
        return float(self.tau0_s)

    @tau.setter
    def tau(self, value: float) -> None:
        self.tau0_s = float(value)
        # Zeit selbst invalidiert Operatoren nicht, aber Energie-Skalierungen/Derived outputs
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
        _invalidate(self._cache, ["L_mix", "spectrum_mix"])

    @property
    def flux(self) -> Tuple[int, int]:
        return (int(self.N_phi_global), int(self.delta_N_phi_walter))

    @flux.setter
    def flux(self, v: Tuple[int, int]) -> None:
        self.N_phi_global = int(v[0])
        self.delta_N_phi_walter = int(v[1])
        _invalidate(self._cache, ["phase_faces", "L_A", "spectrum_A", "chern"])

    # --------------------
    # Geometry derived
    # --------------------
    def edges(self) -> np.ndarray:
        """
        Returns sorted undirected edges array (nE,2).
        Cached.
        """
        if "edges" in self._cache:
            return self._cache["edges"]
        edge_set = set()
        for a,b,c in self.faces.astype(int):
            edge_set.add(tuple(sorted((a,b))))
            edge_set.add(tuple(sorted((b,c))))
            edge_set.add(tuple(sorted((c,a))))
        E = np.array(sorted(edge_set), dtype=int)
        self._cache["edges"] = E
        return E

    def face_areas_chord(self) -> np.ndarray:
        """
        Chord-area on embedding (not exact spherical area).
        Cached; good enough for weighting/partition.
        """
        if "face_area" in self._cache:
            return self._cache["face_area"]
        V = self.vertices
        A = np.zeros(self.nF, dtype=float)
        for i,(a,b,c) in enumerate(self.faces.astype(int)):
            va,vb,vc = V[a],V[b],V[c]
            A[i] = 0.5*np.linalg.norm(np.cross(vb-va, vc-va))
        self._cache["face_area"] = A
        return A

    # --------------------
    # DEC / Laplacian assembly
    # --------------------
    @staticmethod
    def _cotangent(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
        ba = b - a
        ca = c - a
        cross = np.linalg.norm(np.cross(ba, ca))
        dot = float(np.dot(ba, ca))
        if cross <= 0:
            return 0.0
        return dot / cross

    def _assemble_sparse_L_M(self, face_mask: Optional[np.ndarray]=None, phases: Optional[Dict[Tuple[int,int], complex]]=None
                             ) -> Tuple[csr_matrix, np.ndarray]:
        """
        Assemble Laplacian L and lumped mass diag mdiag.
        If phases given: applies U(1) phase on edges (i<j key).
        """
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
            uu,vv = (u,v) if u<v else (v,u)
            phase = 1.0+0j
            if phases is not None and (uu,vv) in phases:
                phase = phases[(uu,vv)]
            I.extend([u,v])
            J.extend([v,u])
            data.extend([-w*phase, -w*np.conjugate(phase)])
            deg[u] += w
            deg[v] += w

        V = self.vertices
        F = self.faces.astype(int)

        for idx,(i,j,k) in enumerate(F):
            if not face_mask[idx]:
                continue
            vi,vj,vk = V[i],V[j],V[k]
            area = 0.5*np.linalg.norm(np.cross(vj-vi, vk-vi))
            mdiag[i] += area/3.0
            mdiag[j] += area/3.0
            mdiag[k] += area/3.0

            cot_i = self._cotangent(vi, vj, vk)
            cot_j = self._cotangent(vj, vk, vi)
            cot_k = self._cotangent(vk, vi, vj)
            add_edge(j,k, 0.5*cot_i)
            add_edge(k,i, 0.5*cot_j)
            add_edge(i,j, 0.5*cot_k)

        # diagonal
        for u in range(nV):
            I.append(u); J.append(u); data.append(deg[u])

        L = coo_matrix((np.array(data), (np.array(I), np.array(J))), shape=(nV,nV)).tocsr()
        return L, mdiag

    def L_full(self) -> Tuple[csr_matrix, np.ndarray]:
        if "L_full" in self._cache:
            return self._cache["L_full"]
        L, m = self._assemble_sparse_L_M()
        self._cache["L_full"] = (L, m)
        return L, m

    def L_channels(self) -> Dict[str, Tuple[csr_matrix, np.ndarray]]:
        """
        Returns dict with L_W, L_R.
        """
        if "L_W" in self._cache and "L_R" in self._cache:
            return {"W": self._cache["L_W"], "R": self._cache["L_R"]}
        mask_W = (self.face_channel.astype(int) == 1)
        mask_R = ~mask_W
        LW, mW = self._assemble_sparse_L_M(face_mask=mask_W)
        LR, mR = self._assemble_sparse_L_M(face_mask=mask_R)
        self._cache["L_W"] = (LW, mW)
        self._cache["L_R"] = (LR, mR)
        return {"W": (LW,mW), "R": (LR,mR)}

    def L_mix_zehntel(self) -> Tuple[csr_matrix, np.ndarray]:
        """
        Operatorisch erzwungener Zehntel-Mix:
          L_mix = w_W * sW*L_W + w_R * sR*L_R
        mit sW,sR so, dass beide Kanäle auf Full-Area skaliert werden.
        """
        if "L_mix" in self._cache:
            return self._cache["L_mix"]
        Lfull, mfull = self.L_full()
        area_full = float(np.sum(mfull))

        chans = self.L_channels()
        LW, mW = chans["W"]
        LR, mR = chans["R"]

        area_W = float(np.sum(mW))
        area_R = float(np.sum(mR))
        sW = area_full / area_W if area_W > 0 else 1.0
        sR = area_full / area_R if area_R > 0 else 1.0

        Lmix = (self.w_walter*sW)*LW + (self.w_rest*sR)*LR
        mmix = (self.w_walter*sW)*mW + (self.w_rest*sR)*mR

        self._cache["L_mix"] = (Lmix, mmix, {"sW":sW, "sR":sR, "area_full":area_full, "area_W":area_W, "area_R":area_R})
        return Lmix, mmix

    # --------------------
    # QHE: phase/flux assignment (Mode C)
    # --------------------
    def phase_faces(self) -> np.ndarray:
        """
        Face flux phases phi_f so that:
          sum_f phi_f = 2π N_phi_global + 2π delta_N_phi_walter  (extra concentrated in Walter channel)
        But we implement mode (C) as:
          global uniform base flux + extra flux distributed uniformly over Walter faces.

        Returns phi_f (radians) per face.
        """
        if "phase_faces" in self._cache:
            return self._cache["phase_faces"]

        F = self.faces
        mask_W = (self.face_channel.astype(int) == 1)
        nW = int(np.sum(mask_W))
        nF = self.nF

        # base uniform flux per face:
        phi_base_total = 2.0*np.pi*int(self.N_phi_global)
        phi_base = phi_base_total / max(nF, 1)

        # extra flux only on Walter faces:
        phi_extra_total = 2.0*np.pi*int(self.delta_N_phi_walter)
        phi_extra = phi_extra_total / max(nW, 1) if nW > 0 else 0.0

        phi = np.full(nF, phi_base, dtype=float)
        phi[mask_W] += phi_extra

        self._cache["phase_faces"] = phi
        return phi

    def phases_on_edges_from_face_flux(self) -> Dict[Tuple[int,int], complex]:
        """
        Build a consistent (simple) edge-phase dictionary from per-face flux.
        For robust Chern/Berry you would do a gauge-fixing step; for now:
          - assign each face's flux equally to its 3 edges (orientation consistent up to gauge)
        This yields a "magnetic" Laplacian that encodes net flux.

        NOTE: This is a pragmatic starter; later we replace by a proper discrete connection solver
              that enforces UijUjkUki = exp(i*phi_f) exactly with consistent edge variables.
        """
        if "edge_phases" in self._cache:
            return self._cache["edge_phases"]

        phi_f = self.phase_faces()
        phases: Dict[Tuple[int,int], complex] = {}

        # naive equal split: each face contributes phi/3 to each of its edges.
        # accumulate contributions additively in the exponent.
        accum: Dict[Tuple[int,int], float] = {}

        for idx, (a,b,c) in enumerate(self.faces.astype(int)):
            d = phi_f[idx] / 3.0
            for u,v in [(a,b),(b,c),(c,a)]:
                uu,vv = (u,v) if u<v else (v,u)
                accum[(uu,vv)] = accum.get((uu,vv), 0.0) + d

        for ekey, theta in accum.items():
            phases[ekey] = np.exp(1j*theta)

        # optional: embed a Morley-like additional phase on Walter edges (as a baseline)
        # This is separate from flux; you may want it as a "spin/texture" field.
        if self.theta_morley is not None:
            # if you want: apply Morley phase only to edges that belong to Walter faces
            mask_W = (self.face_channel.astype(int) == 1)
            walter_edges = set()
            for (f,isW) in zip(self.faces.astype(int), mask_W):
                if not isW:
                    continue
                a,b,c = map(int,f)
                for u,v in [(a,b),(b,c),(c,a)]:
                    uu,vv = (u,v) if u<v else (v,u)
                    walter_edges.add((uu,vv))
            for ekey in walter_edges:
                phases[ekey] *= np.exp(1j*self.theta_morley)

        self._cache["edge_phases"] = phases
        return phases

    def L_magnetic(self) -> Tuple[csr_matrix, np.ndarray]:
        """
        Magnetic Laplacian (nabla - iA)^2 discretized by edge phases.
        Cached.
        """
        if "L_A" in self._cache:
            return self._cache["L_A"]
        phases = self.phases_on_edges_from_face_flux()
        L, m = self._assemble_sparse_L_M(phases=phases)
        self._cache["L_A"] = (L, m)
        return L, m

    # --------------------
    # Spectra
    # --------------------
    @staticmethod
    def _smallest_k_eigs(L: csr_matrix, mdiag: np.ndarray, k: int, want_vecs: bool=False):
        mdiag = np.maximum(mdiag, 1e-18)
        Minv_sqrt = diags(1.0/np.sqrt(mdiag), 0, format="csr")
        A = Minv_sqrt @ (L @ Minv_sqrt)

        # A is Hermitian if phases are used correctly; eigsh works for real/complex Hermitian
        if want_vecs:
            vals, vecs = eigsh(A, k=k, which="SM", tol=1e-8)
            vals = np.real(vals)
            vals = np.maximum(vals, 0.0)
            idx = np.argsort(vals)
            return vals[idx], vecs[:, idx]
        else:
            vals = eigsh(A, k=k, which="SM", return_eigenvectors=False, tol=1e-8)
            vals = np.real(vals)
            vals = np.maximum(vals, 0.0)
            vals.sort()
            return vals, None

    def spectrum_full(self) -> Dict[str, Any]:
        if "spectrum_full" in self._cache:
            return self._cache["spectrum_full"]
        L, m = self.L_full()
        vals, vecs = self._smallest_k_eigs(L, m, k=self.eig_k, want_vecs=self.store_eigenvectors)
        out = {"eigs": vals, "vecs": vecs}
        self._cache["spectrum_full"] = out
        return out

    def spectrum_magnetic(self) -> Dict[str, Any]:
        if "spectrum_A" in self._cache:
            return self._cache["spectrum_A"]
        L, m = self.L_magnetic()
        vals, vecs = self._smallest_k_eigs(L, m, k=self.eig_k, want_vecs=self.store_eigenvectors)
        out = {"eigs": vals, "vecs": vecs}
        self._cache["spectrum_A"] = out
        return out

    def energy_scale(self) -> float:
        """
        Convert dimensionless eigenvalues into energy scale:
          E ~ (ħ^2 / 2m) * (lambda / R^2)
        If radius_m is in meters, E is Joule.
        """
        if "energy_scale" in self._cache:
            return float(self._cache["energy_scale"])
        hbar = self.si.hbar
        scale = (hbar**2) / (2.0*self.si.m_e*(self.radius_m**2))
        self._cache["energy_scale"] = float(scale)
        return float(scale)

    # --------------------
    # (Placeholder) Chern number / QHE observables
    # --------------------
    def chern_number_placeholder(self) -> float:
        """
        Platzhalter: Für echtes QHE brauchst du Berry-Krümmung/Projektor.
        Hier nur ein stub, der zeigt, wo es hingehört.
        """
        # TODO: implement Fukui-Hatsugai-Suzuki (FHS) discrete Chern
        # and allow "filled band" selection (e.g. first q eigenvectors)
        return float("nan")


    # --------------------
    # Serialization
    # --------------------
    def to_dict_meta(self) -> Dict[str, Any]:
        return {
            "model_version": MODEL_VERSION,
            "geometry": {
                "nV": self.nV,
                "nF": self.nF,
                "radius_m": float(self.radius_m),
            },
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

    def save(self, folder: str | Path, include_cache: bool=True) -> Path:
        """
        Serialisiert Modell in einen Ordner:
          - meta.json
          - data.npz (vertices, faces, labels, ggf. caches)
        """
        folder = _as_path(folder)
        _ensure_dir(folder)

        meta_path = folder / "meta.json"
        npz_path = folder / "data.npz"

        meta = self.to_dict_meta()
        _json_dump(meta, meta_path)

        arrays = {
            "vertices": self.vertices.astype(np.float64),
            "faces": self.faces.astype(np.int64),
            "face_channel": self.face_channel.astype(np.int8),
        }

        # optionale Cache-Arrays (nur wenn sinnvoll)
        if include_cache:
            # Spektren speichern ist oft Gold wert:
            if "spectrum_full" in self._cache:
                arrays["eigs_full"] = np.array(self._cache["spectrum_full"]["eigs"])
                if self._cache["spectrum_full"]["vecs"] is not None:
                    arrays["vecs_full"] = np.array(self._cache["spectrum_full"]["vecs"])
            if "spectrum_A" in self._cache:
                arrays["eigs_A"] = np.array(self._cache["spectrum_A"]["eigs"])
                if self._cache["spectrum_A"]["vecs"] is not None:
                    arrays["vecs_A"] = np.array(self._cache["spectrum_A"]["vecs"])
            # face flux phi_f ist klein, aber spart Zeit
            if "phase_faces" in self._cache:
                arrays["phase_faces"] = np.array(self._cache["phase_faces"])

            # edges sind auch sinnvoll:
            if "edges" in self._cache:
                arrays["edges"] = np.array(self._cache["edges"]).astype(np.int64)

        np.savez_compressed(npz_path, **arrays)
        return folder

    @classmethod
    def load(cls, folder: str | Path) -> "BambergQHEModel":
        """
        Lädt Modell aus Ordner (meta.json + data.npz) und rekonstruiert Objekt inkl. wichtiger Caches.
        """
        folder = _as_path(folder)
        meta = _json_load(folder / "meta.json")
        data = np.load(folder / "data.npz", allow_pickle=False)

        vertices = data["vertices"]
        faces = data["faces"]
        face_channel = data["face_channel"]

        si = SIConstants(
            h=float(meta["si"]["h"]),
            c=float(meta["si"]["c"]),
            e=float(meta["si"]["e"]),
            eps0=float(meta["si"]["eps0"]),
            m_e=float(meta["si"]["m_e"]),
        )

        obj = cls(
            vertices=vertices,
            faces=faces,
            face_channel=face_channel,
            radius_m=float(meta["geometry"]["radius_m"]),
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

        # Restore cached arrays if present
        cache = {}
        if "edges" in data.files:
            cache["edges"] = data["edges"]
        if "phase_faces" in data.files:
            cache["phase_faces"] = data["phase_faces"]

        if "eigs_full" in data.files:
            cache["spectrum_full"] = {"eigs": data["eigs_full"], "vecs": data["vecs_full"] if "vecs_full" in data.files else None}
        if "eigs_A" in data.files:
            cache["spectrum_A"] = {"eigs": data["eigs_A"], "vecs": data["vecs_A"] if "vecs_A" in data.files else None}

        obj._cache.update(cache)
        return obj


# ----------------------------
# CLI / Breitenprogramm
# ----------------------------

def _load_mesh_json(mesh_json: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Erwartet JSON mit:
      vertices, faces, face_channel, radius (optional)
    """
    d = json.loads(mesh_json.read_text(encoding="utf-8"))
    V = np.array(d["vertices"], dtype=float)
    F = np.array(d["faces"], dtype=int)
    ch = np.array(d.get("face_channel", np.zeros(len(F), dtype=int)), dtype=int)
    r = float(d.get("radius", 1.0))
    return V, F, ch, r


def main():
    """
    Beispiel:
      python bm_qhe_model.py init --mesh level1_walter_trisection_refined_mesh.json --out run1 --R 1.0
      python bm_qhe_model.py run  --state run1 --Nphi 20 --dNw 5 --k 40
      python bm_qhe_model.py info --state run1
    """
    import argparse

    ap = argparse.ArgumentParser(prog="bm_qhe_model", description="BambergQHEModel: init/run/load/save")
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_init = sub.add_parser("init", help="Initialisiere State aus Mesh JSON und speichere serialisiert")
    ap_init.add_argument("--mesh", required=True, help="Pfad zu Mesh JSON (vertices/faces/face_channel)")
    ap_init.add_argument("--out", required=True, help="Zielordner für serialisierten Zustand")
    ap_init.add_argument("--R", type=float, default=1.0, help="Radius in Meter (oder BM-Skala -> später mappen)")
    ap_init.add_argument("--k", type=int, default=40, help="Anzahl Eigenwerte")
    ap_init.add_argument("--store-vecs", action="store_true", help="Eigenvektoren mit speichern")

    ap_run = sub.add_parser("run", help="Lade State, setze Parameter, berechne Spektren und speichere")
    ap_run.add_argument("--state", required=True, help="Ordner mit serialisiertem Zustand")
    ap_run.add_argument("--Nphi", type=int, default=0, help="Globaler Fluss N_phi")
    ap_run.add_argument("--dNw", type=int, default=0, help="Walter-Exzessfluss delta_N_phi_walter")
    ap_run.add_argument("--theta", type=float, default=2.0*np.pi/3.0, help="Morley-Grundphase")
    ap_run.add_argument("--k", type=int, default=40, help="Anzahl Eigenwerte")
    ap_run.add_argument("--weights", type=str, default="0.1,0.9", help="Kanalgewichte wW,wR (normalisiert)")
    ap_run.add_argument("--save", action="store_true", help="State nach Run wieder speichern (inkl. Cache)")

    ap_info = sub.add_parser("info", help="Zeige Infos aus State")
    ap_info.add_argument("--state", required=True)

    args = ap.parse_args()

    if args.cmd == "init":
        V,F,ch,r = _load_mesh_json(_as_path(args.mesh))
        model = BambergQHEModel(
            vertices=V, faces=F, face_channel=ch,
            radius_m=float(args.R),
            eig_k=int(args.k),
            store_eigenvectors=bool(args.store_vecs)
        )
        # initial caches can be computed optionally:
        model.spectrum_full()
        out = model.save(args.out, include_cache=True)
        print(f"[init] saved state to: {out}")

    elif args.cmd == "run":
        model = BambergQHEModel.load(args.state)
        model.eig_k = int(args.k)
        model.store_eigenvectors = False  # typically no need unless doing Chern now

        wW,wR = args.weights.split(",")
        model.weights = (float(wW), float(wR))

        model.theta_morley = float(args.theta)
        model.flux = (int(args.Nphi), int(args.dNw))

        # compute spectra
        sp0 = model.spectrum_full()["eigs"]
        spA = model.spectrum_magnetic()["eigs"]

        # quick energy view for first few modes
        E_scale = model.energy_scale()
        E0 = E_scale * sp0
        EA = E_scale * spA

        print("[run] params:", model.to_dict_meta()["params"])
        print("[run] first 10 eigs full:", sp0[:10])
        print("[run] first 10 eigs mag :", spA[:10])
        print("[run] first 5 energies (J) full:", E0[:5])
        print("[run] first 5 energies (J) mag :", EA[:5])

        if args.save:
            out = model.save(args.state, include_cache=True)
            print(f"[run] updated state saved to: {out}")

    elif args.cmd == "info":
        model = BambergQHEModel.load(args.state)
        meta = model.to_dict_meta()
        print(json.dumps(meta, indent=2))
        print("\nDerived:")
        print("  tau (s):", model.tau)
        print("  f_C (Hz):", model.si.compton_frequency)
        print("  R_K (Ohm):", model.si.von_klitzing)
        print("  Phi0=h/e (Wb):", model.si.flux_quantum_h_over_e)


if __name__ == "__main__":
    main()