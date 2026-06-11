from __future__ import annotations

import copy
import csv
import math
import random
from dataclasses import dataclass
from math import exp, cos, sin, pi
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import networkx as nx


# =========================================================
# 4D-Kugel / S3-Rand Patch
# =========================================================

COMPARE_EXPANSIONS = True
BOUNDARY_THRESHOLD = 0.60
PRINT_ONLY_FINAL_COMPARE = False
EXPANSION_MODE = "E4"   # "E3" oder "E4"

DEFAULT_SIM_SEED = 42
DEFAULT_SIM_STEPS = 15
DEFAULT_SIM_DT = 0.7

# Konfiguration pro Expansionsmodus (seed, steps, dt überschreibbar)
EXPANSION_MODES: Dict[str, Dict[str, Any]] = {
    "E3": {"seed": DEFAULT_SIM_SEED, "steps": DEFAULT_SIM_STEPS, "dt": DEFAULT_SIM_DT},
    "E4": {"seed": DEFAULT_SIM_SEED, "steps": DEFAULT_SIM_STEPS, "dt": DEFAULT_SIM_DT},
}

# Reihenfolge beim Modusvergleich (main-Schleife)
COMPARE_MODE_ORDER: Tuple[str, ...] = ("E3", "E4")

TRANSFER_ALPHA = 0.08
TRANSFER_SIGMA = 0.03
TRANSFER_GAMMA = 0.05
BOUNDARY_DECAY = 0.00
BULK_REPLENISH = 0.00

CLASS_ORDER = [
    "alpha_like",
    "beta_like",
    "gamma_like",
    "gamma_core",
    "bridge_like",
    "mixed",
    "boundary_like",
    "defect_like",
]


def expansion_params(mode: str) -> Tuple[int, float]:
    """
    Kindzahl und Radiusfaktor für den 4D-Fall.
    E4: vier identische Teilkugeln; E3: drei identische Teilkugeln.
    """
    mode_u = mode.upper().strip()
    if mode_u == "E4":
        child_count = 4
    elif mode_u == "E3":
        child_count = 3
    else:
        raise ValueError(f"Unbekannter Expansionsmodus: {mode}")
    radius_factor = child_count ** (-0.25)
    return child_count, radius_factor


def expansion_modes_for_current_run() -> Dict[str, Dict[str, Any]]:
    """
    Liefert das Dict mode_name -> mode_cfg für die Hauptschleife in main().
    """
    defaults: Dict[str, Any] = {
        "seed": DEFAULT_SIM_SEED,
        "steps": DEFAULT_SIM_STEPS,
        "dt": DEFAULT_SIM_DT,
    }
    if COMPARE_EXPANSIONS:
        return {
            k: {**defaults, **EXPANSION_MODES.get(k, {})}
            for k in COMPARE_MODE_ORDER
        }
    m = EXPANSION_MODE.upper().strip()
    return {m: {**defaults, **EXPANSION_MODES.get(m, {})}}


def _mode_boosts_for_expansion(mode: str) -> Tuple[float, float, float]:
    mode_u = mode.upper().strip()
    if mode_u == "E4":
        return 1.10, 0.95, 0.95
    if mode_u == "E3":
        return 0.95, 1.10, 1.15
    raise ValueError(f"Unbekannter Expansionsmodus: {mode}")


@dataclass
class Node:
    idx: int
    q1: float
    q2: float
    g: float
    rho: int
    ell: int
    q1_star: float
    q2_star: float
    sigma: float = 0.0
    p: float = 0.0
    bulk_info: float = 0.0
    boundary_info: float = 0.0
    transfer_flux: float = 0.0
    phase_tag: str = "bulk"
    sub_class: Optional[str] = None
    boundary_ratio: float = 0.0


def init_bulk_boundary(node: Node) -> None:
    """
    Initialisiert die neuen Bulk-/Randgrößen für einen Knoten.
    """
    node.bulk_info = 0.5 * node.q1 + 1.0 * node.q2 + 0.25 * node.g
    node.boundary_info = 0.0
    node.transfer_flux = 0.0
    node.phase_tag = "bulk"


def compute_boundary_transfer(
    node: Node,
    cls: str,
    mode_alpha_boost: float,
    mode_bridge_boost: float,
    mode_core_boost: float,
) -> float:
    """
    Berechnet den momentanen Bulk->S3-Transfer eines Knotens.
    """
    sigma = max(0.0, node.sigma)
    g = max(0.0, node.g)
    bulk = max(0.0, node.bulk_info)

    class_factor = {
        "alpha_like": 1.00,
        "bridge_like": 0.85,
        "gamma_core": 0.70,
        "gamma_like": 0.55,
        "beta_like": 0.35,
        "mixed": 0.20,
        "boundary_like": 0.50,
        "defect_like": 0.05,
    }.get(cls, 0.20)

    if cls == "alpha_like":
        class_factor *= mode_alpha_boost
    elif cls == "bridge_like":
        class_factor *= mode_bridge_boost
    elif cls == "gamma_core":
        class_factor *= mode_core_boost

    raw = bulk * class_factor * (
        TRANSFER_ALPHA + TRANSFER_SIGMA * sigma + TRANSFER_GAMMA * g
    )
    return max(0.0, raw)


def apply_boundary_projection(
    nodes: List[Node],
    class_labels: List[str],
    mode_alpha_boost: float,
    mode_bridge_boost: float,
    mode_core_boost: float,
) -> None:
    """
    Überträgt pro Schritt einen Teil der Bulk-Information auf den S3-Rand.
    """
    for node, cls in zip(nodes, class_labels):
        flux = compute_boundary_transfer(
            node, cls, mode_alpha_boost, mode_bridge_boost, mode_core_boost
        )
        flux = min(flux, node.bulk_info)

        node.transfer_flux = flux
        node.bulk_info -= flux
        node.boundary_info += flux

        if BOUNDARY_DECAY > 0.0:
            node.boundary_info *= (1.0 - BOUNDARY_DECAY)

        if BULK_REPLENISH > 0.0:
            node.bulk_info += BULK_REPLENISH

        if node.boundary_info > node.bulk_info:
            node.phase_tag = "boundary"
        elif flux > 1e-12:
            node.phase_tag = "transfer"
        else:
            node.phase_tag = "bulk"


def summarize_bulk_boundary(nodes: List[Node]) -> Tuple[float, float, float, float]:
    """
    Liefert globale Bulk-/Rand-Diagnosen.
    """
    bulk_total = sum(max(0.0, n.bulk_info) for n in nodes)
    boundary_total = sum(max(0.0, n.boundary_info) for n in nodes)
    transfer_total = sum(max(0.0, n.transfer_flux) for n in nodes)
    denom = bulk_total + boundary_total
    holo_eff = boundary_total / denom if denom > 0 else 0.0
    return bulk_total, boundary_total, transfer_total, holo_eff


def final_holographic_report(
    nodes: List[Node],
    expansion_mode: str,
    child_count: int,
    radius_factor: float,
) -> None:
    """
    Endbericht zur Bulk-/Rand-Dominanz.
    """
    bulk_nodes = []
    boundary_nodes = []
    transfer_nodes = []

    for i, n in enumerate(nodes):
        if n.phase_tag == "bulk":
            bulk_nodes.append(i)
        elif n.phase_tag == "boundary":
            boundary_nodes.append(i)
        else:
            transfer_nodes.append(i)

    print("\nHolographischer Endbefund")
    print("-" * 90)
    print(f"Expansionsmodus    : {expansion_mode}")
    print(f"Child count        : {child_count}")
    print(f"Radiusfaktor       : {radius_factor:.6f}")
    print(f"Bulk-dominant      : {bulk_nodes}")
    print(f"Boundary-dominant  : {boundary_nodes}")
    print(f"Transfer-aktiv     : {transfer_nodes}")


def holographic_summary(
    history: List[Dict[str, Any]],
    mode: str,
    child_count: int,
    radius_factor: float,
) -> None:
    print("\nHolographische Summary")
    print("-" * 110)
    print(f"Expansionsmodus : {mode}")
    print(f"Child count     : {child_count}")
    print(f"Radiusfaktor    : {radius_factor:.6f}")

    if not history:
        print("(keine Holo-History)")
        return

    h_values = [snap["H"] for snap in history]
    bulk_values = [snap["total_bulk"] for snap in history]
    s3_values = [snap["total_s3"] for snap in history]

    print(f"H(0)            : {h_values[0]:.6f}")
    print(f"H(final)        : {h_values[-1]:.6f}")
    print(f"Bulk(final)     : {bulk_values[-1]:.6f}")
    print(f"S3(final)       : {s3_values[-1]:.6f}")

    monotone = all(h_values[i + 1] >= h_values[i] - 1e-9 for i in range(len(h_values) - 1))
    print(f"H monotone      : {'ja' if monotone else 'nein'}")

    final_nodes = history[-1]["nodes"]
    boundary_dominant = [i for i, n in enumerate(final_nodes) if n.boundary_info > n.bulk_info]
    transfer_active = [i for i, n in enumerate(final_nodes) if n.transfer_flux > 0.05]
    bulk_dominant = [i for i, n in enumerate(final_nodes) if n.bulk_info >= n.boundary_info]

    print(f"Boundary-dominant: {boundary_dominant}")
    print(f"Transfer-aktiv   : {transfer_active}")
    print(f"Bulk-dominant    : {bulk_dominant}")

    if h_values[-1] > 0.66:
        phase = "stark holographische Endphase"
    elif h_values[-1] > 0.50:
        phase = "gemischt randdominierte Endphase"
    else:
        phase = "bulk-dominierte oder unvollständige Randphase"
    print(f"Endphase         : {phase}")


class ReactiveLatticeV2:
    def __init__(
        self,
        nodes: List[Node],
        edges: List[Tuple[int, int]],
        seed: int = 0,
        expansion_mode: Optional[str] = None,
    ):
        self.nodes: Dict[int, Node] = {n.idx: n for n in nodes}
        self.edges = edges
        self.neighbors: Dict[int, List[int]] = {i: [] for i in self.nodes}
        for i, j in edges:
            self.neighbors[i].append(j)
            self.neighbors[j].append(i)

        self.rng = random.Random(seed)

        self.expansion_mode = (expansion_mode or EXPANSION_MODE).upper().strip()
        self.child_count, self.radius_factor = expansion_params(self.expansion_mode)
        mb = _mode_boosts_for_expansion(self.expansion_mode)
        self.mode_alpha_boost, self.mode_bridge_boost, self.mode_core_boost = mb

        # Energien
        self.a1 = 1.0
        self.a2 = 1.4
        self.b_res = 0.6
        self.c_mis = 0.25

        # Spannung / Druck
        self.s1 = 1.0
        self.s2 = 1.2
        self.s3 = 0.7
        self.s4 = 0.20
        self.kappa1 = 1.0
        self.kappa2 = 0.35

        # Flüsse
        self.L1 = 0.10
        self.L2 = 0.07
        self.Lg = 0.22   # etwas stärker

        self.lambda1 = 0.35
        self.lambda2 = 0.55
        self.chi = 0.25

        # Reaktionen
        self.omega0 = 1.0
        self.beta = 1.2
        self.reaction_scale = 0.15

        # -------------------------------------------------
        # Patch V2: neue Stabilisierungsparameter
        # -------------------------------------------------
        self.empty_eps = 0.5        # Leere-Gewichtung der Sollwertspannung
        self.rad_loss = 0.06        # externe Photonabstrahlung pro Schritt
        self.deep_bind_boost = 1.8  # verstärkte Bindung in tiefer Zone

        # -------------------------------------------------
        # Patch V3: strengere Reaktionsbilanz
        # -------------------------------------------------
        self.emit_fraction = 0.85   # Anteil exothermer Energie, der in g landet
        self.absorb_efficiency = 1.0  # benötigte photonische Energie für endotherme Prozesse
        self.max_reaction_step_q1 = 2.0
        self.max_reaction_step_q2 = 1.0

        # -------------------------------------------------
        # Patch V4: Zustandsklassifikation
        # -------------------------------------------------
        self.class_history: List[Dict[str, int]] = []
        self.node_class_snapshots: List[Dict[int, str]] = []

        self.allowed_levels: Dict[int, List[Tuple[int, int]]] = {}
        for i, n in self.nodes.items():
            self.allowed_levels[i] = self._build_local_resonance_levels(n)

        self.history: List[Dict[str, float]] = []
        self.snapshots: List[Dict[int, Dict[str, Any]]] = []
        self.holo_history: List[Dict[str, Any]] = []

        for n in self.nodes.values():
            init_bulk_boundary(n)

    def nodes_in_index_order(self) -> List[Node]:
        return [self.nodes[i] for i in sorted(self.nodes)]

    def apply_holographic_projection(self) -> Tuple[float, float, float, float]:
        nodes_ord = self.nodes_in_index_order()
        labels: List[str] = []
        for n in nodes_ord:
            final_cls, sub, _ = self.classify_node_detailed(n)
            transfer_cls = sub if final_cls == "boundary_like" and sub is not None else final_cls
            labels.append(transfer_cls)
        apply_boundary_projection(
            nodes_ord,
            labels,
            self.mode_alpha_boost,
            self.mode_bridge_boost,
            self.mode_core_boost,
        )
        return summarize_bulk_boundary(nodes_ord)

    def _build_local_resonance_levels(self, node: Node) -> List[Tuple[int, int]]:
        base = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1)]
        if node.ell <= 6:
            base += [(1, 1), (0, 1), (2, 1)]
        elif node.ell <= 12:
            base += [(2, 1), (3, 1), (1, 2)]
        elif node.ell <= 18:
            base += [(2, 2), (3, 2), (1, 3)]
        else:
            base += [(4, 2), (2, 4), (3, 3)]
        return list(dict.fromkeys(base))

    def resonance_penalty(self, node: Node) -> float:
        return min((node.q1 - u) ** 2 + (node.q2 - v) ** 2 for (u, v) in self.allowed_levels[node.idx])

    def local_energy(self, node: Node) -> float:
        e1 = self.a1 * (node.q1 - node.q1_star) ** 2
        e2 = self.a2 * (node.q2 - node.q2_star) ** 2
        eg = self.omega0 * node.g
        eres = self.b_res * self.resonance_penalty(node)

        emis = 0.0
        for j in self.neighbors[node.idx]:
            nj = self.nodes[j]
            emis += self.c_mis * (
                0.8 * (node.q1 - nj.q1) ** 2 +
                1.0 * (node.q2 - nj.q2) ** 2
            )
        return e1 + e2 + eg + eres + emis

    # =====================================================
    # Patch V2: Hilfsfunktionen
    # =====================================================

    def mass_like(self, node: Node) -> float:
        """
        Effektive lokale Besetzungsmasse:
        Einser zählen 1, Zweier zählen 2.
        """
        return node.q1 + 2.0 * node.q2

    def emptiness_weight(self, node: Node) -> float:
        """
        Fast leere Knoten sollen nicht künstlich durch Sollwertabweichung
        hochgespannt bleiben.
        """
        m = self.mass_like(node)
        return m / (m + self.empty_eps)

    def is_deep_sink(self, node: Node) -> bool:
        """
        Tiefe Zone mit bevorzugter Bindung.
        """
        return (node.rho, node.ell) in {(6, 6), (12, 6), (12, 12)}

    # =====================================================
    # Patch V4: heuristische alpha/beta/gamma-Klassifikation
    # =====================================================

    def _classify_base_local(self, node: Node) -> str:
        """Physikalische Basisklasse (ohne boundary_like-Aufsplitting)."""
        m = self.mass_like(node)
        photon_ratio = node.g / (m + 1e-9)
        res_pen = self.resonance_penalty(node)

        if m < 0.35:
            return "defect_like"

        if m < 1.0 and node.sigma > 1.0:
            return "defect_like"

        if (
            node.rho >= 100
            and node.ell >= 100
            and node.g > 2.0
            and node.q2 > 1.0
            and node.sigma < 1.5
        ):
            return "gamma_core"

        if photon_ratio > 1.25 and node.g > 0.8:
            return "gamma_like"

        if node.g > 2.0 and node.sigma < 1.5:
            return "gamma_like"

        if (
            18 <= node.rho <= 30
            and 18 <= node.ell <= 30
            and 0.45 <= node.sigma <= 1.2
            and node.q2 > 1.1
            and node.g < 0.6
        ):
            return "bridge_like"

        if self.is_deep_sink(node) and node.sigma < 0.9 and node.q2 > 0.9:
            return "alpha_like"

        if node.sigma < 0.45 and node.q2 > 1.0 and node.g < 0.6 and res_pen < 0.35:
            return "alpha_like"

        if (node.rho >= 18 or node.ell >= 18) and node.q2 > 1.0 and node.sigma >= 0.45:
            return "beta_like"

        if node.rho >= 100 and node.g > 1.0:
            return "beta_like"

        return "mixed"

    def classify_node_detailed(self, node: Node) -> Tuple[str, Optional[str], float]:
        """
        Klassifikation mit boundary_like-Feinklasse.
        Rückgabe: (node_class, sub_class, boundary_ratio)
        """
        bulk = max(0.0, node.bulk_info)
        s3 = max(0.0, node.boundary_info)
        total_geom = bulk + s3
        boundary_ratio = (s3 / total_geom) if total_geom > 1e-12 else 0.0
        boundary_flag = boundary_ratio >= BOUNDARY_THRESHOLD and s3 > bulk

        base = self._classify_base_local(node)

        if boundary_flag and base in {
            "alpha_like", "bridge_like", "gamma_like", "gamma_core",
        }:
            node.sub_class = base
            node.boundary_ratio = boundary_ratio
            return "boundary_like", base, boundary_ratio

        node.sub_class = None
        node.boundary_ratio = boundary_ratio
        return base, None, boundary_ratio

    def classify_node(self, node: Node) -> str:
        return self.classify_node_detailed(node)[0]

    def classify_all_nodes(self) -> Dict[int, str]:
        return {i: self.classify_node(n) for i, n in self.nodes.items()}

    def class_counts(self) -> Dict[str, int]:
        counts = {k: 0 for k in CLASS_ORDER}
        for cls in self.classify_all_nodes().values():
            counts[cls] += 1
        return counts

    def print_class_timeline(self) -> None:
        """
        Kompakte Zeitreihe der Klassenhäufigkeiten.
        """
        print("\nKlassen-Zeitreihe")
        print("-" * 128)
        print(
            f"{'Schritt':>6} | {'alpha':>5} | {'beta':>4} | {'gamma':>5} | "
            f"{'g_core':>6} | {'bridge':>6} | {'mixed':>5} | {'bnd':>4} | {'defect':>6}"
        )
        print("-" * 128)
        for step_no, counts in enumerate(self.class_history):
            print(
                f"{step_no:6d} | "
                f"{counts['alpha_like']:5d} | "
                f"{counts['beta_like']:4d} | "
                f"{counts['gamma_like']:5d} | "
                f"{counts['gamma_core']:6d} | "
                f"{counts['bridge_like']:6d} | "
                f"{counts['mixed']:5d} | "
                f"{counts['boundary_like']:4d} | "
                f"{counts['defect_like']:6d}"
            )

    def print_class_transition_matrix(self) -> None:
        classes = list(CLASS_ORDER)
        trans = {(a, b): 0 for a in classes for b in classes}

        for t in range(len(self.node_class_snapshots) - 1):
            c0 = self.node_class_snapshots[t]
            c1 = self.node_class_snapshots[t + 1]
            for i in c0:
                trans[(c0[i], c1[i])] += 1

        print("\nKlassen-Übergangsmatrix")
        width = 14 + len(classes) * 15
        print("-" * width)
        header = f"{'von \\ nach':>14}"
        for b in classes:
            header += f" | {b:>12}"
        print(header)
        print("-" * width)

        for a in classes:
            row = f"{a:>14}"
            for b in classes:
                row += f" | {trans[(a, b)]:12d}"
            print(row)

        print("\nNichttriviale Übergänge")
        print("-" * width)
        for a in classes:
            for b in classes:
                if a != b and trans[(a, b)] > 0:
                    print(f"{a:>14} -> {b:<12} : {trans[(a, b)]}")

    def node_role_label(self, idx: int, cls: str) -> str:
        node = self.nodes[idx]

        if cls == "beta_like":
            if node.rho >= 100:
                return "oberer beta-Restkern"
            return "beta-Zone"

        if cls == "gamma_core":
            return "gamma-Kern / Vermittlungskern"

        if cls == "gamma_like":
            return "gamma-Zustand / photonische Anregung"

        if cls == "bridge_like":
            return "Brückenknoten / Übergangszone"

        if cls == "alpha_like":
            if self.is_deep_sink(node):
                return "tiefe alpha-Bindungszone"
            if node.rho <= 30 and node.ell <= 30:
                return "mittlere alpha-Ordnungszone"
            return "alpha-Zone"

        if cls == "mixed":
            return "gemischter Übergangszustand"

        if cls == "defect_like":
            return "Defektzustand"

        if cls == "boundary_like":
            base = self.nodes[idx].sub_class or "?"
            return f"Rand-dominanter Feinkern (Basis: {base})"

        return "unbestimmt"

    def print_node_role_table(self) -> None:
        final_classes = self.classify_all_nodes()

        print("\nKnoten-Rollen-Tabelle (letzter Schritt)")
        print("-" * 118)
        print(
            f"{'Knoten':>6} | {'rho':>4} | {'ell':>4} | {'Klasse':>14} | {'Rolle':>36} | "
            f"{'q1':>6} | {'q2':>6} | {'g':>6} | {'sigma':>7}"
        )
        print("-" * 118)

        for i in sorted(self.nodes):
            n = self.nodes[i]
            cls = final_classes[i]
            role = self.node_role_label(i, cls)
            print(
                f"{i:6d} | {n.rho:4d} | {n.ell:4d} | {cls:>14} | {role:>36} | "
                f"{n.q1:6.3f} | {n.q2:6.3f} | {n.g:6.3f} | {n.sigma:7.3f}"
            )

    def print_structural_report(self) -> None:
        counts = self.class_counts()
        final_classes = self.classify_all_nodes()

        beta_nodes = [i for i, c in final_classes.items() if c == "beta_like"]
        gamma_core_nodes = [i for i, c in final_classes.items() if c == "gamma_core"]
        bridge_nodes = [i for i, c in final_classes.items() if c == "bridge_like"]
        alpha_nodes = [i for i, c in final_classes.items() if c == "alpha_like"]

        print("\nStrukturbefund")
        print("-" * 118)
        print(
            f"Endkonfiguration: alpha={counts['alpha_like']}, beta={counts['beta_like']}, "
            f"gamma={counts['gamma_like']}, gamma_core={counts['gamma_core']}, "
            f"bridge={counts['bridge_like']}, mixed={counts['mixed']}, "
            f"boundary={counts['boundary_like']}, defect={counts['defect_like']}"
        )

        boundary_label_nodes = [i for i, c in final_classes.items() if c == "boundary_like"]
        if boundary_label_nodes:
            print(f"Boundary-like       : {boundary_label_nodes}")

        if beta_nodes:
            print(f"Beta-Restkerne        : {beta_nodes}")
        if gamma_core_nodes:
            print(f"Gamma-Kern            : {gamma_core_nodes}")
        if bridge_nodes:
            print(f"Brückenknoten         : {bridge_nodes}")
        if alpha_nodes:
            print(f"Alpha-Bindungszone    : {alpha_nodes}")

        if counts["mixed"] == 0 and counts["defect_like"] == 0:
            print("Interpretation        : keine diffusen oder defekten Restzustände im Endschritt.")

        if counts["gamma_core"] == 1 and counts["beta_like"] == 2:
            print("Muster                : stabiler einzelner Vermittlungskern mit zwei beta-Restkernen.")

        if counts["alpha_like"] >= 4:
            print("Bindungsphase         : alpha-dominante Endphase ist klar ausgebildet.")

    def update_sigma_and_pressure(self) -> None:
        sigmas: Dict[int, float] = {}
        for i, node in self.nodes.items():
            local_res = self.resonance_penalty(node)
            misfit = 0.0
            for j in self.neighbors[i]:
                nj = self.nodes[j]
                misfit += (node.q1 - nj.q1) ** 2 + (node.q2 - nj.q2) ** 2

            # -------------------------------------------------
            # Patch V2: Leere-Gewichtung der Sollwertspannung
            # -------------------------------------------------
            w_occ = self.emptiness_weight(node)

            sigma = (
                w_occ * (
                    self.s1 * (node.q1 - node.q1_star) ** 2 +
                    self.s2 * (node.q2 - node.q2_star) ** 2
                )
                + self.s3 * local_res
                + self.s4 * misfit
            )

            # -------------------------------------------------
            # Patch V2: tiefe Senken werden bei resonanznaher
            # Besetzung zusätzlich entlastet
            # -------------------------------------------------
            if self.is_deep_sink(node) and local_res < 0.35:
                sigma *= 0.75

            sigmas[i] = sigma

        for i, node in self.nodes.items():
            node.sigma = sigmas[i]

        for i, node in self.nodes.items():
            grad_part = sum(sigmas[i] - sigmas[j] for j in self.neighbors[i])
            node.p = self.kappa1 * node.sigma + self.kappa2 * grad_part

    def resonance_bias_q1(self, i: int, j: int) -> float:
        ni = self.nodes[i]
        nj = self.nodes[j]
        return 1.0 / (1.0 + abs(ni.ell - nj.ell))

    def resonance_bias_q2(self, i: int, j: int) -> float:
        ni = self.nodes[i]
        nj = self.nodes[j]
        dr = abs(ni.rho - nj.rho)
        de = abs(ni.ell - nj.ell)
        return exp(-dr / 40.0) * exp(-de / 12.0)

    def compute_fluxes(self):
        J1, J2, Jg = {}, {}, {}

        for i, j in self.edges:
            ni = self.nodes[i]
            nj = self.nodes[j]

            drive1_ij = (
                (ni.p - nj.p)
                + self.lambda1 * (ni.sigma - nj.sigma)
                + 0.6 * self.resonance_bias_q1(i, j)
            )
            drive1_ji = (
                (nj.p - ni.p)
                + self.lambda1 * (nj.sigma - ni.sigma)
                + 0.6 * self.resonance_bias_q1(j, i)
            )

            drive2_ij = (
                (ni.p - nj.p)
                + self.lambda2 * (ni.sigma - nj.sigma)
                + 1.0 * self.resonance_bias_q2(i, j)
            )
            drive2_ji = (
                (nj.p - ni.p)
                + self.lambda2 * (nj.sigma - ni.sigma)
                + 1.0 * self.resonance_bias_q2(j, i)
            )

            driveg_ij = (ni.sigma - nj.sigma) + self.chi * (ni.p - nj.p)
            driveg_ji = (nj.sigma - ni.sigma) + self.chi * (nj.p - ni.p)

            J1[(i, j)] = max(0.0, self.L1 * drive1_ij)
            J1[(j, i)] = max(0.0, self.L1 * drive1_ji)

            J2[(i, j)] = max(0.0, self.L2 * drive2_ij)
            J2[(j, i)] = max(0.0, self.L2 * drive2_ji)

            Jg[(i, j)] = max(0.0, self.Lg * driveg_ij)
            Jg[(j, i)] = max(0.0, self.Lg * driveg_ji)

        return J1, J2, Jg

    def apply_fluxes(self, dt: float) -> None:
        J1, J2, Jg = self.compute_fluxes()

        dq1 = {i: 0.0 for i in self.nodes}
        dq2 = {i: 0.0 for i in self.nodes}
        dg = {i: 0.0 for i in self.nodes}

        for (i, j), f in J1.items():
            amount = min(dt * f, self.nodes[i].q1)
            dq1[i] -= amount
            dq1[j] += amount

        for (i, j), f in J2.items():
            amount = min(dt * f, self.nodes[i].q2)
            dq2[i] -= amount
            dq2[j] += amount

        for (i, j), f in Jg.items():
            amount = min(dt * f, self.nodes[i].g)
            dq1_i = 0
            dg[i] -= amount
            dg[j] += amount

        for i, node in self.nodes.items():
            node.q1 = max(0.0, node.q1 + dq1[i])
            node.q2 = max(0.0, node.q2 + dq2[i])
            node.g = max(0.0, node.g + dg[i])

        # -------------------------------------------------
        # Patch V2: externe Photonabstrahlung
        # -------------------------------------------------
        for node in self.nodes.values():
            node.g = max(0.0, node.g * (1.0 - self.rad_loss))

    def _reaction_acceptance(self, dE: float) -> float:
        if dE <= 0:
            return 1.0
        return exp(-self.beta * dE)

    def try_local_reaction(self, i: int) -> None:
        node = self.nodes[i]

        candidates = [
            ("q1_decay",        -1.0,  0.0, +0.8),
            ("q2_decay",         0.0, -1.0, +1.6),
            ("bind_2q1_to_q2",  -2.0, +1.0, -0.6),
            ("split_q2_to_2q1", +2.0, -1.0, +0.4),
            ("q1_from_g",       +1.0,  0.0, -1.0),
            ("q2_from_g",        0.0, +1.0, -2.0),
        ]

        self.rng.shuffle(candidates)

        for name, dq1, dq2, dg in candidates:
            if node.q1 + dq1 < 0:
                continue
            if node.q2 + dq2 < 0:
                continue
            if node.g + dg < 0:
                continue

            old_q1, old_q2, old_g = node.q1, node.q2, node.g
            old_sigma, old_p = node.sigma, node.p
            old_pen = self.resonance_penalty(node)
            E_old = self.local_energy(node)

            node.q1 += dq1
            node.q2 += dq2
            node.g += dg

            self.update_sigma_and_pressure()
            E_new = self.local_energy(node)
            new_pen = self.resonance_penalty(node)

            dE = E_new - E_old

            # -------------------------------------------------
            # Patch V3: explizite Energiedeckung
            # -------------------------------------------------
            # endotherm: braucht lokale photonische Energie
            if dE > 0:
                needed = self.absorb_efficiency * dE
                if old_g < needed:
                    node.q1, node.q2, node.g = old_q1, old_q2, old_g
                    node.sigma, node.p = old_sigma, old_p
                    self.update_sigma_and_pressure()
                    continue

            acc = self._reaction_acceptance(dE)

            if name == "q2_decay":
                acc *= min(1.0, 0.3 + 0.08 * old_sigma)

            elif name == "q1_decay":
                acc *= min(1.0, 0.4 + 0.06 * old_sigma)

            elif name == "bind_2q1_to_q2":
                if new_pen < old_pen:
                    acc *= 1.2
                # ---------------------------------------------
                # Patch V2: tiefe Zone bindet bevorzugt
                # ---------------------------------------------
                if self.is_deep_sink(node):
                    acc *= self.deep_bind_boost

            elif name == "split_q2_to_2q1":
                # ---------------------------------------------
                # Patch V2: in tiefer Zone Dissoziation hemmen
                # ---------------------------------------------
                if self.is_deep_sink(node):
                    acc *= 0.55

            elif name in ("q1_from_g", "q2_from_g"):
                acc *= 0.75

            acc = max(0.0, min(1.0, acc * self.reaction_scale * 5.0))

            if self.rng.random() < acc:
                # ---------------------------------------------
                # Patch V3: Energie sauber verbuchen
                # ---------------------------------------------
                if dE > 0:
                    # endotherme Reaktion bezahlt aus lokalem Photonbad
                    node.g = max(0.0, node.g - self.absorb_efficiency * dE)
                elif dE < 0:
                    # exotherme Reaktion speist Photonenergie ein
                    node.g += self.emit_fraction * (-dE)

                self.update_sigma_and_pressure()
                return

            node.q1, node.q2, node.g = old_q1, old_q2, old_g
            node.sigma, node.p = old_sigma, old_p
            self.update_sigma_and_pressure()

    def apply_reactions(self) -> None:
        order = list(self.nodes.keys())
        self.rng.shuffle(order)
        for i in order:
            self.try_local_reaction(i)

    def totals(self) -> Dict[str, float]:
        return {
            "q1": sum(n.q1 for n in self.nodes.values()),
            "q2": sum(n.q2 for n in self.nodes.values()),
            "g": sum(n.g for n in self.nodes.values()),
            "sigma": sum(n.sigma for n in self.nodes.values()),
            "energy": sum(self.local_energy(n) for n in self.nodes.values()),
        }

    def store_snapshot(self) -> None:
        t = self.totals()
        self.history.append(t)
        node_classes = self.classify_all_nodes()
        nodes_ord = self.nodes_in_index_order()
        bulk_total, boundary_total, transfer_total, holo_eff = summarize_bulk_boundary(nodes_ord)
        step_idx = len(self.class_history)
        self.holo_history.append({
            "step": step_idx,
            "total_q1": t["q1"],
            "total_q2": t["q2"],
            "total_g": t["g"],
            "total_sigma": t["sigma"],
            "total_E": t["energy"],
            "total_bulk": bulk_total,
            "total_s3": boundary_total,
            "total_flux": transfer_total,
            "H": holo_eff,
            "nodes": copy.deepcopy(nodes_ord),
        })

        self.snapshots.append({
            i: {
                "q1": n.q1, "q2": n.q2, "g": n.g,
                "sigma": n.sigma, "p": n.p,
                "rho": n.rho, "ell": n.ell,
                "class": node_classes[i],
                "sub_class": n.sub_class,
                "boundary_ratio": n.boundary_ratio,
                "bulk_info": n.bulk_info,
                "boundary_info": n.boundary_info,
                "transfer_flux": n.transfer_flux,
                "phase_tag": n.phase_tag,
            }
            for i, n in self.nodes.items()
        })
        cts: Dict[str, int] = {k: 0 for k in CLASS_ORDER}
        for cl in node_classes.values():
            cts[cl] += 1
        self.class_history.append(cts)
        self.node_class_snapshots.append(node_classes)

    def step(self, dt: float = 1.0) -> None:
        self.update_sigma_and_pressure()
        self.apply_fluxes(dt)
        self.update_sigma_and_pressure()
        self.apply_reactions()
        self.update_sigma_and_pressure()
        self.apply_holographic_projection()
        self.store_snapshot()

    def print_state(self, step_no: int) -> None:
        t = self.totals()
        c = self.class_counts()
        nodes_ord = self.nodes_in_index_order()
        bulk_total, boundary_total, transfer_total, holo_eff = summarize_bulk_boundary(nodes_ord)
        print(f"\n=== Schritt {step_no} ===")
        print(
            f"Total q1={t['q1']:.3f} | q2={t['q2']:.3f} | "
            f"g={t['g']:.3f} | sigma={t['sigma']:.3f} | E={t['energy']:.3f} "
            f"| Bulk={bulk_total:.3f} | S3={boundary_total:.3f} "
            f"| Flux={transfer_total:.3f} | H={holo_eff:.3f}"
        )
        print(
            f"Klassen | alpha={c['alpha_like']} | beta={c['beta_like']} | "
            f"gamma={c['gamma_like']} | gamma_core={c['gamma_core']} | "
            f"bridge={c['bridge_like']} | mixed={c['mixed']} | "
            f"boundary={c['boundary_like']} | defect={c['defect_like']}"
        )
        print("-" * 145)
        for i in sorted(self.nodes):
            n = self.nodes[i]
            cls, sub, br = self.classify_node_detailed(n)
            subtxt = f" ({sub})" if sub is not None else ""
            print(
                f"Knoten {i:2d} | rho={n.rho:3d} ell={n.ell:3d} | "
                f"q1={n.q1:6.3f} q2={n.q2:6.3f} g={n.g:6.3f} | "
                f"sigma={n.sigma:7.3f} p={n.p:7.3f} | "
                f"bulk={n.bulk_info:6.3f} s3={n.boundary_info:6.3f} "
                f"| flux={n.transfer_flux:6.3f} | "
                f"{cls}{subtxt} | br={br:5.3f} | {n.phase_tag}"
            )


# Alias wie im Patch-Diff (Klasse heißt intern ReactiveLatticeV2)
ReactiveLattice = ReactiveLatticeV2

OUTPUT_DIR = Path("holo_output")
OUTPUT_DIR.mkdir(exist_ok=True)


def ensure_mode_store(mode_store: Dict[str, Any], mode_name: str) -> None:
    if mode_name not in mode_store:
        mode_store[mode_name] = {
            "summary_rows": [],
            "node_rows": [],
        }


def append_summary_row(
    mode_store: Dict[str, Any],
    mode_name: str,
    step: int,
    *,
    totals: Dict[str, float],
    class_counts: Dict[str, int],
    holo: Dict[str, float],
) -> None:
    ensure_mode_store(mode_store, mode_name)

    row = {
        "mode": mode_name,
        "step": step,
        "q1_total": totals.get("q1", 0.0),
        "q2_total": totals.get("q2", 0.0),
        "g_total": totals.get("g", 0.0),
        "sigma_total": totals.get("sigma", 0.0),
        "E_total": totals.get("E", totals.get("energy", 0.0)),
        "bulk_total": holo.get("bulk_total", 0.0),
        "s3_total": holo.get("s3_total", 0.0),
        "flux_total": holo.get("flux_total", 0.0),
        "H": holo.get("H", 0.0),
        "alpha": class_counts.get("alpha_like", 0),
        "beta": class_counts.get("beta_like", 0),
        "gamma": class_counts.get("gamma_like", 0),
        "gamma_core": class_counts.get("gamma_core", 0),
        "bridge": class_counts.get("bridge_like", 0),
        "mixed": class_counts.get("mixed", 0),
        "boundary": class_counts.get("boundary_like", 0),
        "defect": class_counts.get("defect_like", 0),
    }
    mode_store[mode_name]["summary_rows"].append(row)


def append_node_rows(
    mode_store: Dict[str, Any],
    mode_name: str,
    step: int,
    nodes: List[Dict[str, Any]],
) -> None:
    """
    nodes: Zeilen ohne mode/step. Kanonisches Schema pro Eintrag (s. build_node_export_rows):

        node = {
            "node": int,      # Knotenindex (zusätzlich zu rho/ell)
            "rho": ...,
            "ell": ...,
            "q1": ...,
            "q2": ...,
            "g": ...,
            "sigma": ...,
            "p": ...,
            "class": ...,
            "role": ...,
            "bulk": ...,
            "s3": ...,
            "flux": ...,
            "br": ...,
            "region": ...,
        }
    """
    ensure_mode_store(mode_store, mode_name)

    for row in nodes:
        mode_store[mode_name]["node_rows"].append(
            {**row, "mode": mode_name, "step": step}
        )


def build_node_export_rows(
    lattice: ReactiveLatticeV2,
    node_list: List[Node],
) -> List[Dict[str, Any]]:
    """
    Eine Exportzeile pro Knoten; Schlüsselreihenfolge wie im kanonischen node-Dict.
    """
    rows: List[Dict[str, Any]] = []
    for node in node_list:
        cls = lattice.classify_node(node)
        role = lattice.node_role_label(node.idx, cls)
        _, _, br = lattice.classify_node_detailed(node)
        row = {
            "node": node.idx,
            "rho": node.rho,
            "ell": node.ell,
            "q1": node.q1,
            "q2": node.q2,
            "g": node.g,
            "sigma": node.sigma,
            "p": node.p,
            "class": cls,
            "role": role,
            "bulk": node.bulk_info,
            "s3": node.boundary_info,
            "flux": node.transfer_flux,
            "br": br,
            "region": node.phase_tag,
        }
        rows.append(row)
    return rows


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def print_holo_output_files() -> None:
    print("\nDateien geschrieben nach:", OUTPUT_DIR.resolve())
    print(" - holo_summary.csv")
    print(" - holo_nodes.csv")
    print(" - plot_H.png")
    print(" - plot_bulk_s3.png")
    print(" - plot_boundary_counts.png")


def export_mode_csvs(mode_store: Dict[str, Any]) -> None:
    all_summary: List[Dict[str, Any]] = []
    all_nodes: List[Dict[str, Any]] = []

    for mode_name, payload in mode_store.items():
        all_summary.extend(payload["summary_rows"])
        all_nodes.extend(payload["node_rows"])

        write_csv(OUTPUT_DIR / f"holo_summary_{mode_name}.csv", payload["summary_rows"])
        write_csv(OUTPUT_DIR / f"holo_nodes_{mode_name}.csv", payload["node_rows"])

    write_csv(OUTPUT_DIR / "holo_summary.csv", all_summary)
    write_csv(OUTPUT_DIR / "holo_nodes.csv", all_nodes)


def make_plots(mode_store: Dict[str, Any]) -> None:
    if not mode_store:
        return
    if not any(payload.get("summary_rows") for payload in mode_store.values()):
        return
    plot_H(mode_store)
    plot_bulk_s3(mode_store)
    plot_boundary_counts(mode_store)


def plot_H(mode_store: Dict[str, Any]) -> None:
    plt.figure(figsize=(10, 6))
    for mode_name, payload in mode_store.items():
        rows = payload["summary_rows"]
        if not rows:
            continue
        xs = [r["step"] for r in rows]
        ys = [r["H"] for r in rows]
        plt.plot(xs, ys, marker="o", label=mode_name)

    plt.xlabel("Schritt")
    plt.ylabel("H")
    plt.title("Holographische Entropie H(t)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "plot_H.png", dpi=180)
    plt.close()


def plot_bulk_s3(mode_store: Dict[str, Any]) -> None:
    for mode_name, payload in mode_store.items():
        rows = payload["summary_rows"]
        if not rows:
            continue
        xs = [r["step"] for r in rows]
        bulk = [r["bulk_total"] for r in rows]
        s3 = [r["s3_total"] for r in rows]

        plt.figure(figsize=(10, 6))
        plt.plot(xs, bulk, marker="o", label=f"{mode_name} Bulk")
        plt.plot(xs, s3, marker="o", label=f"{mode_name} S3")
        plt.xlabel("Schritt")
        plt.ylabel("Anteil")
        plt.title(f"Bulk/S3-Entwicklung ({mode_name})")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / f"plot_bulk_s3_{mode_name}.png", dpi=180)
        plt.close()

    # gemeinsamer Direktvergleich
    plt.figure(figsize=(10, 6))
    for mode_name, payload in mode_store.items():
        rows = payload["summary_rows"]
        if not rows:
            continue
        xs = [r["step"] for r in rows]
        s3 = [r["s3_total"] for r in rows]
        plt.plot(xs, s3, marker="o", label=f"{mode_name} S3")

    plt.xlabel("Schritt")
    plt.ylabel("S3")
    plt.title("Vergleich der S3-Randakkumulation")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "plot_bulk_s3.png", dpi=180)
    plt.close()


def plot_boundary_counts(mode_store: Dict[str, Any]) -> None:
    plt.figure(figsize=(10, 6))
    for mode_name, payload in mode_store.items():
        rows = payload["summary_rows"]
        if not rows:
            continue
        xs = [r["step"] for r in rows]
        ys = [r["boundary"] for r in rows]
        plt.plot(xs, ys, marker="o", label=mode_name)

    plt.xlabel("Schritt")
    plt.ylabel("boundary_like Knoten")
    plt.title("Randdominanz im Zeitverlauf")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "plot_boundary_counts.png", dpi=180)
    plt.close()


def record_holo_csv_step(
    mode_store: Dict[str, Any],
    mode_name: str,
    model: ReactiveLatticeV2,
    step: int,
) -> None:
    t = model.totals()
    totals = {
        "q1": t["q1"],
        "q2": t["q2"],
        "g": t["g"],
        "sigma": t["sigma"],
        "E": t["energy"],
    }
    class_counts = model.class_counts()
    nodes_ord = model.nodes_in_index_order()
    bulk_total, s3_total, flux_total, H_value = summarize_bulk_boundary(nodes_ord)
    holo = {
        "bulk_total": bulk_total,
        "s3_total": s3_total,
        "flux_total": flux_total,
        "H": H_value,
    }
    append_summary_row(
        mode_store,
        mode_name,
        step,
        totals=totals,
        class_counts=class_counts,
        holo=holo,
    )
    nodes = build_node_export_rows(model, nodes_ord)
    append_node_rows(mode_store, mode_name, step, nodes)


def init_nodes_standard_ring(
    expansion_mode: Optional[str] = None,
    seed: int = DEFAULT_SIM_SEED,
) -> ReactiveLatticeV2:
    """
    Standard-Testgitter (Ring); entspricht init_nodes(...) im Pseudocode.
    """
    nodes = [
        Node(0, q1=2.0, q2=3.0, g=1.0, rho=198, ell=102, q1_star=1.0, q2_star=1.0),
        Node(1, q1=1.0, q2=2.5, g=0.5, rho=192, ell=108, q1_star=1.0, q2_star=1.0),
        Node(2, q1=1.0, q2=2.0, g=0.2, rho=102, ell=102, q1_star=1.0, q2_star=1.0),
        Node(3, q1=1.5, q2=1.5, g=0.2, rho=18,  ell=18,  q1_star=1.0, q2_star=1.0),
        Node(4, q1=1.2, q2=1.1, g=0.1, rho=12,  ell=12,  q1_star=1.0, q2_star=1.0),
        Node(5, q1=1.1, q2=0.8, g=0.1, rho=12,  ell=6,   q1_star=1.0, q2_star=1.0),
        Node(6, q1=0.9, q2=0.6, g=0.3, rho=6,   ell=6,   q1_star=1.0, q2_star=1.0),
        Node(7, q1=0.8, q2=0.7, g=0.4, rho=6,   ell=6,   q1_star=1.0, q2_star=1.0),
    ]
    edges = [(i, (i + 1) % len(nodes)) for i in range(len(nodes))]
    model = ReactiveLatticeV2(nodes, edges, seed=seed, expansion_mode=expansion_mode)
    model.update_sigma_and_pressure()
    model.store_snapshot()
    return model


def build_test_ring_v2(
    expansion_mode: Optional[str] = None,
    seed: int = 42,
) -> ReactiveLatticeV2:
    return init_nodes_standard_ring(expansion_mode=expansion_mode, seed=seed)


def run_one_mode_into_store(
    mode_store: Dict[str, Any],
    mode_name: str,
    mode_cfg: Dict[str, Any],
    *,
    do_print: bool,
) -> Dict[str, Any]:
    """
    Ein Expansionsmodus: Schleife über Schritte, pro Schritt
    Dynamik -> Klassen/Holo -> record_holo_csv_step (append_summary_row + append_node_rows).
    """
    seed = int(mode_cfg.get("seed", DEFAULT_SIM_SEED))
    num_steps = int(mode_cfg.get("steps", DEFAULT_SIM_STEPS))
    dt = float(mode_cfg.get("dt", DEFAULT_SIM_DT))

    model = init_nodes_standard_ring(expansion_mode=mode_name, seed=seed)
    mode_key = model.expansion_mode

    record_holo_csv_step(mode_store, mode_key, model, 0)
    if do_print:
        model.print_state(0)

    for step in range(1, num_steps + 1):
        model.step(dt)
        record_holo_csv_step(mode_store, mode_key, model, step)
        if do_print:
            model.print_state(step)

    if do_print:
        model.print_class_timeline()
        model.print_class_transition_matrix()
        model.print_structural_report()
        model.print_node_role_table()

    final_holographic_report(
        model.nodes_in_index_order(),
        model.expansion_mode,
        model.child_count,
        model.radius_factor,
    )
    holographic_summary(
        model.holo_history,
        model.expansion_mode,
        model.child_count,
        model.radius_factor,
    )

    return {
        "mode": model.expansion_mode,
        "child_count": model.child_count,
        "radius_factor": model.radius_factor,
        "history": model.holo_history,
        "final_H": model.holo_history[-1]["H"],
        "final_bulk": model.holo_history[-1]["total_bulk"],
        "final_s3": model.holo_history[-1]["total_s3"],
    }


def run_simulation(
    expansion_mode: str = "E4",
    seed: int = 42,
    steps: int = 15,
    dt: float = 0.7,
    verbose: bool = True,
    mode_store: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    quiet_compare = COMPARE_EXPANSIONS and PRINT_ONLY_FINAL_COMPARE
    do_print = verbose and not quiet_compare

    cfg = {"seed": seed, "steps": steps, "dt": dt}
    if mode_store is None:
        store: Dict[str, Any] = {}
        own_mode_store = True
    else:
        store = mode_store
        own_mode_store = False

    ret = run_one_mode_into_store(
        store,
        expansion_mode,
        cfg,
        do_print=do_print,
    )
    if own_mode_store:
        export_mode_csvs(store)
        make_plots(store)
        print_holo_output_files()
    return ret


def main() -> None:
    # Pseudocode-Struktur:
    #   mode_store = {}
    #   for mode_name, mode_cfg in expansion_modes.items():
    #       init -> für step: Dynamik, Klassen, Holo, Ausgabe, append_* ...
    #   export_mode_csvs(mode_store); make_plots(mode_store)
    mode_store: Dict[str, Any] = {}
    quiet_compare = COMPARE_EXPANSIONS and PRINT_ONLY_FINAL_COMPARE
    do_print = not quiet_compare
    expansion_modes = expansion_modes_for_current_run()

    if COMPARE_EXPANSIONS:
        print("\n" + "=" * 120)
        print("Vergleich E3 / E4")
        print("=" * 120)

    results: List[Dict[str, Any]] = []
    first = True
    for mode_name, mode_cfg in expansion_modes.items():
        if COMPARE_EXPANSIONS and not first:
            print("\n" + "=" * 120 + "\n")
        first = False
        results.append(
            run_one_mode_into_store(
                mode_store,
                mode_name,
                mode_cfg,
                do_print=do_print,
            )
        )

    if COMPARE_EXPANSIONS and len(results) >= 2:
        print("\nVergleichstabelle")
        print("-" * 120)
        print(
            f'{"Modus":>8} | {"Kinder":>6} | {"Radiusfaktor":>12} | '
            f'{"H(final)":>10} | {"Bulk(final)":>12} | {"S3(final)":>10}'
        )
        print("-" * 120)
        for r in results:
            print(
                f'{r["mode"]:>8} | {r["child_count"]:6d} | {r["radius_factor"]:12.6f} | '
                f'{r["final_H"]:10.6f} | {r["final_bulk"]:12.6f} | {r["final_s3"]:10.6f}'
            )

    export_mode_csvs(mode_store)
    make_plots(mode_store)
    print_holo_output_files()


if __name__ == "__main__":
    main()