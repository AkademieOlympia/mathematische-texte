# ============================================================
# Kwant / EABC / 60->420 Achsenmoden-Prototyp
# 42 = 6 Kanäle mod 60 x 7 Phasen mod 7
#
# Ziel:
# H = I_6 ⊗ H_7 + H_axis
#
# Basis:
#   |r, phi>
#   r   in [11,17,29,41,47,59]
#   phi in [0,...,6]
#
# Achsenregel:
#   dist=30 -> d25-Modus
#   dist=18 -> d34-Modus
#   dist=6  -> d16-Modus
#
# Empirische Gewichte beispielhaft aus N=5e10:
#   purity_30 ~ 0.143913
#   purity_18 ~ 0.138114
#   purity_6  ~ 0.131379
#
# Driftmittel:
#   dist=30: d16=-0.019765, d25=+0.023100, d34=-0.001138
#   dist=18: d16=+0.013970, d25=-0.004650, d34=-0.017144
#   dist=6 : d16=-0.019787, d25=+0.005264, d34=-0.006185
# ============================================================

import numpy as np
import pandas as pd
from pathlib import Path


# ------------------------------------------------------------
# 1. Basis
# ------------------------------------------------------------

channels = [11, 17, 29, 41, 47, 59]
phis = list(range(7))

n_r = len(channels)
n_phi = len(phis)
dim = n_r * n_phi


def idx(r, phi):
    """Index der Basis |r,phi>."""
    return channels.index(r) * n_phi + phi


def cyclic_dist(r1, r2, mod=60):
    """Zyklischer Abstand modulo 60."""
    d = abs(r2 - r1) % mod
    return min(d, mod - d)


# ------------------------------------------------------------
# 2. Unperturbierter 7-Phasen-Hamiltonian
# ------------------------------------------------------------

def build_H7_chain(t=1.0):
    """
    Offene 7er-Kette:
    Spektrum: 2t cos(k*pi/8), k=1,...,7.
    """
    H7 = np.zeros((7, 7), dtype=float)
    for p in range(6):
        H7[p, p + 1] = t
        H7[p + 1, p] = t
    return H7


def build_base_H(t_phase=1.0):
    """
    H0 = I_6 ⊗ H7.
    Das erzeugt dein bisheriges Spektrum:
    sieben Eigenwerte, jeweils sechsfach entartet.
    """
    H7 = build_H7_chain(t_phase)
    H0 = np.kron(np.eye(n_r), H7)
    return H0


# ------------------------------------------------------------
# 3. Empirische Achsendaten
# ------------------------------------------------------------

# Purity-Werte aus dem N=5e10-Achsenscan.
# Diese werden als Grundstärke der Achsenkopplung verwendet.
purity_by_dist = {
    30: 0.143913,
    18: 0.138114,
    6: 0.131379,
}

# Driftvektoren nach Distanzklasse.
# Reihenfolge: d16, d25, d34.
drifts_by_dist = {
    30: {"d16": -0.019765, "d25": +0.023100, "d34": -0.001138},
    18: {"d16": +0.013970, "d25": -0.004650, "d34": -0.017144},
    6: {"d16": -0.019787, "d25": +0.005264, "d34": -0.006185},
}

dominant_by_dist = {
    30: "d25",
    18: "d34",
    6: "d16",
}


def phase_pairs_for_mode(mode):
    """
    Gegenphasenpaare.
    d16: phi 1 <-> 6
    d25: phi 2 <-> 5
    d34: phi 3 <-> 4
    """
    if mode == "d16":
        return [(1, 6), (6, 1)]
    if mode == "d25":
        return [(2, 5), (5, 2)]
    if mode == "d34":
        return [(3, 4), (4, 3)]
    raise ValueError(f"Unbekannter Modus: {mode}")


# ------------------------------------------------------------
# 4. Achsenkopplung H_axis
# ------------------------------------------------------------

def build_axis_H(
    alpha=1.0,
    beta=1.0,
    use_signed_drift=True,
    connect_all_phases_weakly=True,
):
    """
    Baut H_axis zwischen unterschiedlichen 60-Kanälen.

    alpha:
        globale Stärke der Achsenkopplung.

    beta:
        Stärke der dominanten Gegenphasen-Selektion.

    use_signed_drift:
        Wenn True, wird das Vorzeichen des empirischen dominanten Drifts
        in die Kopplung eingebaut.
        Wenn False, wird nur der Betrag verwendet.

    connect_all_phases_weakly:
        Wenn True, bekommen alle gleichen phi->phi-Hoppings eine kleine
        isotrope Grundkopplung. Zusätzlich wird der dominante Modus verstärkt.
    """

    H = np.zeros((dim, dim), dtype=float)

    for i, r1 in enumerate(channels):
        for j, r2 in enumerate(channels):
            if j <= i:
                continue

            dist = cyclic_dist(r1, r2)

            # Wir betrachten hier nur die offenen C/A-relevanten Distanzen.
            if dist not in (6, 18, 30):
                continue

            purity = purity_by_dist[dist]
            mode = dominant_by_dist[dist]
            drift = drifts_by_dist[dist][mode]

            if use_signed_drift:
                mode_weight = beta * drift
            else:
                mode_weight = beta * abs(drift)

            # Grundkopplung aus Purity.
            # Skaliere klein, damit H_axis das H7-Spektrum nur aufspaltet,
            # aber nicht völlig dominiert.
            base_weight = alpha * purity

            # 4a. Schwache isotrope phi->phi-Kopplung zwischen Kanälen
            if connect_all_phases_weakly:
                for phi in phis:
                    a = idx(r1, phi)
                    b = idx(r2, phi)
                    H[a, b] += base_weight
                    H[b, a] += base_weight

            # 4b. Dominanter Gegenphasenmodus
            for phi_a, phi_b in phase_pairs_for_mode(mode):
                a = idx(r1, phi_a)
                b = idx(r2, phi_b)

                # base + drift-Selektion
                w = base_weight + mode_weight

                H[a, b] += w
                H[b, a] += w

    return H


# ------------------------------------------------------------
# 5. Diagnostik
# ------------------------------------------------------------

def spectrum_report(H, decimals=8, evals=None):
    if evals is None:
        evals = np.linalg.eigvalsh(H)
    rounded = np.round(evals, decimals)
    unique, counts = np.unique(rounded, return_counts=True)

    print("\n--- Spectrum report ---")
    for u, c in zip(unique, counts):
        print(f"{u: .8f}  multiplicity={c}")

    return evals


def degeneracy_breaking_score(evals, decimals=8):
    """
    Misst grob, wie stark die ursprüngliche sechsfache Entartung
    aufgebrochen wurde.
    """
    rounded = np.round(evals, decimals)
    unique, counts = np.unique(rounded, return_counts=True)

    max_mult = int(np.max(counts))
    n_unique = len(unique)

    return {
        "n_unique": n_unique,
        "max_multiplicity": max_mult,
        "min_e": float(np.min(evals)),
        "max_e": float(np.max(evals)),
        "bandwidth": float(np.max(evals) - np.min(evals)),
    }


def band_splitting_report(evals0, evals):
    """
    Ordnet jeden gestörten Eigenwert dem nächsten ungestörten 7er-Level zu
    und berichtet die Aufspaltung pro Band.
    """
    base_levels = np.unique(np.round(evals0, 12))

    print("\n--- Band splitting report ---")
    for E0 in base_levels:
        # sechs nächste gestörte Eigenwerte zu diesem alten Level
        distances = np.abs(evals - E0)
        inds = np.argsort(distances)[:6]
        band = np.sort(evals[inds])

        print(f"\nBase E0={E0: .8f}")
        print(f"  perturbed: {', '.join(f'{x:.8f}' for x in band)}")
        print(f"  center={np.mean(band): .8f}")
        print(f"  width ={np.max(band)-np.min(band): .8f}")


def axis_weights_of_state(vec):
    """
    Grobe Diagnose:
    Wie viel Gewicht eines Eigenvektors liegt auf Kanalpaaren
    der Distanzen 6,18,30?

    Hier messen wir einfach Kanalgewichte, nicht Kantenströme.
    Für feinere Analyse müsste man bond currents definieren.
    """
    weights_by_channel = {}

    for r in channels:
        s = 0.0
        for phi in phis:
            a = idx(r, phi)
            s += abs(vec[a]) ** 2
        weights_by_channel[r] = s

    # Paargewichte nach Distanzklasse:
    out = {6: 0.0, 18: 0.0, 30: 0.0}

    for i, r1 in enumerate(channels):
        for j, r2 in enumerate(channels):
            if j <= i:
                continue
            dist = cyclic_dist(r1, r2)
            if dist in out:
                out[dist] += weights_by_channel[r1] * weights_by_channel[r2]

    return out


def bond_axis_weights(H, vec):
    """
    Direkte Bond-/Hopping-Gewichte eines Eigenzustands
    nach Achsendistanz 6, 18, 30.

    Misst, welche Hamiltonian-Kanten vom Zustand effektiv getragen werden.
    """
    out = {6: 0.0, 18: 0.0, 30: 0.0}

    for ia, r1 in enumerate(channels):
        for ib, r2 in enumerate(channels):
            if ib <= ia:
                continue

            dist = cyclic_dist(r1, r2)
            if dist not in out:
                continue

            for phi1 in phis:
                for phi2 in phis:
                    a = idx(r1, phi1)
                    b = idx(r2, phi2)

                    hij = H[a, b]
                    if abs(hij) < 1e-14:
                        continue

                    out[dist] += abs(hij) * abs(vec[a]) * abs(vec[b])

    s = sum(out.values())
    if s > 0:
        out = {k: v / s for k, v in out.items()}

    return out


def eigenstate_axis_report(H, n_show=10):
    evals, evecs = np.linalg.eigh(H)

    print("\n--- Eigenstate axis weights (Kanal-Paare), lowest states ---")
    for k in range(min(n_show, len(evals))):
        v = evecs[:, k]
        w = axis_weights_of_state(v)
        print(
            f"k={k:2d} E={evals[k]: .8f} "
            f"w6={w[6]:.4f} w18={w[18]:.4f} w30={w[30]:.4f}"
        )

    print("\n--- Eigenstate axis weights (Kanal-Paare), highest states ---")
    for k in range(len(evals) - n_show, len(evals)):
        v = evecs[:, k]
        w = axis_weights_of_state(v)
        print(
            f"k={k:2d} E={evals[k]: .8f} "
            f"w6={w[6]:.4f} w18={w[18]:.4f} w30={w[30]:.4f}"
        )


def eigenstate_bond_axis_report(H, n_show=10):
    evals, evecs = np.linalg.eigh(H)

    print("\n--- Eigenstate BOND axis weights, lowest states ---")
    for k in range(min(n_show, len(evals))):
        w = bond_axis_weights(H, evecs[:, k])
        print(
            f"k={k:2d} E={evals[k]: .8f} "
            f"B6={w[6]:.4f} B18={w[18]:.4f} B30={w[30]:.4f}"
        )

    print("\n--- Eigenstate BOND axis weights, highest states ---")
    for k in range(len(evals) - n_show, len(evals)):
        w = bond_axis_weights(H, evecs[:, k])
        print(
            f"k={k:2d} E={evals[k]: .8f} "
            f"B6={w[6]:.4f} B18={w[18]:.4f} B30={w[30]:.4f}"
        )


def compute_global_bond_data(H):
    """Ein eigh(H); Mittel über Zustände + dominanter-Bond-Zählungen."""
    evals, evecs = np.linalg.eigh(H)

    rows = []
    for k in range(len(evals)):
        w = bond_axis_weights(H, evecs[:, k])
        rows.append((evals[k], w[6], w[18], w[30]))

    arr = np.array(rows)

    mean6 = float(np.mean(arr[:, 1]))
    mean18 = float(np.mean(arr[:, 2]))
    mean30 = float(np.mean(arr[:, 3]))

    dom_counts = {6: 0, 18: 0, 30: 0}
    for _, b6, b18, b30 in rows:
        vals = {6: b6, 18: b18, 30: b30}
        dom = max(vals, key=vals.get)
        dom_counts[dom] += 1

    return arr, mean6, mean18, mean30, dom_counts


def global_bond_axis_summary_values(H):
    evals, evecs = np.linalg.eigh(H)

    rows = []
    for k in range(len(evals)):
        w = bond_axis_weights(H, evecs[:, k])
        rows.append((evals[k], w[6], w[18], w[30]))

    arr = np.array(rows)

    dom_counts = {6: 0, 18: 0, 30: 0}
    for _, b6, b18, b30 in rows:
        vals = {6: b6, 18: b18, 30: b30}
        dom = max(vals, key=vals.get)
        dom_counts[dom] += 1

    return {
        "mean_B6": float(np.mean(arr[:, 1])),
        "mean_B18": float(np.mean(arr[:, 2])),
        "mean_B30": float(np.mean(arr[:, 3])),
        "dom_counts": dom_counts,
    }


def print_bond_axis_summary_text(bond_vals):
    """Nur Ausgabe, keine neue Diagonalisation."""
    dc = bond_vals["dom_counts"]
    print("\n--- Global BOND axis summary over all eigenstates ---")
    print(f"mean B6  = {bond_vals['mean_B6']:.6f}")
    print(f"mean B18 = {bond_vals['mean_B18']:.6f}")
    print(f"mean B30 = {bond_vals['mean_B30']:.6f}")
    print("dominant bond counts:", dc)


def assemble_variant_h(H0, alpha, beta, use_signed_drift):
    """H = H0 + H_axis für eine Parameterwahl."""
    H_axis = build_axis_H(
        alpha=alpha,
        beta=beta,
        use_signed_drift=use_signed_drift,
        connect_all_phases_weakly=True,
    )
    return H0 + H_axis


def global_bond_axis_summary(H):
    arr, mean6, mean18, mean30, dom_counts = compute_global_bond_data(H)

    print("\n--- Global BOND axis summary over all eigenstates ---")
    print(f"mean B6  = {mean6:.6f}")
    print(f"mean B18 = {mean18:.6f}")
    print(f"mean B30 = {mean30:.6f}")
    print("dominant bond counts:", dom_counts)

    return arr, {
        "mean_B6": mean6,
        "mean_B18": mean18,
        "mean_B30": mean30,
        "dom6": dom_counts[6],
        "dom18": dom_counts[18],
        "dom30": dom_counts[30],
        "dom_counts": dom_counts,
    }


# ------------------------------------------------------------
# 6. Hauptlauf
# ------------------------------------------------------------

VARIANTS = [
    ("purity_only", 0.35, 0.0, False),
    ("abs_drift", 0.35, 4.0, False),
    ("signed_drift", 0.35, 4.0, True),
]


if __name__ == "__main__":
    print("EABC/Kwant-Prototyp: 42 = 6 x 7")
    print("Kanäle:", channels)
    print("Phasen:", phis)

    # 6.1 Basismodell
    H0 = build_base_H(t_phase=1.0)
    evals0 = spectrum_report(H0)

    print("\nBase degeneracy:")
    print(degeneracy_breaking_score(evals0))

    # 6.2 Achsenmodell — Varianten (Name, alpha, beta, use_signed_drift)
    # alpha klein halten; sonst dominiert H_axis die 7-Kette zu stark.
    bond_rows_by_variant = {}

    def run_variant(name, alpha, beta, signed):
        """CSV-Zeile zu H0 (Basishamiltonian aus build_base_H)."""
        H = assemble_variant_h(H0, alpha, beta, signed)
        evals = np.linalg.eigvalsh(H)
        deg = degeneracy_breaking_score(evals)

        bond_vals = global_bond_axis_summary_values(H)

        print(
            f"{name},"
            f"{deg['n_unique']},"
            f"{deg['max_multiplicity']},"
            f"{deg['bandwidth']:.8f},"
            f"{bond_vals['mean_B6']:.6f},"
            f"{bond_vals['mean_B18']:.6f},"
            f"{bond_vals['mean_B30']:.6f},"
            f"{bond_vals['dom_counts'][6]},"
            f"{bond_vals['dom_counts'][18]},"
            f"{bond_vals['dom_counts'][30]}",
        )

        return H, evals, deg, bond_vals

    print(
        "\nvariant,n_unique,max_mult,bandwidth,mean_B6,mean_B18,mean_B30,dom6,dom18,dom30"
    )

    for name, alpha, beta, signed in [
        ("purity_only", 0.35, 0.0, False),
        ("abs_drift", 0.35, 4.0, False),
        ("signed_drift", 0.35, 4.0, True),
    ]:
        run_variant(name, alpha, beta, signed)

    for name, alpha, beta, signed in VARIANTS:
        sep = "=" * 60
        print(
            f"\n{sep}\nVARIANT: {name}  alpha={alpha} beta={beta} signed={signed}\n{sep}"
        )

        H = assemble_variant_h(H0, alpha, beta, signed)
        evals = np.linalg.eigvalsh(H)
        deg = degeneracy_breaking_score(evals)
        arr, mean6, mean18, mean30, dom_counts = compute_global_bond_data(H)
        bond_vals = {
            "mean_B6": mean6,
            "mean_B18": mean18,
            "mean_B30": mean30,
            "dom_counts": dom_counts,
        }
        bond_rows_by_variant[name] = arr

        spectrum_report(H, evals=evals)

        print("\nPerturbed degeneracy:")
        print(deg)

        band_splitting_report(evals0, evals)

        eigenstate_axis_report(H, n_show=8)
        eigenstate_bond_axis_report(H, n_show=8)

        print_bond_axis_summary_text(bond_vals)

    # 6.3 Beta-Sweep bei festem alpha=0.35, signed drift
    print("\nbeta,bandwidth,mean_B6,mean_B18,mean_B30,dom6,dom18,dom30")

    for beta in [0.0, 1.0, 2.0, 4.0, 8.0]:
        H_axis = build_axis_H(
            alpha=0.35,
            beta=beta,
            use_signed_drift=True,
            connect_all_phases_weakly=True,
        )
        H = H0 + H_axis
        evals = np.linalg.eigvalsh(H)
        deg = degeneracy_breaking_score(evals)
        bond = global_bond_axis_summary_values(H)

        print(
            f"{beta:.1f},"
            f"{deg['bandwidth']:.8f},"
            f"{bond['mean_B6']:.6f},"
            f"{bond['mean_B18']:.6f},"
            f"{bond['mean_B30']:.6f},"
            f"{bond['dom_counts'][6]},"
            f"{bond['dom_counts'][18]},"
            f"{bond['dom_counts'][30]}",
        )

    # 6.4 Raster alpha × beta (signed drift) -> kwant_alpha_beta.csv
    _root = Path(__file__).resolve().parent
    alphas = [0.1, 0.2, 0.35, 0.5, 1.0]
    betas = [0.0, 1.0, 2.0, 4.0, 8.0]

    rows = []

    for alpha in alphas:
        for beta in betas:
            H_axis = build_axis_H(
                alpha=alpha,
                beta=beta,
                use_signed_drift=True,
                connect_all_phases_weakly=True,
            )
            H = H0 + H_axis
            evals = np.linalg.eigvalsh(H)
            deg = degeneracy_breaking_score(evals)
            bond = global_bond_axis_summary_values(H)

            rows.append(
                {
                    "alpha": alpha,
                    "beta": beta,
                    "bandwidth": deg["bandwidth"],
                    "mean_B6": bond["mean_B6"],
                    "mean_B18": bond["mean_B18"],
                    "mean_B30": bond["mean_B30"],
                    "dom6": bond["dom_counts"][6],
                    "dom18": bond["dom_counts"][18],
                    "dom30": bond["dom_counts"][30],
                }
            )

    _ab_path = _root / "kwant_alpha_beta.csv"
    pd.DataFrame(rows).to_csv(_ab_path, index=False)
    print(f"\nWrote {_ab_path} ({len(rows)} rows)")
