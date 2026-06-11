#!/usr/bin/env python3
"""Validiert und visualisiert die 2-In-2-Out-Eisregel am EABC-Tetraeder."""

from __future__ import annotations

state_plus = {"katzen": -1.0, "shear_norm": 2.0}   # "In"
state_minus = {"katzen": 1.0, "shear_norm": 2.0}  # "Out"

tetrahedron_states = {
    0: state_plus,
    1: state_plus,
    2: state_minus,
    3: state_minus,
}


def label(state: dict) -> str:
    return "In" if state["katzen"] < 0 else "Out"


def validate_ice_rule(states: dict[int, dict]) -> bool:
    """Eisregel: genau zwei Ecken 'In' und zwei 'Out'."""
    n_in = sum(1 for s in states.values() if s["katzen"] < 0)
    n_out = sum(1 for s in states.values() if s["katzen"] > 0)
    return n_in == 2 and n_out == 2 and len(states) == 4


def bond_energy_local(sigma_a: float, sigma_b: float, j_local: float = 1.0) -> float:
    return j_local * sigma_a * sigma_b


def main() -> None:
    print("=== EABC-Tetraeder: Eisregel-Validierung ===\n")
    for idx, state in tetrahedron_states.items():
        print(f"  Ecke {idx}: {label(state):>3}  (katzen={state['katzen']:+.1f}, "
              f"shear_norm={state['shear_norm']:.1f})")

    ok = validate_ice_rule(tetrahedron_states)
    print(f"\nEisregel (2-In-2-Out) erfüllt: {ok}")

    print("\nLokale Bindungsenergien E_loc = J * sigma_a * sigma_b (J=1):")
    total = 0.0
    for a in range(4):
        for b in range(a + 1, 4):
            sa = tetrahedron_states[a]["katzen"]
            sb = tetrahedron_states[b]["katzen"]
            e = bond_energy_local(sa, sb)
            total += e
            kind = f"{label(tetrahedron_states[a])}--{label(tetrahedron_states[b])}"
            print(f"  {a}<->{b} ({kind:>7}): {e:+.1f}")
    print(f"  Summe E_loc: {total:+.1f}")

    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        # Reguläres Tetraeder
        verts = [
            (1, 1, 1),
            (1, -1, -1),
            (-1, 1, -1),
            (-1, -1, 1),
        ]
        xs, ys, zs = zip(*verts)
        center = (0.0, 0.0, 0.0)

        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_trisurf(xs, ys, zs, color="lightgray", alpha=0.25, edgecolor="gray")

        colors = {"In": "#2166ac", "Out": "#d6604d"}
        for i, (x, y, z) in enumerate(verts):
            lab = label(tetrahedron_states[i])
            ax.scatter([x], [y], [z], c=colors[lab], s=120, depthshade=True)
            ax.text(x, y, z + 0.15, f"{i}: {lab}", fontsize=9, ha="center")
            dx, dy, dz = center[0] - x, center[1] - y, center[2] - z
            if lab == "Out":
                dx, dy, dz = -dx, -dy, -dz
            ax.quiver(x, y, z, dx * 0.35, dy * 0.35, dz * 0.35,
                      color=colors[lab], arrow_length_ratio=0.25, linewidth=1.5)

        ax.set_title("EABC-Tetraeder: 2-In-2-Out (Eisregel)")
        ax.set_axis_off()
        out = "tetrahedron_ice_rule.png"
        plt.tight_layout()
        plt.savefig(out, dpi=150)
        print(f"\nVisualisierung gespeichert: {out}")
    except ImportError:
        print("\n(Hinweis: matplotlib nicht installiert — nur Textausgabe.)")


if __name__ == "__main__":
    main()
