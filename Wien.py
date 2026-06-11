"""
Wien.py — EABC-Resonanz mit optionalen Fenstern.

Ohne Argumente: Animationsfenster öffnet sich; nach Schließen folgt das Peak-Phase-Diagramm
(beides blockiert jeweils mit plt.show()).

  python Wien.py --no-gui           # Keine Fenster: PNGs + npz, Prozess endet sofort
  python Wien.py --skip-animation   # Nur Peak-Phase-Fenster (keine Animation)
  python Wien.py --skip-animation --no-gui   # Nur Phase-PNG + npz

Ctrl+C während ein Fenster offen ist: Abbruch in matplotlibs Hauptschleife —
Tracebacks können dann wie „Zeile beim Zoom/transform“ aussehen, obwohl nur plt.show()
wartete. Besser: Fenster schließen oder gleich --no-gui verwenden.
"""

import sys

CLI_FLAGS = {a for a in sys.argv[1:] if a.startswith("--")}
NO_GUI = "--no-gui" in CLI_FLAGS
SKIP_ANIMATION = "--skip-animation" in CLI_FLAGS

if NO_GUI:
    import matplotlib

    matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

SNAPSHOT_PATH = "wien_resonance_last.png"
PHASE_PLOT_PATH = "wien_peak_phase.png"


def plt_show_blocking():
    """Blockiert bis die Figur geschlossen ist (macOS/Interactive wie gewohnt)."""
    try:
        plt.show()
    except KeyboardInterrupt:
        plt.close("all")
        print(
            "\n[Wien] Abgebrochen, während matplotlib auf geschlossene Fenster wartete "
            "(typisch: Ctrl+C bei plt.show()).\n"
            "Ohne blockierende Fenster:\n"
            "  python Wien.py --no-gui\n",
            file=sys.stderr,
        )
        raise SystemExit(130) from None

# ------------------------------------------------------------
# 1. Datenbasis laden
# ------------------------------------------------------------

zeros = np.load("zeros6.npy")

# Option A: Nullstellen selbst
gamma_raw = zeros[:200]

# Option B: normierte Abstände der Nullstellen
spacings = np.diff(zeros[:201])
gamma_spacings = spacings / np.mean(spacings)

# Umschalten:
USE_SPACINGS = False

if USE_SPACINGS:
    gamma_values = gamma_spacings
    title_mod = "Riemann-Abstände"
else:
    gamma_values = gamma_raw
    title_mod = "Riemann-Nullstellen"


# ------------------------------------------------------------
# 2. EABC / mod-210-Anker
# ------------------------------------------------------------

anchor_classes = [5, 11, 101, 191]

def family210(p):
    r = p % 210
    if r == 5:
        return "Anker 5"
    if r == 11:
        return "Familie 11"
    if r == 101:
        return "Familie 101"
    if r == 191:
        return "Familie 191"
    return "Andere"


# ------------------------------------------------------------
# 3. Resonanzmodell
# ------------------------------------------------------------

# Sichtfenster auf der x-Achse; normierte Peak-Phase := peak_x / X_DOMAIN_MAX
# (bei X_DOMAIN_MAX = 40 entspricht 1/8 genau x = 5, dem ersten mod-210-Anker.)
X_DOMAIN_MAX = 40.0

# Hilfslinien (Restklassen-Skalierung wie im Originalskript)
ANCHOR_LINES_X = {5: 5 / 5, 11: 11 / 5, 101: 101 / 5, 191: 191 / 5}


def add_anchor_vlines(ax):
    for r, xpos in ANCHOR_LINES_X.items():
        if xpos <= X_DOMAIN_MAX:
            ax.axvline(xpos, linestyle="--", alpha=0.25)
            ax.text(xpos, 1.05, str(r), ha="center", va="bottom", fontsize=9)


def resonance_curve(x, T_factor, gamma_values, damping=10.0):
    """
    Riemann-modulierte EABC-Resonanzkurve.

    T_factor skaliert die Frequenzen.
    Größeres T verschiebt die Kopplungsmaxima.
    """
    y = np.zeros_like(x)

    n_modes = min(len(gamma_values), int(T_factor / 2) + 5)

    for g in gamma_values[:n_modes]:
        y += np.cos(g * x / T_factor) * np.exp(-x / damping)

    y2 = y ** 2

    if np.max(y2) > 0:
        y2 = y2 / np.max(y2)

    return y2


# Referenz 1/8 auf der Einheitsstrecke [0,1] (x = 0,125)
X_REF_EIGHTH = 0.125
# Lesart \(x_\ast \Delta_{\mathrm{PV}}\): hier \(\Delta_{\mathrm{PV}}:=8\) zur Lesart „Achtelresonanz“.
DELTA_PV_FACTOR = 8.0


def german_decimal_tex(value, decimals):
    """Dezimalzahl als Mathtext mit Dezimalkomma ({,})."""
    body = f"{value:.{decimals}f}"
    if "." not in body:
        return body
    whole, frac = body.split(".", 1)
    return whole + r"{,}" + frac


def draw_zoom_panel_zero_one(
    axins,
    T_factor,
    peak_x_val,
    peak_y_val,
    gamma_vals,
    *,
    x_domain_max,
):
    r"""Inhalt der Zoom-Achse 0 ≤ x ≤ 1 (axins wird von außen ggf. mit clear() geleert)."""
    xz = np.linspace(0.0, 1.0, 800)
    yz = resonance_curve(xz, T_factor, gamma_vals)
    axins.plot(xz, yz, lw=1.65, color="C0")
    axins.scatter([peak_x_val], [peak_y_val], color="red", s=32, zorder=6)

    axins.axvline(X_REF_EIGHTH, color="C1", linestyle="--", lw=1.35, alpha=0.95)
    axins.axvline(peak_x_val, color="C2", linestyle="-", lw=1.25, alpha=0.9)

    ymax = float(np.max(yz)) if yz.size else 1.0
    axins.set_xlim(0, 1)
    axins.set_ylim(0, ymax * 1.14 + 1e-12)
    axins.set_xlabel("$x$", fontsize=9)
    axins.set_title(r"Zoom: $0 \leq x \leq 1$", fontsize=9, pad=3)
    axins.tick_params(labelsize=8)
    axins.grid(True, alpha=0.28)

    phi_star = peak_x_val / x_domain_max
    x_delta_pv = peak_x_val * DELTA_PV_FACTOR

    tb_text = (
        r"$x_\ast \Delta_{\mathrm{PV}} = "
        + german_decimal_tex(x_delta_pv, 2)
        + r", \qquad \phi_\ast = "
        + german_decimal_tex(phi_star, 4)
        + r"$"
    )
    axins.text(
        0.03,
        0.97,
        tb_text,
        transform=axins.transAxes,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="wheat", edgecolor="0.35", alpha=0.94),
    )

    lbl_dy_ref = 0.08
    lbl_dy_peak = 0.02
    axins.text(
        X_REF_EIGHTH + 0.012,
        lbl_dy_ref,
        r"$1/8$",
        transform=axins.get_xaxis_transform(),
        fontsize=8,
        color="C1",
    )
    peak_lbl_x = min(peak_x_val + 0.012, 0.93)
    axins.text(
        peak_lbl_x,
        lbl_dy_peak,
        r"$x_\ast$",
        transform=axins.get_xaxis_transform(),
        fontsize=8,
        color="C2",
    )


def add_zoom_panel_zero_one(ax, T_factor, peak_x_val, peak_y_val, gamma_vals, *, x_domain_max):
    axins = inset_axes(ax, width="44%", height="44%", loc="upper right", borderpad=0.65)
    draw_zoom_panel_zero_one(
        axins,
        T_factor,
        peak_x_val,
        peak_y_val,
        gamma_vals,
        x_domain_max=x_domain_max,
    )
    return axins


def eabc_resonance_dynamic(T_min=1, T_max=150, frames_count=200):
    x = np.linspace(0.1, X_DOMAIN_MAX, 2000)

    T_values = np.linspace(T_min, T_max, frames_count)

    peak_x = np.empty(frames_count)
    peak_y = np.empty(frames_count)
    for i in range(frames_count):
        y = resonance_curve(x, T_values[i], gamma_values)
        k = int(np.argmax(y))
        peak_x[i] = x[k]
        peak_y[i] = y[k]

    if SKIP_ANIMATION:
        pass
    elif NO_GUI:
        fig, ax = plt.subplots(figsize=(11, 6))
        add_anchor_vlines(ax)
        ax.set_xlim(0, X_DOMAIN_MAX)
        ax.set_ylim(0, 1.15)
        ax.set_xlabel("Strukturabstand / skalierte Wellenlänge")
        ax.set_ylabel("normierte geometrische Starrheit")
        ax.grid(True, alpha=0.3)
        i_last = frames_count - 1
        T_last = T_values[i_last]
        y_last = resonance_curve(x, T_last, gamma_values)
        ax.plot(x, y_last, lw=2, label="Kopplungsintensität")
        ax.plot(peak_x[i_last], peak_y[i_last], "ro", label="Maximum")
        phase_frac = peak_x[i_last] / X_DOMAIN_MAX
        ax.set_title(
            f"EABC-Resonanz ({title_mod}) | T = {T_last:.2f} (letzter Rahmen), "
            f"Peak x = {peak_x[i_last]:.3f}, Phase = {phase_frac:.4f}"
        )
        ax.legend()
        add_zoom_panel_zero_one(
            ax,
            T_last,
            peak_x[i_last],
            peak_y[i_last],
            gamma_values,
            x_domain_max=X_DOMAIN_MAX,
        )
        fig.savefig(SNAPSHOT_PATH, dpi=150)
        plt.close(fig)
        print(f"Letzter Animationsrahmen gespeichert: {SNAPSHOT_PATH}")
    else:
        fig, ax = plt.subplots(figsize=(11, 6))

        line, = ax.plot([], [], lw=2, label="Kopplungsintensität")
        peak_dot, = ax.plot([], [], "ro", label="aktuelles Maximum")

        add_anchor_vlines(ax)

        ax.set_xlim(0, X_DOMAIN_MAX)
        ax.set_ylim(0, 1.15)
        ax.set_title(f"Dynamische Kopplung im EABC-Modell ({title_mod})")
        ax.set_xlabel("Strukturabstand / skalierte Wellenlänge")
        ax.set_ylabel("normierte geometrische Starrheit")
        ax.grid(True, alpha=0.3)
        ax.legend()

        ax_zoom = inset_axes(ax, width="44%", height="44%", loc="upper right", borderpad=0.65)

        def animate(i):
            T_factor = T_values[i]
            y = resonance_curve(x, T_factor, gamma_values)

            line.set_data(x, y)
            px = peak_x[i]
            py = peak_y[i]
            peak_dot.set_data([px], [py])

            phase_frac = px / X_DOMAIN_MAX
            ax.set_title(
                f"EABC-Resonanz ({title_mod}) | T = {T_factor:.2f}, "
                f"Peak x = {px:.3f}, Phase = {phase_frac:.4f}"
            )

            ax_zoom.clear()
            draw_zoom_panel_zero_one(
                ax_zoom,
                T_factor,
                px,
                py,
                gamma_values,
                x_domain_max=X_DOMAIN_MAX,
            )

            return line, peak_dot

        _animation_handle = FuncAnimation(
            fig,
            animate,
            frames=frames_count,
            blit=False,
            interval=80,
        )

        plt_show_blocking()
        plt.close(fig)
        del _animation_handle

    return np.array(T_values), peak_x, peak_y


def plot_peak_phase_vs_T(
    T_values,
    peak_x,
    title_suffix="",
    x_domain_max=X_DOMAIN_MAX,
    *,
    show=True,
    save_path=None,
):
    """Normierte Peak-Phase auf [0,1]; Vergleich mit 1/8 (= erster Anker bei x=5 wenn x_max=40)."""
    peak_phase = peak_x / x_domain_max
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(T_values, peak_phase, lw=1.5, label="Peak-Phase = peak_x / x_max")
    ax.axhline(
        1 / 8,
        color="C1",
        linestyle="--",
        alpha=0.9,
        label=f"1/8  ↔  x = {x_domain_max / 8:.1f} (Anker 5 bei x_max={x_domain_max:g})",
    )
    ax.set_xlabel("T")
    ax.set_ylabel("Peak-Phase")
    ax.set_title(f"Peak-Phase über T ({title_suffix})".strip())
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"Peak-Phase-Diagramm gespeichert: {save_path}")
    if show:
        plt_show_blocking()
    plt.close(fig)
    return peak_phase


# ------------------------------------------------------------
# 4. Lauf
# ------------------------------------------------------------

T_values, peak_x, peak_y = eabc_resonance_dynamic(
    T_min=1,
    T_max=150,
    frames_count=200
)

peak_phase = peak_x / X_DOMAIN_MAX

# Peak-Daten speichern
np.savez(
    "eabc_riemann_resonance_peaks.npz",
    T_values=T_values,
    peak_x=peak_x,
    peak_y=peak_y,
    peak_phase=peak_phase,
    x_domain_max=X_DOMAIN_MAX,
    use_spacings=USE_SPACINGS
)

print("Peak-Daten gespeichert: eabc_riemann_resonance_peaks.npz")
print(
    f"Peak-Phase: min={peak_phase.min():.4f}, max={peak_phase.max():.4f}, "
    f"letzter Wert={peak_phase[-1]:.4f} (Referenz 1/8 = {1/8:.4f})"
)

plot_peak_phase_vs_T(
    T_values,
    peak_x,
    title_suffix=title_mod,
    show=not NO_GUI,
    save_path=PHASE_PLOT_PATH if NO_GUI else None,
)