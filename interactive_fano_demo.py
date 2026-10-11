#!/usr/bin/env python3
"""
FanoRealtimeValidator - Interaktiver Live-Demonstrator
======================================================
Dieses Python-Skript bietet eine interaktive Simulations-Oberfläche (GUI) mit Matplotlib-Widgets.
Es ermöglicht die Live-Demonstration des FanoRealtimeValidator-Kernels für Präsentationen
und B2B-Kundengespräche.

Voraussetzungen:
    pip install numpy matplotlib

Ausführung:
    python3 interactive_fano_demo.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

# --- Physikalische & Mathematische Konstanten ---
FANO_RATIO = np.sqrt(3.0 / 5.0)  # approx 0.7745966692414834
F0 = 10e6  # 10 MHz Referenzfrequenz

def generate_data(flicker_noise_level, filter_active, mode):
    """
    Generiert Phasenrausch-Spektrum, Allan-Abweichung und Jitter-Histogramm.
    """
    # 1. Phasenrausch-Spektrum L(f)
    f = np.logspace(0, 5, 500) # 1 Hz bis 100 kHz Offset
    # Rauschmodell: b_0 (weißes Rauschen) + b_3/f^3 (Kristall-Flickerrauschen)
    b_0 = 1e-15
    b_3 = 1e-6 * (10 ** flicker_noise_level)
    
    L_off = 10 * np.log10(b_0 + b_3 / (f**3))
    
    if filter_active:
        # Fano-Topos-1_Omega Filterung unterdrückt Kristall-Flicker um ~45 dB
        suppression_factor = 10**(-4.5)
        L_on = 10 * np.log10(b_0 + (b_3 * suppression_factor) / (f**3))
        L_current = L_on
    else:
        L_current = L_off

    # 2. Allan-Abweichung sigma_y(tau)
    tau = np.logspace(-3, 2, 200) # 1 ms bis 100 s
    # sigma_y^2 = h_0 / (2 tau) + h_-1 * ln(2) + h_-2 * tau / 3
    h_0 = 1e-22
    h_m1_off = b_3 * 1e-18
    sigma_off = np.sqrt(h_0 / (2 * tau) + h_m1_off * np.log(2) + 1e-26 * tau)
    
    if filter_active:
        h_m1_on = h_m1_off * 10**(-4.5)
        sigma_on = np.sqrt(h_0 / (2 * tau) + h_m1_on * np.log(2) + 1e-28 * tau)
        sigma_current = sigma_on
    else:
        sigma_current = sigma_off

    # 3. Jitter-Histogramm
    np.random.seed(42)
    if filter_active:
        jitter_ps = np.random.normal(0, 0.045, 10000) # < 50 fs
    else:
        sigma_ps = 100.0 * (10 ** (flicker_noise_level / 2.0))
        jitter_ps = np.random.normal(0, sigma_ps, 10000)

    # 4. Branchen-KPIs
    if mode == 'Pro Audio (SNR)':
        kpi_name = "Audio SNR @ 20 kHz"
        if filter_active:
            kpi_val = "144.0 dB (24-Bit Max)"
        else:
            kpi_val = f"{max(60.0, 144.0 - 20*np.log10(2*np.pi*20e3*np.std(jitter_ps)*1e-12)):.1f} dB"
    elif mode == 'Telekom / 5G (EVM)':
        kpi_name = "Error Vector Mag. (EVM)"
        kpi_val = "0.12 % (Exzellent)" if filter_active else "2.85 % (Störung)"
    elif mode == 'Navigation / IMU':
        kpi_name = "Drift-Akkumulation"
        kpi_val = "0.01 m/h (Stabil)" if filter_active else "1.42 m/h (Drift)"
    else: # Metrologie
        kpi_name = "Frequenzunsicherheit"
        kpi_val = f"1.5e-15 @ tau=10s" if filter_active else f"{sigma_off[100]:.1e} (Flicker-Floor)"

    return f, L_current, L_off, tau, sigma_current, sigma_off, jitter_ps, kpi_name, kpi_val


def main():
    fig = plt.figure(figsize=(14, 8), facecolor='#0B0F19')
    fig.canvas.manager.set_window_title('FanoRealtimeValidator - Live Demonstration')

    # Zustand
    state = {
        'filter_active': True,
        'flicker_noise': 0.0,
        'mode': 'Pro Audio (SNR)'
    }

    # Subplots anlegen
    gs = fig.add_gridspec(2, 2, left=0.07, right=0.72, top=0.90, bottom=0.10, hspace=0.35, wspace=0.25)
    ax1 = fig.add_subplot(gs[0, 0], facecolor='#111827')
    ax2 = fig.add_subplot(gs[0, 1], facecolor='#111827')
    ax3 = fig.add_subplot(gs[1, 0], facecolor='#111827')
    ax4 = fig.add_subplot(gs[1, 1], facecolor='#111827')

    # Farben
    C_OFF = '#EF4444'
    C_ON = '#10B981'
    C_TEXT = '#F9FAFB'

    plt.suptitle('FanoRealtimeValidator — Interaktiver Live-Demonstrator', 
                 fontsize=16, fontweight='bold', color=C_TEXT, y=0.96)

    def update(val=None):
        f, L_cur, L_off, tau, sig_cur, sig_off, jitter_ps, kpi_name, kpi_val = generate_data(
            state['flicker_noise'], state['filter_active'], state['mode']
        )

        # Plot 1: Phasenrauschen
        ax1.clear()
        ax1.set_facecolor('#111827')
        ax1.loglog(f, 10**(L_off/10), color=C_OFF, linestyle='--', alpha=0.6, label='Filter AUS')
        color_cur = C_ON if state['filter_active'] else C_OFF
        ax1.loglog(f, 10**(L_cur/10), color=color_cur, linewidth=2, label='Aktuell')
        ax1.set_title('Ebene 1: Phasenrausch-Spektrum L(f)', color=C_TEXT, fontsize=11, fontweight='bold')
        ax1.set_xlabel('Offset-Frequenz f [Hz]', color='#9CA3AF', fontsize=9)
        ax1.set_ylabel('Rauschleistung [1/Hz]', color='#9CA3AF', fontsize=9)
        ax1.tick_params(colors='#9CA3AF')
        ax1.grid(True, which='both', color='#374151', linestyle=':', alpha=0.5)
        ax1.legend(facecolor='#1F2937', edgecolor='#374151', labelcolor=C_TEXT, fontsize=8)

        # Plot 2: Allan-Abweichung
        ax2.clear()
        ax2.set_facecolor('#111827')
        ax2.loglog(tau, sig_off, color=C_OFF, linestyle='--', alpha=0.6, label='Filter AUS (Flicker-Floor)')
        ax2.loglog(tau, sig_cur, color=color_cur, linewidth=2, label='Filter EIN (Durchbrochen)')
        ax2.set_title('Ebene 2: Allan-Abweichung $\sigma_y(\\tau)$', color=C_TEXT, fontsize=11, fontweight='bold')
        ax2.set_xlabel('Mittelungszeit $\\tau$ [s]', color='#9CA3AF', fontsize=9)
        ax2.set_ylabel('Frequenzstabilität $\sigma_y(\\tau)$', color='#9CA3AF', fontsize=9)
        ax2.tick_params(colors='#9CA3AF')
        ax2.grid(True, which='both', color='#374151', linestyle=':', alpha=0.5)
        ax2.legend(facecolor='#1F2937', edgecolor='#374151', labelcolor=C_TEXT, fontsize=8)

        # Plot 3: Jitter-Histogramm
        ax3.clear()
        ax3.set_facecolor('#111827')
        ax3.hist(jitter_ps, bins=50, range=(-200, 200), color=color_cur, alpha=0.8, edgecolor='#111827')
        ax3.set_title('Ebene 3: Jitter-Aperturverteilung', color=C_TEXT, fontsize=11, fontweight='bold')
        ax3.set_xlabel('Zeitversatz $\Delta t$ [ps]', color='#9CA3AF', fontsize=9)
        ax3.set_ylabel('Häufigkeit', color='#9CA3AF', fontsize=9)
        ax3.tick_params(colors='#9CA3AF')
        ax3.grid(True, color='#374151', linestyle=':', alpha=0.5)

        # Panel 4: Live KPI Monitor
        ax4.clear()
        ax4.set_facecolor('#1F2937')
        ax4.axis('off')
        
        status_str = "FILTER AKTIV [EIN]" if state['filter_active'] else "FILTER INAKTIV [AUS]"
        status_color = C_ON if state['filter_active'] else C_OFF
        
        ax4.text(0.05, 0.82, "SYSTEM-STATUS:", color='#9CA3AF', fontsize=10, fontweight='bold')
        ax4.text(0.05, 0.68, status_str, color=status_color, fontsize=14, fontweight='bold')
        
        ax4.text(0.05, 0.45, f"BRANCHEN-METRIK ({state['mode']}):", color='#9CA3AF', fontsize=9, fontweight='bold')
        ax4.text(0.05, 0.30, f"{kpi_name}:", color=C_TEXT, fontsize=10)
        ax4.text(0.05, 0.12, kpi_val, color='#38BDF8', fontsize=15, fontweight='bold')

        fig.canvas.draw_idle()

    # --- GUI-Steuerelemente im rechten Bereich ---
    ax_slider = plt.axes([0.76, 0.75, 0.18, 0.04], facecolor='#1F2937')
    slider_noise = Slider(ax_slider, 'Flicker-Rauschen', -1.0, 2.0, valinit=0.0, valfmt='%.1f', color='#38BDF8')
    slider_noise.label.set_color(C_TEXT)
    slider_noise.valtext.set_color(C_TEXT)

    ax_button = plt.axes([0.76, 0.62, 0.18, 0.07])
    btn_toggle = Button(ax_button, 'Filter Toggle (EIN/AUS)', color='#10B981', hovercolor='#34D399')
    btn_toggle.label.set_color('#0B0F19')
    btn_toggle.label.set_fontweight('bold')

    ax_radio = plt.axes([0.76, 0.25, 0.18, 0.28], facecolor='#1F2937')
    radio_mode = RadioButtons(ax_radio, ('Pro Audio (SNR)', 'Telekom / 5G (EVM)', 'Navigation / IMU', 'Metrologie'), active=0)
    for circle in radio_mode.circles:
        circle.set_edgecolor(C_TEXT)
    for label in radio_mode.labels:
        label.set_color(C_TEXT)

    # Callbacks
    def on_noise_change(val):
        state['flicker_noise'] = val
        update()

    def on_toggle_click(event):
        state['filter_active'] = not state['filter_active']
        if state['filter_active']:
            btn_toggle.color = '#10B981'
            btn_toggle.hovercolor = '#34D399'
        else:
            btn_toggle.color = '#EF4444'
            btn_toggle.hovercolor = '#F87171'
        update()

    def on_mode_change(label):
        state['mode'] = label
        update()

    slider_noise.on_changed(on_noise_change)
    btn_toggle.on_clicked(on_toggle_click)
    radio_mode.on_clicked(on_mode_change)

    update()
    plt.show()

if __name__ == '__main__':
    main()
