#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Energiedoku - Kalibrierungsprogramm der Arithmetischen Dämpfung
Modell: Quaternionisches Primzahlmodell (Bamberger Protokoll)
Lauffähig mit Python 3 (ohne SageMath-Abhängigkeit)
"""

import math
import argparse
from decimal import Decimal, getcontext

def kalibriere_daempfung():
    # 1. Fundamentale Konstanten im Modell
    pi_hochpraezise = math.pi
    phi_golden = (1 + math.sqrt(5)) / 2
    
    # Energieniveaus
    E_urf = pi_hochpraezise / 6          # Urfeld-Energie (~0.523598)
    E_ziel = math.log(phi_golden)        # Ziel-Energie (~0.481211)
    Delta_E_ges = E_ziel - E_urf         # Erforderliche Gesamtdämpfung
    
    print("="*65)
    print("      #ENERGIEDOKU - KALIBRIERUNG DER ARITHMETISCHEN DÄMPFUNG")
    print("="*65)
    print(f"Urfeld-Energie (pi/6):        {E_urf:.15f}")
    print(f"Ziel-Energie (ln(phi)):       {E_ziel:.15f}")
    print(f"Erforderliche Gesamtdämpfung: {Delta_E_ges:.15f}")
    print("-"*65)
    
    # 3. Berechnung der diskreten Kraft-Komponenten (Basiskopplungen)
    # k=1: Tschebyscheff-Bias
    bias_k1 = math.exp(-pi_hochpraezise)
    # k=2: Quadratische Projektion
    proj_k2 = 0.5 * math.exp(-2 * pi_hochpraezise)
    # k=5: Ikosaeder-Fluktuation (Rigiditäts-Anker)
    iko_k5 = math.exp(-5 * pi_hochpraezise)
    
    print("KOPPLUNGS-KOMPONENTEN (BASIS-KALIBRIERUNG):")
    print(f"  1. Tschebyscheff-Bias (e^-pi):       {bias_k1:.15f}")
    print(f"  2. Quadr. Projektion (0.5 * e^-2pi): {proj_k2:.15f}")
    print(f"  3. Ikosaeder-Fluktuation (e^-5pi):   {iko_k5:.15f}")
    print("-"*65)
    
    # 4. Simulation der Skalenkaskade (Iterativer Kollaps)
    print("SIMULATION DER SKALENKASKADE (KOLLAPS GEGEN DIE STRUKTUR):")
    
    E_aktuell = E_urf
    
    # Exakte Interferenz-Kaskade:
    # Delta_k = ln((1 + e^(-5(2k-1)pi)) / (1 + e^(-(2k-1)pi)))
    # Die Vorzeichen sind in dieser Form bereits korrekt codiert.
    # Bis k=5 ist die Zielenergie auf Maschinenpräzision getroffen.
    for k in range(1, 6):
        num = 1 + math.exp(-5 * (2 * k - 1) * pi_hochpraezise)
        den = 1 + math.exp(-(2 * k - 1) * pi_hochpraezise)
        delta_k = math.log(num / den)
        E_aktuell += delta_k
        restfehler = E_aktuell - E_ziel
        
        print(f"\nSchnitt k={k} [Interferenz-Modus]:")
        print(f"  Aktuelle Energie: {E_aktuell:.15f}")
        print(f"  Restfehler:       {restfehler:.6e}")
        
    print("="*65)
    print("KALIBRIERUNG ERGIBT: SÄTTIGUNG BEI k=5 ERREICHT.")
    print("Das System ist im Hurwitz-Gitter eingefroren.")
    print("="*65)


def zeige_high_precision_vergleich(max_k=10, praezision=120):
    getcontext().prec = praezision
    D = Decimal

    # Feste PI-Referenz mit ausreichender Länge für hohe Präzision.
    pi_dec = D(
        "3.14159265358979323846264338327950288419716939937510"
        "58209749445923078164062862089986280348253421170679"
    )
    phi_dec = (D(1) + D(5).sqrt()) / D(2)
    log_phi_dec = phi_dec.ln()

    log_phi_float = math.log((1 + math.sqrt(5)) / 2)

    print("\n" + "=" * 65)
    print("HIGH-PRECISION VERGLEICH (float vs. decimal)")
    print("=" * 65)
    print(f"Kontext-Präzision (decimal): {praezision} Stellen")
    print("k | float_rest            | decimal_rest")
    print("-" * 65)

    for kmax in range(1, max_k + 1):
        # Float-Lauf
        e_float = math.pi / 6
        for k in range(1, kmax + 1):
            num_f = 1 + math.exp(-5 * (2 * k - 1) * math.pi)
            den_f = 1 + math.exp(-(2 * k - 1) * math.pi)
            e_float += math.log(num_f / den_f)
        rest_float = e_float - log_phi_float

        # Decimal-Lauf
        e_dec = pi_dec / D(6)
        for k in range(1, kmax + 1):
            m = D(2 * k - 1)
            num_d = D(1) + (-D(5) * m * pi_dec).exp()
            den_d = D(1) + (-m * pi_dec).exp()
            e_dec += (num_d / den_d).ln()
        rest_dec = e_dec - log_phi_dec

        print(f"{kmax:2d} | {rest_float: .3e} | {rest_dec:.3E}")

    print("-" * 65)
    print("Interpretation: float sättigt am Numerik-Floor, decimal konvergiert weiter.")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Kalibrierung der arithmetischen Dämpfung (#Energiedoku)."
    )
    parser.add_argument(
        "--high-precision",
        action="store_true",
        help="Zeigt zusätzlich einen float-vs-decimal Präzisionsvergleich.",
    )
    parser.add_argument(
        "--max-k",
        type=int,
        default=10,
        help="Anzahl der Schnitte für den High-Precision-Vergleich (default: 10).",
    )
    parser.add_argument(
        "--prec",
        type=int,
        default=120,
        help="Decimal-Präzision in Stellen für den High-Precision-Modus (default: 120).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    kalibriere_daempfung()
    if args.high_precision:
        zeige_high_precision_vergleich(max_k=args.max_k, praezision=args.prec)