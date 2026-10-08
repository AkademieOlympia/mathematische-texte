#!/usr/bin/env python3
"""
Zeta–Keating Provenance Gate
============================

Symbolische Propagation der finite-height-Korrektur der Riemann-Zeta-
Zweipunktfunktion in die Number Variance Sigma^2(L).

Ausgangspunkt (Bogomolny–Bohigas–Leboeuf–Monastra 2006, Eq. 16):
  R2(s;E) = R2_GUE(s)
            - Lambda/(pi^2 rho^2) sin^2(pi s)
            - Q/(2 pi^2 rho^3) s sin(2 pi s)
            + O(rho^-4),

mit
  rho = (1/(2 pi)) log(E/(2 pi)),
  Lambda = gamma_0^2 + 2 gamma_1 + c_0 ~ 1.57314,
  Q = sum_p log(p)^3/(p-1)^2.

Damit, für lambda = log(E/(2 pi)):
  f2(s) = -4 Lambda sin^2(pi s),
  f3(s) = -4 pi Q s sin(2 pi s),

und
  Delta Sigma^2(L;E)
    = 2 int_0^L (L-s) [R2(s;E)-R2_GUE(s)] ds
    = C2(L)/lambda^2 + C3(L)/lambda^3 + O(lambda^-4).

Dieses Modul leitet C2(L), C3(L) symbolisch her und wertet speziell L=1 aus.

WICHTIG:
Die zugrunde liegende finite-height-Zeta-Korrelationsformel ist heuristisch,
nicht als rigoroser Satz bewiesen. Das Modul prüft daher eine
THEORY-CONDITIONAL prediction, keinen mathematischen Beweis.
"""

import argparse
import numpy as np
import pandas as pd
import sympy as sp

s, L = sp.symbols("s L", positive=True, real=True)
Lambda, Q = sp.symbols("Lambda Q", real=True)
pi = sp.pi

# GUE-Limit
R2_gue = 1 - (sp.sin(pi*s)/(pi*s))**2

# Bogomolny et al. finite-height coefficients after rho -> lambda conversion
f2 = -4 * Lambda * sp.sin(pi*s)**2
f3 = -4 * pi * Q * s * sp.sin(2*pi*s)

C2_L = sp.simplify(2 * sp.integrate((L-s)*f2, (s, 0, L)))
C3_L = sp.simplify(2 * sp.integrate((L-s)*f3, (s, 0, L)))

C2_1 = sp.simplify(C2_L.subs(L, 1))
C3_1 = sp.simplify(C3_L.subs(L, 1))

def pred_delta_sigma2(T, lambda_const=1.57314):
    lam = np.log(np.asarray(T, dtype=float)/(2*np.pi))
    return -2.0*lambda_const/(lam**2)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=None,
                    help="Optional CSV from zeta_jitter_scaling_audit.csv")
    ap.add_argument("--Lambda", dest="Lambda_val", type=float, default=1.57314)
    args = ap.parse_args()

    print("="*78)
    print(" ZETA–KEATING PROVENANCE GATE: R2 -> Sigma^2")
    print("="*78)
    print("f2(s) =", f2)
    print("f3(s) =", f3)
    print()
    print("C2(L) =", C2_L)
    print("C3(L) =", C3_L)
    print()
    print("SPECIAL SECTOR L=1")
    print("C2(1) =", C2_1)
    print("C3(1) =", C3_1)
    print()
    print("Therefore:")
    print("  Delta Sigma^2(1;T)")
    print("    = -2 Lambda / log(T/(2*pi))^2 + O(log(T)^-4)")
    print("  because the full O(log^-3) contribution integrates to zero at L=1.")
    print()

    if args.csv:
        df = pd.read_csv(args.csv)
        if "gamma_N" not in df or "Delta" not in df:
            raise ValueError("CSV must contain columns gamma_N and Delta.")
        df = df.copy()
        df["lambda_T"] = np.log(df["gamma_N"]/(2*np.pi))
        df["Delta_BK_leading"] = pred_delta_sigma2(
            df["gamma_N"].to_numpy(), args.Lambda_val
        )
        df["emp_over_theory"] = df["Delta"]/df["Delta_BK_leading"]
        df["residual_emp_minus_BK"] = df["Delta"] - df["Delta_BK_leading"]

        cols = [
            "N", "gamma_N", "Delta",
            "Delta_BK_leading", "emp_over_theory",
            "residual_emp_minus_BK"
        ]
        print(df[cols].to_string(index=False, float_format=lambda x: f"{x:.9f}"))
        print()
        print("CLAIM WALL")
        print("  [EXACT SYMBOLIC] C2(1) = -2 Lambda.")
        print("  [EXACT SYMBOLIC] C3(1) = 0 for the quoted asymptotic R2 expansion.")
        print("  [THEORY-CONDITIONAL] Leading finite-height prediction is log^-2.")
        print("  [DATA TEST] Compare the fixed coefficient -2 Lambda against zeros.")
        print("  [NOT CLAIMED] No alpha or c_eff mapping is used here.")

if __name__ == "__main__":
    main()
