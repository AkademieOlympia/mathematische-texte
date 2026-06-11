"""
Alle Primzahl-Achtlinge: Typ (11/47-Kanäle), R_11/47, Ptolemäus-Energie,
Zusammenfassung, explorative Gruppentests (Schwerpunkt: homogen vs. gemischt).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

DEFAULT_INPUT = "bm_8er_konjugierte_tetraeder_bis_100m.csv"
DEFAULT_OUT_DETAIL = "bm_octet_energy_levels_full.csv"
DEFAULT_OUT_SUMMARY = "bm_octet_energy_levels_full_summary_by_type.csv"
DEFAULT_OUT_TESTS = "bm_octet_energy_levels_full_tests.json"

REQUIRED_COLS = ["p_start", "q_start", "P_E", "P_A", "P_B", "P_C"]

TYPES_ORDER = ["11→11", "11→41", "41→11", "41→41"]


def quad_charge_mod60(r: int) -> int:
    """11-Anteil -1, 47/41-Seite +1 (mod-60-Familie)."""
    r = int(r) % 60
    if r == 11:
        return -1
    if r == 41:
        return 1
    return 0


def sector_label_1147(r: int) -> int | None:
    """Für p,q Start: nur 11 bzw. 41 (Vierlingskanäle in dieser Tabelle)."""
    r = int(r) % 60
    if r == 11:
        return 11
    if r == 41:
        return 41
    return None


def arrow_type(p_start: int, q_start: int) -> str | None:
    a = sector_label_1147(p_start)
    b = sector_label_1147(q_start)
    if a is None or b is None:
        return None
    return f"{a}→{b}"


def mixing_class(arrow: str) -> str | None:
    if arrow in ("11→11", "41→41"):
        return "homogen"
    if arrow in ("11→41", "41→11"):
        return "gemischt"
    return None


def ptolemy_row(row) -> pd.Series:
    uE, uA, uB, uC = row["u_E"], row["u_A"], row["u_B"], row["u_C"]
    d_EA = abs(uE - uA)
    d_AB = abs(uA - uB)
    d_BC = abs(uB - uC)
    d_CE = abs(uC - uE)
    d_EB = abs(uE - uB)
    d_AC = abs(uA - uC)

    E_diag = d_EB * d_AC
    E_1 = d_EA * d_BC
    E_2 = d_AB * d_CE
    Pi = E_diag - E_1 - E_2
    denom = E_diag + E_1 + E_2
    eps = abs(Pi) / denom if denom != 0 else np.nan
    return pd.Series(
        {
            "E_diag": E_diag,
            "E_1": E_1,
            "E_2": E_2,
            "Ptolemy_defect": Pi,
            "Ptolemy_energy_norm": eps,
        }
    )


def mann_whitney_safe(
    a: np.ndarray,
    b: np.ndarray,
    label: str,
) -> dict[str, Any] | None:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 1 or len(b) < 1:
        print(
            f"  {label}: übersprungen (n1={len(a)}, n2={len(b)}).",
            file=sys.stderr,
        )
        return None
    r = mannwhitneyu(a, b, alternative="two-sided")
    d = {
        "test": "Mann-Whitney",
        "label": label,
        "n1": int(len(a)),
        "n2": int(len(b)),
        "U": float(r.statistic),
        "p": float(r.pvalue),
    }
    print(f"  {label}")
    print(f"    n1 = {d['n1']}, n2 = {d['n2']}, U = {d['U']}, p = {d['p']}")
    return d


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Achtlinge: Typ 11/47, R_11/47, Ptolemäus-Energie, Gruppentests."
    )
    parser.add_argument("-i", "--input", default=DEFAULT_INPUT)
    parser.add_argument("-o", "--output", default=DEFAULT_OUT_DETAIL, help="Volldatensatz")
    parser.add_argument(
        "--summary",
        default=DEFAULT_OUT_SUMMARY,
        help="Aggregat nach Pfeiltyp",
    )
    parser.add_argument(
        "--tests-out",
        default=DEFAULT_OUT_TESTS,
        help="JSON mit Testergebnissen (oder 'skip' für kein File)",
    )
    parser.add_argument(
        "--min-p-start",
        type=int,
        default=60,
        help="Zeilen mit p_start < Wert weglassen (Voreinstellung: 60, Anfangsartefakt). 0 = kein Filter.",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Eingabedatei fehlt: {args.input!r}", file=sys.stderr)
        raise SystemExit(1)

    df = pd.read_csv(args.input)
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        print(f"Fehlende Spalten: {missing}", file=sys.stderr)
        raise SystemExit(2)

    n_raw = len(df)
    if args.min_p_start > 0:
        df = df[df["p_start"] >= args.min_p_start].copy()
    n_after = len(df)
    print("============================================================")
    print("bm_octet_energy_levels_full")
    print("Eingabe:", os.path.abspath(args.input))
    print(f"Zeilen: {n_raw} (nach min_p_start>={args.min_p_start}: {n_after})")

    df["p_mod60"] = (df["p_start"] % 60).astype(np.int64)
    df["q_mod60"] = (df["q_start"] % 60).astype(np.int64)
    df["R_11_47"] = df["p_mod60"].map(quad_charge_mod60) + df["q_mod60"].map(
        quad_charge_mod60
    )
    df["typ_11_47_pfeil"] = [
        arrow_type(int(p), int(q)) for p, q in zip(df["p_start"], df["q_start"])
    ]
    df["mischung_homogen_gemischt"] = df["typ_11_47_pfeil"].map(mixing_class)

    bad = df["typ_11_47_pfeil"].isna()
    n_bad = int(bad.sum())
    if n_bad:
        print(
            f"Warnung: {n_bad} Zeilen ohne klaren 11/41-Pfeiltyp (p,q nicht in "
            f"{{11,41}} mod 60) — in Ausgabe markiert, Tests nur auf Gültige.",
        )

    # Logs
    for col in ["P_E", "P_A", "P_B", "P_C"]:
        df[col] = df[col].astype(float)
    df["u_E"] = np.log(df["P_E"])
    df["u_A"] = np.log(df["P_A"])
    df["u_B"] = np.log(df["P_B"])
    df["u_C"] = np.log(df["P_C"])

    extra = df.apply(ptolemy_row, axis=1)
    df = pd.concat([df, extra], axis=1)

    y = "Ptolemy_energy_norm"
    df_test = df[df["typ_11_47_pfeil"].notna()].copy()
    df_test = df_test[np.isfinite(df_test[y])].copy()

    # --- Zusammenfassung nach Pfeiltyp (alle vier Zeilen, fehlend → NaN)
    summary = (
        df_test.groupby("typ_11_47_pfeil", sort=False)[y]
        .agg(n="count", median="median", mean="mean", std="std")
        .reindex(pd.Index(TYPES_ORDER, name="typ_11_47_pfeil"))
    )

    # --- Mischung
    sum_mix = (
        df_test.dropna(subset=["mischung_homogen_gemischt"])
        .groupby("mischung_homogen_gemischt")[y]
        .agg(n="count", median="median", mean="mean", std="std")
    )
    print("\n--- Ptolemäus-Energie (norm) nach Mischung ---")
    print(sum_mix.to_string())
    print("\n--- Ptolemäus-Energie (norm) nach Pfeiltyp ---")
    print(summary.to_string())

    # === Tests
    print("\n" + "=" * 60)
    print("Statistische Tests (Ptolemäus-Energie norm., explorativ)")
    print("=" * 60)

    test_results: list[dict[str, Any]] = []

    # Kruskal-Wallis: vier Pfeiltypen
    groups_4 = [
        df_test.loc[df_test["typ_11_47_pfeil"] == t, y].to_numpy() for t in TYPES_ORDER
    ]
    groups_4 = [g[np.isfinite(g)] for g in groups_4]
    non_empty_4 = [(TYPES_ORDER[i], groups_4[i]) for i in range(4) if len(groups_4[i])]
    print("\nKruskal-Wallis (vier Pfeiltypen):")
    print("  n pro Typ:", {t: len(g) for t, g in non_empty_4})
    if len(non_empty_4) >= 2:
        kw = kruskal(*[g for _, g in non_empty_4])
        rec = {
            "test": "Kruskal-Wallis",
            "groups": "vier Pfeiltypen",
            "n": {t: int(len(g)) for t, g in non_empty_4},
            "H": float(kw.statistic),
            "p": float(kw.pvalue),
        }
        test_results.append(rec)
        print(f"  H = {rec['H']}, p = {rec['p']}")
    else:
        print("  (zu wenige nicht-leere Gruppen)")

    # Homogen vs gemischt (Kernvergleich)
    hom = df_test[df_test["mischung_homogen_gemischt"] == "homogen"][y].to_numpy()
    gem = df_test[df_test["mischung_homogen_gemischt"] == "gemischt"][y].to_numpy()
    print(
        "\n*** Mann-Whitney: homogen {11→11, 41→41} vs. gemischt {11→41, 41→11} ***"
    )
    m1 = mann_whitney_safe(
        hom,
        gem,
        "homogen vs. gemischt",
    )
    if m1:
        m1["groups_detail"] = "11→11 & 41→41  vs.  11→41 & 41→11"
        test_results.append(m1)

    s11 = df_test[df_test["typ_11_47_pfeil"] == "11→11"][y].to_numpy()
    s44 = df_test[df_test["typ_11_47_pfeil"] == "41→41"][y].to_numpy()
    print("\nMann-Whitney: 11→11 vs. 41→41 (innerhalb homogen)")
    m2 = mann_whitney_safe(s11, s44, "11→11 vs. 41→41")
    if m2:
        test_results.append(m2)

    s14 = df_test[df_test["typ_11_47_pfeil"] == "11→41"][y].to_numpy()
    s41 = df_test[df_test["typ_11_47_pfeil"] == "41→11"][y].to_numpy()
    print("\nMann-Whitney: 11→41 vs. 41→11 (innerhalb gemischt)")
    m3 = mann_whitney_safe(s14, s41, "11→41 vs. 41→11")
    if m3:
        test_results.append(m3)

    df.to_csv(args.output, index=False)
    summary.to_csv(args.summary, encoding="utf-8")
    print("\n--- Export ---")
    print("  Voll:"    , os.path.abspath(args.output))
    print("  Summary:" , os.path.abspath(args.summary))

    if args.tests_out and args.tests_out.lower() != "skip":
        payload = {
            "n_raw": n_raw,
            "n_after_min_p": n_after,
            "min_p_start": int(args.min_p_start),
            "n_valid_arrow_type": int(len(df_test)),
            "n_unclassified": n_bad,
            "tests": test_results,
        }
        pth = args.tests_out
        with open(pth, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)
        print("  Tests: ", os.path.abspath(pth))

    print("============================================================")


if __name__ == "__main__":
    main()
