import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, mannwhitneyu, pearsonr, spearmanr


REQUIRED_BASE_COLS = [
    "p_start",
    "P_E",
    "P_A",
    "P_B",
    "P_C",
    "P_EAB",
    "P_CE",
    "log_P_EAB",
    "log_P_CE",
    "projection_distance",
    "start_gap",
    "shadow_center_mu_geomean",
]


def _require_columns(df: pd.DataFrame, cols: list[str]) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        print("Fehlende Spalten in der CSV:", missing, file=sys.stderr)
        print("Vorhandene Spalten (Auszug):", sorted(list(df.columns))[:50], file=sys.stderr)
        raise SystemExit(2)


def run(input_path: str, output_path: str) -> None:
    if not os.path.exists(input_path):
        print(
            f"Eingabedatei fehlt: {input_path!r}.\n"
            "Lege die Datei in diesen Ordner oder rufe das Skript mit\n"
            f"  python3 {os.path.basename(__file__)} --input /pfad/zur/csv\n"
            "auf."
        )
        raise SystemExit(1)

    df = pd.read_csv(input_path)
    _require_columns(df, REQUIRED_BASE_COLS)

    # Anfangsartefakt optional entfernen
    df = df[df["p_start"] > 60].copy()

    # Kanäle
    df["p_mod60"] = df["p_start"] % 60

    # Logs der Familienprodukte
    for col in ["P_E", "P_A", "P_B", "P_C"]:
        df[f"log_{col}"] = np.log(df[col].astype(float))

    # Energie-Asymmetrie
    df["E_EAB_minus_CE"] = df["log_P_EAB"] - df["log_P_CE"]
    df["E_EAB_plus_CE"] = df["log_P_EAB"] + df["log_P_CE"]
    df["A_EAB_CE_norm"] = df["E_EAB_minus_CE"] / df["E_EAB_plus_CE"]

    # Familien-Torsionen
    df["T_EB_minus_AC"] = (df["log_P_E"] + df["log_P_B"]) - (df["log_P_A"] + df["log_P_C"])
    df["T_A_minus_C"] = df["log_P_A"] - df["log_P_C"]
    df["T_E_minus_B"] = df["log_P_E"] - df["log_P_B"]

    # Normierte Projektionsdistanzen
    df["proj_dist_over_gap"] = df["projection_distance"] / df["start_gap"]
    df["proj_dist_over_mu"] = df["projection_distance"] / df["shadow_center_mu_geomean"]

    # Absolutwerte der Energien
    df["abs_A_EAB_CE_norm"] = df["A_EAB_CE_norm"].abs()
    df["abs_T_EB_minus_AC"] = df["T_EB_minus_AC"].abs()
    df["abs_T_A_minus_C"] = df["T_A_minus_C"].abs()
    df["abs_T_E_minus_B"] = df["T_E_minus_B"].abs()

    energy_cols = [
        "E_EAB_minus_CE",
        "A_EAB_CE_norm",
        "abs_A_EAB_CE_norm",
        "T_EB_minus_AC",
        "abs_T_EB_minus_AC",
        "T_A_minus_C",
        "abs_T_A_minus_C",
        "T_E_minus_B",
        "abs_T_E_minus_B",
    ]

    target_cols = [
        "projection_distance",
        "proj_dist_over_gap",
        "proj_dist_over_mu",
    ]

    print("============================================================")
    print("BM Schattenenergie-Analyse")
    print("============================================================")
    print("N rows:", len(df))
    print("Kanäle p mod 60:")
    print(df["p_mod60"].value_counts().sort_index())
    print()

    print("------------------------------------------------------------")
    print("Korrelation Energie -> Projektdistanz")
    print("------------------------------------------------------------")

    rows: list[dict] = []

    for e in energy_cols:
        for t in target_cols:
            x = df[e].replace([np.inf, -np.inf], np.nan)
            y = df[t].replace([np.inf, -np.inf], np.nan)
            mask = x.notna() & y.notna()

            if int(mask.sum()) < 3:
                continue

            pear = pearsonr(x[mask], y[mask])
            spear = spearmanr(x[mask], y[mask])

            rows.append(
                {
                    "energy": e,
                    "target": t,
                    "pearson_r": pear.statistic,
                    "pearson_p": pear.pvalue,
                    "spearman_r": spear.statistic,
                    "spearman_p": spear.pvalue,
                }
            )

    corr_df = pd.DataFrame(rows)
    if corr_df.empty:
        print("Keine brauchbaren Zeilenpaare für Korrelationen.")
    else:
        print(corr_df.sort_values("spearman_p").to_string(index=False))

    print()
    print("------------------------------------------------------------")
    print("Quantiltest: balancierte vs. unbalancierte Schatten")
    print("------------------------------------------------------------")

    q_low = df["abs_A_EAB_CE_norm"].quantile(0.25)
    q_high = df["abs_A_EAB_CE_norm"].quantile(0.75)

    low = df[df["abs_A_EAB_CE_norm"] <= q_low]
    high = df[df["abs_A_EAB_CE_norm"] >= q_high]

    print("abs_A_EAB_CE_norm:")
    print("Q25:", q_low, "Q75:", q_high)
    print("low N:", len(low), "high N:", len(high))
    print("Median proj low:", low["projection_distance"].median())
    print("Median proj high:", high["projection_distance"].median())
    print("Mean proj low:", low["projection_distance"].mean())
    print("Mean proj high:", high["projection_distance"].mean())

    if len(low) and len(high):
        ks = ks_2samp(low["projection_distance"], high["projection_distance"])
        mw = mannwhitneyu(
            low["projection_distance"],
            high["projection_distance"],
            alternative="two-sided",
        )
        print("KS:", ks.statistic, "p:", ks.pvalue)
        print("MW p:", mw.pvalue)
    else:
        print("Zu wenig Daten für KS/MW-Test.")

    print()
    print("------------------------------------------------------------")
    print("Kanalvergleich 11 vs 41 für Energiegrößen")
    print("------------------------------------------------------------")

    d11 = df[df["p_mod60"] == 11]
    d41 = df[df["p_mod60"] == 41]

    for col in ["projection_distance"] + energy_cols:
        if d11[col].notna().sum() and d41[col].notna().sum():
            ks = ks_2samp(d11[col], d41[col])
            mw = mannwhitneyu(d11[col], d41[col], alternative="two-sided")
        else:
            ks = None
            mw = None
        print(
            col,
            "median11=",
            d11[col].median(),
            "median41=",
            d41[col].median(),
            "KS_p=",
            (ks.pvalue if ks is not None else np.nan),
            "MW_p=",
            (mw.pvalue if mw is not None else np.nan),
        )

    df.to_csv(output_path, index=False)
    print()
    print("Geschrieben:", output_path)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="BM Schattenenergie-Analyse")
    p.add_argument(
        "--input",
        default="bm_8er_konjugierte_tetraeder_bis_100m.csv",
        help="Eingabedatei (CSV) mit p_start, Produkte, Projektdistanz, ...",
    )
    p.add_argument(
        "--output",
        default="bm_shadow_energy_analysis.csv",
        help="Ausgabedatei",
    )
    return p


def main() -> None:
    args = build_parser().parse_args()
    run(args.input, args.output)


if __name__ == "__main__":
    main()
