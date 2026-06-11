import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fligner, kruskal, ks_2samp, levene, mannwhitneyu

DEFAULT_INPUT = "bm_8er_konjugierte_tetraeder_bis_100m.csv"
DEFAULT_OUTPUT = "bm_octet_energy_levels_R11_47.csv"

REQUIRED_COLS = ["p_start", "q_start", "P_E", "P_A", "P_B", "P_C"]

# Achtlingstyp (Pfeil p→q in den 11/47-Kanälen); ASCII-Pfeil für Matrizen/CSV
OCTET_TYPE_ORDER = ["11->11", "11->41", "41->11", "41->41", "sonst"]
OCTET_TYPES_2X2 = ["11->11", "11->41", "41->11", "41->41"]

# Ptolemäus-Ziel E je Achtlingstyp (Homogenität 0, Kopplung 1)
PTOL_TARGET_BY_OCTET_TYPE = {
    "11->11": 0.0,
    "41->41": 0.0,
    "11->41": 1.0,
    "41->11": 1.0,
}


def _sector_1147(r: int) -> int | None:
    r = int(r) % 60
    if r == 11:
        return 11
    if r == 41:
        return 41
    return None


def _octet_type_row(p: int, q: int) -> str:
    a, b = _sector_1147(p), _sector_1147(q)
    if a is None or b is None:
        return "sonst"
    return f"{a}->{b}"


def quad_charge_mod60(r):
    """
    Vierlingsstart r mod 60:
    r=11 enthält Twin-Starts 11 und 17 -> 11-Anteil -> -1
    r=41 enthält Twin-Starts 41 und 47 -> 47-Anteil -> +1
    """
    r = int(r) % 60
    if r == 11:
        return -1
    if r == 41:
        return +1
    return 0


def level_name(R):
    if R == -2:
        return "cold_-2"
    if R == 0:
        return "neutral_0"
    if R == 2:
        return "warm_+2"
    return "other"


def ptolemy_energy(row):
    uE = row["u_E"]
    uA = row["u_A"]
    uB = row["u_B"]
    uC = row["u_C"]

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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Primzahl-Achtlinge: Ptolemäus-Energie, R_11/47, Energieniveaus. "
        "Standard: Eingabe vollständig bis Dateiende; Ausgabe-CSV immer vollständig.",
    )
    parser.add_argument(
        "-i",
        "--input",
        default=DEFAULT_INPUT,
        help=f"CSV mit {REQUIRED_COLS} (Standard: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=DEFAULT_OUTPUT,
        help=f"Ausgabe-CSV (Standard: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        metavar="N",
        help="Nur die ersten N Zeilen (Pilot). Standard: gesamte Eingabedatei (bis Dateiende).",
    )
    parser.add_argument(
        "--max-print-rows",
        type=int,
        default=6,
        metavar="K",
        help="Höchstens K Datenzeilen in der Konsole (0 = komplette Tabelle). Standard: 6.",
    )
    args = parser.parse_args()

    input_path = args.input
    if not os.path.isfile(input_path):
        print(
            f"Eingabedatei fehlt: {input_path!r}.\n"
            "Lege die Datei in den Arbeitsordner oder rufe z. B. auf:\n"
            f"  python3 {os.path.basename(__file__)} -i /pfad/zur/bm_8er_....csv",
            file=sys.stderr,
        )
        raise SystemExit(1)

    df = pd.read_csv(input_path)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        print(
            f"Fehlende Spalten: {missing}. Erforderlich: {REQUIRED_COLS}.\n"
            f"Vorhanden (Auszug): {sorted(df.columns)[:40]}",
            file=sys.stderr,
        )
        raise SystemExit(2)

    n_file = len(df)
    if args.limit is not None and args.limit > 0:
        df = df.head(int(args.limit)).copy()
    subset_note = f"erste {len(df)} von {n_file}" if len(df) < n_file else f"alle {len(df)}"
    # ------------------------------------------------------------
    # 1. Resonanzladung R_11/47
    # ------------------------------------------------------------
    df["p_mod60"] = df["p_start"] % 60
    df["q_mod60"] = df["q_start"] % 60

    df["R_11_47"] = df["p_mod60"].apply(quad_charge_mod60) + df["q_mod60"].apply(
        quad_charge_mod60
    )
    df["energy_level"] = df["R_11_47"].apply(level_name)
    df["octet_type"] = [
        _octet_type_row(int(p), int(q)) for p, q in zip(df["p_start"], df["q_start"])
    ]

    # ------------------------------------------------------------
    # 2. Log-Familienenergien
    # ------------------------------------------------------------
    df["u_E"] = np.log(df["P_E"].astype(float))
    df["u_A"] = np.log(df["P_A"].astype(float))
    df["u_B"] = np.log(df["P_B"].astype(float))
    df["u_C"] = np.log(df["P_C"].astype(float))

    # ------------------------------------------------------------
    # 3. Ptolemäus-Eigenenergie im E-A-B-C-Kreis
    # ------------------------------------------------------------
    energy_df = df.apply(ptolemy_energy, axis=1)
    df = pd.concat([df, energy_df], axis=1)

    # ------------------------------------------------------------
    # 4. Zusätzliche BM-Energiegrößen, falls vorhanden
    # ------------------------------------------------------------
    if "shadow_center_mu_geomean" in df.columns:
        df["log_mu"] = np.log(df["shadow_center_mu_geomean"].astype(float))

    if "projection_distance" in df.columns:
        df["log_projection_distance"] = np.log1p(
            df["projection_distance"].astype(float)
        )

    if "log_P_EAB" in df.columns and "log_P_CE" in df.columns:
        df["E_EAB_minus_CE"] = df["log_P_EAB"] - df["log_P_CE"]
        df["A_EAB_CE_norm"] = (df["log_P_EAB"] - df["log_P_CE"]) / (
            df["log_P_EAB"] + df["log_P_CE"]
        )

    # ------------------------------------------------------------
    # 5. Gruppenauswertung
    # ------------------------------------------------------------
    agg_spec = {
        "count": ("R_11_47", "count"),
        "median_ptolemy": ("Ptolemy_energy_norm", "median"),
        "mean_ptolemy": ("Ptolemy_energy_norm", "mean"),
        "std_ptolemy": ("Ptolemy_energy_norm", "std"),
    }
    if "projection_distance" in df.columns:
        agg_spec["median_projection"] = ("projection_distance", "median")
        agg_spec["mean_projection"] = ("projection_distance", "mean")

    summary = df.groupby("energy_level", sort=True).agg(**agg_spec).reset_index()

    type_spec = {
        "count": ("R_11_47", "count"),
        "median_ptolemy": ("Ptolemy_energy_norm", "median"),
        "mean_ptolemy": ("Ptolemy_energy_norm", "mean"),
        "std_ptolemy": ("Ptolemy_energy_norm", "std"),
    }
    if "projection_distance" in df.columns:
        type_spec["median_projection"] = ("projection_distance", "median")
        type_spec["mean_projection"] = ("projection_distance", "mean")
    if "shadow_center_mu_geomean" in df.columns:
        type_spec["median_mu"] = ("shadow_center_mu_geomean", "median")
        type_spec["mean_mu"] = ("shadow_center_mu_geomean", "mean")

    summary_type = (
        df.groupby("octet_type", sort=False)
        .agg(**type_spec)
        .reindex(pd.Index(OCTET_TYPE_ORDER, name="octet_type"))
    )

    n = len(df)
    print(
        f"\n===== Primzahl-Achtlinge: Energieniveaus ({subset_note}, n={n}) ====="
    )
    cols = [
        "p_start",
        "q_start",
        "p_mod60",
        "q_mod60",
        "octet_type",
        "R_11_47",
        "energy_level",
        "Ptolemy_energy_norm",
        "Ptolemy_defect",
        "E_diag",
        "E_1",
        "E_2",
    ]
    if "projection_distance" in df.columns:
        cols.append("projection_distance")
    out = df[cols]
    mpr = args.max_print_rows
    if mpr and mpr > 0 and len(out) > mpr:
        print(out.head(mpr // 2).to_string(index=False))
        print(f"\n[…] {len(out) - mpr} Zeilen ausgelassen […]\n")
        print(out.tail(mpr - mpr // 2).to_string(index=False))
        print(
            f"\n(Tabelle in der Konsole auf {mpr} Zeilen begrenzt; vollständig in CSV.)"
        )
    else:
        print(out.to_string(index=False))

    print("\n===== Zusammenfassung nach Energieniveau =====")
    print(summary.to_string(index=False))

    print("\n===== Ptolemäische Zustandsmatrix nach Achtlingstyp =====")
    print(summary_type.reset_index().to_string(index=False))

    def med_for(t: str, col: str) -> float:
        s = df[df["octet_type"] == t][col].dropna()
        return float(s.median()) if len(s) else float("nan")

    M_ptol = np.array(
        [
            [
                med_for("11->11", "Ptolemy_energy_norm"),
                med_for("11->41", "Ptolemy_energy_norm"),
            ],
            [
                med_for("41->11", "Ptolemy_energy_norm"),
                med_for("41->41", "Ptolemy_energy_norm"),
            ],
        ],
        dtype=float,
    )
    print("\nMedian-Matrizen 2x2 (Zeilen: p-Kanal 11,41; Spalten: q-Kanal 11,41)")
    print("Typ-Liste:", OCTET_TYPES_2X2)
    print("\nPtolemäus-Matrix (Mediane von Ptolemy_energy_norm):")
    print(M_ptol)

    if "projection_distance" in df.columns:
        M_proj = np.array(
            [
                [
                    med_for("11->11", "projection_distance"),
                    med_for("11->41", "projection_distance"),
                ],
                [
                    med_for("41->11", "projection_distance"),
                    med_for("41->41", "projection_distance"),
                ],
            ],
            dtype=float,
        )
        print("\nProjektionsdistanz-Matrix (Mediane):")
        print(M_proj)
    else:
        print(
            "\n(Keine Spalte projection_distance — Projektionsdistanz-Matrix entfällt.)"
        )

    # ------------------------------------------------------------
    # 5b. Diagonal vs. off-diagonal: Strenge der ptolemäischen Zustände
    # ------------------------------------------------------------
    df["state_class"] = np.where(
        df["octet_type"].isin(["11->11", "41->41"]),
        "diagonal_eigenstate",
        np.where(
            df["octet_type"].isin(["11->41", "41->11"]),
            "offdiagonal_coupling",
            "other",
        ),
    )
    df["ptol_target"] = np.where(
        df["state_class"] == "diagonal_eigenstate",
        0.0,
        np.where(df["state_class"] == "offdiagonal_coupling", 1.0, np.nan),
    )
    df["ptol_deviation_from_target"] = (
        df["Ptolemy_energy_norm"] - df["ptol_target"]
    ).abs()

    main = df[
        df["state_class"].isin(["diagonal_eigenstate", "offdiagonal_coupling"])
    ].copy()

    diag = main[main["state_class"] == "diagonal_eigenstate"][
        "ptol_deviation_from_target"
    ].dropna()
    off = main[main["state_class"] == "offdiagonal_coupling"][
        "ptol_deviation_from_target"
    ].dropna()

    print("\n===== Strenge der ptolemäischen Zustände =====")
    print("N diagonal:", len(diag))
    print("N offdiagonal:", len(off))

    print("\nDeviation from target:")
    print("Median diagonal |E-0|:", float(diag.median()))
    print("Mean diagonal   |E-0|:", float(diag.mean()))
    print("Std diagonal    |E-0|:", float(diag.std()))

    print("Median offdiag  |E-1|:", float(off.median()))
    print("Mean offdiag    |E-1|:", float(off.mean()))
    print("Std offdiag     |E-1|:", float(off.std()))

    if len(diag) and len(off):
        mw_dev = mannwhitneyu(diag, off, alternative="two-sided")
        ks_dev = ks_2samp(diag, off)

        print("\nMann-Whitney deviation diagonal vs. offdiagonal:")
        print("U =", mw_dev.statistic, "p =", mw_dev.pvalue)

        print("\nKS deviation diagonal vs. offdiagonal:")
        print("D =", ks_dev.statistic, "p =", ks_dev.pvalue)

        diag_raw = main[main["state_class"] == "diagonal_eigenstate"][
            "Ptolemy_energy_norm"
        ].dropna()
        off_raw = main[main["state_class"] == "offdiagonal_coupling"][
            "Ptolemy_energy_norm"
        ].dropna()

        if len(diag_raw) > 1 and len(off_raw) > 1:
            lev = levene(diag_raw, off_raw, center="median")
            flig = fligner(diag_raw, off_raw)
            print("\nStreuungstest Rohwerte E_Ptol:")
            print("Levene median-centered: stat =", lev.statistic, "p =", lev.pvalue)
            print("Fligner-Killeen:        stat =", flig.statistic, "p =", flig.pvalue)
        else:
            print(
                "\n(Levene/Fligner übersprungen: zu wenige Rohwerte in einer Gruppe.)"
            )
    else:
        print(
            "\n(Mann-Whitney/KS übersprungen: eine Gruppe leer.)",
        )

    # ------------------------------------------------------------
    # 5c. Matrix-Strenge je Achtlingstyp (Abweichung vom typweisen Ziel)
    # ------------------------------------------------------------
    strictness_by_type = df[df["octet_type"].isin(OCTET_TYPES_2X2)].copy()
    strictness_by_type["ptol_target"] = strictness_by_type["octet_type"].map(
        PTOL_TARGET_BY_OCTET_TYPE
    )
    strictness_by_type["target_deviation"] = (
        strictness_by_type["Ptolemy_energy_norm"] - strictness_by_type["ptol_target"]
    ).abs()

    summary_strict = (
        strictness_by_type.groupby("octet_type", sort=False)
        .agg(
            count=("target_deviation", "count"),
            median_E=("Ptolemy_energy_norm", "median"),
            mean_E=("Ptolemy_energy_norm", "mean"),
            std_E=("Ptolemy_energy_norm", "std"),
            median_deviation=("target_deviation", "median"),
            mean_deviation=("target_deviation", "mean"),
            max_deviation=("target_deviation", "max"),
        )
        .reset_index()
    )
    summary_strict = summary_strict.set_index("octet_type").reindex(OCTET_TYPES_2X2).reset_index()

    print("\n===== Matrix-Strenge nach Achtlingstyp =====")
    print(summary_strict.to_string(index=False))

    # ------------------------------------------------------------
    # 6. Explorative Signifikanztests
    # ------------------------------------------------------------
    groups = {
        level: sub["Ptolemy_energy_norm"].dropna().to_numpy()
        for level, sub in df.groupby("energy_level", sort=True)
    }

    print("\n===== Explorative Tests =====")

    if (
        "cold_-2" in groups
        and "warm_+2" in groups
        and len(groups["cold_-2"]) > 0
        and len(groups["warm_+2"]) > 0
    ):
        mw = mannwhitneyu(
            groups["cold_-2"],
            groups["warm_+2"],
            alternative="two-sided",
        )
        print("Mann-Whitney cold_-2 vs warm_+2:")
        print("U =", mw.statistic, "p =", mw.pvalue)

    level_keys_sorted = sorted(groups.keys())
    arrays_ordered = [groups[k] for k in level_keys_sorted]
    if len(arrays_ordered) >= 2 and all(len(v) > 0 for v in arrays_ordered):
        kw = kruskal(*arrays_ordered)
        print("Kruskal-Wallis über Energieniveaus (Reihenfolge):")
        print(level_keys_sorted)
        print("H =", kw.statistic, "p =", kw.pvalue)

    # ------------------------------------------------------------
    # 7. Export
    # ------------------------------------------------------------
    df.to_csv(args.output, index=False)
    print("\nGeschrieben:", os.path.abspath(args.output))


if __name__ == "__main__":
    main()
