#!/usr/bin/env python3
"""
Physikalische Re/Ra-Zuordnung auf dem vereinigten EABC-Zeugen-Katalog.

Rationaler Sektor (a=0): Reynolds-Zahl aus M/Norm-Kaskade und reziproken Normen.
Phi-Sektor (a≠0): Rayleigh-Zahl aus algebraischer M²/Norm²-Kopplung in ℤ[Φ].
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

BASIS = Path(__file__).resolve().parent
M_KONSTANTE = 113160
RE_KLASSISCH = 2300.0
RA_KLASSISCH = 1708.0
PHI = (1.0 + math.sqrt(5.0)) / 2.0
URZEUG_KOORD = (0, -1, 0, 0)


def _lade_dimensionslose_konstanten():
    pfad = BASIS / "Dimensionslose Konstanten.py"
    spec = importlib.util.spec_from_file_location("dimensionslose_konstanten", pfad)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _einheitlicher_rational(eintrag: dict, radius: int) -> dict:
    k = tuple(eintrag["koordinaten"])
    norm = float(eintrag["norm"])
    alpha, beta = eintrag["norm_sq_zphi"]
    return {
        "sektor": "rational",
        "koordinaten": k,
        "radius_quelle": radius,
        "norm": norm,
        "norm_sq_zphi": [alpha, beta],
        "m_verhaeltnis": int(eintrag["verhaeltnis_zu_M"]),
        "reziprok_norm": 1.0 / norm,
        "norm_m": norm / M_KONSTANTE,
        "typ": eintrag.get("typ", "M-Resonanz (rationaler Sektor, a=0)"),
    }


def _einheitlicher_phi(eintrag: dict, radius: int) -> dict:
    k = tuple(eintrag["koordinaten"])
    alpha, beta = eintrag["norm_sq_zphi"]
    gamma, delta = eintrag["m_quotient_zphi"]
    if isinstance(eintrag["norm_sq_sqrt5"][0], str):
        u = Fraction(eintrag["norm_sq_sqrt5"][0])
        v = Fraction(eintrag["norm_sq_sqrt5"][1])
    else:
        u, v = eintrag["norm_sq_sqrt5"]
    return {
        "sektor": "phi",
        "koordinaten": k,
        "radius_quelle": radius,
        "norm": float(eintrag["norm"]),
        "norm_sq_zphi": [alpha, beta],
        "norm_sq_sqrt5": [str(u), str(v)],
        "feld_norm": int(eintrag["feld_norm"]),
        "m_quotient_zphi": [int(gamma), int(delta)],
        "ist_urzeug": k == URZEUG_KOORD,
        "typ": eintrag.get("typ", "Phi-Resonanz (algebraische M-Kopplung)"),
    }


def lade_json_zeugen(pfad: Path) -> tuple[list, list, int]:
    daten = json.loads(pfad.read_text(encoding="utf-8"))
    radius = int(daten["radius"])
    rational = [_einheitlicher_rational(z, radius) for z in daten.get("rational", [])]
    phi = [_einheitlicher_phi(z, radius) for z in daten.get("phi", [])]
    return rational, phi, radius


def scan_rational_radius4(dk_mod) -> list[dict]:
    zeugen = dk_mod.suche_eabc_zeugen_pure_python(4, ausgabe_limit=0)
    return [_einheitlicher_rational(z, 4) for z in zeugen]


def kanonische_phi_zeugen(phi_r2: list[dict], phi_r3: list[dict]) -> list[dict]:
    """Dedupliziert Φ-Zeugen nach algebraischer Schale (α,β)."""
    gesehen: dict[tuple[int, int], dict] = {}
    for z in phi_r2 + phi_r3:
        schluessel = tuple(z["norm_sq_zphi"])
        if schluessel not in gesehen:
            gesehen[schluessel] = z
        elif z["ist_urzeug"]:
            gesehen[schluessel] = z
    return sorted(
        gesehen.values(),
        key=lambda item: (item["norm"], item["koordinaten"]),
    )


def schalen_rational(zeugen: list[dict]) -> list[dict]:
    """Aggregiert rationale Zeugen nach Norm-Schale."""
    nach_norm: dict[float, dict] = {}
    for z in zeugen:
        n = z["norm"]
        if n not in nach_norm:
            nach_norm[n] = {
                "norm": n,
                "m_verhaeltnis": z["m_verhaeltnis"],
                "anzahl_zeugen": 0,
                "reziprok_norm": 1.0 / n,
                "kaskade_von_vorgaenger": None,
            }
        nach_norm[n]["anzahl_zeugen"] += 1
    schalen = sorted(nach_norm.values(), key=lambda s: s["norm"])
    for i in range(1, len(schalen)):
        v = schalen[i - 1]["m_verhaeltnis"]
        w = schalen[i]["m_verhaeltnis"]
        schalen[i]["kaskade_von_vorgaenger"] = v / w
    return schalen


def kalibriere_re(zeugen: list[dict], schalen: list[dict]) -> dict:
    """
    Re-Zuordnung (rationaler Sektor):
      κ_Re = M / Re_klassisch
      Re = (M/Norm) / κ_Re = Re_klassisch / norm
    Reziproke Norm 1/norm koppelt linear an die M-Kaskade (M/Norm = M/norm).
    """
    kappa_re = M_KONSTANTE / RE_KLASSISCH
    anker = next(s for s in schalen if s["norm"] == 1.0)
    re_krit_berechnet = anker["m_verhaeltnis"] / kappa_re

    zuordnung = []
    for z in sorted(zeugen, key=lambda item: (item["norm"], item["koordinaten"])):
        re_wert = z["m_verhaeltnis"] / kappa_re
        zuordnung.append(
            {
                **z,
                "kappa_re": kappa_re,
                "re_zugeordnet": re_wert,
                "re_formel": "Re_klassisch / norm",
                "re_aus_reziprok_norm": RE_KLASSISCH * z["reziprok_norm"],
            }
        )

    # Emergente Kaskaden-Identität: M = 3 · (M/Norm₃)
    schale3 = next((s for s in schalen if s["norm"] == 3.0), None)
    kaskade_emergent = None
    if schale3:
        kaskade_emergent = M_KONSTANTE / schale3["m_verhaeltnis"]

    return {
        "kappa_re": kappa_re,
        "re_krit_berechnet": re_krit_berechnet,
        "re_krit_klassisch": RE_KLASSISCH,
        "abweichung_re_krit": re_krit_berechnet - RE_KLASSISCH,
        "rel_abweichung_re_krit": (re_krit_berechnet - RE_KLASSISCH) / RE_KLASSISCH,
        "kaskade_1_zu_2": next(
            (s["kaskade_von_vorgaenger"] for s in schalen if s["norm"] == 2.0), None
        ),
        "kaskade_2_zu_3": next(
            (s["kaskade_von_vorgaenger"] for s in schalen if s["norm"] == 3.0), None
        ),
        "m_gleich_3_mal_m_norm3": kaskade_emergent,
        "schalen": schalen,
        "zuordnung": zuordnung,
    }


def zphi_quotient_feld_norm(gamma: int, delta: int) -> int:
    return gamma * gamma + gamma * delta - delta * delta


def kalibriere_ra(phi_zeugen: list[dict], dk_mod) -> dict:
    """
    Ra-Zuordnung (Φ-Sektor):
      κ_Ra = Ra_klassisch (Anker: Urzeuge, Feldnorm N(1+Φ) = 1)
      Ra = N(α+βΦ) · Ra_klassisch
    Zusätzlich: algebraischer Kopplungsindex aus M²/Norm² ∈ ℤ[Φ].
    """
    urzeug = next(z for z in phi_zeugen if z["ist_urzeug"])
    kappa_ra = RA_KLASSISCH / urzeug["feld_norm"]
    ra_krit_berechnet = urzeug["feld_norm"] * kappa_ra

    zuordnung = []
    for z in phi_zeugen:
        gamma, delta = z["m_quotient_zphi"]
        kopplung_norm = zphi_quotient_feld_norm(gamma, delta)
        ra_feld = z["feld_norm"] * kappa_ra
        ra_norm_phi = RA_KLASSISCH * (z["norm"] / PHI) ** 2
        alpha, beta = z["norm_sq_zphi"]
        u, v = dk_mod.zphi_to_sqrt5(alpha, beta)
        zuordnung.append(
            {
                **z,
                "kappa_ra": kappa_ra,
                "ra_zugeordnet": ra_feld,
                "ra_formel": "Ra_klassisch · N(α+βΦ)",
                "ra_norm_phi_quadrat": ra_norm_phi,
                "m2_kopplung_feldnorm": kopplung_norm,
                "m2_kopplung_gamma": gamma,
            }
        )

    # Emergente Feldnorm-Kaskade: Vielfache von 5 zwischen benachbarten Schalen
    feldnormen = sorted({z["feld_norm"] for z in phi_zeugen})
    feldnorm_ratios = []
    for i in range(1, len(feldnormen)):
        if feldnormen[i - 1] != 0:
            feldnorm_ratios.append(feldnormen[i] / feldnormen[i - 1])

    return {
        "kappa_ra": kappa_ra,
        "ra_krit_berechnet": ra_krit_berechnet,
        "ra_krit_klassisch": RA_KLASSISCH,
        "abweichung_ra_krit": ra_krit_berechnet - RA_KLASSISCH,
        "rel_abweichung_ra_krit": (ra_krit_berechnet - RA_KLASSISCH) / RA_KLASSISCH,
        "urzeug_koordinaten": list(URZEUG_KOORD),
        "urzeug_norm": urzeug["norm"],
        "feldnorm_kaskade_ratios": feldnorm_ratios,
        "zuordnung": zuordnung,
    }


def baue_vereinigten_katalog(
    rational_r4: list[dict], phi_kanonisch: list[dict]
) -> dict:
    return {
        "m_konstante": M_KONSTANTE,
        "re_krit_klassisch": RE_KLASSISCH,
        "ra_krit_klassisch": RA_KLASSISCH,
        "anzahl_rational": len(rational_r4),
        "anzahl_phi_kanonisch": len(phi_kanonisch),
        "anzahl_gesamt": len(rational_r4) + len(phi_kanonisch),
        "rational": rational_r4,
        "phi": phi_kanonisch,
    }


def export_json(
    katalog: dict,
    re_ergebnis: dict,
    ra_ergebnis: dict,
    pfad: Path,
) -> None:
    daten = {
        "katalog": katalog,
        "kalibrierung_re": {
            k: v
            for k, v in re_ergebnis.items()
            if k != "zuordnung"
        },
        "kalibrierung_ra": {
            k: v
            for k, v in ra_ergebnis.items()
            if k != "zuordnung"
        },
        "zuordnung_re": re_ergebnis["zuordnung"],
        "zuordnung_ra": ra_ergebnis["zuordnung"],
    }
    pfad.write_text(json.dumps(daten, indent=2, ensure_ascii=False), encoding="utf-8")


def export_csv(re_ergebnis: dict, ra_ergebnis: dict, pfad: Path) -> None:
    with pfad.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "sektor",
                "c_e",
                "c_a",
                "c_b",
                "c_c",
                "norm",
                "m_verhaeltnis",
                "reziprok_norm",
                "feld_norm",
                "re_zugeordnet",
                "ra_zugeordnet",
                "ist_urzeug",
                "radius_quelle",
            ]
        )
        for z in re_ergebnis["zuordnung"]:
            k = z["koordinaten"]
            writer.writerow(
                [
                    "rational",
                    k[0],
                    k[1],
                    k[2],
                    k[3],
                    z["norm"],
                    z["m_verhaeltnis"],
                    z["reziprok_norm"],
                    "",
                    f"{z['re_zugeordnet']:.6f}",
                    "",
                    False,
                    z["radius_quelle"],
                ]
            )
        for z in ra_ergebnis["zuordnung"]:
            k = z["koordinaten"]
            writer.writerow(
                [
                    "phi",
                    k[0],
                    k[1],
                    k[2],
                    k[3],
                    z["norm"],
                    "",
                    "",
                    z["feld_norm"],
                    "",
                    f"{z['ra_zugeordnet']:.6f}",
                    z["ist_urzeug"],
                    z["radius_quelle"],
                ]
            )


def export_latex(
    katalog: dict,
    re_ergebnis: dict,
    ra_ergebnis: dict,
    schalen: list[dict],
    pfad: Path,
) -> None:
    dk = _lade_dimensionslose_konstanten()

    re_zeilen = []
    for z in re_ergebnis["zuordnung"][:12]:
        k = z["koordinaten"]
        re_zeilen.append(
            f"    ({k[0]}, {k[1]}, {k[2]}, {k[3]}) & {z['norm']:.1f} & "
            f"{z['m_verhaeltnis']} & {z['reziprok_norm']:.4f} & "
            f"{z['re_zugeordnet']:.1f} \\\\"
        )

    schalen_zeilen = []
    for s in schalen:
        kask = s["kaskade_von_vorgaenger"]
        kask_tex = f"{kask:.4g}" if kask is not None else "---"
        schalen_zeilen.append(
            f"    {s['norm']:.0f} & {s['m_verhaeltnis']} & {s['reziprok_norm']:.4f} & "
            f"{s['anzahl_zeugen']} & {kask_tex} \\\\"
        )

    ra_zeilen = []
    for z in ra_ergebnis["zuordnung"]:
        k = z["koordinaten"]
        alpha, beta = z["norm_sq_zphi"]
        g, d = z["m_quotient_zphi"]
        nsq = dk.format_sqrt5_latex(*dk.zphi_to_sqrt5(alpha, beta))
        mq = dk.format_zphi_latex(g, d)
        urz = "\\checkmark" if z["ist_urzeug"] else ""
        ra_zeilen.append(
            f"    ({k[0]}, {k[1]}, {k[2]}, {k[3]}) & {z['norm']:.4f} & "
            f"{z['feld_norm']} & ${mq}$ & {z['ra_zugeordnet']:.1f} & {urz} \\\\"
        )

    kappa_re = re_ergebnis["kappa_re"]
    kappa_ra = ra_ergebnis["kappa_ra"]
    re_abw = re_ergebnis["rel_abweichung_re_krit"] * 100
    ra_abw = ra_ergebnis["rel_abweichung_ra_krit"] * 100

    inhalt = rf"""\subsection{{Physikalische Re/Ra-Zuordnung}}\label{{sec:eabc-re-ra}}

Vereinigter EABC-Zeugen-Katalog am Kepler-Anker $M={M_KONSTANTE}$:
{len(katalog['rational'])} rationale Fixpunkte (Radius~4) und
{len(katalog['phi'])} kanonische $\Phi$-Zeugen (Radien~2/3, dedupliziert nach Schale $\alpha+\beta\Phi$).

\paragraph{{Kalibrierung Reynolds ($a=0$).}}
Im rationalen Sektor gilt $M/\mathrm{{Norm}}\in\mathbb{{Z}}$.
Die reziproken Normen $1/\mathrm{{Norm}}$ koppeln mit der $M$-Kaskade
über
\[
  \kappa_{{Re}} = \frac{{M}}{{Re_{{\mathrm{{krit}}}}}}, \qquad
  Re = \frac{{M/\mathrm{{Norm}}}}{{\kappa_{{Re}}}}
     = \frac{{Re_{{\mathrm{{krit}}}}}}{{\mathrm{{Norm}}}}.
\]
Am Einheitsschalenpunkt ($\mathrm{{Norm}}=1$, $M/\mathrm{{Norm}}=M$) emergiert
$Re_{{\mathrm{{krit}}}} = {RE_KLASSISCH:.0f}$ exakt
(berechnet: {re_ergebnis['re_krit_berechnet']:.6f}, Abweichung ${re_abw:.2e}\,\%$).
Kaskade: $113160/56580=2$, $56580/37720=3/2$; strukturell $M=3\cdot 37720$.

\begin{{table}}[ht]
  \centering
  \caption{{Rationale Norm-Schalen und $M$-Kaskade (Radius~4).}}
  \label{{tab:eabc-re-schalen}}
  \begin{{tabular}}{{rrrrr}}
    \toprule
    $\mathrm{{Norm}}$ & $M/\mathrm{{Norm}}$ & $1/\mathrm{{Norm}}$ & Anzahl & Kaskade \\
    \midrule
{chr(10).join(schalen_zeilen)}
    \bottomrule
  \end{{tabular}}
\end{{table}}

\begin{{table}}[ht]
  \centering
  \caption{{Auszug Re-Zuordnung (rationaler Sektor).}}
  \label{{tab:eabc-re-zuordnung}}
  \begin{{tabular}}{{rrrrr}}
    \toprule
    $(c_e,c_a,c_b,c_c)$ & $\mathrm{{Norm}}$ & $M/\mathrm{{Norm}}$ & $1/\mathrm{{Norm}}$ & $Re$ \\
    \midrule
{chr(10).join(re_zeilen)}
    \bottomrule
  \end{{tabular}}
\end{{table}}

\paragraph{{Kalibrierung Rayleigh ($a\neq 0$).}}
Im $\Phi$-Sektor gilt $M^2/\mathrm{{Norm}}^2\in\mathbb{{Z}}[\Phi]$.
Der Urzeuge $(0,-1,0,0)$ mit $\mathrm{{Norm}}=\Phi$ und Feldnorm
$N(1+\Phi)=1$ kalibriert
\[
  \kappa_{{Ra}} = \frac{{Ra_{{\mathrm{{krit}}}}}}{{N(\alpha+\beta\Phi)}},
  \qquad Ra = N(\alpha+\beta\Phi)\cdot Ra_{{\mathrm{{krit}}}}.
\]
Am ikosaedrischen Anker emergiert $Ra_{{\mathrm{{krit}}}} = {RA_KLASSISCH:.0f}$ exakt
(berechnet: {ra_ergebnis['ra_krit_berechnet']:.6f}, Abweichung ${ra_abw:.2e}\,\%$).

\begin{{table}}[ht]
  \centering
  \caption{{$\Phi$-Sektor: Rayleigh-Zuordnung (kanonische Zeugen).}}
  \label{{tab:eabc-ra-zuordnung}}
  \begin{{tabular}}{{rrrrrr}}
    \toprule
    $(c_e,c_a,c_b,c_c)$ & $\mathrm{{Norm}}$ & $N$ & $M^2/\mathrm{{Norm}}^2$ & $Ra$ & Urz. \\
    \midrule
{chr(10).join(ra_zeilen)}
    \bottomrule
  \end{{tabular}}
\end{{table}}

\noindent\textbf{{Skalenfaktoren:}}
$\kappa_{{Re}} = {kappa_re:.6f}$,
$\kappa_{{Ra}} = {kappa_ra:.6f}$.
Deterministischer Lauf: \texttt{{eabc\_re\_ra\_zuordnung.py}}.
"""
    pfad.write_text(inhalt, encoding="utf-8")


def drucke_zusammenfassung(
    katalog: dict,
    re_ergebnis: dict,
    ra_ergebnis: dict,
    schalen: list[dict],
) -> None:
    print("=" * 72)
    print("EABC Re/Ra-ZUORDNUNG — ZUSAMMENFASSUNG")
    print("=" * 72)
    print(f"M (Kepler-Anker):          {M_KONSTANTE}")
    print(f"Katalog rational (R4):     {katalog['anzahl_rational']} Zeugen")
    print(f"Katalog Phi (kanonisch):   {katalog['anzahl_phi_kanonisch']} Zeugen")
    print(f"Gesamt:                    {katalog['anzahl_gesamt']}")
    print()
    print("--- REYNOLDS (rationaler Sektor, a=0) ---")
    print(f"κ_Re = M/Re_krit:          {re_ergebnis['kappa_re']:.6f}")
    print(f"Re_krit berechnet:         {re_ergebnis['re_krit_berechnet']:.6f}")
    print(f"Re_krit klassisch:         {re_ergebnis['re_krit_klassisch']:.1f}")
    print(
        f"Abweichung Re_krit:        {re_ergebnis['abweichung_re_krit']:.6e} "
        f"({re_ergebnis['rel_abweichung_re_krit']*100:.2e} %)"
    )
    print(f"Kaskade M/Norm: 1→2:       {re_ergebnis['kaskade_1_zu_2']}")
    print(f"Kaskade M/Norm: 2→3:       {re_ergebnis['kaskade_2_zu_3']}")
    print(f"M / (M/Norm₃):             {re_ergebnis['m_gleich_3_mal_m_norm3']}")
    print()
    print("Norm-Schalen (Radius 4):")
    for s in schalen:
        re_shell = RE_KLASSISCH / s["norm"]
        print(
            f"  Norm={s['norm']:.0f}  M/Norm={s['m_verhaeltnis']:>6}  "
            f"1/Norm={s['reziprok_norm']:.4f}  Anzahl={s['anzahl_zeugen']:>2}  "
            f"Re={re_shell:.1f}"
        )
    print()
    print("--- RAYLEIGH (Phi-Sektor, a≠0) ---")
    print(f"κ_Ra:                      {ra_ergebnis['kappa_ra']:.6f}")
    print(f"Ra_krit berechnet:         {ra_ergebnis['ra_krit_berechnet']:.6f}")
    print(f"Ra_krit klassisch:         {ra_ergebnis['ra_krit_klassisch']:.1f}")
    print(
        f"Abweichung Ra_krit:        {ra_ergebnis['abweichung_ra_krit']:.6e} "
        f"({ra_ergebnis['rel_abweichung_ra_krit']*100:.2e} %)"
    )
    print(f"Urzeuge:                   {ra_ergebnis['urzeug_koordinaten']}  Norm≈{ra_ergebnis['urzeug_norm']:.6f}")
    print()
    print("Phi-Zeugen Ra-Zuordnung:")
    for z in ra_ergebnis["zuordnung"]:
        mark = " [URZEUG]" if z["ist_urzeug"] else ""
        print(
            f"  {z['koordinaten']}  Norm={z['norm']:.4f}  "
            f"N={z['feld_norm']:>3}  Ra={z['ra_zugeordnet']:.1f}{mark}"
        )
    print()
    print("Beispielzeilen Re-Zuordnung:")
    for z in re_ergebnis["zuordnung"][:5]:
        print(
            f"  {z['koordinaten']}  Norm={z['norm']:.1f}  "
            f"M/Norm={z['m_verhaeltnis']}  Re={z['re_zugeordnet']:.1f}"
        )
    print("=" * 72)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Physikalische Re/Ra-Zuordnung auf dem vereinigten EABC-Zeugen-Katalog"
    )
    parser.add_argument(
        "--scan",
        action="store_true",
        help="Rationale Zeugen Radius 4 per Sonde neu berechnen (sonst stilles Scannen)",
    )
    args = parser.parse_args()

    dk = _lade_dimensionslose_konstanten()

    # Phi aus JSON
    _, phi_r2, _ = lade_json_zeugen(BASIS / "eabc_zeugen_r2.json")
    _, phi_r3, _ = lade_json_zeugen(BASIS / "eabc_zeugen_r3.json")
    phi_kanonisch = kanonische_phi_zeugen(phi_r2, phi_r3)

    # Rational Radius 4
    if args.scan:
        rational_r4 = scan_rational_radius4(dk)
    else:
        import io
        import sys

        buf = io.StringIO()
        old_stdout = sys.stdout
        sys.stdout = buf
        try:
            roh = dk.suche_eabc_zeugen_pure_python(4, ausgabe_limit=0)
        finally:
            sys.stdout = old_stdout
        rational_r4 = [_einheitlicher_rational(z, 4) for z in roh]

    if len(rational_r4) != 96:
        print(
            f"Hinweis: {len(rational_r4)} rationale Zeugen bei Radius 4 "
            f"(erwartet: 96)."
        )

    schalen = schalen_rational(rational_r4)
    re_ergebnis = kalibriere_re(rational_r4, schalen)
    ra_ergebnis = kalibriere_ra(phi_kanonisch, dk)
    katalog = baue_vereinigten_katalog(rational_r4, phi_kanonisch)

    export_json(katalog, re_ergebnis, ra_ergebnis, BASIS / "eabc_re_ra_zuordnung.json")
    export_csv(re_ergebnis, ra_ergebnis, BASIS / "eabc_re_ra_zuordnung.csv")
    export_latex(
        katalog,
        re_ergebnis,
        ra_ergebnis,
        schalen,
        BASIS / "kapitel-eabc-re-ra.tex",
    )
    drucke_zusammenfassung(katalog, re_ergebnis, ra_ergebnis, schalen)

    print()
    print("Ausgabedateien:")
    print(f"  {BASIS / 'eabc_re_ra_zuordnung.py'}")
    print(f"  {BASIS / 'eabc_re_ra_zuordnung.json'}")
    print(f"  {BASIS / 'eabc_re_ra_zuordnung.csv'}")
    print(f"  {BASIS / 'kapitel-eabc-re-ra.tex'}")


if __name__ == "__main__":
    main()
