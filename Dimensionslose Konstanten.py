import argparse
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

M_KONSTANTE = 113160


def norm_sq_zphi(c_e, c_a, c_b, c_c):
    """Norm² = c_e² + (c_a·Φ)² + c_b² + c_c² als (α, β) mit α + β·Φ."""
    alpha = c_e**2 + c_b**2 + c_c**2 + c_a**2
    beta = c_a**2
    return alpha, beta


def zphi_to_sqrt5(alpha, beta):
    """Isomorphie ℤ[Φ] ≅ ℚ(√5): α + β·Φ = u + v·√5 mit u = (2α+β)/2, v = β/2."""
    u = Fraction(2 * alpha + beta, 2)
    v = Fraction(beta, 2)
    return u, v


def zphi_mult(a1, b1, a2, b2):
    """(a1 + b1·Φ)(a2 + b2·Φ) unter Φ² = Φ + 1."""
    return a1 * a2 + b1 * b2, a1 * b2 + a2 * b1 + b1 * b2


def zphi_field_norm(alpha, beta):
    """N(α + β·Φ) = α² + αβ − β² ∈ ℤ."""
    return alpha * alpha + alpha * beta - beta * beta


def zphi_divides_zphi(num_a, num_b, den_a, den_b):
    """Prüft, ob (den_a + den_b·Φ) | (num_a + num_b·Φ) in ℤ[Φ]."""
    if den_a == 0 and den_b == 0:
        return False
    norm = zphi_field_norm(den_a, den_b)
    if norm == 0:
        return False
    conj_a, conj_b = den_a + den_b, -den_b
    prod_a, prod_b = zphi_mult(num_a, num_b, conj_a, conj_b)
    return prod_a % norm == 0 and prod_b % norm == 0


def zphi_quotient(num_a, num_b, den_a, den_b):
    """Gibt (γ, δ) ∈ ℤ² mit (num)/(den) = γ + δ·Φ zurück, falls exakt teilbar."""
    norm = zphi_field_norm(den_a, den_b)
    if norm == 0:
        return None
    conj_a, conj_b = den_a + den_b, -den_b
    prod_a, prod_b = zphi_mult(num_a, num_b, conj_a, conj_b)
    if prod_a % norm != 0 or prod_b % norm != 0:
        return None
    return prod_a // norm, prod_b // norm


def norm_float(alpha, beta):
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    return math.sqrt(alpha + beta * phi)


def format_sqrt5(u, v):
    """LaTeX-/Text-Darstellung von u + v·√5."""
    if v == 0:
        return str(u)
    sign = "+" if v > 0 else "-"
    av = abs(v)
    if u == 0:
        if av == 1:
            return f"{sign.strip()}√5" if sign == "+" else "-√5"
        return f"{sign}{av}√5" if sign == "+" else f"-{av}√5"
    if av == 1:
        return f"{u} {sign} √5"
    return f"{u} {sign} {av}√5"


def format_zphi(alpha, beta):
    if beta == 0:
        return str(alpha)
    if alpha == 0:
        return f"{beta}Φ" if beta != 1 else "Φ"
    sign = "+" if beta > 0 else "-"
    ab = abs(beta)
    return f"{alpha} {sign} {ab}Φ" if ab != 1 else f"{alpha} {sign} Φ"


def format_sqrt5_latex(u, v):
    if v == 0:
        return str(u)
    sign = "+" if v > 0 else "-"
    av = abs(v)
    if u == 0:
        coef = "" if av == 1 else str(av)
        return f"{coef}\\sqrt{{5}}" if sign == "+" else f"-{coef}\\sqrt{{5}}"
    coef = "" if av == 1 else str(av)
    return f"{u} {sign} {coef}\\sqrt{{5}}"


def format_zphi_latex(alpha, beta):
    if beta == 0:
        return str(alpha)
    if alpha == 0:
        return f"{beta}\\Phi" if beta != 1 else "\\Phi"
    sign = "+" if beta > 0 else "-"
    ab = abs(beta)
    return f"{alpha} {sign} {ab}\\Phi" if ab != 1 else f"{alpha} {sign} \\Phi"


def suche_eabc_zeugen_pure_python(maximaler_radius, ausgabe_limit=15):
    print("=== STARTE REINES PYTHON-SKRIPT: NUMERISCHE ZEUGENSUCHE ===")
    print(f"Scanne EABC-Gitterkoeffizienten bis Radius {maximaler_radius}...")

    gefundene_zeugen = []

    for c_e in range(-maximaler_radius, maximaler_radius + 1):
        for c_a in range(-maximaler_radius, maximaler_radius + 1):
            for c_b in range(-maximaler_radius, maximaler_radius + 1):
                for c_c in range(-maximaler_radius, maximaler_radius + 1):
                    if c_e == 0 and c_a == 0 and c_b == 0 and c_c == 0:
                        continue

                    alpha, beta = norm_sq_zphi(c_e, c_a, c_b, c_c)
                    norm = norm_float(alpha, beta)

                    if norm <= 1e-9:
                        continue

                    verhaeltnis = M_KONSTANTE / norm
                    rest = verhaeltnis - round(verhaeltnis)

                    if abs(rest) < 1e-7:
                        gefundene_zeugen.append(
                            {
                                "koordinaten": (c_e, c_a, c_b, c_c),
                                "norm": round(norm, 6),
                                "norm_sq_zphi": (alpha, beta),
                                "verhaeltnis_zu_M": int(round(verhaeltnis)),
                                "typ": "M-Resonanz (rationaler Sektor, a=0)",
                            }
                        )

    print(f"\nSuche abgeschlossen. {len(gefundene_zeugen)} rationale Zeugen isoliert.")
    print("-" * 70)
    print(f"{'EABC-Koeffizienten':<22} | {'Norm':<10} | {'M / Norm':<10} | {'Typ'}")
    print("-" * 70)

    for z in gefundene_zeugen[:ausgabe_limit]:
        koord_str = f"({z['koordinaten'][0]}, {z['koordinaten'][1]}, {z['koordinaten'][2]}, {z['koordinaten'][3]})"
        print(
            f"{koord_str:<22} | {z['norm']:<10} | {z['verhaeltnis_zu_M']:<10} | {z['typ']}"
        )
    print("-" * 70)

    normen = Counter(z["norm"] for z in gefundene_zeugen)
    mit_a_achse = sum(1 for z in gefundene_zeugen if z["koordinaten"][1] != 0)
    verhaeltnisse = Counter(z["verhaeltnis_zu_M"] for z in gefundene_zeugen)

    print("Zusammenfassung:")
    print(f"  Normenverteilung: {dict(sorted(normen.items()))}")
    print(f"  M/Norm-Teiler:    {dict(sorted(verhaeltnisse.items()))}")
    print(f"  Zeugen mit a≠0:   {mit_a_achse} / {len(gefundene_zeugen)}")
    print()

    return gefundene_zeugen


def suche_eabc_phizeugen(
    maximaler_radius,
    ausgabe_limit=20,
    m_konstante=M_KONSTANTE,
    max_quotient_norm=None,
):
    """
    Φ-Sektor-Zeugen: c_a ≠ 0, Norm² ∈ ℤ[Φ] ≅ u + v·√5,
    algebraische M-Kopplung M²/(α+βΦ) ∈ ℤ[Φ].
    """
    print("=== ERWEITERTE ZEUGENSUCHE: PHI-RESONANZEN ===")
    print(f"Scanne Φ-Sektor bis Radius {maximaler_radius} (M = {m_konstante})...")

    m_sq = m_konstante * m_konstante
    gefundene_zeugen = []
    gesehen = set()

    for c_e in range(-maximaler_radius, maximaler_radius + 1):
        for c_a in range(-maximaler_radius, maximaler_radius + 1):
            for c_b in range(-maximaler_radius, maximaler_radius + 1):
                for c_c in range(-maximaler_radius, maximaler_radius + 1):
                    if c_e == 0 and c_a == 0 and c_b == 0 and c_c == 0:
                        continue
                    if c_a == 0:
                        continue

                    alpha, beta = norm_sq_zphi(c_e, c_a, c_b, c_c)
                    u, v = zphi_to_sqrt5(alpha, beta)
                    norm = norm_float(alpha, beta)

                    if norm <= 1e-9:
                        continue

                    # Algebraische M-Resonanz: M²/(α+βΦ) ∈ ℤ[Φ]
                    if not zphi_divides_zphi(m_sq, 0, alpha, beta):
                        continue

                    gamma, delta = zphi_quotient(m_sq, 0, alpha, beta)
                    q_norm = max(abs(gamma), abs(delta))

                    if max_quotient_norm is not None and q_norm > max_quotient_norm:
                        continue

                    # Duplikat-Reduktion: gleiche algebraische Schale
                    schluessel = (alpha, beta, gamma, delta)
                    if schluessel in gesehen:
                        continue
                    gesehen.add(schluessel)

                    feld_norm = zphi_field_norm(alpha, beta)

                    gefundene_zeugen.append(
                        {
                            "koordinaten": (c_e, c_a, c_b, c_c),
                            "norm": round(norm, 6),
                            "norm_sq_zphi": (alpha, beta),
                            "norm_sq_sqrt5": (u, v),
                            "feld_norm": feld_norm,
                            "m_quotient_zphi": (gamma, delta),
                            "typ": "Phi-Resonanz (algebraische M-Kopplung)",
                        }
                    )

    gefundene_zeugen.sort(
        key=lambda z: (
            z["norm"],
            z["koordinaten"],
        )
    )

    print(f"\nSuche abgeschlossen. {len(gefundene_zeugen)} Φ-Zeugen isoliert.")
    print("-" * 95)
    print(
        f"{'EABC-Koeffizienten':<22} | {'Norm':<10} | {'Norm² (√5)':<18} | {'M²/Norm² (ℤ[Φ])':<18} | Typ"
    )
    print("-" * 95)

    for z in gefundene_zeugen[:ausgabe_limit]:
        k = z["koordinaten"]
        u, v = z["norm_sq_sqrt5"]
        g, d = z["m_quotient_zphi"]
        koord_str = f"({k[0]}, {k[1]}, {k[2]}, {k[3]})"
        nsq = format_sqrt5(u, v)
        mq = format_zphi(g, d)
        print(
            f"{koord_str:<22} | {z['norm']:<10} | {nsq:<18} | {mq:<18} | {z['typ']}"
        )
    print("-" * 95)

    normen = Counter(z["norm"] for z in gefundene_zeugen)
    schalen = Counter(z["norm_sq_zphi"] for z in gefundene_zeugen)
    print("Zusammenfassung:")
    print(f"  Normenverteilung: {dict(sorted(normen.items()))}")
    print(f"  Algebraische Schalen (α,β): {len(schalen)} verschiedene")
    for schale, count in sorted(
        schalen.items(), key=lambda item: norm_float(item[0][0], item[0][1])
    ):
        a, b = schale
        u, v = zphi_to_sqrt5(a, b)
        print(
            f"    {format_zphi(a, b)} = {format_sqrt5(u, v)}  →  {count} Zeugen, "
            f"N={zphi_field_norm(a, b)}"
        )
    print()

    return gefundene_zeugen


def _json_ready_phi(zeugen):
    out = []
    for z in zeugen:
        u, v = z["norm_sq_sqrt5"]
        out.append(
            {
                **z,
                "norm_sq_sqrt5": [str(u), str(v)],
                "norm_sq_zphi": list(z["norm_sq_zphi"]),
                "m_quotient_zphi": list(z["m_quotient_zphi"]),
            }
        )
    return out


def export_zeugen_json(rational, phi, radius, pfad):
    daten = {
        "m_konstante": M_KONSTANTE,
        "radius": radius,
        "rational": rational,
        "phi": _json_ready_phi(phi),
    }
    pfad.write_text(json.dumps(daten, indent=2, ensure_ascii=False), encoding="utf-8")
    return pfad


def export_zeugen_latex(rational, phi, radius, pfad):
    """\\input{...}-fähige LaTeX-Tabelle für rationale und Φ-Zeugen."""
    zeilen_rational = []
    for z in sorted(rational, key=lambda x: (x["norm"], x["koordinaten"]))[:24]:
        k = z["koordinaten"]
        zeilen_rational.append(
            f"    ({k[0]}, {k[1]}, {k[2]}, {k[3]}) & {z['norm']:.1f} & "
            f"{z['verhaeltnis_zu_M']} & rational \\\\"
        )

    zeilen_phi = []
    for z in phi[:24]:
        k = z["koordinaten"]
        a, b = z["norm_sq_zphi"]
        g, d = z["m_quotient_zphi"]
        u, v = z["norm_sq_sqrt5"]
        nsq_tex = format_sqrt5_latex(u, v)
        mq_tex = format_zphi_latex(g, d)
        zeilen_phi.append(
            f"    ({k[0]}, {k[1]}, {k[2]}, {k[3]}) & {z['norm']:.6f} & "
            f"${nsq_tex}$ & ${mq_tex}$ \\\\"
        )

    inhalt = rf"""\subsection{{EABC-Zeugensonde, Radius {radius}}}\label{{sec:eabc-zeugen-r{radius}}}

Deterministischer Lauf von \texttt{{Dimensionslose Konstanten.py}} am Fixpunkt $M={M_KONSTANTE}$.

\paragraph{{Rationaler Sektor ($a=0$, $M/\mathrm{{Norm}}\in\mathbb{{Z}}$).}}
{len(rational)} Zeugen; Filter \texttt{{abs(M/Norm - round(M/Norm)) < 1e-7}}.

\begin{{table}}[ht]
  \centering
  \caption{{Rationale EABC-Zeugen, Radius {radius}.}}
  \label{{tab:eabc-rational-r{radius}}}
  \begin{{tabular}}{{rrrr}}
    \toprule
    $(c_e,c_a,c_b,c_c)$ & $\mathrm{{Norm}}$ & $M/\mathrm{{Norm}}$ & Sektor \\
    \midrule
{chr(10).join(zeilen_rational)}
    \bottomrule
  \end{{tabular}}
\end{{table}}

\paragraph{{$\Phi$-Sektor ($a\neq 0$, $M^2/\mathrm{{Norm}}^2\in\mathbb{{Z}}[\Phi]$).}}
{len(phi)} Zeugen; algebraische Kopplung statt ganzzahligem $M/\mathrm{{Norm}}$.

\begin{{table}}[ht]
  \centering
  \caption{{$\Phi$-Sektor-Zeugen, Radius {radius}.}}
  \label{{tab:eabc-phi-r{radius}}}
  \begin{{tabular}}{{rrrr}}
    \toprule
    $(c_e,c_a,c_b,c_c)$ & $\mathrm{{Norm}}$ & $\mathrm{{Norm}}^2$ & $M^2/\mathrm{{Norm}}^2$ \\
    \midrule
{chr(10).join(zeilen_phi)}
    \bottomrule
  \end{{tabular}}
\end{{table}}
"""
    pfad.write_text(inhalt, encoding="utf-8")
    return pfad


def run_tiefenscan(radien, modus="beide", ausgabe_limit=15):
    ergebnisse = {}
    for radius in radien:
        print(f"\n{'=' * 70}")
        print(f"TIEFENSCAN: maximaler_radius = {radius}")
        print(f"{'=' * 70}")
        rational = []
        phi = []
        if modus in ("rational", "beide"):
            rational = suche_eabc_zeugen_pure_python(radius, ausgabe_limit=ausgabe_limit)
        if modus in ("phi", "beide"):
            phi = suche_eabc_phizeugen(radius, ausgabe_limit=ausgabe_limit)
        ergebnisse[radius] = {"rational": rational, "phi": phi}

        basis = Path(__file__).resolve().parent
        export_zeugen_json(
            rational,
            phi,
            radius,
            basis / f"eabc_zeugen_r{radius}.json",
        )
        export_zeugen_latex(
            rational,
            phi,
            radius,
            basis / f"kapitel-eabc-zeugen-r{radius}.tex",
        )
    return ergebnisse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="EABC-Zeugensonde im 120-Zeller-Raster")
    parser.add_argument(
        "--radius",
        type=int,
        nargs="+",
        default=[2, 3, 4],
        help="Suchradien (Standard: 2 3 4)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=15,
        help="Anzahl ausgegebener Zeugen pro Scan",
    )
    parser.add_argument(
        "--modus",
        choices=["rational", "phi", "beide"],
        default="beide",
        help="Rationaler Sektor, Phi-Sektor oder beide",
    )
    args = parser.parse_args()

    run_tiefenscan(args.radius, modus=args.modus, ausgabe_limit=args.limit)
