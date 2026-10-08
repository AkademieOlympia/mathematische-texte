#!/usr/bin/env python3
"""Build Primes_EABC_kennzahlen.xlsx from out_eabc/results.json."""

import json
from pathlib import Path

from openpyxl import Workbook
from openpyxl.chart import BarChart, LineChart, Reference
from openpyxl.chart.label import DataLabelList
from openpyxl.chart.series import SeriesLabel
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table, TableStyleInfo
from openpyxl.formatting.rule import CellIsRule, ColorScaleRule
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.drawing.line import LineProperties
from openpyxl.chart.layout import Layout, ManualLayout

ROOT = Path(__file__).resolve().parent
JSON_PATH = ROOT / "results.json"
OUT_PATH = ROOT / "Primes_EABC_kennzahlen.xlsx"

NAVY = "1E3A5F"
NAVY2 = "2C5282"
SURFACE = "F7FAFC"
WHITE = "FFFFFF"
INK = "1A202C"
MUTED = "4A5568"
LINE = "CBD5E0"
WARN = "C05621"
OK = "276749"
GMP = "2B6CB0"
EABC = "C53030"

thin = Border(
    left=Side(style="thin", color=LINE),
    right=Side(style="thin", color=LINE),
    top=Side(style="thin", color=LINE),
    bottom=Side(style="thin", color=LINE),
)
title_font = Font(name="Calibri", size=18, bold=True, color=WHITE)
section_font = Font(name="Calibri", size=12, bold=True, color=WHITE)
header_font = Font(name="Calibri", size=10, bold=True, color=WHITE)
body_font = Font(name="Calibri", size=10, color=INK)
muted_font = Font(name="Calibri", size=9, italic=True, color=MUTED)
kpi_font = Font(name="Calibri", size=16, bold=True, color=NAVY)
label_font = Font(name="Calibri", size=9, color=MUTED)


def fill(hex_color):
    return PatternFill("solid", fgColor=hex_color)


def apply_widths(ws, widths):
    for col, w in widths.items():
        ws.column_dimensions[col].width = w


def style_header_row(ws, row, cols, fill_color=NAVY2):
    for col in range(1, cols + 1):
        cell = ws.cell(row, col)
        cell.font = header_font
        cell.fill = fill(fill_color)
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = thin


def style_body(ws, r1, r2, c1, c2):
    for r in range(r1, r2 + 1):
        bg = WHITE if (r - r1) % 2 == 0 else SURFACE
        for c in range(c1, c2 + 1):
            cell = ws.cell(r, c)
            cell.font = body_font
            cell.fill = fill(bg)
            cell.border = thin
            cell.alignment = Alignment(
                horizontal="right" if c > 1 else "center",
                vertical="center",
            )


def load_runs():
    data = json.loads(JSON_PATH.read_text())
    rows = []
    for run in data["runs"]:
        e = run["engines"]
        rows.append(
            {
                "bits": run["bits"],
                "trial": run["trial"],
                "gmp_ms": e["gmp_stock"]["t_10_ns"] / 1e6,
                "raw_ms": e["openssl_stock"]["t_10_ns"] / 1e6,
                "sieved_ms": e["openssl_sieved"]["t_10_ns"] / 1e6,
                "eabc_ms": e["eabc_mr"]["t_10_ns"] / 1e6,
                "odd": e["openssl_stock"]["candidates_tested"],
                "sieve": e["eabc_mr"]["candidates_tested"],
                "mr": e["eabc_mr"]["mr_rounds_executed"],
                "agree": run["engines_agree"],
                "order": ",".join(run["engine_order"]),
            }
        )
    rows.sort(key=lambda x: (x["bits"], x["trial"]))
    return data, rows


def main():
    meta, rows = load_runs()
    wb = Workbook()

    ws_a = wb.active
    ws_a.title = "Aussage"
    ws_k = wb.create_sheet("Kennzahlen")
    ws_r = wb.create_sheet("Rohdaten")

    # --- Rohdaten first (source of truth) ---
    ws_r.sheet_view.showGridLines = False
    ws_r.merge_cells("A1:K1")
    ws_r["A1"] = "Rohdaten · 15 Läufe · T10 in Millisekunden"
    ws_r["A1"].font = title_font
    ws_r["A1"].fill = fill(NAVY)
    ws_r["A1"].alignment = Alignment(horizontal="left", vertical="center")
    ws_r.row_dimensions[1].height = 28
    ws_r.merge_cells("A2:K2")
    ws_r["A2"] = (
        f"{meta['hardware']['cpu']} · {meta['libraries']['gmp']} · "
        f"{meta['libraries']['openssl']} · isolate={meta['isolate']} · "
        f"failures={meta['failures']} · Seed {meta['master_seed']}"
    )
    ws_r["A2"].font = muted_font
    ws_r["A2"].fill = fill(SURFACE)

    headers = [
        "Bits",
        "Trial",
        "GMP_ms",
        "OpenSSL_raw_ms",
        "OpenSSL_sieved_ms",
        "EABC_ms",
        "Odd_Walk",
        "Sieb",
        "MR_Runden",
        "Agree",
        "Engine_Reihenfolge",
    ]
    for i, h in enumerate(headers, 1):
        ws_r.cell(4, i, h)
    style_header_row(ws_r, 4, 11)

    for i, row in enumerate(rows):
        r = 5 + i
        ws_r.cell(r, 1, row["bits"])
        ws_r.cell(r, 2, row["trial"])
        ws_r.cell(r, 3, round(row["gmp_ms"], 6))
        ws_r.cell(r, 4, round(row["raw_ms"], 6))
        ws_r.cell(r, 5, round(row["sieved_ms"], 6))
        ws_r.cell(r, 6, round(row["eabc_ms"], 6))
        ws_r.cell(r, 7, row["odd"])
        ws_r.cell(r, 8, row["sieve"])
        ws_r.cell(r, 9, row["mr"])
        ws_r.cell(r, 10, row["agree"])
        ws_r.cell(r, 11, row["order"])
    style_body(ws_r, 5, 19, 1, 11)
    for r in range(5, 20):
        for c in (3, 4, 5, 6):
            ws_r.cell(r, c).number_format = "0.000"
        ws_r.cell(r, 11).alignment = Alignment(horizontal="left", vertical="center")
        ws_r.cell(r, 11).font = Font(name="Calibri", size=8, color=INK)

    tab = Table(displayName="RohdatenTab", ref="A4:K19")
    tab.tableStyleInfo = TableStyleInfo(name="TableStyleMedium2", showRowStripes=True)
    ws_r.tables.add(tab)
    apply_widths(ws_r, {
        "A": 10, "B": 10, "C": 12, "D": 18, "E": 20, "F": 12,
        "G": 12, "H": 10, "I": 12, "J": 10, "K": 52,
    })
    ws_r.freeze_panes = "A5"
    ws_r.auto_filter.ref = "A4:K19"

    # blocks: 512=5:9, 1024=10:14, 2048=15:19
    blocks = [(512, 5, 9), (1024, 10, 14), (2048, 15, 19)]

    # --- Kennzahlen (formula-driven) ---
    ws_k.sheet_view.showGridLines = False
    ws_k.merge_cells("A1:H1")
    ws_k["A1"] = "Kennzahlen aus den Medianen der 15 Läufe"
    ws_k["A1"].font = title_font
    ws_k["A1"].fill = fill(NAVY)
    ws_k["A1"].alignment = Alignment(horizontal="left", vertical="center")
    ws_k.row_dimensions[1].height = 28
    ws_k.merge_cells("A2:H2")
    ws_k["A2"] = (
        "Alle abgeleiteten Größen sind Excel-Formeln auf Rohdaten. "
        "µs/Runde = Median(T10_EABC) / Median(MR-Runden), also Verhältnis der Mediane."
    )
    ws_k["A2"].font = muted_font
    ws_k["A2"].fill = fill(SURFACE)

    # Table 1 medians
    ws_k["A4"] = "1 · Mediane T10 (ms)"
    ws_k["A4"].font = section_font
    ws_k["A4"].fill = fill(NAVY2)
    ws_k.merge_cells("A4:E4")
    med_headers = ["Bits", "GMP", "OpenSSL raw", "OpenSSL sieved", "EABC"]
    for i, h in enumerate(med_headers, 1):
        ws_k.cell(5, i, h)
    style_header_row(ws_k, 5, 5)
    for i, (bits, r1, r2) in enumerate(blocks):
        r = 6 + i
        ws_k.cell(r, 1, bits)
        ws_k.cell(r, 2).value = f"=MEDIAN(Rohdaten!C{r1}:C{r2})"
        ws_k.cell(r, 3).value = f"=MEDIAN(Rohdaten!D{r1}:D{r2})"
        ws_k.cell(r, 4).value = f"=MEDIAN(Rohdaten!E{r1}:E{r2})"
        ws_k.cell(r, 5).value = f"=MEDIAN(Rohdaten!F{r1}:F{r2})"
        for c in range(2, 6):
            ws_k.cell(r, c).number_format = "0.000"
    style_body(ws_k, 6, 8, 1, 5)

    # Table 2 decomposition
    ws_k["A10"] = "2 · Zerlegung des Abstands"
    ws_k["A10"].font = section_font
    ws_k["A10"].fill = fill(NAVY2)
    ws_k.merge_cells("A10:F10")
    z_headers = [
        "Bits",
        "Siebgewinn raw/sieved",
        "Testgewinn sieved/EABC",
        "GMP/EABC",
        "GMP-Faktor EABC/GMP",
        "EABC hinter GMP",
    ]
    for i, h in enumerate(z_headers, 1):
        ws_k.cell(11, i, h)
    style_header_row(ws_k, 11, 6)
    for i in range(3):
        r = 12 + i
        src = 6 + i
        ws_k.cell(r, 1).value = f"=A{src}"
        ws_k.cell(r, 2).value = f"=C{src}/D{src}"
        ws_k.cell(r, 3).value = f"=D{src}/E{src}"
        ws_k.cell(r, 4).value = f"=B{src}/E{src}"
        ws_k.cell(r, 5).value = f"=E{src}/B{src}"
        ws_k.cell(r, 6).value = f"=1-D{r}"
        for c in range(2, 7):
            ws_k.cell(r, c).number_format = "0.000"
        ws_k.cell(r, 6).number_format = "0.0%"
    style_body(ws_k, 12, 14, 1, 6)
    ws_k.conditional_formatting.add(
        "D12:D14",
        CellIsRule(operator="lessThan", formula=["1"], fill=fill("FED7D7"), font=Font(color="9B2C2C", bold=True)),
    )

    # Table 3 candidates
    ws_k["A16"] = "3 · Kandidaten (Median)"
    ws_k["A16"].font = section_font
    ws_k["A16"].fill = fill(NAVY2)
    ws_k.merge_cells("A16:D16")
    for i, h in enumerate(["Bits", "Odd-Walk", "Sieb", "Verhältnis odd/Sieb"], 1):
        ws_k.cell(17, i, h)
    style_header_row(ws_k, 17, 4)
    for i, (bits, r1, r2) in enumerate(blocks):
        r = 18 + i
        ws_k.cell(r, 1, bits)
        ws_k.cell(r, 2).value = f"=MEDIAN(Rohdaten!G{r1}:G{r2})"
        ws_k.cell(r, 3).value = f"=MEDIAN(Rohdaten!H{r1}:H{r2})"
        ws_k.cell(r, 4).value = f"=B{r}/C{r}"
        ws_k.cell(r, 2).number_format = "0"
        ws_k.cell(r, 3).number_format = "0"
        ws_k.cell(r, 4).number_format = "0.00"
    style_body(ws_k, 18, 20, 1, 4)

    # Table 4 MR
    ws_k["A22"] = "4 · MR-Runden und Zeit pro Runde"
    ws_k["A22"].font = section_font
    ws_k["A22"].fill = fill(NAVY2)
    ws_k.merge_cells("A22:G22")
    for i, h in enumerate(
        [
            "Bits",
            "MR Median",
            "10×25",
            "Komposite ≈ MR−250",
            "Kand. Sieb",
            "µs/Runde (Med/Med)",
            "relativ zu 512",
        ],
        1,
    ):
        ws_k.cell(23, i, h)
    style_header_row(ws_k, 23, 7)
    for i, (bits, r1, r2) in enumerate(blocks):
        r = 24 + i
        src = 6 + i
        cand = 18 + i
        ws_k.cell(r, 1, bits)
        ws_k.cell(r, 2).value = f"=MEDIAN(Rohdaten!I{r1}:I{r2})"
        ws_k.cell(r, 3, 250)
        ws_k.cell(r, 4).value = f"=B{r}-C{r}"
        ws_k.cell(r, 5).value = f"=C{cand}"
        # µs/round = median eabc ms * 1000 / median mr
        ws_k.cell(r, 6).value = f"=E{src}*1000/B{r}"
        ws_k.cell(r, 7).value = f"=F{r}/$F$24"
        ws_k.cell(r, 2).number_format = "0"
        ws_k.cell(r, 4).number_format = "0"
        ws_k.cell(r, 5).number_format = "0"
        ws_k.cell(r, 6).number_format = "0.0"
        ws_k.cell(r, 7).number_format = "0.00"
    style_body(ws_k, 24, 26, 1, 7)

    ws_k.merge_cells("A28:G29")
    ws_k["A28"] = (
        "Hinweis: µs/Runde verwendet das Verhältnis der Mediane (T10_EABC / MR). "
        "Der Median der Trial-Quotienten liegt bei 2048 etwas niedriger (~1111 µs). "
        "MR-Zerlegung: 10 Treffer × 25 Runden plus 1 Runde je Kompositum (Basis 2). "
        "openssl_sieved und eabc_mr haben identische Kandidatenzahlen."
    )
    ws_k["A28"].font = muted_font
    ws_k["A28"].alignment = Alignment(wrap_text=True, vertical="top")
    ws_k.row_dimensions[28].height = 22
    ws_k.row_dimensions[29].height = 22

    apply_widths(ws_k, {
        "A": 14, "B": 22, "C": 24, "D": 22, "E": 22, "F": 22, "G": 16, "H": 14,
    })
    for r in (5, 11, 17, 23):
        ws_k.row_dimensions[r].height = 28

    # helper chart data on Kennzahlen already in A5:E8
    chart1 = BarChart()
    chart1.type = "col"
    chart1.grouping = "clustered"
    chart1.title = "Median T10 (ms)"
    chart1.y_axis.title = "Millisekunden"
    chart1.x_axis.title = "Bitlänge"
    data = Reference(ws_k, min_col=2, min_row=5, max_col=5, max_row=8)
    cats = Reference(ws_k, min_col=1, min_row=6, max_row=8)
    chart1.add_data(data, titles_from_data=True)
    chart1.set_categories(cats)
    chart1.shape = 4
    chart1.legend.position = "b"
    chart1.width = 15
    chart1.height = 8
    ws_k.add_chart(chart1, "A31")

    chart2 = LineChart()
    chart2.title = "GMP / EABC (Median T10)"
    chart2.y_axis.title = "Verhältnis"
    chart2.x_axis.title = "Bitlänge"
    chart2.y_axis.scaling.min = 0
    chart2.y_axis.scaling.max = 1.2
    data2 = Reference(ws_k, min_col=4, min_row=11, max_col=4, max_row=14)
    cats2 = Reference(ws_k, min_col=1, min_row=12, max_row=14)
    chart2.add_data(data2, titles_from_data=True)
    chart2.set_categories(cats2)
    chart2.width = 15
    chart2.height = 8
    chart2.legend.position = "b"
    ws_k.add_chart(chart2, "E31")

    # --- Aussage ---
    ws_a.sheet_view.showGridLines = False
    ws_a.merge_cells("A1:F1")
    ws_a["A1"] = "Primes_EABC · Außenaussage"
    ws_a["A1"].font = title_font
    ws_a["A1"].fill = fill(NAVY)
    ws_a["A1"].alignment = Alignment(horizontal="left", vertical="center")
    ws_a.row_dimensions[1].height = 32

    ws_a.merge_cells("A2:F3")
    ws_a["A2"] = (
        "Vier Engines finden dieselbe Nachfolgerkette, failures=0, isolate=1. "
        "Der Sieb- und Testvorsprung von EABC existiert, schrumpft aber mit der Bitlänge. "
        "GMP mpz_nextprime bleibt in allen 15 Läufen 1,56–1,73× schneller. "
        "Kein Primzahlbeweis, kein Geschwindigkeitsclaim."
    )
    ws_a["A2"].font = Font(name="Calibri", size=11, color=INK)
    ws_a["A2"].fill = fill("EDF2F7")
    ws_a["A2"].alignment = Alignment(wrap_text=True, vertical="center")
    ws_a.row_dimensions[2].height = 22
    ws_a.row_dimensions[3].height = 22

    kpis = [
        ("A5", "B5", "15 / 15", "GMP schneller als EABC"),
        ("C5", "D5", "0", "failures / Agree-Brüche"),
        ("E5", "F5", "1,56–1,73×", "GMP-Faktor Median T10"),
        ("A8", "B8", "~8×", "Odd-Walk-Kandidaten / Sieb"),
        ("C8", "D8", "1,16 → 1,07", "Siebgewinn, schrumpfend"),
        ("E8", "F8", "1,39 → 1,11", "Testgewinn, schrumpfend"),
    ]
    for a, b, value, label in kpis:
        ws_a.merge_cells(f"{a}:{b}")
        cell = ws_a[a]
        cell.value = value
        cell.font = kpi_font
        cell.fill = fill(SURFACE)
        cell.alignment = Alignment(horizontal="center", vertical="center")
        cell.border = thin
        lab = ws_a.cell(ws_a[a].row + 1, ws_a[a].column, label)
        ws_a.merge_cells(
            start_row=ws_a[a].row + 1,
            start_column=ws_a[a].column,
            end_row=ws_a[a].row + 1,
            end_column=ws_a[a].column + 1,
        )
        lab.font = label_font
        lab.alignment = Alignment(horizontal="center", wrap_text=True)
        lab.fill = fill(SURFACE)
        lab.border = thin
    ws_a.row_dimensions[5].height = 28
    ws_a.row_dimensions[6].height = 28
    ws_a.row_dimensions[8].height = 28
    ws_a.row_dimensions[9].height = 28

    ws_a.merge_cells("A11:F11")
    ws_a["A11"] = "Was die Zahlen trennen — und was sie nicht belegen"
    ws_a["A11"].font = section_font
    ws_a["A11"].fill = fill(NAVY2)

    ws_a.merge_cells("A12:C16")
    ws_a["A12"] = (
        "Belegt\n"
        "• Dieselbe Next-Prime-Kette, probable-prime GMP reps=64 und OpenSSL\n"
        "• Sieb ~8× weniger Kandidaten als Odd-Walk, identisch für sieved und eabc_mr\n"
        "• Beide Gewinne schrumpfen: Sieb 1,16→1,07, Test 1,39→1,11\n"
        "• GMP-Abstand ist strukturell (flach 0,64→0,58), nicht Messrauschen\n"
        "• Trial 3 bei 2048 (EABC 1975 ms) verschiebt den Mittel, nicht den Median"
    )
    ws_a["A12"].font = body_font
    ws_a["A12"].alignment = Alignment(wrap_text=True, vertical="top")
    ws_a["A12"].fill = fill("F0FFF4")

    ws_a.merge_cells("D12:F16")
    ws_a["D12"] = (
        "Nicht belegt\n"
        "• Überlegenheit gegenüber GMP\n"
        "• Neue Primzahltheorie oder geschlossene Formel\n"
        "• Primzahlbeweis (nur probable prime)\n"
        "• FIPS- oder Krypto-Sicherheitsniveau\n"
        "• Publikationsstatistik: n=5 je Länge, rechtsschief; 30+ Trials und p25/p75 fehlen"
    )
    ws_a["D12"].font = body_font
    ws_a["D12"].alignment = Alignment(wrap_text=True, vertical="top")
    ws_a["D12"].fill = fill("FFF5F5")
    for r in range(12, 17):
        ws_a.row_dimensions[r].height = 18

    ws_a.merge_cells("A18:F18")
    ws_a["A18"] = (
        "Quelle: commercial_prime_benchmark/out_eabc/results.json · "
        "Primes_EABC.cpp · Apple M5 · 5 Trials × {512,1024,2048} · count=10"
    )
    ws_a["A18"].font = muted_font

    # live KPI from Kennzahlen
    ws_a["A20"] = "Live aus Kennzahlen"
    ws_a["A20"].font = section_font
    ws_a["A20"].fill = fill(NAVY2)
    ws_a.merge_cells("A20:C20")
    ws_a["A21"] = "GMP/EABC 512"
    ws_a["B21"] = "=Kennzahlen!D12"
    ws_a["C21"] = "GMP-Faktor 512"
    ws_a["D21"] = "=Kennzahlen!E12"
    ws_a["A22"] = "GMP/EABC 1024"
    ws_a["B22"] = "=Kennzahlen!D13"
    ws_a["C22"] = "GMP-Faktor 1024"
    ws_a["D22"] = "=Kennzahlen!E13"
    ws_a["A23"] = "GMP/EABC 2048"
    ws_a["B23"] = "=Kennzahlen!D14"
    ws_a["C23"] = "GMP-Faktor 2048"
    ws_a["D23"] = "=Kennzahlen!E14"
    for r in range(21, 24):
        ws_a.cell(r, 2).number_format = "0.000"
        ws_a.cell(r, 4).number_format = "0.00"
        for c in range(1, 5):
            ws_a.cell(r, c).font = body_font
            ws_a.cell(r, c).border = thin
            ws_a.cell(r, c).fill = fill(WHITE if r % 2 else SURFACE)

    apply_widths(ws_a, {"A": 22, "B": 16, "C": 22, "D": 16, "E": 18, "F": 18})

    # chart data copy for Aussage from Kennzahlen via formulas
    ws_a["A26"] = "Bits"
    ws_a["B26"] = "GMP"
    ws_a["C26"] = "OpenSSL raw"
    ws_a["D26"] = "OpenSSL sieved"
    ws_a["E26"] = "EABC"
    style_header_row(ws_a, 26, 5)
    for i in range(3):
        r = 27 + i
        src = 6 + i
        ws_a.cell(r, 1).value = f"=Kennzahlen!A{src}"
        ws_a.cell(r, 2).value = f"=Kennzahlen!B{src}"
        ws_a.cell(r, 3).value = f"=Kennzahlen!C{src}"
        ws_a.cell(r, 4).value = f"=Kennzahlen!D{src}"
        ws_a.cell(r, 5).value = f"=Kennzahlen!E{src}"
        for c in range(2, 6):
            ws_a.cell(r, c).number_format = "0.0"
    style_body(ws_a, 27, 29, 1, 5)

    ch = BarChart()
    ch.type = "col"
    ch.grouping = "clustered"
    ch.title = "Median T10 nach Engine"
    ch.y_axis.title = "ms"
    dref = Reference(ws_a, min_col=2, min_row=26, max_col=5, max_row=29)
    cref = Reference(ws_a, min_col=1, min_row=27, max_row=29)
    ch.add_data(dref, titles_from_data=True)
    ch.set_categories(cref)
    ch.legend.position = "b"
    ch.width = 18
    ch.height = 9
    ws_a.add_chart(ch, "A31")

    wb.save(OUT_PATH)
    print(OUT_PATH)


if __name__ == "__main__":
    main()
