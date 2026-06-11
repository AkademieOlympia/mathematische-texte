# Pfad B: Hurwitz-Gitter und arithmetische Resonanz

Verknüpfung der **ℚ-Rundweg-Analyse** (`Rundweg.py`, `PAPER_Q_STRUCTURE.md`) mit dem **Gabriel–Hurwitz-Protokoll** (`Gabriel Hurwitz.tex`).

## Identifikation

| EABC | Quaternionenbasis |
|------|-------------------|
| E | 1 |
| A | i |
| B | j |
| C | k |

\(q = e\cdot 1 + a\cdot i + b\cdot j + c\cdot k\), Algebra \(H = \mathrm{QuaternionAlgebra}(\mathbb{Q},-1,-1)\).

## Π\_Γ-Kollaps

Holonomer Rundweg \(\Gamma = R\cdot K\) wie in `Rundweg.py`:

\[
\Pi_\Gamma = \frac{1}{4}(I + \Gamma + \Gamma^2 + \Gamma^3).
\]

Auf Koordinaten \((e,a,b,c)\) (Paper-Notation aus ℚ-Analyse):

\[
\Pi_\Gamma(e,a,b,c) = \bigl(e,\; \tfrac{a+c}{2},\; b,\; \tfrac{a+c}{2}\bigr).
\]

**Eigenraum \(\lambda=0\):** \(V_{\lambda=0} = \mathrm{span}(A-C)\), Basis \((0,1,0,-1)\) ↔ \(q_{\mathrm{asym}} = i-k\).

## Primquaternionen und Ideale

| Prim | \(\pi_p\) | Norm (Sage) |
|------|-----------|-------------|
| 5 | \(2+i\) | 5 |
| 7 | \(2+i+j+k\) | 7 |

**Nicht** verwenden: \((3/2)(1+i+j+k)\) — Norm **9**, nicht 7.

**Linksideal** \(I_p = H\cdot\pi_p = \{h\pi_p \mid h\in H\}\) (Hurwitz-Maximalordnung).

**Resonanzpunkt (Pfad B):** \(0 \neq q \in V_{\lambda=0}\) mit \(q \in I_p\) (nichttrivialer Schnitt mit Primideal-Sphäre).

Operative Prüfung in `rundweg_hurwitz_resonanz.sage`: \(h = q\pi_p^{-1}\); Resonanz genau wenn \(h\) Hurwitz ist. Notwendig: \(N(\pi_p)\mid N(q)\).

## Verifikation

```bash
./run_hurwitz_resonanz.sh
```

## M = 113160

| Quelle | Bedeutung |
|--------|-----------|
| `H32_eabc.py` | Demo-Phasenzentrum `test_H32(113160)`; Mod-30-Anker, 32-Slot-Hülle `H32` |
| `Ptolo Norm.py` / `Ptolo Norm Matrix.py` | Statistik-Zentrum \(10^7 + 113160\) |
| `Jitter Zeit.py` | Bamberger-Zeitwürfel: `n_k=round(α·t_k)`, EABC mod 12, `q∈{1,i,j,k}`, Schwung `u=conj(q_k)q_{k+1}`; CSV `bamberger_zeitwuerfel_active_*.csv` |

`scan_resonanz(M)`: H32-Slot-Zähler → Aggregat-Quaternion. **`zeitwuerfel_resonanz_fuer_M(M)`:** lädt aktive CSV-Zustände, koppelt per `n_k mod 30` an belegte H32-Slots von \(M\), prüft \(I_5,I_7\) (Aggregat + distinkte Schwünge).

**Stand M=113160 (Verifikation `./run_hurwitz_resonanz.sh`):** H32-Aggregat \((4,5,3,5)\), \(N=75\) — kein \(I_5,I_7\). Zeitwürfel-Aggregat aus Simulation (662 aktive Ereignisse) ebenfalls **keine** Linksideal-Resonanz; ganzzahlige Schwung-Typen \(\pm i,\pm j,\pm k\) ebenfalls **nein**. Mod-30-Phasenkopplung verbindet die Module **konzeptionell** (gleiche EABC-/Quaternion-Kodierung), liefert aber **keine** nichttriviale Hurwitz-Ideal-Resonanz in diesem Lauf.

## Empfehlung Einheiten-Scan

**Ja, aber zweistufig:**

1. **Jetzt:** Ideale \(I_5, I_7\) über \(q\pi^{-1}\) (schnell, exakt für maximale Hurwitz-Ordnung in der rationalen Algebra).
2. **Später:** Systematischer Scan der **24** Hurwitz-Einheiten (Eichler–Dupont-Liste) und deren Linksklassen modulo \(\pi_p\) — lohnt sich für Paper-Satz „arithmetischer Schalen-Kollaps“, sobald \(V_{\mathrm{arith}}\) nicht nur \(\mathbb{Q}^4\) ist.

Ohne Einheiten-Scan bleibt Pfad B auf **Einzelzustände** (z. B. \(q_{\mathrm{lokal}}\), \(q_{\mathrm{asym}}\)) und H32-Aggregate beschränkt; das reicht für die erste arXiv-Verifikation.

## Masse vs. Resonanz

| | Eigenraum | Hurwitz | Rolle |
|---|-----------|---------|--------|
| **Fluktuation / Resonanz** | \(V_{\lambda=0}\), \(A-C\) | \(q\in I_p\) (Pfad B) | Antisymmetrisch, \(\Pi_\Gamma q=0\) |
| **Masse-Konstitution** | \(\mathrm{im}(\Pi_\Gamma)\), \(\lambda=1\) | \(N(\Pi_\Gamma q)\ge 1\) | Zwischen ED-sichtbaren EABC-Objekten mit Zweiweg |

Prinzip und Formalisierung: `PAPER_MASSE_KONSTITUTION.md`, LaTeX `paper_masse_axiom.tex`.

## Quadratischer Schwung (Erweiterung)

**Definition:** \(u = \bar{q}\,q'\) zwischen Zustand \(q\) und Referenz/Nachbar \(q'\). Resonanzkriterium (Kapitel 8): \(N(u)\equiv 0 \pmod p\) für \(p\in\{5,7\}\), optional \(u\in I_p\).

| Skript | Inhalt |
|--------|--------|
| `rundweg_schwung_resonanz.sage` | `u = q.conjugate()·q'`, \(N(u)\bmod 5,7\), Linksideal \(I_5,I_7\); Testpaare + `scan_schwung_resonanz(M)` |
| `./run_schwung_resonanz.sh` | Sage-Start (analog `run_hurwitz_resonanz.sh`) |

LaTeX-Kapitel 5 & 8: `kapitel_5_8_schalen_kollaps.tex` (Schalenstruktur \(\mathcal{S}_{\mathrm{halb}}\), erweitertes Kollaps-Theorem).

**Linear vs. quadratisch:** `./run_hurwitz_resonanz.sh` prüft \(q\in I_p\) (negativ für \(q_{\mathrm{asym}}\)); `./run_schwung_resonanz.sh` prüft \(N(\bar{q}q')\) modulo \(\pi_5,\pi_7\).

**Stand Schwung-Scan (Verifikation `./run_schwung_resonanz.sh`):** Testpaare \(q_{\mathrm{lokal}}\times q_{\mathrm{glatt}}\), \(q_{\mathrm{lokal}}\times E\), \(q_{\mathrm{asym}}\times q_{\mathrm{lokal}}\): \(N(u)\equiv 0\pmod 5\) (**ja**, z.\,B. \(N=885,60,30\)), aber **kein** \(u\in I_5,I_7\). M=113160: 5 distinkte CSV-Schwünge und 661 Label-Übergänge — **keine** mod-5/7-Resonanz, **kein** Ideal-Treffer.

Referenz-Schwünge: \(N(u)\equiv 0\pmod 7\) **nie**; Ideal \(I_5\cup I_7\) bei keinem getesteten \(u\).  
\(M=113160\): Schwung-Scan **0** Treffer (mod 5/7 und Ideal).  
LaTeX-Abschnitt „Numerische Verifikation“: `kapitel_5_8_schalen_kollaps.tex` (Proposition + Remark, ohne Overclaim).

## Paper

LaTeX-Satz: `paper_hurwitz_snippet.tex`, Kapitel 5/8: `kapitel_5_8_schalen_kollaps.tex`. ℚ-Tabelle bleibt in `PAPER_Q_STRUCTURE.md` (Pfad A/B-Verweis).
