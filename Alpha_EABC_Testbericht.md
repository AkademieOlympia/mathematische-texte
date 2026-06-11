# Testbericht: EABC/Zeta-Kalibrierung der Feinstrukturkonstante α

**Projekt:** Alpha_EABC — Riemann-Nullstellen als Zeuge für ε = 1/α − 137  
**Skript:** `Alpha_EABC.py :: Zusatztests-v2 :: 2026-06-07`  
**Daten:** `zeros6.npy` — 2 001 052 Riemann-Nullstellen (γ₁ ≈ 14,135 … γ_N ≈ 1,13×10⁶)  
**Referenz:** CODATA 2022, ε_CODATA = 0,035999177 (1/α = 137,035999177)  
**Unsicherheit (σ-Welt):** u(ε) = u(1/α) = 2,1×10⁻⁸ (0,58 ppm auf der ε-Skala)

---

## 1. Executive Summary

Das EABC-Programm kalibriert die elektromagnetische Feinstrukturkonstante über den Median-Defekt normierter Riemann-Lücken:

\[
\varepsilon = 1 - \mathrm{median}(\text{lokale normierte Gaps}), \quad 1/\alpha = 137 + \varepsilon.
\]

**Kernbefund:** Eine rein a-priori festgelegte Fensterregel — center **c = 113**, window **w = 4**, abgeleitet aus {16, 4, 20, 137} via  
c = 137 − (16+4) − 4 = **113**, w = **4** — liefert ohne freie Parameter eine ε-Schätzung, die über eine vierstufige Korrekturkette (Median → Quartil → n(n±1) → σ-Randspannung) auf |Δε| ≈ **6,6×10⁻⁷** gegenüber CODATA konvergiert. Das ist numerisch bemerkenswert (~9×, ~19× und schließlich ~**14,5×** besser als die vorherigen Stufen), aber **statistisch nicht CODATA-kompatibel**: selbst ε_final liegt bei ~**31σ** Entfernung von ε_CODATA (u(ε) ≈ 2,1×10⁻⁸).

**Overfitting-Warnung:** Die σ-Randspannung ist die vierte post-hoc-Korrektur in der Kette. Kontrollfenster (120,5), (105,6), (130,4) liefern |Δ| ~ **0,01–0,1** — nur die 113/114-Kombination trifft ~10⁻⁷. Die numerische Präzision ist damit **stark fenster-spezifisch**, nicht universell.

**Strukturell gesichert:** Die algebraische Umgebung n = 113 — Primzahl p₃₀ mit 30 = 2(16−1), prim 2n+1 = 227, mod-16-Sequenz [113,114,227,228] ≡ [1,2,3,4] — ist ein kohärentes, vorhersagbares Muster unabhängig vom numerischen Fit.

**Schwach / spekulativ:** Sedenion-Interpretation (P₀ = 16, J/4 = 4), physikalische Monopol-BPS-Rücktests (trivial exakt by construction), die „harte Signatur“ mit p_hard ≈ 2,1×10⁻⁵ (3 Treffer ≈ Zufallserwartung bei ~142 000 Kandidaten), und die **σ-Randspannung** als vierte Korrekturstufe (Overfitting-Risiko: Kontrollfenster scheitern).

**Fazit:** Das Projekt liefert eine **intern konsistente, algebraisch motivierte Kalibrierungskette** mit bemerkenswerter numerischer Präzision auf der ε-Skala (ε_final: |Δ| ≈ 6,6×10⁻⁷, ~31σ von CODATA), aber **keinen statistischen Beweis** für physikalische Vorhersagekraft jenseits post-hoc-Selektion. Die stärkste Hypothese betrifft die **mod-16-/Primzahl-Struktur um n = 113**, nicht die absolute CODATA-Übereinstimmung. Die vierte Korrekturstufe (H9) ist **stark fenster-spezifisch**.

---

## 2. Methodik & Testmatrix

### 2.1 Datenpipeline

| Stufe | Definition |
|-------|------------|
| Eingabe | Riemann-Nullstellen γₙ aus `zeros6.npy` |
| Roh-Gaps | Δγ = γₙ₊₁ − γₙ |
| Normierung | gap_norm = Δγ · ln(γ/2π) / (2π) (Montgomery-Odlyzko-Skala) |
| Fenster | center c, Breite w → 2w lokale Gaps [c−w : c+w] |
| Baseline | ε₀ = 1 − median(lokal) |
| EABC-Korrektur | ε_EABC = ε₀ + (q₂₅ + q₇₅ − 2) / 320, mit 320 = 16·20 |
| n(n±1)-Paar | Kern (113,4) + Rand (114,3), Gewichte n/(2n+1), (n+1)/(2n+1) |
| σ-Randspannung (2. Korrektur) | ε_final = ε_pair_EABC − σ̄/(16·20·137), σ̄ = (σ_core+σ_edge)/2 |
| Referenz | ε_CODATA = 0,035999177 |

### 2.4 Korrekturkette (eingefrorene Formeln)

```
Median  →  Quartil  →  n(n±1)  →  σ-Randspannung
 ε₀         ε_EABC      ε_pair      ε_final
```

\[
\boxed{
\begin{aligned}
\varepsilon_0 &= 1 - \mathrm{median}(\text{lokal}) \\
\varepsilon_{\mathrm{EABC}} &= \varepsilon_0 + \frac{q_{25} + q_{75} - 2}{320} \\
\varepsilon_{\mathrm{pair}} &= \frac{n}{2n+1}\,\varepsilon_{\mathrm{EABC}}^{\mathrm{core}} + \frac{n+1}{2n+1}\,\varepsilon_{\mathrm{EABC}}^{\mathrm{edge}} \\
\varepsilon_{\mathrm{final}} &= \varepsilon_{\mathrm{pair}} - \frac{\bar\sigma}{16 \cdot 20 \cdot 137}, \quad \bar\sigma = \frac{\sigma_{\mathrm{core}} + \sigma_{\mathrm{edge}}}{2}
\end{aligned}
\}
\]

Alle Normierungen (320 = 16·20, 137) sind reine Modellzahlen — **kein Fit auf CODATA**.

### 2.2 A-priori-Regel (eingefroren)

| Parameter | Wert | Herleitung |
|-----------|------|------------|
| c | 113 | 137 − (16+4) − 4 |
| w | 4 | Modellzahl {4} |
| Primstützung | 113 = p₃₀ | 30 = 2(16−1) |
| Modellmenge | {16, 4, 20, 137} | Divisionsalgebren / EABC |

**Kein Fit:** center und window werden vor dem Datenvergleich festgelegt.

### 2.3 Testmatrix

| ID | Test | Zweck |
|----|------|-------|
| T1 | Baseline Median (113,4) | A-priori ohne Korrektur |
| T2 | ε_EABC-Quartilkorrektur | Modellinterne Feinjustierung |
| T3 | n(n±1)-Paar (Median / ε_EABC) | Kern-Rand-Kopplung |
| T4 | Default-Kontrolle (138,20) | Vergleich zur alten Heuristik |
| T5 | Fensterrobustheit (w = 1…20) | Stabilität bei c = 113 |
| T6 | Nachbarschaft (c = 100…130, w = 4) | Lokale Eindeutigkeit |
| T7 | Block-Wiederkehr (c = 113 + 1000k) | Skalen-Invarianz |
| T8 | Monte-Carlo-Fensterpermutation (20 000) | Zufalls-Hintergrund (Gaps) |
| T9 | GUE-Surrogat (5 000) | RMT-Kontrolle |
| T10 | Mehrdimensionale Signatur Q | Alpha + Monopol + M_X + Kontrast |
| T11 | Harte Signatur | w ∈ {3,4}, lokaler Rang, Kontrast ≥ 100 |
| T12 | Strukturtest n(n±1) | Prim/mod-16-Diagnostik |
| T13 | Modulo-Profil [113,114,227,228] | EABC-Sequenztest |
| T14 | Kontrollscan 137-Zone | Post-hoc Top-5 (100≤c≤201, 3≤w≤30) |
| T15 | BPS-Monopol-Rücktest | M_X/α = M_M (konstruktiv) |

---

## 3. Ergebnisse

### 3.1 Haupttabelle: ε-Schätzungen vs. CODATA

*(Frisch berechnet, Lauf vom 2026-06-07)*

| Methode | (c, w) | ε | \|Δε\| | σ(CODATA)* | Rel. Fehler |
|---------|--------|---|--------|------------|-------------|
| **CODATA** | — | 0,035999177 | 0 | 0 | — |
| Median A-priori | (113, 4) | 0,036181257 | 1,821×10⁻⁴ | 8 671σ | 0,506 % |
| **ε_EABC** | (113, 4) | 0,036019306 | 2,013×10⁻⁵ | 959σ | 0,056 % |
| ε_EABC Rand | (114, 3) | 0,035998279 | **8,98×10⁻⁷** | **42,8σ** | 0,0025 % |
| **ε_EABC Paar** | (113,4)+(114,3) | 0,036008746 | 9,57×10⁻⁶ | 456σ | 0,026 % |
| **ε_final** | (113,4)+(114,3)+σ̄ | **0,035999839** | **6,62×10⁻⁷** | **~31σ** | **0,0018 %** |
| Default-Kontrolle | (138, 20) | 0,023427151 | 1,257×10⁻² | 5,99×10⁵σ | 34,9 % |
| MC-Best (post-hoc) | start=1 318 537 | ≈0,035996337 | 2,84×10⁻⁶ | 135σ | — |

*\*σ(CODATA) = |Δε| / u(ε), u(ε) = 2,1×10⁻⁸*

**Verbesserungsfaktoren (gegenüber Median A-priori):**

| Vergleich | Faktor |
|-----------|--------|
| ε_EABC vs. ε₀ | **9,05×** |
| ε_EABC-Paar vs. Median-Paar | **19,0×** |
| **ε_final vs. ε_EABC-Paar** | **14,5×** |
| ε_EABC-Paar vs. ε_EABC-Kern allein | 0,48× (Paar schlechter als Rand allein) |

### 3.2 EABC-Quartilkorrektur im Detail (113, 4)

| Größe | Wert |
|-------|------|
| q₂₅ | 0,632604 |
| q₇₅ | 1,315572 |
| A = q₂₅ + q₇₅ − 2 | −0,051824 |
| Korrekturterm | −1,620×10⁻⁴ |
| σ_loc | 0,384040 |

Formel: ε_EABC = ε₀ + (q₂₅ + q₇₅ − 2) / **320**, wobei 320 = 16·20 reine Modellzahl ist.

### 3.3 n(n±1)-Kern-Rand-Kopplung

| Komponente | (c, w) | ε (Median) | ε_EABC | Gewicht |
|------------|--------|------------|--------|---------|
| Kern | (113, 4) | 0,036181257 | 0,036019306 | 113/227 ≈ 0,498 |
| Rand | (114, 3) | 0,036181257 | **0,035998279** | 114/227 ≈ 0,502 |
| **Paar** | gewichtet | 0,036181257 | **0,036008746** | — |

Bemerkung: Median-Kern und Median-Rand sind numerisch identisch (ε₀); die Kopplung wirkt erst über unterschiedliche ε_EABC-Werte. Der Rand allein (114,3) ist der **beste Einzelkandidat**.

Struktur: 2n+1 = **227** (prim), Gewichte n/(2n+1) und (n+1)/(2n+1) summieren zu 1.

### 3.4 Zweite EABC-Randkorrektur (σ-Randspannung)

| Größe | Wert |
|-------|------|
| σ_core (113, 4) | 0,384040 |
| σ_edge (114, 3) | 0,396946 |
| σ̄ = (σ_core + σ_edge)/2 | 0,390493 |
| Korrekturterm σ̄/(16·20·137) | 8,907×10⁻⁶ |
| ε_pair_EABC | 0,036008746 |
| **ε_final** | **0,035999839** |
| \|Δε_final\| vs. CODATA | **6,62×10⁻⁷** |

Formel: ε_final = ε_pair_EABC − σ̄/(16·20·137). Der Term nutzt die lokale Gap-Standardabweichung beider Fenster (Kern und Rand) als „Randspannung“ — normiert durch reine Modellzahlen 16, 20, 137.

**Verbesserung gegenüber ε_EABC-Paar:** Faktor **14,5×** (9,57×10⁻⁶ → 6,62×10⁻⁷). Gegenüber dem Median insgesamt: ~275×.

### 3.5 Overfitting-Diagnose: Kontrollfenster

Um die Spezifität der vierten Korrekturstufe zu prüfen, wurde dieselbe vierstufige Kette auf **alternative Fensterpaare** angewendet (gleiche n/(2n+1)-Gewichtung, gleiche σ-Formel):

| Fensterpaar | ε_final | \|Δε\| vs. CODATA | Urteil |
|-------------|---------|-------------------|--------|
| **(113,4)+(114,3)** | 0,035999839 | **6,62×10⁻⁷** | Treffer ~10⁻⁷ |
| (120,5)+(121,4) | 0,002051634 | 3,39×10⁻² | ~5×10² schlechter |
| (105,6)+(106,5) | 0,136274148 | 1,00×10⁻¹ | ~1,5×10⁵ schlechter |
| (130,4)+(131,3) | −0,015578530 | 5,16×10⁻² | ~7,8×10⁴ schlechter |

**Befund:** Nur die a-priori-Regel (113/114) erreicht CODATA-Nähe auf 10⁻⁷-Skala. Kontrollfenster liefern |Δ| ~ **0,01–0,1** — die vierte Korrektur ist **stark overfitting-anfällig** und darf nicht als universelle Vorhersageformel interpretiert werden.

### 3.6 Robuste Lokalität

**Fensterrobustheit (c = 113):** w = 3 und w = 4 liefern identisch das Minimum (|Δε| = 1,82×10⁻⁴); ab w = 5 verschlechtert sich ε stark (|Δε| > 10⁻²).

| Rang | w | ε | \|Δε\| |
|------|---|-----|--------|
| 1 | 3 | 0,036181257 | 1,82×10⁻⁴ |
| 2 | 4 | 0,036181257 | 1,82×10⁻⁴ |
| 3 | 11 | 0,029167718 | 6,83×10⁻³ |

**Nachbarschaft (w = 4):** c = 113 ist eindeutig bestes Zentrum in [100, 130]; nächster Konkurrent c = 115 mit |Δε| ≈ 1,57×10⁻² (~86× schlechter).

| Rang | c | ε | \|Δε\| |
|------|---|-----|--------|
| 1 | 113 | 0,036181257 | 1,82×10⁻⁴ |
| 2 | 115 | 0,020298236 | 1,57×10⁻² |

**Block-Wiederkehr:** Nur Block 0 (center = 113) trifft CODATA-Nähe; ab k ≥ 1 driftet ε um Größenordnungen — **keine Skalen-Invarianz**.

| Block | center | ε | \|Δε\| |
|-------|--------|-----|--------|
| 0 | 113 | 0,036181257 | 1,82×10⁻⁴ |
| 1 | 1113 | 0,198022413 | 1,62×10⁻¹ |
| 2 | 2113 | −0,176295375 | 2,12×10⁻¹ |

### 3.7 Harte EABC-Signatur

Kriterien: w ∈ {3, 4}, ε_EABC-Rang #1 in c ± 15, Nachbar-Kontrast ≥ 100.

| Kenngröße | Wert |
|-----------|------|
| Kandidaten gesamt | 142 129 |
| Harte Treffer | **3** |
| p_hard | **2,11×10⁻⁵** |
| Erwartung (Zufall) | ≈ **3,0** |

| Rang | (c, w) | \|Δε_EABC\| | Strukturqualität |
|------|--------|-------------|-------------------|
| 1 | (113, 4) | 2,01×10⁻⁵ | n prim, 2n+1 prim, n ≡ 1 (mod 16) |
| 2 | (114, 3) | **8,98×10⁻⁷** | gleiche n-Struktur (n = 113) |
| 3 | (3876, 4) | 1,53×10⁻⁵ | n = 3875 = 5³·31, schwach |

### 3.8 Mehrdimensionale Signatur Q (A-priori)

| Größe | Wert |
|-------|------|
| s1_alpha_rel | 0,00506 |
| s2_monopole_rel | 1,33×10⁻⁶ |
| s3_MX_rel | 1,33×10⁻⁶ |
| s4_local_contrast | 611,2 |
| s5_window_std | 0,00749 |
| Q_target | 9,14×10¹⁸ |

### 3.9 Kontrollscan & Default

**Top-5 im Scan 100 ≤ c ≤ 201, 3 ≤ w ≤ 30:** Alle führenden Treffer liegen bei (113,3/4) und (114,3) — konsistent mit A-priori, aber der Scan-Raum ist klein (~2 800 Paare) und damit anfällig für LEE.

**Default (138, 20):** |Δε| = 1,26×10⁻² — ~**690× schlechter** als unkorrigierter A-priori-Median.

### 3.10 BPS-Monopol-Rücktest

Für jeden Kandidaten α: M_X = α · M_M, BPS-Rückrechnung M_X/α = M_M.

| Beobachtung | Status |
|-------------|--------|
| BPS-relativer Fehler | **exakt 0** (by construction) |
| Physikalische Aussagekraft | **keine** — tautologisch |

Physikalische Ableitungen (ε_EABC, c=113, w=4):

| Größe | Wert |
|-------|------|
| 1/α_EABC | 137,036019306 |
| α_EABC | 0,00729735149 |
| M_X_EABC | 1,992×10¹⁶ GeV |
| α_m | 34,259005 |
| g/e | 68,518010 |

---

## 4. Statistische Bewertung (σ, LEE, p)

### 4.1 σ-Welt (CODATA-Kompatibilität)

**Kein Kandidat liegt innerhalb 1σ von CODATA.** Selbst ε_final (~31σ) ist weit entfernt; der beste Einzeltreffer (114,3, ε_EABC) liegt bei ~43σ. Für echte Übereinstimmung bräuchte |Δε| ≤ 2,1×10⁻⁸ (0,58 ppm). In der σ-Welt ist das Projekt **numerisch interessant, physikalisch nicht bestätigt**.

| Methode | \|Δε\| | rel. ppm | σ vs CODATA |
|---------|--------|----------|-------------|
| Default (138,20) | 1,26×10⁻² | 3,5×10⁵ | ~6×10⁵σ |
| A-priori Median (113,4) | 1,82×10⁻⁴ | 5,1×10³ | ~8,7×10³σ |
| ε_EABC Kern (113,4) | 2,01×10⁻⁵ | 559 | ~959σ |
| ε_EABC-Paar 113/114 | 9,57×10⁻⁶ | 266 | ~456σ |
| Edge ε_EABC (114,3) | 8,98×10⁻⁷ | **25** | **~43σ** |
| **ε_final (113/114+σ̄)** | **6,62×10⁻⁷** | **18** | **~31σ** |

### 4.2 Monte-Carlo-Hintergrund

| Test | trials | hits | p_estimate | Interpretation |
|------|--------|------|------------|----------------|
| Gap-Permutation (113,4), \|Δε₀\| | 20 000 | 29 | 0,00145 | Median-Treffer ~0,15 %-Quantil — **moderat selten** |
| GUE-Surrogat | 5 000 | 7 | 0,0014 | Konsistent mit RMT-typischer Gap-Statistik |
| Signatur Q (mehrdimensional) | 50 000 | 10 016 | 0,200 | Q-Signatur **nicht signifikant** (p ≈ 0,2) |

### 4.3 Look-Elsewhere-Effekt (LEE)

| Suchraum | Größe | Risiko |
|----------|-------|--------|
| Harte Signatur | ~142 000 (c,w)-Paare | 3 Treffer ≈ erwartete Zufallsrate |
| Kontrollscan Top-5 | ~2 856 Paare | Post-hoc; A-priori (113,4) im Top-Cluster |
| MC bestes Fenster | 2×10⁶ mögliche Starts | Best \|Δε\| = 2,84×10⁻⁶ ohne algebraische Begründung |

**LEE-Korrektur:** Bei N ≈ 142 129 Kandidaten und p_hard ≈ 2,1×10⁻⁵ ergibt sich E[Treffer] = N·p ≈ **3,0** — beobachtete 3 Treffer = **kein überschüssiges globales Signal**.

- Lokal (ohne LEE): p ≈ 2×10⁻⁵ ≈ **4,1σ** (1-seitig) — gilt nur für einen vorgegebenen Test
- **Bonferroni** (α = 0,05): p_korr ≈ 3,5×10⁻⁷ → global ~**5σ**-Schwelle; |Δε| müsste < **1,0×10⁻⁷** (~2,9 ppm) sein — Edge (25 ppm) scheitert auch global

**A-priori vs. post-hoc:**

| Regel | LEE-Status | CODATA-σ |
|-------|------------|----------|
| c=113, w=4 (Modellzahlen) | **a-priori** (vor Scan) | ~8700σ |
| ε_EABC-Quartilkorrektur /320 | post-hoc (Modellformel) | ~960σ |
| Paar 113/114, Edge (114,3) | post-hoc (aus hartem Scan) | ~456σ / ~43σ |
| σ-Randspannung (ε_final) | post-hoc (4. Stufe) | ~31σ, **Overfitting-Risiko** |

### 4.4 p-Wert-Zusammenfassung

| Hypothese | p | Urteil |
|-----------|---|--------|
| Median (113,4) zufällig so gut | ~0,0015 | Schwach signifikant |
| Harte Signatur zufällig | p_hard ≈ 2,1×10⁻⁵, 3 Hits ≈ E[N] | **Nicht über Zufall hinaus** |
| Q-Signatur überlegen | ~0,20 | Nicht signifikant |
| CODATA-Übereinstimmung | σ ≫ 1 | **Verworfen** |

---

## 5. Algebraische Strukturbefunde (mod 16 [1,2,3,4])

### 5.1 A-priori-Zahl 113

| Eigenschaft | Wert |
|-------------|------|
| Prim | ja (p₃₀) |
| 30 = 2(16−1) | Primindex-Kopplung an Modellzahl 16 |
| n mod 16 | 1 |
| n mod 4 | 1 |
| Faktoren | 113 (prim) |

### 5.2 EABC-Sequenz [n, n+1, 2n+1, 2n+2] = [113, 114, 227, 228]

| x | prim | x mod 16 |
|---|------|----------|
| 113 | ✓ | **1** |
| 114 | ✗ (2·3·19) | **2** |
| 227 | ✓ | **3** |
| 228 | ✗ (2²·3·19) | **4** |

**Befund:** Die mod-16-Restfolge **[1, 2, 3, 4]** ist exakt und unabhängig vom numerischen ε-Fit reproduzierbar. Dies ist das **stärkste algebraische Muster** des Projekts.

Weitere Relationen:

- 113 + 114 = 227 = 2·113 + 1
- 228 = 2·114 = 2(n+1)
- 320 = 16·20 erscheint als Korrektur-Normierung

### 5.3 Tiefe Modulo-Struktur

| x | mod 4 | mod 16 | mod 20 | mod 137 |
|---|-------|--------|--------|---------|
| 113 | 1 (+1, quarter) | 1 (+1) | 13 | 113 |
| 114 | 2 (half) | 2 | 14 | 114 |
| 227 | 3 (−1, quarter) | 3 | 7 | 90 |
| 228 | 0 | 4 (quarter) | 8 | 91 |

### 5.4 Sedenion-Interpretation (interpretativ)

| Formel | Wert | Rolle |
|--------|------|-------|
| P₀ = 16 | Basis-Divisionsalgebra | Modellanker |
| J/4 = 4 | Window-Parameter w | post-hoc-Deutung |
| 137 − 16 − 4 − 4 = 113 | center | kombinatorische Regel |

**Status:** Plausibilisierende Metaphorik, **kein** unabhängiger Beweis.

### 5.5 Vergleich harte Treffer

Nur **113/114/227** erfüllt gleichzeitig: n prim, 2n+1 prim, n ≡ 1 (mod 16). Der Treffer 3876/4 (n = 3875 = 5³·31) erfüllt diese Struktur **nicht** — spricht für selektive Strukturkorrelation, aber mit nur 3 Hits statistisch schwach.

---

## 6. Gesichert vs. spekulativ

| **Gesichert (reproduzierbar)** | **Spekulativ / nicht gesichert** |
|--------------------------------|----------------------------------|
| A-priori-Regel (113,4) ist vor Fit festgelegt | Physikalische Deutung über Riemann-Hypothese hinaus |
| ε_EABC verbessert ε₀ um Faktor ~9 | ε_EABC trifft CODATA innerhalb Messunsicherheit |
| ε_final verbessert ε_EABC-Paar um Faktor ~14,5 | ε_final trifft CODATA innerhalb u(ε) (~31σ entfernt) |
| Rand (114,3) ist bester Einzelwert (ε_EABC) | n(n±1)-Paar ist global optimal (schlechter als Rand allein) |
| mod-16-Sequenz [1,2,3,4] für n=113 | Sedenion-Formeln als Erklärung |
| Lokales Minimum bei c=113, w∈{3,4} | Globale Eindeutigkeit über alle c |
| Block 0 speziell; höhere Blöcke scheitern | Skalen-unabhängige „Physik“ |
| BPS M_X/α = M_M exakt | Monopol-Observables als unabhängiger Test |
| 3 harte Treffer ≈ Zufallserwartung | Harte Signatur als Discovery-Kriterium |
| Default (138,20) deutlich schlechter | Kontrollscan rechtfertigt allein die Regel |
| CODATA-Referenz mit u(ε) = 2,1×10⁻⁸ | „α innerhalb CODATA-1σ vorhergesagt“ |
| Kontrollfenster (120,5), (105,6), (130,4) scheitern | σ-Randspannung als universelle Korrektur |

**Paper-tauglich (CAN):** A-priori-Spezifikation, numerische Nähe als Beobachtung, BPS als algebraische Folge, globale Zufallserklärung der harten Treffer.

**Nicht behauptbar (CANNOT):** CODATA-1σ-Vorhersage, 5σ-Entdeckung ohne LEE, BPS als empirischer α-Beweis, Edge/Paar/ε_final als a-priori-Vorhersage, σ-Randspannung als universelle Korrekturformel.

---

## 7. Hypothesen H1–H9

### H1 — A-priori-Fensterhypothese

**Aussage:** Das Fenster (c, w) = (137 − 16 − 4 − 4, 4) = (113, 4) ist die einzige a-priori aus {16, 4, 20, 137} gebildete Regel, die in der 137-Zone der normierten Nullstellen-Lücken ein lokales ε-Minimum liefert.

**Testbarkeit:** Vollständiger deterministischer Scan aller kombinatorischen Regeln aus der Modellmenge mit eingefrorener Bewertungsmetrik; Vergleich der Rangverteilung.

**Status:** Teilweise gestützt (lokal bestätigt, global nicht bewiesen).

### H2 — EABC-Quartilkorrektur

**Aussage:** Die Korrektur (q₂₅ + q₇₅ − 2)/320 kompensiert systematische Median-Verschiebung durch lokale Gap-Asymmetrie und ist keine numerische Overfitting-Maßnahme.

**Testbarkeit:** Kreuzvalidierung auf unabhängigen Nullstellen-Dateien; Permutation der lokalen Gaps; Vergleich mit alternativen Modellnormierungen (z. B. 256, 272).

**Status:** Numerisch wirksam (9×), theoretische Begründung offen; **Overfitting-Risiko** steigt mit jeder weiteren Korrekturstufe (vgl. H9).

### H3 — n(n±1)-Kern-Rand-Kopplung

**Aussage:** Physikalisch relevante ε-Schätzung erfordert gewichtete Kombination aus Kern (n, w=4) und Rand (n+1, w=3) mit Gewichten n/(2n+1), (n+1)/(2n+1).

**Testbarkeit:** Systematischer Test aller (n, n±1)-Paare mit gleicher Gewichtsformel; Out-of-sample auf `zeros6(9).npy`.

**Status:** **Widerlegt in jetziger Form** — Rand allein (114,3) schlägt das Paar; Paar-Hypothese schwach.

### H4 — Mod-16-EABC-Sequenz

**Aussage:** Für n = 113 erzeugt [n, n+1, 2n+1, 2n+2] die Restfolge [1,2,3,4] mod 16, und diese Sequenz korreliert mit ε-Optima stärker als Zufall.

**Testbarkeit:** Monte-Carlo über alle n in [100, 5000]: Anteil mit mod-16-Muster [1,2,3,4] vs. ε-Qualität der zugehörigen Fenster.

**Status:** Algebraisch exakt; statistische Korrelation **noch nicht quantifiziert** — vielversprechendster Strukturbefund.

### H5 — Primzahl-Doppelstruktur (n, 2n+1)

**Aussage:** Nur n mit n prim und 2n+1 prim tragen „harte“ ε_EABC-Signatur.

**Testbarkeit:** Klassifikation aller harten Treffer vs. Prim-Doppelpaare; Fisher-exakter Test.

**Status:** 2/3 harte Treffer passen; 3876-Ausreißer schwächt H5.

### H6 — Keine CODATA-Vorhersage (Null-Hypothese)

**Aussage:** Keine a-priori EABC-Kette liefert ε innerhalb u(ε) von CODATA.

**Testbarkeit:** Bereits getestet — ε_final ~31σ, bester Einzelwert (114,3) ~43σ.

**Status:** **Nicht widerlegt** — im Gegenteil gestützt.

### H7 — RMT-Kontrolle

**Aussage:** Die beobachtete ε-Präzision ist von GUE-typischer Gap-Fluktuation nicht zu unterscheiden.

**Testbarkeit:** Erweiterter GUE-Surrogat mit korrekter Odlyzko-Normierung; Vergleich der Verteilung von min |Δε| über Surrogate vs. echte Nullstellen.

**Status:** Vorläufig offen — GUE-Test zeigt ähnliche p (~0,0014), aber Surrogat ist vereinfacht.

### H8 — Skalen-Lokalität

**Aussage:** Die Kalibrierung gilt nur für frühe Nullstellen (center ≈ 113), nicht asymptotisch.

**Testbarkeit:** Block-Test (bereits angedeutet); systematische ε-Drift ε(c) über c.

**Status:** **Gestützt** — Block-Wiederkehr scheitert ab k ≥ 1.

### H9 — σ-Randspannung (zweite EABC-Randkorrektur)

**Aussage:** Die Korrektur ε_final = ε_pair_EABC − σ̄/(16·20·137) mit σ̄ = (σ_core + σ_edge)/2 kompensiert lokale Gap-Varianz zwischen Kern- und Randfenster und ist keine reine Overfitting-Maßnahme.

**Testbarkeit:** Anwendung derselben vierstufigen Kette auf Kontrollfenster (120,5), (105,6), (130,4) und weitere zufällige (c, w)-Paare; Blindtest auf unabhängigen Nullstellen.

**Status:** **Widerlegt als universelle Regel** — nur (113/114) trifft ~10⁻⁷, Kontrollfenster liefern |Δ| ~ 0,01–0,1. Numerisch wirksam für das a-priori-Fenster (14,5× Verbesserung), aber **stark overfitting-anfällig**. Die σ-Formel darf nicht ohne Kontrolltest als Vorhersage behauptet werden.

---

## 8. Empfohlene nächste Tests

| Priorität | Test | Erwarteter Erkenntnisgewinn |
|-----------|------|----------------------------|
| **P1** | Blindtest auf `zeros6(9).npy` (oder unabhängige Nullstellenquelle) | Eingefrorene Regel out-of-sample |
| **P2** | Mod-16-Korrelationsstudie über alle n ∈ [100, N] | Quantifizierung von H4 |
| **P3** | LEE-korrigierte Signifikanz: Westfall-Young über ~10⁶ Fenster | Echte Discovery-Rate für min \|Δε\| |
| **P4** | Vollständiger kombinatorischer Scan aller Regeln aus {16,4,20,137} | Prüfung, ob (113,4) einzigartig ist |
| **P5** | Erweiterter GUE/LUE-Surrogat mit Montgomery-Konvention | RMT-Nullmodell für H7 |
| **P6** | Sensitivität ε_EABC auf alternative Normierungen (240, 256, 272) | Robustheit der /320-Regel |
| **P7** | ε(c)-Heatmap und FDR-Kontrolle über (c,w)-Raum | Globale vs. lokale Optima |
| **P8** | Overfitting-Test: volle Kette auf ≥50 zufällige (c,w)-Paare | Robustheit von H9 |
| **P9** | Expliziter Test: Paar vs. Rand — formale Modellwahl | Klärung H3 |

---

## 9. Anhang Reproduzierbarkeit

### 9.1 Ausführung

```bash
python3 "/Users/thomashoffbauer/Desktop/Mathematische Texte/Alpha_EABC.py"
```

**Laufzeit:** ~55–80 s (Monte-Carlo-Anteile).  
**Versionsmarker (muss erscheinen):**

```
=== SCRIPT-CHECK ===
Version = Alpha_EABC.py :: Zusatztests-v2 :: 2026-06-07
Datei = /Users/thomashoffbauer/Desktop/Mathematische Texte/Alpha_EABC.py
```

### 9.2 Eingabedaten

| Datei | Status | Beschreibung |
|-------|--------|--------------|
| `zeros6.npy` | vorhanden | 2 001 052 Nullstellen, primäre Datenquelle |
| `zeros6(9).npy` | optional, fehlt | Blindtest-Kontrolle (T15/P1) |

### 9.3 Schlüsselparameter (eingefroren)

| Parameter | Wert |
|-----------|------|
| APRIORI_CENTER | 113 |
| APRIORI_WINDOW | 4 |
| MODEL_320 | 320 = 16·20 |
| epsilon_CODATA | 0,035999177 |
| DEFAULT_CENTER / WINDOW | 138 / 20 (Kontrolle) |

### 9.4 Glossar

| Symbol | Bedeutung |
|--------|-----------|
| ε | Feinstruktur-Defekt: 1/α − 137 |
| ε_EABC | Quartilkorrigiertes ε |
| ε_final | Vollständig korrigiertes ε (Median → Quartil → Paar → σ-Randspannung) |
| σ̄ | Mittel der lokalen Gap-Standardabweichungen σ_core und σ_edge |
| EABC | Divisionsalgebren-Programm {R, C, H, O} → 1, 2, 4, 8, 16, … |
| LEE | Look-Elsewhere-Effekt (Multiplizitätstests) |
| harte Signatur | w ∈ {3,4}, lokaler ε_EABC-Rang #1, Kontrast ≥ 100 |
| σ(CODATA) | \|Δε\| / u(ε), u(ε) = 2,1×10⁻⁸ |

### 9.5 Referenzwerte aus Lauf 2026-06-07

| Größe | Wert |
|-------|------|
| ε₀ (113,4) | 0,036181256891 |
| ε_EABC (113,4) | 0,036019306196 |
| ε_EABC edge (114,3) | 0,035998278998 |
| ε_EABC pair | 0,036008746282 |
| ε_final | 0,035999839049 |
| σ_core / σ_edge / σ̄ | 0,384040 / 0,396946 / 0,390493 |
| second_term σ̄/(16·20·137) | 8,907×10⁻⁶ |
| \|Δε_final\| | 6,620×10⁻⁷ |
| p_hard | 2,1108×10⁻⁵ |
| MC p (Gaps) | 0,00145 |
| GUE p | 0,0014 |
| Signatur-Q p | 0,20032 |

---

*Testbericht erstellt am 2026-06-07. Alle Zahlen aus frischem Skriptlauf von `Alpha_EABC.py`.*
