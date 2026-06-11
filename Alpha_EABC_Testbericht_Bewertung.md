# Bewertung des Testberichts: Alpha_EABC / EABC-Zeta-Kalibrierung

**Gegenstand:** Meta-Evaluation von `Alpha_EABC_Testbericht.md` (Lauf `Alpha_EABC.py :: Zusatztests-v2 :: 2026-06-07`)  
**Prüfung:** Inhaltliche Konsistenz mit dem Skript `Alpha_EABC.py` (Stichprobe der Pipeline, Konstanten, Testmatrix, Hypothesen H1–H9)  
**Datum der Bewertung:** 2026-06-07

---

## 1. Qualität des Berichts

Der Testbericht ist **überdurchschnittlich vollständig** für ein exploratives Numerikprojekt. Er folgt einer klaren Architektur — Executive Summary, Methodik, Ergebnisse, Statistik, Algebra, Hypothesenkatalog, Empfehlungen, Reproduzierbarkeit — und trennt durchgängig zwischen Beobachtung, Interpretation und Urteil. Das ist wissenschaftlich ehrlich und selten bei Projekten dieser Art.

**Vollständigkeit:** Die Testmatrix T1–T15 deckt die im Skript implementierte Pipeline ab. Haupttabelle, Korrekturkette, Robuste Lokalität, Monte-Carlo-, GUE- und Signatur-Tests sowie der Hypothesenkatalog H1–H9 sind dokumentiert. Lücken betreffen weniger fehlende Abschnitte als **fehlende durchgeführte Tests** (Blindtest, systematischer kombinatorischer Scan), die der Bericht selbst korrekt als offen markiert.

**Struktur:** Die vierstufige Korrekturkette (Median → Quartil → Paar → σ-Randspannung) wird methodisch und tabellarisch nachvollziehbar geführt. Die Trennung „gesichert vs. spekulativ“ (Abschnitt 6) und die CAN/CANNOT-Liste für Publikationsansprüche sind vorbildlich.

**Wissenschaftliche Ehrlichkeit:** Der Bericht widersteht der Versuchung, numerische Nähe (~6,6×10⁻⁷) als physikalischen Beweis zu verkaufen. CODATA-Inkompatibilität (~31σ), Overfitting bei H9, LEE-Nichtsignifikanz der harten Signatur und die Tautologie des BPS-Rücktests werden explizit benannt. Das erhöht die Glaubwürdigkeit des Dokuments erheblich.

**Reproduzierbarkeit:** Gut bis sehr gut. Ausführungsbefehl, Versionsmarker, eingefrorene Parameter (`APRIORI_CENTER = 113`, `APRIORI_WINDOW = 4`, `MODEL_320 = 320`), Referenzwerte aus dem Lauf vom 2026-06-07 und Datenquelle `zeros6.npy` sind angegeben. Die Konsistenz mit dem Skript (Normierung Montgomery-Odlyzko, Formeln für ε_EABC, Paar-Gewichtung, σ-Randspannung) wurde stichprobenartig verifiziert — keine Abweichungen festgestellt.

**Schwäche der Berichtsqualität:** An einigen Stellen kollidiert die rhetorische Balance: Die Executive Summary und Ergebnistabellen präsentieren ε_final prominent, während H9 und Abschnitt 3.5 dieselbe Stufe als overfitting-anfällig einstufen. H3 wird in Abschnitt 7 korrekt verworfen, aber in Abschnitt 3.3 und der Haupttabelle bleibt das Paar als „primärer“ Korrekturschritt sichtbar präsent — ein **Redaktions-Inkonsistenz**, kein numerischer Fehler.

**Gesamturteil Qualität:** 8/10 — methodisch reif, ehrlich, gut reproduzierbar; kleine narrative Widersprüche und ausstehende Pflichttests verhindern eine höhere Note.

---

## 2. Stärken

**A-priori vs. post-hoc:** Der Bericht klassifiziert jede Stufe der Kette nach LEE-Status (Abschnitt 4.3). Die Regel c = 137 − 16 − 4 − 4 = 113, w = 4 ist vor dem Datenvergleich festgelegt — das ist der stärkste methodische Anker. Post-hoc-Stufen (Quartilkorrektur, Paar 113/114, σ-Randspannung) werden nicht verschleiert.

**σ-Welt:** Die konsequente Verwendung von u(ε) = 2,1×10⁻⁸ (CODATA 2022) verhindert den typischen Fehler, relative ppm-Nähe mit statistischer Übereinstimmung zu verwechseln. Die Tabelle in Abschnitt 4.1 macht deutlich: selbst der beste Wert (ε_final, ~18 ppm) liegt bei ~31σ — weit außerhalb jeder seriösen Vorhersagebehauptung.

**LEE und Bonferroni:** Der Bericht quantifiziert den Look-Elsewhere-Effekt über ~142 000 Kandidaten der harten Signatur (E[Treffer] ≈ 3,0 bei 3 beobachteten Treffern) und skizziert eine Bonferroni-Korrektur (p_korr ≈ 3,5×10⁻⁷). Damit wird die Multiplizität der Suche ernst genommen — auch wenn die Bonferroni-Anwendung noch nicht über alle parallelen Tests hinweg formalisiert ist.

**Mod-16-Befund:** Die exakte Restfolge [113, 114, 227, 228] ≡ [1, 2, 3, 4] (mod 16) ist **algebraisch exakt und unabhängig vom ε-Fit**. Das ist das robusteste Ergebnis des gesamten Projekts: ein reproduzierbares, vorhersagbares Muster mit klarer Prüfbarkeit (H4).

**Overfitting-Warnung (H9):** Abschnitt 3.5 ist ein Musterbeispiel kritischer Selbstkontrolle. Kontrollfenster (120,5), (105,6), (130,4) liefern |Δε| ~ 0,01–0,1, während nur (113/114) ~10⁻⁷ erreicht. Der Bericht zieht die richtige Konsequenz: numerische Präzision ≠ universelle Vorhersageformel.

**Hypothesenkatalog:** H1–H9 mit Status, Testbarkeit und Urteil erleichtern externe Review und vermeiden nachträgliches Rationalisieren.

---

## 3. Schwächen und Lücken

**Kein Blindtest durchgeführt:** Das Skript enthält `blindtest_other_zero_files` und erwartet optional `zeros6(9).npy`; diese Datei fehlt (Abschnitt 9.2). Der Bericht empfiehlt P1 korrekt, führt den Test aber nicht aus. Für jede Kalibrierungskette mit vier post-hoc-Stufen ist ein Out-of-sample-Test **Pflicht**, nicht optional.

**H3 widerlegt, aber noch im Ergebnisnarrativ:** Abschnitt 7 stellt klar: Rand (114,3) schlägt das gewichtete Paar; H3 ist „widerlegt in jetziger Form“. Gleichzeitig bleibt ε_EABC-Paar in der Haupttabelle und in der Korrekturkette als gleichwertige Stufe präsentiert, ohne durchgängige Kennzeichnung als **überholtes Zwischenmodell**. Empfehlung: Paar explizit als historischen/explorativen Schritt markieren oder aus der „primären“ Kette streichen.

**3876 nicht tiefer analysiert:** Der dritte harte Treffer (3876, 4) wird als schwach benannt (n = 3875 = 5³·31, keine Prim-Doppelstruktur, kein mod-16-Muster [1,2,3,4]). Das reicht für eine erste Einordnung, nicht für eine belastbare H5-Bewertung. Offen bleiben: Warum 3876? Ist es ein Artefakt der harten-Signatur-Definition (lokaler Rang, Kontrast ≥ 100)? Korreliert es mit anderen Strukturmerkmen (Abstand zu 137k, Digitensumme)? Ein kurzer Vergleich mit den beiden 113/114-Treffern fehlt.

**Bonferroni nicht durchgängig:** Erwähnt für die harte Signatur, aber nicht für parallele Tests (MC-Gaps p ≈ 0,0015, GUE p ≈ 0,0014, Q p ≈ 0,20, Kontrollscan Top-5). Bei ~15 Tests und α = 0,05 wäre eine Family-wise-Korrektur oder Westfall-Young (wie in P3 empfohlen) überfällig.

**GUE-Surrogat vereinfacht:** Matrixgröße 160, Bulk-Normierung — kein vollständiges RMT-Nullmodell mit Odlyzko-Konvention. H7 bleibt zu Recht „vorläufig offen“.

**Mehrdimensionale Signatur Q:** p ≈ 0,20 — nicht signifikant. Der Bericht sagt das, aber Q erscheint noch in der Testmatrix und könnte in einer Paper-Fassung gekürzt werden, um den Fokus nicht zu verwässern.

**Kombinatorische Eindeutigkeit von (113,4):** H1 ist nur „teilweise gestützt“ (lokal bestätigt, global nicht bewiesen). Ein vollständiger Scan aller Regeln aus {16, 4, 20, 137} fehlt — der Bericht empfiehlt P4, führt ihn nicht aus.

**Physikalische Brücke:** Keine unabhängige Begründung, warum Riemann-Lücken ε = 1/α − 137 kalibrieren sollten. Das ist im Bericht als spekulativ markiert, bleibt aber die größte konzeptuelle Lücke.

---

## 4. Glaubwürdigkeit der Pipeline — Gesamturteil

| Dimension | Note (1–10) | Kurzbegründung |
|-----------|-------------|----------------|
| **Numerisch** | **7** | Reproduzierbar, kohärente Kette, beeindruckende |Δε|-Verbesserung; aber stark fenster- und stufen-spezifisch |
| **Algebraisch** | **8** | mod-16-Sequenz [1,2,3,4] exakt; Primstruktur (113, 227) stimmig; Sedenion-Deutung abwertend |
| **Statistisch** | **6** | Ehrliche σ- und LEE-Diskussion; fehlender Blindtest, unvollständige Multiplizitätskorrektur, vereinfachtes GUE |
| **Physikalisch** | **3** | Kein kausaler Mechanismus; BPS exakt by construction; CODATA ~31σ entfernt |

**Gesamtpipeline (gewichtet):** **6/10** — als **intern konsistentes, exploratives Kalibrierungslabor** glaubwürdig; als **physikalische Vorhersagepipeline für α** derzeit nicht glaubwürdig.

Die Pipeline verdient Vertrauen in ihrer **Rechenintegrität** und **Selbstkritik**, nicht in ihrer **externen Validität**. Das ist genau die Haltung, die der Bericht einnimmt — und sie ist angemessen.

---

## 5. Bewertung der Hypothesen H1–H9

| Hypothese | Kurzinhalt | Urteil | Begründung |
|-----------|------------|--------|------------|
| **H1** | A-priori-Fenster (113,4) einzigartig in 137-Zone | **Teilweise tragfähig** | Lokal bestes Minimum in [100,130], w=3/4; globaler kombinatorischer Beweis fehlt |
| **H2** | Quartilkorrektur /320 systematisch, kein Overfitting | **Offen / vorsichtig tragfähig** | ~9× numerische Verbesserung; theoretische Begründung schwach; Risiko steigt mit Kettentiefe |
| **H3** | n(n±1)-Paar überlegen | **Verworfen** | Rand (114,3) schlägt Paar; Paar verschlechtert gegenüber Einzel-Rand |
| **H4** | mod-16 [1,2,3,4] korreliert mit ε-Optima | **Algebraisch gesichert, statistisch offen** | Muster für n=113 exakt; MC-Korrelation über alle n noch nicht quantifiziert |
| **H5** | Prim-Doppel (n, 2n+1) für harte Signatur | **Schwach / teilweise** | 2/3 Treffer passen; 3876-Ausreißer untergräbt die Regel |
| **H6** | Keine CODATA-Vorhersage (Null-Hypothese) | **Gestützt / nicht widerlegt** | Alle Stufen ≫ 1σ; ε_final ~31σ |
| **H7** | RMT-GUE nicht unterscheidbar | **Offen** | Ähnliche p (~0,0014), aber Surrogat zu grob |
| **H8** | Nur frühe Nullstellen (center ≈ 113) | **Gestützt** | Block-Wiederkehr scheitert ab k ≥ 1 |
| **H9** | σ-Randspannung universell | **Verworfen als universelle Regel** | Nur 113/114 ~10⁻⁷; Kontrollfenster ~10⁻²–10⁻¹; lokal numerisch wirksam (14,5×) |

**Zusammenfassung:** Tragfähig als **explorative Strukturhypothesen** sind vor allem **H4** (algebraisch) und **H8** (Lokalität). **H6** als Null-Hypothese bleibt intakt — ein wichtiges negatives Ergebnis. Verworfen oder stark geschwächt: **H3**, **H9** (universal), **H5** (strikt). Offen mit Nachholbedarf: **H1**, **H2**, **H4** (statistischer Teil), **H7**.

---

## 6. Paper-Tauglichkeit

### Veröffentlichbar (mit klarer Einordnung)

- **A-priori-Spezifikation** der Fensterregel c = 113, w = 4 aus {16, 4, 20, 137} ohne CODATA-Fit.
- **Numerische Beobachtung** der ε-Nähe auf verschiedenen Korrekturstufen, immer mit σ(CODATA)-Spalte.
- **Algebraischer mod-16-Befund** [113, 114, 227, 228] ≡ [1, 2, 3, 4] und Primstruktur um n = 113.
- **Negative Ergebnisse:** CODATA-Inkompatibilität, LEE-Nichtsignifikanz der harten Signatur, Scheitern der Block-Wiederkehr, Overfitting-Diagnose H9.
- **Methodischer Beitrag:** Transparente Trennung a-priori/post-hoc als Vorbild für ähnliche „Zahlenspiel“-Projekte.

### Exploratory / nicht als Discovery verkaufen

- Die **vierstufige Korrekturkette** bis ε_final — insbesondere Stufen 3 und 4 (Paar, σ-Randspannung).
- **Sedenion-Interpretation** und physikalische Monopol-Ableitungen ohne unabhängigen Test.
- **Mehrdimensionale Signatur Q** (p ≈ 0,20).
- Jede Behauptung der Form „α innerhalb CODATA“, „5σ-Entdeckung“ oder „BPS beweist α“.

### Empfohlener Paper-Typ

Ein **kurzer Methods/Exploration-Note** oder **arXiv-Preprint** im Ton „structured numerical coincidence with honest null results“, nicht ein Physics-Letters-Discovery-Paper. Titel und Abstract müssen H6 und H9 prominent tragen, nicht ε_final.

---

## 7. Empfehlung — priorisierte nächste Schritte

1. **Blindtest (P1):** `zeros6(9).npy` oder unabhängige Nullstellenquelle beschaffen; **eingefrorene** Regel (113,4) + dokumentierte Korrekturstufen out-of-sample testen. Ohne diesen Schritt ist keine Eskalation der Behauptungen vertretbar.

2. **Mod-16-Korrelationsstudie (P2):** Monte-Carlo über n ∈ [100, N]: Anteil mit [1,2,3,4]-Muster vs. ε-Qualität der Fenster. Das stärkste algebraische Asset (H4) statistisch quantifizieren.

3. **LEE / Multiplizität (P3):** Westfall-Young oder vergleichbare Korrektur über den vollen (c,w)-Raum (~10⁶ Paare) für min |Δε| und harte Signatur — Bonferroni im Bericht nachziehen und formalisieren.

4. **H3/H9 aufräumen:** Paar-Stufe aus der „primären“ Erzählung streichen oder als verworfen markieren; σ-Randspannung nur als **explorativer numerischer Befund für (113/114)**, nicht als allgemeine Formel. Narrative Konsistenz zwischen Abschnitt 3 und 7 herstellen.

5. **3876 und H5 vertiefen:** Strukturdiagnose des dritten harten Treffers; ggf. harte-Signatur-Kriterien verschärfen (z. B. mod-16-Pflicht) und erneut scannen — klärt, ob H5 Regel oder Artefakt ist.

---

## Schlusswort

Der Testbericht ist **über dem Durchschnitt explorativer Numerik-Arbeit** — nicht wegen der CODATA-Nähe, sondern wegen der **methodischen Disziplin**. Er dokumentiert ein bemerkenswertes numerisches Zusammenspiel aus a-priori-Fenster, Quartilkorrektur und algebraischer mod-16-Struktur, und er sagt offen, dass dies **kein statistischer Beweis für α-Vorhersage** ist. Die größten Lücken sind **fehlende externe Validierung** (Blindtest), **unvollständige Multiplizitätskontrolle** und **narrative Reste verworfener Hypothesen** (H3, Paar in der Hauptkette).

Als Arbeitsgrundlage für die nächste Iteration ist der Bericht **sehr gut geeignet**. Als Grundlage für eine physikalische Publikation reicht er **nur im exploratory/negative-results-Modus** — und genau so sollte er verstanden werden.

---

*Bewertung erstellt am 2026-06-07. Bezugsdokument: `Alpha_EABC_Testbericht.md`; Skript-Stichprobe: `Alpha_EABC.py`.*
