# EABC: Wigner-Funktion-Analogie (signiertes Informationsfeld)

**Status:** Modellabbildung + diskrete Numerik  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59), `collatz/eabc-euklidische-hebung`  
**Tao-Labels:** Definition | Analogie | Hypothese | Modellabbildung | offene Frage | Forschungsfrage

**Querverweise:**
- `collatz_eabc_chirale_polarisation.md` — Helizität $\lambda=\pm 1$, $N_R/N_L$, Phasenkanäle
- `collatz_eabc_holonomie_stufen.md` — Stufen 1–3, $D_E$ als Zirkulationssensor, Wilson-Loop
- `collatz_eabc_potential_geometrie.md` — Bohm/AB/Berry, globale Phasenakkumulation $\oint\omega$
- `collatz_eabc_zirkulationshypothese.md` — kanonisch: $N_\pm$, $C_E$, $D_E$ (5-Block)
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54 / PR #59
- `collatz_eabc_wigner_field.py` — $W_{ij}$, $W_E(i,j;N)$, Informationsüberschuss-Test
- `collatz_eabc_transition_graph.py` — $\Omega_{\mathrm{Pfad}}$, $\Omega_{\mathrm{Hol}}$, Gleitfenster

---

## 0. Boxed Leitfrage

$$\boxed{\;\text{Information als signiertes Feld auf Übergängen } E\!\to\!A\!\to\!B\!\to\!C\!\to\!E\text{ — nicht als Punktzählung }p_n.\;}$$

**4-Block vs. 5-Block (strikt trennen):**

| Observable | Träger | Wörter | Breite | Formel |
|------------|--------|--------|:------:|--------|
| $W_E(N)$ | **Pfad** $P_n^{(4)}$ | ABCE / **CEAB** | 4 | $\#\mathrm{ABCE}-\#\mathrm{CEAB}$ |
| $D_E(N)$ | **Zyklus** $C_n^{(5)}$ | ABCEA / CEABC | 5 | $\#\mathrm{ABCEA}-\#\mathrm{CEABC}$ |

$$\boxed{\;W_E \;\text{misst Pfadorientierung;}\; D_E \;\text{misst geschlossene Holonomie.}\;}$$

**Epistemik:** Wigner-Analogie = **didaktisches Modell** für signierte Quasi-Wahrscheinlichkeit auf dem EABC-Übergangsraum. **Kein** Anspruch auf Quantenphysik in der Primzahlfolge.

---

## 1. Wo die Analogie trägt (**Analogie**)

In der kontinuierlichen QM beschreibt die Wigner-Funktion $W(x,p)$ eine **signierte** Phasenraum-Dichte: positive und negative Regionen kodieren Interferenz, nicht absolute Häufigkeiten.

**EABC-Übertragung:** Das signierte Feld auf dem EABC-Gerüst ist

$$W_E(N) = \#\mathrm{ABCE} - \#\mathrm{CEAB}$$

(4-Pfad-Orientierung; **nicht** mit $D_E$ auf 5-Zyklen verwechseln).

| Wigner (kontinuierlich) | EABC (diskret, 4-Block) |
|-------------------------|-------------------------|
| $W(x,p)$ signiert | $W_E(N)$ signiert |
| negative Quasi-Wahrscheinlichkeit | CEAB-orientierte 4-Pfade ($\omega=-1$) |
| Interferenzregionen | **Vorzeichenwechsel** von $W_E(n)$ |
| Information in Phasenstruktur | Information in **Orientierungsbilanz**, nicht in Einzelprimzahlen |

**Label:** Wigner-Bild = **Analogie**; $W_E$ als signiertes Feld = **Modellabbildung**.

**Verknüpfung:** `collatz_eabc_chirale_polarisation.md` — chiraler Fluss; hier auf 4-Pfaden statt 5-Zyklen.

---

## 2. Phasenraum = Übergänge $E\!\to\!A\!\to\!B\!\to\!C\!\to\!E$, nicht Punkte $p_n$ (**Definition**)

Der relevante Zustandsraum ist **nicht** die Menge der Primzahlen $\{p_n\}$, sondern die **gerichtete Übergangsstruktur** auf $V=\{E,A,B,C\}$:

$$\tau_n = (X_n, X_{n+1}), \qquad X_n = \kappa(p_n) \in V.$$

| Falsch (Punktlesart) | Richtig (Übergangslesart) |
|----------------------|---------------------------|
| Häufigkeit von Prim $p_n$ | Kante $\tau_n: X_n \to X_{n+1}$ |
| Einzelklasse $\chi_c(p_n)$ | Pfad $P_n^{(4)}=(X_n,\ldots,X_{n+3})$ |
| Skalar auf $\mathbb{N}$ | Feld auf $V\times V$ (Kantenpaare) |

**Philosophie:** Die Information liegt in der **Topologie der Übergänge** $E\!\leftrightarrow\!A\!\leftrightarrow\!B\!\leftrightarrow\!C\!\leftrightarrow\!E$, analog zum Phasenraum $(x,p)$ — nicht in der absoluten Primzahldichte.

**Label:** Übergangsraum = **Definition**; Punktzählung als Observable = **irrelevant** für Wigner-Analog.

**Artefakt:** `collatz_eabc_transition_graph.py` — `transition_counts`, `sliding_windows`.

---

## 3. Sagnac-Normalisierung: normiertes chirales Feld (**Definition**)

Analog zur Sagnac-Zirkulation (vgl. `collatz_eabc_zirkulationshypothese.md`, `collatz_eabc_sagnac_circulation.py`) definieren wir das **normalisierte** 4-Block-Feld:

$$S_W(N) := \frac{W_E(N)}{N_+^{(4)} + N_-^{(4)}} = \frac{N_+^{(4)} - N_-^{(4)}}{N_+^{(4)} + N_-^{(4)}} \in [-1,1],$$

wobei $N_+^{(4)}=\#\mathrm{ABCE}$, $N_-^{(4)}=\#\mathrm{CEAB}$ auf Gleitfenstern der Breite 4.

| 4-Block | 5-Block |
|---------|---------|
| $S_W = W_E/(N_+^{(4)}+N_-^{(4)})$ | $S_E = D_E/(N_++N_-)$ |
| Pfad ABCE/CEAB | Zyklus ABCEA/CEABC |
| Offene Orientierung | Geschlossene Holonomie |

**Label:** $S_W$ = **Definition** (normalisiertes chirales Quasi-Wahrscheinlichkeitsfeld).

**Verknüpfung:** `collatz_eabc_holonomie_stufen.md` §1 — $S_E$ auf 5-Zyklen; hier $S_W$ auf 4-Pfaden.

---

## 4. Übergangsgeometrie $W_{ij}=\sum_n \chi_i(n)\chi_j(n+1)$ — **nicht** literale Wigner-Matrix (**Modellabbildung**)

Die diskrete **Übergangskorrelationsmatrix** misst die Geometrie der Klassenfolge:

$$W_{ij} = \sum_n \chi_i(n)\,\chi_j(n+1),$$

wobei $\chi_i(n)\in\{0,1\}$ die Anwesenheit von Klasse $i$ an Position $n$ markiert.

$$
W = \begin{pmatrix}
W_{EE} & W_{EA} & W_{EB} & W_{EC} \\
W_{AE} & W_{AA} & W_{AB} & W_{AC} \\
W_{BE} & W_{BA} & W_{BB} & W_{BC} \\
W_{CE} & W_{CA} & W_{CB} & W_{CC}
\end{pmatrix}.
$$

**Abgrenzung:** Dies ist **arithmetische Übergangsgeometrie** (Kantenfrequenz), **keine** echte Wigner-Transformierte $W(x,p)$ einer Dichte-Matrix. Die Analogie dient der **didaktischen** Übertragung von Korrelationsstruktur — nicht der physikalischen Identität.

**Zusätzlich (signierte Fenster-Korrelation):** $W^{\mathrm{sign}}(a,b)=\sum_n \chi_a(n)\chi_b(n)\,Q(n)$ mit $Q(n)=\Omega_{\mathrm{Pfad}}(P_n^{(4)})\in\{+1,-1,0\}$ — kombiniert Übergangsgeometrie mit Pfadorientierung im gleichen 4-Fenster.

**Label:** $W_{ij}$ = **Definition** (Übergangsgeometrie); Wigner-Vergleich = **Analogie**.

**Implementierung:** `build_w_transition_matrix`, `build_w_matrix_from_windows` in `collatz_eabc_wigner_field.py`.

---

## 5. Wigner lokaler Phasenraum vs. EABC globale Zirkulation $C_E=\oint\omega$ (**Analogie**)

| Wigner (lokal) | EABC (global) |
|----------------|---------------|
| $W(x,p)$ punktweise im Phasenraum | $C_E = \oint_\gamma \omega = D_E$ auf geschlossenen 5-Zyklen |
| lokale Interferenzmuster | globale Zirkulationsstatistik $N_+-N_-$ |
| Vorzeichendomänen im $(x,p)$-Raum | Vorzeichendomänen von $W_E^{\mathrm{cum}}(n)$ auf der Primfolge |

**Stufen-Verknüpfung** (`collatz_eabc_holonomie_stufen.md`):

| Stufe | Objekt | Wigner-Lesart |
|:-----:|--------|---------------|
| 1 | $D_E$ Zirkulationsstatistik | globales Vorzeichen des 5-Block-Feldes |
| 2 | Chirale Phasen $\phi_R,\phi_L$ | kontinuierliche Phase zum diskreten $W_{ij}$ |
| 3 | Wilson-Loop $\chi_{\mathrm{Hol}}$ | geschlossene 5-Zyklen vs. offene 4-Pfade |

**Potentialgeometrie** (`collatz_eabc_potential_geometrie.md`): Aharonov–Bohm-Phasenintegral $\phi_{\mathrm{AB}}=\sum_{e\in\gamma}A(e)$ auf 5-Zyklen — globale Zirkulation, nicht lokale Punktstatistik.

**Label:** Lokaler-vs-globaler Kontrast = **Analogie**; $C_E=D_E$ = **Definition**.

---

## 6. KRITISCH: Arithmetische Wigner-Negativität — Information jenseits der Marginalen (**Hypothese**)

**Quantenmechanik:** Negative Wigner-Werte entstehen aus **Nicht-Kommutativität** ($[\hat{x},\hat{p}]\neq 0$); die Vorzeichenstruktur trägt Information, die in den **Marginalen** allein nicht enthalten ist.

**EABC-Analogie:** Negative Beiträge entstehen aus **orientierter Zählung** ABCE vs. CEAB — nicht aus Nicht-Kommutativität im Hilbert-Raum.

$$\boxed{\;\text{Hypothese (arithmetische Wigner-Negativität): Die Vorzeichenstruktur von }W_E(i,j;N)\text{ trägt Information, die nicht aus den Marginalen allein rekonstruierbar ist.}\;}$$

**Zu beweisen / zu testen:**
1. Berechne $W_E(i,j;N)$ pro Kantenpaar (§7).
2. Rekonstruiere aus Marginalen allein (globales $S_W$, Übergangshäufigkeiten $T_{ij}$).
3. Falls $\|W_E - W_E^{\mathrm{marg}}\| > 0$: **Informationsüberschuss** — analog zu QM-Negativität jenseits klassischer Wahrscheinlichkeiten.

**Label:** Arithmetische Wigner-Negativität = **Hypothese**; Numerik = **Experiment** (`information_excess_test`).

**Verknüpfung:** `collatz_eabc_chirale_polarisation.md` — $C_E=D_E$ als Polarisationsoperator; hier die 4-Pfad-Variante auf Kantenpaaren.

---

## 7. Forschungsprogramm: $W_E(i,j;N)$ pro Kantenpaar (**offene Frage**)

Für jedes gerichtete Kantenpaar $(i,j)$ definieren wir die **chirale Quasi-Wahrscheinlichkeit** entlang $E\!\to\!A\!\to\!B\!\to\!C\!\to\!E$ vs. Umkehr:

$$W_E(i,j;N) = \frac{N_{ij}^{(+)} - N_{ij}^{(-)}}{N_{ij}^{(+)} + N_{ij}^{(-)}},$$

wobei
- $N_{ij}^{(+)}$ = Anzahl 4-Fenster ABCE, die die Kante $i\to j$ enthalten,
- $N_{ij}^{(-)}$ = Anzahl 4-Fenster CEAB, die die Kante $i\to j$ enthalten.

**Forschungsfragen:**
1. Welche Kantenpaare tragen maximale Chiralität $|W_E(i,j)|$?
2. Ist $W_E(i,j)$ von $T_{ij}$ (ungewichtete Übergangshäufigkeit) und globalem $S_W$ **unabhängig** bestimmbar?
3. Korreliert $W_E(i,j)$ mit der Berry-Phasendifferenz $\phi_R-\phi_L$ (`collatz_eabc_potential_geometrie.md`)?

**Label:** $W_E(i,j;N)$ = **Definition**; Abhängigkeiten und Grenzwerte = **offene Frage**.

**Implementierung:** `w_e_edge_pair_field`, `information_excess_test` in `collatz_eabc_wigner_field.py`.

---

## 8. Abgrenzung

$$\boxed{\;\text{KEINE Quantenphysik-Behauptung — diskretes signiertes Korrelationsmodell.}\;}$$

- Kein Hilbert-Raum, keine echte Wigner-Transformierte einer Dichte-Matrix.
- $W_{ij}$ ist **arithmetische Übergangsgeometrie**, kein $W(x,p)$.
- $W_E\neq D_E$ im Allgemeinen (4 vs. 5 Block).
- „Arithmetische Wigner-Negativität" = **Hypothese**, kein bewiesener Satz.

---

## 9. Artefakte und Reproduktion

```bash
python3 collatz_eabc_wigner_field.py --max-p 100000
pytest tests/test_eabc_wigner_field.py -q
```

**JSON:** `collatz_eabc_wigner_field.json` — $W_{ij}$, $W_E(i,j;N)$, $W_E$, $D_E$, Informationsüberschuss-Test, Vorzeichendomänen.
