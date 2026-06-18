# EABC: Wigner-Funktion-Analogie (signiertes Informationsfeld)

**Status:** Modellabbildung + diskrete Numerik  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Analogie | Modellabbildung | Forschungsfrage

**Querverweise:**
- `collatz_eabc_chirale_polarisation.md` — Helizität $\lambda=\pm 1$, $N_R/N_L$, Phasenkanäle
- `collatz_eabc_zirkulationshypothese.md` — kanonisch: $N_\pm$, $C_E$, $D_E$ (5-Block)
- `collatz_eabc_holonomie_stufen.md` — Stufen 1–3, $D_E$ als Zirkulationssensor
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54 / PR #59
- `collatz_eabc_wigner_field.py` — $W_{ab}$, $W_E$, Vorzeichendomänen, Spektral-Stub
- `collatz_eabc_transition_graph.py` — $\Omega_{\mathrm{Pfad}}$, $\Omega_{\mathrm{Hol}}$, Gleitfenster
- `collatz_eabc_graph_laplacian.py` — $\mathrm{Spec}(L_E)$, Near-zero-Moden-Stub
- `collatz_eabc_potential_geometrie.md` — Bohm/AB/Berry, reine Potentialverbindungen

---

## 0. Boxed Leitfrage

$$\boxed{\;\text{Information als signiertes Feld auf } \{E,A,B,C\}\text{ — nicht als Punktzählung.}\;}$$

$$\boxed{\;W_E = \#\mathrm{ABCE} - \#\mathrm{CEAB}\quad\text{(4-Pfad)}\;;\quad D_E = \#\mathrm{ABCEA} - \#\mathrm{CEABC}\quad\text{(5-Zyklus).}\;}$$

**Epistemik:** Wigner-Analogie = **didaktisches Modell** für signierte Quasi-Wahrscheinlichkeit auf dem EABC-Zustandsraum. **Kein** Anspruch auf Quantenphysik in der Primzahlfolge.

---

## 1. Motivation: Wigner-Quasi-Wahrscheinlichkeit als Analogie

In der kontinuierlichen Quantenmechanik beschreibt die **Wigner-Funktion** $W(x,p)$ eine **signierte** Phase-Space-Dichte: positive und negative Regionen kodieren Interferenz, nicht absolute Häufigkeiten.

**Analogie (nicht Identität):** Auf dem diskreten EABC-Gerüst $V=\{E,A,B,C\}$ betrachten wir ein **signiertes Korrelationsfeld** statt bloßer Primzählhäufigkeiten.

| Wigner (Analogie) | EABC-Analog |
|-------------------|-------------|
| $W(x,p)$ signiert | $W(a,b)$ signiert auf $V\times V$ |
| negative „Quasi-Wahrscheinlichkeit“ | CEAB-orientierte 4-Pfade ($\omega=-1$) |
| Interferenzregionen | Vorzeichenwechsel von $W_E(n)$ |
| Zustandsraum $(x,p)$ | Übergangsstruktur $E\!\leftrightarrow\!A\!\leftrightarrow\!B\!\leftrightarrow\!C$ |

**Label:** Wigner-Bild = **Analogie**; EABC-Übertrag = **Modellabbildung**.

---

## 2. Vier-Block $W_E$ vs. Fünf-Block $D_E$

Zentral ist die **Unterscheidung** der Träger (vgl. `collatz_eabc_zyklus_holonomie.md`):

| Observable | Träger | Wörter | Breite | Formel |
|------------|--------|--------|:------:|--------|
| $W_E(N)$ | **Pfad** $P_n^{(4)}$ | ABCE / **CEAB** | 4 | $\#\mathrm{ABCE}-\#\mathrm{CEAB}$ |
| $D_E(N)$ | **Zyklus** $C_n^{(5)}$ | ABCEA / CEABC | 5 | $\#\mathrm{ABCEA}-\#\mathrm{CEABC}$ |

$$\boxed{\;W_E \;\text{misst Pfadorientierung;}\; D_E \;\text{misst geschlossene Holonomie.}\;}$$

**Hinweis (Zirkulationshypothese):** Für 4-Blöcke ist die negative Orientierung **CEAB** (nicht CEABC). Die 5-Block-Holonomie schließt den Zyklus mit fünftem Knoten.

**Normalisierte Varianten:**
$$S_W(N) := \frac{W_E(N)}{N_{\mathrm{ABCE}}+N_{\mathrm{CEAB}}},\qquad
S_E(N) := \frac{D_E(N)}{N_++N_-}.$$

**Label:** $W_E$, $D_E$ = **Definition** (unterschiedliche geometrische Träger).

---

## 3. Diskrete Wigner-Matrix $W(a,b)$

**Vorschlag (dieses Dokument):**

$$W(a,b) = \sum_n \chi_a(n)\,\chi_b(n)\,Q(n),$$

wobei
- $n$ ein Index auf Gleitfenstern der Primklassenfolge oder auf Prim-Vierlingen läuft,
- $\chi_c(n)\in\{0,1,2,3,4\}$ die **Anwesenheit** der Klasse $c$ im 4-Fenster zählt,
- $Q(n)=\Omega_{\mathrm{Pfad}}(P_n^{(4)})\in\{+1,-1,0\}$ den **Vierlings-/Pfadindikator** liefert ($+1$ ABCE, $-1$ CEAB).

**Matrixlayout** über $\{E,A,B,C\}$:

$$
W = \begin{pmatrix}
W_{EE} & W_{EA} & W_{EB} & W_{EC} \\
W_{AE} & W_{AA} & W_{AB} & W_{AC} \\
W_{BE} & W_{BA} & W_{BB} & W_{BC} \\
W_{CE} & W_{CA} & W_{CB} & W_{CC}
\end{pmatrix}.
$$

**Vorzeichenstruktur** = chiral gewichtete Kanalinteraktionen (nicht symmetrisch positiv-definit).

**Zwei Träger in der Numerik:**
1. **Gleitfenster** auf der Primfolge (`sliding_windows`, Breite 4)
2. **Arithmetische Prim-Vierlinge** (`enumerate_quadruplets`)

**Label:** $W_{ab}$ = **Definition** (Modell); Implementierung = **Experiment**.

---

## 4. Vorzeichendomänen und Information in Übergängen

Das kumulative Profil
$$W_E^{\mathrm{cum}}(n) := \sum_{k\le n}\Omega_{\mathrm{Pfad}}(P_k^{(4)})$$
zerlegt die Primfolge in **positive**, **negative** und **Null**-Domänen (Plateaus konstanten Vorzeichens).

**Philosophie:** Die Information liegt in **Vorzeichenwechseln** und **Übergangsstruktur** $E\!\leftrightarrow\!A\!\leftrightarrow\!B\!\leftrightarrow\!C$, nicht in der absoluten Häufigkeit einzelner Primzahlen.

| Domäne | Bedeutung (Modell) |
|--------|-------------------|
| $W_E>0$ | ABCE-dominierte Pfadregion |
| $W_E<0$ | CEAB-dominierte Pfadregion |
| $W_E=0$ | Balance / Null-Linie |
| Vorzeichenwechsel | diskrete „Interferenzkante“ |

**Verknüpfung:** `collatz_eabc_chirale_polarisation.md` — $C_E=D_E$ als Chiralitätsfluss auf 5-Zyklen; hier analog $W_E$ auf 4-Pfaden.

---

## 5. Near-zero-Moden (Spektral-Stub)

**Analogie:** In der Wigner-Darstellung tragen **Interferenzregionen** signifikante Masse nahe Null im Phasenraum.

**EABC-Stub:** Es existiert **kein** Dirac-Operator $D$ im Repo. Stattdessen:

$$L_E^{\mathrm{sym}}\psi_k = \lambda_k \psi_k,\qquad |\lambda_k|\approx 0 \;\Rightarrow\; \text{„POP-Masse“-Analog.}$$

Implementierung: `near_zero_eigenmode_stub` in `collatz_eabc_wigner_field.py`, Daten aus `collatz_eabc_graph_laplacian.py`.

**Label:** Near-zero-Moden = **Analogie** + **Experiment** (Spektral-Stub); kein Dirac-Satz.

---

## 6. Verknüpfung Holonomie-Stufen und Zirkulation

| Stufe | Objekt | Wigner-Lesart |
|:-----:|--------|---------------|
| 1 | $D_E$ Zirkulationsstatistik | globales Vorzeichen des 5-Block-Feldes |
| 2 | Chirale Phasen $\phi_R,\phi_L$ | kontinuierliche Phase zum diskreten $W_{ab}$ |
| 3 | Wilson-Loop $\chi_{\mathrm{Hol}}$ | geschlossene 5-Zyklen vs. offene 4-Pfade |

**Kanonisch:** `collatz_eabc_holonomie_stufen.md`, `collatz_eabc_zirkulationshypothese.md`.

**Strategie:** `collatz_generalangriff_2026.md` — Wigner-Feld als zusätzlicher **Informationskanal** neben $D_E$ und $\mathrm{Spec}(L_E)$.

---

## 7. Abgrenzung

$$\boxed{\;\text{KEINE Quantenphysik-Behauptung — diskretes signiertes Korrelationsmodell.}\;}$$

- Kein Hilbert-Raum, keine echte Wigner-Transformierte einer Dichte-Matrix.
- $W(a,b)$ ist **arithmetische Korrelation** mit Pfadvorzeichen, kein $W(x,p)$.
- $W_E\neq D_E$ im Allgemeinen (4 vs. 5 Block).

---

## 8. Artefakte und Reproduktion

```bash
python3 collatz_eabc_wigner_field.py --max-p 100000
pytest tests/test_eabc_wigner_field.py -q
```

**JSON:** `collatz_eabc_wigner_field.json` — $W$-Matrix, $W_E$, $D_E$, Vorzeichendomänen, Spektral-Stub.
