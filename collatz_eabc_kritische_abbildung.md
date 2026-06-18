# EABC kritische Abbildung — Geschwindigkeitsmodell in der komplexen Ebene

**Status:** Modellabbildung (analytisch) — **kein** Physikanspruch  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Modellabbildung | Analogie

**Querverweise:**
- `collatz_eabc_epistemik_physik.md` — **kanonische drei Ebenen:** Wegfunktion, Nicht-Zwillingsparadoxon, Physik-vs.-EABC-Tabelle, zentrale Holonomie-Frage
- `collatz_eabc_zirkulationshypothese.md` — kanonische Zirkulation, $N_\pm$, $D_E$, Holonomie-Vorzeichen $\pm 1$
- `collatz_eabc_zyklus_holonomie.md` — Hierarchie Klasse→Kante→Zyklus→Holonomie
- `collatz_eabc_sagnac_circulation.py` — Kanten $\omega(e)$, Lückenmuster $(2,4,2,4)$
- `collatz_eabc_holonomie_fehlerterm.py` — ABCEA / CEABC-Zählung
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54 / PR #59
- `collatz_eabc_kritische_abbildung.py` — Numerik, JSON-Ausgabe

---

## 0. Epistemische Abgrenzung

$$\boxed{\;\text{Modellabbildung — keine Behauptung über physikalische Geschwindigkeit oder RH-Beweis.}\;}$$

Die folgende Konstruktion ordnet **EABC-Holonomie-Schaltkreise** (ABCEA / CEABC) einer **analytischen Abbildung** entlang der kritischen Linie $\Re s = \tfrac12$ zu. Sie dient als **Lesartshilfe** für Kantenlängen, Orientierung und Zeta-Nullstellen — analog zur Sagnac-Intuition in `collatz_eabc_sagnac.md`, aber mit expliziter komplexer Parametrisierung.

### 0.1 Abgrenzung SRT

Die Notation **„Geschwindigkeit“** / $v$ ist **rein analytisch** — ein **Skalierungsparameter** für die Abbildung $x \mapsto \tfrac12 + \mathrm{i}\,v(x-\tfrac12)$, **keine** physikalische Geschwindigkeit.

$$\boxed{\;\text{Modellabbildung — keine Physikbehauptung (keine SRT).}\;}$$

| Aspekt | SRT (Spezielle Relativitätstheorie) | EABC (dieses Dokument) |
|--------|-------------------------------------|-------------------------|
| Raumzeit | Minkowski-Metrik, Lorentz-Invarianz | Komplexe Ebene $\Re s=\tfrac12$, euklidische $x$-Achse |
| $v$, $c$ | Relativgeschwindigkeit, Lichtgeschwindigkeit | Skalierungsparameter $v>0$; **kein** $c$ |
| Transformation | Lorentz-Transformationen | Lineare Abbildung $s_v(x)=\tfrac12+\mathrm{i}v(x-\tfrac12)$ |
| Kantengrößen $v_j$ | — | $v_j=\gamma_{\mathrm{ref}}/\ell_j$ aus EABC-Primlücken — **Holonomie-Sensor**, keine relativistischen Geschwindigkeiten |
| Sagnac | Interferometrie, Rotation | **Zirkulationstheorie** (`collatz_eabc_zirkulationshypothese.md`); Sagnac nur didaktische Intuition (`collatz_eabc_sagnac.md`) |

**Kernpunkte:**
- Es gibt **keinen** Minkowski-Raum, **keine** Lorentz-Transformationen und **kein** $c$.
- Die abgeleiteten $v_j=\gamma_{\mathrm{ref}}/\ell_j$ messen **Höheninkremente pro EABC-Lücke** (Holonomie-Sensor), nicht Relativgeschwindigkeiten.
- Das frühere Sagnac-Bild ist der **Zirkulationshypothese** untergeordnet (Prime Race $N_\pm$, $D_E$) — kein Relativitäts- oder Quantenmodell.

**Notation (optional):** Eine Umbenennung $v \to \sigma$ oder $\lambda$ als **Skalierungsparameter** wäre terminologisch klarer; im Code und in älteren Abschnitten bleibt vorerst $v$ (Legacy-Label). Eine triviale Alias-Schicht in Python wäre möglich, ein Massen-Rename ist nicht vorgesehen.

---

## 1. Geschwindigkeitsmodell (Benutzerformeln)

**Quelle:** $P = (\tfrac12, 0)$ auf der kritischen Linie (reeller Teil $\tfrac12$, Höhe $\gamma = 0$).

**Geschwindigkeit:** $v > 0$ (reell, skalar).

**Kritische Parallele:** Für eine Markierung $x$ auf der $x$-Achse:
$$\gamma_v(x) = v\,(x - \tfrac12).$$

**Komplexe Abbildung:**
$$s_v(x) = \tfrac12 + \mathrm{i}\,v\,(x - \tfrac12).$$

$$\boxed{\;x \;\longmapsto\; \tfrac12 + \mathrm{i}\,v\,(x - \tfrac12)\;}$$

**Inverse:**
$$\tfrac12 + \mathrm{i}\,\gamma \;\longmapsto\; x = \tfrac12 + \frac{\gamma}{v}.$$

---

## 2. Zeta-Nullstellen

Für eine nicht-triviale Zeta-Nullstelle $\rho_n = \tfrac12 + \mathrm{i}\,\gamma_n$ bei Geschwindigkeit $v$:
$$x_{n,v} = \tfrac12 + \frac{\gamma_n}{v}.$$

Umkehrung:
$$\rho_n = \tfrac12 + \mathrm{i}\,\gamma_n = \tfrac12 + \mathrm{i}\,v\,(x_{n,v} - \tfrac12).$$

**Beispiel** $\gamma_1 \approx 14{,}134725$:

| $v$ | $x_{1,v}$ |
|-----|-----------|
| $1$ | $\approx 14{,}634725$ |
| $2$ | $\approx 7{,}5673625$ |
| $10$ | $\approx 1{,}9134725$ |

**Skalierung:** Große $v$ **komprimiert** Nullstellen auf der $x$-Achse; kleine $v$ **streckt** sie.

**Label:** $s_v$, $x_{n,v}$ = **Definition** (Modellabbildung).

---

## 3. EABC-$C_4$-Zyklus und Kantenlängen

**Elementarer Zyklus** (Vorwärtsorientierung):
$$A \longrightarrow B \longrightarrow C \longrightarrow E \longrightarrow A.$$

**Restklassen** mod $12$: $E\equiv 1$, $A\equiv 5$, $B\equiv 7$, $C\equiv 11$.

**Gerichtete Kanten** mit Primlücken $(r_{\mathrm{dst}} - r_{\mathrm{src}}) \bmod 12$:

| Kante | Lücke $\ell$ | Rolle im Wort ABCEA |
|-------|-------------|---------------------|
| $A \to B$ | $\ell_{AB} = 2$ | 1. Segment |
| $B \to C$ | $\ell_{BC} = 4$ | 2. Segment |
| $C \to E$ | $\ell_{CE} = 2$ | 3. Segment |
| $E \to A$ | $\ell_{EA} = 4$ | 4. Segment |

$$\boxed{\;\text{Lückenmuster auf } C_4\text{: } (2,4,2,4).\;}$$

**Einheitsmodell:** $\ell_{AB}=\ell_{BC}=\ell_{CE}=\ell_{EA}=1$ (normiert; kanonisch sind die mod-$12$-Lücken).

---

## 4. Geschwindigkeiten pro Kante

### 4.1 Freier Parameter (Legacy)

Pro **Startknoten** der Kante kann eine Kantengeschwindigkeit gesetzt werden:

| Kante | Geschwindigkeit | $\Delta\gamma$ auf Segment |
|-------|-----------------|---------------------------|
| $A \to B$ | $v_A$ | $v_A \cdot \ell_{AB}$ |
| $B \to C$ | $v_B$ | $v_B \cdot \ell_{BC}$ |
| $C \to E$ | $v_C$ | $v_C \cdot \ell_{CE}$ |
| $E \to A$ | $v_E$ | $v_E \cdot \ell_{EA}$ |

**Ein-Parameter-Modell:** $v_A = v_B = v_C = v_E = v$.

### 4.2 Abgeleitete Kantengeschwindigkeiten (Holonomie-Sensor)

**Keine freien $v_j$:** Die Geschwindigkeiten werden aus Primlücken mod $12$ und einer Referenzskala $\gamma_{\mathrm{ref}}$ abgeleitet:

$$\boxed{\;v_j = \frac{\gamma_{\mathrm{ref}}}{\ell_j}\;}$$

**Wahl der Referenz:** $\gamma_{\mathrm{ref}} = \gamma_n$ (typisch $n=1$, $\gamma_1 \approx 14{,}134725$) oder ein Skalenparameter `v_base` mit derselben Rolle. Pro Kante gilt dann

$$\Delta\gamma_j = v_j \cdot \ell_j = \gamma_{\mathrm{ref}},$$

d. h. **konstante Höheninkremente** entlang des Schaltkreises — das Ray-Mapping wird zum **EABC-Holonomie-Sensor** (Modellabbildung, kein Physikanspruch).

**Zwei äquivalente Lückenlistungen** (dieselben vier Kanten):

| Listung | Reihenfolge | Kanonisches Muster |
|---------|-------------|-------------------|
| ABCEA-Traversierung | $(\ell_{AB}, \ell_{BC}, \ell_{CE}, \ell_{EA})$ | $(2, 4, 2, 4)$ |
| EAABC-zyklisch | $(\ell_{EA}, \ell_{AB}, \ell_{BC}, \ell_{CE})$ | $(4, 2, 4, 2)$ |

**Beispiel** $\gamma_{\mathrm{ref}} = \gamma_1$, Muster $(2,4,2,4)$:

| Kante | $\ell_j$ | $v_j = \gamma_1 / \ell_j$ |
|-------|---------|---------------------------|
| $E \to A$ | $4$ | $\approx 3{,}533681$ |
| $A \to B$ | $2$ | $\approx 7{,}067363$ |
| $B \to C$ | $4$ | $\approx 3{,}533681$ |
| $C \to E$ | $2$ | $\approx 7{,}067363$ |

Gesamt-$\gamma$ nach einem Umlauf: $4\,\gamma_{\mathrm{ref}} \approx 56{,}539$.

**Label:** Holonomie-Sensor = **Modellabbildung** / **Sensor** (analog zu $C_E$, $\omega(e)$ in `collatz_eabc_sagnac_circulation.py`).

**Kumulative $x$-Koordinate** entlang ABCEA mit Start $x^{(0)} = \tfrac12 = P$:
$$x^{(k)} = x^{(k-1)} + \ell_k,\qquad
s^{(k)} = s_{v_k}\!\bigl(x^{(k)}\bigr) = \tfrac12 + \mathrm{i}\,v_k\,\bigl(x^{(k)} - \tfrac12\bigr),$$
wobei $\ell_k$ und $v_k$ dem $k$-ten Segment in Orientierungsreihenfolge entsprechen.

---

## 5. ABCEA vs. CEABC — Holonomie-Orientierung

| Orientierung | Wort | Vorzeichen | Segmentreihenfolge |
|--------------|------|------------|-------------------|
| $\gamma^+$ | ABCEA | $+1$ | $AB,\,BC,\,CE,\,EA$ |
| $\gamma^-$ | CEABC | $-1$ | $CE,\,EA,\,AB,\,BC$ |

Beide Wörter tragen dasselbe Lückenmuster $(2,4,2,4)$; sie unterscheiden sich durch **zyklische Verschiebung** und **Orientierung** auf $H_1(C_4,\mathbb{Z})$.

**Modellabbildung:** Bei CEABC wird derselbe Zyklus in umgekehrter Kantenreihenfolge durchlaufen; die kumulative $\gamma$-Summe ändert das Vorzeichen der geschlossenen Zirkulation — konsistent mit $\Omega_{\mathrm{Hol}}(\mathrm{ABCEA}) = +1$, $\Omega_{\mathrm{Hol}}(\mathrm{CEABC}) = -1$ in `collatz_eabc_zyklus_holonomie.md`.

$$\boxed{\;\text{Holonomie } \pm 1 \;\leftrightarrow\; Orientierung } \gamma^\pm \text{ auf demselben } C_4\text{-Gerüst.}\;}$$

**Verknüpfung Zirkulationshypothese:** $N_+(X)$ zählt ABCEA-Fenster ($\gamma^+$), $N_-(X)$ zählt CEABC-Fenster ($\gamma^-$); $D_E = N_+ - N_-$ ist der **Fehlerterm** zwischen den beiden orientierten Schaltkreisen (`collatz_eabc_zirkulationshypothese.md` §3–5).

---

## 6. Beispiel ABCEA — Holonomie-Sensor ($\gamma_{\mathrm{ref}}=\gamma_1$)

Abgeleitete Geschwindigkeiten, Lücken $(2,4,2,4)$:

| Stufe | Kante | $\ell$ | $v_j$ | $\Delta\gamma$ | $x$ | $s_v(x)$ |
|------:|-------|-------|-------|----------------|-----|----------|
| 0 | — | — | — | — | $0{,}5$ | $\tfrac12$ |
| 1 | $A\to B$ | $2$ | $\gamma_1/2$ | $\gamma_1$ | $2{,}5$ | $\tfrac12 + \gamma_1\,\mathrm{i}$ |
| 2 | $B\to C$ | $4$ | $\gamma_1/4$ | $\gamma_1$ | $6{,}5$ | $\tfrac12 + 2\gamma_1\,\mathrm{i}$ |
| 3 | $C\to E$ | $2$ | $\gamma_1/2$ | $\gamma_1$ | $8{,}5$ | $\tfrac12 + 3\gamma_1\,\mathrm{i}$ |
| 4 | $E\to A$ | $4$ | $\gamma_1/4$ | $\gamma_1$ | $12{,}5$ | $\tfrac12 + 4\gamma_1\,\mathrm{i}$ |

Gesamtschritt auf $x$-Achse: $\sum \ell = 12$; Gesamt-$\gamma$-Anstieg: $4\,\gamma_1$.

### 6.1 Legacy: einheitliches $v=1$

Start $x^{(0)} = \tfrac12$, alle $v_k = 1$:

| Stufe | Kante | $\ell$ | $x$ | $s_v(x) = \tfrac12 + \mathrm{i}(x-\tfrac12)$ |
|------:|-------|-------|-----|-----------------------------------------------|
| 0 | — | — | $0{,}5$ | $\tfrac12$ |
| 1 | $A\to B$ | $2$ | $2{,}5$ | $\tfrac12 + 2\mathrm{i}$ |
| 2 | $B\to C$ | $4$ | $6{,}5$ | $\tfrac12 + 6\mathrm{i}$ |
| 3 | $C\to E$ | $2$ | $8{,}5$ | $\tfrac12 + 8\mathrm{i}$ |
| 4 | $E\to A$ | $4$ | $12{,}5$ | $\tfrac12 + 12\mathrm{i}$ |

Gesamtschritt auf $x$-Achse: $\sum \ell = 12$; Gesamt-$\gamma$-Anstieg: $12$ bei $v=1$.

---

## 7. Zwei Wege: Halbkreis vs. gerade

**Epistemische Abgrenzung:** Dieser Abschnitt vergleicht **reine euklidische Weglängen** in der komplexen Ebene — **kein** Einstein-Zwillingsparadoxon, **keine** SRT-Zeitdilatation. Einheitsgeschwindigkeit $|dz/dt| = 1$ ⇒ **Rückkehrzeit = Gesamtweglänge**. Kanonische Drei-Ebenen-Lesart: `collatz_eabc_epistemik_physik.md` §1–2 (Wegparadoxon, Fahrrad-Analogie).

**Quelle:** $P = (\tfrac12, 0)$. Besuchspunkte entlang ABCEA (bzw. CEABC) via Holonomie-Sensor $s_v(x_j)$ mit $v_j = \gamma_{\mathrm{ref}}/\ell_j$.

### 7.1 Gerader Polygonzug („Soldaten“-Photon)

Zwischen aufeinanderfolgenden Besuchspunkten $z_k$ geradlinige Segmente, am Ende Rückkehr $z_n \to P$:

$$L_{\mathrm{lin}} = \sum_{k=0}^{n-1} |z_{k+1} - z_k| + |z_0 - z_n|,
\qquad T_{\mathrm{lin}} = L_{\mathrm{lin}}.$$

### 7.2 Verkettete Halbkreisbögen

Zwischen $z_k$ und $z_{k+1}$: Halbkreis in der **oberen Halbebene** mit Durchmesser = Sehne $z_k z_{k+1}$:

$$L_{\mathrm{semi}}(z_k, z_{k+1}) = \frac{\pi}{2}\,|z_{k+1} - z_k|,
\qquad
T_{\mathrm{semi}} = \sum_{k} L_{\mathrm{semi}}(z_k, z_{k+1}) + L_{\mathrm{semi}}(z_n, z_0).$$

$$\boxed{\;T_{\mathrm{semi}} / T_{\mathrm{lin}} = \pi/2 \quad\text{(wenn alle Segmente dieselbe Sehnenlänge haben).}\;}$$

### 7.3 Beispiel ABCEA, $\gamma_{\mathrm{ref}} = \gamma_1$, Lücken $(2,4,2,4)$

Holonomie-Sensor: konstante Höhenschritte $\Delta\gamma_j = \gamma_1$. Besuchspunkte
$z_k = \tfrac12 + \mathrm{i}\,k\,\gamma_1$ für $k = 0,1,2,3,4$.

| Segment | Sehne $|z_{k+1}-z_k|$ |
|---------|----------------------|
| $P \to$ 1. Knoten | $\gamma_1$ |
| 1. $\to$ 2. | $\gamma_1$ |
| 2. $\to$ 3. | $\gamma_1$ |
| 3. $\to$ 4. | $\gamma_1$ |
| 4. $\to P$ | $4\gamma_1$ |

$$T_{\mathrm{lin}} = 8\,\gamma_1 \approx 113{,}078,\qquad
T_{\mathrm{semi}} = 4\pi\,\gamma_1 \approx 177{,}646,\qquad
\frac{T_{\mathrm{semi}}}{T_{\mathrm{lin}}} = \frac{\pi}{2} \approx 1{,}570796.$$

**ABCEA vs. CEABC:** Bei gleichem Lückenmuster $(2,4,2,4)$ und gleichen $v_j$ sind die **Gesamtweglängen** für beide Orientierungen identisch — es gibt **keine** geometrische Asymmetrie analog zum Zwillingsparadoxon; der Holonomie-Kontrast $\pm 1$ betrifft Orientierung/Zirkulation, nicht die euklidische Pfadlänge.

**Label:** Weglängenvergleich = **Modellabbildung** / **Geometrie** (kein Physikanspruch).

---

## 8. Python-Symbolzuordnung

| LaTeX | Python | Modul |
|-------|--------|-------|
| $s_v(x)$ | `s_v` | `collatz_eabc_kritische_abbildung` |
| $x(\gamma,v)$ | `x_from_gamma` | `collatz_eabc_kritische_abbildung` |
| $x_{n,v}$ | `x_n_v` | `collatz_eabc_kritische_abbildung` |
| ABCEA-Schaltkreis | `eabc_circuit_report(..., orientation="ABCEA")` | `collatz_eabc_kritische_abbildung` |
| CEABC-Schaltkreis | `eabc_circuit_report(..., orientation="CEABC")` | `collatz_eabc_kritische_abbildung` |
| $v_j$ aus Lücken | `edge_velocities_from_gaps(gaps, gamma_ref)` | `collatz_eabc_kritische_abbildung` |
| Holonomie-Sensor | `holonomy_sensor_trajectory(orientation, gaps, gamma_ref=...)` | `collatz_eabc_kritische_abbildung` |
| ABCEA vs CEABC | `compare_holonomy_sensor_trajectories(...)` | `collatz_eabc_kritische_abbildung` |
| Prim-Fenster | `prime_window_gap_samples(max_p, limit)` | `collatz_eabc_kritische_abbildung` |
| Gerader Rundweg | `linear_round_trip_time(trajectory_points)` | `collatz_eabc_kritische_abbildung` |
| Halbkreis-Kette | `semicircle_chain_time(trajectory_points)` | `collatz_eabc_kritische_abbildung` |
| Weglängenvergleich | `compare_path_times(orientation, gamma_ref=...)` | `collatz_eabc_kritische_abbildung` |
| Wegfunktion $T=\sum\ell_j^2/\gamma_{\mathrm{ref}}$ | `path_time_T(gaps, gamma_ref)` | `collatz_eabc_kritische_abbildung` |

---

*Modellabbildung: EABC-Holonomie-Schaltkreise als kritische-Linien-Trajektorie mit skalierbarer Geschwindigkeit $v$. Kein Anspruch auf physikalische oder zahlentheoretische Beweiskraft.*
