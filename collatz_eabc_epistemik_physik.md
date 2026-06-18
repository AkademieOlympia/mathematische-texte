# EABC Epistemik: Physik vs. Modellabbildung (drei Ebenen)

**Status:** Kanonische epistemische Abgrenzung — **kein** Physikanspruch  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Modellabbildung | Analogie | Forschungsfrage

**Querverweise:**
- `collatz_eabc_holonomie_stufen.md` — **drei mathematische Stufen** (Analogie / echte Holonomie / Wilson) + Fall A/B/C in $N$
- `collatz_eabc_kritische_abbildung.md` — §7 Halbkreis vs. gerade; Holonomie-Sensor $v_j=\gamma_{\mathrm{ref}}/\ell_j$
- `collatz_eabc_chirale_polarisation.md` — **Photon-Helizität**, $\phi_R/\phi_L$, Stufe-2-Upgrade
- `collatz_eabc_brachistochrone.md` — $T_R$, $T_L$, Birefringenz-Analogie
- `collatz_eabc_sagnac.md` — Sagnac-Intuition (didaktisch); gegenläufige Zyklen $\gamma^\pm$
- `collatz_eabc_fehlerterm_hypothese.md` — Fehlerterm $D_E$, Prime Race $N_\pm$
- `collatz_eabc_evolution_analytik.md` — Evolution Bell$\to$Sagnac$\to C_E\to\mathrm{Spec}(L_E)$
- `collatz_eabc_zirkulationshypothese.md` — kanonische Zirkulation, $C_E$, $\mathrm{Hol}_E$
- `collatz_eabc_kritische_abbildung.py` — `path_time_T`, `compare_path_times`
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54 / PR #59
- `collatz_eabc_potential_geometrie.md` — Bohm/AB/Berry, $v=f(V)$, $\neq c$

---

## 1. Wegfunktion (ohne Zeitdilatation)

Auf dem EABC-Holonomie-Sensor mit Kantengeschwindigkeiten $v_j = \gamma_{\mathrm{ref}}/\ell_j$ (Primlücken $\ell_j$ mod $12$, Referenzskala $\gamma_{\mathrm{ref}}$) definiert die **Wegfunktion**

$$T = \sum_j \frac{\ell_j}{v_j}
\qquad\text{mit}\qquad
v_j = \frac{\gamma_{\mathrm{ref}}}{\ell_j}
\;\Rightarrow\;
\boxed{\;T = \sum_j \frac{\ell_j^2}{\gamma_{\mathrm{ref}}}\;}$$

**Keine Zeitdilatation:** $T$ ist eine **rein euklidische Laufzeit** bei festgelegten $v_j$ — kein Lorentz-Faktor, kein $c$, kein $ds^2$.

Verschiedene **Geometrien auf denselben Knoten** liefern verschiedene $T$: z. B. gerader Polygonzug vs. verkettete Halbkreisbögen auf dem Holonomie-Sensor ($T_{\mathrm{semi}}/T_{\mathrm{lin}} = \pi/2$ bei gleichen Sehnen, vgl. `collatz_eabc_kritische_abbildung.md` §7).

$$\boxed{\;\text{Wegparadoxon: gleicher Start/Ende/Orientierung, verschiedene Laufzeit — rein euklidisch.}\;}$$

Dasselbe ABCEA- bzw. CEABC-Knotenset, dieselbe Orientierung — aber unterschiedliche Pfadgeometrie ⇒ unterschiedliches $T$. Das ist **kein** Relativitätseffekt, sondern ein **euklidisches Wegparadoxon** in der Modellabbildung.

---

## 2. Nicht das Zwillingsparadoxon

Die **Eigenzeit** in der Speziellen Relativitätstheorie erfordert die Lorentz-Metrik:

$$\tau = \int \sqrt{1 - \frac{v^2}{c^2}}\,\mathrm{d}t.$$

Dafür braucht man $c$, Minkowski-$ds^2$ und relativistische Geschwindigkeiten — **nichts davon ist im EABC-Programm definiert**.

EABC arbeitet mit **Weglängen** $L = \sum_j \ell_j$ auf dem diskreten $C_4$-Gerüst — **ohne** Raumzeit-Metrik, **ohne** $ds^2$.

**Analogie (didaktisch):** Fahrrad auf gerader Strecke vs. Bogen zwischen denselben Endpunkten — gleiche Knoten, verschiedene Weglänge, verschiedene Laufzeit bei Einheitsgeschwindigkeit. Das Fahrrad-Paradoxon ist **geometrisch**, nicht relativistisch.

$$\boxed{\;\text{EABC: kein Zwillingsparadoxon — kein SRT-Anspruch.}\;}$$

---

## 3. Physik vs. EABC (Abgrenzungstabelle)

| Physik (SRT / Interferometrie) | EABC (dieses Repo) |
|--------------------------------|---------------------|
| Zwillingsparadoxon | **nein** |
| Zeitdilatation | **nein** |
| Eigenzeit $\tau$ | **nein** |
| Sagnac | **ja** (didaktisch) |
| Holonomie | **ja** |
| Umlaufsinn | **ja** |
| Zirkulation | **ja** |

**Holonomie / Zirkulation in EABC:** $C_E = N_+ - N_-$ (orientierter Zyklus-Race auf Primfolge-Fenstern). Die Analogie zu Wilson-Schleife / AB-Phase betrifft **Orientierung und Umlaufsinn** — **nicht** die Skalierungsparameter $v_j$ des Holonomie-Sensors (`collatz_eabc_kritische_abbildung.md` §0, §4).

Sagnac in der Physik: Interferenz zweier gegenläufiger Lichtwege in rotierendem Bezugssystem. Sagnac in EABC: **didaktisches Bild** für ABCEA vs. CEABC auf demselben $C_4$-Kreis (`collatz_eabc_sagnac.md`) — Kern ist die Zirkulationshypothese (`collatz_eabc_zirkulationshypothese.md`).

---

## 4. Offene Frage (zentral)

Die **arithmetische Holonomie** (asymptotischer Grenzwert des Prime Race):

$$\mathrm{Hol}_E = \lim_{X\to\infty} \frac{N_+(X) - N_-(X)}{N_+(X) + N_-(X)}.$$

| Grenzfall | Interpretation |
|-----------|----------------|
| $\mathrm{Hol}_E = 0$ | asymptotische Symmetrie der beiden Orientierungen |
| $\mathrm{Hol}_E \neq 0$ | **globale Chiralität** in der Primfolge |

Diese Frage ist **spektakulärer** als jede Zwillings-Analogie: nicht „wer altert langsamer?“, sondern ob die Primfolge einen **bevorzugten Orientierungssinn** im EABC-Fluss trägt.

$$\boxed{\;\text{Sagnac-Frage in EABC: nicht „wer altert langsamer?“, sondern „bevorzugter Umlaufsinn im Primfluss?\"}\;}$$

**Verknüpfungen:**
- Fehlerterm und Prime Race: `collatz_eabc_fehlerterm_hypothese.md`, `collatz_eabc_zirkulationshypothese.md` §5
- Drei Stufen + Wachstumsfälle A/B/C in $N$: `collatz_eabc_holonomie_stufen.md`
- Analytischer Beweisversuch $\mathrm{Hol}_E=0$: `collatz_eabc_holonomie_beweisversuch.md`
- Evolutionspfad und Wachstumsszenarien: `collatz_eabc_evolution_analytik.md`
- Wegparadoxon (euklidisch): `collatz_eabc_kritische_abbildung.md` §7

---

## 5. Chirale Polarisation (Stufe-2-Upgrade)

**Didaktische Physik-Korrektur (kein Behauptungsanspruch):** Photon-Spin $s=1$, Helizität $\lambda=\pm 1$ — nicht spinlos.

**EABC-Lesart:** $\omega(\gamma)=\pm 1$ als **diskrete Helizität**:
- ABCEA $\leftrightarrow$ $\lambda=+1$ ($R$)
- CEABC $\leftrightarrow$ $\lambda=-1$ ($L$)

$$C_E = \sum \omega(\gamma) = N_R - N_L = D_E.$$

**Stufe-2-Transportobjekt** (`collatz_eabc_holonomie_stufen.md` §2):
$$\Psi = (R,L)^\top,\qquad U_E = \mathrm{diag}(e^{i\phi_R}, e^{i\phi_L}),\qquad \text{Observable: }\phi_R - \phi_L.$$

**Brachistochrone-Doppelkanal:** $T_R = \int \mathrm{d}s/v_R$, $T_L = \int \mathrm{d}s/v_L$ — Birefringenz-Analogie (optisch aktive Medien), **nicht** relativistische Eigenzeit.

$$\boxed{\;\text{Chiraler Polarisationsraum — kein Zwillingsparadoxon, kein } \tau = \int\sqrt{1-v^2/c^2}\,\mathrm{d}t.\;}$$

**Artefakte:** `collatz_eabc_chirale_polarisation.md`, `collatz_eabc_chirale_transport.py`, `collatz_eabc_brachistochrone.py`.

**Label:** Helizität/Polarisation = **Modellabbildung**; Photon-Tabelle = **Analogie** (didaktisch).

---

*Drei Ebenen: (1) Wegfunktion $T$ ohne Zeitdilatation, (2) explizite Nicht-SRT-Abgrenzung, (3) Zirkulation/Holonomie als arithmetischer Kern — Physikmetaphern nur didaktisch.*
