# EABC Holonomie: drei mathematische Stufen und Wachstumsfälle A/B/C

**Status:** Kanonische Formalisierung — Analogie → echte Holonomie → Wilson-Schleife  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Analogie | Modellabbildung | Hypothese | Experiment | Forschungsfrage

**Querverweise:**
- `collatz_eabc_epistemik_physik.md` — Physik vs. EABC; $\mathrm{Hol}_E$ als arithmetische Kernfrage
- `collatz_eabc_chirale_polarisation.md` — **Stufe-2-Upgrade:** $\Psi=(R,L)^\top$, $U_E$, Helizität $\lambda=\pm 1$
- `collatz_eabc_zirkulationshypothese.md` — kanonisch: $N_\pm$, $C_E$, $D_E$, $S_E$
- `collatz_eabc_fehlerterm_hypothese.md` — Fehlerterm $D_E$, $\widetilde{D}_E$
- `collatz_eabc_evolution_analytik.md` — Evolution Bell$\to$Sagnac$\to C_E\to\mathrm{Spec}(L_E)$; Wachstum in $X$
- `collatz_eabc_zyklus_holonomie.md` — Hierarchie Klasse→Kante→Pfad→Zyklus→Holonomie
- `collatz_eabc_holonomie_beweisversuch.md` — analytischer Beweisversuch $\mathrm{Hol}_E=0$
- `collatz_eabc_D_growth.py` — Wachstumsdiagnostik: Fall A/B/C in $N$, Legacy-Szenarien in $X$
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54 / PR #59
- `collatz_eabc_wigner_analog.md` — 7 Abschnitte epistemisch: $W_{ij}$, $W_E(i,j;N)$, arithmetische Wigner-Negativität
- `collatz_eabc_diskrete_geometrie.md` — **kanonische Geometrie-Synthese** ($S^1$, Fluss, Zentralvermutung)
- `collatz_eabc_uebergangsraum.md` — $C_4\cong S^1$, Hodge, magnetischer Laplace, Flussdichte-Vermutung
- `collatz_eabc_signierte_massstruktur.md` — signierte Maßstruktur auf $G_E$
- `collatz_eabc_potential_geometrie.md` — Bohm/AB/Berry-Analogien, reine Potentialverbindungen (§2–3)

---

## 0. Boxed Leitfrage

$$\boxed{\;\text{Ist } D_E \text{ asymptotisch nur Rauschen — oder ein stabiler Chiralitätsparameter?}\;}$$

$$\boxed{\;\mathrm{Hol}_E = \lim_{N\to\infty}\frac{D_E}{N}\quad\text{mit}\quad N=N_++N_-.\;}$$

**Philosophie:** Primzahlen $\to$ gerichtete Zyklen $\to$ Zirkulation $\to$ Holonomie. **Keine** versteckte Relativität — gerichteter diskreter Fluss auf dem EABC-Zustandsraum.

---

## 1. Stufe 1 — Analogie (erreicht)

**Status:** **erreicht** — didaktisch und kombinatorisch implementiert.

$D_E$ fungiert als **Zirkulations-/Holonomie-Sensor** auf dem $C_4$-Gerüst:

| Objekt | Inhalt |
|--------|--------|
| Orientierung $\gamma^+$ | ABCEA |
| Orientierung $\gamma^-$ | CEABC |
| Observable | $D_E = N_+ - N_-$, $C_E = \oint_\gamma \alpha$ |
| Kombinatorik | Primlückenmuster $(2,4,2,4)$ mod $12$, zyklische Phase $A\equiv 5$ vs. $C\equiv 11$ |

**Epistemisches Label:** Stufe 1 = **Analogie** + **Definition** (Zirkulationsstatistik). $N_+-N_-$ ist **Zirkulationsstatistik**, noch **keine** echte Holonomie im differentialgeometrischen Sinn.

**Artefakte:** `collatz_eabc_zirkulationshypothese.md`, `collatz_eabc_sagnac.md`, `collatz_eabc_kritische_abbildung.md`, `collatz_eabc_holonomie_fehlerterm.py`.

---

## 2. Stufe 2 — Holonomie im eigentlichen Sinn (offen)

**Status:** **offen** — erfordert expliziten Zusammenhang, Paralleltransport und Transportobjekt.

**Klassische Definition:**
$$\mathrm{Hol}(\gamma) = \mathcal{P}\exp\oint_\gamma A.$$

**Was fehlt im EABC-Programm:**

| Baustein | EABC-Kandidat | Status |
|----------|---------------|--------|
| Raum / Bundel | $G_E=(V,E)$, $H_1(C_4)\cong\mathbb{Z}$ | **Definition** |
| Zusammenhang | diskrete 1-Form $\alpha$, Kantenphase $\omega(e)$ | `collatz_eabc_zirkulation_spektral.md` |
| Paralleltransport | Pfadprodukt entlang $\tau_n$ | `collatz_eabc_transport.md` |
| Transportobjekt | Chiralität? Signatur? Defekt? Quaternion-Phase? **Chirale Polarisation $(R,L)$?** | **Forschungsfrage** → Upgrade: `collatz_eabc_chirale_polarisation.md` |

**Upgrade-Pfad (PR #59):** Transportobjekt $\Psi = (R,L)^\top$ auf Faserbündel; nach Schleife $U_E = \mathrm{diag}(e^{i\phi_R}, e^{i\phi_L})$; Observable $\phi_R - \phi_L$ (Berry-Phase-/Wilson-Analog). Diskrete Phase: $\phi \mathrel{+}= \omega(\gamma)\cdot\theta_{\mathrm{edge}}$. Implementierung: `collatz_eabc_chirale_transport.py`, Stubs: `collatz_eabc_potential_phase.py` (`aharonov_bohm_phase`, `berry_phase_difference`).

**Kernabgrenzung:** $N_+-N_-$ misst heute **Zirkulationsstatistik** (Stufe 1/3), nicht $\mathcal{P}\exp\oint A$ auf einem definierten Bundelwert — Stufe 2 schließt diese Lücke **modellhaft** über chirale Phasenkanäle.

**Label:** Stufe 2 = **Forschungsprogramm** mit konkretem Kandidat (chirale Polarisation); vollständiger Beweis = **offen**.

---

## 3. Stufe 3 — Wilson-Loop-Analogie (diskret, erreicht als Statistik)

**Status:** **diskrete Wilson-Analogie erreicht** — kein Zeitmodell.

**Kontinuierliche Analogie:**
$$W(\gamma) = \mathrm{Tr}\bigl(\mathcal{P}\exp\oint_\gamma A\bigr).$$

**EABC-Diskretisierung:**
$$C_E = \sum_{\text{Fenster } \gamma} \omega(\gamma) = N_+ - N_- = D_E.$$

Das ist eine **diskrete Wilson-Schleife** auf dem orientierten $C_4$-Zyklus: Summe der Orientierungsgewichte über geschlossene 5-Fenster — **kein** Zeitmodell, **keine** SRT.

| Wilson (Physik) | EABC (diskret) |
|-----------------|----------------|
| $\mathrm{Tr}(\mathcal{P}\exp\oint A)$ | $C_E = N_+-N_-$ |
| Eichfeld $A$ | diskrete 1-Form $\alpha$ |
| Schleife $\gamma$ | ABCEA / CEABC auf Primfolge |

**Label:** Stufe 3 = **Analogie** (Wilson) + **Definition** ($C_E$ als diskrete Schleifensumme).

**Verknüpfung Fehlerterm:** $D_E$ ist der **absolut**e Fehlerterm der Wilson-Statistik; $S_E=D_E/N$ die **normierte** Holonomie-Observable (`collatz_eabc_fehlerterm_hypothese.md`).

---

## 4. Wachstumsfälle A/B/C (in $N = N_+ + N_-$)

**Zentrale Normalisierung:** $N(X):=N_+(X)+N_-(X)$ (Anzahl geschlossener 5-Fenster bis Primgrenze $X$).

$$\mathrm{Hol}_E := \lim_{X\to\infty} S_E(X) = \lim_{X\to\infty}\frac{D_E(X)}{N(X)}.$$

| Fall | Asymptotik $D_E$ | Grenzverhalten $\mathrm{Hol}_E$ | Interpretation |
|:----:|------------------|--------------------------------|----------------|
| **A** | $D_E = O(1)$ | $\mathrm{Hol}_E = 0$ | absoluter Effekt stirbt; normiert verschwindet |
| **B** | $D_E = O(\sqrt{N})$ (Random Walk / Nullhypothese) | $\mathrm{Hol}_E \to 0$ | typische Fluktuation, keine stabile Chiralität |
| **C** | $D_E \sim \alpha N$, $\alpha\neq 0$ | $\mathrm{Hol}_E = \alpha$ | **bemerkenswerte asymptotische Orientierung** |

**Äquivalente Lesart:**
- Fall A: $|D_E|$ beschränkt, $N\to\infty$ $\Rightarrow$ $S_E\to 0$.
- Fall B: $D_E/\sqrt{N}$ stabil ($\widetilde{D}_E=O(1)$), $S_E\sim 1/\sqrt{N}\to 0$.
- Fall C: $S_E\to\alpha\neq 0$ — permanente ABCEA- bzw. CEABC-Führung im Hauptterm.

$$\boxed{\;\text{Kritisch: Fall C wäre globale Chiralität der Primfolge im EABC-Fluss.}\;}$$

**Abgrenzung zu Legacy-Szenarien in $X$:** `collatz_eabc_evolution_analytik.md` §3 führt zusätzlich Wachstum $D_E(X)$ relativ zu $X$ (A–D: $O(1)$, $O(\log X)$, $O(\sqrt{X})$, Potenzgesetz). Die **holonomierelevante** Klassifikation ist Fall A/B/C **in $N$** (dieses Dokument).

**Experiment:** `collatz_eabc_D_growth.py::classify_holonomy_growth` — RSS-Fit von $D_E$ gegen $N$, $\sqrt{N}$, $\alpha N$; Diagnose $S_E$-Stabilität.

---

## 5. Chebyshev-Verbindung

**Konservative zahlentheoretische Lesart (Hauptterm):** $\mathrm{Hol}_E = 0$, aber $|D_E|$ kann **große Fluktuationen** tragen — analog zum klassischen Chebyshev-Bias, wo $\pi(x;4,3)-\pi(x;4,1)$ oszilliert, während marginalspezifische Raten sich angleichen.

| Aspekt | Chebyshev mod $4$ | EABC $D_E$ |
|--------|-------------------|------------|
| Hauptterm | Gleichverteilung der Restklassen | $N_+\sim N_-$ $\Rightarrow$ $S_E\to 0$ (**Vermutung**) |
| Fehlerterm | strukturierte Oszillation | $D_E$ mit Bias / $L$-Nullstellen mod $12$ (**Hypothese**) |
| Konservativ | $\mathrm{Hol}_E=0$ erwartet | Fall B als Nullhypothese |

**Label:** Chebyshev-Mechanik = **Theorem** (klassisch, bedingt GRH+LI); Übertrag auf $D_E$ = **Hypothese** (`collatz_eabc_evolution_analytik.md` §7).

---

## 6. Philosophie und Epistemik

**Kette (nicht Relativität):**
$$\text{Primzahlen} \;\to\; \text{gerichtete Zyklen} \;\to\; \text{Zirkulation } C_E \;\to\; \text{Holonomie } \mathrm{Hol}_E.$$

- **Nicht** behauptet: versteckte Lorentz-Symmetrie, Zeitdilatation, Zwillingsparadoxon (`collatz_eabc_epistemik_physik.md`).
- **Behauptet:** gerichteter diskreter Fluss auf dem EABC-Zustandsraum $V=\{E,A,B,C\}$.
- **Offen:** ob $D_E$ asymptotisch nur Rauschen (Fall A/B) oder stabiler Chiralitätsparameter $\alpha$ (Fall C).

---

## 7. Verknüpfungen (Zirkulation, Wilson, Fehlerterm)

| Schicht | Objekt | Dokument |
|---------|--------|----------|
| Zirkulation | $C_E$, $\alpha$, $\omega(e)$ | `collatz_eabc_zirkulationshypothese.md` |
| Wilson (diskret) | $C_E=\sum\omega(\gamma)$ | hier §3; `collatz_eabc_zirkulation_spektral.md` |
| Fehlerterm | $D_E$, $\widetilde{D}_E$, $S_E$ | `collatz_eabc_fehlerterm_hypothese.md` |
| Wachstum $N$ | Fall A/B/C | hier §4; `collatz_eabc_D_growth.py` |
| Wachstum $X$ | Legacy A–D | `collatz_eabc_evolution_analytik.md` §3 |
| Beweisversuch | $\mathrm{Hol}_E=0$ | `collatz_eabc_holonomie_beweisversuch.md` |

$$\boxed{\;D_E \;\leftrightarrow\; \text{diskrete Wilson-Schleife} \;\leftrightarrow\; \text{Fehlerterm} \;\leftrightarrow\; \mathrm{Hol}_E.\;}$$

---

## 8. Python-Symbolzuordnung

| LaTeX | Python | Modul |
|-------|--------|-------|
| $N(X)$ | `N_total` | `collatz_eabc_D_growth`, `collatz_eabc_holonomie_fehlerterm` |
| $S_E$, $\mathrm{Hol}_E$ | `S_E` | `collatz_eabc_holonomie_fehlerterm` |
| $\widetilde{D}_E$ | `D_tilde_E` | `collatz_eabc_holonomie_fehlerterm` |
| Fall A/B/C | `holonomy_fall` | `collatz_eabc_D_growth` |
| Legacy $X$-Szenario | `preferred_scenario` | `collatz_eabc_D_growth` |
| $\phi_R$, $\phi_L$, $U_E$ | `phi_R`, `phi_L`, `holonomy_unitary_phases` | `collatz_eabc_chirale_transport` |
| AB/Berry-Stubs | `aharonov_bohm_phase`, `berry_phase_difference` | `collatz_eabc_potential_phase` |
| Helizität ABCEA$\to R$ | `helicity_channel` | `collatz_eabc_chirale_transport` |
| $T_R$, $T_L$ | `travel_time_birefringent` | `collatz_eabc_brachistochrone` |

---

*Drei Stufen: Analogie (erreicht) → echte Holonomie (offen) → Wilson-Loop-Analogie (diskret). Wachstumsdiagnostik Fall A/B/C in $N$: `collatz_eabc_D_growth.py`.*
