# EABC: Reine Potentialverbindungen (Bohm, Aharonov–Bohm, Berry)

**Status:** Modellabbildung + diskrete Stubs — **kein** Physikanspruch  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Analogie | Modellabbildung | Forschungsfrage

**Querverweise:**
- `collatz_eabc_brachistochrone.md` — $T=\int \mathrm{d}s/v$, fünf Potenzialfamilien, $v=f(V)$
- `collatz_eabc_chirale_polarisation.md` — $\phi_R/\phi_L$, Helizität $\lambda=\pm 1$, Doppelbrechung
- `collatz_eabc_chirale_transport.py` — Phasenakkumulation, `holonomy_phase_difference`
- `collatz_eabc_holonomie_stufen.md` — **§2–3** Stufe-2-Holonomie, Wilson-Loop-Analogie
- `collatz_eabc_epistemik_physik.md` — effektive Geschwindigkeiten $\neq c$, kein SRT
- `collatz_eabc_wigner_analog.md` — $W_{ij}$ Übergangsgeometrie, $W_E(i,j;N)$, arithmetische Wigner-Negativität (§6–7)
- `collatz_eabc_wigner_field.py` — $W_E$, Vorzeichendomänen, Spektral-Stub
- `collatz_eabc_potential_phase.py` — Bohm-/AB-/Berry-Stubs

---

## 0. Boxed Leitfrage

$$\boxed{\;\text{Wirkt reines Potential auf EABC-Bewegung — ohne lokale „Kraft“ und ohne }c\text{?}\;}$$

$$\boxed{\;v = f(V),\quad T = \int \frac{\mathrm{d}s}{v(x)},\quad \phi \mathrel{+}= \omega(\gamma)\,\theta_{\mathrm{edge}},\quad \text{Observable: }\phi_R - \phi_L.\;}$$

**Epistemik:** Bohm-, Aharonov–Bohm- und Berry-Bilder sind **Analogien** (didaktisch). EABC-Übertrag = **Modellabbildung** auf diskretem $C_4$-Gerüst. Effektive Geschwindigkeiten $v_R$, $v_L$ und Kantenphasen $\theta_{\mathrm{edge}}$ sind **Modellparameter** — **nicht** Lichtgeschwindigkeit $c$, **nicht** SRT.

---

## 1. Bohm — Quantenpotential und Führungsgeschwindigkeit (**Analogie**)

**Physik (kontinuierlich):** In der Bohm-de-Broglie-Interpretation beeinflusst das **Quantenpotential** $Q$ die Trajektorie über eine **Führungsgeschwindigkeit** $\mathbf{v}$, die von der Wellenfunktion abhängt und **nicht** mit der „Teilchengeschwindigkeit“ im klassischen Sinn identisch sein muss. Reines Potential kann Bewegung steuern, ohne dass lokal eine klassische Kraft sichtbar ist.

**EABC-Analogie:**

| Physik (Bohm) | EABC (diskret) |
|---------------|----------------|
| Quantenpotential $Q(\mathbf{x})$ | Potenzial $V$ aus Brachistochrone-Familien |
| Führungsgeschwindigkeit $\mathbf{v}\neq$ Teilchengeschw. | Effektive Kantengeschwindigkeit $v_j = f(V_j)$ |
| Potentialgradient $\nabla Q$ | diskreter Gradient $\Delta V$ entlang Kante |
| Kein $c$-Anspruch | $v_j$ Modellparameter (`collatz_eabc_epistemik_physik.md`) |

**Brachistochrone-Kopplung** (`collatz_eabc_brachistochrone.md`):
$$T = \int \frac{\mathrm{d}s}{v(x,\gamma)}, \qquad v = f(V),$$
mit fünf Modellfamilien: $\ln x$, $\zeta$-Summe, $V_E=D_E$, Krümmungshemmung, $-\ln P$.

**Chirale Erweiterung** (`collatz_eabc_chirale_transport.py`, `collatz_eabc_brachistochrone.py`):
$$v_R = v_0 + \alpha\, V_E, \qquad v_L = v_0 - \alpha\, V_E.$$

**Stub:** `bohm_like_velocity(potential_gradient)` — $v = v_0 + \beta\,\mathrm{sign}(\nabla V)\,|\nabla V|$ (Modell, kein QM-Anspruch).

**Label:** Bohm-Bild = **Analogie**; $v=f(V)$ auf EABC-Kanten = **Modellabbildung**.

---

## 2. Aharonov–Bohm — Phasenakkumulation ohne lokale Feldstärke (**Analogie**)

**Physik:** Geladene Quantenobjekte akkumulieren Phase entlang eines Pfades durch das Vektorpotential $\mathbf{A}$:
$$\phi_{\mathrm{AB}} = \oint_\gamma \mathbf{A}\cdot \mathrm{d}\mathbf{l},$$
auch wenn $\mathbf{B}=\nabla\times\mathbf{A}=0$ lokal (feldfreie Region). Die Phase ist **nicht-lokal** im Sinne des geschlossenen Umlaufs — sie hängt vom Pfad und der Topologie ab.

**EABC-Analogie:**

| Physik (AB) | EABC (diskret) |
|-------------|----------------|
| Vektorpotential $\mathbf{A}$ | diskrete 1-Form $\alpha$, Kantenphase $\omega(e)\,\theta_{\mathrm{edge}}$ |
| $\oint \mathbf{A}\cdot\mathrm{d}\mathbf{l}$ | Summe $\sum_{e\in\gamma} A(e)$ auf $C_4$-Kanten |
| Schleife $\gamma$ | ABCEA / CEABC (5-Zyklus) |
| Feldfreie Region | kein lokaler „Kraft“-Term — nur Orientierungsgewicht $\omega(\gamma)=\pm 1$ |

**Diskrete Implementierung** (`collatz_eabc_potential_phase.py`):
$$\phi_{\mathrm{AB}}(\gamma) = \sum_{e\in\gamma} A(e), \qquad A(e) \in \mathbb{R}\;\text{(Stub auf Kanten ``AB'', ``BC'', ``CE'', ``EA'').}$$

Verknüpfung Stufe 2 (`collatz_eabc_holonomie_stufen.md` §2): Phasenakkumulation $\phi \mathrel{+}= \omega(\gamma)\cdot\theta_{\mathrm{edge}}$ entlang Gleitfenster ist **diskrete AB-Analogie** auf dem EABC-Übergangsgraph.

**Label:** AB-Effekt = **Analogie**; Kantensumme auf $G_E$ = **Modellabbildung**.

---

## 3. Berry — Geometrische Phase / Paralleltransport (**Analogie**)

**Physik:** Bei adiabatischem Paralleltransport entlang $\gamma$ in einem Parameter-Hilbertraum akkumuliert ein Zustand eine **geometrische Phase** (Berry-Phase):
$$U(\gamma) = \mathcal{P}\exp\!\left(i\oint_\gamma \mathcal{A}\right), \qquad \gamma_{\mathrm{Berry}} = \arg\langle u|\mathcal{P}u\rangle.$$

Für zwei chiralen Kanäle (rechts-/linkszirkulär): Observable $\phi_R - \phi_L$.

**EABC-Analogie:**

| Physik (Berry) | EABC (diskret) |
|----------------|----------------|
| Paralleltransport $\mathcal{P}$ | Pfadprodukt entlang $\tau_n$ (`collatz_eabc_transport.md`) |
| $U(\gamma)=\mathcal{P}\exp\oint\mathcal{A}$ | $U_E = \mathrm{diag}(e^{i\phi_R}, e^{i\phi_L})$ |
| $\phi_R - \phi_L$ | `holonomy_phase_difference`, `berry_phase_difference` |
| Helizität $\lambda=\pm 1$ | ABCEA$\to R$, CEABC$\to L$ (`collatz_eabc_chirale_polarisation.md`) |

**Stufe 2–3** (`collatz_eabc_holonomie_stufen.md` §2–3): Transportobjekt $\Psi=(R,L)^\top$; nach Schleife unitärer Operator $U_E$; Wilson-Analogie $C_E=\sum\omega(\gamma)=D_E$.

**Stub:** `berry_phase_difference(phi_R, phi_L) = \phi_R - \phi_L` — konsistent mit `collatz_eabc_chirale_transport.py`.

**Label:** Berry-Phase = **Analogie**; $U_E$ auf chiralem Faserbündel = **Modellabbildung** (Stufe 2, offen).

---

## 4. Gesamtkarte: Physik → EABC diskrete Analogie

| Physik-Konzept | EABC diskrete Analogie | Artefakt |
|----------------|------------------------|----------|
| Bohm: Quantenpotential $Q$ | Brachistochrone-Potential $V\in\{\ln x,\zeta,D_E,K,-\ln P\}$ | `collatz_eabc_brachistochrone.py` |
| Bohm: Führungsgeschwindigkeit $\neq$ Teilchengeschw. | $v_j=f(V_j)$, $v_R\neq v_L$ | `bohm_like_velocity`, `birefringent_velocity_pair` |
| AB: $\oint\mathbf{A}\cdot\mathrm{d}\mathbf{l}$ | $\sum_{e\in\gamma}A(e)$, $\phi+=\omega\theta_{\mathrm{edge}}$ | `aharonov_bohm_phase` |
| AB: feldfreie Phase | Orientierung ohne lokale Kraft | $\omega(\gamma)=\pm 1$ auf $C_4$ |
| Berry: $U(\gamma)$ | $U_E=\mathrm{diag}(e^{i\phi_R},e^{i\phi_L})$ | `holonomy_unitary_phases` |
| Berry: $\phi_R-\phi_L$ | Observable Phasendifferenz | `berry_phase_difference` |
| Effektive Geschwindigkeit | Modellparameter $v_j$, **nicht** $c$ | `collatz_eabc_epistemik_physik.md` |
| Brachistochrone | $T=\int\mathrm{d}s/v$ | `travel_time_integral` |
| Chiralität $D_E$ | $C_E=N_R-N_L$, Potenzial $V_E$ | `collatz_eabc_zirkulationshypothese.md` |
| Wigner-Feld | signiertes $W_{ab}$, Informationsphase | `collatz_eabc_wigner_field.py` |

$$\boxed{\;\text{Reines Potential }V \;\Rightarrow\; v=f(V) \;\Rightarrow\; T;\qquad \text{reine Phase }\omega\theta \;\Rightarrow\; \phi_R-\phi_L.\;}$$

---

## 5. Zylinder-Doppelpfad pro Vierling (**geplant**)

**Motivation:** Analog zum AB-Interferenzexperiment mit zwei gegenläufigen Wegen um ein Solenoid soll jeder arithmetische **Prim-Vierling** $(p,p+2,p+6,p+8)$ einen **Zylinder-Doppelpfad** tragen:

| Pfad | Orientierung | Kanal |
|------|--------------|-------|
| äußerer Umlauf | ABCEA | $R$, $\phi_R$ |
| innerer Umlauf | CEABC | $L$, $\phi_L$ |

**Status:** **geplant** — Vierlings-Enumeration existiert (`collatz_eabc_wigner_field.py::enumerate_quadruplets`, `collatz_eabc_holonomie_test.py`); gekoppelte Doppelpfad-Phase $\Delta\phi_{\mathrm{vierling}}$ noch nicht implementiert.

**Forschungsfrage:** Korreliert $\phi_R-\phi_L$ pro Vierling mit $W_E$ (4-Block) und globalem $D_E$ (5-Zyklus)?

---

## 6. Python-Symbolzuordnung

| LaTeX / Konzept | Python | Modul |
|-----------------|--------|-------|
| $v = f(\nabla V)$ (Bohm-Stub) | `bohm_like_velocity` | `collatz_eabc_potential_phase` |
| $\oint\mathbf{A}\cdot\mathrm{d}\mathbf{l}$ | `aharonov_bohm_phase` | `collatz_eabc_potential_phase` |
| $\phi_R - \phi_L$ | `berry_phase_difference` | `collatz_eabc_potential_phase` |
| $U_E$ | `holonomy_unitary_phases` | `collatz_eabc_chirale_transport` |
| $T$, $T_R$, $T_L$ | `travel_time_integral`, `travel_time_birefringent` | `collatz_eabc_brachistochrone` |
| Standard-$A$-Feld | `default_ab_edge_field` | `collatz_eabc_potential_phase` |

---

*Reine Potentialverbindungen: Bohm (Geschwindigkeit aus Potential), AB (Phasenintegral), Berry ($U_E$, $\phi_R-\phi_L$) — durchgängig als **Analogie** gelabelt; effektive $v$ sind Modellparameter, nicht $c$.*
