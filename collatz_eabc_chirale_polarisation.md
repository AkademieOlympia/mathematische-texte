# EABC: Chirale Polarisation und Photon-Helizität (didaktische Korrektur)

**Status:** Modellabbildung + Stufe-2-Upgrade-Pfad  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Analogie | Modellabbildung | Forschungsfrage

**Querverweise:**
- `collatz_eabc_holonomie_stufen.md` — **§2 Upgrade-Pfad** (echte Holonomie, Transportobjekt)
- `collatz_eabc_zirkulationshypothese.md` — $C_E = N_+ - N_- = D_E$
- `collatz_eabc_epistemik_physik.md` — kein SRT-Anspruch; chiraler Polarisationsraum
- `collatz_eabc_brachistochrone.md` — $T_R$, $T_L$, Doppelbrechung
- `collatz_eabc_chirale_transport.py` — Numerik $\phi_R$, $\phi_L$, $\mathrm{Hol}$-Phasendifferenz
- `collatz_eabc_kritische_abbildung.md` — Holonomie-Sensor, Kantenlängen $\ell_j$

---

## 0. Boxed Leitfrage

$$\boxed{\;\omega(\gamma)=\pm 1\;\text{als diskrete Helizität}\;\lambda=\pm 1\;\text{— nicht nur Umlaufzahl.}\;}$$

$$\boxed{\;\Psi = (R,L)^\top,\quad U_E = \mathrm{diag}(e^{i\phi_R}, e^{i\phi_L}),\quad \text{Observable: }\phi_R - \phi_L.\;}$$

**Epistemik:** Didaktische Physik-Korrektur — **kein** Behauptungsanspruch über reale Photonen in der Primzahlfolge.

---

## 1. Physik-Korrektur (pedagogisch, nicht Behauptung)

| Teilchen | Spin $s$ | Helizität $\lambda$ (masselos/transversal) |
|----------|:--------:|:------------------------------------------:|
| Higgs | $0$ | — |
| Photon | $1$ | $\lambda = \pm 1$ |
| Graviton | $2$ | $\lambda = \pm 2$ |

**Korrektur:** Das Photon ist **nicht** spinlos. Transversale Photonen tragen Helizität $\lambda = \pm 1$ (rechts-/linkszirkulare Polarisation).

**Label:** Physik-Tabelle = **Analogie** (didaktisch); Übertrag auf EABC = **Modellabbildung**.

---

## 2. EABC-Interpretation: diskrete Helizität

Die Orientierungsgewichte $\omega(\gamma) = \pm 1$ auf dem $C_4$-Gerüst lesen wir als **diskrete Helizität**, nicht nur als topologische Umlaufzahl:

| Orientierung | Wort | Helizitätskanal | $\lambda$ |
|--------------|------|-----------------|:---------:|
| $\gamma^+$ | ABCEA | $R$ (rechts) | $+1$ |
| $\gamma^-$ | CEABC | $L$ (links) | $-1$ |

**Chiralitätsfluss:**
$$C_E = \sum_{\text{Fenster }\gamma} \omega(\gamma) = N_R - N_L = D_E.$$

Das ist ein **diskreter Polarisationsoperator** auf der Primfolge: zählt rechts- vs. linkszirkulare Zyklusbeiträge.

**Label:** $N_R$, $N_L$, $C_E$ = **Definition** (identisch $N_\pm$, $D_E$ in `collatz_eabc_zirkulationshypothese.md`).

---

## 3. Brachistochrone-Erweiterung: Doppelbrechung

Zwei Geschwindigkeitskanäle entlang desselben Pfads:

$$T_R = \int \frac{\mathrm{d}s}{v_R},\qquad T_L = \int \frac{\mathrm{d}s}{v_L}.$$

**Analogie:** optisch aktive Medien, Berry-Phasen-Systeme mit unterschiedlicher Phasengeschwindigkeit für $R$- und $L$-Polarisation.

Diskret (Holonomie-Sensor):
$$T_R = \sum_j \frac{\Delta s_j}{v_{R,j}},\qquad T_L = \sum_j \frac{\Delta s_j}{v_{L,j}}.$$

Observable: $T_R - T_L$ oder $T_R/T_L$ — **Birefringenz-Analogie**, kein Zwillingsparadoxon.

**Implementierung:** `collatz_eabc_brachistochrone.py` — `travel_time_birefringent`, `birefringent_velocity_pair`.

**Label:** Brachistochrone-Doppelkanal = **Modellabbildung**; Birefringenz = **Analogie**.

---

## 4. Stufe 2 — Holonomie-Fortschritt (Upgrade-Pfad)

**Verknüpfung:** `collatz_eabc_holonomie_stufen.md` §2 — Transportobjekt offen.

**Vorschlag (dieses Dokument):**

| Baustein | EABC-Kandidat | Status |
|----------|---------------|--------|
| Zustandsvektor | $\Psi = (R, L)^\top$ auf Faserbündel | **Modell** |
| Paralleltransport | Pfadprodukt mit Kanalzuordnung | `collatz_eabc_chirale_transport.py` |
| Holonomie-Matrix | $U_E = \mathrm{diag}(e^{i\phi_R}, e^{i\phi_L})$ | **Modell** |
| Observable | $\phi_R - \phi_L$ (Berry-Phase-Analog) | **Definition** |
| Diskrete Phase | $\phi \mathrel{+}= \omega(\gamma)\cdot\theta_{\mathrm{edge}}$ | **Modell** |

Nach geschlossener Schleife $\gamma_{\mathrm{loop}}$:
$$U_E = \begin{pmatrix} e^{i\phi_R} & 0 \\ 0 & e^{i\phi_L} \end{pmatrix}.$$

**Näher an:** Wilson-Schleife, Aharonov–Bohm, Sagnac-Interferometer — **nicht** am Zwillingsparadoxon.

**Label:** Stufe-2-Upgrade = **Forschungsprogramm** mit konkretem Transportobjekt (chirale Polarisation).

---

## 5. Abgrenzung: kein relativistisches Eigenzeit-Modell

$$\boxed{\;\text{NICHT relativistische Eigenzeit — chiraler Polarisationsraum.}\;}$$

- Kein $c$, kein Minkowski-$ds^2$, keine Zeitdilatation (`collatz_eabc_epistemik_physik.md`).
- $T_R$, $T_L$ sind **euklidische Wegzeiten** bei festgelegten $v_R$, $v_L$.
- $\phi_R$, $\phi_L$ sind **Phasenakkumulation** auf dem diskreten $C_4$-Gerüst.

---

## 6. Verknüpfung $D_E \leftrightarrow \phi_R - \phi_L$

**Diskretes Modell:**
$$\phi_R - \phi_L \;\approx\; D_E \cdot \theta_{\mathrm{edge}},$$
wobei $\theta_{\mathrm{edge}}$ die Kantenphase aus Lückenlängen $\ell_j$ mod $12$ ist (`theta_edge_from_gaps`).

| Objekt | Zirkulation | Phasenkanal |
|--------|-------------|-------------|
| $N_R$ | $N_+$ (ABCEA) | Akkumulation $\phi_R$ |
| $N_L$ | $N_-$ (CEABC) | Akkumulation $\phi_L$ |
| $D_E$ | $N_R - N_L$ | $\phi_R - \phi_L$ (skaliert) |

**Experiment:** `collatz_eabc_chirale_transport.py::chiral_transport_report`.

---

## 7. Python-Symbolzuordnung

| LaTeX | Python | Modul |
|-------|--------|-------|
| $\lambda=\pm 1$ | `helicity_map` | `collatz_eabc_chirale_transport` |
| ABCEA $\to R$ | `helicity_channel("ABCEA")` | `collatz_eabc_chirale_transport` |
| CEABC $\to L$ | `helicity_channel("CEABC")` | `collatz_eabc_chirale_transport` |
| $\phi_R$, $\phi_L$ | `phi_R`, `phi_L` | `accumulate_phases_along_windows` |
| $\phi_R - \phi_L$ | `holonomy_phase_difference` | `collatz_eabc_chirale_transport` |
| $U_E$ | `holonomy_unitary_phases` | `collatz_eabc_chirale_transport` |
| $D_E = N_R - N_L$ | `chirality_flux_from_counts` | `collatz_eabc_chirale_transport` |
| $T_R$, $T_L$ | `travel_time_birefringent` | `collatz_eabc_brachistochrone` |

---

*Photon-Helizität als didaktische Korrektur; EABC-Transportobjekt für Stufe 2: chirale Polarisation $(R,L)$ mit Phasen-Holonomie $U_E$. Upgrade-Pfad: `collatz_eabc_holonomie_stufen.md` §2.*
